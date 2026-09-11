"""Recompute review findings in new output directories; never certify scientific closure."""
from __future__ import annotations

import argparse
import json
import platform
import subprocess
import sys
import time
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.linalg import subspace_angles
from sklearn.decomposition import NMF
from sklearn.exceptions import ConvergenceWarning
from sklearn.metrics import roc_auc_score, brier_score_loss

from .clinical_baseline import run_clinical_baseline, _read_cell_composition
from .repair_analysis import run_repair_analysis
from .spatial_control import sha256_file, stable_hash


def write_table(frame, path):
    frame.to_csv(path, sep="\t", index=False)


def bh(p):
    p = np.asarray(p, dtype=float)
    out = np.full(len(p), np.nan)
    valid = np.flatnonzero(np.isfinite(p))
    ordered = valid[np.argsort(p[valid])]
    if len(ordered):
        out[ordered] = np.minimum(1., np.minimum.accumulate(
            (p[ordered] * len(ordered) / np.arange(1, len(ordered)+1))[::-1])[::-1])
    return out


def clinical(output, args):
    result = run_clinical_baseline(output_root=output, response_permutations=199,
                                  bootstrap_draws=1000, seed=args.seed)
    predictions = pd.read_parquet(output / "out_of_fold_predictions.parquet")
    rows = []
    for env, frame in predictions.groupby("environment", sort=True):
        if frame.groupby("patient_row_id").fold.nunique().max() != 1:
            raise ValueError("models used different patient folds")
        wide = frame.pivot(index="patient_row_id", columns="model", values="prediction")
        y = frame.drop_duplicates("patient_row_id").set_index("patient_row_id").response_binary.reindex(wide.index).to_numpy(int)
        folds = frame.drop_duplicates("patient_row_id").set_index("patient_row_id").fold.reindex(wide.index).to_numpy(int)
        wide["null_prevalence"] = np.array([y[folds != f].mean() for f in folds])
        for baseline, model in [("null_prevalence", "abundance_primary"),
                                ("qc", "abundance_primary"),
                                ("composition", "joint_primary"),
                                ("qc", "joint_primary")]:
            a, b = wide[model].to_numpy(), wide[baseline].to_numpy()
            rng = np.random.default_rng(args.seed)
            boot = []
            for _ in range(1000):
                index = rng.integers(0, len(y), len(y))
                if len(np.unique(y[index])) == 2:
                    boot.append(roc_auc_score(y[index], a[index]) - roc_auc_score(y[index], b[index]))
            lo, hi = np.quantile(boot, [.025, .975])
            rows.append(dict(environment=env, model=model, baseline=baseline,
                             n_patients=len(y), auroc_increment=roc_auc_score(y, a)-roc_auc_score(y, b),
                             ci_low=lo, ci_high=hi, brier_improvement=brier_score_loss(y,b)-brier_score_loss(y,a),
                             uncertainty="conditional_on_fixed_oof_predictions",
                             result="EXPLORATORY_POSITIVE" if lo > 0 else "INCONCLUSIVE"))
    write_table(pd.DataFrame(rows), output / "paired_model_increment.tsv")
    table = pd.read_csv(output / "patient_level_table.tsv", sep="\t")
    audit = table.groupby(["analysis_dataset", "treatment_values", "endpoint_values"], dropna=False).agg(
        n_patients=("patient_row_id", "nunique"), n_responder=("response_binary", "sum")).reset_index()
    write_table(audit, output / "environment_audit.tsv")
    result["stage_closure"] = "OPEN_EXTERNAL_VALIDATION_AND_SPATIAL_BRIDGE"
    return result


def repair(output, args):
    result = run_repair_analysis(output_root=output, permutations=199, seed=args.seed)
    paired = pd.read_csv(output / "paired_displacements.tsv", sep="\t")
    composition = _read_cell_composition(
        "results/v7/ontology/scores/expression_unit_scores.parquet", set(paired.patient_key))
    keys = ["analysis_dataset", "patient_key"]
    ccols = [c for c in composition if c.startswith("composition__")]
    composition_rows, effects = [], []
    for (env, target), frame in paired.groupby(["environment", "target_timepoint"], sort=True):
        comp = composition[composition.analysis_dataset.eq(env)]
        c = comp[comp.timepoint.eq("T0")].merge(comp[comp.timepoint.eq(target)], on=keys, suffixes=("_pre", "_post"), validate="one_to_one")
        c = c[c.patient_key.isin(frame.patient_key)].set_index("patient_key")
        cdelta = pd.DataFrame({col:c[col+"_post"]-c[col+"_pre"] for col in ccols})
        labels = frame.drop_duplicates("patient_key").set_index("patient_key").response_binary
        for col in ccols:
            for patient in c.index:
                composition_rows.append(dict(environment=env, target_timepoint=target,
                                             patient_key=patient, feature_id=col,
                                             delta=cdelta.loc[patient,col], response_binary=labels.loc[patient]))
        # Raw displacement contrasts plus composition-adjusted contrasts.
        molecular = frame.pivot(index="patient_key", columns="feature_id", values="delta")
        values = pd.concat([molecular, cdelta.reindex(molecular.index)], axis=1)
        for feature in values:
            y = values[feature]
            label = labels.reindex(y.index)
            valid = y.notna() & label.notna()
            yy, ll = y[valid].to_numpy(float), label[valid].to_numpy(int)
            if min((ll==1).sum(), (ll==0).sum()) < 3:
                continue
            rng = np.random.default_rng(args.seed)
            r, n = yy[ll==1], yy[ll==0]
            boot = np.median(r[rng.integers(0,len(r),(1000,len(r)))],axis=1) - np.median(n[rng.integers(0,len(n),(1000,len(n)))],axis=1)
            adjustment, coef, rank, n_adjusted = "NOT_APPLICABLE_CELL_FRACTION", np.nan, 0, 0
            if feature in molecular:
                x = cdelta.reindex(y.index)
                good = valid & x.notna().all(axis=1)
                design = np.column_stack([np.ones(good.sum()), label[good], x[good]])
                # Remove constant composition columns and one redundant simplex
                # coordinate through a rank check; never replace missing fractions.
                nonconstant = np.std(design[:,2:],axis=0)>1e-10 if good.sum() else np.zeros(len(ccols),bool)
                design = np.column_stack([design[:,:2], design[:,2:][:,nonconstant]])
                if design.shape[1] > 2:
                    design = design[:,:-1]  # compositional sum-to-zero constraint
                rank = np.linalg.matrix_rank(design) if len(design) else 0
                n_adjusted = len(design)
                if rank == design.shape[1] and n_adjusted >= rank + 5:
                    coef = np.linalg.lstsq(design,y[good],rcond=None)[0][1]
                    adjustment = "EXPLORATORY_COMPOSITION_ADJUSTED_OLS"
                else:
                    adjustment = "NOT_TESTABLE_RANK_OR_RESIDUAL_DF"
            lo, hi = np.quantile(boot,[.025,.975])
            effects.append(dict(environment=env,target_timepoint=target,feature_id=feature,
                                n_responder=len(r),n_non_responder=len(n),median_contrast=np.median(r)-np.median(n),
                                ci_low=lo,ci_high=hi,adjusted_response_coefficient=coef,
                                adjustment_status=adjustment,adjustment_n=n_adjusted,design_rank=rank,
                                repair_direction="NOT_ASSIGNED_NO_INDEPENDENT_DIRECTION_VALIDATION"))
    write_table(pd.DataFrame(composition_rows), output / "cell_fraction_displacements.tsv")
    write_table(pd.DataFrame(effects), output / "displacement_uncertainty_and_composition.tsv")
    result["stage_closure"] = "OPEN_DIRECTION_EXTERNAL_VALIDATION_SPATIAL_BRIDGE"
    result["pairing_null"] = "whole_patient_post_vector_within_environment_timepoint"
    return result


def spatial_factors(output, args):
    root = Path("results/v7/spatial_discovery")
    frames, audits = [], []
    for path in sorted((root / "scores").glob("*.parquet")):
        frame = pd.read_parquet(path)
        columns = [c for c in frame if c.endswith("__score_standardized")]
        finite = np.isfinite(frame[columns].to_numpy(float)).all(axis=1)
        audits.append(dict(unit_id=frame.unit_id.iloc[0],dataset_id=frame.dataset_id.iloc[0],
                           patient_id=frame.patient_id.iloc[0],n_observations=len(frame),
                           n_complete_observations=int(finite.sum()),
                           status="ELIGIBLE" if finite.sum()>=20 else "NOT_TESTABLE_INCOMPLETE_PANEL"))
        if finite.sum() >= 20:
            frames.append(frame.loc[finite,["patient_id","dataset_id",*columns]])
    write_table(pd.DataFrame(audits), output / "factor_input_eligibility.tsv")
    full = pd.concat(frames,ignore_index=True)
    rows, stability = [], []
    for held_out in sorted(full.dataset_id.unique()):
        train, test = full[full.dataset_id.ne(held_out)], full[full.dataset_id.eq(held_out)]
        # Balanced fitting cap, fixed before evaluation; all complete test rows retained.
        sampled = [f.sample(min(len(f),1500),random_state=args.seed) for _,f in train.groupby("patient_id",sort=True)]
        xtrain = pd.concat(sampled)[columns].to_numpy(float)
        offset = xtrain.min(axis=0)
        xtrain = xtrain-offset+1e-6
        mean = xtrain.mean(axis=0)
        components = {}
        for k in [0,2,4,8]:
            for seed in ([args.seed] if k==0 else [args.seed,args.seed+1,args.seed+2]):
                model, converged, fit_seconds = None, True, 0.
                if k:
                    model = NMF(n_components=k,init="random",random_state=seed,max_iter=1000,tol=1e-4)
                    start = time.monotonic()
                    with warnings.catch_warnings(record=True) as caught:
                        warnings.simplefilter("always",ConvergenceWarning)
                        model.fit(xtrain)
                    converged = not any(issubclass(w.category,ConvergenceWarning) for w in caught)
                    fit_seconds = time.monotonic()-start
                    components[(k,seed)] = model.components_.copy()
                for patient, f in test.groupby("patient_id",sort=True):
                    raw = f[columns].to_numpy(float)-offset+1e-6
                    clipped_fraction = float((raw<0).mean())
                    x = np.maximum(raw,0)
                    transform_converged = True
                    if model is None:
                        prediction = np.broadcast_to(mean,x.shape)
                    else:
                        with warnings.catch_warnings(record=True) as caught:
                            warnings.simplefilter("always",ConvergenceWarning)
                            prediction = model.transform(x) @ model.components_
                        transform_converged = not any(issubclass(w.category,ConvergenceWarning) for w in caught)
                    rows.append(dict(held_out_dataset=held_out,patient_id=patient,k=k,seed=seed,
                                     n_train_patients=train.patient_id.nunique(),n_train_rows=len(xtrain),
                                     n_test_rows=len(x),rmse=float(np.sqrt(np.mean((x-prediction)**2))),
                                     below_training_min_fraction=clipped_fraction,
                                     fit_converged=converged,transform_converged=transform_converged,
                                     n_iter=0 if model is None else model.n_iter_,fit_seconds=fit_seconds,
                                     claim_scope="held_out_module_compressibility_not_open_spatial_field"))
                print(f"factor {held_out} k={k} seed={seed} fit_converged={converged}",flush=True)
        for k in [2,4,8]:
            for a,b in [(args.seed,args.seed+1),(args.seed,args.seed+2),(args.seed+1,args.seed+2)]:
                angle = subspace_angles(components[(k,a)].T,components[(k,b)].T)
                stability.append(dict(held_out_dataset=held_out,k=k,seed_a=a,seed_b=b,
                                      max_angle_degrees=float(np.degrees(angle).max()),
                                      mean_cosine_squared=float(np.mean(np.cos(angle)**2))))
    write_table(pd.DataFrame(rows),output / "patient_held_out_nmf.tsv")
    write_table(pd.DataFrame(stability),output / "nmf_subspace_stability.tsv")
    primitive = pd.read_csv(root / "niche_primitives.tsv",sep="\t")
    patient = primitive.groupby(["dataset_id","patient_id","primitive_id"],as_index=False).estimate.median()
    write_table(patient,output / "patient_balanced_primitives.tsv")
    leave = []
    for primitive_id, f in patient.groupby("primitive_id",sort=True):
        for dataset in sorted(f.dataset_id.unique()):
            train, test = f[f.dataset_id.ne(dataset)], f[f.dataset_id.eq(dataset)]
            estimate, reference = test.estimate.median(), train.estimate.median()
            leave.append(dict(primitive_id=primitive_id,held_out_dataset=dataset,n_test_patients=len(test),
                              test_median=estimate,train_median=reference,sign_agreement=np.sign(estimate)==np.sign(reference),
                              role="development_dataset_transfer_not_independent_validation"))
    write_table(pd.DataFrame(leave),output / "held_out_primitive_results.tsv")
    return dict(stage_closure="OPEN_TOPOLOGY_GT_AND_OPEN_FIELD",n_input_units=len(audits),
                fitting_cap_per_patient=1500,nmf_max_iter=1000,seeds=3,
                evaluation="full_complete_test_rows_with_patient_outer_unit")


def spatial_permutation(output, args):
    from .spatial_discovery import build_stage4_selection, _load_unit
    from .spatial_stats.graph import build_section_graph
    from .moran_matrix import matrix_moran
    root = Path("results/v7/spatial_discovery")
    selection = build_stage4_selection()
    sanitized = output / "sanitized"
    sanitized.mkdir()
    rows = []
    completed = set()
    if args.resume_permutation:
        previous = Path(args.resume_permutation)
        receipt = json.loads((previous / "RUN_RECEIPT.json").read_text())
        if receipt["mode"] != "spatial_permutation" or receipt["seed"] != args.seed:
            raise ValueError("resume mode/seed mismatch")
        for path, expected in receipt["input_hashes"].items():
            if sha256_file(path) != expected:
                raise ValueError(f"resume input changed: {path}")
        prior = pd.read_csv(previous / "moran_progress.tsv",sep="\t")
        for unit_id, frame in prior.groupby("unit_id"):
            if len(frame)!=39 or frame.feature_id.nunique()!=39:
                raise ValueError("partial capture cannot be silently resumed")
            if not frame.loc[frame.status.eq("ESTIMABLE"),"n_permutations"].eq(999).all():
                raise ValueError("resume permutation count differs")
            completed.add(unit_id)
        rows = prior.to_dict("records")
        for row in rows:
            row.setdefault("null_kernel", "independent_feature_permutation_previous_run")
        print(f"reusing {len(completed)} completed captures after input hash checks",flush=True)
    for ordinal,resolved in enumerate(selection.units,1):
        if resolved.unit_id in completed:
            continue
        unit, _ = _load_unit(resolved,sanitized)
        scores_path = root / "scores" / f"{resolved.unit_id}.parquet"
        scores = pd.read_parquet(scores_path)
        if scores.observation_id.astype(str).tolist() != unit.observations.observation_id.astype(str).tolist():
            raise ValueError("score/coordinate observation identity mismatch")
        sections = unit.observations.section_id.to_numpy()
        graph = build_section_graph(unit.observations[["analysis_x","analysis_y"]].to_numpy(float),
                                    sections,method="knn",k=6,weight_mode="binary")
        columns = [c for c in scores if c.endswith("__score_standardized")]
        statistics = matrix_moran(scores[columns].to_numpy(float), graph, sections, permutations=999,
                                  seed=args.seed ^ int(stable_hash({"unit":resolved.unit_id})[:8],16))
        for col, stat in zip(columns, statistics, strict=True):
            feature = col.removesuffix("__score_standardized")
            rows.append(dict(unit_id=resolved.unit_id,dataset_id=resolved.dataset_id,
                             patient_id=resolved.opaque_patient_id,feature_id=feature,
                             observed=stat["observed"],p_value=stat["p_value"],
                             n_permutations=stat["n_permutations"],status=stat["status"],
                             null_kernel="shared_mask_vector_permutation_same_marginal_null"))
        write_table(pd.DataFrame(rows),output / "moran_progress.tsv")
        print(f"permutation capture {ordinal}/{len(selection.units)} {resolved.unit_id}",flush=True)
        del unit, scores, graph
    table = pd.DataFrame(rows)
    table["q_within_capture_39_features"] = table.groupby("unit_id").p_value.transform(lambda s:bh(s))
    table["claim_scope"] = "within_capture_spatial_autocorrelation_not_patient_replication_or_topology_increment"
    write_table(table,output / "moran_999_permutations.tsv")
    return dict(stage_closure="OPEN_TOPOLOGY_INCREMENT_GT",n_records=len(rows),permutations=999,
                multiple_testing_family="39_features_per_capture",patient_inference="NOT_RUN")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("mode",choices=["clinical","repair","spatial_factors","spatial_permutation"])
    parser.add_argument("--output",required=True)
    parser.add_argument("--seed",type=int,default=20260911)
    parser.add_argument("--resume-permutation")
    args = parser.parse_args()
    output = Path(args.output)
    output.mkdir(parents=True,exist_ok=False)
    started = time.time()
    inputs = [Path("results/v7/registry/treatment_response_completeness.tsv"),
              Path("results/v7/ontology/measurement_reliability.tsv")]
    if args.mode in ["clinical","repair"]:
        inputs += [Path("results/v7/ontology/scores/patient_timepoint_scores.parquet"),
                   Path("results/v7/ontology/scores/expression_unit_scores.parquet")]
    else:
        inputs += sorted(Path("results/v7/spatial_discovery/scores").glob("*.parquet"))
        inputs += [Path("results/v7/spatial_discovery/model_safe_input_manifest.tsv"),
                   Path("results/v7/spatial_discovery/niche_primitives.tsv")]
    receipt = dict(status="RUNNING",mode=args.mode,argv=sys.argv,seed=args.seed,started_unix=started,
                   python=sys.version,platform=platform.platform(),
                   git_head=subprocess.check_output(["git","rev-parse","HEAD"],text=True).strip(),
                   input_hashes={str(p):sha256_file(p) for p in inputs},
                   code_hashes={str(p):sha256_file(p) for p in sorted(Path("src/v7").rglob("*.py"))})
    (output / "RUN_RECEIPT.json").write_text(json.dumps(receipt,indent=2)+"\n")
    try:
        result = globals()[args.mode](output,args)
    except Exception as exc:
        receipt.update(status="FAILED",error=f"{type(exc).__name__}: {exc}",elapsed_seconds=time.time()-started)
        (output / "RUN_RECEIPT.json").write_text(json.dumps(receipt,indent=2)+"\n")
        raise
    receipt.update(status="COMPUTATION_COMPLETE",result=result,elapsed_seconds=time.time()-started,
                   outputs={str(p):sha256_file(p) for p in sorted(output.iterdir()) if p.is_file() and p.name!="RUN_RECEIPT.json"})
    (output / "RUN_RECEIPT.json").write_text(json.dumps(receipt,indent=2)+"\n")
    print(json.dumps(receipt["result"],ensure_ascii=False),flush=True)


if __name__ == "__main__":
    main()
