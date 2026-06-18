#!/usr/bin/env python3
"""Phase7A: frozen-topic semantic and measurement validity."""

from __future__ import annotations

import argparse
import hashlib
import os
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

from phase7_common import (
    MODULES,
    ROOT,
    bh_fdr,
    build_design,
    ensure_dirs,
    file_record,
    forbidden_columns,
    hypergeom_sf,
    load_config,
    normalize_timepoint,
    parse_gmt,
    precision_weight,
    rank_stability,
    robust_z,
    sha256,
    spearman_pair,
    weighted_group_mean,
    weighted_residual,
    write_yaml,
)


WHITELIST_SAMPLE = ["sample_key", "cohort_id", "patient_key", "tissue_source", "lesion_context", "timepoint", "treatment_context"]
WHITELIST_COHORT = ["cohort_id", "cancer_type", "platform"]
WHITELIST_QC = [
    "sample_key", "total_cells", "coarse_label_coverage", "mid_label_coverage", "fine_label_coverage",
    "low_quality_fraction", "detected_genes", "umi", "mitochondrial_fraction", "ribosomal_fraction",
    "tissue_source", "treatment_context", "expression_layer", "integration_status",
]


def read_phase6(cfg: dict) -> tuple[dict, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    handoff = yaml.safe_load(cfg["inputs"]["phase6_handoff"].read_text())
    if handoff.get("verdict") != "GO_TO_PHASE7_WITH_FROZEN_MODULE_ENTRY":
        raise RuntimeError(f"Phase6 handoff not open: {handoff.get('verdict')}")
    score_path = Path(handoff["default_module_score_matrix"])
    member_path = Path(handoff["default_module_membership"])
    dictionary_path = Path(handoff["module_dictionary"])
    return handoff, pd.read_parquet(score_path), pd.read_csv(member_path), pd.read_csv(dictionary_path)


def preflight(cfg: dict, handoff: dict, score: pd.DataFrame) -> None:
    out = cfg["out"]
    required = {"phase6_handoff": cfg["inputs"]["phase6_handoff"], **cfg["inputs"]}
    rows = []
    for role, path in required.items():
        rec = file_record(path)
        rec.update(role=role, required=role not in {"tf_regulon"}, status="PASS" if rec["exists"] else "MISSING")
        rows.append(rec)
    pd.DataFrame(rows).to_csv(out / "preflight/phase7_input_audit.csv", index=False)

    key_rows = []
    for key in ("expression_unit_id", "sample_key", "patient_key"):
        present = key in score
        key_rows.append({"key": key, "present": present, "n_missing": int(score[key].isna().sum()) if present else len(score),
                         "n_unique": int(score[key].nunique()) if present else 0})
    pd.DataFrame(key_rows).to_csv(out / "preflight/phase7_key_join_audit.csv", index=False)

    leakage = []
    for source, cols, permitted in [
        ("frozen_module_score", score.columns, False),
        ("sample_metadata_source", pd.read_csv(cfg["inputs"]["sample_metadata"], nrows=1).columns, True),
    ]:
        bad = forbidden_columns(cols, cfg["prohibited_tokens"])
        leakage.append({"source": source, "prohibited_columns": "|".join(bad), "n_prohibited": len(bad),
                        "read_into_model": not permitted, "status": "FAIL" if bad and not permitted else "PASS_DROPPED_AT_SOURCE" if bad else "PASS"})
    pd.DataFrame(leakage).to_csv(out / "preflight/phase7_leakage_preflight.csv", index=False)

    paths = [Path(handoff[k]) for k in ("default_module_score_matrix", "default_module_membership", "module_dictionary", "module_method_support_map")]
    write_yaml({
        "phase": "phase6_closeout_v6_2_1",
        "immutable": True,
        "algorithm_selection_reopened": False,
        "primary_run": handoff["primary_backbone_run"],
        "n_frozen_topics": len(MODULES),
        "assets": [file_record(p) for p in paths],
        "phase7_config": file_record(Path(os.environ.get(
            "PHASE7_CONFIG_PATH", ROOT / "scripts/v6_2/phase7_v6_2_1_config.yaml"))),
        "phase7a_code": file_record(ROOT / "scripts/v6_2/run_phase7a_measurement.py"),
        "phase7_allowed_entry": "frozen_assets_only",
    }, out / "preflight/phase6_closeout_manifest_v6_2_1.yaml")


def prepare_units(cfg: dict, score: pd.DataFrame) -> pd.DataFrame:
    registry = pd.read_csv(cfg["inputs"]["phase6_unit_registry"])
    keep_reg = [c for c in ("expression_unit_id", "n_genes", "n_cells_used", "library_size", "inclusion_status") if c in registry]
    shared = [c for c in ("n_cells_used", "library_size") if c in score and c in registry]
    reg = registry[keep_reg].rename(columns={c: f"{c}__registry" for c in shared})
    units = score.merge(reg, on="expression_unit_id", how="left", validate="one_to_one")
    for column in shared:
        right = pd.to_numeric(units.pop(f"{column}__registry"), errors="coerce")
        left = pd.to_numeric(units[column], errors="coerce")
        if not np.allclose(left, right, equal_nan=True):
            raise ValueError(f"phase6 score/registry mismatch for {column}")
    fractions = pd.read_csv(cfg["inputs"]["sample_fractions"])
    sample_meta = pd.read_csv(cfg["inputs"]["sample_metadata"], usecols=lambda c: c in WHITELIST_SAMPLE)
    cohort = pd.read_csv(cfg["inputs"]["cohort_registry"], usecols=lambda c: c in WHITELIST_COHORT)
    qc = pd.read_csv(cfg["inputs"]["qc_covariates"], usecols=lambda c: c in WHITELIST_QC)
    fcols = [c for c in fractions if c in {"sample_key", "patient_key", "timepoint", "total_cells"} or c.startswith("frac_mid__") or c.startswith("frac_coarse__")]
    units = units.merge(fractions[fcols], on="sample_key", how="left", suffixes=("", "_fraction"), validate="many_to_one")
    units = units.merge(sample_meta.drop_duplicates("sample_key"), on="sample_key", how="left", suffixes=("", "_metadata"), validate="many_to_one")
    units = units.merge(cohort.drop_duplicates("cohort_id"), on="cohort_id", how="left", validate="many_to_one")
    units = units.merge(qc.drop_duplicates("sample_key"), on="sample_key", how="left", suffixes=("", "_qc"), validate="many_to_one")
    source_timepoint = units.get("timepoint", pd.Series("unknown", index=units.index)).fillna("unknown").astype(str)
    metadata_timepoint = units.get("timepoint_metadata", pd.Series("unknown", index=units.index)).fillna("unknown").astype(str)
    unresolved_timepoint = source_timepoint.str.strip().str.lower().isin({"", "unknown", "nan", "none", "na"})
    source_timepoint = source_timepoint.where(~unresolved_timepoint, metadata_timepoint)
    units["timepoint"] = source_timepoint.map(normalize_timepoint)
    units["patient_key_raw"] = units["patient_key"]
    null_like = units["patient_key"].isna() | units["patient_key"].astype(str).str.lower().isin({"", "unknown", "nan", "none"})
    fallback = units.get("patient_key_fraction")
    units.loc[null_like & fallback.notna(), "patient_key"] = fallback[null_like & fallback.notna()]
    units["patient_binding_source"] = np.where(null_like & fallback.notna(), "phase4b_fraction_fallback", "phase6_frozen_key")
    tissue = units.get("tissue_source", pd.Series("unknown", index=units.index)).fillna("unknown").astype(str)
    lesion = units.get("lesion_context", pd.Series("unknown", index=units.index)).fillna("unknown").astype(str)
    arm = units.get("treatment_context", pd.Series("unknown", index=units.index)).fillna("unknown").astype(str)
    units["tissue_context"] = tissue + "::lesion=" + lesion + "::arm=" + arm
    invalid_patient = units["patient_key"].isna() | units["patient_key"].astype(str).str.lower().isin({"", "unknown", "nan", "none"})
    units["binding_status"] = np.where(invalid_patient, "unresolved_patient_or_sample", "bound_for_patient_context")
    units["patient_timepoint_context_id"] = np.where(
        invalid_patient, "",
        units["cohort_id"].astype(str) + "::" + units["patient_key"].astype(str) + "::"
        + units["timepoint"].astype(str) + "::" + units["tissue_context"],
    )
    units["state_fraction"] = [row.get(f"frac_{row.cell_state_level}__{row.cell_state}", np.nan) for _, row in units.iterrows()]
    state_sets = registry.groupby("cell_state_level")["cell_state"].apply(lambda x: sorted(set(x.astype(str)))).to_dict()
    denominator = {}
    for level, states in state_sets.items():
        cols = [f"frac_{level}__{state}" for state in states if f"frac_{level}__{state}" in fractions]
        if cols:
            denominator[level] = fractions.set_index("sample_key")[cols].sum(axis=1)
    units["eligible_fraction_denominator"] = [denominator.get(row.cell_state_level, pd.Series(dtype=float)).get(row.sample_key, np.nan) for _, row in units.iterrows()]
    units["precision_weight"] = precision_weight(units["n_cells_used"], cfg["measurement"]["precision_weight_reference_cells"])
    units["log_library"] = np.log1p(pd.to_numeric(units["library_size"], errors="coerce"))
    units["log_cells"] = np.log1p(pd.to_numeric(units["n_cells_used"], errors="coerce"))
    units["annotation_coverage"] = np.where(units["cell_state_level"].eq("mid"), units.get("mid_label_coverage"), units.get("coarse_label_coverage"))
    return units


def make_activity(cfg: dict, units: pd.DataFrame) -> pd.DataFrame:
    id_cols = ["expression_unit_id", "cohort_id", "sample_key", "patient_key", "patient_timepoint_context_id", "binding_status", "cell_state_level", "cell_state", "timepoint", "tissue_context",
               "n_cells_used", "library_size", "state_fraction", "precision_weight", "cancer_type", "platform", "tissue_source", "treatment_context"]
    long = units[id_cols + MODULES].melt(id_vars=id_cols, value_vars=MODULES, var_name="module_id", value_name="topic_mass")
    long["activity_lograte"] = np.log1p(1e6 * pd.to_numeric(long["topic_mass"], errors="coerce") / pd.to_numeric(long["library_size"], errors="coerce").replace(0, np.nan))
    long["activity_global_robust_z"] = long.groupby("module_id", group_keys=False)["activity_lograte"].transform(
        lambda x: robust_z(x, pd.Series("all", index=x.index), min_n=10)
    )
    group = long["cohort_id"].astype(str) + "::" + long["cell_state_level"].astype(str) + "::" + long["cell_state"].astype(str)
    long["activity_robust_z"] = long.groupby("module_id", group_keys=False).apply(
        lambda x: robust_z(x["activity_lograte"], group.loc[x.index], cfg["measurement"]["minimum_group_size_for_robust_z"]),
        include_groups=False,
    ).reset_index(level=0, drop=True)
    long["scoring_level"] = "cell_state_expression_unit"
    long["score_provenance"] = "phase6_frozen_topic_mass_per_library_then_robust_z"
    return long


def residualize(cfg: dict, units: pd.DataFrame, activity: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    continuous = ["state_fraction", "log_library", "log_cells", "annotation_coverage", "low_quality_fraction", "mitochondrial_fraction"]
    categorical = ["cell_state", "cancer_type", "platform", "tissue_source", "timepoint", "treatment_context", "cohort_id"]
    fit = units.copy()
    x_full_raw, design_names = build_design(fit, continuous, categorical)
    _, r = np.linalg.qr(x_full_raw, mode="reduced")
    keep = np.abs(np.diag(r)) > max(x_full_raw.shape) * np.finfo(float).eps * np.abs(np.diag(r)).max()
    x_full = x_full_raw[:, keep]
    design_names = [name for name, selected in zip(design_names, keep) if selected]
    design_rank = int(np.linalg.matrix_rank(x_full))
    blocks = {
        "composition": ["state_fraction"],
        "lineage": ["cell_state"],
        "context": ["cancer_type", "platform", "tissue_source", "timepoint", "treatment_context", "cohort_id"],
        "qc": ["log_library", "log_cells", "annotation_coverage", "low_quality_fraction", "mitochondrial_fraction"],
    }
    outputs, variance = [], []
    for module in MODULES:
        y = activity.loc[activity["module_id"].eq(module), "activity_global_robust_z"].to_numpy()
        resid, full_r2 = weighted_residual(y, x_full, fit["precision_weight"].to_numpy())
        tmp = fit[["expression_unit_id", "cohort_id", "sample_key", "patient_key", "patient_timepoint_context_id", "binding_status", "patient_binding_source", "cell_state_level", "cell_state", "timepoint", "tissue_context", "state_fraction", "eligible_fraction_denominator", "precision_weight", "n_cells_used"]].copy()
        tmp["module_id"] = module
        tmp["observed_activity"] = y
        tmp["route_a_residual"] = resid
        outputs.append(tmp)
        valid_resid = np.isfinite(resid)
        technical_rhos = []
        for col in ("log_library", "log_cells"):
            if valid_resid.sum() >= 20:
                technical_rhos.append(abs(spearman_pair(pd.Series(resid[valid_resid]), fit.loc[valid_resid, col].reset_index(drop=True))))
        row = {"module_id": module, "full_model_r2": full_r2, "within_state_residual_fraction": 1 - full_r2 if np.isfinite(full_r2) else np.nan,
               "design_n_columns": x_full.shape[1], "design_rank": design_rank, "design_full_rank": design_rank == x_full.shape[1],
               "design_redundant_columns_removed": int(x_full_raw.shape[1] - x_full.shape[1]),
               "max_abs_residual_technical_spearman": np.nanmax(technical_rhos) if technical_rhos else np.nan}
        for block, names in blocks.items():
            c = [v for v in continuous if v not in names]
            k = [v for v in categorical if v not in names]
            x_reduced, _ = build_design(fit, c, k)
            _, reduced_r2 = weighted_residual(y, x_reduced, fit["precision_weight"].to_numpy())
            row[f"incremental_r2_{block}"] = max(0.0, full_r2 - reduced_r2) if np.isfinite(full_r2) and np.isfinite(reduced_r2) else np.nan
        variance.append(row)
    return pd.concat(outputs, ignore_index=True), pd.DataFrame(variance)


def aggregate_scores(cfg: dict, residual: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    primary = residual[residual["cell_state_level"].eq(cfg["measurement"]["primary_level"])].copy()
    primary = primary[primary["binding_status"].eq("bound_for_patient_context")].copy()
    primary["abundance_weight"] = primary["state_fraction"].fillna(0) * primary["precision_weight"]
    keys = ["cohort_id", "patient_key", "timepoint", "tissue_context", "patient_timepoint_context_id", "sample_key", "module_id"]
    observed = weighted_group_mean(primary, "observed_activity", "abundance_weight", keys)
    coverage = primary.groupby(keys, dropna=False).agg(represented_fraction=("state_fraction", "sum"), eligible_fraction_denominator=("eligible_fraction_denominator", "max"), n_states=("cell_state", "nunique"), n_cells=("n_cells_used", "sum")).reset_index()
    coverage["available_fraction"] = coverage.represented_fraction / coverage.eligible_fraction_denominator.replace(0, np.nan)
    observed = observed.merge(coverage, on=keys, how="left")
    observed["scoring_status"] = np.where((observed.available_fraction >= cfg["measurement"]["minimum_available_fraction"]) &
                                           (observed.n_states >= cfg["measurement"]["minimum_states_per_sample"]), "primary", "low_coverage")
    observed.to_csv(cfg["out"] / "measurement/sample_module_coverage_gate.csv", index=False)
    pkeys = ["cohort_id", "patient_key", "timepoint", "tissue_context", "patient_timepoint_context_id", "module_id"]
    observed["sample_equal_weight"] = 1.0
    observed_valid = observed[observed.scoring_status.eq("primary")].copy()
    observed_patient = weighted_group_mean(observed_valid, "observed_activity", "sample_equal_weight", pkeys)

    valid_samples = observed_valid[keys].drop_duplicates()
    primary_valid = primary.merge(valid_samples.assign(_coverage_gate=True), on=keys, how="inner")
    route_a_sample = weighted_group_mean(primary_valid, "route_a_residual", "precision_weight", keys)
    route_a_sample = route_a_sample.rename(columns={"_w": "sample_measurement_precision"})
    route_a_sample["sample_equal_weight"] = 1.0
    route_a = weighted_group_mean(route_a_sample, "route_a_residual", "sample_equal_weight", pkeys)
    route_a = route_a.rename(columns={"route_a_residual": "residual_activity_route_a"})
    context_precision = route_a_sample.groupby(pkeys, dropna=False)["sample_measurement_precision"].mean().reset_index()
    route_a = route_a.merge(context_precision, on=pkeys, how="left", validate="one_to_one")
    route_a["route_b_prior_mean"] = route_a.groupby(["cohort_id", "module_id"])["residual_activity_route_a"].transform("mean")
    route_a["route_b_prior_variance"] = route_a.groupby(["cohort_id", "module_id"])["residual_activity_route_a"].transform("var").fillna(1.0)
    route_a["sampling_variance"] = 1.0 / route_a["sample_measurement_precision"].clip(lower=0.05)
    shrink = route_a["route_b_prior_variance"] / (route_a["route_b_prior_variance"] + route_a["sampling_variance"])
    route_a["residual_activity_route_b"] = shrink * route_a["residual_activity_route_a"] + (1 - shrink) * route_a["route_b_prior_mean"]
    route_a["route_b_posterior_se"] = np.sqrt((route_a["route_b_prior_variance"] * route_a["sampling_variance"]) /
                                               (route_a["route_b_prior_variance"] + route_a["sampling_variance"]))
    route_a["measurement_contract"] = "route_a_wls_residual_primary;route_b_empirical_bayes_sensitivity"
    wide = route_a.pivot_table(index=["cohort_id", "patient_key", "timepoint", "tissue_context", "patient_timepoint_context_id"], columns="module_id", values="residual_activity_route_a").reset_index()
    return observed_patient, route_a, wide


def aggregate_sensitivity_level(cfg: dict, residual: pd.DataFrame) -> pd.DataFrame:
    data = residual[
        residual.cell_state_level.eq(cfg["measurement"]["sensitivity_level"])
        & residual.binding_status.eq("bound_for_patient_context")
    ].copy()
    sample_keys = ["cohort_id", "patient_key", "timepoint", "tissue_context", "patient_timepoint_context_id", "sample_key", "module_id"]
    coverage = data.groupby(sample_keys, dropna=False).agg(represented_fraction=("state_fraction", "sum"), denominator=("eligible_fraction_denominator", "max"), n_states=("cell_state", "nunique")).reset_index()
    coverage["available_fraction"] = coverage.represented_fraction / coverage.denominator.replace(0, np.nan)
    valid = coverage[(coverage.available_fraction >= cfg["measurement"]["minimum_available_fraction"]) & (coverage.n_states >= cfg["measurement"]["minimum_states_per_sample"])][sample_keys]
    data = data.merge(valid.assign(_coverage_gate=True), on=sample_keys, how="inner")
    sample = weighted_group_mean(data, "route_a_residual", "precision_weight", sample_keys)
    sample["sample_equal_weight"] = 1.0
    context_keys = ["cohort_id", "patient_key", "timepoint", "tissue_context", "patient_timepoint_context_id", "module_id"]
    context = weighted_group_mean(sample, "route_a_residual", "sample_equal_weight", context_keys)
    return context.pivot_table(index=context_keys[:-1], columns="module_id", values="route_a_residual").reset_index()


def bootstrap_stability(cfg: dict, residual: pd.DataFrame, mode: str, allowed_contexts: set[str]) -> pd.DataFrame:
    n_boot = cfg["measurement"][f"bootstrap_n_{mode}"]
    rng = np.random.default_rng(cfg["seed"])
    data = residual[
        residual.binding_status.eq("bound_for_patient_context")
        & residual.cell_state_level.eq(cfg["measurement"]["primary_level"])
        & residual.patient_timepoint_context_id.isin(allowed_contexts)
    ].copy()
    index = ["expression_unit_id", "patient_timepoint_context_id", "precision_weight"]
    unit = data.pivot_table(index=index, columns="module_id", values="route_a_residual").reset_index()
    values = unit[MODULES].to_numpy(dtype=float)
    valid = np.isfinite(values)
    values = np.nan_to_num(values)
    base_weight = unit.precision_weight.to_numpy(dtype=float)
    context_code, contexts = pd.factorize(unit.patient_timepoint_context_id, sort=True)
    n_context = len(contexts)

    def aggregate(weight: np.ndarray) -> np.ndarray:
        numerator = np.zeros((n_context, len(MODULES)))
        denominator = np.zeros((n_context, len(MODULES)))
        np.add.at(numerator, context_code, values * weight[:, None])
        np.add.at(denominator, context_code, valid * weight[:, None])
        return np.divide(numerator, denominator, out=np.full_like(numerator, np.nan), where=denominator > 0)

    base = aggregate(base_weight)
    records = []
    for b in range(n_boot):
        candidate = aggregate(base_weight * rng.exponential(1.0, len(unit)))
        for j, module in enumerate(MODULES):
            records.append({"bootstrap_id": b + 1, "module_id": module,
                            "patient_context_score_rank_stability": rank_stability(pd.Series(base[:, j]), pd.Series(candidate[:, j]))})
    return pd.DataFrame(records)


def semantic_annotation(cfg: dict, membership: pd.DataFrame, residual: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    universe = {g.upper() for g in membership.gene.astype(str)}
    gene_sets = {**parse_gmt(cfg["inputs"]["hallmark_gmt"]), **parse_gmt(cfg["inputs"]["reactome_gmt"])}
    carrier = residual.groupby(["module_id", "cell_state"])["route_a_residual"].median().reset_index()
    carrier = carrier.sort_values(["module_id", "route_a_residual"], ascending=[True, False]).groupby("module_id").head(3)
    rows, enrichment = [], []
    for module in MODULES:
        genes = membership[membership.frozen_module_id.eq(module)].sort_values("membership_weight", ascending=False).gene.astype(str).str.upper().head(50).tolist()
        module_sets = []
        for name, gs in gene_sets.items():
            overlap = set(genes) & gs & universe
            if len(overlap) < 2:
                continue
            p = hypergeom_sf(len(overlap) - 1, len(universe), len(gs & universe), len(set(genes)))
            module_sets.append((name, len(overlap), p, ";".join(sorted(overlap))))
        q = bh_fdr([x[2] for x in module_sets])
        for item, fdr in zip(module_sets, q):
            enrichment.append({"module_id": module, "gene_set": item[0], "overlap_n": item[1], "p_value": item[2], "fdr": fdr, "overlap_genes": item[3]})
        top = sorted(zip(module_sets, q), key=lambda x: x[1] if np.isfinite(x[1]) else 1)[:3]
        pathways = [x[0][0] for x in top if x[1] <= 0.1]
        states = carrier.loc[carrier.module_id.eq(module), "cell_state"].tolist()
        confidence = "medium" if pathways and states else "low"
        name = pathways[0].replace("HALLMARK_", "").replace("REACTOME_", "") if pathways else f"{module}_provisional_expression_topic"
        rows.append({"module_id": module, "provisional_name": name, "name_confidence": confidence, "top_genes": ";".join(genes[:15]),
                     "top_pathways": ";".join(pathways), "carrier_states": ";".join(states), "n_cohorts": residual.loc[residual.module_id.eq(module), "cohort_id"].nunique(),
                     "claim_boundary": "response_blind_candidate_program_not_failure_or_resistance"})
    return pd.DataFrame(rows), pd.DataFrame(enrichment)


def sentinel_scores(cfg: dict, units: pd.DataFrame, mode: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    expr = pd.read_parquet(cfg["inputs"]["phase6_logcpm"])
    id_cols = [c for c in ("expression_unit_id", "cohort_id", "sample_key", "patient_key", "cell_state_level", "cell_state") if c in expr]
    universe = {c.upper(): c for c in expr.columns if c not in id_cols and c not in {"n_cells_used", "library_size", "object_id", "layer_family"}}
    score_rows, eligibility = [], []
    for program, spec in cfg["sentinel_programs"].items():
        requested = [g.upper() for g in spec["genes"]]
        found = [universe[g] for g in requested if g in universe]
        coverage = len(found) / len(requested)
        if coverage >= 0.5 and len(found) >= 4:
            scoreability = "scoreable"
        elif coverage >= 0.4 and len(found) >= 3:
            scoreability = "limited_scoreability"
        else:
            scoreability = "unsupported"
        stability, n_cohorts = np.nan, 0
        if found:
            tmp = expr[id_cols].copy()
            tmp["program_id"] = program
            tmp["activity"] = expr[found].mean(axis=1)
            tmp["activity_robust_z"] = tmp.groupby(["cohort_id", "cell_state_level", "cell_state"], dropna=False)["activity"].transform(lambda x: robust_z(x, pd.Series("g", index=x.index), min_n=5))
            binding = units[["expression_unit_id", "patient_timepoint_context_id", "binding_status", "precision_weight"]].drop_duplicates("expression_unit_id")
            tmp = tmp.merge(binding, on="expression_unit_id", how="left")
            audit = tmp[tmp.cell_state_level.eq(cfg["measurement"]["primary_level"]) & tmp.binding_status.eq("bound_for_patient_context")].copy()
            n_cohorts = audit.cohort_id.nunique()
            if not audit.empty:
                code, contexts = pd.factorize(audit.patient_timepoint_context_id, sort=True)
                values = audit.activity_robust_z.to_numpy(dtype=float)
                weights = audit.precision_weight.to_numpy(dtype=float)
                def agg(w):
                    num, den = np.zeros(len(contexts)), np.zeros(len(contexts))
                    valid = np.isfinite(values)
                    np.add.at(num, code, np.nan_to_num(values) * w)
                    np.add.at(den, code, valid * w)
                    return np.divide(num, den, out=np.full(len(contexts), np.nan), where=den > 0)
                base = agg(weights)
                rng = np.random.default_rng(cfg["seed"] + sum(map(ord, program)))
                reps = cfg["measurement"][f"bootstrap_n_{mode}"]
                stability = float(np.nanmedian([rank_stability(pd.Series(base), pd.Series(agg(weights * rng.exponential(1, len(audit))))) for _ in range(reps)]))
            status = "mainline_candidate" if scoreability == "scoreable" and n_cohorts >= 3 and stability >= cfg["measurement_gate"]["residual_stability_mainline_min"] else "sensitivity_only" if scoreability != "unsupported" else "unsupported"
            tmp["eligibility_status"] = status
            score_rows.append(tmp)
        else:
            status = "unsupported"
        eligibility.append({"program_id": program, "requested_genes": len(requested), "available_genes": len(found), "gene_coverage": coverage,
                            "available_gene_symbols": ";".join(found), "scoreability_status": scoreability,
                            "n_independent_cohorts": n_cohorts, "bootstrap_stability": stability, "eligibility_status": status,
                            "boundary": "prespecified_program_not_frozen_module"})
    return (pd.concat(score_rows, ignore_index=True) if score_rows else pd.DataFrame()), pd.DataFrame(eligibility)


def top40_projection_audit(cfg: dict, membership: pd.DataFrame, activity: pd.DataFrame) -> pd.DataFrame:
    expr = pd.read_parquet(cfg["inputs"]["phase6_logcpm"])
    expr = expr.set_index("expression_unit_id")
    rows = []
    for module in MODULES:
        mem = membership[membership.frozen_module_id.eq(module)].sort_values("membership_weight", ascending=False)
        top = mem.head(40)
        genes = [g for g in top.gene.astype(str) if g in expr.columns]
        weights = top.set_index("gene").reindex(genes).membership_weight.to_numpy(dtype=float)
        projected = expr[genes].to_numpy() @ (weights / weights.sum()) if genes else np.full(len(expr), np.nan)
        projected = pd.Series(projected, index=expr.index, name="top40_projection")
        frozen = activity[activity.module_id.eq(module)].set_index("expression_unit_id")["activity_lograte"]
        rows.append({"module_id": module, "top40_genes_available": len(genes),
                     "top40_membership_weight_fraction": float(top.membership_weight.sum() / mem.membership_weight.sum()),
                     "top40_vs_full_frozen_score_spearman": spearman_pair(frozen, projected)})
    return pd.DataFrame(rows)


def route_agreement(residualized: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for module, data in residualized.groupby("module_id"):
        a, b = data.residual_activity_route_a, data.residual_activity_route_b
        n_top = max(1, int(len(data) * 0.1))
        top_a = set(data.nlargest(n_top, "residual_activity_route_a").index)
        top_b = set(data.nlargest(n_top, "residual_activity_route_b").index)
        rows.append({"module_id": module, "route_ab_spearman": spearman_pair(a, b),
                     "route_ab_sign_agreement": float((np.sign(a) == np.sign(b)).mean()),
                     "route_ab_top_decile_jaccard": len(top_a & top_b) / len(top_a | top_b)})
    return pd.DataFrame(rows)


def decide(cfg: dict, variance: pd.DataFrame, stability: pd.DataFrame, wide: pd.DataFrame, annotation: pd.DataFrame, agreement: pd.DataFrame) -> pd.DataFrame:
    stab = stability.groupby("module_id")["patient_context_score_rank_stability"].median().rename("bootstrap_stability")
    coverage = pd.read_csv(cfg["out"] / "measurement/sample_module_coverage_gate.csv")
    total_contexts = coverage[["cohort_id", "patient_key", "timepoint", "tissue_context", "patient_timepoint_context_id"]].drop_duplicates().shape[0]
    complete = (wide[MODULES].notna().sum() / max(total_contexts, 1)).rename("score_completeness")
    rows = variance.merge(stab, on="module_id", how="left").merge(complete, left_on="module_id", right_index=True, how="left").merge(annotation[["module_id", "name_confidence"]], on="module_id").merge(agreement, on="module_id", how="left")
    rows["n_patient_timepoint_contexts"] = len(wide)
    rows["n_independent_cohorts"] = wide.cohort_id.nunique()
    rows["maximum_single_cohort_fraction"] = wide.cohort_id.value_counts(normalize=True).max()
    gate = cfg["measurement_gate"]
    def status(row):
        if row.score_completeness < gate["completeness_block_min"] or row.bootstrap_stability < gate["residual_stability_block_max"]:
            return "blocked_measurement_invalid"
        if row.max_abs_residual_technical_spearman > gate["technical_residual_abs_spearman_max"]:
            return "sensitivity_only"
        if row.route_ab_spearman < cfg["measurement"]["route_ab_spearman_min"] or row.route_ab_sign_agreement < cfg["measurement"]["route_ab_sign_agreement_min"] or row.route_ab_top_decile_jaccard < cfg["measurement"]["route_ab_top_decile_jaccard_min"]:
            return "sensitivity_only"
        if row.n_independent_cohorts < gate["minimum_independent_cohorts"] or row.n_patient_timepoint_contexts < gate["minimum_patient_timepoint_contexts"] or row.maximum_single_cohort_fraction > gate["maximum_single_cohort_fraction"]:
            return "sensitivity_only"
        if row.incremental_r2_context >= gate["context_variance_context_only_min"]:
            return "context_only"
        if row.incremental_r2_qc >= gate["qc_variance_sensitivity_min"]:
            return "sensitivity_only"
        if row.within_state_residual_fraction >= gate["residual_variance_mainline_min"] and row.bootstrap_stability >= gate["residual_stability_mainline_min"]:
            return "eligible_for_identifiability"
        return "sensitivity_only"
    rows["phase7a_status"] = rows.apply(status, axis=1)
    rows["phase7b_allowed"] = ~rows.phase7a_status.eq("blocked_measurement_invalid")
    rows["decision_rule_version"] = "v6.2.1_preregistered"
    return rows


def run(mode: str = "full") -> dict:
    cfg = load_config()
    ensure_dirs(cfg["out"])
    handoff, score, membership, dictionary = read_phase6(cfg)
    preflight(cfg, handoff, score)

    # The full run can exceed an interactive shell lifetime. Cache deterministic
    # response-blind stages under an input-version token so interruption never
    # requires lowering bootstrap counts or reusing stale metadata.
    cache_inputs = [Path(path) for path in cfg["inputs"].values() if Path(path).exists()]
    cache_payload = "|".join(
        f"{path}:{path.stat().st_size}:{path.stat().st_mtime_ns}" for path in sorted(cache_inputs)
    ) + f"|mode={mode}"
    cache_token = hashlib.sha256(cache_payload.encode()).hexdigest()[:16]
    cache = cfg["out"] / "preflight" / f"phase7a_stage_cache_{cache_token}"
    cache.mkdir(parents=True, exist_ok=True)

    def frame(name: str, builder) -> pd.DataFrame:
        path = cache / f"{name}.parquet"
        if path.exists():
            return pd.read_parquet(path)
        result = builder()
        result.to_parquet(path, index=False)
        return result

    units = frame("units", lambda: prepare_units(cfg, score))
    activity = frame("activity", lambda: make_activity(cfg, units))
    residual_path, variance_path = cache / "residual.parquet", cache / "variance.parquet"
    if residual_path.exists() and variance_path.exists():
        residual, variance = pd.read_parquet(residual_path), pd.read_parquet(variance_path)
    else:
        residual, variance = residualize(cfg, units, activity)
        residual.to_parquet(residual_path, index=False)
        variance.to_parquet(variance_path, index=False)
    observed_path, residualized_path, wide_path = (
        cache / "observed.parquet", cache / "residualized.parquet", cache / "wide.parquet"
    )
    if observed_path.exists() and residualized_path.exists() and wide_path.exists():
        observed = pd.read_parquet(observed_path)
        residualized = pd.read_parquet(residualized_path)
        wide = pd.read_parquet(wide_path)
    else:
        observed, residualized, wide = aggregate_scores(cfg, residual)
        observed.to_parquet(observed_path, index=False)
        residualized.to_parquet(residualized_path, index=False)
        wide.to_parquet(wide_path, index=False)
    coarse_wide = frame("coarse_wide", lambda: aggregate_sensitivity_level(cfg, residual))
    stability = frame(
        "bootstrap_stability",
        lambda: bootstrap_stability(cfg, residual, mode, set(wide.patient_timepoint_context_id)),
    )
    annotation_path, enrichment_path = cache / "annotation.parquet", cache / "enrichment.parquet"
    if annotation_path.exists() and enrichment_path.exists():
        annotation = pd.read_parquet(annotation_path)
        enrichment = pd.read_parquet(enrichment_path)
    else:
        annotation, enrichment = semantic_annotation(cfg, membership, residual)
        annotation.to_parquet(annotation_path, index=False)
        enrichment.to_parquet(enrichment_path, index=False)
    top40 = frame("top40", lambda: top40_projection_audit(cfg, membership, activity))
    annotation = annotation.merge(top40, on="module_id", how="left")
    annotation.loc[annotation.top40_vs_full_frozen_score_spearman < 0.5, "name_confidence"] = "low"
    sent_scores_path, sent_eligibility_path = cache / "sentinel_scores.parquet", cache / "sentinel_eligibility.parquet"
    if sent_scores_path.exists() and sent_eligibility_path.exists():
        sent_scores = pd.read_parquet(sent_scores_path)
        sent_eligibility = pd.read_parquet(sent_eligibility_path)
    else:
        sent_scores, sent_eligibility = sentinel_scores(cfg, units, mode)
        sent_scores.to_parquet(sent_scores_path, index=False)
        sent_eligibility.to_parquet(sent_eligibility_path, index=False)
    agreement = route_agreement(residualized)
    eligibility = decide(cfg, variance, stability, wide, annotation, agreement)

    out = cfg["out"]
    activity.to_csv(out / "measurement/module_activity_by_cell_state.csv", index=False)
    units[["expression_unit_id", "cohort_id", "sample_key", "patient_key_raw", "patient_key", "patient_binding_source", "timepoint", "tissue_context", "patient_timepoint_context_id", "binding_status"]].to_csv(out / "measurement/expression_unit_binding_table.csv", index=False)
    observed.to_csv(out / "measurement/module_activity_by_patient_timepoint_observed.csv", index=False)
    residualized.to_csv(out / "measurement/module_activity_by_patient_timepoint_residualized.csv", index=False)
    wide.to_parquet(out / "measurement/module_activity_by_patient_timepoint_residualized.wide.parquet", index=False)
    coarse_wide.to_parquet(out / "measurement/module_activity_by_patient_timepoint_residualized.coarse_sensitivity.parquet", index=False)
    variance.to_csv(out / "measurement/module_variance_decomposition.csv", index=False)
    stability.to_csv(out / "measurement/module_bootstrap_stability.csv", index=False)
    annotation.to_csv(out / "annotation/module_semantic_annotation_v6_2_1.csv", index=False)
    enrichment.to_csv(out / "annotation/module_pathway_enrichment_v6_2_1.csv", index=False)
    top40.to_csv(out / "annotation/module_top40_projection_sensitivity.csv", index=False)
    annotation[["module_id", "provisional_name", "name_confidence", "claim_boundary"]].to_csv(out / "annotation/module_name_confidence_and_claim_boundary.csv", index=False)
    sent_scores.to_csv(out / "sentinel/prespecified_sentinel_program_scores.csv", index=False)
    sent_eligibility.to_csv(out / "sentinel/prespecified_sentinel_program_eligibility.csv", index=False)
    eligibility.to_csv(out / "handoff/module_eligibility_for_phase7b.csv", index=False)
    corr = wide[MODULES].rank().corr()
    overlap_rows = [{"module_a": a, "module_b": b, "score_spearman": corr.loc[a, b]} for i, a in enumerate(MODULES) for b in MODULES[i + 1:]]
    pd.DataFrame(overlap_rows).to_csv(out / "measurement/module_redundancy_and_overlap_audit.csv", index=False)

    hard = []
    if eligibility.phase7b_allowed.sum() == 0:
        hard.append("no_module_passed_measurement_gate")
    verdict = "BLOCKED" if hard else "GO_TO_PHASE7B" if eligibility.phase7a_status.eq("eligible_for_identifiability").all() else "CONDITIONAL_GO_TO_PHASE7B"
    manifest = {
        "phase": "phase7a_frozen_module_semantic_and_measurement_validity", "mode": mode, "verdict": verdict,
        "response_or_outcome_used": False, "module_discovery_reopened": False, "full_integration_claimed": False,
        "hard_blockers": hard, "conditional_items": sorted(eligibility.loc[~eligibility.phase7a_status.eq("eligible_for_identifiability"), "module_id"].tolist()),
        "phase7b_unique_inputs": {
            "residualized_module_matrix": str((out / "measurement/module_activity_by_patient_timepoint_residualized.wide.parquet").relative_to(ROOT)),
            "module_eligibility": str((out / "handoff/module_eligibility_for_phase7b.csv").relative_to(ROOT)),
            "coarse_sensitivity_matrix": str((out / "measurement/module_activity_by_patient_timepoint_residualized.coarse_sensitivity.parquet").relative_to(ROOT)),
            "semantic_annotation": str((out / "annotation/module_semantic_annotation_v6_2_1.csv").relative_to(ROOT)),
            "sentinel_eligibility": str((out / "sentinel/prespecified_sentinel_program_eligibility.csv").relative_to(ROOT)),
        },
        "phase6_asset_hashes": {k: sha256(Path(handoff[k])) for k in ("default_module_score_matrix", "default_module_membership", "module_dictionary")},
    }
    write_yaml(manifest, out / "handoff/phase7a_to_phase7b_handoff.yaml")
    report = f"""# Phase7A Module Measurement Validity Report

## Executive Verdict

**{verdict}**。Phase6 的 8 个 frozen candidate expression topics 未重训；本阶段未读取 response、outcome 或 split 构造分数。

## Measurement Contract

- Cell-state activity: topic mass / library size 后取 log rate，并在 cohort × resolution × cell state 内稳健标准化。
- Patient-timepoint observed: mid-state abundance × precision 加权汇总。
- Route A: 显式调整 composition、lineage、cohort/cancer/platform、tissue/timepoint/treatment 与 QC。
- Route B: 对 Route A patient-timepoint 结果做 cohort 内 empirical-Bayes partial pooling，仅作增强路线。

## Results

- expression units: {len(units):,}
- patient-timepoints: {len(wide):,}
- Phase7B fully measured topics: {int(eligibility.phase7a_status.eq('eligible_for_identifiability').sum())}/8
- Phase7B sensitivity/context audit topics: {int(eligibility.phase7b_allowed.sum() - eligibility.phase7a_status.eq('eligible_for_identifiability').sum())}/8
- Coverage-qualified patient-context completeness: {eligibility.score_completeness.median():.3f}
- sentinel programs: {len(sent_eligibility)}，其中 mainline candidate {int(sent_eligibility.eligibility_status.eq('mainline_candidate').sum())}。

## Boundaries

这些 topic 仍不是 failure/resistance barrier。本阶段不做 response association、不重新发现 module，也不把 sentinel 自动补入 frozen module set。
"""
    (out / "PHASE7A_MODULE_MEASUREMENT_VALIDITY_REPORT.md").write_text(report)
    return manifest


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--mode", choices=["smoke", "full"], default="full")
    args = parser.parse_args()
    print(yaml.safe_dump(run(args.mode), sort_keys=False))
