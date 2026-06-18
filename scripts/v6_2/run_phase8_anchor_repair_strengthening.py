#!/usr/bin/env python3
"""Phase8 anchor repair and bounded statistical strengthening.

The runner preserves FM01/FM04/FM07 membership, applies source-proven metadata
overlays, and publishes atomically into a new result directory. Response is
never used to construct or select barrier scores.
"""
from __future__ import annotations

import hashlib
import math
import os
import shutil
import atexit
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

from phase7_common import ROOT, normalize_timepoint, sha256
from run_phase8_anchor_context import auc_rank, bayesian_partial_pool, hedges_g, load_config, random_effects
from run_phase8_bounded_repair import (
    PRIMARY,
    cliffs_delta,
    coverage_audit,
    repaired_environment,
    sensitivity_and_sentinels,
)

OUT = Path(os.environ.get(
    "PHASE8_OUTPUT_DIR",
    ROOT / "results/v6_2/phase8_anchor_repair_and_statistical_strengthening",
))
STAGING = OUT.parent / f".{OUT.name}.staging"
P7 = Path(os.environ.get(
    "PHASE7_RESULT_DIR",
    ROOT / "results/v6_2/phase7_module_measurement_and_barrier_identifiability",
))
P6 = Path(os.environ.get(
    "PHASE6_RESULT_DIR",
    ROOT / "results/v6_2/phase6_module_algorithm_benchmark_strengthened",
))
OLD = Path(os.environ.get(
    "PHASE8_REFERENCE_DIR",
    ROOT / "results/v6_2/phase8_anchor_context_adjudication",
))
PHASE4B = Path(os.environ.get(
    "PHASE4B_RESULT_DIR",
    ROOT / "results/v6_2/phase4b_immune_state_feature_construction",
))
PHASE4A = Path(os.environ.get(
    "PHASE4A_RESULT_DIR",
    ROOT / "results/v6_2/phase4a_cell_state_harmonization",
))
SEED = 1729
MIN_CLASS = 3

GSE120_SUPP = ROOT / "data/combo/GSE120575/NIHMS1510803-supplement-10.xlsx"
GSE120_H5AD = ROOT / "data/processed/srt/raw/gse120575.h5ad"
GSE120_ADDENDUM = ROOT / "results/v6_2/phase2_5_data_onboarding/01_metadata_freeze/sample_metadata_addendum.GSE120575_rescue.tsv"
GSE123_SUPP = ROOT / "data/imm/GSE123813/GSE123813_Yost2019_Supplementary_Tables.xlsx"
GSE243_META = ROOT / "data/imm/GSE243013/GSE243013_NSCLC_immune_scRNA_metadata.csv.gz"


def now_iso() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat()


def write_yaml(path: Path, obj: dict) -> None:
    path.write_text(yaml.safe_dump(obj, sort_keys=False, allow_unicode=True), encoding="utf-8")


def normalized_counts_audit(path: Path) -> dict:
    import anndata as ad

    a = ad.read_h5ad(path, backed="r")
    layer = "counts" if "counts" in a.layers else "X"
    x = a.layers[layer][: min(200, a.n_obs), : min(200, a.n_vars)] if layer != "X" else a.X[: min(200, a.n_obs), : min(200, a.n_vars)]
    x = x.toarray() if hasattr(x, "toarray") else np.asarray(x)
    finite = x[np.isfinite(x)]
    noninteger_fraction = float(np.mean(np.abs(finite - np.round(finite)) > 1e-6)) if finite.size else np.nan
    result = {
        "layer": layer,
        "n_obs": int(a.n_obs),
        "n_vars": int(a.n_vars),
        "sample_noninteger_fraction": noninteger_fraction,
        "sample_max": float(np.max(finite)) if finite.size else np.nan,
        "raw_count_semantics_pass": bool(noninteger_fraction <= 0.01),
    }
    a.file.close()
    return result


def source_overlays() -> tuple[pd.DataFrame, pd.DataFrame, dict]:
    """Build patient-level response/treatment overlays without score access."""
    p120 = pd.read_excel(GSE120_SUPP, sheet_name="Patient-scRNA data", header=1)
    response_col = "Clinical response (RECIST; R=CR, PR; NR=SD, PD)"
    p120 = p120[["Patient ID", "Therapy", response_col]].dropna(subset=["Patient ID"])
    p120["patient_key"] = "GSE120575::gse120575_" + p120["Patient ID"].astype(str)
    p120["response_binary"] = np.where(p120[response_col].astype(str).str.strip().eq("R"), 1, 0)
    p120["response_raw_repaired"] = p120[response_col].astype(str)
    p120["endpoint_type_repaired"] = "RECIST"
    p120["treatment_arm_repaired"] = np.where(p120.Therapy.eq("PD1"), "monotherapy", "combination_or_sequential")
    p120["overlay_source"] = str(GSE120_SUPP.relative_to(ROOT)) + "::Patient-scRNA data"

    p123 = pd.read_excel(GSE123_SUPP, sheet_name="SuppTable1", header=3).dropna(subset=["Patient"])
    p123 = p123[p123["Tumor Type"].isin(["BCC", "SCC"])].copy()
    p123["source_patient"] = p123.Patient.astype(str)
    p123["cohort_id"] = np.where(p123["Tumor Type"].eq("BCC"), "GSE123813_bcc", "GSE123813_scc")
    p123["patient_id_join"] = p123.source_patient.str.replace("-S", "", regex=False)
    # BCC and SCC are lesion contexts from one study subject. Cohort remains
    # available for lesion-specific response binding and is not part of the
    # independent patient key.
    p123["patient_key"] = "GSE123813::" + p123.patient_id_join
    response_text = p123.Response.astype(str).str.strip()
    p123["response_binary"] = np.where(response_text.str.startswith("Yes"), 1, np.where(response_text.str.startswith("No"), 0, np.nan))
    p123["response_raw_repaired"] = p123.Response.astype(str)
    p123["endpoint_type_repaired"] = "RECIST"
    p123["treatment_arm_repaired"] = np.where(
        p123["Ongoing Vismodegib treatment"].astype(str).str.strip().eq("+"),
        "combination_PD1_plus_vismodegib",
        "monotherapy",
    )
    p123["overlay_source"] = str(GSE123_SUPP.relative_to(ROOT)) + "::SuppTable1"
    audit120 = normalized_counts_audit(GSE120_H5AD)
    return p120, p123, audit120


def repair_arm_classifier(env: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Correct GSE243013 arm assignment using explicit source columns."""
    meta = pd.read_csv(GSE243_META)
    cols = [c for c in ["patient", "Patient", "patient_id", "sampleID", "anti-PD1_therapy", "chemotherapy", "targeted_therapy", "pathological_response", "radiological_response"] if c in meta]
    meta = meta[cols].copy()
    patient_col = next(c for c in ["patient", "Patient", "patient_id", "sampleID"] if c in meta)
    patient = meta.groupby(patient_col, as_index=False).first()
    patient["patient_key"] = "GSE243013::" + patient[patient_col].astype(str)
    clean = (
        patient.get("anti-PD1_therapy", pd.Series("", index=patient.index)).astype(str).str.strip().ne("")
        & patient.get("chemotherapy", pd.Series("", index=patient.index)).astype(str).str.lower().isin(["no", "n", "0", "false"])
        & patient.get("targeted_therapy", pd.Series("", index=patient.index)).astype(str).str.lower().isin(["no", "n", "0", "false"])
    )
    patient["source_clean_monotherapy"] = clean
    env = env.merge(patient[["patient_key", "source_clean_monotherapy"]], on="patient_key", how="left")
    mask = env.cohort_id.eq("GSE243013")
    env.loc[mask, "treatment_arm"] = np.where(env.loc[mask, "source_clean_monotherapy"].fillna(False), "monotherapy", "combination_or_other")
    env.loc[mask, "treatment_cleanliness"] = np.where(env.loc[mask, "source_clean_monotherapy"].fillna(False), "source_clean_but_frozen_support_only", "high_confounding")
    ledger = patient.assign(
        record_type="GSE243013_arm_repair",
        status=np.where(clean, "clean_monotherapy_source_verified", "combination_or_other"),
        detail="explicit anti-PD1/chemotherapy/targeted-therapy fields; no string-substring inference",
    )
    return env, ledger


def build_anchor_master(cfg: dict) -> tuple[pd.DataFrame, pd.DataFrame, dict]:
    env, _ = repaired_environment(cfg)
    env["all_primary_scores_available"] = env[PRIMARY].notna().all(axis=1)
    env["score_source"] = env.score_source.fillna("phase7_frozen_mid_barrier_matrix")
    env["measurement_route"] = env.measurement_route.fillna("direct_mid")
    fallback_coverage = pd.Series(np.where(env.all_primary_scores_available, "phase7_primary", "missing"), index=env.index)
    env["coverage_status"] = env.coverage_status.fillna(fallback_coverage)
    env["primary_role"] = env.primary_role.fillna("not_primary")
    env["sensitivity_role"] = env.sensitivity_role.fillna("registry_or_support")
    env["endpoint_quality"] = env.endpoint_quality.fillna("unknown")

    p120, p123, audit120 = source_overlays()
    overlay_rows = []
    map120 = p120.set_index("patient_key")
    mask120 = env.cohort_id.eq("GSE120575") & env.patient_key.isin(map120.index)
    for idx in env.index[mask120]:
        x = map120.loc[env.at[idx, "patient_key"]]
        env.at[idx, "response_binary"] = int(x.response_binary)
        env.at[idx, "response_raw"] = x.response_raw_repaired
        env.at[idx, "endpoint_type"] = "RECIST"
        env.at[idx, "endpoint_quality"] = "high_source_supplement"
        env.at[idx, "treatment_arm"] = x.treatment_arm_repaired
        env.at[idx, "treatment_cleanliness"] = "clean" if x.treatment_arm_repaired == "monotherapy" else "confounded"
        env.at[idx, "primary_role"] = "not_primary"
        env.at[idx, "sensitivity_role"] = "normalized_expression_projection_support"
        env.at[idx, "score_source"] = "phase7_score_invalid_for_raw_count_claim"
        env.at[idx, "measurement_route"] = "normalized_expression_support_only"
    overlay_rows.append({
        "record_type": "GSE120575_expression_semantics",
        "dataset_id": "GSE120575",
        "status": "support_only_normalized_expression",
        "n_patients": int(p120[p120.treatment_arm_repaired.eq("monotherapy")].patient_key.nunique()),
        "detail": f"counts-layer noninteger fraction={audit120['sample_noninteger_fraction']:.4f}; raw-count semantics pass={audit120['raw_count_semantics_pass']}",
    })
    add120 = pd.read_csv(GSE120_ADDENDUM, sep="\t")
    baseline_mono = add120[add120.timepoint.astype(str).str.lower().eq("pre") & add120.treatment_context.eq("anti-PD1")]
    cov120 = pd.read_csv(P7 / "measurement/sample_module_coverage_gate.csv")
    cov120 = cov120[cov120.cohort_id.eq("GSE120575") & cov120.sample_key.isin("GSE120575::" + baseline_mono.sample_id.astype(str))]
    cov_patient = cov120.groupby("patient_key").scoring_status.agg(lambda x: "primary" if (x == "primary").all() else "low_coverage")
    overlay_rows.append({
        "record_type": "GSE120575_baseline_monotherapy_inventory", "dataset_id": "GSE120575",
        "status": "normalized_expression_support_only", "n_patients": int(baseline_mono.patient_id.nunique()),
        "detail": f"source baseline mono=12; Phase7 direct-score rows primary={int((cov_patient == 'primary').sum())}; low_coverage={int((cov_patient == 'low_coverage').sum())}; raw-count semantics failed",
    })

    map123 = p123.set_index(["cohort_id", "patient_key"])
    env123_keys = pd.MultiIndex.from_frame(env[["cohort_id", "patient_key"]])
    mask123 = env123_keys.isin(map123.index)
    for idx in env.index[mask123]:
        x = map123.loc[(env.at[idx, "cohort_id"], env.at[idx, "patient_key"])]
        env.at[idx, "response_binary"] = int(x.response_binary)
        env.at[idx, "response_raw"] = x.response_raw_repaired
        env.at[idx, "endpoint_type"] = "RECIST"
        env.at[idx, "endpoint_quality"] = "high_source_supplement_RECIST_v1_1"
        env.at[idx, "treatment_arm"] = x.treatment_arm_repaired
        clean = x.treatment_arm_repaired == "monotherapy"
        env.at[idx, "treatment_cleanliness"] = "clean" if clean else "PD1_plus_vismodegib"
        lesion_discordant = env.at[idx, "patient_key"] == "GSE123813::su010"
        env.at[idx, "primary_role"] = (
            "conditional_primary"
            if clean and env.at[idx, "normalized_timepoint"] == "baseline" and not lesion_discordant
            else "not_primary"
        )
        env.at[idx, "sensitivity_role"] = (
            "lesion_discordant_clustered_sensitivity" if lesion_discordant
            else "treatment_context_sensitivity" if not clean
            else "clean_anchor"
        )
        env.at[idx, "score_source"] = "phase7_frozen_mid_barrier_matrix_with_source_RECIST_overlay"
    overlay_rows.append({
        "record_type": "GSE123813_response_rescue",
        "dataset_id": "GSE123813_bcc+scc",
        "status": "source_RECIST_restored",
        "n_patients": int(p123.patient_key.nunique()),
        "detail": f"public Supplementary Table 1; clean mono={int(p123.treatment_arm_repaired.eq('monotherapy').sum())}; vismodegib combination={int(p123.treatment_arm_repaired.ne('monotherapy').sum())}",
    })
    env, arm_ledger = repair_arm_classifier(env)

    env["baseline_eligible"] = env.normalized_timepoint.eq("baseline")
    env["response_known"] = pd.to_numeric(env.response_binary, errors="coerce").isin([0, 1])
    env["statistical_estimability"] = "not_eligible"
    strict = env.primary_role.eq("conditional_primary") & env.baseline_eligible & env.response_known & env[PRIMARY].notna().all(axis=1)
    counts = env[strict].groupby("cohort_id").response_binary.agg([lambda x: int((x == 1).sum()), lambda x: int((x == 0).sum())])
    counts.columns = ["n_R", "n_NR"]
    for cohort, row in counts.iterrows():
        status = "estimable" if min(row.n_R, row.n_NR) >= MIN_CLASS else "weakly_estimable"
        env.loc[strict & env.cohort_id.eq(cohort), "statistical_estimability"] = status
    env.loc[strict & env.statistical_estimability.eq("not_eligible"), "statistical_estimability"] = "weakly_estimable"
    env["exclusion_reason"] = np.select(
        [~env.baseline_eligible, ~env.response_known, env.primary_role.ne("conditional_primary"), ~env[PRIMARY].notna().all(axis=1)],
        ["NONBASELINE", "RESPONSE_UNKNOWN", "NOT_CLEAN_PRIMARY_ROLE", "MISSING_PRIMARY_BARRIER"],
        default="",
    )
    env = env.sort_values(["cohort_id", "patient_key", "normalized_timepoint", "tissue_context", "resolution"]).reset_index(drop=True)
    bridge = pd.concat([pd.DataFrame(overlay_rows), arm_ledger], ignore_index=True, sort=False)
    return env, bridge, audit120


def gse301_loss_trace(cfg: dict) -> pd.DataFrame:
    sample = pd.read_csv(cfg["inputs"]["sample_metadata"])
    sample = sample[sample.cohort_id.eq("GSE301741")][
        ["sample_key", "legacy_sample_keys", "patient_key", "timepoint", "response_binary_harmonized"]
    ].copy()
    sample["sample_id"] = (
        sample["legacy_sample_keys"].fillna("").astype(str)
        .str.split("|").str[0].str.rsplit("::", n=1).str[-1]
    )
    elig = pd.read_csv(PHASE4A / "handoff/phase4a_to_phase4b_aggregation_eligibility.csv")
    elig = elig[(elig.cohort_id.eq("GSE301741")) & elig.allowed_phase4b_pseudobulk.eq("yes")]
    ec = elig.groupby("sample_key").agg(
        n_frozen_eligible_states=("harmonized_mid_label", "nunique"),
        n_frozen_eligible_cells=("n_cells", "sum"),
    ).reset_index()
    from phase7_common import load_config as load_phase7_config
    cfg7 = load_phase7_config()
    phase6_handoff = yaml.safe_load(cfg7["inputs"]["phase6_handoff"].read_text())
    score = pd.read_parquet(phase6_handoff["default_module_score_matrix"])
    sc = score[score.cohort_id.eq("GSE301741")].groupby("sample_key").agg(n_frozen_mid_units=("expression_unit_id", "nunique")).reset_index()
    cov = pd.read_csv(P7 / "measurement/sample_module_coverage_gate.csv")
    cov = cov[cov.cohort_id.eq("GSE301741")].groupby("sample_key").agg(
        phase7_scoring_status=("scoring_status", lambda x: "primary" if (x == "primary").all() else "low_coverage"),
        available_fraction=("available_fraction", "mean"), n_states=("n_states", "max"),
    ).reset_index()
    direct = pd.read_parquet(cfg["inputs"]["barrier_matrix"])
    direct = set(direct.loc[direct.cohort_id.eq("GSE301741"), "patient_timepoint_context_id"])
    out = sample.merge(ec, on="sample_key", how="left").merge(sc, on="sample_key", how="left").merge(cov, on="sample_key", how="left")
    out["normalized_timepoint"] = out.timepoint.map(normalize_timepoint)
    out["patient_timepoint_context_id"] = out.patient_key + "::" + out.normalized_timepoint + "::tumor"
    out["direct_mid_available"] = out.patient_timepoint_context_id.isin(direct)
    out["record_type"] = "GSE301741_frozen_direct_mid_reconstruction"
    out["dataset_id"] = "GSE301741"
    out["status"] = np.select(
        [out.direct_mid_available, out.n_frozen_eligible_states.fillna(0).eq(0), out.phase7_scoring_status.eq("low_coverage")],
        ["direct_mid_available", "BRIDGE_FAILED_NO_FROZEN_ELIGIBLE_MID_STATE", "LOW_COVERAGE"],
        default="frozen_interface_missing_or_unresolved",
    )
    out["detail"] = "response-blind replay of frozen Phase4A eligibility; thresholds not relaxed"
    return out


def firth_logistic(data: pd.DataFrame, feature: str) -> dict:
    d = data[["response_binary", feature]].apply(pd.to_numeric, errors="coerce").dropna()
    y = d.response_binary.to_numpy(float)
    x = d[feature].to_numpy(float)
    if len(d) < 6 or len(np.unique(y)) < 2 or np.std(x) < 1e-10:
        return {"firth_log_odds": np.nan, "firth_se": np.nan, "firth_converged": False}
    z = (x - x.mean()) / x.std(ddof=0)
    X = np.column_stack([np.ones(len(z)), z])
    beta = np.zeros(2)
    converged = False
    for _ in range(100):
        p = 1 / (1 + np.exp(-np.clip(X @ beta, -25, 25)))
        w = np.maximum(p * (1 - p), 1e-8)
        info = X.T @ (w[:, None] * X)
        inv = np.linalg.pinv(info)
        h = np.sum((X @ inv) * X, axis=1) * w
        score = X.T @ (y - p + h * (0.5 - p))
        step = inv @ score
        beta += step
        if np.max(np.abs(step)) < 1e-8:
            converged = True
            break
    se = math.sqrt(max(np.linalg.pinv(info)[1, 1], 0))
    return {"firth_log_odds": float(beta[1]), "firth_se": float(se), "firth_converged": converged}


def within_effects(anchor: pd.DataFrame) -> pd.DataFrame:
    primary = anchor[
        anchor.primary_role.eq("conditional_primary")
        & anchor.baseline_eligible
        & anchor.response_known
        & anchor[PRIMARY].notna().all(axis=1)
    ].copy()
    if primary.duplicated(["cohort_id", "patient_key"]).any():
        raise RuntimeError("Primary response surface has repeated patient contributions across context/resolution.")
    primary = primary.groupby(
        ["cohort_id", "patient_key", "endpoint_type", "treatment_arm", "resolution", "cancer_type", "statistical_estimability"],
        as_index=False, dropna=False,
    ).agg({**{fm: "mean" for fm in PRIMARY}, "response_binary": "first"})
    rows = []
    for keys, d in primary.groupby(["cohort_id", "endpoint_type", "treatment_arm", "resolution", "cancer_type", "statistical_estimability"], dropna=False):
        for fm in PRIMARY:
            h = hedges_g(d, fm)
            c = cliffs_delta(d, fm, seed=SEED + PRIMARY.index(fm))
            f = firth_logistic(d, fm)
            signs = [np.sign(v) for v in [h.get("effect"), c.get("cliffs_delta"), f.get("firth_log_odds")] if np.isfinite(v)]
            rows.append(dict(zip(["cohort_id", "endpoint_type", "treatment_arm", "resolution", "cancer_type", "statistical_estimability"], keys)) | {
                "record_type": "within_cohort_effect", "barrier_id": fm, **h, **c, **f,
                "method_direction_agreement": bool(len(signs) >= 2 and len(set(signs)) == 1),
                "abstention_code": "" if keys[-1] == "estimable" else "NOT_ESTIMABLE_ONE_CLASS_OR_SMALL_CLASS",
            })
    return pd.DataFrame(rows)


def endpoint_meta(effects: pd.DataFrame) -> pd.DataFrame:
    valid = effects[(effects.statistical_estimability == "estimable")].dropna(subset=["effect", "se"])
    valid = valid.copy()
    valid["arm_family"] = np.where(valid.treatment_arm.astype(str).str.contains("monotherapy", case=False), "monotherapy", "other")
    rows = []
    for (fm, endpoint, resolution, arm_family), d in valid.groupby(["barrier_id", "endpoint_type", "resolution", "arm_family"]):
        m = random_effects(d)
        m.update({"record_type": "endpoint_meta", "barrier_id": fm, "endpoint_type": endpoint,
                  "resolution": resolution, "arm_family": arm_family, "n_independent_cohorts": int(d.cohort_id.nunique())})
        if len(d) >= 2:
            w = 1 / (d.se.to_numpy() ** 2 + m["tau2"])
            qstar = float(np.sum(w * (d.effect.to_numpy() - m["pooled_effect"]) ** 2) / max(len(d) - 1, 1))
            hk_se = math.sqrt(max(qstar, 1e-12) / w.sum())
            tcrit = 12.706 if len(d) == 2 else 4.303 if len(d) == 3 else 3.182
            m.update(hartung_knapp_se=hk_se, hartung_knapp_ci_low=m["pooled_effect"] - tcrit * hk_se,
                     hartung_knapp_ci_high=m["pooled_effect"] + tcrit * hk_se)
        m["direction_consistency"] = float((np.sign(d.effect) == np.sign(d.effect.mean())).mean()) if len(d) else np.nan
        m["abstention_code"] = "" if d.cohort_id.nunique() >= 2 else "ENDPOINT_UNREPLICATED"
        rows.append(m)
    for fm, d in valid.groupby("barrier_id"):
        m = random_effects(d)
        m.update({"record_type": "cross_endpoint_direction_only", "barrier_id": fm, "endpoint_type": "mixed_not_pooled_for_claim",
                  "direction_consistency": float((np.sign(d.effect) == np.sign(d.effect.mean())).mean()),
                  "abstention_code": "CONTEXT_CONFOUNDED"})
        rows.append(m)
    return pd.DataFrame(rows)


def bayes_sensitivity(effects: pd.DataFrame, cfg: dict) -> pd.DataFrame:
    valid = effects[(effects.statistical_estimability == "estimable")].dropna(subset=["effect", "se"])
    rows = []
    for fm, d in valid.groupby("barrier_id"):
        replicated_endpoint = d.groupby("endpoint_type").size().max() >= 2
        identifiable = len(d) >= 4 and d.endpoint_type.nunique() >= 2 and replicated_endpoint
        for mu_sd in (0.5, 1.0, 2.0, 4.0):
            for tau_sd in (0.5, 1.0, 2.0):
                b = bayesian_partial_pool(d, cfg, mu_sd, tau_sd)
                rows.append({"record_type": "bayesian_prior_sensitivity", "barrier_id": fm, "endpoint_type": "mixed_hierarchical_sensitivity",
                             "identifiable": identifiable, "can_trigger_route": False, **b,
                             "abstention_code": "" if identifiable else "CONTEXT_CONFOUNDED"})
    return pd.DataFrame(rows)


def crossfit_effects(cfg: dict, anchor: pd.DataFrame, estimable_cohorts: list[str]) -> pd.DataFrame:
    from phase7_common import build_design, load_config as load_phase7_config
    from run_phase7a_measurement import make_activity, prepare_units, read_phase6

    cfg7 = load_phase7_config()
    _, score, _, _ = read_phase6(cfg7)
    units = prepare_units(cfg7, score)
    activity = make_activity(cfg7, units)
    continuous = ["state_fraction", "log_library", "log_cells", "annotation_coverage", "low_quality_fraction", "mitochondrial_fraction"]
    categorical = ["cell_state", "cancer_type", "platform", "tissue_source", "timepoint", "treatment_context"]
    xraw, _ = build_design(units, continuous, categorical)
    _, r = np.linalg.qr(xraw, mode="reduced")
    diag = np.abs(np.diag(r)); keep = diag > max(xraw.shape) * np.finfo(float).eps * max(diag.max(), 1e-12)
    x = xraw[:, keep]
    labels = anchor[["cohort_id", "patient_key", "response_binary", "primary_role", "endpoint_type", "treatment_arm", "resolution"]].drop_duplicates(["cohort_id", "patient_key"])
    rows = []
    for target in estimable_cohorts:
        train = units.cohort_id.ne(target).to_numpy(); held = units.cohort_id.eq(target).to_numpy()
        weights = units.precision_weight.to_numpy(float)
        module_context = []
        for fm in PRIMARY:
            y = activity.loc[activity.module_id.eq(fm), "activity_global_robust_z"].to_numpy(float)
            valid = train & np.isfinite(y) & (weights > 0) & np.isfinite(x).all(axis=1)
            rw = np.sqrt(weights[valid]); xv = x[valid]
            gram = (xv * rw[:, None]).T @ (xv * rw[:, None])
            penalty = np.eye(gram.shape[0]) * 1e-6; penalty[0, 0] = 0
            beta = np.linalg.solve(gram + penalty, (xv * rw[:, None]).T @ (y[valid] * rw))
            idx = np.flatnonzero(held)
            z = units.iloc[idx][["cohort_id", "patient_key", "timepoint", "tissue_context", "patient_timepoint_context_id", "cell_state_level", "state_fraction", "eligible_fraction_denominator", "precision_weight", "cell_state"]].copy()
            z["barrier_id"] = fm; z["crossfit_residual"] = y[idx] - x[idx] @ beta
            module_context.append(z)
        long = pd.concat(module_context, ignore_index=True)
        long = long[long.cell_state_level.eq("mid")]
        keys = ["cohort_id", "patient_key", "timepoint", "tissue_context", "patient_timepoint_context_id", "barrier_id"]
        cov = long.groupby(keys, as_index=False).agg(represented=("state_fraction", "sum"), denominator=("eligible_fraction_denominator", "max"), n_states=("cell_state", "nunique"))
        cov["coverage"] = cov.represented / cov.denominator.replace(0, np.nan)
        valid_keys = cov[(cov.coverage >= .5) & (cov.n_states >= 2)][keys]
        long = long.merge(valid_keys.assign(ok=True), on=keys, how="inner")
        long["wx"] = long.crossfit_residual * long.precision_weight
        ctx = long.groupby(keys, as_index=False).agg(wx=("wx", "sum"), w=("precision_weight", "sum"))
        ctx["score"] = ctx.wx / ctx.w
        wide = ctx.pivot_table(index=keys[:-1], columns="barrier_id", values="score").reset_index()
        wide["timepoint"] = [
            normalize_timepoint(value, cohort)
            for value, cohort in zip(wide["timepoint"], wide["cohort_id"])
        ]
        wide = wide.merge(labels, on=["cohort_id", "patient_key"], how="left")
        wide = wide[(wide.timepoint.eq("baseline")) & wide.primary_role.eq("conditional_primary")]
        wide = wide.groupby(["cohort_id", "patient_key", "endpoint_type", "treatment_arm", "resolution", "response_binary"], as_index=False)[PRIMARY].mean()
        for fm in PRIMARY:
            h = hedges_g(wide, fm)
            rows.append({"record_type": "outer_cohort_crossfit_effect", "cohort_id": target, "barrier_id": fm,
                         "endpoint_type": wide.endpoint_type.iloc[0] if len(wide) else "unknown", **h,
                         "abstention_code": "" if np.isfinite(h["effect"]) else "NOT_ESTIMABLE_ONE_CLASS"})
    return pd.DataFrame(rows)


def coverage_and_ipw(cfg: dict, anchor: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    cov, _ = coverage_audit(cfg, anchor.rename(columns={"endpoint_type": "response_endpoint_type"}))
    label = anchor[["cohort_id", "patient_key", "response_binary"]].drop_duplicates(["cohort_id", "patient_key"])
    cov = cov.drop(columns=["response_binary"], errors="ignore").merge(label, on=["cohort_id", "patient_key"], how="left")
    features = ["available_fraction", "n_states", "n_cells", "total_cells", "mid_label_coverage", "low_quality_fraction"]
    xnum = cov[features].apply(pd.to_numeric, errors="coerce").fillna(cov[features].apply(pd.to_numeric, errors="coerce").median()).fillna(0)
    xd = pd.get_dummies(cov.cohort_id, prefix="cohort", dtype=float)
    X = pd.concat([xnum, xd], axis=1).to_numpy(float); y = cov.included.to_numpy(int)
    pred = np.full(len(cov), np.nan)
    for cohort in sorted(cov.cohort_id.unique()):
        test = cov.cohort_id.eq(cohort).to_numpy(); train = ~test
        if len(np.unique(y[train])) < 2:
            continue
        mean = X[train].mean(axis=0); sd = X[train].std(axis=0); sd[sd < 1e-8] = 1.0
        xtr = (X[train] - mean) / sd; xte = (X[test] - mean) / sd
        design = np.column_stack([np.ones(len(xtr)), xtr]); beta = np.zeros(design.shape[1])
        penalty = np.eye(design.shape[1]) * .1; penalty[0, 0] = 0
        for _ in range(100):
            prob = 1 / (1 + np.exp(-np.clip(design @ beta, -20, 20)))
            weight = np.maximum(prob * (1 - prob), 1e-5)
            hessian = design.T @ (weight[:, None] * design) + penalty
            step = np.linalg.solve(hessian, design.T @ (y[train] - prob) - penalty @ beta)
            beta += step
            if np.max(np.abs(step)) < 1e-7:
                break
        pred[test] = 1 / (1 + np.exp(-np.clip(np.column_stack([np.ones(len(xte)), xte]) @ beta, -20, 20)))
    cov["oof_inclusion_probability"] = pred
    prevalence = float(np.mean(y))
    cov["stabilized_ipw"] = np.where(cov.included.eq(1), prevalence / np.clip(pred, .05, .95), np.nan)
    w = cov.loc[cov.included.eq(1), "stabilized_ipw"].dropna()
    ess = float(w.sum() ** 2 / w.pow(2).sum()) if len(w) else np.nan
    positivity = bool(((pred < .05) | (pred > .95)).mean() <= .10) if np.isfinite(pred).any() else False
    summary = pd.DataFrame([{
        "record_type": "coverage_selection_oof_IPW", "dataset_id": "all", "status": "informative_sensitivity" if positivity and ess >= .7 * len(w) else "POSITIVITY_FAILURE",
        "n_patients": int(cov.patient_key.nunique()), "detail": f"OOF inclusion AUC={auc_rank(y, pred):.4f}; ESS={ess:.2f}/{len(w)}; max_weight={w.max() if len(w) else np.nan:.3f}; positivity={positivity}",
        "effective_sample_size": ess, "positivity_pass": positivity,
    }])
    return cov, summary


def negative_controls(anchor: pd.DataFrame, effects: pd.DataFrame, n_perm: int = 500) -> pd.DataFrame:
    primary = anchor[(anchor.primary_role.eq("conditional_primary")) & anchor.baseline_eligible & anchor.response_known].copy()
    primary = primary.groupby(["cohort_id", "patient_key", "response_binary"], as_index=False)[PRIMARY].mean()
    rng = np.random.default_rng(SEED)
    rows = []
    observed = effects[effects.statistical_estimability.eq("estimable")].groupby("barrier_id").effect.mean()
    null = {fm: [] for fm in PRIMARY}
    for _ in range(n_perm):
        shuffled = primary.copy()
        shuffled["response_binary"] = shuffled.groupby("cohort_id").response_binary.transform(lambda x: rng.permutation(x.to_numpy()))
        for fm in PRIMARY:
            vals = []
            for _, d in shuffled.groupby("cohort_id"):
                h = hedges_g(d, fm)
                if np.isfinite(h["effect"]): vals.append(h["effect"])
            null[fm].append(float(np.mean(vals)) if vals else np.nan)
    for fm in PRIMARY:
        arr = np.asarray(null[fm], float); obs = float(observed.get(fm, np.nan))
        p = float((1 + np.sum(np.abs(arr[np.isfinite(arr)]) >= abs(obs))) / (1 + np.isfinite(arr).sum())) if np.isfinite(obs) else np.nan
        rows.append({"control_id": "within_cohort_response_permutation", "barrier_id": fm, "observed_statistic": obs,
                     "null_mean": float(np.nanmean(arr)), "empirical_p": p, "n_permutations": n_perm,
                     "status": "negative_control_collapses" if p <= .10 else "negative_control_not_separated"})
    return pd.DataFrame(rows)


def measurement_error_sensitivity(anchor: pd.DataFrame, effects: pd.DataFrame) -> pd.DataFrame:
    rows = []
    primary = anchor[(anchor.primary_role.eq("conditional_primary")) & anchor.baseline_eligible & anchor.response_known]
    for cohort, d in primary.groupby("cohort_id"):
        for fm in PRIMARY:
            se_col = f"{fm}_measurement_se"
            finite = pd.to_numeric(d.get(se_col), errors="coerce").notna().mean() if se_col in d else 0.0
            score_var = float(pd.to_numeric(d[fm], errors="coerce").var(ddof=1))
            err_var = float(np.nanmean(pd.to_numeric(d.get(se_col), errors="coerce") ** 2)) if se_col in d else np.nan
            reliability = max(0.0, min(1.0, (score_var - err_var) / score_var)) if np.isfinite(score_var) and score_var > 0 and np.isfinite(err_var) else np.nan
            base = effects[(effects.cohort_id.eq(cohort)) & effects.barrier_id.eq(fm)]
            effect = float(base.effect.iloc[0]) if len(base) else np.nan
            corrected = effect / math.sqrt(reliability) if finite >= .8 and np.isfinite(reliability) and reliability > .1 else np.nan
            rows.append({"record_type": "measurement_error_sensitivity", "cohort_id": cohort, "barrier_id": fm,
                         "finite_measurement_se_fraction": finite, "reliability_ratio": reliability, "error_corrected_effect": corrected,
                         "can_trigger_route": False, "abstention_code": "" if np.isfinite(corrected) else "MEASUREMENT_ERROR_UNCALIBRATED"})
    return pd.DataFrame(rows)


def route_verdict(effects: pd.DataFrame, meta: pd.DataFrame, crossfit: pd.DataFrame, controls: pd.DataFrame, bridge: pd.DataFrame) -> dict:
    barrier_rows = []
    for fm in PRIMARY:
        eligible = effects[(effects.barrier_id.eq(fm)) & effects.statistical_estimability.eq("estimable")]
        eligible = eligible.copy()
        eligible["arm_family"] = np.where(eligible.treatment_arm.astype(str).str.contains("monotherapy", case=False), "monotherapy", "other")
        strata_cols = ["endpoint_type", "resolution", "arm_family"]
        strata_sizes = eligible.groupby(strata_cols).cohort_id.nunique().sort_values(ascending=False)
        selected_stratum = sorted(strata_sizes[strata_sizes.eq(strata_sizes.max())].index)[0] if len(strata_sizes) else ("unknown", "unknown", "unknown")
        selected_endpoint, selected_resolution, selected_arm = selected_stratum
        selected = eligible[
            eligible.endpoint_type.eq(selected_endpoint) & eligible.resolution.eq(selected_resolution) & eligible.arm_family.eq(selected_arm)
        ]
        endpoint_row = meta[(meta.barrier_id.eq(fm)) & meta.endpoint_type.eq(selected_endpoint) & meta.resolution.eq(selected_resolution)
                            & meta.arm_family.eq(selected_arm) & meta.record_type.eq("endpoint_meta")]
        n_independent = int(selected.cohort_id.nunique())
        direction_methods_ok = len(selected) > 0 and selected.effect.apply(np.sign).nunique() == 1 and selected.method_direction_agreement.all()
        same_direction = n_independent >= 2 and direction_methods_ok
        pi_excludes_zero = False
        if len(endpoint_row):
            lo, hi = endpoint_row.iloc[0].get("prediction_low"), endpoint_row.iloc[0].get("prediction_high")
            pi_excludes_zero = bool(np.isfinite(lo) and np.isfinite(hi) and lo * hi > 0)
        cf = crossfit[(crossfit.barrier_id.eq(fm)) & crossfit.cohort_id.isin(selected.cohort_id)]
        merged = selected[["cohort_id", "effect"]].merge(cf[["cohort_id", "effect"]], on="cohort_id", suffixes=("_frozen", "_crossfit"))
        crossfit_stable = len(merged) == len(selected) and len(merged) > 0 and (np.sign(merged.effect_frozen) == np.sign(merged.effect_crossfit)).all()
        control = controls[(controls.barrier_id.eq(fm)) & controls.control_id.eq("within_cohort_response_permutation")]
        negative_ok = bool(len(control) and control.empirical_p.iloc[0] <= .10)
        g301_complete = not bridge[(bridge.get("dataset_id") == "GSE301741") & bridge.status.astype(str).str.contains("BRIDGE_FAILED")].shape[0]
        # GSE301741 is an optional projection-support lane by contract. Its
        # bridge status is reported but cannot veto a replicated direct-mid R1.
        hcc = eligible[eligible.cancer_type.astype(str).str.lower().eq("hcc")].copy()
        hcc_direction_ok = (
            hcc.cohort_id.nunique() >= 2
            and hcc.effect.apply(np.sign).nunique() == 1
            and hcc.method_direction_agreement.fillna(False).all()
        )
        if same_direction and pi_excludes_zero and crossfit_stable and negative_ok:
            route = "R1"
        elif same_direction and crossfit_stable and negative_ok:
            route = "R1-lite"
        elif hcc_direction_ok and negative_ok:
            route = "R2"
        else:
            route = "Unresolved"
        if route == "Unresolved":
            reasons = []
            if n_independent < 2: reasons.append("ENDPOINT_UNREPLICATED")
            if n_independent >= 2 and not direction_methods_ok: reasons.append("DIRECTION_OR_METHOD_DISAGREEMENT")
            if not crossfit_stable: reasons.append("CROSSFIT_FLIP_OR_NOT_ESTIMABLE")
            if not negative_ok: reasons.append("NEGATIVE_CONTROL_NOT_SEPARATED")
            if not pi_excludes_zero: reasons.append("PREDICTION_INTERVAL_CROSSES_ZERO")
            if not g301_complete: reasons.append("GSE301741_PROJECTION_SUPPORT_BRIDGE_FAILED")
        else:
            reasons = []
        direction_source = hcc if route == "R2" else selected
        barrier_rows.append({"barrier_id": fm, "route": route, "direction": "responder_higher" if len(direction_source) and direction_source.effect.mean() > 0 else "nonresponder_higher" if len(direction_source) else "undetermined",
                             "selected_endpoint": selected_endpoint, "selected_resolution": selected_resolution, "selected_arm_family": selected_arm,
                             "reasons": reasons, "same_endpoint_cohorts": n_independent,
                             "hcc_independent_cohorts": int(hcc.cohort_id.nunique()),
                             "hcc_direction_ok": bool(hcc_direction_ok),
                             "gse301741_projection_bridge_pass": bool(g301_complete),
                             "crossfit_stable": bool(crossfit_stable), "negative_control_pass": bool(negative_ok)})
    eligible = [x for x in barrier_rows if x["route"] in {"R1", "R1-lite", "R2"}]
    verdict = "GO_TO_PHASE9_LIMITED" if eligible else "UNRESOLVED_PHASE9_BLOCKED"
    return {
        "phase": "v6.2.1_phase8_anchor_repair_and_statistical_strengthening",
        "verdict": verdict,
        "phase9_entry_allowed": bool(eligible),
        "supervised_SRB_allowed": False,
        "counterfactual_repair_allowed": False,
        "X_class_ranking_allowed": False,
        "barriers": barrier_rows,
        "calibration_status": "NOT_APPLICABLE_NO_PREDICTIVE_MODEL",
        "claim_boundary": "response_association_candidate_only; no causal resistance or repair claim",
    }


def report(anchor: pd.DataFrame, effects: pd.DataFrame, meta: pd.DataFrame, verdict: dict, audit120: dict, bridge: pd.DataFrame) -> str:
    strict = anchor[(anchor.primary_role.eq("conditional_primary")) & anchor.baseline_eligible & anchor.response_known & anchor[PRIMARY].notna().all(axis=1)]
    cohort = strict.groupby(["cohort_id", "endpoint_type", "statistical_estimability"]).response_binary.agg(n="size", R="sum").reset_index()
    cohort["NR"] = cohort.n - cohort.R
    lines = [
        "# Phase8 Anchor Repair and Statistical Strengthening Final Report", "", f"Created: {now_iso()}", "",
        "## Executive Verdict", f"**{verdict['verdict']}**", "",
        f"真实 strict baseline surface：{strict.patient_key.nunique()} patients，R={int((strict.response_binary == 1).sum())}，NR={int((strict.response_binary == 0).sum())}。",
        "", "## Anchor Surface", "", cohort.to_markdown(index=False), "",
        "## Data Repairs", "",
        "- GSE123813：公开 Supplementary Table 1 恢复 RECIST 1.1、PD1 药物和 ongoing vismodegib；后者只进组合治疗 sensitivity。",
        f"- GSE120575：患者级 RECIST 已恢复，但 `counts` 抽样非整数比例为 {audit120['sample_noninteger_fraction']:.4f}，因此降级为 normalized-expression support，不进入 strict direct-mid 主线。",
        "- GSE301741：严格重放冻结 Phase4A eligibility；未放宽 cell-state 规则，无法把 projection 伪装成 direct-mid。",
        "- GSE243013：改用显式 anti-PD1/chemotherapy/targeted-therapy 字段拆臂，不再按字符串误判 212 名患者为单药。",
        "", "## Statistical Strengthening", "",
        "- 主结果按 cohort × endpoint × treatment arm × resolution 分层；执行 Hedges g、Cliff's delta、Firth logistic。",
        "- 只在同 endpoint 至少两个可估环境时做 Hartung-Knapp random-effects meta；跨 endpoint 仅比较方向。",
        "- Bayesian partial pooling 为 prior-sensitivity，不可独立触发 route。",
        "- coverage 使用 response-blind OOF inclusion model 和 stabilized IPW；outer-cohort residualization 单独检查方向翻转。",
        "- uncertainty/abstention 使用结构化原因；未训练预测模型，因此 calibration 为 NOT_APPLICABLE。",
        "", "## Route Decision", "",
    ]
    for row in verdict["barriers"]:
        lines.append(f"- {row['barrier_id']}: {row['route']}; reasons={','.join(row['reasons']) or 'none'}")
    lines += ["", "## Boundaries", "", "- 未修改 FM01/FM04/FM07 membership。", "- 未训练 supervised SRB。", "- 未执行 counterfactual repair 或 X-class ranking。", "- Unresolved 不解释为生物学阴性。", ""]
    return "\n".join(lines)


def publish(staging: Path) -> None:
    if OUT.exists():
        archive = OUT.parent / f"{OUT.name}.previous"
        if archive.exists(): shutil.rmtree(archive)
        os.replace(OUT, archive)
    os.replace(staging, OUT)


def analysis_input_paths(cfg: dict) -> list[Path]:
    """Enumerate every direct and crossfit-indirect input consumed by the run."""
    from phase7_common import load_config as load_phase7_config

    cfg7 = load_phase7_config()
    paths = set(cfg["inputs"].values()) | set(cfg7["inputs"].values())
    paths.update({
        Path(os.environ.get("PHASE8_CONFIG_PATH", ROOT / "scripts/v6_2/phase8_v6_2_1_config.yaml")),
        Path(os.environ.get("PHASE7_CONFIG_PATH", ROOT / "scripts/v6_2/phase7_v6_2_1_config.yaml")),
        P7 / "measurement/module_activity_by_patient_timepoint_residualized.csv",
        P7 / "measurement/sample_module_coverage_gate.csv",
        P7 / "sentinel/prespecified_sentinel_program_scores.csv",
        PHASE4A / "handoff/phase4a_to_phase4b_aggregation_eligibility.csv",
        PHASE4B / "qc_covariates/sample_patient_timepoint_qc_covariates.csv",
        OLD / "hcc_sensitivity/gse206325_coarse_response_sensitivity.csv",
        GSE120_SUPP, GSE120_H5AD, GSE120_ADDENDUM, GSE123_SUPP, GSE243_META,
    })
    handoff = yaml.safe_load(cfg7["inputs"]["phase6_handoff"].read_text())
    for key in ("default_module_score_matrix", "default_module_membership", "module_dictionary", "module_method_support_map"):
        if handoff.get(key): paths.add(Path(handoff[key]))
    missing = [str(p) for p in paths if not p.exists()]
    if missing:
        raise FileNotFoundError(f"Analysis input manifest contains missing paths: {missing}")
    return sorted(paths)


def main() -> None:
    attempt = OUT.parent / f"{OUT.name}.last_attempt.yaml"
    write_yaml(attempt, {"status": "RUNNING", "started_at": now_iso(), "output": str(OUT)})
    def mark_interrupted() -> None:
        state = yaml.safe_load(attempt.read_text()) if attempt.exists() else {}
        if state.get("status") == "RUNNING":
            write_yaml(attempt, {**state, "status": "FAILED_OR_INTERRUPTED", "finished_at": now_iso()})
    atexit.register(mark_interrupted)
    if OUT.exists():
        previous = OUT.parent / f"{OUT.name}.previous"
        if previous.exists(): shutil.rmtree(previous)
        os.replace(OUT, previous)
    if STAGING.exists(): shutil.rmtree(STAGING)
    STAGING.mkdir(parents=True)
    cfg = load_config()
    anchor, bridge, audit120 = build_anchor_master(cfg)
    loss = gse301_loss_trace(cfg)
    cov, cov_summary = coverage_and_ipw(cfg, anchor)
    bridge = pd.concat([bridge, loss, cov_summary], ignore_index=True, sort=False)
    g206 = pd.read_csv(OLD / "hcc_sensitivity/gse206325_coarse_response_sensitivity.csv")
    for _, x in g206.iterrows():
        bridge = pd.concat([bridge, pd.DataFrame([{"record_type": "GSE206325_coarse_HCC_sensitivity", "dataset_id": "GSE206325", "barrier_id": x.get("barrier_id"), "status": "coarse_post_treatment_support_only", "detail": f"effect={x.get('effect')}; baseline unavailable"}])], ignore_index=True, sort=False)

    effects = within_effects(anchor)
    meta = endpoint_meta(effects)
    bayes = bayes_sensitivity(effects, cfg)
    estimable = sorted(effects.loc[effects.statistical_estimability.eq("estimable"), "cohort_id"].unique())
    crossfit = crossfit_effects(cfg, anchor, estimable)
    measurement = measurement_error_sensitivity(anchor, effects)
    effects_master = pd.concat([effects, meta, bayes, crossfit, measurement], ignore_index=True, sort=False)
    controls = negative_controls(anchor, effects)
    sensitivity = sensitivity_and_sentinels(anchor)
    verdict = route_verdict(effects, meta, crossfit, controls, bridge)

    anchor.to_parquet(STAGING / "phase8_anchor_master.parquet", index=False)
    effects_master.to_csv(STAGING / "phase8_effects_master.csv", index=False)
    bridge.to_csv(STAGING / "phase8_coverage_and_bridge.csv", index=False)
    sensitivity.to_csv(STAGING / "phase8_sensitivity_and_sentinel.csv", index=False)
    controls.to_csv(STAGING / "phase8_negative_controls.csv", index=False)
    write_yaml(STAGING / "phase8_route_gate_and_handoff.yaml", verdict)
    (STAGING / "PHASE8_FINAL_REPORT.md").write_text(report(anchor, effects, meta, verdict, audit120, bridge), encoding="utf-8")

    manifest = {
        "phase": "v6.2.1_phase8_anchor_repair_and_statistical_strengthening",
        "created_at": now_iso(), "seed": SEED, "run_status": "COMPLETE", "verdict": verdict["verdict"],
        "primary_barriers": PRIMARY,
        "input_hashes": {str(p.relative_to(ROOT)): sha256(p) for p in analysis_input_paths(cfg)},
        "runner_sha256": sha256(Path(__file__)),
        "rules": {"module_membership_modified": False, "response_used_for_score_construction": False, "endpoint_pooling_for_claim": False,
                  "projection_as_direct_mid": False, "raw_fastq_downloaded": False},
        "outputs": {},
    }
    for p in sorted(STAGING.iterdir()):
        if p.name == "phase8_contract_and_run_manifest.yaml": continue
        manifest["outputs"][p.name] = {"sha256": sha256(p), "size_bytes": p.stat().st_size}
    write_yaml(STAGING / "phase8_contract_and_run_manifest.yaml", manifest)
    publish(STAGING)
    write_yaml(attempt, {"status": "COMPLETE", "started_at": yaml.safe_load(attempt.read_text())["started_at"], "finished_at": now_iso(), "output": str(OUT), "verdict": verdict["verdict"]})
    print(yaml.safe_dump({"output": str(OUT), "verdict": verdict["verdict"], "n_outputs": len(list(OUT.iterdir()))}, sort_keys=False))


if __name__ == "__main__":
    main()
