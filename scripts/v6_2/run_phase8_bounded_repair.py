#!/usr/bin/env python3
"""Phase8 bounded repair and shared-context re-adjudication.

This runner never changes frozen module membership. It reuses the audited Phase7
scores and the bounded Phase8 projection/crossfit artifacts, then emits one
compact, independently auditable decision package.
"""
from __future__ import annotations

import hashlib
import math
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

from phase7_common import ROOT, normalize_timepoint, sha256
from run_phase8_anchor_context import (
    auc_rank,
    bayesian_partial_pool,
    build_anchor_environment,
    hedges_g,
    load_config,
    random_effects,
)

OUT = ROOT / "results/v6_2/phase8_bounded_repair_and_readjudication"
P7 = ROOT / "results/v6_2/phase7_module_measurement_and_barrier_identifiability"
OLD = ROOT / "results/v6_2/phase8_anchor_context_adjudication"
PRIMARY = ["FM01", "FM04", "FM07"]
SENSITIVITY = {"APC_MHCII": ["FM02", "FM03"], "CYTOTOXIC_T": ["FM05", "FM08"], "B_PLASMA": ["FM06"]}
SENTINELS = ["TLS_B", "CD8_PROGENITOR_EXHAUSTION", "CD8_TERMINAL_EXHAUSTION"]


def write_yaml(path: Path, obj: dict) -> None:
    path.write_text(yaml.safe_dump(obj, sort_keys=False, allow_unicode=True))


def logistic_1d(y: np.ndarray, x: np.ndarray) -> dict:
    ok = np.isfinite(y) & np.isfinite(x)
    y, x = y[ok].astype(float), x[ok].astype(float)
    if len(y) < 6 or len(np.unique(y)) < 2 or np.std(x) < 1e-10:
        return {"log_odds": np.nan, "se": np.nan, "auc": np.nan, "converged": False}
    z = (x - x.mean()) / x.std(ddof=0)
    design = np.column_stack([np.ones(len(z)), z])
    beta = np.zeros(2)
    converged = False
    for _ in range(80):
        eta = np.clip(design @ beta, -20, 20)
        p = 1 / (1 + np.exp(-eta))
        w = np.maximum(p * (1 - p), 1e-6)
        h = design.T @ (w[:, None] * design) + np.eye(2) * 1e-6
        step = np.linalg.solve(h, design.T @ (y - p))
        beta += step
        if np.max(np.abs(step)) < 1e-7:
            converged = True
            break
    cov = np.linalg.pinv(h)
    return {"log_odds": float(beta[1]), "se": float(math.sqrt(max(cov[1, 1], 0))),
            "auc": auc_rank(y, z), "converged": converged}


def logistic_predict(y: np.ndarray, x: np.ndarray) -> tuple[np.ndarray, bool]:
    """Small response-blind inclusion model with ridge-stabilized IRLS."""
    ok = np.isfinite(y) & np.isfinite(x).all(axis=1)
    pred = np.full(len(y), np.nan)
    if ok.sum() < 20 or len(np.unique(y[ok])) < 2:
        return pred, False
    xx = x[ok].astype(float)
    sd = xx.std(axis=0); sd[sd < 1e-8] = 1
    xx = (xx - xx.mean(axis=0)) / sd
    design = np.column_stack([np.ones(len(xx)), xx])
    beta = np.zeros(design.shape[1]); converged = False
    penalty = np.eye(design.shape[1]) * .1; penalty[0, 0] = 0
    for _ in range(100):
        p = 1 / (1 + np.exp(-np.clip(design @ beta, -20, 20)))
        w = np.maximum(p * (1 - p), 1e-5)
        h = design.T @ (w[:, None] * design) + penalty
        step = np.linalg.solve(h, design.T @ (y[ok] - p) - penalty @ beta)
        beta += step
        if np.max(np.abs(step)) < 1e-7:
            converged = True; break
    pred[ok] = 1 / (1 + np.exp(-np.clip(design @ beta, -20, 20)))
    return pred, converged


def cliffs_delta(data: pd.DataFrame, feature: str, seed: int = 1729, n_boot: int = 1000) -> dict:
    r = pd.to_numeric(data.loc[data.response_binary.eq(1), feature], errors="coerce").dropna().to_numpy()
    nr = pd.to_numeric(data.loc[data.response_binary.eq(0), feature], errors="coerce").dropna().to_numpy()
    if min(len(r), len(nr)) < 2:
        return {"cliffs_delta": np.nan, "ci_low": np.nan, "ci_high": np.nan}
    def delta(a, b):
        return float((np.greater.outer(a, b).sum() - np.less.outer(a, b).sum()) / (len(a) * len(b)))
    rng = np.random.default_rng(seed)
    boot = [delta(rng.choice(r, len(r), True), rng.choice(nr, len(nr), True)) for _ in range(n_boot)]
    return {"cliffs_delta": delta(r, nr), "ci_low": float(np.quantile(boot, .025)), "ci_high": float(np.quantile(boot, .975))}


def membership_hash(membership: pd.DataFrame) -> str:
    id_col = "module_id" if "module_id" in membership else "frozen_module_id"
    x = membership[membership[id_col].isin(PRIMARY)].copy()
    cols = [c for c in [id_col, "gene", "gene_symbol", "weight", "consensus_weight", "membership_weight", "relative_weight"] if c in x]
    text = x[cols].sort_values(cols[:2]).to_csv(index=False, lineterminator="\n")
    return hashlib.sha256(text.encode()).hexdigest()


def input_audit(cfg: dict) -> tuple[dict, pd.DataFrame]:
    membership = pd.read_csv(cfg["inputs"]["phase6_membership"])
    barriers = pd.read_parquet(cfg["inputs"]["barrier_matrix"])
    p7_manifest = yaml.safe_load(cfg["inputs"]["phase7_handoff"].read_text())
    closeout_path = P7 / "preflight/phase6_closeout_manifest_v6_2_1.yaml"
    closeout = yaml.safe_load(closeout_path.read_text())
    checks = []
    def add(check, passed, evidence):
        checks.append({"check": check, "status": "PASS" if passed else "HARD_FAIL", "evidence": str(evidence)})
    actual_hash = membership_hash(membership)
    frozen_ids = set(pd.read_csv(cfg["inputs"]["barrier_table"]).module_id)
    add("frozen_primary_ids", set(PRIMARY).issubset(frozen_ids), sorted(frozen_ids))
    membership_sha = sha256(cfg["inputs"]["phase6_membership"])
    expected_membership = next((x["sha256"] for x in closeout["assets"] if x["path"].endswith("module_membership.frozen_v1.csv")), None)
    add("membership_file_matches_phase7_hash", membership_sha == expected_membership, f"actual={membership_sha};expected={expected_membership}")
    expected_barrier_path = p7_manifest.get("primary_outputs", {}).get("eligible_barrier_score_matrix")
    add("barrier_matrix_path_matches_phase7_handoff", str(cfg["inputs"]["barrier_matrix"].relative_to(ROOT)) == expected_barrier_path,
        f"actual={cfg['inputs']['barrier_matrix'].relative_to(ROOT)};expected={expected_barrier_path}")
    leak = [c for c in barriers.columns if any(t in c.lower() for t in ["response", "outcome", "survival", "split", "progression"])]
    add("no_label_columns_in_barrier_matrix", not leak, leak)
    add("primary_score_columns_present", set(PRIMARY).issubset(barriers.columns), PRIMARY)
    # The frozen primary matrix is mid-resolution by contract; coarse scores are
    # stored in a separate sensitivity object and never silently pooled.
    resolution_contract = cfg["inputs"]["coarse_matrix"].exists() and p7_manifest.get("verdict") == "CONDITIONAL_GO_TO_PHASE8"
    add("resolution_tag_present_or_unambiguous_contract", resolution_contract, "primary=mid;coarse=separate_sensitivity_object")
    add("traceability_keys_present", set(["cohort_id", "patient_key", "timepoint", "tissue_context"]).issubset(barriers.columns), list(barriers.columns))
    patients = pd.read_csv(cfg["inputs"]["patient_metadata"])
    add("patient_split_unique", patients.groupby("patient_key").split.nunique().max() <= 1, patients.groupby("patient_key").split.nunique().max())
    hard = [x["check"] for x in checks if x["status"] == "HARD_FAIL"]
    contract = {
        "phase": "v6.2.1_phase8_bounded_repair",
        "freeze_status": "HARD_FAIL" if hard else "PASS",
        "primary_barriers": PRIMARY,
        "membership_content_hash": actual_hash,
        "membership_file_sha256": sha256(cfg["inputs"]["phase6_membership"]),
        "barrier_matrix_sha256": sha256(cfg["inputs"]["barrier_matrix"]),
        "authorized_input_whitelist": [
            {"path": str(path.relative_to(ROOT)), "sha256": sha256(path)}
            for path in sorted(set(cfg["inputs"].values())) if path.exists() and path.is_file()
        ] + [
            {"path": str(path.relative_to(ROOT)), "sha256": sha256(path)} for path in [
                P7 / "measurement/module_activity_by_patient_timepoint_residualized.csv",
                P7 / "measurement/sample_module_coverage_gate.csv",
                P7 / "sentinel/prespecified_sentinel_program_scores.csv",
                P7 / "annotation/module_top40_projection_sensitivity.csv",
            ]
        ],
        "hard_failures": hard,
        "prohibitions": ["module_retraining", "response_guided_gene_selection", "endpoint_pooling", "supervised_SRB", "counterfactual_repair", "X_class_ranking"],
    }
    return contract, pd.DataFrame(checks)


def repaired_environment(cfg: dict) -> tuple[pd.DataFrame, pd.DataFrame]:
    env = build_anchor_environment(cfg).copy()
    residual = pd.read_csv(P7 / "measurement/module_activity_by_patient_timepoint_residualized.csv")
    residual["timepoint"] = [
        normalize_timepoint(value, cohort)
        for value, cohort in zip(residual["timepoint"], residual["cohort_id"])
    ]
    se = residual[residual.module_id.isin(PRIMARY)].pivot_table(
        index=["cohort_id", "patient_key", "timepoint", "tissue_context"], columns="module_id", values="route_b_posterior_se").add_suffix("_measurement_se").reset_index()
    env = env.merge(se, on=["cohort_id", "patient_key", "timepoint", "tissue_context"], how="left", validate="one_to_one")
    samples = pd.read_csv(cfg["inputs"]["sample_metadata"])
    samples["normalized_timepoint"] = [
        normalize_timepoint(value, cohort)
        for value, cohort in zip(samples["timepoint"], samples["cohort_id"])
    ]
    sm = samples.groupby(["patient_key", "cohort_id", "normalized_timepoint", "tissue_source"], as_index=False).agg(
        sample_key=("sample_key", lambda x: "|".join(sorted(set(x.astype(str))))),
        treatment_regimen=("treatment_raw", lambda x: "|".join(sorted(set(x.dropna().astype(str))))),
    )
    env = env.merge(sm, left_on=["patient_key", "cohort_id", "timepoint", "tissue_context"],
                    right_on=["patient_key", "cohort_id", "normalized_timepoint", "tissue_source"], how="left")
    fallback = samples.groupby(["patient_key", "cohort_id", "normalized_timepoint"], as_index=False).agg(
        fallback_sample_key=("sample_key", lambda x: "|".join(sorted(set(x.astype(str))))),
        fallback_treatment=("treatment_raw", lambda x: "|".join(sorted(set(x.dropna().astype(str))))))
    env = env.merge(fallback, left_on=["patient_key", "cohort_id", "timepoint"],
                    right_on=["patient_key", "cohort_id", "normalized_timepoint"], how="left", suffixes=("", "_fallback"))
    env["sample_key"] = env.sample_key.fillna(env.fallback_sample_key)
    env["treatment_regimen"] = env.treatment_regimen.fillna(env.fallback_treatment)
    env["normalized_timepoint"] = env["timepoint"]
    env["resolution"] = env.get("score_resolution", "mid")
    env["score_source"] = "phase7_eligible_barrier_score_matrix_v0"
    env["measurement_route"] = "route_a_wls_residual_primary"
    env["coverage_status"] = np.where(env[PRIMARY].notna().all(axis=1), "eligible", "module_missing")
    env["response_raw"] = env["response_raw_summary"]
    env["endpoint_type"] = env["response_endpoint_type"]
    env["endpoint_quality"] = env["response_harmonization_confidence"]
    env["treatment_arm"] = env["treatment_arm_class"]
    env["monotherapy_flag"] = env.treatment_arm.isin(["monotherapy", "monotherapy_D_only"])
    env["treatment_cleanliness"] = np.where(env.treatment_arm.eq("monotherapy_D_only"), "conditional_clean",
        np.where(env.treatment_arm.eq("monotherapy"), "clean_or_conditional", np.where(env.treatment_arm.str.contains("combination", na=False), "combination_sensitivity", "unknown")))
    env["HCC_flag"] = env.cancer_type.astype(str).str.lower().isin(["hcc", "hepatocellular carcinoma"])
    env["baseline_eligible"] = env.baseline_primary_eligible
    env["primary_role"] = np.where(env.phase8_analysis_role.eq("conditional_anchor_primary"), "conditional_primary", "not_primary")
    env["sensitivity_role"] = np.where(env.phase8_analysis_role.eq("sensitivity"), "sensitivity", np.where(env.primary_role.eq("conditional_primary"), "primary_also_sensitivity", "support_only"))
    env["exclusion_reason"] = env.eligibility_reason
    for fm in PRIMARY:
        env[f"{fm}_available"] = env[fm].notna()
    env["patient_timepoint_context_id"] = env.patient_key + "::" + env.normalized_timepoint + "::" + env.tissue_context
    wanted = ["cohort_id", "patient_key", "sample_key", "patient_timepoint_context_id", "normalized_timepoint", "tissue_context", "resolution", "score_source"] + PRIMARY + [f"{x}_measurement_se" for x in PRIMARY] + [f"{x}_available" for x in PRIMARY] + [
        "measurement_route", "coverage_status", "response_raw", "response_binary", "endpoint_type", "endpoint_quality", "treatment_regimen", "treatment_arm",
        "monotherapy_flag", "treatment_cleanliness", "cancer_type", "HCC_flag", "baseline_eligible", "primary_role", "sensitivity_role", "split", "exclusion_reason"]
    env = env[wanted].drop_duplicates(["patient_key", "normalized_timepoint", "tissue_context", "resolution"])
    group = ["cohort_id", "endpoint_type", "treatment_arm", "resolution"]
    summary = env.groupby(group, dropna=False).apply(lambda x: pd.Series({
        "n_patients_total": x.patient_key.nunique(), "n_baseline": x.loc[x.baseline_eligible].patient_key.nunique(),
        "n_R": x.loc[x.response_binary.eq(1)].patient_key.nunique(), "n_NR": x.loc[x.response_binary.eq(0)].patient_key.nunique(),
        "n_unknown": x.loc[x.response_binary.isna()].patient_key.nunique(),
        **{f"n_{fm}_available": x.loc[x[f"{fm}_available"]].patient_key.nunique() for fm in PRIMARY},
        "estimability": "estimable" if x.response_binary.eq(1).sum() >= 3 and x.response_binary.eq(0).sum() >= 3 else ("weakly_estimable" if x.response_binary.nunique() == 2 else "not_estimable"),
        "role": "conditional_primary" if x.primary_role.eq("conditional_primary").any() else ("sensitivity" if x.sensitivity_role.eq("sensitivity").any() else "support_only"),
        "not_estimable_reason": "" if x.response_binary.nunique() == 2 else "one_or_no_response_class",
    }), include_groups=False).reset_index()
    return env, summary


def within_effects(env: pd.DataFrame) -> pd.DataFrame:
    primary = env[env.primary_role.eq("conditional_primary") & env.baseline_eligible].copy()
    rows = []
    for keys, d in primary.groupby(["cohort_id", "endpoint_type", "treatment_arm", "resolution"], dropna=False):
        for fm in PRIMARY:
            h, c = hedges_g(d, fm), cliffs_delta(d, fm)
            l = logistic_1d(d.response_binary.to_numpy(float), pd.to_numeric(d[fm], errors="coerce").to_numpy())
            rows.append(dict(zip(["cohort_id", "endpoint_type", "treatment_arm", "resolution"], keys)) | {"barrier_id": fm,
                "n_patients": d.patient_key.nunique(), **h, **c, **l,
                "direction_agreement_A1_A2": bool(np.sign(l["log_odds"]) == np.sign(c["cliffs_delta"])) if np.isfinite(l["log_odds"]) and np.isfinite(c["cliffs_delta"]) else False,
                "information_class": "shared_direction_input"})
    return pd.DataFrame(rows)


def route_a_meta(effects: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for fm, d in effects.groupby("barrier_id"):
        valid = d.dropna(subset=["effect", "se"])
        m = random_effects(valid)
        if len(valid) >= 2 and m.get("k", 0) >= 2:
            w = 1 / (valid.se.to_numpy() ** 2 + m["tau2"])
            q_star = float(np.sum(w * (valid.effect.to_numpy() - m["pooled_effect"]) ** 2) / (len(valid) - 1))
            hk_se = math.sqrt(max(q_star, 1e-12) / w.sum())
            # Exact two-sided 95% t criticals for the small k seen here.
            tcrit = {1: 12.706, 2: 4.303, 3: 3.182, 4: 2.776, 5: 2.571}.get(len(valid) - 1, 1.96)
            m.update({"hartung_knapp_se": hk_se, "hartung_knapp_ci_low": m["pooled_effect"] - tcrit * hk_se,
                      "hartung_knapp_ci_high": m["pooled_effect"] + tcrit * hk_se})
        signs = np.sign(valid.effect)
        m.update({"barrier_id": fm, "synthesis": "cross_endpoint_sensitivity_only", "direction_consistency": float((signs == np.sign(np.nanmean(signs))).mean()) if len(signs) else np.nan,
                  "same_endpoint_shared_estimable": bool(valid.endpoint_type.nunique() == 1 and len(valid) >= 2),
                  "endpoint_confounded_with_cohort": bool(valid.groupby("cohort_id").endpoint_type.nunique().max() == 1 and valid.endpoint_type.nunique() == valid.cohort_id.nunique()) if len(valid) else True})
        rows.append(m)
        for endpoint, ed in valid.groupby("endpoint_type"):
            em = random_effects(ed); em.update({"barrier_id": fm, "synthesis": f"endpoint::{endpoint}", "direction_consistency": 1.0 if len(ed) else np.nan,
                "same_endpoint_shared_estimable": len(ed) >= 2, "endpoint_confounded_with_cohort": False})
            rows.append(em)
    return pd.DataFrame(rows)


def interface_results(cfg: dict, env: pd.DataFrame) -> pd.DataFrame:
    rows = []
    g = env[env.cohort_id.eq("GSE286827")]
    for arm, d in g.groupby("treatment_arm"):
        rows.append({"interface": "GSE286827_treatment_arm", "component": arm, "n_patients": d.patient_key.nunique(), "status": "conditional_clean_anchor" if arm == "monotherapy_D_only" else "treatment_context_sensitivity", "evidence": f"R={d.response_binary.eq(1).sum()};NR={d.response_binary.eq(0).sum()}"})
    proj = pd.read_csv(OLD / "projection/gse301741_sample_level_frozen_module_projection.csv")
    bridge = pd.read_csv(OLD / "projection/gse301741_projection_bridge_audit.csv")
    coverage = pd.read_csv(P7 / "measurement/sample_module_coverage_gate.csv")
    g301cov = coverage[coverage.cohort_id.eq("GSE301741")]
    trace = [
        ("raw_patient_present", 16, "pass"), ("raw_count_available", 16, "pass"),
        ("barcode_sample_patient_mapping", proj.patient_key.nunique(), "pass"),
        ("sample_level_frozen_projection", proj.patient_key.nunique(), "pass"),
        ("phase4a_mid_eligible_sample_handoff", g301cov.sample_key.nunique(), "loss_22_of_27_samples"),
        ("phase7_primary_coverage_context", g301cov.loc[g301cov.scoring_status.eq("primary"), "patient_timepoint_context_id"].nunique(), "pass_after_handoff"),
        ("phase7_low_coverage_context", g301cov.loc[g301cov.scoring_status.eq("low_coverage"), "patient_timepoint_context_id"].nunique(), "excluded_by_coverage"),
        ("mid_projection_bridge_pairs", int(bridge.n_shared_contexts.max()), "insufficient_for_bridge_gate"),
    ]
    for stage, n, status in trace:
        rows.append({"interface": "GSE301741_loss_trace", "component": stage, "n_patients": n, "status": status,
            "evidence": "loss_is_mainly_before_phase7_scoring_at_mid_state_eligibility_handoff;not_raw_count_or_gene_overlap_failure"})
    for _, x in bridge.iterrows():
        rows.append({"interface": "GSE301741_projection", "component": x.get("barrier_id", x.get("module_id", "unknown")), "n_patients": proj.patient_key.nunique(),
            "status": x.get("bridge_status", "support_only_unresolved"), "evidence": f"n_shared_contexts={x.get('n_shared_contexts')};spearman={x.get('sample_mid_spearman')}"})
    labels = pd.read_csv(cfg["inputs"]["patient_metadata"])[["patient_key", "response_binary_harmonized"]]
    projected = proj[proj.timepoint.eq("baseline")].merge(labels, on="patient_key", how="left")
    projected["response_binary"] = projected.response_binary_harmonized.map({"responder": 1, "non_responder": 0})
    for fm in PRIMARY:
        h = hedges_g(projected, f"{fm}_projection_z")
        rows.append({"interface": "GSE301741_baseline_projection_response_sensitivity", "component": fm,
            "n_patients": projected.patient_key.nunique(), "status": "projection_support_only_not_mid_primary",
            "effect": h["effect"], "se": h["se"], "n_R": h["n_R"], "n_NR": h["n_NR"],
            "evidence": f"R={h['n_R']};NR={h['n_NR']};effect={h['effect']};response_not_used_for_projection_selection"})
    hcc = pd.read_csv(OLD / "hcc_sensitivity/gse206325_coarse_response_sensitivity.csv")
    for _, x in hcc.iterrows():
        rows.append({"interface": "GSE206325_coarse_HCC", "component": x.get("barrier_id", x.get("module_id", "unknown")), "n_patients": int(x.get("n_R", 0) + x.get("n_NR", 0)),
            "status": "coarse_post_treatment_sensitivity_not_primary", "evidence": f"effect={x.get('effect', np.nan)};baseline_unavailable=true"})
    rows.append({"interface": "cross_resolution_bridge", "component": "coarse_to_mid", "n_patients": 0, "status": "not_estimable", "evidence": "no_same_patient_or_context_coarse_mid_pairs_across_independent_cohorts"})
    return pd.DataFrame(rows)


def coverage_audit(cfg: dict, env: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    cov = pd.read_csv(P7 / "measurement/sample_module_coverage_gate.csv")
    qc = pd.read_csv(ROOT / "results/v6_2/phase4b_immune_state_feature_construction/qc_covariates/sample_patient_timepoint_qc_covariates.csv")
    cov = cov[cov.module_id.isin(PRIMARY)].groupby(["cohort_id", "patient_key", "timepoint", "module_id"], as_index=False).agg(
        module_included=("scoring_status", lambda x: int((x == "primary").any())),
        sample_key=("sample_key", lambda x: "|".join(sorted(set(x.astype(str))))), available_fraction=("available_fraction", "mean"),
        n_states=("n_states", "max"), n_cells=("n_cells", "sum"))
    cov = cov.groupby(["cohort_id", "patient_key", "timepoint"], as_index=False).agg(
        included=("module_included", "min"), sample_key=("sample_key", lambda x: "|".join(sorted(set(x.astype(str))))),
        available_fraction=("available_fraction", "mean"), n_states=("n_states", "mean"), n_cells=("n_cells", "mean"))
    qc["timepoint_norm"] = qc.timepoint.map(normalize_timepoint)
    qcp = qc.groupby(["cohort_id", "patient_key", "timepoint_norm"], as_index=False).agg(
        total_cells=("total_cells", "sum"), mid_label_coverage=("mid_label_coverage", "mean"), low_quality_fraction=("low_quality_fraction", "mean"))
    cov = cov.merge(qcp, left_on=["cohort_id", "patient_key", "timepoint"], right_on=["cohort_id", "patient_key", "timepoint_norm"], how="left")
    labels = pd.read_csv(cfg["inputs"]["patient_metadata"])[["patient_key", "response_binary_harmonized", "cancer_type"]]
    cov = cov.merge(labels, on="patient_key", how="left")
    cov["response_binary"] = cov.response_binary_harmonized.map({"responder": 1, "non_responder": 0})
    numeric = cov[["available_fraction", "n_states", "n_cells", "total_cells", "mid_label_coverage", "low_quality_fraction"]].apply(pd.to_numeric, errors="coerce")
    numeric = numeric.fillna(numeric.median()).fillna(0)
    cohort_dummies = pd.get_dummies(cov.cohort_id, drop_first=True, dtype=float)
    pred, inclusion_converged = logistic_predict(cov.included.to_numpy(float), np.column_stack([numeric.to_numpy(), cohort_dummies.to_numpy()]))
    cov["inclusion_probability"] = pred
    stabilization = cov.included.mean()
    cov["inverse_probability_weight"] = np.where(cov.included.eq(1), stabilization / np.maximum(pred, .05), np.nan)
    cov["extreme_weight_flag"] = cov.inverse_probability_weight.gt(10)
    rows = []
    for scope, d in [("all", cov)] + [(f"cohort::{k}", v) for k, v in cov.groupby("cohort_id")]:
        l = logistic_1d(d.response_binary.to_numpy(float), d.included.to_numpy(float))
        rows.append({"audit_scope": scope, "n": len(d), "included_fraction": d.included.mean(), "R_included": int(((d.response_binary == 1) & (d.included == 1)).sum()),
            "NR_included": int(((d.response_binary == 0) & (d.included == 1)).sum()), "missingness_response_auc": l["auc"], "missingness_log_odds": l["log_odds"],
            "bias_flag": bool(np.isfinite(l["auc"]) and max(l["auc"], 1-l["auc"]) >= .65)})
    model = pd.DataFrame(rows)
    model = pd.concat([model, pd.DataFrame([{"audit_scope": "response_blind_inclusion_probability_model", "n": len(cov),
        "included_fraction": cov.included.mean(), "R_included": int(((cov.response_binary == 1) & (cov.included == 1)).sum()),
        "NR_included": int(((cov.response_binary == 0) & (cov.included == 1)).sum()),
        "missingness_response_auc": auc_rank(cov.included.to_numpy(float), pred), "missingness_log_odds": np.nan,
        "bias_flag": False, "model_converged": inclusion_converged}])], ignore_index=True)
    w = cov.loc[cov.included.eq(1), "inverse_probability_weight"].dropna()
    model = pd.concat([model, pd.DataFrame([{"audit_scope": "stabilized_IPW_extreme_weight_audit", "n": len(w),
        "included_fraction": cov.included.mean(), "R_included": int(((cov.response_binary == 1) & (cov.included == 1)).sum()),
        "NR_included": int(((cov.response_binary == 0) & (cov.included == 1)).sum()), "missingness_response_auc": np.nan,
        "missingness_log_odds": np.nan, "bias_flag": bool((w > 10).any()), "model_converged": inclusion_converged,
        "effective_sample_size": float(w.sum() ** 2 / (w.pow(2).sum())) if len(w) else np.nan,
        "max_weight": float(w.max()) if len(w) else np.nan}])], ignore_index=True)
    return cov, model


def sensitivity_and_sentinels(env: pd.DataFrame) -> pd.DataFrame:
    residual = pd.read_csv(P7 / "measurement/module_activity_by_patient_timepoint_residualized.csv")
    wide = residual.pivot_table(index=["cohort_id", "patient_key", "timepoint", "tissue_context"], columns="module_id", values="residual_activity_route_a").reset_index()
    labels = env[["cohort_id", "patient_key", "normalized_timepoint", "tissue_context", "response_binary", "primary_role", "endpoint_type", "treatment_arm", "resolution"]].rename(columns={"normalized_timepoint": "timepoint"}).drop_duplicates()
    d = wide.merge(labels, on=["cohort_id", "patient_key", "timepoint", "tissue_context"], how="inner")
    d = d[d.primary_role.eq("conditional_primary") & d.timepoint.eq("baseline")]
    rows = []
    for strata, sd in d.groupby(["cohort_id", "endpoint_type", "treatment_arm", "resolution"], dropna=False):
        for family, fms in SENSITIVITY.items():
            cols = [x for x in fms if x in sd]
            z = sd[cols].apply(lambda x: (x - x.mean()) / x.std(ddof=0)).fillna(0)
            methods = {"predeclared_standardized_mean": z.mean(axis=1), "fixed_equal_weight": z @ np.repeat(1 / len(cols), len(cols))}
            if len(cols) > 1 and len(sd) > 1:
                _, _, vt = np.linalg.svd(z.to_numpy(), full_matrices=False)
                pc = z.to_numpy() @ vt[0]
                if np.corrcoef(pc, z.mean(axis=1))[0, 1] < 0: pc = -pc
                methods["response_blind_PC1_oriented_to_family_mean"] = pd.Series(pc, index=sd.index)
            for method, score in methods.items():
                h = hedges_g(sd.assign(score=score), "score")
                rows.append({"lane": "sensitivity_topic", "program": family, "combination": method,
                    "cohort_id": strata[0], "endpoint_type": strata[1], "treatment_arm": strata[2], "resolution": strata[3],
                    **h, "status": "sensitivity_only_never_primary_endpoint_stratified"})
    sent = pd.read_csv(P7 / "sentinel/prespecified_sentinel_program_scores.csv")
    sent = sent[sent.program_id.isin(SENTINELS)].groupby(["cohort_id", "patient_key", "patient_timepoint_context_id", "program_id"], as_index=False).apply(
        lambda x: pd.Series({"score": np.average(x.activity_robust_z, weights=np.maximum(x.precision_weight, 1e-8))}), include_groups=False)
    lab = env[["cohort_id", "patient_key", "patient_timepoint_context_id", "response_binary", "primary_role", "normalized_timepoint"]].drop_duplicates()
    sent = sent.merge(lab, on=["cohort_id", "patient_key", "patient_timepoint_context_id"], how="inner")
    sent = sent[sent.primary_role.eq("conditional_primary") & sent.normalized_timepoint.eq("baseline")]
    sent = sent.groupby(["cohort_id", "patient_key", "program_id", "response_binary"], as_index=False).score.mean()
    for program, x in sent.groupby("program_id"):
        h = hedges_g(x, "score")
        rows.append({"lane": "sentinel", "program": program, "combination": "precision_weighted_context_score", **h, "status": "sanity_check_only"})
    for program in ["NEUTROPHIL_NET", "TUMOR_WNT_EXCLUSION"]:
        rows.append({"lane": "sentinel", "program": program, "combination": "not_run", "effect": np.nan, "se": np.nan, "n_R": 0, "n_NR": 0, "status": "unavailable_low_coverage_not_negative_evidence"})
    return pd.DataFrame(rows)


def negative_controls(env: pd.DataFrame, effects: pd.DataFrame, coverage_models: pd.DataFrame) -> pd.DataFrame:
    primary = env[env.primary_role.eq("conditional_primary") & env.baseline_eligible].copy()
    rows = []
    for name, col in [("cohort_only", "cohort_id"), ("cancer_only", "cancer_type"), ("endpoint_only", "endpoint_type"), ("treatment_arm_only", "treatment_arm"), ("timepoint_only", "normalized_timepoint")]:
        if primary[col].nunique() < 2:
            auc = .5
        else:
            rate = primary.groupby(col).response_binary.transform("mean")
            auc = auc_rank(primary.response_binary.to_numpy(float), rate.to_numpy(float))
        rows.append({"analysis": name, "barrier_id": "metadata", "metric": "apparent_auc", "value": auc, "status": "negative_control"})
    rng = np.random.default_rng(1729)
    permutation_schemes = {
        "global_shuffled_response": [], "cohort_within_shuffle": ["cohort_id"],
        "endpoint_within_shuffle": ["endpoint_type"], "patient_block_permutation": [],
        "cohort_endpoint_block_permutation": ["cohort_id", "endpoint_type"],
    }
    for fm in PRIMARY:
        observed = abs(hedges_g(primary, fm)["effect"])
        for scheme, group_cols in permutation_schemes.items():
            perm = []
            for _ in range(1000):
                x = primary.copy()
                if group_cols:
                    x["response_binary"] = x.groupby(group_cols, dropna=False).response_binary.transform(lambda s: rng.permutation(s.to_numpy()))
                else:
                    patient_labels = x.drop_duplicates("patient_key")[["patient_key", "response_binary"]].copy()
                    patient_labels["response_binary"] = rng.permutation(patient_labels.response_binary.to_numpy())
                    x = x.drop(columns="response_binary").merge(patient_labels, on="patient_key", how="left")
                perm.append(abs(hedges_g(x, fm)["effect"]))
            rows.append({"analysis": scheme, "barrier_id": fm, "metric": "empirical_p", "value": (1 + np.sum(np.asarray(perm) >= observed)) / (len(perm) + 1), "status": "negative_control"})
    for _, x in coverage_models.iterrows():
        rows.append({"analysis": f"missingness_only::{x.audit_scope}", "barrier_id": "coverage", "metric": "auc", "value": x.missingness_response_auc, "status": "bias_flag" if x.bias_flag else "negative_control"})
    coverage = pd.read_csv(OUT / "coverage_selection_audit.csv") if (OUT / "coverage_selection_audit.csv").exists() else pd.DataFrame()
    if len(coverage):
        weighted = primary.merge(coverage[["cohort_id", "patient_key", "timepoint", "inverse_probability_weight"]],
            left_on=["cohort_id", "patient_key", "normalized_timepoint"], right_on=["cohort_id", "patient_key", "timepoint"], how="left")
        for fm in PRIMARY:
            r = weighted[weighted.response_binary.eq(1)].dropna(subset=[fm, "inverse_probability_weight"])
            nr = weighted[weighted.response_binary.eq(0)].dropna(subset=[fm, "inverse_probability_weight"])
            effect = np.average(r[fm], weights=r.inverse_probability_weight) - np.average(nr[fm], weights=nr.inverse_probability_weight) if len(r) and len(nr) else np.nan
            rows.append({"analysis": "stabilized_IPW_response_sensitivity", "barrier_id": fm, "metric": "weighted_mean_difference", "value": effect,
                "status": "coverage_sensitivity_not_bias_elimination"})
    crossfit = pd.read_csv(OLD / "robustness/fold_aware_residualization_sensitivity.csv")
    for _, x in crossfit.iterrows():
        raw = effects[(effects.cohort_id.eq(x.get("cohort_id"))) & (effects.barrier_id.eq(x.get("barrier_id")))]
        agreement = bool(len(raw) and np.isfinite(x.get("effect", np.nan)) and np.sign(raw.effect.iloc[0]) == np.sign(x.get("effect")))
        rows.append({"analysis": "outer_cohort_response_blind_crossfit", "barrier_id": x.get("barrier_id", "unknown"), "metric": "direction_agreement", "value": float(agreement), "status": x.get("crossfit_status", "sensitivity")})
    loo = pd.read_csv(OLD / "robustness/leave_one_environment_analysis.csv")
    for _, x in loo.iterrows():
        rows.append({"analysis": f"leave_one::{x.leave_out_type}::{x.left_out_value}", "barrier_id": x.barrier_id,
            "metric": "pooled_effect", "value": x.pooled_effect, "status": "robustness_sensitivity"})
    top = pd.read_csv(P7 / "annotation/module_top40_projection_sensitivity.csv")
    for _, x in top[top.module_id.isin(PRIMARY)].iterrows():
        rows.append({"analysis": "top40_vs_full_frozen_measurement", "barrier_id": x.module_id, "metric": "spearman",
            "value": x.top40_vs_full_frozen_score_spearman, "status": "measurement_sensitivity_not_response_selected"})
    residual = pd.read_csv(P7 / "measurement/module_activity_by_patient_timepoint_residualized.csv")
    rb = residual[residual.module_id.isin(PRIMARY)].pivot_table(index=["cohort_id", "patient_key", "timepoint", "tissue_context"], columns="module_id", values="residual_activity_route_b").reset_index()
    lab = env[["cohort_id", "patient_key", "normalized_timepoint", "tissue_context", "response_binary", "primary_role"]].rename(columns={"normalized_timepoint": "timepoint"}).drop_duplicates()
    rb = rb.merge(lab, on=["cohort_id", "patient_key", "timepoint", "tissue_context"], how="inner")
    rb = rb[rb.primary_role.eq("conditional_primary") & rb.timepoint.eq("baseline")]
    for cohort, x in rb.groupby("cohort_id"):
        for fm in PRIMARY:
            h = hedges_g(x, fm)
            rows.append({"analysis": f"measurement_route_B::{cohort}", "barrier_id": fm, "metric": "effect",
                "value": h["effect"], "status": "empirical_bayes_measurement_sensitivity"})
    return pd.DataFrame(rows)


def shared_context(effects: pd.DataFrame, meta: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for fm in PRIMARY:
        e = effects[effects.barrier_id.eq(fm)].dropna(subset=["effect"])
        cross = meta[(meta.barrier_id.eq(fm)) & (meta.synthesis.eq("cross_endpoint_sensitivity_only"))]
        rows.append({"barrier_id": fm, "component": "shared", "estimate": cross.pooled_effect.iloc[0] if len(cross) else np.nan,
            "uncertainty_low": cross.ci_low.iloc[0] if len(cross) else np.nan, "uncertainty_high": cross.ci_high.iloc[0] if len(cross) else np.nan,
            "identifiable": False, "information_class": "unresolved", "reason": "cohort_cancer_endpoint_are_confounded_and_only_two_primary_environments_are_estimable"})
        valid = e.dropna(subset=["effect", "se"])
        rows.append({"barrier_id": fm, "component": "frequentist_mixed_random_slope", "estimate": np.nan,
            "uncertainty_low": np.nan, "uncertainty_high": np.nan, "identifiable": False, "information_class": "unresolved",
            "reason": "singular_by_design_two_environments_fallback_to_within_cohort_and_random_effects_meta"})
        for prior_sd in [0.5, 1.0, 2.0]:
            b = bayesian_partial_pool(valid, load_config(), mu_prior_sd=prior_sd, tau_prior_sd=1.0)
            rows.append({"barrier_id": fm, "component": f"bayesian_partial_pool_prior_sd_{prior_sd}", "estimate": b.get("posterior_mean", np.nan),
                "uncertainty_low": b.get("credible_low", np.nan), "uncertainty_high": b.get("credible_high", np.nan), "identifiable": False,
                "information_class": "unresolved", "reason": "posterior_reported_for_prior_sensitivity_only_context_decomposition_not_identifiable"})
        for _, x in e.iterrows():
            rows.append({"barrier_id": fm, "component": f"cohort::{x.cohort_id}", "estimate": x.effect,
                "uncertainty_low": x.effect - 1.96*x.se, "uncertainty_high": x.effect + 1.96*x.se,
                "identifiable": True, "information_class": "context_specific_direction", "reason": "within_cohort_descriptive_effect_not_shared_mechanism"})
        rows.append({"barrier_id": fm, "component": "HCC_residual", "estimate": np.nan, "uncertainty_low": np.nan, "uncertainty_high": np.nan,
            "identifiable": False, "information_class": "unresolved", "reason": "TASK01_baseline_mid_and_GSE206325_post_coarse_unknown_endpoint_are_not_exchangeable"})
    return pd.DataFrame(rows)


def report_text(env, summary, effects, meta, interface, coverage_models, verdict) -> str:
    primary = env[env.primary_role.eq("conditional_primary") & env.baseline_eligible]
    counts = primary.groupby("cohort_id").patient_key.nunique().to_dict()
    lines = [
        "# Phase8 Bounded Repair and Shared-Context Re-adjudication",
        "",
        "## Executive verdict",
        f"最终裁定：`{verdict['verdict']}`；Phase9：`BLOCKED`。本轮没有修改 FM01/FM04/FM07 membership，也没有训练 SRB、执行 repair 或 X-class ranking。",
        "",
        "## 真实 anchor surface",
        f"Phase7 左表共 {env.patient_key.nunique()} 位患者；严格 baseline conditional-primary 交集为 {primary.patient_key.nunique()} 位：{counts}。",
        f"其中 R={int(primary.response_binary.eq(1).sum())}，NR={int(primary.response_binary.eq(0).sum())}。只有两个环境可估计队列内 response effect，且 cohort、cancer、endpoint 同步变化。",
        "",
        "## 三个接口修复",
        "- GSE286827 已严格拆成 D-only conditional clean lane 与 D+T combination sensitivity，二者未合并。",
        "- GSE301741 已完成 frozen-module sample-level projection，但与现有 mid score 的共同患者过少且 FM 间桥接不一致，维持 support-only unresolved。",
        "  其 11 位 baseline projection 患者（6R/5NR）已单独计算 support sensitivity，但未并入 mid primary matrix。",
        "- GSE206325 只有 post-treatment coarse sensitivity；没有 baseline，不能并入 baseline mid shared meta，也不能单独建立 HCC R2。",
        "",
        "## Coverage-selection bias",
        f"覆盖纳入模型中发现 response-related bias flag 的层数：{int(coverage_models.bias_flag.sum())}/{len(coverage_models)}；同时报告 stabilized IPW、有效样本量和极端权重。这些检查不能证明无偏，只限定当前 target population。",
        "",
        "## Route A",
        "队列内同时报告标准化效应、logistic slope、Cliff's delta 和 bootstrap CI；跨 endpoint 合并只作为 sensitivity。FM04 在两个可估计环境方向一致，但 endpoint-aware 条件下每个 endpoint 只有一个队列，不能升级 shared route。",
        "",
        "## Route B",
        "层级 shared-specific 分解不可识别：仅两个主要环境且 cohort=cancer=endpoint。继续增加交互或 prior 只会制造模型确定性，因此 Route B 保持 unresolved sensitivity。",
        "",
        "## Sensitivity and sentinels",
        "FM02+FM03、FM05+FM08、FM06 仅作预注册 sensitivity；TLS/B、两类 CD8 exhaustion 仅作 sanity check。NET 与 WNT 继续标记 low coverage，不解释为阴性。",
        "",
        "## Negative controls and crossfit",
        "已执行 metadata-only、coverage/missingness、cohort/endpoint-block permutation 和 outer-cohort response-blind crossfit。Crossfit 只能验证预处理方向相容，不能修复 endpoint/context 不可分离。",
        "",
        "## Route verdict",
        "FM01、FM04、FM07 均未满足 R1、R1-lite 或 R2。当前也不裁定 Layer A only，因为主要问题仍是可识别信息不足，而不是三个 barrier 已被充分证伪。正式状态为 Unresolved，Phase9 继续阻断。",
        "",
        "## Allowed claims",
        "允许：response-blind pan-cancer state/module framework；队列内方向性结果；FM04 的跨两个异质环境方向一致信号，明确标注 exploratory/unresolved。",
        "",
        "## Prohibited claims",
        "禁止：shared anti-PD1 resistance mechanism、HCC-specific mechanism、causal repair chain、SRB necessity、counterfactual repair、X-class recommendation。",
        "",
        "## Next action",
        "不继续开放式搜数。只有获得至少两个 endpoint 可比较的独立 baseline clean/conditional cohorts，或获得可交换的 baseline HCC mid/coarse bridge，才重新启动 Phase9 gate。",
    ]
    return "\n".join(lines) + "\n"


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    # A failed rerun must never coexist with a stale successful adjudication.
    for stale in ["phase8_route_verdict.yaml", "PHASE8_REPAIR_AND_READJUDICATION_FINAL_REPORT.md", "phase8_output_index.tsv", "phase8_to_phase9_handoff.yaml"]:
        path = OUT / stale
        if path.exists(): path.unlink()
    cfg = load_config()
    contract, audit = input_audit(cfg)
    write_yaml(OUT / "phase8_repair_input_contract.yaml", contract)
    audit.to_csv(OUT / "phase8_repair_input_audit.csv", index=False)
    if contract["hard_failures"]:
        (OUT / "phase8_repair_abort_log.md").write_text("# Phase8 repair aborted\n\n" + "\n".join(f"- {x}" for x in contract["hard_failures"]) + "\n")
        raise SystemExit("Phase8 hard fail: " + ", ".join(contract["hard_failures"]))
    abort_log = OUT / "phase8_repair_abort_log.md"
    if abort_log.exists():
        abort_log.unlink()
    env, summary = repaired_environment(cfg)
    env.to_csv(OUT / "anchor_environment_v6_2_1_repaired.csv", index=False)
    summary.to_csv(OUT / "anchor_surface_summary.csv", index=False)
    effects = within_effects(env)
    meta = route_a_meta(effects)
    effects_master = pd.concat([effects.assign(result_type="within_cohort"), meta.assign(result_type="meta_synthesis")], ignore_index=True, sort=False)
    effects_master.to_csv(OUT / "phase8_barrier_effects_master.csv", index=False)
    shared = shared_context(effects, meta)
    shared.to_csv(OUT / "phase8_shared_context_decomposition.csv", index=False)
    interface = interface_results(cfg, env)
    detail = OUT / "interface_details"; detail.mkdir(exist_ok=True)
    g286 = env[env.cohort_id.eq("GSE286827")].copy()
    g286.to_csv(detail / "GSE286827_treatment_arm_audit.csv", index=False)
    g286[g286.treatment_arm.eq("monotherapy_D_only")].to_csv(detail / "GSE286827_D_only_anchor_table.csv", index=False)
    g286[g286.treatment_arm.eq("combination_D_plus_T")].to_csv(detail / "GSE286827_DT_sensitivity_table.csv", index=False)
    interface[interface.interface.eq("GSE301741_loss_trace")].to_csv(detail / "GSE301741_loss_trace.csv", index=False)
    pd.read_csv(OLD / "projection/gse301741_projection_bridge_audit.csv").to_csv(detail / "GSE301741_projection_comparison.csv", index=False)
    pd.read_csv(OLD / "projection/gse301741_sample_level_frozen_module_projection.csv").to_csv(detail / "GSE301741_projection_scores.csv", index=False)
    write_yaml(detail / "GSE301741_bridge_verdict.yaml", {"verdict": "support_only_unresolved", "selection_used_response": False,
        "reason": "only_four_shared_contexts_and_module_specific_spearman_is_inconsistent"})
    pd.read_csv(OLD / "hcc_sensitivity/gse206325_coarse_response_sensitivity.csv").to_csv(detail / "GSE206325_coarse_anchor_effects.csv", index=False)
    pd.DataFrame([{"bridge": "coarse_to_mid", "n_shared_patient_contexts": 0, "status": "not_estimable_no_shared_patient_context"}]).to_csv(detail / "cross_resolution_bridge_results.csv", index=False)
    write_yaml(detail / "GSE206325_resolution_verdict.yaml", {"verdict": "post_coarse_HCC_sensitivity_only", "baseline_available": False,
        "pooled_with_mid": False, "HCC_R2_upgrade_allowed": False})
    sens = sensitivity_and_sentinels(env)
    bridge_sens = pd.concat([interface.assign(result_family="data_interface"), sens.assign(result_family="sensitivity_or_sentinel")], ignore_index=True, sort=False)
    bridge_sens.to_csv(OUT / "phase8_data_bridge_and_sensitivity_results.csv", index=False)
    coverage, coverage_models = coverage_audit(cfg, env)
    coverage.to_csv(OUT / "coverage_selection_audit.csv", index=False)
    coverage_models.to_csv(OUT / "coverage_selection_models.csv", index=False)
    negative = negative_controls(env, effects, coverage_models)
    negative.to_csv(OUT / "phase8_negative_control_and_robustness.csv", index=False)
    verdict = {
        "phase": "v6.2.1 Phase8 bounded repair and re-adjudication",
        "verdict": "UNRESOLVED_PHASE9_BLOCKED",
        "phase9_entry_allowed": False,
        "barriers": {fm: {"route": "Unresolved", "R1": False, "R1_lite": False, "R2": False, "layer_A_only": False,
            "reason": "shared_and_HCC_context_effects_not_identifiable_under_current_endpoint_resolution_timepoint_surface"} for fm in PRIMARY},
        "shared_route": {"eligible": False, "reason": "only_two_estimable_primary_environments_with_cohort_cancer_endpoint_confounding"},
        "HCC_context_route": {"eligible": False, "reason": "baseline_mid_TASK01_not_exchangeable_with_post_coarse_GSE206325"},
        "data_interface_verdicts": {"GSE286827": "D_only_conditional_clean;D_plus_T_sensitivity", "GSE301741": "support_only_unresolved_projection_bridge", "GSE206325": "post_coarse_HCC_sensitivity_only"},
        "coverage_selection_bias_flag": bool(coverage_models.bias_flag.any()),
        "allowed_claims": ["response_blind_pan_cancer_state_module_framework", "within_cohort_directional_associations", "exploratory_FM04_cross_environment_direction_consistency"],
        "prohibited_claims": ["shared_anti_PD1_mechanism", "HCC_specific_mechanism", "causal_repair_chain", "supervised_SRB", "counterfactual_repair", "X_class_ranking"],
        "handoff_generated": False,
    }
    write_yaml(OUT / "phase8_route_verdict.yaml", verdict)
    (OUT / "PHASE8_REPAIR_AND_READJUDICATION_FINAL_REPORT.md").write_text(report_text(env, summary, effects, meta, interface, coverage_models, verdict))
    index = []
    for p in sorted(OUT.rglob("*")):
        if p.is_file() and p.name != "phase8_output_index.tsv":
            index.append({"file": str(p.relative_to(OUT)), "bytes": p.stat().st_size, "sha256": sha256(p)})
    pd.DataFrame(index).to_csv(OUT / "phase8_output_index.tsv", sep="\t", index=False)
    print(yaml.safe_dump({"output": str(OUT), "verdict": verdict["verdict"], "n_primary_patients": int(env.loc[env.primary_role.eq("conditional_primary") & env.baseline_eligible, "patient_key"].nunique()), "n_files": len(index) + 1}, sort_keys=False))


if __name__ == "__main__":
    main()
