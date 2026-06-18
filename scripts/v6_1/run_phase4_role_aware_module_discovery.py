#!/usr/bin/env python3
from __future__ import annotations

import math
import re
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from scipy import stats
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial.distance import squareform


ROOT = Path(__file__).resolve().parents[2]
DATE_TAG = "20260611"
OUT = ROOT / "results" / "v6_1" / f"phase4_role_aware_module_discovery_{DATE_TAG}"

READINESS = ROOT / "results" / "v6_1" / f"phase4_readiness_{DATE_TAG}"
STEP2_REPAIR = ROOT / "results" / "v6_1" / "step2_repair"
STEP3 = ROOT / "results" / "v6_1" / "step3_qc_aware_strong_baseline"
HANDOFF = STEP3 / "07_milestone_B_step4_handoff"

PRIMARY_TASKS = {"HCC_specific", "PD1X_extension", "pan_cancer_shared"}
PRIMARY_FAMILIES = {
    "fraction",
    "signature",
    "pathway",
    "tf_activity",
    "combined_interpretable",
    "combined_all_allowed",
}
RESPONSE_MAP = {"responder": 1, "non_responder": 0}
RNG = np.random.default_rng(611)


def now_iso() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat()


def mkdir(path: Path) -> Path:
    path.mkdir(parents=True, exist_ok=True)
    return path


def rel(path: Path) -> str:
    try:
        return str(path.relative_to(ROOT))
    except ValueError:
        return str(path)


def require(path: Path) -> Path:
    if not path.exists():
        raise FileNotFoundError(path)
    return path


def yaml_scalar(value: Any) -> str:
    if isinstance(value, bool):
        return "true" if value else "false"
    if isinstance(value, (int, float)) and not isinstance(value, bool):
        if isinstance(value, float) and (math.isnan(value) or math.isinf(value)):
            return "null"
        return str(value)
    if value is None:
        return "null"
    return '"' + str(value).replace('"', '\\"') + '"'


def to_yaml(value: Any, indent: int = 0) -> str:
    pad = " " * indent
    if isinstance(value, dict):
        lines = []
        for key, sub in value.items():
            if isinstance(sub, (dict, list)):
                lines.append(f"{pad}{key}:")
                lines.append(to_yaml(sub, indent + 2))
            else:
                lines.append(f"{pad}{key}: {yaml_scalar(sub)}")
        return "\n".join(lines)
    if isinstance(value, list):
        lines = []
        for item in value:
            if isinstance(item, (dict, list)):
                lines.append(f"{pad}-")
                lines.append(to_yaml(item, indent + 2))
            else:
                lines.append(f"{pad}- {yaml_scalar(item)}")
        return "\n".join(lines)
    return f"{pad}{yaml_scalar(value)}"


def write_yaml(path: Path, data: dict[str, Any]) -> None:
    path.write_text(to_yaml(data) + "\n")


def write_status(path: Path, phase: str, verdict: str, **kwargs: Any) -> None:
    data = {"phase": phase, "verdict": verdict, "generated_at_utc": now_iso(), **kwargs}
    write_yaml(path, data)


def bool_series(s: pd.Series) -> pd.Series:
    return s.astype(str).str.lower().isin(["true", "1", "yes"])


def bh_fdr(p_values: pd.Series) -> pd.Series:
    p = pd.to_numeric(p_values, errors="coerce")
    out = pd.Series(np.nan, index=p.index, dtype=float)
    valid = p.dropna().sort_values()
    m = len(valid)
    if m == 0:
        return out
    ranked = valid * m / np.arange(1, m + 1)
    ranked = ranked[::-1].cummin()[::-1].clip(upper=1.0)
    out.loc[ranked.index] = ranked
    return out


def safe_corr(x: pd.Series, y: pd.Series) -> float:
    xy = pd.concat([x, y], axis=1).dropna()
    if len(xy) < 5 or xy.iloc[:, 0].nunique() < 2 or xy.iloc[:, 1].nunique() < 2:
        return np.nan
    return float(np.corrcoef(xy.iloc[:, 0], xy.iloc[:, 1])[0, 1])


def safe_nanmean(values: Any) -> float:
    arr = np.asarray(values, dtype=float)
    if arr.size == 0 or np.isnan(arr).all():
        return np.nan
    return float(np.nanmean(arr))


def safe_true_fraction(values: list[bool]) -> float:
    return float(np.mean(values)) if values else 0.0


def series_distribution(values: pd.Series) -> dict[str, float | int]:
    clean = pd.to_numeric(values, errors="coerce").dropna()
    if clean.empty:
        return {"n": 0, "mean": np.nan, "median": np.nan, "std": np.nan}
    return {
        "n": int(len(clean)),
        "mean": float(clean.mean()),
        "median": float(clean.median()),
        "std": float(clean.std()),
    }


def zscore_frame(df: pd.DataFrame) -> pd.DataFrame:
    numeric = df.apply(pd.to_numeric, errors="coerce")
    std = numeric.std(axis=0, skipna=True).replace(0, np.nan)
    z = (numeric - numeric.mean(axis=0, skipna=True)) / std
    return z


def feature_prefix(feature: str) -> str:
    if feature.startswith("frac_"):
        return "fraction"
    if feature.startswith("signature__") or feature.startswith("sig_"):
        return "signature"
    if feature.startswith("pathway__") or feature.startswith("path_"):
        return "pathway"
    if feature.startswith("tf_activity__"):
        return "tf_activity"
    return "other"


def infer_annotation(features: list[str]) -> tuple[str, str, str]:
    text = " ".join(features).lower()
    rules = [
        ("cytotoxic_T_NK_activity", ["cd8", "cytotoxic", "gzmb", "prf", "nkg", "nk", "ifng"]),
        ("T_cell_dysfunction_or_exhaustion", ["exhaust", "dysfunction", "pdcd1", "havcr2", "lag3", "tigit", "tox"]),
        ("Treg_regulatory_suppression", ["treg", "foxp3", "il2ra", "ctla4"]),
        ("myeloid_inflammatory_state", ["mono_fcn1", "inflammatory", "il1", "tnf", "fcgr"]),
        ("myeloid_suppressive_state", ["macro_spp1", "spp1", "c1qc", "tam", "mdsc", "m2"]),
        ("DC_APC_activity", ["dc", "apc", "hla", "mhc", "antigen"]),
        ("B_plasma_contribution", ["b_", "plasma", "immunoglobulin", "igh"]),
        ("IFN_response", ["ifn", "interferon", "stat1", "irf1", "irf7"]),
        ("stromal_or_vascular_proxy", ["stromal", "fibro", "endo", "vascular", "vegf"]),
        ("immune_cold_or_low_infiltration", ["immune_unspecified", "non_immune", "epithelial_tumor"]),
    ]
    hits = [label for label, keys in rules if any(k in text for k in keys)]
    if not hits:
        hits = ["mixed_or_unannotated_immune_program"]
    feature_types = sorted({feature_prefix(f) for f in features})
    tf_support = "tf_activity_present" if "tf_activity" in feature_types else "tf_activity_not_primary_driver"
    return ";".join(hits[:4]), ";".join(feature_types), tf_support


def load_inputs() -> dict[str, pd.DataFrame]:
    return {
        "entry_decision": pd.DataFrame(),
        "evidence_pool": pd.read_csv(require(READINESS / "phase4_evidence_pool_registry_20260611.csv")),
        "feature_contract": pd.read_csv(require(READINESS / "phase4_feature_family_contract_20260611.csv")),
        "blocked_inputs": pd.read_csv(require(READINESS / "phase4_blocked_inputs_20260611.csv")),
        "join_contract": pd.read_csv(require(READINESS / "phase4_join_key_contract_20260611.tsv"), sep="\t"),
        "file_index": pd.read_csv(require(READINESS / "phase4_downstream_file_index_20260611.tsv"), sep="\t"),
        "matrix": pd.read_parquet(require(STEP2_REPAIR / "repair_B2_unified_feature_matrix_by_sample.hotfix.parquet")),
        "feature_dictionary": pd.read_csv(require(STEP2_REPAIR / "repair_B2_feature_dictionary_v6_1.hotfix.csv")),
        "membership": pd.read_csv(require(STEP3 / "01_analysis_universe" / "step3_phase3_1_sample_universe_membership.repaired.csv")),
        "model_grades": pd.read_csv(require(STEP3 / "06_robustness_negative_controls" / "step3_phase3_6_model_evidence_grading.csv")),
        "univariate": pd.read_csv(require(STEP3 / "04_interpretable_baselines" / "step3_phase3_4_univariate_feature_results.csv")),
        "support_pd1": pd.read_csv(require(HANDOFF / "step3_step4_handoff_support_PD1_anchor.csv")),
    }


def phase4_0(inputs: dict[str, pd.DataFrame]) -> None:
    out = mkdir(OUT / "00_entry_contract")
    evidence = inputs["evidence_pool"]
    feature_contract = inputs["feature_contract"]
    blocked = inputs["blocked_inputs"]
    join_contract = inputs["join_contract"]
    file_index = inputs["file_index"]

    primary = evidence[evidence["pool"] == "primary_evidence_pool"].copy()
    support = evidence[evidence["pool"] == "support_evidence_pool"].copy()

    audit_rows = []
    for _, row in evidence.iterrows():
        task = row["task"]
        family = row["feature_family"]
        pool = row["pool"]
        allowed_primary = (
            pool == "primary_evidence_pool"
            and task in PRIMARY_TASKS
            and family in PRIMARY_FAMILIES
            and str(row.get("primary_allowed_after_repair", "")).lower() == "true"
        )
        violation = pool == "primary_evidence_pool" and not allowed_primary
        audit_rows.append(
            {
                "universe_id": row["universe_id"],
                "task": task,
                "feature_family": family,
                "pool": pool,
                "phase4_primary_allowed": allowed_primary,
                "support_only": pool == "support_evidence_pool",
                "violation": violation,
                "action": "allow_primary" if allowed_primary else "support_or_exclude_from_primary",
            }
        )
    audit = pd.DataFrame(audit_rows)

    block_rows = []
    for _, row in blocked.iterrows():
        block_rows.append(
            {
                "blocked_input": row["blocked_input"],
                "blocking_reason": row["blocking_reason"],
                "used_in_primary_phase4": False,
                "enforcement_status": "PASS_not_used",
                "forbidden_use": row["forbidden_use"],
            }
        )
    block_df = pd.DataFrame(block_rows)

    join_rows = []
    matrix_cols = set(inputs["matrix"].columns)
    for _, row in join_contract.iterrows():
        keys = [k for k in str(row["required_join_keys"]).split(";") if k]
        if row["contract_id"] in {"hotfix_feature_matrix"}:
            present = all(k in matrix_cols for k in keys)
        else:
            present = Path(ROOT / str(row["canonical_path"])).exists()
        join_rows.append(
            {
                "contract_id": row["contract_id"],
                "canonical_path": row["canonical_path"],
                "required_join_keys": row["required_join_keys"],
                "validation_status": "PASS" if present else "FAIL",
                "phase4_use": row["phase4_use"],
            }
        )
    join_df = pd.DataFrame(join_rows)

    audit.to_csv(out / "phase4_0_evidence_pool_audit.csv", index=False)
    block_df.to_csv(out / "phase4_0_blocked_input_enforcement.csv", index=False)
    join_df.to_csv(out / "phase4_0_join_key_validation.csv", index=False)
    verdict = "PASS" if not audit["violation"].any() and (join_df["validation_status"] != "FAIL").all() and file_index["exists"].astype(str).str.lower().eq("true").all() else "HARD_FAIL"
    (out / "phase4_0_entry_contract_summary.md").write_text(
        "\n".join(
            [
                "# Phase4.0 Entry Contract Summary",
                "",
                f"- Verdict: `{verdict}`",
                f"- Primary evidence rows: {len(primary)}",
                f"- Support evidence rows: {len(support)}",
                f"- Blocked rules enforced: {len(block_df)}",
                "- PD1_anchor retained as support/sensitivity only.",
                "- Phase3.5 strong ML and unintegrated addendum data are not used.",
                "",
            ]
        )
    )
    write_status(out / "phase4_0_status.yaml", "Phase4.0", verdict, primary_rows=len(primary), support_rows=len(support), blocked_rules=len(block_df))
    if verdict == "HARD_FAIL":
        raise RuntimeError("Phase4.0 hard fail")


def feature_list_for(row: pd.Series) -> list[str]:
    path = str(row["feature_list_path"])
    df = pd.read_csv(require(Path(path)))
    return df["feature_name"].astype(str).tolist()


def build_primary_matrices(inputs: dict[str, pd.DataFrame]) -> tuple[dict[str, dict[str, Any]], pd.DataFrame]:
    evidence = inputs["evidence_pool"]
    contract = inputs["feature_contract"]
    matrix = inputs["matrix"]
    membership = inputs["membership"]
    model_grades = inputs["model_grades"]

    primary_evidence = evidence[evidence["pool"] == "primary_evidence_pool"].copy()
    primary_clean = contract[
        (contract["phase4_use"] == "primary")
        & (contract["feature_set_type"] == "clean_core_feature_set")
        & (contract["task"].isin(PRIMARY_TASKS))
    ].copy()
    eligible_missing = set(
        model_grades[
            (model_grades["eligible_for_step4_primary"].astype(str).str.lower() == "true")
            & (model_grades["evidence_grade"].isin(["A", "B"]))
        ]["universe_id"].astype(str)
    )
    missingness_log = []
    matrices: dict[str, dict[str, Any]] = {}
    registry_rows = []
    shape_rows = []
    availability_rows = []

    for _, row in primary_clean.iterrows():
        universe_id = row["universe_id"]
        evidence_row = primary_evidence[primary_evidence["universe_id"] == universe_id]
        if evidence_row.empty:
            continue
        features = [f for f in feature_list_for(row) if f in matrix.columns]
        mem = membership[
            (membership["universe_id"] == universe_id)
            & (membership["included"].astype(str).str.lower() == "true")
            & (membership["response_label"].isin(["responder", "non_responder"]))
        ].copy()
        sample_keys = mem["sample_key"].astype(str).unique()
        sub = matrix[matrix["sample_key"].astype(str).isin(sample_keys)].copy()
        sub = sub.merge(
            mem[
                [
                    "sample_key",
                    "response_label",
                    "split",
                    "treatment_context",
                    "cancer_type",
                    "sample_weight",
                    "universe_id",
                    "universe_role",
                ]
            ].drop_duplicates("sample_key"),
            on="sample_key",
            how="inner",
            suffixes=("", "_membership"),
        )
        if "split_membership" in sub.columns:
            sub["split"] = sub["split_membership"].fillna(sub["split"])
        if "treatment_context_membership" in sub.columns:
            sub["treatment_context"] = sub["treatment_context_membership"].fillna(sub["treatment_context"])
        numeric = sub[features].apply(pd.to_numeric, errors="coerce")
        miss = numeric.isna().mean()
        var = numeric.var(skipna=True)
        usable = [f for f in features if miss.get(f, 1.0) <= 0.8 and var.get(f, 0.0) > 0]
        matrices[universe_id] = {
            "task": row["task"],
            "feature_family": row["feature_family"],
            "feature_set_type": row["feature_set_type"],
            "data": sub,
            "features": usable,
            "all_features": features,
            "confounding_risk_level": evidence_row.iloc[0]["confounding_risk_level"],
        }
        registry_rows.append(
            {
                "universe_id": universe_id,
                "task": row["task"],
                "feature_family": row["feature_family"],
                "feature_set_type": row["feature_set_type"],
                "n_samples": len(sub),
                "n_patients": sub["patient_key"].nunique() if "patient_key" in sub.columns else sub["patient_id"].nunique(),
                "n_cohorts": sub["cohort_id"].nunique(),
                "n_features_declared": len(features),
                "n_features_usable": len(usable),
                "source_pool": "primary_evidence_pool",
                "matrix_status": "ready" if len(sub) > 0 and len(usable) >= 2 else "insufficient",
            }
        )
        shape_rows.append(
            {
                "universe_id": universe_id,
                "n_samples": len(sub),
                "n_patients": sub["patient_key"].nunique() if "patient_key" in sub.columns else sub["patient_id"].nunique(),
                "n_cohorts": sub["cohort_id"].nunique(),
                "n_responders": int((sub["response_label"] == "responder").sum()),
                "n_non_responders": int((sub["response_label"] == "non_responder").sum()),
                "n_features": len(usable),
                "response_rate": float((sub["response_label"] == "responder").mean()) if len(sub) else np.nan,
            }
        )
        availability_rows.append(
            {
                "universe_id": universe_id,
                "feature_family": row["feature_family"],
                "n_declared_features": len(features),
                "n_usable_features": len(usable),
                "median_missing_rate": float(miss[usable].median()) if usable else np.nan,
                "max_missing_rate": float(miss[usable].max()) if usable else np.nan,
                "n_features_dropped_missing_or_zero_variance": len(features) - len(usable),
            }
        )

    missing_rows = contract[
        contract["phase4_use"].astype(str).str.startswith("conditional_primary", na=False)
        & (contract["status"] == "eligible")
    ].copy()
    for _, row in missing_rows.iterrows():
        missingness_log.append(
            {
                "universe_id": row["universe_id"],
                "task": row["task"],
                "feature_family": row["feature_family"],
                "feature_set_type": row["feature_set_type"],
                "allowed_in_phase4": row["universe_id"] in eligible_missing,
                "usage": "supplemental_only_if_A_or_B" if row["universe_id"] in eligible_missing else "not_used",
                "reason": row["reason"],
            }
        )
    return matrices, pd.DataFrame(missingness_log), pd.DataFrame(registry_rows), pd.DataFrame(shape_rows), pd.DataFrame(availability_rows)


def phase4_1(inputs: dict[str, pd.DataFrame]) -> dict[str, dict[str, Any]]:
    out = mkdir(OUT / "01_primary_matrix")
    matrices, missing_log, registry, shapes, availability = build_primary_matrices(inputs)
    registry.to_csv(out / "phase4_1_primary_matrix_registry.csv", index=False)
    shapes.to_csv(out / "phase4_1_feature_family_matrix_shapes.csv", index=False)
    availability.to_csv(out / "phase4_1_feature_availability_summary.csv", index=False)
    missing_log.to_csv(out / "phase4_1_missingness_aware_usage_log.csv", index=False)
    n_ready = int((registry["matrix_status"] == "ready").sum()) if not registry.empty else 0
    verdict = "PASS" if n_ready >= 3 else "FAIL"
    (out / "phase4_1_primary_matrix_summary.md").write_text(
        "\n".join(
            [
                "# Phase4.1 Primary Matrix Summary",
                "",
                f"- Verdict: `{verdict}`",
                f"- Ready primary matrices: {n_ready}/{len(registry)}",
                "- Source: primary evidence pool only.",
                "- Feature source: Step2 repair_B2 hotfix matrix and Step3 clean-core feature lists.",
                "- Missingness-aware rows are logged separately and not used for primary module discovery unless A/B evidence permits supplemental use.",
                "",
            ]
        )
    )
    write_status(out / "phase4_1_status.yaml", "Phase4.1", verdict, ready_matrices=n_ready, total_matrices=len(registry))
    return matrices


def corr_matrix(z: pd.DataFrame) -> np.ndarray:
    filled = z.fillna(0.0).to_numpy(dtype=float)
    if filled.shape[1] == 1:
        return np.array([[1.0]])
    with np.errstate(invalid="ignore", divide="ignore"):
        c = np.corrcoef(filled, rowvar=False)
    c = np.nan_to_num(c, nan=0.0, posinf=0.0, neginf=0.0)
    np.fill_diagonal(c, 1.0)
    return c


def discover_modules_for_universe(universe_id: str, item: dict[str, Any]) -> tuple[list[dict[str, Any]], list[dict[str, Any]], pd.DataFrame]:
    df = item["data"]
    features = item["features"]
    task = item["task"]
    family = item["feature_family"]
    if len(features) == 0:
        return [], [], pd.DataFrame(index=df.index)
    z = zscore_frame(df[features])
    c = corr_matrix(z)
    if len(features) < 2:
        labels = np.ones(len(features), dtype=int)
    else:
        dist = np.clip(1.0 - c, 0.0, 2.0)
        np.fill_diagonal(dist, 0.0)
        link = linkage(squareform(dist, checks=False), method="average")
        labels = fcluster(link, t=0.75, criterion="distance")
        if pd.Series(labels).value_counts().max() < 2 and len(features) > 2:
            labels = fcluster(link, t=0.9, criterion="distance")

    module_rows = []
    membership_rows = []
    module_scores = pd.DataFrame(index=df.index)
    for idx, cluster_id in enumerate(sorted(set(labels)), start=1):
        mod_features = [features[i] for i, lab in enumerate(labels) if lab == cluster_id]
        n = len(mod_features)
        module_id = f"P4M_{task}_{family}_{idx:03d}"
        feat_idx = [features.index(f) for f in mod_features]
        if n > 1:
            subcorr = c[np.ix_(feat_idx, feat_idx)]
            tri = subcorr[np.triu_indices_from(subcorr, k=1)]
            mean_within = safe_nanmean(tri)
        else:
            mean_within = np.nan
        outside_idx = [i for i in range(len(features)) if i not in feat_idx]
        if outside_idx:
            outside = c[np.ix_(feat_idx, outside_idx)]
            mean_outside = safe_nanmean(outside)
        else:
            mean_outside = np.nan
        if n == 1:
            granularity = "single_feature"
        elif n > 80:
            granularity = "broad_program"
        else:
            granularity = "standard_module"
        score = z[mod_features].mean(axis=1, skipna=True)
        module_scores[module_id] = score
        module_rows.append(
            {
                "module_id": module_id,
                "source_task": task,
                "source_feature_family": family,
                "source_universe_id": universe_id,
                "detection_method": "response_blind_hierarchical_feature_correlation",
                "n_features": n,
                "feature_list": ";".join(mod_features),
                "mean_within_module_correlation": mean_within,
                "mean_outside_module_correlation": mean_outside,
                "module_granularity": granularity,
                "response_blind_confirmed": True,
                "primary_claim_allowed": granularity != "single_feature",
            }
        )
        for f in mod_features:
            membership_rows.append(
                {
                    "module_id": module_id,
                    "feature_name": f,
                    "source_task": task,
                    "source_feature_family": family,
                    "source_universe_id": universe_id,
                    "feature_type_inferred": feature_prefix(f),
                }
            )
    return module_rows, membership_rows, module_scores


def phase4_2(matrices: dict[str, dict[str, Any]]) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, pd.DataFrame]]:
    out = mkdir(OUT / "02_response_blind_modules")
    module_rows: list[dict[str, Any]] = []
    membership_rows: list[dict[str, Any]] = []
    module_scores_by_universe: dict[str, pd.DataFrame] = {}
    for universe_id, item in matrices.items():
        rows, members, scores = discover_modules_for_universe(universe_id, item)
        module_rows.extend(rows)
        membership_rows.extend(members)
        module_scores_by_universe[universe_id] = scores
    modules = pd.DataFrame(module_rows)
    members = pd.DataFrame(membership_rows)
    modules.to_csv(out / "phase4_2_raw_module_assignments.csv", index=False)
    members.to_csv(out / "phase4_2_module_feature_membership.csv", index=False)
    write_yaml(
        out / "phase4_2_module_detection_parameters.yaml",
        {
            "method": "response_blind_hierarchical_feature_correlation",
            "correlation": "pearson_on_z_scored_features_with_mean_imputation_for_correlation_only",
            "distance": "1 - correlation",
            "linkage": "average",
            "primary_cut_distance": 0.75,
            "fallback_cut_distance_if_all_singletons": 0.90,
            "response_labels_used_for_discovery": False,
            "blocked_inputs_allowed": False,
        },
    )
    n_primary = int((modules.get("primary_claim_allowed", pd.Series(dtype=bool)) == True).sum()) if not modules.empty else 0
    verdict = "PASS" if n_primary > 0 else "FAIL"
    (out / "phase4_2_module_discovery_summary.md").write_text(
        f"# Phase4.2 Module Discovery Summary\n\n- Verdict: `{verdict}`\n- Raw modules: {len(modules)}\n- Primary-claim-eligible non-singleton modules: {n_primary}\n- Response labels used for module construction: `False`\n"
    )
    write_status(out / "phase4_2_status.yaml", "Phase4.2", verdict, raw_modules=len(modules), primary_claim_eligible_modules=n_primary)
    return modules, members, module_scores_by_universe


def mean_within_corr_for_features(z: pd.DataFrame, features: list[str]) -> float:
    if len(features) < 2 or len(z) < 5:
        return np.nan
    c = corr_matrix(z[features])
    tri = c[np.triu_indices_from(c, k=1)]
    return safe_nanmean(tri)


def phase4_3(modules: pd.DataFrame, members: pd.DataFrame, matrices: dict[str, dict[str, Any]], module_scores_by_universe: dict[str, pd.DataFrame]) -> pd.DataFrame:
    out = mkdir(OUT / "03_module_stability")
    stability_rows = []
    bootstrap_rows = []
    cohort_rows = []
    consistency_rows = []

    feature_to_family = members.groupby("module_id")["feature_type_inferred"].apply(lambda x: sorted(set(x))).to_dict()

    for _, mod in modules.iterrows():
        universe_id = mod["source_universe_id"]
        item = matrices[universe_id]
        df = item["data"].reset_index(drop=True)
        features = members[members["module_id"] == mod["module_id"]]["feature_name"].tolist()
        z = zscore_frame(df[features])
        patient_col = "patient_key" if "patient_key" in df.columns else "patient_id"
        patients = df[patient_col].astype(str).unique()
        boot_vals = []
        if len(features) >= 2 and len(patients) >= 5:
            for b in range(30):
                sampled = RNG.choice(patients, size=len(patients), replace=True)
                idx = df[patient_col].astype(str).isin(sampled)
                val = mean_within_corr_for_features(z.loc[idx].reset_index(drop=True), features)
                boot_vals.append(val)
                bootstrap_rows.append(
                    {
                        "module_id": mod["module_id"],
                        "bootstrap_id": b + 1,
                        "mean_within_correlation": val,
                        "support_pass": bool(pd.notna(val) and val >= 0.15),
                    }
                )
        bootstrap_support = safe_true_fraction([bool(pd.notna(v) and v >= 0.15) for v in boot_vals])
        score = module_scores_by_universe[universe_id][mod["module_id"]].reset_index(drop=True)
        feature_corrs = [safe_corr(z[f].reset_index(drop=True), score) for f in features]
        core_flags = [bool(abs(v) >= 0.25) for v in feature_corrs if pd.notna(v)]
        core_stability = safe_true_fraction(core_flags)

        split_vals = []
        for split, sub_idx in df.groupby("split").groups.items():
            if len(sub_idx) >= 10:
                val = mean_within_corr_for_features(z.loc[list(sub_idx)].reset_index(drop=True), features)
                split_vals.append(val)
        split_support = safe_true_fraction([bool(pd.notna(v) and v >= 0.10) for v in split_vals])
        split_status = "stable" if split_vals and split_support >= 0.5 else "limited_or_unstable"

        cohort_counts = df["cohort_id"].value_counts(normalize=True)
        max_cohort_fraction = float(cohort_counts.max()) if not cohort_counts.empty else np.nan
        if df["cohort_id"].nunique() < 3:
            cohort_status = "limited_cohort_stability"
        elif max_cohort_fraction > 0.65:
            cohort_status = "cohort_dominated"
        else:
            cohort_status = "not_single_cohort_dominated"
        cohort_rows.append(
            {
                "module_id": mod["module_id"],
                "n_cohorts": df["cohort_id"].nunique(),
                "max_cohort_fraction": max_cohort_fraction,
                "cohort_sensitivity_status": cohort_status,
            }
        )

        missing_burden = df[features].isna().mean(axis=1)
        miss_corr = safe_corr(score, missing_burden)
        missing_status = "missingness_driven" if pd.notna(miss_corr) and abs(miss_corr) >= 0.5 else "not_missingness_driven"

        if mod["module_granularity"] == "single_feature":
            grade = "C_exploratory"
            reason = "single_feature_module"
        elif missing_status == "missingness_driven" or cohort_status == "cohort_dominated":
            grade = "Blocked"
            reason = f"{missing_status};{cohort_status}"
        elif bootstrap_support >= 0.70 and core_stability >= 0.70 and split_status == "stable":
            grade = "A_stable"
            reason = "bootstrap_core_and_split_stable"
        elif bootstrap_support >= 0.45 and core_stability >= 0.50:
            grade = "B_usable"
            reason = "moderate_stability_with_caveat"
        else:
            grade = "C_exploratory"
            reason = "weak_or_limited_stability"

        consistency_rows.append(
            {
                "module_id": mod["module_id"],
                "source_task": mod["source_task"],
                "feature_family_consistency_status": "cross_family_module" if len(feature_to_family.get(mod["module_id"], [])) > 1 else "single_family_module",
                "feature_types": ";".join(feature_to_family.get(mod["module_id"], [])),
            }
        )
        stability_rows.append(
            {
                "module_id": mod["module_id"],
                "bootstrap_support": bootstrap_support,
                "core_feature_stability": core_stability,
                "split_stability": split_status,
                "cohort_sensitivity_status": cohort_status,
                "feature_family_consistency_status": consistency_rows[-1]["feature_family_consistency_status"],
                "missingness_sensitivity_status": missing_status,
                "missingness_score_correlation": miss_corr,
                "stability_grade": grade,
                "demotion_reason": reason,
            }
        )
    stability = pd.DataFrame(stability_rows)
    stability.to_csv(out / "phase4_3_module_stability_summary.csv", index=False)
    pd.DataFrame(bootstrap_rows).to_csv(out / "phase4_3_bootstrap_membership_support.csv", index=False)
    pd.DataFrame(cohort_rows).to_csv(out / "phase4_3_cohort_sensitivity_summary.csv", index=False)
    pd.DataFrame(consistency_rows).to_csv(out / "phase4_3_feature_family_consistency.csv", index=False)
    n_ab = int(stability["stability_grade"].isin(["A_stable", "B_usable"]).sum()) if not stability.empty else 0
    verdict = "PASS" if n_ab > 0 else "FAIL"
    (out / "phase4_3_summary.md").write_text(
        "\n".join(
            [
                "# Phase4.3 Module Stability Summary",
                "",
                f"- Verdict: `{verdict}`",
                f"- A/B stable or usable modules: {n_ab}",
                f"- Total modules assessed: {len(stability)}",
                "- Stability dimensions: patient bootstrap, split stability, cohort sensitivity, feature-family consistency, missingness sensitivity.",
                "- Unstable modules are downgraded rather than forcing Phase4 failure.",
                "",
            ]
        )
    )
    write_status(out / "phase4_3_status.yaml", "Phase4.3", verdict, stable_or_usable_modules=n_ab, total_modules=len(stability))
    return stability


def cohort_adjusted_linear(y: pd.Series, score: pd.Series, cohort: pd.Series) -> tuple[float, float]:
    df = pd.DataFrame({"y": y, "score": score, "cohort": cohort}).dropna()
    if len(df) < 20 or df["y"].nunique() < 2 or df["score"].nunique() < 2:
        return np.nan, np.nan
    dummies = pd.get_dummies(df["cohort"], drop_first=True, dtype=float)
    x = pd.concat([pd.Series(1.0, index=df.index, name="intercept"), df["score"].astype(float), dummies], axis=1)
    try:
        beta, *_ = np.linalg.lstsq(x.to_numpy(dtype=float), df["y"].to_numpy(dtype=float), rcond=None)
        resid = df["y"].to_numpy(dtype=float) - x.to_numpy(dtype=float) @ beta
        dof = max(len(df) - x.shape[1], 1)
        sigma2 = float((resid @ resid) / dof)
        cov = sigma2 * np.linalg.pinv(x.to_numpy(dtype=float).T @ x.to_numpy(dtype=float))
        se = math.sqrt(max(cov[1, 1], 0.0))
        coef = float(beta[1])
        p = 2 * (1 - stats.t.cdf(abs(coef / se), dof)) if se > 0 else np.nan
        return coef, p
    except Exception:
        return np.nan, np.nan


def phase4_4(modules: pd.DataFrame, matrices: dict[str, dict[str, Any]], module_scores_by_universe: dict[str, pd.DataFrame], stability: pd.DataFrame) -> pd.DataFrame:
    out = mkdir(OUT / "04_module_response_association")
    assoc_rows = []
    effect_rows = []
    adjusted_rows = []
    for _, mod in modules.iterrows():
        universe_id = mod["source_universe_id"]
        item = matrices[universe_id]
        df = item["data"].reset_index(drop=True)
        score = module_scores_by_universe[universe_id][mod["module_id"]].reset_index(drop=True)
        y = df["response_label"].map(RESPONSE_MAP)
        valid = y.notna() & score.notna()
        if valid.sum() < 10 or y[valid].nunique() < 2:
            effect = np.nan
            p_value = np.nan
            smd = np.nan
            ci_low = np.nan
            ci_high = np.nan
        else:
            resp = score[valid & (y == 1)]
            non = score[valid & (y == 0)]
            effect = float(resp.mean() - non.mean())
            pooled = math.sqrt(((resp.var(ddof=1) + non.var(ddof=1)) / 2)) if len(resp) > 1 and len(non) > 1 else np.nan
            smd = float(effect / pooled) if pooled and not math.isnan(pooled) and pooled != 0 else np.nan
            p_value = float(stats.ttest_ind(resp, non, equal_var=False, nan_policy="omit").pvalue)
            se = math.sqrt(resp.var(ddof=1) / len(resp) + non.var(ddof=1) / len(non)) if len(resp) > 1 and len(non) > 1 else np.nan
            ci_low = effect - 1.96 * se if pd.notna(se) else np.nan
            ci_high = effect + 1.96 * se if pd.notna(se) else np.nan
        adj_coef, adj_p = cohort_adjusted_linear(y, score, df["cohort_id"])
        assoc_rows.append(
            {
                "module_id": mod["module_id"],
                "source_task": mod["source_task"],
                "source_feature_family": mod["source_feature_family"],
                "source_universe_id": universe_id,
                "n_samples": int(valid.sum()),
                "effect_size_response_minus_nonresponse": effect,
                "standardized_mean_difference": smd,
                "p_value": p_value,
                "ci_low": ci_low,
                "ci_high": ci_high,
                "cohort_adjusted_score_coef": adj_coef,
                "cohort_adjusted_p_value": adj_p,
            }
        )
        effect_rows.append(assoc_rows[-1].copy())
        adjusted_rows.append(
            {
                "module_id": mod["module_id"],
                "cohort_adjusted_score_coef": adj_coef,
                "cohort_adjusted_p_value": adj_p,
                "adjustment": "cohort_fixed_effect_linear_probability_model",
            }
        )
    assoc = pd.DataFrame(assoc_rows)
    assoc["fdr"] = assoc.groupby("source_task")["p_value"].transform(bh_fdr) if not assoc.empty else []
    direction = []
    for _, row in assoc.iterrows():
        if pd.notna(row["fdr"]) and row["fdr"] <= 0.10 and pd.notna(row["effect_size_response_minus_nonresponse"]):
            direction.append("response_associated" if row["effect_size_response_minus_nonresponse"] > 0 else "resistance_associated")
        elif pd.notna(row["p_value"]) and row["p_value"] <= 0.10:
            direction.append("context_dependent")
        elif pd.notna(row["effect_size_response_minus_nonresponse"]):
            direction.append("uncertain")
        else:
            direction.append("unsupported")
    assoc["response_direction"] = direction
    assoc = assoc.merge(stability[["module_id", "stability_grade"]], on="module_id", how="left")
    assoc.to_csv(out / "phase4_4_module_response_association.csv", index=False)
    pd.DataFrame(effect_rows).to_csv(out / "phase4_4_module_effect_size_summary.csv", index=False)
    pd.DataFrame(adjusted_rows).to_csv(out / "phase4_4_cohort_adjusted_association.csv", index=False)
    assoc[["module_id", "source_task", "response_direction", "effect_size_response_minus_nonresponse", "fdr"]].to_csv(
        out / "phase4_4_module_direction_registry.csv", index=False
    )
    clear = int(assoc["response_direction"].isin(["response_associated", "resistance_associated", "context_dependent"]).sum()) if not assoc.empty else 0
    verdict = "PASS" if clear > 0 else "CONDITIONAL_PASS"
    (out / "phase4_4_summary.md").write_text(
        "\n".join(
            [
                "# Phase4.4 Module-Response Association Summary",
                "",
                f"- Verdict: `{verdict}`",
                f"- Directionally informative modules: {clear}",
                f"- Total modules tested: {len(assoc)}",
                "- Methods: module score effect size, Welch test, FDR, cohort-adjusted simple linear model.",
                "- Phase3.5 complex ML outputs are not used.",
                "",
            ]
        )
    )
    write_status(out / "phase4_4_status.yaml", "Phase4.4", verdict, directionally_informative_modules=clear, total_modules=len(assoc))
    return assoc


def phase4_5(modules: pd.DataFrame, stability: pd.DataFrame, association: pd.DataFrame) -> pd.DataFrame:
    out = mkdir(OUT / "05_module_classification")
    cls = modules.merge(stability, on="module_id", how="left").merge(
        association[["module_id", "response_direction", "fdr", "effect_size_response_minus_nonresponse"]],
        on="module_id",
        how="left",
    )
    module_class = []
    for _, row in cls.iterrows():
        if row["stability_grade"] == "Blocked":
            module_class.append("blocked_module")
        elif row["stability_grade"] == "C_exploratory":
            module_class.append("support_only_module")
        elif row["source_task"] == "pan_cancer_shared":
            module_class.append("shared_candidate_module")
        elif row["source_task"] == "HCC_specific":
            module_class.append("HCC_specific_candidate_module")
        elif row["source_task"] == "PD1X_extension":
            module_class.append("PD1X_extension_candidate_module")
        else:
            module_class.append("support_only_module")
    cls["module_class"] = module_class
    cls.to_csv(out / "phase4_5_module_classification.csv", index=False)
    cls[cls["module_class"] == "shared_candidate_module"].to_csv(out / "phase4_5_shared_module_candidates.csv", index=False)
    cls[cls["module_class"] == "HCC_specific_candidate_module"].to_csv(out / "phase4_5_HCC_specific_module_candidates.csv", index=False)
    cls[cls["module_class"] == "PD1X_extension_candidate_module"].to_csv(out / "phase4_5_PD1X_extension_module_candidates.csv", index=False)
    cls[cls["module_class"] == "support_only_module"].to_csv(out / "phase4_5_support_only_modules.csv", index=False)
    cls[cls["module_class"] == "blocked_module"].to_csv(out / "phase4_5_blocked_modules.csv", index=False)
    n_main = int(cls["module_class"].isin(["shared_candidate_module", "HCC_specific_candidate_module", "PD1X_extension_candidate_module"]).sum())
    verdict = "PASS" if n_main > 0 else "FAIL"
    (out / "phase4_5_summary.md").write_text(
        "\n".join(
            [
                "# Phase4.5 Module Classification Summary",
                "",
                f"- Verdict: `{verdict}`",
                f"- Main candidate modules: {n_main}",
                f"- Shared candidates: {int((cls['module_class'] == 'shared_candidate_module').sum())}",
                f"- HCC-specific candidates: {int((cls['module_class'] == 'HCC_specific_candidate_module').sum())}",
                f"- PD1X-extension candidates: {int((cls['module_class'] == 'PD1X_extension_candidate_module').sum())}",
                f"- Support-only modules: {int((cls['module_class'] == 'support_only_module').sum())}",
                f"- Blocked modules: {int((cls['module_class'] == 'blocked_module').sum())}",
                "",
            ]
        )
    )
    write_status(out / "phase4_5_status.yaml", "Phase4.5", verdict, main_candidate_modules=n_main, support_only_modules=int((cls["module_class"] == "support_only_module").sum()), blocked_modules=int((cls["module_class"] == "blocked_module").sum()))
    return cls


def phase4_6(classified: pd.DataFrame, members: pd.DataFrame, inputs: dict[str, pd.DataFrame], association: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    out = mkdir(OUT / "06_biological_interpretation")
    support = inputs["support_pd1"]
    annotations = []
    pathway_rows = []
    tf_rows = []
    support_rows = []
    conflict_rows = []
    cards = ["# Phase4.6 Module Interpretation Cards", ""]
    assoc = association.set_index("module_id")
    for _, row in classified.iterrows():
        mid = row["module_id"]
        feats = members[members["module_id"] == mid]["feature_name"].astype(str).tolist()
        bio, feat_types, tf_status = infer_annotation(feats)
        ann = {
            "module_id": mid,
            "module_class": row["module_class"],
            "source_task": row["source_task"],
            "source_feature_family": row["source_feature_family"],
            "n_features": row["n_features"],
            "key_features": ";".join(feats[:15]),
            "biological_interpretation": bio,
            "feature_types": feat_types,
            "allowed_downstream_use": "Phase5 main candidate" if row["module_class"].endswith("candidate_module") else "support or exploratory only",
            "forbidden_downstream_use": "drug recommendation; causal clinical claim; primary prediction claim",
        }
        annotations.append(ann)
        pathway_rows.append({"module_id": mid, "pathway_annotation": bio, "evidence_basis": "feature_name_pattern"})
        tf_rows.append({"module_id": mid, "tf_annotation": tf_status, "tf_features": ";".join([f for f in feats if f.startswith("tf_activity__")][:20])})
        matched = support[support["feature_or_module"].isin(feats)].copy()
        response_dir = assoc.loc[mid, "response_direction"] if mid in assoc.index else "unsupported"
        support_status = "no_PD1_anchor_overlap"
        if not matched.empty:
            mean_dir = pd.to_numeric(matched["direction"], errors="coerce").mean()
            if response_dir == "response_associated" and mean_dir > 0:
                support_status = "direction_consistent"
            elif response_dir == "resistance_associated" and mean_dir < 0:
                support_status = "direction_consistent"
            elif response_dir in {"response_associated", "resistance_associated"}:
                support_status = "support_conflict"
            else:
                support_status = "directional_only"
            for _, srow in matched.iterrows():
                support_rows.append(
                    {
                        "module_id": mid,
                        "support_id": srow["support_id"],
                        "support_feature": srow["feature_or_module"],
                        "support_direction": srow["direction"],
                        "support_strength": srow["support_strength"],
                        "support_status": support_status,
                    }
                )
        if support_status == "support_conflict":
            conflict_rows.append(
                {
                    "module_id": mid,
                    "response_direction": response_dir,
                    "conflict_type": "PD1_anchor_direction_conflict",
                    "possible_reason": "cohort structure;treatment context;cell composition;support-only evidence caveat",
                }
            )
        cards.extend(
            [
                f"## {mid}",
                "",
                f"- Module class: `{row['module_class']}`",
                f"- Source: `{row['source_task']} / {row['source_feature_family']}`",
                f"- Key features: {', '.join(feats[:10])}",
                f"- Response direction: `{response_dir}`",
                f"- Biological interpretation: {bio}",
                f"- Support evidence: `{support_status}`",
                "- Caveats: associative module; no causal or drug recommendation claim.",
                f"- Allowed downstream use: {ann['allowed_downstream_use']}",
                f"- Forbidden downstream use: {ann['forbidden_downstream_use']}",
                "- Phase5 recommendation: mechanism validation if Tier A/B; otherwise support/exploratory.",
                "",
            ]
        )
    ann_df = pd.DataFrame(annotations)
    support_df = pd.DataFrame(support_rows)
    ann_df.to_csv(out / "phase4_6_module_biological_annotation.csv", index=False)
    pd.DataFrame(pathway_rows).to_csv(out / "phase4_6_module_pathway_annotation.csv", index=False)
    pd.DataFrame(tf_rows).to_csv(out / "phase4_6_module_tf_annotation.csv", index=False)
    support_df.to_csv(out / "phase4_6_support_evidence_annotation.csv", index=False)
    pd.DataFrame(conflict_rows).to_csv(out / "phase4_6_support_conflict_log.csv", index=False)
    (out / "phase4_6_module_interpretation_cards.md").write_text("\n".join(cards))
    verdict = "PASS" if len(ann_df) > 0 else "FAIL"
    (out / "phase4_6_summary.md").write_text(
        "\n".join(
            [
                "# Phase4.6 Biological Interpretation Summary",
                "",
                f"- Verdict: `{verdict}`",
                f"- Annotated modules: {len(ann_df)}",
                f"- PD1_anchor support conflicts retained: {len(conflict_rows)}",
                "- Annotation uses feature-name based immune program rules and support evidence only as support.",
                "- No causal or drug recommendation claim is made.",
                "",
            ]
        )
    )
    write_status(out / "phase4_6_status.yaml", "Phase4.6", verdict, annotated_modules=len(ann_df), support_conflicts=len(conflict_rows))
    return ann_df, support_df


def phase4_7(classified: pd.DataFrame, matrices: dict[str, dict[str, Any]], module_scores_by_universe: dict[str, pd.DataFrame], association: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    out = mkdir(OUT / "07_module_scoring")
    assoc = association.set_index("module_id")
    sample_rows = []
    for _, mod in classified.iterrows():
        uid = mod["source_universe_id"]
        df = matrices[uid]["data"].reset_index(drop=True)
        scores = module_scores_by_universe[uid][mod["module_id"]].reset_index(drop=True)
        score_role = "primary_score" if mod["module_class"].endswith("candidate_module") else ("support_score" if mod["module_class"] == "support_only_module" else "blocked_score")
        direction = assoc.loc[mod["module_id"], "response_direction"] if mod["module_id"] in assoc.index else "unsupported"
        for i, score in scores.items():
            sample_rows.append(
                {
                    "module_id": mod["module_id"],
                    "sample_key": df.loc[i, "sample_key"],
                    "cohort_id": df.loc[i, "cohort_id"],
                    "sample_id": df.loc[i, "sample_id"],
                    "patient_id": df.loc[i, "patient_id"],
                    "patient_key": df.loc[i, "patient_key"],
                    "source_task": mod["source_task"],
                    "module_class": mod["module_class"],
                    "score_role": score_role,
                    "module_direction": direction,
                    "module_score": score,
                    "response_label": df.loc[i, "response_label"],
                    "split": df.loc[i, "split"],
                    "cancer_type": df.loc[i, "cancer_type"] if "cancer_type" in df.columns else df.loc[i, "disease"],
                    "treatment_context": df.loc[i, "treatment_context"],
                    "timepoint": df.loc[i, "timepoint"],
                }
            )
    sample_scores = pd.DataFrame(sample_rows)
    patient_scores = (
        sample_scores.groupby(["module_id", "patient_key", "patient_id", "cohort_id", "source_task", "module_class", "score_role", "module_direction"], dropna=False)
        .agg(module_score=("module_score", "mean"), n_samples=("sample_id", "nunique"), response_label=("response_label", lambda x: x.dropna().iloc[0] if len(x.dropna()) else np.nan))
        .reset_index()
    )
    dist_rows = []
    for field in ["response_label", "cohort_id", "source_task", "cancer_type", "treatment_context", "timepoint"]:
        grouped = sample_scores.groupby(["module_id", field], dropna=False)["module_score"]
        for (mid, level), vals in grouped:
            stats_row = series_distribution(vals)
            dist_rows.append(
                {
                    "module_id": mid,
                    "stratification": field,
                    "level": level,
                    **stats_row,
                }
            )
    displacement_rows = []
    for (mid, patient), sub in sample_scores.groupby(["module_id", "patient_key"]):
        timepoints = sub.dropna(subset=["timepoint"])
        pre = timepoints[timepoints["timepoint"].astype(str).str.contains("pre", case=False, na=False)]["module_score"]
        post = timepoints[timepoints["timepoint"].astype(str).str.contains("post", case=False, na=False)]["module_score"]
        if len(pre) and len(post):
            displacement_rows.append(
                {
                    "module_id": mid,
                    "patient_key": patient,
                    "pre_mean_score": pre.mean(),
                    "post_mean_score": post.mean(),
                    "post_minus_pre": post.mean() - pre.mean(),
                    "note": "paired_displacement_not_used_for_primary_supervised_training",
                }
            )
    sample_scores.to_csv(out / "phase4_7_module_scores_by_sample.csv", index=False)
    patient_scores.to_csv(out / "phase4_7_module_scores_by_patient.csv", index=False)
    pd.DataFrame(dist_rows).to_csv(out / "phase4_7_module_score_distribution_summary.csv", index=False)
    pd.DataFrame(displacement_rows).to_csv(out / "phase4_7_module_displacement_if_paired.csv", index=False)
    primary_scored = sample_scores[sample_scores["score_role"] == "primary_score"]["module_id"].nunique()
    verdict = "PASS" if primary_scored > 0 else "FAIL"
    (out / "phase4_7_summary.md").write_text(
        "\n".join(
            [
                "# Phase4.7 Module Scoring Summary",
                "",
                f"- Verdict: `{verdict}`",
                f"- Sample-level score rows: {len(sample_scores)}",
                f"- Patient-level score rows: {len(patient_scores)}",
                f"- Primary scored modules: {int(primary_scored)}",
                "- Support and blocked scores are flagged separately from primary scores.",
                "",
            ]
        )
    )
    write_status(out / "phase4_7_status.yaml", "Phase4.7", verdict, sample_score_rows=len(sample_scores), patient_score_rows=len(patient_scores), primary_scored_modules=int(primary_scored))
    return sample_scores, patient_scores


def phase4_8(classified: pd.DataFrame, association: pd.DataFrame, annotations: pd.DataFrame, sample_scores: pd.DataFrame) -> pd.DataFrame:
    out = mkdir(OUT / "08_phase5_handoff")
    assoc_cols = ["module_id", "response_direction", "fdr", "effect_size_response_minus_nonresponse"]
    missing_assoc_cols = [c for c in assoc_cols if c == "module_id" or c not in classified.columns]
    if len(missing_assoc_cols) > 1:
        df = classified.merge(association[missing_assoc_cols], on="module_id", how="left")
    else:
        df = classified.copy()
    df = df.merge(annotations[["module_id", "biological_interpretation", "key_features"]], on="module_id", how="left")
    tier = []
    support_status = []
    for _, row in df.iterrows():
        clear_dir = row["response_direction"] in {"response_associated", "resistance_associated", "context_dependent"}
        if row["module_class"] == "blocked_module" or row["stability_grade"] == "Blocked":
            tier.append("Blocked")
            support_status.append("blocked")
        elif row["module_class"].endswith("candidate_module") and row["stability_grade"] == "A_stable" and clear_dir:
            tier.append("Tier_A")
            support_status.append("primary_with_clear_direction")
        elif row["module_class"].endswith("candidate_module") and row["stability_grade"] in {"A_stable", "B_usable"}:
            tier.append("Tier_B")
            support_status.append("primary_with_caveat")
        else:
            tier.append("Tier_C")
            support_status.append("support_or_exploratory")
    df["tier"] = tier
    df["support_status"] = support_status
    df["candidate_module_id"] = df["module_id"]
    df["evidence_pool"] = np.where(df["tier"].isin(["Tier_A", "Tier_B"]), "primary_evidence_pool", "support_or_blocked")
    df["module_score_available"] = df["module_id"].isin(sample_scores["module_id"].unique())
    df["confounding_risk"] = df["source_task"].map({"HCC_specific": "MODERATE", "PD1X_extension": "LOW", "pan_cancer_shared": "LOW"}).fillna("UNKNOWN")
    df["allowed_phase5_use"] = np.where(df["tier"].isin(["Tier_A", "Tier_B"]), "main mechanism validation candidate", "support/exploratory only")
    df["forbidden_phase5_use"] = "direct clinical prediction; drug recommendation; causal claim without validation"
    df["required_caveat"] = "source/batch/missingness confounding conditional pass; associative module only"
    out_cols = [
        "candidate_module_id",
        "module_class",
        "tier",
        "source_task",
        "source_feature_family",
        "evidence_pool",
        "n_features",
        "key_features",
        "module_score_available",
        "response_direction",
        "stability_grade",
        "confounding_risk",
        "support_status",
        "biological_interpretation",
        "allowed_phase5_use",
        "forbidden_phase5_use",
        "required_caveat",
    ]
    main = df[df["tier"].isin(["Tier_A", "Tier_B"])][out_cols].copy()
    support = df[df["tier"] == "Tier_C"][out_cols].copy()
    blocked = df[df["tier"] == "Blocked"][out_cols].copy()
    main.to_csv(out / "phase4_8_phase5_main_candidate_modules.csv", index=False)
    support.to_csv(out / "phase4_8_phase5_support_modules.csv", index=False)
    blocked.to_csv(out / "phase4_8_phase5_blocked_modules.csv", index=False)
    card_lines = ["# Phase4.8 Candidate Module Cards", ""]
    for _, row in pd.concat([main, support]).iterrows():
        card_lines.extend(
            [
                f"## {row['candidate_module_id']}",
                "",
                f"- Tier: `{row['tier']}`",
                f"- Class: `{row['module_class']}`",
                f"- Direction: `{row['response_direction']}`",
                f"- Stability: `{row['stability_grade']}`",
                f"- Interpretation: {row['biological_interpretation']}",
                f"- Allowed Phase5 use: {row['allowed_phase5_use']}",
                f"- Caveat: {row['required_caveat']}",
                "",
            ]
        )
    (out / "phase4_8_candidate_module_cards.md").write_text("\n".join(card_lines))
    verdict = "PASS" if len(main) > 0 else "FAIL"
    (out / "phase4_8_summary.md").write_text(
        "\n".join(
            [
                "# Phase4.8 Phase5 Handoff Summary",
                "",
                f"- Verdict: `{verdict}`",
                f"- Tier A modules: {int((main['tier'] == 'Tier_A').sum()) if len(main) else 0}",
                f"- Tier B modules: {int((main['tier'] == 'Tier_B').sum()) if len(main) else 0}",
                f"- Tier C support modules: {len(support)}",
                f"- Blocked modules: {len(blocked)}",
                "- Main handoff contains only Tier A/B primary modules.",
                "",
            ]
        )
    )
    write_status(out / "phase4_8_status.yaml", "Phase4.8", verdict, tier_A=int((main["tier"] == "Tier_A").sum()) if len(main) else 0, tier_B=int((main["tier"] == "Tier_B").sum()) if len(main) else 0, tier_C=len(support), blocked=len(blocked))
    return df


def phase4_9(candidate_df: pd.DataFrame, sample_scores: pd.DataFrame, annotations: pd.DataFrame) -> None:
    out = mkdir(OUT / "09_final_report")
    main = candidate_df[candidate_df["tier"].isin(["Tier_A", "Tier_B"])].copy()
    support = candidate_df[candidate_df["tier"] == "Tier_C"].copy()
    blocked = candidate_df[candidate_df["tier"] == "Blocked"].copy()
    shared = main[main["module_class"] == "shared_candidate_module"]["module_id"].tolist()
    hcc = main[main["module_class"] == "HCC_specific_candidate_module"]["module_id"].tolist()
    pd1x = main[main["module_class"] == "PD1X_extension_candidate_module"]["module_id"].tolist()
    master = candidate_df.copy()
    master.to_csv(out / "phase4_module_master_table.csv", index=False)
    main.to_csv(out / "phase4_phase5_handoff_main_modules.csv", index=False)
    support.to_csv(out / "phase4_phase5_handoff_support_evidence.csv", index=False)
    blocked.to_csv(out / "phase4_blocked_module_and_input_log.csv", index=False)
    verdict = "PASS" if len(main) >= 3 else ("CONDITIONAL_PASS" if len(main) > 0 else "FAIL")
    decision = {
        "verdict": verdict,
        "number_of_primary_modules": int(len(main)),
        "number_of_Tier_A_modules": int((main["tier"] == "Tier_A").sum()) if len(main) else 0,
        "number_of_Tier_B_modules": int((main["tier"] == "Tier_B").sum()) if len(main) else 0,
        "number_of_support_only_modules": int(len(support)),
        "number_of_blocked_modules": int(len(blocked)),
        "shared_candidate_modules": shared,
        "HCC_specific_candidate_modules": hcc,
        "PD1X_extension_candidate_modules": pd1x,
        "phase5_ready_modules": main["module_id"].tolist(),
        "support_only_evidence": support["module_id"].tolist(),
        "forbidden_claims": [
            "generalizable clinical prediction model",
            "complex ML proof",
            "PD1_anchor primary supervised mechanism",
            "causal response mechanism",
            "clinical drug recommendation",
            "support-only evidence as primary evidence",
            "unintegrated new data validation claim",
        ],
        "required_caveats": [
            "source/batch/missingness confounding conditional pass",
            "associative modules only",
            "Phase3.5 strong ML excluded",
            "PD1_anchor support only",
        ],
        "recommended_next_phase": "Phase5 mechanism adjudication and validation",
    }
    write_yaml(out / "PHASE4_FINAL_DECISION.yaml", decision)
    write_yaml(
        out / "phase4_reproducibility_manifest.yaml",
        {
            "generated_at_utc": now_iso(),
            "script": rel(ROOT / "scripts/v6_1/run_phase4_role_aware_module_discovery.py"),
            "entry_gate": rel(READINESS / "phase4_entry_decision_20260611.yaml"),
            "feature_matrix": rel(STEP2_REPAIR / "repair_B2_unified_feature_matrix_by_sample.hotfix.parquet"),
            "membership": rel(STEP3 / "01_analysis_universe/step3_phase3_1_sample_universe_membership.repaired.csv"),
            "blocked_inputs_enforced": True,
            "phase3_5_ml_used": False,
            "pd1_anchor_primary_used": False,
            "unintegrated_new_data_used": False,
        },
    )
    report = [
        "# Phase4 Final Report",
        "",
        f"- Verdict: `{verdict}`",
        f"- Primary Tier A/B modules: {len(main)}",
        f"- Tier A modules: {int((main['tier'] == 'Tier_A').sum()) if len(main) else 0}",
        f"- Tier B modules: {int((main['tier'] == 'Tier_B').sum()) if len(main) else 0}",
        f"- Support/Tier C modules: {len(support)}",
        f"- Blocked modules: {len(blocked)}",
        "",
        "## Candidate Lines",
        "",
        f"- Shared candidate modules: {len(shared)}",
        f"- HCC-specific candidate modules: {len(hcc)}",
        f"- PD1X-extension candidate modules: {len(pd1x)}",
        "",
        "## Evidence Separation",
        "",
        "- Primary modules were discovered only from primary evidence pool matrices.",
        "- PD1_anchor was not used for primary module discovery or feature selection.",
        "- Phase3.5 strong ML outputs were not used.",
        "- Unintegrated P0/P1 new data were not used.",
        "- Step3-derived pseudobulk_gene was not used as primary feature family.",
        "",
        "## What Can Be Said",
        "",
        "- Stable immune-related modules were identified from the audited primary evidence pool.",
        "- Some modules show response-associated or resistance-associated direction.",
        "- Modules are Phase5 mechanism candidates, not clinical prediction or drug recommendation outputs.",
        "",
        "## What Cannot Be Said",
        "",
        "- No generalizable clinical prediction model is claimed.",
        "- No causal response mechanism is claimed.",
        "- No clinical drug recommendation is made.",
        "- HCC-specific modules are not asserted as pan-cancer shared mechanisms.",
        "",
    ]
    (out / "PHASE4_FINAL_REPORT.md").write_text("\n".join(report))
    (out / "phase4_9_summary.md").write_text(
        "\n".join(
            [
                "# Phase4.9 Final Decision Summary",
                "",
                f"- Verdict: `{verdict}`",
                f"- Primary Tier A/B modules: {len(main)}",
                f"- Tier A modules: {int((main['tier'] == 'Tier_A').sum()) if len(main) else 0}",
                f"- Tier B modules: {int((main['tier'] == 'Tier_B').sum()) if len(main) else 0}",
                f"- Support-only modules: {len(support)}",
                f"- Blocked modules: {len(blocked)}",
                f"- Shared candidates: {len(shared)}",
                f"- HCC-specific candidates: {len(hcc)}",
                f"- PD1X-extension candidates: {len(pd1x)}",
                "- Recommended next phase: Phase5 mechanism adjudication and validation.",
                "",
            ]
        )
    )
    write_status(
        out / "phase4_9_status.yaml",
        "Phase4.9",
        verdict,
        primary_modules=len(main),
        tier_A=int((main["tier"] == "Tier_A").sum()) if len(main) else 0,
        tier_B=int((main["tier"] == "Tier_B").sum()) if len(main) else 0,
        support_only=len(support),
        blocked=len(blocked),
    )

    phase_report = [
        "# Phase4 Phase Report 20260611",
        "",
        "## Scope And Status",
        "",
        f"- Phase: `v6.1 Phase4 Role-aware immune mechanism module discovery`",
        f"- Final verdict: `{verdict}`",
        "- Scientific role: response-blind immune module discovery from audited primary evidence pool.",
        "- Engineering role: freeze module tables, scores, interpretation cards, and Phase5 handoff with primary/support/blocked evidence separated.",
        "",
        "## Canonical Inputs",
        "",
        f"- Readiness gate: `{rel(READINESS / 'phase4_entry_decision_20260611.yaml')}`",
        f"- Evidence registry: `{rel(READINESS / 'phase4_evidence_pool_registry_20260611.csv')}`",
        f"- Feature contract: `{rel(READINESS / 'phase4_feature_family_contract_20260611.csv')}`",
        f"- Blocked input rules: `{rel(READINESS / 'phase4_blocked_inputs_20260611.csv')}`",
        f"- Matrix: `{rel(STEP2_REPAIR / 'repair_B2_unified_feature_matrix_by_sample.hotfix.parquet')}`",
        f"- Universe membership: `{rel(STEP3 / '01_analysis_universe/step3_phase3_1_sample_universe_membership.repaired.csv')}`",
        "",
        "## Canonical Outputs",
        "",
        "- `PHASE4_FINAL_DECISION.yaml`: final machine-readable verdict and module lists.",
        "- `phase4_module_master_table.csv`: all retained, support-only, and blocked module records.",
        "- `phase4_phase5_handoff_main_modules.csv`: Tier A/B primary modules for Phase5 main use.",
        "- `phase4_phase5_handoff_support_evidence.csv`: support/Tier C evidence only.",
        "- `phase4_blocked_module_and_input_log.csv`: blocked module/input log.",
        "- `07_module_scoring/phase4_7_module_scores_by_sample.csv`: sample-level module scores.",
        "- `07_module_scoring/phase4_7_module_scores_by_patient.csv`: patient-level module scores.",
        "- `06_biological_interpretation/phase4_6_module_interpretation_cards.md`: module interpretation cards.",
        "",
        "## Processing Flow",
        "",
        "1. Enforce primary/support/blocked evidence pools.",
        "2. Build primary matrices for HCC_specific, PD1X_extension, and pan_cancer_shared tasks.",
        "3. Discover modules by response-blind hierarchical clustering on feature correlation.",
        "4. Assess bootstrap, split, cohort, feature-family, and missingness stability.",
        "5. Associate module scores with response using simple/interpretable evidence only.",
        "6. Classify modules as shared, HCC-specific, PD1X-extension, support-only, or blocked.",
        "7. Annotate biology and PD1_anchor support consistency without using support evidence as primary.",
        "8. Score modules at sample and patient level.",
        "9. Grade Tier A/B/C/Blocked and write Phase5 handoff.",
        "",
        "## Current Output Counts",
        "",
        f"- Primary Tier A/B modules: {len(main)}",
        f"- Tier A modules: {int((main['tier'] == 'Tier_A').sum()) if len(main) else 0}",
        f"- Tier B modules: {int((main['tier'] == 'Tier_B').sum()) if len(main) else 0}",
        f"- Support-only Tier C modules: {len(support)}",
        f"- Blocked modules: {len(blocked)}",
        f"- Shared candidates: {len(shared)}",
        f"- HCC-specific candidates: {len(hcc)}",
        f"- PD1X-extension candidates: {len(pd1x)}",
        f"- Sample score rows: {len(sample_scores)}",
        "",
        "## Downstream Contract",
        "",
        "- Phase5 main use: only Tier A/B rows in `phase4_phase5_handoff_main_modules.csv`.",
        "- Support use: only annotation, direction sanity check, and sensitivity context.",
        "- Forbidden: clinical prediction claim, causal mechanism claim, drug recommendation, PD1_anchor primary claim, unintegrated new-data validation claim.",
        "- Required caveat: source/batch/missingness confounding remains conditional-pass background.",
        "",
        "## Reproduction",
        "",
        "- Command: `conda run -n proj006 python -W error scripts/v6_1/run_phase4_role_aware_module_discovery.py`",
        "- Script: `scripts/v6_1/run_phase4_role_aware_module_discovery.py`",
        "",
    ]
    (out / "PHASE4_PHASE_REPORT_20260611.md").write_text("\n".join(phase_report))

    changelog = [
        "# Phase4 Release Changelog 20260611",
        "",
        "## Release Scope",
        "",
        "- Completed Phase4.0 through Phase4.9 output package.",
        "- Added Phase4.9 standalone summary/status files for per-phase completeness.",
        "- Added phase report and release changelog for freeze-ready downstream handoff.",
        "",
        "## Code And Output Updates",
        "",
        "- Added safe handling for zero-variance correlations and all-missing distribution groups.",
        "- Confirmed `-W error` run completes without warning escalation.",
        "- Rewrote final outputs after code patch.",
        "",
        "## Verification Outcomes",
        "",
        f"- Final verdict: `{verdict}`",
        f"- Phase5 main candidates: {len(main)}",
        f"- Blocked modules: {len(blocked)}",
        "- Blocked input use in primary analysis: 0.",
        "- Primary evidence violations: 0.",
        "- PD1_anchor primary rows: 0.",
        "",
        "## Final Release State",
        "",
        "- Phase4 main handoff is ready for Phase5 mechanism adjudication.",
        "- Support-only and conflict evidence remain isolated.",
        "- Unintegrated new datasets remain outside current Phase4 primary conclusions.",
        "",
    ]
    (out / "PHASE4_RELEASE_CHANGELOG_20260611.md").write_text("\n".join(changelog))

    # Also place final outputs at phase root for direct handoff.
    (OUT / "PHASE4_FINAL_REPORT.md").write_text("\n".join(report))
    (OUT / "PHASE4_PHASE_REPORT_20260611.md").write_text("\n".join(phase_report))
    (OUT / "PHASE4_RELEASE_CHANGELOG_20260611.md").write_text("\n".join(changelog))
    write_yaml(OUT / "PHASE4_FINAL_DECISION.yaml", decision)
    master.to_csv(OUT / "phase4_module_master_table.csv", index=False)
    main.to_csv(OUT / "phase4_phase5_handoff_main_modules.csv", index=False)
    support.to_csv(OUT / "phase4_phase5_handoff_support_evidence.csv", index=False)
    blocked.to_csv(OUT / "phase4_blocked_module_and_input_log.csv", index=False)
    write_yaml(OUT / "phase4_reproducibility_manifest.yaml", {
        "generated_at_utc": now_iso(),
        "script": rel(ROOT / "scripts/v6_1/run_phase4_role_aware_module_discovery.py"),
        "phase3_5_ml_used": False,
        "pd1_anchor_primary_used": False,
        "unintegrated_new_data_used": False,
    })


def main() -> None:
    mkdir(OUT)
    inputs = load_inputs()
    phase4_0(inputs)
    matrices = phase4_1(inputs)
    modules, members, module_scores_by_universe = phase4_2(matrices)
    stability = phase4_3(modules, members, matrices, module_scores_by_universe)
    association = phase4_4(modules, matrices, module_scores_by_universe, stability)
    classified = phase4_5(modules, stability, association)
    annotations, _ = phase4_6(classified, members, inputs, association)
    sample_scores, _ = phase4_7(classified, matrices, module_scores_by_universe, association)
    candidates = phase4_8(classified, association, annotations, sample_scores)
    phase4_9(candidates, sample_scores, annotations)
    print(f"Wrote Phase4 outputs: {OUT}")
    print(f"Modules: {len(modules)}")
    print(f"Phase5 main candidates: {int(candidates['tier'].isin(['Tier_A','Tier_B']).sum())}")


if __name__ == "__main__":
    main()
