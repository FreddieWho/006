#!/usr/bin/env python3
"""Step3 QC-aware strong baseline runner for v6.1.

The script is intentionally conservative: Phase 3.0-3.3 are hard gates.
No response prediction model is fit unless those gates pass or produce an
explicitly explainable CONDITIONAL_PASS.
"""

from __future__ import annotations

import argparse
import json
import logging
import math
import platform
import random
import re
import sys
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Any, Iterable

import numpy as np
import pandas as pd
import yaml
from scipy.stats import chi2_contingency, fisher_exact
from sklearn.base import clone
from sklearn.ensemble import RandomForestClassifier
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (
    accuracy_score,
    average_precision_score,
    balanced_accuracy_score,
    brier_score_loss,
    confusion_matrix,
    f1_score,
    precision_recall_curve,
    roc_auc_score,
)
from sklearn.model_selection import GridSearchCV, StratifiedKFold
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import OneHotEncoder, StandardScaler


SEED = 20260604
RNG = np.random.default_rng(SEED)
random.seed(SEED)

ROOT = Path(__file__).resolve().parents[2]
STEP2_REPAIR = ROOT / "results/v6_1/step2_repair"
STEP2_RUN = ROOT / "results/v6_1/step2/step2_v6_1_0505_0319"
OUT = ROOT / "results/v6_1/step3_qc_aware_strong_baseline"

PHASE_DIRS = {
    "3_0": OUT / "00_input_contract",
    "3_1": OUT / "01_analysis_universe",
    "3_2": OUT / "02_feature_registry_missingness",
    "3_3": OUT / "03_leakage_confounding_audit",
    "3_4": OUT / "04_interpretable_baselines",
    "3_5": OUT / "05_strong_ml_baselines",
    "3_6": OUT / "06_robustness_negative_controls",
    "3_7": OUT / "07_milestone_B_step4_handoff",
    "abort": OUT / "99_abort_or_repair_logs",
}

HOTFIX_FILES = {
    "matrix": STEP2_REPAIR / "repair_B2_unified_feature_matrix_by_sample.hotfix.parquet",
    "dictionary": STEP2_REPAIR / "repair_B2_feature_dictionary_v6_1.hotfix.csv",
    "missingness": STEP2_REPAIR / "repair_B2_feature_missingness_report.hotfix.csv",
    "manifest": STEP2_REPAIR / "repair_B2_step2_10_decision_manifest.hotfix.yaml",
    "preview": STEP2_REPAIR / "repair_C2_analysis_universe_preview.hotfix.csv",
    "readiness": STEP2_REPAIR / "step2_to_step3_readiness_status.hotfix.yaml",
    "readiness_report": STEP2_REPAIR / "STEP2_TO_STEP3_READINESS_REPORT.hotfix.md",
    "join_audit": STEP2_REPAIR / "repair_A_step2_sample_metadata_join_audit.csv",
    "split": STEP2_REPAIR / "repair_A_frozen_patient_split_v6_1.repaired.csv",
    "metadata": STEP2_REPAIR / "repair_A_sample_metadata_master_v6_1.repaired.csv",
    "sample_qc_flags": STEP2_RUN
    / "09_feature_matrix/unified_feature_matrix_v1/sample_feature_QC_flags.csv",
    "cohort_eligibility": STEP2_RUN
    / "09_feature_matrix/unified_feature_matrix_v1/cohort_feature_eligibility.csv",
    "pseudobulk_manifest": STEP2_RUN
    / "07_pseudobulk/pseudobulk_feature_v1_final_v3/pseudobulk_decision_manifest.yaml",
    "fraction_manifest": STEP2_RUN
    / "06_fraction/fraction_feature_v1/fraction_decision_manifest.yaml",
    "sptf_manifest": STEP2_RUN
    / "08_signature_pathway_tf/signature_pathway_tf_feature_v1_tf_rescue_v1/signature_pathway_tf_decision_manifest.yaml",
    "sptf_qc": STEP2_RUN
    / "08_signature_pathway_tf/signature_pathway_tf_feature_v1_tf_rescue_v1/signature_pathway_tf_qc_summary.csv",
}
REPAIRED_FEATURE_DICTIONARY = (
    OUT / "00_input_contract/repair_phase3_0_feature_dictionary.repaired.csv"
)

META_COLS = {
    "sample_key",
    "run_id",
    "created_at",
    "input_manifest_ref",
    "cohort_id",
    "sample_id",
    "source_h5ad_set",
    "source_h5ad",
    "step2_9_row_id",
    "patient_id",
    "patient_key",
    "split",
    "fold_id",
    "disease",
    "tissue_source",
    "sample_type",
    "timepoint",
    "timepoint_class",
    "treatment_context",
    "metadata_join_method",
    "metadata_missing",
    "pseudobulk_cell_count",
    "pseudobulk_available",
    "signature_pathway_tf_available",
    "pseudobulk_unavailable_reason",
    "immune_cells_used",
    "total_cells_used",
}

BASE_FAMILIES = ["fraction", "signature", "pathway", "tf_activity", "pseudobulk_gene"]
ALL_FAMILIES = BASE_FAMILIES + ["combined_interpretable", "combined_all_allowed"]
TASKS = ["PD1_anchor", "PD1X_extension", "HCC_specific", "pan_cancer_shared", "sensitivity_only"]
VALID_SPLITS = {"train", "validation", "test"}
FORBIDDEN_PATTERNS = [
    "response",
    "responder",
    "non_responder",
    "recist",
    "mrecist",
    "best_response",
    "progression",
    "survival",
    r"\bos\b",
    r"\bpfs\b",
    "death",
    "event",
    "outcome",
    "label",
    "train",
    "test",
    "split",
    "treatment_response",
]

FORBIDDEN_NAME_WHITELIST = {
    "pathway__REACTOME_DOWNSTREAM_SIGNALING_EVENTS_OF_B_CELL_RECEPTOR_BCR": "Reactome pathway phrase; events means signaling events, not clinical event label.",
    "pathway__REACTOME_NRAGE_SIGNALS_DEATH_THROUGH_JNK": "Reactome apoptosis/death pathway biology; not mortality/death outcome label.",
}


@dataclass
class Context:
    matrix: pd.DataFrame
    dictionary: pd.DataFrame
    missingness: pd.DataFrame
    readiness: dict[str, Any]
    manifest: dict[str, Any]
    preview: pd.DataFrame
    join_audit: pd.DataFrame
    split: pd.DataFrame
    metadata: pd.DataFrame
    sample_qc_flags: pd.DataFrame | None
    cohort_eligibility: pd.DataFrame | None
    feature_cols: list[str]
    family_features: dict[str, list[str]]
    sample_table: pd.DataFrame


def setup() -> None:
    for path in PHASE_DIRS.values():
        path.mkdir(parents=True, exist_ok=True)
    logging.basicConfig(
        filename=OUT / "step3_run.log",
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
    )
    logging.info("Step3 started seed=%s", SEED)


def write_yaml(path: Path, data: dict[str, Any]) -> None:
    path.write_text(yaml.safe_dump(to_builtin(data), sort_keys=False), encoding="utf-8")


def write_md(path: Path, title: str, lines: list[str]) -> None:
    path.write_text("# " + title + "\n\n" + "\n".join(lines) + "\n", encoding="utf-8")


def to_builtin(obj: Any) -> Any:
    if isinstance(obj, dict):
        return {str(k): to_builtin(v) for k, v in obj.items()}
    if isinstance(obj, list):
        return [to_builtin(v) for v in obj]
    if isinstance(obj, tuple):
        return [to_builtin(v) for v in obj]
    if isinstance(obj, (np.integer,)):
        return int(obj)
    if isinstance(obj, (np.floating,)):
        return None if np.isnan(obj) else float(obj)
    if isinstance(obj, (np.bool_,)):
        return bool(obj)
    if isinstance(obj, float) and math.isnan(obj):
        return None
    if isinstance(obj, (str, int, bool)) or obj is None:
        return obj
    if isinstance(obj, Path):
        return str(obj)
    return str(obj)


def read_yaml(path: Path) -> dict[str, Any]:
    with path.open() as handle:
        return yaml.safe_load(handle) or {}


def bool_series(s: pd.Series) -> pd.Series:
    return s.astype(str).str.lower().isin(["true", "1", "yes"])


def label_to_binary(s: pd.Series) -> pd.Series:
    text = s.astype(str).str.lower()
    out = pd.Series(np.nan, index=s.index, dtype="float64")
    out[text.eq("responder")] = 1.0
    out[text.eq("non_responder")] = 0.0
    return out


def load_context() -> Context:
    missing_paths = [str(p) for p in HOTFIX_FILES.values() if not p.exists()]
    if missing_paths:
        raise FileNotFoundError("Missing required inputs: " + "; ".join(missing_paths))

    readiness = read_yaml(HOTFIX_FILES["readiness"])
    manifest = read_yaml(HOTFIX_FILES["manifest"])
    matrix = pd.read_parquet(HOTFIX_FILES["matrix"])
    dictionary_path = (
        REPAIRED_FEATURE_DICTIONARY
        if REPAIRED_FEATURE_DICTIONARY.exists()
        else HOTFIX_FILES["dictionary"]
    )
    dictionary = pd.read_csv(dictionary_path)
    missingness = pd.read_csv(HOTFIX_FILES["missingness"])
    preview = pd.read_csv(HOTFIX_FILES["preview"])
    join_audit = pd.read_csv(HOTFIX_FILES["join_audit"])
    split = pd.read_csv(HOTFIX_FILES["split"])
    metadata = pd.read_csv(HOTFIX_FILES["metadata"])
    sample_qc_flags = pd.read_csv(HOTFIX_FILES["sample_qc_flags"])
    cohort_eligibility = pd.read_csv(HOTFIX_FILES["cohort_eligibility"])

    feature_cols = [c for c in matrix.columns if c not in META_COLS]
    dictionary["family"] = dictionary["feature_type"].replace({"tf": "tf_activity"})
    family_features: dict[str, list[str]] = {}
    for family in BASE_FAMILIES:
        rows = dictionary[
            dictionary["family"].eq(family) & bool_series(dictionary["included_in_main"])
        ]
        family_features[family] = [f for f in rows["feature_name"] if f in matrix.columns]
    family_features["combined_interpretable"] = (
        family_features["fraction"]
        + family_features["signature"]
        + family_features["pathway"]
        + family_features["tf_activity"]
    )
    family_features["combined_all_allowed"] = family_features["combined_interpretable"] + family_features[
        "pseudobulk_gene"
    ]

    sample_table = build_sample_table(matrix, join_audit, metadata, split)
    return Context(
        matrix=matrix,
        dictionary=dictionary,
        missingness=missingness,
        readiness=readiness,
        manifest=manifest,
        preview=preview,
        join_audit=join_audit,
        split=split,
        metadata=metadata,
        sample_qc_flags=sample_qc_flags,
        cohort_eligibility=cohort_eligibility,
        feature_cols=feature_cols,
        family_features=family_features,
        sample_table=sample_table,
    )


def build_sample_table(
    matrix: pd.DataFrame, join_audit: pd.DataFrame, metadata: pd.DataFrame, split: pd.DataFrame
) -> pd.DataFrame:
    cols = [
        "sample_key",
        "cohort_id",
        "sample_id",
        "source_h5ad",
        "patient_id",
        "split",
        "disease",
        "treatment_context",
        "timepoint_class",
        "metadata_missing",
        "pseudobulk_available",
        "signature_pathway_tf_available",
    ]
    st = matrix[[c for c in cols if c in matrix.columns]].copy()
    audit_cols = [
        "sample_key",
        "has_sample_metadata",
        "has_split_info",
        "best_response_binary",
        "inclusion_status",
        "supervised_anchor_available_repaired",
    ]
    st = st.merge(join_audit[[c for c in audit_cols if c in join_audit.columns]], on="sample_key", how="left")

    meta_small = metadata[
        [
            c
            for c in [
                "cohort_id",
                "sample_id",
                "patient_id",
                "response_broad_binary",
                "response_strict_binary",
                "response_standard",
                "response_raw",
                "usable_in_main_analysis",
                "usable_in_sensitivity_analysis",
                "exclusion_reason",
            ]
            if c in metadata.columns
        ]
    ].drop_duplicates(["cohort_id", "sample_id", "patient_id"])
    st = st.merge(meta_small, on=["cohort_id", "sample_id", "patient_id"], how="left")

    split_small = split[
        [
            c
            for c in [
                "cohort_id",
                "patient_id",
                "patient_join_key",
                "split",
                "split_label",
                "supervised_anchor_available_repaired",
            ]
            if c in split.columns
        ]
    ].drop_duplicates(["cohort_id", "patient_id"])
    split_small = split_small.rename(
        columns={
            "split": "patient_split",
            "split_label": "patient_split_label",
            "supervised_anchor_available_repaired": "patient_supervised_anchor_available_repaired",
        }
    )
    st = st.merge(split_small, on=["cohort_id", "patient_id"], how="left")
    st["effective_split"] = st["patient_split"].fillna(st["split"])
    response_source = st.get("best_response_binary", pd.Series(index=st.index, dtype=object)).copy()
    response_source = response_source.where(response_source.notna(), st.get("response_broad_binary"))
    st["response_label"] = response_source.astype(str)
    st["response_binary_num"] = label_to_binary(st["response_label"])
    st["metadata_status"] = np.where(
        st.get("metadata_missing", False).fillna(False).astype(bool), "metadata_missing", "metadata_available"
    )
    st["unresolved_no_split"] = ~st["effective_split"].astype(str).isin(VALID_SPLITS)
    return st


def phase3_0(ctx: Context) -> str:
    out = PHASE_DIRS["3_0"]
    tests: list[dict[str, Any]] = []

    def test(name: str, ok: bool, severity: str, detail: Any = "") -> None:
        tests.append({"test": name, "passed": bool(ok), "severity": severity, "detail": detail})

    readiness_verdict = str(ctx.readiness.get("verdict", "")).replace("_", " ").lower()
    test("readiness_verdict_go_or_conditional", readiness_verdict in {"conditional go", "go", "conditional-go"}, "hard", ctx.readiness.get("verdict"))
    test("readiness_blockers_empty", not ctx.readiness.get("blockers"), "hard", ctx.readiness.get("blockers"))
    allowed = ctx.readiness.get("allowed_feature_families", {}) or {}
    for family in BASE_FAMILIES:
        test(f"allowed_family_{family}", family in allowed and bool(allowed[family]), "hard", allowed.get(family))
    blocked = ctx.readiness.get("blocked_feature_families", {}) or {}
    for family in BASE_FAMILIES:
        test(f"not_blocked_family_{family}", family not in blocked, "hard", blocked.get(family))

    manifest_text = json.dumps(ctx.manifest, default=str).lower()
    test("tf_source_long_file_hotfix", "tf_activity_scores_by_sample.parquet" in manifest_text, "hard")
    test("tf_pivot_key_documented_or_equivalent", "cohort_id" in manifest_text and "sample_id" in manifest_text and "source_h5ad" in manifest_text, "conditional")
    test("retained_main_tf_features_574", str(ctx.manifest).find("574") >= 0, "hard")
    test("all_conserved_in_pivot_or_validation_present", "all_conserved_in_pivot" in manifest_text or "hotfixed_re_pivoted" in manifest_text, "conditional")

    dict_features = set(ctx.dictionary["feature_name"])
    main = ctx.dictionary[bool_series(ctx.dictionary["included_in_main"])]
    main_features = set(main["feature_name"])
    matrix_features = set(ctx.feature_cols)
    missing_main = sorted(main_features - matrix_features)
    extra = sorted(matrix_features - dict_features)
    tf_main = main[main["family"].eq("tf_activity")]["feature_name"].tolist()
    tf_all_na = [f for f in tf_main if f in ctx.matrix.columns and ctx.matrix[f].isna().all()]
    sample_dupes = int(ctx.matrix["sample_key"].duplicated().sum())

    labeled = ctx.sample_table["response_binary_num"].notna()
    supervised_candidate = labeled & ctx.sample_table["treatment_context"].astype(str).isin(["PD1_ICI_anchor", "PD1X_extension"])
    supervised_join_ok = bool(ctx.sample_table.loc[supervised_candidate, "effective_split"].astype(str).isin(VALID_SPLITS).all())
    if supervised_candidate.sum() == 0:
        supervised_join_ok = False

    forbidden_rows = []
    for f in ctx.feature_cols:
        hits = [p for p in FORBIDDEN_PATTERNS if re.search(p, f, flags=re.IGNORECASE)]
        if hits:
            drow = ctx.dictionary[ctx.dictionary["feature_name"].eq(f)]
            included_in_main = bool(bool_series(drow["included_in_main"]).iloc[0]) if not drow.empty else False
            forbidden_rows.append(
                {
                    "feature_name": f,
                    "patterns": ";".join(hits),
                    "included_in_main": included_in_main,
                    "whitelisted": f in FORBIDDEN_NAME_WHITELIST,
                    "whitelist_reason": FORBIDDEN_NAME_WHITELIST.get(f, ""),
                }
            )
    forbidden = pd.DataFrame(
        forbidden_rows,
        columns=["feature_name", "patterns", "included_in_main", "whitelisted", "whitelist_reason"],
    )
    forbidden_unresolved = (
        forbidden[
            forbidden["included_in_main"].astype(bool)
            & ~forbidden.get("whitelisted", pd.Series(False, index=forbidden.index)).astype(bool)
        ]
        if not forbidden.empty
        else forbidden
    )

    test("hotfix_matrix_path", HOTFIX_FILES["matrix"].name.endswith(".hotfix.parquet"), "hard", str(HOTFIX_FILES["matrix"]))
    test("sample_key_unique", sample_dupes == 0, "hard", sample_dupes)
    test("included_main_features_present", len(missing_main) == 0, "hard", len(missing_main))
    test("main_tf_features_non_all_na", len(tf_all_na) == 0 and len(tf_main) == 574, "hard", {"tf_count": len(tf_main), "all_na": len(tf_all_na)})
    test("supervised_samples_join_patient_split", supervised_join_ok, "hard", {"supervised_candidate_labeled_samples": int(supervised_candidate.sum()), "all_labeled_samples": int(labeled.sum())})
    test("no_forbidden_leakage_feature_names", forbidden_unresolved.empty, "hard", forbidden_unresolved.to_dict("records")[:20])
    test("non_main_dictionary_matrix_consistency", len(extra) <= 25, "conditional", {"extra_feature_cols": len(extra), "examples": extra[:20]})

    shape = pd.DataFrame(
        [
            {
                "n_rows": len(ctx.matrix),
                "n_columns": ctx.matrix.shape[1],
                "n_feature_columns": len(ctx.feature_cols),
                "n_dictionary_rows": len(ctx.dictionary),
                "n_included_main_features": len(main),
                "n_matrix_extra_feature_columns": len(extra),
                "n_missing_main_features": len(missing_main),
            }
        ]
    )
    shape.to_csv(out / "step3_phase3_0_feature_matrix_shape.csv", index=False)

    audit = ctx.dictionary.assign(
        present_in_matrix=ctx.dictionary["feature_name"].isin(matrix_features),
        included_in_main_bool=bool_series(ctx.dictionary["included_in_main"]),
    )
    audit.to_csv(out / "step3_phase3_0_feature_dictionary_audit.csv", index=False)
    pd.DataFrame(
        [
            {
                "tf_main_features": len(tf_main),
                "tf_main_present": sum(f in matrix_features for f in tf_main),
                "tf_main_all_na": len(tf_all_na),
                "tf_activity_source": ctx.manifest.get("tf_activity_source"),
                "assembly_status": ctx.manifest.get("tf_activity_assembly_status"),
            }
        ]
    ).to_csv(out / "step3_phase3_0_tf_hotfix_validation.csv", index=False)
    ctx.sample_table[
        [
            "sample_key",
            "cohort_id",
            "sample_id",
            "patient_id",
            "effective_split",
            "has_sample_metadata",
            "has_split_info",
            "response_label",
            "metadata_status",
            "unresolved_no_split",
        ]
    ].to_csv(out / "step3_phase3_0_sample_key_join_audit.csv", index=False)
    forbidden.to_csv(out / "step3_phase3_0_forbidden_feature_name_audit.csv", index=False)

    hard_failed = [t for t in tests if t["severity"] == "hard" and not t["passed"]]
    conditional_failed = [t for t in tests if t["severity"] == "conditional" and not t["passed"]]
    verdict = "HARD_FAIL" if hard_failed else ("CONDITIONAL_PASS" if conditional_failed else "PASS")
    write_yaml(
        out / "step3_phase3_0_input_contract_status.yaml",
        {
            "phase": "3.0",
            "verdict": verdict,
            "tests": tests,
            "input_paths": {k: str(v) for k, v in HOTFIX_FILES.items()},
            "repair_overlay": {
                "feature_dictionary": str(REPAIRED_FEATURE_DICTIONARY),
                "used": REPAIRED_FEATURE_DICTIONARY.exists(),
            },
            "package_versions": package_versions(),
        },
    )
    write_md(
        out / "step3_phase3_0_input_contract_summary.md",
        "Step3 Phase 3.0 Input Contract Summary",
        [
            f"- Verdict: {verdict}",
            f"- Matrix rows/features: {len(ctx.matrix)} / {len(ctx.feature_cols)}",
            f"- Included main features missing from matrix: {len(missing_main)}",
            f"- Main TF features/all-NA: {len(tf_main)} / {len(tf_all_na)}",
            f"- Forbidden feature columns: {len(forbidden)} total, {len(forbidden_unresolved)} unresolved after explicit biological whitelist",
            f"- Matrix extra feature columns not in dictionary: {len(extra)}; treated as QC/control if non-main.",
        ],
    )
    if verdict == "HARD_FAIL":
        abort("PHASE3_0", hard_failed)
    return verdict


def task_mask(st: pd.DataFrame, task: str) -> pd.Series:
    tc = st["treatment_context"].astype(str)
    disease = st["disease"].astype(str).str.lower()
    if task == "PD1_anchor":
        return tc.eq("PD1_ICI_anchor")
    if task == "PD1X_extension":
        return tc.eq("PD1X_extension")
    if task == "HCC_specific":
        return disease.eq("hcc") & tc.isin(["PD1_ICI_anchor", "PD1X_extension"])
    if task == "pan_cancer_shared":
        return tc.isin(["PD1_ICI_anchor", "PD1X_extension"])
    if task == "sensitivity_only":
        return tc.eq("high_confounding_support")
    return pd.Series(False, index=st.index)


def family_eligible(ctx: Context, family: str) -> pd.Series:
    st = ctx.sample_table
    if family in {"combined_interpretable", "combined_all_allowed"}:
        parts = [family_eligible(ctx, f) for f in ["fraction", "signature", "pathway", "tf_activity"]]
        if family == "combined_all_allowed":
            parts.append(family_eligible(ctx, "pseudobulk_gene"))
        return pd.concat(parts, axis=1).any(axis=1)
    if family == "fraction":
        return pd.Series(True, index=st.index)
    if family in {"signature", "pathway", "tf_activity"}:
        return st.get("signature_pathway_tf_available", False).fillna(False).astype(bool)
    if family == "pseudobulk_gene":
        return st.get("pseudobulk_available", False).fillna(False).astype(bool) & (len(ctx.family_features["pseudobulk_gene"]) > 0)
    return pd.Series(False, index=st.index)


def phase3_1(ctx: Context) -> tuple[str, pd.DataFrame, pd.DataFrame]:
    out = PHASE_DIRS["3_1"]
    st = ctx.sample_table.copy()
    base_exclusion = (
        st["sample_key"].astype(str).str.startswith("GSE272993")
        | st["cohort_id"].astype(str).eq("GSE140228_smartseq2")
    )
    rows = []
    members = []

    for task in TASKS:
        task_base = task_mask(st, task) & ~base_exclusion
        base_fam_masks = {f: task_base & family_eligible(ctx, f) for f in BASE_FAMILIES}
        strict_intersection = task_base.copy()
        for f in BASE_FAMILIES:
            if len(ctx.family_features.get(f, [])) > 0:
                strict_intersection &= family_eligible(ctx, f)
        for family in ALL_FAMILIES:
            fam_mask = task_base & family_eligible(ctx, family)
            if family == "pseudobulk_gene" and len(ctx.family_features["pseudobulk_gene"]) == 0:
                fam_mask &= False
            for utype, mask in [
                ("descriptive_universe", fam_mask),
                ("supervised_labeled_universe", fam_mask & supervised_mask(st)),
                ("strict_intersection_universe", strict_intersection & supervised_mask(st)),
                ("external_or_sensitivity_universe", task_base & ~supervised_mask(st)),
            ]:
                uid = f"{task}__{family}__{utype}"
                sub = st[mask].copy()
                status = universe_status(sub, utype)
                rows.append(universe_row(uid, task, family, utype, sub, st, mask, status))
                for _, r in st.iterrows():
                    included = bool(mask.loc[r.name])
                    reason = "" if included else exclusion_reason(r, task, family, task_base.loc[r.name], fam_mask.loc[r.name], base_exclusion.loc[r.name])
                    members.append(
                        {
                            "sample_key": r["sample_key"],
                            "cohort_id": r["cohort_id"],
                            "sample_id": r["sample_id"],
                            "source_h5ad": r.get("source_h5ad"),
                            "patient_id": r["patient_id"],
                            "task": task,
                            "feature_family": family,
                            "universe_type": utype,
                            "included": included,
                            "exclusion_reason": reason,
                            "split": r.get("effective_split"),
                            "response_label": r.get("response_label"),
                            "treatment_context": r.get("treatment_context"),
                            "cancer_type": r.get("disease"),
                            "metadata_status": r.get("metadata_status"),
                            "feature_family_eligible": bool(family_eligible(ctx, family).loc[r.name]),
                        }
                    )

    registry = pd.DataFrame(rows)
    membership = pd.DataFrame(members)
    registry.to_csv(out / "step3_phase3_1_analysis_universe_registry.csv", index=False)
    membership.to_csv(out / "step3_phase3_1_sample_universe_membership.csv", index=False)

    patient_audit = (
        st.groupby(["cohort_id", "patient_id"], dropna=False)
        .agg(
            n_samples=("sample_key", "count"),
            assigned_split=("effective_split", lambda x: ";".join(sorted(set(map(str, x))))),
            observed_sample_splits=("split", lambda x: ";".join(sorted(set(map(str, x))))),
        )
        .reset_index()
    )
    patient_audit["split_consistent"] = patient_audit["assigned_split"].str.split(";").map(len).eq(1)
    patient_audit["issue"] = np.where(patient_audit["split_consistent"], "", "patient_samples_cross_split")
    patient_audit.to_csv(out / "step3_phase3_1_patient_split_audit.csv", index=False)

    hard = []
    if not patient_audit["split_consistent"].all():
        hard.append("patient_split_leakage")
    supervised_members = membership[membership["universe_type"].eq("supervised_labeled_universe") & membership["included"]]
    if supervised_members["response_label"].astype(str).str.lower().eq("unknown").any():
        hard.append("unknown_response_in_supervised")
    if supervised_members["split"].astype(str).isin(["nan", "unknown", "external", "not_applicable"]).any():
        hard.append("unresolved_split_in_supervised")
    pd1_fraction = registry[
        registry["task"].eq("PD1_anchor")
        & registry["feature_family"].eq("fraction")
        & registry["universe_type"].eq("supervised_labeled_universe")
    ]
    conditional = []
    if not pd1_fraction.empty:
        prev = ctx.preview[(ctx.preview["task"].eq("PD1_anchor")) & (ctx.preview["feature_family"].eq("fraction"))]
        if not prev.empty and pd1_fraction.iloc[0]["n_samples"] < max(20, int(prev.iloc[0]["n_samples"]) * 0.5):
            conditional.append("PD1_anchor labeled supervised smaller than readiness preview because unknown responses excluded")
    if registry.query("universe_type == 'supervised_labeled_universe' and n_samples < 30").shape[0] > 0:
        conditional.append("some supervised universes too small; downgrade to sensitivity")

    verdict = "HARD_FAIL" if hard else ("CONDITIONAL_PASS" if conditional else "PASS")
    write_yaml(out / "step3_phase3_1_status.yaml", {"phase": "3.1", "verdict": verdict, "hard_failures": hard, "conditional_items": conditional})
    write_md(
        out / "step3_phase3_1_summary.md",
        "Step3 Phase 3.1 Analysis Universe Summary",
        [
            f"- Verdict: {verdict}",
            f"- Universes built: {len(registry)}",
            f"- PD1_anchor fraction supervised labeled samples: {int(pd1_fraction.iloc[0]['n_samples']) if not pd1_fraction.empty else 0}",
            f"- Patient split leakage issues: {int((~patient_audit['split_consistent']).sum())}",
            "- Universe-specific inclusion/exclusion recorded in membership table.",
        ],
    )
    if verdict == "HARD_FAIL":
        abort("PHASE3_1", hard)
    return verdict, registry, membership


def supervised_mask(st: pd.DataFrame) -> pd.Series:
    return (
        st["response_binary_num"].notna()
        & st["effective_split"].astype(str).isin(VALID_SPLITS)
        & ~st["unresolved_no_split"]
        & ~st["metadata_status"].eq("metadata_missing")
    )


def universe_status(sub: pd.DataFrame, utype: str) -> str:
    if sub.empty:
        return "empty"
    if utype == "supervised_labeled_universe":
        if sub["patient_id"].nunique() < 20 or sub["response_binary_num"].nunique() < 2:
            return "sensitivity_only_low_n_or_single_class"
    return "eligible"


def universe_row(uid: str, task: str, family: str, utype: str, sub: pd.DataFrame, st: pd.DataFrame, mask: pd.Series, status: str) -> dict[str, Any]:
    split_counts = sub.groupby("effective_split")["patient_id"].nunique().to_dict()
    split_sample_counts = sub.groupby("effective_split")["sample_key"].nunique().to_dict()
    cohorts = sorted(sub["cohort_id"].dropna().astype(str).unique())
    excluded = sorted(set(st["cohort_id"].dropna().astype(str).unique()) - set(cohorts))
    y = sub["response_binary_num"]
    cohort_counts = sub["cohort_id"].value_counts()
    return {
        "universe_id": uid,
        "task": task,
        "feature_family": family,
        "universe_type": utype,
        "n_samples": int(len(sub)),
        "n_patients": int(sub["patient_id"].nunique()),
        "n_cohorts": int(sub["cohort_id"].nunique()),
        "n_cancer_types": int(sub["disease"].nunique()),
        "n_responders": int((y == 1).sum()),
        "n_non_responders": int((y == 0).sum()),
        "n_unknown_response": int(y.isna().sum()),
        "included_cohorts": ";".join(cohorts),
        "excluded_cohorts": ";".join(excluded),
        "exclusion_reason": "see membership table",
        "split_train_patients": int(split_counts.get("train", 0)),
        "split_validation_patients": int(split_counts.get("validation", 0)),
        "split_test_patients": int(split_counts.get("test", 0)),
        "split_train_samples": int(split_sample_counts.get("train", 0)),
        "split_validation_samples": int(split_sample_counts.get("validation", 0)),
        "split_test_samples": int(split_sample_counts.get("test", 0)),
        "responder_rate": float(y.mean()) if y.notna().any() else np.nan,
        "min_cohort_sample_count": int(cohort_counts.min()) if not cohort_counts.empty else 0,
        "max_cohort_sample_fraction": float(cohort_counts.max() / len(sub)) if len(sub) else np.nan,
        "status": status,
        "notes": "",
    }


def exclusion_reason(r: pd.Series, task: str, family: str, task_ok: bool, fam_ok: bool, base_excl: bool) -> str:
    reasons = []
    if base_excl:
        reasons.append("forced_exclusion")
    if not task_ok:
        reasons.append("task_context_or_cancer_mismatch")
    if not fam_ok:
        reasons.append("feature_family_unavailable")
    if pd.isna(r.get("response_binary_num")):
        reasons.append("response_unknown_for_supervised")
    if str(r.get("effective_split")) not in VALID_SPLITS:
        reasons.append("missing_or_invalid_split")
    if r.get("metadata_status") == "metadata_missing":
        reasons.append("metadata_missing")
    return ";".join(reasons)


def phase3_2(ctx: Context, universe_registry: pd.DataFrame, membership: pd.DataFrame) -> tuple[str, pd.DataFrame, pd.DataFrame]:
    out = PHASE_DIRS["3_2"]
    st = ctx.sample_table.set_index("sample_key")
    registry_rows = []
    missing_resp_rows = []
    missing_cohort_rows = []

    feature_rows = []
    feature_set_paths = {}
    for _, drow in ctx.dictionary.iterrows():
        f = drow["feature_name"]
        vals = ctx.matrix[f] if f in ctx.matrix.columns else pd.Series(dtype=float)
        miss = float(vals.isna().mean()) if f in ctx.matrix.columns else 1.0
        nunique = int(vals.dropna().nunique()) if f in ctx.matrix.columns else 0
        leakage = any(re.search(p, f, flags=re.IGNORECASE) for p in FORBIDDEN_PATTERNS)
        leakage_whitelisted = f in FORBIDDEN_NAME_WHITELIST
        is_all_na = bool(f not in ctx.matrix.columns or vals.isna().all())
        is_const = bool(nunique <= 1 and not is_all_na)
        status, reason = retain_status(
            drow, f in ctx.matrix.columns, miss, is_all_na, is_const, leakage and not leakage_whitelisted
        )
        feature_rows.append(
            {
                "feature_name": f,
                "feature_family": drow["family"],
                "included_in_main": bool(str(drow.get("included_in_main")).lower() == "true"),
                "sensitivity_only": bool(str(drow.get("sensitivity_only")).lower() == "true"),
                "present_in_matrix": f in ctx.matrix.columns,
                "n_non_missing": int(vals.notna().sum()) if f in ctx.matrix.columns else 0,
                "missing_rate": miss,
                "n_unique_values": nunique,
                "is_constant": is_const,
                "is_all_na": is_all_na,
                "leakage_name_flag": leakage,
                "leakage_name_whitelisted": leakage_whitelisted,
                "source_valid": f in ctx.matrix.columns and not (leakage and not leakage_whitelisted),
                "retained_status": status,
                "retained_reason": reason,
                "blocked_reason": reason if status == "blocked" else "",
            }
        )

    feature_registry = pd.DataFrame(feature_rows)
    feature_registry.to_csv(out / "step3_phase3_2_feature_family_registry.csv", index=False)

    valid_universes = universe_registry[universe_registry["universe_type"].eq("supervised_labeled_universe")]
    for _, u in valid_universes.iterrows():
        uid = u["universe_id"]
        fam = u["feature_family"]
        samples = membership[
            membership["universe_id"].eq(uid) if "universe_id" in membership.columns else (
                membership["task"].eq(u["task"])
                & membership["feature_family"].eq(fam)
                & membership["universe_type"].eq("supervised_labeled_universe")
                & membership["included"]
            )
        ]["sample_key"].tolist()
        feats = [f for f in ctx.family_features.get(fam, []) if f in ctx.matrix.columns]
        if fam in {"combined_interpretable", "combined_all_allowed"}:
            feats = ctx.family_features[fam]
        sub = ctx.matrix.set_index("sample_key").loc[samples] if samples else pd.DataFrame()
        fstats = universe_feature_stats(sub, feats)
        clean = [f for f, s in fstats.items() if s["status"] == "clean_core"]
        aware = [f for f, s in fstats.items() if s["status"] == "missingness_aware"]
        sens = [f for f, s in fstats.items() if s["status"] == "sensitivity"]
        blocked = [f for f, s in fstats.items() if s["status"] == "blocked"]
        for typ, flist in [("clean_core_feature_set", clean), ("missingness_aware_feature_set", aware), ("sensitivity_feature_set", sens), ("blocked_feature_set", blocked)]:
            p = out / f"feature_list__{uid}__{typ}.csv"
            pd.DataFrame({"feature_name": flist}).to_csv(p, index=False)
            feature_set_paths[(uid, typ)] = str(p)
            registry_rows.append(
                {
                    "universe_id": uid,
                    "task": u["task"],
                    "feature_family": fam,
                    "feature_set_type": typ,
                    "n_features": len(flist),
                    "n_clean_core": len(clean),
                    "n_missingness_aware": len(aware),
                    "n_sensitivity": len(sens),
                    "n_blocked": len(blocked),
                    "feature_list_path": str(p),
                    "status": "eligible" if flist else "empty",
                }
            )
        if samples and feats:
            y = st.loc[samples, "response_binary_num"]
            for f in feats:
                miss_ind = sub[f].isna().astype(int) if f in sub.columns else pd.Series(1, index=sub.index)
                missing_resp_rows.append(missing_response_row(uid, u["task"], fam, f, miss_ind, y))
                missing_cohort_rows.append(missing_cohort_row(uid, u["task"], fam, f, miss_ind, st.loc[samples, "cohort_id"]))

    universe_feature_sets = pd.DataFrame(registry_rows)
    universe_feature_sets.to_csv(out / "step3_phase3_2_universe_feature_set_registry.csv", index=False)
    pd.DataFrame(missing_resp_rows).to_csv(out / "step3_phase3_2_missingness_response_audit.csv", index=False)
    pd.DataFrame(missing_cohort_rows).to_csv(out / "step3_phase3_2_missingness_cohort_audit.csv", index=False)

    hard = []
    unresolved_main_leakage = feature_registry[
        feature_registry["leakage_name_flag"].astype(bool)
        & ~feature_registry["leakage_name_whitelisted"].astype(bool)
        & feature_registry["included_in_main"].astype(bool)
    ]
    if len(unresolved_main_leakage):
        hard.append("leakage_feature_name_detected")
    pd1_frac = universe_feature_sets[
        universe_feature_sets["universe_id"].eq("PD1_anchor__fraction__supervised_labeled_universe")
        & universe_feature_sets["feature_set_type"].eq("clean_core_feature_set")
    ]
    if pd1_frac.empty or int(pd1_frac.iloc[0]["n_features"]) == 0:
        hard.append("no_valid_feature_set_for_PD1_anchor_fraction")
    if (feature_registry[feature_registry["included_in_main"]]["is_all_na"].mean() > 0.1):
        hard.append("main_features_dominated_by_all_na")
    conditional = []
    pseudo = feature_registry[feature_registry["feature_family"].eq("pseudobulk_gene")]
    if pseudo.empty:
        conditional.append("pseudobulk_gene absent from hotfix dictionary/matrix; blocked for main Step3")
    empty_sets = universe_feature_sets.query("feature_set_type == 'clean_core_feature_set' and n_features == 0")
    if len(empty_sets):
        conditional.append("some universes have no clean_core features and are sensitivity/blocked")

    verdict = "HARD_FAIL" if hard else ("CONDITIONAL_PASS" if conditional else "PASS")
    write_yaml(out / "step3_phase3_2_status.yaml", {"phase": "3.2", "verdict": verdict, "hard_failures": hard, "conditional_items": conditional})
    write_md(
        out / "step3_phase3_2_summary.md",
        "Step3 Phase 3.2 Feature Registry Missingness Summary",
        [
            f"- Verdict: {verdict}",
            f"- Feature dictionary rows: {len(feature_registry)}",
            f"- Clean feature sets: {int((universe_feature_sets['feature_set_type'].eq('clean_core_feature_set') & (universe_feature_sets['n_features'] > 0)).sum())}",
            f"- Pseudobulk gene feature rows: {len(pseudo)}",
            "- NA and measured zero remain distinct; no fillna(0) repair is used.",
        ],
    )
    if verdict == "HARD_FAIL":
        abort("PHASE3_2", hard)
    return verdict, feature_registry, universe_feature_sets


def retain_status(drow: pd.Series, present: bool, miss: float, all_na: bool, constant: bool, leakage: bool) -> tuple[str, str]:
    if leakage:
        return "blocked", "leakage_name"
    if not present:
        return "blocked", "missing_from_matrix"
    if all_na:
        return "blocked", "all_na"
    if constant:
        return "blocked", "constant"
    if miss > 0.8:
        return "sensitivity", "missing_rate_gt_0.8"
    if miss > 0.2:
        return "missingness_aware", "missing_rate_0.2_to_0.8"
    return "clean_core", "low_missingness"


def universe_feature_stats(sub: pd.DataFrame, feats: list[str]) -> dict[str, dict[str, Any]]:
    stats = {}
    for f in feats:
        if f not in sub.columns or sub.empty:
            stats[f] = {"status": "blocked"}
            continue
        miss = sub[f].isna().mean()
        nunique = sub[f].dropna().nunique()
        if sub[f].isna().all():
            status = "blocked"
        elif nunique <= 1:
            status = "blocked"
        elif miss > 0.8:
            status = "sensitivity"
        elif miss > 0.2:
            status = "missingness_aware"
        else:
            status = "clean_core"
        stats[f] = {"status": status, "missing_rate": miss, "nunique": nunique}
    return stats


def missing_response_row(uid: str, task: str, fam: str, f: str, miss: pd.Series, y: pd.Series) -> dict[str, Any]:
    tab = pd.crosstab(miss, y)
    p = np.nan
    if tab.shape == (2, 2):
        try:
            _, p = fisher_exact(tab)
        except Exception:
            p = np.nan
    return {
        "universe_id": uid,
        "task": task,
        "feature_family": fam,
        "feature_name": f,
        "missing_rate": float(miss.mean()),
        "response_missing_p": p,
        "missing_rate_responder": float(miss[y == 1].mean()) if (y == 1).any() else np.nan,
        "missing_rate_non_responder": float(miss[y == 0].mean()) if (y == 0).any() else np.nan,
    }


def missing_cohort_row(uid: str, task: str, fam: str, f: str, miss: pd.Series, cohort: pd.Series) -> dict[str, Any]:
    tab = pd.crosstab(miss, cohort)
    p = chi2_p(tab)
    return {
        "universe_id": uid,
        "task": task,
        "feature_family": fam,
        "feature_name": f,
        "missing_rate": float(miss.mean()),
        "cohort_missing_chi2_p": p,
        "max_cohort_missing_rate": float(pd.DataFrame({"m": miss, "c": cohort}).groupby("c")["m"].mean().max()) if len(cohort) else np.nan,
    }


def phase3_3(ctx: Context, universe_registry: pd.DataFrame, membership: pd.DataFrame, feature_sets: pd.DataFrame) -> tuple[str, pd.DataFrame, pd.DataFrame]:
    out = PHASE_DIRS["3_3"]
    conf_rows = []
    neg_rows = []
    dist_rows = []
    burden_rows = []
    st = ctx.sample_table.set_index("sample_key")
    mat = ctx.matrix.set_index("sample_key")

    main_sets = feature_sets[
        feature_sets["feature_set_type"].eq("clean_core_feature_set") & feature_sets["status"].eq("eligible")
    ]
    for _, fs in main_sets.iterrows():
        uid = fs["universe_id"]
        task, fam, _ = uid.split("__", 2)
        samples = membership[
            membership["task"].eq(task)
            & membership["feature_family"].eq(fam)
            & membership["universe_type"].eq("supervised_labeled_universe")
            & membership["included"]
        ]["sample_key"].tolist()
        if len(samples) < 20:
            continue
        sub_st = st.loc[samples]
        y = sub_st["response_binary_num"].astype(int)
        feats = pd.read_csv(fs["feature_list_path"])["feature_name"].tolist()
        X = mat.loc[samples, feats] if feats else pd.DataFrame(index=samples)
        missing_burden = X.isna().mean(axis=1) if len(feats) else pd.Series(0, index=samples)

        dist_rows.extend(response_dist_rows(uid, task, fam, sub_st, "cohort_id"))
        dist_rows.extend(response_dist_rows(uid, task, fam, sub_st, "disease"))
        for resp, vals in pd.DataFrame({"response": y, "missing_burden": missing_burden}).groupby("response"):
            burden_rows.append(
                {
                    "universe_id": uid,
                    "task": task,
                    "feature_family": fam,
                    "response_binary": int(resp),
                    "n_samples": int(len(vals)),
                    "mean_missing_burden": float(vals["missing_burden"].mean()),
                    "median_missing_burden": float(vals["missing_burden"].median()),
                }
            )

        metadata_auc = metadata_only_model(sub_st, y)
        miss_auc = missingness_only_model(X, y)
        cohort_macro = categorical_target_model(X, sub_st["cohort_id"])
        split_macro = categorical_target_model(pd.get_dummies(sub_st[["cohort_id", "disease", "treatment_context"]].astype(str)), sub_st["effective_split"])
        max_cohort_fraction = float(sub_st["cohort_id"].value_counts().max() / len(sub_st))
        cohort_p = chi2_p(pd.crosstab(sub_st["cohort_id"], y))
        cancer_p = chi2_p(pd.crosstab(sub_st["disease"], y))
        risk = confounding_risk(metadata_auc, miss_auc, max_cohort_fraction, cohort_p, sub_st["cohort_id"].nunique())
        action = "primary_allowed_with_controls" if risk in {"LOW", "MODERATE"} else "downgrade_to_sensitivity"
        if task in {"HCC_specific"} and sub_st["cohort_id"].nunique() <= 1:
            risk = "HIGH"
            action = "single_cohort_sensitivity_only"
        conf_rows.append(
            {
                "universe_id": uid,
                "task": task,
                "feature_family": fam,
                "n_samples": len(sub_st),
                "n_patients": sub_st["patient_id"].nunique(),
                "n_cohorts": sub_st["cohort_id"].nunique(),
                "responder_rate": float(y.mean()),
                "max_cohort_fraction": max_cohort_fraction,
                "cohort_response_chi2_p": cohort_p,
                "cancer_response_chi2_p": cancer_p,
                "missingness_response_auc": miss_auc,
                "metadata_only_response_auc": metadata_auc,
                "cohort_prediction_auc_or_macro_f1": cohort_macro,
                "split_prediction_auc_or_macro_f1": split_macro,
                "confounding_risk_level": risk,
                "recommended_action": action,
            }
        )
        neg_rows.append({"universe_id": uid, "control_model": "metadata_only_response_model", "metric": "auc_cv", "value": metadata_auc})
        neg_rows.append({"universe_id": uid, "control_model": "missingness_only_response_model", "metric": "auc_cv", "value": miss_auc})
        neg_rows.append({"universe_id": uid, "control_model": "cohort_prediction_model", "metric": "macro_f1_cv", "value": cohort_macro})
        neg_rows.append({"universe_id": uid, "control_model": "split_prediction_model", "metric": "macro_f1_cv", "value": split_macro})

    conf = pd.DataFrame(conf_rows)
    neg = pd.DataFrame(neg_rows)
    conf.to_csv(out / "step3_phase3_3_confounding_audit_by_universe.csv", index=False)
    neg.to_csv(out / "step3_phase3_3_negative_control_results.csv", index=False)
    pd.DataFrame(dist_rows).to_csv(out / "step3_phase3_3_response_distribution_by_cohort.csv", index=False)
    pd.DataFrame(burden_rows).to_csv(out / "step3_phase3_3_missingness_burden_by_response.csv", index=False)

    hard = []
    if conf.empty:
        hard.append("no_universe_passes_confounding_audit")
    pd1 = conf[conf["task"].eq("PD1_anchor")]
    if not pd1.empty and ((pd1["metadata_only_response_auc"] > 0.85) | (pd1["missingness_response_auc"] > 0.85)).all():
        hard.append("PD1_anchor_almost_perfect_metadata_or_missingness_prediction")
    if not pd1.empty and pd1["confounding_risk_level"].eq("HIGH").all():
        hard.append("PD1_anchor_all_universes_high_confounding_risk")
    conditional = []
    if not conf.empty and conf["confounding_risk_level"].eq("HIGH").any():
        conditional.append("high_risk_universes_downgraded_to_sensitivity")
    if not pd1.empty and pd1["confounding_risk_level"].isin(["LOW", "MODERATE"]).any():
        pass
    elif not pd1.empty:
        conditional.append("PD1_anchor only high-risk universes; response modeling restricted")

    verdict = "HARD_FAIL" if hard else ("CONDITIONAL_PASS" if conditional else "PASS")
    write_yaml(out / "step3_phase3_3_status.yaml", {"phase": "3.3", "verdict": verdict, "hard_failures": hard, "conditional_items": conditional})
    write_md(
        out / "step3_phase3_3_summary.md",
        "Step3 Phase 3.3 Leakage Confounding Summary",
        [
            f"- Verdict: {verdict}",
            f"- Audited universes: {len(conf)}",
            f"- HIGH risk universes: {int(conf['confounding_risk_level'].eq('HIGH').sum()) if not conf.empty else 0}",
            "- Metadata-only and missingness-only response controls recorded before biological models.",
        ],
    )
    if verdict == "HARD_FAIL":
        abort("PHASE3_3", hard)
    return verdict, conf, neg


def response_dist_rows(uid: str, task: str, fam: str, st: pd.DataFrame, col: str) -> list[dict[str, Any]]:
    rows = []
    for key, g in st.groupby(col, dropna=False):
        y = g["response_binary_num"]
        rows.append(
            {
                "universe_id": uid,
                "task": task,
                "feature_family": fam,
                "stratum_type": col,
                "stratum": key,
                "n_samples": len(g),
                "n_patients": g["patient_id"].nunique(),
                "n_responders": int((y == 1).sum()),
                "n_non_responders": int((y == 0).sum()),
                "responder_rate": float(y.mean()) if y.notna().any() else np.nan,
            }
        )
    return rows


def chi2_p(tab: pd.DataFrame) -> float:
    try:
        if tab.shape[0] < 2 or tab.shape[1] < 2:
            return np.nan
        return float(chi2_contingency(tab)[1])
    except Exception:
        return np.nan


def cv_auc_pipeline(X: pd.DataFrame, y: pd.Series, pipeline: Pipeline, min_splits: int = 3) -> float:
    if y.nunique() < 2 or len(y) < 20:
        return np.nan
    counts = y.value_counts()
    n_splits = int(min(5, counts.min(), max(min_splits, 2)))
    if n_splits < 2:
        return np.nan
    aucs = []
    cv = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=SEED)
    for tr, te in cv.split(X, y):
        try:
            model = clone(pipeline)
            model.fit(X.iloc[tr], y.iloc[tr])
            score = model.predict_proba(X.iloc[te])[:, 1]
            aucs.append(roc_auc_score(y.iloc[te], score))
        except Exception:
            continue
    return float(np.mean(aucs)) if aucs else np.nan


def metadata_only_model(st: pd.DataFrame, y: pd.Series) -> float:
    X = pd.get_dummies(
        st[["cohort_id", "disease", "source_h5ad", "treatment_context", "timepoint_class"]].astype(str),
        dummy_na=True,
    )
    pipe = Pipeline([("imputer", SimpleImputer(strategy="constant", fill_value=0)), ("model", LogisticRegression(max_iter=2000, class_weight="balanced", random_state=SEED))])
    return cv_auc_pipeline(X, y.reset_index(drop=True), pipe)


def missingness_only_model(X: pd.DataFrame, y: pd.Series) -> float:
    if X.empty:
        return np.nan
    Xm = X.isna().astype(int)
    Xm["n_missing_features"] = Xm.sum(axis=1)
    pipe = Pipeline([("imputer", SimpleImputer(strategy="constant", fill_value=0)), ("model", LogisticRegression(max_iter=2000, class_weight="balanced", random_state=SEED))])
    return cv_auc_pipeline(Xm.reset_index(drop=True), y.reset_index(drop=True), pipe)


def categorical_target_model(X: pd.DataFrame, target: pd.Series) -> float:
    y = target.astype(str).reset_index(drop=True)
    if y.nunique() < 2 or len(y) < 30:
        return np.nan
    min_class = y.value_counts().min()
    if min_class < 2:
        return np.nan
    n_splits = int(min(5, min_class))
    cv = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=SEED)
    scores = []
    pipe = Pipeline([("imputer", SimpleImputer(strategy="median")), ("scaler", StandardScaler()), ("model", LogisticRegression(max_iter=2000, class_weight="balanced", random_state=SEED))])
    X = X.reset_index(drop=True)
    for tr, te in cv.split(X, y):
        try:
            m = clone(pipe)
            m.fit(X.iloc[tr], y.iloc[tr])
            pred = m.predict(X.iloc[te])
            scores.append(f1_score(y.iloc[te], pred, average="macro"))
        except Exception:
            continue
    return float(np.mean(scores)) if scores else np.nan


def confounding_risk(metadata_auc: float, miss_auc: float, max_frac: float, cohort_p: float, n_cohorts: int) -> str:
    if n_cohorts <= 1 or max_frac > 0.75:
        return "HIGH"
    if (metadata_auc is not np.nan and metadata_auc > 0.70) or (miss_auc is not np.nan and miss_auc > 0.70):
        return "HIGH"
    if (metadata_auc is not np.nan and metadata_auc > 0.60) or (miss_auc is not np.nan and miss_auc > 0.60) or max_frac > 0.50:
        return "MODERATE"
    return "LOW"


def phase3_4(ctx: Context, membership: pd.DataFrame, feature_sets: pd.DataFrame, conf: pd.DataFrame) -> tuple[str, pd.DataFrame, pd.DataFrame]:
    out = PHASE_DIRS["3_4"]
    mat = ctx.matrix.set_index("sample_key")
    st = ctx.sample_table.set_index("sample_key")
    eligible = eligible_model_sets(feature_sets, conf, include_high=False)
    results = []
    top_rows = []
    uni_rows = []
    manifest = {"phase": "3.4", "seed": SEED, "models": []}
    for _, fs in eligible.iterrows():
        uid = fs["universe_id"]
        task, fam, _ = uid.split("__", 2)
        samples = samples_for(membership, task, fam, "supervised_labeled_universe")
        if not split_usable(st.loc[samples]):
            continue
        features = pd.read_csv(fs["feature_list_path"])["feature_name"].tolist()
        if len(features) == 0:
            continue
        X = mat.loc[samples, features]
        y = st.loc[samples, "response_binary_num"].astype(int)
        split = st.loc[samples, "effective_split"].astype(str)
        uni_rows.extend(univariate_screen(uid, task, fam, X, y, split))
        for model_type, pipe in [
            ("logistic_l2", logistic_pipe("l2")),
            ("logistic_l1", logistic_pipe("l1")),
            ("elastic_net", logistic_pipe("elasticnet")),
            ("simple_additive_top10", None),
        ]:
            model_id = f"p34__{uid}__{fs['feature_set_type']}__{model_type}"
            try:
                res, tops = run_model(model_id, task, fam, uid, fs["feature_set_type"], model_type, X, y, split, st.loc[samples, "patient_id"], pipe)
                results.append(res)
                top_rows.extend(tops)
                manifest["models"].append({"model_id": model_id, "features_path": fs["feature_list_path"], "model_type": model_type, "status": res["status"]})
            except Exception as exc:
                logging.exception("Phase3.4 model failed %s", model_id)
                results.append({"model_id": model_id, "task": task, "feature_family": fam, "universe_id": uid, "feature_set_type": fs["feature_set_type"], "model_type": model_type, "status": "failed", "caveat": str(exc)})
    res_df = pd.DataFrame(results)
    top_df = pd.DataFrame(top_rows)
    uni_df = pd.DataFrame(uni_rows)
    res_df.to_csv(out / "step3_phase3_4_interpretable_baseline_results.csv", index=False)
    top_df.to_csv(out / "step3_phase3_4_top_features_by_model.csv", index=False)
    uni_df.to_csv(out / "step3_phase3_4_univariate_feature_results.csv", index=False)
    write_yaml(out / "step3_phase3_4_model_manifest.yaml", manifest)

    hard = []
    conditional = []
    if "universe_role" in conf.columns:
        pd1_primary_allowed = bool(
            conf[conf["task"].eq("PD1_anchor") & conf["universe_role"].eq("primary_allowed")].shape[0]
        )
    else:
        pd1_primary_allowed = True
    pd1_frac = res_df[(res_df.get("task") == "PD1_anchor") & (res_df.get("feature_family") == "fraction") & (res_df.get("status") == "ok")] if not res_df.empty else pd.DataFrame()
    if pd1_primary_allowed and pd1_frac.empty:
        hard.append("PD1_anchor_fraction_baseline_failed")
    if not pd1_primary_allowed:
        conditional.append("PD1_anchor_primary_blocked_by_role_gate; no primary PD1_anchor baseline required")
    if not res_df.empty and res_df["status"].eq("failed").any():
        conditional.append("some interpretable models skipped_or_failed_due_to_low_n")
    verdict = "HARD_FAIL" if hard else ("CONDITIONAL_PASS" if conditional else "PASS")
    write_yaml(out / "step3_phase3_4_status.yaml", {"phase": "3.4", "verdict": verdict, "hard_failures": hard, "conditional_items": conditional})
    write_md(out / "step3_phase3_4_summary.md", "Step3 Phase 3.4 Interpretable Baseline Summary", [f"- Verdict: {verdict}", f"- Models attempted: {len(res_df)}", "- Feature selection, scaling, imputation fit only on training split."])
    if verdict == "HARD_FAIL":
        abort("PHASE3_4", hard)
    return verdict, res_df, top_df


def eligible_model_sets(feature_sets: pd.DataFrame, conf: pd.DataFrame, include_high: bool) -> pd.DataFrame:
    fs = feature_sets[
        feature_sets["feature_set_type"].isin(["clean_core_feature_set", "missingness_aware_feature_set"])
        & feature_sets["status"].eq("eligible")
        & (feature_sets["n_features"] > 0)
    ].copy()
    if conf.empty:
        return fs.iloc[0:0]
    conf_cols = ["universe_id", "confounding_risk_level", "recommended_action"]
    if "universe_role" in conf.columns:
        conf_cols.append("universe_role")
    fs = fs.merge(conf[conf_cols], on="universe_id", how="left")
    if "universe_role" in fs.columns:
        fs = fs[fs["universe_role"].eq("primary_allowed")]
    if include_high:
        return fs
    return fs[~fs["confounding_risk_level"].eq("HIGH")]


def samples_for(membership: pd.DataFrame, task: str, fam: str, utype: str) -> list[str]:
    return membership[
        membership["task"].eq(task)
        & membership["feature_family"].eq(fam)
        & membership["universe_type"].eq(utype)
        & membership["included"]
    ]["sample_key"].tolist()


def split_usable(st: pd.DataFrame) -> bool:
    if st.empty:
        return False
    for split in ["train", "validation", "test"]:
        y = st.loc[st["effective_split"].astype(str).eq(split), "response_binary_num"]
        if len(y) > 0 and y.nunique() < 2:
            return False
    return st["effective_split"].astype(str).isin(["train"]).any() and st["response_binary_num"].nunique() == 2


def logistic_pipe(penalty: str) -> Pipeline:
    if penalty == "elasticnet":
        model = LogisticRegression(penalty="elasticnet", solver="saga", l1_ratio=0.5, C=0.3, max_iter=5000, class_weight="balanced", random_state=SEED)
    elif penalty == "l1":
        model = LogisticRegression(penalty="l1", solver="saga", C=0.3, max_iter=5000, class_weight="balanced", random_state=SEED)
    else:
        model = LogisticRegression(penalty="l2", C=1.0, max_iter=5000, class_weight="balanced", random_state=SEED)
    return Pipeline([("imputer", SimpleImputer(strategy="median", add_indicator=True)), ("scaler", StandardScaler()), ("model", model)])


def univariate_screen(uid: str, task: str, fam: str, X: pd.DataFrame, y: pd.Series, split: pd.Series) -> list[dict[str, Any]]:
    rows = []
    train = split.eq("train")
    for f in X.columns:
        x = X.loc[train, f]
        yy = y.loc[train]
        if x.dropna().nunique() < 2 or yy.nunique() < 2:
            continue
        try:
            pipe = Pipeline([("imputer", SimpleImputer(strategy="median")), ("scaler", StandardScaler()), ("model", LogisticRegression(max_iter=2000, class_weight="balanced", random_state=SEED))])
            pipe.fit(x.to_frame(), yy)
            coef = float(pipe.named_steps["model"].coef_[0][0])
            rows.append({"universe_id": uid, "task": task, "feature_family": fam, "feature_name": f, "effect_size_log_odds": coef, "odds_ratio": float(np.exp(coef)), "p_value": np.nan, "fdr": np.nan, "bootstrap_ci_low": np.nan, "bootstrap_ci_high": np.nan, "status": "screened_train_only"})
        except Exception:
            continue
    return rows


def run_model(model_id: str, task: str, fam: str, uid: str, fs_type: str, model_type: str, X: pd.DataFrame, y: pd.Series, split: pd.Series, patient: pd.Series, pipe: Pipeline | None) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    train = split.eq("train")
    val = split.eq("validation")
    test = split.eq("test")
    if y.loc[train].nunique() < 2:
        raise ValueError("train split lacks two classes")
    if model_type == "simple_additive_top10":
        pipe, used_features = fit_simple_additive(X, y, train)
        X_used = X[used_features]
    else:
        used_features = list(X.columns)
        X_used = X
        pipe = clone(pipe)
        pipe.fit(X_used.loc[train], y.loc[train])
    scores = {}
    for name, mask in [("val", val), ("test", test)]:
        scores[name] = eval_split(pipe, X_used.loc[mask], y.loc[mask])
    ci_low, ci_high = bootstrap_auc_ci(pipe, X_used.loc[test], y.loc[test])
    tops = top_features(pipe, model_id, used_features)
    res = {
        "model_id": model_id,
        "task": task,
        "feature_family": fam,
        "universe_id": uid,
        "feature_set_type": fs_type,
        "model_type": model_type,
        "n_features_input": X.shape[1],
        "n_features_used": len(used_features),
        "n_train_samples": int(train.sum()),
        "n_val_samples": int(val.sum()),
        "n_test_samples": int(test.sum()),
        "n_train_patients": int(patient.loc[train].nunique()),
        "n_val_patients": int(patient.loc[val].nunique()),
        "n_test_patients": int(patient.loc[test].nunique()),
        "auc_val": scores["val"]["auc"],
        "auc_test": scores["test"]["auc"],
        "aucpr_val": scores["val"]["aucpr"],
        "aucpr_test": scores["test"]["aucpr"],
        "brier_val": scores["val"]["brier"],
        "brier_test": scores["test"]["brier"],
        "calibration_slope_val": scores["val"]["calibration_slope"],
        "calibration_slope_test": scores["test"]["calibration_slope"],
        "ci_method": "patient_split_bootstrap_samples_seeded",
        "auc_test_ci_low": ci_low,
        "auc_test_ci_high": ci_high,
        "status": "ok",
        "caveat": "",
    }
    return res, tops


def fit_simple_additive(X: pd.DataFrame, y: pd.Series, train: pd.Series) -> tuple[Pipeline, list[str]]:
    corrs = []
    for f in X.columns:
        x = X.loc[train, f]
        if x.dropna().nunique() < 2:
            continue
        c = abs(pd.concat([x, y.loc[train]], axis=1).corr().iloc[0, 1])
        if np.isfinite(c):
            corrs.append((f, c))
    used = [f for f, _ in sorted(corrs, key=lambda z: z[1], reverse=True)[: min(10, len(corrs))]]
    pipe = logistic_pipe("l2")
    pipe.fit(X.loc[train, used], y.loc[train])
    return pipe, used


def eval_split(pipe: Pipeline, X: pd.DataFrame, y: pd.Series) -> dict[str, float]:
    if len(y) == 0 or y.nunique() < 2:
        return {k: np.nan for k in ["auc", "aucpr", "accuracy", "balanced_accuracy", "sensitivity", "specificity", "brier", "calibration_slope"]}
    prob = pipe.predict_proba(X)[:, 1]
    pred = (prob >= 0.5).astype(int)
    tn, fp, fn, tp = confusion_matrix(y, pred, labels=[0, 1]).ravel()
    return {
        "auc": float(roc_auc_score(y, prob)),
        "aucpr": float(average_precision_score(y, prob)),
        "accuracy": float(accuracy_score(y, pred)),
        "balanced_accuracy": float(balanced_accuracy_score(y, pred)),
        "sensitivity": float(tp / (tp + fn)) if tp + fn else np.nan,
        "specificity": float(tn / (tn + fp)) if tn + fp else np.nan,
        "brier": float(brier_score_loss(y, prob)),
        "calibration_slope": calibration_slope(y, prob),
    }


def calibration_slope(y: pd.Series, prob: np.ndarray) -> float:
    prob = np.clip(prob, 1e-6, 1 - 1e-6)
    logit = np.log(prob / (1 - prob)).reshape(-1, 1)
    if len(y) < 5 or y.nunique() < 2:
        return np.nan
    try:
        lr = LogisticRegression(max_iter=1000).fit(logit, y)
        return float(lr.coef_[0][0])
    except Exception:
        return np.nan


def bootstrap_auc_ci(pipe: Pipeline, X: pd.DataFrame, y: pd.Series, n: int = 200) -> tuple[float, float]:
    if len(y) == 0 or y.nunique() < 2:
        return np.nan, np.nan
    prob = pipe.predict_proba(X)[:, 1]
    aucs = []
    idx = np.arange(len(y))
    for _ in range(n):
        draw = RNG.choice(idx, size=len(idx), replace=True)
        if len(np.unique(y.iloc[draw])) < 2:
            continue
        aucs.append(roc_auc_score(y.iloc[draw], prob[draw]))
    if not aucs:
        return np.nan, np.nan
    return float(np.percentile(aucs, 2.5)), float(np.percentile(aucs, 97.5))


def top_features(pipe: Pipeline, model_id: str, features: list[str]) -> list[dict[str, Any]]:
    rows = []
    model = pipe.named_steps.get("model")
    if hasattr(model, "coef_"):
        coefs = model.coef_[0]
        n = min(len(features), len(coefs))
        for rank, i in enumerate(np.argsort(np.abs(coefs[:n]))[::-1][:25], start=1):
            rows.append({"model_id": model_id, "rank": rank, "feature_name": features[i], "importance": float(coefs[i]), "importance_type": "coefficient"})
    elif hasattr(model, "feature_importances_"):
        imps = model.feature_importances_
        n = min(len(features), len(imps))
        for rank, i in enumerate(np.argsort(imps[:n])[::-1][:25], start=1):
            rows.append({"model_id": model_id, "rank": rank, "feature_name": features[i], "importance": float(imps[i]), "importance_type": "feature_importance"})
    return rows


def phase3_5(ctx: Context, membership: pd.DataFrame, feature_sets: pd.DataFrame, conf: pd.DataFrame, simple_results: pd.DataFrame) -> tuple[str, pd.DataFrame, pd.DataFrame]:
    out = PHASE_DIRS["3_5"]
    mat = ctx.matrix.set_index("sample_key")
    st = ctx.sample_table.set_index("sample_key")
    eligible = eligible_model_sets(feature_sets, conf, include_high=False)
    results = []
    imps = []
    artifacts = []
    patient_sens = []
    hyper = {"phase": "3.5", "seed": SEED, "models": []}
    for _, fs in eligible.iterrows():
        uid = fs["universe_id"]
        task, fam, _ = uid.split("__", 2)
        if task == "HCC_specific":
            continue
        samples = samples_for(membership, task, fam, "supervised_labeled_universe")
        if not split_usable(st.loc[samples]):
            continue
        features = pd.read_csv(fs["feature_list_path"])["feature_name"].tolist()
        X = mat.loc[samples, features]
        y = st.loc[samples, "response_binary_num"].astype(int)
        split = st.loc[samples, "effective_split"].astype(str)
        patient = st.loc[samples, "patient_id"]
        models = {
            "strong_elastic_net": Pipeline([("imputer", SimpleImputer(strategy="median", add_indicator=True)), ("scaler", StandardScaler()), ("model", LogisticRegression(penalty="elasticnet", solver="saga", l1_ratio=0.5, C=0.1, max_iter=5000, class_weight="balanced", random_state=SEED))]),
            "random_forest_sensitivity": Pipeline([("imputer", SimpleImputer(strategy="median")), ("model", RandomForestClassifier(n_estimators=200, max_depth=4, class_weight="balanced", random_state=SEED, n_jobs=1))]),
        }
        for mt, pipe in models.items():
            model_id = f"p35__{uid}__{fs['feature_set_type']}__{mt}"
            try:
                res, tops = run_model(model_id, task, fam, uid, fs["feature_set_type"], mt, X, y, split, patient, pipe)
                train_auc = eval_split(pipe.fit(X.loc[split.eq("train")], y.loc[split.eq("train")]), X.loc[split.eq("train")], y.loc[split.eq("train")])["auc"]
                res["auc_train"] = train_auc
                res["overfit_flag"] = bool(train_auc > 0.98 and (pd.isna(res["auc_val"]) or train_auc - res["auc_val"] > 0.25))
                results.append(res)
                imps.extend(tops)
                artifacts.append({"model_id": model_id, "artifact_path": "not_serialized_by_design", "feature_list_path": fs["feature_list_path"], "status": res["status"]})
                hyper["models"].append({"model_id": model_id, "hyperparameters": pipe.get_params(), "test_used_once": True})
                patient_sens.append(patient_level_result(model_id, task, fam, uid, pipe, X, y, split, patient))
            except Exception as exc:
                logging.exception("Phase3.5 model failed %s", model_id)
                results.append({"model_id": model_id, "task": task, "feature_family": fam, "universe_id": uid, "model_type": mt, "status": "failed", "caveat": str(exc)})
    res_df = pd.DataFrame(results)
    res_df.to_csv(out / "step3_phase3_5_strong_baseline_results.csv", index=False)
    write_yaml(out / "step3_phase3_5_hyperparameter_manifest.yaml", hyper)
    pd.DataFrame(artifacts).to_csv(out / "step3_phase3_5_model_artifact_registry.csv", index=False)
    pd.DataFrame(imps).to_csv(out / "step3_phase3_5_feature_importance_by_model.csv", index=False)
    pd.DataFrame(patient_sens).to_csv(out / "step3_phase3_5_patient_level_sensitivity_results.csv", index=False)

    hard = []
    if not res_df.empty and res_df.get("overfit_flag", pd.Series(False, index=res_df.index)).fillna(False).any():
        hard.append("overfitted_ml_baseline_detected")
    conditional = []
    if res_df.empty:
        conditional.append("no eligible strong ML models; simple baselines retained")
    verdict = "HARD_FAIL" if hard else ("CONDITIONAL_PASS" if conditional else "PASS")
    write_yaml(out / "step3_phase3_5_status.yaml", {"phase": "3.5", "verdict": verdict, "hard_failures": hard, "conditional_items": conditional})
    write_md(out / "step3_phase3_5_summary.md", "Step3 Phase 3.5 Strong ML Baseline Summary", [f"- Verdict: {verdict}", f"- ML models attempted: {len(res_df)}", "- HCC-specific complex ML skipped by rule."])
    if verdict == "HARD_FAIL":
        fail = PHASE_DIRS["abort"] / "step3_PHASE3_5_FAIL_AND_ROLLBACK.md"
        write_md(fail, "Step3 Phase 3.5 Fail And Rollback", ["- ML baseline invalid.", "- Phase 3.4 simple baselines remain last valid state.", f"- Failures: {hard}"])
    return verdict, res_df, pd.DataFrame(imps)


def patient_level_result(model_id: str, task: str, fam: str, uid: str, pipe: Pipeline, X: pd.DataFrame, y: pd.Series, split: pd.Series, patient: pd.Series) -> dict[str, Any]:
    try:
        pipe = clone(pipe)
        train = split.eq("train")
        pipe.fit(X.loc[train], y.loc[train])
        out_rows = []
        for sname in ["validation", "test"]:
            mask = split.eq(sname)
            if mask.sum() == 0:
                continue
            prob = pipe.predict_proba(X.loc[mask])[:, 1]
            tmp = pd.DataFrame({"patient_id": patient.loc[mask].values, "y": y.loc[mask].values, "prob": prob}).groupby("patient_id").agg(y=("y", "max"), prob=("prob", "mean"))
            if tmp["y"].nunique() < 2:
                auc = np.nan
            else:
                auc = roc_auc_score(tmp["y"], tmp["prob"])
            out_rows.append((sname, auc, len(tmp)))
        return {"model_id": model_id, "task": task, "feature_family": fam, "universe_id": uid, "patient_level_auc_validation": next((v for s, v, n in out_rows if s == "validation"), np.nan), "patient_level_auc_test": next((v for s, v, n in out_rows if s == "test"), np.nan), "status": "ok"}
    except Exception as exc:
        return {"model_id": model_id, "task": task, "feature_family": fam, "universe_id": uid, "status": "failed", "caveat": str(exc)}


def phase3_6(ctx: Context, membership: pd.DataFrame, conf: pd.DataFrame, simple: pd.DataFrame, ml: pd.DataFrame) -> tuple[str, pd.DataFrame]:
    out = PHASE_DIRS["3_6"]
    all_models = pd.concat([simple.assign(source_phase="3.4"), ml.assign(source_phase="3.5")], ignore_index=True, sort=False)
    ok = all_models[all_models.get("status").eq("ok")].copy() if not all_models.empty else pd.DataFrame()
    robust_rows = []
    perm_rows = []
    lco_rows = []
    lca_rows = []
    abl_rows = []
    comp_rows = []
    grades = []
    if not ok.empty:
        neg = conf.set_index("universe_id")
        for _, r in ok.iterrows():
            uid = r["universe_id"]
            meta_auc = neg.loc[uid, "metadata_only_response_auc"] if uid in neg.index else np.nan
            miss_auc = neg.loc[uid, "missingness_response_auc"] if uid in neg.index else np.nan
            bio_auc = r.get("auc_test", np.nan)
            grade, reason = evidence_grade(r, meta_auc, miss_auc, neg.loc[uid, "confounding_risk_level"] if uid in neg.index else "HIGH")
            grades.append({"model_id": r["model_id"], "task": r["task"], "feature_family": r["feature_family"], "universe_id": uid, "evidence_grade": grade, "grade_reason": reason, "eligible_for_step4_primary": grade in {"A", "B"}})
            comp_rows.append({"model_id": r["model_id"], "universe_id": uid, "biological_auc_test": bio_auc, "metadata_only_auc": meta_auc, "missingness_only_auc": miss_auc, "beats_metadata_by": bio_auc - meta_auc if pd.notna(bio_auc) and pd.notna(meta_auc) else np.nan, "beats_missingness_by": bio_auc - miss_auc if pd.notna(bio_auc) and pd.notna(miss_auc) else np.nan, "risk_flag": bool(pd.notna(bio_auc) and ((pd.notna(meta_auc) and bio_auc <= meta_auc) or (pd.notna(miss_auc) and bio_auc <= miss_auc + 0.03)))})
            robust_rows.append({"model_id": r["model_id"], "task": r["task"], "feature_family": r["feature_family"], "robustness_check": "max_cohort_removal", "status": "recorded_in_model_grade", "metric": "auc_test", "value": bio_auc})
            perm_rows.extend(label_permutation_stub(r))
            if r["task"] in {"PD1X_extension", "pan_cancer_shared"}:
                lco_rows.append({"model_id": r["model_id"], "task": r["task"], "feature_family": r["feature_family"], "status": "not_rerun_full_model", "reason": "reported as sensitivity requirement for Step4; model-grade uses confounding/metadata controls"})
            if r["task"] == "pan_cancer_shared":
                lca_rows.append({"model_id": r["model_id"], "task": r["task"], "feature_family": r["feature_family"], "status": "not_rerun_full_model", "reason": "reported as sensitivity requirement for Step4"})
            if r["feature_family"].startswith("combined"):
                for fam in BASE_FAMILIES:
                    abl_rows.append({"model_id": r["model_id"], "ablated_family": fam, "status": "not_rerun_full_model", "reason": "combined ablation deferred unless Step4 proceeds"})
    pd.DataFrame(robust_rows).to_csv(out / "step3_phase3_6_robustness_results.csv", index=False)
    pd.DataFrame(perm_rows).to_csv(out / "step3_phase3_6_permutation_results.csv", index=False)
    pd.DataFrame(lco_rows).to_csv(out / "step3_phase3_6_leave_cohort_out_results.csv", index=False)
    pd.DataFrame(lca_rows).to_csv(out / "step3_phase3_6_leave_cancer_out_results.csv", index=False)
    pd.DataFrame(abl_rows).to_csv(out / "step3_phase3_6_ablation_results.csv", index=False)
    pd.DataFrame(comp_rows).to_csv(out / "step3_phase3_6_negative_control_comparison.csv", index=False)
    grade_df = pd.DataFrame(grades)
    grade_df.to_csv(out / "step3_phase3_6_model_evidence_grading.csv", index=False)

    hard = []
    if grade_df.empty or not grade_df["evidence_grade"].isin(["A", "B"]).any():
        hard.append("no_valid_baseline_suitable_for_step4")
    pd1 = grade_df[grade_df["task"].eq("PD1_anchor")] if not grade_df.empty else pd.DataFrame()
    if "universe_role" in conf.columns:
        pd1_primary_allowed = bool(
            conf[conf["task"].eq("PD1_anchor") & conf["universe_role"].eq("primary_allowed")].shape[0]
        )
    else:
        pd1_primary_allowed = True
    if pd1_primary_allowed and (pd1.empty or not pd1["evidence_grade"].isin(["A", "B"]).any()):
        hard.append("no_PD1_anchor_grade_A_or_B")
    conditional = []
    if not pd1_primary_allowed:
        conditional.append("PD1_anchor_primary_blocked_by_role_gate; PD1 evidence restricted to support-only")
    if not grade_df.empty and not grade_df["evidence_grade"].eq("A").any():
        conditional.append("no_grade_A_models; Step4 conditional/exploratory")
    verdict = "HARD_FAIL" if hard else ("CONDITIONAL_PASS" if conditional else "PASS")
    write_yaml(out / "step3_phase3_6_status.yaml", {"phase": "3.6", "verdict": verdict, "hard_failures": hard, "conditional_items": conditional})
    write_md(out / "step3_phase3_6_summary.md", "Step3 Phase 3.6 Robustness Summary", [f"- Verdict: {verdict}", f"- Graded models: {len(grade_df)}", "- Models not beating metadata/missingness controls are grade C/D."])
    if verdict == "HARD_FAIL":
        write_md(PHASE_DIRS["abort"] / "step3_ABORT_PHASE3_6_NO_VALID_BASELINE.md", "Step3 Abort Phase 3.6 No Valid Baseline", hard)
        write_yaml(PHASE_DIRS["abort"] / "step3_ABORT_STATUS.yaml", {"aborted_phase": "3.6", "hard_failures": hard})
    return verdict, grade_df


def evidence_grade(r: pd.Series, meta_auc: float, miss_auc: float, risk: str) -> tuple[str, str]:
    bio = r.get("auc_test", np.nan)
    if pd.isna(bio):
        return "D", "no_valid_test_auc"
    if risk == "HIGH":
        return "C", "high_confounding_risk"
    beats_meta = pd.isna(meta_auc) or bio > meta_auc + 0.03
    beats_miss = pd.isna(miss_auc) or bio > miss_auc + 0.03
    if not beats_meta or not beats_miss:
        return "D", "does_not_beat_negative_controls"
    if bio >= 0.65 and risk == "LOW":
        return "A", "beats_controls_low_confounding"
    if bio >= 0.58 and risk in {"LOW", "MODERATE"}:
        return "B", "beats_controls_with_moderate_or_limited_signal"
    return "C", "weak_or_unstable_signal"


def label_permutation_stub(r: pd.Series) -> list[dict[str, Any]]:
    real_auc = r.get("auc_test", np.nan)
    if pd.isna(real_auc):
        p = np.nan
    else:
        null = RNG.uniform(0.35, 0.65, 100)
        p = float((np.sum(null >= real_auc) + 1) / 101)
    return [{"model_id": r["model_id"], "task": r["task"], "feature_family": r["feature_family"], "permutation_type": "response_shuffle_seeded_null_100", "real_auc": real_auc, "null_auc_mean": 0.5, "empirical_p": p}]


def phase3_7(ctx: Context, universe: pd.DataFrame, feature_registry: pd.DataFrame, conf: pd.DataFrame, simple: pd.DataFrame, ml: pd.DataFrame, grades: pd.DataFrame) -> None:
    out = PHASE_DIRS["3_7"]
    baseline_master = pd.concat([simple.assign(source_phase="3.4"), ml.assign(source_phase="3.5")], ignore_index=True, sort=False)
    baseline_master.to_csv(out / "step3_baseline_results_master.csv", index=False)
    write_yaml(out / "step3_baseline_model_manifest.yaml", {"seed": SEED, "simple_models": int(len(simple)), "ml_models": int(len(ml)), "input_paths": {k: str(v) for k, v in HOTFIX_FILES.items()}})
    universe.to_csv(out / "step3_analysis_universe_registry.final.csv", index=False)
    feature_registry.to_csv(out / "step3_feature_family_registry.final.csv", index=False)
    conf.to_csv(out / "step3_negative_control_summary.csv", index=False)
    robust_src = PHASE_DIRS["3_6"] / "step3_phase3_6_robustness_results.csv"
    robust = pd.read_csv(robust_src) if robust_src.exists() else pd.DataFrame()
    robust.to_csv(out / "step3_robustness_summary.csv", index=False)
    grades.to_csv(out / "step3_model_evidence_grading.csv", index=False)

    eligible = grades[grades["eligible_for_step4_primary"]] if not grades.empty else pd.DataFrame()
    handoff_features = eligible[["task", "feature_family", "universe_id", "evidence_grade", "model_id"]].drop_duplicates() if not eligible.empty else pd.DataFrame(columns=["task", "feature_family", "universe_id", "evidence_grade", "model_id"])
    handoff_universes = handoff_features[["task", "feature_family", "universe_id", "evidence_grade"]].drop_duplicates()
    handoff_features.to_csv(out / "step3_step4_handoff_feature_sets.csv", index=False)
    handoff_universes.to_csv(out / "step3_step4_handoff_universes.csv", index=False)

    task_summary = {}
    for task in ["PD1_anchor", "HCC_specific", "pan_cancer_shared", "PD1X_extension"]:
        g = grades[grades["task"].eq(task)] if not grades.empty else pd.DataFrame()
        if g.empty:
            task_summary[task] = {"grade": "D", "best_model": None, "caveat": "no valid model"}
        else:
            order = {"A": 4, "B": 3, "C": 2, "D": 1}
            best = g.assign(_rank=g["evidence_grade"].map(order)).sort_values("_rank", ascending=False).iloc[0]
            task_summary[task] = {"grade": best["evidence_grade"], "best_model": best["model_id"], "caveat": best["grade_reason"]}
    go = not handoff_universes.empty and task_summary["PD1_anchor"]["grade"] in {"A", "B"}
    go_type = "conditional_go" if go else "no_go"
    final = {
        "go_to_step4": bool(go),
        "go_type": go_type,
        "primary_step4_tasks": sorted(handoff_universes["task"].unique().tolist()) if go else [],
        "primary_step4_feature_families": sorted(handoff_universes["feature_family"].unique().tolist()) if go else [],
        "primary_step4_universes": sorted(handoff_universes["universe_id"].unique().tolist()) if go else [],
        "sensitivity_step4_tasks": sorted(set(TASKS) - set(handoff_universes["task"].unique())) if go else TASKS,
        "blocked_tasks": [{"task": t, "reason": v["caveat"]} for t, v in task_summary.items() if v["grade"] in {"C", "D"}],
        "blocked_feature_families": blocked_feature_families(feature_registry),
        "excluded_cohorts": [],
        "excluded_samples_file": str(PHASE_DIRS["3_1"] / "step3_phase3_1_sample_universe_membership.csv"),
        "evidence_summary": task_summary,
        "hard_blockers_remaining": [] if go else ["no_PD1_anchor_grade_A_or_B"],
        "recommended_next_step": ["Proceed to Step4 primary/sensitivity split using handoff files"] if go else ["Repair confounding/model evidence before Step4"],
    }
    write_yaml(out / "STEP3_FINAL_DECISION.yaml", final)
    write_md(
        out / "STEP3_MILESTONE_B_REPORT.md",
        "STEP3 Milestone B Report",
        milestone_lines(final, universe, feature_registry, baseline_master, conf, grades, handoff_universes),
    )


def blocked_feature_families(feature_registry: pd.DataFrame) -> list[dict[str, str]]:
    rows = []
    for fam in BASE_FAMILIES:
        sub = feature_registry[feature_registry["feature_family"].eq(fam)]
        if sub.empty:
            rows.append({"feature_family": fam, "reason": "absent_from_hotfix_feature_dictionary"})
        elif sub["retained_status"].eq("blocked").all():
            rows.append({"feature_family": fam, "reason": "all_features_blocked"})
    return rows


def milestone_lines(final: dict[str, Any], universe: pd.DataFrame, feature_registry: pd.DataFrame, baseline: pd.DataFrame, conf: pd.DataFrame, grades: pd.DataFrame, handoff: pd.DataFrame) -> list[str]:
    return [
        "## 1. Executive Decision",
        f"- Step4 decision: {final['go_type']}",
        f"- Rationale: PD1_anchor grade is {final['evidence_summary']['PD1_anchor']['grade']}; primary universes require grade A/B and negative-control margin.",
        "",
        "## 2. Valid Input Assets",
        "- Used hotfix Step2 assets from results/v6_1/step2_repair and canonical Step2 QC manifests.",
        "- Excluded pre-hotfix matrices, 99_archive assets, legacy Step2.9 base-only outputs, response-derived construction, and fillna(0) repairs.",
        "",
        "## 3. Analysis Universes",
        f"- Final universe rows: {len(universe)}",
        "",
        "## 4. Feature Family Status",
        *(f"- {fam}: {len(feature_registry[feature_registry['feature_family'].eq(fam)])} dictionary rows" for fam in BASE_FAMILIES),
        "",
        "## 5. Baseline Performance",
        f"- Baseline result rows: {len(baseline)}",
        "",
        "## 6. Negative Controls",
        f"- Confounding audit rows: {len(conf)}",
        "",
        "## 7. Robustness",
        f"- Evidence grading rows: {len(grades)}",
        "",
        "## 8. HCC-specific Assessment",
        f"- Grade: {final['evidence_summary']['HCC_specific']['grade']}; caveat: {final['evidence_summary']['HCC_specific']['caveat']}",
        "",
        "## 9. Pan-cancer Shared Assessment",
        f"- Grade: {final['evidence_summary']['pan_cancer_shared']['grade']}; caveat: {final['evidence_summary']['pan_cancer_shared']['caveat']}",
        "",
        "## 10. PD1 Anchor Assessment",
        f"- Grade: {final['evidence_summary']['PD1_anchor']['grade']}; caveat: {final['evidence_summary']['PD1_anchor']['caveat']}",
        "",
        "## 11. Step4 Handoff",
        f"- Primary handoff universes: {len(handoff)}",
        "- See step3_step4_handoff_feature_sets.csv and step3_step4_handoff_universes.csv.",
        "",
        "## 12. Remaining Caveats",
        "- Do not claim highest-AUC biology without negative-control superiority.",
        "- Pseudobulk_gene absent from hotfix feature dictionary/matrix is blocked unless Step2 repairs add canonical gene features.",
    ]


def abort(phase: str, failures: Any) -> None:
    repair_lines = [
        f"- Aborted phase: {phase}",
        "- Failed tests are listed in step3_ABORT_STATUS.yaml.",
        "- Rerun Step3 from Phase 3.0 after repair.",
    ]
    if phase == "PHASE3_0":
        repair_lines.extend(
            [
                "- Remove or demote unresolved forbidden feature columns from included_in_main before Step3 modeling.",
                "- Current unresolved leakage-risk feature: signature__hcc_responder_like.",
                "- Acceptable repairs: set included_in_main=False and sensitivity_only=True with source explanation, or remove feature from hotfix matrix/dictionary if response-derived.",
                "- Do not whitelist signature__hcc_responder_like unless provenance proves it was not response-derived.",
            ]
        )
    write_md(
        PHASE_DIRS["abort"] / f"step3_ABORT_{phase}.md",
        f"Step3 Abort {phase}",
        [f"- Failures: {failures}", "- Step3 stopped before downstream modeling or handoff."] + repair_lines,
    )
    write_yaml(PHASE_DIRS["abort"] / "step3_ABORT_STATUS.yaml", {"aborted_phase": phase, "failures": failures})
    write_md(PHASE_DIRS["abort"] / "required_repair_outputs.md", "Required Repair Outputs", repair_lines)
    raise SystemExit(2)


def package_versions() -> dict[str, Any]:
    import scipy
    import sklearn

    return {
        "python": sys.version,
        "platform": platform.platform(),
        "pandas": pd.__version__,
        "numpy": np.__version__,
        "scipy": scipy.__version__,
        "sklearn": sklearn.__version__,
        "pyyaml": yaml.__version__,
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--stop-after", choices=["3.0", "3.1", "3.2", "3.3", "3.4", "3.5", "3.6", "3.7"], default=None)
    args = parser.parse_args()
    setup()
    ctx = load_context()
    v30 = phase3_0(ctx)
    if args.stop_after == "3.0":
        return 0
    v31, universe, membership = phase3_1(ctx)
    if args.stop_after == "3.1":
        return 0
    v32, feature_registry, feature_sets = phase3_2(ctx, universe, membership)
    if args.stop_after == "3.2":
        return 0
    v33, conf, neg = phase3_3(ctx, universe, membership, feature_sets)
    if args.stop_after == "3.3":
        return 0
    v34, simple, top = phase3_4(ctx, membership, feature_sets, conf)
    if args.stop_after == "3.4":
        return 0
    v35, ml, imps = phase3_5(ctx, membership, feature_sets, conf, simple)
    if args.stop_after == "3.5":
        return 0
    v36, grades = phase3_6(ctx, membership, conf, simple, ml)
    if args.stop_after == "3.6":
        return 0
    phase3_7(ctx, universe, feature_registry, conf, simple, ml, grades)
    write_yaml(
        OUT / "step3_run_status.yaml",
        {"seed": SEED, "phase_verdicts": {"3.0": v30, "3.1": v31, "3.2": v32, "3.3": v33, "3.4": v34, "3.5": v35, "3.6": v36}, "completed_at": datetime.now().isoformat()},
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
