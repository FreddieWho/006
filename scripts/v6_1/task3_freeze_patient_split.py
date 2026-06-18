from __future__ import annotations

import hashlib
from pathlib import Path
from typing import Dict, List

import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
STEP1_DIR = ROOT / "results" / "v6_1" / "step1"

PATIENT_PATH = STEP1_DIR / "patient_metadata_master_v6_1.csv"
SAMPLE_PATH = STEP1_DIR / "sample_metadata_master_v6_1.csv"
COHORT_SPLIT_PATH = STEP1_DIR / "frozen_cohort_inclusion_v6_1.csv"
FLAGS_PATH = STEP1_DIR / "treatment_context_flags.csv"

OUTPUT_CSV = STEP1_DIR / "frozen_patient_split_v6_1.csv"
OUTPUT_MD = STEP1_DIR / "split_design_report.md"


def clean_str(value: object) -> str:
    if pd.isna(value):
        return ""
    text = str(value).strip()
    if text.lower() in {"", "nan", "na", "n/a", "none", "null"}:
        return ""
    return text


def truthy(value: object) -> bool:
    return str(value).strip().lower() in {"true", "1", "yes", "y", "t"}


def stable_hash_key(cohort_id: str, patient_id: str) -> str:
    return hashlib.md5(f"{cohort_id}::{patient_id}".encode("utf-8")).hexdigest()


def count_anchor_supervision(patient: pd.DataFrame, cohort_freeze: pd.DataFrame) -> Dict[str, object]:
    anchor_main = cohort_freeze[
        (cohort_freeze["treatment_context"] == "PD1_ICI_anchor")
        & (cohort_freeze["inclusion_status"] == "main")
    ]["cohort_id"].tolist()
    sub = patient[patient["cohort_id"].isin(anchor_main)].copy()
    pre = sub[sub["has_pre_sample"].map(truthy)].copy()
    pre_known = pre[pre["best_response_binary"].isin(["responder", "non_responder"])].copy()
    available = bool(pre_known["cohort_id"].nunique() >= 2 and len(pre_known) >= 40)
    return {
        "anchor_main_cohorts": anchor_main,
        "anchor_patient_total": int(len(sub)),
        "anchor_pre_total": int(len(pre)),
        "anchor_pre_known_total": int(len(pre_known)),
        "anchor_pre_known_cohorts": int(pre_known["cohort_id"].nunique()),
        "supervised_anchor_available": available,
    }


def summarize_flags(flags: pd.DataFrame) -> pd.DataFrame:
    rank = {"high": 3, "medium": 2, "low": 1}

    def patient_confidence(series: pd.Series) -> str:
        values = [clean_str(x).lower() for x in series.tolist() if clean_str(x)]
        if not values:
            return "unknown"
        return min(values, key=lambda x: rank.get(x, 0))

    def collect_unique(series: pd.Series) -> str:
        vals = sorted({clean_str(x) for x in series.tolist() if clean_str(x)})
        return ";".join(vals) if vals else "none"

    out = (
        flags.groupby(["cohort_id", "patient_id"])
        .agg(
            patient_confidence_level=("confidence_level", patient_confidence),
            treatment_context_from_flags=("treatment_context", collect_unique),
            flag_reason_summary=("reason_for_assignment", collect_unique),
            flag_ambiguity_summary=("ambiguity_note", collect_unique),
        )
        .reset_index()
    )
    return out


def assign_main_splits(main_df: pd.DataFrame) -> pd.DataFrame:
    rows: List[pd.DataFrame] = []
    for cohort_id, sub in main_df.groupby("cohort_id", sort=True):
        sub = sub.copy()
        sub["hash_key"] = [stable_hash_key(cohort_id, pid) for pid in sub["patient_id"].astype(str)]
        sub = sub.sort_values("hash_key").reset_index(drop=True)
        n = len(sub)

        if n < 8:
            n_test = 0
            n_val = 0
            split_mode = "small_cohort_all_train"
        elif n < 15:
            n_test = 1
            n_val = 1
            split_mode = "small_cohort_minimal_val_test"
        else:
            n_test = max(1, round(n * 0.2))
            n_val = max(1, round(n * 0.2))
            if n - n_test - n_val < 3:
                n_test = max(1, n_test - 1)
            split_mode = "cohort_hash_60_20_20"

        splits: List[str] = []
        for i in range(n):
            if i < n_test:
                splits.append("test")
            elif i < n_test + n_val:
                splits.append("validation")
            else:
                splits.append("train")

        sub["split_label"] = splits
        sub["fold_id"] = [f"fold_{i % 5}" for i in range(n)]
        sub["split_reason"] = [
            f"main_executable_cohort;mode={split_mode};cohort_patient_n={n};paired_locked_by_patient"
            for _ in range(n)
        ]
        rows.append(sub)

    return pd.concat(rows, ignore_index=True) if rows else main_df.copy()


def build_report(df: pd.DataFrame, anchor_stats: Dict[str, object]) -> str:
    lines: List[str] = []
    lines.append("# split_design_report (v6.1 Task3)")
    lines.append("")
    lines.append("## 设计原则")
    lines.append("- 运行时 split 键采用 `cohort_id + patient_id`，因为 `patient_id` 在跨 cohort 情况下不是全局唯一。")
    lines.append("- 同一 patient 的所有 sample 强制共享同一 split/fold。")
    lines.append("- paired pre/post 通过 patient-level 冻结自然绑定，不允许跨 split。")
    lines.append("- `external_anchor_only` 不进入主训练，统一冻结到 `external`。")
    lines.append("- `high_confounding_support` 默认只进入 `sensitivity`。")
    lines.append("- `excluded_pending_metadata` 的已入库 patient 也保守放入 `sensitivity`，避免未清理 cohort 混入主训练。")
    lines.append("- `PD1_ICI_anchor` 当前采用 unsupervised split design；`supervised_anchor_available=False`。")
    lines.append("")

    lines.append("## Supervised Anchor Availability")
    lines.append(f"- anchor main cohorts: `{';'.join(anchor_stats['anchor_main_cohorts']) or 'none'}`")
    lines.append(f"- anchor patients total: `{anchor_stats['anchor_patient_total']}`")
    lines.append(f"- anchor pre-treatment patients: `{anchor_stats['anchor_pre_total']}`")
    lines.append(f"- anchor pre-treatment patients with known binary response: `{anchor_stats['anchor_pre_known_total']}`")
    lines.append(f"- anchor pre-treatment known-label cohorts: `{anchor_stats['anchor_pre_known_cohorts']}`")
    lines.append(f"- supervised_anchor_available: `{anchor_stats['supervised_anchor_available']}`")
    lines.append("- judgment: current anchor supervision is too narrow because known pre-treatment labels come from only one executable cohort (`TASK01`).")
    lines.append("")

    lines.append("## Split Summary")
    for key, count in df["split_label"].value_counts().sort_index().items():
        lines.append(f"- `{key}`: {int(count)} patients")
    lines.append("")

    lines.append("## Split by Treatment Context")
    context_counts = df.groupby(["split_label", "treatment_context_primary"]).size()
    for (split_label, context), count in context_counts.items():
        lines.append(f"- `{split_label}` / `{context}`: {int(count)}")
    lines.append("")

    lines.append("## Main Cohort Allocation")
    main = df[df["split_label"].isin(["train", "validation", "test"])].copy()
    for cohort_id, sub in main.groupby("cohort_id", sort=True):
        counts = sub["split_label"].value_counts().to_dict()
        lines.append(f"- `{cohort_id}`: n=`{len(sub)}`, train/validation/test=`{counts}`")
    lines.append("")

    lines.append("## Integrity Checks")
    lines.append(f"- unique patient keys in output: `{df['patient_key'].nunique()}`")
    lines.append(f"- rows in output: `{len(df)}`")
    lines.append(f"- patients with paired pre/post: `{int(df['has_paired_pre_post'].map(truthy).sum())}`")
    lines.append("- patient-level assignment guarantees no same-patient leakage across splits.")
    lines.append("- external and sensitivity partitions are isolated from main train/validation/test.")
    lines.append("")

    lines.append("## Follow-up Constraint")
    lines.append("- Any later sample-level matrix must join this file by `cohort_id + patient_id`; direct sample-wise random split is prohibited.")
    return "\n".join(lines) + "\n"


def add_split_alias_column(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    if "split_label" in out.columns:
        out["split"] = out["split_label"].fillna("")
    return out


def main() -> None:
    patient = pd.read_csv(PATIENT_PATH)
    sample = pd.read_csv(SAMPLE_PATH)
    cohort_freeze = pd.read_csv(COHORT_SPLIT_PATH)
    flags = pd.read_csv(FLAGS_PATH)

    sample_patient_counts = (
        sample.groupby(["cohort_id", "patient_id"])["sample_id"].nunique().reset_index(name="sample_count_from_sample_table")
    )
    flag_summary = summarize_flags(flags)
    anchor_stats = count_anchor_supervision(patient, cohort_freeze)

    merge_cols = ["cohort_id", "inclusion_status", "executable_role", "treatment_context"]
    out = patient.merge(cohort_freeze[merge_cols], on="cohort_id", how="left")
    out = out.merge(sample_patient_counts, on=["cohort_id", "patient_id"], how="left")
    out = out.merge(flag_summary, on=["cohort_id", "patient_id"], how="left")
    out["patient_key"] = out["cohort_id"].astype(str) + "::" + out["patient_id"].astype(str)
    out["sample_count_from_sample_table"] = out["sample_count_from_sample_table"].fillna(0).astype(int)
    out["patient_confidence_level"] = out["patient_confidence_level"].fillna("unknown")
    out["treatment_context_from_flags"] = out["treatment_context_from_flags"].fillna("none")
    out["flag_reason_summary"] = out["flag_reason_summary"].fillna("none")
    out["flag_ambiguity_summary"] = out["flag_ambiguity_summary"].fillna("none")
    out["supervised_anchor_available"] = anchor_stats["supervised_anchor_available"]

    out["split_label"] = "sensitivity"
    out["fold_id"] = "not_applicable"
    out["split_reason"] = "default_sensitivity_holding_partition"

    main_mask = out["inclusion_status"] == "main"
    external_mask = out["inclusion_status"] == "external"
    sensitivity_mask = out["inclusion_status"].isin(["sensitivity", "spatial_only", "excluded_pending_metadata"])

    if main_mask.any():
        assigned_main = assign_main_splits(out.loc[main_mask].copy())
        assigned_main = assigned_main[["cohort_id", "patient_id", "split_label", "fold_id", "split_reason"]].copy()
        assigned_main = assigned_main.rename(
            columns={
                "split_label": "split_label_main",
                "fold_id": "fold_id_main",
                "split_reason": "split_reason_main",
            }
        )
        out = out.merge(assigned_main, on=["cohort_id", "patient_id"], how="left")
        out.loc[main_mask, "split_label"] = out.loc[main_mask, "split_label_main"]
        out.loc[main_mask, "fold_id"] = out.loc[main_mask, "fold_id_main"]
        out.loc[main_mask, "split_reason"] = out.loc[main_mask, "split_reason_main"]
        out = out.drop(columns=["split_label_main", "fold_id_main", "split_reason_main"])

    out.loc[external_mask, "split_label"] = "external"
    out.loc[external_mask, "fold_id"] = "not_applicable"
    out.loc[external_mask, "split_reason"] = "external_anchor_only_partition;excluded_from_main_training"

    out.loc[sensitivity_mask, "split_label"] = "sensitivity"
    out.loc[sensitivity_mask, "fold_id"] = "not_applicable"
    out.loc[sensitivity_mask, "split_reason"] = out.loc[sensitivity_mask].apply(
        lambda r: (
            "sensitivity_partition;"
            f"inclusion_status={clean_str(r['inclusion_status'])};"
            f"treatment_context={clean_str(r['treatment_context_primary'])};"
            "not_eligible_for_main_training"
        ),
        axis=1,
    )

    out = out[
        [
            "cohort_id",
            "patient_id",
            "patient_key",
            "disease",
            "treatment_context_primary",
            "inclusion_status",
            "executable_role",
            "sample_count_from_sample_table",
            "has_pre_sample",
            "has_on_treatment_sample",
            "has_post_sample",
            "has_paired_pre_post",
            "best_response_raw",
            "best_response_binary",
            "response_standard_primary",
            "survival_available",
            "number_of_samples",
            "number_of_modalities",
            "patient_confidence_level",
            "treatment_context_from_flags",
            "flag_reason_summary",
            "flag_ambiguity_summary",
            "split_label",
            "fold_id",
            "supervised_anchor_available",
            "split_reason",
            "notes",
        ]
    ].sort_values(["split_label", "cohort_id", "patient_id"]).reset_index(drop=True)
    out = add_split_alias_column(out)

    out.to_csv(OUTPUT_CSV, index=False)
    OUTPUT_MD.write_text(build_report(out, anchor_stats), encoding="utf-8")


if __name__ == "__main__":
    main()
