#!/usr/bin/env python3
"""Response-blind, cohort-level anchor feasibility audit (WP6a).

This stage reads only cohort/design metadata.  It must not load patient-level
response values, outcome values, or responder/failure assignments.
"""

from __future__ import annotations

import argparse
import hashlib
from pathlib import Path
from typing import Iterable

import pandas as pd
import yaml


REGISTRY_COLS = [
    "cohort_id",
    "disease",
    "cancer_group",
    "data_modality",
    "response_label_type",
    "response_label_quality",
    "paired_available",
    "spatial_available",
    "tcr_raw_available",
    "perturb_available",
    "bulk_available",
    "etiology_metadata",
    "b1_etiology_applicable",
    "reaudit_status",
]
ENV_COLS = [
    "cohort_id",
    "cancer_group",
    "treatment_context",
    "response_label_type",
    "response_label_quality",
    "timepoint_schema",
    "paired_available",
    "environment_role",
    "etiology_env_B1",
]

# These are value-bearing fields, not design metadata.  The allowed metadata
# fields above deliberately include response_label_type/quality as schemas.
FORBIDDEN_VALUE_TOKENS = (
    "response_raw",
    "response_binary",
    "outcome_value",
    "survival_time",
    "failure_value",
    "responder_value",
    "endpoint_value",
    "response_direction",
)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def require_columns(path: Path, expected: Iterable[str]) -> None:
    got = list(pd.read_csv(path, nrows=0).columns)
    missing = sorted(set(expected) - set(got))
    if missing:
        raise ValueError(f"{path}: missing required metadata columns: {missing}")
    lowered = [x.lower() for x in got]
    bad = [x for x in lowered if any(t in x for t in FORBIDDEN_VALUE_TOKENS)]
    if bad:
        raise ValueError(f"{path}: forbidden value-bearing columns present: {bad}")


def value_counts(df: pd.DataFrame, column: str) -> dict[str, int]:
    return {str(k): int(v) for k, v in df[column].fillna("<NA>").value_counts(dropna=False).items()}


def run(registry_path: Path, environment_path: Path, output_dir: Path, run_id: str) -> None:
    require_columns(registry_path, REGISTRY_COLS)
    require_columns(environment_path, ENV_COLS)
    registry = pd.read_csv(registry_path, usecols=REGISTRY_COLS)
    environment = pd.read_csv(environment_path, usecols=ENV_COLS)
    if registry.cohort_id.duplicated().any() or environment.cohort_id.duplicated().any():
        raise ValueError("cohort-level inputs must have one row per cohort")
    if set(registry.cohort_id) != set(environment.cohort_id):
        raise ValueError("registry and environment metadata cohort sets differ")

    df = registry.merge(environment, on="cohort_id", how="inner", suffixes=("_registry", "_environment"), validate="one_to_one")
    same_cancer = df["cancer_group_registry"].fillna("") == df["cancer_group_environment"].fillna("")
    same_label_type = df["response_label_type_registry"].fillna("") == df["response_label_type_environment"].fillna("")
    same_quality = df["response_label_quality_registry"].fillna("") == df["response_label_quality_environment"].fillna("")

    known_endpoint = ~df["response_label_type_environment"].fillna("unknown").str.lower().eq("unknown")
    known_timepoint = ~df["timepoint_schema"].fillna("unknown").str.lower().eq("unknown")
    high_quality = df["response_label_quality_environment"].fillna("").str.lower().eq("high")
    unknown_context = df["treatment_context"].fillna("unknown").str.lower().eq("unknown")
    confounded_context = df["treatment_context"].fillna("").str.lower().eq("high_confounding_support")
    excluded_role = df["environment_role"].fillna("").str.lower().eq("stratification_or_excluded")
    unknown_cancer = df["cancer_group_environment"].fillna("unknown").str.lower().eq("unknown")

    audit = pd.DataFrame(
        {
            "cohort_id": df.cohort_id,
            "disease": df.disease,
            "cancer_group": df.cancer_group_environment,
            "treatment_context": df.treatment_context,
            "environment_role": df.environment_role,
            "response_label_type_metadata": df.response_label_type_environment,
            "response_label_quality_metadata": df.response_label_quality_environment,
            "timepoint_schema_metadata": df.timepoint_schema,
            "paired_available": df.paired_available_environment,
            "endpoint_metadata_consistent": known_endpoint & same_label_type,
            "label_quality_metadata_consistent": same_quality,
            "cohort_metadata_consistent": same_cancer,
            "metadata_only": True,
            "response_values_read": False,
            "numeric_contrast_run": False,
        }
    )
    audit["design_warning"] = ""
    audit.loc[~known_endpoint, "design_warning"] += "endpoint_schema_unknown;"
    audit.loc[~known_timepoint, "design_warning"] += "timepoint_schema_unknown;"
    audit.loc[~high_quality, "design_warning"] += "label_quality_below_high;"
    audit.loc[confounded_context, "design_warning"] += "context_marked_high_confounding;"
    audit.loc[unknown_cancer, "design_warning"] += "cancer_group_unknown;"
    audit.loc[excluded_role, "design_warning"] += "stratification_or_excluded_role;"
    audit.loc[unknown_context, "design_warning"] += "treatment_context_unknown;"
    audit["design_warning"] = audit.design_warning.str.rstrip(";")

    context_cancer = pd.crosstab(df.treatment_context, df.cancer_group_environment)
    context_summary = []
    for context, row in context_cancer.iterrows():
        present = [str(c) for c, n in row.items() if n > 0]
        context_summary.append(
            {
                "treatment_context": str(context),
                "cancer_groups_present": ";".join(present),
                "n_cancer_groups": len(present),
                "n_cohorts": int(row.sum()),
                "shared_context_metadata_supported": len(present) >= 2 and "unknown" not in present,
            }
        )
    context_summary_df = pd.DataFrame(context_summary)

    known_contexts = int((~unknown_context).sum())
    summary = {
        "schema_version": "v1",
        "record_type": "wp6a_anchor_metadata_summary",
        "status": "COMPLETE_METADATA_ONLY",
        "run_id": run_id,
        "response_lock": "LOCKED",
        "response_values_read": False,
        "outcome_values_read": False,
        "numeric_contrast_run": False,
        "n_cohorts": int(len(df)),
        "n_independent_cohorts": int(df.cohort_id.nunique()),
        "duplicate_cohort_rows": int(len(df) - df.cohort_id.nunique()),
        "known_treatment_context_cohorts": known_contexts,
        "known_endpoint_metadata_cohorts": int(known_endpoint.sum()),
        "high_quality_label_metadata_cohorts": int(high_quality.sum()),
        "known_timepoint_schema_cohorts": int(known_timepoint.sum()),
        "label_quality_counts": value_counts(df, "response_label_quality_environment"),
        "endpoint_schema_counts": value_counts(df, "response_label_type_environment"),
        "treatment_context_counts": value_counts(df, "treatment_context"),
        "environment_role_counts": value_counts(df, "environment_role"),
        "timepoint_schema_counts": value_counts(df, "timepoint_schema"),
        "context_cancer_table": {
            str(k): {str(c): int(v) for c, v in row.items() if int(v) > 0}
            for k, row in context_cancer.iterrows()
        },
        "metadata_findings": [
            "Cohort rows are independent at the registry level; no duplicate cohort_id was found.",
            "Treatment contexts and cancer groups are not interchangeable metadata environments; context-stratified design is required.",
            "Shared-direction transport is not identified under response lock because this stage contains no numeric response contrast.",
            "Unknown endpoints/timepoints, low-quality labels, high-confounding contexts, and stratification/excluded roles require support-only or exclusion handling.",
        ],
        "proposed_estimand_id": "anchor_metadata_stratified_direction_v0",
        "estimand_status": "NOT_ESTIMABLE_UNDER_RESPONSE_LOCK",
        "allowed_claim_layer_now": "metadata_design_only",
        "wp6b_required": True,
        "wp6b_unlock": "separate_response_unlock_required",
    }
    manifest = {
        "schema_version": "v1",
        "record_type": "wp6a_anchor_metadata_manifest",
        "status": "COMPLETE_METADATA_ONLY",
        "run_id": run_id,
        "response_lock": "LOCKED",
        "response_values_read": False,
        "outcome_values_read": False,
        "numeric_contrast_run": False,
        "input_files": [
            {"path": str(registry_path), "sha256": sha256(registry_path), "columns_read": REGISTRY_COLS},
            {"path": str(environment_path), "sha256": sha256(environment_path), "columns_read": ENV_COLS},
        ],
        "forbidden_inputs_not_read": ["patient_metadata_master.csv", "response_raw", "response_binary_harmonized", "outcome_value", "failure_value"],
        "outputs": ["anchor_metadata_audit.csv", "anchor_context_summary.csv", "anchor_metadata_summary.yaml"],
        "next_gate": "WP6b requires separate response unlock and an explicit estimand decision",
    }

    output_dir.mkdir(parents=True, exist_ok=True)
    audit.to_csv(output_dir / "anchor_metadata_audit.csv", index=False)
    context_summary_df.to_csv(output_dir / "anchor_context_summary.csv", index=False)
    (output_dir / "anchor_metadata_summary.yaml").write_text(yaml.safe_dump(summary, sort_keys=False, allow_unicode=True))
    (output_dir / "anchor_metadata_manifest.yaml").write_text(yaml.safe_dump(manifest, sort_keys=False, allow_unicode=True))


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--registry", type=Path, required=True)
    ap.add_argument("--environment", type=Path, required=True)
    ap.add_argument("--output-dir", type=Path, required=True)
    ap.add_argument("--run-id", required=True)
    args = ap.parse_args()
    run(args.registry, args.environment, args.output_dir, args.run_id)


if __name__ == "__main__":
    main()
