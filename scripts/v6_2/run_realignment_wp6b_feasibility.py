#!/usr/bin/env python3
"""WP6b-0 numeric anchor feasibility audit.

The unlock validator is called before any response-bearing column is read.
This stage audits coverage and design structure only; it never fits an effect
model and cannot produce a shared-direction or barrier verdict.
"""

from __future__ import annotations

import argparse
import hashlib
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

from validate_wp6b_unlock import validate as validate_unlock


MODULES = [f"FM{i:02d}" for i in range(1, 9)]
WIDE_META = ["cohort_id", "patient_key", "timepoint", "tissue_context", "patient_timepoint_context_id"]
PATIENT_COLS = [
    "patient_key",
    "cancer_type",
    "response_binary_harmonized",
    "response_harmonization_confidence",
    "supervised_use_allowed",
    "support_use_allowed",
    "timepoint_schema",
    "treatment_context_summary",
    "response_endpoint_type",
]
SAMPLE_COLS = [
    "patient_key",
    "cohort_id",
    "timepoint",
    "treatment_arm",
    "treatment_context",
    "response_endpoint_type",
    "response_binary_harmonized",
    "response_harmonization_confidence",
    "supervised_use_allowed",
    "support_use_allowed",
    "timepoint_use_boundary",
]
REGISTRY_COLS = ["cohort_id", "cancer_type", "treatment_context"]
ENVIRONMENT_COLS = [
    "cohort_id",
    "endpoint_type",
    "treatment_context",
    "clean_anchor_eligible",
    "supervised_eligible",
    "support_only",
    "confidence",
    "use_boundary",
]
FAILURE_VALUES = {"non_responder", "failure", "nr", "non-responder"}
RESPONSE_VALUES = {"responder", "response", "r"}


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def as_bool(value: object) -> bool:
    return str(value).strip().lower() in {"true", "1", "yes", "y"}


def endpoint_family(value: object) -> str:
    text = str(value).strip().lower()
    if not text or text in {"unknown", "nan", "none"}:
        return "unknown"
    return text.replace(" ", "_").replace("/", "_")


def deterministic_value(series: pd.Series) -> str:
    values = sorted({str(x) for x in series.dropna() if str(x).strip() not in {"", "nan", "None"}})
    return values[0] if len(values) == 1 else "__CONFLICT__" if values else "unknown"


def summarize_sample_metadata(sample: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for (patient_key, cohort_id), group in sample.groupby(["patient_key", "cohort_id"], dropna=False):
        row = {"patient_key": patient_key, "cohort_id": cohort_id}
        for col in SAMPLE_COLS:
            if col not in {"patient_key", "cohort_id", "timepoint"}:
                row[col] = deterministic_value(group[col])
        rows.append(row)
    return pd.DataFrame(rows)


def run(
    unlock_contract: Path,
    wide_matrix: Path,
    patient_metadata: Path,
    sample_metadata: Path,
    cohort_registry: Path,
    environment_table: Path,
    output_dir: Path,
    run_id: str,
) -> None:
    existing_manifest = output_dir / "wp6b_numeric_feasibility_manifest.yaml"
    if existing_manifest.exists():
        existing = yaml.safe_load(existing_manifest.read_text()) or {}
        if existing.get("run_id") == run_id and existing.get("response_values_read") is True:
            raise RuntimeError("existing_consumed_wp6b_run_output")
    # This is intentionally before any response-bearing usecols/read.
    unlock_result = validate_unlock(unlock_contract)
    if not unlock_result["valid"]:
        raise RuntimeError("WP6b unlock validation failed")

    wide = pd.read_parquet(wide_matrix, columns=WIDE_META + MODULES)
    patient = pd.read_csv(patient_metadata, usecols=PATIENT_COLS)
    sample = pd.read_csv(sample_metadata, usecols=SAMPLE_COLS)
    registry = pd.read_csv(cohort_registry, usecols=REGISTRY_COLS)
    environment = pd.read_csv(environment_table, usecols=ENVIRONMENT_COLS)

    baseline_timepoints = {"baseline", "pre", "t0", "b1"}
    wide["timepoint_norm"] = wide["timepoint"].astype(str).str.strip().str.lower()
    wide = wide[wide["timepoint_norm"].isin(baseline_timepoints)].copy()
    if wide.empty:
        raise RuntimeError("no baseline module rows available")

    # Patient-level baseline aggregation is deterministic and response-blind.
    agg = wide.groupby(["cohort_id", "patient_key"], as_index=False)[MODULES].mean()
    observed = wide.groupby(["cohort_id", "patient_key"], as_index=False)[MODULES].apply(lambda x: x.notna().any())
    observed = observed.rename(columns={m: f"{m}_observed" for m in MODULES})
    agg = agg.merge(observed, on=["cohort_id", "patient_key"], validate="one_to_one")
    agg["observed_fraction"] = agg[[f"{m}_observed" for m in MODULES]].mean(axis=1)

    sample["timepoint_norm"] = sample["timepoint"].astype(str).str.strip().str.lower()
    sample = sample[sample["timepoint_norm"].isin(baseline_timepoints)].copy()
    sample_summary = summarize_sample_metadata(sample)
    table = agg.merge(patient, on="patient_key", how="left", validate="many_to_one")
    table = table.merge(sample_summary, on=["patient_key", "cohort_id"], how="left", suffixes=("_patient", "_sample"), validate="many_to_one")
    table = table.merge(registry, on="cohort_id", how="left", suffixes=("", "_registry"), validate="many_to_one")
    table = table.merge(environment, on="cohort_id", how="left", suffixes=("", "_environment"), validate="many_to_one")

    table["cancer_type"] = table["cancer_type"].where(
        table["cancer_type"].notna() & ~table["cancer_type"].astype(str).str.lower().isin({"unknown", "nan"}),
        table["cancer_type_registry"],
    )
    table["cancer_group"] = table["cancer_type"].astype(str).str.lower()
    table["response_label"] = table["response_binary_harmonized_sample"].where(
        table["response_binary_harmonized_sample"].notna(), table["response_binary_harmonized_patient"]
    ).astype(str).str.strip().str.lower()
    table["known_response"] = table["response_label"].isin(RESPONSE_VALUES | FAILURE_VALUES)
    table["endpoint_family"] = table["response_endpoint_type_sample"].where(
        table["response_endpoint_type_sample"].notna(), table["response_endpoint_type_patient"]
    ).map(endpoint_family)
    table["treatment_axis"] = table["treatment_arm"].where(
        table["treatment_arm"].notna() & ~table["treatment_arm"].astype(str).str.lower().isin({"unknown", "nan", "__conflict__"}),
        table["treatment_context_registry"],
    ).astype(str).str.strip().str.lower()
    table["environment_signature"] = (
        table["cancer_type"].astype(str).str.lower()
        + "__" + table["treatment_axis"]
        + "__" + table["endpoint_family"]
        + "__baseline"
    )
    table["analysis_environment_id"] = table["environment_signature"] + "__" + table["cohort_id"].astype(str)
    table["metadata_conflict"] = table[[c for c in table.columns if c.endswith("_sample")]].eq("__CONFLICT__").any(axis=1)
    table["high_confidence"] = table["response_harmonization_confidence_sample"].where(
        table["response_harmonization_confidence_sample"].notna(), table["response_harmonization_confidence_patient"]
    ).astype(str).str.lower().eq("high")
    table["supervised_allowed"] = table["supervised_use_allowed_sample"].where(
        table["supervised_use_allowed_sample"].notna(), table["supervised_use_allowed_patient"]
    ).map(as_bool)
    table["clean_anchor"] = table["clean_anchor_eligible"].map(as_bool)
    table["primary_treatment_axis"] = table["treatment_context_registry"].astype(str).str.lower().eq("pd1_ici_anchor")
    table["primary_candidate"] = (
        table["known_response"]
        & table["high_confidence"]
        & table["supervised_allowed"]
        & table["clean_anchor"]
        & table["primary_treatment_axis"]
        & table["endpoint_family"].ne("unknown")
        & ~table["metadata_conflict"]
        & (table["observed_fraction"] >= 0.80)
    )

    rows = []
    for env_id, group in table.groupby("analysis_environment_id", dropna=False):
        # Counts must be computed on rows that pass the relevant row-level
        # quality gate.  Using the full environment after `.any()` could let
        # one eligible row promote an environment whose remaining rows are
        # unknown, low-confidence, or support-only.
        primary_group = group.loc[group["primary_candidate"]]
        primary = not primary_group.empty
        if primary:
            count_group = primary_group
        else:
            count_group = group.loc[
                group["known_response"]
                & group["endpoint_family"].ne("unknown")
                & ~group["metadata_conflict"]
            ]
        labels = count_group.loc[count_group["known_response"], "response_label"]
        n_r = int(labels.isin(RESPONSE_VALUES).sum())
        n_f = int(labels.isin(FAILURE_VALUES).sum())
        rows.append(
            {
                "analysis_environment_id": str(env_id),
                "environment_signature": str(group["environment_signature"].iloc[0]),
                "cohort_id": str(group["cohort_id"].iloc[0]),
                "cancer_type": str(group["cancer_type"].iloc[0]),
                "cancer_group": str(group["cancer_group"].iloc[0]),
                "treatment_axis": str(group["treatment_axis"].iloc[0]),
                "treatment_context": str(group["treatment_context_registry"].iloc[0]),
                "endpoint_family": str(group["endpoint_family"].iloc[0]),
                "n_baseline_patients": int(len(count_group)),
                "n_known_response": int(n_r + n_f),
                "n_responder": n_r,
                "n_failure": n_f,
                "min_class_n": int(min(n_r, n_f)),
                "median_observed_fraction": float(count_group["observed_fraction"].median()) if not count_group.empty else 0.0,
                "primary_candidate": primary,
                "metadata_conflict_n": int(group["metadata_conflict"].sum()),
            }
        )
    env = pd.DataFrame(rows).sort_values(["primary_candidate", "cancer_type", "cohort_id"], ascending=[False, True, True])
    env["primary_eligible"] = env["primary_candidate"] & (env["n_baseline_patients"] >= 15) & (env["min_class_n"] >= 5)
    env["validation_candidate"] = (
        ~env["primary_eligible"]
        & (env["n_baseline_patients"] >= 10)
        & (env["min_class_n"] >= 3)
        & env["endpoint_family"].ne("unknown")
    )

    primary = env[env["primary_eligible"]].copy()
    validation = env[env["validation_candidate"]].copy()
    cancer_groups = sorted(primary["cancer_group"].dropna().unique().tolist())
    endpoint_axes = primary.groupby(["endpoint_family", "treatment_axis"], dropna=False).size()
    compatible_axis_n = int((endpoint_axes >= 2).sum())
    go_wp7 = bool(
        len(primary) >= 2
        and len(cancer_groups) >= 2
        and (len(primary) >= 3 or (len(primary) >= 2 and len(validation) >= 1))
        and compatible_axis_n >= 1
    )
    stop_reasons = []
    if len(primary) < 2:
        stop_reasons.append("fewer_than_two_primary_eligible_environments")
    if len(cancer_groups) < 2:
        stop_reasons.append("fewer_than_two_cancer_groups")
    if len(primary) == 2 and len(validation) < 1:
        stop_reasons.append("no_independent_validation_environment")
    if compatible_axis_n < 1:
        stop_reasons.append("no_repeated_endpoint_and_treatment_axis")
    status = "WP6B_FEASIBILITY_PASS_WP7_PENDING_G2" if go_wp7 else "WP6B_STOPPED_FEASIBILITY"
    summary = {
        "schema_version": "v1",
        "record_type": "wp6b_numeric_feasibility_summary",
        "status": status,
        "run_id": run_id,
        "response_lock": "UNLOCKED_FOR_WP6B",
        "response_values_read": True,
        "outcome_values_read": False,
        "effect_model_fit": False,
        "shared_direction_fit": False,
        "n_baseline_rows_after_aggregation": int(len(table)),
        "n_analysis_environments": int(len(env)),
        "n_primary_candidate_environments": int(env["primary_candidate"].sum()),
        "n_primary_eligible_environments": int(len(primary)),
        "n_validation_candidate_environments": int(len(validation)),
        "primary_cancer_groups": cancer_groups,
        "repeated_endpoint_treatment_axes": compatible_axis_n,
        "go_wp7": go_wp7,
        "stop_reasons": stop_reasons,
        "pre_registered_thresholds": {"environment_n": 15, "minimum_class_n": 5, "observed_fraction": 0.80},
        "next_gate": "G2_WP7_CONTINUATION" if go_wp7 else "DATA_GAP_SPECIFICATION",
    }
    manifest = {
        "schema_version": "v1",
        "record_type": "wp6b_numeric_feasibility_manifest",
        "status": status,
        "run_id": run_id,
        "unlock_contract": {"path": str(unlock_contract), "sha256": sha256(unlock_contract)},
        "response_values_read": True,
        "outcome_values_read": False,
        "effect_model_fit": False,
        "input_files": [
            {"path": str(wide_matrix), "sha256": sha256(wide_matrix), "columns_read": WIDE_META + MODULES},
            {"path": str(patient_metadata), "sha256": sha256(patient_metadata), "columns_read": PATIENT_COLS},
            {"path": str(sample_metadata), "sha256": sha256(sample_metadata), "columns_read": SAMPLE_COLS},
            {"path": str(cohort_registry), "sha256": sha256(cohort_registry), "columns_read": REGISTRY_COLS},
            {"path": str(environment_table), "sha256": sha256(environment_table), "columns_read": ENVIRONMENT_COLS},
        ],
        "outputs": ["anchor_environment_feasibility.csv", "wp6b_feasibility_summary.yaml"],
        "forbidden_next_steps_before_g2": ["shared_direction_fit", "context_modulation_fit", "barrier_attribution", "phase8"],
    }
    output_dir.mkdir(parents=True, exist_ok=True)
    env.to_csv(output_dir / "anchor_environment_feasibility.csv", index=False)
    (output_dir / "wp6b_feasibility_summary.yaml").write_text(yaml.safe_dump(summary, sort_keys=False, allow_unicode=True))
    (output_dir / "wp6b_numeric_feasibility_manifest.yaml").write_text(yaml.safe_dump(manifest, sort_keys=False, allow_unicode=True))


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--unlock-contract", type=Path, required=True)
    ap.add_argument("--wide-matrix", type=Path, required=True)
    ap.add_argument("--patient-metadata", type=Path, required=True)
    ap.add_argument("--sample-metadata", type=Path, required=True)
    ap.add_argument("--cohort-registry", type=Path, required=True)
    ap.add_argument("--environment-table", type=Path, required=True)
    ap.add_argument("--output-dir", type=Path, required=True)
    ap.add_argument("--run-id", required=True)
    args = ap.parse_args()
    run(**vars(args))


if __name__ == "__main__":
    main()
