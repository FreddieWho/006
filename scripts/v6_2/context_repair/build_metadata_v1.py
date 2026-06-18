#!/usr/bin/env python3
"""Build context-repair v1 metadata without modifying the frozen_v0 inputs.

The identity sidecar is a CSV/TSV with ``study_subject_key`` and either
``sample_key`` or the pair ``cohort_id,sample_id``.  It is intentionally used
only to resolve identity: all descriptive metadata originates in frozen_v0.
"""

from __future__ import annotations

import argparse
import csv
import sys
from collections import defaultdict
from pathlib import Path
from typing import Iterable

try:  # Supports both direct execution and package-style test imports.
    from .cohort_rules import (
        COHORT_CORRECTIONS,
        FROZEN_VERSION,
        STRICT_GSE272993_SUBSET,
        correction_for,
        default_study_subject_key,
        gse272993_strict_subset,
        response_binary,
    )
except ImportError:  # pragma: no cover - exercised by the CLI invocation.
    from cohort_rules import (
        COHORT_CORRECTIONS,
        FROZEN_VERSION,
        STRICT_GSE272993_SUBSET,
        correction_for,
        default_study_subject_key,
        gse272993_strict_subset,
        response_binary,
    )


TABLES = {
    "cohort": ("cohort_registry.frozen_v0.csv", "cohort_id"),
    "sample": ("sample_metadata_master.frozen_v0.csv", "sample_key"),
    "patient": ("patient_metadata_master.frozen_v0.csv", "patient_key"),
    "response": ("response_label_environment.frozen_v0.csv", "response_environment_id"),
    "dataset_role": ("dataset_role_and_feature_eligibility.frozen_v0.csv", "cohort_id"),
    "patient_split": ("patient_split.frozen_v0.csv", "patient_key"),
}

EXTRA_COLUMNS = {
    "cohort": ["metadata_revision_reason"],
    "sample": ["study_subject_key", "aggregation_key", "strict_subset_eligible", "metadata_revision_reason"],
    "patient": ["study_subject_key", "aggregation_key", "metadata_revision_reason"],
    "response": ["response_scale", "metadata_revision_reason"],
    "dataset_role": ["metadata_revision_reason"],
    "patient_split": ["study_subject_key", "aggregation_key", "metadata_revision_reason"],
}


def read_csv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        if not reader.fieldnames:
            raise ValueError(f"empty CSV: {path}")
        return list(reader.fieldnames), [dict(row) for row in reader]


def write_csv(path: Path, fieldnames: list[str], rows: Iterable[dict[str, str]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def append_reason(existing: str, reason: str) -> str:
    parts = [part for part in (existing or "").split("; ") if part]
    if reason and reason not in parts:
        parts.append(reason)
    return "; ".join(parts)


def read_identity_sidecar(path: Path) -> tuple[dict[str, str], dict[tuple[str, str], str]]:
    if not path.is_file():
        raise FileNotFoundError(f"identity sidecar not found: {path}")
    delimiter = "\t" if path.suffix.lower() == ".tsv" else ","
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter=delimiter)
        fields = set(reader.fieldnames or [])
        if "study_subject_key" not in fields:
            raise ValueError("identity sidecar requires study_subject_key")
        if "sample_key" not in fields and not {"cohort_id", "sample_id"}.issubset(fields):
            raise ValueError("identity sidecar requires sample_key or cohort_id and sample_id")
        by_sample_key: dict[str, str] = {}
        by_cohort_sample: dict[tuple[str, str], str] = {}
        for line_number, row in enumerate(reader, start=2):
            subject_key = (row.get("study_subject_key") or "").strip()
            if not subject_key:
                raise ValueError(f"identity sidecar line {line_number} has an empty study_subject_key")
            sample_key = (row.get("sample_key") or "").strip()
            cohort_id = (row.get("cohort_id") or "").strip()
            sample_id = (row.get("sample_id") or "").strip()
            if sample_key:
                _insert_identity(by_sample_key, sample_key, subject_key, line_number)
            if cohort_id and sample_id:
                _insert_identity(by_cohort_sample, (cohort_id, sample_id), subject_key, line_number)
            if not sample_key and not (cohort_id and sample_id):
                raise ValueError(f"identity sidecar line {line_number} has no usable sample identity")
    if not by_sample_key and not by_cohort_sample:
        raise ValueError("identity sidecar contains no usable rows")
    return by_sample_key, by_cohort_sample


def _insert_identity(mapping: dict, key: object, value: str, line_number: int) -> None:
    previous = mapping.setdefault(key, value)
    if previous != value:
        raise ValueError(f"identity sidecar line {line_number} conflicts for {key!r}: {previous!r} vs {value!r}")


def sidecar_subject_key(
    row: dict[str, str], by_sample_key: dict[str, str], by_cohort_sample: dict[tuple[str, str], str]
) -> str:
    return (
        by_sample_key.get(row.get("sample_key", ""))
        or by_cohort_sample.get((row.get("cohort_id", ""), row.get("sample_id", "")))
        or default_study_subject_key(row)
    )


def normalize_sample(
    row: dict[str, str], by_sample_key: dict[str, str], by_cohort_sample: dict[tuple[str, str], str]
) -> dict[str, str]:
    result = dict(row)
    cohort_id = result["cohort_id"]
    correction = correction_for(cohort_id)
    if correction:
        result["modality"] = correction["modality"]
        result["dataset_role"] = correction["dataset_role"]
        if cohort_id == "GSE200996":
            result["response_endpoint_type"] = correction["response_endpoint_type"]
        result["exclusion_or_downgrade_reason"] = append_reason(
            result.get("exclusion_or_downgrade_reason", ""), correction["revision_reason"]
        )
    binary, confidence, _ = response_binary(result.get("response_raw", ""), cohort_id)
    current_binary = result.get("response_binary_harmonized", "").strip().lower()
    if binary != "unknown" or current_binary in {"", "unknown", "na", "nan"}:
        result["response_binary_harmonized"] = binary
        result["response_harmonization_confidence"] = confidence
    if cohort_id == "TASK01":
        baseline = result.get("timepoint", "").strip().lower() in {"pre", "baseline"}
        pd1 = result.get("treatment_raw", "").strip().upper() == "PD1"
        result["supervised_use_allowed"] = "yes" if binary != "unknown" and baseline and pd1 else "no"
        result["support_use_allowed"] = "yes"
        result["metadata_revision_reason"] = (
            "v1: source mRECIST CR/PR->responder, SD/PD->non_responder; "
            "baseline PD1 samples supervised, later visits temporal sensitivity only"
        )
    result["study_subject_key"] = sidecar_subject_key(result, by_sample_key, by_cohort_sample)
    result["aggregation_key"] = result["study_subject_key"]
    result["strict_subset_eligible"] = "yes" if gse272993_strict_subset(result) else "no"
    if cohort_id == "GSE272993":
        result["metadata_revision_reason"] = (
            f"v1: strict analysis subset is {STRICT_GSE272993_SUBSET}; all other arms/timepoints remain non-strict"
        )
    elif cohort_id == "TASK01":
        pass
    elif correction:
        result["metadata_revision_reason"] = correction["revision_reason"]
    elif result.get("response_raw", "").strip().upper() == "NE":
        result["metadata_revision_reason"] = "v1: NE preserved as unknown; never mapped to non_responder"
    elif result["study_subject_key"] != default_study_subject_key(result):
        result["metadata_revision_reason"] = "v1: identity sidecar study_subject_key"
    else:
        result["metadata_revision_reason"] = ""
    result["input_version"] = FROZEN_VERSION
    return result


def normalize_cohort(row: dict[str, str]) -> dict[str, str]:
    result = dict(row)
    correction = correction_for(result["cohort_id"])
    if correction:
        result["modality"] = correction["modality"]
        result["dataset_role"] = correction["dataset_role"]
        if "response_endpoint_type" in correction:
            result["response_endpoint_type"] = correction["response_endpoint_type"]
        result["downgrade_reason"] = append_reason(result.get("downgrade_reason", ""), correction["revision_reason"])
        result["metadata_revision_reason"] = correction["revision_reason"]
    else:
        result["metadata_revision_reason"] = ""
    result["input_version"] = FROZEN_VERSION
    return result


def normalize_dataset_role(row: dict[str, str]) -> dict[str, str]:
    result = dict(row)
    correction = correction_for(result["cohort_id"])
    if correction:
        result["dataset_role"] = correction["dataset_role"]
        result["feature_family"] = correction["feature_family"]
        result["supervised_use_allowed"] = "no"
        result["support_use_allowed"] = "yes"
        if correction["modality"] == "spatial":
            result["spatial_lane"] = "yes"
            result["spatial_coordinates_available"] = "yes"
        if correction["modality"] == "bulkRNA":
            result["bulk_expression_available"] = "yes"
            result["external_validation_lane"] = "yes"
        result["downgrade_reason"] = append_reason(result.get("downgrade_reason", ""), correction["revision_reason"])
        result["metadata_revision_reason"] = correction["revision_reason"]
    else:
        result["metadata_revision_reason"] = ""
    return result


def normalize_patient(row: dict[str, str], samples_by_patient: dict[str, list[dict[str, str]]]) -> dict[str, str]:
    result = dict(row)
    cohort_id = result["cohort_id"]
    samples = samples_by_patient.get(result["patient_key"], [])
    correction = correction_for(cohort_id)
    if samples:
        subjects = {sample["study_subject_key"] for sample in samples}
        if len(subjects) != 1:
            raise ValueError(f"patient {result['patient_key']} resolves to multiple study_subject_keys")
        result["study_subject_key"] = subjects.pop()
        result["aggregation_key"] = result["study_subject_key"]
        result["n_samples"] = str(len(samples))
        result["available_modalities"] = ";".join(sorted({sample["modality"] for sample in samples}))
        timepoints = {sample["timepoint"] for sample in samples if sample.get("timepoint")}
        result["timepoint_schema"] = ";".join(sorted(timepoints)) or result.get("timepoint_schema", "")
        result["paired_pre_post_available"] = "yes" if {"pre", "post"}.issubset(timepoints) else "no"
        raw_values = sorted({sample["response_raw"] for sample in samples if sample.get("response_raw")})
        result["response_raw_summary"] = ";".join(raw_values) or "unknown"
        binaries = {sample["response_binary_harmonized"] for sample in samples}
        result["response_binary_harmonized"] = binaries.pop() if len(binaries) == 1 else "unknown"
    else:
        result["study_subject_key"] = result["patient_key"]
        result["aggregation_key"] = result["study_subject_key"]
    if correction:
        result["available_modalities"] = correction["modality"]
        result["supervised_use_allowed"] = "no"
        result["support_use_allowed"] = "yes"
        result["exclusion_or_downgrade_reason"] = append_reason(
            result.get("exclusion_or_downgrade_reason", ""), correction["revision_reason"]
        )
        result["metadata_revision_reason"] = correction["revision_reason"]
    elif cohort_id == "GSE272993":
        result["metadata_revision_reason"] = f"v1: strict analysis subset is {STRICT_GSE272993_SUBSET}"
    else:
        result["metadata_revision_reason"] = ""
    result["input_version"] = FROZEN_VERSION
    return result


def normalize_split(row: dict[str, str], patients_by_key: dict[str, dict[str, str]]) -> dict[str, str]:
    result = dict(row)
    patient = patients_by_key.get(result["patient_key"])
    if patient:
        result["study_subject_key"] = patient["study_subject_key"]
        result["aggregation_key"] = patient["aggregation_key"]
    else:
        result["study_subject_key"] = result["patient_key"]
        result["aggregation_key"] = result["patient_key"]
    correction = correction_for(result["cohort_id"])
    if correction:
        result["split"] = "support_only"
        result["split_source"] = "context_repair_v1_role_based_exclusion"
        result["split_eligibility"] = "no"
        result["reason_if_excluded"] = correction["revision_reason"]
        result["metadata_revision_reason"] = correction["revision_reason"]
    elif result["cohort_id"] == "GSE272993":
        result["metadata_revision_reason"] = f"v1: strict analysis subset is {STRICT_GSE272993_SUBSET}"
    else:
        result["metadata_revision_reason"] = ""
    result["input_version"] = FROZEN_VERSION
    return result


def normalize_response(rows: list[dict[str, str]], samples: list[dict[str, str]]) -> list[dict[str, str]]:
    samples_by_cohort: dict[str, list[dict[str, str]]] = defaultdict(list)
    for sample in samples:
        samples_by_cohort[sample["cohort_id"]].append(sample)
    normalized: list[dict[str, str]] = []
    for row in rows:
        result = dict(row)
        cohort_samples = samples_by_cohort.get(result["cohort_id"], [])
        if cohort_samples:
            labels_by_patient: dict[str, set[str]] = defaultdict(set)
            sample_counts = defaultdict(int)
            for sample in cohort_samples:
                label, _, mapping_rule = response_binary(sample.get("response_raw", ""), sample["cohort_id"])
                labels_by_patient[sample["patient_key"]].add(label)
                sample_counts[label] += 1
            patient_counts = defaultdict(int)
            for labels in labels_by_patient.values():
                patient_counts[labels.pop() if len(labels) == 1 else "unknown"] += 1
            result["n_patients_R"] = str(patient_counts["responder"])
            result["n_patients_NR"] = str(patient_counts["non_responder"])
            result["n_patients_unknown"] = str(patient_counts["unknown"])
            result["n_samples_R"] = str(sample_counts["responder"])
            result["n_samples_NR"] = str(sample_counts["non_responder"])
            result["n_samples_unknown"] = str(sample_counts["unknown"])
            if result["cohort_id"] == "GSE200996":
                result["endpoint_type"] = "ordinal_clinical_benefit"
                result["response_mapping_rule"] = "preserve_ordinal_response_no_binary_mapping"
                result["response_scale"] = "ordinal:Low<Medium<High"
                result["metadata_revision_reason"] = correction_for("GSE200996")["revision_reason"]
            else:
                result["response_scale"] = "binary_harmonized_with_unknown"
                result["metadata_revision_reason"] = (
                    "v1: NE preserved as unknown; never mapped to non_responder"
                    if result["cohort_id"] == "GSE266919"
                    else ""
                )
        else:
            result["response_scale"] = "registry_only"
            result["metadata_revision_reason"] = correction_for(result["cohort_id"]).get("revision_reason", "")
        normalized.append(result)
    return normalized


def is_changed(table: str, original: dict[str, str], updated: dict[str, str]) -> bool:
    relevant = set(original)
    if table in {"sample", "patient", "patient_split"}:
        relevant.add("study_subject_key")
    if table == "sample":
        relevant.add("strict_subset_eligible")
    if table == "response":
        relevant.add("response_scale")
    return any(original.get(column, "") != updated.get(column, "") for column in relevant)


def aggregation_spec() -> list[dict[str, str]]:
    return [
        {
            "table_name": "sample_metadata_master.effective_v1",
            "entity_level": "sample",
            "source_key": "sample_key",
            "aggregation_key": "study_subject_key",
            "uniqueness_rule": "sample_key is unique; multiple samples may share a study_subject_key",
            "use_boundary": "aggregate repeated samples only after the analysis declares its sample/timepoint rule",
        },
        {
            "table_name": "patient_metadata_master.effective_v1",
            "entity_level": "registry_patient",
            "source_key": "patient_key",
            "aggregation_key": "study_subject_key",
            "uniqueness_rule": "patient_key is unique; cross-subcohort aliases may share a study_subject_key",
            "use_boundary": "split and leakage control use study_subject_key",
        },
        {
            "table_name": "GSE272993",
            "entity_level": "strict_analysis_sample",
            "source_key": "sample_key",
            "aggregation_key": "study_subject_key",
            "uniqueness_rule": "strict_subset_eligible=yes requires treatment_raw=aPD1 and timepoint_raw=Baseline",
            "use_boundary": f"strict subset: {STRICT_GSE272993_SUBSET}",
        },
        {
            "table_name": "GSE200996",
            "entity_level": "response",
            "source_key": "patient_key",
            "aggregation_key": "study_subject_key",
            "uniqueness_rule": "High/Medium/Low remain ordinal values",
            "use_boundary": "no binary response aggregation",
        },
    ]


def build(frozen_root: Path, identity_sidecar: Path, output_dir: Path) -> None:
    raw: dict[str, tuple[list[str], list[dict[str, str]]]] = {}
    for name, (filename, _) in TABLES.items():
        raw[name] = read_csv(frozen_root / filename)
    by_sample_key, by_cohort_sample = read_identity_sidecar(identity_sidecar)

    sample_rows = [normalize_sample(row, by_sample_key, by_cohort_sample) for row in raw["sample"][1]]
    samples_by_patient: dict[str, list[dict[str, str]]] = defaultdict(list)
    for sample in sample_rows:
        samples_by_patient[sample["patient_key"]].append(sample)
    effective = {
        "cohort": [normalize_cohort(row) for row in raw["cohort"][1]],
        "sample": sample_rows,
        "patient": [normalize_patient(row, samples_by_patient) for row in raw["patient"][1]],
        "dataset_role": [normalize_dataset_role(row) for row in raw["dataset_role"][1]],
    }
    patients_by_key = {row["patient_key"]: row for row in effective["patient"]}
    effective["patient_split"] = [normalize_split(row, patients_by_key) for row in raw["patient_split"][1]]
    effective["response"] = normalize_response(raw["response"][1], sample_rows)

    output_dir.mkdir(parents=True, exist_ok=True)
    for name, (_, key_column) in TABLES.items():
        original_rows = {row[key_column]: row for row in raw[name][1]}
        fields = raw[name][0] + [column for column in EXTRA_COLUMNS[name] if column not in raw[name][0]]
        effective_rows = effective[name]
        addendum_rows = [row for row in effective_rows if is_changed(name, original_rows[row[key_column]], row)]
        source_stem = TABLES[name][0].removesuffix(".frozen_v0.csv")
        write_csv(output_dir / f"{source_stem}.addendum_v1.csv", fields, addendum_rows)
        write_csv(output_dir / f"{source_stem}.effective_v1.csv", fields, effective_rows)
    spec_fields = ["table_name", "entity_level", "source_key", "aggregation_key", "uniqueness_rule", "use_boundary"]
    write_csv(output_dir / "aggregation_key_spec_v1.csv", spec_fields, aggregation_spec())


def parse_args(argv: list[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--frozen-root",
        type=Path,
        default=Path("results/v6_2"),
        help="directory containing the six frozen_v0 CSV registries",
    )
    parser.add_argument("--identity-sidecar", type=Path, required=True, help="CSV/TSV identity mapping sidecar")
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("results/v6_2/context_repair_v1"),
        help="new output directory; frozen_v0 is never modified",
    )
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv or sys.argv[1:])
    try:
        build(args.frozen_root, args.identity_sidecar, args.output_dir)
    except (FileNotFoundError, ValueError) as error:
        print(f"context_repair_v1: {error}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
