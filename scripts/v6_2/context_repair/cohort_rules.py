"""Cohort-specific, reviewable rules for the context-repair v1 metadata freeze.

This module deliberately contains no file I/O.  Keeping the exceptional rules
here makes the metadata builder deterministic and prevents downstream code from
silently reinterpreting response labels or cohort roles.
"""

from __future__ import annotations

from typing import Mapping


FROZEN_VERSION = "frozen_v1"

# These corrections are intentionally narrow.  They are applied to the cohort,
# sample, patient, feature-eligibility, and split views by the builder.
COHORT_CORRECTIONS: dict[str, dict[str, str]] = {
    "GSE200996": {
        "modality": "scRNA",
        "dataset_role": "conditional_scRNA_support",
        "response_endpoint_type": "ordinal_clinical_benefit",
        "feature_family": "scRNA_expression",
        "revision_reason": "v1: local 220141-cell scRNA object; retain High/Medium/Low as ordinal response",
    },
    "GSE179994": {
        "modality": "scRNA",
        "dataset_role": "support_only",
        "feature_family": "scRNA_expression",
        "revision_reason": "v1: local scRNA object restored from erroneous bulk registration; response remains unavailable",
    },
    "mendeley_skrx2fz79n": {
        "modality": "scRNA",
        "dataset_role": "localization_support_only",
        "feature_family": "scRNA_expression",
        "revision_reason": "v1: local h5ad is scRNA, not spatial; retain only localization/support role because response is unavailable",
    },
    "bi_2021_rcc": {
        "modality": "scRNA",
        "dataset_role": "sensitivity_only",
        "feature_family": "scRNA_expression",
        "revision_reason": "v1: scRNA RCC sensitivity cohort; not a primary training cohort",
    },
}

STRICT_GSE272993_ARM = "apd1"
STRICT_GSE272993_TIMEPOINT = "baseline"
STRICT_GSE272993_SUBSET = "aPD1_baseline_only"


def correction_for(cohort_id: str) -> dict[str, str]:
    """Return a copy so callers cannot mutate the frozen rule registry."""

    return dict(COHORT_CORRECTIONS.get(cohort_id, {}))


def response_binary(raw_response: str, cohort_id: str) -> tuple[str, str, str]:
    """Return (binary_label, confidence, mapping_rule) without inventing labels.

    GSE200996 is explicitly ordinal.  GSE266919 ``NE`` is explicitly unknown;
    it must never fall through to a non-responder mapping.
    """

    raw = (raw_response or "").strip()
    normalized = raw.upper()
    if cohort_id == "GSE200996":
        return "unknown", "high", "preserve_ordinal_response_no_binary_mapping"
    if normalized in {"CR", "PR", "OR", "RESPONDER", "R"}:
        return "responder", "high", "recist_cr_pr_or_responder"
    if normalized in {"PD", "SD", "NR", "NON_RESPONDER", "NON-RESPONDER"}:
        return "non_responder", "high", "recist_pd_sd_or_non_responder"
    if normalized in {"", "UNKNOWN", "NA", "N/A", "NE", "NOT_EVALUABLE"}:
        return "unknown", "high" if normalized in {"NE", "NOT_EVALUABLE"} else "low", "preserve_unknown_or_ne"
    return "unknown", "low", "unrecognized_raw_response_preserved"


def gse272993_strict_subset(row: Mapping[str, str]) -> bool:
    """Accept only the pre-treatment anti-PD1 arm for strict GSE272993 use."""

    if row.get("cohort_id") != "GSE272993":
        return False
    arm = (row.get("treatment_raw") or "").strip().casefold()
    timepoint = (row.get("timepoint_raw") or "").strip().casefold()
    return arm == STRICT_GSE272993_ARM and timepoint == STRICT_GSE272993_TIMEPOINT


def default_study_subject_key(row: Mapping[str, str]) -> str:
    """Supply deterministic keys when the sidecar has no row for a sample."""

    cohort_id = row.get("cohort_id", "")
    patient_id = row.get("patient_id", "")
    if cohort_id in {"GSE123813_bcc", "GSE123813_scc"} and patient_id == "su010":
        return "GSE123813::su010"
    return row.get("patient_key") or f"{cohort_id}::{patient_id}"
