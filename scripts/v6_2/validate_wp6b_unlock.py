#!/usr/bin/env python3
"""Fail-closed validator for the one-time WP6b response unlock contract."""

from __future__ import annotations

import argparse
import hashlib
import sys
from datetime import datetime
from pathlib import Path

import yaml


EXPECTED_READS = {
    "response_binary_harmonized",
    "response_harmonization_confidence",
    "response_endpoint_type",
    "treatment_arm",
    "timepoint_use_boundary",
    "supervised_use_allowed",
    "support_use_allowed",
}
EXPECTED_ALLOWED_OPERATIONS = {
    "canonical_anchor_join",
    "pre_registered_class_coverage_audit",
    "environment_specific_joint_contrast",
    "patient_bootstrap",
    "within_environment_permutation",
}
EXPECTED_FORBIDDEN_OPERATIONS = {
    "response_driven_module_selection",
    "response_driven_environment_definition",
    "endpoint_pooling_outside_contract",
    "response_driven_split_selection",
    "shared_direction_or_barrier_claim_before_G2",
}
EXPECTED_SEPARATE_FROM = {"WP7", "Phase8", "repair", "spatial", "perturbation", "X_class"}


def file_sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def validate(path: Path) -> dict:
    data = yaml.safe_load(path.read_text())
    errors: list[str] = []
    if data.get("contract_id") != "wp6b_response_unlock_contract_v1":
        errors.append("wrong_contract_id")
    approval = data.get("approval", {})
    if data.get("status") != "APPROVED_FOR_WP6B":
        errors.append("status_not_APPROVED_FOR_WP6B")
    if approval.get("required") is not True:
        errors.append("approval_required_must_be_true")
    if approval.get("approver_role") != "mentor":
        errors.append("approver_role_must_be_mentor")
    if not approval.get("approved_by"):
        errors.append("approved_by_missing")
    if not approval.get("approved_at_utc"):
        errors.append("approved_at_utc_missing")
    if approval.get("approval_basis") != "explicit_mentor_separate_unlock":
        errors.append("approval_basis_missing_or_wrong")
    approved_at = approval.get("approved_at_utc")
    if approved_at:
        try:
            datetime.fromisoformat(str(approved_at).replace("Z", "+00:00"))
        except ValueError:
            errors.append("approved_at_utc_not_iso8601")
    lock = data.get("lock", {})
    if lock.get("current_status") != "UNLOCKED_FOR_WP6B":
        errors.append("current_status_must_be_UNLOCKED_FOR_WP6B")
    if lock.get("unlock_scope") != "WP6b_numeric_anchor_audit_only":
        errors.append("unlock_scope_must_be_WP6B_only")
    if lock.get("one_time_run") is not True:
        errors.append("one_time_run_must_be_true")
    if lock.get("representation_retraining") != "forbidden":
        errors.append("representation_retraining_must_remain_forbidden")
    if lock.get("response_model_retraining") != "forbidden":
        errors.append("response_model_retraining_must_remain_forbidden")
    if set(lock.get("separate_from", [])) != EXPECTED_SEPARATE_FROM:
        errors.append("separate_from_scope_changed")
    if data.get("required_lock_state_after_approval") != "UNLOCKED_FOR_WP6B":
        errors.append("required_lock_state_after_approval_missing")
    if data.get("allowed_reads_exact") is not True or set(data.get("allowed_reads_after_approval", [])) != EXPECTED_READS:
        errors.append("allowed_reads_whitelist_changed")
    if data.get("forbidden_operations_exact") is not True or set(data.get("forbidden_operations", [])) != EXPECTED_FORBIDDEN_OPERATIONS:
        errors.append("forbidden_operations_whitelist_changed")
    consumption = data.get("consumption", {})
    if consumption.get("status") != "NOT_CONSUMED" or consumption.get("consumed_at_utc") is not None or consumption.get("consumed_run_id") is not None:
        errors.append("unlock_already_consumed_or_consumption_state_invalid")

    lineage = data.get("lineage", {})
    for key in ("parent_snapshot_path", "g0_freeze_manifest_path"):
        p = Path(lineage.get(key, ""))
        if not p.exists():
            errors.append(f"lineage_path_missing:{key}")
    if lineage.get("parent_snapshot_path") and Path(lineage["parent_snapshot_path"]).exists():
        if file_sha256(Path(lineage["parent_snapshot_path"])) != lineage.get("parent_snapshot_sha256"):
            errors.append("parent_snapshot_hash_mismatch")
    if lineage.get("g0_freeze_manifest_path") and Path(lineage["g0_freeze_manifest_path"]).exists():
        if file_sha256(Path(lineage["g0_freeze_manifest_path"])) != lineage.get("g0_freeze_manifest_sha256"):
            errors.append("g0_freeze_manifest_hash_mismatch")

    g1 = data.get("g1_estimand_approval", {})
    if g1.get("required") is not True or g1.get("atomic_with_wp6b_unlock") is not True:
        errors.append("g1_estimand_approval_must_be_atomic")
    g1_path = Path(g1.get("contract_path", ""))
    if not g1_path.exists():
        errors.append("g1_estimand_contract_missing")
    elif file_sha256(g1_path) != g1.get("contract_sha256"):
        errors.append("g1_estimand_contract_hash_mismatch")
    else:
        g1_data = yaml.safe_load(g1_path.read_text())
        if g1_data.get("status") != "FROZEN_APPROVED_FOR_WP6B":
            errors.append("g1_estimand_status_not_FROZEN_APPROVED_FOR_WP6B")
        if g1_data.get("contract_id") != "g1_shared_direction_estimand_contract_v1":
            errors.append("g1_estimand_contract_id_mismatch")
        if g1_data.get("model", {}).get("WP6b_shared_direction_or_context_verdict") != "forbidden_before_G2_WP7":
            errors.append("g1_wp6b_wp7_boundary_changed")
        if g1_data.get("endpoint_policy", {}).get("cross_endpoint_raw_pooling") != "forbidden":
            errors.append("cross_endpoint_pooling_policy_changed")
    result = {"valid": not errors, "status": data.get("status"), "errors": errors}
    if errors:
        raise RuntimeError(yaml.safe_dump(result, sort_keys=False))
    return result


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--contract", type=Path, required=True)
    args = ap.parse_args()
    try:
        print(yaml.safe_dump(validate(args.contract), sort_keys=False))
    except RuntimeError as exc:
        print(str(exc), file=sys.stderr)
        raise SystemExit(3)


if __name__ == "__main__":
    main()
