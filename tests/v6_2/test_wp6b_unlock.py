from __future__ import annotations

import importlib.util
from pathlib import Path

import pytest
import yaml


SCRIPT = Path(__file__).parents[2] / "scripts" / "v6_2" / "validate_wp6b_unlock.py"
spec = importlib.util.spec_from_file_location("validate_wp6b_unlock", SCRIPT)
unlock = importlib.util.module_from_spec(spec)
assert spec.loader is not None
spec.loader.exec_module(unlock)


def test_locked_contract_fails_closed(tmp_path):
    source = Path(__file__).parents[2] / "results" / "v6_2" / "scientific_engineering_realignment_v1" / "wp6b_response_unlock_contract_v1.yaml"
    data = yaml.safe_load(source.read_text())
    data["status"] = "PENDING_MENTOR_SEPARATE_UNLOCK"
    data["lock"]["current_status"] = "LOCKED"
    data["approval"] = {
        "required": True,
        "approver_role": "mentor",
        "approved_by": None,
        "approved_at_utc": None,
        "approval_basis": None,
    }
    path = tmp_path / "locked_unlock_contract.yaml"
    path.write_text(yaml.safe_dump(data, sort_keys=False))
    with pytest.raises(RuntimeError, match="status_not_APPROVED_FOR_WP6B"):
        unlock.validate(path)


def test_approved_contract_requires_explicit_separate_unlock(tmp_path):
    data = yaml.safe_load(
        (Path(__file__).parents[2] / "results" / "v6_2" / "scientific_engineering_realignment_v1" / "wp6b_response_unlock_contract_v1.yaml").read_text()
    )
    data["status"] = "APPROVED_FOR_WP6B"
    data["approval"] = {
        "required": True,
        "approver_role": "mentor",
        "approved_by": "mentor_user",
        "approved_at_utc": "2026-08-22T00:00:00Z",
        "approval_basis": "explicit_mentor_separate_unlock",
    }
    data["lock"]["current_status"] = "UNLOCKED_FOR_WP6B"
    data["consumption"] = {"status": "NOT_CONSUMED", "consumed_at_utc": None, "consumed_run_id": None}
    path = tmp_path / "approved.yaml"
    path.write_text(yaml.safe_dump(data, sort_keys=False))
    assert unlock.validate(path)["valid"] is True


@pytest.mark.parametrize(
    "mutation,expected",
    [
        (lambda d: d["lock"].update(current_status="LOCKED"), "current_status_must_be_UNLOCKED_FOR_WP6B"),
        (lambda d: d["lineage"].update(parent_snapshot_sha256="0" * 64), "parent_snapshot_hash_mismatch"),
        (lambda d: d["forbidden_operations"].clear(), "forbidden_operations_whitelist_changed"),
        (lambda d: d["approval"].update(approver_role="not_mentor"), "approver_role_must_be_mentor"),
        (lambda d: d["approval"].update(approved_at_utc="not-a-timestamp"), "approved_at_utc_not_iso8601"),
    ],
)
def test_approved_contract_rejects_tampering(tmp_path, mutation, expected):
    source = Path(__file__).parents[2] / "results" / "v6_2" / "scientific_engineering_realignment_v1" / "wp6b_response_unlock_contract_v1.yaml"
    data = yaml.safe_load(source.read_text())
    data["status"] = "APPROVED_FOR_WP6B"
    data["approval"] = {
        "required": True,
        "approver_role": "mentor",
        "approved_by": "mentor_user",
        "approved_at_utc": "2026-08-22T00:00:00Z",
        "approval_basis": "explicit_mentor_separate_unlock",
    }
    data["lock"]["current_status"] = "UNLOCKED_FOR_WP6B"
    data["consumption"] = {"status": "NOT_CONSUMED", "consumed_at_utc": None, "consumed_run_id": None}
    mutation(data)
    path = tmp_path / "tampered.yaml"
    path.write_text(yaml.safe_dump(data, sort_keys=False))
    with pytest.raises(RuntimeError, match=expected):
        unlock.validate(path)


def test_consumed_contract_cannot_be_reused(tmp_path):
    source = Path(__file__).parents[2] / "results" / "v6_2" / "scientific_engineering_realignment_v1" / "wp6b_response_unlock_contract_v1.yaml"
    data = yaml.safe_load(source.read_text())
    data["status"] = "APPROVED_FOR_WP6B"
    data["lock"]["current_status"] = "UNLOCKED_FOR_WP6B"
    data["approval"].update(
        approver_role="mentor",
        approved_by="mentor_user",
        approved_at_utc="2026-08-22T00:00:00Z",
        approval_basis="explicit_mentor_separate_unlock",
    )
    data["consumption"] = {
        "status": "CONSUMED",
        "consumed_at_utc": "2026-08-22T00:01:00Z",
        "consumed_run_id": "prior_run",
    }
    path = tmp_path / "consumed.yaml"
    path.write_text(yaml.safe_dump(data, sort_keys=False))
    with pytest.raises(RuntimeError, match="unlock_already_consumed_or_consumption_state_invalid"):
        unlock.validate(path)
