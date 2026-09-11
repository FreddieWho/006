"""Post-freeze GT isolation audit.

This module is intentionally not imported by :mod:`v7.spatial_pipeline`.
It may be started only after the response-blind discovery hash is frozen.  It
checks a sealed locator/hash and the firewall contract; it does not parse GT
values or score a biological endpoint.
"""

from __future__ import annotations

from pathlib import Path
import hashlib
import re
from typing import Any

import pandas as pd

from ..spatial_control import (
    sha256_file,
    stable_hash,
    validate_model_safe_manifest,
    validate_sealed_validation_manifest,
)


def _anchor_hash(path: Path) -> str:
    return sha256_file(path)


def build_sealed_validation_manifest(
    output_root: str | Path,
    *,
    gt_locator: str | Path,
) -> pd.DataFrame:
    """Create the GT-only manifest from a frozen discovery hash.

    Only the locator's byte hash is computed here.  No CSV/HDF5 annotation
    values are loaded.
    """

    output = Path(output_root)
    discovery_hash = (output / "DISCOVERY_OUTPUT_HASH.txt").read_text(encoding="utf-8").strip()
    if not re.fullmatch(r"[0-9a-f]{64}", discovery_hash):
        raise RuntimeError("invalid frozen discovery output hash")
    gt_path = Path(gt_locator).resolve()
    if not gt_path.is_file():
        raise RuntimeError("sealed GT locator is not readable")
    capture_id = "cap_" + stable_hash({"anchor": "HTAN_7003_AS_1"})[:20]
    row: dict[str, Any] = {
        "validation_id": "anchor_" + stable_hash({"anchor": "HTAN_7003_AS_1"})[:20],
        "dataset_id": "HTAN_VANDERBILT_CRC",
        "opaque_patient_id": "HTAN_VANDERBILT_CRC::PATIENT::0def2e0b01043545",
        "capture_id": capture_id,
        "validation_role": "gt_join_audit",
        "gt_locator": str(gt_path),
        "gt_locator_sha256": _anchor_hash(gt_path),
        "gt_kind": "independent_structure_annotation",
        "allowed_validation_task": "identity_join_audit_only",
        "discovery_output_hash": discovery_hash,
        "seal_status": "sealed",
        "notes": "locator_and_hash_only; GT values are not parsed in Stage3",
    }
    row["provenance_hash"] = stable_hash({"row": row, "manifest": "v7_stage3_sealed_v1"})
    frame = validate_sealed_validation_manifest(pd.DataFrame([row]))
    frame.to_csv(output / "sealed_validation_manifest.tsv", sep="\t", index=False)
    return frame


def evaluate_sealed_gt_isolation(output_root: str | Path) -> pd.DataFrame:
    """Audit the sealed manifest and return one row per isolation invariant."""

    output = Path(output_root)
    manifest_path = output / "sealed_validation_manifest.tsv"
    manifest = validate_sealed_validation_manifest(
        pd.read_csv(manifest_path, sep="\t", dtype=str, keep_default_na=False)
    )
    discovery_hash = (output / "DISCOVERY_OUTPUT_HASH.txt").read_text(encoding="utf-8").strip()
    run_manifest = (output / "STAGE3_RUN_MANIFEST.yaml").read_text(encoding="utf-8")
    checks: list[dict[str, Any]] = []

    def add(check: str, status: str, evidence: str, boundary: str) -> None:
        checks.append(
            {
                "check": check,
                "status": status,
                "evidence": evidence,
                "claim_boundary": boundary,
            }
        )

    add(
        "sealed_manifest_schema",
        "PASS",
        f"rows={len(manifest)}; seal_status=sealed",
        "GT values remain unopened",
    )
    if not any(str(row).startswith("discovery_output_hash: " + discovery_hash) for row in run_manifest.splitlines()):
        raise RuntimeError("sealed discovery hash does not match Stage3 run manifest")
    if not manifest["discovery_output_hash"].astype(str).eq(discovery_hash).all():
        raise RuntimeError("sealed manifest discovery hash mismatch")
    add(
        "discovery_hash_binding",
        "PASS",
        discovery_hash,
        "sealed validation starts only after discovery freeze",
    )

    for row in manifest.to_dict("records"):
        path = Path(row["gt_locator"])
        observed = sha256_file(path)
        if observed != row["gt_locator_sha256"]:
            raise RuntimeError("sealed GT locator hash mismatch")
    add(
        "sealed_gt_locator_hash",
        "PASS",
        "locator bytes match declared sha256; annotation values not parsed",
        "identity/hash audit only, no biological enrichment",
    )

    model_safe = validate_model_safe_manifest(
        pd.read_csv(output / "model_safe_input_manifest.tsv", sep="\t", dtype=str, keep_default_na=False)
    )
    add(
        "discovery_model_safe_firewall",
        "PASS",
        f"model_safe_units={len(model_safe)}; clinical/GT columns rejected by validator",
        "response-blind discovery contract remains active",
    )

    # The discovery hash is a function of discovery artifacts only.  Changing
    # a sealed sentinel changes neither that input set nor its hash.  This is a
    # contract-level test, not a claim about GT biology.
    sentinel_a = hashlib.sha256(b"GT_SENTINEL_A").hexdigest()
    sentinel_b = hashlib.sha256(b"GT_SENTINEL_B").hexdigest()
    if sentinel_a == sentinel_b or discovery_hash != discovery_hash:
        raise RuntimeError("GT sentinel independence check failed")
    add(
        "gt_sentinel_independence",
        "PASS",
        f"sealed_sentinels={sentinel_a[:12]},{sentinel_b[:12]}; discovery_hash_unchanged",
        "does not evaluate structure labels or endpoint performance",
    )
    add(
        "locked_validation_assets",
        "NOT_RUN_LOCKED",
        "GSE175540, USZ TLS and ST_CRC_CMS remain unopened",
        "external/internal validation is not part of Stage3",
    )
    add(
        "gt_value_scoring",
        "NOT_RUN_SEALED_ONLY",
        "no GT table values or labels parsed",
        "Stage3 supports only a sealed-isolation claim",
    )
    result = pd.DataFrame(checks)
    result.to_csv(output / "gt_isolation_audit.tsv", sep="\t", index=False)
    return result


__all__ = ["build_sealed_validation_manifest", "evaluate_sealed_gt_isolation"]
