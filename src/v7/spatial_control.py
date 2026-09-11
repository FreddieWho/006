"""Outcome-blind control plane for the v7 Stage 3 spatial pilot.

This module deliberately does not read spatial expression.  It filters the
Stage 1 registry through an explicit technical allowlist, verifies the
read-only R-04 locator manifest, creates a reproducible pilot selection, and
validates the discovery/validation manifest firewall.
"""

from __future__ import annotations

from dataclasses import dataclass
import hashlib
import json
import math
from pathlib import Path
import re
from typing import Any, Iterable, Mapping, Sequence

import pandas as pd


class Stage3ControlError(RuntimeError):
    """Fail-closed Stage 3 control-plane violation."""


CONTROL_SCHEMA = "v7.stage3.spatial_control.v1"
R04_SCHEMA = "r04.locator_manifest.v1"
PILOT_SELECTION_NAMESPACE = "v7_stage3_pilot_patient_v1"
FIXED_HTAN_CAPTURES = ("7003_AS_1", "7003_AS_2", "7003_AS_3", "6723_KL_1")
PATIENT_QUOTAS = {"006_GSE238264": 2, "006_GSE291246": 2}

STAGE1_TECHNICAL_COLUMNS = (
    "physical_unit_id",
    "source_id",
    "dataset_id",
    "patient_id",
    "block_id",
    "modality",
    "platform",
    "resolution",
    "coordinate_system",
    "counts_layer",
    "image_available",
    "segmentation_available",
    "identity_confidence",
    "duplicate_lineage",
    "leakage_group_id",
    "record_status",
    "n_observations",
    "audit_status",
)

R04_INPUT_COLUMNS = (
    "section_id",
    "patient_id",
    "block_id",
    "lineage",
    "matrix_locator",
    "coordinate_locator",
    "matrix_kind",
    "read_authorization",
    "count_key",
    "outer_fold",
    "leakage_group_id",
    "primary_role",
)

PILOT_COLUMNS = (
    "pilot_unit_id",
    "source_id",
    "dataset_id",
    "unit_kind",
    "opaque_patient_id",
    "opaque_block_id",
    "opaque_section_id",
    "capture_id",
    "platform",
    "modality",
    "resolution",
    "counts_layer",
    "coordinate_system",
    "identity_status",
    "leakage_group_id",
    "member_capture_count",
    "member_set_hash",
    "selection_mode",
    "selection_hash",
    "technical_eligibility",
    "stage1_registry_sha256",
    "r04_manifest_sha256",
    "r04_payload_hash",
)

MODEL_SAFE_REQUIRED_COLUMNS = frozenset(
    {
        "unit_id",
        "dataset_id",
        "opaque_patient_id",
        "capture_id",
        "modality",
        "platform",
        "matrix_locator",
        "matrix_locator_sha256",
        "matrix_kind",
        "count_key",
        "counts_semantics",
        "coordinate_locator",
        "coordinate_locator_sha256",
        "coordinate_system",
        "read_authorization",
        "leakage_group_id",
        "selection_hash",
        "provenance_hash",
    }
)
MODEL_SAFE_ALLOWED_COLUMNS = MODEL_SAFE_REQUIRED_COLUMNS | frozenset(
    {
        "opaque_block_id",
        "opaque_section_id",
        "scale_status",
        "scalefactor_locator",
        "scalefactor_sha256",
        "image_locator",
        "image_sha256",
        "image_usage",
        "segmentation_locator",
        "segmentation_sha256",
        "tissue_mask_locator",
        "tissue_mask_sha256",
        "tissue_mask_source",
        "source_registry_row_hash",
        "qc_status",
    }
)

SEALED_VALIDATION_REQUIRED_COLUMNS = frozenset(
    {
        "validation_id",
        "dataset_id",
        "opaque_patient_id",
        "capture_id",
        "validation_role",
        "gt_locator",
        "gt_locator_sha256",
        "gt_kind",
        "allowed_validation_task",
        "discovery_output_hash",
        "seal_status",
        "provenance_hash",
    }
)
SEALED_VALIDATION_ALLOWED_COLUMNS = SEALED_VALIDATION_REQUIRED_COLUMNS | frozenset(
    {"opaque_block_id", "opaque_section_id", "notes"}
)

_PROHIBITED_CLINICAL_COLUMN_TOKENS = (
    "response",
    "responder",
    "treatment",
    "timepoint",
    "endpoint",
    "outcome",
    "survival",
    "recist",
    "progression",
    "recurrence",
    "relapse",
)
_PROHIBITED_DISCOVERY_COLUMN_TOKENS = _PROHIBITED_CLINICAL_COLUMN_TOKENS + (
    "ground_truth",
    "structure_gt",
    "gt_distance",
    "gt_mask",
    "boundary_annotation",
)
_SHA256_PATTERN = re.compile(r"^[0-9a-f]{64}$")


@dataclass(frozen=True)
class PilotSelection:
    """Deterministic pilot plus complete technical and rank audits."""

    pilot: pd.DataFrame
    technical_audit: pd.DataFrame
    selection_audit: pd.DataFrame
    provenance: Mapping[str, Any]


def _jsonable(value: Any) -> Any:
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, Mapping):
        return {str(key): _jsonable(item) for key, item in value.items()}
    if isinstance(value, (list, tuple, set, frozenset)):
        values = [_jsonable(item) for item in value]
        return sorted(values, key=lambda item: json.dumps(item, sort_keys=True)) if isinstance(value, (set, frozenset)) else values
    if hasattr(value, "item") and callable(value.item):
        value = value.item()
    if value is None or pd.isna(value):
        return None
    if isinstance(value, float) and not math.isfinite(value):
        raise Stage3ControlError("STAGE3_BLOCKED_HASH_INPUT: non-finite value")
    return value


def stable_hash(value: Any) -> str:
    """SHA256 of canonical, compact JSON."""

    payload = json.dumps(
        _jsonable(value),
        ensure_ascii=False,
        sort_keys=True,
        separators=(",", ":"),
        allow_nan=False,
    ).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def sha256_file(path: str | Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def file_provenance(path: str | Path) -> dict[str, Any]:
    """Return content provenance without exposing the filesystem locator."""

    source = Path(path)
    if not source.is_file():
        raise Stage3ControlError("STAGE3_BLOCKED_MISSING_INPUT")
    return {
        "locator_sha256": stable_hash({"locator": str(source.resolve())}),
        "content_sha256": sha256_file(source),
        "bytes": source.stat().st_size,
        "exists": True,
    }


def stable_dataframe_hash(frame: pd.DataFrame, columns: Sequence[str] | None = None) -> str:
    """Hash a table independent of its current row order."""

    selected = list(columns) if columns is not None else sorted(map(str, frame.columns))
    missing = sorted(set(selected) - set(frame.columns))
    if missing:
        raise Stage3ControlError(
            "STAGE3_BLOCKED_HASH_SCHEMA: missing columns " + ",".join(missing)
        )
    records = [
        {column: _jsonable(row[column]) for column in selected}
        for row in frame[selected].to_dict("records")
    ]
    records.sort(key=lambda row: json.dumps(row, ensure_ascii=False, sort_keys=True, separators=(",", ":")))
    return stable_hash({"columns": selected, "records": records})


def _forbidden_column_hits(columns: Iterable[str], tokens: Sequence[str]) -> list[str]:
    hits: list[str] = []
    for column in columns:
        normalized = str(column).lower()
        for token in tokens:
            pattern = rf"(^|[^a-z0-9]){re.escape(token)}([^a-z0-9]|$)"
            if re.search(pattern, normalized):
                hits.append(str(column))
                break
    return sorted(set(hits))


def _forbidden_value_columns(
    frame: pd.DataFrame, columns: Iterable[str], tokens: Sequence[str]
) -> list[str]:
    """Return semantic columns containing forbidden whole tokens, never values."""

    hits: list[str] = []
    for column in columns:
        if column not in frame:
            continue
        for value in frame[column].astype(str):
            normalized = value.lower()
            if any(
                re.search(rf"(^|[^a-z0-9]){re.escape(token)}([^a-z0-9]|$)", normalized)
                for token in tokens
            ):
                hits.append(column)
                break
    return sorted(set(hits))


def _require_columns(columns: Iterable[str], required: Iterable[str], context: str) -> None:
    missing = sorted(set(required) - set(columns))
    if missing:
        raise Stage3ControlError(
            f"STAGE3_BLOCKED_{context}_SCHEMA: missing " + ",".join(missing)
        )


def _canonical_frame(frame: pd.DataFrame, columns: Iterable[str], sort_by: Sequence[str]) -> pd.DataFrame:
    order = [column for column in columns if column in frame.columns]
    result = frame.reindex(columns=order).copy()
    if not result.empty:
        result = result.sort_values(list(sort_by), kind="mergesort", na_position="last")
    return result.reset_index(drop=True)


def load_stage1_technical_registry(path: str | Path) -> pd.DataFrame:
    """Read only the technical Stage 1 allowlist; clinical columns never materialize."""

    source = Path(path)
    header = pd.read_csv(source, sep="\t", nrows=0)
    _require_columns(header.columns, STAGE1_TECHNICAL_COLUMNS, "STAGE1_REGISTRY")
    frame = pd.read_csv(
        source,
        sep="\t",
        usecols=list(STAGE1_TECHNICAL_COLUMNS),
        dtype=str,
        keep_default_na=False,
    )
    if tuple(frame.columns) != STAGE1_TECHNICAL_COLUMNS:
        frame = frame.reindex(columns=STAGE1_TECHNICAL_COLUMNS)
    return _canonical_frame(frame, STAGE1_TECHNICAL_COLUMNS, ["source_id", "physical_unit_id"])


def load_r04_technical_manifest(path: str | Path) -> tuple[pd.DataFrame, dict[str, Any]]:
    """Verify R-04 provenance and return hash-only locator records."""

    source = Path(path)
    payload = json.loads(source.read_text(encoding="utf-8"))
    if payload.get("schema") != R04_SCHEMA or payload.get("status") != "READY":
        raise Stage3ControlError("STAGE3_BLOCKED_R04_MANIFEST_STATUS")
    rows = payload.get("rows")
    if not isinstance(rows, list) or not rows:
        raise Stage3ControlError("STAGE3_BLOCKED_R04_MANIFEST_ROWS")
    declared_hash = str(payload.get("input_manifest_hash", ""))
    observed_hash = stable_hash(rows)
    if declared_hash != observed_hash:
        raise Stage3ControlError("STAGE3_BLOCKED_R04_PAYLOAD_HASH")
    for row in rows:
        if not isinstance(row, Mapping):
            raise Stage3ControlError("STAGE3_BLOCKED_R04_MANIFEST_ROWS")
        _require_columns(row, R04_INPUT_COLUMNS, "R04_MANIFEST")
        hits = _forbidden_column_hits(row, _PROHIBITED_CLINICAL_COLUMN_TOKENS)
        if hits:
            raise Stage3ControlError(
                "STAGE3_BLOCKED_R04_PROHIBITED_COLUMNS: " + ",".join(hits)
            )

    records: list[dict[str, Any]] = []
    for row in rows:
        safe = {key: str(row.get(key, "")) for key in R04_INPUT_COLUMNS if "locator" not in key}
        safe.update(
            {
                "physical_unit_id": str(row["section_id"]),
                "matrix_locator_hash": stable_hash({"locator": str(row["matrix_locator"])}),
                "coordinate_locator_hash": stable_hash({"locator": str(row["coordinate_locator"])}),
                "r04_row_hash": stable_hash(row),
                "r04_payload_hash": declared_hash,
            }
        )
        records.append(safe)
    frame = pd.DataFrame(records)
    if frame["physical_unit_id"].duplicated().any():
        raise Stage3ControlError("STAGE3_BLOCKED_R04_DUPLICATE_SECTION")
    frame = frame.sort_values("physical_unit_id", kind="mergesort").reset_index(drop=True)
    provenance = file_provenance(source)
    provenance["payload_hash"] = declared_hash
    provenance["row_count"] = len(frame)
    return frame, provenance


def technical_eligibility_audit(registry: pd.DataFrame) -> pd.DataFrame:
    """Evaluate only platform, count, coordinate, identity and read-QC fields."""

    _require_columns(registry.columns, STAGE1_TECHNICAL_COLUMNS, "TECHNICAL_REGISTRY")
    records: list[dict[str, Any]] = []
    allowed_record_status = {"EFFECTIVE_SPATIAL_MANIFEST", "RESOLVED_INCLUDED_CANDIDATE"}
    allowed_counts = {"h5ad_or_spot_matrix", "replay_indexed"}
    unavailable = {"", "unknown", "not_available", "not_applicable", "none"}
    for row in registry.to_dict("records"):
        reasons: list[str] = []
        if str(row["modality"]).strip().lower() != "spatial":
            reasons.append("modality_not_spatial")
        if str(row["counts_layer"]).strip() not in allowed_counts:
            reasons.append("counts_not_technically_available")
        if str(row["coordinate_system"]).strip().lower() in unavailable:
            reasons.append("coordinates_unavailable")
        if not str(row["patient_id"]).strip():
            reasons.append("patient_identity_missing")
        if not str(row["physical_unit_id"]).strip():
            reasons.append("physical_unit_missing")
        if str(row["record_status"]).strip() not in allowed_record_status:
            reasons.append("record_status_not_current")
        if str(row["identity_confidence"]).strip().lower() in unavailable:
            reasons.append("identity_confidence_unavailable")
        try:
            observations = int(float(str(row["n_observations"])))
            if observations <= 0:
                reasons.append("nonpositive_observation_count")
        except (TypeError, ValueError):
            reasons.append("observation_count_unavailable")
        records.append(
            {
                "physical_unit_id": str(row["physical_unit_id"]),
                "source_id": str(row["source_id"]),
                "dataset_id": str(row["dataset_id"]),
                "opaque_patient_id": str(row["patient_id"]),
                "technical_eligibility": "eligible" if not reasons else "ineligible",
                "technical_reasons": "|".join(sorted(reasons)),
                "technical_row_hash": stable_hash(
                    {column: row[column] for column in STAGE1_TECHNICAL_COLUMNS}
                ),
            }
        )
    return _canonical_frame(
        pd.DataFrame(records),
        records[0].keys() if records else (),
        ["source_id", "opaque_patient_id", "physical_unit_id"],
    )


def _capture_id(section_id: str) -> str:
    value = section_id.removeprefix("HTAN::")
    for suffix in ("_filtered_trimmed.h5ad", "_filtered.h5ad", ".h5ad"):
        if value.endswith(suffix):
            return value[: -len(suffix)]
    return value


def _pilot_row_hashes(rows: pd.DataFrame) -> tuple[int, str]:
    ids = sorted(rows["physical_unit_id"].astype(str))
    return len(ids), stable_hash({"physical_unit_ids": ids})


def select_stage3_pilot(
    registry: pd.DataFrame,
    r04: pd.DataFrame,
    *,
    stage1_registry_sha256: str,
    r04_manifest_sha256: str,
    r04_payload_hash: str,
    selection_namespace: str = PILOT_SELECTION_NAMESPACE,
    fixed_htan_captures: Sequence[str] = FIXED_HTAN_CAPTURES,
    patient_quotas: Mapping[str, int] = PATIENT_QUOTAS,
) -> PilotSelection:
    """Select the eight-unit pilot without clinical fields or filename semantics."""

    technical = technical_eligibility_audit(registry)
    eligible_ids = set(
        technical.loc[technical["technical_eligibility"].eq("eligible"), "physical_unit_id"]
    )
    eligible = registry[registry["physical_unit_id"].isin(eligible_ids)].copy()
    pilot_rows: list[dict[str, Any]] = []
    audit_rows: list[dict[str, Any]] = []

    r04_material = r04.copy()
    r04_material["capture_id"] = r04_material["section_id"].astype(str).map(_capture_id)
    if r04_material["capture_id"].duplicated().any():
        raise Stage3ControlError("STAGE3_BLOCKED_HTAN_CAPTURE_ALIAS_COLLISION")
    r04_by_capture = r04_material.set_index("capture_id", drop=False)
    registry_by_id = eligible.set_index("physical_unit_id", drop=False)
    for capture in fixed_htan_captures:
        if capture not in r04_by_capture.index:
            raise Stage3ControlError("STAGE3_BLOCKED_HTAN_FIXED_CAPTURE_MISSING")
        source = r04_by_capture.loc[capture]
        physical_id = str(source["physical_unit_id"])
        if physical_id not in registry_by_id.index:
            raise Stage3ControlError("STAGE3_BLOCKED_HTAN_FIXED_CAPTURE_INELIGIBLE")
        registry_row = registry_by_id.loc[physical_id]
        selection_hash = stable_hash(
            {"namespace": selection_namespace, "mode": "fixed_capture", "capture_id": capture}
        )
        pilot_rows.append(
            {
                "pilot_unit_id": f"013_HTAN_CRC::capture::{capture}",
                "source_id": "013_HTAN_CRC",
                "dataset_id": str(registry_row["dataset_id"]),
                "unit_kind": "technical_capture",
                "opaque_patient_id": str(source["patient_id"]),
                "opaque_block_id": str(source["block_id"]),
                "opaque_section_id": physical_id,
                "capture_id": capture,
                "platform": str(registry_row["platform"]),
                "modality": str(registry_row["modality"]),
                "resolution": str(registry_row["resolution"]),
                "counts_layer": str(registry_row["counts_layer"]),
                "coordinate_system": str(registry_row["coordinate_system"]),
                "identity_status": str(registry_row["identity_confidence"]),
                "leakage_group_id": str(source["leakage_group_id"]),
                "member_capture_count": 1,
                "member_set_hash": stable_hash({"physical_unit_ids": [physical_id]}),
                "selection_mode": "fixed_predeclared_technical_capture",
                "selection_hash": selection_hash,
                "technical_eligibility": "eligible",
                "stage1_registry_sha256": stage1_registry_sha256,
                "r04_manifest_sha256": r04_manifest_sha256,
                "r04_payload_hash": r04_payload_hash,
            }
        )
        audit_rows.append(
            {
                "source_id": "013_HTAN_CRC",
                "candidate_kind": "technical_capture",
                "candidate_id": capture,
                "selection_hash": selection_hash,
                "selection_rank": list(fixed_htan_captures).index(capture) + 1,
                "selected": True,
                "selection_reason": "fixed_in_stage3_technical_plan",
            }
        )

    for source_id, quota in sorted(patient_quotas.items()):
        candidates = eligible[eligible["source_id"].eq(source_id)].copy()
        if candidates.empty:
            raise Stage3ControlError("STAGE3_BLOCKED_PATIENT_PILOT_NO_ELIGIBLE_UNITS")
        patient_groups = list(candidates.groupby("patient_id", sort=True))
        if len(patient_groups) < quota:
            raise Stage3ControlError("STAGE3_BLOCKED_PATIENT_PILOT_QUOTA")
        ranked: list[tuple[str, str, pd.DataFrame]] = []
        for patient_id, group in patient_groups:
            rank_hash = stable_hash(
                {
                    "namespace": selection_namespace,
                    "parts": [source_id, str(patient_id)],
                }
            )
            ranked.append((rank_hash, str(patient_id), group))
        ranked.sort(key=lambda item: (item[0], item[1]))
        selected_patients = {patient for _, patient, _ in ranked[:quota]}
        for rank, (rank_hash, patient_id, group) in enumerate(ranked, start=1):
            is_selected = patient_id in selected_patients
            audit_rows.append(
                {
                    "source_id": source_id,
                    "candidate_kind": "patient_unit",
                    "candidate_id": patient_id,
                    "selection_hash": rank_hash,
                    "selection_rank": rank,
                    "selected": is_selected,
                    "selection_reason": "stable_patient_hash_rank" if is_selected else "outside_fixed_quota",
                }
            )
            if not is_selected:
                continue
            member_count, member_hash = _pilot_row_hashes(group)
            first = group.sort_values("physical_unit_id", kind="mergesort").iloc[0]
            unique = lambda column: sorted(set(group[column].astype(str)))
            for column in ("dataset_id", "platform", "modality", "resolution", "counts_layer", "coordinate_system", "identity_confidence", "leakage_group_id"):
                if len(unique(column)) != 1:
                    raise Stage3ControlError(
                        f"STAGE3_BLOCKED_PATIENT_UNIT_HETEROGENEITY: {source_id}:{column}"
                    )
            pilot_rows.append(
                {
                    "pilot_unit_id": f"{source_id}::patient::{patient_id}",
                    "source_id": source_id,
                    "dataset_id": str(first["dataset_id"]),
                    "unit_kind": "patient_all_eligible_captures",
                    "opaque_patient_id": patient_id,
                    "opaque_block_id": "",
                    "opaque_section_id": "",
                    "capture_id": "",
                    "platform": str(first["platform"]),
                    "modality": str(first["modality"]),
                    "resolution": str(first["resolution"]),
                    "counts_layer": str(first["counts_layer"]),
                    "coordinate_system": str(first["coordinate_system"]),
                    "identity_status": str(first["identity_confidence"]),
                    "leakage_group_id": str(first["leakage_group_id"]),
                    "member_capture_count": member_count,
                    "member_set_hash": member_hash,
                    "selection_mode": "stable_hash_of_source_and_patient_key",
                    "selection_hash": rank_hash,
                    "technical_eligibility": "eligible",
                    "stage1_registry_sha256": stage1_registry_sha256,
                    "r04_manifest_sha256": r04_manifest_sha256,
                    "r04_payload_hash": r04_payload_hash,
                }
            )

    pilot = _canonical_frame(
        pd.DataFrame(pilot_rows), PILOT_COLUMNS, ["source_id", "pilot_unit_id"]
    )
    if len(pilot) != len(fixed_htan_captures) + sum(patient_quotas.values()):
        raise Stage3ControlError("STAGE3_BLOCKED_PILOT_UNIT_COUNT")
    if pilot["pilot_unit_id"].duplicated().any():
        raise Stage3ControlError("STAGE3_BLOCKED_PILOT_DUPLICATE_UNIT")
    selection_audit = _canonical_frame(
        pd.DataFrame(audit_rows),
        [
            "source_id",
            "candidate_kind",
            "candidate_id",
            "selection_hash",
            "selection_rank",
            "selected",
            "selection_reason",
        ],
        ["source_id", "selection_rank", "candidate_id"],
    )
    return PilotSelection(
        pilot=pilot,
        technical_audit=technical,
        selection_audit=selection_audit,
        provenance={
            "selection_namespace": selection_namespace,
            "stage1_registry_sha256": stage1_registry_sha256,
            "r04_manifest_sha256": r04_manifest_sha256,
            "r04_payload_hash": r04_payload_hash,
            "pilot_hash": stable_dataframe_hash(pilot, PILOT_COLUMNS),
        },
    )


def build_stage3_pilot(
    stage1_registry_path: str | Path,
    r04_manifest_path: str | Path,
    *,
    selection_namespace: str = PILOT_SELECTION_NAMESPACE,
) -> PilotSelection:
    """Load the two audited sources and create the deterministic pilot."""

    registry_path = Path(stage1_registry_path)
    r04_path = Path(r04_manifest_path)
    registry = load_stage1_technical_registry(registry_path)
    r04, r04_provenance = load_r04_technical_manifest(r04_path)
    return select_stage3_pilot(
        registry,
        r04,
        stage1_registry_sha256=sha256_file(registry_path),
        r04_manifest_sha256=str(r04_provenance["content_sha256"]),
        r04_payload_hash=str(r04_provenance["payload_hash"]),
        selection_namespace=selection_namespace,
    )


def validate_pilot_manifest(frame: pd.DataFrame) -> pd.DataFrame:
    _require_columns(frame.columns, PILOT_COLUMNS, "PILOT_MANIFEST")
    extras = sorted(set(frame.columns) - set(PILOT_COLUMNS))
    if extras:
        raise Stage3ControlError("STAGE3_BLOCKED_PILOT_EXTRA_COLUMNS: " + ",".join(extras))
    hits = _forbidden_column_hits(frame.columns, _PROHIBITED_CLINICAL_COLUMN_TOKENS)
    if hits:
        raise Stage3ControlError("STAGE3_BLOCKED_PILOT_PROHIBITED_COLUMNS: " + ",".join(hits))
    value_hits = _forbidden_value_columns(
        frame,
        ["selection_mode", "technical_eligibility", "identity_status"],
        _PROHIBITED_CLINICAL_COLUMN_TOKENS,
    )
    if value_hits:
        raise Stage3ControlError(
            "STAGE3_BLOCKED_PILOT_PROHIBITED_SELECTION_SIGNAL: " + ",".join(value_hits)
        )
    if len(frame) != 8 or frame["pilot_unit_id"].duplicated().any():
        raise Stage3ControlError("STAGE3_BLOCKED_PILOT_CARDINALITY")
    for column in ("selection_hash", "member_set_hash", "stage1_registry_sha256", "r04_manifest_sha256", "r04_payload_hash"):
        if not frame[column].astype(str).map(lambda value: bool(_SHA256_PATTERN.fullmatch(value))).all():
            raise Stage3ControlError(f"STAGE3_BLOCKED_PILOT_HASH: {column}")
    return _canonical_frame(frame, PILOT_COLUMNS, ["source_id", "pilot_unit_id"])


def _validate_manifest_schema(
    frame: pd.DataFrame,
    *,
    required: frozenset[str],
    allowed: frozenset[str],
    context: str,
    forbidden_tokens: Sequence[str],
    unique_column: str,
) -> pd.DataFrame:
    _require_columns(frame.columns, required, context)
    extras = sorted(set(frame.columns) - allowed)
    if extras:
        raise Stage3ControlError(f"STAGE3_BLOCKED_{context}_EXTRA_COLUMNS: " + ",".join(extras))
    hits = _forbidden_column_hits(frame.columns, forbidden_tokens)
    if hits:
        raise Stage3ControlError(f"STAGE3_BLOCKED_{context}_PROHIBITED_COLUMNS: " + ",".join(hits))
    if frame.empty or frame[unique_column].astype(str).str.strip().eq("").any():
        raise Stage3ControlError(f"STAGE3_BLOCKED_{context}_EMPTY_ID")
    if frame[unique_column].duplicated().any():
        raise Stage3ControlError(f"STAGE3_BLOCKED_{context}_DUPLICATE_ID")
    blank_required = [
        column
        for column in sorted(required)
        if frame[column].astype(str).str.strip().eq("").any()
    ]
    if blank_required:
        raise Stage3ControlError(
            f"STAGE3_BLOCKED_{context}_BLANK_REQUIRED: " + ",".join(blank_required)
        )
    return _canonical_frame(frame, sorted(frame.columns), [unique_column])


def validate_model_safe_manifest(frame: pd.DataFrame) -> pd.DataFrame:
    """Validate the only manifest accepted by the discovery pipeline."""

    result = _validate_manifest_schema(
        frame,
        required=MODEL_SAFE_REQUIRED_COLUMNS,
        allowed=MODEL_SAFE_ALLOWED_COLUMNS,
        context="MODEL_SAFE_MANIFEST",
        forbidden_tokens=_PROHIBITED_DISCOVERY_COLUMN_TOKENS,
        unique_column="unit_id",
    )
    value_hits = _forbidden_value_columns(
        result,
        [
            "modality",
            "platform",
            "matrix_kind",
            "count_key",
            "counts_semantics",
            "coordinate_system",
            "read_authorization",
            "tissue_mask_source",
            "image_usage",
            "qc_status",
        ],
        _PROHIBITED_CLINICAL_COLUMN_TOKENS,
    )
    if value_hits:
        raise Stage3ControlError(
            "STAGE3_BLOCKED_MODEL_SAFE_PROHIBITED_SEMANTICS: " + ",".join(value_hits)
        )
    if not result["counts_semantics"].astype(str).eq("integer_raw_counts").all():
        raise Stage3ControlError("STAGE3_BLOCKED_MODEL_SAFE_COUNTS_SEMANTICS")
    for column in ("matrix_locator_sha256", "coordinate_locator_sha256", "selection_hash", "provenance_hash"):
        if not result[column].astype(str).map(lambda value: bool(_SHA256_PATTERN.fullmatch(value))).all():
            raise Stage3ControlError(f"STAGE3_BLOCKED_MODEL_SAFE_HASH: {column}")
    if "tissue_mask_locator" in result:
        present = result["tissue_mask_locator"].astype(str).str.strip().ne("")
        if present.any():
            if "tissue_mask_source" not in result:
                raise Stage3ControlError("STAGE3_BLOCKED_MODEL_SAFE_TISSUE_MASK_SOURCE")
            valid = result.loc[present, "tissue_mask_source"].astype(str).eq("non_gt_tissue_mask")
            if not valid.all():
                raise Stage3ControlError("STAGE3_BLOCKED_MODEL_SAFE_GT_DERIVED_MASK")
    if "image_locator" in result:
        present = result["image_locator"].astype(str).str.strip().ne("")
        if present.any() and (
            "image_usage" not in result
            or not result.loc[present, "image_usage"].astype(str).eq("locator_only_not_feature").all()
        ):
            raise Stage3ControlError("STAGE3_BLOCKED_MODEL_SAFE_IMAGE_USAGE")
    return result


def validate_sealed_validation_manifest(frame: pd.DataFrame) -> pd.DataFrame:
    """Validate the GT-only manifest kept outside discovery imports."""

    result = _validate_manifest_schema(
        frame,
        required=SEALED_VALIDATION_REQUIRED_COLUMNS,
        allowed=SEALED_VALIDATION_ALLOWED_COLUMNS,
        context="SEALED_VALIDATION_MANIFEST",
        forbidden_tokens=_PROHIBITED_CLINICAL_COLUMN_TOKENS,
        unique_column="validation_id",
    )
    value_hits = _forbidden_value_columns(
        result,
        ["validation_role", "gt_kind", "allowed_validation_task", "seal_status", "notes"],
        _PROHIBITED_CLINICAL_COLUMN_TOKENS,
    )
    if value_hits:
        raise Stage3ControlError(
            "STAGE3_BLOCKED_SEALED_VALIDATION_PROHIBITED_SEMANTICS: "
            + ",".join(value_hits)
        )
    if not result["seal_status"].astype(str).eq("sealed").all():
        raise Stage3ControlError("STAGE3_BLOCKED_VALIDATION_MANIFEST_NOT_SEALED")
    allowed_roles = {"internal_validation", "external_validation", "gt_join_audit"}
    if not result["validation_role"].astype(str).isin(allowed_roles).all():
        raise Stage3ControlError("STAGE3_BLOCKED_VALIDATION_ROLE")
    for column in ("gt_locator_sha256", "discovery_output_hash", "provenance_hash"):
        if not result[column].astype(str).map(lambda value: bool(_SHA256_PATTERN.fullmatch(value))).all():
            raise Stage3ControlError(f"STAGE3_BLOCKED_SEALED_VALIDATION_HASH: {column}")
    return result


__all__ = [
    "CONTROL_SCHEMA",
    "FIXED_HTAN_CAPTURES",
    "MODEL_SAFE_ALLOWED_COLUMNS",
    "MODEL_SAFE_REQUIRED_COLUMNS",
    "PATIENT_QUOTAS",
    "PILOT_COLUMNS",
    "PILOT_SELECTION_NAMESPACE",
    "PilotSelection",
    "SEALED_VALIDATION_ALLOWED_COLUMNS",
    "SEALED_VALIDATION_REQUIRED_COLUMNS",
    "STAGE1_TECHNICAL_COLUMNS",
    "Stage3ControlError",
    "build_stage3_pilot",
    "file_provenance",
    "load_r04_technical_manifest",
    "load_stage1_technical_registry",
    "select_stage3_pilot",
    "sha256_file",
    "stable_dataframe_hash",
    "stable_hash",
    "technical_eligibility_audit",
    "validate_model_safe_manifest",
    "validate_pilot_manifest",
    "validate_sealed_validation_manifest",
]
