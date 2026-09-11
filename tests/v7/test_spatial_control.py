from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
import pytest
import yaml

from src.v7.spatial_control import (
    FIXED_HTAN_CAPTURES,
    MODEL_SAFE_REQUIRED_COLUMNS,
    PILOT_COLUMNS,
    R04_INPUT_COLUMNS,
    SEALED_VALIDATION_REQUIRED_COLUMNS,
    STAGE1_TECHNICAL_COLUMNS,
    Stage3ControlError,
    file_provenance,
    load_r04_technical_manifest,
    load_stage1_technical_registry,
    select_stage3_pilot,
    stable_dataframe_hash,
    stable_hash,
    validate_model_safe_manifest,
    validate_pilot_manifest,
    validate_sealed_validation_manifest,
)


ROOT = Path(__file__).resolve().parents[2]


def _technical_row(
    physical_unit_id: str,
    source_id: str,
    patient_id: str,
    *,
    platform: str = "Visium",
    resolution: str = "spot",
    counts_layer: str = "h5ad_or_spot_matrix",
) -> dict[str, str]:
    return {
        "physical_unit_id": physical_unit_id,
        "source_id": source_id,
        "dataset_id": source_id.removeprefix("006_"),
        "patient_id": patient_id,
        "block_id": "BLOCK",
        "modality": "spatial",
        "platform": platform,
        "resolution": resolution,
        "coordinate_system": "native_coordinates",
        "counts_layer": counts_layer,
        "image_available": "unknown",
        "segmentation_available": "unknown",
        "identity_confidence": "explicit_manifest_patient",
        "duplicate_lineage": patient_id,
        "leakage_group_id": patient_id,
        "record_status": (
            "RESOLVED_INCLUDED_CANDIDATE"
            if source_id == "013_HTAN_CRC"
            else "EFFECTIVE_SPATIAL_MANIFEST"
        ),
        "n_observations": "100",
        "audit_status": "technical_qc_linked",
    }


def _write_stage1(path: Path, clinical_marker: str) -> None:
    rows: list[dict[str, str]] = []
    for capture in FIXED_HTAN_CAPTURES:
        rows.append(
            _technical_row(
                f"HTAN::{capture}_filtered_trimmed.h5ad",
                "013_HTAN_CRC",
                "HTAN_PATIENT_A" if capture in {"7003_AS_1", "7003_AS_2"} else f"HTAN_{capture}",
                counts_layer="replay_indexed",
            )
        )
    for index in range(1, 5):
        rows.append(
            _technical_row(
                f"HCC_UNIT_{index}", "006_GSE238264", f"GSE238264::HCC{index}"
            )
        )
        rows.append(
            _technical_row(
                f"BCC_UNIT_{index}",
                "006_GSE291246",
                f"GSE291246::su{index:03d}",
                platform="Xenium",
                resolution="cell",
            )
        )
    rows.append(
        _technical_row(
            "BCC_INELIGIBLE", "006_GSE291246", "GSE291246::su999", counts_layer="manifest_only"
        )
    )
    for row in rows:
        # These columns exist in the Stage 1 source but must never materialize
        # in the technical loader or affect selection.
        row["response"] = clinical_marker
        row["treatment"] = f"therapy_{clinical_marker}"
        row["timepoint"] = f"visit_{clinical_marker}"
        row["endpoint"] = f"endpoint_{clinical_marker}"
        row["source_path"] = f"/data/{clinical_marker}_encoded_file.h5ad"
    pd.DataFrame(rows).to_csv(path, sep="\t", index=False)


def _write_r04(path: Path) -> None:
    rows = []
    for index, capture in enumerate(FIXED_HTAN_CAPTURES):
        section = f"HTAN::{capture}_filtered_trimmed.h5ad"
        rows.append(
            {
                "section_id": section,
                "patient_id": "HTAN_PATIENT_A" if capture in {"7003_AS_1", "7003_AS_2"} else f"HTAN_{capture}",
                "block_id": "BLOCK_A" if capture in {"7003_AS_1", "7003_AS_2"} else f"BLOCK_{index}",
                "lineage": "HTAN_VANDERBILT_CRC",
                "matrix_locator": f"/read-only/matrix_{index}.h5ad",
                "coordinate_locator": f"/read-only/coordinates_{index}.csv",
                "matrix_kind": "h5ad_counts",
                "read_authorization": "HTAN_COUNTS_AND_COORDINATES",
                "count_key": "X",
                "outer_fold": f"training::{index}",
                "leakage_group_id": "HTAN_LINEAGE",
                "primary_role": "training",
            }
        )
    path.write_text(
        json.dumps(
            {
                "schema": "r04.locator_manifest.v1",
                "status": "READY",
                "input_manifest_hash": stable_hash(rows),
                "rows": rows,
            },
            sort_keys=True,
        ),
        encoding="utf-8",
    )


def test_technical_loader_never_materializes_clinical_or_filename_columns(tmp_path: Path) -> None:
    first_path = tmp_path / "registry_a.tsv"
    second_path = tmp_path / "registry_b.tsv"
    _write_stage1(first_path, "R")
    _write_stage1(second_path, "NR")

    first = load_stage1_technical_registry(first_path)
    second = load_stage1_technical_registry(second_path)
    assert tuple(first.columns) == STAGE1_TECHNICAL_COLUMNS
    assert tuple(second.columns) == STAGE1_TECHNICAL_COLUMNS
    pd.testing.assert_frame_equal(first, second)
    assert not {"response", "treatment", "timepoint", "endpoint", "source_path"} & set(first)


def test_r04_locator_values_are_verified_then_reduced_to_hashes(tmp_path: Path) -> None:
    path = tmp_path / "r04.json"
    _write_r04(path)
    frame, provenance = load_r04_technical_manifest(path)

    assert len(frame) == 4
    assert "matrix_locator" not in frame
    assert "coordinate_locator" not in frame
    assert frame["matrix_locator_hash"].str.fullmatch(r"[0-9a-f]{64}").all()
    assert provenance["payload_hash"] == json.loads(path.read_text())["input_manifest_hash"]

    payload = json.loads(path.read_text())
    payload["input_manifest_hash"] = "0" * 64
    path.write_text(json.dumps(payload), encoding="utf-8")
    with pytest.raises(Stage3ControlError, match="STAGE3_BLOCKED_R04_PAYLOAD_HASH"):
        load_r04_technical_manifest(path)


def test_stable_selection_uses_patient_key_and_technical_eligibility_only(tmp_path: Path) -> None:
    stage1_path = tmp_path / "registry.tsv"
    r04_path = tmp_path / "r04.json"
    _write_stage1(stage1_path, "R")
    _write_r04(r04_path)
    registry = load_stage1_technical_registry(stage1_path)
    r04, provenance = load_r04_technical_manifest(r04_path)
    kwargs = {
        "stage1_registry_sha256": "a" * 64,
        "r04_manifest_sha256": str(provenance["content_sha256"]),
        "r04_payload_hash": str(provenance["payload_hash"]),
    }
    first = select_stage3_pilot(registry, r04, **kwargs)
    second = select_stage3_pilot(registry.sample(frac=1, random_state=4), r04.sample(frac=1, random_state=5), **kwargs)

    pd.testing.assert_frame_equal(first.pilot, second.pilot)
    assert len(first.pilot) == 8
    assert first.pilot.groupby("source_id").size().to_dict() == {
        "006_GSE238264": 2,
        "006_GSE291246": 2,
        "013_HTAN_CRC": 4,
    }
    assert "GSE291246::su999" not in set(first.pilot["opaque_patient_id"])
    assert first.pilot.loc[first.pilot["source_id"].eq("013_HTAN_CRC"), "opaque_patient_id"].nunique() == 3
    assert stable_dataframe_hash(first.pilot, PILOT_COLUMNS) == first.provenance["pilot_hash"]
    validate_pilot_manifest(first.pilot)
    leaked_signal = first.pilot.copy()
    leaked_signal.loc[0, "selection_mode"] = "rank_by_response"
    with pytest.raises(Stage3ControlError, match="PROHIBITED_SELECTION_SIGNAL"):
        validate_pilot_manifest(leaked_signal)


def _safe_manifest() -> pd.DataFrame:
    row = {column: "value" for column in MODEL_SAFE_REQUIRED_COLUMNS}
    row.update(
        {
            "unit_id": "U1",
            "matrix_locator": "/data/R_encoded_but_not_a_selection_signal.h5",
            "coordinate_locator": "/data/NR_coordinates.csv",
            "matrix_locator_sha256": "a" * 64,
            "coordinate_locator_sha256": "b" * 64,
            "selection_hash": "c" * 64,
            "provenance_hash": "d" * 64,
            "counts_semantics": "integer_raw_counts",
        }
    )
    return pd.DataFrame([row])


def test_model_safe_manifest_is_strict_allowlist() -> None:
    safe = _safe_manifest()
    validated = validate_model_safe_manifest(safe)
    assert validated.loc[0, "unit_id"] == "U1"

    leaked = safe.assign(response="R")
    with pytest.raises(Stage3ControlError, match="EXTRA_COLUMNS"):
        validate_model_safe_manifest(leaked)

    normalized = safe.assign(counts_semantics="log_normalized")
    with pytest.raises(Stage3ControlError, match="COUNTS_SEMANTICS"):
        validate_model_safe_manifest(normalized)

    masked = safe.assign(tissue_mask_locator="/data/gt.tsv", tissue_mask_source="structure_gt")
    with pytest.raises(Stage3ControlError, match="GT_DERIVED_MASK"):
        validate_model_safe_manifest(masked)


def test_sealed_manifest_is_gt_only_but_clinical_columns_are_rejected() -> None:
    row = {column: "value" for column in SEALED_VALIDATION_REQUIRED_COLUMNS}
    row.update(
        {
            "validation_id": "V1",
            "validation_role": "external_validation",
            "gt_locator_sha256": "a" * 64,
            "discovery_output_hash": "b" * 64,
            "provenance_hash": "c" * 64,
            "seal_status": "sealed",
        }
    )
    safe = pd.DataFrame([row])
    assert validate_sealed_validation_manifest(safe).loc[0, "validation_id"] == "V1"

    leaked = safe.assign(endpoint="forbidden")
    with pytest.raises(Stage3ControlError, match="EXTRA_COLUMNS"):
        validate_sealed_validation_manifest(leaked)


def test_provenance_does_not_emit_locator(tmp_path: Path) -> None:
    source = tmp_path / "source.txt"
    source.write_text("payload", encoding="utf-8")
    provenance = file_provenance(source)
    assert set(provenance) == {"locator_sha256", "content_sha256", "bytes", "exists"}
    assert str(source) not in json.dumps(provenance)


def test_checked_in_control_config_and_pilot_are_consistent() -> None:
    config = yaml.safe_load((ROOT / "config/v7/stage3.yaml").read_text(encoding="utf-8"))
    assert config["schema"] == "v7.stage3.spatial_control.v1"
    assert config["selection"]["technical_columns_only"] is True
    assert config["firewall"]["filename_semantics_used_for_hash_ranking"] is False
    assert config["firewall"]["clinical_fields_used_for_selection"] is False
    assert config["firewall"]["discovery_pipeline_may_read_sealed_manifest"] is False

    pilot = pd.read_csv(
        ROOT / "config/v7/stage3_pilot.tsv", sep="\t", dtype=str, keep_default_na=False
    )
    validated = validate_pilot_manifest(pilot)
    assert validated.groupby("source_id").size().to_dict() == {
        "006_GSE238264": 2,
        "006_GSE291246": 2,
        "013_HTAN_CRC": 4,
    }
    assert set(validated.loc[validated["source_id"].eq("013_HTAN_CRC"), "capture_id"]) == set(
        FIXED_HTAN_CAPTURES
    )
