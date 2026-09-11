"""Stage 3 response-blind spatial foundation pipeline.

The pipeline resolves the frozen technical pilot, loads each physical capture
through the fail-closed adapters, projects the frozen Stage 2 vocabulary, and
computes section-local spatial statistics.  It never imports the sealed GT
evaluator and never reads clinical columns.
"""

from __future__ import annotations

from dataclasses import dataclass
import json
from pathlib import Path
import subprocess
from typing import Any, Iterable, Mapping

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

from .spatial_control import (
    PILOT_COLUMNS,
    build_stage3_pilot,
    file_provenance,
    load_r04_technical_manifest,
    load_stage1_technical_registry,
    sha256_file,
    stable_dataframe_hash,
    stable_hash,
    validate_model_safe_manifest,
    validate_pilot_manifest,
)
from .spatial_io import SpatialIdentity
from .spatial_io.adapters import load_h5ad_counts, load_visium_10x, load_xenium_h5ad
from .spatial_io.sanitization import sanitize_h5ad_for_discovery
from .spatial_projection import project_stage2_features
from .spatial_stats import (
    build_section_graph,
    edge_association,
    moran_permutation_test,
    moran_statistic,
    patient_block_split,
    patient_nested_summary,
    sparse_variogram,
)


@dataclass(frozen=True)
class ResolvedSpatialUnit:
    """One physical capture selected by the outcome-blind pilot."""

    unit_id: str
    pilot_unit_id: str
    source_id: str
    dataset_id: str
    physical_unit_id: str
    opaque_patient_id: str
    opaque_block_id: str
    opaque_section_id: str
    capture_id: str
    platform: str
    modality: str
    resolution: str
    counts_layer: str
    coordinate_system: str
    leakage_group_id: str
    selection_hash: str
    identity_status: str
    matrix_kind: str
    matrix_path: Path
    coordinate_path: Path | None
    scalefactors_path: Path | None
    image_path: Path | None
    source_record_id: str
    source_registry_row_hash: str
    read_authorization: str


def _safe_id(value: str) -> str:
    return "u_" + stable_hash({"physical_unit_id": value})[:16]


def _opaque_id(prefix: str, value: str) -> str:
    return f"{prefix}_{stable_hash({prefix: value})[:20]}"


def _file_hash(path: Path | None) -> str:
    return sha256_file(path) if path is not None and path.is_file() else ""


def _source_row_hash(row: Mapping[str, Any], columns: Iterable[str]) -> str:
    return stable_hash({column: str(row.get(column, "")) for column in columns})


def _metadata_free_gsm_map(path: Path) -> dict[str, str]:
    """Read only technical sample/GSM keys; clinical columns never load."""

    frame = pd.read_csv(
        path,
        sep="\t",
        usecols=["gsm_id", "sample_id"],
        dtype=str,
        keep_default_na=False,
    )
    if frame["sample_id"].duplicated().any() or frame["gsm_id"].duplicated().any():
        raise RuntimeError("GSE291246 technical sample map is not one-to-one")
    return dict(zip(frame["sample_id"], frame["gsm_id"]))


def resolve_stage3_units(
    pilot: pd.DataFrame,
    registry: pd.DataFrame,
    r04: pd.DataFrame,
    *,
    r04_manifest_path: str | Path = "/home/huyudi/013_spatial/infra/r04/input_manifest.json",
    metadata_manifest: str | Path = (
        "results/v6_2/phase2_5_data_onboarding/01_metadata_freeze/"
        "spatial_sample_manifest.GSE291246.tsv"
    ),
) -> list[ResolvedSpatialUnit]:
    """Expand the eight logical pilot rows to their selected physical captures."""

    registry = registry.copy()
    r04 = r04.copy()
    registry["physical_unit_id"] = registry["physical_unit_id"].astype(str)
    r04["capture_id"] = r04["section_id"].astype(str).str.replace(
        "HTAN::", "", regex=False
    )
    for suffix in ("_filtered_trimmed.h5ad", "_filtered.h5ad", ".h5ad"):
        r04["capture_id"] = r04["capture_id"].str.removesuffix(suffix)
    r04_by_capture = r04.set_index("capture_id", drop=False)
    # The control-plane frame intentionally contains locator hashes only.  A
    # post-selection materialisation step may resolve the original read-only
    # locator, but never uses it for patient/capture ranking.
    raw_r04 = json.loads(Path(r04_manifest_path).read_text(encoding="utf-8"))
    raw_r04_by_capture: dict[str, Mapping[str, Any]] = {}
    for raw_row in raw_r04.get("rows", []):
        raw_capture = str(raw_row["section_id"]).removeprefix("HTAN::")
        for suffix in ("_filtered_trimmed.h5ad", "_filtered.h5ad", ".h5ad"):
            raw_capture = raw_capture.removesuffix(suffix)
        raw_r04_by_capture[raw_capture] = raw_row
    reg_by_id = registry.set_index("physical_unit_id", drop=False)
    gsm_map = _metadata_free_gsm_map(Path(metadata_manifest))
    resolved: list[ResolvedSpatialUnit] = []

    for pilot_row in pilot.to_dict("records"):
        source_id = str(pilot_row["source_id"])
        pilot_id = str(pilot_row["pilot_unit_id"])
        if source_id == "013_HTAN_CRC":
            capture = str(pilot_row["capture_id"])
            source = r04_by_capture.loc[capture]
            raw_source = raw_r04_by_capture.get(capture)
            if raw_source is None:
                raise RuntimeError(f"R04 materialization locator missing for {capture}")
            physical_id = str(source["physical_unit_id"])
            reg_row = reg_by_id.loc[physical_id]
            members = [physical_id]
            source_record = capture
            matrix_path = Path(str(raw_source["matrix_locator"]))
            coordinate_path = None
            scale_path = None
            image_path = None
            matrix_kind = str(raw_source["matrix_kind"])
            count_key = str(raw_source["count_key"])
            read_authorization = str(raw_source["read_authorization"])
            block_id = str(raw_source["block_id"])
            patient_id = str(raw_source["patient_id"])
            leakage = str(raw_source["leakage_group_id"])
        else:
            rows = registry[
                (registry["source_id"].astype(str) == source_id)
                & (registry["patient_id"].astype(str) == str(pilot_row["opaque_patient_id"]))
                & (registry["counts_layer"].astype(str) == "h5ad_or_spot_matrix")
            ].sort_values("physical_unit_id", kind="mergesort")
            if rows.empty:
                raise RuntimeError(f"no materializable technical captures for {pilot_id}")
            observed_hash = stable_hash(
                {"physical_unit_ids": sorted(rows["physical_unit_id"].astype(str))}
            )
            if observed_hash != str(pilot_row["member_set_hash"]):
                raise RuntimeError(f"pilot membership hash mismatch for {pilot_id}")
            members = list(rows["physical_unit_id"].astype(str))
            for physical_id in members:
                reg_row = reg_by_id.loc[physical_id]
                sample_id = physical_id.split("::section::", 1)[-1]
                if source_id == "006_GSE238264":
                    sample_root = Path(
                        "results/v6_2/phase2_5_data_onboarding/_scratch/"
                        "GSE238264_verify"
                    ) / sample_id
                    matrix_path = sample_root / "filtered_feature_bc_matrix.h5"
                    coordinate_path = sample_root / "spatial/tissue_positions_list.csv"
                    scale_path = sample_root / "spatial/scalefactors_json.json"
                    image_path = sample_root / "spatial/tissue_hires_image.png"
                    matrix_kind = "10x_h5"
                    count_key = "matrix/raw_umi"
                    read_authorization = "LOCAL_QC_VISIUM_COUNTS_COORDINATES"
                elif source_id == "006_GSE291246":
                    gsm = gsm_map.get(sample_id)
                    if not gsm:
                        raise RuntimeError(f"technical GSM mapping missing for {sample_id}")
                    matrix_path = (
                        Path(
                            "results/v6_2/phase2_5_data_onboarding/02_conversion_qc/h5ad"
                        )
                        / f"GSE291246.{gsm}.h5ad"
                    )
                    coordinate_path = None
                    scale_path = None
                    image_path = None
                    matrix_kind = "h5ad_counts"
                    count_key = "X"
                    read_authorization = "LOCAL_QC_XENIUM_COUNTS_COORDINATES"
                else:
                    raise RuntimeError(f"unsupported Stage 3 source {source_id}")
                if not matrix_path.is_file():
                    raise RuntimeError(f"materialization asset missing for {pilot_id}")
                source_record = sample_id
                patient_id = str(reg_row["patient_id"])
                block_id = str(reg_row["block_id"]) or _opaque_id("BLOCK_UNKNOWN", patient_id)
                leakage = str(reg_row["leakage_group_id"])
                section_opaque = _opaque_id("section", physical_id)
                capture_opaque = _opaque_id("capture", physical_id)
                resolved.append(
                    ResolvedSpatialUnit(
                        unit_id=_safe_id(physical_id),
                        pilot_unit_id=pilot_id,
                        source_id=source_id,
                        dataset_id=str(pilot_row["dataset_id"]),
                        physical_unit_id=physical_id,
                        opaque_patient_id=patient_id,
                        opaque_block_id=block_id,
                        opaque_section_id=section_opaque,
                        capture_id=capture_opaque,
                        platform=str(reg_row["platform"]),
                        modality=str(reg_row["modality"]),
                        resolution=str(reg_row["resolution"]),
                        counts_layer=count_key,
                        coordinate_system=str(reg_row["coordinate_system"]),
                        leakage_group_id=leakage,
                        selection_hash=str(pilot_row["selection_hash"]),
                        identity_status=str(reg_row["identity_confidence"]),
                        matrix_kind=matrix_kind,
                        matrix_path=matrix_path,
                        coordinate_path=coordinate_path,
                        scalefactors_path=scale_path,
                        image_path=image_path if image_path and image_path.is_file() else None,
                        source_record_id=source_record,
                        source_registry_row_hash=_source_row_hash(
                            reg_row,
                            (
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
                                "identity_confidence",
                                "leakage_group_id",
                                "record_status",
                                "n_observations",
                                "audit_status",
                            ),
                        ),
                        read_authorization=read_authorization,
                    )
                )
            continue

        section_opaque = _opaque_id("section", physical_id)
        capture_opaque = _opaque_id("capture", physical_id)
        resolved.append(
            ResolvedSpatialUnit(
                unit_id=_safe_id(physical_id),
                pilot_unit_id=pilot_id,
                source_id=source_id,
                dataset_id=str(pilot_row["dataset_id"]),
                physical_unit_id=physical_id,
                opaque_patient_id=patient_id,
                opaque_block_id=block_id,
                opaque_section_id=section_opaque,
                capture_id=capture_opaque,
                platform=str(reg_row["platform"]),
                modality=str(reg_row["modality"]),
                resolution=str(reg_row["resolution"]),
                counts_layer=count_key,
                coordinate_system=str(reg_row["coordinate_system"]),
                leakage_group_id=leakage,
                selection_hash=str(pilot_row["selection_hash"]),
                identity_status=str(reg_row["identity_confidence"]),
                matrix_kind=matrix_kind,
                matrix_path=matrix_path,
                coordinate_path=coordinate_path,
                scalefactors_path=scale_path,
                image_path=image_path,
                source_record_id=source_record,
                source_registry_row_hash=_source_row_hash(
                    reg_row,
                    (
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
                        "identity_confidence",
                        "leakage_group_id",
                        "record_status",
                        "n_observations",
                        "audit_status",
                    ),
                ),
                read_authorization=read_authorization,
            )
        )
    return sorted(resolved, key=lambda item: (item.source_id, item.physical_unit_id))


def _identity(unit: ResolvedSpatialUnit) -> SpatialIdentity:
    return SpatialIdentity(
        dataset_id=unit.dataset_id,
        opaque_patient_id=unit.opaque_patient_id,
        opaque_block_id=unit.opaque_block_id,
        opaque_section_id=unit.opaque_section_id,
        capture_id=unit.capture_id,
        leakage_group_id=unit.leakage_group_id,
        relationship_provenance="v7_stage3_technical_registry",
        relationship_confidence=unit.identity_status,
    )


def build_model_safe_manifest(
    units: Iterable[ResolvedSpatialUnit],
    *,
    config_hash: str,
) -> pd.DataFrame:
    """Build a locator-opaque manifest consumed by the discovery pipeline."""

    rows: list[dict[str, Any]] = []
    for unit in units:
        matrix_hash = _file_hash(unit.matrix_path)
        coordinate_hash = _file_hash(unit.coordinate_path) or matrix_hash
        scale_hash = _file_hash(unit.scalefactors_path)
        image_hash = _file_hash(unit.image_path)
        scale_status = (
            "pixel_scalefactors_available_no_micron_claim"
            if scale_hash
            else ("native_micrometer_coordinates" if unit.platform == "Xenium" else "not_provided")
        )
        row: dict[str, Any] = {
            "unit_id": unit.unit_id,
            "dataset_id": unit.dataset_id,
            "opaque_patient_id": unit.opaque_patient_id,
            "capture_id": unit.capture_id,
            "opaque_block_id": unit.opaque_block_id,
            "opaque_section_id": unit.opaque_section_id,
            "modality": "spatial",
            "platform": unit.platform,
            "matrix_locator": f"asset://{matrix_hash}",
            "matrix_locator_sha256": matrix_hash,
            "matrix_kind": unit.matrix_kind,
            "count_key": unit.counts_layer,
            "counts_semantics": "integer_raw_counts",
            "coordinate_locator": f"asset://{coordinate_hash}",
            "coordinate_locator_sha256": coordinate_hash,
            "coordinate_system": unit.coordinate_system,
            "read_authorization": unit.read_authorization,
            "leakage_group_id": unit.leakage_group_id,
            "selection_hash": unit.selection_hash,
            "scale_status": scale_status,
            "scalefactor_locator": f"asset://{scale_hash}" if scale_hash else "",
            "scalefactor_sha256": scale_hash,
            "image_locator": f"asset://{image_hash}" if image_hash else "",
            "image_sha256": image_hash,
            "image_usage": "locator_only_not_feature" if image_hash else "not_available",
            "segmentation_locator": (
                "native_observation_id" if unit.platform == "Xenium" else ""
            ),
            "segmentation_sha256": (
                stable_hash({"unit_id": unit.unit_id, "segmentation": "native_observation_id"})
                if unit.platform == "Xenium"
                else ""
            ),
            "tissue_mask_locator": "visium_positions_in_tissue" if unit.platform == "Visium" else "",
            "tissue_mask_sha256": coordinate_hash if unit.platform == "Visium" else "",
            "tissue_mask_source": (
                "non_gt_tissue_mask" if unit.platform == "Visium" else "not_provided"
            ),
            "source_registry_row_hash": unit.source_registry_row_hash,
            "qc_status": "PENDING",
        }
        row["provenance_hash"] = stable_hash(
            {"row": row, "config_hash": config_hash, "source": unit.source_record_id}
        )
        rows.append(row)
    frame = pd.DataFrame(rows)
    return validate_model_safe_manifest(frame)


def _nearest_scale(coords: np.ndarray) -> float:
    if len(coords) < 2:
        return float("nan")
    distances, _ = cKDTree(coords).query(coords, k=2)
    positive = np.asarray(distances[:, 1], dtype=float)
    positive = positive[np.isfinite(positive) & (positive > 0)]
    return float(np.median(positive)) if len(positive) else float("nan")


def _bootstrap_patient_ci(
    patient_frame: pd.DataFrame,
    *,
    group_columns: tuple[str, ...],
    seed: int,
    draws: int = 1000,
) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    keys = ["cohort_id", *group_columns]
    for key, frame in patient_frame.groupby(keys, dropna=False, sort=True):
        if not isinstance(key, tuple):
            key = (key,)
        valid = pd.to_numeric(frame.loc[frame["status"].eq("ESTIMABLE"), "estimate"], errors="coerce").dropna().to_numpy(float)
        row = dict(zip(keys, key))
        row["n_patients"] = int(len(valid))
        if len(valid) < 2:
            row.update({"ci_low": np.nan, "ci_high": np.nan, "status": "NOT_ESTIMABLE_TOO_FEW_PATIENTS"})
        else:
            rng = np.random.default_rng(seed ^ int(stable_hash({"key": list(key)})[:8], 16))
            sampled = rng.integers(0, len(valid), size=(draws, len(valid)))
            means = valid[sampled].mean(axis=1)
            row.update(
                {
                    "ci_low": float(np.quantile(means, 0.025)),
                    "ci_high": float(np.quantile(means, 0.975)),
                    "status": "ESTIMABLE",
                }
            )
        rows.append(row)
    return pd.DataFrame(rows)


def run_spatial_foundation(
    *,
    output_root: str | Path = "results/v7/spatial_foundation",
    stage1_registry_path: str | Path = "results/v7/registry/spatial_physical_units.tsv",
    r04_manifest_path: str | Path = "/home/huyudi/013_spatial/infra/r04/input_manifest.json",
    pilot_path: str | Path = "config/v7/stage3_pilot.tsv",
    module_dictionary_path: str | Path = "results/v7/ontology/module_dictionary.tsv",
    module_membership_path: str | Path = "results/v7/ontology/module_membership.tsv",
    measurement_reliability_path: str | Path = "results/v7/ontology/measurement_reliability.tsv",
    metadata_manifest_path: str | Path = (
        "results/v6_2/phase2_5_data_onboarding/01_metadata_freeze/"
        "spatial_sample_manifest.GSE291246.tsv"
    ),
    config_path: str | Path = "config/v7/stage3.yaml",
    permutations: int = 99,
    bootstrap_draws: int = 1000,
    seed: int = 20260904,
) -> dict[str, Any]:
    """Run the complete response-blind Stage 3 physical replay."""

    output = Path(output_root)
    output.mkdir(parents=True, exist_ok=True)
    scores_dir = output / "scores"
    scores_dir.mkdir(parents=True, exist_ok=True)

    registry = load_stage1_technical_registry(stage1_registry_path)
    r04, r04_provenance = load_r04_technical_manifest(r04_manifest_path)
    dynamic = build_stage3_pilot(stage1_registry_path, r04_manifest_path)
    static = validate_pilot_manifest(
        pd.read_csv(
            pilot_path,
            sep="\t",
            dtype=str,
            keep_default_na=False,
            comment="#",
        )
    )
    static_for_hash = static[list(PILOT_COLUMNS)].astype(str)
    dynamic_for_hash = dynamic.pilot[list(PILOT_COLUMNS)].astype(str)
    if stable_dataframe_hash(static_for_hash, PILOT_COLUMNS) != stable_dataframe_hash(
        dynamic_for_hash, PILOT_COLUMNS
    ):
        raise RuntimeError("static Stage 3 pilot differs from current technical selection")
    static.to_csv(output / "pilot_manifest.tsv", sep="\t", index=False)
    dynamic.technical_audit.to_csv(output / "technical_eligibility_audit.tsv", sep="\t", index=False)
    dynamic.selection_audit.to_csv(output / "pilot_selection_audit.tsv", sep="\t", index=False)
    (output / "pilot_provenance.json").write_text(
        json.dumps(
            {
                **dynamic.provenance,
                "pilot_file_sha256": sha256_file(pilot_path),
                "stage1_registry_file": file_provenance(stage1_registry_path),
                "r04_manifest_file": r04_provenance,
                "module_dictionary_sha256": sha256_file(module_dictionary_path),
                "module_membership_sha256": sha256_file(module_membership_path),
                "measurement_reliability_sha256": sha256_file(measurement_reliability_path),
                "config_sha256": sha256_file(config_path),
            },
            ensure_ascii=False,
            indent=2,
            sort_keys=True,
        ),
        encoding="utf-8",
    )
    config_hash = sha256_file(config_path)
    code_paths = [
        Path(__file__),
        Path(__file__).parent / "spatial_control.py",
        Path(__file__).parent / "spatial_projection.py",
        Path(__file__).parent / "spatial_stats/statistics.py",
        Path(__file__).parent / "spatial_stats/graph.py",
        Path(__file__).parent / "spatial_io/contracts.py",
        Path(__file__).parent / "spatial_io/sanitization.py",
    ]
    code_hashes = {
        str(path.relative_to(Path.cwd())): sha256_file(path)
        for path in code_paths
        if path.is_file()
    }
    try:
        git_commit = subprocess.run(
            ["git", "rev-parse", "HEAD"], capture_output=True, text=True, check=True
        ).stdout.strip()
        git_dirty = bool(
            subprocess.run(
                ["git", "status", "--porcelain"],
                capture_output=True,
                text=True,
                check=True,
            ).stdout.strip()
        )
    except (OSError, subprocess.CalledProcessError):
        git_commit, git_dirty = "unknown", True
    units = resolve_stage3_units(
        static,
        registry,
        r04,
        r04_manifest_path=r04_manifest_path,
        metadata_manifest=metadata_manifest_path,
    )
    if len(units) < 6 or len({unit.platform for unit in units}) < 2 or len({unit.source_id for unit in units}) < 3:
        raise RuntimeError("Stage 3 D3 minimum panel is not met")

    model_safe = build_model_safe_manifest(units, config_hash=config_hash)
    model_safe.to_csv(output / "model_safe_input_manifest.tsv", sep="\t", index=False)

    dictionary = pd.read_csv(module_dictionary_path, sep="\t")
    membership = pd.read_csv(module_membership_path, sep="\t")
    reliability = pd.read_csv(measurement_reliability_path, sep="\t")

    adapter_rows: list[dict[str, Any]] = []
    conservation_rows: list[dict[str, Any]] = []
    scale_rows: list[dict[str, Any]] = []
    identity_rows: list[dict[str, Any]] = []
    coverage_rows: list[pd.DataFrame] = []
    stat_rows: list[dict[str, Any]] = []
    null_rows: list[dict[str, Any]] = []
    variogram_rows: list[dict[str, Any]] = []
    composition_rows: list[dict[str, Any]] = []
    materialization_rows: list[dict[str, Any]] = []
    sanitization_rows: list[dict[str, Any]] = []
    successful_units = 0
    sanitization_dir = Path("/tmp/v7_stage3_sanitized")
    sanitization_dir.mkdir(parents=True, exist_ok=True)

    for ordinal, resolved in enumerate(units, start=1):
        identity = _identity(resolved)
        sanitizer_audit: dict[str, Any] = {
            "status": "NOT_APPLICABLE",
            "source_sha256": _file_hash(resolved.matrix_path),
        }
        load_path = resolved.matrix_path
        if resolved.platform == "Xenium":
            load_path = sanitization_dir / f"{resolved.unit_id}.h5ad"
            sanitizer_audit = sanitize_h5ad_for_discovery(resolved.matrix_path, load_path)
        removed_fields = [str(field) for field in sanitizer_audit.get("removed_fields", [])]
        sanitization_rows.append(
            {
                "unit_id": resolved.unit_id,
                "source_sha256": sanitizer_audit.get("source_sha256", ""),
                "sanitized_sha256": sanitizer_audit.get("destination_sha256", ""),
                "removed_field_count": len(removed_fields),
                "removed_field_names_hash": stable_hash({"fields": removed_fields}),
                "status": sanitizer_audit.get("status", "NOT_APPLICABLE"),
            }
        )
        if resolved.platform == "Visium" and resolved.matrix_kind == "h5ad_counts":
            spatial_unit = load_h5ad_counts(
                load_path,
                identity,
                counts_layer="X",
                coordinate_key="spatial",
                coordinate_unit="array_unit",
                platform="Visium",
                resolution="spot",
            )
        elif resolved.platform == "Visium":
            spatial_unit = load_visium_10x(
                load_path,
                resolved.coordinate_path,
                identity,
                scalefactors_path=resolved.scalefactors_path,
            )
        elif resolved.platform == "Xenium":
            spatial_unit = load_xenium_h5ad(
                load_path,
                identity,
                counts_layer="X",
                coordinate_key="spatial",
            )
        else:
            raise RuntimeError(f"unsupported adapter platform: {resolved.platform}")

        successful_units += 1
        observation_ids = spatial_unit.observations["observation_id"].astype(str).to_numpy()
        coordinates = spatial_unit.observations[["analysis_x", "analysis_y"]].to_numpy(float)
        sections = spatial_unit.observations["section_id"].astype(str).to_numpy()
        graph = build_section_graph(coordinates, sections, method="knn", k=6, weight_mode="binary")
        nearest_scale = _nearest_scale(coordinates)
        scaled_coordinates = coordinates / nearest_scale if np.isfinite(nearest_scale) and nearest_scale > 0 else coordinates
        scores, coverage = project_stage2_features(
            spatial_unit.counts,
            spatial_unit.features["gene_symbol"],
            dictionary,
            membership,
            section_ids=sections,
            measurement_reliability=reliability,
        )
        scores.insert(0, "observation_id", observation_ids)
        scores.insert(1, "capture_id", resolved.capture_id)
        scores.insert(2, "section_id", resolved.opaque_section_id)
        scores.insert(3, "opaque_patient_id", resolved.opaque_patient_id)
        scores.insert(4, "dataset_id", resolved.dataset_id)
        scores.insert(5, "platform", resolved.platform)
        scores.insert(6, "unit_id", resolved.unit_id)
        scores.to_parquet(scores_dir / f"{resolved.unit_id}.parquet", index=False)
        coverage = coverage.assign(
            unit_id=resolved.unit_id,
            dataset_id=resolved.dataset_id,
            platform=resolved.platform,
            opaque_patient_id=resolved.opaque_patient_id,
            section_id=resolved.opaque_section_id,
        )
        coverage_rows.append(coverage)

        audit = spatial_unit.audit
        adapter_rows.append(
            {
                "unit_id": resolved.unit_id,
                "dataset_id": resolved.dataset_id,
                "platform": resolved.platform,
                "resolution": resolved.resolution,
                "n_observations": spatial_unit.n_observations,
                "n_features": spatial_unit.n_features,
                "nnz": int(spatial_unit.counts.nnz),
                "count_total": spatial_unit.count_total,
                "counts_semantics": spatial_unit.counts_semantics,
                "coordinate_system": spatial_unit.coordinate_system,
                "adapter_status": audit.get("adapter_status", "PASS"),
                "position_join_status": audit.get("position_join_status", "not_applicable"),
                "scale_status": audit.get("scale_status", "not_provided"),
                "image_status": audit.get("image_status", "not_loaded"),
                "segmentation_status": audit.get("segmentation_status", "not_provided"),
                "sanitization_status": sanitizer_audit.get("status", "NOT_APPLICABLE"),
                "sanitization_removed_fields": len(sanitizer_audit.get("removed_fields", [])),
                "qc_status": "PASS",
            }
        )
        conservation_rows.append(
            {
                "unit_id": resolved.unit_id,
                "count_total": spatial_unit.count_total,
                "nnz": int(spatial_unit.counts.nnz),
                "n_observations": spatial_unit.n_observations,
                "n_features": spatial_unit.n_features,
                "row_total_finite": bool(np.isfinite(np.asarray(spatial_unit.counts.sum(axis=1))).all()),
                "feature_total_finite": bool(np.isfinite(np.asarray(spatial_unit.counts.sum(axis=0))).all()),
                "status": "PASS",
            }
        )
        scale_rows.append(
            {
                "unit_id": resolved.unit_id,
                "platform": resolved.platform,
                "coordinate_system": spatial_unit.coordinate_system,
                "native_coordinate_unit": str(spatial_unit.observations["native_coordinate_unit"].iloc[0]),
                "analysis_coordinate_unit": str(spatial_unit.observations["analysis_coordinate_unit"].iloc[0]),
                "nearest_neighbor_median_native": nearest_scale,
                "distance_scale_used": "median_nearest_neighbor",
                "physical_distance_status": audit.get("scale_status", "not_provided"),
            }
        )
        identity_rows.append(
            {
                "unit_id": resolved.unit_id,
                "opaque_patient_id": resolved.opaque_patient_id,
                "opaque_block_id": resolved.opaque_block_id,
                "opaque_section_id": resolved.opaque_section_id,
                "capture_id": resolved.capture_id,
                "leakage_group_id": resolved.leakage_group_id,
                "relationship_provenance": identity.relationship_provenance,
                "relationship_confidence": identity.relationship_confidence,
                "n_observations": spatial_unit.n_observations,
                "observation_id_unique": not spatial_unit.observations["observation_id"].duplicated().any(),
                "status": "PASS",
            }
        )

        feature_columns = [column for column in scores.columns if column.endswith("__score_standardized")]
        for feature_column in feature_columns:
            feature_id = feature_column.removesuffix("__score_standardized")
            values = scores[feature_column].to_numpy(float)
            observed = moran_statistic(values, graph, section_ids=sections)
            stat_rows.append(
                {
                    "unit_id": resolved.unit_id,
                    "dataset_id": resolved.dataset_id,
                    "platform": resolved.platform,
                    "patient_id": resolved.opaque_patient_id,
                    "section_id": resolved.opaque_section_id,
                    "cohort_id": resolved.dataset_id,
                    "feature_id": feature_id,
                    "statistic": "moran",
                    "status": observed["status"],
                    "estimate": observed["estimate"],
                    "n_observations": observed["n_observations"],
                    "n_directed_edges": observed["n_directed_edges"],
                    "distance_scale": nearest_scale,
                }
            )
            permutation = moran_permutation_test(
                values,
                graph,
                sections,
                permutations=permutations,
                seed=seed ^ int(stable_hash({"unit": resolved.unit_id, "feature": feature_id})[:8], 16),
            )
            null = np.asarray(permutation.get("permuted_statistics", []), dtype=float)
            null_rows.append(
                {
                    "unit_id": resolved.unit_id,
                    "dataset_id": resolved.dataset_id,
                    "platform": resolved.platform,
                    "patient_id": resolved.opaque_patient_id,
                    "section_id": resolved.opaque_section_id,
                    "cohort_id": resolved.dataset_id,
                    "feature_id": feature_id,
                    "statistic": "moran",
                    "status": permutation["status"],
                    "observed": permutation["observed"],
                    "p_value": permutation["p_value"],
                    "n_permutations": permutation["n_permutations"],
                    "null_median": float(np.median(null)) if len(null) else np.nan,
                    "null_q025": float(np.quantile(null, 0.025)) if len(null) else np.nan,
                    "null_q975": float(np.quantile(null, 0.975)) if len(null) else np.nan,
                }
            )
            variogram = sparse_variogram(
                values,
                graph,
                scaled_coordinates,
                section_ids=sections,
                distance_bins=8,
            )
            for record in variogram.to_dict("records"):
                variogram_rows.append(
                    {
                        "unit_id": resolved.unit_id,
                        "dataset_id": resolved.dataset_id,
                        "platform": resolved.platform,
                        "patient_id": resolved.opaque_patient_id,
                        "section_id": resolved.opaque_section_id,
                        "cohort_id": resolved.dataset_id,
                        "feature_id": feature_id,
                        "statistic": "variogram",
                        "distance_scale": nearest_scale,
                        **record,
                    }
                )
        composition_rows.append(
            {
                "unit_id": resolved.unit_id,
                "dataset_id": resolved.dataset_id,
                "platform": resolved.platform,
                "patient_id": resolved.opaque_patient_id,
                "section_id": resolved.opaque_section_id,
                "cohort_id": resolved.dataset_id,
                "source": "T",
                "target": "Myeloid",
                **edge_association(
                    graph,
                    None,
                    source="T",
                    target="Myeloid",
                    section_ids=sections,
                ),
            }
        )
        materialization_rows.append(
            {
                "unit_id": resolved.unit_id,
                "source_id": resolved.source_id,
                "source_record_id_hash": stable_hash(
                    {"source_record_id": resolved.source_record_id}
                ),
                "physical_unit_id_hash": stable_hash(
                    {"physical_unit_id": resolved.physical_unit_id}
                ),
                "matrix_locator": f"asset://{_file_hash(resolved.matrix_path)}",
                "matrix_sha256": _file_hash(resolved.matrix_path),
                "coordinate_locator": (
                    f"asset://{_file_hash(resolved.coordinate_path)}"
                    if resolved.coordinate_path
                    else "in_matrix_obsm/spatial"
                ),
                "coordinate_sha256": _file_hash(resolved.coordinate_path) or _file_hash(resolved.matrix_path),
                "status": "PASS",
            }
        )

    if successful_units != len(units):
        raise RuntimeError("one or more selected Stage 3 units failed core replay")

    coverage_frame = pd.concat(coverage_rows, ignore_index=True)
    coverage_frame.to_csv(output / "feature_coverage.tsv", sep="\t", index=False)
    projection_summary = (
        coverage_frame.groupby(
            ["dataset_id", "platform", "feature_id", "spatial_projection_status"],
            dropna=False,
            sort=True,
        )
        .agg(
            n_units=("unit_id", "nunique"),
            median_coverage=("coverage", "median"),
            min_coverage=("coverage", "min"),
            max_coverage=("coverage", "max"),
            n_finite_score_rows=("n_finite_scores", "sum"),
        )
        .reset_index()
    )
    projection_summary.to_csv(
        output / "spatial_projection_summary.tsv", sep="\t", index=False
    )
    adapter_frame = pd.DataFrame(adapter_rows)
    adapter_frame.to_csv(output / "adapter_qc.tsv", sep="\t", index=False)
    pd.DataFrame(conservation_rows).to_csv(output / "count_conservation.tsv", sep="\t", index=False)
    pd.DataFrame(scale_rows).to_csv(output / "coordinate_scale_audit.tsv", sep="\t", index=False)
    identity_frame = pd.DataFrame(identity_rows)
    identity_frame.to_csv(output / "identity_leakage_audit.tsv", sep="\t", index=False)
    pd.DataFrame(materialization_rows).to_csv(output / "source_asset_audit.tsv", sep="\t", index=False)
    stats_frame = pd.DataFrame(stat_rows)
    null_frame = pd.DataFrame(null_rows)
    variogram_frame = pd.DataFrame(variogram_rows)
    composition_frame = pd.DataFrame(composition_rows)
    stats_frame.to_parquet(output / "graph_statistics.parquet", index=False)
    null_frame.to_parquet(output / "spatial_null_statistics.parquet", index=False)
    variogram_frame.to_parquet(output / "spatial_variogram.parquet", index=False)
    composition_frame.to_csv(output / "composition_edge_status.tsv", sep="\t", index=False)
    pd.DataFrame(sanitization_rows).to_csv(
        output / "sanitization_audit.tsv", sep="\t", index=False
    )

    patient_section, cohort = patient_nested_summary(
        stats_frame,
        value_column="estimate",
        status_column="status",
        patient_column="patient_id",
        section_column="section_id",
        cohort_column="cohort_id",
        group_columns=("feature_id", "statistic"),
    )
    patient_section["level"] = "patient"
    cohort["level"] = "cohort"
    patient_summary = pd.concat([patient_section, cohort], ignore_index=True, sort=False)
    patient_summary.to_parquet(output / "patient_level_summary.parquet", index=False)
    bootstrap = _bootstrap_patient_ci(
        patient_section,
        group_columns=("feature_id", "statistic"),
        seed=seed,
        draws=bootstrap_draws,
    )
    bootstrap.to_csv(output / "patient_bootstrap_ci.tsv", sep="\t", index=False)

    split_source = identity_frame[
        ["unit_id", "opaque_patient_id", "leakage_group_id"]
    ].drop_duplicates("unit_id").rename(columns={"opaque_patient_id": "patient_id"})
    split_source["dataset_id"] = [
        next(unit.dataset_id for unit in units if unit.unit_id == unit_id)
        for unit_id in split_source["unit_id"]
    ]
    split = patient_block_split(
        split_source,
        patient_column="patient_id",
        n_splits=3,
        seed=seed,
        stratify_columns=("dataset_id",),
        fold_column="fold",
    )
    split.to_csv(output / "split_audit.tsv", sep="\t", index=False)

    replay_lines = [
        "schema: v7.stage3.platform_replay_manifest.v1",
        "status: COMPLETE",
        "response_blind: true",
        "units:",
    ]
    for unit in units:
        replay_lines.extend(
            [
                f"  - unit_id: {unit.unit_id}",
                f"    source_id: {unit.source_id}",
                f"    platform: {unit.platform}",
                f"    matrix_kind: {unit.matrix_kind}",
                f"    matrix_sha256: {_file_hash(unit.matrix_path)}",
                f"    coordinate_sha256: {_file_hash(unit.coordinate_path) or _file_hash(unit.matrix_path)}",
                "    image_usage: locator_only_not_feature",
                "    replay_status: PASS",
            ]
        )
    (output / "platform_replay_manifest.yaml").write_text(
        "\n".join(replay_lines) + "\n", encoding="utf-8"
    )

    post_freeze_files = {
        "DISCOVERY_OUTPUT_HASH.txt",
        "STAGE3_RUN_MANIFEST.yaml",
        "D3_GATE.json",
        "SPATIAL_FOUNDATION_REPORT.md",
        "sealed_validation_manifest.tsv",
        "gt_isolation_audit.tsv",
    }
    discovery_files = sorted(
        path
        for path in output.iterdir()
        if path.is_file() and path.name not in post_freeze_files
    )
    discovery_hash = stable_hash(
        {
            "files": [
                {"name": path.name, "sha256": sha256_file(path), "bytes": path.stat().st_size}
                for path in discovery_files
            ]
        }
    )
    (output / "DISCOVERY_OUTPUT_HASH.txt").write_text(discovery_hash + "\n", encoding="utf-8")
    run_manifest = {
        "schema": "v7.stage3.spatial_foundation.run.v1",
        "status": "S3_REPLAY_COMPLETE_PENDING_GT_ISOLATION",
        "n_logical_pilot_rows": int(len(static)),
        "n_physical_captures": int(len(units)),
        "n_patients": int(identity_frame["opaque_patient_id"].nunique()),
        "platforms": sorted(identity_frame["unit_id"].map({u.unit_id: u.platform for u in units}).dropna().unique()),
        "datasets": sorted(identity_frame["unit_id"].map({u.unit_id: u.dataset_id for u in units}).dropna().unique()),
        "permutations": permutations,
        "bootstrap_draws": bootstrap_draws,
        "seed": seed,
        "discovery_output_hash": discovery_hash,
        "model_safe_manifest_hash": sha256_file(output / "model_safe_input_manifest.tsv"),
        "outputs": [path.name for path in discovery_files],
        "response_blind": True,
        "sealed_gt_read": False,
        "git_commit": git_commit,
        "git_dirty": git_dirty,
        "code_sha256": code_hashes,
    }
    (output / "STAGE3_RUN_MANIFEST.yaml").write_text(
        "\n".join(
            [
                f"schema: {run_manifest['schema']}",
                f"status: {run_manifest['status']}",
                f"n_logical_pilot_rows: {run_manifest['n_logical_pilot_rows']}",
                f"n_physical_captures: {run_manifest['n_physical_captures']}",
                f"n_patients: {run_manifest['n_patients']}",
                "platforms: [" + ", ".join(run_manifest["platforms"]) + "]",
                "datasets: [" + ", ".join(run_manifest["datasets"]) + "]",
                f"permutations: {permutations}",
                f"bootstrap_draws: {bootstrap_draws}",
                f"seed: {seed}",
                f"discovery_output_hash: {discovery_hash}",
                f"model_safe_manifest_hash: {run_manifest['model_safe_manifest_hash']}",
                "response_blind: true",
                "sealed_gt_read: false",
                f"git_commit: {git_commit}",
                f"git_dirty: {str(git_dirty).lower()}",
                "code_sha256:",
                *[f"  {name}: {digest}" for name, digest in sorted(code_hashes.items())],
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    return {
        **run_manifest,
        "output_root": str(output.resolve()),
        "adapter_qc": adapter_frame,
        "feature_coverage": coverage_frame,
        "stats": stats_frame,
        "null": null_frame,
    }


__all__ = ["ResolvedSpatialUnit", "resolve_stage3_units", "build_model_safe_manifest", "run_spatial_foundation"]
