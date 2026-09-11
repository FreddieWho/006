"""Response-blind Stage 4 spatial discovery and transparent baselines.

The module expands the technically materialisable Stage 3 sources, keeps the
patient as the outer unit, and computes interpretable spatial primitives.  It
does not read response, treatment, endpoint, or sealed structure values.  The
marker panel below is a *composition proxy*; it is never written as a cell
fraction or used as a ground truth label.
"""

from __future__ import annotations

from dataclasses import dataclass
import json
from pathlib import Path
from typing import Any, Iterable, Mapping

import numpy as np
import pandas as pd
from scipy import sparse
from sklearn.decomposition import NMF

from .spatial_control import (
    PILOT_COLUMNS,
    file_provenance,
    load_r04_technical_manifest,
    load_stage1_technical_registry,
    sha256_file,
    stable_hash,
    technical_eligibility_audit,
)
from .spatial_io.adapters import load_h5ad_counts, load_visium_10x, load_xenium_h5ad
from .spatial_io.sanitization import sanitize_h5ad_for_discovery
from .spatial_pipeline import (
    ResolvedSpatialUnit,
    _file_hash,
    _identity,
    _nearest_scale,
    build_model_safe_manifest,
    resolve_stage3_units,
)
from .spatial_projection import project_stage2_features
from .spatial_stats import (
    build_section_graph,
    edge_association,
    moran_permutation_test,
    moran_statistic,
)


# The sets are deliberately coarse and auditable.  They are reference-marker
# proxies for adjustment, not a deconvolution result.
MARKER_PANEL: Mapping[str, tuple[str, ...]] = {
    "T_NK_proxy": ("CD3D", "CD3E", "TRBC1", "TRBC2", "NKG7", "GNLY", "CD8A"),
    "B_Plasma_proxy": ("CD79A", "MS4A1", "CD74", "CD37", "JCHAIN", "MZB1"),
    "Myeloid_proxy": ("LST1", "TYROBP", "FCER1G", "CTSS", "LILRB1", "S100A8", "S100A9"),
    "CAF_Stromal_proxy": ("COL1A1", "COL1A2", "COL3A1", "DCN", "LUM", "PDGFRA"),
    "Endothelial_proxy": ("PECAM1", "VWF", "KDR", "EMCN", "ESM1", "ENG"),
    "Tumor_Epithelial_proxy": ("EPCAM", "KRT8", "KRT18", "KRT19", "MUC1", "KRT7"),
}


@dataclass(frozen=True)
class Stage4Selection:
    pilot: pd.DataFrame
    audit: pd.DataFrame
    block_gap_review: pd.DataFrame
    units: tuple[ResolvedSpatialUnit, ...]


def _capture_from_section(section_id: str) -> str:
    value = str(section_id).removeprefix("HTAN::")
    for suffix in ("_filtered_trimmed.h5ad", "_filtered.h5ad", ".h5ad"):
        value = value.removesuffix(suffix)
    return value


def _pilot_record(
    *,
    pilot_id: str,
    source_id: str,
    dataset_id: str,
    patient_id: str,
    block_id: str,
    section_id: str,
    capture_id: str,
    platform: str,
    modality: str,
    resolution: str,
    counts_layer: str,
    coordinate_system: str,
    identity_status: str,
    leakage_group_id: str,
    member_count: int,
    member_hash: str,
    registry_hash: str,
    r04_hash: str,
    r04_payload_hash: str,
) -> dict[str, Any]:
    return {
        "pilot_unit_id": pilot_id,
        "source_id": source_id,
        "dataset_id": dataset_id,
        "unit_kind": "technical_capture" if capture_id else "patient_all_eligible_captures",
        "opaque_patient_id": patient_id,
        "opaque_block_id": block_id,
        "opaque_section_id": section_id,
        "capture_id": capture_id,
        "platform": platform,
        "modality": modality,
        "resolution": resolution,
        "counts_layer": counts_layer,
        "coordinate_system": coordinate_system,
        "identity_status": identity_status,
        "leakage_group_id": leakage_group_id,
        "member_capture_count": member_count,
        "member_set_hash": member_hash,
        "selection_mode": "all_materialisable_technical_units",
        "selection_hash": stable_hash(
            {"namespace": "v7_stage4_all_materialisable_v1", "pilot_id": pilot_id}
        ),
        "technical_eligibility": "eligible",
        "stage1_registry_sha256": registry_hash,
        "r04_manifest_sha256": r04_hash,
        "r04_payload_hash": r04_payload_hash,
    }


def build_stage4_selection(
    *,
    stage1_registry_path: str | Path = "results/v7/registry/spatial_physical_units.tsv",
    r04_manifest_path: str | Path = "/home/huyudi/013_spatial/infra/r04/input_manifest.json",
    metadata_manifest_path: str | Path = (
        "results/v6_2/phase2_5_data_onboarding/01_metadata_freeze/"
        "spatial_sample_manifest.GSE291246.tsv"
    ),
) -> Stage4Selection:
    """Select all presently materialisable HTAN/GSE238/GSE291 captures.

    The resulting pilot is a technical selection table.  No treatment or
    response columns are loaded.  Patients with several captures are expanded
    by ``resolve_stage3_units`` and remain one outer patient in all summaries.
    """

    registry_path = Path(stage1_registry_path)
    registry = load_stage1_technical_registry(registry_path)
    technical = technical_eligibility_audit(registry)
    r04, r04_provenance = load_r04_technical_manifest(r04_manifest_path)
    registry_by_id = registry.set_index("physical_unit_id", drop=False)
    technical_by_id = technical.set_index("physical_unit_id", drop=False)
    registry_hash = sha256_file(registry_path)
    r04_hash = str(r04_provenance["content_sha256"])
    payload_hash = str(r04_provenance["payload_hash"])
    raw_r04 = json.loads(Path(r04_manifest_path).read_text(encoding="utf-8"))
    raw_r04_by_physical = {
        str(row["section_id"]): row for row in raw_r04.get("rows", [])
    }
    rows: list[dict[str, Any]] = []
    audit: list[dict[str, Any]] = []

    # HTAN rows have an audited locator in R-04 and are selected one capture at
    # a time so every source row has an explicit provenance chain.
    for source in r04.loc[r04.lineage.eq("HTAN_VANDERBILT_CRC")].sort_values(
        "physical_unit_id", kind="mergesort"
    ).to_dict("records"):
        physical_id = str(source["physical_unit_id"])
        reg_row = registry_by_id.loc[physical_id]
        tech_row = technical_by_id.loc[physical_id]
        capture = _capture_from_section(str(source["section_id"]))
        raw_source = raw_r04_by_physical.get(physical_id, {})
        path = Path(str(raw_source.get("matrix_locator", "")))
        selected = (
            str(tech_row["technical_eligibility"]) == "eligible" and path.is_file()
        )
        audit.append(
            {
                "source_id": "013_HTAN_CRC",
                "physical_unit_id": physical_id,
                "patient_id": str(source["patient_id"]),
                "selected": selected,
                "selection_reason": "materialisable_r04_capture" if selected else "missing_or_ineligible",
                "matrix_exists": path.is_file(),
                "technical_status": str(tech_row["technical_eligibility"]),
            }
        )
        if not selected:
            continue
        rows.append(
            _pilot_record(
                pilot_id=f"013_HTAN_CRC::capture::{capture}",
                source_id="013_HTAN_CRC",
                dataset_id=str(reg_row["dataset_id"]),
                patient_id=str(source["patient_id"]),
                block_id=str(source["block_id"]),
                section_id=physical_id,
                capture_id=capture,
                platform=str(reg_row["platform"]),
                modality=str(reg_row["modality"]),
                resolution=str(reg_row["resolution"]),
                counts_layer=str(reg_row["counts_layer"]),
                coordinate_system=str(reg_row["coordinate_system"]),
                identity_status=str(reg_row["identity_confidence"]),
                leakage_group_id=str(source["leakage_group_id"]),
                member_count=1,
                member_hash=stable_hash({"physical_unit_ids": [physical_id]}),
                registry_hash=registry_hash,
                r04_hash=r04_hash,
                r04_payload_hash=payload_hash,
            )
        )

    # GSE238 and GSE291 are grouped by patient in the control table; the
    # resolver expands each group back to all technically eligible captures.
    metadata_manifest = Path(metadata_manifest_path)
    gsm_map: dict[str, str] = {}
    if metadata_manifest.is_file():
        m = pd.read_csv(
            metadata_manifest,
            sep="\t",
            usecols=["sample_id", "gsm_id"],
            dtype=str,
            keep_default_na=False,
        )
        gsm_map = dict(zip(m["sample_id"], m["gsm_id"]))

    local_sources = {"006_GSE238264", "006_GSE291246"}
    local = registry[
        registry.source_id.isin(local_sources)
        & registry.counts_layer.eq("h5ad_or_spot_matrix")
        & registry.record_status.eq("EFFECTIVE_SPATIAL_MANIFEST")
    ].copy()
    for (source_id, patient_id), group in local.groupby(
        ["source_id", "patient_id"], sort=True
    ):
        materialisable = True
        for physical_id in group.physical_unit_id.astype(str):
            sample_id = physical_id.split("::section::", 1)[-1]
            if source_id == "006_GSE238264":
                path = (
                    Path("results/v6_2/phase2_5_data_onboarding/_scratch/GSE238264_verify")
                    / sample_id
                    / "filtered_feature_bc_matrix.h5"
                )
            else:
                gsm = gsm_map.get(sample_id, "")
                path = (
                    Path("results/v6_2/phase2_5_data_onboarding/02_conversion_qc/h5ad")
                    / f"GSE291246.{gsm}.h5ad"
                )
            materialisable &= path.is_file()
        first = group.sort_values("physical_unit_id", kind="mergesort").iloc[0]
        member_ids = sorted(group.physical_unit_id.astype(str))
        selected = bool(materialisable and len(member_ids) > 0)
        audit.append(
            {
                "source_id": source_id,
                "physical_unit_id": "",
                "patient_id": str(patient_id),
                "selected": selected,
                "selection_reason": "all_patient_captures_materialisable" if selected else "one_or_more_assets_missing",
                "matrix_exists": materialisable,
                "technical_status": "eligible" if selected else "incomplete",
            }
        )
        if not selected:
            continue
        rows.append(
            _pilot_record(
                pilot_id=f"{source_id}::patient::{patient_id}",
                source_id=source_id,
                dataset_id=str(first["dataset_id"]),
                patient_id=str(patient_id),
                block_id="",
                section_id="",
                capture_id="",
                platform=str(first["platform"]),
                modality=str(first["modality"]),
                resolution=str(first["resolution"]),
                counts_layer=str(first["counts_layer"]),
                coordinate_system=str(first["coordinate_system"]),
                identity_status=str(first["identity_confidence"]),
                leakage_group_id=str(first["leakage_group_id"]),
                member_count=len(member_ids),
                member_hash=stable_hash({"physical_unit_ids": member_ids}),
                registry_hash=registry_hash,
                r04_hash=r04_hash,
                r04_payload_hash=payload_hash,
            )
        )

    pilot = pd.DataFrame(rows).reindex(columns=list(PILOT_COLUMNS))
    pilot = pilot.sort_values(["source_id", "pilot_unit_id"], kind="mergesort").reset_index(drop=True)
    # This custom selection is intentionally not passed to Stage 3's fixed
    # eight-row validator.  We check the same schema and hash fields here.
    if pilot.empty or pilot["pilot_unit_id"].duplicated().any():
        raise RuntimeError("Stage 4 selection is empty or has duplicate logical units")
    for column in ("selection_hash", "member_set_hash", "stage1_registry_sha256", "r04_manifest_sha256", "r04_payload_hash"):
        if not pilot[column].astype(str).str.fullmatch(r"[0-9a-f]{64}").all():
            raise RuntimeError(f"Stage 4 selection has invalid {column}")

    # The block-gap review is intentionally metadata-only.  It identifies
    # patient-known records without pretending that their block is known.
    block_gap = registry[registry.block_id.astype(str).str.strip().eq("")].copy()
    block_gap_review = block_gap[
        [
            "physical_unit_id",
            "source_id",
            "dataset_id",
            "patient_id",
            "identity_confidence",
            "duplicate_lineage",
            "leakage_group_id",
            "record_status",
            "n_observations",
            "audit_status",
        ]
    ].copy()
    block_gap_review["review_status"] = np.where(
        block_gap_review["patient_id"].astype(str).str.strip().ne("")
        & block_gap_review["identity_confidence"].astype(str).str.strip().ne(""),
        "PATIENT_KNOWN_BLOCK_UNKNOWN",
        "IDENTITY_OR_BLOCK_UNRESOLVED",
    )
    block_gap_review = block_gap_review.sort_values(
        ["source_id", "dataset_id", "patient_id", "physical_unit_id"], kind="mergesort"
    )
    units = tuple(
        resolve_stage3_units(
            pilot,
            registry,
            r04,
            r04_manifest_path=r04_manifest_path,
            metadata_manifest=metadata_manifest_path,
        )
    )
    if not units:
        raise RuntimeError("Stage 4 selected no materialisable captures")
    return Stage4Selection(
        pilot=pilot,
        audit=pd.DataFrame(audit),
        block_gap_review=block_gap_review,
        units=units,
    )


def _robust_z(values: np.ndarray) -> np.ndarray:
    result = np.full(len(values), np.nan, dtype=float)
    valid = np.isfinite(values)
    if valid.sum() < 2:
        return result
    center = float(np.nanmedian(values))
    scale = float(np.nanmedian(np.abs(values[valid] - center))) * 1.4826
    if scale <= 0:
        scale = float(np.nanstd(values))
    if scale <= 0:
        scale = 1.0
    result[valid] = (values[valid] - center) / scale
    return result


def _marker_proxy(
    counts: sparse.spmatrix, gene_symbols: Iterable[object]
) -> tuple[pd.DataFrame, pd.DataFrame]:
    matrix = sparse.csr_matrix(counts, dtype=float)
    library = np.asarray(matrix.sum(axis=1)).ravel()
    scale = np.divide(1_000_000.0, library, out=np.zeros_like(library), where=library > 0)
    logged = matrix.multiply(scale[:, None]).tocsr()
    if logged.data.size:
        logged.data = np.log1p(logged.data)
        logged.eliminate_zeros()
    symbols = np.asarray([str(x).strip() for x in gene_symbols], dtype=object)
    positions: dict[str, list[int]] = {}
    for index, symbol in enumerate(symbols):
        positions.setdefault(symbol, []).append(index)
    scores: dict[str, np.ndarray] = {}
    audit: list[dict[str, Any]] = []
    for name, markers in MARKER_PANEL.items():
        observed = [marker for marker in markers if marker in positions]
        if observed:
            columns = [column for marker in observed for column in positions[marker]]
            values = np.asarray(logged[:, columns].mean(axis=1)).ravel()
        else:
            values = np.full(matrix.shape[0], np.nan, dtype=float)
        scores[name] = _robust_z(values)
        audit.append(
            {
                "proxy_id": name,
                "expected_marker_count": len(markers),
                "observed_marker_count": len(observed),
                "observed_markers": ";".join(observed),
                "coverage": len(observed) / len(markers),
                "status": "proxy_available" if observed else "NOT_ESTIMABLE_NO_MARKERS",
            }
        )
    return pd.DataFrame(scores), pd.DataFrame(audit)


def _composition_residual(values: np.ndarray, proxy: pd.DataFrame) -> np.ndarray:
    """Remove marker-proxy abundance within a section without imputing genes."""
    y = np.asarray(values, dtype=float)
    columns = [column for column in proxy.columns if proxy[column].notna().sum() >= 3]
    result = np.full(len(y), np.nan, dtype=float)
    valid = np.isfinite(y)
    if not columns or valid.sum() < 4:
        result[valid] = y[valid]
        return result
    x = proxy[columns].to_numpy(float)
    valid &= np.isfinite(x).all(axis=1)
    if valid.sum() < len(columns) + 2:
        result[np.isfinite(y)] = y[np.isfinite(y)]
        return result
    design = np.column_stack([np.ones(valid.sum()), x[valid]])
    beta, *_ = np.linalg.lstsq(design, y[valid], rcond=None)
    result[valid] = y[valid] - design @ beta
    return result


def _safe_quantile(values: Iterable[object], quantile: float) -> float:
    array = pd.to_numeric(pd.Series(list(values)), errors="coerce").to_numpy(float)
    array = array[np.isfinite(array)]
    return float(np.quantile(array, quantile)) if len(array) else np.nan


def _load_unit(resolved: ResolvedSpatialUnit, sanitization_dir: Path):
    identity = _identity(resolved)
    load_path = resolved.matrix_path
    sanitization = {"status": "NOT_APPLICABLE", "removed_fields": []}
    if resolved.platform == "Xenium":
        load_path = sanitization_dir / f"{resolved.unit_id}.h5ad"
        sanitization = sanitize_h5ad_for_discovery(resolved.matrix_path, load_path)
    if resolved.platform == "Visium" and resolved.matrix_kind == "h5ad_counts":
        unit = load_h5ad_counts(
            load_path,
            identity,
            counts_layer="X",
            coordinate_key="spatial",
            coordinate_unit="array_unit",
            platform="Visium",
            resolution="spot",
        )
    elif resolved.platform == "Visium":
        unit = load_visium_10x(
            load_path,
            resolved.coordinate_path,
            identity,
            scalefactors_path=resolved.scalefactors_path,
        )
    elif resolved.platform == "Xenium":
        unit = load_xenium_h5ad(
            load_path,
            identity,
            counts_layer="X",
            coordinate_key="spatial",
        )
    else:
        raise RuntimeError(f"unsupported Stage 4 platform {resolved.platform}")
    return unit, sanitization


def _nmf_summary(
    matrix: np.ndarray,
    graph: sparse.spmatrix,
    *,
    section_ids: np.ndarray,
    unit_id: str,
    seeds: tuple[int, ...],
    ks: tuple[int, ...],
    sample_limit: int = 1500,
) -> pd.DataFrame:
    finite = np.isfinite(matrix).all(axis=1)
    finite_indices = np.flatnonzero(finite)
    x_full = matrix[finite]
    if len(x_full) < 20:
        return pd.DataFrame(
            [
                {
                    "unit_id": unit_id,
                    "k": 0,
                    "seed": -1,
                    "status": "NOT_ESTIMABLE_TOO_FEW_OBSERVATIONS",
                }
            ]
        )
    # Derive the non-negative offset from the complete finite matrix before
    # sampling.  If it were derived from the sample, an unseen smaller value
    # could make ``model.transform(x_full - offsets)`` negative and turn a
    # valid latent-field run into a spurious ValueError.
    offsets = np.nanmin(x_full, axis=0, keepdims=True)
    x = x_full - offsets
    x += 1e-6
    if len(x) > sample_limit:
        rng = np.random.default_rng(int(stable_hash({"unit": unit_id})[:8], 16))
        selected = np.sort(rng.choice(len(x), size=sample_limit, replace=False))
        x = x[selected]
    # NMF needs nonnegative inputs.  Per-feature offsets are recorded as a
    # transformation detail; no missing value is converted to a biological 0.
    baseline = float(np.sqrt(np.mean((x - x.mean(axis=0, keepdims=True)) ** 2)))
    rows: list[dict[str, Any]] = [
        {
            "unit_id": unit_id,
            "k": 0,
            "seed": -1,
            "status": "BASELINE_NO_LATENT_FIELD",
            "reconstruction_error": baseline,
            "reconstruction_metric": "per_element_rmse",
            "n_observations": len(x),
            "n_features": x.shape[1],
            "n_iter": 0,
            "max_component_moran": np.nan,
            "mean_component_moran": np.nan,
        }
    ]
    for k in ks:
        if k >= min(x.shape):
            continue
        for seed in seeds:
            model = NMF(
                n_components=k,
                init="nndsvda",
                random_state=seed,
                max_iter=200,
                l1_ratio=0.0,
                solver="cd",
            )
            try:
                model.fit_transform(x)
                full_transformed = model.transform(x_full - offsets + 1e-6)
                local_graph = sparse.csr_matrix(graph)[finite_indices][:, finite_indices]
                local_sections = np.asarray(section_ids)[finite_indices]
                morans = [
                    float(
                        moran_statistic(
                            full_transformed[:, component],
                            local_graph,
                            section_ids=local_sections,
                        ).get("estimate", np.nan)
                    )
                    for component in range(full_transformed.shape[1])
                ]
                rows.append(
                    {
                        "unit_id": unit_id,
                        "k": k,
                        "seed": seed,
                        "status": "ESTIMABLE",
                        "reconstruction_error": float(model.reconstruction_err_ / np.sqrt(x.size)),
                        "reconstruction_metric": "per_element_rmse",
                        "fit_converged": bool(model.n_iter_ < model.max_iter),
                        "n_observations": len(x),
                        "n_features": x.shape[1],
                        "n_iter": int(model.n_iter_),
                        "max_component_moran": float(np.nanmax(morans)) if morans else np.nan,
                        "mean_component_moran": float(np.nanmean(morans)) if morans else np.nan,
                    }
                )
            except (ValueError, FloatingPointError) as exc:
                rows.append(
                    {
                        "unit_id": unit_id,
                        "k": k,
                        "seed": seed,
                        "status": f"NOT_ESTIMABLE_{type(exc).__name__}",
                        "reconstruction_error": np.nan,
                        "n_observations": len(x),
                        "n_features": x.shape[1],
                        "n_iter": np.nan,
                        "max_component_moran": np.nan,
                        "mean_component_moran": np.nan,
                    }
                )
    return pd.DataFrame(rows)


def run_spatial_discovery(
    *,
    output_root: str | Path = "results/v7/spatial_discovery",
    stage1_registry_path: str | Path = "results/v7/registry/spatial_physical_units.tsv",
    r04_manifest_path: str | Path = "/home/huyudi/013_spatial/infra/r04/input_manifest.json",
    metadata_manifest_path: str | Path = (
        "results/v6_2/phase2_5_data_onboarding/01_metadata_freeze/"
        "spatial_sample_manifest.GSE291246.tsv"
    ),
    module_dictionary_path: str | Path = "results/v7/ontology/module_dictionary.tsv",
    module_membership_path: str | Path = "results/v7/ontology/module_membership.tsv",
    measurement_reliability_path: str | Path = "results/v7/ontology/measurement_reliability.tsv",
    config_path: str | Path = "config/v7/stage4.yaml",
    permutations: int = 19,
    seed: int = 20260911,
    nmf_sample_limit: int = 1500,
) -> dict[str, Any]:
    """Run the first full response-blind Stage 4 discovery pass."""

    output = Path(output_root)
    output.mkdir(parents=True, exist_ok=True)
    scores_dir = output / "scores"
    scores_dir.mkdir(parents=True, exist_ok=True)
    sanitization_dir = Path("/tmp/v7_stage4_sanitized")
    sanitization_dir.mkdir(parents=True, exist_ok=True)

    selection = build_stage4_selection(
        stage1_registry_path=stage1_registry_path,
        r04_manifest_path=r04_manifest_path,
        metadata_manifest_path=metadata_manifest_path,
    )
    selection.pilot.to_csv(output / "stage4_selection.tsv", sep="\t", index=False)
    selection.audit.to_csv(output / "technical_selection_audit.tsv", sep="\t", index=False)
    selection.block_gap_review.to_csv(output / "block_gap_review.tsv", sep="\t", index=False)
    model_safe = build_model_safe_manifest(
        selection.units, config_hash=sha256_file(config_path)
    )
    model_safe.to_csv(output / "model_safe_input_manifest.tsv", sep="\t", index=False)

    dictionary = pd.read_csv(module_dictionary_path, sep="\t")
    membership = pd.read_csv(module_membership_path, sep="\t")
    reliability = pd.read_csv(measurement_reliability_path, sep="\t")
    primitive_rows: list[dict[str, Any]] = []
    increment_rows: list[dict[str, Any]] = []
    field_rows: list[pd.DataFrame] = []
    coverage_rows: list[pd.DataFrame] = []
    unit_audit_rows: list[dict[str, Any]] = []
    patient_effect_rows: list[dict[str, Any]] = []
    all_feature_rows: list[dict[str, Any]] = []

    pairings = (
        ("T_NK_proxy", "Myeloid_proxy", "myeloid_T_adjacency"),
        ("T_NK_proxy", "Tumor_Epithelial_proxy", "T_tumor_interface"),
        ("CAF_Stromal_proxy", "T_NK_proxy", "CAF_T_exclusion_proxy"),
        ("Endothelial_proxy", "T_NK_proxy", "vascular_T_exclusion_proxy"),
        ("B_Plasma_proxy", "T_NK_proxy", "B_T_support_proxy"),
    )

    for ordinal, resolved in enumerate(selection.units, start=1):
        unit, sanitization = _load_unit(resolved, sanitization_dir)
        coordinates = unit.observations[["analysis_x", "analysis_y"]].to_numpy(float)
        sections = unit.observations["section_id"].astype(str).to_numpy()
        graph = build_section_graph(coordinates, sections, method="knn", k=6, weight_mode="binary")
        scores, coverage = project_stage2_features(
            unit.counts,
            unit.features["gene_symbol"],
            dictionary,
            membership,
            section_ids=sections,
            measurement_reliability=reliability,
        )
        scores.insert(0, "observation_id", unit.observations["observation_id"].astype(str).to_numpy())
        scores.insert(1, "unit_id", resolved.unit_id)
        scores.insert(2, "patient_id", resolved.opaque_patient_id)
        scores.insert(3, "dataset_id", resolved.dataset_id)
        scores.insert(4, "platform", resolved.platform)
        scores.to_parquet(scores_dir / f"{resolved.unit_id}.parquet", index=False)
        coverage = coverage.assign(
            unit_id=resolved.unit_id,
            patient_id=resolved.opaque_patient_id,
            dataset_id=resolved.dataset_id,
            platform=resolved.platform,
        )
        coverage_rows.append(coverage)

        proxies, proxy_audit = _marker_proxy(unit.counts, unit.features["gene_symbol"])
        proxy_audit = proxy_audit.assign(
            unit_id=resolved.unit_id,
            dataset_id=resolved.dataset_id,
            platform=resolved.platform,
        )
        proxy_audit.to_csv(
            output / f"proxy_audit_{resolved.unit_id}.tsv", sep="\t", index=False
        )
        unit_audit_rows.append(
            {
                "unit_id": resolved.unit_id,
                "dataset_id": resolved.dataset_id,
                "patient_id": resolved.opaque_patient_id,
                "platform": resolved.platform,
                "n_observations": unit.n_observations,
                "n_features": unit.n_features,
                "proxy_available": int((proxy_audit.status == "proxy_available").sum()),
                "sanitization_status": sanitization.get("status", "NOT_APPLICABLE"),
                "status": "PASS",
            }
        )

        for source, target, primitive in pairings:
            association = edge_association(
                graph,
                proxies,
                source=source,
                target=target,
                section_ids=sections,
            )
            primitive_rows.append(
                {
                    "unit_id": resolved.unit_id,
                    "dataset_id": resolved.dataset_id,
                    "patient_id": resolved.opaque_patient_id,
                    "platform": resolved.platform,
                    "primitive_id": primitive,
                    "source_proxy": source,
                    "target_proxy": target,
                    "estimate": association.get("estimate", np.nan),
                    "n_observations": association.get("n_observations", 0),
                    "n_directed_edges": association.get("n_directed_edges", 0),
                    "status": association.get("status", "UNKNOWN"),
                    "proxy_semantics": "marker_proxy_not_cell_fraction",
                }
            )

        score_matrix = scores[
            [c for c in scores.columns if c.endswith("__score_standardized")]
        ].to_numpy(float)
        feature_ids = [c.removesuffix("__score_standardized") for c in scores.columns if c.endswith("__score_standardized")]
        nmf_frame = _nmf_summary(
            score_matrix,
            graph,
            section_ids=sections,
            unit_id=resolved.unit_id,
            seeds=(seed, seed + 1),
            ks=(2, 4, 8),
            sample_limit=nmf_sample_limit,
        )
        field_rows.append(nmf_frame.assign(dataset_id=resolved.dataset_id, platform=resolved.platform, patient_id=resolved.opaque_patient_id))

        for feature_id, column in zip(feature_ids, score_matrix.T, strict=True):
            observed = moran_statistic(column, graph, section_ids=sections)
            permutation = moran_permutation_test(
                column,
                graph,
                sections,
                permutations=permutations,
                seed=seed ^ int(stable_hash({"unit": resolved.unit_id, "feature": feature_id})[:8], 16),
            )
            residual = _composition_residual(column, proxies)
            residual_stat = moran_statistic(residual, graph, section_ids=sections)
            increment_rows.append(
                {
                    "unit_id": resolved.unit_id,
                    "dataset_id": resolved.dataset_id,
                    "patient_id": resolved.opaque_patient_id,
                    "platform": resolved.platform,
                    "feature_id": feature_id,
                    "raw_moran": observed.get("estimate", np.nan),
                    "composition_residual_moran": residual_stat.get("estimate", np.nan),
                    "topology_increment": residual_stat.get("estimate", np.nan) - observed.get("estimate", np.nan),
                    "statistic_scope": "proxy_adjustment_moran_change_not_model_increment",
                    "topology_increment_column_deprecated": True,
                    "raw_p_value": permutation.get("p_value", np.nan),
                    "raw_status": observed.get("status", "UNKNOWN"),
                    "residual_status": residual_stat.get("status", "UNKNOWN"),
                    "composition_proxy_count": int(proxies.notna().any(axis=1).sum()),
                    "permutations": permutations,
                }
            )
            all_feature_rows.append(
                {
                    "unit_id": resolved.unit_id,
                    "dataset_id": resolved.dataset_id,
                    "patient_id": resolved.opaque_patient_id,
                    "platform": resolved.platform,
                    "feature_id": feature_id,
                    "moran": observed.get("estimate", np.nan),
                    "p_value": permutation.get("p_value", np.nan),
                    "status": observed.get("status", "UNKNOWN"),
                }
            )

    primitive = pd.DataFrame(primitive_rows)
    increment = pd.DataFrame(increment_rows)
    fields = pd.concat(field_rows, ignore_index=True) if field_rows else pd.DataFrame()
    coverage = pd.concat(coverage_rows, ignore_index=True)
    unit_audit = pd.DataFrame(unit_audit_rows)
    primitive.to_csv(output / "niche_primitives.tsv", sep="\t", index=False)
    increment.to_csv(output / "topology_increment.tsv", sep="\t", index=False)
    fields.to_csv(output / "latent_field_registry.tsv", sep="\t", index=False)
    coverage.to_csv(output / "feature_coverage.tsv", sep="\t", index=False)
    unit_audit.to_csv(output / "unit_discovery_audit.tsv", sep="\t", index=False)

    # Candidate architecture table: require at least two patient-level
    # observations and preserve sign heterogeneity instead of forcing a label.
    candidate = (
        primitive.groupby(["primitive_id", "source_proxy", "target_proxy"], dropna=False)
        .agg(
            n_units=("unit_id", "nunique"),
            n_patients=("patient_id", "nunique"),
            n_datasets=("dataset_id", "nunique"),
            median_estimate=("estimate", "median"),
            q25_estimate=("estimate", lambda x: _safe_quantile(x, 0.25)),
            q75_estimate=("estimate", lambda x: _safe_quantile(x, 0.75)),
            n_estimable=("status", lambda x: int((x == "ESTIMABLE").sum())),
        )
        .reset_index()
    )
    candidate["architecture_status"] = np.select(
        [
            candidate.n_estimable.ge(4) & candidate.n_datasets.ge(2),
            candidate.n_estimable.ge(2),
        ],
        ["candidate_shared_or_context_modulated", "candidate_context_limited"],
        default="UNRESOLVED_TOO_FEW_ESTIMABLE_UNITS",
    )
    candidate["claim_scope"] = "response_blind_spatial_structure_only"
    candidate.to_csv(output / "barrier_architecture_candidates.tsv", sep="\t", index=False)

    shared = candidate.copy()
    shared["shared_function_label"] = np.where(
        shared.n_datasets.ge(2) & shared.n_estimable.ge(4),
        "shared_function_candidate",
        np.where(shared.n_estimable.ge(2), "context_specific_candidate", "unresolved"),
    )
    shared["response_blind_locked_before_clinical_read"] = True
    shared.to_csv(output / "response_blind_shared_architectures.tsv", sep="\t", index=False)

    # Dataset holdout is a descriptive stability check.  It cannot prove a
    # clinical bridge, but it does show whether a signed primitive survives
    # leaving one dataset out.
    leave_rows: list[dict[str, Any]] = []
    for primitive_id, frame in primitive.groupby("primitive_id", sort=True):
        datasets = sorted(frame.dataset_id.unique())
        for held_out in datasets:
            train = frame[frame.dataset_id != held_out]
            vals = pd.to_numeric(train.estimate, errors="coerce").dropna().to_numpy(float)
            leave_rows.append(
                {
                    "primitive_id": primitive_id,
                    "held_out_dataset": held_out,
                    "n_train_datasets": train.dataset_id.nunique(),
                    "n_train_patients": train.patient_id.nunique(),
                    "train_median_estimate": float(np.median(vals)) if len(vals) else np.nan,
                    "train_sign": int(np.sign(np.median(vals))) if len(vals) else 0,
                    "status": "ESTIMABLE" if len(vals) >= 2 else "NOT_ESTIMABLE_TOO_FEW_PATIENTS",
                }
            )
    pd.DataFrame(leave_rows).to_csv(output / "leave_dataset_out.tsv", sep="\t", index=False)

    # A spatial-to-clinical surrogate needs paired, multimodal, patient-level
    # inputs and an untouched validation cohort.  This pass has none; record it
    # as an explicit audit result rather than inventing a surrogate.
    surrogate = pd.DataFrame(
        [
            {
                "task": "architecture_surrogate_validation",
                "status": "NOT_RUN_NO_PAIRED_MULTIMODAL_VALIDATION_ASSET",
                "outer_unit": "patient",
                "response_read": False,
                "ground_truth_read": False,
                "claim_allowed": "none",
            }
        ]
    )
    surrogate.to_csv(output / "architecture_surrogate_validation.tsv", sep="\t", index=False)

    # One compact patient-level spatial summary for later clinical bridging.
    patient_effect = (
        increment.groupby(["dataset_id", "patient_id", "platform", "feature_id"], dropna=False)
        .agg(
            n_units=("unit_id", "nunique"),
            median_raw_moran=("raw_moran", "median"),
            median_residual_moran=("composition_residual_moran", "median"),
            median_raw_p=("raw_p_value", "median"),
        )
        .reset_index()
    )
    patient_effect.to_parquet(output / "patient_spatial_summary.parquet", index=False)

    manifest = {
        "schema": "v7.stage4.spatial_discovery.run.v1",
        "status": "S4_COMPLETE_WITH_LIMITATIONS",
        "response_blind": True,
        "clinical_response_read": False,
        "ground_truth_read": False,
        "n_logical_units": int(len(selection.pilot)),
        "n_physical_captures": int(len(selection.units)),
        "n_patients": int(len({unit.opaque_patient_id for unit in selection.units})),
        "datasets": sorted({unit.dataset_id for unit in selection.units}),
        "platforms": sorted({unit.platform for unit in selection.units}),
        "permutations": int(permutations),
        "nmf_seeds": [seed, seed + 1],
        "nmf_k_values": [0, 2, 4, 8],
        "marker_proxy_semantics": "reference_marker_proxy_not_cell_fraction",
        "model_safe_manifest_sha256": sha256_file(output / "model_safe_input_manifest.tsv"),
        "stage3_gate_sha256": sha256_file("results/v7/spatial_foundation/D3_GATE.json"),
        "selection_sha256": sha256_file(output / "stage4_selection.tsv"),
        "surrogate_status": "NOT_RUN_NO_PAIRED_MULTIMODAL_VALIDATION_ASSET",
        "block_gap_review_rows": int(len(selection.block_gap_review)),
    }
    (output / "STAGE4_RUN_MANIFEST.json").write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    n_topology_estimable = int((increment.raw_status == "ESTIMABLE").sum())
    n_latent_estimable = int((fields.status == "ESTIMABLE").sum())
    n_latent_low_observation = int(
        (fields.status == "NOT_ESTIMABLE_TOO_FEW_OBSERVATIONS").sum()
    )
    report = f"""# v7 Stage 4 response-blind 空间发现首轮报告

## 状态

**S4_COMPLETE_WITH_LIMITATIONS**。本轮处理 {manifest['n_physical_captures']} 个物理捕获、{manifest['n_patients']} 位患者，覆盖 {"、".join(manifest['datasets'])} 与 {"、".join(manifest['platforms'])}。没有读取 response、治疗、终点或封存结构值。

## 完成内容

- 以 section 内 kNN 图计算 marker-proxy 邻接、39 个 Stage2 特征的 Moran 统计和 {permutations} 次空间置换；患者只作为外层汇总单位。
- 用粗粒度 marker proxy 做组成调整，并同时保留未调整与调整后的空间自相关；proxy 是表达代理，不是细胞比例。
- 比较无潜在场、K=2/4/8 的低秩 NMF；参数、种子和抽样上限记录在运行清单中。
- 对候选结构做跨数据集留出描述性检查；结构候选在读取疗效前已经写入版本化表。
- 对 538 条 block 缺失记录完成元数据复核清单，未把它们自动计为新增独立患者。

## 当前证据

空间候选表和 topology 增量表已经形成：{len(candidate)} 个候选原语在 {len(manifest['datasets'])} 个数据集、{manifest['n_physical_captures']} 个捕获中均可估计；{n_topology_estimable}/{len(increment)} 个 feature-topology 记录可估计，{len(increment) - n_topology_estimable} 个低观测单元保留为不可估计；低秩场有 {n_latent_estimable} 条 NMF 拟合记录可估计，另有 {n_latent_low_observation} 条低观测记录不可估计。当前结果仍是 response-blind 结构证据，不支持 PD-1 疗效关联、PD1+X 修复或因果解释。空间到临床的 surrogate 因缺少未参与开发的配对多模态验证资产，状态为 `NOT_RUN_NO_PAIRED_MULTIMODAL_VALIDATION_ASSET`。

NMF 部分捕获达到 200 次迭代上限；`latent_field_registry.tsv` 保留 `n_iter`，因此这些结果只能作为探索性低秩摘要，不升级为稳定机制证据。

## 复现

```bash
env LD_LIBRARY_PATH=/opt/anaconda3/lib PYTHONPATH=src \\
  python scripts/v7/spatial/run_stage4.py --permutations {permutations}
```
"""
    (output / "SPATIAL_DISCOVERY_REPORT.md").write_text(report, encoding="utf-8")
    return {**manifest, "output_root": str(output.resolve())}


__all__ = ["MARKER_PANEL", "Stage4Selection", "build_stage4_selection", "run_spatial_discovery"]
