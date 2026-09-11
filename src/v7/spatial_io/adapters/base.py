"""Shared constructors for Stage 3 spatial adapters."""

from __future__ import annotations

from pathlib import Path
from typing import Any, Sequence

import numpy as np
import pandas as pd
from scipy import sparse

from ..contracts import (
    ADAPTER_VERSION,
    RAW_INTEGER_COUNTS,
    SpatialContractError,
    SpatialIdentity,
    SpatialUnitData,
    audit_raw_integer_counts,
    validate_spatial_unit,
)
from ..identity import namespace_observation_ids


def materialize_matrix(matrix: Any) -> Any:
    """Materialize an AnnData-backed matrix without changing its values."""

    if hasattr(matrix, "to_memory"):
        return matrix.to_memory()
    if sparse.issparse(matrix):
        return matrix.copy()
    try:
        return matrix[:]
    except (TypeError, IndexError):
        return np.asarray(matrix)


def normalize_in_tissue(values: Sequence[object] | None, n_observations: int) -> pd.Series:
    """Normalize a technical tissue mask while rejecting unknown encodings."""

    if values is None:
        return pd.Series([pd.NA] * n_observations, dtype="boolean")
    series = pd.Series(list(values))
    if len(series) != n_observations:
        raise SpatialContractError("in_tissue length does not match the count matrix")
    if pd.api.types.is_bool_dtype(series):
        return series.astype("boolean")
    numeric = pd.to_numeric(series, errors="coerce")
    if numeric.isna().any() or not numeric.isin([0, 1]).all():
        raise SpatialContractError("in_tissue must contain only boolean or 0/1 values")
    return numeric.astype(bool).astype("boolean")


def build_observation_table(
    *,
    native_ids: Sequence[object] | pd.Index,
    coordinates: Any,
    identity: SpatialIdentity,
    native_coordinate_unit: str,
    analysis_coordinate_unit: str | None = None,
    in_tissue: Sequence[object] | None = None,
    segmentation_ids: Sequence[object] | None = None,
) -> pd.DataFrame:
    """Build the shared observation table without discarding native IDs."""

    ids = namespace_observation_ids(native_ids, identity.capture_id)
    xy = np.asarray(coordinates)
    if xy.ndim != 2 or xy.shape[0] != len(ids) or xy.shape[1] < 2:
        raise SpatialContractError("coordinates must have shape n_observations x >=2")
    xy = np.asarray(xy[:, :2], dtype=float)
    if not np.isfinite(xy).all():
        raise SpatialContractError("coordinates contain non-finite values")
    if not native_coordinate_unit.strip():
        raise SpatialContractError("native_coordinate_unit must be explicit")
    analysis_unit = analysis_coordinate_unit or native_coordinate_unit

    if segmentation_ids is None:
        segmentation = pd.Series([pd.NA] * len(ids), dtype="string")
    else:
        segmentation = pd.Series(list(segmentation_ids), dtype="string")
        if len(segmentation) != len(ids):
            raise SpatialContractError("segmentation_id length does not match the count matrix")

    observations = ids.assign(
        capture_id=identity.capture_id,
        section_id=identity.opaque_section_id,
        native_x=xy[:, 0],
        native_y=xy[:, 1],
        native_coordinate_unit=native_coordinate_unit,
        analysis_x=xy[:, 0],
        analysis_y=xy[:, 1],
        analysis_coordinate_unit=analysis_unit,
        in_tissue=normalize_in_tissue(in_tissue, len(ids)),
        segmentation_id=segmentation,
        observation_qc="PASS",
    )
    return observations


def build_feature_table(
    *,
    native_feature_ids: Sequence[object] | pd.Index,
    gene_symbols: Sequence[object] | pd.Index | None = None,
    feature_types: Sequence[object] | pd.Index | None = None,
) -> pd.DataFrame:
    """Build the measured feature universe; absent genes are not synthesized."""

    native = pd.Series(list(native_feature_ids), dtype="string")
    if native.isna().any() or native.str.strip().eq("").any():
        raise SpatialContractError("native feature IDs contain missing or empty values")
    if native.duplicated().any():
        raise SpatialContractError("native feature IDs must be unique")
    symbols = native.copy() if gene_symbols is None else pd.Series(list(gene_symbols), dtype="string")
    types = (
        pd.Series(["Gene Expression"] * len(native), dtype="string")
        if feature_types is None
        else pd.Series(list(feature_types), dtype="string")
    )
    if len(symbols) != len(native) or len(types) != len(native):
        raise SpatialContractError("feature annotations do not align to the count matrix")
    return pd.DataFrame(
        {
            "feature_id": native.astype(str),
            "native_feature_id": native.astype(str),
            "gene_symbol": symbols,
            "feature_type": types,
            "is_measured": True,
        }
    )


def make_spatial_unit(
    *,
    counts: Any,
    observations: pd.DataFrame,
    features: pd.DataFrame,
    identity: SpatialIdentity,
    platform: str,
    resolution: str,
    coordinate_system: str,
    source_assets: Sequence[str | Path],
    audit: dict[str, Any],
) -> SpatialUnitData:
    """Construct and validate one canonical unit."""

    canonical_counts = audit_raw_integer_counts(counts, f"{platform} adapter")
    complete_audit = {
        "adapter_version": ADAPTER_VERSION,
        "adapter_status": "PASS",
        "n_observations": int(canonical_counts.shape[0]),
        "n_features": int(canonical_counts.shape[1]),
        "count_total": int(canonical_counts.sum()),
        "nnz": int(canonical_counts.nnz),
        **audit,
    }
    unit = SpatialUnitData(
        counts=canonical_counts,
        observations=observations.reset_index(drop=True),
        features=features.reset_index(drop=True),
        identity=identity,
        platform=platform,
        modality="spatial_transcriptomics",
        resolution=resolution,
        counts_semantics=RAW_INTEGER_COUNTS,
        coordinate_system=coordinate_system,
        source_assets=tuple(str(Path(path).resolve()) for path in source_assets),
        audit=complete_audit,
    )
    return validate_spatial_unit(unit)
