"""Adapter for spatial AnnData with an explicit raw-count source."""

from __future__ import annotations

from pathlib import Path
from collections.abc import Mapping
from typing import Any, Iterable

import anndata as ad
import h5py
import numpy as np

from ..contracts import SpatialContractError, SpatialIdentity, SpatialUnitData, assert_model_safe_fields
from .base import (
    build_feature_table,
    build_observation_table,
    make_spatial_unit,
    materialize_matrix,
)


def _nested_field_names(value: Any, prefix: str = "") -> list[str]:
    names: list[str] = []
    if isinstance(value, Mapping):
        for key, nested in value.items():
            path = f"{prefix}/{key}" if prefix else str(key)
            names.append(path)
            names.extend(_nested_field_names(nested, path))
    elif isinstance(value, np.ndarray) and value.dtype.names:
        names.extend(f"{prefix}/{name}" for name in value.dtype.names)
    return names


def _audit_h5ad_field_names(adata: ad.AnnData) -> None:
    namespaces: Iterable[tuple[str, Iterable[object]]] = (
        ("h5ad.obs", list(adata.obs.columns) + [adata.obs.index.name]),
        ("h5ad.var", list(adata.var.columns) + [adata.var.index.name]),
        ("h5ad.obsm", adata.obsm.keys()),
        ("h5ad.varm", adata.varm.keys()),
        ("h5ad.obsp", adata.obsp.keys()),
        ("h5ad.varp", adata.varp.keys()),
        ("h5ad.layers", adata.layers.keys()),
        ("h5ad.uns", _nested_field_names(adata.uns)),
    )
    for context, names in namespaces:
        assert_model_safe_fields((name for name in names if name is not None), context)


def _audit_h5ad_storage_fields(path: Path) -> None:
    """Inspect HDF5 object names before AnnData can materialize table values."""

    with h5py.File(path, "r") as handle:
        names: list[str] = []
        handle.visit(names.append)
    assert_model_safe_fields(names, "h5ad storage fields")


def load_h5ad_counts(
    path: str | Path,
    identity: SpatialIdentity,
    *,
    counts_layer: str,
    coordinate_key: str = "spatial",
    coordinate_unit: str,
    platform: str = "spatial_h5ad",
    resolution: str = "spot",
    feature_id_column: str | None = None,
    gene_symbol_column: str | None = None,
    feature_type_column: str | None = None,
    in_tissue_column: str | None = "in_tissue",
    segmentation_column: str | None = None,
) -> SpatialUnitData:
    """Load a spatial h5ad without guessing whether X or a layer is raw counts.

    ``counts_layer`` must be exactly ``"X"`` or an existing ``adata.layers``
    key.  Integer-valued storage is audited after loading.
    """

    source = Path(path)
    if not isinstance(counts_layer, str) or not counts_layer.strip():
        raise SpatialContractError("counts_layer must explicitly name X or an AnnData layer")
    if not isinstance(coordinate_key, str) or not coordinate_key.strip():
        raise SpatialContractError("coordinate_key must be explicit")
    if not source.is_file():
        raise SpatialContractError(f"h5ad source is not a readable file: {source}")

    _audit_h5ad_storage_fields(source)
    adata = ad.read_h5ad(source, backed="r")
    try:
        _audit_h5ad_field_names(adata)
        if coordinate_key not in adata.obsm:
            raise SpatialContractError(f"h5ad is missing obsm[{coordinate_key!r}]")
        if counts_layer == "X":
            matrix_source = adata.X
        elif counts_layer in adata.layers:
            matrix_source = adata.layers[counts_layer]
        else:
            raise SpatialContractError(f"h5ad is missing counts layer {counts_layer!r}")

        counts = materialize_matrix(matrix_source)
        coordinates = np.asarray(adata.obsm[coordinate_key]).copy()
        native_ids = adata.obs_names.astype(str).copy()
        if feature_id_column is not None:
            if feature_id_column not in adata.var:
                raise SpatialContractError(
                    f"h5ad is missing feature ID column {feature_id_column!r}"
                )
            native_features = adata.var[feature_id_column].astype("string").copy()
        else:
            native_features = adata.var_names.astype(str).copy()
        if gene_symbol_column is not None:
            if gene_symbol_column not in adata.var:
                raise SpatialContractError(
                    f"h5ad is missing gene symbol column {gene_symbol_column!r}"
                )
            gene_symbols = adata.var[gene_symbol_column].astype("string").copy()
        else:
            gene_symbols = native_features
        if feature_type_column is not None:
            if feature_type_column not in adata.var:
                raise SpatialContractError(
                    f"h5ad is missing feature type column {feature_type_column!r}"
                )
            feature_types = adata.var[feature_type_column].astype("string").copy()
        else:
            feature_types = None
        if in_tissue_column is not None and in_tissue_column in adata.obs:
            in_tissue = adata.obs[in_tissue_column].copy()
        else:
            in_tissue = None
        if segmentation_column is not None:
            if segmentation_column not in adata.obs:
                raise SpatialContractError(
                    f"h5ad is missing segmentation column {segmentation_column!r}"
                )
            segmentation_ids = adata.obs[segmentation_column].copy()
        else:
            segmentation_ids = None
    finally:
        adata.file.close()

    observations = build_observation_table(
        native_ids=native_ids,
        coordinates=coordinates,
        identity=identity,
        native_coordinate_unit=coordinate_unit,
        in_tissue=in_tissue,
        segmentation_ids=segmentation_ids,
    )
    features = build_feature_table(
        native_feature_ids=native_features,
        gene_symbols=gene_symbols,
        feature_types=feature_types,
    )
    return make_spatial_unit(
        counts=counts,
        observations=observations,
        features=features,
        identity=identity,
        platform=platform,
        resolution=resolution,
        coordinate_system=coordinate_unit,
        source_assets=[source],
        audit={
            "counts_source": counts_layer,
            "coordinate_source": f"obsm/{coordinate_key}",
            "storage_dtype": str(counts.dtype),
            "in_tissue_status": "provided" if in_tissue is not None else "not_provided",
            "segmentation_status": "provided" if segmentation_ids is not None else "not_provided",
        },
    )
