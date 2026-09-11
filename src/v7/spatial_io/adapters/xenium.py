"""Cell-level Xenium h5ad adapter with explicit panel semantics."""

from __future__ import annotations

from pathlib import Path

from ..contracts import SpatialIdentity, SpatialUnitData, validate_spatial_unit
from .h5ad_counts import load_h5ad_counts


def load_xenium_h5ad(
    path: str | Path,
    identity: SpatialIdentity,
    *,
    counts_layer: str,
    coordinate_key: str = "spatial",
    segmentation_column: str | None = None,
    feature_id_column: str | None = None,
    gene_symbol_column: str | None = None,
) -> SpatialUnitData:
    """Load a Xenium cell-feature h5ad without inferring unmeasured genes.

    When no separate segmentation column is supplied, the native Xenium cell
    ID is retained as the segmentation ID; no transcript-level reconstruction
    or morphology-image feature extraction is performed.
    """

    unit = load_h5ad_counts(
        path,
        identity,
        counts_layer=counts_layer,
        coordinate_key=coordinate_key,
        coordinate_unit="micrometer",
        platform="Xenium",
        resolution="cell",
        feature_id_column=feature_id_column,
        gene_symbol_column=gene_symbol_column,
        segmentation_column=segmentation_column,
    )
    if segmentation_column is None:
        unit.observations["segmentation_id"] = unit.observations["native_observation_id"].astype(
            "string"
        )
        unit.audit = {
            **unit.audit,
            "segmentation_status": "native_cell_id",
            "panel_semantics": "targeted_measured_gene_universe",
            "unmeasured_genes": "absent_not_zero",
        }
    else:
        unit.audit = {
            **unit.audit,
            "panel_semantics": "targeted_measured_gene_universe",
            "unmeasured_genes": "absent_not_zero",
        }
    return validate_spatial_unit(unit)
