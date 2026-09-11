"""Canonical response-blind spatial input layer for v7 Stage 3."""

from .contracts import (
    ADAPTER_VERSION,
    RAW_INTEGER_COUNTS,
    SpatialContractError,
    SpatialIdentity,
    SpatialUnitData,
    assert_model_safe_fields,
    audit_raw_integer_counts,
    validate_spatial_unit,
)
from .identity import namespace_observation_ids, validate_identity_collection
from .adapters import load_h5ad_counts, load_visium_10x, load_xenium_h5ad

__all__ = [
    "ADAPTER_VERSION",
    "RAW_INTEGER_COUNTS",
    "SpatialContractError",
    "SpatialIdentity",
    "SpatialUnitData",
    "assert_model_safe_fields",
    "audit_raw_integer_counts",
    "validate_spatial_unit",
    "namespace_observation_ids",
    "validate_identity_collection",
    "load_h5ad_counts",
    "load_visium_10x",
    "load_xenium_h5ad",
]
