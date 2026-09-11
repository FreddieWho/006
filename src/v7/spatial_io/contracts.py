"""Fail-closed data contracts shared by v7 spatial input adapters."""

from __future__ import annotations

from dataclasses import dataclass, field
import re
from typing import Any, Iterable, Mapping

import numpy as np
import pandas as pd
from scipy import sparse


ADAPTER_VERSION = "v7-stage3-spatial-io-1"
RAW_INTEGER_COUNTS = "raw_integer_counts"

OBSERVATION_COLUMNS = (
    "observation_id",
    "native_observation_id",
    "capture_id",
    "section_id",
    "native_x",
    "native_y",
    "native_coordinate_unit",
    "analysis_x",
    "analysis_y",
    "analysis_coordinate_unit",
    "in_tissue",
    "segmentation_id",
    "observation_qc",
)

FEATURE_COLUMNS = (
    "feature_id",
    "native_feature_id",
    "gene_symbol",
    "feature_type",
    "is_measured",
)

_FORBIDDEN_FIELD_PATTERN = re.compile(
    r"(?:^|[^a-z0-9])(?:"
    r"response|responder|non[_ -]?responder|outcome|endpoint|recist|survival|"
    r"ground[_ -]?truth|structure[_ -]?gt|gt|gt[_ -]?(?:mask|distance|label)|"
    r"(?:tls|tumou?r[_ -]?boundary)[_ -]?(?:gt|mask|distance|label)?"
    r")(?:$|[^a-z0-9])",
    flags=re.IGNORECASE,
)


class SpatialContractError(ValueError):
    """Raised when an input cannot satisfy the Stage 3 discovery contract."""


@dataclass(frozen=True)
class SpatialIdentity:
    """Opaque patient hierarchy carried by one technical capture.

    Values are identifiers only.  Treatment and clinical endpoint fields are
    intentionally absent from this object.
    """

    dataset_id: str
    opaque_patient_id: str
    opaque_block_id: str
    opaque_section_id: str
    capture_id: str
    leakage_group_id: str
    relationship_provenance: str
    relationship_confidence: str

    def __post_init__(self) -> None:
        for name, value in self.__dict__.items():
            if not isinstance(value, str) or not value.strip():
                raise SpatialContractError(f"SpatialIdentity.{name} must be a non-empty string")


@dataclass
class SpatialUnitData:
    """Canonical in-memory object returned by every spatial adapter."""

    counts: sparse.csr_matrix
    observations: pd.DataFrame
    features: pd.DataFrame
    identity: SpatialIdentity
    platform: str
    modality: str
    resolution: str
    counts_semantics: str
    coordinate_system: str
    source_assets: tuple[str, ...]
    audit: Mapping[str, Any] = field(default_factory=dict)

    @property
    def n_observations(self) -> int:
        return self.counts.shape[0]

    @property
    def n_features(self) -> int:
        return self.counts.shape[1]

    @property
    def count_total(self) -> int:
        return int(self.counts.sum())


def forbidden_field_names(field_names: Iterable[object]) -> list[str]:
    """Return discovery-input field names that encode endpoints or sealed GT."""

    return sorted(
        {
            str(name)
            for name in field_names
            if _FORBIDDEN_FIELD_PATTERN.search(str(name).strip()) is not None
        }
    )


def assert_model_safe_fields(field_names: Iterable[object], context: str) -> None:
    """Fail without printing any values when forbidden field names are present."""

    forbidden = forbidden_field_names(field_names)
    if forbidden:
        raise SpatialContractError(
            f"{context} contains forbidden response/outcome/ground-truth fields: {forbidden}"
        )


def audit_raw_integer_counts(matrix: Any, context: str) -> sparse.csr_matrix:
    """Validate count semantics without rounding, clipping, or imputing values."""

    if sparse.issparse(matrix):
        counts = matrix.tocsr(copy=True)
        if not counts.has_canonical_format:
            raise SpatialContractError(f"{context} sparse matrix has duplicate or unsorted entries")
        values = counts.data
    else:
        array = np.asarray(matrix)
        if array.ndim != 2:
            raise SpatialContractError(f"{context} counts must be a two-dimensional matrix")
        values = array.ravel()
        counts = sparse.csr_matrix(array)

    if values.size and not np.isfinite(values).all():
        raise SpatialContractError(f"{context} counts contain non-finite values")
    if values.size and (values < 0).any():
        raise SpatialContractError(f"{context} counts contain negative values")
    if values.size and not np.equal(values, np.floor(values)).all():
        raise SpatialContractError(
            f"{context} is not raw integer counts; normalized or inferred expression is forbidden"
        )
    if counts.shape[0] == 0 or counts.shape[1] == 0:
        raise SpatialContractError(f"{context} counts matrix is empty")
    return counts


def validate_spatial_unit(unit: SpatialUnitData) -> SpatialUnitData:
    """Validate alignment, identity, coordinates and count semantics."""

    if unit.counts_semantics != RAW_INTEGER_COUNTS:
        raise SpatialContractError(
            f"unsupported counts_semantics={unit.counts_semantics!r}; raw integer counts are required"
        )
    unit.counts = audit_raw_integer_counts(unit.counts, "SpatialUnitData")
    if unit.counts.shape != (len(unit.observations), len(unit.features)):
        raise SpatialContractError(
            "counts shape does not match observation/feature tables: "
            f"{unit.counts.shape} != {(len(unit.observations), len(unit.features))}"
        )

    missing_obs = [column for column in OBSERVATION_COLUMNS if column not in unit.observations]
    missing_features = [column for column in FEATURE_COLUMNS if column not in unit.features]
    if missing_obs:
        raise SpatialContractError(f"observation table is missing columns: {missing_obs}")
    if missing_features:
        raise SpatialContractError(f"feature table is missing columns: {missing_features}")
    assert_model_safe_fields(unit.observations.columns, "observation table")
    assert_model_safe_fields(unit.features.columns, "feature table")
    assert_model_safe_fields(unit.audit.keys(), "adapter audit")

    observations = unit.observations
    for column in ("observation_id", "native_observation_id"):
        values = observations[column].astype("string")
        if values.isna().any() or values.str.strip().eq("").any():
            raise SpatialContractError(f"{column} contains missing or empty values")
        if values.duplicated().any():
            raise SpatialContractError(f"{column} must be unique within a capture")
    if not observations["capture_id"].astype(str).eq(unit.identity.capture_id).all():
        raise SpatialContractError("observation capture_id does not match SpatialIdentity")
    if not observations["section_id"].astype(str).eq(unit.identity.opaque_section_id).all():
        raise SpatialContractError("observation section_id does not match SpatialIdentity")

    coordinates = observations.loc[
        :, ["native_x", "native_y", "analysis_x", "analysis_y"]
    ].apply(pd.to_numeric, errors="coerce")
    if not np.isfinite(coordinates.to_numpy(dtype=float)).all():
        raise SpatialContractError("spatial coordinates must be finite numeric values")
    if len(observations) > 1 and coordinates[["analysis_x", "analysis_y"]].drop_duplicates().shape[0] < 2:
        raise SpatialContractError("analysis coordinates are degenerate within the capture")

    features = unit.features
    for column in ("feature_id", "native_feature_id"):
        values = features[column].astype("string")
        if values.isna().any() or values.str.strip().eq("").any():
            raise SpatialContractError(f"{column} contains missing or empty values")
        if values.duplicated().any():
            raise SpatialContractError(f"{column} must be unique within a capture")
    if not features["is_measured"].astype(bool).all():
        raise SpatialContractError(
            "feature table may list only the measured panel; unavailable genes stay absent, not zero-filled"
        )
    if not unit.platform.strip() or not unit.coordinate_system.strip():
        raise SpatialContractError("platform and coordinate_system must be explicit")
    if not unit.source_assets or any(not str(path).strip() for path in unit.source_assets):
        raise SpatialContractError("at least one non-empty source asset locator is required")
    return unit
