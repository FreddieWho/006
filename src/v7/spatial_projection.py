"""Response-blind projection of the frozen v7 Stage-2 vocabulary to space.

This module deliberately keeps the spatial layer independent from the Stage-2
scRNA diagnostics.  It reports what a platform actually measured and retains
missing panel members as missing; no gene or observation is imputed.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd
from scipy import sparse

from .ontology.scoring import score_legacy_fm


@dataclass(frozen=True)
class SpatialFeatureResult:
    feature_id: str
    score_native: np.ndarray
    score_standardized: np.ndarray
    expected_genes: tuple[str, ...]
    observed_genes: tuple[str, ...]
    coverage: float
    spatial_projection_status: str
    scoring_semantics: str
    source_scRNA_D2_class: str


def _robust_standardize(values: np.ndarray) -> np.ndarray:
    result = np.full(values.shape, np.nan, dtype=float)
    finite = np.isfinite(values)
    if not finite.any():
        return result
    center = float(np.median(values[finite]))
    mad = float(np.median(np.abs(values[finite] - center))) * 1.4826
    if mad <= 0:
        q25, q75 = np.percentile(values[finite], [25, 75])
        mad = float((q75 - q25) / 1.349)
    if mad <= 0:
        mad = float(np.std(values[finite]))
    if mad <= 0:
        mad = 1.0
    result[finite] = (values[finite] - center) / mad
    return result


def _collapse_symbols(
    counts: sparse.spmatrix, symbols: Iterable[object]
) -> tuple[sparse.csr_matrix, tuple[str, ...], dict[str, tuple[int, ...]]]:
    """Collapse duplicate gene symbols in count space with an audit mapping."""

    raw = tuple(str(value).strip() for value in symbols)
    if len(raw) != counts.shape[1] or any(not value for value in raw):
        raise ValueError("gene symbols must be non-empty and aligned to counts")
    unique = tuple(sorted(set(raw)))
    index = {symbol: i for i, symbol in enumerate(unique)}
    columns: dict[str, list[int]] = {symbol: [] for symbol in unique}
    for column, symbol in enumerate(raw):
        columns[symbol].append(column)
    rows = np.arange(len(raw), dtype=np.int64)
    target = np.asarray([index[symbol] for symbol in raw], dtype=np.int64)
    aggregator = sparse.csr_matrix(
        (np.ones(len(raw), dtype=float), (rows, target)),
        shape=(len(raw), len(unique)),
    )
    collapsed = sparse.csr_matrix(counts, dtype=float) @ aggregator
    return collapsed.tocsr(), unique, {key: tuple(value) for key, value in columns.items()}


def _log1p_cpm(counts: sparse.csr_matrix) -> tuple[sparse.csr_matrix, np.ndarray]:
    library = np.asarray(counts.sum(axis=1)).ravel().astype(float)
    valid = library > 0
    scale = np.divide(1_000_000.0, library, out=np.zeros_like(library), where=valid)
    logged = counts.multiply(scale[:, None]).tocsr()
    if logged.data.size:
        logged.data = np.log1p(logged.data)
        logged.eliminate_zeros()
    return logged, valid


def _weighted_component(
    logged: sparse.csr_matrix,
    positions: dict[str, int],
    weights: pd.DataFrame,
    direction: int,
) -> tuple[np.ndarray | None, tuple[str, ...], float]:
    subset = weights.loc[weights.direction.eq(direction)].copy()
    expected = tuple(sorted(set(subset.gene_symbol.astype(str))))
    observed = tuple(gene for gene in expected if gene in positions)
    if not observed:
        return None, observed, 0.0
    by_gene = subset.groupby("gene_symbol", sort=True).weight.sum().to_dict()
    raw_weights = np.asarray([float(by_gene[gene]) for gene in observed], dtype=float)
    if not np.isfinite(raw_weights).all() or np.any(raw_weights <= 0):
        raise ValueError("module membership contains invalid positive weights")
    normalized = raw_weights / raw_weights.sum()
    columns = [positions[gene] for gene in observed]
    values = np.asarray(logged[:, columns] @ normalized).ravel().astype(float)
    return values, observed, float(
        sum(by_gene[gene] for gene in observed) / max(sum(by_gene.values()), np.finfo(float).eps)
    )


def project_stage2_features(
    counts: sparse.spmatrix,
    gene_symbols: Iterable[object],
    module_dictionary: pd.DataFrame,
    module_membership: pd.DataFrame,
    *,
    section_ids: Iterable[object] | None = None,
    measurement_reliability: pd.DataFrame | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Project scoreable Stage-2 features and return wide scores plus an audit.

    The score table has one row per spatial observation and paired native /
    section-standardized columns for every scoreable feature.  ``audit`` keeps
    expected and observed gene sets and the independent spatial status.
    """

    matrix = sparse.csr_matrix(counts, dtype=float)
    matrix.sum_duplicates()
    matrix.eliminate_zeros()
    matrix, symbols, _ = _collapse_symbols(matrix, gene_symbols)
    logged, valid_rows = _log1p_cpm(matrix)
    positions = {symbol: i for i, symbol in enumerate(symbols)}
    dictionary = module_dictionary.copy()
    membership = module_membership.copy()
    required_dict = {"feature_id", "feature_family", "scoring_semantics", "scoreable"}
    required_members = {"feature_id", "gene_symbol", "direction", "weight"}
    if missing := sorted(required_dict - set(dictionary.columns)):
        raise ValueError(f"module dictionary missing columns: {missing}")
    if missing := sorted(required_members - set(membership.columns)):
        raise ValueError(f"module membership missing columns: {missing}")
    scoreable = dictionary.loc[
        dictionary.scoreable.astype(str).str.lower().eq("true")
        & dictionary.feature_family.isin(["legacy_fm", "mechanism_component"])
    ].copy()
    if scoreable.empty:
        raise ValueError("no scoreable Stage-2 features available")
    d2_map: dict[str, str] = {}
    if measurement_reliability is not None and {"feature_id", "D2_class"}.issubset(
        measurement_reliability.columns
    ):
        d2_map = {
            str(row.feature_id): str(row.D2_class)
            for row in measurement_reliability[["feature_id", "D2_class"]].itertuples()
        }
    scores: dict[str, np.ndarray] = {}
    audits: list[dict[str, object]] = []
    for row in scoreable.sort_values("feature_id", kind="mergesort").itertuples():
        feature_id = str(row.feature_id)
        members = membership.loc[membership.feature_id.eq(feature_id)].copy()
        members["gene_symbol"] = members.gene_symbol.astype(str).str.strip()
        members = members.loc[members.gene_symbol.ne("")]
        expected = tuple(sorted(set(members.gene_symbol)))
        if row.feature_family == "legacy_fm":
            weights = members.groupby("gene_symbol", sort=True).weight.sum().to_dict()
            legacy = score_legacy_fm(matrix, symbols, weights, layer="counts")
            native = np.asarray(legacy.log1p_topic_mass_per_million, dtype=float)
            observed = tuple(legacy.observed_genes)
            coverage = float(legacy.weight_coverage)
        else:
            positive, observed_positive, positive_coverage = _weighted_component(
                logged, positions, members, 1
            )
            negative, observed_negative, negative_coverage = _weighted_component(
                logged, positions, members, -1
            )
            if positive is not None and negative is not None:
                native = positive - negative
            elif positive is not None:
                native = positive
            elif negative is not None:
                native = -negative
            else:
                native = np.full(matrix.shape[0], np.nan, dtype=float)
            observed = tuple(sorted(set(observed_positive + observed_negative)))
            total_weight = float(members.groupby("gene_symbol").weight.sum().sum())
            coverage = float(
                (positive_coverage * members.loc[members.direction.eq(1), "weight"].sum()
                 + negative_coverage * members.loc[members.direction.eq(-1), "weight"].sum())
                / max(total_weight, np.finfo(float).eps)
            )
        native[~valid_rows] = np.nan
        standardized = _robust_standardize(native)
        scores[f"{feature_id}__score_native"] = native
        scores[f"{feature_id}__score_standardized"] = standardized
        if not observed:
            status = "not_estimable"
        elif len(observed) == len(expected):
            status = "measurable"
        else:
            status = "limited_with_uncertainty"
        audits.append(
            {
                "feature_id": feature_id,
                "feature_family": str(row.feature_family),
                "scoring_semantics": str(row.scoring_semantics),
                "expected_genes": ";".join(expected),
                "observed_genes": ";".join(observed),
                "expected_gene_count": len(expected),
                "observed_gene_count": len(observed),
                "coverage": coverage,
                "spatial_projection_status": status,
                "source_scRNA_D2_class": d2_map.get(feature_id, "not_provided"),
                "D2_scope": "scRNA_measurement_only",
                "n_finite_scores": int(np.isfinite(native).sum()),
                "n_observations": int(len(native)),
            }
        )
    score_frame = pd.DataFrame(scores)
    audit_frame = pd.DataFrame(audits).sort_values("feature_id").reset_index(drop=True)
    if section_ids is not None:
        sections = np.asarray(list(section_ids), dtype=object)
        if len(sections) != len(score_frame):
            raise ValueError("section_ids must align with counts")
        for feature_id in audit_frame.feature_id:
            key = f"{feature_id}__score_native"
            standard = np.full(len(score_frame), np.nan, dtype=float)
            for section in pd.unique(sections):
                indices = np.flatnonzero(sections == section)
                standard[indices] = _robust_standardize(score_frame.loc[indices, key].to_numpy(float))
            score_frame[f"{feature_id}__score_standardized"] = standard
    return score_frame, audit_frame


__all__ = ["SpatialFeatureResult", "project_stage2_features"]
