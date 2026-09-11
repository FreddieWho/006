"""Deterministic, layer-aware gene resolution and module scoring for v7.

The functions in this module deliberately fail closed when expression-layer
semantics are unknown.  Missing signature genes are excluded from the score;
they are never represented by artificial zero-expression columns.
"""

from __future__ import annotations

from dataclasses import dataclass
from collections import Counter
from enum import Enum
import math
import re
from typing import Iterable, Mapping, Sequence

import numpy as np
from scipy import sparse
from scipy.stats import rankdata


ArrayLike = np.ndarray | sparse.spmatrix
_ENSEMBL_GENE = re.compile(r"^(ENSG\d+)(?:\.\d+)?$")


class ScoringError(ValueError):
    """Base class for invalid or unscoreable measurement inputs."""


class UnknownExpressionLayerError(ScoringError):
    """Raised when an input layer is not one of the explicitly supported layers."""


class GeneResolutionError(ScoringError):
    """Raised when a gene namespace or reference mapping is invalid."""


class InvalidExpressionError(ScoringError):
    """Raised when the expression matrix violates its declared layer semantics."""


class InsufficientCoverageError(ScoringError):
    """Raised when no positive-weight member of a module is measurable."""


class ExpressionLayer(str, Enum):
    """Supported expression layers; spellings are intentionally exact."""

    COUNTS = "counts"
    CPM = "CPM"
    TPM = "TPM"
    LOG1P_CPM = "log1p_CPM"
    LOG2_TPM = "log2_TPM"


def parse_expression_layer(layer: ExpressionLayer | str) -> ExpressionLayer:
    """Return a supported layer or fail closed without guessing aliases."""

    if isinstance(layer, ExpressionLayer):
        return layer
    try:
        return ExpressionLayer(layer)
    except (TypeError, ValueError) as exc:
        allowed = ", ".join(item.value for item in ExpressionLayer)
        raise UnknownExpressionLayerError(
            f"Unsupported expression layer {layer!r}; expected exactly one of: {allowed}"
        ) from exc


def strip_ensembl_version(gene_id: str) -> str:
    """Remove a numeric version suffix from a valid human Ensembl gene ID."""

    value = str(gene_id).strip()
    match = _ENSEMBL_GENE.fullmatch(value)
    return match.group(1) if match else value


def _candidate_tuple(value: str | Iterable[str] | None) -> tuple[str, ...]:
    if value is None:
        return ()
    values = (value,) if isinstance(value, str) else tuple(value)
    return tuple(sorted({str(item).strip() for item in values if str(item).strip()}))


def _normalise_reference(
    reference: Mapping[str, str | Iterable[str] | None],
) -> dict[str, tuple[str, ...]]:
    normalised: dict[str, set[str]] = {}
    for raw_key, raw_value in reference.items():
        key = strip_ensembl_version(str(raw_key))
        if not key:
            continue
        normalised.setdefault(key, set()).update(_candidate_tuple(raw_value))
    return {key: tuple(sorted(values)) for key, values in normalised.items()}


@dataclass(frozen=True)
class GeneMappingRecord:
    """One auditable gene-identifier resolution decision."""

    input_id: str
    normalized_id: str
    resolved_symbol: str | None
    source: str
    status: str
    candidates: tuple[str, ...]


class GeneResolver:
    """Resolve Ensembl IDs with GENCODE v49 priority and rank-1 fallback.

    ``gencode_v49`` and ``org_hs_eg_db_rank1`` may map one ID to one symbol or
    to an iterable. Multiple distinct candidates are audited as ambiguous and
    excluded. A GENCODE entry, including an ambiguous one, always takes
    precedence over the fallback reference.

    Non-Ensembl identifiers are treated as exact symbols. Supplying
    ``known_symbols`` makes that path fail closed for symbols outside the set.
    """

    def __init__(
        self,
        gencode_v49: Mapping[str, str | Iterable[str] | None],
        org_hs_eg_db_rank1: Mapping[str, str | Iterable[str] | None],
        *,
        known_symbols: Iterable[str] | None = None,
    ) -> None:
        self._gencode = _normalise_reference(gencode_v49)
        self._fallback = _normalise_reference(org_hs_eg_db_rank1)
        self._known_symbols = (
            None
            if known_symbols is None
            else frozenset(str(symbol).strip() for symbol in known_symbols if str(symbol).strip())
        )

    @staticmethod
    def _from_candidates(
        input_id: str,
        normalized_id: str,
        candidates: tuple[str, ...],
        source: str,
    ) -> GeneMappingRecord:
        if len(candidates) == 1:
            return GeneMappingRecord(
                input_id, normalized_id, candidates[0], source, "resolved", candidates
            )
        if len(candidates) > 1:
            return GeneMappingRecord(
                input_id, normalized_id, None, source, "ambiguous", candidates
            )
        return GeneMappingRecord(input_id, normalized_id, None, source, "unmapped", ())

    def resolve(self, gene_id: str) -> GeneMappingRecord:
        """Resolve one gene ID without silently selecting an ambiguous symbol."""

        input_id = str(gene_id)
        value = input_id.strip()
        if not value:
            return GeneMappingRecord(input_id, "", None, "none", "empty", ())

        normalized = strip_ensembl_version(value)
        if value.startswith("ENSG"):
            if not _ENSEMBL_GENE.fullmatch(value):
                return GeneMappingRecord(
                    input_id, normalized, None, "none", "invalid_ensembl", ()
                )
            if normalized in self._gencode:
                return self._from_candidates(
                    input_id, normalized, self._gencode[normalized], "gencode_v49"
                )
            if normalized in self._fallback:
                return self._from_candidates(
                    input_id,
                    normalized,
                    self._fallback[normalized],
                    "org_hs_eg_db_rank1",
                )
            return GeneMappingRecord(input_id, normalized, None, "none", "unmapped", ())

        if self._known_symbols is not None and value not in self._known_symbols:
            return GeneMappingRecord(input_id, value, None, "exact_symbol", "unmapped", ())
        return GeneMappingRecord(
            input_id, value, value, "exact_symbol", "resolved", (value,)
        )

    def resolve_many(self, gene_ids: Sequence[str]) -> tuple[GeneMappingRecord, ...]:
        """Resolve IDs in input order; records are suitable for an audit table."""

        return tuple(self.resolve(gene_id) for gene_id in gene_ids)


@dataclass(frozen=True)
class CollapsedExpression:
    """Expression matrix after unresolved features are dropped and duplicates merged."""

    matrix: ArrayLike
    symbols: tuple[str, ...]
    source_column_counts: tuple[int, ...]


def _shape_and_validate(matrix: ArrayLike, *, require_nonnegative: bool = True) -> None:
    if getattr(matrix, "ndim", None) != 2:
        raise InvalidExpressionError("Expression matrix must be two-dimensional")
    values = matrix.data if sparse.issparse(matrix) else np.asarray(matrix)
    if not np.issubdtype(values.dtype, np.number):
        raise InvalidExpressionError("Expression matrix must be numeric")
    if not np.all(np.isfinite(values)):
        raise InvalidExpressionError("Expression matrix contains NaN or infinite values")
    if require_nonnegative and np.any(values < 0):
        raise InvalidExpressionError("Expression layers must contain non-negative values")


def _to_linear(matrix: ArrayLike, layer: ExpressionLayer) -> ArrayLike:
    if sparse.issparse(matrix):
        result = sparse.csr_matrix(matrix, dtype=float, copy=True)
        if layer is ExpressionLayer.LOG1P_CPM:
            result.data = np.expm1(result.data)
        elif layer is ExpressionLayer.LOG2_TPM:
            result.data = np.exp2(result.data) - 1.0
        result.eliminate_zeros()
        return result

    result = np.asarray(matrix, dtype=float).copy()
    if layer is ExpressionLayer.LOG1P_CPM:
        return np.expm1(result)
    if layer is ExpressionLayer.LOG2_TPM:
        return np.exp2(result) - 1.0
    return result


def _from_linear(matrix: ArrayLike, layer: ExpressionLayer) -> ArrayLike:
    if layer not in {ExpressionLayer.LOG1P_CPM, ExpressionLayer.LOG2_TPM}:
        return matrix
    if sparse.issparse(matrix):
        result = sparse.csr_matrix(matrix, dtype=float, copy=True)
        if layer is ExpressionLayer.LOG1P_CPM:
            result.data = np.log1p(result.data)
        else:
            result.data = np.log2(1.0 + result.data)
        result.eliminate_zeros()
        return result
    if layer is ExpressionLayer.LOG1P_CPM:
        return np.log1p(np.asarray(matrix))
    return np.log2(1.0 + np.asarray(matrix))


def collapse_duplicate_symbols(
    matrix: ArrayLike,
    resolved_symbols: Sequence[str | None],
    layer: ExpressionLayer | str,
) -> CollapsedExpression:
    """Drop unresolved genes and sum duplicate symbols in linear space.

    Output symbols are sorted, so downstream results do not depend on feature
    order. Known log layers are inverted before summing and restored after it.
    """

    parsed_layer = parse_expression_layer(layer)
    _shape_and_validate(matrix)
    if matrix.shape[1] != len(resolved_symbols):
        raise InvalidExpressionError(
            "Number of resolved symbols must equal the number of matrix columns"
        )

    cleaned = tuple(
        None if symbol is None or not str(symbol).strip() else str(symbol).strip()
        for symbol in resolved_symbols
    )
    symbols = tuple(sorted({symbol for symbol in cleaned if symbol is not None}))
    symbol_index = {symbol: index for index, symbol in enumerate(symbols)}
    source_counts = Counter(symbol for symbol in cleaned if symbol is not None)
    counts = tuple(source_counts[symbol] for symbol in symbols)
    linear = _to_linear(matrix, parsed_layer)

    if sparse.issparse(linear):
        rows = [column for column, symbol in enumerate(cleaned) if symbol is not None]
        cols = [symbol_index[cleaned[column]] for column in rows]
        aggregator = sparse.csr_matrix(
            (np.ones(len(rows), dtype=float), (rows, cols)),
            shape=(matrix.shape[1], len(symbols)),
        )
        collapsed: ArrayLike = sparse.csr_matrix(linear) @ aggregator
    else:
        collapsed = np.zeros((matrix.shape[0], len(symbols)), dtype=float)
        for source_index, symbol in enumerate(cleaned):
            if symbol is not None:
                collapsed[:, symbol_index[symbol]] += linear[:, source_index]

    return CollapsedExpression(
        matrix=_from_linear(collapsed, parsed_layer),
        symbols=symbols,
        source_column_counts=counts,
    )


def _validate_symbols(matrix: ArrayLike, symbols: Sequence[str]) -> tuple[str, ...]:
    _shape_and_validate(matrix)
    cleaned = tuple(str(symbol).strip() for symbol in symbols)
    if matrix.shape[1] != len(cleaned):
        raise InvalidExpressionError("Number of symbols must equal matrix columns")
    if any(not symbol for symbol in cleaned):
        raise InvalidExpressionError("Resolved symbol names cannot be empty")
    if len(set(cleaned)) != len(cleaned):
        raise InvalidExpressionError(
            "Duplicate symbols must be collapsed before module scoring"
        )
    return cleaned


def _row_vector(value: ArrayLike | np.ndarray) -> np.ndarray:
    return np.asarray(value).reshape(-1).astype(float, copy=False)


@dataclass(frozen=True)
class LegacyModuleScore:
    """Legacy FM score and explicit gene-coverage diagnostics."""

    topic_mass_per_million: np.ndarray
    log1p_topic_mass_per_million: np.ndarray
    topic_mass: np.ndarray
    observed_genes: tuple[str, ...]
    missing_genes: tuple[str, ...]
    expected_gene_count: int
    observed_gene_count: int
    weight_coverage: float
    scorable: np.ndarray


def score_legacy_fm(
    matrix: ArrayLike,
    symbols: Sequence[str],
    weights: Mapping[str, float],
    *,
    layer: ExpressionLayer | str = ExpressionLayer.COUNTS,
    library_size: Sequence[float] | np.ndarray | None = None,
) -> LegacyModuleScore:
    """Score a frozen non-negative FM from raw counts.

    Observed member weights are renormalized to sum to one. The native score is
    topic mass divided by the complete row library size, multiplied by one
    million. Rows with zero library size are returned as NaN and ``scorable``
    is false.
    """

    parsed_layer = parse_expression_layer(layer)
    if parsed_layer is not ExpressionLayer.COUNTS:
        raise InvalidExpressionError("Legacy FM topic-mass scoring requires counts")
    cleaned_symbols = _validate_symbols(matrix, symbols)

    cleaned_weights: dict[str, float] = {}
    for raw_gene, raw_weight in weights.items():
        gene = str(raw_gene).strip()
        weight = float(raw_weight)
        if not gene or not math.isfinite(weight) or weight < 0:
            raise InvalidExpressionError(
                "Legacy FM weights require non-empty genes and finite non-negative values"
            )
        if weight > 0:
            cleaned_weights[gene] = weight
    if not cleaned_weights:
        raise InsufficientCoverageError("Legacy FM has no positive-weight genes")

    column_index = {symbol: index for index, symbol in enumerate(cleaned_symbols)}
    expected = tuple(sorted(cleaned_weights))
    observed = tuple(gene for gene in expected if gene in column_index)
    missing = tuple(gene for gene in expected if gene not in column_index)
    if not observed:
        raise InsufficientCoverageError("No legacy FM member is present in the matrix")

    observed_weight = sum(cleaned_weights[gene] for gene in observed)
    total_weight = sum(cleaned_weights.values())
    normalized = np.array(
        [cleaned_weights[gene] / observed_weight for gene in observed], dtype=float
    )
    subset = matrix[:, [column_index[gene] for gene in observed]]
    topic_mass = _row_vector(subset @ normalized)

    if library_size is None:
        library = _row_vector(matrix.sum(axis=1))
    else:
        library = np.asarray(library_size, dtype=float).reshape(-1)
        if library.shape != (matrix.shape[0],):
            raise InvalidExpressionError("library_size must contain one value per row")
        if not np.all(np.isfinite(library)) or np.any(library < 0):
            raise InvalidExpressionError("library_size must be finite and non-negative")

    scorable = library > 0
    topic_mass_per_million = np.full(matrix.shape[0], np.nan, dtype=float)
    np.divide(
        topic_mass * 1_000_000.0,
        library,
        out=topic_mass_per_million,
        where=scorable,
    )
    return LegacyModuleScore(
        topic_mass_per_million=topic_mass_per_million,
        log1p_topic_mass_per_million=np.log1p(topic_mass_per_million),
        topic_mass=topic_mass,
        observed_genes=observed,
        missing_genes=missing,
        expected_gene_count=len(expected),
        observed_gene_count=len(observed),
        weight_coverage=observed_weight / total_weight,
        scorable=scorable,
    )


def _as_log1p_expression(
    matrix: ArrayLike, layer: ExpressionLayer
) -> tuple[ArrayLike, np.ndarray]:
    """Convert supported layers to natural-log ``log1p`` expression."""

    _shape_and_validate(matrix)
    valid_rows = np.ones(matrix.shape[0], dtype=bool)
    if layer is ExpressionLayer.COUNTS:
        library = _row_vector(matrix.sum(axis=1))
        valid_rows = library > 0
        scale = np.zeros_like(library)
        np.divide(1_000_000.0, library, out=scale, where=valid_rows)
        if sparse.issparse(matrix):
            transformed = sparse.csr_matrix(matrix, dtype=float).multiply(scale[:, None])
            transformed = transformed.tocsr()
            transformed.data = np.log1p(transformed.data)
            transformed.eliminate_zeros()
            return transformed, valid_rows
        return np.log1p(np.asarray(matrix, dtype=float) * scale[:, None]), valid_rows

    if layer in {ExpressionLayer.CPM, ExpressionLayer.TPM}:
        if sparse.issparse(matrix):
            transformed = sparse.csr_matrix(matrix, dtype=float, copy=True)
            transformed.data = np.log1p(transformed.data)
            transformed.eliminate_zeros()
            return transformed, valid_rows
        return np.log1p(np.asarray(matrix, dtype=float)), valid_rows

    if layer is ExpressionLayer.LOG1P_CPM:
        return (
            sparse.csr_matrix(matrix, dtype=float, copy=True)
            if sparse.issparse(matrix)
            else np.asarray(matrix, dtype=float).copy(),
            valid_rows,
        )

    # log2(TPM + 1) -> ln(TPM + 1), preserving the same underlying abundance.
    factor = math.log(2.0)
    return (
        sparse.csr_matrix(matrix, dtype=float, copy=True) * factor
        if sparse.issparse(matrix)
        else np.asarray(matrix, dtype=float) * factor,
        valid_rows,
    )


@dataclass(frozen=True)
class CuratedSignatureScore:
    """Separate positive/negative component scores and optional signed contrast."""

    positive_score: np.ndarray | None
    negative_score: np.ndarray | None
    contrast: np.ndarray | None
    observed_positive_genes: tuple[str, ...]
    missing_positive_genes: tuple[str, ...]
    observed_negative_genes: tuple[str, ...]
    missing_negative_genes: tuple[str, ...]
    positive_coverage: float
    negative_coverage: float
    scorable: np.ndarray


def _clean_gene_set(genes: Iterable[str]) -> tuple[str, ...]:
    return tuple(sorted({str(gene).strip() for gene in genes if str(gene).strip()}))


def score_curated_signature(
    matrix: ArrayLike,
    symbols: Sequence[str],
    positive_genes: Iterable[str],
    negative_genes: Iterable[str] = (),
    *,
    layer: ExpressionLayer | str,
) -> CuratedSignatureScore:
    """Score observed signature genes without imputing absent genes as zero.

    Counts are library-normalized to CPM before ``log1p``. Linear CPM/TPM are
    transformed with ``log1p``; declared log layers are converted to the same
    natural-log scale. A signed contrast is emitted only when at least one gene
    from each declared component is observed.
    """

    parsed_layer = parse_expression_layer(layer)
    cleaned_symbols = _validate_symbols(matrix, symbols)
    positive = _clean_gene_set(positive_genes)
    negative = _clean_gene_set(negative_genes)
    if not positive and not negative:
        raise InsufficientCoverageError("Curated signature has no genes")
    overlap = set(positive).intersection(negative)
    if overlap:
        raise InvalidExpressionError(
            "Positive and negative components overlap: " + ",".join(sorted(overlap))
        )

    column_index = {symbol: index for index, symbol in enumerate(cleaned_symbols)}
    observed_positive = tuple(gene for gene in positive if gene in column_index)
    missing_positive = tuple(gene for gene in positive if gene not in column_index)
    observed_negative = tuple(gene for gene in negative if gene in column_index)
    missing_negative = tuple(gene for gene in negative if gene not in column_index)
    logged, valid_rows = _as_log1p_expression(matrix, parsed_layer)

    def component(observed: tuple[str, ...]) -> np.ndarray | None:
        if not observed:
            return None
        values = _row_vector(
            logged[:, [column_index[gene] for gene in observed]].mean(axis=1)
        )
        values[~valid_rows] = np.nan
        return values

    positive_score = component(observed_positive)
    negative_score = component(observed_negative)
    contrast = (
        positive_score - negative_score
        if positive_score is not None and negative_score is not None
        else None
    )
    scorable = valid_rows.copy()
    if positive and not observed_positive:
        scorable[:] = False
    if negative and not observed_negative:
        scorable[:] = False

    return CuratedSignatureScore(
        positive_score=positive_score,
        negative_score=negative_score,
        contrast=contrast,
        observed_positive_genes=observed_positive,
        missing_positive_genes=missing_positive,
        observed_negative_genes=observed_negative,
        missing_negative_genes=missing_negative,
        positive_coverage=(len(observed_positive) / len(positive) if positive else math.nan),
        negative_coverage=(len(observed_negative) / len(negative) if negative else math.nan),
        scorable=scorable,
    )


def robust_zscore(values: Sequence[float] | np.ndarray, *, axis: int = 0) -> np.ndarray:
    """Median/MAD z-score with deterministic dispersion fallbacks.

    If MAD is zero, the scale falls back to IQR, then standard deviation, then
    one for a constant vector. NaNs remain NaN and are excluded from estimates.
    """

    array = np.asarray(values, dtype=float)
    if array.ndim not in {1, 2}:
        raise InvalidExpressionError("robust_zscore supports one- or two-dimensional arrays")
    if axis not in {0, 1, -1}:
        raise InvalidExpressionError("axis must be 0 or 1")
    axis = axis % array.ndim
    median = np.nanmedian(array, axis=axis, keepdims=True)
    mad = np.nanmedian(np.abs(array - median), axis=axis, keepdims=True)
    scale = 1.4826 * mad
    q75 = np.nanpercentile(array, 75, axis=axis, keepdims=True)
    q25 = np.nanpercentile(array, 25, axis=axis, keepdims=True)
    iqr_scale = (q75 - q25) / 1.349
    std_scale = np.nanstd(array, axis=axis, keepdims=True)
    scale = np.where(scale > 0, scale, iqr_scale)
    scale = np.where(scale > 0, scale, std_scale)
    scale = np.where(scale > 0, scale, 1.0)
    return (array - median) / scale


def rank_percentile(
    values: Sequence[float] | np.ndarray, *, axis: int = 0
) -> np.ndarray:
    """Average-tie percentile ranks in ``(0, 1)`` with NaNs preserved."""

    array = np.asarray(values, dtype=float)
    if array.ndim not in {1, 2}:
        raise InvalidExpressionError("rank_percentile supports one- or two-dimensional arrays")
    if axis not in {0, 1, -1}:
        raise InvalidExpressionError("axis must be 0 or 1")
    axis = axis % array.ndim

    def rank_one(vector: np.ndarray) -> np.ndarray:
        result = np.full(vector.shape, np.nan, dtype=float)
        finite = np.isfinite(vector)
        count = int(finite.sum())
        if count:
            result[finite] = (rankdata(vector[finite], method="average") - 0.5) / count
        return result

    if array.ndim == 1:
        return rank_one(array)
    return np.apply_along_axis(rank_one, axis, array)


__all__ = [
    "CollapsedExpression",
    "CuratedSignatureScore",
    "ExpressionLayer",
    "GeneMappingRecord",
    "GeneResolver",
    "GeneResolutionError",
    "InsufficientCoverageError",
    "InvalidExpressionError",
    "LegacyModuleScore",
    "ScoringError",
    "UnknownExpressionLayerError",
    "collapse_duplicate_symbols",
    "parse_expression_layer",
    "rank_percentile",
    "robust_zscore",
    "score_curated_signature",
    "score_legacy_fm",
    "strip_ensembl_version",
]
