"""Streaming cell-level technical validation for v7 Stage 2.

The production measurement is count pseudobulk.  This module independently
scores individual cells for a small, named technical panel and compares the
mean cell score within each expression unit with the score obtained after
count aggregation.  The two estimands are intentionally not treated as
numerically equivalent because log normalisation and aggregation do not
commute.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable

import numpy as np
import pandas as pd
from scipy import sparse

from .aggregate import ObjectSpec, _identity_table, _raw_gene_ids
from .contracts import Stage2Error
from .diagnostics import spearman_correlation
from .scoring import GeneResolver
from .vocabulary import Stage2Vocabulary


@dataclass(frozen=True)
class _FeaturePlan:
    feature_family: str
    feature_id: str
    legacy_indices: tuple[int, ...]
    legacy_weights: tuple[float, ...]
    positive_indices: tuple[int, ...]
    negative_indices: tuple[int, ...]
    coverage: float
    scoreable: bool


def _resolved_aggregator(
    raw_gene_ids: Iterable[str], resolver: GeneResolver
) -> tuple[sparse.csr_matrix, tuple[str, ...]]:
    records = resolver.resolve_many(tuple(map(str, raw_gene_ids)))
    resolved = tuple(record.resolved_symbol for record in records)
    symbols = tuple(sorted({symbol for symbol in resolved if symbol}))
    if not symbols:
        raise Stage2Error("STAGE2_BLOCKED_TECHNICAL_PANEL_NO_RESOLVED_GENES")
    position = {symbol: index for index, symbol in enumerate(symbols)}
    rows = [index for index, symbol in enumerate(resolved) if symbol]
    columns = [position[resolved[index]] for index in rows]
    transform = sparse.csr_matrix(
        (np.ones(len(rows), dtype=float), (rows, columns)),
        shape=(len(resolved), len(symbols)),
    )
    return transform, symbols


def _feature_plans(
    vocabulary: Stage2Vocabulary, symbols: tuple[str, ...]
) -> tuple[_FeaturePlan, ...]:
    position = {symbol: index for index, symbol in enumerate(symbols)}
    scoreable = vocabulary.module_dictionary.loc[
        vocabulary.module_dictionary.scoreable.astype(str).str.lower().isin(
            {"true", "1", "yes"}
        )
    ]
    plans: list[_FeaturePlan] = []
    for feature in scoreable.itertuples(index=False):
        membership = vocabulary.module_membership.loc[
            vocabulary.module_membership.feature_id.eq(feature.feature_id)
        ].copy()
        if membership.empty:
            raise Stage2Error(
                f"STAGE2_BLOCKED_EMPTY_MEMBERSHIP: {feature.feature_id}"
            )
        if feature.feature_family == "legacy_fm":
            weights = dict(
                zip(
                    membership.gene_symbol.astype(str),
                    pd.to_numeric(membership.relative_weight, errors="raise"),
                )
            )
            observed = [gene for gene in sorted(weights) if gene in position]
            observed_weight = float(sum(weights[gene] for gene in observed))
            total_weight = float(sum(weights.values()))
            normalized = (
                tuple(weights[gene] / observed_weight for gene in observed)
                if observed_weight > 0
                else ()
            )
            plans.append(
                _FeaturePlan(
                    feature_family=str(feature.feature_family),
                    feature_id=str(feature.feature_id),
                    legacy_indices=tuple(position[gene] for gene in observed),
                    legacy_weights=normalized,
                    positive_indices=(),
                    negative_indices=(),
                    coverage=(observed_weight / total_weight if total_weight > 0 else 0.0),
                    scoreable=bool(observed),
                )
            )
            continue

        positive = sorted(
            set(
                membership.loc[
                    pd.to_numeric(membership.direction, errors="coerce").gt(0),
                    "gene_symbol",
                ].astype(str)
            )
        )
        negative = sorted(
            set(
                membership.loc[
                    pd.to_numeric(membership.direction, errors="coerce").lt(0),
                    "gene_symbol",
                ].astype(str)
            )
        )
        if set(positive).intersection(negative):
            raise Stage2Error(
                f"STAGE2_BLOCKED_OVERLAPPING_SIGNATURE: {feature.feature_id}"
            )
        observed_positive = [gene for gene in positive if gene in position]
        observed_negative = [gene for gene in negative if gene in position]
        expected = len(positive) + len(negative)
        observed = len(observed_positive) + len(observed_negative)
        valid = (not positive or bool(observed_positive)) and (
            not negative or bool(observed_negative)
        )
        plans.append(
            _FeaturePlan(
                feature_family=str(feature.feature_family),
                feature_id=str(feature.feature_id),
                legacy_indices=(),
                legacy_weights=(),
                positive_indices=tuple(position[gene] for gene in observed_positive),
                negative_indices=tuple(position[gene] for gene in observed_negative),
                coverage=(observed / expected if expected else 0.0),
                scoreable=valid,
            )
        )
    return tuple(plans)


def _group_codes(
    identity: pd.DataFrame, spec: ObjectSpec, level: str
) -> tuple[np.ndarray, pd.DataFrame]:
    eligible = (
        identity._resolved
        & identity.analysis_unit_key.notna()
        & identity.patient_key.astype(str).str.strip().ne("")
        & identity.analysis_unit_key.astype(str).str.strip().ne("")
        & identity[level].astype(str).str.strip().ne("")
    )
    key_columns = ["analysis_unit_key", "patient_key", "timepoint", "tissue_context", level]
    keys = (
        identity.loc[eligible, key_columns]
        .astype(str)
        .drop_duplicates()
        .sort_values(key_columns)
        .reset_index(drop=True)
    )
    if keys.empty:
        raise Stage2Error(
            f"STAGE2_BLOCKED_TECHNICAL_PANEL_NO_ELIGIBLE_CELLS: {spec.object_id}:{level}"
        )
    key_to_code = {tuple(row): index for index, row in keys.iterrows()}
    codes = np.full(len(identity), -1, dtype=np.int64)
    eligible_rows = identity.index[eligible].to_numpy()
    codes[eligible_rows] = np.fromiter(
        (
            key_to_code[tuple(row)]
            for _, row in identity.loc[eligible, key_columns].astype(str).iterrows()
        ),
        dtype=np.int64,
        count=len(eligible_rows),
    )
    keys = keys.rename(columns={level: "cell_state"})
    keys["cell_state_level"] = level
    keys["expression_unit_id"] = [
        f"{spec.cohort_id}::{spec.object_id}::{row.analysis_unit_key}::{level}::{row.cell_state}::count_pseudobulk"
        for row in keys.itertuples()
    ]
    return codes, keys


def _log1p_cpm(matrix: sparse.csr_matrix) -> tuple[sparse.csr_matrix, np.ndarray]:
    library = np.asarray(matrix.sum(axis=1)).reshape(-1).astype(float)
    scale = np.divide(
        1_000_000.0,
        library,
        out=np.zeros_like(library),
        where=library > 0,
    )
    logged = matrix.astype(float, copy=False).multiply(scale[:, None]).tocsr()
    logged.data = np.log1p(logged.data)
    logged.eliminate_zeros()
    return logged, library > 0


def _score_chunk(
    collapsed: sparse.csr_matrix,
    raw_library: np.ndarray,
    plans: tuple[_FeaturePlan, ...],
) -> np.ndarray:
    logged, resolved_library_valid = _log1p_cpm(collapsed)
    values = np.full((collapsed.shape[0], len(plans)), np.nan, dtype=float)
    for feature_index, plan in enumerate(plans):
        if not plan.scoreable:
            continue
        if plan.feature_family == "legacy_fm":
            topic_mass = np.asarray(
                collapsed[:, list(plan.legacy_indices)]
                @ np.asarray(plan.legacy_weights, dtype=float)
            ).reshape(-1)
            valid = raw_library > 0
            per_million = np.divide(
                topic_mass * 1_000_000.0,
                raw_library,
                out=np.full_like(topic_mass, np.nan, dtype=float),
                where=valid,
            )
            values[:, feature_index] = np.log1p(per_million)
            continue

        positive = (
            np.asarray(logged[:, list(plan.positive_indices)].mean(axis=1)).reshape(-1)
            if plan.positive_indices
            else None
        )
        negative = (
            np.asarray(logged[:, list(plan.negative_indices)].mean(axis=1)).reshape(-1)
            if plan.negative_indices
            else None
        )
        if positive is not None and negative is not None:
            score = positive - negative
        elif positive is not None:
            score = positive
        elif negative is not None:
            score = -negative
        else:
            continue
        score[~resolved_library_valid] = np.nan
        values[:, feature_index] = score
    return values


def summarize_cell_vs_pseudobulk(
    unit_comparison: pd.DataFrame,
) -> pd.DataFrame:
    """Summarize rank agreement without asserting estimator equivalence."""

    rows: list[dict[str, Any]] = []
    group_columns = [
        "cohort_id",
        "object_id",
        "cell_state_level",
        "feature_family",
        "feature_id",
    ]
    for key, frame in unit_comparison.groupby(group_columns, sort=True, dropna=False):
        comparable = frame.dropna(subset=["cell_mean_score", "pseudobulk_score"])
        correlation = spearman_correlation(
            comparable.cell_mean_score, comparable.pseudobulk_score
        )
        row = dict(zip(group_columns, key))
        row.update(
            {
                "n_expression_units": int(len(frame)),
                "n_comparable_units": int(len(comparable)),
                "n_cells_total": int(frame.n_cells_total.sum()),
                "n_cells_scoreable": int(frame.n_cells_scoreable.sum()),
                "cell_scoreable_fraction": (
                    float(frame.n_cells_scoreable.sum() / frame.n_cells_total.sum())
                    if frame.n_cells_total.sum() > 0
                    else np.nan
                ),
                "cell_mean_vs_pseudobulk_spearman": correlation,
                "comparison_status": (
                    "OBSERVED_RANK_CONCORDANCE"
                    if np.isfinite(correlation)
                    else "NOT_TESTABLE_CONSTANT_OR_FEW_UNITS"
                ),
                "estimand_relation": "EXPECTED_NON_EQUIVALENCE_LOG_NORMALIZATION_BEFORE_VS_AFTER_AGGREGATION",
            }
        )
        rows.append(row)
    return pd.DataFrame(rows)


def score_object_cell_level_panel(
    spec: ObjectSpec,
    vocabulary: Stage2Vocabulary,
    resolver: GeneResolver,
    pseudobulk_scores: pd.DataFrame,
    *,
    chunk_size: int = 50_000,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Run an all-cell, no-cell-export technical comparison for one object."""

    import anndata as ad

    adata = ad.read_h5ad(spec.h5ad, backed="r")
    try:
        if spec.count_layer != "X" and spec.count_layer not in adata.layers:
            raise Stage2Error(
                f"STAGE2_BLOCKED_MISSING_COUNT_LAYER: {spec.object_id}"
            )
        source = adata.X if spec.count_layer == "X" else adata.layers[spec.count_layer]
        identity, _ = _identity_table(spec, adata.obs)
        transform, symbols = _resolved_aggregator(
            _raw_gene_ids(spec, adata).astype(str), resolver
        )
        plans = _feature_plans(vocabulary, symbols)
        group_data = {
            level: _group_codes(identity, spec, level) for level in ("coarse", "mid")
        }
        selected_rows = np.flatnonzero(
            np.logical_or(
                group_data["coarse"][0] >= 0,
                group_data["mid"][0] >= 0,
            )
        )
        accumulators: dict[str, tuple[np.ndarray, np.ndarray]] = {}
        for level, (_, keys) in group_data.items():
            shape = (len(keys), len(plans))
            accumulators[level] = (np.zeros(shape, dtype=float), np.zeros(shape, dtype=np.int64))

        for start in range(0, len(selected_rows), chunk_size):
            rows = selected_rows[start : start + chunk_size]
            raw = source[rows]
            raw = raw.tocsr() if sparse.issparse(raw) else sparse.csr_matrix(raw)
            raw = raw.astype(float, copy=False)
            if raw.data.size:
                if np.any(raw.data < 0) or not np.all(np.isfinite(raw.data)):
                    raise Stage2Error(
                        f"STAGE2_BLOCKED_INVALID_COUNTS: {spec.object_id}"
                    )
                if not np.allclose(raw.data, np.rint(raw.data), rtol=0, atol=1e-6):
                    raise Stage2Error(
                        f"STAGE2_BLOCKED_LAYER_SEMANTICS: {spec.object_id}"
                    )
            raw_library = np.asarray(raw.sum(axis=1)).reshape(-1).astype(float)
            collapsed = (raw @ transform).tocsr()
            chunk_scores = _score_chunk(collapsed, raw_library, plans)
            finite = np.isfinite(chunk_scores)
            for level, (all_codes, keys) in group_data.items():
                codes = all_codes[rows]
                eligible = codes >= 0
                if not eligible.any():
                    continue
                local_codes = codes[eligible]
                selector = sparse.csr_matrix(
                    (
                        np.ones(int(eligible.sum()), dtype=float),
                        (local_codes, np.arange(int(eligible.sum()))),
                    ),
                    shape=(len(keys), int(eligible.sum())),
                )
                sums, counts = accumulators[level]
                observed = chunk_scores[eligible]
                observed_finite = finite[eligible]
                sums += np.asarray(selector @ np.where(observed_finite, observed, 0.0))
                counts += np.asarray(selector @ observed_finite.astype(np.int64)).astype(
                    np.int64
                )

        comparison_frames: list[pd.DataFrame] = []
        for level, (codes, keys) in group_data.items():
            sums, counts = accumulators[level]
            total_cells = np.bincount(codes[codes >= 0], minlength=len(keys))
            means = np.divide(
                sums,
                counts,
                out=np.full_like(sums, np.nan),
                where=counts > 0,
            )
            for feature_index, plan in enumerate(plans):
                frame = keys[["expression_unit_id"]].copy()
                frame["cohort_id"] = spec.cohort_id
                frame["object_id"] = spec.object_id
                frame["cell_state_level"] = level
                frame["feature_family"] = plan.feature_family
                frame["feature_id"] = plan.feature_id
                frame["coverage"] = plan.coverage
                frame["n_cells_total"] = total_cells
                frame["n_cells_scoreable"] = counts[:, feature_index]
                frame["cell_mean_score"] = means[:, feature_index]
                comparison_frames.append(frame)

        comparison = pd.concat(comparison_frames, ignore_index=True)
        direct = pseudobulk_scores.loc[
            pseudobulk_scores.cohort_id.eq(spec.cohort_id)
            & pseudobulk_scores.object_id.eq(spec.object_id)
            & pseudobulk_scores.modality.eq("scRNA"),
            ["expression_unit_id", "feature_id", "score_native"],
        ].rename(columns={"score_native": "pseudobulk_score"})
        if direct.duplicated(["expression_unit_id", "feature_id"]).any():
            raise Stage2Error(
                f"STAGE2_BLOCKED_TECHNICAL_PANEL_DUPLICATE_DIRECT_SCORE: {spec.object_id}"
            )
        comparison = comparison.merge(
            direct,
            on=["expression_unit_id", "feature_id"],
            how="left",
            validate="one_to_one",
        )
        if comparison.pseudobulk_score.notna().sum() == 0:
            raise Stage2Error(
                f"STAGE2_BLOCKED_TECHNICAL_PANEL_NO_DIRECT_MATCH: {spec.object_id}"
            )
        return summarize_cell_vs_pseudobulk(comparison), comparison
    finally:
        if getattr(adata, "file", None) is not None:
            adata.file.close()


def run_cell_level_technical_panel(
    specs: Iterable[ObjectSpec],
    panel_cohorts: Iterable[str],
    vocabulary: Stage2Vocabulary,
    resolver: GeneResolver,
    pseudobulk_scores: pd.DataFrame,
    *,
    chunk_size: int = 50_000,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    requested = tuple(map(str, panel_cohorts))
    selected = [spec for spec in specs if spec.cohort_id in requested]
    observed = {spec.cohort_id for spec in selected}
    missing = sorted(set(requested) - observed)
    if missing:
        raise Stage2Error(
            "STAGE2_BLOCKED_TECHNICAL_PANEL_COHORTS: " + ",".join(missing)
        )
    summaries: list[pd.DataFrame] = []
    comparisons: list[pd.DataFrame] = []
    for spec in selected:
        summary, comparison = score_object_cell_level_panel(
            spec,
            vocabulary,
            resolver,
            pseudobulk_scores,
            chunk_size=chunk_size,
        )
        summaries.append(summary)
        comparisons.append(comparison)
    return (
        pd.concat(summaries, ignore_index=True).sort_values(
            ["cohort_id", "object_id", "cell_state_level", "feature_id"]
        ).reset_index(drop=True),
        pd.concat(comparisons, ignore_index=True).sort_values(
            ["cohort_id", "object_id", "cell_state_level", "feature_id", "expression_unit_id"]
        ).reset_index(drop=True),
    )


__all__ = [
    "run_cell_level_technical_panel",
    "score_object_cell_level_panel",
    "summarize_cell_vs_pseudobulk",
]
