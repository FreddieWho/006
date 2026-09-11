"""Orchestration helpers for v7 Stage 2."""

from __future__ import annotations

from concurrent.futures import ProcessPoolExecutor, as_completed
import json
import math
from pathlib import Path
from typing import Any, Iterable

import numpy as np
import pandas as pd
from scipy import sparse

from .aggregate import AggregatedShard, load_shards
from .bridges import (
    build_gse207422_crosswalk,
    read_gse193736_bulk_counts,
    read_gse207422_bulk_log2tpm,
    read_gse207422_metadata,
)
from .contracts import Stage2Error, sha256, write_tsv
from .diagnostics import robust_standardize_scores
from .scoring import (
    ExpressionLayer,
    GeneResolver,
    collapse_duplicate_symbols,
    score_curated_signature,
    score_legacy_fm,
)
from .vocabulary import Stage2Vocabulary, compile_stage2_vocabulary


def compile_vocabulary(config: dict[str, Any]) -> Stage2Vocabulary:
    paths = config["inputs"]
    return compile_stage2_vocabulary(
        legacy_membership_path=paths["legacy_membership"],
        legacy_dictionary_path=paths["legacy_dictionary"],
        signature_registry_path=paths["curated_signature_registry"],
        hallmark_gmt_path=paths["hallmark_gmt"],
        reactome_gmt_path=paths["reactome_gmt"],
        sentinel_config_path=paths["sentinel_config"],
        ontology_path=paths["cell_state_ontology"],
        original_mapping_path=paths["original_label_mapping"],
        shadow_mapping_paths=[
            paths["marker_coarse_shadow"],
            paths["marker_mid_shadow"],
            paths["marker_fine_shadow"],
        ],
    )


def write_vocabulary(vocabulary: Stage2Vocabulary, output_root: Path) -> None:
    write_tsv(vocabulary.module_dictionary, output_root / "module_dictionary.tsv", ["feature_family", "parent_axis_id", "feature_id"])
    write_tsv(vocabulary.module_membership, output_root / "module_membership.tsv", ["feature_family", "parent_axis_id", "feature_id", "gene_symbol"])
    write_tsv(vocabulary.cell_state_ontology, output_root / "cell_state_ontology.tsv", ["node_kind", "level", "ontology_id"])
    write_tsv(vocabulary.label_mapping, output_root / "cell_state_label_mapping.tsv", ["mapping_status", "cohort_id", "annotation_field", "original_label"])


def load_written_vocabulary(output_root: Path) -> Stage2Vocabulary:
    return Stage2Vocabulary(
        module_dictionary=pd.read_csv(output_root / "module_dictionary.tsv", sep="\t"),
        module_membership=pd.read_csv(output_root / "module_membership.tsv", sep="\t"),
        cell_state_ontology=pd.read_csv(output_root / "cell_state_ontology.tsv", sep="\t"),
        label_mapping=pd.read_csv(output_root / "cell_state_label_mapping.tsv", sep="\t"),
    )


def _mapping_dictionary(frame: pd.DataFrame, id_column: str, symbol_column: str) -> dict[str, tuple[str, ...]]:
    source = frame[[id_column, symbol_column]].dropna().copy()
    source[id_column] = source[id_column].astype(str).str.replace(r"\.\d+$", "", regex=True)
    source[symbol_column] = source[symbol_column].astype(str).str.strip()
    return {
        key: tuple(sorted(set(group[symbol_column])))
        for key, group in source.groupby(id_column, sort=True)
        if key and group[symbol_column].ne("").any()
    }


def build_gene_resolver(config: dict[str, Any]) -> GeneResolver:
    gencode = pd.read_csv(
        config["inputs"]["gencode_gene_map"],
        header=None,
        names=["ENSEMBL", "SYMBOL"],
        dtype=str,
    )
    fallback = pd.read_csv(config["inputs"]["org_hs_gene_map"], dtype=str)
    fallback["mapping_rank"] = pd.to_numeric(fallback.mapping_rank, errors="coerce")
    fallback = fallback.loc[fallback.mapping_rank.eq(1)]
    return GeneResolver(
        _mapping_dictionary(gencode, "ENSEMBL", "SYMBOL"),
        _mapping_dictionary(fallback, "ENSEMBL", "SYMBOL"),
    )


def _natural_log_expression(matrix: Any, layer: ExpressionLayer) -> Any:
    if layer is ExpressionLayer.COUNTS:
        library = np.asarray(matrix.sum(axis=1)).reshape(-1).astype(float)
        scale = np.divide(1_000_000.0, library, out=np.zeros_like(library), where=library > 0)
        if sparse.issparse(matrix):
            result = sparse.csr_matrix(matrix, dtype=float).multiply(scale[:, None]).tocsr()
            result.data = np.log1p(result.data)
            return result
        return np.log1p(np.asarray(matrix, dtype=float) * scale[:, None])
    if layer in {ExpressionLayer.CPM, ExpressionLayer.TPM}:
        if sparse.issparse(matrix):
            result = sparse.csr_matrix(matrix, dtype=float, copy=True)
            result.data = np.log1p(result.data)
            return result
        return np.log1p(np.asarray(matrix, dtype=float))
    if layer is ExpressionLayer.LOG1P_CPM:
        return sparse.csr_matrix(matrix, dtype=float, copy=True) if sparse.issparse(matrix) else np.asarray(matrix, dtype=float)
    if layer is ExpressionLayer.LOG2_TPM:
        return matrix * math.log(2.0)
    raise Stage2Error(f"STAGE2_BLOCKED_LAYER: {layer}")


def _independent_legacy_count_formula(
    matrix: Any,
    symbols: tuple[str, ...],
    weights: dict[str, float],
    library_size: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Replay the frozen count formula without calling the production scorer."""

    position = {gene: index for index, gene in enumerate(symbols)}
    expected = {
        str(gene): float(weight)
        for gene, weight in weights.items()
        if np.isfinite(float(weight)) and float(weight) > 0
    }
    observed = [gene for gene in sorted(expected) if gene in position]
    if not observed:
        empty = np.full(matrix.shape[0], np.nan)
        return empty, empty.copy()
    normalized = np.asarray([expected[gene] for gene in observed], dtype=float)
    normalized /= normalized.sum()
    values = matrix[:, [position[gene] for gene in observed]]
    topic_mass = np.asarray(values @ normalized).reshape(-1).astype(float)
    per_million = np.divide(
        topic_mass * 1_000_000.0,
        np.asarray(library_size, dtype=float),
        out=np.full(matrix.shape[0], np.nan),
        where=np.asarray(library_size, dtype=float) > 0,
    )
    return topic_mass, np.log1p(per_million)


def _weighted_log_score(
    matrix: Any,
    symbols: tuple[str, ...],
    weights: dict[str, float],
    layer: ExpressionLayer,
    *,
    logged: Any | None = None,
) -> tuple[np.ndarray, float, int, int]:
    position = {gene: index for index, gene in enumerate(symbols)}
    expected = {str(gene): float(weight) for gene, weight in weights.items() if float(weight) > 0}
    observed = [gene for gene in sorted(expected) if gene in position]
    if not observed:
        return np.full(matrix.shape[0], np.nan), 0.0, len(expected), 0
    observed_total = sum(expected[gene] for gene in observed)
    normalized = np.asarray([expected[gene] / observed_total for gene in observed])
    if logged is None:
        logged = _natural_log_expression(matrix, layer)
    values = np.asarray(logged[:, [position[gene] for gene in observed]] @ normalized).reshape(-1)
    return values, observed_total / sum(expected.values()), len(expected), len(observed)


def _signature_median_sensitivity(
    matrix: Any,
    symbols: tuple[str, ...],
    positive: list[str],
    negative: list[str],
    layer: ExpressionLayer,
    *,
    logged: Any | None = None,
) -> np.ndarray:
    position = {gene: index for index, gene in enumerate(symbols)}
    if logged is None:
        logged = _natural_log_expression(matrix, layer)

    def median(genes: list[str]) -> np.ndarray | None:
        observed = [gene for gene in sorted(set(genes)) if gene in position]
        if not observed:
            return None
        subset = logged[:, [position[gene] for gene in observed]]
        values = subset.toarray() if sparse.issparse(subset) else np.asarray(subset)
        return np.median(values, axis=1)

    pos, neg = median(positive), median(negative)
    if pos is not None and neg is not None:
        return pos - neg
    if pos is not None:
        return pos
    if neg is not None:
        return -neg
    return np.full(matrix.shape[0], np.nan)


def score_expression_matrix(
    matrix: Any,
    raw_gene_ids: Iterable[str],
    metadata: pd.DataFrame,
    *,
    layer: ExpressionLayer,
    vocabulary: Stage2Vocabulary,
    resolver: GeneResolver,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    raw_gene_ids = tuple(map(str, raw_gene_ids))
    raw_library_size = (
        np.asarray(matrix.sum(axis=1)).reshape(-1).astype(float)
        if layer is ExpressionLayer.COUNTS
        else np.full(matrix.shape[0], np.nan)
    )
    mapping = resolver.resolve_many(raw_gene_ids)
    collapsed = collapse_duplicate_symbols(
        matrix,
        [record.resolved_symbol for record in mapping],
        layer,
    )
    symbols = collapsed.symbols
    matrix = collapsed.matrix
    if len(metadata) != matrix.shape[0]:
        raise Stage2Error("STAGE2_BLOCKED_SCORE_METADATA_SHAPE")
    audit = pd.DataFrame(
        [
            {
                "input_id": record.input_id,
                "normalized_id": record.normalized_id,
                "resolved_symbol": record.resolved_symbol or "",
                "mapping_source": record.source,
                "mapping_status": record.status,
                "candidates": "|".join(record.candidates),
            }
            for record in mapping
        ]
    )
    base_columns = [
        "expression_unit_id",
        "cohort_id",
        "study_family",
        "object_id",
        "sample_key",
        "patient_key",
        "timepoint",
        "tissue_context",
        "cell_state_level",
        "cell_state",
        "n_cells_used",
        "library_size",
        "modality",
        "platform",
        "aggregation_level",
    ]
    base = metadata.copy()
    for column in base_columns:
        if column not in base:
            base[column] = "" if column not in {"n_cells_used", "library_size"} else np.nan
    rows: list[pd.DataFrame] = []
    scoreable = vocabulary.module_dictionary.loc[vocabulary.module_dictionary.scoreable.astype(str).str.lower().isin({"true", "1", "yes"})]
    logged = _natural_log_expression(matrix, layer)
    valid_log_rows = (
        np.asarray(matrix.sum(axis=1)).reshape(-1).astype(float) > 0
        if layer is ExpressionLayer.COUNTS
        else np.ones(matrix.shape[0], dtype=bool)
    )
    symbol_position = {gene: index for index, gene in enumerate(symbols)}
    for feature in scoreable.itertuples(index=False):
        membership = vocabulary.module_membership.loc[vocabulary.module_membership.feature_id.eq(feature.feature_id)].copy()
        if membership.empty:
            raise Stage2Error(f"STAGE2_BLOCKED_EMPTY_MEMBERSHIP: {feature.feature_id}")
        expected = int(membership.gene_symbol.nunique())
        legacy_topic_mass = np.full(matrix.shape[0], np.nan)
        formula_replay_pass = np.full(matrix.shape[0], np.nan)
        formula_replay_topic_mass_pass = np.full(matrix.shape[0], np.nan)
        formula_replay_native_pass = np.full(matrix.shape[0], np.nan)
        formula_replay_sensitivity_pass = np.full(matrix.shape[0], np.nan)
        library_size_replay_pass = np.full(matrix.shape[0], np.nan)
        positive_coverage = np.nan
        negative_coverage = np.nan
        if feature.feature_family == "legacy_fm":
            weights = dict(zip(membership.gene_symbol.astype(str), pd.to_numeric(membership.relative_weight, errors="raise")))
            if layer is ExpressionLayer.COUNTS:
                library_size = pd.to_numeric(
                    base.library_size, errors="coerce"
                ).to_numpy(float)
                scored = score_legacy_fm(
                    matrix,
                    symbols,
                    weights,
                    layer=layer,
                    library_size=library_size,
                )
                native = scored.log1p_topic_mass_per_million
                legacy_topic_mass = scored.topic_mass
                direct_topic, direct_native = _independent_legacy_count_formula(
                    matrix, symbols, weights, library_size
                )
                formula_replay_topic_mass_pass = np.isclose(
                    direct_topic, scored.topic_mass, rtol=1e-12, atol=1e-12,
                    equal_nan=True,
                )
                formula_replay_native_pass = np.isclose(
                    direct_native, native, rtol=1e-12, atol=1e-12,
                    equal_nan=True,
                )
                library_size_replay_pass = np.isclose(
                    raw_library_size, library_size, rtol=1e-12, atol=1e-12,
                    equal_nan=False,
                )
                coverage = scored.weight_coverage
                observed = scored.observed_gene_count
                valid = scored.scorable & np.isfinite(native)
            else:
                native, coverage, expected, observed = _weighted_log_score(
                    matrix, symbols, weights, layer, logged=logged
                )
                valid = np.isfinite(native)
            top = membership.loc[
                membership.is_top_gene.astype(str)
                .str.lower()
                .isin({"yes", "true", "1"})
            ]
            top_weights = dict(zip(top.gene_symbol.astype(str), pd.to_numeric(top.relative_weight, errors="raise")))
            if top_weights:
                if layer is ExpressionLayer.COUNTS:
                    sensitivity = score_legacy_fm(
                        matrix,
                        symbols,
                        top_weights,
                        layer=layer,
                        library_size=library_size,
                    ).log1p_topic_mass_per_million
                    _, direct_sensitivity = _independent_legacy_count_formula(
                        matrix, symbols, top_weights, library_size
                    )
                    formula_replay_sensitivity_pass = np.isclose(
                        direct_sensitivity,
                        sensitivity,
                        rtol=1e-12,
                        atol=1e-12,
                        equal_nan=True,
                    )
                else:
                    sensitivity = _weighted_log_score(matrix, symbols, top_weights, layer)[0]
            else:
                sensitivity = np.full(matrix.shape[0], np.nan)
                if layer is ExpressionLayer.COUNTS:
                    formula_replay_sensitivity_pass = np.ones(
                        matrix.shape[0], dtype=bool
                    )
            if layer is ExpressionLayer.COUNTS:
                formula_replay_pass = (
                    formula_replay_topic_mass_pass.astype(bool)
                    & formula_replay_native_pass.astype(bool)
                    & formula_replay_sensitivity_pass.astype(bool)
                    & library_size_replay_pass.astype(bool)
                )
        else:
            positive = membership.loc[pd.to_numeric(membership.direction, errors="coerce").gt(0), "gene_symbol"].astype(str).tolist()
            negative = membership.loc[pd.to_numeric(membership.direction, errors="coerce").lt(0), "gene_symbol"].astype(str).tolist()
            observed_positive = sorted(set(positive).intersection(symbol_position))
            observed_negative = sorted(set(negative).intersection(symbol_position))

            def mean_component(genes: list[str]) -> np.ndarray | None:
                if not genes:
                    return None
                return np.asarray(
                    logged[:, [symbol_position[gene] for gene in genes]].mean(axis=1)
                ).reshape(-1)

            positive_score = mean_component(observed_positive)
            negative_score = mean_component(observed_negative)
            if positive_score is not None and negative_score is not None:
                native = positive_score - negative_score
            elif positive_score is not None:
                native = positive_score
            elif negative_score is not None:
                native = -negative_score
            else:
                native = np.full(matrix.shape[0], np.nan)
            observed = len(observed_positive) + len(observed_negative)
            coverage = observed / expected if expected else 0.0
            positive_coverage = len(observed_positive) / len(set(positive)) if positive else np.nan
            negative_coverage = len(observed_negative) / len(set(negative)) if negative else np.nan
            components_observed = (not positive or bool(observed_positive)) and (
                not negative or bool(observed_negative)
            )
            valid = valid_log_rows & components_observed & np.isfinite(native)
            sensitivity = _signature_median_sensitivity(
                matrix, symbols, positive, negative, layer, logged=logged
            )
        frame = base[base_columns].reset_index(drop=True).copy()
        frame["feature_family"] = feature.feature_family
        frame["feature_id"] = feature.feature_id
        frame["parent_axis_id"] = getattr(feature, "parent_axis_id", "")
        frame["score_native"] = np.asarray(native, dtype=float)
        frame["score_sensitivity"] = np.asarray(sensitivity, dtype=float)
        frame["legacy_topic_mass"] = legacy_topic_mass
        frame["formula_replay_pass"] = formula_replay_pass
        frame["formula_replay_topic_mass_pass"] = formula_replay_topic_mass_pass
        frame["formula_replay_native_pass"] = formula_replay_native_pass
        frame["formula_replay_sensitivity_pass"] = formula_replay_sensitivity_pass
        frame["library_size_replay_pass"] = library_size_replay_pass
        frame["expected_genes"] = expected
        frame["observed_genes"] = observed
        frame["coverage"] = coverage
        frame["positive_coverage"] = positive_coverage
        frame["negative_coverage"] = negative_coverage
        frame["technical_status"] = np.where(valid, "valid", "reject")
        frame["qc_flags"] = np.where(valid, "", "unscoreable_or_zero_library")
        rows.append(frame)
    scores = pd.concat(rows, ignore_index=True)
    return scores, audit


def _read_shard(shard: AggregatedShard) -> tuple[Any, list[str], pd.DataFrame]:
    matrix = sparse.load_npz(shard.matrix_path)
    genes = json.loads(shard.genes_path.read_text(encoding="utf-8"))["genes"]
    units = pd.read_parquet(shard.units_path)
    units["aggregation_level"] = "expression_unit"
    return matrix, genes, units


def _score_pseudobulk_shard(
    shard: AggregatedShard,
    vocabulary: Stage2Vocabulary,
    resolver: GeneResolver,
    shard_root: Path,
) -> tuple[str, Path, pd.DataFrame]:
    matrix, genes, units = _read_shard(shard)
    scores, audit = score_expression_matrix(
        matrix,
        genes,
        units,
        layer=ExpressionLayer.COUNTS,
        vocabulary=vocabulary,
        resolver=resolver,
    )
    safe = shard.object_id.replace("::", "__").replace("/", "_")
    path = shard_root / f"{safe}__{shard.level}.parquet"
    temporary = path.with_suffix(".parquet.tmp")
    scores.to_parquet(temporary, index=False)
    temporary.replace(path)
    audit["cohort_id"] = shard.cohort_id
    audit["object_id"] = shard.object_id
    return f"{shard.object_id}::{shard.level}", path, audit


def score_pseudobulk_shards(
    output_root: Path,
    vocabulary: Stage2Vocabulary,
    resolver: GeneResolver,
    *,
    workers: int = 1,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    shard_root = output_root / "cache" / "scores"
    shard_root.mkdir(parents=True, exist_ok=True)
    shards = load_shards(output_root)
    if workers < 1:
        raise Stage2Error("STAGE2_BLOCKED_SCORE_WORKERS: workers must be positive")
    results: list[tuple[str, Path, pd.DataFrame]] = []
    if workers == 1:
        results = [
            _score_pseudobulk_shard(shard, vocabulary, resolver, shard_root)
            for shard in shards
        ]
    else:
        with ProcessPoolExecutor(max_workers=workers) as executor:
            futures = {
                executor.submit(
                    _score_pseudobulk_shard,
                    shard,
                    vocabulary,
                    resolver,
                    shard_root,
                ): shard
                for shard in shards
            }
            for future in as_completed(futures):
                results.append(future.result())
    results.sort(key=lambda item: item[0])
    score_frames = [pd.read_parquet(path) for _, path, _ in results]
    audit_frames = [audit for _, _, audit in results]
    if not score_frames:
        raise Stage2Error("STAGE2_BLOCKED_NO_PSEUDOBULK_SHARDS")
    scores = robust_standardize_scores(pd.concat(score_frames, ignore_index=True))
    audit = pd.concat(audit_frames, ignore_index=True)
    grouped = (
        audit.groupby(
            ["input_id", "normalized_id", "resolved_symbol", "mapping_source", "mapping_status", "candidates"],
            as_index=False,
            dropna=False,
        )
        .agg(n_objects=("object_id", "nunique"), n_cohorts=("cohort_id", "nunique"))
    )
    return scores, grouped


def score_gse207422_whole_sample(
    output_root: Path,
    vocabulary: Stage2Vocabulary,
    resolver: GeneResolver,
    *,
    cohort_id: str = "GSE207422_sc",
) -> pd.DataFrame:
    """Resum cell-state count pseudobulks into whole-sample scRNA profiles."""

    frames: list[pd.DataFrame] = []
    for shard in load_shards(output_root):
        if shard.cohort_id != cohort_id or shard.level != "coarse":
            continue
        matrix, genes, units = _read_shard(shard)
        keys = units[["patient_key", "timepoint"]].astype(str)
        unique = keys.drop_duplicates().sort_values(
            ["patient_key", "timepoint"]
        ).reset_index(drop=True)
        key_to_code = {tuple(row): index for index, row in unique.iterrows()}
        codes = np.fromiter(
            (key_to_code[tuple(row)] for _, row in keys.iterrows()),
            dtype=np.int64,
            count=len(keys),
        )
        selector = sparse.csr_matrix(
            (np.ones(len(keys)), (codes, np.arange(len(keys)))),
            shape=(len(unique), len(keys)),
        )
        whole = sparse.csr_matrix(selector @ matrix)
        unique["expression_unit_id"] = (
            cohort_id
            + "::whole_sample::"
            + unique.patient_key.astype(str)
            + "::"
            + unique.timepoint.astype(str)
        )
        unique["cohort_id"] = cohort_id
        unique["study_family"] = "GSE207422"
        unique["object_id"] = shard.object_id
        unique["sample_key"] = unique.patient_key.astype(str) + "::" + unique.timepoint.astype(str)
        unique["tissue_context"] = "NSCLC_scRNA_whole_sample"
        unique["cell_state_level"] = "whole_sample"
        unique["cell_state"] = "whole_sample"
        unique["n_cells_used"] = np.bincount(
            codes,
            weights=pd.to_numeric(units.n_cells_used, errors="coerce").fillna(0).to_numpy(float),
            minlength=len(unique),
        ).astype(int)
        unique["library_size"] = np.asarray(whole.sum(axis=1)).reshape(-1)
        unique["modality"] = "scRNA_whole_sample"
        unique["platform"] = "scRNA"
        unique["aggregation_level"] = "whole_sample"
        scored, _ = score_expression_matrix(
            whole,
            genes,
            unique,
            layer=ExpressionLayer.COUNTS,
            vocabulary=vocabulary,
            resolver=resolver,
        )
        frames.append(scored)
    if not frames:
        raise Stage2Error("STAGE2_BLOCKED_GSE207422_WHOLE_SAMPLE_MISSING")
    result = pd.concat(frames, ignore_index=True)
    if result.duplicated(["patient_key", "timepoint", "feature_id"]).any():
        raise Stage2Error("STAGE2_BLOCKED_GSE207422_WHOLE_SAMPLE_DUPLICATE")
    return robust_standardize_scores(
        result,
        group_columns=("modality", "aggregation_level", "feature_id"),
    )


def score_bridges(
    config: dict[str, Any],
    vocabulary: Stage2Vocabulary,
    resolver: GeneResolver,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    inputs = config["inputs"]
    bulk_meta = read_gse207422_metadata(inputs["gse207422_bulk_metadata"])
    sc_meta = read_gse207422_metadata(inputs["gse207422_sc_metadata"])
    crosswalk = build_gse207422_crosswalk(sc_meta, bulk_meta)
    bulk = read_gse207422_bulk_log2tpm(inputs["gse207422_bulk_expression"])
    bulk_units = bulk_meta.set_index("sample_id").reindex(bulk.index).reset_index()
    if bulk_units.patient_id.isna().any():
        raise Stage2Error("STAGE2_BLOCKED_GSE207422_BULK_IDENTITY")
    bulk_units = bulk_units.rename(columns={"sample_id": "sample_key"})
    bulk_units["expression_unit_id"] = "GSE207422_bulk::" + bulk_units.sample_key.astype(str)
    bulk_units["cohort_id"] = "GSE207422_bulk"
    bulk_units["study_family"] = "GSE207422"
    bulk_units["object_id"] = "GSE207422_bulk::matrix"
    bulk_units["patient_key"] = "GSE207422::" + bulk_units.patient_id.astype(str)
    bulk_units["tissue_context"] = "NSCLC_bulk"
    bulk_units["cell_state_level"] = "bulk"
    bulk_units["cell_state"] = "bulk_tissue"
    bulk_units["modality"] = "bulkRNA"
    bulk_units["platform"] = "bulk_log2TPM"
    bulk_units["aggregation_level"] = "bulk_sample"
    bulk_units["n_cells_used"] = np.nan
    bulk_units["library_size"] = np.nan
    bulk_scores, bulk_audit = score_expression_matrix(
        bulk.to_numpy(float),
        bulk.columns.astype(str),
        bulk_units,
        layer=ExpressionLayer.LOG2_TPM,
        vocabulary=vocabulary,
        resolver=resolver,
    )

    perturb, design = read_gse193736_bulk_counts(inputs["gse193736_bulk_counts"])
    perturb_units = design.rename(columns={"sample_id": "sample_key"}).copy()
    perturb_units["expression_unit_id"] = "GSE193736::" + perturb_units.sample_key.astype(str)
    perturb_units["cohort_id"] = "GSE193736"
    perturb_units["study_family"] = "GSE193736"
    perturb_units["object_id"] = "GSE193736::bulkRNAseq_counts"
    perturb_units["patient_key"] = ""
    perturb_units["timepoint"] = perturb_units.rest_stim
    perturb_units["tissue_context"] = "primary_T_cell_perturbation"
    perturb_units["cell_state_level"] = "lineage"
    perturb_units["cell_state"] = perturb_units.lineage
    perturb_units["modality"] = "bulk_perturbation"
    perturb_units["platform"] = "bulkRNA_counts"
    perturb_units["aggregation_level"] = "perturbation_sample"
    perturb_units["n_cells_used"] = np.nan
    perturb_units["library_size"] = perturb.sum(axis=1).to_numpy(float)
    perturb_scores, perturb_audit = score_expression_matrix(
        sparse.csr_matrix(perturb.to_numpy(float)),
        perturb.columns.astype(str),
        perturb_units,
        layer=ExpressionLayer.COUNTS,
        vocabulary=vocabulary,
        resolver=resolver,
    )
    scores = robust_standardize_scores(pd.concat([bulk_scores, perturb_scores], ignore_index=True))
    audit = pd.concat(
        [bulk_audit.assign(cohort_id="GSE207422_bulk", object_id="GSE207422_bulk::matrix"), perturb_audit.assign(cohort_id="GSE193736", object_id="GSE193736::bulkRNAseq_counts")],
        ignore_index=True,
    )
    return scores, audit, crosswalk, design


def compatibility_audit(scores: pd.DataFrame, reference_path: Path) -> pd.DataFrame:
    reference = pd.read_parquet(reference_path)
    legacy = scores.loc[scores.feature_family.eq("legacy_fm")].pivot_table(
        index="expression_unit_id", columns="feature_id", values="legacy_topic_mass", aggfunc="first"
    )
    reference = reference.set_index("expression_unit_id")
    rows = []
    for feature_id in [f"FM{index:02d}" for index in range(1, 9)]:
        common = legacy.index.intersection(reference.index)
        left = pd.to_numeric(legacy.loc[common, feature_id], errors="coerce").to_numpy(float)
        right = pd.to_numeric(reference.loc[common, feature_id], errors="coerce").to_numpy(float)
        finite = np.isfinite(left) & np.isfinite(right)
        close = np.isclose(left[finite], right[finite], rtol=1e-8, atol=1e-10)
        feature_scores = scores.loc[
            scores.feature_id.eq(feature_id) & scores.modality.eq("scRNA")
        ]

        def all_true(column: str) -> bool:
            values = feature_scores[column].dropna()
            return bool(len(values) and values.astype(bool).all())

        topic_mass_pass = all_true("formula_replay_topic_mass_pass")
        native_pass = all_true("formula_replay_native_pass")
        sensitivity_pass = all_true("formula_replay_sensitivity_pass")
        library_size_pass = all_true("library_size_replay_pass")
        formula_pass = all_true("formula_replay_pass") and all(
            (topic_mass_pass, native_pass, sensitivity_pass, library_size_pass)
        )
        historical_exact = bool(finite.any() and close.all())
        rows.append(
            {
                "feature_id": feature_id,
                "n_reference_units": len(reference),
                "n_new_units": len(legacy),
                "n_common_units": len(common),
                "n_finite_compared": int(finite.sum()),
                "n_not_close": int((~close).sum()),
                "max_absolute_difference": float(np.max(np.abs(left[finite] - right[finite]))) if finite.any() else np.nan,
                "topic_mass_replay_pass": topic_mass_pass,
                "native_score_replay_pass": native_pass,
                "top_gene_sensitivity_replay_pass": sensitivity_pass,
                "library_size_replay_pass": library_size_pass,
                "formula_replay_pass": formula_pass,
                "historical_reference_exact": historical_exact,
                "context_mapping_status": (
                    "PASS_COMMON_CONTEXT_PRESENT"
                    if len(common) > 0
                    else "FAIL_NO_COMMON_CONTEXT"
                ),
                "status": (
                    "PASS_EXACT"
                    if formula_pass and historical_exact
                    else "PASS_FORMULA_REPLAY_HISTORICAL_CONTEXT_DRIFT_REMEASURED"
                    if formula_pass
                    else "FAIL_FORMULA_REPLAY"
                ),
            }
        )
    return pd.DataFrame(rows)


def family_state_gene_profiles(
    output_root: Path,
    resolver: GeneResolver,
) -> tuple[pd.DataFrame, np.ndarray, list[str], np.ndarray]:
    """Build family × mid-state profiles without encoding panel absence as zero."""

    pieces: list[tuple[pd.DataFrame, sparse.csr_matrix, tuple[str, ...]]] = []
    union: set[str] = set()
    for shard in load_shards(output_root):
        if shard.level != "mid":
            continue
        matrix, genes, units = _read_shard(shard)
        mapping = resolver.resolve_many(tuple(map(str, genes)))
        collapsed = collapse_duplicate_symbols(
            matrix,
            [record.resolved_symbol for record in mapping],
            ExpressionLayer.COUNTS,
        )
        logged = _natural_log_expression(collapsed.matrix, ExpressionLayer.COUNTS)
        logged = sparse.csr_matrix(logged)
        keys = units[["study_family", "cell_state"]].astype(str)
        unique = keys.drop_duplicates().sort_values(["study_family", "cell_state"]).reset_index(drop=True)
        key_to_code = {tuple(row): index for index, row in unique.iterrows()}
        codes = np.fromiter((key_to_code[tuple(row)] for _, row in keys.iterrows()), dtype=np.int64, count=len(keys))
        selector = sparse.csr_matrix(
            (np.ones(len(keys)), (codes, np.arange(len(keys)))),
            shape=(len(unique), len(keys)),
        )
        sums = sparse.csr_matrix(selector @ logged)
        counts = np.bincount(codes, minlength=len(unique)).astype(float)
        unique["n_expression_units"] = counts.astype(int)
        pieces.append((unique, sums, collapsed.symbols))
        union.update(collapsed.symbols)
    if not pieces:
        raise Stage2Error("STAGE2_BLOCKED_NO_MID_STATE_PROFILES")
    symbols = sorted(union)
    position = {gene: index for index, gene in enumerate(symbols)}
    metadata_parts: list[pd.DataFrame] = []
    numerator_parts: list[sparse.csr_matrix] = []
    denominator_parts: list[sparse.csr_matrix] = []
    for metadata, matrix, local_symbols in pieces:
        mapper = sparse.csr_matrix(
            (
                np.ones(len(local_symbols)),
                (np.arange(len(local_symbols)), [position[gene] for gene in local_symbols]),
            ),
            shape=(len(local_symbols), len(symbols)),
        )
        metadata_parts.append(metadata)
        numerator_parts.append(matrix @ mapper)
        local_observed = sparse.csr_matrix(
            np.repeat(
                metadata.n_expression_units.to_numpy(float)[:, None],
                len(local_symbols),
                axis=1,
            )
        )
        denominator_parts.append(local_observed @ mapper)
    metadata = pd.concat(metadata_parts, ignore_index=True)
    numerator_stack = sparse.vstack(numerator_parts, format="csr")
    denominator_stack = sparse.vstack(denominator_parts, format="csr")
    keys = metadata[["study_family", "cell_state"]].astype(str)
    unique = keys.drop_duplicates().sort_values(["study_family", "cell_state"]).reset_index(drop=True)
    key_to_code = {tuple(row): index for index, row in unique.iterrows()}
    codes = np.fromiter((key_to_code[tuple(row)] for _, row in keys.iterrows()), dtype=np.int64, count=len(keys))
    selector = sparse.csr_matrix(
        (np.ones(len(metadata)), (codes, np.arange(len(metadata)))),
        shape=(len(unique), len(metadata)),
    )
    numerator = (selector @ numerator_stack).toarray()
    denominator = (selector @ denominator_stack).toarray()
    profile = np.divide(
        numerator,
        denominator,
        out=np.full(numerator.shape, np.nan, dtype=float),
        where=denominator > 0,
    )
    total_units = np.bincount(codes, weights=metadata.n_expression_units.to_numpy(float), minlength=len(unique))
    unique["n_expression_units"] = total_units.astype(int)
    return unique, profile, symbols, denominator > 0
