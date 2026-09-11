"""Unit tests for the fail-closed Stage 2 measurement primitives."""

from pathlib import Path
import sys
import unittest

import numpy as np
from scipy import sparse


ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

from v7.ontology.scoring import (  # noqa: E402
    ExpressionLayer,
    GeneResolver,
    InsufficientCoverageError,
    InvalidExpressionError,
    UnknownExpressionLayerError,
    collapse_duplicate_symbols,
    parse_expression_layer,
    rank_percentile,
    robust_zscore,
    score_curated_signature,
    score_legacy_fm,
)


class ExpressionLayerTest(unittest.TestCase):
    def test_only_explicit_layer_spellings_are_accepted(self) -> None:
        self.assertIs(parse_expression_layer("counts"), ExpressionLayer.COUNTS)
        self.assertIs(parse_expression_layer("CPM"), ExpressionLayer.CPM)
        with self.assertRaises(UnknownExpressionLayerError):
            parse_expression_layer("cpm")
        with self.assertRaises(UnknownExpressionLayerError):
            parse_expression_layer("normalized")


class GeneResolverTest(unittest.TestCase):
    def test_priority_version_stripping_fallback_and_ambiguity(self) -> None:
        resolver = GeneResolver(
            {
                "ENSG000001.9": "TP53",
                "ENSG000002": ("B", "A"),
                "ENSG000003": ("GENCODE",),
            },
            {
                "ENSG000003": "FALLBACK",
                "ENSG000004": "CD274",
            },
            known_symbols={"CXCL13"},
        )
        records = resolver.resolve_many(
            [
                "ENSG000001.2",
                "ENSG000002",
                "ENSG000003",
                "ENSG000004.7",
                "CXCL13",
                "NOT_A_KNOWN_SYMBOL",
                "ENSGbad",
            ]
        )
        self.assertEqual(records[0].resolved_symbol, "TP53")
        self.assertEqual(records[0].normalized_id, "ENSG000001")
        self.assertEqual(records[1].status, "ambiguous")
        self.assertIsNone(records[1].resolved_symbol)
        self.assertEqual(records[1].candidates, ("A", "B"))
        self.assertEqual(records[2].resolved_symbol, "GENCODE")
        self.assertEqual(records[2].source, "gencode_v49")
        self.assertEqual(records[3].resolved_symbol, "CD274")
        self.assertEqual(records[3].source, "org_hs_eg_db_rank1")
        self.assertEqual(records[4].source, "exact_symbol")
        self.assertEqual(records[5].status, "unmapped")
        self.assertEqual(records[6].status, "invalid_ensembl")


class DuplicateCollapseTest(unittest.TestCase):
    def test_linear_duplicates_are_summed_and_unresolved_columns_dropped(self) -> None:
        matrix = np.array([[1.0, 2.0, 50.0, 3.0], [4.0, 5.0, 60.0, 6.0]])
        result = collapse_duplicate_symbols(
            matrix, ["B", "A", None, "A"], ExpressionLayer.COUNTS
        )
        self.assertEqual(result.symbols, ("A", "B"))
        self.assertEqual(result.source_column_counts, (2, 1))
        np.testing.assert_allclose(result.matrix, [[5.0, 1.0], [11.0, 4.0]])

    def test_log_layers_are_collapsed_in_linear_space(self) -> None:
        linear = np.array([[2.0, 3.0, 7.0]])
        log1p = np.log1p(linear)
        collapsed = collapse_duplicate_symbols(
            log1p, ["A", "A", "B"], ExpressionLayer.LOG1P_CPM
        )
        np.testing.assert_allclose(collapsed.matrix, np.log1p([[5.0, 7.0]]))

        log2 = np.log2(linear + 1.0)
        collapsed_log2 = collapse_duplicate_symbols(
            sparse.csr_matrix(log2),
            ["A", "A", "B"],
            ExpressionLayer.LOG2_TPM,
        )
        self.assertTrue(sparse.issparse(collapsed_log2.matrix))
        np.testing.assert_allclose(
            collapsed_log2.matrix.toarray(), np.log2(np.array([[5.0, 7.0]]) + 1.0)
        )

    def test_invalid_matrix_and_unknown_layer_fail_closed(self) -> None:
        with self.assertRaises(InvalidExpressionError):
            collapse_duplicate_symbols(np.array([[1.0, -1.0]]), ["A", "B"], "counts")
        with self.assertRaises(UnknownExpressionLayerError):
            collapse_duplicate_symbols(np.ones((1, 1)), ["A"], "logcounts")


class LegacyFmTest(unittest.TestCase):
    def test_observed_weights_are_renormalized_and_missing_genes_audited(self) -> None:
        matrix = np.array([[10.0, 30.0, 60.0], [0.0, 20.0, 80.0], [0.0, 0.0, 0.0]])
        score = score_legacy_fm(
            matrix,
            ["A", "B", "OTHER"],
            {"A": 1.0, "B": 3.0, "MISSING": 6.0},
            layer="counts",
        )
        # Renormalized observed score: A*0.25 + B*0.75.
        np.testing.assert_allclose(score.topic_mass[:2], [25.0, 15.0])
        np.testing.assert_allclose(score.topic_mass_per_million[:2], [250_000.0, 150_000.0])
        self.assertTrue(np.isnan(score.topic_mass_per_million[2]))
        np.testing.assert_array_equal(score.scorable, [True, True, False])
        self.assertEqual(score.observed_genes, ("A", "B"))
        self.assertEqual(score.missing_genes, ("MISSING",))
        self.assertAlmostEqual(score.weight_coverage, 0.4)

    def test_dense_sparse_and_feature_order_are_equivalent(self) -> None:
        matrix = np.array([[2.0, 4.0, 10.0], [8.0, 1.0, 11.0]])
        dense = score_legacy_fm(matrix, ["A", "B", "X"], {"A": 1, "B": 2})
        sparse_score = score_legacy_fm(
            sparse.csr_matrix(matrix[:, [2, 1, 0]]),
            ["X", "B", "A"],
            {"B": 2, "A": 1},
        )
        np.testing.assert_allclose(
            dense.topic_mass_per_million, sparse_score.topic_mass_per_million
        )

    def test_noncount_or_zero_coverage_fails(self) -> None:
        with self.assertRaises(InvalidExpressionError):
            score_legacy_fm(np.ones((1, 1)), ["A"], {"A": 1}, layer="CPM")
        with self.assertRaises(InsufficientCoverageError):
            score_legacy_fm(np.ones((1, 1)), ["A"], {"B": 1})


class CuratedSignatureTest(unittest.TestCase):
    def test_missing_genes_are_excluded_and_contrast_requires_both_components(self) -> None:
        counts = np.array([[10.0, 30.0, 60.0], [0.0, 20.0, 80.0]])
        result = score_curated_signature(
            counts,
            ["A", "B", "OTHER"],
            ["A", "NOT_MEASURED"],
            ["B"],
            layer="counts",
        )
        expected_positive = np.log1p([100_000.0, 0.0])
        expected_negative = np.log1p([300_000.0, 200_000.0])
        np.testing.assert_allclose(result.positive_score, expected_positive)
        np.testing.assert_allclose(result.negative_score, expected_negative)
        np.testing.assert_allclose(result.contrast, expected_positive - expected_negative)
        self.assertEqual(result.observed_positive_genes, ("A",))
        self.assertEqual(result.missing_positive_genes, ("NOT_MEASURED",))
        self.assertEqual(result.positive_coverage, 0.5)

        one_sided = score_curated_signature(
            counts, ["A", "B", "OTHER"], ["A"], ["ABSENT"], layer="counts"
        )
        self.assertIsNone(one_sided.negative_score)
        self.assertIsNone(one_sided.contrast)
        np.testing.assert_array_equal(one_sided.scorable, [False, False])

    def test_dense_sparse_and_declared_log_layers_agree(self) -> None:
        tpm = np.array([[4.0, 8.0], [2.0, 0.0]])
        linear = score_curated_signature(
            tpm, ["A", "B"], ["A", "B"], layer="TPM"
        )
        log2_sparse = score_curated_signature(
            sparse.csr_matrix(np.log2(tpm + 1.0)),
            ["A", "B"],
            ["A", "B"],
            layer="log2_TPM",
        )
        np.testing.assert_allclose(linear.positive_score, log2_sparse.positive_score)

    def test_overlapping_components_and_duplicates_fail(self) -> None:
        with self.assertRaises(InvalidExpressionError):
            score_curated_signature(
                np.ones((1, 2)), ["A", "B"], ["A"], ["A"], layer="CPM"
            )
        with self.assertRaises(InvalidExpressionError):
            score_curated_signature(
                np.ones((1, 2)), ["A", "A"], ["A"], layer="CPM"
            )


class SensitivityTransformTest(unittest.TestCase):
    def test_robust_zscore_is_stable_for_constant_and_outlier_vectors(self) -> None:
        np.testing.assert_allclose(robust_zscore([2.0, 2.0, 2.0]), [0.0, 0.0, 0.0])
        scores = robust_zscore([1.0, 2.0, 3.0, 100.0])
        self.assertTrue(np.all(np.isfinite(scores)))
        self.assertLess(scores[0], scores[1])
        self.assertLess(scores[1], scores[2])
        self.assertLess(scores[2], scores[3])

    def test_rank_percentile_is_deterministic_with_ties_and_nans(self) -> None:
        values = np.array([5.0, 1.0, 5.0, np.nan])
        expected = np.array([2.0 / 3.0, 0.5 / 3.0, 2.0 / 3.0, np.nan])
        np.testing.assert_allclose(rank_percentile(values), expected, equal_nan=True)
        matrix = np.column_stack([values, values[::-1]])
        ranked = rank_percentile(matrix, axis=0)
        self.assertEqual(ranked.shape, matrix.shape)
        self.assertTrue(np.isnan(ranked[3, 0]))


if __name__ == "__main__":
    unittest.main()
