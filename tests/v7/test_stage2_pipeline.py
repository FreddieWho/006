from __future__ import annotations

import json

import numpy as np
import pandas as pd
from scipy import sparse

from src.v7.ontology.pipeline import (
    family_state_gene_profiles,
    score_expression_matrix,
)
from src.v7.ontology.scoring import ExpressionLayer, GeneResolver, score_legacy_fm
from src.v7.ontology.vocabulary import Stage2Vocabulary


def test_pipeline_recognizes_boolean_top_genes_for_legacy_sensitivity() -> None:
    matrix = sparse.csr_matrix(
        np.asarray([[10, 0, 5], [1, 4, 8], [3, 2, 1]], dtype=float)
    )
    dictionary = pd.DataFrame(
        {
            "feature_family": ["legacy_fm"],
            "feature_id": ["FM01"],
            "parent_axis_id": [""],
            "scoreable": [True],
        }
    )
    membership = pd.DataFrame(
        {
            "feature_family": ["legacy_fm"] * 3,
            "feature_id": ["FM01"] * 3,
            "gene_symbol": ["A", "B", "C"],
            "relative_weight": [0.5, 0.3, 0.2],
            "weight": [0.5, 0.3, 0.2],
            "direction": [1, 1, 1],
            "is_top_gene": [True, True, False],
        }
    )
    vocabulary = Stage2Vocabulary(
        dictionary, membership, pd.DataFrame(), pd.DataFrame()
    )
    metadata = pd.DataFrame(
        {
            "expression_unit_id": ["u1", "u2", "u3"],
            "library_size": np.asarray(matrix.sum(axis=1)).reshape(-1),
        }
    )

    scores, _ = score_expression_matrix(
        matrix,
        ["A", "B", "C"],
        metadata,
        layer=ExpressionLayer.COUNTS,
        vocabulary=vocabulary,
        resolver=GeneResolver({}, {}),
    )
    expected = score_legacy_fm(
        matrix,
        ["A", "B", "C"],
        {"A": 0.5, "B": 0.3},
        layer=ExpressionLayer.COUNTS,
        library_size=metadata.library_size,
    ).log1p_topic_mass_per_million

    assert np.allclose(scores.score_sensitivity, expected)


def test_legacy_formula_replay_covers_native_sensitivity_and_library_size() -> None:
    matrix = sparse.csr_matrix(
        np.asarray([[10, 0, 5], [1, 4, 8], [3, 2, 1]], dtype=float)
    )
    dictionary = pd.DataFrame(
        {
            "feature_family": ["legacy_fm"],
            "feature_id": ["FM01"],
            "parent_axis_id": [""],
            "scoreable": [True],
        }
    )
    membership = pd.DataFrame(
        {
            "feature_family": ["legacy_fm"] * 3,
            "feature_id": ["FM01"] * 3,
            "gene_symbol": ["A", "B", "C"],
            "relative_weight": [0.5, 0.3, 0.2],
            "weight": [0.5, 0.3, 0.2],
            "direction": [1, 1, 1],
            "is_top_gene": [True, True, False],
        }
    )
    vocabulary = Stage2Vocabulary(
        dictionary, membership, pd.DataFrame(), pd.DataFrame()
    )
    true_library = np.asarray(matrix.sum(axis=1)).reshape(-1)
    metadata = pd.DataFrame(
        {
            "expression_unit_id": ["u1", "u2", "u3"],
            "library_size": true_library,
        }
    )
    scores, _ = score_expression_matrix(
        matrix,
        ["A", "B", "C"],
        metadata,
        layer=ExpressionLayer.COUNTS,
        vocabulary=vocabulary,
        resolver=GeneResolver({}, {}),
    )
    replay_columns = [
        "formula_replay_topic_mass_pass",
        "formula_replay_native_pass",
        "formula_replay_sensitivity_pass",
        "library_size_replay_pass",
        "formula_replay_pass",
    ]
    assert scores[replay_columns].astype(bool).all().all()

    corrupted = metadata.assign(library_size=true_library + 1)
    corrupted_scores, _ = score_expression_matrix(
        matrix,
        ["A", "B", "C"],
        corrupted,
        layer=ExpressionLayer.COUNTS,
        vocabulary=vocabulary,
        resolver=GeneResolver({}, {}),
    )
    assert not corrupted_scores.library_size_replay_pass.astype(bool).any()
    assert not corrupted_scores.formula_replay_pass.astype(bool).any()


def test_family_state_profiles_preserve_unmeasured_genes_as_nan(tmp_path) -> None:
    cache = tmp_path / "cache" / "pseudobulk"
    cache.mkdir(parents=True)
    specifications = [
        ("C1", "O1", "F1", ["A", "B"], [[10, 1], [5, 4]]),
        ("C2", "O2", "F2", ["A", "C"], [[8, 3], [2, 9]]),
    ]
    for cohort, object_id, family, genes, values in specifications:
        base = cache / f"{object_id}__mid"
        matrix_path = base.with_suffix(".counts.npz")
        units_path = base.with_suffix(".units.parquet")
        genes_path = base.with_suffix(".genes.json")
        audit_path = base.with_suffix(".audit.json")
        sparse.save_npz(matrix_path, sparse.csr_matrix(values))
        pd.DataFrame(
            {
                "study_family": [family, family],
                "cell_state": ["S1", "S2"],
            }
        ).to_parquet(units_path, index=False)
        genes_path.write_text(json.dumps({"genes": genes}), encoding="utf-8")
        audit_path.write_text(
            json.dumps(
                {
                    "cohort_id": cohort,
                    "object_id": object_id,
                    "cell_state_level": "mid",
                    "matrix_path": str(matrix_path),
                    "units_path": str(units_path),
                    "genes_path": str(genes_path),
                }
            ),
            encoding="utf-8",
        )

    metadata, profile, symbols, observed = family_state_gene_profiles(
        tmp_path, GeneResolver({}, {})
    )
    rows = {
        (row.study_family, row.cell_state): index
        for index, row in metadata.iterrows()
    }
    columns = {gene: index for index, gene in enumerate(symbols)}
    assert np.isnan(profile[rows[("F1", "S1")], columns["C"]])
    assert np.isnan(profile[rows[("F2", "S1")], columns["B"]])
    assert not observed[rows[("F1", "S1")], columns["C"]]
    assert not observed[rows[("F2", "S1")], columns["B"]]
    assert np.isfinite(profile[rows[("F1", "S1")], columns["A"]])
