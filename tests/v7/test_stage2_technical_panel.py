from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import sparse

from src.v7.ontology.aggregate import ObjectSpec, aggregate_object
from src.v7.ontology.pipeline import score_pseudobulk_shards
from src.v7.ontology.scoring import GeneResolver
from src.v7.ontology.technical_panel import (
    score_object_cell_level_panel,
    summarize_cell_vs_pseudobulk,
)
from src.v7.ontology.vocabulary import Stage2Vocabulary


def test_cell_level_summary_reports_rank_agreement_without_equivalence_claim() -> None:
    frame = pd.DataFrame(
        {
            "cohort_id": ["C"] * 3,
            "object_id": ["C::object_1"] * 3,
            "cell_state_level": ["coarse"] * 3,
            "feature_family": ["legacy_fm"] * 3,
            "feature_id": ["FM01"] * 3,
            "n_cells_total": [10, 20, 30],
            "n_cells_scoreable": [9, 18, 27],
            "cell_mean_score": [1.0, 2.0, 3.0],
            "pseudobulk_score": [4.0, 8.0, 10.0],
        }
    )

    result = summarize_cell_vs_pseudobulk(frame)

    assert len(result) == 1
    row = result.iloc[0]
    assert row.n_expression_units == 3
    assert row.n_comparable_units == 3
    assert row.n_cells_total == 60
    assert row.n_cells_scoreable == 54
    assert np.isclose(row.cell_scoreable_fraction, 0.9)
    assert np.isclose(row.cell_mean_vs_pseudobulk_spearman, 1.0)
    assert row.comparison_status == "OBSERVED_RANK_CONCORDANCE"
    assert row.estimand_relation.startswith("EXPECTED_NON_EQUIVALENCE")


def test_cell_level_summary_marks_constant_comparison_not_testable() -> None:
    frame = pd.DataFrame(
        {
            "cohort_id": ["C"] * 3,
            "object_id": ["O"] * 3,
            "cell_state_level": ["mid"] * 3,
            "feature_family": ["mechanism_component"] * 3,
            "feature_id": ["M01"] * 3,
            "n_cells_total": [2, 2, 2],
            "n_cells_scoreable": [2, 2, 2],
            "cell_mean_score": [1.0, 1.0, 1.0],
            "pseudobulk_score": [1.0, 2.0, 3.0],
        }
    )

    result = summarize_cell_vs_pseudobulk(frame)

    assert np.isnan(result.iloc[0].cell_mean_vs_pseudobulk_spearman)
    assert result.iloc[0].comparison_status == "NOT_TESTABLE_CONSTANT_OR_FEW_UNITS"


def test_full_cell_panel_path_matches_real_pseudobulk_keys(tmp_path) -> None:
    import anndata as ad

    counts = sparse.csr_matrix(
        np.asarray(
            [
                [10, 0, 1],
                [8, 1, 0],
                [0, 9, 2],
                [1, 7, 1],
                [4, 2, 6],
                [3, 1, 8],
            ],
            dtype=float,
        )
    )
    cells = [f"c{index}" for index in range(6)]
    h5ad = tmp_path / "object.h5ad"
    ad.AnnData(
        X=counts,
        obs=pd.DataFrame(index=cells),
        var=pd.DataFrame(index=["A", "B", "C"]),
    ).write_h5ad(h5ad)
    sidecar = tmp_path / "identity.parquet"
    pd.DataFrame(
        {
            "source_cell_id": cells,
            "analysis_unit_key": ["S1"] * 3 + ["S2"] * 3,
            "patient_key": ["P1"] * 3 + ["P2"] * 3,
            "normalized_timepoint": ["pre"] * 6,
            "tissue_context": ["tumor"] * 6,
            "assignment_status": ["resolved"] * 6,
        }
    ).to_parquet(sidecar, index=False)
    annotation = tmp_path / "annotations"
    annotation.mkdir()
    pd.DataFrame(
        {
            "source_object_id": ["C::object_1"] * 6,
            "source_cell_id": cells,
            "analysis_unit_key": ["S1"] * 3 + ["S2"] * 3,
            "study_subject_key": ["P1"] * 3 + ["P2"] * 3,
            "patient_key": ["P1"] * 3 + ["P2"] * 3,
            "patient_id": ["P1"] * 3 + ["P2"] * 3,
            "normalized_timepoint": ["pre"] * 6,
            "tissue_context": ["tumor"] * 6,
            "assignment_status": ["resolved"] * 6,
            "harmonized_coarse_label": ["T_NK", "T_NK", "Myeloid"] * 2,
            "harmonized_mid_label": ["CD8_T", "CD4_T", "Macrophage"] * 2,
            "harmonized_fine_label": ["Cytotoxic", "Helper", "Macro"] * 2,
            "original_annotation": ["x"] * 6,
            "mapping_confidence": ["high"] * 6,
            "low_quality_flag": [0] * 6,
            "ambient_rna_risk_flag": [0] * 6,
        }
    ).to_parquet(annotation / "cohort_id=C.parquet", index=False)
    spec = ObjectSpec(
        cohort_id="C",
        object_id="C::object_1",
        h5ad=h5ad,
        count_layer="X",
        identity_sidecar=sidecar,
        annotation_path=annotation,
    )
    vocabulary = Stage2Vocabulary(
        pd.DataFrame(
            {
                "feature_family": ["legacy_fm"],
                "feature_id": ["FM01"],
                "parent_axis_id": [""],
                "scoreable": [True],
            }
        ),
        pd.DataFrame(
            {
                "feature_family": ["legacy_fm"] * 3,
                "feature_id": ["FM01"] * 3,
                "gene_symbol": ["A", "B", "C"],
                "relative_weight": [0.5, 0.3, 0.2],
                "weight": [0.5, 0.3, 0.2],
                "direction": [1, 1, 1],
                "is_top_gene": [True, True, False],
            }
        ),
        pd.DataFrame(),
        pd.DataFrame(),
    )
    resolver = GeneResolver({}, {})
    output_root = tmp_path / "results"
    aggregate_object(spec, output_root, chunk_size=2, study_family="F")
    pseudobulk_scores, _ = score_pseudobulk_shards(
        output_root, vocabulary, resolver, workers=1
    )

    summary, comparison = score_object_cell_level_panel(
        spec,
        vocabulary,
        resolver,
        pseudobulk_scores,
        chunk_size=2,
    )

    assert set(comparison.cell_state_level) == {"coarse", "mid"}
    assert comparison.pseudobulk_score.notna().all()
    assert comparison.expression_unit_id.isin(
        pseudobulk_scores.expression_unit_id
    ).all()
    assert summary.n_cells_total.sum() == 12
