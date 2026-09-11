from __future__ import annotations

import numpy as np
import pandas as pd
from pandas.testing import assert_frame_equal

from src.v7.ontology.diagnostics import (
    _profile_score,
    ambient_diagnostics,
    confounding_diagnostics,
    leave_family_out_state_profiles,
    matched_null_coherence,
    patient_timepoint_scores,
    robust_standardize_scores,
    summarize_reliability,
    technical_panel_reliability,
)


def _reference_patient_timepoint_scores(expression_scores: pd.DataFrame) -> pd.DataFrame:
    keys = [
        "cohort_id",
        "study_family",
        "patient_key",
        "timepoint",
        "cell_state_level",
        "feature_family",
        "feature_id",
    ]
    rows: list[dict[str, object]] = []
    source = expression_scores.loc[
        expression_scores.patient_key.fillna("").astype(str).str.strip().ne("")
    ]
    for key, frame in source.groupby(keys, dropna=False, sort=True):
        values = pd.to_numeric(frame.score_native, errors="coerce").to_numpy(float)
        weights = pd.to_numeric(frame.n_cells_used, errors="coerce").fillna(0).to_numpy(float)
        valid = np.isfinite(values)
        for method in ("abundance_weighted", "equal_state"):
            if not valid.any():
                score = float("nan")
            elif method == "abundance_weighted" and weights[valid].sum() > 0:
                score = float(np.average(values[valid], weights=weights[valid]))
            else:
                score = float(np.mean(values[valid]))
            row = dict(zip(keys, key))
            row.update(
                {
                    "patient_timepoint_id": f"{row['cohort_id']}::{row['patient_key']}::{row['timepoint']}",
                    "aggregation_level": "patient_timepoint",
                    "aggregation_method": method,
                    "modality": "scRNA",
                    "cell_state": "ALL_STATES",
                    "score_native": score,
                    "n_states": int(frame.loc[valid, "cell_state"].nunique()),
                    "n_cells_used": int(weights[valid].sum()),
                    "coverage": float(pd.to_numeric(frame.coverage, errors="coerce").median()),
                    "technical_status": "valid" if np.isfinite(score) else "reject",
                }
            )
            rows.append(row)
    return pd.DataFrame(rows)


def test_vectorized_patient_aggregation_matches_reference_and_keeps_feature_level() -> None:
    frame = pd.DataFrame(
        {
            "cohort_id": ["C"] * 8,
            "study_family": ["F"] * 8,
            "patient_key": ["P"] * 8,
            "timepoint": ["pre"] * 8,
            "cell_state_level": ["coarse"] * 4 + ["mid"] * 4,
            "cell_state": ["A", "B", "A", "B", "A1", "A2", "A1", "A2"],
            "feature_family": ["legacy_fm"] * 8,
            "feature_id": ["FM01", "FM01", "FM02", "FM02"] * 2,
            "score_native": [1.0, 3.0, np.nan, np.nan, 2.0, 4.0, 8.0, 10.0],
            "n_cells_used": [1, 3, 2, 2, 0, 0, 1, 1],
            "coverage": [0.9, 0.8, 0.7, 0.6, 0.9, 0.9, 0.5, 0.7],
        }
    )

    observed = patient_timepoint_scores(frame)
    expected = _reference_patient_timepoint_scores(frame)
    sort = ["cell_state_level", "feature_id", "aggregation_method"]
    observed = observed.sort_values(sort).reset_index(drop=True)
    expected = expected.sort_values(sort).reset_index(drop=True)

    assert set(observed.feature_id) == {"FM01", "FM02"}
    assert set(observed.cell_state_level) == {"coarse", "mid"}
    common = sorted(set(expected.columns).intersection(observed.columns))
    assert_frame_equal(
        observed[common], expected[common], check_dtype=False, check_like=True
    )


def test_vectorized_robust_standardization_matches_scalar_reference() -> None:
    frame = pd.DataFrame(
        {
            "cohort_id": ["A"] * 5 + ["B"] * 4 + ["C"] * 3,
            "feature_id": ["F1"] * 5 + ["F1"] * 4 + ["F2"] * 3,
            "score_native": [1.0, 2.0, 4.0, 10.0, np.nan, 3.0, 3.0, 3.0, 3.0, 1.0, 1.0, 2.0],
        }
    )

    def reference(series: pd.Series) -> pd.Series:
        values = pd.to_numeric(series, errors="coerce")
        median = values.median()
        scale = 1.4826 * (values - median).abs().median()
        if not np.isfinite(scale) or scale == 0:
            scale = (values.quantile(0.75) - values.quantile(0.25)) / 1.349
        if not np.isfinite(scale) or scale == 0:
            scale = values.std(ddof=0)
        if not np.isfinite(scale) or scale == 0:
            scale = 1.0
        return (values - median) / scale

    expected = frame.groupby(
        ["cohort_id", "feature_id"], dropna=False, sort=False
    ).score_native.transform(reference)
    observed = robust_standardize_scores(
        frame, group_columns=("cohort_id", "feature_id")
    ).score_standardized

    assert np.allclose(observed, expected, equal_nan=True, rtol=0, atol=1e-12)


def test_technical_panel_reliability_preserves_mixed_route_evidence() -> None:
    panel = pd.DataFrame(
        {
            "feature_id": ["stable"] * 6 + ["mixed"] * 6,
            "cell_mean_vs_pseudobulk_spearman": [
                0.7,
                0.8,
                0.9,
                0.75,
                0.85,
                0.95,
                -0.4,
                -0.2,
                0.1,
                0.2,
                0.3,
                0.4,
            ],
        }
    )

    result = technical_panel_reliability(
        panel, seed=17, bootstrap_replicates=500
    ).set_index("feature_id")

    assert result.loc["stable", "technical_panel_status"] == "PASS_DIRECTIONAL_STABILITY"
    assert result.loc["stable", "n_negative_technical_panel_comparisons"] == 0
    assert result.loc["mixed", "n_negative_technical_panel_comparisons"] == 2
    assert result.loc["mixed", "technical_panel_status"] == "MIXED_OR_NOT_TESTABLE"


def test_matched_null_uses_boolean_legacy_top_genes() -> None:
    matrix = np.arange(1, 61, dtype=float).reshape(5, 12)
    membership = pd.DataFrame(
        {
            "feature_id": ["FM01"] * 6,
            "feature_family": ["legacy_fm"] * 6,
            "gene_symbol": list("ABCDEF"),
            "direction": [1] * 6,
            "weight": [1.0] * 6,
            "is_top_gene": [True, True, True, True, False, False],
        }
    )
    dictionary = pd.DataFrame(
        {"feature_id": ["FM01"], "scoreable": [True]}
    )

    summary, records = matched_null_coherence(
        matrix,
        list("ABCDEFGHIJKL"),
        membership,
        dictionary,
        replicates=2,
        seed=3,
    )

    assert summary.loc[0, "n_profile_genes"] == 4
    assert len(records) == 2


def test_d2_measurable_requires_positive_cell_vs_pseudobulk_panel_interval() -> None:
    expression = pd.DataFrame(
        {
            "feature_id": ["stable", "mixed"],
            "feature_family": ["mechanism_component"] * 2,
            "expression_unit_id": ["u1", "u2"],
            "study_family": ["S1", "S2"],
            "coverage": [1.0, 1.0],
            "technical_status": ["valid", "valid"],
        }
    )
    leave_out = pd.DataFrame(
        {
            "feature_id": ["stable"] * 3 + ["mixed"] * 3,
            "state_profile_spearman": [0.5, 0.6, 0.7] * 2,
        }
    )
    method = pd.DataFrame(
        {
            "feature_id": ["stable", "mixed"],
            "method_agreement_ci_low": [0.3, 0.3],
        }
    )
    nulls = pd.DataFrame(
        {
            "feature_id": ["stable", "mixed"],
            "matched_null_ci_high": [0.2, 0.2],
            "observed_coherence_ci_low": [0.8, 0.8],
        }
    )
    panel = pd.DataFrame(
        {
            "feature_id": ["stable", "mixed"],
            "technical_panel_ci_low": [0.25, -0.1],
        }
    )

    result = summarize_reliability(
        expression,
        leave_out,
        pd.DataFrame(),
        pd.DataFrame(),
        pd.DataFrame(),
        method,
        null_intervals=nulls,
        technical_panel=panel,
        seed=4,
        bootstrap_replicates=200,
    ).set_index("feature_id")

    assert result.loc["stable", "D2_class"] == "measurable"
    assert result.loc["mixed", "D2_class"] == "joint_representation_with_uncertainty"


def test_leave_family_out_compares_within_family_state_profiles() -> None:
    rows = []
    for family, offset in (("S1", 100.0), ("S2", 10.0), ("S3", -20.0)):
        for state, value in (("A", 1.0), ("B", 2.0), ("C", 4.0), ("D", 8.0)):
            rows.append(
                {
                    "feature_id": "F1",
                    "study_family": family,
                    "cell_state": state,
                    "score_native": offset + value,
                    "score_standardized": 0.0,
                    "technical_status": "valid",
                }
            )
    result = leave_family_out_state_profiles(pd.DataFrame(rows))

    assert len(result) == 3
    assert result.diagnostic_status.eq("PASS").all()
    assert np.allclose(result.state_profile_spearman, 1.0)


def test_profile_score_renormalizes_finite_genes_instead_of_filling_nan_with_zero() -> None:
    matrix = np.asarray([[2.0, np.nan], [2.0, 0.0], [np.nan, np.nan]])

    score = _profile_score(matrix, [0, 1], [1, 1], [1.0, 1.0])

    assert score[0] == 2.0
    assert score[1] == 1.0
    assert np.isnan(score[2])


def test_matched_null_requires_identical_observation_masks() -> None:
    matrix = np.asarray(
        [
            [1.0, 2.0, 3.0, 4.0, 1.2, 2.2, 3.2, 4.2, 9.0],
            [2.0, np.nan, 4.0, np.nan, 2.2, np.nan, 4.2, np.nan, 8.0],
            [3.0, 4.0, 5.0, 6.0, 3.2, 4.2, 5.2, 6.2, 7.0],
            [4.0, np.nan, 6.0, np.nan, 4.2, np.nan, 6.2, np.nan, 6.0],
            [5.0, 6.0, 7.0, 8.0, 5.2, 6.2, 7.2, 8.2, 5.0],
        ]
    )
    symbols = list("ABCDEFGHI")
    membership = pd.DataFrame(
        {
            "feature_id": ["F"] * 4,
            "feature_family": ["mechanism_component"] * 4,
            "gene_symbol": list("ABCD"),
            "direction": [1, 1, 1, 1],
            "weight": [1.0] * 4,
            "is_top_gene": [True] * 4,
        }
    )
    dictionary = pd.DataFrame({"feature_id": ["F"], "scoreable": [True]})

    summary, records = matched_null_coherence(
        matrix, symbols, membership, dictionary, replicates=4, seed=11
    )

    assert summary.loc[0, "diagnostic_status"] == "PASS_COMPUTED"
    assert records.observation_mask_match.all()
    for row in records.itertuples(index=False):
        targets = row.target_gene_symbols.split("|")
        nulls = row.null_gene_symbols.split("|")
        assert len(targets) == len(nulls) == 4
        assert len(set(nulls)) == 4
        assert set(targets).isdisjoint(nulls)
        for target, null in zip(targets, nulls):
            assert np.array_equal(
                np.isfinite(matrix[:, symbols.index(target)]),
                np.isfinite(matrix[:, symbols.index(null)]),
            )


def test_matched_null_is_explicitly_not_testable_when_exact_masks_are_insufficient() -> None:
    matrix = np.asarray(
        [
            [1.0, 2.0, 3.0, 4.0, 5.0],
            [2.0, np.nan, 4.0, np.nan, 6.0],
            [3.0, 4.0, 5.0, 6.0, 7.0],
            [4.0, np.nan, 6.0, np.nan, 8.0],
        ]
    )
    membership = pd.DataFrame(
        {
            "feature_id": ["F"] * 4,
            "feature_family": ["mechanism_component"] * 4,
            "gene_symbol": list("ABCD"),
            "direction": [1] * 4,
            "weight": [1.0] * 4,
            "is_top_gene": [True] * 4,
        }
    )
    dictionary = pd.DataFrame({"feature_id": ["F"], "scoreable": [True]})

    summary, records = matched_null_coherence(
        matrix, list("ABCDE"), membership, dictionary, replicates=2, seed=2
    )

    assert summary.loc[0, "diagnostic_status"] == "NOT_TESTABLE_INSUFFICIENT_EXACT_MASK_NULLS"
    assert records.empty


def test_ambient_diagnostic_uses_native_scores_within_matched_context() -> None:
    scores = pd.DataFrame(
        {
            "feature_id": ["F"] * 6,
            "cohort_id": ["C1"] * 3 + ["C2"] * 3,
            "patient_key": ["P1"] * 3 + ["P2"] * 3,
            "timepoint": ["pre"] * 6,
            "cell_state_level": ["coarse"] * 6,
            "cell_state": ["Low_quality_or_ambient", "T_NK", "Myeloid"] * 2,
            "score_native": [10.0, 2.0, 4.0, 20.0, 8.0, 12.0],
            "score_standardized": [0.0] * 6,
        }
    )

    result = ambient_diagnostics(scores).iloc[0]

    assert result["ambient_minus_other_median_native"] == 8.5
    assert result["n_matched_contexts"] == 2
    assert result["ambient_status"] == "OBSERVED_LOW_QUALITY_STATE_PROXY_AVAILABLE"


def test_confounding_keeps_cohort_shift_after_global_state_residualization() -> None:
    rows = []
    for cohort, shift in (("C1", 0.0), ("C2", 10.0)):
        for patient in ("P1", "P2"):
            for state, baseline in (("T_NK", 1.0), ("Myeloid", 4.0)):
                rows.append(
                    {
                        "feature_id": "F",
                        "cohort_id": cohort,
                        "patient_key": f"{cohort}_{patient}",
                        "timepoint": "pre",
                        "cell_state_level": "coarse",
                        "cell_state": state,
                        "platform": "p1" if patient == "P1" else "p2",
                        "tissue_context": "tumor",
                        "score_native": baseline + shift,
                        "score_standardized": 0.0,
                        "library_size": 100.0,
                        "n_cells_used": 10,
                    }
                )

    result = confounding_diagnostics(pd.DataFrame(rows)).iloc[0]

    assert result["cohort_eta_squared"] == 1.0
    assert result["cohort_status"] == "COMPUTED_PATIENT_TIMEPOINT_EQUAL_WEIGHT"
    assert result["n_patient_timepoints"] == 4


def test_confounding_marks_single_platform_not_testable() -> None:
    frame = pd.DataFrame(
        {
            "feature_id": ["F"] * 4,
            "cohort_id": ["C1", "C1", "C2", "C2"],
            "patient_key": ["P1", "P1", "P2", "P2"],
            "timepoint": ["pre"] * 4,
            "cell_state_level": ["coarse"] * 4,
            "cell_state": ["T_NK", "Myeloid"] * 2,
            "platform": ["single"] * 4,
            "tissue_context": ["tumor"] * 4,
            "score_native": [1.0, 2.0, 3.0, 4.0],
            "library_size": [1.0, 2.0, 3.0, 4.0],
            "n_cells_used": [10] * 4,
        }
    )

    result = confounding_diagnostics(frame).iloc[0]

    assert np.isnan(result["platform_eta_squared"])
    assert result["platform_status"] == "NOT_TESTABLE_SINGLE_PLATFORM"
