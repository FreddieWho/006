"""Regression tests for patient-level comparisons and vector-preserving nulls."""
import numpy as np
import pandas as pd
import pytest
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import brier_score_loss

from v7.clinical_baseline import _evaluate_model, _null_metric, _environment_groups
from v7.repair_analysis import _permute_patient_vectors


def test_common_patient_folds_survive_model_and_seed_changes():
    frame = pd.DataFrame({
        "patient_row_id": [f"p{i}" for i in range(20)],
        "response_binary": [0, 1] * 10,
        "x": np.arange(20), "z": np.arange(20) ** 2,
    })
    splits = list(StratifiedKFold(4, shuffle=True, random_state=13).split(frame, frame.response_binary))
    outputs = []
    for seed, features in [(1, ["x"]), (2, ["x", "z"])]:
        _, prediction = _evaluate_model(frame, features, environment="test", model_name=str(seed),
                                         seed=seed, bootstrap_draws=0, cv_splits=splits)
        outputs.append(prediction)
    assert outputs[0].fold.tolist() == outputs[1].fold.tolist()
    assert outputs[0].fold.nunique() == 4


def test_one_accession_with_two_treatments_is_two_environments():
    frame = pd.DataFrame({"analysis_dataset": ["g", "g"],
                          "treatment_values": ["mono", "combo"],
                          "endpoint_values": ["RECIST", "RECIST"]})
    groups = list(_environment_groups(frame))
    assert len(groups) == 2
    assert all(len(rows) == 1 for _, rows in groups)


def test_prevalence_baseline_is_trained_without_test_labels():
    frame = pd.DataFrame({"response_binary": [0, 0, 0, 1, 1, 1]})
    splits = [(np.array([2, 3, 4, 5]), np.array([0, 1])),
              (np.array([0, 1, 4, 5]), np.array([2, 3])),
              (np.array([0, 1, 2, 3]), np.array([4, 5]))]
    result = _null_metric(frame, "test", splits)
    assert result["brier"] == pytest.approx(brier_score_loss(frame.response_binary, [.75, .75, .5, .5, .25, .25]))


def test_pair_permutation_preserves_patient_feature_covariance():
    frame = pd.DataFrame([
        {"patient_key": str(i), "feature_id": f, "score_standardized_pre": 0.,
         "score_standardized_post": float(i * multiplier)}
        for i in range(10) for f, multiplier in [("a", 1), ("b", 7)]
    ])
    result = _permute_patient_vectors(frame, np.random.default_rng(11))
    vectors = result.pivot(index="patient_key", columns="feature_id", values="delta")
    np.testing.assert_array_equal(vectors.b, vectors.a * 7)
    assert sorted(vectors.a) == list(range(10))


def test_pair_permutation_rejects_duplicate_patient_features():
    frame = pd.DataFrame({"patient_key": ["p", "p"], "feature_id": ["a", "a"]})
    with pytest.raises(ValueError, match="duplicate"):
        _permute_patient_vectors(frame, np.random.default_rng(1))


def test_batched_moran_matches_reference_null_with_missing_and_sections():
    from scipy import sparse
    from v7.spatial_stats.statistics import moran_permutation_test, moran_statistic, permute_within_sections
    sections=np.array(["a"]*5+["b"]*5)
    graph=sparse.block_diag([np.ones((5,5))-np.eye(5)]*2,format="csr")
    values=np.array([1.,3.,np.nan,2.,8.,101.,109.,104.,np.nan,102.])
    rng=np.random.default_rng(41)
    reference=[moran_statistic(permute_within_sections(values,sections,rng),graph,section_ids=sections)["estimate"] for _ in range(67)]
    actual=moran_permutation_test(values,graph,sections,permutations=67,seed=41)
    np.testing.assert_allclose(actual["permuted_statistics"],reference,atol=1e-12,rtol=1e-12)
