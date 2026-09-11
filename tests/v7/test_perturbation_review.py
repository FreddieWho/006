import pytest
from v7.perturbation_review import condition_contrast


def test_exact_unpaired_null_enumerates_all_sample_assignments():
    result=condition_contrast([4,5,6],[1,2,3])
    assert result["difference"] == 3
    assert result["n_label_assignments"] == 20
    assert result["exact_sample_label_p"] == pytest.approx(.1)


def test_missing_conditions_do_not_create_a_direction():
    assert condition_contrast([1],[2,3])["status"].startswith("NOT_TESTABLE")
