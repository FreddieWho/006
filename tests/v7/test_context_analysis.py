import numpy as np

from v7.context_analysis import _sign_label


def test_context_sign_label_keeps_mixed_direction_as_context_modulated():
    label, positive, negative = _sign_label([0.2, 0.1, -0.02, 0.05])
    assert label == "CONTEXT_MODULATED"
    assert np.isclose(positive, 0.75)
    assert np.isclose(negative, 0.25)


def test_context_sign_label_is_unresolved_without_finite_estimates():
    label, positive, negative = _sign_label([np.nan, np.inf])
    assert label == "UNRESOLVED_NO_ESTIMATES"
    assert np.isnan(positive)
    assert np.isnan(negative)
