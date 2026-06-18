import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "scripts/v6_2"))

from phase7_common import ROOT, forbidden_columns, hypergeom_sf, normalize_timepoint, precision_weight, robust_z, weighted_residual


def test_normalize_timepoint():
    assert normalize_timepoint("Baseline") == "baseline"
    assert normalize_timepoint("week 6") == "week_6"
    assert normalize_timepoint("post-treatment") == "post_treatment"


def test_precision_weight_is_capped():
    got = precision_weight(pd.Series([0, 25, 100, 4_000_000]), 100)
    assert np.allclose(got, [0, 0.5, 1, 1])


def test_robust_z_falls_back_to_global():
    values = pd.Series([1.0, 2.0, 3.0, 4.0])
    groups = pd.Series(["a", "a", "b", "b"])
    got = robust_z(values, groups, min_n=3)
    assert got.notna().all()
    assert abs(got.median()) < 1e-12


def test_weighted_residual_removes_linear_signal():
    x = np.column_stack([np.ones(20), np.arange(20)])
    y = 2 + 3 * np.arange(20, dtype=float)
    residual, r2 = weighted_residual(y, x, np.ones(20))
    assert np.nanmax(np.abs(residual)) < 1e-5
    assert r2 > 0.999


def test_leakage_token_boundaries():
    bad = forbidden_columns(["FM01", "response_binary", "sample_key", "train_split"], ["response", "split", "train"])
    assert bad == ["response_binary", "train_split"]


def test_hypergeometric_tail_without_scipy():
    assert np.isclose(hypergeom_sf(0, 10, 2, 2), 17 / 45)


def test_full_phase7_gate_artifacts_are_self_consistent():
    out = ROOT / "results/v6_2/phase7_module_measurement_and_barrier_identifiability"
    if not out.exists():
        return
    coverage = pd.read_csv(out / "measurement/sample_module_coverage_gate.csv")
    valid = coverage[coverage.scoring_status.eq("primary")]
    assert not valid.empty
    assert (valid.available_fraction >= 0.5).all()
    assert (valid.n_states >= 2).all()
    binding = pd.read_csv(out / "measurement/expression_unit_binding_table.csv")
    assert len(binding) == binding.expression_unit_id.nunique() == 10947
    assert not binding.binding_status.eq("unresolved_patient_or_sample").any()


def test_full_phase7_identifiability_contract():
    out = ROOT / "results/v6_2/phase7_module_measurement_and_barrier_identifiability"
    if not out.exists():
        return
    pairs = pd.read_csv(out / "identifiability/barrier_pairwise_identifiability.csv")
    assert len(pairs) == 28
    assert pairs.partial_correlation.notna().all()
    bootstrap = pd.read_csv(out / "identifiability/bootstrap_rank_stability.csv")
    assert bootstrap.bootstrap_id.nunique() in {100, 1000}
    assert bootstrap.rank_stability.median() >= 0.70
    barriers = pd.read_csv(out / "identifiability/barrier_identifiability_table_v6_2_1.csv")
    matrix = pd.read_parquet(out / "handoff/eligible_barrier_score_matrix_v0.parquet")
    expected = set(barriers.loc[barriers.allowed_for_phase8_conditional, "phase8_barrier_id"])
    identity = {"cohort_id", "patient_key", "timepoint", "tissue_context", "patient_timepoint_context_id"}
    assert set(matrix.columns) - identity == expected
