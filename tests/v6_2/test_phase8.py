import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "scripts/v6_2"))

from run_phase8_anchor_context import auc_rank, hedges_g, random_effects, route_flags, same_endpoint_same_direction_replicated


OUT = Path(__file__).resolve().parents[2] / "results/v6_2/phase8_anchor_context_adjudication"


def test_auc_rank():
    assert auc_rank(np.array([0, 0, 1, 1]), np.array([0.1, 0.2, 0.8, 0.9])) == 1.0


def test_hedges_g_direction():
    d = pd.DataFrame({"response_binary": [0, 0, 0, 1, 1, 1], "x": [0, 1, 2, 3, 4, 5]})
    assert hedges_g(d, "x")["effect"] > 0


def test_random_effects():
    d = pd.DataFrame({"effect": [0.2, 0.4, 0.3], "se": [0.1, 0.1, 0.1]})
    r = random_effects(d)
    assert r["k"] == 3
    assert 0.2 < r["pooled_effect"] < 0.4


def test_phase8_anchor_contract():
    if not OUT.exists():
        return
    env = pd.read_csv(OUT / "preflight/anchor_environment_v6_2_1.csv")
    assert len(env) == 1036
    assert env.patient_timepoint_context_id.is_unique
    primary = pd.read_csv(OUT / "association/phase8_primary_patient_table.csv")
    assert len(primary) == 39
    assert (primary.response_binary == 1).sum() == 17
    assert (primary.response_binary == 0).sum() == 22
    assert "GSE243013" not in set(primary.cohort_id)
    g243 = env[env.cohort_id.eq("GSE243013")]
    assert not g243.eligibility_reason.eq("eligible_conditional_anchor").any()
    g286 = env[(env.cohort_id == "GSE286827") & (env.phase8_analysis_role == "conditional_anchor_primary")]
    assert len(g286) == 13
    assert g286.Neoadj_type.eq("D").all()


def test_phase8_bounded_strengthening_contract():
    if not OUT.exists():
        return
    projection = pd.read_csv(OUT / "projection/gse301741_sample_level_frozen_module_projection.csv")
    bridge = pd.read_csv(OUT / "projection/gse301741_projection_bridge_audit.csv")
    assert len(projection) == 27 and projection.patient_key.nunique() == 16
    assert bridge.n_shared_contexts.eq(4).all()
    assert not bridge.bridge_status.eq("projection_support").any()
    hcc = pd.read_csv(OUT / "hcc_sensitivity/gse206325_coarse_response_sensitivity.csv")
    hcc = hcc[hcc.analysis_scope.eq("post_treatment_patient_mean")]
    assert hcc.n_R.eq(7).all() and hcc.n_NR.eq(17).all()


def test_phase8_gate_and_audit_contract():
    if not OUT.exists():
        return
    gate = pd.read_csv(OUT / "handoff/phase9_entry_gate_by_barrier.csv")
    assert set(gate.barrier_id) == {"FM01", "FM04", "FM07"}
    assert not gate.phase9_gate_pass.any()
    assert gate.fold_aware_residualization_complete.all()
    audit = pd.read_csv(OUT / "audit/phase8_strict_audit.csv")
    assert audit.passed.all()


def test_synthetic_route_remains_blocked_without_context_and_fold_evidence():
    gate = pd.DataFrame({"statistical_screen_pass": [True], "shared_context_separation_estimable": [False],
                         "hcc_residual_estimable": [False], "fold_aware_residualization_complete": [False]})
    assert not (gate.statistical_screen_pass & gate.shared_context_separation_estimable & gate.fold_aware_residualization_complete).any()
    assert not (gate.statistical_screen_pass & gate.hcc_residual_estimable & gate.fold_aware_residualization_complete).any()


def test_r1_r1lite_r2_are_distinct_routes():
    r1 = route_flags(True, True, False, True, True, True)
    assert r1["R1_pass"] and not r1["R1_lite_pass"] and not r1["R2_pass"]
    lite = route_flags(True, True, False, True, False, True)
    assert lite["R1_lite_pass"] and not lite["R1_pass"]
    r2 = route_flags(True, False, True, True, False, False)
    assert r2["R2_pass"] and not r2["R1_pass"] and not r2["R1_lite_pass"]
    incomplete = route_flags(False, False, False, True, False, False)
    assert not incomplete["route_adjudication_complete"]


def test_same_endpoint_replication_requires_same_direction():
    conflicting = pd.DataFrame({"cohort_id": ["A", "B", "C"], "response_endpoint_type": ["RECIST", "RECIST", "mRECIST"], "effect": [1.0, -1.0, 1.0]})
    assert not same_endpoint_same_direction_replicated(conflicting)
    concordant = pd.DataFrame({"cohort_id": ["A", "B", "C"], "response_endpoint_type": ["RECIST", "RECIST", "mRECIST"], "effect": [1.0, 0.5, -1.0]})
    assert same_endpoint_same_direction_replicated(concordant)


def test_low_power_cohort_cannot_trigger_r1_replication():
    low_power = pd.DataFrame({"cohort_id": ["A", "B"], "response_endpoint_type": ["RECIST", "RECIST"],
                              "effect": [1.0, 0.5], "n_R": [10, 2], "n_NR": [10, 2]})
    assert not same_endpoint_same_direction_replicated(low_power, minimum_cohorts=2, minimum_patients_per_class=3)
