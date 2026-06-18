from pathlib import Path
import sys

import numpy as np
import pandas as pd
import pytest
import yaml

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "results/v6_2/phase8_anchor_repair_and_statistical_strengthening_repair_v1"
RUN_MANIFEST = yaml.safe_load((OUT / "phase8_contract_and_run_manifest.yaml").read_text())
FULL_RUN = RUN_MANIFEST.get("run_status") == "COMPLETE"


def test_current_phase8_honors_phase7_primary_barrier_gate():
    route = yaml.safe_load((OUT / "phase8_route_gate_and_handoff.yaml").read_text())
    assert RUN_MANIFEST["run_status"] == "ABORTED_PRE_MODEL"
    assert route["verdict"] == "PHASE8_NOT_RUN_PHASE7_GATE_BLOCKED"
    assert route["eligible_primary_barriers"] == ["FM07"]
    assert set(route["blocked_primary_barriers"]) == {"FM01", "FM04"}
    assert RUN_MANIFEST["rules"]["phase7_gate_bypassed"] is False
    assert RUN_MANIFEST["rules"]["response_model_run"] is False
    assert not (OUT / "phase8_anchor_master.parquet").exists()


@pytest.mark.skipif(not FULL_RUN, reason="current Phase8 stopped at the Phase7 primary-barrier gate")
def test_source_repairs_and_expression_semantics_are_honest():
    bridge = pd.read_csv(OUT / "phase8_coverage_and_bridge.csv")
    g120 = bridge[bridge.record_type.eq("GSE120575_expression_semantics")].iloc[0]
    assert g120.status == "support_only_normalized_expression"
    g123 = bridge[bridge.record_type.eq("GSE123813_response_rescue")].iloc[0]
    assert g123.status == "source_RECIST_restored"
    inventory = bridge[bridge.record_type.eq("GSE120575_baseline_monotherapy_inventory")].iloc[0]
    assert inventory.n_patients == 12
    assert "primary=0" in inventory.detail and "raw-count semantics failed" in inventory.detail
    arm = bridge[bridge.record_type.eq("GSE243013_arm_repair")]
    assert arm.status.eq("clean_monotherapy_source_verified").sum() == 1


@pytest.mark.skipif(not FULL_RUN, reason="current Phase8 stopped at the Phase7 primary-barrier gate")
def test_strict_surface_is_patient_level_and_endpoint_preserving():
    x = pd.read_parquet(OUT / "phase8_anchor_master.parquet")
    strict = x[x.primary_role.eq("conditional_primary") & x.baseline_eligible & x.response_known & x[["FM01", "FM04", "FM07"]].notna().all(axis=1)]
    assert not strict.duplicated(["cohort_id", "patient_key", "normalized_timepoint", "tissue_context", "resolution"]).any()
    assert not strict.duplicated(["cohort_id", "patient_key"]).any()
    assert set(strict.endpoint_type).issuperset({"RECIST", "mRECIST", "pathologic_response"})
    assert not strict.cohort_id.eq("GSE120575").any()
    assert strict.cohort_id.str.startswith("GSE123813").any()


@pytest.mark.skipif(not FULL_RUN, reason="current Phase8 stopped at the Phase7 primary-barrier gate")
def test_gse301_projection_is_not_promoted_to_direct_mid():
    bridge = pd.read_csv(OUT / "phase8_coverage_and_bridge.csv")
    x = bridge[bridge.record_type.eq("GSE301741_frozen_direct_mid_reconstruction")]
    assert x.sample_id.nunique() == 27
    assert x[x.normalized_timepoint.eq("baseline")].patient_key.nunique() == 11
    assert x.status.str.contains("BRIDGE_FAILED").any()


@pytest.mark.skipif(not FULL_RUN, reason="current Phase8 stopped at the Phase7 primary-barrier gate")
def test_effects_are_endpoint_aware_and_cross_endpoint_is_not_claim_meta():
    x = pd.read_csv(OUT / "phase8_effects_master.csv")
    endpoint = x[x.record_type.eq("endpoint_meta")]
    assert endpoint.endpoint_type.ne("mixed_not_pooled_for_claim").all()
    mixed = x[x.record_type.eq("cross_endpoint_direction_only")]
    assert mixed.endpoint_type.eq("mixed_not_pooled_for_claim").all()
    within = x[x.record_type.eq("within_cohort_effect")]
    assert {"cliffs_delta", "firth_log_odds", "method_direction_agreement"}.issubset(within.columns)
    for _, row in endpoint.iterrows():
        matching = within[(within.barrier_id.eq(row.barrier_id)) & within.endpoint_type.eq(row.endpoint_type)
                          & within.resolution.eq(row.resolution) & within.statistical_estimability.eq("estimable")]
        matching = matching[np.where(matching.treatment_arm.str.contains("monotherapy", case=False), "monotherapy", "other") == row.arm_family]
        assert row.n_independent_cohorts == matching.cohort_id.nunique()


def test_route_and_manifest_do_not_overclaim():
    route = yaml.safe_load((OUT / "phase8_route_gate_and_handoff.yaml").read_text())
    manifest = yaml.safe_load((OUT / "phase8_contract_and_run_manifest.yaml").read_text())
    assert route["supervised_SRB_allowed"] is False
    assert route["counterfactual_repair_allowed"] is False
    assert route["X_class_ranking_allowed"] is False
    assert manifest["rules"]["module_membership_modified"] is False
    assert manifest["rules"]["response_used_for_score_construction"] is False
    if manifest["run_status"] == "ABORTED_PRE_MODEL":
        assert manifest["rules"]["response_model_run"] is False
        assert manifest["rules"]["phase7_gate_bypassed"] is False
    else:
        assert route["calibration_status"] == "NOT_APPLICABLE_NO_PREDICTIVE_MODEL"
    import hashlib
    for relative, expected in manifest["input_hashes"].items():
        path = ROOT / relative
        assert path.exists()
        assert hashlib.sha256(path.read_bytes()).hexdigest() == expected
    if manifest["run_status"] == "COMPLETE":
        attempt = yaml.safe_load((OUT.parent / f"{OUT.name}.last_attempt.yaml").read_text())
        assert attempt["status"] == "COMPLETE"


@pytest.mark.skipif(not FULL_RUN, reason="current Phase8 stopped at the Phase7 primary-barrier gate")
def test_sentinel_counts_do_not_exceed_strict_surface():
    x = pd.read_csv(OUT / "phase8_sensitivity_and_sentinel.csv")
    anchor = pd.read_parquet(OUT / "phase8_anchor_master.parquet")
    strict = anchor[anchor.primary_role.eq("conditional_primary") & anchor.baseline_eligible & anchor.response_known]
    sent = x[x.lane.eq("sentinel") & x.program.isin(["TLS_B", "CD8_PROGENITOR_EXHAUSTION", "CD8_TERMINAL_EXHAUSTION"])]
    assert ((sent.n_R + sent.n_NR) <= strict.patient_key.nunique()).all()


def test_r2_state_is_reachable_only_with_two_hcc_cohorts():
    module_dir = ROOT / "scripts/v6_2"
    sys.path.insert(0, str(module_dir))
    try:
        import run_phase8_anchor_repair_strengthening as repair
    finally:
        sys.path.pop(0)
    effects = pd.DataFrame({
        "barrier_id": ["FM01", "FM01"], "statistical_estimability": ["estimable"] * 2,
        "treatment_arm": ["monotherapy"] * 2, "endpoint_type": ["RECIST"] * 2,
        "resolution": ["coarse"] * 2, "cohort_id": ["HCC_A", "HCC_B"],
        "effect": [0.4, 0.2], "method_direction_agreement": [True, True],
        "cancer_type": ["hcc", "HCC"],
    })
    meta = pd.DataFrame(columns=["barrier_id", "endpoint_type", "resolution", "arm_family", "record_type"])
    crossfit = pd.DataFrame(columns=["barrier_id", "cohort_id", "effect"])
    controls = pd.DataFrame({
        "barrier_id": ["FM01"], "control_id": ["within_cohort_response_permutation"],
        "empirical_p": [0.05],
    })
    bridge = pd.DataFrame(columns=["dataset_id", "status"])
    verdict = repair.route_verdict(effects, meta, crossfit, controls, bridge)
    fm01 = next(row for row in verdict["barriers"] if row["barrier_id"] == "FM01")
    assert fm01["route"] == "R2"
    assert verdict["phase9_entry_allowed"] is True
