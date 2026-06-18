from pathlib import Path
import pandas as pd
import yaml

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "results/v6_2/phase8_bounded_repair_and_readjudication"


def test_contract_and_verdict_are_honest():
    contract = yaml.safe_load((OUT / "phase8_repair_input_contract.yaml").read_text())
    verdict = yaml.safe_load((OUT / "phase8_route_verdict.yaml").read_text())
    assert contract["freeze_status"] == "PASS"
    assert contract["primary_barriers"] == ["FM01", "FM04", "FM07"]
    assert verdict["verdict"] == "UNRESOLVED_PHASE9_BLOCKED"
    assert verdict["phase9_entry_allowed"] is False
    assert not (OUT / "phase8_to_phase9_handoff.yaml").exists()


def test_anchor_environment_has_no_endpoint_pooling_or_duplicate_contexts():
    x = pd.read_csv(OUT / "anchor_environment_v6_2_1_repaired.csv")
    required = {"patient_key", "normalized_timepoint", "tissue_context", "resolution", "endpoint_type", "treatment_arm", "FM01_available", "FM04_available", "FM07_available"}
    assert required.issubset(x.columns)
    assert not x.duplicated(["patient_key", "normalized_timepoint", "tissue_context", "resolution"]).any()
    g = x[x.cohort_id.eq("GSE286827")]
    assert set(g.treatment_arm.dropna()).issubset({"monotherapy_D_only", "combination_D_plus_T"})


def test_shared_decomposition_does_not_overclaim():
    x = pd.read_csv(OUT / "phase8_shared_context_decomposition.csv")
    shared = x[x.component.eq("shared")]
    hcc = x[x.component.eq("HCC_residual")]
    assert set(shared.barrier_id) == {"FM01", "FM04", "FM07"}
    assert not shared.identifiable.any()
    assert not hcc.identifiable.any()


def test_coverage_gate_uses_real_phase7_status_values():
    x = pd.read_csv(OUT / "coverage_selection_audit.csv")
    assert 0 < x.included.mean() < 1
    models = pd.read_csv(OUT / "coverage_selection_models.csv")
    assert models.loc[models.audit_scope.eq("all"), "missingness_response_auc"].notna().all()


def test_no_response_features_and_sentinel_limits():
    audit = pd.read_csv(OUT / "phase8_repair_input_audit.csv")
    assert audit.loc[audit.check.eq("no_label_columns_in_barrier_matrix"), "status"].iloc[0] == "PASS"
    x = pd.read_csv(OUT / "phase8_data_bridge_and_sensitivity_results.csv")
    low = x[x.program.isin(["NEUTROPHIL_NET", "TUMOR_WNT_EXCLUSION"])]
    assert low.status.str.contains("not_negative_evidence").all()
    sentinel = x[x.lane.eq("sentinel") & x.program.isin(["TLS_B", "CD8_PROGENITOR_EXHAUSTION", "CD8_TERMINAL_EXHAUSTION"])]
    assert ((sentinel.n_R + sentinel.n_NR) <= 39).all()
    topics = x[x.lane.eq("sensitivity_topic")]
    assert topics.endpoint_type.notna().all()
