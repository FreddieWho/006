#!/usr/bin/env python3
"""Hostile final audit for v6.2.1 Phase7."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd
import yaml

from phase7_common import MODULES, ROOT, ensure_dirs, forbidden_columns, load_config, sha256, write_yaml


def run() -> dict:
    cfg = load_config()
    out = cfg["out"]
    ensure_dirs(out)
    checks = []

    def check(name: str, passed: bool, severity: str, evidence: str) -> None:
        checks.append({"check_id": name, "passed": bool(passed), "severity": severity, "evidence": evidence})

    closeout = yaml.safe_load((out / "preflight/phase6_closeout_manifest_v6_2_1.yaml").read_text())
    for asset in closeout["assets"]:
        path = ROOT / asset["path"]
        check(f"hash::{path.name}", path.exists() and sha256(path) == asset["sha256"], "hard", asset["path"])

    phase7a = yaml.safe_load((out / "handoff/phase7a_to_phase7b_handoff.yaml").read_text())
    phase7b = yaml.safe_load((out / "handoff/phase7b_to_phase8_handoff.yaml").read_text())
    eligibility = pd.read_csv(out / "handoff/module_eligibility_for_phase7b.csv")
    barriers = pd.read_csv(out / "identifiability/barrier_identifiability_table_v6_2_1.csv")
    sentinels = pd.read_csv(out / "sentinel/prespecified_sentinel_program_eligibility.csv")
    coverage = pd.read_csv(out / "measurement/sample_module_coverage_gate.csv")
    stability_a = pd.read_csv(out / "measurement/module_bootstrap_stability.csv")
    stability_b = pd.read_csv(out / "identifiability/bootstrap_rank_stability.csv")
    pairwise = pd.read_csv(out / "identifiability/barrier_pairwise_identifiability.csv")
    binding = pd.read_csv(out / "measurement/expression_unit_binding_table.csv")
    unit_registry = pd.read_csv(cfg["inputs"]["phase6_unit_registry"])
    matrix = pd.read_parquet(out / "handoff/eligible_barrier_score_matrix_v0.parquet")
    merge = pd.read_csv(out / "identifiability/coarse_barrier_merge_map_v6_2_1.csv")

    check("all_8_topics_phase7a", set(eligibility.module_id) == set(MODULES), "hard", str(sorted(eligibility.module_id.unique())))
    check("all_8_topics_phase7b", set(barriers.module_id) == set(MODULES), "hard", str(sorted(barriers.module_id.unique())))
    check("all_5_sentinels", len(sentinels) == 5 and sentinels.program_id.nunique() == 5, "hard", str(sentinels.program_id.tolist()))
    check("no_response_used_phase7a", phase7a.get("response_or_outcome_used") is False, "hard", str(phase7a.get("response_or_outcome_used")))
    check("no_response_direction_phase7b", phase7b.get("response_direction_evaluated") is False, "hard", str(phase7b.get("response_direction_evaluated")))
    check("module_discovery_closed", phase7a.get("module_discovery_reopened") is False and phase7b.get("module_discovery_reopened") is False, "hard", "both manifests")
    check("phase8_matrix_no_leakage_columns", not forbidden_columns(matrix.columns, cfg["prohibited_tokens"]), "hard", "|".join(forbidden_columns(matrix.columns, cfg["prohibited_tokens"])))
    identity = {"cohort_id", "patient_key", "timepoint", "tissue_context", "patient_timepoint_context_id"}
    barrier_cols = [c for c in matrix if c not in identity]
    expected = set(barriers.loc[barriers.allowed_for_phase8_conditional.astype(bool), "phase8_barrier_id"].dropna())
    check("phase8_matrix_only_eligible", set(barrier_cols) == expected, "hard", f"matrix={barrier_cols}; expected={sorted(expected)}")
    key = [c for c in ("cohort_id", "patient_key", "timepoint", "tissue_context", "patient_timepoint_context_id") if c in matrix]
    check("patient_timepoint_context_unique", not matrix.duplicated(key).any(), "hard", str(matrix.duplicated(key).sum()))
    check("no_unknown_patient_in_phase8", not matrix.patient_key.astype(str).str.lower().isin({"unknown", "nan", "none", ""}).any(), "hard", "patient_key")
    check("merge_map_unique_source", not merge.source_module_id.duplicated().any(), "hard", str(merge.source_module_id.duplicated().sum()))
    check("phase8_has_barrier", len(barrier_cols) >= 1, "hard", str(barrier_cols))
    check("full_atlas_not_claimed", phase7a.get("full_integration_claimed") is False, "soft", str(phase7a.get("full_integration_claimed")))
    valid_coverage = coverage[coverage.scoring_status.eq("primary")]
    check("coverage_gate_enforced", not valid_coverage.empty and (valid_coverage.available_fraction >= cfg["measurement"]["minimum_available_fraction"]).all() and (valid_coverage.n_states >= cfg["measurement"]["minimum_states_per_sample"]).all(), "hard", f"primary={len(valid_coverage)} low={int(coverage.scoring_status.eq('low_coverage').sum())}")
    expected_units = set(unit_registry.expression_unit_id.astype(str))
    bound_units = set(binding.expression_unit_id.astype(str))
    check(
        "all_expression_units_bound_or_audited",
        len(binding) == len(expected_units) and bound_units == expected_units,
        "hard",
        f"registry={len(expected_units)} rows={len(binding)} unique={len(bound_units)} missing={len(expected_units - bound_units)} extra={len(bound_units - expected_units)}",
    )
    check("recoverable_unknown_patient_fallback_applied", not binding.binding_status.eq("unresolved_patient_or_sample").any(), "hard", str(binding.binding_status.value_counts().to_dict()))
    expected_a = cfg["measurement"][f"bootstrap_n_{phase7a.get('mode', 'full')}"]
    expected_b = cfg["identifiability"][f"bootstrap_n_{phase7b.get('mode', 'full')}"]
    check("phase7a_bootstrap_complete", stability_a.bootstrap_id.nunique() == expected_a and stability_a.module_id.nunique() == 8, "hard", f"B={stability_a.bootstrap_id.nunique()} expected={expected_a}")
    geometry = stability_b.rank_stability.median()
    check("phase7b_geometry_rank_stability", stability_b.bootstrap_id.nunique() == expected_b and geometry >= cfg["identifiability"]["rank_stability_pass_min"], "hard", f"B={stability_b.bootstrap_id.nunique()} expected={expected_b} median={geometry}")
    check("partial_correlation_audited", {"partial_correlation", "abs_partial_correlation"}.issubset(pairwise.columns) and len(pairwise) == 28, "hard", str(pairwise.columns.tolist()))
    sentinel_valid = sentinels[sentinels.eligibility_status.eq("mainline_candidate")]
    check("sentinel_stability_and_coverage_gate", ((sentinel_valid.n_independent_cohorts >= 3) & (sentinel_valid.bootstrap_stability >= cfg["measurement_gate"]["residual_stability_mainline_min"])).all(), "hard", sentinel_valid.program_id.tolist())

    table = pd.DataFrame(checks)
    table.to_csv(out / "audit/phase7_strict_audit.csv", index=False)
    hard_fail = table[(table.severity == "hard") & (~table.passed)]
    verdict = "PASS" if hard_fail.empty else "FAIL"
    summary = {"phase": "v6.2.1_phase7_hostile_audit", "verdict": verdict, "n_checks": len(table),
               "n_hard_failures": len(hard_fail), "hard_failures": hard_fail.check_id.tolist()}
    write_yaml(summary, out / "audit/phase7_strict_audit_manifest.yaml")
    lines = ["# Phase7 Hostile Audit", "", f"**{verdict}**", "", f"- Checks: {len(table)}", f"- Hard failures: {len(hard_fail)}", "",
             "Phase7A 未使用疗效构造分数；Phase7B 未判断疗效方向；Phase8 仅允许读取冻结后的 eligible barrier matrix。"]
    if not hard_fail.empty:
        lines += ["", "## Hard Failures", ""] + [f"- {r.check_id}: {r.evidence}" for r in hard_fail.itertuples()]
    (out / "audit/PHASE7_HOSTILE_AUDIT_REPORT.md").write_text("\n".join(lines) + "\n")
    return summary


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.parse_args()
    result = run()
    print(yaml.safe_dump(result, sort_keys=False))
    raise SystemExit(0 if result["verdict"] == "PASS" else 1)
