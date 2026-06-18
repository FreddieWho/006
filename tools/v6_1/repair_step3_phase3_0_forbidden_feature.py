#!/usr/bin/env python3
"""Repair Step3 Phase 3.0 forbidden feature blocker.

This repair demotes one unresolved leakage-risk feature from the main feature
set. It does not modify the hotfix matrix, sample universe, split, response
labels, treatment context, or TF hotfix assets.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd
import yaml


ROOT = Path(__file__).resolve().parents[2]
STEP2_REPAIR = ROOT / "results/v6_1/step2_repair"
STEP2_SPTF = (
    ROOT
    / "results/v6_1/step2/step2_v6_1_0505_0319/08_signature_pathway_tf/signature_pathway_tf_feature_v1_tf_rescue_v1"
)
OUT = ROOT / "results/v6_1/step3_qc_aware_strong_baseline/00_input_contract"
TARGET = "signature__hcc_responder_like"
BLOCKED_REASON = "forbidden_name_and_unresolved_response_derivation_risk"


def bool_value(value: object) -> bool:
    return str(value).strip().lower() in {"true", "1", "yes"}


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    dictionary_path = STEP2_REPAIR / "repair_B2_feature_dictionary_v6_1.hotfix.csv"
    matrix_path = STEP2_REPAIR / "repair_B2_unified_feature_matrix_by_sample.hotfix.parquet"
    registry_path = STEP2_SPTF / "signature_registry_v6_1.csv"
    program_path = ROOT / "mvp/outputs/program_signatures.json"

    dictionary = pd.read_csv(dictionary_path)
    matrix_columns = set(pd.read_parquet(matrix_path).columns)
    registry = pd.read_csv(registry_path)
    target_rows = dictionary[dictionary["feature_name"].eq(TARGET)]
    if len(target_rows) != 1:
        raise SystemExit(f"Expected exactly one dictionary row for {TARGET}, found {len(target_rows)}")
    row = target_rows.iloc[0]
    registry_target = registry[
        (registry.get("feature_id", pd.Series(dtype=str)).astype(str).eq(TARGET))
        | (registry.get("feature_name", pd.Series(dtype=str)).astype(str).eq("hcc_responder_like"))
    ]

    included_before = bool_value(row["included_in_main"])
    sensitivity_before = bool_value(row["sensitivity_only"])
    repaired = dictionary.copy()
    if "blocked_reason" not in repaired.columns:
        repaired["blocked_reason"] = ""
    if "downstream_allowed" not in repaired.columns:
        repaired["downstream_allowed"] = repaired["included_in_main"].map(bool_value)
    mask = repaired["feature_name"].eq(TARGET)
    repaired.loc[mask, "included_in_main"] = False
    repaired.loc[mask, "sensitivity_only"] = True
    repaired.loc[mask, "blocked_reason"] = BLOCKED_REASON
    repaired.loc[mask, "downstream_allowed"] = False

    audit = pd.DataFrame(
        [
            {
                "feature_name": TARGET,
                "forbidden_pattern": "responder",
                "present_in_matrix": TARGET in matrix_columns,
                "present_in_dictionary": True,
                "feature_type": row.get("feature_type", ""),
                "source_step": row.get("source_step", ""),
                "calculation_method": row.get("calculation_method", ""),
                "gene_set_source": f"{row.get('gene_set', '')}; registry_source={';'.join(sorted(registry_target.get('source_resource', pd.Series(dtype=str)).dropna().astype(str).unique()))}; program_json={program_path}",
                "provenance_status": "insufficient_to_prove_frozen_external_non_response_derived",
                "response_derived_risk": "unresolved_high_risk_due_name_and_mvp_response_workflow_provenance",
                "repair_action": "demote_from_main_keep_matrix_column",
                "included_in_main_before": included_before,
                "included_in_main_after": False,
                "sensitivity_only_before": sensitivity_before,
                "sensitivity_only_after": True,
                "blocked_reason": BLOCKED_REASON,
                "notes": "Step2.9 report says response usage forbidden_not_used, but source traces to MVP program_signatures and responder-like naming; no evidence proves signature was frozen before response labels. Default to demotion.",
            }
        ]
    )

    blocked = audit[
        [
            "feature_name",
            "forbidden_pattern",
            "present_in_matrix",
            "present_in_dictionary",
            "feature_type",
            "source_step",
            "blocked_reason",
            "repair_action",
            "response_derived_risk",
            "notes",
        ]
    ].copy()
    blocked["downstream_allowed"] = False

    audit.to_csv(OUT / "repair_phase3_0_forbidden_feature_audit.csv", index=False)
    repaired.to_csv(OUT / "repair_phase3_0_feature_dictionary.repaired.csv", index=False)
    blocked.to_csv(OUT / "repair_phase3_0_blocked_features.csv", index=False)

    summary = "\n".join(
        [
            "# Repair Phase3.0 Forbidden Feature Summary",
            "",
            "## Why Phase 3.0 Failed",
            f"- `{TARGET}` matched forbidden pattern `responder` while `included_in_main=True`.",
            "",
            "## Provenance Decision",
            "- Feature dictionary: signature, Step2.9 signature/pathway/TF factory, mean signed gene z-score, gene set `hcc_responder_like`.",
            "- Upstream registry source: `mvp/outputs/program_signatures.json`.",
            "- Step2.9 report records `response usage: forbidden_not_used`, but available provenance does not prove this gene set was frozen before response labels or independent of responder/non-responder contrasts.",
            "- Therefore provenance is insufficient for whitelist.",
            "",
            "## Repair Action",
            "- Action: demote from main, keep matrix column unchanged.",
            "- `included_in_main`: True -> False.",
            "- `sensitivity_only`: False -> True.",
            f"- `blocked_reason`: `{BLOCKED_REASON}`.",
            "",
            "## Non-Changes",
            "- sample universe changed: false.",
            "- patient split changed: false.",
            "- response label changed: false.",
            "- treatment_context changed: false.",
            "- TF hotfix changed: false.",
            "- feature matrix changed: false.",
            "- Step3 downstream workflow changed: false; rerun starts at Phase 3.0 input contract only.",
        ]
    )
    (OUT / "repair_phase3_0_summary.md").write_text(summary + "\n", encoding="utf-8")

    status = {
        "repair_target": TARGET,
        "repair_type": "demote_or_remove_forbidden_feature",
        "sample_universe_changed": False,
        "patient_split_changed": False,
        "response_label_changed": False,
        "tf_hotfix_changed": False,
        "feature_matrix_changed": False,
        "feature_dictionary_changed": True,
        "rerun_from_phase": "Phase3_0",
        "allowed_to_continue_step3": False,
        "next_action": "rerun Phase3_0 input contract validation",
    }
    (OUT / "step3_phase3_0_input_contract_status.repair_ready.yaml").write_text(
        yaml.safe_dump(status, sort_keys=False), encoding="utf-8"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
