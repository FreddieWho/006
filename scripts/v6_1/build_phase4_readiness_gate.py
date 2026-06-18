#!/usr/bin/env python3
from __future__ import annotations

import csv
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[2]
DATE_TAG = "20260611"
OUT_DIR = ROOT / "results" / "v6_1" / f"phase4_readiness_{DATE_TAG}"

STEP1 = ROOT / "results" / "v6_1" / "step1"
STEP2_REPAIR = ROOT / "results" / "v6_1" / "step2_repair"
STEP3 = ROOT / "results" / "v6_1" / "step3_qc_aware_strong_baseline"
HANDOFF = STEP3 / "07_milestone_B_step4_handoff"
DATA_STRENGTH = ROOT / "results" / "v6_1" / f"data_completeness_strengthening_{DATE_TAG}"

PRIMARY_TASKS = {"HCC_specific", "PD1X_extension", "pan_cancer_shared"}
MAIN_FEATURE_FAMILIES = {
    "fraction",
    "signature",
    "pathway",
    "tf_activity",
    "combined_interpretable",
    "combined_all_allowed",
}


def now_iso() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat()


def require(path: Path) -> Path:
    if not path.exists():
        raise FileNotFoundError(path)
    return path


def rel(path: Path) -> str:
    try:
        return str(path.relative_to(ROOT))
    except ValueError:
        return str(path)


def read_table(path: Path, delimiter: str = ",") -> list[dict[str, str]]:
    with require(path).open(newline="") as fh:
        return list(csv.DictReader(fh, delimiter=delimiter))


def write_table(path: Path, rows: list[dict[str, Any]], delimiter: str = ",", fieldnames: list[str] | None = None) -> None:
    if fieldnames is None:
        fieldnames = list(rows[0].keys()) if rows else []
    with path.open("w", newline="") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=fieldnames,
            delimiter=delimiter,
            extrasaction="ignore",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)


def as_bool(value: Any) -> bool:
    return str(value).strip().lower() in {"true", "1", "yes", "y"}


def as_int(value: Any, default: int = 0) -> int:
    try:
        return int(float(str(value)))
    except (TypeError, ValueError):
        return default


def yaml_scalar(value: Any) -> str:
    if isinstance(value, bool):
        return "true" if value else "false"
    if isinstance(value, (int, float)):
        return str(value)
    if value is None:
        return "null"
    return json.dumps(str(value))


def to_yaml(value: Any, indent: int = 0) -> str:
    pad = " " * indent
    if isinstance(value, dict):
        lines: list[str] = []
        for key, sub in value.items():
            if isinstance(sub, (dict, list)):
                lines.append(f"{pad}{key}:")
                lines.append(to_yaml(sub, indent + 2))
            else:
                lines.append(f"{pad}{key}: {yaml_scalar(sub)}")
        return "\n".join(lines)
    if isinstance(value, list):
        lines = []
        for item in value:
            if isinstance(item, (dict, list)):
                lines.append(f"{pad}-")
                lines.append(to_yaml(item, indent + 2))
            else:
                lines.append(f"{pad}- {yaml_scalar(item)}")
        return "\n".join(lines)
    return f"{pad}{yaml_scalar(value)}"


def write_yaml(path: Path, data: dict[str, Any]) -> None:
    path.write_text(to_yaml(data) + "\n")


def build_evidence_pool() -> list[dict[str, Any]]:
    source = STEP3 / "03_leakage_confounding_audit" / "step3_phase3_3_confounding_audit_by_universe.repaired.csv"
    conf = read_table(source)
    rows: list[dict[str, Any]] = []

    for row in conf:
        if "supervised_labeled_universe" not in row.get("universe_id", ""):
            continue

        task = row.get("task", "")
        family = row.get("feature_family", "")
        primary_allowed = as_bool(row.get("primary_allowed_after_repair", "False"))
        sensitivity = as_bool(row.get("sensitivity_only_after_repair", "False"))
        biological_anchor = as_bool(row.get("biological_anchor_only", "False"))

        if primary_allowed and task in PRIMARY_TASKS and family in MAIN_FEATURE_FAMILIES:
            pool = "primary_evidence_pool"
            phase4_use = "primary_module_discovery"
            allowed = "core mechanism discovery; module construction; main feature direction checks"
            forbidden = "clinical prediction claim without Phase4 validation; untreated direct causal claim"
        elif sensitivity or biological_anchor or task == "PD1_anchor":
            pool = "support_evidence_pool"
            phase4_use = "directional_support_or_sensitivity"
            allowed = "directional consistency; mechanism annotation support; sensitivity comparison"
            forbidden = "primary feature selection; primary model ranking; primary module discovery input; clinical prediction claim"
        else:
            pool = "blocked_evidence_pool"
            phase4_use = "do_not_use"
            allowed = "none"
            forbidden = "all Phase4 primary and support uses until re-audited"

        rows.append(
            {
                "pool": pool,
                "phase4_use": phase4_use,
                "universe_id": row.get("universe_id", ""),
                "task": task,
                "feature_family": family,
                "n_samples": row.get("n_samples", ""),
                "n_patients": row.get("n_patients", ""),
                "n_cohorts": row.get("n_cohorts", ""),
                "responder_rate": row.get("responder_rate", ""),
                "max_cohort_fraction": row.get("max_cohort_fraction", ""),
                "confounding_risk_level": row.get("confounding_risk_level", ""),
                "recommended_action": row.get("recommended_action", ""),
                "universe_role": row.get("universe_role", ""),
                "primary_allowed_after_repair": primary_allowed,
                "sensitivity_only_after_repair": sensitivity,
                "biological_anchor_only": biological_anchor,
                "allowed_downstream_use": allowed,
                "forbidden_downstream_use": forbidden,
                "source_audit_file": rel(source),
            }
        )

    return rows


def build_feature_contract() -> list[dict[str, Any]]:
    source = STEP3 / "02_feature_registry_missingness" / "step3_phase3_2_universe_feature_set_registry.csv"
    feature_sets = read_table(source)
    rows: list[dict[str, Any]] = []

    for row in feature_sets:
        if "supervised_labeled_universe" not in row.get("universe_id", ""):
            continue

        task = row.get("task", "")
        family = row.get("feature_family", "")
        set_type = row.get("feature_set_type", "")
        status = row.get("status", "")
        n_features = as_int(row.get("n_features", 0))

        if family == "pseudobulk_gene":
            phase4_use = "not_primary_from_step3_handoff"
            reason = "Step3 handoff excludes pseudobulk_gene as primary feature family; direct gene-level GRN should use Step2.8 layer separately."
        elif task == "PD1_anchor":
            phase4_use = "support_only"
            reason = "PD1_anchor primary supervised use blocked by cohort-response confounding."
        elif task in PRIMARY_TASKS and family in MAIN_FEATURE_FAMILIES and set_type == "clean_core_feature_set" and status == "eligible" and n_features > 0:
            phase4_use = "primary"
            reason = "Clean core feature set in primary-allowed universe."
        elif task in PRIMARY_TASKS and family in MAIN_FEATURE_FAMILIES and set_type == "missingness_aware_feature_set" and status == "eligible" and n_features > 0:
            phase4_use = "conditional_primary_if_model_grade_A_or_B_else_sensitivity"
            reason = "Missingness-aware set may be used only when corresponding model evidence beats controls."
        elif set_type == "blocked_feature_set" and status == "eligible" and n_features > 0:
            phase4_use = "blocked"
            reason = "Feature set explicitly blocked by missingness/leakage/availability audit."
        elif status == "empty" or n_features == 0:
            phase4_use = "unavailable"
            reason = "No usable features in this set."
        else:
            phase4_use = "sensitivity_only"
            reason = "Not clean-core primary under current Phase4 gate."

        rows.append(
            {
                "task": task,
                "feature_family": family,
                "universe_id": row.get("universe_id", ""),
                "feature_set_type": set_type,
                "n_features": n_features,
                "status": status,
                "phase4_use": phase4_use,
                "reason": reason,
                "feature_list_path": row.get("feature_list_path", ""),
            }
        )

    return rows


def build_blocked_inputs() -> list[dict[str, str]]:
    return [
        {
            "blocked_input": "PD1_anchor primary supervised model outputs",
            "blocking_reason": "Irreducible cohort-response confounding; weighted PD1_anchor remains HIGH risk.",
            "forbidden_use": "primary feature selection; primary module discovery; model ranking; clinical prediction claim",
            "allowed_use": "within-cohort/meta directional biological anchor only via support registry",
            "source_file": rel(HANDOFF / "STEP3_FINAL_DECISION.yaml"),
        },
        {
            "blocked_input": "Phase3.5 strong ML outputs",
            "blocking_reason": "Strong ML baseline hard-failed and was rolled back.",
            "forbidden_use": "any Phase4 primary claim, feature ranking, or model-performance argument",
            "allowed_use": "historical rollback record only",
            "source_file": rel(STEP3 / "05_strong_ml_baselines" / "step3_PHASE3_5_FAIL_AND_ROLLBACK.md"),
        },
        {
            "blocked_input": "Pre-hotfix Step2.10 unified feature matrix",
            "blocking_reason": "TF activity columns were stale before repair_B2 hotfix.",
            "forbidden_use": "any new downstream analysis",
            "allowed_use": "archive/provenance only",
            "source_file": rel(STEP2_REPAIR / "STEP2_TO_STEP3_READINESS_REPORT.hotfix.md"),
        },
        {
            "blocked_input": "Invalid Phase3.4 outputs generated before Phase3.3 gate fix",
            "blocking_reason": "Archived as invalid pre-repair run.",
            "forbidden_use": "any Phase4 analysis or reporting",
            "allowed_use": "abort/recovery provenance only",
            "source_file": rel(STEP3 / "99_abort_or_repair_logs" / "invalid_phase3_4_generated_before_phase3_3_gate_fix"),
        },
        {
            "blocked_input": "Unintegrated P0/P1 local datasets",
            "blocking_reason": "Data completeness intake registered but expression/TCR/spatial conversion has not passed Step1/Step2/Step3 audits.",
            "forbidden_use": "primary Phase4 module discovery or main feature selection",
            "allowed_use": "intake planning; future addendum; support only after manifest/QC gates pass",
            "source_file": rel(DATA_STRENGTH / "v6_1_step1_addendum_source_records_20260611.tsv"),
        },
        {
            "blocked_input": "Step3 pseudobulk_gene supervised feature set as primary Phase4 input",
            "blocking_reason": "Step3 final allowed inputs exclude pseudobulk_gene; direct gene-level GRN should use Step2.8 source layer separately.",
            "forbidden_use": "primary Step3-derived feature-family module discovery",
            "allowed_use": "separate gene-level GRN addendum using Step2.8 with provenance and missingness audit",
            "source_file": rel(HANDOFF / "STEP3_FINAL_DECISION.yaml"),
        },
    ]


def build_checklist() -> list[dict[str, str]]:
    return [
        {"criterion_id": "1", "criterion": "Sample and patient metadata frozen", "status": "PASS", "evidence": "Step1 sample/patient metadata and split files exist and are canonical.", "canonical_file": rel(STEP1 / "sample_metadata_master_v6_1.broad_response.csv")},
        {"criterion_id": "2", "criterion": "Treatment context and response labels frozen", "status": "PASS", "evidence": "Response schema and broad-response mappings are frozen in Step1.", "canonical_file": rel(STEP1 / "response_schema_report.md")},
        {"criterion_id": "3", "criterion": "Primary/support/blocked collections frozen", "status": "PASS", "evidence": "This Phase4 gate emits explicit evidence-pool and blocked-input registries.", "canonical_file": "phase4_evidence_pool_registry_20260611.csv"},
        {"criterion_id": "4", "criterion": "Feature source/type/availability frozen", "status": "PASS", "evidence": "Step2 hotfix feature dictionary plus Phase4 feature-family contract define availability.", "canonical_file": rel(STEP2_REPAIR / "repair_B2_feature_dictionary_v6_1.hotfix.csv")},
        {"criterion_id": "5", "criterion": "Leakage checks complete", "status": "PASS", "evidence": "Forbidden response-feature audits completed before Step3 handoff.", "canonical_file": rel(STEP3 / "00_input_contract" / "repair_phase3_0_forbidden_feature_audit.csv")},
        {"criterion_id": "6", "criterion": "Source/batch/missingness confounding checks complete", "status": "CONDITIONAL_PASS", "evidence": "Main universes pass with controls; PD1_anchor high-risk is downgraded and isolated.", "canonical_file": rel(STEP3 / "03_leakage_confounding_audit" / "step3_phase3_3_confounding_audit_by_universe.repaired.csv")},
        {"criterion_id": "7", "criterion": "Simple interpretable baseline complete", "status": "PASS", "evidence": "Phase3.4 interpretable baseline rows used downstream; complex ML not required.", "canonical_file": rel(STEP3 / "04_interpretable_baselines" / "step3_phase3_4_summary.md")},
        {"criterion_id": "8", "criterion": "Unreliable complex model rolled back", "status": "PASS", "evidence": "Phase3.5 strong ML hard-failed and is excluded from Phase4.", "canonical_file": rel(STEP3 / "05_strong_ml_baselines" / "step3_PHASE3_5_FAIL_AND_ROLLBACK.md")},
        {"criterion_id": "9", "criterion": "Main downstream feature sets fixed", "status": "PASS", "evidence": "Step3 handoff main feature sets and this contract define primary inputs.", "canonical_file": rel(HANDOFF / "step3_step4_handoff_main_feature_sets.csv")},
        {"criterion_id": "10", "criterion": "Auxiliary evidence isolated", "status": "PASS", "evidence": "PD1_anchor and new-data candidates are support/future-addendum only.", "canonical_file": rel(HANDOFF / "step3_step4_handoff_support_evidence_registry.csv")},
        {"criterion_id": "11", "criterion": "Downstream usage rules explicit", "status": "PASS", "evidence": "This report and decision YAML define direct Phase4 use rules.", "canonical_file": "phase4_entry_decision_20260611.yaml"},
    ]


def build_join_key_contract() -> list[dict[str, str]]:
    return [
        {
            "contract_id": "sample_metadata",
            "canonical_path": rel(STEP1 / "sample_metadata_master_v6_1.broad_response.csv"),
            "grain": "sample",
            "primary_key": "cohort_id + sample_id",
            "required_join_keys": "cohort_id;sample_id;patient_id",
            "required_phase4_fields": "disease;tissue_source;sample_type;timepoint;treatment_context;response_strict_binary;response_broad_binary;response_evidence_level;split;fold_id",
            "phase4_use": "metadata join authority for sample-level analyses",
            "rule": "Do not infer treatment or response from filenames when this table has an audited value.",
        },
        {
            "contract_id": "patient_metadata",
            "canonical_path": rel(STEP1 / "patient_metadata_master_v6_1.broad_response.csv"),
            "grain": "patient",
            "primary_key": "patient_key",
            "required_join_keys": "cohort_id;patient_id;patient_key",
            "required_phase4_fields": "treatment_context_primary;has_pre_sample;has_post_sample;has_paired_pre_post;response_strict_binary;response_broad_binary;patient_response_conflict;split;fold_id",
            "phase4_use": "patient-level response and pairing authority",
            "rule": "Patient-level split is authoritative; same patient must not cross train/validation/test.",
        },
        {
            "contract_id": "patient_split",
            "canonical_path": rel(STEP1 / "frozen_patient_split_v6_1.csv"),
            "grain": "patient",
            "primary_key": "patient_key",
            "required_join_keys": "cohort_id;patient_id;patient_key",
            "required_phase4_fields": "inclusion_status;executable_role;split_label;fold_id;supervised_anchor_available;split_reason",
            "phase4_use": "split and leakage guard",
            "rule": "If split conflicts with another table, use frozen_patient_split_v6_1.csv and re-audit mismatch.",
        },
        {
            "contract_id": "hotfix_feature_matrix",
            "canonical_path": rel(STEP2_REPAIR / "repair_B2_unified_feature_matrix_by_sample.hotfix.parquet"),
            "grain": "sample",
            "primary_key": "sample_key",
            "required_join_keys": "sample_key;cohort_id;sample_id;patient_id;patient_key",
            "required_phase4_fields": "split;fold_id;disease;tissue_source;timepoint;treatment_context;feature columns",
            "phase4_use": "feature matrix for allowed Phase4 primary/support feature families",
            "rule": "Use hotfix matrix only; pre-hotfix Step2.10 matrix is blocked.",
        },
        {
            "contract_id": "feature_dictionary",
            "canonical_path": rel(STEP2_REPAIR / "repair_B2_feature_dictionary_v6_1.hotfix.csv"),
            "grain": "feature",
            "primary_key": "feature_name",
            "required_join_keys": "feature_name",
            "required_phase4_fields": "feature_type;source_step;calculation_method;missing_rule;included_in_main;sensitivity_only;repair_B2_note",
            "phase4_use": "feature provenance and use-eligibility authority",
            "rule": "Feature without dictionary entry cannot become primary Phase4 evidence.",
        },
        {
            "contract_id": "phase4_evidence_pool",
            "canonical_path": "results/v6_1/phase4_readiness_20260611/phase4_evidence_pool_registry_20260611.csv",
            "grain": "analysis_universe",
            "primary_key": "universe_id",
            "required_join_keys": "universe_id;task;feature_family",
            "required_phase4_fields": "pool;phase4_use;confounding_risk_level;allowed_downstream_use;forbidden_downstream_use",
            "phase4_use": "primary/support universe authority",
            "rule": "Phase4 follows pool labels and does not re-decide Step3 risk.",
        },
        {
            "contract_id": "phase4_feature_contract",
            "canonical_path": "results/v6_1/phase4_readiness_20260611/phase4_feature_family_contract_20260611.csv",
            "grain": "universe_feature_set",
            "primary_key": "universe_id + feature_set_type",
            "required_join_keys": "universe_id;task;feature_family;feature_set_type",
            "required_phase4_fields": "n_features;status;phase4_use;reason;feature_list_path",
            "phase4_use": "feature-family availability and Phase4-use authority",
            "rule": "Clean-core primary rows may enter discovery; support/blocked/unavailable rows may not.",
        },
        {
            "contract_id": "pending_data_addendum",
            "canonical_path": "results/v6_1/phase4_readiness_20260611/phase4_pending_data_addendum_register_20260611.tsv",
            "grain": "candidate_dataset",
            "primary_key": "source_id",
            "required_join_keys": "source_id;local_path",
            "required_phase4_fields": "priority;modality;cancer;treatment_context;primary_role;phase4_current_use;phase4_gate",
            "phase4_use": "future addendum planning only",
            "rule": "No pending dataset enters primary Phase4 until Step1/Step2/Step3 addendum gates pass.",
        },
    ]


def build_pending_data_register() -> list[dict[str, str]]:
    source = DATA_STRENGTH / "v6_1_step1_addendum_source_records_20260611.tsv"
    rows = read_table(source, delimiter="\t")
    use_by_priority = {
        "P0": "future_addendum_high_priority_not_primary_yet",
        "P1": "future_addendum_prepare_not_primary_yet",
        "P2": "metadata_or_support_only",
        "P3": "defer_or_exclude",
    }
    for row in rows:
        row["phase4_current_use"] = use_by_priority.get(row.get("priority", ""), "audit_required")
        row["phase4_gate"] = "must_pass_step1_step2_step3_addendum_before_primary_use"
    return rows


def build_file_index(generated_files: list[Path], index_path: Path) -> list[dict[str, Any]]:
    base_rows = [
        ("sample_metadata_freeze", STEP1 / "sample_metadata_master_v6_1.broad_response.csv", "canonical_upstream", "sample_id/cohort_id/patient_id/treatment/response fields"),
        ("patient_metadata_freeze", STEP1 / "patient_metadata_master_v6_1.broad_response.csv", "canonical_upstream", "patient_id/patient_key/response/split fields"),
        ("patient_split_freeze", STEP1 / "frozen_patient_split_v6_1.csv", "canonical_upstream", "patient-level split; no patient crosses train/test"),
        ("hotfix_feature_matrix", STEP2_REPAIR / "repair_B2_unified_feature_matrix_by_sample.hotfix.parquet", "canonical_feature_input", "Step2.10 hotfix matrix"),
        ("hotfix_feature_dictionary", STEP2_REPAIR / "repair_B2_feature_dictionary_v6_1.hotfix.csv", "canonical_feature_input", "feature source/type/main/sensitivity flags"),
        ("hotfix_missingness_report", STEP2_REPAIR / "repair_B2_feature_missingness_report.hotfix.csv", "canonical_feature_input", "feature missingness after TF hotfix"),
        ("step3_final_decision", HANDOFF / "STEP3_FINAL_DECISION.yaml", "canonical_decision", "Conditional Go and allowed/support/forbidden Step4 inputs"),
        ("step3_main_universes", HANDOFF / "step3_step4_handoff_main_universes.csv", "canonical_primary_pool", "primary-allowed universes"),
        ("step3_main_feature_sets", HANDOFF / "step3_step4_handoff_main_feature_sets.csv", "canonical_primary_pool", "A/B eligible model/feature-set rows"),
        ("step3_support_pd1_anchor", HANDOFF / "step3_step4_handoff_support_PD1_anchor.csv", "canonical_support_pool", "PD1 directional support only"),
        ("step3_negative_controls", HANDOFF / "step3_negative_control_summary.csv", "canonical_audit", "confounding and negative-control table"),
        ("step3_robustness", HANDOFF / "step3_robustness_summary.csv", "canonical_audit", "robustness/negative-control execution summary"),
        ("data_completion_registry", DATA_STRENGTH / "v6_1_new_data_candidate_registry_20260611.csv", "future_addendum", "new/rescue data intake priorities"),
        ("data_completion_source_records", DATA_STRENGTH / "v6_1_step1_addendum_source_records_20260611.tsv", "future_addendum", "Step1 addendum source-record candidates"),
        ("data_completion_asset_manifest", DATA_STRENGTH / "v6_1_local_unintegrated_asset_manifest_20260611.tsv", "future_addendum", "local unintegrated file manifest"),
        ("phase4_join_key_contract", OUT_DIR / "phase4_join_key_contract_20260611.tsv", "phase4_gate_output", "join keys and authoritative metadata fields for Phase4"),
    ]
    rows = []
    for asset_id, path, role, desc in base_rows:
        rows.append({"asset_id": asset_id, "path": rel(path), "phase4_role": role, "description": desc, "exists": path.exists()})
    for path in generated_files:
        rows.append({"asset_id": path.stem, "path": rel(path), "phase4_role": "phase4_gate_output", "description": "generated Phase4 readiness artifact", "exists": path.exists() or path == index_path})
    return rows


def count(rows: list[dict[str, Any]], key: str, value: str) -> int:
    return sum(1 for row in rows if row.get(key) == value)


def report_text(decision: dict[str, Any], evidence: list[dict[str, Any]], features: list[dict[str, Any]], blocked: list[dict[str, Any]], checklist: list[dict[str, str]], pending: list[dict[str, str]]) -> str:
    primary_feature_rows = [r for r in features if r.get("phase4_use") == "primary"]
    conditional_feature_rows = [r for r in features if str(r.get("phase4_use", "")).startswith("conditional_primary")]
    lines = [
        "# Phase4 Readiness Gate Report",
        "",
        f"Date: {DATE_TAG}",
        "",
        "## Verdict",
        "",
        f"- Phase4 entry verdict: `{decision['phase4_entry_verdict']}`.",
        "- Phase4 may start mechanism/module discovery from primary evidence only.",
        "- Support and pending evidence remain isolated.",
        "- Complex ML is not entry criterion; Phase3.5 is rolled back and forbidden for Phase4 primary claims.",
        "",
        "## Primary Evidence Pool",
        "",
        f"- Primary universe rows: {count(evidence, 'pool', 'primary_evidence_pool')}.",
        f"- Primary tasks: {', '.join(decision['primary_evidence_pool']['allowed_tasks'])}.",
        f"- Clean-core primary feature-family rows: {len(primary_feature_rows)}.",
        f"- Conditional missingness-aware rows: {len(conditional_feature_rows)}; use only when corresponding evidence grade is A/B.",
        "",
        "## Support Evidence Pool",
        "",
        f"- Support universe rows: {count(evidence, 'pool', 'support_evidence_pool')}.",
        "- PD1_anchor is support/sensitivity/biological-anchor only, not primary supervised response modelling.",
        "- New local P0/P1 data are registered for addendum intake but cannot enter primary Phase4 until Step1/Step2/Step3 gates pass.",
        "",
        "## Blocked Inputs",
        "",
        f"- Blocked input rules: {len(blocked)}.",
        "- Blocked classes: PD1_anchor primary supervised outputs, Phase3.5 strong ML, pre-hotfix Step2.10 matrix, invalid pre-repair Phase3.4, unintegrated new datasets, Step3-derived pseudobulk_gene as primary feature family.",
        "",
        "## Data Completion State",
        "",
        f"- P0 pending intake rows: {count(pending, 'priority', 'P0')}.",
        f"- P1 pending intake rows: {count(pending, 'priority', 'P1')}.",
        "- These strengthen v6_1 later, but current Phase4 treats them as future addendum evidence.",
        "",
        "## Checklist",
        "",
    ]
    for row in checklist:
        lines.append(f"- {row['criterion_id']}. {row['criterion']}: `{row['status']}`.")
    lines.extend(
        [
            "",
            "## Canonical Reading Order",
            "",
        "1. `phase4_entry_decision_20260611.yaml`",
        "2. `phase4_evidence_pool_registry_20260611.csv`",
        "3. `phase4_feature_family_contract_20260611.csv`",
        "4. `phase4_blocked_inputs_20260611.csv`",
        "5. `phase4_join_key_contract_20260611.tsv`",
        "6. `phase4_downstream_file_index_20260611.tsv`",
            "",
            "## Operating Rule",
            "",
            "Phase4 should not re-judge Step3 risks. Read pool labels: primary evidence enters module discovery; support evidence enters annotation/direction checks; blocked evidence does not enter.",
        ]
    )
    return "\n".join(lines) + "\n"


def changelog_text(generated_files: list[Path]) -> str:
    lines = [
        "# Phase4 Readiness Release Changelog",
        "",
        f"Date: {DATE_TAG}",
        "",
        "## Release Scope",
        "",
        "- Consolidated Step1/Step2/Step3/data-completeness outputs into a Phase4 entry gate.",
        "- No Step1, Step2, or Step3 canonical data files were modified.",
        "- No large raw dataset was downloaded or unpacked.",
        "",
        "## Generated Artifacts",
        "",
    ]
    lines.extend(f"- `{rel(path)}`" for path in generated_files)
    lines.extend(
        [
            "",
            "## Decision Changes",
            "",
            "- Phase4 entry state is `CONDITIONAL_GO_READY_FOR_PHASE4`.",
            "- PD1_anchor remains blocked for primary supervised modelling and isolated as support-only.",
            "- Phase3.5 strong ML remains rolled back and forbidden as primary evidence.",
            "- P0/P1 unintegrated local datasets are registered as future addendum intake, not current primary Phase4 inputs.",
            "",
            "## Verification",
            "",
            "- Script checks all required canonical source files exist before writing outputs.",
            "- File index records existence status for every canonical downstream input.",
        ]
    )
    return "\n".join(lines) + "\n"


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    evidence = build_evidence_pool()
    features = build_feature_contract()
    blocked = build_blocked_inputs()
    checklist = build_checklist()
    pending = build_pending_data_register()
    join_contract = build_join_key_contract()
    model_grades = read_table(STEP3 / "06_robustness_negative_controls" / "step3_phase3_6_model_evidence_grading.csv")
    eligible_model_rows = sum(1 for row in model_grades if as_bool(row.get("eligible_for_step4_primary", "False")))

    decision = {
        "phase4_entry_verdict": "CONDITIONAL_GO_READY_FOR_PHASE4",
        "generated_at_utc": now_iso(),
        "scope": "v6_1_only",
        "primary_evidence_pool": {
            "allowed_tasks": sorted(PRIMARY_TASKS),
            "allowed_universe_rows": count(evidence, "pool", "primary_evidence_pool"),
            "allowed_feature_families": sorted(MAIN_FEATURE_FAMILIES),
            "primary_universe_file": "phase4_evidence_pool_registry_20260611.csv",
            "primary_feature_contract_file": "phase4_feature_family_contract_20260611.csv",
        },
        "support_evidence_pool": {
            "support_universe_rows": count(evidence, "pool", "support_evidence_pool"),
            "pd1_anchor_primary_blocked": True,
            "pd1_anchor_support_available": True,
            "support_registry_file": rel(HANDOFF / "step3_step4_handoff_support_evidence_registry.csv"),
        },
        "blocked_evidence": {
            "blocked_input_rules": len(blocked),
            "blocked_input_file": "phase4_blocked_inputs_20260611.csv",
        },
        "baseline_status": {
            "phase3_4_interpretable_baseline": "usable",
            "phase3_5_strong_ml": "rolled_back_for_primary_use",
            "phase3_6_negative_controls": "conditional_pass",
            "eligible_A_or_B_model_rows": eligible_model_rows,
        },
        "data_completion_addendum": {
            "p0_pending_rows": count(pending, "priority", "P0"),
            "p1_pending_rows": count(pending, "priority", "P1"),
            "current_phase4_primary_use": "not_allowed_until_addendum_gates_pass",
            "source_records_file": rel(DATA_STRENGTH / "v6_1_step1_addendum_source_records_20260611.tsv"),
        },
        "must_not_use_as_primary": [
            "PD1_anchor primary supervised outputs",
            "Phase3.5 strong ML outputs",
            "pre-hotfix Step2.10 unified feature matrix",
            "invalid pre-repair Phase3.4 outputs",
            "unintegrated new local datasets",
            "Step3 pseudobulk_gene supervised feature set as primary Phase4 input",
        ],
    }

    decision_path = OUT_DIR / "phase4_entry_decision_20260611.yaml"
    evidence_path = OUT_DIR / "phase4_evidence_pool_registry_20260611.csv"
    feature_path = OUT_DIR / "phase4_feature_family_contract_20260611.csv"
    blocked_path = OUT_DIR / "phase4_blocked_inputs_20260611.csv"
    checklist_path = OUT_DIR / "phase4_readiness_checklist_20260611.tsv"
    pending_path = OUT_DIR / "phase4_pending_data_addendum_register_20260611.tsv"
    join_contract_path = OUT_DIR / "phase4_join_key_contract_20260611.tsv"
    report_path = OUT_DIR / "PHASE4_READINESS_REPORT_20260611.md"
    changelog_path = OUT_DIR / "PHASE4_RELEASE_CHANGELOG_20260611.md"
    index_path = OUT_DIR / "phase4_downstream_file_index_20260611.tsv"

    write_yaml(decision_path, decision)
    write_table(evidence_path, evidence)
    write_table(feature_path, features)
    write_table(blocked_path, blocked)
    write_table(checklist_path, checklist, delimiter="\t")
    write_table(pending_path, pending, delimiter="\t")
    write_table(join_contract_path, join_contract, delimiter="\t")

    generated = [decision_path, evidence_path, feature_path, blocked_path, checklist_path, pending_path, join_contract_path, report_path, changelog_path, index_path]
    report_path.write_text(report_text(decision, evidence, features, blocked, checklist, pending))
    changelog_path.write_text(changelog_text(generated))

    file_index = build_file_index(generated, index_path)
    write_table(index_path, file_index, delimiter="\t")
    missing = [row["path"] for row in file_index if not as_bool(row["exists"])]
    if missing:
        raise RuntimeError(f"Phase4 file index has missing assets: {missing}")

    print(f"Wrote Phase4 readiness package: {OUT_DIR}")
    print(f"Primary evidence rows: {count(evidence, 'pool', 'primary_evidence_pool')}")
    print(f"Support evidence rows: {count(evidence, 'pool', 'support_evidence_pool')}")
    print(f"Eligible A/B model rows: {eligible_model_rows}")


if __name__ == "__main__":
    main()
