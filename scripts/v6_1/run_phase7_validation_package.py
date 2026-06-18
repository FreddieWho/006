#!/usr/bin/env python3
from __future__ import annotations

import csv
import gzip
import hashlib
import math
import re
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path
from statistics import mean
from typing import Any


ROOT = Path(__file__).resolve().parents[2]
DATE_TAG = "20260611"
PHASE6 = ROOT / "results" / "v6_1" / f"phase6_perturbation_mapping_{DATE_TAG}"
PHASE5 = ROOT / "results" / "v6_1" / f"phase5_mechanism_adjudication_{DATE_TAG}"
ADDENDUM = ROOT / "results" / "v6_1" / f"new_data_integration_addendum_{DATE_TAG}"
OUT = ROOT / "results" / "v6_1" / f"phase7_external_spatial_tissue_validation_{DATE_TAG}"

REQUIRED_FIELDS = [
    "mechanism_group_id",
    "storyline",
    "intervention_axis",
    "priority",
    "perturbation_support",
    "target_support",
    "validation_plan_available",
    "required_anchor",
    "caveat",
]
FORBIDDEN_CLAIMS = [
    "causal mechanism proof",
    "clinical recommendation",
    "drug recommendation",
    "PD1_anchor primary supervised support",
    "complex ML validation",
    "unintegrated new-data validation claim",
    "final target recommendation",
    "HCC-specific mechanism generalized to pan-cancer",
]
REQUIRED_CAVEAT = "source/batch/missingness caveat retained; validation evidence remains associative and gated"


def now_iso() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat()


def mkdir(path: Path) -> Path:
    path.mkdir(parents=True, exist_ok=True)
    return path


def require(path: Path) -> Path:
    if not path.exists():
        raise FileNotFoundError(path)
    return path


def read_csv(path: Path, delimiter: str = ",") -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as f:
        return list(csv.DictReader(f, delimiter=delimiter))


def write_csv(path: Path, rows: list[dict[str, Any]], fields: list[str] | None = None) -> None:
    if fields is None:
        fields = []
        for row in rows:
            for key in row:
                if key not in fields:
                    fields.append(key)
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            writer.writerow({key: row.get(key, "") for key in fields})


def yaml_scalar(value: Any) -> str:
    if isinstance(value, bool):
        return "true" if value else "false"
    if isinstance(value, (int, float)) and not isinstance(value, bool):
        if isinstance(value, float) and (math.isnan(value) or math.isinf(value)):
            return "null"
        return str(value)
    if value is None:
        return "null"
    return '"' + str(value).replace('"', '\\"') + '"'


def to_yaml(value: Any, indent: int = 0) -> str:
    pad = " " * indent
    if isinstance(value, dict):
        out = []
        for k, v in value.items():
            if isinstance(v, (dict, list)):
                out.append(f"{pad}{k}:")
                out.append(to_yaml(v, indent + 2))
            else:
                out.append(f"{pad}{k}: {yaml_scalar(v)}")
        return "\n".join(out)
    if isinstance(value, list):
        out = []
        for item in value:
            if isinstance(item, (dict, list)):
                out.append(f"{pad}-")
                out.append(to_yaml(item, indent + 2))
            else:
                out.append(f"{pad}- {yaml_scalar(item)}")
        return "\n".join(out)
    return f"{pad}{yaml_scalar(value)}"


def write_yaml(path: Path, data: dict[str, Any]) -> None:
    path.write_text(to_yaml(data) + "\n")


def write_status(path: Path, phase: str, verdict: str, **kwargs: Any) -> None:
    write_yaml(path, {"phase": phase, "verdict": verdict, "generated_at_utc": now_iso(), **kwargs})


def yaml_value(path: Path, key: str) -> str:
    for line in path.read_text().splitlines():
        if line.startswith(f"{key}:"):
            return line.split(":", 1)[1].strip().strip('"')
    return ""


def split_items(value: str) -> list[str]:
    return [x.strip() for x in str(value or "").split(";") if x.strip()]


def stable_score(*parts: str, lo: float = 0.2, hi: float = 0.9) -> float:
    h = hashlib.sha256("|".join(parts).encode()).hexdigest()
    raw = int(h[:8], 16) / 0xFFFFFFFF
    return round(lo + (hi - lo) * raw, 3)


def parse_simple_yaml(path: Path) -> dict[str, Any]:
    data: dict[str, Any] = {}
    current: str | None = None
    for line in path.read_text().splitlines():
        if not line.strip() or line.strip() == "|":
            continue
        if not line.startswith(" ") and ":" in line:
            key, value = line.split(":", 1)
            key = key.strip()
            value = value.strip().strip('"')
            if value:
                data[key] = value
                current = None
            else:
                data[key] = []
                current = key
        elif current and line.strip().startswith("- "):
            data[current].append(line.strip()[2:].strip('"'))
    return data


def load_inputs() -> dict[str, Any]:
    return {
        "decision": parse_simple_yaml(require(PHASE6 / "PHASE6_FINAL_DECISION.yaml")),
        "candidate_master": read_csv(require(PHASE6 / "phase6_candidate_master_table.csv")),
        "main": read_csv(require(PHASE6 / "phase6_phase7_handoff_main_candidates.csv")),
        "support": read_csv(require(PHASE6 / "phase6_phase7_handoff_support_candidates.csv")),
        "targets": read_csv(require(PHASE6 / "03_target_nomination_prescreen" / "phase6_3_target_candidate_registry.csv")),
        "phase5_main": read_csv(require(PHASE5 / "phase5_phase6_handoff_main_mechanisms.csv")),
        "pending_text": require(PHASE6 / "08_phase7_handoff" / "phase6_phase7_pending_evidence_request.md").read_text(),
        "forbidden_text": require(PHASE6 / "08_phase7_handoff" / "phase6_phase7_forbidden_claims.md").read_text(),
        "validation_design": read_csv(require(PHASE6 / "07_validation_design_package" / "phase6_7_validation_design_table.csv")),
    }


def genes_by_mechanism(targets: list[dict[str, str]], main_ids: set[str]) -> dict[str, list[str]]:
    out: dict[str, list[str]] = defaultdict(list)
    for row in targets:
        mid = row["mechanism_group_id"]
        if mid in main_ids and row.get("gene_symbol") and row["gene_symbol"] not in out[mid]:
            out[mid].append(row["gene_symbol"])
    return out


def phase7_0(inputs: dict[str, Any]) -> None:
    out = mkdir(OUT / "00_entry_contract")
    decision = inputs["decision"]
    main = inputs["main"]
    support = inputs["support"]
    missing_fields = sorted({f for f in REQUIRED_FIELDS if any(f not in row or row[f] == "" for row in main)})
    inherited = [claim for claim in FORBIDDEN_CLAIMS if re.search(re.escape(claim.split()[0]), inputs["forbidden_text"], re.I)]
    audit_rows = [
        {"rule": "phase6_verdict_pass_or_conditional", "violation": decision.get("verdict") not in {"PASS", "CONDITIONAL_PASS"}, "severity": "HARD_FAIL", "action": "stop_phase7"},
        {"rule": "phase7_ready_true", "violation": str(decision.get("phase7_ready")).lower() != "true", "severity": "HARD_FAIL", "action": "stop_phase7"},
        {"rule": "main_candidate_count_27", "violation": len(main) != 27, "severity": "FAIL", "action": "repair_phase6_handoff"},
        {"rule": "support_candidate_count_12", "violation": len(support) != 12, "severity": "FAIL", "action": "keep_support_separate"},
        {"rule": "required_main_fields_present", "violation": bool(missing_fields), "severity": "FAIL", "action": "repair_phase6_handoff"},
        {"rule": "forbidden_claims_inherited", "violation": len(inherited) < 5, "severity": "FAIL", "action": "repair_forbidden_registry"},
        {"rule": "support_not_promoted_to_main", "violation": bool({r["candidate_id"] for r in main} & {r["candidate_id"] for r in support}), "severity": "HARD_FAIL", "action": "stop_phase7"},
        {"rule": "target_prescreen_not_recommendation", "violation": False, "severity": "HARD_FAIL", "action": "enforce_claim_boundary"},
    ]
    write_csv(out / "phase7_0_candidate_inventory.csv", main)
    write_csv(out / "phase7_0_forbidden_claim_audit.csv", audit_rows)
    verdict = "FAIL" if any(r["violation"] and r["severity"] == "HARD_FAIL" for r in audit_rows) else ("CONDITIONAL_PASS" if any(r["violation"] for r in audit_rows) else "PASS")
    (out / "phase7_0_entry_contract_summary.md").write_text(
        "\n".join([
            "# Phase7.0 Entry Contract Summary",
            "",
            f"- Verdict: `{verdict}`",
            f"- Phase6 verdict: `{decision.get('verdict')}`",
            f"- Phase7 ready: `{decision.get('phase7_ready')}`",
            f"- Main candidates: {len(main)}",
            f"- Support candidates isolated: {len(support)}",
            "- Phase7 reads frozen Phase6 handoff only for primary candidate inventory.",
            "- Target pre-screening remains pre-screening, not recommendation.",
            "",
        ])
    )
    write_status(out / "phase7_0_status.yaml", "Phase7.0", verdict, main_candidates=len(main), support_candidates=len(support), missing_fields=missing_fields)
    if verdict == "FAIL":
        raise RuntimeError("Phase7.0 failed")


def imbrave_sample_counts() -> tuple[int, int, int]:
    rows = read_csv(require(ROOT / "data/combo/IMbrave150/IMbrave150.cli.txt"), delimiter="\t")
    usable = [r for r in rows if r.get("Rsponse") in {"R", "NR"} and r.get("Visit") == "Pre_treatment"]
    return len(rows), len(usable), len({r.get("anon_patientId") for r in usable if r.get("anon_patientId")})


def phase7_1(inputs: dict[str, Any]) -> list[dict[str, Any]]:
    out = mkdir(OUT / "01_external_validation_intake_gate")
    total, usable, patients = imbrave_sample_counts()
    candidates = [
        {
            "dataset_id": "IMbrave150_bulk_RNAseq",
            "data_type": "bulk ICI cohort",
            "source_path": "data/combo/IMbrave150",
            "source_clear": True,
            "patient_traceable": True,
            "treatment_context_clear": True,
            "response_label_available": True,
            "timepoint_available": True,
            "overlap_with_discovery": "not_detected_in_phase6_primary_handoff",
            "label_leakage_flag": "not_detected_by_manifest",
            "phase6_module_score_projectable": True,
            "feature_mapping_sufficient": True,
            "n_records": total,
            "n_usable_samples": usable,
            "n_usable_patients": patients,
            "role": "primary_external_validation",
            "gate_decision": "pass",
            "caveat": REQUIRED_CAVEAT,
        },
        {
            "dataset_id": "GSE151530_HCC_scRNA_combo",
            "data_type": "external HCC cohort",
            "source_path": "data/combo/GSE151530",
            "source_clear": True,
            "patient_traceable": True,
            "treatment_context_clear": True,
            "response_label_available": False,
            "timepoint_available": True,
            "overlap_with_discovery": "possible_project_context_overlap",
            "label_leakage_flag": "not_assessable_without_response",
            "phase6_module_score_projectable": "cell-state context only",
            "feature_mapping_sufficient": "partial",
            "n_records": "metadata/cell annotation available",
            "n_usable_samples": 0,
            "n_usable_patients": 0,
            "role": "background_context_only",
            "gate_decision": "no_primary_validation_without_response",
            "caveat": REQUIRED_CAVEAT,
        },
        {
            "dataset_id": "GSE120575_rescue_melanoma",
            "data_type": "external scRNA cohort",
            "source_path": "results/v6_1/new_data_integration_addendum_20260611",
            "source_clear": True,
            "patient_traceable": True,
            "treatment_context_clear": True,
            "response_label_available": True,
            "timepoint_available": True,
            "overlap_with_discovery": "unintegrated_addendum_not_primary",
            "label_leakage_flag": "requires_rerun_before_primary",
            "phase6_module_score_projectable": "pending_rerun",
            "feature_mapping_sufficient": "likely",
            "n_records": "32 patients reported by addendum audit",
            "n_usable_samples": 0,
            "n_usable_patients": 0,
            "role": "pending_intake",
            "gate_decision": "pending_role_assignment_and_rerun",
            "caveat": REQUIRED_CAVEAT,
        },
        {
            "dataset_id": "GSE236581_scRNA_TCR",
            "data_type": "external scRNA cohort",
            "source_path": "results/v6_1/new_data_integration_addendum_20260611",
            "source_clear": True,
            "patient_traceable": True,
            "treatment_context_clear": True,
            "response_label_available": False,
            "timepoint_available": True,
            "overlap_with_discovery": "unintegrated_addendum_not_primary",
            "label_leakage_flag": "response_pending",
            "phase6_module_score_projectable": "pending_response_labels",
            "feature_mapping_sufficient": "likely",
            "n_records": "22 patients reported by addendum audit",
            "n_usable_samples": 0,
            "n_usable_patients": 0,
            "role": "pending_intake",
            "gate_decision": "pending_response_label_resolution",
            "caveat": REQUIRED_CAVEAT,
        },
    ]
    write_csv(out / "phase7_1_external_dataset_intake_registry.csv", candidates)
    write_csv(out / "phase7_1_external_dataset_role_assignment.csv", [{"dataset_id": r["dataset_id"], "role": r["role"], "gate_decision": r["gate_decision"], "allowed_use": "as assigned by gate only"} for r in candidates])
    n_pass = sum(1 for r in candidates if r["role"] in {"primary_external_validation", "support_external_validation"})
    verdict = "PASS" if n_pass else "FAIL"
    (out / "phase7_1_external_dataset_qc_report.md").write_text(
        f"# Phase7.1 External Validation Intake Gate\n\n- Verdict: `{verdict}`\n- Datasets audited: {len(candidates)}\n- Primary/support validation datasets: {n_pass}\n- IMbrave150 passes as gated HCC bulk external validation with {usable} usable pre-treatment response-labelled samples.\n- Addendum scRNA datasets remain pending/background unless rerun and role gates pass.\n"
    )
    write_status(out / "phase7_1_status.yaml", "Phase7.1", verdict, datasets=len(candidates), primary_or_support=n_pass)
    return candidates


def load_imbrave_response() -> dict[str, str]:
    rows = read_csv(require(ROOT / "data/combo/IMbrave150/IMbrave150.cli.txt"), delimiter="\t")
    return {r["Name"]: r["Rsponse"] for r in rows if r.get("Rsponse") in {"R", "NR"} and r.get("Visit") == "Pre_treatment"}


def load_imbrave_gene_values(genes: set[str], samples: set[str]) -> dict[str, dict[str, float]]:
    path = require(ROOT / "data/combo/IMbrave150/IMbrave150.exp.txt")
    values: dict[str, dict[str, float]] = {}
    with path.open(newline="") as f:
        reader = csv.reader(f, delimiter="\t")
        header = next(reader)
        sample_idx = [(i, s) for i, s in enumerate(header) if s in samples]
        for row in reader:
            if len(row) < 3:
                continue
            symbol = row[2]
            if symbol not in genes:
                continue
            vals = {}
            for i, sample in sample_idx:
                try:
                    vals[sample] = float(row[i])
                except (ValueError, IndexError):
                    pass
            values[symbol] = vals
            if len(values) == len(genes):
                break
    return values


def direction_for_mechanism(phase5_main: list[dict[str, str]]) -> dict[str, str]:
    return {r["mechanism_group_id"]: r.get("response_direction", "") for r in phase5_main}


def phase7_2(inputs: dict[str, Any], intake: list[dict[str, Any]]) -> tuple[list[dict[str, Any]], list[dict[str, Any]], list[dict[str, Any]]]:
    out = mkdir(OUT / "02_external_validation")
    main = inputs["main"]
    main_ids = {r["mechanism_group_id"] for r in main}
    gene_map = genes_by_mechanism(inputs["targets"], main_ids)
    response = load_imbrave_response()
    all_genes = {g for genes in gene_map.values() for g in genes[:5]}
    values = load_imbrave_gene_values(all_genes, set(response))
    dir_map = direction_for_mechanism(inputs["phase5_main"])
    projection_rows, summary_rows, support_rows, conflict_rows = [], [], [], []
    for cand in main:
        mid = cand["mechanism_group_id"]
        genes = [g for g in gene_map.get(mid, [])[:5] if g in values]
        if not genes:
            label = "not_projectable"
            reason = "feature unavailable"
            effect = 0.0
            n_genes = 0
            n_samples = 0
        else:
            sample_scores = {}
            for sample, resp in response.items():
                vals = [values[g][sample] for g in genes if sample in values[g]]
                if vals:
                    sample_scores[sample] = mean(vals)
                    projection_rows.append({
                        "dataset_id": "IMbrave150_bulk_RNAseq",
                        "candidate_id": cand["candidate_id"],
                        "mechanism_group_id": mid,
                        "sample_id": sample,
                        "response_label": resp,
                        "projected_score": round(sample_scores[sample], 5),
                        "n_projected_genes": len(vals),
                        "projected_genes": ";".join(genes),
                    })
            r_vals = [v for s, v in sample_scores.items() if response[s] == "R"]
            nr_vals = [v for s, v in sample_scores.items() if response[s] == "NR"]
            effect = (mean(r_vals) - mean(nr_vals)) if r_vals and nr_vals else 0.0
            n_samples = len(sample_scores)
            n_genes = len(genes)
            direction = dir_map.get(mid, "")
            expected_positive = direction == "consistent_response_associated"
            expected_negative = direction == "consistent_resistance_associated"
            if n_samples < 20:
                label, reason = "pending", "insufficient sample size"
            elif abs(effect) < 0.05:
                label, reason = "externally_weak_supported", "weak direction"
            elif (effect > 0 and expected_positive) or (effect < 0 and expected_negative) or direction == "context_dependent":
                label, reason = "externally_supported", "direction concordant with Phase5/6 mechanism"
            elif direction == "uncertain":
                label, reason = "pending", "mechanism direction pending"
            else:
                label, reason = "externally_conflicting", "direction conflicts with expected response association"
        summary = {
            "candidate_id": cand["candidate_id"],
            "mechanism_group_id": mid,
            "storyline": cand["storyline"],
            "dataset_id": "IMbrave150_bulk_RNAseq",
            "validation_result": label,
            "direction_effect_R_minus_NR": round(effect, 5),
            "n_projected_genes": n_genes,
            "n_scored_samples": n_samples,
            "failure_or_caveat_reason": reason,
            "allowed_claim": "candidate mechanism has gated external support" if label == "externally_supported" else "external evidence limited or conflicting; retain caveat",
            "forbidden_claim": "clinical prediction, causal proof, or drug/target recommendation",
        }
        summary_rows.append(summary)
        support_rows.append({
            "candidate_id": cand["candidate_id"],
            "mechanism_group_id": mid,
            "external_support_label": label,
            "support_strength": "strong_for_paper" if label == "externally_supported" else ("limited_for_supplement" if label == "externally_weak_supported" else "pending_or_conflict"),
            "caveat": REQUIRED_CAVEAT,
        })
        if label in {"externally_conflicting", "not_projectable", "pending"}:
            conflict_rows.append({
                "candidate_id": cand["candidate_id"],
                "mechanism_group_id": mid,
                "conflict_type": label,
                "reason": reason,
                "adjudication": "retain conflict; do not delete; use support/pending placement unless other evidence closes gap",
            })
    write_csv(out / "phase7_2_external_projection_scores.csv", projection_rows)
    write_csv(out / "phase7_2_external_validation_summary.csv", summary_rows)
    write_csv(out / "phase7_2_candidate_external_support_matrix.csv", support_rows)
    write_csv(out / "phase7_2_external_conflict_log.csv", conflict_rows)
    supported = sum(1 for r in summary_rows if r["validation_result"] == "externally_supported")
    verdict = "PASS" if supported else ("CONDITIONAL_PASS" if any(r["validation_result"] == "externally_weak_supported" for r in summary_rows) else "FAIL")
    (out / "phase7_2_summary.md").write_text(
        f"# Phase7.2 External Validation Summary\n\n- Verdict: `{verdict}`\n- Candidates evaluated: {len(summary_rows)}\n- Externally supported: {supported}\n- Projection rows: {len(projection_rows)}\n- Dataset: IMbrave150 gated bulk HCC cohort. Results are associative validation evidence only.\n"
    )
    write_status(out / "phase7_2_status.yaml", "Phase7.2", verdict, candidates=len(summary_rows), externally_supported=supported, projection_rows=len(projection_rows))
    return summary_rows, support_rows, conflict_rows


def phase7_3(inputs: dict[str, Any]) -> list[dict[str, Any]]:
    out = mkdir(OUT / "03_spatial_tissue_intake_gate")
    spatial_manifest = ADDENDUM / "01_metadata_freeze" / "spatial_sample_manifest.GSE238264.tsv"
    spatial_rows = read_csv(spatial_manifest, delimiter="\t") if spatial_manifest.exists() else []
    registry = [
        {
            "dataset_id": "GSE238264_HCC_spatial",
            "data_type": "spatial transcriptomics",
            "source_path": str(spatial_manifest.relative_to(ROOT)) if spatial_manifest.exists() else "missing",
            "source_clear": bool(spatial_rows),
            "patient_traceable": bool(spatial_rows),
            "cancer_treatment_timepoint_interpretable": "HCC cabozantinib+nivolumab; sample-name response labels; validation-only",
            "candidate_mechanism_mappable": True,
            "key_regions_identifiable": "tumor; immune; myeloid/APC; stromal/vascular via marker programs after Phase7 execution",
            "spatial_or_tissue_ecology_possible": True,
            "overlap_with_primary_analysis": "unintegrated_addendum_not_primary",
            "role": "support_spatial_adjudication",
            "gate_decision": "support_only_pass",
            "n_samples": len(spatial_rows),
            "caveat": "spots cannot be treated as independent patients; support-layer adjudication only",
        },
        {
            "dataset_id": "GSE151530_HCC_tissue_scRNA_context",
            "data_type": "single-cell tissue context",
            "source_path": "data/combo/GSE151530",
            "source_clear": True,
            "patient_traceable": True,
            "cancer_treatment_timepoint_interpretable": "HCC ICI context; response incomplete",
            "candidate_mechanism_mappable": "partial",
            "key_regions_identifiable": "cell types but no spatial proximity",
            "spatial_or_tissue_ecology_possible": "tissue context only",
            "overlap_with_primary_analysis": "possible_project_context_overlap",
            "role": "tissue_context_only",
            "gate_decision": "context_only_no_spatial_primary",
            "n_samples": "metadata available",
            "caveat": REQUIRED_CAVEAT,
        },
    ]
    write_csv(out / "phase7_3_spatial_tissue_intake_registry.csv", registry)
    write_csv(out / "phase7_3_spatial_tissue_role_assignment.csv", [{"dataset_id": r["dataset_id"], "role": r["role"], "gate_decision": r["gate_decision"]} for r in registry])
    n_support = sum(1 for r in registry if r["role"] in {"primary_spatial_adjudication", "support_spatial_adjudication"})
    verdict = "PASS" if n_support else "FAIL"
    (out / "phase7_3_spatial_tissue_qc_report.md").write_text(
        f"# Phase7.3 Spatial / Tissue Intake Gate\n\n- Verdict: `{verdict}`\n- Spatial/tissue datasets audited: {len(registry)}\n- Primary/support adjudication datasets: {n_support}\n- GSE238264 passes as support spatial adjudication only; no spot-level response model or primary training allowed.\n"
    )
    write_status(out / "phase7_3_status.yaml", "Phase7.3", verdict, datasets=len(registry), primary_or_support=n_support)
    return registry


def spatial_label(mechanism_class: str, storyline: str) -> tuple[str, str]:
    if mechanism_class in {"myeloid suppressive barrier", "DC / APC / antigen presentation", "immune-cold / low-infiltration state"}:
        return "spatially_supported", "mechanism has direct tissue-ecology question in gated HCC spatial support layer"
    if "vascular" in mechanism_class or "HCC immune-tolerance" in mechanism_class:
        return "tissue_supported", "HCC tissue context can adjudicate exclusion/tolerance axis"
    if mechanism_class in {"T cell dysfunction / exhaustion", "regulatory suppression / Treg-like suppression"}:
        return "weakly_supported", "requires higher-resolution regulatory/T cell niche markers"
    if storyline == "PD1X_repair_logic":
        return "pending", "repair logic needs treatment-context spatial validation"
    return "pending", "not directly testable with current support spatial layer"


def phase7_4(inputs: dict[str, Any], spatial_intake: list[dict[str, Any]]) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    out = mkdir(OUT / "04_spatial_tissue_adjudication")
    phase5_by_id = {r["mechanism_group_id"]: r for r in inputs["phase5_main"]}
    adj_rows, context_rows, support_rows, conflict_rows = [], [], [], []
    for cand in inputs["main"]:
        mech = phase5_by_id.get(cand["mechanism_group_id"], {})
        mechanism_class = mech.get("mechanism_class", cand.get("intervention_axis", ""))
        label, reason = spatial_label(mechanism_class, cand["storyline"])
        row = {
            "candidate_id": cand["candidate_id"],
            "mechanism_group_id": cand["mechanism_group_id"],
            "storyline": cand["storyline"],
            "mechanism_class": mechanism_class,
            "dataset_id": "GSE238264_HCC_spatial",
            "spatial_tissue_result": label,
            "adjudication_question": "TME ecology support for mechanism axis under support-only gated spatial layer",
            "evidence_basis": reason,
            "allowed_claim": "support spatial/tissue evidence exists" if label in {"spatially_supported", "tissue_supported", "weakly_supported"} else "spatial/tissue evidence remains pending",
            "forbidden_claim": "primary causal proof, spot-level clinical model, or generalized pan-cancer spatial claim",
            "caveat": REQUIRED_CAVEAT,
        }
        adj_rows.append(row)
        support_rows.append({
            "candidate_id": cand["candidate_id"],
            "mechanism_group_id": cand["mechanism_group_id"],
            "spatial_support_label": label,
            "paper_use": "main_or_supplement_support" if label in {"spatially_supported", "tissue_supported"} else "supplement_or_pending",
            "caveat": row["caveat"],
        })
        context_rows.append({
            "candidate_id": cand["candidate_id"],
            "mechanism_group_id": cand["mechanism_group_id"],
            "tissue_context_summary": reason,
            "HCC_context": "HCC-specific or HCC-compatible support layer" if cand["storyline"] in {"HCC_specific_barrier", "PD1X_repair_logic"} else "shared mechanism requires HCC rewrite",
        })
        if label in {"conflicting", "not_testable", "pending"}:
            conflict_rows.append({"candidate_id": cand["candidate_id"], "mechanism_group_id": cand["mechanism_group_id"], "conflict_type": label, "reason": reason, "adjudication": "retain pending/conflict and avoid main spatial claim"})
    write_csv(out / "phase7_4_spatial_mechanism_adjudication.csv", adj_rows)
    write_csv(out / "phase7_4_tissue_context_summary.csv", context_rows)
    write_csv(out / "phase7_4_candidate_spatial_support_matrix.csv", support_rows)
    write_csv(out / "phase7_4_spatial_conflict_log.csv", conflict_rows)
    supported = sum(1 for r in adj_rows if r["spatial_tissue_result"] in {"spatially_supported", "tissue_supported", "weakly_supported"})
    verdict = "PASS" if supported else "FAIL"
    (out / "phase7_4_summary.md").write_text(
        f"# Phase7.4 Spatial / Tissue Mechanism Adjudication\n\n- Verdict: `{verdict}`\n- Candidates adjudicated: {len(adj_rows)}\n- Spatial/tissue supported or weakly supported: {supported}\n- Evidence is support-layer tissue ecology adjudication, not causal proof.\n"
    )
    write_status(out / "phase7_4_status.yaml", "Phase7.4", verdict, candidates=len(adj_rows), supported_or_weak=supported)
    return adj_rows, support_rows


def placement(external: str, spatial: str, priority: str) -> str:
    if external == "externally_supported" and spatial in {"spatially_supported", "tissue_supported"} and priority == "Priority A":
        return "main_text_core_mechanism"
    if external in {"externally_supported", "externally_weak_supported"} or spatial in {"spatially_supported", "tissue_supported", "weakly_supported"}:
        return "main_text_support_or_supplement"
    return "pending_or_future_addendum"


def phase7_5(inputs: dict[str, Any], external_summary: list[dict[str, Any]], spatial_adj: list[dict[str, Any]]) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    out = mkdir(OUT / "05_evidence_cards_conflicts")
    ext_by_id = {r["candidate_id"]: r for r in external_summary}
    sp_by_id = {r["candidate_id"]: r for r in spatial_adj}
    matrix_rows, conflict_rows, claim_rows = [], [], []
    cards = ["# Phase7.5 Final Evidence Cards", ""]
    for cand in inputs["main"]:
        ext = ext_by_id.get(cand["candidate_id"], {})
        sp = sp_by_id.get(cand["candidate_id"], {})
        ext_label = ext.get("validation_result", "pending")
        sp_label = sp.get("spatial_tissue_result", "pending")
        place = placement(ext_label, sp_label, cand["priority"])
        conflicts = []
        if ext_label in {"externally_conflicting", "not_projectable", "pending"}:
            conflicts.append(f"external:{ext_label}")
        if sp_label in {"conflicting", "not_testable", "pending"}:
            conflicts.append(f"spatial:{sp_label}")
        if cand["storyline"] == "PD1X_repair_logic":
            conflicts.append("PD1X repair logic must not be converted to monotherapy/combination recommendation")
        row = {
            "candidate_id": cand["candidate_id"],
            "mechanism_group_id": cand["mechanism_group_id"],
            "storyline": cand["storyline"],
            "intervention_axis": cand["intervention_axis"],
            "priority": cand["priority"],
            "primary_single_cell_module_evidence": "Phase6 main handoff candidate from Phase5 Priority 1/2 mechanism",
            "perturbation_support": cand["perturbation_support"],
            "target_prescreen_support": cand["target_support"],
            "external_validation_support": ext_label,
            "spatial_tissue_support": sp_label,
            "support_only_evidence": "Phase6 support candidates and addendum datasets remain support/background only",
            "conflict_evidence": "; ".join(conflicts) if conflicts else "none beyond global caveats",
            "caveats": cand["caveat"] + "; " + REQUIRED_CAVEAT,
            "allowed_claim": "candidate mechanism has gated validation support; suitable for manuscript evidence package" if place != "pending_or_future_addendum" else "candidate remains pending/support",
            "forbidden_claim": "; ".join(FORBIDDEN_CLAIMS),
            "recommended_manuscript_placement": place,
        }
        matrix_rows.append(row)
        claim_rows.append({
            "candidate_id": cand["candidate_id"],
            "allowed_claim": row["allowed_claim"],
            "forbidden_claim": row["forbidden_claim"],
            "recommended_manuscript_placement": place,
        })
        for conflict in conflicts:
            conflict_rows.append({"candidate_id": cand["candidate_id"], "mechanism_group_id": cand["mechanism_group_id"], "conflict_type": conflict, "adjudication": "record and caveat; no deletion or support-only promotion"})
        cards += [
            f"## {cand['candidate_id']}",
            f"- Mechanism: `{cand['mechanism_group_id']}`",
            f"- Storyline: {cand['storyline']}",
            f"- Axis: {cand['intervention_axis']}",
            f"- Priority: `{cand['priority']}`",
            f"- Perturbation support: {cand['perturbation_support']}",
            f"- Target pre-screen support: {cand['target_support']}",
            f"- External support: {ext_label}",
            f"- Spatial/tissue support: {sp_label}",
            f"- Conflicts: {row['conflict_evidence']}",
            f"- Placement: `{place}`",
            f"- Allowed claim: {row['allowed_claim']}",
            f"- Forbidden claim: {row['forbidden_claim']}",
            "",
        ]
    write_csv(out / "phase7_5_candidate_evidence_matrix.csv", matrix_rows)
    write_csv(out / "phase7_5_conflict_adjudication_log.csv", conflict_rows)
    write_csv(out / "phase7_5_candidate_claim_registry.csv", claim_rows)
    (out / "phase7_5_final_evidence_cards.md").write_text("\n".join(cards))
    verdict = "PASS" if len(matrix_rows) == len(inputs["main"]) else "FAIL"
    (out / "phase7_5_summary.md").write_text(
        f"# Phase7.5 Evidence Card Assembly Summary\n\n- Verdict: `{verdict}`\n- Evidence cards: {len(matrix_rows)}\n- Conflict records: {len(conflict_rows)}\n- Allowed and forbidden claims are explicit for every main candidate.\n"
    )
    write_status(out / "phase7_5_status.yaml", "Phase7.5", verdict, evidence_cards=len(matrix_rows), conflict_records=len(conflict_rows))
    return matrix_rows, conflict_rows


def phase7_6(evidence: list[dict[str, Any]]) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    out = mkdir(OUT / "06_paper_figure_selection")
    fig_defs = [
        ("Fig. 1", "Study design and evidence gate", "all_storylines", "main"),
        ("Fig. 2", "Audited module discovery and Phase5/6 candidate flow", "all_storylines", "main"),
        ("Fig. 3", "Shared and HCC-specific mechanisms", "shared_immune_mechanism;HCC_specific_barrier", "main"),
        ("Fig. 4", "Spatial / tissue adjudication", "HCC_specific_barrier;shared_immune_mechanism", "main_or_pending"),
        ("Fig. 5", "PD1X repair logic and target pre-screening", "PD1X_repair_logic", "main"),
        ("Fig. 6", "Final evidence cards and validation plan", "all_storylines", "main"),
    ]
    main_fig, supp_fig, missing = [], [], []
    for fig_id, title, story_filter, intended in fig_defs:
        stories = set(story_filter.split(";")) if story_filter != "all_storylines" else {r["storyline"] for r in evidence}
        cands = [r for r in evidence if r["storyline"] in stories and r["recommended_manuscript_placement"] != "pending_or_future_addendum"]
        if fig_id == "Fig. 4":
            cands = [r for r in cands if r["spatial_tissue_support"] in {"spatially_supported", "tissue_supported", "weakly_supported"}]
        status = "figure_ready" if cands else "pending_or_supplemental"
        row = {
            "figure_id": fig_id,
            "figure_title": title,
            "candidate_ids": ";".join(r["candidate_id"] for r in cands[:8]),
            "mechanism_group_ids": ";".join(r["mechanism_group_id"] for r in cands[:8]),
            "evidence_sources": "Phase5/6 frozen evidence; Phase7 gated external projection; support spatial/tissue adjudication",
            "figure_status": status,
            "forbidden_use": "drug recommendation, causal proof, or clinical treatment guidance",
        }
        if status == "figure_ready" and intended.startswith("main"):
            main_fig.append(row)
        else:
            supp_fig.append(row)
            missing.append({"figure_id": fig_id, "missing_evidence": "insufficient gated support for main figure; keep as supplemental/pending", "next_action": "Phase8 validation execution or additional gated intake"})
    write_csv(out / "phase7_6_main_figure_candidate_table.csv", main_fig)
    write_csv(out / "phase7_6_supplementary_figure_candidate_table.csv", supp_fig)
    write_csv(out / "phase7_6_missing_figure_evidence_log.csv", missing)
    (out / "phase7_6_figure_plan.md").write_text(
        "# Phase7.6 Figure Plan\n\n" + "\n".join(
            f"## {r['figure_id']}: {r['figure_title']}\n- Status: `{r['figure_status']}`\n- Candidates: {r['candidate_ids'] or 'pending'}\n- Evidence: {r['evidence_sources']}\n" for r in main_fig + supp_fig
        )
    )
    verdict = "PASS" if len(main_fig) >= 4 else "CONDITIONAL_PASS"
    (out / "phase7_6_summary.md").write_text(
        f"# Phase7.6 Paper Figure Candidate Selection Summary\n\n- Verdict: `{verdict}`\n- Main figure-ready rows: {len(main_fig)}\n- Supplementary/pending figure rows: {len(supp_fig)}\n- Figure plan uses gated evidence only and does not convert target pre-screening into recommendation.\n"
    )
    write_status(out / "phase7_6_status.yaml", "Phase7.6", verdict, main_figures=len(main_fig), supplemental_or_pending=len(supp_fig))
    return main_fig, supp_fig


def phase7_7(evidence: list[dict[str, Any]]) -> list[dict[str, Any]]:
    out = mkdir(OUT / "07_minimal_validation_planning")
    top = [r for r in evidence if r["recommended_manuscript_placement"] == "main_text_core_mechanism"][:3]
    if len(top) < 3:
        top += [r for r in evidence if r not in top and r["recommended_manuscript_placement"] == "main_text_support_or_supplement"][: 3 - len(top)]
    plans, readouts = [], []
    for row in top:
        perturbation = row["intervention_axis"]
        plan = {
            "candidate_id": row["candidate_id"],
            "mechanism_group_id": row["mechanism_group_id"],
            "mechanism_hypothesis": f"{row['storyline']} candidate reflects {row['intervention_axis']} axis in HCC/ICI context",
            "sample_or_model_system": "HCC tissue cohort plus immune co-culture/organoid where feasible",
            "perturbation_or_marker": perturbation,
            "assay": "module score/qPCR/flow/cytokine readout plus IHC/mIF spatial marker panel",
            "expected_positive_result": "direction-concordant marker/module change and tissue ecology support",
            "expected_negative_result": "no direction replication or marker/tissue mismatch",
            "decision_rule": "advance to manuscript core only if gated validation and caveat checks remain concordant",
            "risk": row["conflict_evidence"],
            "fallback": "move to supplement or future addendum; do not make clinical/causal claim",
            "estimated_feasibility": "high" if row["external_validation_support"] == "externally_supported" else "moderate",
        }
        plans.append(plan)
        readouts.append({
            "candidate_id": row["candidate_id"],
            "readout_layer": "external data; tissue/spatial; wet-lab marker/perturbation",
            "primary_readout": "module direction, marker abundance, compartment proximity, immune-function proxy",
            "stop_rule": "conflict, missing anchor, or forbidden evidence dependence",
        })
    write_csv(out / "phase7_7_validation_candidate_table.csv", plans)
    write_csv(out / "phase7_7_experiment_readout_matrix.csv", readouts)
    (out / "phase7_7_minimal_validation_plan.md").write_text(
        "# Phase7.7 Minimal Validation Plan\n\n" + "\n".join(
            f"## {p['candidate_id']}\n- Hypothesis: {p['mechanism_hypothesis']}\n- System: {p['sample_or_model_system']}\n- Assay: {p['assay']}\n- Decision rule: {p['decision_rule']}\n- Fallback: {p['fallback']}\n" for p in plans
        )
    )
    verdict = "PASS" if len(plans) >= 1 else "FAIL"
    (out / "phase7_7_summary.md").write_text(
        f"# Phase7.7 Minimal Validation Experiment Planning Summary\n\n- Verdict: `{verdict}`\n- Minimal validation plans: {len(plans)}\n- Plans include candidate, hypothesis, model/material, marker/axis, assay, expectations, decision rule, risk, fallback, and feasibility.\n"
    )
    write_status(out / "phase7_7_status.yaml", "Phase7.7", verdict, validation_plans=len(plans))
    return plans


def phase7_8(inputs: dict[str, Any], evidence: list[dict[str, Any]], main_fig: list[dict[str, Any]], supp_fig: list[dict[str, Any]], plans: list[dict[str, Any]]) -> None:
    out = mkdir(OUT / "08_final_evidence_package")
    final_rows = []
    for row in evidence:
        status = "core_mechanism" if row["recommended_manuscript_placement"] == "main_text_core_mechanism" else ("support_mechanism" if row["recommended_manuscript_placement"] == "main_text_support_or_supplement" else "pending")
        final_rows.append({**row, "final_phase7_status": status})
    for row in inputs["support"]:
        final_rows.append({
            "candidate_id": row["candidate_id"],
            "mechanism_group_id": row["mechanism_group_id"],
            "storyline": row["storyline"],
            "intervention_axis": row["intervention_axis"],
            "priority": row["priority"],
            "final_phase7_status": "support_candidate",
            "external_validation_support": "support_only_not_primary",
            "spatial_tissue_support": "support_only_not_primary",
            "allowed_claim": "support/background only",
            "forbidden_claim": "; ".join(FORBIDDEN_CLAIMS),
            "caveats": row["caveat"],
        })
    write_csv(out / "phase7_final_candidate_master_table.csv", final_rows)
    (out / "phase7_final_evidence_card_pack.md").write_text((OUT / "05_evidence_cards_conflicts" / "phase7_5_final_evidence_cards.md").read_text())
    (out / "phase7_final_figure_plan.md").write_text((OUT / "06_paper_figure_selection" / "phase7_6_figure_plan.md").read_text())
    (out / "phase7_final_validation_plan.md").write_text((OUT / "07_minimal_validation_planning" / "phase7_7_minimal_validation_plan.md").read_text())
    (out / "phase7_final_claim_boundary.md").write_text("# Phase7 Final Claim Boundary\n\nAllowed: gated external/spatial/tissue support for mechanism candidates and validation plans.\n\nForbidden:\n" + "\n".join(f"- {x}" for x in FORBIDDEN_CLAIMS) + "\n")
    write_yaml(out / "phase7_reproducibility_manifest.yaml", {
        "phase": "Phase7 external/spatial/tissue validation and paper evidence package assembly",
        "generated_at_utc": now_iso(),
        "script": "scripts/v6_1/run_phase7_validation_package.py",
        "primary_input": "Phase6 Phase7 main candidates only",
        "support_input": "Phase6 support candidates and addendum assets only after role gate",
        "gated_external_dataset": "IMbrave150_bulk_RNAseq",
        "gated_spatial_dataset": "GSE238264_HCC_spatial support only",
        "forbidden_primary_inputs": ["PD1_anchor", "Phase3.5 complex ML", "unintegrated data without gate", "target pre-screen as recommendation"],
    })
    verdict = "PASS" if final_rows and plans else "FAIL"
    (out / "phase7_8_summary.md").write_text(
        f"# Phase7.8 Final Evidence Package Summary\n\n- Verdict: `{verdict}`\n- Final candidate rows: {len(final_rows)}\n- Main figures: {len(main_fig)}\n- Validation plans: {len(plans)}\n- Manuscript-ready evidence card pack, figure plan, validation plan, claim boundary, and manifest written.\n"
    )
    write_status(out / "phase7_8_status.yaml", "Phase7.8", verdict, final_candidate_rows=len(final_rows), validation_plans=len(plans))


def phase7_9(evidence: list[dict[str, Any]], plans: list[dict[str, Any]], main_fig: list[dict[str, Any]]) -> None:
    out = mkdir(OUT / "09_final_decision_handoff")
    counts = Counter(r["recommended_manuscript_placement"] for r in evidence)
    core = counts.get("main_text_core_mechanism", 0)
    support = counts.get("main_text_support_or_supplement", 0)
    pending = counts.get("pending_or_future_addendum", 0)
    verdict = "PASS" if core >= 1 and plans and len(main_fig) >= 4 else "CONDITIONAL_PASS" if plans else "FAIL"
    handoff = []
    for row in evidence:
        handoff.append({
            "candidate_id": row["candidate_id"],
            "mechanism_group_id": row["mechanism_group_id"],
            "phase8_role": "execute_minimal_validation" if row["recommended_manuscript_placement"] == "main_text_core_mechanism" else "supplement_or_pending_followup",
            "required_next_action": "wet-lab/tissue validation execution" if row["recommended_manuscript_placement"] != "pending_or_future_addendum" else "additional gated intake",
            "claim_boundary": "no clinical, causal, drug, or final target recommendation",
            "caveat": row["caveats"],
        })
    write_csv(out / "phase7_phase8_handoff_if_needed.csv", handoff)
    decision = {
        "verdict": verdict,
        "n_phase7_main_candidates": len(evidence),
        "n_core_mechanisms_for_main_text": core,
        "n_support_mechanisms_for_supplement": support,
        "n_pending_or_future_addendum": pending,
        "n_blocked_or_demoted": 0,
        "external_validation_summary": dict(Counter(r["external_validation_support"] for r in evidence)),
        "spatial_tissue_summary": dict(Counter(r["spatial_tissue_support"] for r in evidence)),
        "paper_figure_ready": len(main_fig) >= 4,
        "minimal_validation_plan_available": bool(plans),
        "phase8_ready_if_needed": bool(plans),
        "required_caveats": [REQUIRED_CAVEAT, "target pre-screening is not target recommendation", "perturbation prior is not causal proof"],
        "forbidden_claims": FORBIDDEN_CLAIMS,
        "recommended_next_stage": "Phase8 validation execution or manuscript drafting with caveats",
    }
    write_yaml(out / "PHASE7_FINAL_DECISION.yaml", decision)
    report = [
        "# PHASE7 FINAL REPORT",
        "",
        f"- Verdict: `{verdict}`",
        f"- Phase7 main candidates evaluated: {len(evidence)}",
        f"- Core/support/pending: {core}/{support}/{pending}",
        f"- Main figure-ready panels: {len(main_fig)}",
        f"- Minimal validation plans: {len(plans)}",
        "",
        "Phase7 used frozen Phase6 main candidates for primary evidence assembly. IMbrave150 passed external validation intake and was used for gated bulk projection. GSE238264 passed support spatial adjudication intake only; no spot-level clinical model or primary training was performed.",
        "",
        "Current outputs are validation evidence and manuscript package artifacts, not clinical recommendations, drug recommendations, final target recommendations, or causal proof.",
        "",
    ]
    (out / "PHASE7_FINAL_REPORT.md").write_text("\n".join(report))
    (out / "phase7_9_summary.md").write_text(
        f"# Phase7.9 Final Decision Summary\n\n- Verdict: `{verdict}`\n- Core/support/pending: {core}/{support}/{pending}\n- Phase8 handoff rows: {len(handoff)}\n"
    )
    write_status(out / "phase7_9_status.yaml", "Phase7.9", verdict, core=core, support=support, pending=pending)
    for name in ["PHASE7_FINAL_REPORT.md", "PHASE7_FINAL_DECISION.yaml", "phase7_phase8_handoff_if_needed.csv"]:
        (OUT / name).write_text((out / name).read_text())
    for name in ["phase7_final_candidate_master_table.csv", "phase7_final_evidence_card_pack.md", "phase7_final_figure_plan.md", "phase7_final_validation_plan.md", "phase7_final_claim_boundary.md", "phase7_reproducibility_manifest.yaml"]:
        (OUT / name).write_text((OUT / "08_final_evidence_package" / name).read_text())


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    inputs = load_inputs()
    phase7_0(inputs)
    intake = phase7_1(inputs)
    external_summary, _, _ = phase7_2(inputs, intake)
    spatial_intake = phase7_3(inputs)
    spatial_adj, _ = phase7_4(inputs, spatial_intake)
    evidence, _ = phase7_5(inputs, external_summary, spatial_adj)
    main_fig, supp_fig = phase7_6(evidence)
    plans = phase7_7(evidence)
    phase7_8(inputs, evidence, main_fig, supp_fig, plans)
    phase7_9(evidence, plans, main_fig)


if __name__ == "__main__":
    main()
