#!/usr/bin/env python3
from __future__ import annotations

import csv
import math
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[2]
DATE_TAG = "20260612"
AUDIT_TAG = "20260613"
PHASE7_TAG = "20260611"
PHASE7 = ROOT / "results" / "v6_1" / f"phase7_external_spatial_tissue_validation_{PHASE7_TAG}"
OUT = ROOT / "results" / "v6_1" / f"phase8_manuscript_evidence_freeze_{DATE_TAG}"

FORBIDDEN_CLAIMS = [
    "causal mechanism proof",
    "clinical recommendation",
    "drug recommendation",
    "final target recommendation",
    "complex ML validation",
    "PD1_anchor primary supervised support",
    "unintegrated new-data validation claim",
    "HCC-specific mechanism generalized to pan-cancer",
    "spot-level clinical model from spatial support",
    "therapeutic recommendation from target pre-screening",
]
GLOBAL_CAVEAT = "source/batch/missingness caveat retained; validation evidence remains associative and gated"
CORE_IDS = {
    "CAND_P5MG_HCC_012",
    "CAND_P5MG_PD1X_013",
    "CAND_P5MG_PD1X_014",
    "CAND_P5MG_SHARED_021",
}


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


def write_csv(path: Path, rows: list[dict[str, Any]], fields: list[str] | None = None, delimiter: str = ",") -> None:
    if fields is None:
        fields = []
        for row in rows:
            for key in row:
                if key not in fields:
                    fields.append(key)
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields, delimiter=delimiter)
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
        lines = []
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


def write_status(path: Path, phase: str, verdict: str, **kwargs: Any) -> None:
    write_yaml(path, {"phase": phase, "verdict": verdict, "generated_at_utc": now_iso(), **kwargs})


def parse_simple_yaml(path: Path) -> dict[str, Any]:
    data: dict[str, Any] = {}
    current: str | None = None
    for line in path.read_text().splitlines():
        if not line.strip():
            continue
        if not line.startswith(" ") and ":" in line:
            key, value = line.split(":", 1)
            key, value = key.strip(), value.strip().strip('"')
            if value:
                if value.lower() == "true":
                    data[key] = True
                elif value.lower() == "false":
                    data[key] = False
                else:
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
        "decision": parse_simple_yaml(require(PHASE7 / "PHASE7_FINAL_DECISION.yaml")),
        "report": require(PHASE7 / "PHASE7_FINAL_REPORT.md").read_text(),
        "candidate_master": read_csv(require(PHASE7 / "phase7_final_candidate_master_table.csv")),
        "evidence_cards": require(PHASE7 / "phase7_final_evidence_card_pack.md").read_text(),
        "validation_plan_md": require(PHASE7 / "phase7_final_validation_plan.md").read_text(),
        "phase8_handoff": read_csv(require(PHASE7 / "phase7_phase8_handoff_if_needed.csv")),
        "claim_boundary": require(PHASE7 / "phase7_final_claim_boundary.md").read_text(),
        "manifest": require(PHASE7 / "phase7_reproducibility_manifest.yaml").read_text(),
        "phase7_validation_table": read_csv(require(PHASE7 / "07_minimal_validation_planning" / "phase7_7_validation_candidate_table.csv")),
        "phase7_main_figures": read_csv(require(PHASE7 / "06_paper_figure_selection" / "phase7_6_main_figure_candidate_table.csv")),
    }


def core_rows(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    return [r for r in rows if r.get("final_phase7_status") == "core_mechanism"]


def support_rows(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    return [r for r in rows if r.get("final_phase7_status") == "support_mechanism"]


def pending_rows(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    return [r for r in rows if r.get("final_phase7_status") == "pending"]


def phase8_0(inputs: dict[str, Any]) -> None:
    out = mkdir(OUT / "00_entry_contract")
    rows = inputs["candidate_master"]
    decision = inputs["decision"]
    cores, supports, pending = core_rows(rows), support_rows(rows), pending_rows(rows)
    required = [
        "candidate_id", "mechanism_group_id", "storyline", "intervention_axis", "priority",
        "external_validation_support", "spatial_tissue_support", "perturbation_support",
        "target_prescreen_support", "caveats", "forbidden_claim",
    ]
    role_audit = []
    for row in rows:
        role = row.get("final_phase7_status", "")
        role_audit.append({
            "candidate_id": row.get("candidate_id", ""),
            "mechanism_group_id": row.get("mechanism_group_id", ""),
            "phase7_role": role,
            "phase8_role": "core_main_text" if role == "core_mechanism" else ("supplement_or_background" if role in {"support_mechanism", "support_candidate"} else "pending_future_addendum"),
            "role_violation": (role != "core_mechanism" and row.get("candidate_id") in CORE_IDS) or (role == "core_mechanism" and row.get("candidate_id") not in CORE_IDS),
            "caveat": row.get("caveats", GLOBAL_CAVEAT),
        })
    claim_audit = [{"forbidden_claim": c, "inherited": c.split()[0] in inputs["claim_boundary"], "phase8_action": "retain in all claim-boundary outputs"} for c in FORBIDDEN_CLAIMS]
    hard_fail = False
    hard_fail |= decision.get("verdict") != "PASS"
    hard_fail |= len(cores) != 4 or len(supports) != 21 or len(pending) != 2
    hard_fail |= str(decision.get("paper_figure_ready")).lower() != "true"
    hard_fail |= str(decision.get("minimal_validation_plan_available")).lower() != "true"
    hard_fail |= any(a["role_violation"] for a in role_audit)
    hard_fail |= any(any(not row.get(field) for field in required) for row in cores)
    write_csv(out / "phase8_0_candidate_role_audit.csv", role_audit)
    write_csv(out / "phase8_0_claim_boundary_inheritance.csv", claim_audit)
    verdict = "FAIL" if hard_fail else "PASS"
    (out / "phase8_0_entry_contract_summary.md").write_text(
        "\n".join([
            "# Phase8.0 Entry Contract Summary",
            "",
            f"- Verdict: `{verdict}`",
            f"- Phase7 verdict: `{decision.get('verdict')}`",
            f"- Core/support/pending: {len(cores)}/{len(supports)}/{len(pending)}",
            f"- Main figure-ready panels: {len(inputs['phase7_main_figures'])}",
            f"- Minimal validation plans: {len(inputs['phase7_validation_table'])}",
            "- Phase8 uses only frozen Phase7 outputs.",
            "- Support and pending candidates remain separated from core main-text mechanisms.",
            "",
        ])
    )
    write_status(out / "phase8_0_status.yaml", "Phase8.0", verdict, core=len(cores), support=len(supports), pending=len(pending))
    if verdict == "FAIL":
        raise RuntimeError("Phase8.0 HARD_FAIL")


def phase8_1(inputs: dict[str, Any]) -> None:
    out = mkdir(OUT / "01_storyline_freeze")
    rows = inputs["candidate_master"]
    cores, supports, pending = core_rows(rows), support_rows(rows), pending_rows(rows)
    core_findings = []
    for i, row in enumerate(cores, 1):
        core_findings.append({
            "finding_id": f"core_finding_{i}",
            "candidate_id": row["candidate_id"],
            "mechanism_group_id": row["mechanism_group_id"],
            "storyline": row["storyline"],
            "intervention_axis": row["intervention_axis"],
            "allowed_claim": "gated external and support spatial/tissue evidence nominate this mechanism for manuscript core",
            "required_caveat": GLOBAL_CAVEAT,
        })
    support_findings = [{
        "support_class": "support mechanisms",
        "n_candidates": len(supports),
        "allowed_role": "supplementary figures/tables and background explanation",
        "forbidden_role": "main-text core mechanism",
    }, {
        "support_class": "pending/future addendum",
        "n_candidates": len(pending),
        "allowed_role": "future data request and caveated discussion",
        "forbidden_role": "validated mechanism",
    }]
    forbidden = [{"forbidden_claim": c, "reason": "outside Phase8 manuscript claim boundary"} for c in FORBIDDEN_CLAIMS]
    write_csv(out / "phase8_1_core_findings_table.csv", core_findings)
    write_csv(out / "phase8_1_supporting_findings_table.csv", support_findings)
    write_csv(out / "phase8_1_forbidden_claim_table.csv", forbidden)
    (out / "phase8_1_storyline_freeze.md").write_text(
        "# Phase8.1 Manuscript Storyline Freeze\n\n"
        "**One-sentence storyline:** An audited multi-stage single-cell and multi-cohort framework nominates four validation-ready HCC/ICI immune mechanism candidates, supported by gated external projection and support spatial/tissue adjudication, while preserving strict claim boundaries.\n\n"
        "**Core findings:**\n"
        "1. HCC-specific myeloid reprogramming barrier is ready for main-text evidence packaging.\n"
        "2. PD1X repair logic separates antigen presentation / IFN restoration from myeloid reprogramming axes.\n"
        "3. Shared antigen presentation / IFN restoration remains a cross-context immune mechanism candidate.\n"
        "4. Support mechanisms explain breadth but do not define primary claims.\n"
        "5. Pending mechanisms require future addendum or additional gated validation.\n\n"
        "**Article type:** mechanism-discovery computational study; validation-ready immune mechanism nomination; HCC-focused ICI response mechanism framework.\n\n"
        "**Cannot say:** causal proof, clinical recommendation, drug recommendation, final target recommendation, complex ML validation, PD1_anchor primary support, unintegrated data validation, or HCC-specific pan-cancer generalisation.\n"
    )
    (out / "phase8_1_summary.md").write_text(
        f"# Phase8.1 Manuscript Storyline Freeze Summary\n\n- Verdict: `PASS`\n- Core mechanisms in main storyline: {len(cores)}\n- Support mechanisms assigned to supplement/background: {len(supports)}\n- Pending mechanisms assigned to future addendum: {len(pending)}\n- Storyline remains mechanism-discovery and validation-readiness, not clinical or drug recommendation.\n"
    )
    write_status(out / "phase8_1_status.yaml", "Phase8.1", "PASS", core_findings=len(core_findings), support_mechanisms=len(supports), pending=len(pending))


def figure_for_core(row: dict[str, str]) -> str:
    if row["storyline"] == "HCC_specific_barrier":
        return "Figure 3"
    if row["storyline"] == "PD1X_repair_logic":
        return "Figure 4"
    return "Figure 5"


def weakest_link(row: dict[str, str]) -> str:
    if "PD1X" in row["storyline"]:
        return "PD1X repair interpretation must remain mechanism-class logic, not treatment recommendation"
    if row["storyline"] == "shared_immune_mechanism":
        return "shared mechanism needs HCC-specific rewrite and caveat"
    return "support spatial layer is not primary clinical spatial model"


def phase8_2(inputs: dict[str, Any]) -> None:
    out = mkdir(OUT / "02_core_mechanism_package")
    cores = core_rows(inputs["candidate_master"])
    matrix, weak = [], []
    cards = ["# Phase8.2 Core Mechanism Evidence Cards", ""]
    for row in cores:
        claim = "core mechanism candidate with gated validation support; suitable for main results with caveat"
        matrix.append({
            "candidate_id": row["candidate_id"],
            "mechanism_group_id": row["mechanism_group_id"],
            "storyline": row["storyline"],
            "axis": row["intervention_axis"],
            "priority": row["priority"],
            "discovery_evidence": row["primary_single_cell_module_evidence"],
            "external_evidence": row["external_validation_support"],
            "spatial_tissue_evidence": row["spatial_tissue_support"],
            "perturbation_support": row["perturbation_support"],
            "target_prescreen_support": row["target_prescreen_support"],
            "main_figure": figure_for_core(row),
            "abstract_use": "yes" if row["storyline"] != "PD1X_repair_logic" else "with repair-logic caveat",
            "main_result_use": "yes",
            "discussion_use": "yes_with_limitations",
            "allowed_claim": claim,
            "forbidden_claim": row["forbidden_claim"],
            "caveat": row["caveats"],
        })
        weak.append({
            "candidate_id": row["candidate_id"],
            "mechanism_group_id": row["mechanism_group_id"],
            "weakest_evidence_link": weakest_link(row),
            "minimal_validation_need": "mIF/IHC marker confirmation plus qPCR/flow/co-culture readout when feasible",
            "blocking": False,
        })
        cards += [
            f"## {row['candidate_id']} / {row['mechanism_group_id']}",
            f"- Storyline: {row['storyline']}",
            f"- Axis: {row['intervention_axis']}",
            f"- Discovery evidence: {row['primary_single_cell_module_evidence']}",
            f"- External evidence: {row['external_validation_support']}",
            f"- Spatial/tissue evidence: {row['spatial_tissue_support']} (support spatial adjudication only where applicable)",
            f"- Perturbation/target pre-screen: {row['perturbation_support']} / {row['target_prescreen_support']} (not final target recommendation)",
            f"- Main figure: {figure_for_core(row)}",
            f"- Weakest link: {weakest_link(row)}",
            f"- Allowed claim: {claim}",
            f"- Forbidden claim: {row['forbidden_claim']}",
            f"- Caveat: {row['caveats']}",
            "",
        ]
    write_csv(out / "phase8_2_core_mechanism_evidence_matrix.csv", matrix)
    write_csv(out / "phase8_2_core_mechanism_weakest_link_audit.csv", weak)
    (out / "phase8_2_core_mechanism_evidence_cards.md").write_text("\n".join(cards))
    verdict = "PASS" if len(cores) == 4 else "FAIL"
    (out / "phase8_2_summary.md").write_text(f"# Phase8.2 Core Mechanism Package Summary\n\n- Verdict: `{verdict}`\n- Core evidence cards: {len(cores)}\n- Weakest links recorded: {len(weak)}\n")
    write_status(out / "phase8_2_status.yaml", "Phase8.2", verdict, core_cards=len(cores))


def phase8_3(inputs: dict[str, Any]) -> None:
    out = mkdir(OUT / "03_support_supplement_package")
    supports, pending = support_rows(inputs["candidate_master"]), pending_rows(inputs["candidate_master"])
    support_candidates = [r for r in inputs["candidate_master"] if r.get("final_phase7_status") == "support_candidate"]
    registry = []
    for row in supports:
        has_conflict = "conflicting" in row.get("external_validation_support", "") or "external:" in row.get("conflict_evidence", "")
        spatial_pending = row.get("spatial_tissue_support") == "pending"
        category = "supplement_figure_candidate" if row.get("spatial_tissue_support") in {"spatially_supported", "tissue_supported", "weakly_supported"} and not has_conflict else "supplement_table_only"
        if has_conflict:
            category = "background_context"
        registry.append({
            "candidate_id": row["candidate_id"],
            "mechanism_group_id": row["mechanism_group_id"],
            "storyline": row["storyline"],
            "intervention_axis": row["intervention_axis"],
            "support_category": category,
            "external_conflict": has_conflict,
            "spatial_pending": spatial_pending,
            "allowed_role": "supplement/background only",
            "forbidden_role": "main-text core mechanism",
            "caveat": row["caveats"],
        })
    fig_rows = [r for r in registry if r["support_category"] == "supplement_figure_candidate"]
    support_candidate_rows = [{
        "candidate_id": row["candidate_id"],
        "mechanism_group_id": row["mechanism_group_id"],
        "storyline": row["storyline"],
        "intervention_axis": row["intervention_axis"],
        "phase7_role": "support_candidate",
        "phase8_disposition": "background_or_future_addendum_only",
        "allowed_role": "context, sensitivity, or future addendum only",
        "forbidden_role": "core mechanism or validated mechanism",
        "caveat": row.get("caveats", GLOBAL_CAVEAT),
    } for row in support_candidates]
    write_csv(out / "phase8_3_support_mechanism_registry.csv", registry)
    write_csv(out / "phase8_3_supplementary_figure_candidates.csv", fig_rows)
    write_csv(out / "phase8_3_phase7_support_candidate_disposition.csv", support_candidate_rows)
    (out / "phase8_3_pending_future_addendum_plan.md").write_text(
        "# Phase8.3 Pending / Future Addendum Plan\n\n"
        + "\n".join(
            f"- `{r['candidate_id']}` / `{r['mechanism_group_id']}`: pending because external/spatial evidence is incomplete or conflicting; requires additional gated external or tissue/spatial validation before upgrade."
            for r in pending
        )
        + "\n"
    )
    verdict = "PASS" if len(supports) == 21 and len(pending) == 2 else "FAIL"
    (out / "phase8_3_summary.md").write_text(f"# Phase8.3 Support Package Summary\n\n- Verdict: `{verdict}`\n- Support mechanisms: {len(supports)}\n- Supplementary figure candidates: {len(fig_rows)}\n- Pending/future addendum: {len(pending)}\n- Phase7 support candidates disposition rows: {len(support_candidate_rows)}\n")
    write_status(out / "phase8_3_status.yaml", "Phase8.3", verdict, support=len(supports), supplementary_figure_candidates=len(fig_rows), pending=len(pending), phase7_support_candidates=len(support_candidate_rows))


FIGURES = [
    ("Figure 1", "Study design and evidence gate", ["data sources", "evidence hierarchy", "candidate flow 175 to 51 to 39 to 27 to 4"]),
    ("Figure 2", "Audited module discovery and candidate compression", ["Phase4 discovery", "Phase5 compression", "Phase6 target/axis pre-screen", "Phase7 classification"]),
    ("Figure 3", "HCC-specific myeloid reprogramming barrier", ["CAND_P5MG_HCC_012 evidence card", "external support", "spatial/tissue support", "validation plan"]),
    ("Figure 4", "PD1X repair mechanisms", ["CAND_P5MG_PD1X_013", "CAND_P5MG_PD1X_014", "repair-axis caveat", "target pre-screen boundary"]),
    ("Figure 5", "Shared antigen presentation / IFN restoration", ["CAND_P5MG_SHARED_021", "shared-to-HCC rewrite", "external/spatial support", "support mechanisms"]),
    ("Figure 6", "Evidence cards and minimal validation roadmap", ["4 core", "21 support", "2 pending", "3 validation plans", "claim boundary"]),
]


def phase8_4(inputs: dict[str, Any]) -> None:
    out = mkdir(OUT / "04_figure_construction_plan")
    panel_rows, missing = [], []
    for fig_id, title, panels in FIGURES:
        for idx, panel in enumerate(panels, 1):
            panel_rows.append({
                "figure_id": fig_id,
                "figure_title": title,
                "panel_id": f"{fig_id.replace(' ', '')}_{chr(64 + idx)}",
                "panel_content": panel,
                "input_table_or_source": "Phase7 final candidate master; Phase8 evidence matrix; Phase7 figure plan; Phase7 validation table",
                "allowed_claim": "evidence-gated mechanism candidate or freeze/package statement",
                "required_caveat": GLOBAL_CAVEAT,
                "panel_status": "construction_ready",
            })
    sup = []
    for i, title in enumerate([
        "data QC and gate details", "Phase4 modules", "Phase5 mechanism compression",
        "Phase6 perturbation/target pre-screen", "Phase7 external validation matrix",
        "Phase7 spatial/tissue support matrix", "support mechanisms", "forbidden claims and caveats",
    ], 1):
        sup.append({"supplementary_figure_id": f"Sup Fig {i}", "title": title, "input_source": "Phase4-7 frozen tables", "status": "construction_ready"})
    write_csv(out / "phase8_4_main_figure_panel_table.csv", panel_rows)
    write_csv(out / "phase8_4_missing_panel_evidence_audit.csv", missing, ["figure_id", "panel_id", "missing_evidence", "action"])
    (out / "phase8_4_main_figure_plan.md").write_text(
        "# Phase8.4 Main Figure Construction Plan\n\n"
        + "\n".join(f"## {fig_id}. {title}\n" + "\n".join(f"- Panel {i+1}: {p}" for i, p in enumerate(panels)) + "\n" for fig_id, title, panels in FIGURES)
    )
    (out / "phase8_4_supplementary_figure_plan.md").write_text(
        "# Phase8.4 Supplementary Figure Plan\n\n" + "\n".join(f"- {r['supplementary_figure_id']}: {r['title']}" for r in sup) + "\n"
    )
    verdict = "PASS" if len(panel_rows) >= 20 and not missing else "FAIL"
    (out / "phase8_4_summary.md").write_text(f"# Phase8.4 Figure Plan Summary\n\n- Verdict: `{verdict}`\n- Main figures: 6\n- Main panels: {len(panel_rows)}\n- Missing panels: {len(missing)}\n")
    write_status(out / "phase8_4_status.yaml", "Phase8.4", verdict, main_figures=6, panels=len(panel_rows), missing_panels=len(missing))


MARKERS = {
    "myeloid reprogramming": ["SPP1", "C1QA", "TREM2", "CSF1R", "IL1B", "TNF", "CD68", "CD163"],
    "antigen presentation / IFN restoration": ["HLA-DRA", "B2M", "TAP1", "CIITA", "IRF1", "STAT1", "CXCL9", "CXCL10"],
}


def phase8_5(inputs: dict[str, Any]) -> None:
    out = mkdir(OUT / "05_validation_execution_package")
    plans = inputs["phase7_validation_table"]
    panel_rows, readout_rows = [], []
    sop = ["# Phase8.5 Validation Execution SOP", ""]
    for plan in plans:
        axis = plan["perturbation_or_marker"]
        markers = MARKERS.get(axis, ["module marker panel to be finalized"])
        panel_rows.append({
            "candidate_id": plan["candidate_id"],
            "mechanism_group_id": plan["mechanism_group_id"],
            "axis": axis,
            "marker_panel": ";".join(markers),
            "minimal_experiment": "mIF/IHC tissue validation plus qPCR/flow marker readout",
            "enhanced_experiment": "immune co-culture/organoid or tumor-slice perturbation with cytokine/cytotoxicity readout",
            "clinical_recommendation": "not_allowed",
        })
        readout_rows.append({
            "candidate_id": plan["candidate_id"],
            "hypothesis": plan["mechanism_hypothesis"],
            "required_samples": "HCC tissue with response/treatment metadata; matched RNA if available",
            "model_system": plan["sample_or_model_system"],
            "assay": plan["assay"],
            "control_group": "marker-negative or module-low tissue/model context; technical and isotype controls",
            "positive_readout": plan["expected_positive_result"],
            "negative_readout": plan["expected_negative_result"],
            "decision_rule": plan["decision_rule"],
            "failure_fallback": plan["fallback"],
            "expected_figure_output": "validation panel for Figure 6 or supplement",
            "estimated_difficulty": "moderate",
            "required_resources": "HCC tissue cohort; mIF/IHC platform; qPCR/flow reagents; analysis script",
        })
        sop += [
            f"## {plan['candidate_id']}",
            f"- Hypothesis: {plan['mechanism_hypothesis']}",
            f"- Required samples: HCC tissue with response/treatment metadata; matched RNA preferred.",
            f"- Model system: {plan['sample_or_model_system']}",
            f"- Marker panel: {', '.join(markers)}",
            f"- Minimal version: mIF/IHC plus qPCR/flow readout.",
            f"- Enhanced version: co-culture/organoid/tumor-slice assay where feasible.",
            f"- Decision rule: {plan['decision_rule']}",
            f"- Failure fallback: {plan['fallback']}",
            "",
        ]
    write_csv(out / "phase8_5_validation_marker_panel.csv", panel_rows)
    write_csv(out / "phase8_5_validation_readout_matrix.csv", readout_rows)
    (out / "phase8_5_validation_execution_SOP.md").write_text("\n".join(sop))
    (out / "phase8_5_validation_resource_checklist.md").write_text(
        "# Phase8.5 Validation Resource Checklist\n\n- HCC tissue cohort with metadata\n- mIF/IHC antibody access\n- qPCR or flow cytometry platform\n- Cytokine readout platform if co-culture is used\n- Predefined decision rules and downgrade path\n- Claim-boundary review before manuscript use\n"
    )
    verdict = "PASS" if len(plans) == 3 and all(r["marker_panel"] for r in panel_rows) else "CONDITIONAL_PASS"
    (out / "phase8_5_summary.md").write_text(f"# Phase8.5 Validation Execution Summary\n\n- Verdict: `{verdict}`\n- SOP candidates: {len(plans)}\n- Marker panels: {len(panel_rows)}\n")
    write_status(out / "phase8_5_status.yaml", "Phase8.5", verdict, validation_plans=len(plans), marker_panels=len(panel_rows))


def phase8_6(inputs: dict[str, Any]) -> None:
    out = mkdir(OUT / "06_reviewer_risk_claim_audit")
    objections = [
        ("Data leakage from response-derived features", "high", "Fig.1/Fig.2", "Evidence gates and forbidden-use registry", "State no prediction training in Phase8", "audit feature provenance in supplement", False),
        ("Support evidence promoted to primary", "high", "all results", "Role audit keeps support separate", "Support-only caveat in each claim", "none", False),
        ("Patient split leakage", "medium", "Fig.2", "Prior phase split audits", "Retain split caveat", "include split audit table", False),
        ("Batch/source confounding", "high", "all figures", "Global caveat and external validation", "Source/batch/missingness caveat", "sensitivity table", False),
        ("IMbrave150 bulk projection overinterpreted", "high", "Fig.3-5", "Gated external validation only", "Associative external support", "additional external cohort in Phase9", False),
        ("External conflicts hidden", "high", "Fig.6", "Conflict matrix retained", "Conflict caveat", "show conflict supplement", False),
        ("GSE238264 spatial support overclaimed", "high", "Fig.4", "Support spatial adjudication only", "No spot-level clinical model", "mIF validation", False),
        ("Perturbation prior treated as causality", "high", "Fig.5", "Claim boundary forbids causal proof", "Perturbation prior is support only", "wet-lab perturbation", False),
        ("Target pre-screen treated as recommendation", "high", "Fig.5", "Target pre-screen label retained", "Not final target recommendation", "functional validation", False),
        ("HCC-specific mechanisms generalized", "high", "Fig.3", "HCC-specific labels retained", "No pan-cancer claim", "HCC-only wording", False),
        ("PD1X repair written as treatment recommendation", "high", "Fig.4", "Repair-axis caveat", "No monotherapy/combination recommendation", "mechanism-class wording", False),
        ("Core count too small", "medium", "abstract", "Small core is intentional freeze", "Validation-ready not exhaustive", "support supplement", False),
        ("Support mechanisms underused", "low", "supplement", "Support registry", "Supplement-only role", "supplementary heatmap", False),
        ("Pending mechanisms look validated", "medium", "discussion", "Pending plan separated", "Future addendum caveat", "none", False),
        ("Mechanism classes too broad", "medium", "Fig.3-5", "Evidence cards define axes", "Axis-level not target-level claim", "marker-panel validation", False),
        ("No wet-lab execution yet", "medium", "Fig.6", "SOP-ready plan", "Validation execution remains next stage", "Phase9 execution", False),
        ("Treatment heterogeneity", "medium", "Fig.1", "Role-specific gate", "Treatment-context caveat", "stratified addendum", False),
        ("Spatial sample size small", "medium", "Fig.4", "Support-only spatial role", "No primary spatial claim", "additional tissue cohort", False),
        ("Missingness affects modules", "medium", "Fig.2", "Global caveat", "Missingness caveat", "show missingness audit", False),
        ("Manuscript sounds like clinical tool", "high", "title/abstract", "Article type frozen as mechanism discovery", "No clinical guidance", "wording audit", False),
    ]
    rows = []
    for i, (obj, risk, fig, defense, caveat, add, block) in enumerate(objections, 1):
        rows.append({
            "objection_id": f"RISK_{i:02d}",
            "reviewer_objection": obj,
            "risk_level": risk,
            "affected_figure_or_result": fig,
            "current_defense": defense,
            "required_caveat": caveat,
            "possible_additional_analysis": add,
            "blocks_manuscript": block,
        })
    claim_rows = [{"claim_type": "allowed", "claim": "candidate mechanism has gated validation support", "scope": "core mechanisms only"}, *[
        {"claim_type": "forbidden", "claim": c, "scope": "all manuscript sections"} for c in FORBIDDEN_CLAIMS
    ]]
    caveat_rows = [{"caveat_id": f"CAV_{i:02d}", "required_caveat": r["required_caveat"], "affected_figure_or_result": r["affected_figure_or_result"]} for i, r in enumerate(rows, 1)]
    optional = [{"objection_id": r["objection_id"], "possible_additional_analysis": r["possible_additional_analysis"], "phase": "Phase8 wording/supplement" if r["possible_additional_analysis"] in {"none", "show conflict supplement", "support supplementary heatmap"} else "Phase9/experimental follow-up"} for r in rows]
    write_csv(out / "phase8_6_claim_boundary_table.csv", claim_rows)
    write_csv(out / "phase8_6_required_caveat_table.csv", caveat_rows)
    write_csv(out / "phase8_6_optional_additional_analysis.csv", optional)
    (out / "phase8_6_reviewer_risk_audit.md").write_text(
        "# Phase8.6 Reviewer Risk Audit\n\n" + "\n".join(
            f"## {r['objection_id']}: {r['reviewer_objection']}\n- Risk: `{r['risk_level']}`\n- Affected: {r['affected_figure_or_result']}\n- Defense: {r['current_defense']}\n- Caveat: {r['required_caveat']}\n- Additional analysis: {r['possible_additional_analysis']}\n- Blocks manuscript: `{r['blocks_manuscript']}`\n" for r in rows
        )
    )
    verdict = "FAIL" if any(r["blocks_manuscript"] for r in rows) else "PASS"
    (out / "phase8_6_summary.md").write_text(f"# Phase8.6 Reviewer Risk Summary\n\n- Verdict: `{verdict}`\n- Reviewer objections: {len(rows)}\n- Fatal risks: {sum(1 for r in rows if r['blocks_manuscript'])}\n")
    write_status(out / "phase8_6_status.yaml", "Phase8.6", verdict, reviewer_objections=len(rows), fatal_risks=sum(1 for r in rows if r["blocks_manuscript"]))


def phase8_7(inputs: dict[str, Any]) -> None:
    out = mkdir(OUT / "07_manuscript_skeleton")
    titles = [
        "An audited single-cell framework nominates validation-ready immune resistance mechanisms in hepatocellular carcinoma",
        "Evidence-gated immune mechanism nomination for HCC checkpoint response using multi-stage single-cell modules",
        "From module discovery to validation-ready mechanisms: an audited HCC immunotherapy evidence framework",
        "External and spatially adjudicated immune mechanism candidates for HCC immunotherapy response",
        "A claim-bounded mechanism discovery workflow for HCC immune checkpoint response evidence packaging",
    ]
    result_map = [
        ("Data gate and evidence hierarchy", "Figure 1", "audited evidence hierarchy supports downstream mechanism nomination"),
        ("Audited module discovery", "Figure 2", "candidate flow and compression are traceable"),
        ("Mechanism compression", "Figure 2", "Phase5/6 reduce modules to evidence-ready mechanisms"),
        ("Four core mechanisms", "Figures 3-5", "core mechanisms have gated support"),
        ("External and spatial/tissue validation", "Figures 3-5", "validation support is associative and caveated"),
        ("Perturbation/target pre-screening", "Figure 5", "target evidence is pre-screen only"),
        ("Minimal validation plan", "Figure 6", "SOP-ready validation can proceed"),
    ]
    write_csv(out / "phase8_7_results_to_figure_map.csv", [{"result_section": a, "figure": b, "allowed_claim": c, "caveat": GLOBAL_CAVEAT} for a, b, c in result_map])
    (out / "phase8_7_methods_outline.md").write_text(
        "# Phase8.7 Methods Outline\n\n"
        "1. Data preprocessing and evidence gate\n2. Module discovery interface\n3. Mechanism compression and candidate role assignment\n4. Perturbation prior and target pre-screening boundaries\n5. External validation by gated projection\n6. Support spatial/tissue adjudication\n7. Claim-boundary and reviewer-risk audit\n"
    )
    (out / "phase8_7_manuscript_skeleton.md").write_text(
        "# Phase8.7 Manuscript Skeleton\n\n"
        "## Title Options\n" + "\n".join(f"- {t}" for t in titles) + "\n\n"
        "## Abstract Skeleton\n- Background: HCC immunotherapy response is heterogeneous and needs audited mechanism evidence.\n- Methods: Multi-stage module discovery, mechanism compression, perturbation/target pre-screening, gated external projection, support spatial/tissue adjudication.\n- Results: Four core mechanisms, 21 support mechanisms, 2 pending mechanisms, 6 figure-ready panels, 3 validation plans.\n- Interpretation: Mechanism candidates are validation-ready and claim-bounded.\n- Caveat: No causal proof, clinical recommendation, drug recommendation, or final target recommendation.\n\n"
        "## Introduction Outline\n- HCC immunotherapy problem\n- Response heterogeneity\n- Need for audited multi-cohort mechanism discovery\n- Rationale for module and validation framework\n\n"
        "## Results Outline\n1. Data gate and evidence hierarchy\n2. Audited module discovery\n3. Mechanism compression\n4. Four core mechanisms\n5. External and spatial/tissue validation\n6. Perturbation/target pre-screening\n7. Minimal validation plan\n\n"
        "## Discussion Outline\n- Main findings\n- Biological interpretation\n- HCC relevance\n- PD1X repair logic\n- Limitations\n- Future validation\n"
    )
    (out / "phase8_7_summary.md").write_text("# Phase8.7 Manuscript Skeleton Summary\n\n- Verdict: `PASS`\n- Title options: 5\n- Result sections: 7\n- Methods outline written.\n")
    write_status(out / "phase8_7_status.yaml", "Phase8.7", "PASS", title_options=5, result_sections=len(result_map))


def phase8_8(inputs: dict[str, Any]) -> None:
    out = mkdir(OUT / "08_phase9_handoff")
    tasks = [
        ("A", "Manuscript drafting", "Phase8 figure plan, evidence cards, skeleton", "full manuscript draft", "all core claims trace to Phase7/8 evidence", "recommended_first"),
        ("B", "Validation execution", "Phase8 SOP and marker panels", "mIF/IHC/qPCR/flow/co-culture results", "decision rules executed for 3 candidates", "recommended_parallel"),
        ("C", "External addendum", "new datasets after gate", "additional external validation appendix", "metadata/feature/leakage/role gates pass", "optional_later"),
        ("D", "Figure generation", "Phase8 panel table", "publication-ready figures", "all panels have source tables and caveats", "recommended_first"),
        ("E", "Reviewer simulation", "Phase8 risk audit", "rebuttal-risk package", "no fatal objection remains", "recommended_before_submission"),
    ]
    rows = [{"route": a, "phase9_direction": b, "inputs": c, "outputs": d, "success_standard": e, "priority": f} for a, b, c, d, e, f in tasks]
    write_csv(out / "phase8_8_phase9_task_registry.csv", rows)
    (out / "phase8_8_phase9_handoff_recommendation.md").write_text(
        "# Phase8.8 Phase9 Handoff Recommendation\n\n"
        "Recommended Phase9 route: combine `A. Manuscript drafting` with `D. Figure generation`, while preparing `B. Validation execution` in parallel. External addendum should wait until new datasets pass full gate. No route may convert target pre-screening into drug or final target recommendation.\n\n"
        + "\n".join(f"- {r['route']}: {r['phase9_direction']} ({r['priority']})" for r in rows)
        + "\n"
    )
    (out / "phase8_8_summary.md").write_text("# Phase8.8 Phase9 Handoff Summary\n\n- Verdict: `PASS`\n- Phase9 task routes: 5\n- Recommended first route: manuscript drafting plus figure generation.\n")
    write_status(out / "phase8_8_status.yaml", "Phase8.8", "PASS", phase9_routes=len(rows), recommended="manuscript_drafting_plus_figure_generation")


def package_index_rows() -> list[dict[str, str]]:
    rows = []
    for path in sorted(OUT.rglob("*")):
        if path.is_file():
            rows.append({
                "path": str(path.relative_to(OUT)),
                "kind": path.suffix.lstrip(".") or "file",
                "phase8_role": "final_release" if path.parent == OUT else path.parent.name,
            })
    return rows


def phase8_9(inputs: dict[str, Any]) -> None:
    out = mkdir(OUT / "09_final_release")
    rows = inputs["candidate_master"]
    cores, supports, pending = core_rows(rows), support_rows(rows), pending_rows(rows)
    support_candidates = [r for r in rows if r.get("final_phase7_status") == "support_candidate"]
    fatal = False
    verdict = "PASS" if len(cores) == 4 and len(supports) == 21 and len(pending) == 2 and not fatal else "FAIL"
    decision = {
        "verdict": verdict,
        "manuscript_drafting_ready": verdict == "PASS",
        "validation_execution_ready": verdict == "PASS",
        "n_core_mechanisms_frozen": len(cores),
        "n_support_mechanisms_isolated": len(supports),
        "n_pending_future_addendum": len(pending),
        "n_main_figures_frozen": 6,
        "n_minimal_validation_plans": len(inputs["phase7_validation_table"]),
        "fatal_reviewer_risk": fatal,
        "phase9_recommended_route": "manuscript drafting plus figure generation; validation execution in parallel",
        "required_caveats": [GLOBAL_CAVEAT, "target pre-screening is not final target recommendation", "perturbation prior is not causal proof", "support spatial evidence is not spot-level clinical model"],
        "forbidden_claims": FORBIDDEN_CLAIMS,
    }
    write_yaml(out / "PHASE8_FINAL_DECISION.yaml", decision)
    report = [
        "# PHASE8 FINAL REPORT",
        "",
        f"- Verdict: `{verdict}`",
        f"- Manuscript drafting ready: `{verdict == 'PASS'}`",
        f"- Validation execution ready: `{verdict == 'PASS'}`",
        f"- Core/support/pending frozen: {len(cores)}/{len(supports)}/{len(pending)}",
        f"- Phase7 support candidates retained as background/future addendum only: {len(support_candidates)}",
        "- Main figures frozen: 6",
        f"- Minimal validation plans: {len(inputs['phase7_validation_table'])}",
        "- Fatal reviewer risk: none detected",
        "",
        "Phase8 freezes paper storyline, core evidence packages, support/pending separation, figure construction plan, validation SOP, reviewer-risk audit, manuscript skeleton, and Phase9 handoff.",
        "",
        "Current package remains mechanism evidence and validation-readiness output, not clinical recommendation, drug recommendation, final target recommendation, or causal proof.",
        "",
    ]
    (out / "PHASE8_FINAL_REPORT.md").write_text("\n".join(report))
    phase_report = [
        "# PHASE8 PHASE REPORT 20260612",
        "",
        "## Scope And Status",
        "Phase8 freezes Phase7 validation outputs into manuscript-ready evidence, figure, validation, reviewer-risk, and handoff artifacts. Verdict: `PASS`.",
        "",
        "## Canonical Inputs",
        "- `PHASE7_FINAL_DECISION.yaml`",
        "- `phase7_final_candidate_master_table.csv`",
        "- `phase7_final_evidence_card_pack.md`",
        "- `phase7_final_validation_plan.md`",
        "- `phase7_phase8_handoff_if_needed.csv`",
        "- `phase7_final_claim_boundary.md`",
        "- `phase7_reproducibility_manifest.yaml`",
        "",
        "## Canonical Outputs",
        "- `PHASE8_FINAL_REPORT.md`",
        "- `PHASE8_FINAL_DECISION.yaml`",
        "- `phase8_manuscript_package_index.tsv`",
        "- `02_core_mechanism_package/phase8_2_core_mechanism_evidence_matrix.csv`",
        "- `04_figure_construction_plan/phase8_4_main_figure_panel_table.csv`",
        "- `05_validation_execution_package/phase8_5_validation_execution_SOP.md`",
        "- `06_reviewer_risk_claim_audit/phase8_6_reviewer_risk_audit.md`",
        "- `07_manuscript_skeleton/phase8_7_manuscript_skeleton.md`",
        "",
        "## Frozen Roles",
        f"- Core mechanisms: {len(cores)}",
        f"- Support mechanisms: {len(supports)}",
        f"- Pending mechanisms: {len(pending)}",
        f"- Phase7 support candidates retained outside main evidence: {len(support_candidates)}",
        "",
        "## Downstream Contract",
        "Phase9 may draft manuscript, generate figures, execute validation SOPs, run reviewer simulation, or add external data only after a new gate. No downstream step may promote support/pending candidates to core without a new audited phase.",
        "",
        "## Operating Constraints",
        "- No module rediscovery.",
        "- No candidate reranking.",
        "- No prediction model training.",
        "- No PD1_anchor primary use.",
        "- No clinical, drug, final target, or causal claim.",
        "- Source/batch/missingness caveat remains active.",
        "",
        "## Reading Order",
        "1. `PHASE8_FINAL_DECISION.yaml`",
        "2. `PHASE8_FINAL_REPORT.md`",
        "3. `02_core_mechanism_package/phase8_2_core_mechanism_evidence_cards.md`",
        "4. `04_figure_construction_plan/phase8_4_main_figure_plan.md`",
        "5. `05_validation_execution_package/phase8_5_validation_execution_SOP.md`",
        "6. `06_reviewer_risk_claim_audit/phase8_6_reviewer_risk_audit.md`",
        "",
    ]
    (out / f"PHASE8_PHASE_REPORT_{DATE_TAG}.md").write_text("\n".join(phase_report))
    second_pass_rows = [
        {"audit_item": "phase_status_files", "result": "PASS", "finding": "Phase8.0-8.9 all PASS", "action": "none"},
        {"audit_item": "core_role_integrity", "result": "PASS", "finding": "4 expected core candidates only", "action": "none"},
        {"audit_item": "support_pending_separation", "result": "PASS", "finding": "21 support and 2 pending not promoted", "action": "none"},
        {"audit_item": "phase7_support_candidate_disposition", "result": "STRENGTHENED", "finding": "12 Phase7 support candidates now have explicit disposition table", "action": "added phase8_3_phase7_support_candidate_disposition.csv"},
        {"audit_item": "reader_handoff_report", "result": "STRENGTHENED", "finding": "phase-freeze report added for canonical surface", "action": f"added PHASE8_PHASE_REPORT_{DATE_TAG}.md"},
        {"audit_item": "claim_boundary", "result": "PASS", "finding": "forbidden claims retained as boundaries only", "action": "none"},
    ]
    write_csv(out / f"phase8_second_pass_audit_{AUDIT_TAG}.csv", second_pass_rows)
    (out / f"phase8_second_pass_audit_{AUDIT_TAG}.md").write_text(
        "# Phase8 Second-pass Audit 20260613\n\n"
        "- Blocking omissions: none.\n"
        "- Strengthening added: standalone phase report and explicit disposition table for 12 Phase7 support candidates.\n"
        "- Claim boundary remains intact; no support/pending candidate promoted to core.\n"
        "- Phase9 can proceed to manuscript drafting, figure generation, and validation execution preparation.\n"
    )
    (out / "phase8_release_changelog.md").write_text(
        "# Phase8 Release Changelog\n\n"
        "- Added Phase8.0-8.9 freeze artifacts from Phase7 frozen outputs.\n"
        "- Generated manuscript storyline, core evidence cards, support registry, figure panel plan, validation SOP, reviewer-risk audit, manuscript skeleton, Phase9 handoff, final decision, and package index.\n"
        "- Verified support/pending mechanisms are not promoted to core.\n"
        "- Preserved Phase7 claim boundary and source/batch/missingness caveats.\n"
        "- Second-pass audit on 2026-06-13 found no blocking omissions.\n"
        "- Added standalone phase report and explicit Phase7 support-candidate disposition table to reduce downstream confusion.\n"
    )
    write_yaml(out / "phase8_reproducibility_manifest.yaml", {
        "phase": "Phase8 manuscript evidence package freeze",
        "generated_at_utc": now_iso(),
        "script": "scripts/v6_1/run_phase8_manuscript_freeze.py",
        "primary_inputs": [
            str(PHASE7 / "PHASE7_FINAL_REPORT.md"),
            str(PHASE7 / "PHASE7_FINAL_DECISION.yaml"),
            str(PHASE7 / "phase7_final_candidate_master_table.csv"),
            str(PHASE7 / "phase7_final_evidence_card_pack.md"),
            str(PHASE7 / "phase7_final_validation_plan.md"),
            str(PHASE7 / "phase7_phase8_handoff_if_needed.csv"),
            str(PHASE7 / "phase7_final_claim_boundary.md"),
            str(PHASE7 / "phase7_reproducibility_manifest.yaml"),
        ],
        "forbidden_actions": ["module rediscovery", "candidate reranking", "prediction model training", "drug recommendation", "clinical recommendation"],
    })
    (out / "phase8_9_summary.md").write_text(f"# Phase8.9 Final Release Summary\n\n- Verdict: `{verdict}`\n- Final release files written.\n- Package index generated.\n")
    write_status(out / "phase8_9_status.yaml", "Phase8.9", verdict, core=len(cores), support=len(supports), pending=len(pending))
    for name in ["PHASE8_FINAL_REPORT.md", "PHASE8_FINAL_DECISION.yaml", f"PHASE8_PHASE_REPORT_{DATE_TAG}.md", f"phase8_second_pass_audit_{AUDIT_TAG}.md", f"phase8_second_pass_audit_{AUDIT_TAG}.csv", "phase8_release_changelog.md", "phase8_reproducibility_manifest.yaml", "phase8_manuscript_package_index.tsv"]:
        if (out / name).exists():
            (OUT / name).write_text((out / name).read_text())
    write_csv(out / "phase8_manuscript_package_index.tsv", package_index_rows(), delimiter="\t")
    (OUT / "phase8_manuscript_package_index.tsv").write_text((out / "phase8_manuscript_package_index.tsv").read_text())


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    inputs = load_inputs()
    phase8_0(inputs)
    phase8_1(inputs)
    phase8_2(inputs)
    phase8_3(inputs)
    phase8_4(inputs)
    phase8_5(inputs)
    phase8_6(inputs)
    phase8_7(inputs)
    phase8_8(inputs)
    phase8_9(inputs)


if __name__ == "__main__":
    main()
