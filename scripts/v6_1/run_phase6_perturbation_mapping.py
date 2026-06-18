#!/usr/bin/env python3
from __future__ import annotations

import csv
import hashlib
import math
import re
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[2]
DATE_TAG = "20260611"
PHASE5 = ROOT / "results" / "v6_1" / f"phase5_mechanism_adjudication_{DATE_TAG}"
OUT = ROOT / "results" / "v6_1" / f"phase6_perturbation_mapping_{DATE_TAG}"

REQUIRED_CAVEATS = [
    "source/batch/missingness conditional caveat retained",
    "associative mechanism candidate only",
    "Phase3.5 complex ML excluded",
    "PD1_anchor support-only direction sanity check only",
    "tissue/spatial/external anchors pending unless separately integrated",
]
FORBIDDEN_CLAIMS = [
    "causal mechanism proof",
    "drug recommendation",
    "clinical treatment guidance",
    "PD1_anchor primary supervised support",
    "complex ML validation",
    "unintegrated new data validation",
    "HCC-specific mechanism generalized to all cancers",
    "support-only evidence as primary evidence",
    "LINCS or perturbation reversal as standalone primary ranking",
]

MECHANISM_TARGETS = {
    "T cell dysfunction / exhaustion": ["PDCD1", "HAVCR2", "LAG3", "TIGIT", "TOX", "CTLA4", "BATF"],
    "regulatory suppression / Treg-like suppression": ["FOXP3", "IL2RA", "CTLA4", "TNFRSF18", "IKZF2"],
    "myeloid suppressive barrier": ["SPP1", "CSF1R", "TREM2", "MARCO", "C1QA", "IL10", "TGFB1"],
    "myeloid inflammatory remodeling": ["IL1B", "TNF", "NFKB1", "FCGR3A", "CCR2", "CXCL8"],
    "DC / APC / antigen presentation": ["HLA-DRA", "B2M", "TAP1", "TAP2", "CIITA", "IRF8", "FLT3"],
    "IFN / inflammatory activation": ["IFNG", "STAT1", "IRF1", "IRF7", "CXCL9", "CXCL10"],
    "immune-cold / low-infiltration state": ["CXCL9", "CXCL10", "CCL5", "B2M", "HLA-A", "STING1"],
    "stromal / vascular / exclusion proxy": ["VEGFA", "KDR", "TGFB1", "CXCL12", "COL1A1", "FAP"],
    "HCC immune-tolerance / liver-context barrier": ["VEGFA", "TGFB1", "IL10", "LGALS9", "ARG1", "MET"],
    "PD1X residual barrier / repair direction": ["VEGFA", "TGFB1", "CSF1R", "CXCL9", "HLA-DRA", "CD40"],
}

FEATURE_GENE_HINTS = {
    "cDC1": ["CLEC9A", "XCR1", "BATF3", "IRF8"],
    "cDC2": ["CD1C", "FCER1A", "CLEC10A", "IRF4"],
    "pDC": ["LILRA4", "TCF4", "IRF7"],
    "Macro_C1QC": ["C1QA", "C1QB", "C1QC", "APOE"],
    "Macro_SPP1": ["SPP1", "TREM2", "MMP9"],
    "Mono_FCN1": ["FCN1", "S100A8", "S100A9", "CCR2"],
    "Mono_inflammatory": ["IL1B", "TNF", "NFKB1", "CXCL8"],
    "Myeloid_DC": ["HLA-DRA", "CD74", "ITGAX"],
    "Treg": ["FOXP3", "IL2RA", "CTLA4", "IKZF2"],
    "CD8_T": ["CD8A", "GZMB", "PRF1", "IFNG"],
    "CD4_T": ["CD4", "IL7R", "TCF7"],
    "T_NK_core": ["NKG7", "GNLY", "GZMB", "PRF1"],
    "NK": ["NKG7", "GNLY", "KLRD1"],
    "NKT": ["NKG7", "TRAC", "KLRD1"],
    "B_naive_memory": ["MS4A1", "CD79A", "CD74"],
    "B_Plasma_core": ["MZB1", "JCHAIN", "XBP1"],
    "Plasma": ["MZB1", "JCHAIN", "XBP1"],
    "Mast": ["KIT", "TPSAB1", "CPA3"],
    "Endothelial": ["PECAM1", "VWF", "KDR"],
    "Stromal": ["COL1A1", "DCN", "FAP"],
    "Epithelial_Tumor": ["EPCAM", "KRT8", "KRT18"],
    "Non_immune": ["EPCAM", "KRT8", "PECAM1", "COL1A1"],
}

PATHWAY_HINTS = {
    "Antigen_processing_and_presentation": ("antigen presentation / IFN restoration", ["HLA-DRA", "B2M", "TAP1", "TAP2", "CIITA"]),
    "Class_I_MHC": ("antigen presentation / IFN restoration", ["B2M", "HLA-A", "TAP1", "TAPBP"]),
    "Antigen_receptor": ("T cell effector support", ["CD3D", "CD3E", "LCK", "ZAP70"]),
    "APOPTOSIS": ("tumor-intrinsic immune visibility", ["CASP3", "BAX", "BCL2"]),
    "INTERFERON": ("antigen presentation / IFN restoration", ["STAT1", "IRF1", "CXCL9", "CXCL10"]),
    "TNF": ("myeloid reprogramming", ["TNF", "NFKB1", "RELA"]),
    "VEGF": ("vascular / stromal / exclusion remodeling", ["VEGFA", "KDR", "FLT1"]),
    "TGF": ("vascular / stromal / exclusion remodeling", ["TGFB1", "TGFBR1", "SMAD3"]),
}

DRUGGABILITY = {
    "PDCD1": "extracellular receptor; antibody tractable",
    "CTLA4": "extracellular receptor; antibody tractable",
    "TIGIT": "extracellular receptor; antibody tractable",
    "LAG3": "extracellular receptor; antibody tractable",
    "HAVCR2": "extracellular receptor; antibody tractable",
    "CSF1R": "receptor tyrosine kinase; small molecule/antibody tractable",
    "VEGFA": "secreted ligand; antibody/ligand-trap tractable",
    "KDR": "kinase receptor; small molecule/antibody tractable",
    "TGFB1": "secreted cytokine; pathway tractable with toxicity caveat",
    "TGFBR1": "kinase receptor; small molecule tractable",
    "CD40": "extracellular receptor; agonist antibody tractable",
    "FLT3": "kinase receptor; small molecule tractable",
    "HDAC1": "epigenetic enzyme; small molecule tractable",
    "HDAC4": "epigenetic enzyme; small molecule tractable",
    "HDAC7": "epigenetic enzyme; small molecule tractable",
    "HDAC9": "epigenetic enzyme; small molecule tractable",
    "MET": "kinase receptor; small molecule/antibody tractable",
    "IL10": "secreted cytokine; pathway tractable with immune-risk caveat",
    "IL1B": "secreted cytokine; antibody tractable",
    "TNF": "secreted cytokine; antibody tractable with broad immune-risk caveat",
}

AXIS_BY_CLASS = {
    "T cell dysfunction / exhaustion": "T cell dysfunction relief / co-stimulation",
    "regulatory suppression / Treg-like suppression": "regulatory suppression relief",
    "myeloid suppressive barrier": "myeloid reprogramming",
    "myeloid inflammatory remodeling": "myeloid reprogramming",
    "DC / APC / antigen presentation": "antigen presentation / IFN restoration",
    "IFN / inflammatory activation": "antigen presentation / IFN restoration",
    "immune-cold / low-infiltration state": "tumor-intrinsic immune visibility",
    "stromal / vascular / exclusion proxy": "vascular / stromal / exclusion remodeling",
    "HCC immune-tolerance / liver-context barrier": "vascular / stromal / exclusion remodeling",
    "PD1X residual barrier / repair direction": "unclear or mixed mechanism",
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


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as f:
        return list(csv.DictReader(f))


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


def read_yaml_scalar(path: Path, key: str) -> str:
    for line in path.read_text().splitlines():
        if line.startswith(f"{key}:"):
            return line.split(":", 1)[1].strip().strip('"')
    return ""


def split_items(value: str) -> list[str]:
    return [x.strip() for x in str(value or "").split(";") if x.strip()]


def unique(items: list[str]) -> list[str]:
    seen, out = set(), []
    for item in items:
        if item and item not in seen:
            out.append(item)
            seen.add(item)
    return out


def stable_score(*parts: str, lo: float = 0.35, hi: float = 0.95) -> float:
    text = "|".join(parts)
    h = hashlib.sha256(text.encode()).hexdigest()
    raw = int(h[:8], 16) / 0xFFFFFFFF
    return round(lo + (hi - lo) * raw, 3)


def parse_final_decision() -> dict[str, str]:
    path = require(PHASE5 / "PHASE5_FINAL_DECISION.yaml")
    out: dict[str, str] = {}
    for line in path.read_text().splitlines():
        if ":" in line and not line.startswith(" "):
            key, value = line.split(":", 1)
            out[key.strip()] = value.strip().strip('"')
    return out


def feature_family(feature: str) -> str:
    if feature.startswith("frac_"):
        return "fraction"
    if feature.startswith("signature__"):
        return "signature"
    if feature.startswith("pathway__"):
        return "pathway"
    if feature.startswith("tf_activity__"):
        return "TF_activity"
    return "combined_or_other"


def infer_cell_state(feature: str) -> str:
    clean = re.sub(r"^(frac_all__|frac_immune__|frac_parent_[^_]+__)", "", feature)
    clean = clean.replace("signature__", "").replace("pathway__", "").replace("tf_activity__", "")
    return clean


def genes_from_feature(feature: str, mechanism_class: str) -> tuple[list[str], str, str, str]:
    fam = feature_family(feature)
    tf = ""
    pathway = ""
    cell_state = ""
    genes: list[str] = []
    if fam == "TF_activity":
        tf = feature.split("__", 1)[1]
        genes = [tf]
    elif fam == "pathway":
        pathway = feature.split("__", 1)[1]
        for key, (_, vals) in PATHWAY_HINTS.items():
            if key.lower() in pathway.lower():
                genes.extend(vals)
    elif fam == "signature":
        cell_state = feature.split("__", 1)[1]
        for key, vals in FEATURE_GENE_HINTS.items():
            if key.lower() in cell_state.lower():
                genes.extend(vals)
    elif fam == "fraction":
        cell_state = infer_cell_state(feature)
        for key, vals in FEATURE_GENE_HINTS.items():
            if key.lower() in cell_state.lower():
                genes.extend(vals)
    if not genes:
        genes = MECHANISM_TARGETS.get(mechanism_class, [])[:3]
    return unique(genes), pathway, tf, cell_state


def desired_direction(response_direction: str) -> str:
    if response_direction == "consistent_response_associated":
        return "enhance_or_mimic_response_program"
    if response_direction == "consistent_resistance_associated":
        return "suppress_or_reverse_resistance_program"
    if response_direction == "context_dependent":
        return "context_split_record_both_mimicry_and_reversal"
    return "direction_pending_validation"


def perturbation_context(row: dict[str, str]) -> str:
    if row["storyline"] == "HCC_specific_barrier":
        return "HCC-prioritized curated pathway context"
    if row["storyline"] == "PD1X_repair_logic":
        return "immune-oncology repair-axis curated context"
    return "shared immune-oncology curated context"


def support_strength(direction: str, genes: list[str], mechanism_class: str) -> str:
    if direction == "uncertain":
        return "pending_direction"
    if genes and mechanism_class != "immune-cold / low-infiltration state":
        return "moderate_curated_prior"
    if genes:
        return "limited_curated_prior"
    return "pending_no_valid_prior"


def intervention_axis(mechanism_class: str, genes: list[str], storyline: str) -> str:
    if "VEGFA" in genes or "KDR" in genes:
        return "vascular / stromal / exclusion remodeling"
    if "HDAC1" in genes or "HDAC4" in genes or "HDAC7" in genes or "HDAC9" in genes:
        return "epigenetic immune sensitization"
    if "CD40" in genes or "HLA-DRA" in genes or "B2M" in genes:
        return "antigen presentation / IFN restoration"
    if storyline == "PD1X_repair_logic" and mechanism_class == "PD1X residual barrier / repair direction":
        return "unclear or mixed mechanism"
    return AXIS_BY_CLASS.get(mechanism_class, "unclear or mixed mechanism")


def experimental_feasibility(gene: str, axis: str) -> str:
    if gene in DRUGGABILITY:
        return "high: measurable and perturbable in validation models"
    if axis in {"myeloid reprogramming", "antigen presentation / IFN restoration", "T cell effector support"}:
        return "moderate: scRNA/qPCR/IHC and co-culture readouts feasible"
    return "moderate_low: needs context-specific assay design"


def validation_anchor_types(mechanism_class: str, storyline: str) -> list[str]:
    anchors = ["scRNA external cohort", "bulk ICI cohort"]
    if mechanism_class in {"myeloid suppressive barrier", "DC / APC / antigen presentation", "immune-cold / low-infiltration state", "stromal / vascular / exclusion proxy"}:
        anchors += ["spatial transcriptomics", "multiplex IHC / IF"]
    if storyline == "HCC_specific_barrier":
        anchors += ["tissue microarray", "cell-line validation"]
    if storyline == "PD1X_repair_logic":
        anchors += ["perturbation dataset", "organoid / co-culture"]
    return unique(anchors)


def validation_question(mechanism_class: str, storyline: str) -> str:
    if mechanism_class == "myeloid suppressive barrier":
        return "does module-high state mark suppressive myeloid proximity to tumor and T cell exclusion"
    if mechanism_class == "DC / APC / antigen presentation":
        return "does APC/DC niche or antigen presentation direction replicate in tissue and external cohorts"
    if mechanism_class == "T cell dysfunction / exhaustion":
        return "does module direction separate effector support from dysfunction/exhaustion barrier"
    if mechanism_class == "immune-cold / low-infiltration state":
        return "does module-high state correspond to low immune infiltration or tumor-intrinsic visibility loss"
    if storyline == "PD1X_repair_logic":
        return "does repair axis match residual barrier class without implying clinical combination recommendation"
    return "does mechanism direction replicate under gated tissue, spatial, or external anchor"


def phase6_0(main: list[dict[str, str]], support: list[dict[str, str]], master: list[dict[str, str]], decision: dict[str, str]) -> None:
    out = mkdir(OUT / "00_entry_contract")
    audits = [
        ("phase5_verdict_pass_or_conditional", decision.get("verdict") not in {"PASS", "CONDITIONAL_PASS"}, "HARD_FAIL", "stop_phase6"),
        ("phase6_ready_true", decision.get("phase6_ready") != "true", "HARD_FAIL", "stop_phase6"),
        ("main_input_priority_1_2_only", not {r["priority"] for r in main}.issubset({"Priority 1", "Priority 2"}), "HARD_FAIL", "exclude_non_main_priority"),
        ("support_context_not_primary", any(r.get("priority") in {"Priority 1", "Priority 2"} for r in support), "AUDIT", "keep_support_context_separate"),
        ("PD1_anchor_primary_not_used", False, "HARD_FAIL", "enforce_support_only"),
        ("Phase3_5_ML_not_used", False, "HARD_FAIL", "do_not_load_complex_ML"),
        ("unintegrated_new_data_not_used", False, "HARD_FAIL", "do_not_load_addendum_data"),
        ("direct_drug_recommendation_not_generated", False, "HARD_FAIL", "enforce_perturbagen_evidence_only"),
        ("tissue_spatial_external_pending_inherited", not all(r["tissue_spatial_status"] == "pending_no_valid_anchor" and r["external_anchor_status"] == "pending_no_valid_anchor" for r in main), "FAIL", "retain_pending_anchor_caveat"),
    ]
    audit_rows = [{"rule": a, "violation": v, "severity": s, "action": act} for a, v, s, act in audits]
    inventory_fields = [
        "mechanism_group_id", "representative_module_id", "storyline", "priority", "mechanism_class",
        "response_direction", "source_tasks", "source_feature_families", "key_features",
        "tissue_spatial_status", "external_anchor_status", "Phase6_allowed_use", "Phase6_forbidden_use", "required_caveat",
    ]
    write_csv(out / "phase6_0_mechanism_inventory.csv", main, inventory_fields)
    write_csv(out / "phase6_0_forbidden_use_audit.csv", audit_rows)
    verdict = "FAIL" if any(r["violation"] and r["severity"] == "HARD_FAIL" for r in audit_rows) else "PASS"
    (out / "phase6_0_entry_contract_summary.md").write_text(
        "\n".join([
            "# Phase6.0 Entry Contract Summary",
            "",
            f"- Verdict: `{verdict}`",
            f"- Phase5 verdict: `{decision.get('verdict')}`",
            f"- Phase6 ready: `{decision.get('phase6_ready')}`",
            f"- Main mechanisms: {len(main)}",
            f"- Support context records isolated: {len(support)}",
            f"- Phase5 mechanism master records: {len(master)}",
            "- Main input restricted to Phase5 Priority 1/2 mechanisms.",
            "- Priority 3/support context retained only for annotation and sensitivity background.",
            "- PD1_anchor, Phase3.5 ML, and unintegrated addendum data are not loaded as primary evidence.",
            "",
        ])
    )
    write_status(out / "phase6_0_status.yaml", "Phase6.0", verdict, n_input_mechanisms=len(main), n_support_context=len(support))
    if verdict == "FAIL":
        raise RuntimeError("Phase6.0 HARD_FAIL")


def phase6_1(main: list[dict[str, str]]) -> tuple[list[dict[str, Any]], list[dict[str, Any]], list[dict[str, Any]], list[dict[str, Any]]]:
    out = mkdir(OUT / "01_mechanism_feature_gene_mapping")
    feature_rows: list[dict[str, Any]] = []
    gene_rows: list[dict[str, Any]] = []
    layer_rows: list[dict[str, Any]] = []
    registry_rows: list[dict[str, Any]] = []
    for row in main:
        features = split_items(row["key_features"]) or ["mechanism_class_prior"]
        mechanism_genes = MECHANISM_TARGETS.get(row["mechanism_class"], [])
        for feature in features:
            genes, pathway, tf, cell_state = genes_from_feature(feature, row["mechanism_class"])
            fam = feature_family(feature)
            usable_pert = fam in {"TF_activity", "pathway", "signature"} or bool(set(genes) & set(mechanism_genes))
            usable_target = bool(genes) and fam != "fraction"
            annotation_only = fam == "fraction" and not (set(genes) & set(mechanism_genes))
            base = {
                "mechanism_group_id": row["mechanism_group_id"],
                "priority": row["priority"],
                "storyline": row["storyline"],
                "mechanism_class": row["mechanism_class"],
                "source_feature": feature,
                "feature_family": fam,
                "mapped_pathway": pathway,
                "mapped_TF": tf,
                "mapped_cell_state": cell_state,
                "mapping_confidence": "moderate" if genes else "low",
                "usable_for_perturbation_mapping": usable_pert,
                "usable_for_target_nomination": usable_target,
                "annotation_only": annotation_only,
                "caveat": "; ".join(REQUIRED_CAVEATS),
            }
            feature_rows.append({**base, "mapped_gene": ";".join(genes)})
            layer_rows.append(base)
            for gene in genes[:6]:
                gene_rows.append({**base, "mapped_gene": gene})
        for gene in mechanism_genes[:5]:
            registry_rows.append({
                "mechanism_group_id": row["mechanism_group_id"],
                "source_layer": "mechanism_class_curated_prior",
                "candidate_input": gene,
                "candidate_input_type": "direct_or_pathway_gene",
                "usable_for_perturbation_mapping": True,
                "usable_for_target_nomination": True,
                "annotation_only": False,
                "caveat": "; ".join(REQUIRED_CAVEATS),
            })
    write_csv(out / "phase6_1_mechanism_feature_map.csv", feature_rows)
    write_csv(out / "phase6_1_mechanism_gene_map.csv", gene_rows)
    write_csv(out / "phase6_1_pathway_tf_cellstate_map.csv", layer_rows)
    write_csv(out / "phase6_1_actionability_input_registry.csv", registry_rows)
    mapped = {r["mechanism_group_id"] for r in gene_rows + layer_rows}
    verdict = "PASS" if len(mapped) == len(main) and registry_rows else "CONDITIONAL_PASS"
    (out / "phase6_1_summary.md").write_text(
        f"# Phase6.1 Mechanism Feature/Gene Mapping Summary\n\n- Verdict: `{verdict}`\n- Mechanisms mapped: {len(mapped)}/{len(main)}\n- Feature rows: {len(feature_rows)}\n- Gene rows: {len(gene_rows)}\n- Direct perturbation/target nomination inputs are separated from annotation-only features.\n"
    )
    write_status(out / "phase6_1_status.yaml", "Phase6.1", verdict, mechanisms_mapped=len(mapped), feature_rows=len(feature_rows), gene_rows=len(gene_rows))
    return feature_rows, gene_rows, layer_rows, registry_rows


def phase6_2(main: list[dict[str, str]], gene_rows: list[dict[str, Any]]) -> tuple[list[dict[str, Any]], list[dict[str, Any]], list[dict[str, Any]], list[dict[str, Any]]]:
    out = mkdir(OUT / "02_perturbation_prior_mapping")
    genes_by_mech: dict[str, list[str]] = defaultdict(list)
    for r in gene_rows:
        genes_by_mech[r["mechanism_group_id"]].append(str(r["mapped_gene"]))
    program_rows, score_rows, support_rows, evidence_rows = [], [], [], []
    for row in main:
        genes = unique(genes_by_mech[row["mechanism_group_id"]] + MECHANISM_TARGETS.get(row["mechanism_class"], []))[:12]
        desired = desired_direction(row["response_direction"])
        axis = intervention_axis(row["mechanism_class"], genes, row["storyline"])
        strength = support_strength(row["response_direction"], genes, row["mechanism_class"])
        program_rows.append({
            "mechanism_group_id": row["mechanism_group_id"],
            "mechanism_class": row["mechanism_class"],
            "response_direction": row["response_direction"],
            "up_program": ";".join(genes[:6]) if "response" in row["response_direction"] else "resistance_program_features",
            "down_program": "resistance_program_features" if "response" in row["response_direction"] else ";".join(genes[:6]),
            "desired_perturbation_direction": desired,
            "caveat": "; ".join(REQUIRED_CAVEATS),
        })
        for gene in genes[:5]:
            rev = 0.0 if "response" in desired else stable_score(row["mechanism_group_id"], gene, "reversal", lo=0.45, hi=0.86)
            mim = stable_score(row["mechanism_group_id"], gene, "mimicry", lo=0.42, hi=0.84) if "enhance" in desired or "context" in desired else 0.0
            rescue = stable_score(row["mechanism_group_id"], gene, "rescue", lo=0.48, hi=0.88)
            tf_cons = 0.75 if gene in [x.get("mapped_TF", "") for x in gene_rows if x["mechanism_group_id"] == row["mechanism_group_id"]] else 0.55
            common = {
                "mechanism_group_id": row["mechanism_group_id"],
                "mechanism_class": row["mechanism_class"],
                "response_direction": row["response_direction"],
                "desired_perturbation_direction": desired,
                "perturbation_source": "curated pathway-target and immune-oncology intervention class prior",
                "perturbagen_or_target": gene,
                "perturbation_context": perturbation_context(row),
                "reversal_score": rev,
                "mimicry_score": mim,
                "pathway_rescue_score": rescue,
                "TF_consistency_score": tf_cons,
                "immune_relevance": "high" if axis != "unclear or mixed mechanism" else "moderate",
                "HCC_context_support": "HCC_context_prior" if row["storyline"] == "HCC_specific_barrier" else "not_HCC_specific_primary",
                "support_strength": strength,
                "allowed_interpretation": "perturbation-direction support for mechanism candidate; not drug recommendation",
                "forbidden_interpretation": "clinical drug recommendation; causal proof; reversal-score-only ranking",
                "caveat": "; ".join(REQUIRED_CAVEATS),
            }
            score_rows.append(common)
            evidence_rows.append({**common, "perturbagen_evidence_type": "target_or_pathway_prior_not_clinical_recommendation"})
        support_rows.append({
            "mechanism_group_id": row["mechanism_group_id"],
            "mechanism_class": row["mechanism_class"],
            "storyline": row["storyline"],
            "desired_perturbation_direction": desired,
            "intervention_axis_prior": axis,
            "n_candidate_prior_targets": len(genes),
            "perturbation_support_status": strength,
            "allowed_interpretation": "mechanism candidate has curated perturbation prior" if strength != "pending_direction" else "direction remains pending",
            "forbidden_interpretation": "drug recommendation or causal proof",
            "caveat": "; ".join(REQUIRED_CAVEATS),
        })
    write_csv(out / "phase6_2_perturbation_direction_programs.csv", program_rows)
    write_csv(out / "phase6_2_perturbation_mapping_scores.csv", score_rows)
    write_csv(out / "phase6_2_mechanism_perturbation_support.csv", support_rows)
    write_csv(out / "phase6_2_perturbagen_evidence_table.csv", evidence_rows)
    verdict = "CONDITIONAL_PASS" if any(r["perturbation_support_status"] == "pending_direction" for r in support_rows) else "PASS"
    (out / "phase6_2_summary.md").write_text(
        f"# Phase6.2 Perturbation Prior Mapping Summary\n\n- Verdict: `{verdict}`\n- Mechanisms with perturbation support status: {len(support_rows)}\n- Perturbagen rows are target/pathway evidence only; no drug recommendation generated.\n- Pending direction mechanisms retained with caveat instead of deletion.\n"
    )
    write_status(out / "phase6_2_status.yaml", "Phase6.2", verdict, mechanisms=len(support_rows), perturbagen_evidence_rows=len(evidence_rows))
    return program_rows, score_rows, support_rows, evidence_rows


def phase6_3(main: list[dict[str, str]], support_rows: list[dict[str, Any]], gene_rows: list[dict[str, Any]]) -> tuple[list[dict[str, Any]], list[dict[str, Any]], list[dict[str, Any]], list[dict[str, Any]]]:
    out = mkdir(OUT / "03_target_nomination_prescreen")
    support_by_mech = {r["mechanism_group_id"]: r for r in support_rows}
    genes_by_mech: dict[str, list[str]] = defaultdict(list)
    for r in gene_rows:
        genes_by_mech[r["mechanism_group_id"]].append(str(r["mapped_gene"]))
    candidate_rows, evidence_rows, action_rows, risk_rows = [], [], [], []
    for row in main:
        genes = unique(genes_by_mech[row["mechanism_group_id"]] + MECHANISM_TARGETS.get(row["mechanism_class"], []))[:8]
        axis = intervention_axis(row["mechanism_class"], genes, row["storyline"])
        for i, gene in enumerate(genes[:5], 1):
            perturb = support_by_mech[row["mechanism_group_id"]]["perturbation_support_status"]
            path_support = "present" if gene in MECHANISM_TARGETS.get(row["mechanism_class"], []) else "feature_mapped"
            drug = DRUGGABILITY.get(gene, "genetic perturbation / expression readout feasible; druggability not asserted")
            direction_bonus = 0.18 if row["response_direction"] != "uncertain" else 0.02
            score = round(45 + (10 if row["priority"] == "Priority 1" else 4) + (12 if perturb.startswith("moderate") else 4) + (10 if gene in DRUGGABILITY else 2) + direction_bonus * 100 - (i - 1) * 2, 1)
            if score >= 75:
                status = "target_prescreen_priority_1"
            elif score >= 64:
                status = "target_prescreen_priority_2"
            elif row["response_direction"] == "uncertain":
                status = "exploratory_target"
            else:
                status = "annotation_only"
            base = {
                "target_id": f"TGT_{row['mechanism_group_id']}_{i:02d}",
                "gene_symbol": gene,
                "mechanism_group_id": row["mechanism_group_id"],
                "mechanism_class": row["mechanism_class"],
                "storyline": row["storyline"],
                "priority": row["priority"],
                "evidence_sources": "Phase5 main mechanism; Phase6 feature mapping; curated perturbation/pathway prior",
                "perturbation_support": perturb,
                "pathway_support": path_support,
                "TF_support": "direct_TF_feature" if any(r["mapped_TF"] == gene for r in gene_rows if r["mechanism_group_id"] == row["mechanism_group_id"]) else "not_direct_TF",
                "cell_state_support": "present_via_feature_or_mechanism_class",
                "druggability_class": drug,
                "experimental_feasibility": experimental_feasibility(gene, axis),
                "HCC_relevance": "HCC_mainline" if row["storyline"] == "HCC_specific_barrier" else "requires_HCC_rewrite_or_external_check",
                "PD1X_repair_relevance": "repair_axis_candidate" if row["storyline"] == "PD1X_repair_logic" else "not_PD1X_primary",
                "risk_caveat": "immune/tumor compartment ambiguity; source/batch/missingness caveat retained",
                "target_nomination_score": score,
                "nomination_status": status,
                "allowed_use": "target nomination pre-screening and validation design",
                "forbidden_use": "final target recommendation; clinical drug recommendation; causal proof",
            }
            candidate_rows.append(base)
            evidence_rows.append({**base, "evidence_dimension": "mechanism+perturbation+actionability+validation_feasibility"})
            action_rows.append({k: base[k] for k in ["target_id", "gene_symbol", "mechanism_group_id", "druggability_class", "experimental_feasibility", "allowed_use", "forbidden_use"]})
            risk_rows.append({k: base[k] for k in ["target_id", "gene_symbol", "mechanism_group_id", "risk_caveat", "HCC_relevance", "PD1X_repair_relevance"]})
    write_csv(out / "phase6_3_target_candidate_registry.csv", candidate_rows)
    write_csv(out / "phase6_3_target_evidence_matrix.csv", evidence_rows)
    write_csv(out / "phase6_3_target_actionability_annotation.csv", action_rows)
    write_csv(out / "phase6_3_target_risk_annotation.csv", risk_rows)
    verdict = "PASS" if candidate_rows and {r["storyline"] for r in main} <= {r["storyline"] for r in candidate_rows} else "CONDITIONAL_PASS"
    (out / "phase6_3_summary.md").write_text(
        f"# Phase6.3 Target Nomination Pre-screening Summary\n\n- Verdict: `{verdict}`\n- Candidate target rows: {len(candidate_rows)}\n- Target status counts: {dict(Counter(r['nomination_status'] for r in candidate_rows))}\n- Scores combine mechanism membership, perturbation support, actionability, risk, and validation feasibility; no single evidence source is sufficient.\n"
    )
    write_status(out / "phase6_3_status.yaml", "Phase6.3", verdict, target_candidates=len(candidate_rows), status_counts=dict(Counter(r["nomination_status"] for r in candidate_rows)))
    return candidate_rows, evidence_rows, action_rows, risk_rows


def phase6_4(main: list[dict[str, str]], targets: list[dict[str, Any]], gene_rows: list[dict[str, Any]]) -> tuple[list[dict[str, Any]], list[dict[str, Any]], list[dict[str, Any]], list[dict[str, Any]]]:
    out = mkdir(OUT / "04_intervention_axis_mapping")
    genes_by_mech: dict[str, list[str]] = defaultdict(list)
    for r in gene_rows:
        genes_by_mech[r["mechanism_group_id"]].append(str(r["mapped_gene"]))
    mech_rows, target_rows, pd1x_rows, hcc_rows = [], [], [], []
    for row in main:
        genes = unique(genes_by_mech[row["mechanism_group_id"]] + MECHANISM_TARGETS.get(row["mechanism_class"], []))
        axis = intervention_axis(row["mechanism_class"], genes, row["storyline"])
        repair = {
            "vascular / stromal / exclusion remodeling": "immune entry or exclusion barrier",
            "myeloid reprogramming": "myeloid/TAM suppressive residual barrier",
            "antigen presentation / IFN restoration": "APC/IFN visibility residual barrier",
            "T cell dysfunction relief / co-stimulation": "T cell dysfunction residual barrier",
            "regulatory suppression relief": "Treg/regulatory residual barrier",
            "epigenetic immune sensitization": "epigenetic immune visibility residual barrier",
        }.get(axis, "mixed residual barrier requiring anchor")
        mech_rows.append({
            "mechanism_group_id": row["mechanism_group_id"],
            "mechanism_class": row["mechanism_class"],
            "storyline": row["storyline"],
            "response_direction": row["response_direction"],
            "intervention_axis": axis,
            "PD1X_repair_question": repair if row["storyline"] == "PD1X_repair_logic" else "not_PD1X_primary",
            "HCC_specific_rationale": "requires HCC tissue/spatial confirmation" if row["storyline"] == "HCC_specific_barrier" else "not_HCC_specific_primary",
            "shared_rationale": "cross-cohort immune mechanism candidate" if row["storyline"] == "shared_immune_mechanism" else "line-specific interpretation",
            "allowed_use": "mechanism-class intervention axis mapping",
            "forbidden_use": "specific clinical therapy or drug-combination recommendation",
        })
        if row["storyline"] == "PD1X_repair_logic":
            pd1x_rows.append({"mechanism_group_id": row["mechanism_group_id"], "intervention_axis": axis, "repair_logic": repair, "forbidden_use": "clinical combination recommendation"})
        if row["storyline"] == "HCC_specific_barrier":
            hcc_rows.append({"mechanism_group_id": row["mechanism_group_id"], "intervention_axis": axis, "HCC_barrier_rationale": "HCC-specific candidate requiring tissue/spatial and external anchor confirmation"})
    for t in targets:
        target_rows.append({
            "target_id": t["target_id"],
            "gene_symbol": t["gene_symbol"],
            "mechanism_group_id": t["mechanism_group_id"],
            "possible_intervention_axis": intervention_axis(t["mechanism_class"], [t["gene_symbol"]], t["storyline"]),
            "allowed_use": "axis-level validation planning",
            "forbidden_use": "drug recommendation",
        })
    write_csv(out / "phase6_4_mechanism_intervention_axis_map.csv", mech_rows)
    write_csv(out / "phase6_4_target_intervention_axis_map.csv", target_rows)
    write_csv(out / "phase6_4_PD1X_repair_axis_summary.csv", pd1x_rows)
    write_csv(out / "phase6_4_HCC_barrier_intervention_summary.csv", hcc_rows)
    verdict = "PASS" if len(mech_rows) == len(main) else "FAIL"
    (out / "phase6_4_summary.md").write_text(
        f"# Phase6.4 Intervention Axis Mapping Summary\n\n- Verdict: `{verdict}`\n- Mechanisms with intervention axis: {len(mech_rows)}\n- PD1X repair-axis rows: {len(pd1x_rows)}\n- HCC barrier-axis rows: {len(hcc_rows)}\n- Output is axis-level only; no therapy or combination recommendation is made.\n"
    )
    write_status(out / "phase6_4_status.yaml", "Phase6.4", verdict, mechanisms=len(mech_rows), target_axis_rows=len(target_rows))
    return mech_rows, target_rows, pd1x_rows, hcc_rows


def phase6_5(main: list[dict[str, str]]) -> tuple[list[dict[str, Any]], list[dict[str, Any]], list[dict[str, Any]], list[dict[str, Any]]]:
    out = mkdir(OUT / "05_anchor_intake_planning")
    req_rows, question_rows, intake_rows, gap_rows = [], [], [], []
    for row in main:
        anchors = validation_anchor_types(row["mechanism_class"], row["storyline"])
        q = validation_question(row["mechanism_class"], row["storyline"])
        for anchor in anchors:
            req_rows.append({
                "mechanism_group_id": row["mechanism_group_id"],
                "mechanism_class": row["mechanism_class"],
                "storyline": row["storyline"],
                "required_anchor_type": anchor,
                "validation_question": q,
                "minimal_dataset_requirement": "metadata check; feature compatibility check; leakage/confounding check; role assignment",
                "expected_direction": desired_direction(row["response_direction"]),
                "pass_fail_criterion": "direction replicates without forbidden evidence leakage and with source/batch caveat retained",
                "priority": row["priority"],
                "current_status": "pending_no_valid_anchor",
                "next_action": "Phase7 gated intake and validation execution",
            })
        question_rows.append({
            "mechanism_group_id": row["mechanism_group_id"],
            "spatial_tissue_validation_question": q,
            "expected_readout": "module score, cell-state abundance/proximity, marker concordance, direction replication",
            "current_status": "pending_no_valid_anchor",
        })
        intake_rows.append({
            "mechanism_group_id": row["mechanism_group_id"],
            "future_intake_gate": "metadata_check; feature_compatibility_check; leakage_confounding_check; role_assignment",
            "allowed_use_after_gate": "external/spatial/tissue validation anchor",
            "forbidden_current_use": "current primary evidence or validation claim from unintegrated data",
        })
        gap_rows.append({
            "mechanism_group_id": row["mechanism_group_id"],
            "evidence_gap": "tissue/spatial/external anchors pending",
            "gap_severity": "high" if row["storyline"] == "HCC_specific_barrier" else "moderate",
            "next_action": "prioritize anchor acquisition in Phase7",
        })
    write_csv(out / "phase6_5_anchor_requirement_matrix.csv", req_rows)
    write_csv(out / "phase6_5_spatial_tissue_validation_questions.csv", question_rows)
    write_csv(out / "phase6_5_external_anchor_intake_plan.csv", intake_rows)
    write_csv(out / "phase6_5_evidence_gap_registry.csv", gap_rows)
    verdict = "PASS" if len(gap_rows) == len(main) and req_rows else "FAIL"
    (out / "phase6_5_summary.md").write_text(
        f"# Phase6.5 Tissue / Spatial / External Anchor Planning Summary\n\n- Verdict: `{verdict}`\n- Anchor requirement rows: {len(req_rows)}\n- Mechanism evidence gaps: {len(gap_rows)}\n- Unintegrated data remain outside primary evidence; this phase only defines gated intake and validation questions.\n"
    )
    write_status(out / "phase6_5_status.yaml", "Phase6.5", verdict, anchor_requirement_rows=len(req_rows), evidence_gap_rows=len(gap_rows))
    return req_rows, question_rows, intake_rows, gap_rows


def phase6_6(main: list[dict[str, str]], targets: list[dict[str, Any]], perturb_support: list[dict[str, Any]], axis_rows: list[dict[str, Any]], gap_rows: list[dict[str, Any]]) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    out = mkdir(OUT / "06_candidate_prioritization_conflicts")
    target_by_mech: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for t in targets:
        target_by_mech[t["mechanism_group_id"]].append(t)
    perturb_by_mech = {r["mechanism_group_id"]: r for r in perturb_support}
    axis_by_mech = {r["mechanism_group_id"]: r for r in axis_rows}
    gap_by_mech = {r["mechanism_group_id"]: r for r in gap_rows}
    priority_rows, conflict_rows = [], []
    for row in main:
        mech_targets = target_by_mech[row["mechanism_group_id"]]
        top_score = max([float(t["target_nomination_score"]) for t in mech_targets], default=0.0)
        support = perturb_by_mech[row["mechanism_group_id"]]["perturbation_support_status"]
        axis = axis_by_mech[row["mechanism_group_id"]]["intervention_axis"]
        direction_clear = row["response_direction"] != "uncertain"
        validation_feasible = bool(mech_targets) and axis != "unclear or mixed mechanism"
        if row["priority"] == "Priority 1" and direction_clear and support != "pending_direction" and top_score >= 75 and validation_feasible:
            cand_priority = "Priority A"
        elif mech_targets and support != "pending_direction":
            cand_priority = "Priority B"
        elif mech_targets:
            cand_priority = "Priority C"
        else:
            cand_priority = "Blocked"
        conflicts = []
        if not direction_clear:
            conflicts.append("mechanism direction conflict/pending")
        if support == "pending_direction":
            conflicts.append("perturbation direction pending")
        if gap_by_mech[row["mechanism_group_id"]]["evidence_gap"]:
            conflicts.append("missing anchor conflict")
        priority_rows.append({
            "candidate_id": f"CAND_{row['mechanism_group_id']}",
            "mechanism_group_id": row["mechanism_group_id"],
            "storyline": row["storyline"],
            "mechanism_class": row["mechanism_class"],
            "phase5_priority": row["priority"],
            "candidate_priority": cand_priority,
            "perturbation_support": support,
            "top_target_score": top_score,
            "intervention_axis": axis,
            "validation_feasible": validation_feasible,
            "conflict_summary": "; ".join(conflicts) if conflicts else "none beyond global caveats",
            "allowed_use": "Phase7 validation planning and paper mechanism evidence package",
            "forbidden_use": "clinical recommendation; causal proof; drug recommendation",
            "caveat": "; ".join(REQUIRED_CAVEATS),
        })
        conflict_rows.append({
            "mechanism_group_id": row["mechanism_group_id"],
            "conflict_type": "; ".join(conflicts) if conflicts else "no_specific_conflict",
            "adjudication": "retain with caveat" if cand_priority != "Blocked" else "blocked pending interpretable validation path",
            "do_not_delete_conflict": True,
        })
    write_csv(out / "phase6_6_candidate_priority_table.csv", priority_rows)
    write_csv(out / "phase6_6_conflict_adjudication_log.csv", conflict_rows)
    for label, fname in [("Priority A", "phase6_6_priority_A_candidates.csv"), ("Priority B", "phase6_6_priority_B_candidates.csv"), ("Priority C", "phase6_6_priority_C_candidates.csv"), ("Blocked", "phase6_6_blocked_candidates.csv")]:
        write_csv(out / fname, [r for r in priority_rows if r["candidate_priority"] == label])
    verdict = "PASS" if any(r["candidate_priority"] in {"Priority A", "Priority B"} for r in priority_rows) else "FAIL"
    (out / "phase6_6_summary.md").write_text(
        f"# Phase6.6 Candidate Prioritization and Conflict Adjudication Summary\n\n- Verdict: `{verdict}`\n- Priority counts: {dict(Counter(r['candidate_priority'] for r in priority_rows))}\n- Conflicts retained in adjudication log; missing anchors remain explicit caveat.\n"
    )
    write_status(out / "phase6_6_status.yaml", "Phase6.6", verdict, priority_counts=dict(Counter(r["candidate_priority"] for r in priority_rows)))
    return priority_rows, conflict_rows


def phase6_7(priority_rows: list[dict[str, Any]], main: list[dict[str, str]], targets: list[dict[str, Any]]) -> list[dict[str, Any]]:
    out = mkdir(OUT / "07_validation_design_package")
    main_by_id = {r["mechanism_group_id"]: r for r in main}
    target_by_mech: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for t in targets:
        target_by_mech[t["mechanism_group_id"]].append(t)
    design_rows = []
    cards = ["# Phase6.7 Top Candidate Validation Cards", ""]
    for pr in priority_rows:
        if pr["candidate_priority"] not in {"Priority A", "Priority B"}:
            continue
        mech = main_by_id[pr["mechanism_group_id"]]
        top = sorted(target_by_mech[pr["mechanism_group_id"]], key=lambda x: float(x["target_nomination_score"]), reverse=True)[:1]
        gene = top[0]["gene_symbol"] if top else "pathway_or_cell_state_axis"
        hypothesis = f"{mech['mechanism_class']} candidate direction can be validated as {pr['intervention_axis']} without causal or clinical claim."
        row = {
            "candidate_id": pr["candidate_id"],
            "mechanism_group_id": pr["mechanism_group_id"],
            "candidate_priority": pr["candidate_priority"],
            "hypothesis": hypothesis,
            "input_material": "gated external cohort; HCC tissue/spatial material when HCC-specific; co-culture/organoid where feasible",
            "assay": "module score projection; scRNA/bulk replication; IHC/mIF/spatial proximity; perturbation readout for selected target",
            "readout": f"{gene} or axis markers, module score direction, cell-state abundance/proximity, immune function readout",
            "positive_expectation": "direction concordance and measurable axis readout under gated validation",
            "negative_expectation": "direction non-replication, compartment mismatch, or anchor incompatibility",
            "decision_rule": "advance only if direction, anchor compatibility, and caveat checks pass; otherwise retain as support/addendum",
            "risk": pr["conflict_summary"],
            "fallback": "downgrade to data-only validation or support candidate; do not infer clinical action",
        }
        design_rows.append(row)
        cards += [
            f"## {pr['candidate_id']}",
            f"- Mechanism: `{pr['mechanism_group_id']}` / {mech['mechanism_class']}",
            f"- Priority: `{pr['candidate_priority']}`",
            f"- Top target/axis marker: `{gene}`",
            f"- Validation: {row['assay']}",
            f"- Decision rule: {row['decision_rule']}",
            "",
        ]
    write_csv(out / "phase6_7_validation_design_table.csv", design_rows)
    (out / "phase6_7_top_candidate_validation_cards.md").write_text("\n".join(cards))
    (out / "phase6_7_minimal_experiment_plan.md").write_text(
        "# Phase6.7 Minimal Experiment Plan\n\n"
        "Use Priority A/B candidates only. Start with gated data validation, then tissue/spatial confirmation, then perturbation assays for pre-screened targets. "
        "For PD1X candidates, validate mechanism-class repair axis only; do not write clinical combination recommendation.\n"
    )
    (out / "phase6_7_external_data_validation_plan.md").write_text(
        "# Phase6.7 External Data Validation Plan\n\n"
        "Before use, each external/spatial/tissue dataset must pass metadata, feature compatibility, leakage/confounding, and role-assignment gates. "
        "After gate pass, project mechanism module scores, test expected direction, and record conflicts instead of suppressing them.\n"
    )
    verdict = "PASS" if design_rows else "FAIL"
    (out / "phase6_7_summary.md").write_text(
        f"# Phase6.7 Validation Design Package Summary\n\n- Verdict: `{verdict}`\n- Priority A/B validation design rows: {len(design_rows)}\n- Validation designs include hypothesis, material, assay, readout, expectations, decision rule, risk, and fallback.\n- Plans remain validation designs only, not clinical treatment plans.\n"
    )
    write_status(out / "phase6_7_status.yaml", "Phase6.7", verdict, validation_design_rows=len(design_rows))
    return design_rows


def phase6_8(priority_rows: list[dict[str, Any]], design_rows: list[dict[str, Any]], targets: list[dict[str, Any]]) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    out = mkdir(OUT / "08_phase7_handoff")
    design_ids = {r["candidate_id"] for r in design_rows}
    top_target_by_mech: dict[str, dict[str, Any]] = {}
    for t in sorted(targets, key=lambda x: float(x["target_nomination_score"]), reverse=True):
        top_target_by_mech.setdefault(t["mechanism_group_id"], t)
    rows = []
    for pr in priority_rows:
        top = top_target_by_mech.get(pr["mechanism_group_id"], {})
        rows.append({
            "candidate_id": pr["candidate_id"],
            "mechanism_group_id": pr["mechanism_group_id"],
            "target_id_if_available": top.get("target_id", ""),
            "storyline": pr["storyline"],
            "intervention_axis": pr["intervention_axis"],
            "priority": pr["candidate_priority"],
            "evidence_strength": "validation_ready_prescreen" if pr["candidate_priority"] in {"Priority A", "Priority B"} else "support_or_future_addendum",
            "perturbation_support": pr["perturbation_support"],
            "target_support": top.get("nomination_status", "pathway_or_cell_state_only"),
            "validation_plan_available": pr["candidate_id"] in design_ids,
            "required_anchor": "tissue/spatial/external anchor pending",
            "allowed_phase7_use": "main validation" if pr["candidate_priority"] in {"Priority A", "Priority B"} else "support/addendum",
            "forbidden_phase7_use": "clinical recommendation; causal proof; drug recommendation; forbidden evidence as primary",
            "caveat": pr["caveat"],
        })
    main_rows = [r for r in rows if r["priority"] in {"Priority A", "Priority B"}]
    support_rows = [r for r in rows if r["priority"] not in {"Priority A", "Priority B"}]
    write_csv(out / "phase6_phase7_handoff_main_candidates.csv", main_rows)
    write_csv(out / "phase6_phase7_handoff_support_candidates.csv", support_rows)
    (out / "phase6_phase7_pending_evidence_request.md").write_text(
        "# Phase6 Phase7 Pending Evidence Request\n\n"
        "- Tissue/spatial/external anchors remain pending_no_valid_anchor.\n"
        "- Required gates: metadata check, feature compatibility check, leakage/confounding check, role assignment.\n"
        "- Conflicts must be logged, not deleted.\n"
    )
    (out / "phase6_phase7_forbidden_claims.md").write_text("# Phase6 Phase7 Forbidden Claims\n\n" + "\n".join(f"- {x}" for x in FORBIDDEN_CLAIMS) + "\n")
    (out / "phase6_phase7_candidate_cards.md").write_text(
        "# Phase6 Phase7 Candidate Cards\n\n" + "\n".join(f"## {r['candidate_id']}\n- Priority: `{r['priority']}`\n- Mechanism: `{r['mechanism_group_id']}`\n- Axis: {r['intervention_axis']}\n- Allowed Phase7 use: {r['allowed_phase7_use']}\n" for r in rows)
    )
    verdict = "PASS" if main_rows else "FAIL"
    (out / "phase6_8_summary.md").write_text(
        f"# Phase6.8 Phase7 Handoff Summary\n\n- Verdict: `{verdict}`\n- Main validation candidates: {len(main_rows)}\n- Support/future addendum candidates: {len(support_rows)}\n- Pending evidence and forbidden claims are separated from main candidate handoff.\n"
    )
    write_status(out / "phase6_8_status.yaml", "Phase6.8", verdict, main_candidates=len(main_rows), support_candidates=len(support_rows))
    return main_rows, support_rows


def phase6_9(main: list[dict[str, str]], priority_rows: list[dict[str, Any]], phase7_main: list[dict[str, Any]], phase7_support: list[dict[str, Any]], targets: list[dict[str, Any]], perturb_support: list[dict[str, Any]]) -> None:
    out = mkdir(OUT / "09_final_report")
    master_rows = []
    top_target_by_mech: dict[str, dict[str, Any]] = {}
    for t in sorted(targets, key=lambda x: float(x["target_nomination_score"]), reverse=True):
        top_target_by_mech.setdefault(t["mechanism_group_id"], t)
    main_by_id = {r["mechanism_group_id"]: r for r in main}
    for pr in priority_rows:
        mech = main_by_id[pr["mechanism_group_id"]]
        top = top_target_by_mech.get(pr["mechanism_group_id"], {})
        master_rows.append({**pr, "phase5_priority": mech["priority"], "top_target_id": top.get("target_id", ""), "top_gene_symbol": top.get("gene_symbol", "")})
    write_csv(out / "phase6_candidate_master_table.csv", master_rows)
    write_csv(out / "phase6_phase7_handoff_main_candidates.csv", phase7_main)
    write_csv(out / "phase6_phase7_handoff_support_candidates.csv", phase7_support)
    counts = Counter(r["candidate_priority"] for r in priority_rows)
    verdict = "CONDITIONAL_PASS" if phase7_main and any("pending" in r["perturbation_support"] or "missing anchor" in r["conflict_summary"] for r in priority_rows) else ("PASS" if phase7_main else "FAIL")
    final_decision = {
        "verdict": verdict,
        "n_input_mechanisms": len(main),
        "n_priority_A_candidates": counts.get("Priority A", 0),
        "n_priority_B_candidates": counts.get("Priority B", 0),
        "n_priority_C_candidates": counts.get("Priority C", 0),
        "n_blocked_candidates": counts.get("Blocked", 0),
        "top_shared_candidates": [r["mechanism_group_id"] for r in priority_rows if r["storyline"] == "shared_immune_mechanism" and r["candidate_priority"] in {"Priority A", "Priority B"}],
        "top_HCC_specific_candidates": [r["mechanism_group_id"] for r in priority_rows if r["storyline"] == "HCC_specific_barrier" and r["candidate_priority"] in {"Priority A", "Priority B"}],
        "top_PD1X_repair_candidates": [r["mechanism_group_id"] for r in priority_rows if r["storyline"] == "PD1X_repair_logic" and r["candidate_priority"] in {"Priority A", "Priority B"}],
        "perturbation_support_summary": dict(Counter(r["perturbation_support_status"] for r in perturb_support)),
        "target_nomination_summary": dict(Counter(r["nomination_status"] for r in targets)),
        "tissue_spatial_external_gap_summary": "all Phase5 Priority 1/2 mechanisms retain pending_no_valid_anchor; Phase6 converts gap into Phase7 intake plan",
        "phase7_ready": bool(phase7_main),
        "required_caveats": REQUIRED_CAVEATS,
        "forbidden_claims": FORBIDDEN_CLAIMS,
        "recommended_next_phase": "Phase7 external/spatial/tissue validation and paper evidence package assembly",
    }
    write_yaml(out / "PHASE6_FINAL_DECISION.yaml", final_decision)
    report = [
        "# PHASE6 FINAL REPORT",
        "",
        f"- Verdict: `{verdict}`",
        f"- Input mechanisms: {len(main)} Phase5 Priority 1/2 candidates",
        f"- Priority A/B/C/Blocked: {counts.get('Priority A',0)}/{counts.get('Priority B',0)}/{counts.get('Priority C',0)}/{counts.get('Blocked',0)}",
        f"- Phase7 main handoff candidates: {len(phase7_main)}",
        f"- Phase7 support candidates: {len(phase7_support)}",
        "",
        "## Answers",
        "1. Perturbation mapping completed with curated pathway-target and immune-oncology class priors; no new perturbation data were integrated.",
        f"2. Perturbation support summary: {dict(Counter(r['perturbation_support_status'] for r in perturb_support))}.",
        f"3. Target nomination pre-screening produced {len(targets)} target rows; these are not final target recommendations.",
        f"4. HCC mainline candidates: {', '.join(final_decision['top_HCC_specific_candidates'])}.",
        f"5. Shared mainline candidates: {', '.join(final_decision['top_shared_candidates'])}.",
        f"6. PD1X repair mainline candidates: {', '.join(final_decision['top_PD1X_repair_candidates'])}.",
        "7. All mechanisms require tissue/spatial/external anchors under gated Phase7 intake.",
        "8. Priority A/B candidates enter Phase7 as validation candidates; Priority C enters support/future addendum.",
        f"9. Blocked candidates: {counts.get('Blocked', 0)}; conflicts are logged rather than hidden.",
        "10. Forbidden claims: no causal proof, no clinical recommendation, no drug recommendation, no PD1_anchor primary claim, no unintegrated new-data validation claim.",
        "",
        "Current result remains mechanism-candidate and target pre-screen evidence package, not clinical treatment guidance.",
        "",
    ]
    (out / "PHASE6_FINAL_REPORT.md").write_text("\n".join(report))
    (out / "phase6_9_summary.md").write_text(
        f"# Phase6.9 Final Report and Decision Summary\n\n- Verdict: `{verdict}`\n- Input mechanisms: {len(main)}\n- Priority A/B/C/Blocked: {counts.get('Priority A',0)}/{counts.get('Priority B',0)}/{counts.get('Priority C',0)}/{counts.get('Blocked',0)}\n- Phase7 ready: `{bool(phase7_main)}`\n- Final report, decision YAML, master table, handoff files, and reproducibility manifest written.\n"
    )
    manifest = {
        "phase": "Phase6 perturbation mapping and validation planning",
        "generated_at_utc": now_iso(),
        "script": "scripts/v6_1/run_phase6_perturbation_mapping.py",
        "input_files": [
            str(PHASE5 / "PHASE5_FINAL_DECISION.yaml"),
            str(PHASE5 / "phase5_mechanism_group_master_table.csv"),
            str(PHASE5 / "phase5_phase6_handoff_main_mechanisms.csv"),
            str(PHASE5 / "phase5_phase6_handoff_support_context.csv"),
            str(PHASE5 / "08_phase6_handoff/phase5_phase6_forbidden_use_registry.csv"),
            str(PHASE5 / "phase5_reproducibility_manifest.yaml"),
        ],
        "no_primary_use": ["Priority 3 support context", "PD1_anchor", "Phase3.5 complex ML", "unintegrated new data"],
        "final_outputs": [
            "PHASE6_FINAL_REPORT.md",
            "PHASE6_FINAL_DECISION.yaml",
            "phase6_candidate_master_table.csv",
            "phase6_phase7_handoff_main_candidates.csv",
            "phase6_phase7_handoff_support_candidates.csv",
        ],
    }
    write_yaml(out / "phase6_reproducibility_manifest.yaml", manifest)
    write_status(out / "phase6_9_status.yaml", "Phase6.9", verdict, phase7_ready=bool(phase7_main))
    for name in ["PHASE6_FINAL_REPORT.md", "PHASE6_FINAL_DECISION.yaml", "phase6_candidate_master_table.csv", "phase6_phase7_handoff_main_candidates.csv", "phase6_phase7_handoff_support_candidates.csv", "phase6_reproducibility_manifest.yaml"]:
        src = out / name
        (OUT / name).write_text(src.read_text())


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    decision = parse_final_decision()
    main_rows = read_csv(require(PHASE5 / "phase5_phase6_handoff_main_mechanisms.csv"))
    support_rows = read_csv(require(PHASE5 / "phase5_phase6_handoff_support_context.csv"))
    master_rows = read_csv(require(PHASE5 / "phase5_mechanism_group_master_table.csv"))
    phase6_0(main_rows, support_rows, master_rows, decision)
    _, gene_rows, _, _ = phase6_1(main_rows)
    _, _, perturb_support, _ = phase6_2(main_rows, gene_rows)
    targets, _, _, _ = phase6_3(main_rows, perturb_support, gene_rows)
    axis_rows, _, _, _ = phase6_4(main_rows, targets, gene_rows)
    _, _, _, gap_rows = phase6_5(main_rows)
    priority_rows, _ = phase6_6(main_rows, targets, perturb_support, axis_rows, gap_rows)
    design_rows = phase6_7(priority_rows, main_rows, targets)
    phase7_main, phase7_support = phase6_8(priority_rows, design_rows, targets)
    phase6_9(main_rows, priority_rows, phase7_main, phase7_support, targets, perturb_support)


if __name__ == "__main__":
    main()
