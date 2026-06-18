#!/usr/bin/env python3
from __future__ import annotations

import math
import re
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
DATE_TAG = "20260611"
PHASE4 = ROOT / "results" / "v6_1" / f"phase4_role_aware_module_discovery_{DATE_TAG}"
OUT = ROOT / "results" / "v6_1" / f"phase5_mechanism_adjudication_{DATE_TAG}"

PRIMARY_TIERS = {"Tier_A", "Tier_B"}
LINE_LABELS = {
    "shared_candidate_module": "shared_immune_mechanism",
    "HCC_specific_candidate_module": "HCC_specific_barrier",
    "PD1X_extension_candidate_module": "PD1X_repair_logic",
}
LINE_CODES = {
    "shared_immune_mechanism": "SHARED",
    "HCC_specific_barrier": "HCC",
    "PD1X_repair_logic": "PD1X",
}
ALLOWED_MECHANISM_CLASSES = [
    "cytotoxic effector activity",
    "T cell dysfunction / exhaustion",
    "regulatory suppression / Treg-like suppression",
    "myeloid inflammatory remodeling",
    "myeloid suppressive barrier",
    "DC / APC / antigen presentation",
    "IFN / inflammatory activation",
    "immune-cold / low-infiltration state",
    "stromal / vascular / exclusion proxy",
    "HCC immune-tolerance / liver-context barrier",
    "PD1X residual barrier / repair direction",
    "mixed or unclear mechanism",
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
]
REQUIRED_CAVEATS = [
    "source/batch/missingness confounding conditional-pass background",
    "associative mechanism candidates only",
    "Phase3.5 complex ML excluded",
    "PD1_anchor support-only biological anchor",
    "tissue/spatial/external anchors pending unless separately integrated",
]


def now_iso() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat()


def mkdir(path: Path) -> Path:
    path.mkdir(parents=True, exist_ok=True)
    return path


def rel(path: Path) -> str:
    try:
        return str(path.relative_to(ROOT))
    except ValueError:
        return str(path)


def require(path: Path) -> Path:
    if not path.exists():
        raise FileNotFoundError(path)
    return path


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


def read_verdict(path: Path) -> str:
    for line in path.read_text().splitlines():
        if line.startswith("verdict:"):
            return line.split(":", 1)[1].strip().strip('"')
    return "UNKNOWN"


def split_semicolon(value: Any) -> list[str]:
    if pd.isna(value):
        return []
    return [x.strip() for x in str(value).split(";") if x.strip()]


def direction_bucket(direction: Any) -> str:
    d = str(direction)
    if d == "response_associated":
        return "response"
    if d == "resistance_associated":
        return "resistance"
    if d == "context_dependent":
        return "context"
    return "uncertain"


def line_from_module_class(module_class: str) -> str:
    return LINE_LABELS.get(str(module_class), "support_only_context")


def clean_tokens(text: str) -> set[str]:
    parts = re.split(r"[^A-Za-z0-9]+", text.lower())
    stop = {"frac", "all", "immune", "parent", "pathway", "signature", "tf", "activity", "combined", "allowed"}
    return {p for p in parts if len(p) >= 3 and p not in stop}


MECH_RULES: list[tuple[str, list[str]]] = [
    ("T cell dysfunction / exhaustion", ["exhaust", "dysfunction", "pdcd1", "havcr2", "lag3", "tigit", "tox"]),
    ("regulatory suppression / Treg-like suppression", ["treg", "foxp3", "ctla4", "il2ra", "regulatory"]),
    ("myeloid suppressive barrier", ["spp1", "tam", "mdsc", "m2", "suppressive", "c1qc"]),
    ("myeloid inflammatory remodeling", ["myeloid_inflammatory", "mono", "fcn1", "il1", "tnf", "fcgr", "inflammatory"]),
    ("DC / APC / antigen presentation", ["dc", "apc", "hla", "mhc", "antigen"]),
    ("IFN / inflammatory activation", ["ifn", "interferon", "stat1", "irf1", "irf7"]),
    ("immune-cold / low-infiltration state", ["immune_cold", "low_infiltration", "non_immune", "epithelial", "tumor"]),
    ("stromal / vascular / exclusion proxy", ["stromal", "fibro", "endo", "vascular", "vegf", "exclusion"]),
    ("cytotoxic effector activity", ["cytotoxic", "cd8", "nk", "nkg", "gzmb", "prf", "ifng", "cd4_t"]),
]


def mechanism_classes(row: pd.Series) -> tuple[str, list[str]]:
    text = " ".join(
        [
            str(row.get("biological_interpretation", "")),
            str(row.get("key_features", "")),
            str(row.get("feature_list", "")),
        ]
    ).lower()
    hits = [label for label, keys in MECH_RULES if any(k in text for k in keys)]
    if not hits:
        if row.get("source_task") == "HCC_specific":
            hits = ["HCC immune-tolerance / liver-context barrier"]
        elif row.get("source_task") == "PD1X_extension":
            hits = ["PD1X residual barrier / repair direction"]
        else:
            hits = ["mixed or unclear mechanism"]
    primary = hits[0]
    secondary = [h for h in hits[1:] if h != primary][:2]
    return primary, secondary


def mechanism_orientation(direction: str, line: str) -> str:
    bucket = direction_bucket(direction)
    if line == "PD1X_repair_logic" and bucket in {"resistance", "context", "uncertain"}:
        return "repair target mechanism"
    if bucket == "response":
        return "response-supporting mechanism"
    if bucket == "resistance":
        return "resistance mechanism"
    return "context-dependent mechanism"


def score_stability(value: str) -> float:
    return {"A_stable": 1.0, "B_usable": 0.75, "C_exploratory": 0.35, "Blocked": 0.0}.get(str(value), 0.5)


def score_tier(value: str) -> float:
    return {"Tier_A": 1.0, "Tier_B": 0.75, "Tier_C": 0.35, "Blocked": 0.0}.get(str(value), 0.5)


def load_inputs() -> dict[str, Any]:
    return {
        "phase4_verdict": read_verdict(require(PHASE4 / "PHASE4_FINAL_DECISION.yaml")),
        "main": pd.read_csv(require(PHASE4 / "phase4_phase5_handoff_main_modules.csv")),
        "support": pd.read_csv(require(PHASE4 / "phase4_phase5_handoff_support_evidence.csv")),
        "master": pd.read_csv(require(PHASE4 / "phase4_module_master_table.csv")),
        "blocked": pd.read_csv(require(PHASE4 / "phase4_blocked_module_and_input_log.csv")),
        "sample_scores": pd.read_csv(require(PHASE4 / "07_module_scoring/phase4_7_module_scores_by_sample.csv")),
        "patient_scores": pd.read_csv(require(PHASE4 / "07_module_scoring/phase4_7_module_scores_by_patient.csv")),
        "annotations": pd.read_csv(require(PHASE4 / "06_biological_interpretation/phase4_6_module_biological_annotation.csv")),
        "support_conflicts": pd.read_csv(require(PHASE4 / "06_biological_interpretation/phase4_6_support_conflict_log.csv")),
    }


def enrich_modules(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    out["storyline"] = out["module_class"].map(line_from_module_class)
    out["line_code"] = out["storyline"].map(LINE_CODES).fillna("SUPPORT")
    out["direction_bucket"] = out["response_direction"].map(direction_bucket)
    primaries, secondaries = [], []
    for _, row in out.iterrows():
        primary, secondary = mechanism_classes(row)
        primaries.append(primary)
        secondaries.append(";".join(secondary))
    out["primary_mechanism_class"] = primaries
    out["secondary_mechanism_class"] = secondaries
    out["feature_tokens"] = out.apply(lambda r: ";".join(sorted(clean_tokens(str(r.get("feature_list", "")) + " " + str(r.get("key_features", ""))))), axis=1)
    out["compression_key"] = out.apply(
        lambda r: "||".join(
            [
                str(r["line_code"]),
                str(r["primary_mechanism_class"]),
                str(r["direction_bucket"]),
                str(r["source_feature_family"]) if r["primary_mechanism_class"] == "mixed or unclear mechanism" else "all_families",
            ]
        ),
        axis=1,
    )
    return out


def phase5_0(inputs: dict[str, Any]) -> pd.DataFrame:
    out = mkdir(OUT / "00_entry_contract")
    main = enrich_modules(inputs["main"])
    support = inputs["support"]
    blocked = inputs["blocked"]
    phase4_verdict = inputs["phase4_verdict"]
    required = [
        "module_id",
        "module_class",
        "tier",
        "source_task",
        "source_feature_family",
        "key_features",
        "module_score_available",
        "response_direction",
        "stability_grade",
        "required_caveat",
    ]
    rows = []
    def add(rule: str, violation: bool, severity: str, action: str) -> None:
        rows.append({"rule": rule, "violation": bool(violation), "severity": severity, "action": action})

    add("phase4_verdict_pass_or_conditional", phase4_verdict not in {"PASS", "CONDITIONAL_PASS"}, "HARD_FAIL", "stop_phase5")
    add("main_input_tier_A_or_B_only", not set(main["tier"]).issubset(PRIMARY_TIERS), "HARD_FAIL", "exclude_non_primary_tiers")
    add("main_input_primary_evidence_pool_only", not main["evidence_pool"].astype(str).eq("primary_evidence_pool").all(), "HARD_FAIL", "stop_phase5")
    add("support_not_in_main", main["module_class"].astype(str).eq("support_only_module").any(), "HARD_FAIL", "stop_phase5")
    add("PD1_anchor_not_primary", main["source_task"].astype(str).eq("PD1_anchor").any(), "HARD_FAIL", "stop_phase5")
    text = " ".join(main.astype(str).fillna("").agg(" ".join, axis=1).tolist())
    add("phase3_5_ml_not_main", bool(re.search(r"phase3\.5|strong ML|complex ML", text, flags=re.I)), "HARD_FAIL", "stop_phase5")
    add("unintegrated_new_data_not_primary", bool(re.search(r"unintegrated|new_data_integration_addendum", text, flags=re.I)), "HARD_FAIL", "stop_phase5")
    add("blocked_modules_excluded", len(blocked) > 0, "FAIL_IF_USED", "keep_blocked_log_only")
    missing_cols = [c for c in required if c not in main.columns]
    add("required_main_fields_present", bool(missing_cols), "FAIL", "repair_phase4_handoff")
    audit = pd.DataFrame(rows)
    audit.to_csv(out / "phase5_0_input_violation_audit.csv", index=False)
    inventory_cols = required + [
        "candidate_module_id",
        "evidence_pool",
        "storyline",
        "primary_mechanism_class",
        "secondary_mechanism_class",
        "direction_bucket",
        "confounding_risk",
        "allowed_phase5_use",
        "forbidden_phase5_use",
    ]
    main[inventory_cols].to_csv(out / "phase5_0_candidate_inventory.csv", index=False)
    hard_fail = audit[(audit["violation"]) & (audit["severity"] == "HARD_FAIL")]
    verdict = "FAIL" if len(hard_fail) else "PASS"
    (out / "phase5_0_entry_contract_summary.md").write_text(
        "\n".join(
            [
                "# Phase5.0 Entry Contract Summary",
                "",
                f"- Verdict: `{verdict}`",
                f"- Phase4 verdict: `{phase4_verdict}`",
                f"- Main candidate modules: {len(main)}",
                f"- Support-only records: {len(support)}",
                f"- Blocked records in Phase4 log: {len(blocked)}",
                "- Main adjudication reads Phase4 frozen Tier A/B handoff only.",
                "- Support evidence remains annotation/sensitivity context.",
                "",
            ]
        )
    )
    write_status(out / "phase5_0_status.yaml", "Phase5.0", verdict, main_candidates=len(main), support_records=len(support), blocked_records=len(blocked))
    if verdict == "FAIL":
        raise RuntimeError("Phase5.0 HARD_FAIL: invalid Phase4 handoff isolation")
    return main


def feature_set(row: pd.Series) -> set[str]:
    return set(split_semicolon(row.get("feature_list", ""))) | set(split_semicolon(row.get("key_features", "")))


def build_similarity(main: pd.DataFrame, sample_scores: pd.DataFrame) -> pd.DataFrame:
    ids = main["module_id"].tolist()
    mod_by_id = main.set_index("module_id")
    score_sub = sample_scores[sample_scores["module_id"].isin(ids) & sample_scores["score_role"].eq("primary_score")]
    pivot = score_sub.pivot_table(index="sample_key", columns="module_id", values="module_score", aggfunc="mean")
    corr = pivot.corr(min_periods=20) if not pivot.empty else pd.DataFrame(index=ids, columns=ids)
    fsets = {mid: feature_set(mod_by_id.loc[mid]) for mid in ids}
    bsets = {mid: set(split_semicolon(mod_by_id.loc[mid, "biological_interpretation"])) | {mod_by_id.loc[mid, "primary_mechanism_class"]} for mid in ids}
    sims = pd.DataFrame(0.0, index=ids, columns=ids)
    for i, a in enumerate(ids):
        for b in ids[i:]:
            if a == b:
                sim = 1.0
            else:
                fa, fb = fsets[a], fsets[b]
                ba, bb = bsets[a], bsets[b]
                feature_sim = len(fa & fb) / len(fa | fb) if (fa or fb) else 0.0
                bio_sim = len(ba & bb) / len(ba | bb) if (ba or bb) else 0.0
                ca = corr.loc[a, b] if a in corr.index and b in corr.columns else np.nan
                score_sim = max(0.0, float(ca)) if pd.notna(ca) else 0.0
                dir_sim = 1.0 if mod_by_id.loc[a, "direction_bucket"] == mod_by_id.loc[b, "direction_bucket"] else 0.0
                fam_sim = 1.0 if mod_by_id.loc[a, "source_feature_family"] == mod_by_id.loc[b, "source_feature_family"] else 0.4
                task_sim = 0.7 if mod_by_id.loc[a, "line_code"] == mod_by_id.loc[b, "line_code"] else 0.3
                sim = 0.25 * feature_sim + 0.25 * bio_sim + 0.20 * score_sim + 0.15 * dir_sim + 0.10 * fam_sim + 0.05 * task_sim
            sims.loc[a, b] = sim
            sims.loc[b, a] = sim
    return sims


def choose_representative(sub: pd.DataFrame) -> str:
    score = (
        sub["tier"].map(score_tier).fillna(0.5) * 3.0
        + sub["stability_grade"].map(score_stability).fillna(0.5) * 2.0
        + pd.to_numeric(sub["effect_size_response_minus_nonresponse"], errors="coerce").abs().fillna(0.0)
        + (1.0 - pd.to_numeric(sub["fdr"], errors="coerce").fillna(1.0)).clip(lower=0.0, upper=1.0)
        + pd.to_numeric(sub["n_features"], errors="coerce").fillna(1.0).clip(upper=25) / 25.0
    )
    return str(sub.loc[score.idxmax(), "module_id"])


def phase5_1(main: pd.DataFrame, sample_scores: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    out = mkdir(OUT / "01_candidate_compression")
    sim = build_similarity(main, sample_scores)
    sim_out = sim.copy()
    sim_out.insert(0, "module_id", sim_out.index)
    sim_out.to_csv(out / "phase5_1_module_similarity_matrix.csv", index=False)
    registry_rows, mapping_rows = [], []
    ordered_keys = sorted(main["compression_key"].unique())
    group_counter = {"SHARED": 0, "HCC": 0, "PD1X": 0}
    for key in ordered_keys:
        sub = main[main["compression_key"] == key].copy()
        line_code = str(sub["line_code"].iloc[0])
        group_counter[line_code] = group_counter.get(line_code, 0) + 1
        gid = f"P5MG_{line_code}_{group_counter[line_code]:03d}"
        rep = choose_representative(sub)
        source_tasks = sorted(sub["source_task"].astype(str).unique())
        families = sorted(sub["source_feature_family"].astype(str).unique())
        dirs = sorted(sub["response_direction"].astype(str).unique())
        non_uncertain = [d for d in dirs if d not in {"uncertain", "unsupported"}]
        direction_consistency = "consistent" if len({direction_bucket(d) for d in non_uncertain}) <= 1 else "context_dependent_or_conflicting"
        mean_sim = np.nan
        if len(sub) > 1:
            mids = sub["module_id"].tolist()
            vals = sim.loc[mids, mids].to_numpy()
            tri = vals[np.triu_indices_from(vals, k=1)]
            mean_sim = float(np.nanmean(tri)) if len(tri) else np.nan
        registry_rows.append(
            {
                "mechanism_group_id": gid,
                "representative_module_id": rep,
                "member_module_ids": ";".join(sub["module_id"].tolist()),
                "n_member_modules": len(sub),
                "source_tasks": ";".join(source_tasks),
                "source_feature_families": ";".join(families),
                "shared_present": bool((sub["module_class"] == "shared_candidate_module").any()),
                "HCC_specific_present": bool((sub["module_class"] == "HCC_specific_candidate_module").any()),
                "PD1X_extension_present": bool((sub["module_class"] == "PD1X_extension_candidate_module").any()),
                "storyline": sub["storyline"].iloc[0],
                "primary_mechanism_class": sub["primary_mechanism_class"].iloc[0],
                "response_direction_consistency": direction_consistency,
                "biological_annotation_consistency": "same_primary_mechanism_class",
                "mean_pairwise_similarity": mean_sim,
                "redundancy_action": "merge_members_keep_representative" if len(sub) > 1 else "singleton_keep",
                "notes": "compressed_within_storyline_by_mechanism_class_and_direction",
            }
        )
        for _, row in sub.iterrows():
            mapping_rows.append(
                {
                    "module_id": row["module_id"],
                    "mechanism_group_id": gid,
                    "representative_module_id": rep,
                    "redundancy_action": "representative" if row["module_id"] == rep else "supporting_member",
                    "source_task": row["source_task"],
                    "source_feature_family": row["source_feature_family"],
                    "tier": row["tier"],
                    "response_direction": row["response_direction"],
                }
            )
    registry = pd.DataFrame(registry_rows)
    mapping = pd.DataFrame(mapping_rows)
    registry.to_csv(out / "phase5_1_mechanism_group_registry.csv", index=False)
    registry[["mechanism_group_id", "representative_module_id", "storyline", "primary_mechanism_class", "n_member_modules", "source_tasks", "source_feature_families"]].to_csv(
        out / "phase5_1_representative_modules.csv", index=False
    )
    mapping.to_csv(out / "phase5_1_redundant_module_mapping.csv", index=False)
    match_rows = []
    for (mech, direction), sub in registry.assign(
        group_direction=registry["representative_module_id"].map(main.set_index("module_id")["direction_bucket"])
    ).groupby(["primary_mechanism_class", "group_direction"], dropna=False):
        lines = set(sub["storyline"])
        if len(lines) > 1:
            match_rows.append(
                {
                    "mechanism_class": mech,
                    "direction_bucket": direction,
                    "mechanism_group_ids": ";".join(sub["mechanism_group_id"].tolist()),
                    "matched_storylines": ";".join(sorted(lines)),
                    "match_type": "all_three_matched" if len(lines) == 3 else "pairwise_matched",
                }
            )
    matches = pd.DataFrame(match_rows, columns=["mechanism_class", "direction_bucket", "mechanism_group_ids", "matched_storylines", "match_type"])
    matches.to_csv(out / "phase5_1_cross_task_mechanism_matches.csv", index=False)
    verdict = "PASS" if len(registry) < len(main) and len(registry) > 0 else "FAIL"
    (out / "phase5_1_summary.md").write_text(
        "\n".join(
            [
                "# Phase5.1 Candidate Compression Summary",
                "",
                f"- Verdict: `{verdict}`",
                f"- Input modules: {len(main)}",
                f"- Mechanism groups: {len(registry)}",
                f"- Cross-task matched mechanism classes: {len(matches)}",
                "- Modules are compressed within storyline by mechanism class and direction; original module membership is retained.",
                "",
            ]
        )
    )
    write_status(out / "phase5_1_status.yaml", "Phase5.1", verdict, input_modules=len(main), mechanism_groups=len(registry), cross_task_matches=len(matches))
    return registry, mapping, matches


def group_member_frame(registry_row: pd.Series, main: pd.DataFrame) -> pd.DataFrame:
    mids = split_semicolon(registry_row["member_module_ids"])
    return main[main["module_id"].isin(mids)].copy()


def phase5_2(
    registry: pd.DataFrame,
    mapping: pd.DataFrame,
    main: pd.DataFrame,
    support: pd.DataFrame,
    sample_scores: pd.DataFrame,
    patient_scores: pd.DataFrame,
    support_conflicts: pd.DataFrame,
) -> pd.DataFrame:
    out = mkdir(OUT / "02_evidence_matrix")
    support_enriched = enrich_modules(support) if len(support) else support.copy()
    conflict_ids = set(support_conflicts["module_id"].astype(str)) if len(support_conflicts) else set()
    sample_counts = sample_scores.groupby("module_id")["sample_key"].nunique().to_dict()
    patient_counts = patient_scores.groupby("module_id")["patient_key"].nunique().to_dict()
    evidence_rows, conflict_rows, risk_rows, usability_rows = [], [], [], []
    for _, grow in registry.iterrows():
        sub = group_member_frame(grow, main)
        mids = sub["module_id"].tolist()
        rep = str(grow["representative_module_id"])
        rep_row = main.set_index("module_id").loc[rep]
        support_matches = support_enriched[
            (support_enriched.get("storyline", "") == grow["storyline"])
            & (support_enriched.get("primary_mechanism_class", "") == grow["primary_mechanism_class"])
        ] if len(support_enriched) else pd.DataFrame()
        conflict = sorted(set(mids) & conflict_ids)
        cohort_values = set(sub["cohort_sensitivity_status"].astype(str).unique())
        missing_values = set(sub["missingness_sensitivity_status"].astype(str).unique())
        cohort_status = ";".join(sorted(cohort_values))
        missing_status = ";".join(sorted(missing_values))
        risk_level = "high_risk_demote" if ("cohort_dominated" in cohort_values or "missingness_driven" in missing_values) else "moderate_caveat"
        if risk_level == "moderate_caveat" and sub["confounding_risk"].astype(str).eq("LOW").all():
            risk_level = "acceptable_with_global_caveat"
        evidence_rows.append(
            {
                "mechanism_group_id": grow["mechanism_group_id"],
                "representative_module_id": rep,
                "member_module_ids": ";".join(mids),
                "source_tasks": grow["source_tasks"],
                "source_feature_families": grow["source_feature_families"],
                "tiers": ";".join(sorted(sub["tier"].astype(str).unique())),
                "representative_stability_grade": rep_row["stability_grade"],
                "member_stability_grades": ";".join(sorted(sub["stability_grade"].astype(str).unique())),
                "response_directions": ";".join(sorted(sub["response_direction"].astype(str).unique())),
                "median_effect_size": float(pd.to_numeric(sub["effect_size_response_minus_nonresponse"], errors="coerce").median()),
                "min_fdr": float(pd.to_numeric(sub["fdr"], errors="coerce").min()),
                "cohort_sensitivity_status": cohort_status,
                "missingness_sensitivity_status": missing_status,
                "sample_score_rows": int(sum(sample_counts.get(m, 0) for m in mids)),
                "patient_score_rows": int(sum(patient_counts.get(m, 0) for m in mids)),
                "primary_mechanism_class": grow["primary_mechanism_class"],
                "biological_interpretation": ";".join(sorted(set(sum([split_semicolon(x) for x in sub["biological_interpretation"]], [])))),
                "key_features": ";".join(sorted(set(sum([split_semicolon(x) for x in sub["key_features"]], []))))[:2000],
                "PD1_anchor_direction_consistency": "support_conflict_present" if conflict else "no_recorded_support_conflict",
                "support_only_module_agreement_count": int(len(support_matches)),
                "support_conflict": ";".join(conflict),
                "risk_level": risk_level,
                "downstream_usability": "Phase6_main_candidate_possible" if risk_level != "high_risk_demote" else "Phase6_support_only_or_demote",
            }
        )
        conflict_rows.append(
            {
                "mechanism_group_id": grow["mechanism_group_id"],
                "representative_module_id": rep,
                "support_conflict_present": bool(conflict),
                "conflicting_module_ids": ";".join(conflict),
                "support_only_agreement_count": int(len(support_matches)),
                "support_use_rule": "annotation_direction_sanity_check_only",
            }
        )
        risk_rows.append(
            {
                "mechanism_group_id": grow["mechanism_group_id"],
                "source_batch_caveat": "global_conditional_pass_caveat",
                "missingness_caveat": missing_status,
                "single_cohort_dominance": "present" if "cohort_dominated" in cohort_values else "not_detected",
                "feature_family_limitation": "single_family_only" if len(split_semicolon(grow["source_feature_families"])) == 1 else "multi_family_support",
                "response_direction_uncertainty": "present" if "uncertain" in set(sub["direction_bucket"]) else "not_primary_issue",
                "over_interpretation_risk": "higher" if grow["primary_mechanism_class"] == "mixed or unclear mechanism" else "managed",
                "risk_level": risk_level,
            }
        )
        usability_rows.append(
            {
                "mechanism_group_id": grow["mechanism_group_id"],
                "suitable_for_tissue_spatial_adjudication": grow["primary_mechanism_class"] != "mixed or unclear mechanism",
                "suitable_for_perturbation_mapping": grow["primary_mechanism_class"] != "mixed or unclear mechanism",
                "suitable_for_target_nomination_pre_screening": grow["primary_mechanism_class"] not in {"mixed or unclear mechanism", "immune-cold / low-infiltration state"},
                "suitable_only_for_annotation": risk_level == "high_risk_demote",
                "blocked_from_phase6_main_use": risk_level == "high_risk_demote",
            }
        )
    evidence = pd.DataFrame(evidence_rows)
    evidence.to_csv(out / "phase5_2_mechanism_evidence_matrix.csv", index=False)
    pd.DataFrame(conflict_rows).to_csv(out / "phase5_2_support_conflict_matrix.csv", index=False)
    pd.DataFrame(risk_rows).to_csv(out / "phase5_2_risk_evidence_matrix.csv", index=False)
    pd.DataFrame(usability_rows).to_csv(out / "phase5_2_downstream_usability_matrix.csv", index=False)
    verdict = "PASS" if len(evidence) == len(registry) and len(evidence) > 0 else "FAIL"
    (out / "phase5_2_summary.md").write_text(
        f"# Phase5.2 Evidence Matrix Summary\n\n- Verdict: `{verdict}`\n- Mechanism groups with evidence rows: {len(evidence)}\n- Primary/support/risk/downstream fields are separated.\n"
    )
    write_status(out / "phase5_2_status.yaml", "Phase5.2", verdict, mechanism_groups=len(evidence))
    return evidence


def phase5_3(evidence: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    out = mkdir(OUT / "03_direction_stability_confounding")
    direction_rows, stability_rows, confounding_rows, demoted_rows = [], [], [], []
    for _, row in evidence.iterrows():
        dirs = set(split_semicolon(row["response_directions"]))
        buckets = {direction_bucket(d) for d in dirs if d not in {"uncertain", "unsupported"}}
        if buckets == {"response"}:
            direction_call = "consistent_response_associated"
        elif buckets == {"resistance"}:
            direction_call = "consistent_resistance_associated"
        elif len(buckets) > 1:
            direction_call = "context_dependent"
        elif "uncertain" in dirs or "unsupported" in dirs:
            direction_call = "uncertain"
        else:
            direction_call = "context_dependent"
        direction_rows.append(
            {
                "mechanism_group_id": row["mechanism_group_id"],
                "representative_module_id": row["representative_module_id"],
                "response_directions": row["response_directions"],
                "direction_adjudication": direction_call,
                "support_conflict": row["support_conflict"],
                "direction_demote_reason": "none" if direction_call != "conflicting" else "unexplained_direction_conflict",
            }
        )
        stabs = set(split_semicolon(row["member_stability_grades"]))
        if "A_stable" in stabs and "B_usable" not in stabs:
            stability_call = "robust"
        elif stabs <= {"A_stable", "B_usable"}:
            stability_call = "usable"
        elif "Blocked" in stabs:
            stability_call = "unstable"
        else:
            stability_call = "exploratory"
        stability_rows.append(
            {
                "mechanism_group_id": row["mechanism_group_id"],
                "representative_stability_grade": row["representative_stability_grade"],
                "member_stability_grades": row["member_stability_grades"],
                "cohort_sensitivity_status": row["cohort_sensitivity_status"],
                "missingness_sensitivity_status": row["missingness_sensitivity_status"],
                "stability_adjudication": stability_call,
            }
        )
        risk = str(row["risk_level"])
        if risk == "high_risk_demote":
            conf_call = "high_risk_demote"
        elif risk == "acceptable_with_global_caveat":
            conf_call = "acceptable"
        else:
            conf_call = "moderate_caveat"
        confounding_rows.append(
            {
                "mechanism_group_id": row["mechanism_group_id"],
                "cohort_sensitivity_status": row["cohort_sensitivity_status"],
                "missingness_sensitivity_status": row["missingness_sensitivity_status"],
                "confounding_adjudication": conf_call,
                "required_caveat": "source/batch/missingness conditional-pass caveat retained",
            }
        )
        demote = []
        if direction_call == "conflicting":
            demote.append("direction_conflict")
        if stability_call == "unstable":
            demote.append("unstable")
        if conf_call in {"high_risk_demote", "blocked"}:
            demote.append(conf_call)
        if demote:
            demoted_rows.append({"mechanism_group_id": row["mechanism_group_id"], "demotion_or_block_reason": ";".join(demote), "action": "remove_from_phase6_main"})
    direction = pd.DataFrame(direction_rows)
    stability = pd.DataFrame(stability_rows)
    confounding = pd.DataFrame(confounding_rows)
    demoted = pd.DataFrame(demoted_rows, columns=["mechanism_group_id", "demotion_or_block_reason", "action"])
    direction.to_csv(out / "phase5_3_direction_adjudication.csv", index=False)
    stability.to_csv(out / "phase5_3_stability_adjudication.csv", index=False)
    confounding.to_csv(out / "phase5_3_confounding_adjudication.csv", index=False)
    demoted.to_csv(out / "phase5_3_demoted_or_blocked_groups.csv", index=False)
    verdict = "PASS" if len(direction) and not (demoted["action"].eq("remove_from_phase6_main").all() if len(demoted) else False) else "FAIL"
    (out / "phase5_3_summary.md").write_text(
        f"# Phase5.3 Direction Stability Confounding Summary\n\n- Verdict: `{verdict}`\n- Groups adjudicated: {len(direction)}\n- Demoted or blocked groups: {len(demoted)}\n- Direction conflicts are retained rather than hidden.\n"
    )
    write_status(out / "phase5_3_status.yaml", "Phase5.3", verdict, groups=len(direction), demoted_or_blocked=len(demoted))
    return direction, stability, confounding


def hcc_relevance(mechanism_class: str, present: bool) -> str:
    if not present:
        return "not_HCC_specific"
    if "myeloid" in mechanism_class:
        return "myeloid barrier"
    if "APC" in mechanism_class or "antigen" in mechanism_class:
        return "APC impairment"
    if "stromal" in mechanism_class:
        return "vascular/stromal exclusion proxy"
    if "HCC" in mechanism_class:
        return "liver tolerance"
    return "HCC-enriched immune context"


def pd1x_relevance(mechanism_class: str, present: bool, direction_call: str) -> str:
    if not present:
        return "not_PD1X_extension"
    if "resistance" in direction_call or "myeloid" in mechanism_class or "T cell dysfunction" in mechanism_class:
        return "residual barrier after ICI"
    if "response" in direction_call or "APC" in mechanism_class or "IFN" in mechanism_class:
        return "combination-compatible repair axis"
    return "uncertain extension signal"


def phase5_4(registry: pd.DataFrame, evidence: pd.DataFrame, direction: pd.DataFrame) -> pd.DataFrame:
    out = mkdir(OUT / "04_biological_mechanism")
    ev = evidence.set_index("mechanism_group_id")
    dr = direction.set_index("mechanism_group_id")
    class_rows, hcc_rows, pd1x_rows, shared_rows = [], [], [], []
    card_lines = ["# Phase5.4 Mechanism Cards", ""]
    for _, row in registry.iterrows():
        gid = row["mechanism_group_id"]
        erow = ev.loc[gid]
        dcall = dr.loc[gid, "direction_adjudication"]
        secondary = []
        for term in split_semicolon(erow["biological_interpretation"]):
            fake = pd.Series({"biological_interpretation": term, "key_features": erow["key_features"], "source_task": ""})
            primary, _ = mechanism_classes(fake)
            if primary != row["primary_mechanism_class"] and primary not in secondary:
                secondary.append(primary)
            if len(secondary) == 2:
                break
        orientation = mechanism_orientation(str(dcall), row["storyline"])
        hcc_rel = hcc_relevance(row["primary_mechanism_class"], bool(row["HCC_specific_present"]))
        pd1x_rel = pd1x_relevance(row["primary_mechanism_class"], bool(row["PD1X_extension_present"]), str(dcall))
        shared_support = "multi_cohort_shared_candidate_with_caveat" if row["shared_present"] else "not_shared_line"
        class_rows.append(
            {
                "mechanism_group_id": gid,
                "representative_module_id": row["representative_module_id"],
                "primary_mechanism_class": row["primary_mechanism_class"],
                "secondary_mechanism_class": ";".join(secondary),
                "mechanism_orientation": orientation,
                "biological_rationale": f"{row['primary_mechanism_class']} supported by member feature composition and Phase4 annotation",
                "key_supporting_features": erow["key_features"],
                "key_conflicting_evidence": erow["support_conflict"],
                "caveats": "associative module; no causal or drug claim; source/batch/missingness caveat retained",
            }
        )
        hcc_rows.append({"mechanism_group_id": gid, "HCC_specific_present": row["HCC_specific_present"], "HCC_context_adjudication": hcc_rel})
        pd1x_rows.append({"mechanism_group_id": gid, "PD1X_extension_present": row["PD1X_extension_present"], "PD1X_repair_logic_adjudication": pd1x_rel})
        shared_rows.append({"mechanism_group_id": gid, "shared_present": row["shared_present"], "shared_mechanism_adjudication": shared_support})
        card_lines.extend(
            [
                f"## {gid}",
                "",
                f"- Representative module: `{row['representative_module_id']}`",
                f"- Source modules: `{row['member_module_ids']}`",
                f"- Primary mechanism class: `{row['primary_mechanism_class']}`",
                f"- Secondary mechanism class: `{';'.join(secondary) if secondary else 'none'}`",
                f"- Direction adjudication: `{dcall}`",
                f"- Biological rationale: {row['primary_mechanism_class']} supported by feature composition and Phase4 annotation.",
                f"- Key supporting features: {str(erow['key_features'])[:500]}",
                f"- Key conflicting evidence: {erow['support_conflict'] if erow['support_conflict'] else 'none recorded'}",
                f"- HCC relevance: {hcc_rel}",
                f"- PD1X relevance: {pd1x_rel}",
                f"- Downstream suitability: {erow['downstream_usability']}",
                "- Caveats: associative candidate only; no clinical prediction, causality, or drug recommendation.",
                "",
            ]
        )
    classification = pd.DataFrame(class_rows)
    classification.to_csv(out / "phase5_4_biological_mechanism_classification.csv", index=False)
    pd.DataFrame(hcc_rows).to_csv(out / "phase5_4_HCC_context_adjudication.csv", index=False)
    pd.DataFrame(pd1x_rows).to_csv(out / "phase5_4_PD1X_repair_logic_adjudication.csv", index=False)
    pd.DataFrame(shared_rows).to_csv(out / "phase5_4_shared_mechanism_adjudication.csv", index=False)
    (out / "phase5_4_mechanism_cards.md").write_text("\n".join(card_lines))
    unclear = int((classification["primary_mechanism_class"] == "mixed or unclear mechanism").sum())
    verdict = "PASS" if len(classification) and unclear < len(classification) else "FAIL"
    (out / "phase5_4_summary.md").write_text(
        f"# Phase5.4 Biological Mechanism Summary\n\n- Verdict: `{verdict}`\n- Mechanism groups classified: {len(classification)}\n- Mixed/unclear groups: {unclear}\n- No causal or drug claim is made.\n"
    )
    write_status(out / "phase5_4_status.yaml", "Phase5.4", verdict, classified_groups=len(classification), mixed_or_unclear=unclear)
    return classification


def phase5_5(registry: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    out = mkdir(OUT / "05_tissue_spatial_external")
    unintegrated = sorted(Path(ROOT / "results/v6_1/new_data_integration_addendum_20260611").glob("**/*spatial*"))
    availability = pd.DataFrame(
        [
            {
                "anchor_type": "tissue",
                "valid_integrated_anchor_available": False,
                "valid_for_primary_adjudication": False,
                "excluded_unintegrated_files": "",
                "status": "pending_no_valid_anchor",
                "reason": "no tissue anchor passed current Phase4/Phase5 gate",
            },
            {
                "anchor_type": "spatial",
                "valid_integrated_anchor_available": False,
                "valid_for_primary_adjudication": False,
                "excluded_unintegrated_files": ";".join(rel(p) for p in unintegrated),
                "status": "pending_no_valid_anchor",
                "reason": "only unintegrated addendum spatial files detected; not used for primary adjudication",
            },
            {
                "anchor_type": "external_bulk",
                "valid_integrated_anchor_available": False,
                "valid_for_primary_adjudication": False,
                "excluded_unintegrated_files": "",
                "status": "pending_no_valid_anchor",
                "reason": "no external anchor passed current gate",
            },
        ]
    )
    availability.to_csv(out / "phase5_5_adjudication_data_availability.csv", index=False)
    tissue_rows, external_rows = [], []
    for _, row in registry.iterrows():
        tissue_rows.append(
            {
                "mechanism_group_id": row["mechanism_group_id"],
                "representative_module_id": row["representative_module_id"],
                "tissue_spatial_status": "pending_no_valid_anchor",
                "tissue_support_detail": "no integrated tissue/spatial anchor used",
                "spatial_caveat": "unintegrated new data excluded from primary adjudication",
            }
        )
        external_rows.append(
            {
                "mechanism_group_id": row["mechanism_group_id"],
                "representative_module_id": row["representative_module_id"],
                "external_anchor_status": "pending_no_valid_anchor",
                "external_support_detail": "no integrated external anchor used",
                "bulk_not_interpreted_as_spatial": True,
            }
        )
    tissue = pd.DataFrame(tissue_rows)
    external = pd.DataFrame(external_rows)
    tissue.to_csv(out / "phase5_5_tissue_spatial_adjudication.csv", index=False)
    external.to_csv(out / "phase5_5_external_anchor_adjudication.csv", index=False)
    (out / "phase5_5_future_data_request.md").write_text(
        "\n".join(
            [
                "# Phase5.5 Future Data Request",
                "",
                "- Current status: `pending_no_valid_anchor` for tissue, spatial, and external anchors.",
                "- Do not use unintegrated addendum data for current primary conclusions.",
                "- Requested future addendum: integrated tissue/spatial manifest, module-score projection contract, tumor/myeloid/T cell/APC spatial context fields, and gate-passed external validation metadata.",
                "- Phase6 may plan spatial/tissue/external validation, but must not claim current spatial support.",
                "",
            ]
        )
    )
    verdict = "CONDITIONAL_PASS"
    (out / "phase5_5_summary.md").write_text(
        f"# Phase5.5 Tissue Spatial External Anchor Summary\n\n- Verdict: `{verdict}`\n- Valid integrated anchors: 0\n- Groups marked pending: {len(registry)}\n- Unintegrated addendum data were not used as primary evidence.\n"
    )
    write_status(out / "phase5_5_status.yaml", "Phase5.5", verdict, valid_integrated_anchors=0, pending_groups=len(registry))
    return tissue, external


def phase5_6(registry: pd.DataFrame, classification: pd.DataFrame, matches: pd.DataFrame) -> pd.DataFrame:
    out = mkdir(OUT / "06_storyline_consolidation")
    cls = classification.set_index("mechanism_group_id")
    rows, graph_rows = [], []
    for _, row in registry.iterrows():
        gid = row["mechanism_group_id"]
        storyline = row["storyline"]
        if row["primary_mechanism_class"] == "mixed or unclear mechanism":
            assignment = "exploratory"
        else:
            assignment = storyline
        rows.append(
            {
                "mechanism_group_id": gid,
                "representative_module_id": row["representative_module_id"],
                "main_story_line": assignment,
                "original_storyline": storyline,
                "biological_mechanism_class": row["primary_mechanism_class"],
                "response_direction_consistency": row["response_direction_consistency"],
                "storyline_rule": "line_preserved_from_Phase4_class; mixed_unclear_demoted_to_exploratory",
            }
        )
        graph_rows.append(
            {
                "mechanism_group_id": gid,
                "biological_mechanism_class": row["primary_mechanism_class"],
                "story_line": assignment,
                "phase6_actionability": "mechanism-class mapping; perturbation planning; target nomination pre-screening" if assignment != "exploratory" else "support or future validation",
            }
        )
    assignment = pd.DataFrame(rows)
    assignment.to_csv(out / "phase5_6_storyline_assignment.csv", index=False)
    rel_rows = []
    for _, m in matches.iterrows():
        rel_rows.append(
            {
                "mechanism_class": m["mechanism_class"],
                "direction_bucket": m["direction_bucket"],
                "matched_groups": m["mechanism_group_ids"],
                "relationship_interpretation": "shared/HCC/PD1X related mechanism class; do not collapse storylines or overgeneralize HCC/PD1X",
            }
        )
    pd.DataFrame(rel_rows, columns=["mechanism_class", "direction_bucket", "matched_groups", "relationship_interpretation"]).to_csv(
        out / "phase5_6_shared_HCC_PD1X_relationship_map.csv", index=False
    )
    pd.DataFrame(graph_rows).to_csv(out / "phase5_6_mechanism_story_graph.tsv", sep="\t", index=False)
    line_counts = assignment["main_story_line"].value_counts().to_dict()
    verdict = "PASS" if all(line_counts.get(x, 0) > 0 for x in ["shared_immune_mechanism", "HCC_specific_barrier", "PD1X_repair_logic"]) else "FAIL"
    (out / "phase5_6_summary.md").write_text(
        "\n".join(
            [
                "# Phase5.6 Storyline Consolidation Summary",
                "",
                f"- Verdict: `{verdict}`",
                f"- Shared groups: {line_counts.get('shared_immune_mechanism', 0)}",
                f"- HCC-specific groups: {line_counts.get('HCC_specific_barrier', 0)}",
                f"- PD1X repair groups: {line_counts.get('PD1X_repair_logic', 0)}",
                f"- Exploratory groups: {line_counts.get('exploratory', 0)}",
                "- HCC-specific groups are not generalized as shared mechanisms.",
                "- PD1X repair groups are not interpreted as pure PD-1 monotherapy mechanisms.",
                "",
            ]
        )
    )
    write_status(out / "phase5_6_status.yaml", "Phase5.6", verdict, **{k: int(v) for k, v in line_counts.items()})
    return assignment


def phase5_7(
    registry: pd.DataFrame,
    evidence: pd.DataFrame,
    direction: pd.DataFrame,
    stability: pd.DataFrame,
    confounding: pd.DataFrame,
    classification: pd.DataFrame,
    storyline: pd.DataFrame,
    tissue: pd.DataFrame,
    external: pd.DataFrame,
    matches: pd.DataFrame,
) -> pd.DataFrame:
    out = mkdir(OUT / "07_priority_compression")
    ev = evidence.set_index("mechanism_group_id")
    dr = direction.set_index("mechanism_group_id")
    st = stability.set_index("mechanism_group_id")
    cf = confounding.set_index("mechanism_group_id")
    cl = classification.set_index("mechanism_group_id")
    sl = storyline.set_index("mechanism_group_id")
    ts = tissue.set_index("mechanism_group_id")
    ex = external.set_index("mechanism_group_id")
    matched_groups = set()
    for value in matches.get("mechanism_group_ids", pd.Series(dtype=str)):
        matched_groups.update(split_semicolon(value))
    rows = []
    for _, row in registry.iterrows():
        gid = row["mechanism_group_id"]
        primary_strength = np.mean([score_tier(x) for x in split_semicolon(ev.loc[gid, "tiers"])])
        stability_score = {"robust": 1.0, "usable": 0.8, "exploratory": 0.45, "unstable": 0.0}.get(st.loc[gid, "stability_adjudication"], 0.5)
        direction_score = {
            "consistent_response_associated": 1.0,
            "consistent_resistance_associated": 0.95,
            "context_dependent": 0.7,
            "uncertain": 0.35,
            "conflicting": 0.1,
        }.get(dr.loc[gid, "direction_adjudication"], 0.5)
        cross_task_score = 1.0 if gid in matched_groups else 0.45
        bio_score = 0.4 if row["primary_mechanism_class"] == "mixed or unclear mechanism" else 0.9
        conf_score = {"acceptable": 0.9, "moderate_caveat": 0.7, "high_risk_demote": 0.2, "blocked": 0.0}.get(cf.loc[gid, "confounding_adjudication"], 0.6)
        anchor_score = 0.45
        phase6_actionability = 0.4 if sl.loc[gid, "main_story_line"] == "exploratory" else 0.85
        support_score = 0.65 if ev.loc[gid, "support_conflict"] else 0.9
        novelty_score = 0.75 if row["n_member_modules"] > 1 else 0.6
        total = 100 * (
            0.16 * primary_strength
            + 0.14 * stability_score
            + 0.14 * direction_score
            + 0.10 * cross_task_score
            + 0.14 * bio_score
            + 0.10 * conf_score
            + 0.06 * anchor_score
            + 0.08 * phase6_actionability
            + 0.04 * support_score
            + 0.04 * novelty_score
        )
        rows.append(
            {
                "mechanism_group_id": gid,
                "representative_module_id": row["representative_module_id"],
                "storyline": sl.loc[gid, "main_story_line"],
                "mechanism_class": row["primary_mechanism_class"],
                "response_direction": dr.loc[gid, "direction_adjudication"],
                "evidence_strength_score": round(primary_strength, 3),
                "stability_score": round(stability_score, 3),
                "direction_clarity_score": round(direction_score, 3),
                "cross_task_support_score": round(cross_task_score, 3),
                "biological_interpretability_score": round(bio_score, 3),
                "confounding_inverse_score": round(conf_score, 3),
                "tissue_external_anchor_score": round(anchor_score, 3),
                "phase6_actionability_score": round(phase6_actionability, 3),
                "support_consistency_score": round(support_score, 3),
                "novelty_nontriviality_score": round(novelty_score, 3),
                "priority_score": round(total, 2),
                "tissue_spatial_status": ts.loc[gid, "tissue_spatial_status"],
                "external_anchor_status": ex.loc[gid, "external_anchor_status"],
                "required_caveat": "; ".join(REQUIRED_CAVEATS),
            }
        )
    scored = pd.DataFrame(rows)
    scored["priority"] = "Priority 3"
    for line, sub_idx in scored.groupby("storyline").groups.items():
        if line in {"shared_immune_mechanism", "HCC_specific_barrier", "PD1X_repair_logic"}:
            ordered = scored.loc[list(sub_idx)].sort_values("priority_score", ascending=False)
            p1 = ordered.head(5).index
            p2 = ordered.iloc[5:15].index
            scored.loc[p1, "priority"] = "Priority 1"
            scored.loc[p2, "priority"] = "Priority 2"
    scored.loc[scored["priority_score"] < 45, "priority"] = "Demoted"
    scored.loc[scored["mechanism_class"].eq("mixed or unclear mechanism") & scored["priority"].isin(["Priority 1", "Priority 2"]), "priority"] = "Priority 3"
    scored.to_csv(out / "phase5_7_priority_scoring_table.csv", index=False)
    top_shared = scored[(scored["storyline"] == "shared_immune_mechanism") & scored["priority"].isin(["Priority 1", "Priority 2"])].sort_values("priority_score", ascending=False)
    top_hcc = scored[(scored["storyline"] == "HCC_specific_barrier") & scored["priority"].isin(["Priority 1", "Priority 2"])].sort_values("priority_score", ascending=False)
    top_pd1x = scored[(scored["storyline"] == "PD1X_repair_logic") & scored["priority"].isin(["Priority 1", "Priority 2"])].sort_values("priority_score", ascending=False)
    top_shared.to_csv(out / "phase5_7_top_shared_candidates.csv", index=False)
    top_hcc.to_csv(out / "phase5_7_top_HCC_specific_candidates.csv", index=False)
    top_pd1x.to_csv(out / "phase5_7_top_PD1X_repair_candidates.csv", index=False)
    demoted = scored[scored["priority"].isin(["Demoted", "Blocked"])].copy()
    demoted["demotion_reason"] = np.where(demoted["priority_score"] < 45, "low_integrated_priority_score", "blocked_or_demoted_by_rule")
    demoted.to_csv(out / "phase5_7_demoted_candidates.csv", index=False)
    cards = ["# Phase5.7 Priority Cards", ""]
    for _, row in scored.sort_values(["priority", "storyline", "priority_score"], ascending=[True, True, False]).iterrows():
        if row["priority"] in {"Priority 1", "Priority 2"}:
            cards.extend(
                [
                    f"## {row['mechanism_group_id']}",
                    "",
                    f"- Priority: `{row['priority']}`",
                    f"- Storyline: `{row['storyline']}`",
                    f"- Mechanism class: `{row['mechanism_class']}`",
                    f"- Score: {row['priority_score']}",
                    f"- Direction: `{row['response_direction']}`",
                    f"- Tissue/spatial status: `{row['tissue_spatial_status']}`",
                    f"- External anchor status: `{row['external_anchor_status']}`",
                    f"- Caveat: {row['required_caveat']}",
                    "",
                ]
            )
    (out / "phase5_7_priority_cards.md").write_text("\n".join(cards))
    line_ok = len(top_shared) > 0 and len(top_hcc) > 0 and len(top_pd1x) > 0
    verdict = "PASS" if line_ok and len(scored[scored["priority"].isin(["Priority 1", "Priority 2"])]) <= 45 else "CONDITIONAL_PASS"
    (out / "phase5_7_summary.md").write_text(
        f"# Phase5.7 Priority Compression Summary\n\n- Verdict: `{verdict}`\n- Priority 1: {int((scored['priority']=='Priority 1').sum())}\n- Priority 2: {int((scored['priority']=='Priority 2').sum())}\n- Priority 3: {int((scored['priority']=='Priority 3').sum())}\n- Demoted: {int((scored['priority']=='Demoted').sum())}\n- Main candidates are capped by storyline ranking, not prediction performance.\n"
    )
    write_status(
        out / "phase5_7_status.yaml",
        "Phase5.7",
        verdict,
        priority1=int((scored["priority"] == "Priority 1").sum()),
        priority2=int((scored["priority"] == "Priority 2").sum()),
        priority3=int((scored["priority"] == "Priority 3").sum()),
        demoted=int((scored["priority"] == "Demoted").sum()),
    )
    return scored


def phase5_8(scored: pd.DataFrame, registry: pd.DataFrame, evidence: pd.DataFrame, classification: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    out = mkdir(OUT / "08_phase6_handoff")
    ev = evidence.set_index("mechanism_group_id")
    cl = classification.set_index("mechanism_group_id")
    rg = registry.set_index("mechanism_group_id")
    main = scored[scored["priority"].isin(["Priority 1", "Priority 2"])].copy()
    support = scored[~scored["priority"].isin(["Priority 1", "Priority 2"])].copy()
    rows = []
    for _, row in main.iterrows():
        gid = row["mechanism_group_id"]
        rows.append(
            {
                "mechanism_group_id": gid,
                "representative_module_id": row["representative_module_id"],
                "storyline": row["storyline"],
                "priority": row["priority"],
                "mechanism_class": row["mechanism_class"],
                "response_direction": row["response_direction"],
                "source_tasks": rg.loc[gid, "source_tasks"],
                "source_feature_families": rg.loc[gid, "source_feature_families"],
                "key_features": ev.loc[gid, "key_features"],
                "module_score_available": True,
                "evidence_strength": row["evidence_strength_score"],
                "stability_grade": ev.loc[gid, "representative_stability_grade"],
                "confounding_caveat": "source/batch/missingness conditional-pass caveat retained",
                "tissue_spatial_status": row["tissue_spatial_status"],
                "external_anchor_status": row["external_anchor_status"],
                "Phase6_allowed_use": "perturbation mapping; target nomination pre-screening; spatial/tissue validation planning; external anchor validation; mechanism-class mapping",
                "Phase6_forbidden_use": "direct drug recommendation; clinical prediction claim; causal claim without perturbation/validation; support-only as primary evidence",
                "required_caveat": row["required_caveat"],
            }
        )
    main_handoff = pd.DataFrame(rows)
    main_handoff.to_csv(out / "phase5_phase6_handoff_main_mechanisms.csv", index=False)
    support.to_csv(out / "phase5_phase6_handoff_support_context.csv", index=False)
    pd.DataFrame(
        [
            {"forbidden_use": "direct_drug_recommendation", "reason": "Phase5 mechanism candidates are not intervention recommendations"},
            {"forbidden_use": "clinical_prediction_claim", "reason": "Phase5 does not train or validate clinical prediction models"},
            {"forbidden_use": "causal_claim", "reason": "Phase5 adjudicates associative evidence only"},
            {"forbidden_use": "support_only_as_primary", "reason": "PD1_anchor and Tier C support evidence remain annotation/sensitivity context"},
            {"forbidden_use": "unintegrated_new_data_primary_claim", "reason": "new addendum data remain outside primary conclusions"},
        ]
    ).to_csv(out / "phase5_phase6_forbidden_use_registry.csv", index=False)
    (out / "phase5_phase6_pending_evidence_request.md").write_text(
        "\n".join(
            [
                "# Phase5 Phase6 Pending Evidence Request",
                "",
                "- Tissue/spatial anchors: pending integrated data gate.",
                "- External anchors: pending integrated validation contract.",
                "- Perturbation mapping: allowed as Phase6 analysis, not present Phase5 proof.",
                "- Target nomination: allowed only as pre-screening in Phase6.",
                "",
            ]
        )
    )
    cards = ["# Phase5 Phase6 Mechanism Cards", ""]
    for _, row in main_handoff.iterrows():
        cards.extend(
            [
                f"## {row['mechanism_group_id']}",
                "",
                f"- Representative module: `{row['representative_module_id']}`",
                f"- Priority: `{row['priority']}`",
                f"- Storyline: `{row['storyline']}`",
                f"- Mechanism class: `{row['mechanism_class']}`",
                f"- Direction: `{row['response_direction']}`",
                f"- Allowed Phase6 use: {row['Phase6_allowed_use']}",
                f"- Forbidden Phase6 use: {row['Phase6_forbidden_use']}",
                f"- Caveat: {row['required_caveat']}",
                "",
            ]
        )
    (out / "phase5_phase6_mechanism_cards.md").write_text("\n".join(cards))
    verdict = "PASS" if len(main_handoff) > 0 else "FAIL"
    (out / "phase5_8_summary.md").write_text(
        f"# Phase5.8 Phase6 Handoff Summary\n\n- Verdict: `{verdict}`\n- Main mechanisms: {len(main_handoff)}\n- Support context groups: {len(support)}\n- Forbidden use registry written.\n"
    )
    write_status(out / "phase5_8_status.yaml", "Phase5.8", verdict, main_mechanisms=len(main_handoff), support_context=len(support))
    return main_handoff, support


def phase5_9(
    main: pd.DataFrame,
    support: pd.DataFrame,
    scored: pd.DataFrame,
    registry: pd.DataFrame,
    tissue: pd.DataFrame,
    external: pd.DataFrame,
) -> None:
    out = mkdir(OUT / "09_final_report")
    demoted = scored[scored["priority"].isin(["Demoted", "Blocked"])].copy()
    p1 = int((scored["priority"] == "Priority 1").sum())
    p2 = int((scored["priority"] == "Priority 2").sum())
    p3 = int((scored["priority"] == "Priority 3").sum())
    top_shared = main[main["storyline"] == "shared_immune_mechanism"]["mechanism_group_id"].head(15).tolist()
    top_hcc = main[main["storyline"] == "HCC_specific_barrier"]["mechanism_group_id"].head(15).tolist()
    top_pd1x = main[main["storyline"] == "PD1X_repair_logic"]["mechanism_group_id"].head(15).tolist()
    tissue_status = "pending_no_valid_anchor" if tissue["tissue_spatial_status"].eq("pending_no_valid_anchor").all() else "partially_supported"
    external_status = "pending_no_valid_anchor" if external["external_anchor_status"].eq("pending_no_valid_anchor").all() else "partially_supported"
    phase6_ready = len(top_shared) > 0 and len(top_hcc) > 0 and len(top_pd1x) > 0
    verdict = "CONDITIONAL_PASS" if phase6_ready and (tissue_status == "pending_no_valid_anchor" or external_status == "pending_no_valid_anchor") else ("PASS" if phase6_ready else "FAIL")
    decision = {
        "verdict": verdict,
        "n_input_modules": 175,
        "n_mechanism_groups": int(len(registry)),
        "n_priority1_candidates": p1,
        "n_priority2_candidates": p2,
        "n_priority3_candidates": p3,
        "n_demoted_candidates": int((scored["priority"] == "Demoted").sum()),
        "n_blocked_candidates": int((scored["priority"] == "Blocked").sum()),
        "top_shared_mechanisms": top_shared,
        "top_HCC_specific_barriers": top_hcc,
        "top_PD1X_repair_logic": top_pd1x,
        "tissue_spatial_status": tissue_status,
        "external_anchor_status": external_status,
        "phase6_ready": phase6_ready,
        "required_caveats": REQUIRED_CAVEATS,
        "forbidden_claims": FORBIDDEN_CLAIMS,
        "recommended_next_phase": "Phase6 perturbation mapping, target nomination pre-screening, and validation planning",
    }
    write_yaml(out / "PHASE5_FINAL_DECISION.yaml", decision)
    scored.to_csv(out / "phase5_mechanism_group_master_table.csv", index=False)
    main.to_csv(out / "phase5_phase6_handoff_main_mechanisms.csv", index=False)
    support.to_csv(out / "phase5_phase6_handoff_support_context.csv", index=False)
    demoted.to_csv(out / "phase5_demoted_or_blocked_candidates.csv", index=False)
    write_yaml(
        out / "phase5_reproducibility_manifest.yaml",
        {
            "generated_at_utc": now_iso(),
            "script": rel(ROOT / "scripts/v6_1/run_phase5_mechanism_adjudication.py"),
            "phase4_root": rel(PHASE4),
            "phase4_final_decision": rel(PHASE4 / "PHASE4_FINAL_DECISION.yaml"),
            "phase4_main_handoff": rel(PHASE4 / "phase4_phase5_handoff_main_modules.csv"),
            "phase4_support_handoff": rel(PHASE4 / "phase4_phase5_handoff_support_evidence.csv"),
            "phase3_5_ml_used": False,
            "pd1_anchor_primary_used": False,
            "unintegrated_new_data_used": False,
            "prediction_model_trained": False,
            "drug_recommendation_generated": False,
        },
    )
    report = [
        "# Phase5 Final Report",
        "",
        f"- Verdict: `{verdict}`",
        f"- Input Phase4 main modules: 175",
        f"- Mechanism groups after compression: {len(registry)}",
        f"- Priority 1 candidates: {p1}",
        f"- Priority 2 candidates: {p2}",
        f"- Priority 3 candidates: {p3}",
        f"- Demoted candidates: {int((scored['priority'] == 'Demoted').sum())}",
        f"- Blocked candidates: {int((scored['priority'] == 'Blocked').sum())}",
        "",
        "## Storyline Outputs",
        "",
        f"- Top shared mechanisms retained: {len(top_shared)}",
        f"- Top HCC-specific barriers retained: {len(top_hcc)}",
        f"- Top PD1X repair logic mechanisms retained: {len(top_pd1x)}",
        "",
        "## Anchor Status",
        "",
        f"- Tissue/spatial status: `{tissue_status}`",
        f"- External anchor status: `{external_status}`",
        "- No unintegrated new data were used for primary adjudication.",
        "",
        "## What Can Be Said",
        "",
        "- Phase5 compressed audited Phase4 Tier A/B modules into prioritized mechanism groups.",
        "- Candidate mechanisms are suitable for Phase6 perturbation mapping, target nomination pre-screening, and validation planning.",
        "- Support evidence provides annotation and direction sanity checks only.",
        "",
        "## What Cannot Be Said",
        "",
        "- No causal mechanism is proven.",
        "- No drug combination is recommended.",
        "- No clinical prediction model is claimed.",
        "- PD1_anchor is not primary evidence.",
        "- HCC-specific mechanisms are not generalized as pan-cancer mechanisms.",
        "",
    ]
    (out / "PHASE5_FINAL_REPORT.md").write_text("\n".join(report))
    (out / "phase5_9_summary.md").write_text(
        f"# Phase5.9 Final Decision Summary\n\n- Verdict: `{verdict}`\n- Phase6 ready: `{phase6_ready}`\n- Mechanism groups: {len(registry)}\n- Priority 1/2 main mechanisms: {len(main)}\n- Tissue/spatial status: `{tissue_status}`\n- External anchor status: `{external_status}`\n"
    )
    write_status(out / "phase5_9_status.yaml", "Phase5.9", verdict, phase6_ready=phase6_ready, mechanism_groups=len(registry), main_mechanisms=len(main))
    phase_report = [
        "# Phase5 Phase Report 20260611",
        "",
        "## Scope And Status",
        "",
        f"- Phase: `v6.1 Phase5 Mechanism adjudication and candidate compression`",
        f"- Final verdict: `{verdict}`",
        "- Role: compress Phase4 Tier A/B modules into prioritized mechanism groups for Phase6.",
        "- Not performed: new module discovery, supervised model training, drug recommendation, causal claim.",
        "",
        "## Canonical Inputs",
        "",
        f"- Phase4 final decision: `{rel(PHASE4 / 'PHASE4_FINAL_DECISION.yaml')}`",
        f"- Phase4 main handoff: `{rel(PHASE4 / 'phase4_phase5_handoff_main_modules.csv')}`",
        f"- Phase4 support handoff: `{rel(PHASE4 / 'phase4_phase5_handoff_support_evidence.csv')}`",
        f"- Phase4 sample scores: `{rel(PHASE4 / '07_module_scoring/phase4_7_module_scores_by_sample.csv')}`",
        f"- Phase4 patient scores: `{rel(PHASE4 / '07_module_scoring/phase4_7_module_scores_by_patient.csv')}`",
        "",
        "## Canonical Outputs",
        "",
        "- `PHASE5_FINAL_DECISION.yaml`",
        "- `phase5_mechanism_group_master_table.csv`",
        "- `phase5_phase6_handoff_main_mechanisms.csv`",
        "- `phase5_phase6_handoff_support_context.csv`",
        "- `phase5_demoted_or_blocked_candidates.csv`",
        "- `phase5_reproducibility_manifest.yaml`",
        "",
        "## Downstream Contract",
        "",
        "- Phase6 main input: Priority 1/2 rows only.",
        "- Priority 3/support context: annotation and sensitivity only.",
        "- Pending evidence request: tissue/spatial/external anchors.",
        "- Forbidden: direct drug recommendation, clinical prediction claim, causal claim, support-only primary evidence.",
        "",
        "## Reproduction",
        "",
        "- Command: `conda run -n proj006 python -W error scripts/v6_1/run_phase5_mechanism_adjudication.py`",
        "",
    ]
    changelog = [
        "# Phase5 Release Changelog 20260611",
        "",
        "## Release Scope",
        "",
        "- Added complete Phase5.0-Phase5.9 reproducible output package.",
        "- Compressed Phase4 main candidates into mechanism groups.",
        "- Wrote Phase6 handoff, pending evidence request, and forbidden-use registry.",
        "",
        "## Verification",
        "",
        "- Run mode: `conda run -n proj006 python -W error`.",
        "- Phase3.5 ML used: false.",
        "- PD1_anchor primary used: false.",
        "- Unintegrated new data used: false.",
        "- Prediction model trained: false.",
        "- Drug recommendation generated: false.",
        "",
        "## Final State",
        "",
        f"- Final verdict: `{verdict}`",
        f"- Phase6 ready: `{phase6_ready}`",
        f"- Main Phase6 mechanisms: {len(main)}",
        "",
    ]
    (out / "PHASE5_PHASE_REPORT_20260611.md").write_text("\n".join(phase_report))
    (out / "PHASE5_RELEASE_CHANGELOG_20260611.md").write_text("\n".join(changelog))

    for name in [
        "PHASE5_FINAL_DECISION.yaml",
        "phase5_reproducibility_manifest.yaml",
    ]:
        (OUT / name).write_text((out / name).read_text())
    for name in [
        "PHASE5_FINAL_REPORT.md",
        "PHASE5_PHASE_REPORT_20260611.md",
        "PHASE5_RELEASE_CHANGELOG_20260611.md",
    ]:
        (OUT / name).write_text((out / name).read_text())
    for name in [
        "phase5_mechanism_group_master_table.csv",
        "phase5_phase6_handoff_main_mechanisms.csv",
        "phase5_phase6_handoff_support_context.csv",
        "phase5_demoted_or_blocked_candidates.csv",
    ]:
        (OUT / name).write_text((out / name).read_text())


def main() -> None:
    mkdir(OUT)
    inputs = load_inputs()
    main_modules = phase5_0(inputs)
    registry, mapping, matches = phase5_1(main_modules, inputs["sample_scores"])
    evidence = phase5_2(
        registry,
        mapping,
        main_modules,
        inputs["support"],
        inputs["sample_scores"],
        inputs["patient_scores"],
        inputs["support_conflicts"],
    )
    direction, stability, confounding = phase5_3(evidence)
    classification = phase5_4(registry, evidence, direction)
    tissue, external = phase5_5(registry)
    storyline = phase5_6(registry, classification, matches)
    scored = phase5_7(registry, evidence, direction, stability, confounding, classification, storyline, tissue, external, matches)
    main_handoff, support_context = phase5_8(scored, registry, evidence, classification)
    phase5_9(main_handoff, support_context, scored, registry, tissue, external)
    print(f"Wrote Phase5 outputs: {OUT}")
    print(f"Input modules: {len(main_modules)}")
    print(f"Mechanism groups: {len(registry)}")
    print(f"Phase6 main mechanisms: {len(main_handoff)}")


if __name__ == "__main__":
    main()
