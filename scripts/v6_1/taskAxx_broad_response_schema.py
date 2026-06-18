from __future__ import annotations

from pathlib import Path
from typing import Dict, List, Tuple

import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
STEP1_DIR = ROOT / "results" / "v6_1" / "step1"

SAMPLE_PATH = STEP1_DIR / "sample_metadata_master_v6_1.csv"
PATIENT_PATH = STEP1_DIR / "patient_metadata_master_v6_1.csv"
MAPPING_PATH = STEP1_DIR / "response_mapping_dictionary_v6_1.csv"
COHORT_FREEZE_PATH = STEP1_DIR / "frozen_cohort_inclusion_v6_1.csv"
PATIENT_SPLIT_PATH = STEP1_DIR / "frozen_patient_split_v6_1.csv"

DICT_OUT = STEP1_DIR / "broad_response_mapping_dictionary_v6_1.csv"
SAMPLE_OUT = STEP1_DIR / "sample_metadata_master_v6_1.broad_response.csv"
PATIENT_OUT = STEP1_DIR / "patient_metadata_master_v6_1.broad_response.csv"
REPORT_OUT = STEP1_DIR / "response_schema_report.md"
POLICY_OUT = STEP1_DIR / "response_usage_policy.md"
SOURCE_OUT = STEP1_DIR / "source_records_broad_response.tsv"


EVIDENCE_SCORE = {
    "A_clinical_radiologic": 4,
    "B_pathologic_or_study_defined": 3,
    "C_biological_proxy": 2,
    "D_unknown_or_not_evaluable": 1,
}

TIMEPOINT_PRIORITY = {
    "pre": 1,
    "on_treatment": 2,
    "post": 3,
    "long_post": 4,
    "unknown": 5,
}


def clean_str(value: object) -> str:
    if pd.isna(value):
        return ""
    text = str(value).strip()
    if text.lower() in {"", "nan", "na", "n/a", "none", "null"}:
        return ""
    return text


def norm(value: object) -> str:
    return clean_str(value).upper()


def truthy(value: object) -> bool:
    return str(value).strip().lower() in {"true", "1", "yes", "y", "t"}


def make_record(
    response_label_type_original: str,
    response_strict_binary: str,
    response_broad_binary: str,
    response_ordered: str,
    response_response_type: str,
    response_evidence_level: str,
    mapping_confidence: str,
    mapping_rule: str,
    should_use_for_strict_baseline: bool,
    should_use_for_broad_baseline: bool,
    should_use_for_proxy_baseline: bool,
    notes: str,
) -> Dict[str, object]:
    return {
        "response_label_type_original": response_label_type_original,
        "response_strict_binary": response_strict_binary,
        "response_broad_binary": response_broad_binary,
        "response_ordered": response_ordered,
        "response_response_type": response_response_type,
        "response_evidence_level": response_evidence_level,
        "mapping_confidence": mapping_confidence,
        "mapping_rule": mapping_rule,
        "should_use_for_strict_baseline": should_use_for_strict_baseline,
        "should_use_for_broad_baseline": should_use_for_broad_baseline,
        "should_use_for_proxy_baseline": should_use_for_proxy_baseline,
        "notes": notes,
    }


def map_broad_response(cohort_id: str, response_standard: str, response_raw: str) -> Dict[str, object]:
    raw = norm(response_raw)
    standard = norm(response_standard) or "unknown"

    unknown = make_record(
        response_label_type_original=standard,
        response_strict_binary="unknown",
        response_broad_binary="unknown",
        response_ordered="not_evaluable",
        response_response_type="unknown",
        response_evidence_level="D_unknown_or_not_evaluable",
        mapping_confidence="high",
        mapping_rule="unknown_or_not_evaluable",
        should_use_for_strict_baseline=False,
        should_use_for_broad_baseline=False,
        should_use_for_proxy_baseline=False,
        notes="Unknown / not evaluable label kept outside all supervised baselines.",
    )

    if raw in {"", "UNKNOWN", "NE", "UN", "-", "NA", "N/A", "NOT EVALUABLE", "NOT AVAILABLE"}:
        return unknown

    if standard in {"RECIST", "MRECIST"}:
        if raw in {"CR"}:
            return make_record(
                standard, "responder", "responder", "complete_response", "RECIST_or_mRECIST",
                "A_clinical_radiologic", "high", "direct_RECIST_complete_response",
                True, True, False,
                "Canonical RECIST/mRECIST complete response.",
            )
        if raw in {"PR", "OR"}:
            return make_record(
                standard, "responder", "responder", "partial_or_objective_response", "RECIST_or_mRECIST",
                "A_clinical_radiologic", "high", "direct_RECIST_partial_or_objective_response",
                True, True, False,
                "Canonical RECIST/mRECIST partial or objective response.",
            )
        if raw in {"PD"}:
            return make_record(
                standard, "non_responder", "non_responder", "progressive_disease", "RECIST_or_mRECIST",
                "A_clinical_radiologic", "high", "direct_RECIST_progressive_disease",
                True, True, False,
                "Canonical RECIST/mRECIST progressive disease.",
            )
        if raw in {"SD"}:
            return make_record(
                standard, "stable_disease", "stable_disease", "stable_disease", "disease_control",
                "A_clinical_radiologic", "high", "RECIST_stable_disease_preserved",
                False, False, False,
                "Stable disease preserved; may enter clinical-benefit sensitivity only.",
            )
        if raw in {"NON CR/PD"}:
            return make_record(
                standard, "unknown", "disease_control_uncertain", "non_CR_PD", "disease_control",
                "A_clinical_radiologic", "medium", "RECIST_non_CR_PD_preserved",
                False, False, False,
                "NON CR/PD preserved as disease-control-uncertain; not collapsed into responder.",
            )
        if raw in {"RESPONDER"}:
            return make_record(
                standard, "responder", "responder", "study_defined_responder", "RECIST_or_mRECIST",
                "A_clinical_radiologic", "medium", "RECIST_literal_responder",
                True, True, False,
                "Literal responder label under RECIST context treated as clinical radiologic response.",
            )
        if raw in {"NON-RESPONDER", "NON_RESPONDER", "NONRESPONDER"}:
            return make_record(
                standard, "non_responder", "non_responder", "study_defined_non_responder", "RECIST_or_mRECIST",
                "A_clinical_radiologic", "medium", "RECIST_literal_non_responder",
                True, True, False,
                "Literal non-responder label under RECIST context treated as clinical radiologic non-response.",
            )
        if raw in {"R"}:
            return make_record(
                standard, "responder", "responder", "study_defined_responder", "RECIST_or_mRECIST",
                "A_clinical_radiologic", "medium", "RECIST_R_shorthand_responder",
                True, True, False,
                "Responder shorthand under explicit RECIST context treated as clinical response.",
            )
        if raw in {"NR"}:
            return make_record(
                standard, "non_responder", "non_responder", "study_defined_non_responder", "RECIST_or_mRECIST",
                "A_clinical_radiologic", "medium", "RECIST_NR_shorthand_non_responder",
                True, True, False,
                "Non-responder shorthand under explicit RECIST context treated as clinical non-response.",
            )

    if standard == "RECIST_LIKE":
        if raw in {"PR"}:
            return make_record(
                standard, "unknown", "responder", "partial_response", "study_defined_R_NR",
                "B_pathologic_or_study_defined", "medium", "RECIST_like_partial_response",
                False, True, False,
                "RECIST-like partial response used only in broad baseline.",
            )
        if raw in {"PD"}:
            return make_record(
                standard, "unknown", "non_responder", "progressive_disease", "study_defined_R_NR",
                "B_pathologic_or_study_defined", "medium", "RECIST_like_progressive_disease",
                False, True, False,
                "RECIST-like progression used only in broad baseline.",
            )
        if raw in {"SD"}:
            return make_record(
                standard, "unknown", "stable_disease", "stable_disease", "disease_control",
                "B_pathologic_or_study_defined", "medium", "RECIST_like_stable_disease",
                False, False, False,
                "RECIST-like stable disease preserved for sensitivity only.",
            )
        return unknown | {"response_label_type_original": standard}

    if standard.startswith("PATHOLOGICAL_RESPONSE:PCR") or standard == "PCR":
        if raw in {"PR", "R", "R1", "R2", "R;R1", "R;R2"}:
            ordered = "pathologic_complete_response" if raw == "PR" else "pathologic_response_positive"
            rtype = "pathologic_complete_response"
            return make_record(
                standard, "unknown", "responder", ordered, rtype,
                "B_pathologic_or_study_defined", "medium", "pathologic_complete_or_positive_response",
                False, True, False,
                "Pathologic complete/positive response enters broad baseline only.",
            )
        if raw in {"NR"}:
            return make_record(
                standard, "unknown", "non_responder", "pathologic_non_response", "pathologic_complete_response",
                "B_pathologic_or_study_defined", "medium", "pathologic_non_response",
                False, True, False,
                "Pathologic non-response enters broad baseline only.",
            )

    if standard.startswith("PATHOLOGICAL_RESPONSE:NON-PCR"):
        if raw in {"NE"}:
            return unknown | {"response_label_type_original": standard}

    if standard == "PATHOLOGICAL_RESPONSE":
        if raw in {"R", "NR"}:
            return make_record(
                standard,
                "unknown",
                "responder" if raw == "R" else "non_responder",
                "study_defined_responder" if raw == "R" else "study_defined_non_responder",
                "study_defined_R_NR",
                "B_pathologic_or_study_defined",
                "medium",
                "pathologic_study_defined_R_NR",
                False,
                True,
                False,
                "Study-defined R/NR under pathological response context used only in broad baseline.",
            )
        if raw in {"HIGH", "LOW", "MEDIUM"}:
            broad = "responder" if raw == "HIGH" else "non_responder" if raw == "LOW" else "intermediate"
            return make_record(
                standard, "unknown", broad, f"T_cell_expansion_{raw.lower()}", "T_cell_expansion",
                "C_biological_proxy", "high", "T_cell_expansion_proxy",
                False, False, raw in {"HIGH", "LOW"},
                "T-cell expansion proxy; High/Low may enter proxy baseline, Medium remains intermediate only.",
            )

    if standard == "TUMOR_SIZE_CHANGE":
        if raw == "PR":
            return make_record(
                standard, "unknown", "responder", "tumor_size_reduction", "biological_proxy",
                "C_biological_proxy", "medium", "study_defined_tumor_size_partial_response",
                False, False, True,
                "Study-specific tumor-size reduction retained as proxy response, not strict clinical response.",
            )
        if raw == "NE":
            return unknown | {"response_label_type_original": standard}

    if standard == "TUMOR_SIZE":
        return unknown | {"response_label_type_original": standard}

    if standard == "UNKNOWN":
        return unknown | {"response_label_type_original": standard}

    if cohort_id == "lambrecht_hcc" and raw == "DBI":
        return make_record(
            standard, "unknown", "responder", "durable_benefit", "clinical_benefit",
            "B_pathologic_or_study_defined", "medium", "DBI_durable_benefit_inferred_from_study_context",
            False, False, False,
            "DBI interpreted as durable clinical benefit in atezo/bev HCC study context; kept outside default broad baseline pending explicit metadata confirmation. Source: Nature Communications 2023 article defines response as objective response at 3 months or disease control for >=6 months.",
        )

    return unknown | {"response_label_type_original": standard}


def build_dictionary(sample: pd.DataFrame) -> pd.DataFrame:
    rows: List[Dict[str, object]] = []
    grouped = sample.groupby(["cohort_id", "response_standard", "response_raw"], dropna=False)
    for (cohort_id, response_standard, response_raw), sub in grouped:
        rec = map_broad_response(str(cohort_id), clean_str(response_standard) or "unknown", clean_str(response_raw) or "unknown")
        row = {
            "cohort_id": cohort_id,
            "response_raw": clean_str(response_raw) or "unknown",
            **rec,
        }
        rows.append(row)
    return pd.DataFrame(rows).sort_values(["cohort_id", "response_label_type_original", "response_raw"]).reset_index(drop=True)


def attach_sample_schema(sample: pd.DataFrame, dictionary: pd.DataFrame) -> pd.DataFrame:
    lookup = dictionary.set_index(["cohort_id", "response_raw", "response_label_type_original"]).to_dict("index")

    out_rows: List[Dict[str, object]] = []
    for _, row in sample.iterrows():
        cohort_id = row["cohort_id"]
        raw = clean_str(row["response_raw"]) or "unknown"
        standard = clean_str(row["response_standard"]).upper() or "unknown"
        rec = lookup[(cohort_id, raw, standard)]
        new_row = row.to_dict()
        new_row["response_strict_binary"] = rec["response_strict_binary"]
        new_row["response_broad_binary"] = rec["response_broad_binary"]
        new_row["response_ordered_broad"] = rec["response_ordered"]
        new_row["response_response_type"] = rec["response_response_type"]
        new_row["response_evidence_level"] = rec["response_evidence_level"]
        new_row["response_use_strict"] = rec["should_use_for_strict_baseline"]
        new_row["response_use_broad"] = rec["should_use_for_broad_baseline"]
        new_row["response_use_proxy"] = rec["should_use_for_proxy_baseline"]
        new_row["response_mapping_confidence_broad"] = rec["mapping_confidence"]
        new_row["response_mapping_note_broad"] = rec["notes"]
        out_rows.append(new_row)
    return pd.DataFrame(out_rows)


def ordered_unique_join(values: pd.Series, priority_map: Dict[str, int] | None = None) -> str:
    uniq = {clean_str(v) for v in values.tolist() if clean_str(v)}
    if not uniq:
        return ""
    if priority_map:
        ordered = sorted(uniq, key=lambda x: (priority_map.get(x, 999), x))
    else:
        ordered = sorted(uniq)
    return ";".join(ordered)


def aggregate_patient(sample_broad: pd.DataFrame, patient: pd.DataFrame, patient_split: pd.DataFrame) -> pd.DataFrame:
    split_cols = ["cohort_id", "patient_id", "patient_key", "split_label", "fold_id", "supervised_anchor_available"]
    split = patient_split[split_cols].drop_duplicates()
    rows: List[Dict[str, object]] = []

    for _, prow in patient.iterrows():
        cohort_id = prow["cohort_id"]
        patient_id = prow["patient_id"]
        sub = sample_broad[(sample_broad["cohort_id"] == cohort_id) & (sample_broad["patient_id"] == patient_id)].copy()

        usable = sub[
            sub["response_evidence_level"].isin(EVIDENCE_SCORE.keys())
            & (
                sub["response_use_strict"].map(truthy)
                | sub["response_use_broad"].map(truthy)
                | sub["response_use_proxy"].map(truthy)
                | sub["response_strict_binary"].isin(["stable_disease"])
                | sub["response_broad_binary"].isin(["stable_disease", "disease_control_uncertain", "intermediate"])
            )
        ].copy()

        if usable.empty:
            selected = None
            conflict_note = "no_usable_response_label"
            patient_response_conflict = False
        else:
            usable["evidence_score"] = usable["response_evidence_level"].map(EVIDENCE_SCORE)
            usable["timepoint_score"] = usable["timepoint_normalized"].map(TIMEPOINT_PRIORITY).fillna(9)
            usable["use_priority"] = usable.apply(
                lambda r: 3 if truthy(r["response_use_strict"]) else 2 if truthy(r["response_use_broad"]) else 1 if truthy(r["response_use_proxy"]) else 0,
                axis=1,
            )
            usable = usable.sort_values(
                ["evidence_score", "use_priority", "timepoint_score", "sample_id"],
                ascending=[False, False, True, True],
            ).reset_index(drop=True)
            selected = usable.iloc[0]
            unique_labels = usable[["response_strict_binary", "response_broad_binary", "response_evidence_level"]].astype(str).drop_duplicates()
            patient_response_conflict = len(unique_labels) > 1
            conflict_note = "multiple_response_labels_across_samples" if patient_response_conflict else "single_consistent_response_label"

        out = prow.to_dict()
        if selected is None:
            out["response_source_sample_id"] = "none"
            out["response_source_timepoint"] = "unknown"
            out["response_strict_binary"] = "unknown"
            out["response_broad_binary"] = "unknown"
            out["response_ordered_broad"] = "not_evaluable"
            out["response_response_type"] = "unknown"
            out["response_evidence_level"] = "D_unknown_or_not_evaluable"
            out["response_use_strict"] = False
            out["response_use_broad"] = False
            out["response_use_proxy"] = False
            out["response_mapping_confidence_broad"] = "high"
            out["response_mapping_note_broad"] = "No usable broad-response label at patient level."
        else:
            out["response_source_sample_id"] = selected["sample_id"]
            out["response_source_timepoint"] = selected["timepoint_normalized"]
            out["response_strict_binary"] = selected["response_strict_binary"]
            out["response_broad_binary"] = selected["response_broad_binary"]
            out["response_ordered_broad"] = selected["response_ordered_broad"]
            out["response_response_type"] = selected["response_response_type"]
            out["response_evidence_level"] = selected["response_evidence_level"]
            out["response_use_strict"] = bool(selected["response_use_strict"])
            out["response_use_broad"] = bool(selected["response_use_broad"])
            out["response_use_proxy"] = bool(selected["response_use_proxy"])
            out["response_mapping_confidence_broad"] = selected["response_mapping_confidence_broad"]
            out["response_mapping_note_broad"] = selected["response_mapping_note_broad"]

        pre_strict = sub[
            (sub["timepoint_normalized"] == "pre")
            & (sub["response_use_strict"].map(truthy))
            & (sub["response_strict_binary"].isin(["responder", "non_responder"]))
        ]
        out["patient_response_conflict"] = patient_response_conflict
        out["patient_response_conflict_note"] = conflict_note
        out["clean_pre_treatment_strict_anchor"] = not pre_strict.empty
        rows.append(out)

    out_df = pd.DataFrame(rows).merge(split, on=["cohort_id", "patient_id"], how="left")
    return out_df


def add_step2_compatibility_aliases(
    sample_broad: pd.DataFrame,
    patient_broad: pd.DataFrame,
    patient_split: pd.DataFrame,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    sample_out = sample_broad.copy()
    split_cols = ["cohort_id", "patient_id", "split_label", "fold_id"]
    split = patient_split[[c for c in split_cols if c in patient_split.columns]].drop_duplicates()
    if not split.empty:
        sample_out = sample_out.merge(split, on=["cohort_id", "patient_id"], how="left")
    else:
        sample_out["split_label"] = ""
        sample_out["fold_id"] = ""

    sample_out["timepoint"] = sample_out.get("timepoint_normalized", "").fillna("")
    sample_out["split"] = sample_out.get("split_label", "").fillna("")

    patient_out = patient_broad.copy()
    agg_spec = {
        "sample_id": ("sample_id", "first"),
        "timepoint": ("timepoint_normalized", lambda s: ordered_unique_join(s, TIMEPOINT_PRIORITY)),
        "tissue_source": ("tissue_source", ordered_unique_join),
    }
    if "treatment_context" in sample_broad.columns:
        agg_spec["treatment_context"] = ("treatment_context", "first")
    grouped = sample_broad.groupby(["cohort_id", "patient_id"], dropna=False).agg(**agg_spec).reset_index()
    patient_out = patient_out.merge(grouped, on=["cohort_id", "patient_id"], how="left")
    response_source_sample = patient_out.get("response_source_sample_id", "").fillna("")
    patient_out["sample_id"] = response_source_sample.where(
        ~response_source_sample.astype(str).str.strip().str.lower().isin(["", "none", "unknown"]),
        patient_out.get("sample_id", "").fillna(""),
    )
    patient_out["timepoint"] = patient_out.get("timepoint", "").fillna("")
    patient_out["tissue_source"] = patient_out.get("tissue_source", "").fillna("")
    patient_out["split"] = patient_out.get("split_label", "").fillna("")
    return sample_out, patient_out


def write_report(sample_broad: pd.DataFrame, patient_broad: pd.DataFrame, cohort_freeze: pd.DataFrame) -> None:
    merged = sample_broad.merge(
        cohort_freeze[["cohort_id", "inclusion_status", "executable_role"]],
        on="cohort_id",
        how="left",
    )
    cohort_lines: List[str] = []
    cohort_summary = (
        merged.groupby("cohort_id")
        .agg(
            samples_total=("sample_id", "size"),
            strict_samples=("response_use_strict", lambda x: int(pd.Series(x).map(truthy).sum())),
            broad_samples=("response_use_broad", lambda x: int(pd.Series(x).map(truthy).sum())),
            proxy_samples=("response_use_proxy", lambda x: int(pd.Series(x).map(truthy).sum())),
            inclusion_status=("inclusion_status", "first"),
            treatment_context=("treatment_context", "first"),
        )
        .reset_index()
        .sort_values(["samples_total", "cohort_id"], ascending=[False, True])
    )
    for _, r in cohort_summary.iterrows():
        cohort_lines.append(
            f"- `{r['cohort_id']}`: total=`{r['samples_total']}`, strict=`{r['strict_samples']}`, broad=`{r['broad_samples']}`, proxy=`{r['proxy_samples']}`, "
            f"inclusion=`{r['inclusion_status']}`, context=`{r['treatment_context']}`"
        )

    sample_evidence = sample_broad["response_evidence_level"].value_counts().sort_index()
    patient_evidence = patient_broad["response_evidence_level"].value_counts().sort_index()

    strict_cohorts = cohort_summary[cohort_summary["strict_samples"] > 0]["cohort_id"].tolist()
    broad_cohorts = cohort_summary[cohort_summary["broad_samples"] > 0]["cohort_id"].tolist()
    proxy_cohorts = cohort_summary[cohort_summary["proxy_samples"] > 0]["cohort_id"].tolist()

    support_only = cohort_summary[cohort_summary["inclusion_status"].isin(["external", "sensitivity", "spatial_only", "excluded_pending_metadata"])]["cohort_id"].tolist()

    lines: List[str] = []
    lines.append("# response_schema_report (v6.1 Task A++)")
    lines.append("")
    lines.append("## Input note")
    lines.append("- Requested `sample_metadata_master_v6_1.response_repaired.csv` / `patient_metadata_master_v6_1.response_repaired.csv` were not present.")
    lines.append("- Reconstruction used the current repaired Step1 masters plus the broad mapping dictionary.")
    lines.append("")
    lines.append("## Cohort-level strict / broad / proxy availability")
    lines.extend(cohort_lines)
    lines.append("")
    lines.append("## Evidence-level counts")
    lines.append("### Sample-level")
    for key, value in sample_evidence.items():
        lines.append(f"- `{key}`: `{int(value)}` samples")
    lines.append("### Patient-level")
    for key, value in patient_evidence.items():
        lines.append(f"- `{key}`: `{int(value)}` patients")
    lines.append("")
    lines.append("## Baseline-usable cohort groups")
    lines.append(f"- strict response baseline cohorts: `{';'.join(strict_cohorts) or 'none'}`")
    lines.append(f"- broad response baseline cohorts: `{';'.join(broad_cohorts) or 'none'}`")
    lines.append(f"- proxy response baseline cohorts: `{';'.join(proxy_cohorts) or 'none'}`")
    lines.append("")
    lines.append("## Support-only cohorts")
    lines.append(f"- support-only or non-main cohorts with any response-related information: `{';'.join(support_only) or 'none'}`")
    lines.append("")
    lines.append("## Impact on downstream branches")
    lines.append("- PD1 anchor: strict baseline remains narrow; only `TASK01` provides clean pre-treatment strict anchor labels, while `GSE206325` contributes broad study-defined/pathology response support rather than clean baseline supervision.")
    lines.append("- PD1X extension: broad schema materially expands usable supervision by admitting pathologic/study-defined response cohorts such as `GSE246613`, `GSE205506`, and study-defined RECIST-like cohorts.")
    lines.append("- HCC-specific module: HCC cohorts now retain more clinically relevant but heterogeneous response information, especially `lambrecht_hcc` durable-benefit context and `GSE206325` pathology-linked R/NR.")
    lines.append("- shared module: shared-module discovery can now test module association against strict, broad, and proxy labels separately instead of collapsing everything into one binary response.")
    lines.append("")
    lines.append("## Caution")
    lines.append("- `GSE200996` T-cell expansion labels are proxy only; they must not be written as clinical response.")
    lines.append("- `SD`, `NON CR/PD`, and `DBI` are preserved as sensitivity/support labels rather than default responders.")
    REPORT_OUT.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_policy() -> None:
    lines: List[str] = []
    lines.append("# response_usage_policy (v6.1 Task A++)")
    lines.append("")
    lines.append("## 1. Strict clinical baseline")
    lines.append("- Use only samples/patients with `response_use_strict=True`.")
    lines.append("- Valid strict labels are restricted to `A_clinical_radiologic` contexts derived from explicit `RECIST` / `mRECIST` style endpoints.")
    lines.append("- `CR`, `PR`, `OR` -> strict responder; `PD` -> strict non_responder.")
    lines.append("- `SD` is not a strict responder and must not be merged into responder in the default strict baseline.")
    lines.append("- For clean pre-treatment PD1 anchor, require in addition: executable main cohort, pre-treatment sample, and non-post-only label source.")
    lines.append("")
    lines.append("## 2. Broad treatment-response baseline")
    lines.append("- Use only samples/patients with `response_use_broad=True`.")
    lines.append("- Broad baseline may include explicit clinical radiologic labels plus `B_pathologic_or_study_defined` labels such as pathologic response and study-defined R/NR.")
    lines.append("- Broad baseline supports treatment-response association analyses, but resulting claims must be phrased as broad treatment-response association rather than strict RECIST prediction unless the contributing labels are exclusively evidence level A.")
    lines.append("")
    lines.append("## 3. Proxy biological-response baseline")
    lines.append("- Use only samples/patients with `response_use_proxy=True`.")
    lines.append("- Proxy baseline includes biological surrogate labels such as `T_cell_expansion` High/Low.")
    lines.append("- Proxy results must be described as biological proxy association and cannot be written as clinical efficacy prediction.")
    lines.append("")
    lines.append("## 4. Claim restrictions")
    lines.append("- Labels with `C_biological_proxy` evidence cannot support claims of clinical response prediction.")
    lines.append("- Pathologic/study-defined labels can support mechanistic or broad-response claims, but not be overstated as strict RECIST efficacy unless explicitly validated.")
    lines.append("- Post-only response labels can support broad/proxy association analyses but cannot define a clean pre-treatment supervised anchor.")
    lines.append("")
    lines.append("## 5. Sensitivity handling")
    lines.append("- `SD`: preserve as `stable_disease`; optionally merge into responder only in an explicitly named clinical-benefit sensitivity analysis.")
    lines.append("- `Medium` T-cell expansion: preserve as `intermediate`; do not force into responder/non_responder.")
    lines.append("- `DBI`: preserve as broad `clinical_benefit` annotation with medium confidence, but exclude from default broad baseline unless explicit metadata confirmation is added.")
    lines.append("- `NON CR/PD`: preserve as `disease_control_uncertain`; exclude from default strict and broad baselines, but allow disease-control sensitivity analyses.")
    lines.append("- `unknown`, `NE`, `UN`, `-`: exclude from all supervised baselines.")
    POLICY_OUT.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_source_record() -> None:
    df = pd.DataFrame(
        [
            {
                "cohort_id": "lambrecht_hcc",
                "url": "https://www.nature.com/articles/s41467-023-43381-1",
                "source_checked_on": pd.Timestamp.today().strftime("%Y-%m-%d"),
                "verified_fields": "response_definition;mRECIST;durable_response_context",
                "evidence_note": "Nature article states modified RECIST assessment and defines response as objective response at 3 months or disease control for >=6 months; used to interpret DBI as a durable clinical-benefit context with medium confidence.",
                "applied_to": "broad_response_schema",
            }
        ]
    )
    df.to_csv(SOURCE_OUT, sep="\t", index=False)


def main() -> None:
    sample = pd.read_csv(SAMPLE_PATH)
    patient = pd.read_csv(PATIENT_PATH)
    cohort_freeze = pd.read_csv(COHORT_FREEZE_PATH)
    patient_split = pd.read_csv(PATIENT_SPLIT_PATH)

    dictionary = build_dictionary(sample)
    dictionary.to_csv(DICT_OUT, index=False)

    sample_broad = attach_sample_schema(sample, dictionary)
    patient_broad = aggregate_patient(sample_broad, patient, patient_split)
    sample_broad, patient_broad = add_step2_compatibility_aliases(sample_broad, patient_broad, patient_split)
    sample_broad.to_csv(SAMPLE_OUT, index=False)
    patient_broad.to_csv(PATIENT_OUT, index=False)

    write_report(sample_broad, patient_broad, cohort_freeze)
    write_policy()
    write_source_record()


if __name__ == "__main__":
    main()
