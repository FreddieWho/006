from __future__ import annotations

from pathlib import Path
from typing import Dict, List, Tuple

import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
STEP1_DIR = ROOT / "results" / "v6_1" / "step1"

SAMPLE_PATH = STEP1_DIR / "sample_metadata_master_v6_1.csv"
PATIENT_PATH = STEP1_DIR / "patient_metadata_master_v6_1.csv"
FLAGS_PATH = STEP1_DIR / "treatment_context_flags.csv"
COHORT_PATH = STEP1_DIR / "cohort_registry_v6_1.csv"
SOURCE_LOCAL_PATH = STEP1_DIR / "source_records_local.tsv"
SOURCE_WEB_PATH = STEP1_DIR / "source_records_web_verified.tsv"

DICT_PATH = STEP1_DIR / "response_mapping_dictionary_v6_1.csv"
ANCHOR_REPORT_PATH = STEP1_DIR / "PD1_anchor_label_availability_report.md"


def clean_str(value: object) -> str:
    if pd.isna(value):
        return ""
    text = str(value).strip()
    if text.lower() in {"", "nan", "na", "n/a", "none", "null"}:
        return ""
    return text


def normalize_token(value: str) -> str:
    return clean_str(value).upper()


def truthy_mask(series: pd.Series) -> pd.Series:
    return series.astype(str).str.lower().isin(["true", "1", "yes", "y", "t"])


def infer_mapping(response_raw: str, response_standard: str, cohort_id: str) -> Tuple[str, str, str, str]:
    raw = normalize_token(response_raw)
    standard = normalize_token(response_standard)

    if raw in {"", "UNKNOWN", "NA", "N/A", "UN", "NE", "-", "NOT AVAILABLE"}:
        return (
            "unknown",
            "raw_missing_or_not_evaluable",
            "local_h5ad_obs",
            "Label absent or explicitly not evaluable; kept as unknown.",
        )

    if raw in {"CR", "PR", "OR"}:
        return (
            "responder",
            "objective_response_direct",
            "local_h5ad_obs;rule:objective_response",
            "Canonical objective response label mapped directly to responder.",
        )

    if raw in {"PD", "NR"}:
        return (
            "non_responder",
            "progression_or_nonresponse_direct",
            "local_h5ad_obs;rule:progression_nonresponse",
            "Canonical progression/non-response label mapped directly to non_responder.",
        )

    if raw in {"RESPONDER", "RESPONSE"}:
        return (
            "responder",
            "literal_responder_direct",
            "local_h5ad_obs;rule:literal_responder",
            "Literal responder text label mapped directly to responder.",
        )

    if raw in {"NON-RESPONDER", "NON_RESPONDER", "NONRESPONDER"}:
        return (
            "non_responder",
            "literal_nonresponder_direct",
            "local_h5ad_obs;rule:literal_nonresponder",
            "Literal non-responder text label mapped directly to non_responder.",
        )

    if raw == "SD":
        return (
            "non_responder",
            "stable_disease_binary_collapse",
            "local_h5ad_obs;rule:SD_to_non_responder",
            "Stable disease collapsed to non_responder for binary response analysis.",
        )

    if raw in {"R", "R1", "R2", "R;R1", "R;R2"}:
        if "PCR" in standard or "PATHOLOGICAL_RESPONSE" in standard:
            return (
                "responder",
                "pathological_response_direct",
                "local_h5ad_obs;rule:pathological_R_family",
                "Pathological response family label treated as responder.",
            )
        if "RECIST" in standard or "MRECIST" in standard:
            return (
                "responder",
                "response_shorthand_direct",
                "local_h5ad_obs;rule:R_to_responder",
                "Responder shorthand under RECIST-like endpoint treated as responder.",
            )
        return (
            "responder",
            "response_shorthand_direct",
            "local_h5ad_obs;rule:R_family_to_responder",
            "Responder shorthand treated as responder.",
        )

    if raw in {"HIGH", "MEDIUM", "LOW"}:
        return (
            "unknown",
            "ordinal_endpoint_not_binarized",
            "local_h5ad_obs;manual_review_required",
            "Ordinal pathological response scale preserved as unknown; no forced binarization.",
        )

    if raw == "DBI":
        return (
            "unknown",
            "ambiguous_benefit_acronym",
            "local_h5ad_obs;manual_review_required",
            "DBI meaning is not explicit enough in current source package; preserved as unknown.",
        )

    if raw == "NON CR/PD":
        return (
            "unknown",
            "recist_intermediate_state_ambiguous",
            "local_h5ad_obs;manual_review_required",
            "NON CR/PD was not forcibly collapsed because its binary interpretation is cohort-dependent.",
        )

    return (
        "unknown",
        "unrecognized_label",
        "local_h5ad_obs;manual_review_required",
        "Response label not recognized confidently; preserved as unknown.",
    )


def build_mapping_dictionary(sample: pd.DataFrame, cohort: pd.DataFrame) -> pd.DataFrame:
    rows: List[Dict[str, object]] = []
    group_cols = ["response_standard", "response_raw"]
    grouped = sample.groupby(group_cols, dropna=False)

    cohort_context_map = cohort.set_index("cohort_id")["treatment_context"].to_dict()

    for (response_standard, response_raw), sub in grouped:
        mapped_binary, rule_type, evidence_source, notes = infer_mapping(
            response_raw=response_raw,
            response_standard=response_standard,
            cohort_id="",
        )
        cohort_ids = sorted(sub["cohort_id"].astype(str).unique().tolist())
        context_values = sorted({cohort_context_map.get(cid, "unknown") for cid in cohort_ids})
        rows.append(
            {
                "response_standard": clean_str(response_standard) or "unknown",
                "response_raw": clean_str(response_raw) or "unknown",
                "mapped_response_binary": mapped_binary,
                "mapping_rule_type": rule_type,
                "n_samples": int(len(sub)),
                "n_patients": int(sub["patient_id"].astype(str).nunique()),
                "n_cohorts": int(len(cohort_ids)),
                "cohort_ids": ";".join(cohort_ids),
                "treatment_contexts": ";".join(context_values),
                "evidence_source": evidence_source,
                "notes": notes,
            }
        )

    out = pd.DataFrame(rows).sort_values(["response_standard", "response_raw"]).reset_index(drop=True)
    return out


def apply_mapping(sample: pd.DataFrame, mapping: pd.DataFrame) -> pd.DataFrame:
    lookup = {
        (clean_str(r["response_standard"]) or "unknown", clean_str(r["response_raw"]) or "unknown"): (
            r["mapped_response_binary"],
            r["mapping_rule_type"],
            r["evidence_source"],
            r["notes"],
        )
        for _, r in mapping.iterrows()
    }

    mapped_binary: List[str] = []
    mapping_rule_type: List[str] = []
    mapping_evidence: List[str] = []
    mapping_notes: List[str] = []
    binary_before: List[str] = []

    for _, row in sample.iterrows():
        key = (
            clean_str(row["response_standard"]) or "unknown",
            clean_str(row["response_raw"]) or "unknown",
        )
        repaired_binary, rule_type, evidence_source, notes = lookup.get(
            key,
            ("unknown", "missing_lookup", "manual_review_required", "No dictionary entry found."),
        )
        binary_before.append(clean_str(row["response_binary"]) or "unknown")
        mapped_binary.append(repaired_binary)
        mapping_rule_type.append(rule_type)
        mapping_evidence.append(evidence_source)
        mapping_notes.append(notes)

    sample = sample.copy()
    sample["response_binary_before_repair"] = binary_before
    sample["response_binary"] = mapped_binary
    sample["response_mapping_rule_type"] = mapping_rule_type
    sample["response_mapping_evidence_source"] = mapping_evidence
    sample["response_mapping_notes"] = mapping_notes
    return sample


def repair_patient_table(patient: pd.DataFrame, sample: pd.DataFrame) -> pd.DataFrame:
    repaired_rows: List[Dict[str, object]] = []

    for _, prow in patient.iterrows():
        cohort_id = clean_str(prow["cohort_id"])
        patient_id = clean_str(prow["patient_id"])
        sub = sample[(sample["cohort_id"] == cohort_id) & (sample["patient_id"] == patient_id)].copy()

        non_na = sub[~sub["response_binary"].isin(["unknown", "not_applicable"])].copy()
        if (non_na["response_binary"] == "responder").any():
            best_binary = "responder"
            notes = "At least one sample-level repaired responder label found."
        elif (non_na["response_binary"] == "non_responder").any():
            best_binary = "non_responder"
            notes = "Only non_responder labels remained after repair."
        elif not sub.empty and (sub["response_binary"] == "not_applicable").all():
            best_binary = "not_applicable"
            notes = "All sample-level labels are not_applicable."
        else:
            best_binary = "unknown"
            notes = "No sample-level repaired binary label available."

        response_standard_primary = "unknown"
        non_unknown_standard = [clean_str(x) for x in sub["response_standard"].tolist() if clean_str(x) and clean_str(x).lower() != "unknown"]
        if non_unknown_standard:
            response_standard_primary = pd.Series(non_unknown_standard).mode().iloc[0]

        evidence_source = ";".join(sorted({clean_str(x) for x in sub["response_mapping_evidence_source"].tolist() if clean_str(x)}))

        new_row = prow.to_dict()
        new_row["best_response_binary_before_repair"] = clean_str(prow.get("best_response_binary", "")) or "unknown"
        new_row["best_response_binary"] = best_binary
        new_row["response_standard_primary"] = response_standard_primary
        new_row["best_response_mapping_evidence_source"] = evidence_source or "none"
        new_row["best_response_mapping_notes"] = notes
        repaired_rows.append(new_row)

    return pd.DataFrame(repaired_rows)


def build_anchor_report(
    sample: pd.DataFrame,
    cohort: pd.DataFrame,
    mapping: pd.DataFrame,
) -> str:
    anchor_cohorts = cohort[cohort["treatment_context"] == "PD1_ICI_anchor"]["cohort_id"].tolist()
    lines: List[str] = []
    lines.append("# PD1_anchor_label_availability_report (v6.1 Task1)")
    lines.append("")
    lines.append("## Summary")

    anchor_samples = sample[sample["treatment_context"] == "PD1_ICI_anchor"].copy()
    lines.append(f"- Anchor cohorts in registry: `{len(anchor_cohorts)}`")
    lines.append(f"- Anchor cohorts with ingested samples: `{anchor_samples['cohort_id'].nunique()}`")
    lines.append(f"- Anchor samples total: `{len(anchor_samples)}`")
    lines.append(
        f"- Repaired anchor binary counts: `{anchor_samples['response_binary'].value_counts().to_dict()}`"
    )
    lines.append("")
    lines.append("## Cohort Assessment")

    ingested_anchor = set(anchor_samples["cohort_id"].unique().tolist())
    registry_only_anchor = [cid for cid in anchor_cohorts if cid not in ingested_anchor]

    for cid in sorted(anchor_cohorts):
        lines.append(f"### {cid}")
        if cid not in ingested_anchor:
            lines.append("- status: `registry_only`")
            lines.append("- finding: cohort is assigned to PD1 anchor in registry but no sample-level record is currently ingested.")
            lines.append("- interpretation: cannot support supervised anchor analysis until sample metadata are loaded.")
            lines.append("")
            continue

        sub = anchor_samples[anchor_samples["cohort_id"] == cid].copy()
        repaired_counts = sub["response_binary"].value_counts().to_dict()
        raw_counts = sub["response_raw"].value_counts().to_dict()
        timepoint_counts = sub["timepoint_normalized"].value_counts().to_dict()
        endpoint_counts = sub["response_standard"].value_counts().to_dict()

        issue_notes: List[str] = []
        if (sub["response_raw"].replace("unknown", "") == "").all():
            issue_notes.append("original_response_missing")
        elif (sub["response_binary_before_repair"] == "unknown").all() and (sub["response_binary"] != "unknown").any():
            issue_notes.append("response_existed_but_was_not_parsed")
        if (sub["response_binary"] == "unknown").all():
            issue_notes.append("endpoint_not_currently_supervisable")
        if set(sub["timepoint_normalized"].unique()) == {"post"}:
            issue_notes.append("post_only_not_clean_baseline_anchor")

        lines.append(f"- response_raw counts: `{raw_counts}`")
        lines.append(f"- repaired response_binary counts: `{repaired_counts}`")
        lines.append(f"- response_standard counts: `{endpoint_counts}`")
        lines.append(f"- timepoint counts: `{timepoint_counts}`")
        lines.append(f"- issues: `{';'.join(issue_notes) if issue_notes else 'none'}`")

        if cid == "GSE206325":
            lines.append("- judgment: labels exist and are now parsed (`R/NR`), but all samples are `post`; this cohort is not suitable as a clean pre-treatment supervised anchor by itself.")
            lines.append("- evidence_source: `results/v6_1/step1/source_records_web_verified.tsv` + local h5ad obs.")
        elif cid == "TASK01":
            lines.append("- judgment: labels exist and are now parsed (`CR/PR/SD/PD`) under `mRECIST`; pre/on/post are all present, so pre-treatment samples can support supervised anchor analysis.")
            lines.append("- evidence_source: local h5ad obs.")
        else:
            lines.append("- judgment: see counts above.")

        lines.append("")

    lines.append("## Supervision Decision")
    if "TASK01" in ingested_anchor:
        lines.append("- `TASK01` can support supervised anchor analysis on its pre-treatment subset after patient-level split control.")
    else:
        lines.append("- No ingested anchor cohort currently supports clean supervised analysis.")

    if "GSE206325" in ingested_anchor:
        lines.append("- `GSE206325` supports response-associated analysis, but mainly in post-treatment/pathological-response context rather than clean baseline supervision.")

    if registry_only_anchor:
        lines.append(f"- Registry-only anchor cohorts still missing from sample table: `{';'.join(sorted(registry_only_anchor))}`")

    lines.append("- Overall conclusion: `PD1 anchor is partially supervisable, but not yet strong enough as a multi-cohort baseline supervision package.`")
    return "\n".join(lines) + "\n"


def main() -> None:
    sample = pd.read_csv(SAMPLE_PATH).fillna("")
    patient = pd.read_csv(PATIENT_PATH).fillna("")
    cohort = pd.read_csv(COHORT_PATH).fillna("")

    mapping = build_mapping_dictionary(sample=sample, cohort=cohort)
    mapping.to_csv(DICT_PATH, index=False)

    repaired_sample = apply_mapping(sample=sample, mapping=mapping)
    repaired_sample.to_csv(SAMPLE_PATH, index=False)

    repaired_patient = repair_patient_table(patient=patient, sample=repaired_sample)
    repaired_patient.to_csv(PATIENT_PATH, index=False)

    report_text = build_anchor_report(sample=repaired_sample, cohort=cohort, mapping=mapping)
    ANCHOR_REPORT_PATH.write_text(report_text, encoding="utf-8")


if __name__ == "__main__":
    main()
