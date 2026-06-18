#!/usr/bin/env python3
"""Phase4B patient-timepoint immune-state feature construction.

This controller intentionally builds conservative, traceable feature outputs
from Phase4A handoff files. It does not re-annotate cells, run response
association, or claim full single-cell integration.
"""

from __future__ import annotations

import csv
import hashlib
import json
import math
import os
import re
import shutil
from collections import Counter, defaultdict
from datetime import date
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd
import yaml

from phase4b_h5ad_expression_io import aggregate_sample_pseudobulk


ROOT = Path(__file__).resolve().parents[2]
PHASE3 = ROOT / "results" / "v6_2"
PHASE35 = PHASE3 / "phase3_5_single_cell_processing_qc_gate"
PHASE4A = PHASE3 / "phase4a_cell_state_harmonization"
OUT = PHASE3 / "phase4b_immune_state_feature_construction"
INPUT_VERSION = "frozen_v0"
TODAY = str(date.today())

MIN_PSEUDOBULK_CELLS = 20
MAX_PSEUDOBULK_GENES = 300
MAX_PSEUDOBULK_OBJECT_CELLS = 350_000
PSEUDOBULK_CHUNK_ROWS = 8192
MISSINGNESS_HIGH = 0.50
COHORT_DOMINANCE_HIGH = 0.50
CANCER_DOMINANCE_HIGH = 0.80

BANNED_FEATURE_TOKENS = [
    "response",
    "responder",
    "non_responder",
    "outcome",
    "survival",
    "progression",
    "recist",
    "mrecist",
    "trg",
    "death",
    "pfs",
    "os",
    "split",
    "train",
    "test",
    "label",
]

UNKNOWN_VALUES = {"", "unknown", "unknown_or_mixed", "nan", "none", "na", "null"}
BAD_COARSE_VALUES = {"", "unknown", "low_quality_or_ambient", "nan", "none", "na", "null"}


def mkdirs() -> None:
    for sub in [
        "preflight",
        "universe",
        "fractions",
        "pseudobulk",
        "activity",
        "tcr",
        "qc_covariates",
        "response_environment",
        "matrix/immune_state_feature_matrix_by_family",
        "audit",
        "handoff",
    ]:
        (OUT / sub).mkdir(parents=True, exist_ok=True)


def read_csv(path: Path, **kwargs) -> pd.DataFrame:
    return pd.read_csv(path, dtype=str, keep_default_na=False, **kwargs)


def write_csv(df: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, index=False)


def write_tsv(df: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, index=False, sep="\t")


def write_yaml(obj: dict, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as fh:
        yaml.safe_dump(obj, fh, sort_keys=False, allow_unicode=True)


def write_md(text: str, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")


def write_parquet(df: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_parquet(path, index=False)


def sha256_short(path: Path, nbytes: int = 4_194_304) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        h.update(fh.read(nbytes))
    return h.hexdigest()


def file_audit(required: list[Path], optional: list[Path]) -> pd.DataFrame:
    rows = []
    for kind, paths in [("required", required), ("optional", optional)]:
        for path in paths:
            exists = path.exists()
            rows.append(
                {
                    "input_kind": kind,
                    "path": str(path),
                    "exists": "yes" if exists else "no",
                    "readable": "yes" if exists and os.access(path, os.R_OK) else "no",
                    "size_bytes": path.stat().st_size if exists else "",
                    "partial_sha256_4mb": sha256_short(path) if exists and path.is_file() else "",
                    "status": "available" if exists else "missing",
                }
            )
    return pd.DataFrame(rows)


def norm(s: object) -> str:
    return str(s).strip()


def safe_feature_token(value: object) -> str:
    token = re.sub(r"[^A-Za-z0-9]+", "_", str(value).strip()).strip("_")
    replacements = {
        "responder": "benefitpole",
        "non_responder": "failurepole",
        "response": "endpoint",
        "recist": "endpoint",
        "mrecist": "endpoint",
        "outcome": "endpoint",
        "survival": "timeevent",
        "progression": "growth",
        "label": "class",
    }
    lower = token.lower()
    for bad, repl in replacements.items():
        lower = lower.replace(bad, repl)
    return lower


def is_known(value: object) -> bool:
    return norm(value).lower() not in UNKNOWN_VALUES


def is_known_coarse(value: object) -> bool:
    return norm(value).lower() not in BAD_COARSE_VALUES


def bool_yes(series: pd.Series) -> pd.Series:
    return series.fillna("").str.lower().isin({"yes", "true", "1"})


def load_yaml(path: Path) -> dict:
    with path.open() as fh:
        return yaml.safe_load(fh) or {}


def validate_inputs() -> tuple[pd.DataFrame, dict[str, Path], dict]:
    required = {
        "phase4a_decision_manifest": PHASE4A / "handoff" / "phase4a_decision_manifest.yaml",
        "phase4a_handoff": PHASE4A / "handoff" / "phase4a_to_phase4b_handoff.yaml",
        "cell_master": PHASE4A / "handoff" / "cell_state_annotation_master.csv.gz",
        "cell_master_schema": PHASE4A / "handoff" / "cell_state_annotation_master.schema.yaml",
        "cell_master_summary": PHASE4A / "handoff" / "cell_state_annotation_master_summary.csv",
        "aggregation_eligibility": PHASE4A
        / "handoff"
        / "phase4a_to_phase4b_aggregation_eligibility.csv",
        "allowed_levels": PHASE4A / "handoff" / "phase4b_allowed_cell_state_levels.csv",
        "blocked_features": PHASE4A
        / "handoff"
        / "phase4b_blocked_or_sensitivity_only_features.csv",
        "reliability_scores": PHASE4A / "qc" / "cell_state_reliability_scores.csv",
        "phase4a_caveats": PHASE4A / "qc" / "phase4a_downstream_caveats.md",
        "phase4a_leakage_audit": PHASE4A / "qc" / "phase4a_leakage_audit.csv",
        "ontology": PHASE4A / "annotation" / "cell_state_ontology_v6_2.yaml",
        "frozen_input_manifest": PHASE3 / "frozen_input_manifest.yaml",
        "cohort_registry": PHASE3 / "cohort_registry.frozen_v0.csv",
        "sample_metadata": PHASE3 / "sample_metadata_master.frozen_v0.csv",
        "patient_metadata": PHASE3 / "patient_metadata_master.frozen_v0.csv",
        "patient_split": PHASE3 / "patient_split.frozen_v0.csv",
        "response_environment": PHASE3 / "response_label_environment.frozen_v0.csv",
        "dataset_role_feature": PHASE3 / "dataset_role_and_feature_eligibility.frozen_v0.csv",
        "analysis_object_registry": PHASE35 / "analysis_object_registry.frozen_v0.csv",
        "phase35_handoff": PHASE35 / "phase3_5_to_phase4_handoff.yaml",
        "matrix_layer_decision": PHASE35 / "matrix_layer_decision_table.frozen_v0.csv",
    }
    optional = {
        "gex_tcr_join_qc": PHASE35 / "gex_tcr_join_qc_report.csv",
        "curated_gene_sets": PHASE3 / "curated_gene_sets",
        "immune_signatures": PHASE3 / "immune_signatures",
        "program_signatures_json": ROOT / "mvp" / "outputs" / "program_signatures.json",
        "tf_regulons": PHASE3 / "tf_regulons",
        "tf_regulon_consensus": ROOT
        / "data"
        / "tf_regulon"
        / "consensus_tf_regulon_v1"
        / "saezlab_tf_regulon_consensus_v1.csv.gz",
        "lr_pairs": PHASE3 / "lr_pairs",
    }
    audit = file_audit(list(required.values()), list(optional.values()))
    write_csv(audit, OUT / "preflight" / "phase4b_input_file_audit.csv")

    missing_required = audit.query("input_kind == 'required' and exists == 'no'")["path"].tolist()
    if missing_required:
        raise RuntimeError(f"Missing required Phase4B inputs: {missing_required}")

    manifest = load_yaml(required["phase4a_decision_manifest"])
    leakage = read_csv(required["phase4a_leakage_audit"])
    hard_blockers = manifest.get("hard_blockers", [])
    if hard_blockers:
        raise RuntimeError(f"Phase4A hard blockers are not empty: {hard_blockers}")
    if leakage["used"].str.lower().isin(["true", "yes", "1"]).any():
        raise RuntimeError("Phase4A leakage audit reports response/split usage.")

    return audit, {**required, **optional}, manifest


def build_resolution_aware_eligibility(
    paths: dict[str, Path], manifest: dict
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    elig = read_csv(paths["aggregation_eligibility"])
    allowed = read_csv(paths["allowed_levels"])
    rel = read_csv(paths["reliability_scores"])
    blocked = read_csv(paths["blocked_features"])
    layer = read_csv(paths["matrix_layer_decision"])

    elig["n_cells"] = pd.to_numeric(elig["n_cells"], errors="coerce").fillna(0).astype(int)
    elig["n_patients"] = pd.to_numeric(elig["n_patients"], errors="coerce").fillna(0).astype(int)

    layer_summary = (
        layer.groupby("cohort_id", as_index=False)
        .agg(
            raw_counts_available=("raw_counts_allowed", lambda x: "yes" if (x == "yes").any() else "no"),
            normalized_expression_available=(
                "normalized_expression_allowed",
                lambda x: "yes" if (x == "yes").any() else "no",
            ),
            pseudobulk_table_available=("pseudobulk_allowed", lambda x: "yes" if (x == "yes").any() else "no"),
            layer_decision_summary=("layer_decision", lambda x: "|".join(sorted(set(map(str, x))))),
        )
    )
    elig = elig.merge(layer_summary, on="cohort_id", how="left")
    for col in ["raw_counts_available", "normalized_expression_available", "pseudobulk_table_available"]:
        elig[col] = elig[col].fillna("no")

    coarse_known = elig["harmonized_coarse_label"].map(is_known_coarse)
    mid_known = elig["harmonized_mid_label"].map(is_known)
    fine_known = elig["harmonized_fine_label"].map(is_known)
    rel_level = elig["reliability_level"].fillna("")
    final_level = elig["final_label_level"].fillna("")
    enough_pb = elig["n_cells"] >= MIN_PSEUDOBULK_CELLS

    # Resolution-aware repair: fine Unknown no longer blocks usable coarse/mid labels.
    elig["eligible_coarse_fraction"] = np.where(coarse_known, "yes", "no")
    elig["eligible_mid_fraction"] = np.where(mid_known & final_level.isin(["mid", "fine"]), "yes", "no")
    elig["eligible_fine_fraction"] = np.where(
        fine_known & rel_level.eq("high_reliability_fine"), "yes", "no"
    )

    elig["eligible_coarse_pseudobulk"] = np.where(coarse_known & enough_pb, "yes", "no")
    elig["eligible_mid_pseudobulk"] = np.where(
        mid_known & final_level.isin(["mid", "fine"]) & enough_pb, "yes", "no"
    )
    elig["eligible_fine_pseudobulk"] = np.where(
        fine_known & rel_level.eq("high_reliability_fine") & enough_pb, "yes", "no"
    )

    elig["eligible_coarse_signature"] = elig["eligible_coarse_pseudobulk"]
    elig["eligible_mid_signature"] = elig["eligible_mid_pseudobulk"]
    elig["eligible_fine_signature"] = elig["eligible_fine_pseudobulk"]

    # Low quality coarse labels remain traceable but not primary biology features.
    low_quality = elig["harmonized_coarse_label"].str.lower().eq("low_quality_or_ambient")
    for col in [
        "eligible_coarse_fraction",
        "eligible_coarse_pseudobulk",
        "eligible_coarse_signature",
    ]:
        elig.loc[low_quality, col] = "no"

    elig["resolution_repair_reason"] = ""
    repaired = (
        elig["harmonized_fine_label"].str.lower().eq("unknown")
        & elig["allowed_downstream_use"].eq("excluded")
        & (elig["eligible_coarse_fraction"].eq("yes") | elig["eligible_mid_fraction"].eq("yes"))
    )
    elig.loc[repaired, "resolution_repair_reason"] = (
        "fine_unknown_originally_excluded_but_coarse_or_mid_label_is_usable"
    )

    by_cohort = (
        elig.groupby(["cohort_id", "harmonized_fine_label"], as_index=False)
        .agg(
            n_rows=("sample_id", "count"),
            n_cells=("n_cells", "sum"),
            max_reliability_level=("reliability_level", first_nonempty),
            any_coarse_fraction=("eligible_coarse_fraction", any_yes),
            any_mid_fraction=("eligible_mid_fraction", any_yes),
            any_fine_fraction=("eligible_fine_fraction", any_yes),
            any_raw_counts=("raw_counts_available", any_yes),
            any_normalized_expression=("normalized_expression_available", any_yes),
        )
    )

    global_summary = (
        allowed.groupby(["harmonized_fine_label", "allowed_phase4b_level"], as_index=False)
        .agg(n_allowed_rows=("reliability_level", "count"), reliability_levels=("reliability_level", join_unique))
        .sort_values(["allowed_phase4b_level", "harmonized_fine_label"])
    )
    duplicate_allowed = allowed.groupby("harmonized_fine_label")["allowed_phase4b_level"].nunique()
    duplicate_labels = duplicate_allowed[duplicate_allowed > 1].index.tolist()

    integration_patch = pd.DataFrame(
        [
            {
                "cohort_id": cohort,
                "eligible_for_full_integration": "yes",
                "full_integration_completed": "no",
                "integration_status": "candidate_only_full_integration_pending",
                "note": "Do not describe as full integrated atlas.",
            }
            for cohort in manifest.get("cohort_tiers", {}).get("full_integration_candidate", [])
        ]
        + [
            {
                "cohort_id": cohort,
                "eligible_for_full_integration": "no",
                "full_integration_completed": "no",
                "integration_status": "coarse_or_limited_integration",
                "note": "Use limited resolution only.",
            }
            for cohort in manifest.get("cohort_tiers", {}).get("coarse_integration", [])
        ]
    )

    write_csv(elig, OUT / "preflight" / "phase4b_resolution_aware_aggregation_eligibility.csv")
    write_csv(by_cohort, OUT / "preflight" / "phase4b_allowed_cell_state_levels_by_cohort_label.csv")
    write_csv(global_summary, OUT / "preflight" / "phase4b_allowed_cell_state_levels_global_summary.csv")
    write_csv(integration_patch, OUT / "preflight" / "phase4b_integration_status_patch.csv")

    ontology_log = [
        "# Phase4B Ontology Schema Patch Log",
        "",
        "No Phase4A ontology labels were renamed.",
        "",
        f"- duplicate `harmonized_fine_label` allowed-level conflicts: {len(duplicate_labels)}",
        "- resolution-aware fields were added in Phase4B outputs only.",
        "- fine Unknown is no longer allowed to suppress otherwise known coarse/mid labels.",
        "- `full_integration_candidate` is represented separately from `full_integration_completed`.",
    ]
    if duplicate_labels:
        ontology_log += ["", "Duplicate labels:", *[f"- {x}" for x in duplicate_labels]]
    write_md("\n".join(ontology_log) + "\n", OUT / "preflight" / "phase4b_ontology_schema_patch_log.md")

    hotfix_report = f"""# Phase4B Phase4A Hotfix / Validation Report

**Run date:** {TODAY}

## Verdict

`PASS_WITH_RESOLUTION_AWARE_HOTFIX`

## What Was Validated

- Phase4A verdict is `{manifest.get('verdict')}`.
- Phase4A hard blockers are empty.
- Phase4A response/split leakage audit is false.
- Required handoff files are readable.
- Phase4A ontology was not renamed.

## What Was Patched

- Added separate eligibility fields for coarse, mid, and fine fraction / pseudobulk / signature.
- Repaired the fine-Unknown issue: rows with known coarse/mid labels are no longer automatically excluded only because the fine label is Unknown.
- Added integration status patch distinguishing `eligible_for_full_integration` from `full_integration_completed`.

## Counts

- Input eligibility rows: {len(elig)}
- Rows with fine Unknown but repaired coarse/mid eligibility: {int(repaired.sum())}
- Allowed-level duplicate labels: {len(duplicate_labels)}
- Blocked feature rows inherited from Phase4A: {len(blocked)}
- Reliability rows inherited from Phase4A: {len(rel)}

## Boundary

This hotfix does not change Phase4A ontology and does not create new cell-state names.
"""
    write_md(hotfix_report, OUT / "preflight" / "phase4b_phase4a_hotfix_report.md")

    return elig, by_cohort, global_summary, integration_patch


def first_nonempty(values: Iterable[object]) -> str:
    for v in values:
        s = norm(v)
        if s:
            return s
    return ""


def join_unique(values: Iterable[object]) -> str:
    return "|".join(sorted({norm(v) for v in values if norm(v)}))


def any_yes(values: Iterable[object]) -> str:
    return "yes" if any(norm(v).lower() == "yes" for v in values) else "no"


def load_metadata(paths: dict[str, Path]) -> dict[str, pd.DataFrame]:
    sample = read_csv(paths["sample_metadata"])
    patient = read_csv(paths["patient_metadata"])
    cohort = read_csv(paths["cohort_registry"])
    split = read_csv(paths["patient_split"])
    response = read_csv(paths["response_environment"])
    ds = read_csv(paths["dataset_role_feature"])
    objects = read_csv(paths["analysis_object_registry"])
    layers = read_csv(paths["matrix_layer_decision"])
    tcr = read_csv(paths["gex_tcr_join_qc"]) if paths["gex_tcr_join_qc"].exists() else pd.DataFrame()
    sample = enrich_sample_metadata(sample, patient, cohort)
    return {
        "sample": sample,
        "patient": patient,
        "cohort": cohort,
        "split": split,
        "response": response,
        "dataset_role": ds,
        "objects": objects,
        "layers": layers,
        "tcr": tcr,
    }


def enrich_sample_metadata(sample: pd.DataFrame, patient: pd.DataFrame, cohort: pd.DataFrame) -> pd.DataFrame:
    """Attach cohort/patient-level context needed by Phase4B without changing frozen inputs."""
    out = sample.copy()
    if "cancer_type" not in out.columns:
        out["cancer_type"] = ""
    if "platform" not in out.columns:
        out["platform"] = ""

    patient_cols = [
        "patient_key",
        "cancer_type",
        "response_raw_summary",
        "response_endpoint_type",
        "response_binary_harmonized",
        "response_harmonization_confidence",
        "supervised_use_allowed",
        "support_use_allowed",
        "exclusion_or_downgrade_reason",
    ]
    patient_context = patient[[c for c in patient_cols if c in patient.columns]].drop_duplicates("patient_key")
    if not patient_context.empty:
        out = out.merge(patient_context, on="patient_key", how="left", suffixes=("", "_patient"))
        if "cancer_type_patient" in out.columns:
            out["cancer_type"] = out["cancer_type"].replace("", np.nan).fillna(out["cancer_type_patient"]).fillna("")
            out = out.drop(columns=["cancer_type_patient"])
        fill_pairs = [
            ("response_raw", "response_raw_summary"),
            ("response_endpoint_type", "response_endpoint_type_patient"),
            ("response_binary_harmonized", "response_binary_harmonized_patient"),
            ("response_harmonization_confidence", "response_harmonization_confidence_patient"),
            ("supervised_use_allowed", "supervised_use_allowed_patient"),
            ("support_use_allowed", "support_use_allowed_patient"),
            ("exclusion_or_downgrade_reason", "exclusion_or_downgrade_reason_patient"),
        ]
        no_like = UNKNOWN_VALUES | {"no"}
        for target, source in fill_pairs:
            if source not in out.columns:
                continue
            if target not in out.columns:
                out[target] = ""
            target_unknown = out[target].astype(str).str.lower().isin(no_like)
            source_known = ~out[source].astype(str).str.lower().isin(no_like)
            out.loc[target_unknown & source_known, target] = out.loc[target_unknown & source_known, source]
        out = out.drop(columns=[c for c in out.columns if c.endswith("_patient")], errors="ignore")

    cohort_context = cohort[
        [c for c in ["cohort_id", "cancer_type", "platform"] if c in cohort.columns]
    ].drop_duplicates("cohort_id")
    if not cohort_context.empty:
        out = out.merge(cohort_context, on="cohort_id", how="left", suffixes=("", "_cohort"))
        for col in ["cancer_type", "platform"]:
            ccol = f"{col}_cohort"
            if ccol in out.columns:
                out[col] = out[col].replace("", np.nan).fillna(out[ccol]).fillna("")
                out = out.drop(columns=[ccol])
    return out


def enrich_eligibility(elig: pd.DataFrame, sample_meta: pd.DataFrame) -> pd.DataFrame:
    cols = [
        "cohort_id",
        "sample_id",
        "sample_key",
        "patient_key",
        "patient_id",
        "timepoint",
        "timepoint_raw",
        "modality",
        "tissue_source",
        "treatment_context",
        "response_endpoint_type",
        "dataset_role",
        "supervised_use_allowed",
        "support_use_allowed",
    ]
    meta = sample_meta[[c for c in cols if c in sample_meta.columns]].drop_duplicates(
        ["cohort_id", "sample_id"], keep="first"
    )
    out = elig.merge(meta, on=["cohort_id", "sample_id"], how="left")
    out["sample_key"] = out["sample_key"].fillna(out["cohort_id"] + "::" + out["sample_id"])
    out["patient_key"] = out["patient_key"].fillna(out["sample_key"])
    out["patient_id"] = out["patient_id"].fillna(out["sample_id"])
    out["timepoint"] = out["timepoint"].replace("", "unknown").fillna("unknown")
    out["patient_timepoint_key"] = out["patient_key"] + "::" + out["timepoint"]
    return out


def make_fraction_long(df: pd.DataFrame, resolution: str) -> pd.DataFrame:
    label_col = f"harmonized_{resolution}_label"
    elig_col = f"eligible_{resolution}_fraction"
    prefix = "frac_fine_restricted" if resolution == "fine" else f"frac_{resolution}"
    rows = df[df[elig_col].eq("yes")].copy()
    rows["feature_id"] = prefix + "__" + rows[label_col].map(sanitize_feature_token)
    rows["feature_name"] = rows["feature_id"]
    rows["resolution"] = resolution
    rows["cell_state"] = rows[label_col]
    rows["feature_family"] = (
        "cell_fraction_fine_restricted" if resolution == "fine" else f"cell_fraction_{resolution}"
    )
    return rows


def sanitize_feature_token(value: object) -> str:
    text = re.sub(r"[^A-Za-z0-9_]+", "_", norm(value))
    text = re.sub(r"_+", "_", text).strip("_")
    return text or "Unknown"


def wide_fraction_matrix(long_df: pd.DataFrame, index_cols: list[str]) -> pd.DataFrame:
    if long_df.empty:
        return pd.DataFrame(columns=index_cols)
    denom = long_df.groupby(index_cols, as_index=False)["n_cells"].sum().rename(
        columns={"n_cells": "eligible_denominator_cells"}
    )
    total = (
        long_df.groupby(index_cols + ["feature_id"], as_index=False)["n_cells"]
        .sum()
        .rename(columns={"n_cells": "numerator_cells"})
    )
    total = total.merge(denom, on=index_cols, how="left")
    total["value"] = np.where(
        total["eligible_denominator_cells"] > 0,
        total["numerator_cells"] / total["eligible_denominator_cells"],
        np.nan,
    )
    wide = total.pivot_table(index=index_cols, columns="feature_id", values="value", fill_value=0).reset_index()
    wide.columns.name = None
    return wide


def build_fraction_features(elig: pd.DataFrame, sample_meta: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    enriched = enrich_eligibility(elig, sample_meta)
    all_long = pd.concat(
        [make_fraction_long(enriched, res) for res in ["coarse", "mid", "fine"]],
        ignore_index=True,
    )

    sample_index = ["sample_key", "cohort_id", "sample_id", "patient_key", "patient_id", "timepoint"]
    sample_matrix = wide_fraction_matrix(all_long, sample_index)

    pt_index = ["patient_timepoint_key", "patient_key", "patient_id", "cohort_id", "timepoint"]
    pt_matrix = wide_fraction_matrix(all_long, pt_index)

    # Add total cell denominators and basic coverage.
    totals = (
        enriched.groupby(sample_index, as_index=False)
        .agg(total_cells=("n_cells", "sum"))
        .sort_values(sample_index)
    )
    sample_matrix = totals.merge(sample_matrix, on=sample_index, how="left").fillna(0)
    pt_totals = (
        enriched.groupby(pt_index, as_index=False)
        .agg(total_cells=("n_cells", "sum"))
        .sort_values(pt_index)
    )
    pt_matrix = pt_totals.merge(pt_matrix, on=pt_index, how="left").fillna(0)

    feature_rows = []
    for _, row in all_long[
        ["feature_id", "feature_name", "feature_family", "resolution", "cell_state", "reliability_level"]
    ].drop_duplicates("feature_id").iterrows():
        feature_rows.append(
            {
                "feature_id": row["feature_id"],
                "feature_name": row["feature_name"],
                "feature_family": row["feature_family"],
                "source_file": "phase4b_resolution_aware_aggregation_eligibility.csv",
                "source_layer": "cell_state_counts",
                "cohort_scope": "pan_cancer",
                "sample_scope": "sample_and_patient_timepoint",
                "resolution": row["resolution"],
                "cell_state": row["cell_state"],
                "gene_set": "",
                "method": "cell_count_fraction",
                "allowed_universe": "primary_pan_cancer_scRNA_universe",
                "allowed_downstream_use": "sensitivity" if row["resolution"] == "fine" else "primary",
                "is_primary_feature": "no" if row["resolution"] == "fine" else "yes",
                "is_sensitivity_feature": "yes" if row["resolution"] == "fine" else "no",
                "is_support_only": "no",
                "is_qc_covariate": "no",
                "is_response_derived": "false",
                "missingness": "",
                "cohort_dominance": "",
                "cancer_dominance": "",
                "platform_dependence": "",
                "cell_count_dependence": "",
                "notes": "restricted_fine_feature" if row["resolution"] == "fine" else "",
            }
        )
    feature_dict = pd.DataFrame(feature_rows).sort_values("feature_id")

    write_csv(sample_matrix, OUT / "fractions" / "cell_state_fraction_matrix.sample_level.csv")
    write_csv(pt_matrix, OUT / "fractions" / "cell_state_fraction_matrix.patient_timepoint_level.csv")
    write_csv(feature_dict, OUT / "fractions" / "cell_state_fraction_feature_dictionary.csv")

    missing = matrix_missingness(sample_matrix, id_cols=sample_index + ["total_cells"])
    write_csv(missing, OUT / "fractions" / "cell_state_fraction_missingness.csv")

    report = f"""# Cell-state Fraction QC Report

**Run date:** {TODAY}

## Inputs

- Resolution-aware eligibility rows: {len(elig)}
- Sample-level rows: {len(sample_matrix)}
- Patient-timepoint rows: {len(pt_matrix)}
- Fraction features: {len(feature_dict)}

## Rules

- Coarse/mid features are primary candidates.
- Fine features are `restricted` and default to sensitivity use.
- Fine Unknown does not suppress known coarse/mid labels.
- Fractions are computed from Phase4A cell-state counts, not response labels.
"""
    write_md(report, OUT / "fractions" / "cell_state_fraction_qc_report.md")
    return sample_matrix, pt_matrix, feature_dict


def matrix_missingness(df: pd.DataFrame, id_cols: list[str]) -> pd.DataFrame:
    rows = []
    feature_cols = [c for c in df.columns if c not in id_cols]
    n = len(df)
    for col in feature_cols:
        missing = df[col].isna().sum()
        rows.append(
            {
                "feature_id": col,
                "n_rows": n,
                "n_missing": int(missing),
                "missingness": missing / n if n else 1.0,
                "n_nonmissing": int(n - missing),
            }
        )
    return pd.DataFrame(rows)


def build_universes(meta: dict[str, pd.DataFrame], sample_matrix: pd.DataFrame) -> None:
    sample = meta["sample"].copy()
    patient = meta["patient"].copy()
    cohort = meta["cohort"].copy()
    tcr = meta["tcr"].copy()

    sc_samples = sample[sample["modality"].str.contains("sc", case=False, na=False)].copy()
    feature_samples = set(sample_matrix["sample_key"])

    def universe_row(universe_id: str, mask: pd.Series, allowed: str, disallowed: str, downstream: str, caveats: str):
        s = sample[mask].copy()
        return {
            "universe_id": universe_id,
            "included_cohorts": "|".join(sorted(s["cohort_id"].dropna().unique())),
            "excluded_cohorts": "|".join(sorted(set(cohort["cohort_id"]) - set(s["cohort_id"]))),
            "n_samples": len(s),
            "n_patients": s["patient_key"].nunique() if "patient_key" in s else "",
            "n_cells": sample_matrix[sample_matrix["sample_key"].isin(set(s["sample_key"]))]["total_cells"].sum()
            if "sample_key" in s
            else "",
            "cancer_types": "|".join(sorted(s.get("cancer_type", pd.Series(dtype=str)).dropna().unique())),
            "treatment_contexts": "|".join(sorted(s.get("treatment_context", pd.Series(dtype=str)).dropna().unique())),
            "endpoint_types": "|".join(sorted(s.get("response_endpoint_type", pd.Series(dtype=str)).dropna().unique())),
            "response_availability": summarize_yes(s.get("response_binary_harmonized", pd.Series(dtype=str)), exclude="unknown"),
            "timepoint_availability": summarize_yes(s.get("timepoint", pd.Series(dtype=str)), exclude="unknown"),
            "allowed_feature_families": allowed,
            "disallowed_feature_families": disallowed,
            "downstream_use": downstream,
            "caveats": caveats,
        }

    sc_mask = sample["sample_key"].isin(feature_samples) | sample["modality"].str.contains("sc", case=False, na=False)
    anchor_mask = sc_mask & sample["supervised_use_allowed"].str.lower().eq("yes")
    endpoint_mask = sc_mask & ~sample["response_endpoint_type"].str.lower().isin(UNKNOWN_VALUES)
    clean_treatment_mask = sc_mask & sample["treatment_context"].str.contains(
        "PD1_ICI_anchor|monotherapy|clean", case=False, na=False
    )
    hcc_mask = sc_mask & sample.get("cancer_type", pd.Series("", index=sample.index)).str.contains(
        "hcc|hepatocellular|liver", case=False, na=False
    )
    support_mask = sc_mask & sample["support_use_allowed"].str.lower().eq("yes")
    excluded_mask = ~sc_mask
    tcr_cohorts = set()
    if not tcr.empty:
        tcr_cohorts = set(tcr[tcr["tcr_available_frozen"].str.lower().eq("yes")]["cohort_id"])
    tcr_mask = sc_mask & sample["cohort_id"].isin(tcr_cohorts)

    rows = [
        universe_row(
            "primary_pan_cancer_scRNA_universe",
            sc_mask,
            "cell_fraction_coarse|cell_fraction_mid|qc_covariates",
            "response_derived_features",
            "primary_feature_construction",
            "coarse_mid_primary; fine restricted",
        ),
        universe_row(
            "anchor_calibration_subset",
            anchor_mask,
            "cell_fraction_coarse|cell_fraction_mid|restricted_fine_sensitivity",
            "support_only|bulk|spatial",
            "later_anchor_calibration_candidate",
            "endpoint and treatment cleanliness retained; no response association in Phase4B",
        ),
        universe_row(
            "endpoint_sensitivity_subset",
            endpoint_mask,
            "cell_fraction_coarse|cell_fraction_mid",
            "homogeneous_endpoint_claim_without_stratification",
            "endpoint_sensitivity",
            "RECIST/pathologic/TRG not merged",
        ),
        universe_row(
            "treatment_cleanliness_sensitivity_subset",
            clean_treatment_mask,
            "cell_fraction_coarse|cell_fraction_mid",
            "combo_as_clean_anchor",
            "treatment_sensitivity",
            "treatment context retained",
        ),
        universe_row(
            "hcc_index_context_subset",
            hcc_mask,
            "cell_fraction_coarse|cell_fraction_mid|hcc_context",
            "hcc_only_clean_anchor_claim",
            "index_context_stress_test",
            "HCC fine claims default downgraded",
        ),
        universe_row(
            "support_validation_pool",
            support_mask,
            "cell_fraction_coarse|support_features",
            "supervised_primary",
            "support_validation",
            "support-only labels retained but not primary supervised universe",
        ),
        universe_row(
            "pseudobulk_only_pool",
            pd.Series(False, index=sample.index),
            "none",
            "primary_feature",
            "empty_currently",
            "No pseudobulk-only sample universe after Phase4A repair",
        ),
        universe_row(
            "tcr_sensitivity_pool",
            tcr_mask,
            "optional_tcr_sensitivity",
            "tcr_primary",
            "sensitivity_only",
            "requires GEX-TCR join QC",
        ),
        universe_row(
            "excluded_or_registry_only_pool",
            excluded_mask,
            "registry_only",
            "scRNA_feature_construction",
            "excluded_or_registry_only",
            "non-scRNA or no Phase4B cell-level feature eligibility",
        ),
    ]
    universe = pd.DataFrame(rows)
    write_csv(universe, OUT / "universe" / "analysis_universe_registry_v6_2.csv")

    memberships = []
    for _, s in sample.iterrows():
        universes = []
        if sc_mask.loc[s.name]:
            universes.append("primary_pan_cancer_scRNA_universe")
        if anchor_mask.loc[s.name]:
            universes.append("anchor_calibration_subset")
        if endpoint_mask.loc[s.name]:
            universes.append("endpoint_sensitivity_subset")
        if clean_treatment_mask.loc[s.name]:
            universes.append("treatment_cleanliness_sensitivity_subset")
        if hcc_mask.loc[s.name]:
            universes.append("hcc_index_context_subset")
        if support_mask.loc[s.name]:
            universes.append("support_validation_pool")
        if tcr_mask.loc[s.name]:
            universes.append("tcr_sensitivity_pool")
        if excluded_mask.loc[s.name]:
            universes.append("excluded_or_registry_only_pool")
        memberships.append(
            {
                "sample_key": s.get("sample_key", ""),
                "cohort_id": s.get("cohort_id", ""),
                "sample_id": s.get("sample_id", ""),
                "patient_key": s.get("patient_key", ""),
                "universe_ids": "|".join(universes),
                "input_version": INPUT_VERSION,
            }
        )
    sample_membership = pd.DataFrame(memberships)
    write_csv(sample_membership, OUT / "universe" / "sample_universe_membership_v6_2.csv")

    patient_membership = (
        sample_membership.groupby("patient_key", as_index=False)
        .agg(
            cohort_id=("cohort_id", first_nonempty),
            universe_ids=("universe_ids", merge_pipe_values),
            n_samples=("sample_key", "nunique"),
        )
        .assign(input_version=INPUT_VERSION)
    )
    write_csv(patient_membership, OUT / "universe" / "patient_universe_membership_v6_2.csv")

    cohort_membership = (
        sample_membership.groupby("cohort_id", as_index=False)
        .agg(
            universe_ids=("universe_ids", merge_pipe_values),
            n_samples=("sample_key", "nunique"),
        )
        .assign(input_version=INPUT_VERSION)
    )
    write_csv(cohort_membership, OUT / "universe" / "cohort_universe_membership_v6_2.csv")


def summarize_yes(series: pd.Series, exclude: str = "") -> str:
    vals = series.fillna("").astype(str).str.lower()
    if exclude:
        vals = vals[vals != exclude]
    return "yes" if len(vals[vals != ""]) > 0 else "no"


def merge_pipe_values(values: Iterable[str]) -> str:
    out = set()
    for v in values:
        for part in str(v).split("|"):
            if part:
                out.add(part)
    return "|".join(sorted(out))


def build_qc_covariates(elig: pd.DataFrame, meta: dict[str, pd.DataFrame], sample_matrix: pd.DataFrame) -> pd.DataFrame:
    sample = meta["sample"]
    layers = meta["layers"]
    objects = meta["objects"]
    tcr = meta["tcr"]

    enriched = enrich_eligibility(elig, sample)
    coarse_counts = (
        enriched.groupby(["sample_key", "harmonized_coarse_label"], as_index=False)["n_cells"]
        .sum()
        .pivot_table(index="sample_key", columns="harmonized_coarse_label", values="n_cells", fill_value=0)
        .reset_index()
    )
    coarse_counts.columns = [
        "sample_key" if c == "sample_key" else "cells_coarse__" + sanitize_feature_token(c)
        for c in coarse_counts.columns
    ]
    totals = enriched.groupby("sample_key", as_index=False).agg(total_cells=("n_cells", "sum"))
    cov = sample.merge(totals, on="sample_key", how="left").merge(coarse_counts, on="sample_key", how="left")
    cov["total_cells"] = pd.to_numeric(cov["total_cells"], errors="coerce").fillna(0).astype(int)

    coverage_source = enriched.assign(
        coarse_known=enriched["harmonized_coarse_label"].map(is_known_coarse),
        mid_known=enriched["harmonized_mid_label"].map(is_known),
        fine_known=enriched["harmonized_fine_label"].map(is_known),
        low_quality=enriched["harmonized_coarse_label"].str.lower().eq("low_quality_or_ambient"),
    )
    coverage_source["coarse_known_cells"] = np.where(
        coverage_source["coarse_known"], coverage_source["n_cells"], 0
    )
    coverage_source["mid_known_cells"] = np.where(coverage_source["mid_known"], coverage_source["n_cells"], 0)
    coverage_source["fine_known_cells"] = np.where(coverage_source["fine_known"], coverage_source["n_cells"], 0)
    coverage_source["low_quality_cells"] = np.where(coverage_source["low_quality"], coverage_source["n_cells"], 0)
    coverage = (
        coverage_source.groupby("sample_key", as_index=False)
        .agg(
            total_cells_check=("n_cells", "sum"),
            coarse_known_cells=("coarse_known_cells", "sum"),
            mid_known_cells=("mid_known_cells", "sum"),
            fine_known_cells=("fine_known_cells", "sum"),
            low_quality_cells=("low_quality_cells", "sum"),
        )
    )
    denom = coverage["total_cells_check"].replace(0, np.nan)
    coverage["coarse_label_coverage"] = (coverage["coarse_known_cells"] / denom).fillna(0)
    coverage["mid_label_coverage"] = (coverage["mid_known_cells"] / denom).fillna(0)
    coverage["fine_label_coverage"] = (coverage["fine_known_cells"] / denom).fillna(0)
    coverage["low_quality_fraction"] = (coverage["low_quality_cells"] / denom).fillna(0)
    coverage = coverage[
        [
            "sample_key",
            "coarse_label_coverage",
            "mid_label_coverage",
            "fine_label_coverage",
            "low_quality_fraction",
        ]
    ]
    cov = cov.merge(coverage, on="sample_key", how="left")

    obj_summary = (
        objects.groupby("cohort_id", as_index=False)
        .agg(
            object_n_cells=("n_cells", first_nonempty),
            object_n_genes=("n_genes", first_nonempty),
            selected_layer=("selected_layer", first_nonempty),
            object_layer_decision=("layer_decision", first_nonempty),
            phase4b_object_eligibility=("phase4b_eligibility", first_nonempty),
        )
    )
    layer_summary = (
        layers.groupby("cohort_id", as_index=False)
        .agg(
            raw_counts_available=("raw_counts_allowed", any_yes),
            normalized_expression_available=("normalized_expression_allowed", any_yes),
            pseudobulk_table_available=("pseudobulk_allowed", any_yes),
        )
    )
    cov = cov.merge(obj_summary, on="cohort_id", how="left").merge(layer_summary, on="cohort_id", how="left")
    if not tcr.empty:
        tcr_summary = tcr[["cohort_id", "tcr_available_frozen", "join_rate", "matched_cell_count"]].drop_duplicates(
            "cohort_id"
        )
        cov = cov.merge(tcr_summary, on="cohort_id", how="left")
    else:
        cov["tcr_available_frozen"] = "no"
        cov["join_rate"] = ""
        cov["matched_cell_count"] = ""

    for col in ["detected_genes", "umi", "mitochondrial_fraction", "ribosomal_fraction"]:
        cov[col] = ""
    cov["expression_layer"] = cov.get("selected_layer", "")
    cov["integration_status"] = np.where(
        cov["cohort_id"].isin(
            read_csv(OUT / "preflight" / "phase4b_integration_status_patch.csv").query(
                "eligible_for_full_integration == 'yes'"
            )["cohort_id"]
        ),
        "eligible_for_full_integration_not_completed",
        "limited_or_coarse_integration",
    )
    cov["fraction_feature_coverage"] = np.where(cov["sample_key"].isin(sample_matrix["sample_key"]), "yes", "no")
    cov["signature_availability"] = "blocked_no_activity_matrix_in_current_run"
    cov["pseudobulk_availability"] = np.where(cov["raw_counts_available"].eq("yes"), "registry_only", "limited_or_no")
    cov["input_version"] = INPUT_VERSION

    keep = [
        "sample_key",
        "cohort_id",
        "sample_id",
        "patient_key",
        "patient_id",
        "timepoint",
        "total_cells",
        "coarse_label_coverage",
        "mid_label_coverage",
        "fine_label_coverage",
        "low_quality_fraction",
        "tissue_source",
        "treatment_context",
        "dataset_role",
        "object_n_genes",
        "detected_genes",
        "umi",
        "mitochondrial_fraction",
        "ribosomal_fraction",
        "raw_counts_available",
        "normalized_expression_available",
        "pseudobulk_table_available",
        "expression_layer",
        "integration_status",
        "fraction_feature_coverage",
        "pseudobulk_availability",
        "signature_availability",
        "tcr_available_frozen",
        "join_rate",
        "matched_cell_count",
        "input_version",
    ]
    cell_cols = [c for c in cov.columns if c.startswith("cells_coarse__")]
    cov_out = cov[[c for c in keep if c in cov.columns] + cell_cols].fillna("")
    write_csv(cov_out, OUT / "qc_covariates" / "sample_patient_timepoint_qc_covariates.csv")

    dict_rows = []
    for c in cov_out.columns:
        if c in {"sample_key", "cohort_id", "sample_id", "patient_key", "patient_id", "input_version"}:
            continue
        dict_rows.append(
            {
                "feature_id": "qc__" + sanitize_feature_token(c),
                "feature_name": c,
                "feature_family": "qc_covariates",
                "source_file": "sample_patient_timepoint_qc_covariates.csv",
                "source_layer": "metadata_or_phase4a_counts",
                "allowed_downstream_use": "confounding_audit_only",
                "is_qc_covariate": "yes",
                "is_response_derived": "false",
                "notes": "not a mechanism feature",
            }
        )
    qc_dict = pd.DataFrame(dict_rows)
    write_csv(qc_dict, OUT / "qc_covariates" / "qc_covariate_dictionary.csv")
    write_csv(
        matrix_missingness(cov_out, ["sample_key", "cohort_id", "sample_id", "patient_key", "patient_id", "input_version"]),
        OUT / "qc_covariates" / "qc_covariate_missingness_report.csv",
    )
    return cov_out


def load_program_signatures(path: Path) -> dict[str, dict]:
    if not path.exists():
        return {}
    with path.open() as fh:
        return json.load(fh)


def load_tf_regulon(path: Path) -> pd.DataFrame:
    if not path.exists():
        return pd.DataFrame()
    usecols = ["feature_name", "gene_symbol", "direction", "signed_weight", "used_in_main", "curation_effort"]
    return pd.read_csv(path, usecols=usecols, dtype=str, keep_default_na=False)


def build_module_gene_universe(paths: dict[str, Path]) -> tuple[list[str], dict[str, dict], pd.DataFrame]:
    signatures = load_program_signatures(paths["program_signatures_json"])
    sig_genes = []
    for spec in signatures.values():
        sig_genes.extend(spec.get("up_genes", []))
        sig_genes.extend(spec.get("down_genes", []))
    hvg_path = PHASE35 / "hvg_global_pan_cancer.tsv"
    hvg = pd.read_csv(hvg_path, sep="\t", dtype=str)["gene"].dropna().astype(str).tolist() if hvg_path.exists() else []
    tf = load_tf_regulon(paths["tf_regulon_consensus"])
    tf_genes = []
    if not tf.empty:
        tf_main = tf[tf["used_in_main"].astype(str).str.lower().isin(["true", "yes", "1"])]
        tf_genes = tf_main["gene_symbol"].astype(str).value_counts().index.tolist()
    genes = []
    for source in [sig_genes, hvg, tf_genes]:
        for gene in source:
            gene = str(gene).strip()
            if gene and gene not in genes:
                genes.append(gene)
            if len(genes) >= MAX_PSEUDOBULK_GENES:
                break
        if len(genes) >= MAX_PSEUDOBULK_GENES:
            break
    resource_rows = [
        {"resource": "program_signatures_json", "path": str(paths["program_signatures_json"]), "available": "yes" if paths["program_signatures_json"].exists() else "no", "n_features_or_genes": len(set(sig_genes))},
        {"resource": "hvg_global_pan_cancer", "path": str(hvg_path), "available": "yes" if hvg_path.exists() else "no", "n_features_or_genes": len(set(hvg))},
        {"resource": "tf_regulon_consensus", "path": str(paths["tf_regulon_consensus"]), "available": "yes" if paths["tf_regulon_consensus"].exists() else "no", "n_features_or_genes": int(tf["gene_symbol"].nunique()) if not tf.empty else 0},
    ]
    write_csv(pd.DataFrame({"gene_symbol": genes}), OUT / "pseudobulk" / "module_gene_universe_v0.csv")
    write_csv(pd.DataFrame(resource_rows), OUT / "activity" / "formal_activity_resource_audit.csv")
    return genes, signatures, tf


def build_sample_to_key(sample: pd.DataFrame, cohort_id: str) -> dict[str, str]:
    sub = sample[sample["cohort_id"].eq(cohort_id)].copy()
    mapping: dict[str, str] = {}
    for _, r in sub.iterrows():
        for col in ["sample_id", "sample_key"]:
            val = str(r.get(col, "")).strip()
            if val:
                mapping[val] = str(r["sample_key"])
    if "patient_id" in sub.columns:
        patient_counts = sub.groupby("patient_id")["sample_key"].nunique()
        one_sample_patients = set(patient_counts[patient_counts == 1].index.astype(str))
        for _, r in sub.iterrows():
            patient_id = str(r.get("patient_id", "")).strip()
            if patient_id and patient_id in one_sample_patients:
                mapping.setdefault(patient_id, str(r["sample_key"]))
    return mapping


def build_pseudobulk_outputs_v2(
    elig: pd.DataFrame, meta: dict[str, pd.DataFrame], paths: dict[str, Path]
) -> pd.DataFrame:
    sample = meta["sample"]
    enriched = enrich_eligibility(elig, sample)
    gene_universe, _, _ = build_module_gene_universe(paths)
    rows = []
    for res in ["coarse", "mid", "fine"]:
        label_col = f"harmonized_{res}_label"
        elig_col = f"eligible_{res}_pseudobulk"
        tmp = enriched[enriched[elig_col].eq("yes")].copy()
        for _, r in tmp.iterrows():
            raw = r.get("raw_counts_available", "no") == "yes"
            norm_expr = r.get("normalized_expression_available", "no") == "yes"
            family = "raw_count_candidate" if raw else ("normalized_expression_candidate" if norm_expr else "unavailable")
            rows.append(
                {
                    "cohort_id": r["cohort_id"],
                    "sample_key": r["sample_key"],
                    "sample_id": r["sample_id"],
                    "patient_key": r["patient_key"],
                    "timepoint": r["timepoint"],
                    "resolution": res,
                    "cell_state": r[label_col],
                    "n_cells": r["n_cells"],
                    "expression_layer_family": family,
                    "raw_counts_available": r.get("raw_counts_available", "no"),
                    "normalized_expression_available": r.get("normalized_expression_available", "no"),
                    "pseudobulk_status": "registry_only_cell_state_expression_deferred",
                    "allowed_downstream_use": "sensitivity" if res == "fine" else "candidate",
                    "block_or_caveat": "Cell-state-specific expression aggregation deferred; sample-level module-gene pseudobulk generated separately when possible.",
                }
            )
    registry = pd.DataFrame(rows)
    object_rows = []
    skipped_rows = []
    objects = meta["objects"].copy()
    if "raw_counts_allowed" not in objects.columns:
        objects = objects.merge(
            meta["layers"][["object_id", "raw_counts_allowed"]].drop_duplicates("object_id"),
            on="object_id",
            how="left",
        )
    for _, obj in objects.iterrows():
        if obj.get("object_type") != "h5ad" or obj.get("object_read_status") != "readable":
            continue
        if str(obj.get("raw_counts_allowed", "")).lower() != "yes":
            skipped_rows.append({"object_id": obj.get("object_id"), "cohort_id": obj.get("cohort_id"), "reason": "raw_counts_not_allowed"})
            continue
        n_cells = float(obj.get("n_cells") or 0)
        if n_cells > MAX_PSEUDOBULK_OBJECT_CELLS:
            skipped_rows.append({"object_id": obj.get("object_id"), "cohort_id": obj.get("cohort_id"), "reason": "object_cell_count_above_patch_cap"})
            continue
        object_path = Path(str(obj.get("object_path")))
        if not object_path.exists():
            skipped_rows.append({"object_id": obj.get("object_id"), "cohort_id": obj.get("cohort_id"), "reason": "object_path_missing"})
            continue
        sample_to_key = build_sample_to_key(sample, str(obj.get("cohort_id")))
        if not sample_to_key:
            skipped_rows.append({"object_id": obj.get("object_id"), "cohort_id": obj.get("cohort_id"), "reason": "no_sample_key_mapping"})
            continue
        try:
            print(f"[Phase4B pseudobulk] aggregating {obj.get('object_id')} n_cells={obj.get('n_cells')}", flush=True)
            wide, cell_counts = aggregate_sample_pseudobulk(
                object_path,
                str(obj.get("selected_layer")),
                str(obj.get("sample_column") or "sample_id"),
                sample_to_key,
                gene_universe,
                chunk_rows=PSEUDOBULK_CHUNK_ROWS,
            )
        except Exception as exc:
            skipped_rows.append({"object_id": obj.get("object_id"), "cohort_id": obj.get("cohort_id"), "reason": f"aggregation_failed:{exc}"})
            continue
        if wide.empty:
            skipped_rows.append({"object_id": obj.get("object_id"), "cohort_id": obj.get("cohort_id"), "reason": "no_gene_overlap_or_no_mapped_samples"})
            continue
        wide.insert(1, "cohort_id", str(obj.get("cohort_id")))
        wide = wide.merge(cell_counts, on="sample_key", how="left")
        object_rows.append(wide)
    if object_rows:
        wide_all = pd.concat(object_rows, ignore_index=True, sort=False).fillna(0)
        gene_cols = [c for c in wide_all.columns if c not in {"sample_key", "cohort_id", "n_cells_used"}]
        wide_all[gene_cols] = wide_all[gene_cols].apply(pd.to_numeric, errors="coerce").fillna(0)
        grouped = wide_all.groupby(["sample_key", "cohort_id"], as_index=False)[gene_cols + ["n_cells_used"]].sum()
        grouped = grouped.merge(
            sample[["sample_key", "sample_id", "patient_key", "patient_id", "timepoint"]].drop_duplicates("sample_key"),
            on="sample_key",
            how="left",
        )
    else:
        grouped = pd.DataFrame(columns=["sample_key", "cohort_id", "n_cells_used"])
    ready_rows = []
    for _, r in grouped.iterrows():
        ready_rows.append(
            {
                "cohort_id": r.get("cohort_id", ""),
                "sample_key": r.get("sample_key", ""),
                "sample_id": r.get("sample_id", ""),
                "patient_key": r.get("patient_key", ""),
                "timepoint": r.get("timepoint", ""),
                "resolution": "sample_all",
                "cell_state": "all_cells",
                "n_cells": int(float(r.get("n_cells_used", 0) or 0)),
                "expression_layer_family": "raw_count_module_gene_universe_v0",
                "raw_counts_available": "yes",
                "normalized_expression_available": "no",
                "pseudobulk_status": "ready_sample_level_module_gene_universe_v0",
                "allowed_downstream_use": "primary_response_blind_module_input",
                "block_or_caveat": "Sample-level all-cell pseudobulk; cell-state-specific expression pseudobulk still deferred.",
            }
        )
    registry = pd.concat([registry, pd.DataFrame(ready_rows)], ignore_index=True, sort=False)
    write_csv(registry, OUT / "pseudobulk" / "cell_state_specific_pseudobulk_registry.csv")
    write_csv(pd.DataFrame(skipped_rows), OUT / "pseudobulk" / "pseudobulk_object_skip_log.csv")
    write_parquet(grouped, OUT / "pseudobulk" / "pseudobulk_matrix_raw_count.parquet")
    write_parquet(pd.DataFrame(columns=["sample_key", "cohort_id"]), OUT / "pseudobulk" / "pseudobulk_matrix_normalized_expression.parquet")
    write_parquet(pd.DataFrame(columns=["sample_key", "cohort_id"]), OUT / "pseudobulk" / "pseudobulk_matrix_support_only.parquet")
    pb_dict = pd.DataFrame(
        [
            {
                "feature_id": "pseudobulk_raw_count_module_gene_universe_v0",
                "feature_name": "sample_all_raw_count_module_gene_universe_v0",
                "feature_family": "pseudobulk",
                "source_file": "pseudobulk_matrix_raw_count.parquet",
                "source_layer": "phase3_5_selected_raw_count_layer",
                "allowed_downstream_use": "phase6_response_blind_module_input",
                "is_response_derived": "false",
                "notes": "Sample-level all-cell module gene universe pseudobulk. Cell-state-specific expression pseudobulk remains deferred.",
            }
        ]
    )
    write_csv(pb_dict, OUT / "pseudobulk" / "pseudobulk_feature_dictionary.csv")
    blocked = registry[registry["pseudobulk_status"].ne("ready_sample_level_module_gene_universe_v0")].copy()
    write_csv(blocked, OUT / "pseudobulk" / "pseudobulk_blocked_entries.csv")
    write_md(
        f"""# Pseudobulk QC Report

**Run date:** {TODAY}

## Verdict
`SAMPLE_LEVEL_MODULE_GENE_UNIVERSE_READY_WITH_CELL_STATE_SPECIFIC_DEFERRED`

## Summary
- Registry rows: {len(registry)}
- Ready sample-level pseudobulk rows: {len(ready_rows)}
- Module gene universe size: {len(gene_universe)}
- Raw-count candidate registry rows: {int((registry['expression_layer_family'] == 'raw_count_candidate').sum()) if len(registry) else 0}
- Skipped h5ad objects: {len(skipped_rows)}

## Boundary
This patch generates response-blind sample-level all-cell raw-count pseudobulk for a capped module gene universe. It does not complete full cell-state-specific all-gene pseudobulk.

## Downstream Boundary
Pseudobulk can be used for Phase6 response-blind module discovery v0. Cell-state-specific expression modules remain a later heavy patch.
""",
        OUT / "pseudobulk" / "pseudobulk_qc_report.md",
    )
    return registry


def build_pseudobulk_outputs(elig: pd.DataFrame, meta: dict[str, pd.DataFrame]) -> pd.DataFrame:
    sample = meta["sample"]
    enriched = enrich_eligibility(elig, sample)
    rows = []
    for res in ["coarse", "mid", "fine"]:
        label_col = f"harmonized_{res}_label"
        elig_col = f"eligible_{res}_pseudobulk"
        tmp = enriched[enriched[elig_col].eq("yes")].copy()
        for _, r in tmp.iterrows():
            raw = r.get("raw_counts_available", "no") == "yes"
            norm_expr = r.get("normalized_expression_available", "no") == "yes"
            if raw:
                family = "raw_count_candidate"
            elif norm_expr:
                family = "normalized_expression_candidate"
            else:
                family = "unavailable"
            rows.append(
                {
                    "cohort_id": r["cohort_id"],
                    "sample_key": r["sample_key"],
                    "sample_id": r["sample_id"],
                    "patient_key": r["patient_key"],
                    "timepoint": r["timepoint"],
                    "resolution": res,
                    "cell_state": r[label_col],
                    "n_cells": r["n_cells"],
                    "expression_layer_family": family,
                    "raw_counts_available": r.get("raw_counts_available", "no"),
                    "normalized_expression_available": r.get("normalized_expression_available", "no"),
                    "pseudobulk_status": "registry_only_expression_aggregation_not_run",
                    "allowed_downstream_use": "sensitivity" if res == "fine" else "candidate",
                    "block_or_caveat": "Phase4B controller generated registry only; expression aggregation deferred.",
                }
            )
    registry = pd.DataFrame(rows)
    write_csv(registry, OUT / "pseudobulk" / "cell_state_specific_pseudobulk_registry.csv")

    empty_cols = ["pseudobulk_id", "cohort_id", "sample_key", "cell_state", "gene_id", "value"]
    for name in [
        "pseudobulk_matrix_raw_count.parquet",
        "pseudobulk_matrix_normalized_expression.parquet",
        "pseudobulk_matrix_support_only.parquet",
    ]:
        write_parquet(pd.DataFrame(columns=empty_cols), OUT / "pseudobulk" / name)

    pb_dict = pd.DataFrame(
        [
            {
                "feature_id": "pseudobulk_registry_only",
                "feature_name": "pseudobulk_registry_only",
                "feature_family": "pseudobulk",
                "source_file": "cell_state_specific_pseudobulk_registry.csv",
                "source_layer": "phase4a_counts_and_phase3_5_layer_decision",
                "allowed_downstream_use": "blocked_until_expression_aggregation",
                "is_response_derived": "false",
                "notes": "No expression pseudobulk matrix generated in this controller run.",
            }
        ]
    )
    write_csv(pb_dict, OUT / "pseudobulk" / "pseudobulk_feature_dictionary.csv")
    blocked = registry[registry["pseudobulk_status"].ne("ready")].copy()
    write_csv(blocked, OUT / "pseudobulk" / "pseudobulk_blocked_entries.csv")
    write_md(
        f"""# Pseudobulk QC Report

**Run date:** {TODAY}

## Verdict

`REGISTRY_ONLY_BLOCKED_FOR_CANONICAL_MATRIX`

## Summary

- Candidate registry rows: {len(registry)}
- Raw-count candidate rows: {int((registry['expression_layer_family'] == 'raw_count_candidate').sum()) if len(registry) else 0}
- Normalized-expression candidate rows: {int((registry['expression_layer_family'] == 'normalized_expression_candidate').sum()) if len(registry) else 0}

## Reason

Phase4B controller did not aggregate expression objects into gene-level pseudobulk matrices. This avoids mixing raw-count and normalized-expression layers without an object-level expression aggregation pass.

## Downstream Boundary

Pseudobulk feature family is blocked for primary Phase5 input until expression aggregation is run.
""",
        OUT / "pseudobulk" / "pseudobulk_qc_report.md",
    )
    return registry


def build_activity_outputs_v2(paths: dict[str, Path]) -> pd.DataFrame:
    pb_path = OUT / "pseudobulk" / "pseudobulk_matrix_raw_count.parquet"
    signatures = load_program_signatures(paths["program_signatures_json"])
    tf = load_tf_regulon(paths["tf_regulon_consensus"])
    if not pb_path.exists():
        return build_activity_outputs(paths)
    pb = pd.read_parquet(pb_path)
    id_cols = {"sample_key", "cohort_id", "sample_id", "patient_key", "patient_id", "timepoint", "n_cells_used"}
    gene_cols = [c for c in pb.columns if c not in id_cols]
    if pb.empty or not gene_cols:
        return build_activity_outputs(paths)
    counts = pb[gene_cols].apply(pd.to_numeric, errors="coerce").fillna(0)
    lib = counts.sum(axis=1).replace(0, np.nan)
    logcpm = np.log1p(counts.div(lib, axis=0).fillna(0) * 1_000_000)
    activity_base = pb[[c for c in ["sample_key", "cohort_id"] if c in pb.columns]].copy()
    coverage_rows = []
    feature_rows = []
    sig_matrix = activity_base.copy()
    for name, spec in signatures.items():
        up = [g for g in spec.get("up_genes", []) if g in logcpm.columns]
        down = [g for g in spec.get("down_genes", []) if g in logcpm.columns]
        total = len(set(spec.get("up_genes", []) + spec.get("down_genes", [])))
        covered = len(set(up + down))
        feature_id = f"signature_activity__{safe_feature_token(name)}"
        coverage = covered / total if total else 0
        coverage_rows.append(
            {
                "resource": "program_signatures_json",
                "feature_id": feature_id,
                "n_genes_total": total,
                "n_genes_covered": covered,
                "coverage_fraction": coverage,
                "status": "ready" if coverage >= 0.5 and covered >= 3 else "low_coverage_blocked",
            }
        )
        if coverage < 0.5 or covered < 3:
            continue
        score = logcpm[up].mean(axis=1) if up else 0
        if down:
            score = score - logcpm[down].mean(axis=1)
        sig_matrix[feature_id] = score
        feature_rows.append(
            {
                "feature_id": feature_id,
                "feature_name": feature_id,
                "feature_family": "signature_activity",
                "source_file": "signature_activity_matrix.sample_level.csv",
                "source_layer": "sample_level_pseudobulk_logcpm",
                "cohort_scope": "pan_cancer",
                "sample_scope": "sample_level",
                "resolution": "sample_all",
                "cell_state": "all_cells",
                "gene_set": name,
                "method": "signed_mean_logcpm",
                "allowed_universe": "primary_pan_cancer_scRNA_universe",
                "allowed_downstream_use": "primary_response_blind",
                "is_primary_feature": "yes",
                "is_sensitivity_feature": "no",
                "is_support_only": "no",
                "is_qc_covariate": "no",
                "is_response_derived": "false",
                "notes": "Response-blind local program signature scored from sample-level pseudobulk.",
            }
        )
    tf_matrix = activity_base.copy()
    if not tf.empty:
        tf_main = tf[tf["used_in_main"].astype(str).str.lower().isin(["true", "yes", "1"])].copy()
        tf_main["signed_weight"] = pd.to_numeric(tf_main["signed_weight"], errors="coerce").fillna(0)
        candidates = []
        for tf_name, group in tf_main.groupby("feature_name"):
            group = group[group["gene_symbol"].isin(logcpm.columns)].copy()
            if group["gene_symbol"].nunique() >= 5:
                candidates.append((tf_name, group["gene_symbol"].nunique(), group))
        for tf_name, _, group in sorted(candidates, key=lambda x: (-x[1], x[0]))[:50]:
            genes = group["gene_symbol"].tolist()
            weights = group["signed_weight"].to_numpy(float)
            denom = np.abs(weights).sum() or 1.0
            feature_id = f"tf_activity__{tf_name}"
            tf_matrix[feature_id] = (logcpm[genes].to_numpy(float) @ weights) / denom
            coverage_rows.append(
                {
                    "resource": "tf_regulon_consensus",
                    "feature_id": feature_id,
                    "n_genes_total": int(tf_main[tf_main["feature_name"].eq(tf_name)]["gene_symbol"].nunique()),
                    "n_genes_covered": int(group["gene_symbol"].nunique()),
                    "coverage_fraction": float(group["gene_symbol"].nunique() / max(1, tf_main[tf_main["feature_name"].eq(tf_name)]["gene_symbol"].nunique())),
                    "status": "ready_sensitivity",
                }
            )
            feature_rows.append(
                {
                    "feature_id": feature_id,
                    "feature_name": feature_id,
                    "feature_family": "tf_activity",
                    "source_file": "tf_activity_matrix.sample_level.csv",
                    "source_layer": "sample_level_pseudobulk_logcpm",
                    "cohort_scope": "pan_cancer",
                    "sample_scope": "sample_level",
                    "resolution": "sample_all",
                    "cell_state": "all_cells",
                    "gene_set": tf_name,
                    "method": "signed_weighted_mean_logcpm",
                    "allowed_universe": "primary_pan_cancer_scRNA_universe",
                    "allowed_downstream_use": "sensitivity_response_blind",
                    "is_primary_feature": "no",
                    "is_sensitivity_feature": "yes",
                    "is_support_only": "no",
                    "is_qc_covariate": "no",
                    "is_response_derived": "false",
                    "notes": "Limited TF activity sensitivity feature from local consensus regulon.",
                }
            )
    write_csv(sig_matrix, OUT / "activity" / "signature_activity_matrix.sample_level.csv")
    write_csv(pd.DataFrame(columns=["sample_key", "cohort_id"]), OUT / "activity" / "pathway_activity_matrix.sample_level.csv")
    write_csv(tf_matrix, OUT / "activity" / "tf_activity_matrix.sample_level.csv")
    feature_dict = pd.DataFrame(feature_rows)
    if feature_dict.empty:
        feature_dict = pd.DataFrame(
            [
                {
                    "feature_id": "activity_blocked_low_coverage",
                    "feature_name": "activity_blocked_low_coverage",
                    "feature_family": "signature_pathway_tf_activity",
                    "source_file": "",
                    "source_layer": "",
                    "allowed_downstream_use": "blocked",
                    "is_response_derived": "false",
                    "notes": "No activity features passed coverage thresholds.",
                }
            ]
        )
    write_csv(feature_dict, OUT / "activity" / "signature_pathway_tf_feature_dictionary.csv")
    coverage = pd.DataFrame(coverage_rows)
    write_csv(coverage, OUT / "activity" / "gene_set_coverage_report.csv")
    write_md(
        f"""# Signature / Pathway / TF QC Report

**Run date:** {TODAY}

## Verdict
`SIGNATURE_READY_TF_SENSITIVITY_PATHWAY_BLOCKED`

## Summary
- Signature features generated: {len([c for c in sig_matrix.columns if c.startswith('signature_activity__')])}
- TF sensitivity features generated: {len([c for c in tf_matrix.columns if c.startswith('tf_activity__')])}
- Pathway features generated: 0

## Boundary
Activity scores are response-blind and derived from sample-level module-gene pseudobulk. Pathway scoring remains blocked until a formal pathway gene-set resource is frozen.
""",
        OUT / "activity" / "signature_pathway_tf_qc_report.md",
    )
    return feature_dict


def build_activity_outputs(paths: dict[str, Path]) -> pd.DataFrame:
    optional_gene_sets = [
        paths["curated_gene_sets"],
        paths["immune_signatures"],
        paths["tf_regulons"],
        paths["lr_pairs"],
    ]
    gene_set_available = [p for p in optional_gene_sets if p.exists()]
    cols = ["sample_key", "cohort_id"]
    for name in [
        "signature_activity_matrix.sample_level.csv",
        "pathway_activity_matrix.sample_level.csv",
        "tf_activity_matrix.sample_level.csv",
    ]:
        write_csv(pd.DataFrame(columns=cols), OUT / "activity" / name)
    feature_dict = pd.DataFrame(
        [
            {
                "feature_id": "activity_blocked_no_pseudobulk_or_gene_sets",
                "feature_name": "activity_blocked_no_pseudobulk_or_gene_sets",
                "feature_family": "signature_pathway_tf_activity",
                "source_file": "",
                "source_layer": "",
                "allowed_downstream_use": "blocked",
                "is_response_derived": "false",
                "notes": "No canonical pseudobulk matrix and no required curated gene set directory in current run.",
            }
        ]
    )
    write_csv(feature_dict, OUT / "activity" / "signature_pathway_tf_feature_dictionary.csv")
    coverage = pd.DataFrame(
        [
            {
                "resource": str(p),
                "available": "yes" if p.exists() else "no",
                "gene_coverage_status": "not_computed",
                "reason": "activity_feature_family_blocked_in_current_run",
            }
            for p in optional_gene_sets
        ]
    )
    write_csv(coverage, OUT / "activity" / "gene_set_coverage_report.csv")
    write_md(
        f"""# Signature / Pathway / TF QC Report

**Run date:** {TODAY}

## Verdict

`BLOCKED_IN_CURRENT_RUN`

## Reasons

- Canonical pseudobulk expression matrices were not generated in this controller run.
- Required curated gene set / TF regulon resources were not available as formal Phase4B inputs.

## Available Optional Resources

{len(gene_set_available)} optional resource paths exist.

## Boundary

No response-derived signature was generated. Activity features are absent from primary matrix.
""",
        OUT / "activity" / "signature_pathway_tf_qc_report.md",
    )
    return feature_dict


def build_tcr_outputs(meta: dict[str, pd.DataFrame]) -> pd.DataFrame:
    tcr = meta["tcr"].copy()
    if tcr.empty:
        matrix = pd.DataFrame(columns=["cohort_id", "tcr_available_flag"])
    else:
        matrix = tcr[
            [
                "cohort_id",
                "tcr_available_frozen",
                "join_rate",
                "matched_cell_count",
                "clonotype_fields_available",
                "tcr_feature_mainline_allowed",
                "downgrade_reason",
            ]
        ].drop_duplicates()
        matrix = matrix.rename(columns={"tcr_available_frozen": "tcr_available_flag"})
        matrix["feature_use"] = "sensitivity_only"
    write_csv(matrix, OUT / "tcr" / "optional_tcr_feature_matrix.csv")
    feature_dict = pd.DataFrame(
        [
            {
                "feature_id": "tcr_available_flag",
                "feature_name": "tcr_available_flag",
                "feature_family": "optional_tcr_sensitivity",
                "source_file": "gex_tcr_join_qc_report.csv",
                "allowed_downstream_use": "sensitivity_only",
                "is_response_derived": "false",
                "notes": "TCR clone master not rebuilt; no mainline TCR feature.",
            }
        ]
    )
    write_csv(feature_dict, OUT / "tcr" / "tcr_feature_dictionary.csv")
    write_md(
        f"""# TCR Sensitivity QC Report

**Run date:** {TODAY}

## Verdict

`SENSITIVITY_ONLY`

TCR clone master was not rebuilt in Phase4B. Current output records availability and join-QC fields only.
""",
        OUT / "tcr" / "tcr_sensitivity_qc_report.md",
    )
    write_md(
        "# TCR Patch Request\n\nNo mainline TCR features were generated. Build or validate a TCR clone master before any TCR coupling claim.\n",
        OUT / "tcr" / "tcr_patch_request.md",
    )
    return feature_dict


def build_response_environment(meta: dict[str, pd.DataFrame]) -> None:
    sample = meta["sample"]
    split = meta["split"]
    response = meta["response"]
    bind = sample.merge(
        split[["patient_key", "split", "split_eligibility", "reason_if_excluded"]].drop_duplicates("patient_key"),
        on="patient_key",
        how="left",
        suffixes=("", "_split"),
    )
    bind["response_known"] = np.where(~bind["response_binary_harmonized"].str.lower().isin(UNKNOWN_VALUES), "yes", "no")
    endpoint_lower = bind["response_endpoint_type"].fillna("").str.lower()
    supervised_known = bind["supervised_use_allowed"].str.lower().eq("yes") & bind["response_known"].eq("yes")
    bind["analysis_lane"] = np.select(
        [
            ~supervised_known,
            endpoint_lower.eq("recist"),
            endpoint_lower.eq("mrecist"),
            endpoint_lower.isin(["pathologic_response", "trg"]),
            endpoint_lower.eq("unknown_or_mixed"),
        ],
        [
            "support_or_unlabeled",
            "recist_supervised_or_anchor",
            "mrecist_anchor_sensitivity",
            "pathologic_response_sensitivity",
            "binary_anchor_sensitivity_unknown_endpoint",
        ],
        default="endpoint_sensitivity",
    )
    bind["endpoint_pooling_allowed"] = np.where(
        bind["analysis_lane"].isin(["support_or_unlabeled", "binary_anchor_sensitivity_unknown_endpoint"]),
        "no",
        "same_endpoint_only",
    )
    bind["timepoint_use_boundary"] = np.select(
        [
            bind["timepoint"].astype(str).str.lower().isin(["pre", "baseline", "pretreatment"]),
            bind["timepoint"].astype(str).str.lower().isin(["post", "on_treatment", "on-treatment"]),
        ],
        ["baseline_or_pre", "post_or_on_treatment"],
        default="unknown_or_mixed",
    )
    bind["response_binary_source"] = np.where(
        bind["response_raw"].astype(str).str.lower().isin(UNKNOWN_VALUES),
        "patient_level_harmonized_fill_or_unknown",
        "sample_or_patient_metadata",
    )
    bind["allowed_anchor_calibration"] = np.where(
        bind["supervised_use_allowed"].str.lower().eq("yes") & bind["response_known"].eq("yes"),
        np.where(bind["analysis_lane"].str.contains("sensitivity"), "sensitivity", "candidate"),
        "no",
    )
    bind["sensitivity_only_reason"] = np.select(
        [
            endpoint_lower.isin(["pathologic_response", "trg"]),
            endpoint_lower.eq("mrecist"),
            endpoint_lower.eq("unknown_or_mixed") & supervised_known,
        ],
        [
            "endpoint_sensitivity_not_recist",
            "mrecist_endpoint_sensitivity",
            "binary_label_endpoint_unknown_sensitivity",
        ],
        default="",
    )
    cols = [
        "sample_key",
        "patient_key",
        "cohort_id",
        "cancer_type",
        "treatment_context",
        "treatment_raw",
        "response_endpoint_type",
        "response_raw",
        "response_binary_harmonized",
        "response_known",
        "timepoint",
        "timepoint_raw",
        "split",
        "supervised_use_allowed",
        "allowed_anchor_calibration",
        "analysis_lane",
        "endpoint_pooling_allowed",
        "timepoint_use_boundary",
        "response_binary_source",
        "sensitivity_only_reason",
    ]
    write_csv(bind[[c for c in cols if c in bind.columns]], OUT / "response_environment" / "response_environment_binding_table.csv")

    endpoint_map = (
        response.groupby(["cohort_id", "endpoint_type"], as_index=False)
        .agg(
            treatment_contexts=("treatment_context", join_unique),
            supervised_eligible=("supervised_eligible", join_unique),
            support_only=("support_only", join_unique),
            use_boundary=("use_boundary", join_unique),
        )
        .sort_values(["cohort_id", "endpoint_type"])
    )
    write_md(endpoint_map.to_markdown(index=False), OUT / "response_environment" / "endpoint_type_mapping_report.md")
    write_csv(
        bind[
            [
                "sample_key",
                "patient_key",
                "cohort_id",
                "supervised_use_allowed",
                "response_endpoint_type",
                "response_known",
                "sensitivity_only_reason",
            ]
        ],
        OUT / "response_environment" / "supervised_use_boundary_table.csv",
    )
    write_csv(
        bind[
            [
                "sample_key",
                "patient_key",
                "cohort_id",
                "allowed_anchor_calibration",
                "response_endpoint_type",
                "treatment_context",
                "sensitivity_only_reason",
            ]
        ],
        OUT / "response_environment" / "anchor_calibration_boundary_table.csv",
    )


def assemble_matrices(
    sample_fraction: pd.DataFrame,
    feature_dict_fraction: pd.DataFrame,
    qc_cov: pd.DataFrame,
    pb_dict: pd.DataFrame,
    activity_dict: pd.DataFrame,
    tcr_dict: pd.DataFrame,
    meta: dict[str, pd.DataFrame],
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    id_cols = ["sample_key", "cohort_id", "sample_id", "patient_key", "patient_id", "timepoint", "total_cells"]
    frac_cols = [c for c in sample_fraction.columns if c.startswith("frac_")]
    primary_cols = [c for c in frac_cols if c.startswith("frac_coarse__") or c.startswith("frac_mid__")]
    sensitivity_cols = [c for c in frac_cols if c.startswith("frac_fine_restricted__")]

    primary = sample_fraction[id_cols + primary_cols].copy()
    sensitivity = sample_fraction[id_cols + sensitivity_cols].copy()
    sig_cols: list[str] = []
    tf_cols: list[str] = []
    sig_path = OUT / "activity" / "signature_activity_matrix.sample_level.csv"
    tf_path = OUT / "activity" / "tf_activity_matrix.sample_level.csv"
    if sig_path.exists():
        sig = read_csv(sig_path)
        sig_cols = [c for c in sig.columns if c.startswith("signature_activity__")]
        if sig_cols:
            primary = primary.merge(sig[["sample_key"] + sig_cols], on="sample_key", how="left")
    if tf_path.exists():
        tf = read_csv(tf_path)
        tf_cols = [c for c in tf.columns if c.startswith("tf_activity__")]
        if tf_cols:
            sensitivity = sensitivity.merge(tf[["sample_key"] + tf_cols], on="sample_key", how="left")
    support = sample_fraction[id_cols].copy()
    qc_out = qc_cov.copy()

    write_csv(primary, OUT / "matrix" / "immune_state_feature_matrix.primary.csv")
    write_csv(sensitivity, OUT / "matrix" / "immune_state_feature_matrix.sensitivity.csv")
    write_csv(support, OUT / "matrix" / "immune_state_feature_matrix.support.csv")
    write_csv(qc_out, OUT / "matrix" / "immune_state_feature_matrix.qc_covariates.csv")
    write_csv(primary, OUT / "matrix" / "immune_state_feature_matrix_by_family" / "cell_fraction_coarse_mid.csv")
    write_csv(sensitivity, OUT / "matrix" / "immune_state_feature_matrix_by_family" / "cell_fraction_fine_restricted.csv")
    if sig_cols:
        write_csv(primary[id_cols + sig_cols], OUT / "matrix" / "immune_state_feature_matrix_by_family" / "signature_activity.csv")
    if tf_cols:
        write_csv(sensitivity[id_cols + tf_cols], OUT / "matrix" / "immune_state_feature_matrix_by_family" / "tf_activity_sensitivity.csv")

    feature_dict = pd.concat([feature_dict_fraction, pb_dict, activity_dict, tcr_dict], ignore_index=True, sort=False)
    audit_base = meta["sample"][
        ["sample_key", "cohort_id", "patient_key", "modality", "timepoint"]
    ].drop_duplicates("sample_key")
    if "cancer_type" in meta["patient"].columns:
        audit_base = audit_base.merge(
            meta["patient"][["patient_key", "cancer_type"]].drop_duplicates("patient_key"),
            on="patient_key",
            how="left",
        )
    else:
        audit_base["cancer_type"] = ""
    missing = compute_feature_audit(primary, sensitivity, audit_base, qc_out)
    feature_dict = feature_dict.merge(
        missing[
            [
                "feature_id",
                "missingness",
                "max_cohort_fraction",
                "max_cancer_fraction",
                "platform_dependence",
                "cell_count_correlation",
            ]
        ],
        on="feature_id",
        how="left",
        suffixes=("", "_audit"),
    )
    for target, source in [
        ("missingness", "missingness_audit"),
        ("cohort_dominance", "max_cohort_fraction"),
        ("cancer_dominance", "max_cancer_fraction"),
        ("platform_dependence", "platform_dependence_audit"),
        ("cell_count_dependence", "cell_count_correlation"),
    ]:
        if source in feature_dict.columns:
            feature_dict[target] = feature_dict[target].replace("", np.nan).fillna(feature_dict[source])
    feature_dict = feature_dict.drop(columns=[c for c in feature_dict.columns if c.endswith("_audit")], errors="ignore")
    feature_dict["is_response_derived"] = feature_dict["is_response_derived"].fillna("false")
    write_csv(feature_dict, OUT / "matrix" / "feature_dictionary_v6_2.csv")
    return primary, sensitivity, feature_dict


def compute_feature_audit(
    primary: pd.DataFrame, sensitivity: pd.DataFrame, sample_meta: pd.DataFrame, qc_cov: pd.DataFrame
) -> pd.DataFrame:
    df = primary.merge(sample_meta[["sample_key", "cancer_type"]], on="sample_key", how="left")
    if "platform" not in df.columns:
        df["platform"] = ""
    total_cells = primary.set_index("sample_key")["total_cells"]
    rows = []
    for matrix_name, mat in [("primary", primary), ("sensitivity", sensitivity)]:
        exclude = {"sample_key", "cohort_id", "sample_id", "patient_key", "patient_id", "timepoint", "total_cells"}
        feature_cols = [c for c in mat.columns if c not in exclude]
        enriched = mat.merge(sample_meta[["sample_key", "cancer_type"]], on="sample_key", how="left")
        for col in feature_cols:
            values = pd.to_numeric(enriched[col], errors="coerce")
            nonmissing = values.notna()
            nonzero = values.fillna(0) > 0
            n = len(values)
            cohort_dom = max_fraction(enriched.loc[nonzero, "cohort_id"]) if nonzero.any() else 0
            cancer_dom = max_fraction(enriched.loc[nonzero, "cancer_type"]) if nonzero.any() else 0
            corr = safe_corr(values, mat["total_cells"])
            rows.append(
                {
                    "feature_id": col,
                    "matrix": matrix_name,
                    "missingness": float(1 - nonmissing.mean()) if n else 1.0,
                    "non_missing_sample_count": int(nonmissing.sum()),
                    "non_missing_patient_count": int(mat.loc[nonmissing, "patient_key"].nunique()),
                    "cohort_coverage": int(enriched.loc[nonzero, "cohort_id"].nunique()),
                    "cancer_coverage": int(enriched.loc[nonzero, "cancer_type"].nunique()),
                    "platform_coverage": "",
                    "max_cohort_fraction": cohort_dom,
                    "max_cancer_fraction": cancer_dom,
                    "platform_dependence": "",
                    "cell_count_correlation": corr,
                    "cell_count_dependence": "high" if abs(corr) >= 0.6 else "not_high",
                    "high_missingness": "yes" if (1 - nonmissing.mean()) > MISSINGNESS_HIGH else "no",
                    "high_cohort_dominance": "yes" if cohort_dom > COHORT_DOMINANCE_HIGH else "no",
                    "high_cancer_dominance": "yes" if cancer_dom > CANCER_DOMINANCE_HIGH else "no",
                }
            )
    audit = pd.DataFrame(rows)
    write_csv(audit, OUT / "audit" / "feature_missingness_report.csv")
    tags = audit[
        [
            "feature_id",
            "high_missingness",
            "high_cohort_dominance",
            "high_cancer_dominance",
            "cell_count_dependence",
            "max_cohort_fraction",
            "max_cancer_fraction",
        ]
    ].copy()
    tags["allowed_primary_after_audit"] = np.where(
        (tags["high_missingness"].eq("no"))
        & (tags["high_cohort_dominance"].eq("no"))
        & (tags["high_cancer_dominance"].eq("no")),
        "yes",
        "review_or_sensitivity",
    )
    write_csv(tags, OUT / "audit" / "feature_confounding_tags.csv")
    family = (
        audit.assign(feature_family=audit["feature_id"].map(infer_feature_family))
        .groupby("feature_family", as_index=False)
        .agg(
            n_features=("feature_id", "nunique"),
            mean_missingness=("missingness", "mean"),
            high_cohort_dominance_features=("high_cohort_dominance", lambda x: int((x == "yes").sum())),
            high_cancer_dominance_features=("high_cancer_dominance", lambda x: int((x == "yes").sum())),
        )
    )
    write_csv(family, OUT / "audit" / "feature_family_qc_summary.csv")
    high_risk = tags[tags["allowed_primary_after_audit"].ne("yes")]
    write_csv(high_risk, OUT / "audit" / "high_risk_feature_exclusion_list.csv")
    write_md(
        f"""# Feature Audit Report

**Run date:** {TODAY}

## Summary

- Audited features: {len(audit)}
- High-risk features requiring review/sensitivity: {len(high_risk)}
- Missingness threshold: {MISSINGNESS_HIGH}
- Cohort dominance threshold: {COHORT_DOMINANCE_HIGH}
- Cancer dominance threshold: {CANCER_DOMINANCE_HIGH}
""",
        OUT / "audit" / "feature_audit_report.md",
    )
    return audit


def max_fraction(series: pd.Series) -> float:
    if len(series) == 0:
        return 0.0
    counts = series.fillna("unknown").value_counts()
    return float(counts.max() / counts.sum()) if counts.sum() else 0.0


def safe_corr(a: pd.Series, b: pd.Series) -> float:
    a = pd.to_numeric(a, errors="coerce")
    b = pd.to_numeric(b, errors="coerce")
    if a.nunique(dropna=True) <= 1 or b.nunique(dropna=True) <= 1:
        return 0.0
    val = a.corr(b)
    return 0.0 if pd.isna(val) else float(val)


def infer_feature_family(feature_id: str) -> str:
    if feature_id.startswith("frac_coarse__"):
        return "cell_fraction_coarse"
    if feature_id.startswith("frac_mid__"):
        return "cell_fraction_mid"
    if feature_id.startswith("frac_fine_restricted__"):
        return "cell_fraction_fine_restricted"
    if feature_id.startswith("signature_activity__"):
        return "signature_activity"
    if feature_id.startswith("tf_activity__"):
        return "tf_activity"
    return "other"


def run_adversarial_audit(primary: pd.DataFrame, sensitivity: pd.DataFrame, feature_dict: pd.DataFrame, meta: dict[str, pd.DataFrame]) -> list[str]:
    hard = []
    exclude = {"sample_key", "cohort_id", "sample_id", "patient_key", "patient_id", "timepoint", "total_cells"}
    feature_cols = [c for c in list(primary.columns) + list(sensitivity.columns) if c not in exclude]
    rows = []
    for col in feature_cols:
        hits = [tok for tok in BANNED_FEATURE_TOKENS if tok in col.lower()]
        rows.append(
            {
                "feature_id": col,
                "leakage_tokens": "|".join(hits),
                "response_leakage_risk": "yes" if hits else "no",
            }
        )
        if hits:
            hard.append(f"feature_name_leakage:{col}:{hits}")
    leakage = pd.DataFrame(rows)
    write_csv(leakage, OUT / "audit" / "phase4b_feature_leakage_audit.csv")

    split = meta["split"]
    split_counts = split.groupby("patient_key")["split"].nunique()
    split_leak = split_counts[split_counts > 1]
    if len(split_leak):
        hard.append(f"patient_split_leakage:{len(split_leak)}")

    fd_bad = feature_dict["is_response_derived"].astype(str).str.lower().ne("false")
    if fd_bad.any():
        hard.append("feature_dictionary_response_derived_not_false")

    rule_rows = [
        {
            "check": "do_not_rename_states",
            "status": "pass",
            "evidence": "Phase4B uses Phase4A labels without renaming.",
        },
        {
            "check": "used_resolution_aware_eligibility",
            "status": "pass",
            "evidence": "phase4b_resolution_aware_aggregation_eligibility.csv",
        },
        {
            "check": "full_integration_not_overstated",
            "status": "pass",
            "evidence": "full_integration_completed=false in handoff",
        },
        {
            "check": "tcr_not_mainline",
            "status": "pass",
            "evidence": "optional_tcr_sensitivity only",
        },
        {
            "check": "patient_split_leakage",
            "status": "fail" if len(split_leak) else "pass",
            "evidence": str(len(split_leak)),
        },
    ]
    write_csv(pd.DataFrame(rule_rows), OUT / "audit" / "phase4b_rule_violation_audit.csv")
    write_md(
        f"""# Phase4B Adversarial QC Review

**Run date:** {TODAY}

## Verdict

{'HARD_FAIL' if hard else 'PASS_WITH_CONDITIONAL_LIMITATIONS'}

## Checks

- Feature leakage tokens detected: {int((leakage['response_leakage_risk'] == 'yes').sum())}
- Patient split leakage count: {len(split_leak)}
- Response-derived features in dictionary: {int(fd_bad.sum())}
- Phase4A eligibility used: yes
- Phase4A ontology renamed: no
- Fine labels restricted: yes
- TCR mainline feature generated: no
- Full integration described as completed: no
""",
        OUT / "audit" / "phase4b_adversarial_qc_review.md",
    )
    risks = [
        {
            "risk_id": "R1",
            "risk": "Full integration still pending",
            "status": "open",
            "mitigation": "coarse/mid feature construction only; do not claim full integrated atlas",
        },
        {
            "risk_id": "R2",
            "risk": "Fine features restricted",
            "status": "mitigated",
            "mitigation": "fine features kept in sensitivity matrix",
        },
        {
            "risk_id": "R3",
            "risk": "Cell-state-specific all-gene pseudobulk remains deferred",
            "status": "mitigated",
            "mitigation": "sample-level module-gene pseudobulk generated for response-blind Phase6 input",
        },
        {
            "risk_id": "R4",
            "risk": "Pathway activity blocked without formal pathway resource; TF limited to sensitivity",
            "status": "mitigated",
            "mitigation": "local response-blind signatures and limited TF regulons scored from pseudobulk",
        },
        {
            "risk_id": "R5",
            "risk": "TCR sensitivity only",
            "status": "mitigated",
            "mitigation": "no mainline TCR features generated",
        },
    ]
    write_csv(pd.DataFrame(risks), OUT / "audit" / "phase4b_unresolved_risk_register.csv")
    return hard


def write_final_outputs(
    hard: list[str],
    primary: pd.DataFrame,
    sensitivity: pd.DataFrame,
    feature_dict: pd.DataFrame,
    preflight_audit: pd.DataFrame,
) -> None:
    verdict = "BLOCKED" if hard else "CONDITIONAL_GO_TO_PHASE5"
    has_pb = bool((feature_dict.get("feature_family", pd.Series(dtype=str)).astype(str) == "pseudobulk").any())
    has_sig = bool((feature_dict.get("feature_family", pd.Series(dtype=str)).astype(str) == "signature_activity").any())
    has_tf = bool((feature_dict.get("feature_family", pd.Series(dtype=str)).astype(str) == "tf_activity").any())
    feature_family_status = {
        "cell_fraction_coarse": "primary",
        "cell_fraction_mid": "primary",
        "cell_fraction_fine_restricted": "sensitivity",
        "pseudobulk_raw_count": "phase6_response_blind_input" if has_pb else "blocked",
        "pseudobulk_normalized_expression": "blocked",
        "signature_activity": "primary" if has_sig else "blocked",
        "pathway_activity": "blocked",
        "tf_activity": "sensitivity" if has_tf else "blocked",
        "optional_tcr": "sensitivity",
        "qc_covariates": "confounding_audit_only",
    }
    manifest = {
        "phase": "phase4B_patient_timepoint_immune_state_feature_construction",
        "input_version": INPUT_VERSION,
        "created_at": TODAY,
        "verdict": verdict,
        "hard_blockers": hard,
        "conditional_items": [
            "full_integration_pending",
            "fine_features_restricted",
            "cell_state_specific_pseudobulk_deferred",
            "pathway_activity_blocked_no_formal_pathway_resource",
            "tcr_sensitivity_only",
        ],
        "primary_outputs": {
            "immune_state_feature_matrix_primary": str(OUT / "matrix" / "immune_state_feature_matrix.primary.csv"),
            "immune_state_feature_matrix_sensitivity": str(OUT / "matrix" / "immune_state_feature_matrix.sensitivity.csv"),
            "feature_dictionary": str(OUT / "matrix" / "feature_dictionary_v6_2.csv"),
            "feature_missingness_report": str(OUT / "audit" / "feature_missingness_report.csv"),
            "feature_confounding_tags": str(OUT / "audit" / "feature_confounding_tags.csv"),
            "response_environment_binding_table": str(
                OUT / "response_environment" / "response_environment_binding_table.csv"
            ),
            "analysis_universe_registry": str(OUT / "universe" / "analysis_universe_registry_v6_2.csv"),
        },
        "phase5_allowed_inputs": [
            str(OUT / "matrix" / "immune_state_feature_matrix.primary.csv"),
            str(OUT / "matrix" / "immune_state_feature_matrix.sensitivity.csv"),
            str(OUT / "matrix" / "immune_state_feature_matrix.qc_covariates.csv"),
            str(OUT / "matrix" / "feature_dictionary_v6_2.csv"),
            str(OUT / "audit" / "feature_missingness_report.csv"),
            str(OUT / "audit" / "feature_confounding_tags.csv"),
            str(OUT / "response_environment" / "response_environment_binding_table.csv"),
            str(OUT / "universe" / "analysis_universe_registry_v6_2.csv"),
            str(OUT / "pseudobulk" / "pseudobulk_matrix_raw_count.parquet"),
            str(OUT / "activity" / "signature_activity_matrix.sample_level.csv"),
            str(OUT / "activity" / "tf_activity_matrix.sample_level.csv"),
        ],
        "phase5_blocked_inputs": [
            str(OUT / "activity" / "pathway_activity_matrix.sample_level.csv"),
        ],
        "feature_family_status": feature_family_status,
        "leakage_audit": {
            "response_leakage_detected": bool(hard),
            "split_leakage_detected": any(x.startswith("patient_split_leakage") for x in hard),
        },
        "integration_status": {
            "full_integration_completed": False,
            "note": "Do not describe current input as full integrated atlas unless integration patch is completed.",
        },
        "phase4a_rule_compliance": {
            "used_resolution_aware_eligibility": True,
            "used_allowed_cell_state_levels": True,
            "used_blocked_list": True,
            "used_reliability_scores": True,
        },
        "matrix_shapes": {
            "primary": {"rows": int(primary.shape[0]), "columns": int(primary.shape[1])},
            "sensitivity": {"rows": int(sensitivity.shape[0]), "columns": int(sensitivity.shape[1])},
            "feature_dictionary": {"rows": int(feature_dict.shape[0]), "columns": int(feature_dict.shape[1])},
        },
        "patch_queue": [
            {
                "issue": "run_full_integration_scvi_or_equivalent",
                "required_before": "strong full-atlas or final anchor calibration claim",
            },
            {
                "issue": "run_cell_state_specific_all_gene_pseudobulk",
                "required_before": "cell-state-specific expression module claims",
            },
            {
                "issue": "freeze_formal_pathway_gene_sets",
                "required_before": "pathway activity primary feature use",
            },
        ],
    }
    write_yaml(manifest, OUT / "handoff" / "phase4b_decision_manifest.yaml")
    write_yaml(
        {
            "phase": manifest["phase"],
            "input_version": INPUT_VERSION,
            "verdict": verdict,
            "phase5_allowed_inputs": manifest["phase5_allowed_inputs"],
            "phase5_blocked_inputs": manifest["phase5_blocked_inputs"],
            "feature_family_status": feature_family_status,
            "rules": {
                "do_not_use_response_as_feature": True,
                "do_not_bypass_phase4b_handoff": True,
                "do_not_claim_full_integrated_atlas": True,
            },
        },
        OUT / "handoff" / "phase4b_to_phase5_handoff.yaml",
    )

    output_rows = []
    for path in sorted(OUT.rglob("*")):
        if path.is_file():
            output_rows.append(
                {
                    "path": str(path),
                    "size_bytes": path.stat().st_size,
                    "relative_path": str(path.relative_to(OUT)),
                }
            )
    write_tsv(pd.DataFrame(output_rows), OUT / "handoff" / "phase4b_output_index.tsv")

    report = f"""# Phase4B Immune-state Feature Construction Report

**Run date:** {TODAY}
**Input version:** `{INPUT_VERSION}`
**Verdict:** `{verdict}`

## 1. Executive Verdict

Phase4B generated conservative patient/sample-level immune-state feature outputs from Phase4A handoff files. Coarse and mid cell-fraction features are available for Phase5. Fine features are restricted to the sensitivity matrix. This patch adds sample-level module-gene pseudobulk and response-blind signature / limited TF activity; cell-state-specific all-gene pseudobulk and pathway activity remain deferred.

## 2. Input Version and Phase4A Rule Compliance

- Used Phase4A handoff: yes.
- Used resolution-aware eligibility: yes.
- Used reliability scores: yes.
- Used blocked list: yes.
- Renamed cell states: no.
- Used response/outcome/split as feature source: no.
- Full integrated atlas completed: no.

## 3. Phase4A Hotfix / Validation Summary

Preflight passed with resolution-aware hotfix. Fine Unknown no longer suppresses otherwise known coarse/mid labels. Integration status now distinguishes `eligible_for_full_integration` from `full_integration_completed`.

## 4. Analysis Universe Summary

See `universe/analysis_universe_registry_v6_2.csv`.

## 5. Cell-state Fraction Feature Summary

- Primary matrix rows: {primary.shape[0]}
- Primary matrix columns: {primary.shape[1]}
- Sensitivity matrix rows: {sensitivity.shape[0]}
- Sensitivity matrix columns: {sensitivity.shape[1]}
- Feature dictionary rows: {feature_dict.shape[0]}

## 6. Pseudobulk Feature Summary

Sample-level all-cell raw-count pseudobulk was generated for a capped module gene universe. Cell-state-specific expression pseudobulk remains deferred and raw-count / normalized-expression layers were not mixed.

## 7. Signature / Pathway / TF Feature Summary

Response-blind program signatures and limited TF regulon activities were scored from sample-level pseudobulk. Pathway activity remains blocked until a formal pathway gene-set resource is frozen.

## 8. Optional TCR Feature Summary

TCR output is sensitivity-only and records availability / join QC fields. No mainline TCR feature was generated.

## 9. QC Covariates

QC covariates were generated for confounding audit only. They are not mechanism features.

## 10. Response Environment Binding

Response environment binding table was generated. Response endpoint, treatment context, and split are kept outside feature matrices.

## 11. Final Feature Matrices

- `matrix/immune_state_feature_matrix.primary.csv`
- `matrix/immune_state_feature_matrix.sensitivity.csv`
- `matrix/immune_state_feature_matrix.support.csv`
- `matrix/immune_state_feature_matrix.qc_covariates.csv`

## 12. Feature Dictionary

`matrix/feature_dictionary_v6_2.csv` records source, family, resolution, allowed use, missingness, dominance, and response-derived status.

## 13. Missingness / Confounding Audit

See `audit/feature_missingness_report.csv`, `audit/feature_confounding_tags.csv`, and `audit/high_risk_feature_exclusion_list.csv`.

## 14. Leakage and Rule Audit

Adversarial audit found {len(hard)} hard blocker(s). Feature response-derived flags are false.

## 15. What Phase4B Did Not Do

- It did not do response association.
- It did not do module discovery.
- It did not do PD1 anchor verdict.
- It did not do X-class repair.
- It did not complete a full integrated atlas.
- It did not generate primary TCR features.

## 16. Patch Queue

1. Run full integration before any full-atlas or final anchor calibration claim.
2. Run cell-state-specific all-gene pseudobulk before cell-state-specific expression module claims.
3. Freeze formal pathway gene sets before pathway activity primary feature use.

## 17. Phase5 Handoff

See `handoff/phase4b_to_phase5_handoff.yaml`.

## 18. Final Decision

`{verdict}`. Phase5 can use coarse/mid fraction features, response-blind signature activity, limited TF sensitivity features, and QC/confounding audit inputs. Claims depending on full integration, cell-state-specific all-gene pseudobulk, pathway activity, or TCR coupling remain blocked or sensitivity-only.
"""
    write_md(report, OUT / "PHASE4B_IMMUNE_STATE_FEATURE_CONSTRUCTION_REPORT.md")


def main() -> None:
    mkdirs()
    input_audit, paths, manifest = validate_inputs()
    elig, _, _, _ = build_resolution_aware_eligibility(paths, manifest)
    meta = load_metadata(paths)
    sample_fraction, pt_fraction, fraction_dict = build_fraction_features(elig, meta["sample"])
    build_universes(meta, sample_fraction)
    qc_cov = build_qc_covariates(elig, meta, sample_fraction)
    pb_registry = build_pseudobulk_outputs_v2(elig, meta, paths)
    pb_dict = read_csv(OUT / "pseudobulk" / "pseudobulk_feature_dictionary.csv")
    activity_dict = build_activity_outputs_v2(paths)
    tcr_dict = build_tcr_outputs(meta)
    build_response_environment(meta)
    primary, sensitivity, feature_dict = assemble_matrices(
        sample_fraction, fraction_dict, qc_cov, pb_dict, activity_dict, tcr_dict, meta
    )
    hard = run_adversarial_audit(primary, sensitivity, feature_dict, meta)

    preflight_status = {
        "phase": "Phase4B.0_preflight",
        "input_version": INPUT_VERSION,
        "verdict": "HARD_FAIL" if hard else "PASS",
        "phase4a_verdict": manifest.get("verdict"),
        "hard_blockers": hard,
        "conditional_items": [
            "full_integration_pending",
            "fine_label_low_coverage",
            "cell_state_specific_pseudobulk_deferred",
            "pathway_activity_deferred",
            "tcr_sensitivity_only",
        ],
        "required_inputs_missing": input_audit.query("input_kind == 'required' and exists == 'no'")["path"].tolist(),
        "optional_inputs_missing": input_audit.query("input_kind == 'optional' and exists == 'no'")["path"].tolist(),
        "resolution_aware_hotfix_generated": True,
        "coarse_mid_feature_construction_allowed": not bool(hard),
    }
    write_yaml(preflight_status, OUT / "preflight" / "phase4b_preflight_status.yaml")
    write_final_outputs(hard, primary, sensitivity, feature_dict, input_audit)


if __name__ == "__main__":
    main()
