#!/usr/bin/env python3
from __future__ import annotations

import argparse
import hashlib
import json
import re
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import yaml

RUN_ID = "step2_v6_1_0505_0319"
STEP2_ROOT = Path("/home/huyudi/006/results/v6_1/step2") / RUN_ID
CONS_DIR = STEP2_ROOT / "04_annotation" / "consensus_annotation_v1"
MYELOID_DIR = STEP2_ROOT / "05_myeloid_qc" / "myeloid_adjudication_v1"
OUT_DIR = STEP2_ROOT / "06_fraction" / "fraction_feature_v1"

CONSENSUS_PATH = CONS_DIR / "cell_state_annotation_consensus_v1.parquet"
CONSENSUS_MANIFEST_PATH = CONS_DIR / "consensus_annotation_manifest.yaml"
MYELOID_PATH = MYELOID_DIR / "myeloid_annotation_adjudication.parquet"
MYELOID_MANIFEST_PATH = MYELOID_DIR / "myeloid_decision_manifest.yaml"
FLAGS_PATH = STEP2_ROOT / "03_qc" / "final_cell_inclusion_flags.qc0513_balanced.parquet"
INPUT_MANIFEST_REF = MYELOID_MANIFEST_PATH

KEY_COLS = ["source_h5ad", "cell_barcode", "sample_id"]
CELL_KEY_COLS = ["source_h5ad", "cohort_id", "cell_barcode", "sample_id"]
SAMPLE_KEY = ["cohort_id", "sample_id"]
PROVENANCE_COLS = ["run_id", "created_at", "input_manifest_ref"]
FORBIDDEN_RESPONSE_TERMS = ["response_strict", "response_broad", "response_binary", "response_ordered", "response_raw", "response_"]

IMMUNE_MAJORS = {"T_NK", "B_Plasma", "Myeloid_DC", "Mast_or_minor_immune", "Immune"}
NON_MYELOID_IMMUNE = {"T_NK", "B_Plasma", "Mast_or_minor_immune", "Immune"}
MYELOID_MAIN_LABELS = {
    "Mono_FCN1",
    "Mono_inflammatory",
    "Macro_C1QC",
    "Macro_SPP1_TAM_like",
    "cDC1",
    "cDC2",
    "pDC",
    "Mast_or_minor_immune",
    "Myeloid_DC_core",
    "Myeloid_unspecified",
}
MISSING_MID = {"", "NA", "nan", "None", "Unknown"}

THRESHOLDS = {
    "low_total_cells": 200,
    "low_immune_cells": 30,
    "low_confidence_denominator": 30,
    "missing_denominator": 10,
    "high_unknown_fraction": 0.20,
    "high_immune_unspecified_fraction": 0.30,
    "high_conflict_fraction": 0.05,
    "myeloid_uncertain_fraction": 0.20,
    "myeloid_sensitivity_fraction": 0.20,
}

PARENT_CHILDREN = {
    "T_NK": ["CD8_T", "CD4_T", "Treg", "NK", "NKT_like", "T_NK_core"],
    "CD4_T": ["Treg", "CD4_T"],
    "B_Plasma": ["B_naive_memory", "Plasma", "B_Plasma_core"],
    "Myeloid_DC": [
        "Mono_FCN1",
        "Mono_inflammatory",
        "Macro_C1QC",
        "Macro_SPP1_TAM_like",
        "cDC1",
        "cDC2",
        "pDC",
        "Myeloid_DC_core",
        "Myeloid_unspecified",
    ],
    "Mast_or_minor_immune": ["Mast_or_minor_immune", "Mast_minor_immune"],
    "Immune": ["Immune_unspecified"],
}


def now_iso() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat()


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def sanitize_label(label: Any) -> str:
    text = "Unknown" if pd.isna(label) else str(label).strip()
    if text in MISSING_MID:
        text = "Unknown"
    text = re.sub(r"[^0-9A-Za-z]+", "_", text).strip("_")
    return text or "Unknown"


def feature_name(family: str, label: str, parent: str | None = None) -> str:
    clean_label = sanitize_label(label)
    if family == "all":
        return f"frac_all__{clean_label}"
    if family == "immune":
        return f"frac_immune__{clean_label}"
    if family == "parent":
        if parent is None:
            raise ValueError("parent feature requires parent label")
        return f"frac_parent_{sanitize_label(parent)}__{clean_label}"
    raise ValueError(f"unknown feature family: {family}")


def fraction_value(numerator: int | float, denominator: int | float) -> float:
    if pd.isna(denominator) or float(denominator) < THRESHOLDS["missing_denominator"]:
        return np.nan
    return float(numerator) / float(denominator)


def add_provenance(df: pd.DataFrame, created_at: str) -> pd.DataFrame:
    out = df.copy()
    out.insert(0, "input_manifest_ref", str(INPUT_MANIFEST_REF))
    out.insert(0, "created_at", created_at)
    out.insert(0, "run_id", RUN_ID)
    return out


def load_step2_6_low_coverage_sources(path: Path) -> set[str]:
    with path.open("r", encoding="utf-8") as f:
        manifest = yaml.safe_load(f)
    validation = manifest.get("validation", {}) if isinstance(manifest, dict) else {}
    out = set()
    for row in validation.get("low_coverage_sources", []) or []:
        if isinstance(row, dict) and row.get("source"):
            out.add(str(row["source"]))
    return out


def validate_consensus_against_flags(consensus: pd.DataFrame, flags_path: Path) -> dict[str, Any]:
    flags = pd.read_parquet(flags_path, columns=CELL_KEY_COLS + ["final_inclusion"])
    flags_main = flags[flags["final_inclusion"].astype(str).eq("include_main")][CELL_KEY_COLS].copy()
    consensus_main = consensus[consensus["included_in_main_annotation"].astype(bool)][CELL_KEY_COLS].copy()

    if len(flags_main) != len(consensus_main):
        raise RuntimeError(
            f"include_main row count mismatch between flags ({len(flags_main)}) "
            f"and consensus ({len(consensus_main)})"
        )
    if flags_main.duplicated(CELL_KEY_COLS).any():
        raise RuntimeError("final inclusion flags include_main keys are not unique")
    if consensus_main.duplicated(CELL_KEY_COLS).any():
        raise RuntimeError("consensus include_main keys are not unique")

    left = consensus_main.merge(flags_main, on=CELL_KEY_COLS, how="left", indicator=True)
    missing_in_flags = int((left["_merge"] != "both").sum())
    right = flags_main.merge(consensus_main, on=CELL_KEY_COLS, how="left", indicator=True)
    missing_in_consensus = int((right["_merge"] != "both").sum())
    if missing_in_flags or missing_in_consensus:
        raise RuntimeError(
            "include_main key mismatch between consensus and flags: "
            f"missing_in_flags={missing_in_flags}, missing_in_consensus={missing_in_consensus}"
        )
    return {
        "flags_path": str(flags_path),
        "flags_rows": int(len(flags)),
        "flags_include_main_rows": int(len(flags_main)),
        "consensus_include_main_rows": int(len(consensus_main)),
        "include_main_keys_match": True,
    }


def load_inputs() -> tuple[pd.DataFrame, pd.DataFrame, set[str]]:
    consensus_cols = [
        "source_h5ad",
        "cohort_id",
        "cell_barcode",
        "sample_id",
        "included_in_main_annotation",
        "final_major_lineage",
        "final_immune_mid_state",
        "final_confidence",
        "decision_rule_id",
        "conflict_type",
    ]
    myeloid_cols = [
        "source_h5ad",
        "cohort_id",
        "cell_barcode",
        "sample_id",
        "myeloid_adjudication_status",
        "final_myeloid_label",
        "included_in_myeloid_main",
        "myeloid_sensitivity_only",
    ]
    consensus_all = pd.read_parquet(CONSENSUS_PATH, columns=consensus_cols)
    flags_check = validate_consensus_against_flags(consensus_all, FLAGS_PATH)
    consensus = consensus_all[consensus_all["included_in_main_annotation"].astype(bool)].copy()
    consensus.attrs["flags_check"] = flags_check
    myeloid = pd.read_parquet(MYELOID_PATH, columns=myeloid_cols)
    low_cov_sources = load_step2_6_low_coverage_sources(MYELOID_MANIFEST_PATH)
    return consensus, myeloid, low_cov_sources


def apply_myeloid_overlay(consensus: pd.DataFrame, myeloid: pd.DataFrame) -> pd.DataFrame:
    my_cols = KEY_COLS + [
        "myeloid_adjudication_status",
        "final_myeloid_label",
        "included_in_myeloid_main",
        "myeloid_sensitivity_only",
    ]
    df = consensus.merge(myeloid[my_cols], on=KEY_COLS, how="left", validate="one_to_one")
    df["final_major_lineage"] = df["final_major_lineage"].fillna("Unknown").astype(str)
    df["final_immune_mid_state"] = df["final_immune_mid_state"].fillna("NA").astype(str)
    df["conflict_type"] = df["conflict_type"].fillna("").astype(str)
    df["final_confidence"] = df["final_confidence"].fillna("low").astype(str)
    df["myeloid_adjudication_status"] = df["myeloid_adjudication_status"].fillna("not_myeloid_candidate").astype(str)
    df["final_myeloid_label"] = df["final_myeloid_label"].fillna("NA").astype(str)
    df["included_in_myeloid_main"] = df["included_in_myeloid_main"].fillna(False).astype(bool)
    df["myeloid_sensitivity_only"] = df["myeloid_sensitivity_only"].fillna(False).astype(bool)

    state = df["final_major_lineage"].copy()
    parent = df["final_major_lineage"].copy()

    immune_mask = df["final_major_lineage"].isin(IMMUNE_MAJORS)
    mid = df["final_immune_mid_state"].where(~df["final_immune_mid_state"].isin(MISSING_MID), "Immune_unspecified")
    state.loc[immune_mask] = mid.loc[immune_mask]

    myeloid_main = df["included_in_myeloid_main"] & df["final_myeloid_label"].isin(MYELOID_MAIN_LABELS)
    state.loc[myeloid_main] = df.loc[myeloid_main, "final_myeloid_label"]
    parent.loc[myeloid_main] = "Myeloid_DC"

    myeloid_uncertain = df["myeloid_adjudication_status"].isin(
        ["sensitivity_only", "exclude_from_myeloid_main", "downgrade_to_myeloid_unspecified"]
    )
    myeloid_uncertain = myeloid_uncertain & df["final_major_lineage"].isin(["Myeloid_DC", "Immune"])
    state.loc[myeloid_uncertain & (~myeloid_main)] = "Myeloid_unspecified"
    parent.loc[myeloid_uncertain & (~myeloid_main)] = "Myeloid_DC"

    state.loc[df["final_major_lineage"].eq("Unknown")] = "Unknown"
    parent.loc[df["final_major_lineage"].eq("Unknown")] = "Unknown"

    df["fraction_state_label"] = state.map(sanitize_label)
    df["fraction_parent_label"] = parent.map(sanitize_label)
    df["is_immune_fraction_cell"] = df["final_major_lineage"].isin(IMMUNE_MAJORS) | df["fraction_parent_label"].isin(IMMUNE_MAJORS)
    df["is_myeloid_uncertain"] = myeloid_uncertain
    return df


def count_by_sample(df: pd.DataFrame, label_col: str, count_family: str) -> pd.DataFrame:
    out = df.groupby(SAMPLE_KEY + [label_col], dropna=False).size().reset_index(name="n_cells")
    out = out.rename(columns={label_col: "label"})
    out.insert(2, "count_family", count_family)
    return out


def build_counts(df: pd.DataFrame) -> pd.DataFrame:
    parts = [
        count_by_sample(df, "final_major_lineage", "major_lineage"),
        count_by_sample(df, "fraction_state_label", "fraction_state"),
        count_by_sample(df[df["is_immune_fraction_cell"]], "fraction_state_label", "immune_state"),
        count_by_sample(df, "fraction_parent_label", "parent_lineage"),
    ]
    return pd.concat(parts, ignore_index=True)


def build_denominators(df: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    sample_base = df.groupby(SAMPLE_KEY, dropna=False).agg(
        total_cells_used=("cell_barcode", "size"),
        immune_cells_used=("is_immune_fraction_cell", "sum"),
        source_h5ad_set=("source_h5ad", lambda s: ";".join(sorted(set(map(str, s))))),
    ).reset_index()
    parent_denoms = df.groupby(SAMPLE_KEY + ["fraction_parent_label"], dropna=False).size().reset_index(name="parent_cells_used")
    state_counts = df.groupby(SAMPLE_KEY + ["fraction_state_label"], dropna=False).size().reset_index(name="state_cells")
    immune_state_counts = (
        df[df["is_immune_fraction_cell"]]
        .groupby(SAMPLE_KEY + ["fraction_state_label"], dropna=False)
        .size()
        .reset_index(name="immune_state_cells")
    )
    return sample_base, parent_denoms, pd.merge(state_counts, immune_state_counts, on=SAMPLE_KEY + ["fraction_state_label"], how="outer").fillna(0)


def feature_dict_row(name: str, family: str, label: str, denom_label: str, main_or_sensitivity: str = "main") -> dict[str, str]:
    denom_def = {
        "all": "all included Step2.5 consensus main cells in sample",
        "immune": "immune included Step2.5 consensus main cells in sample",
        "parent": f"included main cells in parent lineage {denom_label} within sample",
    }[family]
    return {
        "feature_name": name,
        "feature_family": family,
        "numerator_label": label,
        "denominator_label": denom_label,
        "denominator_definition": denom_def,
        "missing_rule": "NA when denominator < 10; zero numerator with valid denominator is 0",
        "low_confidence_rule": "low_confidence_fraction when denominator < 30",
        "source_step": "Step2.5_consensus+Step2.6_myeloid" if "Myeloid" in label or label in MYELOID_MAIN_LABELS else "Step2.5_consensus",
        "main_or_sensitivity": main_or_sensitivity,
    }


def build_fraction_outputs(df: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    sample_base, parent_denoms, state_counts = build_denominators(df)
    samples = sample_base[SAMPLE_KEY].copy()
    matrix = sample_base.copy()
    long_rows: list[dict[str, Any]] = []
    dict_rows: list[dict[str, str]] = []

    all_labels = sorted(df["fraction_state_label"].dropna().unique().tolist())
    immune_labels = sorted(df.loc[df["is_immune_fraction_cell"], "fraction_state_label"].dropna().unique().tolist())
    state_lookup = state_counts.set_index(SAMPLE_KEY + ["fraction_state_label"]).to_dict("index")
    sample_lookup = sample_base.set_index(SAMPLE_KEY).to_dict("index")

    for label in all_labels:
        fname = feature_name("all", label)
        dict_rows.append(feature_dict_row(fname, "all", label, "all_cells"))
        vals = []
        for _, srow in samples.iterrows():
            key = (srow["cohort_id"], srow["sample_id"])
            denom = sample_lookup[key]["total_cells_used"]
            numerator = state_lookup.get((key[0], key[1], label), {}).get("state_cells", 0)
            value = fraction_value(numerator, denom)
            vals.append(value)
            long_rows.append({
                "cohort_id": key[0], "sample_id": key[1], "feature_name": fname, "feature_family": "all",
                "numerator_label": label, "denominator_label": "all_cells", "numerator_count": int(numerator),
                "denominator_count": int(denom), "fraction_value": value,
                "low_confidence_fraction": bool(denom < THRESHOLDS["low_confidence_denominator"]),
                "missing_fraction": bool(denom < THRESHOLDS["missing_denominator"]),
            })
        matrix[fname] = vals

    for label in immune_labels:
        fname = feature_name("immune", label)
        dict_rows.append(feature_dict_row(fname, "immune", label, "immune_cells"))
        vals = []
        for _, srow in samples.iterrows():
            key = (srow["cohort_id"], srow["sample_id"])
            denom = sample_lookup[key]["immune_cells_used"]
            numerator = state_lookup.get((key[0], key[1], label), {}).get("immune_state_cells", 0)
            value = fraction_value(numerator, denom)
            vals.append(value)
            long_rows.append({
                "cohort_id": key[0], "sample_id": key[1], "feature_name": fname, "feature_family": "immune",
                "numerator_label": label, "denominator_label": "immune_cells", "numerator_count": int(numerator),
                "denominator_count": int(denom), "fraction_value": value,
                "low_confidence_fraction": bool(denom < THRESHOLDS["low_confidence_denominator"]),
                "missing_fraction": bool(denom < THRESHOLDS["missing_denominator"]),
            })
        matrix[fname] = vals

    parent_count_lookup = parent_denoms.set_index(SAMPLE_KEY + ["fraction_parent_label"])["parent_cells_used"].to_dict()
    for parent, children in PARENT_CHILDREN.items():
        for child in children:
            if child not in all_labels:
                continue
            fname = feature_name("parent", child, parent)
            dict_rows.append(feature_dict_row(fname, "parent", child, parent))
            vals = []
            for _, srow in samples.iterrows():
                key = (srow["cohort_id"], srow["sample_id"])
                denom = parent_count_lookup.get((key[0], key[1], sanitize_label(parent)), 0)
                numerator = state_lookup.get((key[0], key[1], sanitize_label(child)), {}).get("state_cells", 0)
                value = fraction_value(numerator, denom)
                vals.append(value)
                long_rows.append({
                    "cohort_id": key[0], "sample_id": key[1], "feature_name": fname, "feature_family": "parent",
                    "numerator_label": child, "denominator_label": parent, "numerator_count": int(numerator),
                    "denominator_count": int(denom), "fraction_value": value,
                    "low_confidence_fraction": bool(denom < THRESHOLDS["low_confidence_denominator"]),
                    "missing_fraction": bool(denom < THRESHOLDS["missing_denominator"]),
                })
            matrix[fname] = vals

    return pd.DataFrame(long_rows), matrix, pd.DataFrame(dict_rows).drop_duplicates("feature_name")


def build_qc_flags(df: pd.DataFrame, sample_base: pd.DataFrame, low_cov_sources: set[str]) -> pd.DataFrame:
    total_by_sample = sample_base.set_index(SAMPLE_KEY).copy()
    unknown = df[df["fraction_state_label"].eq("Unknown")].groupby(SAMPLE_KEY).size().rename("unknown_cells")
    immune_unspec = df[df["fraction_state_label"].eq("Immune_unspecified")].groupby(SAMPLE_KEY).size().rename("immune_unspecified_cells")
    conflict = df[df["conflict_type"].astype(str).ne("")].groupby(SAMPLE_KEY).size().rename("conflict_cells")
    myeloid_candidates = df[df["myeloid_adjudication_status"].ne("not_myeloid_candidate")].groupby(SAMPLE_KEY).size().rename("myeloid_candidate_cells")
    myeloid_uncertain = df[df["is_myeloid_uncertain"]].groupby(SAMPLE_KEY).size().rename("myeloid_uncertain_cells")
    myeloid_sens = df[df["myeloid_sensitivity_only"]].groupby(SAMPLE_KEY).size().rename("myeloid_sensitivity_cells")

    flags = total_by_sample.join([unknown, immune_unspec, conflict, myeloid_candidates, myeloid_uncertain, myeloid_sens], how="left").fillna(0).reset_index()
    for c in ["unknown_cells", "immune_unspecified_cells", "conflict_cells", "myeloid_candidate_cells", "myeloid_uncertain_cells", "myeloid_sensitivity_cells"]:
        flags[c] = flags[c].astype(int)
    flags["unknown_fraction"] = flags["unknown_cells"] / flags["total_cells_used"].replace(0, np.nan)
    flags["immune_unspecified_fraction"] = flags["immune_unspecified_cells"] / flags["immune_cells_used"].replace(0, np.nan)
    flags["conflict_fraction"] = flags["conflict_cells"] / flags["total_cells_used"].replace(0, np.nan)
    flags["myeloid_uncertain_fraction"] = flags["myeloid_uncertain_cells"] / flags["myeloid_candidate_cells"].replace(0, np.nan)
    flags["myeloid_sensitivity_fraction"] = flags["myeloid_sensitivity_cells"] / flags["myeloid_candidate_cells"].replace(0, np.nan)
    flags["low_total_cells"] = flags["total_cells_used"] < THRESHOLDS["low_total_cells"]
    flags["low_immune_cells"] = flags["immune_cells_used"] < THRESHOLDS["low_immune_cells"]
    flags["high_unknown_fraction"] = flags["unknown_fraction"].fillna(0) > THRESHOLDS["high_unknown_fraction"]
    flags["high_immune_unspecified_fraction"] = flags["immune_unspecified_fraction"].fillna(0) > THRESHOLDS["high_immune_unspecified_fraction"]
    flags["high_conflict_fraction"] = flags["conflict_fraction"].fillna(0) > THRESHOLDS["high_conflict_fraction"]
    flags["myeloid_annotation_uncertain"] = flags["myeloid_uncertain_fraction"].fillna(0) > THRESHOLDS["myeloid_uncertain_fraction"]
    flags["myeloid_sensitivity_fraction_high"] = flags["myeloid_sensitivity_fraction"].fillna(0) > THRESHOLDS["myeloid_sensitivity_fraction"]
    flags["has_step2_6_coverage_warning"] = flags["source_h5ad_set"].map(lambda x: any(s in low_cov_sources for s in str(x).split(";")))
    flags["myeloid_low_coverage_source"] = flags["has_step2_6_coverage_warning"]
    return flags


def validate_outputs(matrix: pd.DataFrame, long_df: pd.DataFrame, dictionary: pd.DataFrame, flags: pd.DataFrame, expected_samples: int) -> dict[str, Any]:
    feature_cols = [c for c in matrix.columns if c.startswith("frac_")]
    dict_features = set(dictionary["feature_name"].astype(str))
    response_cols = [c for frame in [matrix, long_df, dictionary, flags] for c in frame.columns if any(term in c.lower() for term in FORBIDDEN_RESPONSE_TERMS)]
    all_cols = [c for c in feature_cols if c.startswith("frac_all__")]
    immune_cols = [c for c in feature_cols if c.startswith("frac_immune__")]
    all_valid = matrix["total_cells_used"] >= THRESHOLDS["missing_denominator"]
    all_sum = matrix.loc[all_valid, all_cols].sum(axis=1, skipna=True) if all_cols else pd.Series(dtype=float)
    immune_valid = matrix["immune_cells_used"] >= THRESHOLDS["missing_denominator"]
    immune_sum = matrix.loc[immune_valid, immune_cols].sum(axis=1, skipna=True) if immune_cols else pd.Series(dtype=float)
    validation = {
        "matrix_rows": int(len(matrix)),
        "expected_samples": int(expected_samples),
        "long_fraction_rows": int(len(long_df)),
        "dictionary_rows": int(len(dictionary)),
        "feature_columns": int(len(feature_cols)),
        "sample_count_matches": len(matrix) == expected_samples,
        "all_features_have_dictionary": set(feature_cols).issubset(dict_features),
        "all_fraction_sum_max_abs_error_valid_denominators": float((all_sum - 1.0).abs().max()) if len(all_sum) else np.nan,
        "immune_fraction_sum_max_abs_error": float((immune_sum - 1.0).abs().max()) if len(immune_sum) else np.nan,
        "forbidden_response_cols": sorted(set(response_cols)),
        "provenance_present": all(c in matrix.columns for c in PROVENANCE_COLS) and all(c in flags.columns for c in PROVENANCE_COLS),
    }
    validation["passed"] = (
        validation["sample_count_matches"]
        and validation["all_features_have_dictionary"]
        and (validation["all_fraction_sum_max_abs_error_valid_denominators"] <= 1e-6)
        and (pd.isna(validation["immune_fraction_sum_max_abs_error"]) or validation["immune_fraction_sum_max_abs_error"] <= 1e-6)
        and len(validation["forbidden_response_cols"]) == 0
        and validation["provenance_present"]
    )
    return validation


def write_report(path: Path, validation: dict[str, Any], flags: pd.DataFrame, created_at: str) -> None:
    flag_cols = [
        "low_total_cells",
        "low_immune_cells",
        "high_unknown_fraction",
        "high_immune_unspecified_fraction",
        "high_conflict_fraction",
        "myeloid_annotation_uncertain",
        "myeloid_low_coverage_source",
        "myeloid_sensitivity_fraction_high",
        "has_step2_6_coverage_warning",
    ]
    lines = [
        "# Step2.7 Cell Fraction Feature Factory Report",
        "",
        f"- run_id: `{RUN_ID}`",
        f"- created_at: `{created_at}`",
        "- branch: `fraction_feature_v1`",
        f"- samples: `{validation['matrix_rows']}`",
        f"- feature_columns: `{validation['feature_columns']}`",
        f"- validation_passed: `{validation['passed']}`",
        "",
        "## Compliance",
        "- response fields were not used or propagated.",
        "- no responder/non-responder comparison performed.",
        "- no differential analysis performed.",
        "- no model training performed.",
        "",
        "## QC Flag Counts",
        "",
    ]
    for c in flag_cols:
        lines.append(f"- `{c}`: `{int(flags[c].sum())}`")
    lines += ["", "## Validation", ""]
    for k, v in validation.items():
        lines.append(f"- `{k}`: `{v}`")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--dry-run", action="store_true", help="load inputs and print expected sample count without writing outputs")
    args = parser.parse_args()

    created_at = now_iso()
    consensus, myeloid, low_cov_sources = load_inputs()
    if args.dry_run:
        print(json.dumps({
            "main_rows": len(consensus),
            "samples": consensus[SAMPLE_KEY].drop_duplicates().shape[0],
            "flags_check": consensus.attrs.get("flags_check", {}),
        }, indent=2))
        return

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    annotated = apply_myeloid_overlay(consensus, myeloid)
    counts = build_counts(annotated)
    long_df, matrix, dictionary = build_fraction_outputs(annotated)
    sample_base, _, _ = build_denominators(annotated)
    flags = build_qc_flags(annotated, sample_base, low_cov_sources)

    expected_samples = annotated[SAMPLE_KEY].drop_duplicates().shape[0]
    matrix = add_provenance(matrix, created_at)
    long_df = add_provenance(long_df, created_at)
    counts = add_provenance(counts, created_at)
    flags = add_provenance(flags, created_at)
    dictionary = add_provenance(dictionary, created_at)

    validation = validate_outputs(matrix, long_df, dictionary, flags, expected_samples)

    counts.to_csv(OUT_DIR / "cell_counts_by_sample.csv", index=False)
    long_df.to_csv(OUT_DIR / "cell_fraction_by_sample.csv", index=False)
    matrix.to_parquet(OUT_DIR / "cell_fraction_feature_matrix.parquet", index=False)
    matrix.to_csv(OUT_DIR / "cell_fraction_feature_matrix.csv.gz", index=False, compression="gzip")
    dictionary.to_csv(OUT_DIR / "cell_fraction_feature_dictionary.csv", index=False)
    flags.to_csv(OUT_DIR / "fraction_qc_flags_by_sample.csv", index=False)

    summary = flags.groupby("cohort_id", dropna=False).agg(
        n_samples=("sample_id", "nunique"),
        median_total_cells_used=("total_cells_used", "median"),
        median_immune_cells_used=("immune_cells_used", "median"),
        low_total_cells_samples=("low_total_cells", "sum"),
        low_immune_cells_samples=("low_immune_cells", "sum"),
        high_unknown_fraction_samples=("high_unknown_fraction", "sum"),
        high_conflict_fraction_samples=("high_conflict_fraction", "sum"),
        myeloid_annotation_uncertain_samples=("myeloid_annotation_uncertain", "sum"),
        has_step2_6_coverage_warning_samples=("has_step2_6_coverage_warning", "sum"),
    ).reset_index()
    summary = add_provenance(summary, created_at)
    summary.to_csv(OUT_DIR / "fraction_summary_by_cohort.csv", index=False)

    manifest = {
        "run_id": RUN_ID,
        "created_at": created_at,
        "branch": "fraction_feature_v1",
        "input_manifest_ref": str(INPUT_MANIFEST_REF),
        "input_paths": {
            "consensus_annotation": str(CONSENSUS_PATH),
            "consensus_manifest": str(CONSENSUS_MANIFEST_PATH),
            "myeloid_adjudication": str(MYELOID_PATH),
            "myeloid_manifest": str(MYELOID_MANIFEST_PATH),
            "inclusion_flags": str(FLAGS_PATH),
        },
        "input_hashes": {p.name: sha256_file(p) for p in [CONSENSUS_PATH, CONSENSUS_MANIFEST_PATH, MYELOID_PATH, MYELOID_MANIFEST_PATH, FLAGS_PATH]},
        "thresholds": THRESHOLDS,
        "low_coverage_sources_from_step2_6": sorted(low_cov_sources),
        "flags_consensus_crosscheck": consensus.attrs.get("flags_check", {}),
        "excluded_labels_by_design": {
            "Macro_inflammatory": "Step2.6 myeloid_adjudication_v1 does not emit this as a main myeloid label; Step2.7 excludes it from parent-lineage main fraction features unless Step2.6 is revised.",
        },
        "validation": validation,
        "output_files": {
            "counts": str(OUT_DIR / "cell_counts_by_sample.csv"),
            "long_fractions": str(OUT_DIR / "cell_fraction_by_sample.csv"),
            "feature_matrix_parquet": str(OUT_DIR / "cell_fraction_feature_matrix.parquet"),
            "feature_matrix_csv_gz": str(OUT_DIR / "cell_fraction_feature_matrix.csv.gz"),
            "dictionary": str(OUT_DIR / "cell_fraction_feature_dictionary.csv"),
            "qc_flags": str(OUT_DIR / "fraction_qc_flags_by_sample.csv"),
            "summary_by_cohort": str(OUT_DIR / "fraction_summary_by_cohort.csv"),
            "report": str(OUT_DIR / "fraction_summary_report.md"),
        },
        "response_used": False,
        "responder_nonresponder_comparison_used": False,
        "global_clustering_used": False,
        "differential_analysis_used": False,
        "model_training_used": False,
    }
    with (OUT_DIR / "fraction_decision_manifest.yaml").open("w", encoding="utf-8") as f:
        yaml.safe_dump(manifest, f, sort_keys=False, allow_unicode=True)
    write_report(OUT_DIR / "fraction_summary_report.md", validation, flags, created_at)

    if not validation["passed"]:
        raise RuntimeError(f"Step2.7 validation failed: {validation}")
    print(json.dumps(validation, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
