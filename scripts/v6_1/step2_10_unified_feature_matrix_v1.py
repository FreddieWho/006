#!/usr/bin/env python3
from __future__ import annotations

import argparse
import hashlib
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import pyarrow.parquet as pq
import yaml

RUN_ID = "step2_v6_1_0505_0319"
ROOT = Path("/home/huyudi/006")
STEP2_ROOT = ROOT / "results/v6_1/step2" / RUN_ID
FRACTION_DIR = STEP2_ROOT / "06_fraction/fraction_feature_v1"
PSEUDO_DIR = STEP2_ROOT / "07_pseudobulk/pseudobulk_feature_v1_final_v3"
SIGTF_DIR = STEP2_ROOT / "08_signature_pathway_tf/signature_pathway_tf_feature_v1_tf_rescue_v1"
OUT_DIR = STEP2_ROOT / "09_feature_matrix/unified_feature_matrix_v1"
STEP1_DIR = ROOT / "results/v6_1/step1"
PLAN_PATH = ROOT / "results/v6_1/step2/plans/09_feature_matrix/STEP2_10_unified_feature_matrix_and_final_audit_plan.md"

SAMPLE_META = STEP1_DIR / "sample_metadata_master_v6_1.broad_response.csv"
PATIENT_META = STEP1_DIR / "patient_metadata_master_v6_1.broad_response.csv"
PATIENT_SPLIT = STEP1_DIR / "frozen_patient_split_v6_1.csv"

FORBIDDEN_RESPONSE_TERMS = ["response_strict", "response_broad", "response_binary", "response_ordered", "response_raw", "response_"]
PROVENANCE_COLS = ["run_id", "created_at", "input_manifest_ref"]
SIG_META_COLS = ["run_id", "created_at", "input_manifest_ref", "row_id", "cohort_id", "sample_id", "source_h5ad", "aggregation_scope", "cell_state_label", "parent_lineage", "cell_count", "gene_universe", "matrix_role"]
BASE_SAMPLE_COLS = ["run_id", "created_at", "input_manifest_ref", "sample_key", "cohort_id", "sample_id", "source_h5ad_set", "source_h5ad", "step2_9_row_id", "patient_id", "patient_key", "split", "fold_id", "disease", "tissue_source", "sample_type", "timepoint", "timepoint_class", "treatment_context"]


def now_iso() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat()


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def forbidden_response_cols(cols: list[str]) -> list[str]:
    return [c for c in cols if any(term in c.lower() for term in FORBIDDEN_RESPONSE_TERMS)]


def assert_no_response_cols(df: pd.DataFrame, name: str) -> None:
    bad = forbidden_response_cols(df.columns.tolist())
    if bad:
        raise RuntimeError(f"forbidden response columns in {name}: {bad}")


def sample_key(cohort: pd.Series, sample: pd.Series, source: pd.Series) -> pd.Series:
    return cohort.astype(str) + "||" + sample.astype(str) + "||" + source.astype(str)


def normalize_source(x: pd.Series) -> pd.Series:
    return x.fillna("").astype(str).str.strip().str.lower()


def safe_read_manifest(path: Path) -> dict[str, Any]:
    return yaml.safe_load(path.read_text()) if path.exists() else {}


def load_sample_metadata() -> pd.DataFrame:
    meta = pd.read_csv(SAMPLE_META, dtype=str)
    keep = [c for c in meta.columns if c in {
        "cohort_id", "sample_id", "patient_id", "disease", "tissue_source", "sample_type",
        "timepoint", "timepoint_normalized", "treatment_context", "split", "split_label", "fold_id",
        "usable_in_main_analysis", "usable_in_sensitivity_analysis", "exclusion_reason",
    } and not forbidden_response_cols([c])]
    meta = meta[keep].copy()
    if "timepoint" not in meta.columns and "timepoint_normalized" in meta.columns:
        meta["timepoint"] = meta["timepoint_normalized"]
    if "split" not in meta.columns and "split_label" in meta.columns:
        meta["split"] = meta["split_label"]
    return meta.drop_duplicates()


def attach_metadata(frac: pd.DataFrame, meta: pd.DataFrame) -> pd.DataFrame:
    direct_cols = [c for c in meta.columns if c not in {"cohort_id", "sample_id"}]
    meta_direct = meta.drop_duplicates(["cohort_id", "sample_id"])
    out = frac.merge(meta_direct, on=["cohort_id", "sample_id"], how="left")
    out["metadata_join_method"] = np.where(out["patient_id"].notna(), "cohort_sample", "missing")

    missing = out["patient_id"].isna()
    if missing.any():
        counts = meta.groupby("sample_id").size().rename("n")
        unique_ids = counts[counts == 1].index
        meta_unique = meta[meta["sample_id"].isin(unique_ids)].drop_duplicates("sample_id")
        fallback = out.loc[missing, ["sample_id"]].merge(meta_unique, on="sample_id", how="left")
        for c in direct_cols:
            if c in out.columns and c in fallback.columns:
                vals = fallback[c].to_numpy()
                idx = out.index[missing]
                fill_mask = pd.isna(out.loc[idx, c].to_numpy()) & pd.notna(vals)
                out.loc[idx[fill_mask], c] = vals[fill_mask]
        idx = out.index[missing]
        recovered = out.loc[idx, "patient_id"].notna()
        out.loc[idx[recovered], "metadata_join_method"] = "sample_unique_fallback"
    out["metadata_missing"] = out["patient_id"].isna()
    out["patient_id"] = out["patient_id"].fillna(out["cohort_id"].astype(str) + "::" + out["sample_id"].astype(str))
    for c in ["split", "fold_id", "disease", "tissue_source", "sample_type", "timepoint", "treatment_context"]:
        if c not in out.columns:
            out[c] = "unknown"
        out[c] = out[c].fillna("unknown")
    out["patient_key"] = out["cohort_id"].astype(str) + "::" + out["patient_id"].astype(str)
    return out


def classify_timepoint(x: object) -> str:
    s = str(x).strip().lower()
    if any(k in s for k in ["pre", "baseline", "prior", "before"]):
        return "pre"
    if "post" in s or "after" in s:
        return "post"
    if "on" in s or "during" in s or "treat" in s:
        return "on_treatment"
    return "unknown"


def load_fraction(created_at: str) -> tuple[pd.DataFrame, list[str]]:
    frac = pd.read_parquet(FRACTION_DIR / "cell_fraction_feature_matrix.parquet")
    frac = frac.copy()
    frac["source_h5ad_set"] = normalize_source(frac["source_h5ad_set"])
    frac["source_h5ad"] = frac["source_h5ad_set"]
    frac["sample_key"] = sample_key(frac["cohort_id"], frac["sample_id"], frac["source_h5ad_set"])
    frac_feature_cols = [c for c in frac.columns if c.startswith("frac_")]
    frac["run_id"] = RUN_ID
    frac["created_at"] = created_at
    return frac, frac_feature_cols


def load_sigtf_sample() -> tuple[pd.DataFrame, list[str]]:
    pf = pq.ParquetFile(SIGTF_DIR / "signature_pathway_tf_feature_matrix.parquet")
    cols = pf.schema_arrow.names
    bad = forbidden_response_cols(cols)
    if bad:
        raise RuntimeError(f"forbidden response columns in Step2.9 feature matrix: {bad}")
    df = pd.read_parquet(SIGTF_DIR / "signature_pathway_tf_feature_matrix.parquet")
    df = df[df["matrix_role"] == "sample_core"].copy()
    df["source_h5ad"] = normalize_source(df["source_h5ad"])
    df["sample_key"] = sample_key(df["cohort_id"], df["sample_id"], df["source_h5ad"])
    feature_cols = [c for c in df.columns if c not in SIG_META_COLS + ["sample_key"]]
    keep = ["sample_key", "row_id", "source_h5ad", "cell_count"] + feature_cols
    df = df[keep].rename(columns={"row_id": "step2_9_row_id", "cell_count": "pseudobulk_cell_count"})
    return df, feature_cols


def build_sample_matrix(created_at: str) -> tuple[pd.DataFrame, list[str], list[str]]:
    frac, frac_feature_cols = load_fraction(created_at)
    meta = load_sample_metadata()
    sample = attach_metadata(frac, meta)
    sigtf, sigtf_feature_cols = load_sigtf_sample()
    if sigtf["sample_key"].duplicated().any():
        raise RuntimeError("Step2.9 sample_core has duplicated composite sample_key")
    sample = sample.merge(sigtf, on=["sample_key", "source_h5ad"], how="left")
    sample["pseudobulk_available"] = sample["step2_9_row_id"].notna()
    sample["signature_pathway_tf_available"] = sample["pseudobulk_available"]
    sample["pseudobulk_unavailable_reason"] = np.where(
        sample["pseudobulk_available"],
        "available",
        np.where(sample["source_h5ad"].isin(["gse120575.h5ad", "gse229772.h5ad"]), "step2_8_unclear_or_normalized_raw_count_exclusion", "not_matched_to_step2_9_sample_core"),
    )
    sample["timepoint_class"] = sample["timepoint"].map(classify_timepoint)
    qc = pd.read_csv(FRACTION_DIR / "fraction_qc_flags_by_sample.csv")
    qc["source_h5ad_set"] = normalize_source(qc["source_h5ad_set"])
    qc["sample_key"] = sample_key(qc["cohort_id"], qc["sample_id"], qc["source_h5ad_set"])
    qc_cols = [c for c in qc.columns if c not in set(sample.columns) and not forbidden_response_cols([c])]
    sample = sample.merge(qc[["sample_key"] + qc_cols], on="sample_key", how="left")
    feature_cols = frac_feature_cols + sigtf_feature_cols
    metadata_cols = [c for c in BASE_SAMPLE_COLS + ["metadata_join_method", "metadata_missing", "pseudobulk_cell_count", "pseudobulk_available", "signature_pathway_tf_available", "pseudobulk_unavailable_reason"] if c in sample.columns]
    qc_out_cols = [c for c in sample.columns if c.startswith("low_") or c.startswith("high_") or c.startswith("myeloid_") or c in ["has_step2_6_coverage_warning"]]
    ordered = metadata_cols + [c for c in feature_cols if c in sample.columns] + [c for c in qc_out_cols if c not in metadata_cols and c in sample.columns]
    # Keep any useful cell-count columns near the front if not already included.
    for c in ["total_cells_used", "immune_cells_used"]:
        if c in sample.columns and c not in ordered:
            ordered.insert(min(len(metadata_cols), len(ordered)), c)
    sample = sample[ordered].copy()
    assert_no_response_cols(sample, "sample_matrix")
    return sample, feature_cols, qc_out_cols


def build_feature_dictionary(frac_feature_cols: list[str], sigtf_feature_cols: list[str], sample_feature_cols: list[str]) -> pd.DataFrame:
    frac_dict = pd.read_csv(FRACTION_DIR / "cell_fraction_feature_dictionary.csv")
    frac_rows = []
    for r in frac_dict.itertuples(index=False):
        name = getattr(r, "feature_name")
        if name in sample_feature_cols:
            frac_rows.append({
                "feature_name": name,
                "feature_type": "fraction",
                "source_step": "Step2.7_fraction_feature_v1",
                "calculation_method": "cell_count_fraction",
                "denominator": getattr(r, "denominator_label", "sample_cells"),
                "gene_set": "not_applicable",
                "gene_coverage_rule": "not_applicable",
                "missing_rule": getattr(r, "missing_rule", "NA when denominator invalid"),
                "included_in_main": getattr(r, "main_or_sensitivity", "main") == "main",
                "sensitivity_only": getattr(r, "main_or_sensitivity", "main") != "main",
            })
    sig = pd.read_csv(SIGTF_DIR / "signature_pathway_tf_feature_dictionary.csv")
    sig = sig.drop_duplicates("feature_id")
    sig_rows = []
    for r in sig.itertuples(index=False):
        fid = getattr(r, "feature_id")
        if fid in sample_feature_cols:
            fam = getattr(r, "feature_family")
            sig_rows.append({
                "feature_name": fid,
                "feature_type": fam,
                "source_step": "Step2.9_signature_pathway_tf_feature_v1_tf_rescue_v1",
                "calculation_method": "mean_signed_gene_zscore_within_matrix" if fam != "tf_activity" else "saezlab_omnipath_regulon_weighted_signed_mean_gene_zscore",
                "denominator": "core_gene_universe",
                "gene_set": getattr(r, "feature_name"),
                "gene_coverage_rule": f"available {getattr(r, 'n_genes_available')} / expected {getattr(r, 'n_genes_expected')}; coverage {getattr(r, 'gene_coverage')}",
                "missing_rule": "NA if source lacks gene overlap or pseudobulk unavailable; no expression imputation",
                "included_in_main": bool(getattr(r, "included_in_main_feature_matrix")),
                "sensitivity_only": not bool(getattr(r, "included_in_main_feature_matrix")),
            })
    out = pd.DataFrame(frac_rows + sig_rows)
    missing = sorted(set(sample_feature_cols) - set(out["feature_name"]))
    if missing:
        raise RuntimeError(f"feature dictionary missing {len(missing)} sample feature columns, first={missing[:10]}")
    return out.sort_values(["feature_type", "feature_name"]).reset_index(drop=True)


def weighted_mean(df: pd.DataFrame, cols: list[str], weights: pd.Series) -> pd.Series:
    vals = df[cols].astype(float)
    w = pd.to_numeric(weights, errors="coerce").fillna(0).astype(float)
    out = {}
    for c in cols:
        x = vals[c]
        mask = x.notna() & w.gt(0)
        out[c] = np.nan if not mask.any() else float((x[mask] * w[mask]).sum() / w[mask].sum())
    return pd.Series(out)


def build_patient_matrix(sample: pd.DataFrame, feature_cols: list[str]) -> pd.DataFrame:
    records = []
    weights_col = "total_cells_used" if "total_cells_used" in sample.columns else "pseudobulk_cell_count"
    for (cohort_id, patient_id), g in sample.groupby(["cohort_id", "patient_id"], dropna=False):
        rec: dict[str, Any] = {
            "run_id": RUN_ID,
            "created_at": sample["created_at"].iloc[0],
            "patient_key": f"{cohort_id}::{patient_id}",
            "cohort_id": cohort_id,
            "patient_id": patient_id,
            "n_samples_total": int(len(g)),
            "n_samples_pre": int((g["timepoint_class"] == "pre").sum()),
            "n_samples_post": int((g["timepoint_class"] == "post").sum()),
            "has_pre_sample": bool((g["timepoint_class"] == "pre").any()),
            "has_post_sample": bool((g["timepoint_class"] == "post").any()),
            "has_paired_pre_post": bool((g["timepoint_class"] == "pre").any() and (g["timepoint_class"] == "post").any()),
            "split": g["split"].dropna().iloc[0] if "split" in g and g["split"].notna().any() else "unknown",
            "fold_id": g["fold_id"].dropna().iloc[0] if "fold_id" in g and g["fold_id"].notna().any() else "unknown",
            "disease": g["disease"].dropna().iloc[0] if "disease" in g and g["disease"].notna().any() else "unknown",
            "treatment_context": g["treatment_context"].dropna().iloc[0] if "treatment_context" in g and g["treatment_context"].notna().any() else "unknown",
        }
        summaries: dict[str, pd.Series] = {}
        for tp in ["pre", "post"]:
            sub = g[g["timepoint_class"] == tp]
            if len(sub):
                summaries[f"{tp}_unweighted_mean"] = sub[feature_cols].astype(float).mean(axis=0, skipna=True)
                summaries[f"{tp}_cell_count_weighted_mean"] = weighted_mean(sub, feature_cols, sub[weights_col])
            else:
                summaries[f"{tp}_unweighted_mean"] = pd.Series({c: np.nan for c in feature_cols})
                summaries[f"{tp}_cell_count_weighted_mean"] = pd.Series({c: np.nan for c in feature_cols})
        for c in feature_cols:
            pre_u = summaries["pre_unweighted_mean"][c]
            post_u = summaries["post_unweighted_mean"][c]
            pre_w = summaries["pre_cell_count_weighted_mean"][c]
            post_w = summaries["post_cell_count_weighted_mean"][c]
            rec[f"pre_unweighted_mean__{c}"] = pre_u
            rec[f"pre_cell_count_weighted_mean__{c}"] = pre_w
            rec[f"post_unweighted_mean__{c}"] = post_u
            rec[f"post_cell_count_weighted_mean__{c}"] = post_w
            rec[f"delta_post_minus_pre_unweighted_mean__{c}"] = post_u - pre_u if pd.notna(post_u) and pd.notna(pre_u) else np.nan
            rec[f"delta_post_minus_pre_cell_count_weighted_mean__{c}"] = post_w - pre_w if pd.notna(post_w) and pd.notna(pre_w) else np.nan
        records.append(rec)
    patient = pd.DataFrame(records)
    assert_no_response_cols(patient, "patient_matrix")
    return patient


def build_missingness(sample: pd.DataFrame, feature_cols: list[str]) -> pd.DataFrame:
    rows = []
    for c in feature_cols:
        rows.append({
            "scope": "overall",
            "cohort_id": "all",
            "feature_name": c,
            "n_rows": int(len(sample)),
            "n_missing": int(sample[c].isna().sum()),
            "missing_rate": float(sample[c].isna().mean()),
        })
    for cohort, g in sample.groupby("cohort_id", dropna=False):
        for c in feature_cols:
            rows.append({
                "scope": "cohort",
                "cohort_id": cohort,
                "feature_name": c,
                "n_rows": int(len(g)),
                "n_missing": int(g[c].isna().sum()),
                "missing_rate": float(g[c].isna().mean()),
            })
    return pd.DataFrame(rows)


def build_sample_qc(sample: pd.DataFrame, feature_cols: list[str]) -> pd.DataFrame:
    qc_cols = [c for c in sample.columns if c.startswith("low_") or c.startswith("high_") or c.startswith("myeloid_") or c in ["has_step2_6_coverage_warning"]]
    out_cols = [c for c in ["sample_key", "cohort_id", "sample_id", "source_h5ad", "patient_id", "pseudobulk_available", "signature_pathway_tf_available", "pseudobulk_unavailable_reason", "metadata_join_method", "metadata_missing"] if c in sample.columns]
    qc = sample[out_cols + qc_cols].copy()
    qc["feature_missing_rate"] = sample[feature_cols].isna().mean(axis=1)
    return qc


def build_cohort_eligibility(sample: pd.DataFrame, feature_cols: list[str]) -> pd.DataFrame:
    rows = []
    for cohort, g in sample.groupby("cohort_id", dropna=False):
        n = len(g)
        n_pb = int(g["pseudobulk_available"].sum())
        mean_missing = float(g[feature_cols].isna().mean(axis=1).mean())
        if n_pb == 0:
            status = "annotation_qc_only"
            reason = "no_pseudobulk_expression_features_available"
        elif n_pb < n:
            status = "sensitivity_only"
            reason = "partial_pseudobulk_expression_feature_availability"
        elif n < 3:
            status = "sensitivity_only"
            reason = "low_sample_count"
        else:
            status = "main"
            reason = "fraction_and_expression_features_available"
        rows.append({
            "cohort_id": cohort,
            "eligibility": status,
            "reason": reason,
            "n_sample_rows": int(n),
            "n_pseudobulk_available": int(n_pb),
            "n_pseudobulk_unavailable": int(n - n_pb),
            "mean_feature_missing_rate": mean_missing,
        })
    return pd.DataFrame(rows).sort_values(["eligibility", "cohort_id"])


def write_outputs(sample: pd.DataFrame, patient: pd.DataFrame, feature_dict: pd.DataFrame, missing: pd.DataFrame, qc: pd.DataFrame, eligibility: pd.DataFrame, feature_cols: list[str], created_at: str) -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    sample.to_parquet(OUT_DIR / "unified_feature_matrix_by_sample.parquet", index=False)
    sample.to_csv(OUT_DIR / "unified_feature_matrix_by_sample.csv.gz", index=False)
    patient.to_parquet(OUT_DIR / "unified_feature_matrix_by_patient.parquet", index=False)
    patient.to_csv(OUT_DIR / "unified_feature_matrix_by_patient.csv.gz", index=False)
    feature_dict.to_csv(OUT_DIR / "feature_dictionary_v6_1.csv", index=False)
    missing.to_csv(OUT_DIR / "feature_missingness_report.csv", index=False)
    qc.to_csv(OUT_DIR / "sample_feature_QC_flags.csv", index=False)
    eligibility.to_csv(OUT_DIR / "cohort_feature_eligibility.csv", index=False)

    inputs = [
        FRACTION_DIR / "cell_fraction_feature_matrix.parquet",
        FRACTION_DIR / "cell_fraction_feature_dictionary.csv",
        FRACTION_DIR / "fraction_qc_flags_by_sample.csv",
        FRACTION_DIR / "fraction_decision_manifest.yaml",
        PSEUDO_DIR / "pseudobulk_decision_manifest.yaml",
        PSEUDO_DIR / "pseudobulk_sample_metadata.csv",
        SIGTF_DIR / "signature_pathway_tf_feature_matrix.parquet",
        SIGTF_DIR / "signature_pathway_tf_feature_dictionary.csv",
        SIGTF_DIR / "signature_pathway_tf_decision_manifest.yaml",
        SAMPLE_META,
        PATIENT_META,
        PATIENT_SPLIT,
        PLAN_PATH,
    ]
    frac_manifest = safe_read_manifest(FRACTION_DIR / "fraction_decision_manifest.yaml")
    pseudo_manifest = safe_read_manifest(PSEUDO_DIR / "pseudobulk_decision_manifest.yaml")
    sig_manifest = safe_read_manifest(SIGTF_DIR / "signature_pathway_tf_decision_manifest.yaml")
    sample_bad = forbidden_response_cols(sample.columns.tolist())
    patient_bad = forbidden_response_cols(patient.columns.tolist())
    validation = {
        "sample_rows": int(len(sample)),
        "sample_rows_expected": 2292,
        "sample_rows_match": int(len(sample)) == 2292,
        "patient_rows": int(len(patient)),
        "feature_columns": int(len(feature_cols)),
        "feature_dictionary_rows": int(len(feature_dict)),
        "feature_dictionary_covers_features": set(feature_cols).issubset(set(feature_dict["feature_name"])),
        "pseudobulk_unavailable_rows": int((~sample["pseudobulk_available"]).sum()),
        "pseudobulk_unavailable_rows_expected": 79,
        "forbidden_response_cols_sample": sample_bad,
        "forbidden_response_cols_patient": patient_bad,
        "metadata_missing_rows": int(sample["metadata_missing"].sum()),
        "has_patient_pre_post_delta_columns": any(c.startswith("pre_unweighted_mean__") for c in patient.columns) and any(c.startswith("post_unweighted_mean__") for c in patient.columns) and any(c.startswith("delta_post_minus_pre_unweighted_mean__") for c in patient.columns),
    }
    pass_status = all([
        validation["sample_rows_match"],
        validation["feature_dictionary_covers_features"],
        validation["pseudobulk_unavailable_rows"] == 79,
        not sample_bad,
        not patient_bad,
        validation["has_patient_pre_post_delta_columns"],
    ])
    gate = "pass_with_expected_pseudobulk_exclusions_and_metadata_fallbacks" if pass_status else "blocked"
    manifest = {
        "run_id": RUN_ID,
        "step_id": "Step2.10_unified_feature_matrix_v1",
        "created_at": created_at,
        "gate_status": gate,
        "input_hashes": {str(p): sha256_file(p) for p in inputs if p.exists()},
        "upstream_gates": {
            "step2_7_fraction": frac_manifest.get("validation", {}).get("passed"),
            "step2_8_pseudobulk": pseudo_manifest.get("gate_status"),
            "step2_9_signature_pathway_tf": sig_manifest.get("gate_status"),
        },
        "sample_universe_policy": "Step2.7 fraction matrix left-joined to Step2.9 sample_core by cohort_id+sample_id+source_h5ad",
        "pseudobulk_unavailable_sources": ["gse120575.h5ad", "gse229772.h5ad"],
        "response_usage": "forbidden_not_used",
        "forbidden_operations": ["response_leakage", "differential_analysis", "model_training", "global_clustering"],
        "validation": validation,
        "outputs": {
            "sample_matrix": str(OUT_DIR / "unified_feature_matrix_by_sample.parquet"),
            "patient_matrix": str(OUT_DIR / "unified_feature_matrix_by_patient.parquet"),
            "feature_dictionary": str(OUT_DIR / "feature_dictionary_v6_1.csv"),
            "missingness_report": str(OUT_DIR / "feature_missingness_report.csv"),
            "cohort_eligibility": str(OUT_DIR / "cohort_feature_eligibility.csv"),
            "sample_qc_flags": str(OUT_DIR / "sample_feature_QC_flags.csv"),
        },
    }
    (OUT_DIR / "step2_10_decision_manifest.yaml").write_text(yaml.safe_dump(manifest, sort_keys=False), encoding="utf-8")
    report = [
        "# Step2.10 Unified Feature Matrix Summary", "",
        f"- created_at: `{created_at}`",
        f"- gate_status: `{gate}`",
        f"- sample matrix: `{len(sample)}` rows x `{len(sample.columns)}` columns",
        f"- patient matrix: `{len(patient)}` rows x `{len(patient.columns)}` columns",
        f"- feature columns: `{len(feature_cols)}`",
        f"- pseudobulk unavailable rows retained: `{validation['pseudobulk_unavailable_rows']}`",
        f"- metadata missing rows after fallback: `{validation['metadata_missing_rows']}`",
        "- response usage: `forbidden_not_used`",
        "- differential analysis: `not_run`",
        "- global clustering: `not_run`",
        "- model training: `not_run`", "",
        "## Required Conditions", "",
        "- Sample universe is Step2.7 fraction matrix; no fraction rows were dropped.",
        "- Step2.9 features were left-joined by composite sample key, not sample_id alone.",
        "- 79 expected pseudobulk-unavailable samples from `gse120575.h5ad` and `gse229772.h5ad` are retained and flagged.",
    ]
    text = "\n".join(report) + "\n"
    (OUT_DIR / "step2_feature_summary.md").write_text(text, encoding="utf-8")
    (OUT_DIR / "step2_final_audit_report.md").write_text(text, encoding="utf-8")
    reports_dir = STEP2_ROOT / "10_reports"
    reports_dir.mkdir(parents=True, exist_ok=True)
    (reports_dir / "step2_feature_summary.md").write_text(text, encoding="utf-8")
    (reports_dir / "step2_final_audit_report.md").write_text(text, encoding="utf-8")


def run() -> None:
    created_at = now_iso()
    sample, feature_cols, _ = build_sample_matrix(created_at)
    # Feature columns can be duplicated if upstream dictionaries contain role duplicates; enforce unique order.
    feature_cols = [c for c in dict.fromkeys(feature_cols) if c in sample.columns]
    feature_dict = build_feature_dictionary([c for c in feature_cols if c.startswith("frac_")], [c for c in feature_cols if c.startswith(("signature__", "pathway__", "tf_activity__"))], feature_cols)
    patient = build_patient_matrix(sample, feature_cols)
    missing = build_missingness(sample, feature_cols)
    qc = build_sample_qc(sample, feature_cols)
    eligibility = build_cohort_eligibility(sample, feature_cols)
    write_outputs(sample, patient, feature_dict, missing, qc, eligibility, feature_cols, created_at)
    print("Step2.10 complete")
    print(f"out_dir={OUT_DIR}")
    print(f"sample_shape={sample.shape}")
    print(f"patient_shape={patient.shape}")
    print(f"feature_cols={len(feature_cols)}")
    print(f"pseudobulk_unavailable={(~sample['pseudobulk_available']).sum()}")


def parse_args() -> argparse.Namespace:
    return argparse.ArgumentParser(description="Step2.10 unified feature matrix v1").parse_args()


if __name__ == "__main__":
    parse_args()
    run()
