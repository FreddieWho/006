#!/usr/bin/env python3
"""Targeted Phase4B blocker patch from local rescue/addendum pseudobulk tables.

This script does not rebuild cell-state annotation or raw-count pseudobulk.
It standardizes already-local pseudobulk addendum tables into support/sensitivity
normalized-expression assets, then scores response-blind signatures and limited
TF activities for downstream support use.
"""
from __future__ import annotations

import hashlib
import json
import re
import shutil
from datetime import date
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

ROOT = Path(__file__).resolve().parents[2]
PHASE = ROOT / "results" / "v6_2" / "phase4b_immune_state_feature_construction"
PSEUDOBULK = PHASE / "pseudobulk"
ACTIVITY = PHASE / "activity"
MATRIX = PHASE / "matrix"
AUDIT = PHASE / "audit"
HANDOFF = PHASE / "handoff"
TODAY = str(date.today())
INPUT_VERSION = "frozen_v0"
MAX_TF_FEATURES = 50
RUN_TF_ACTIVITY = False

SOURCES = [
    {
        "cohort_id": "GSE120575_rescue",
        "path": ROOT
        / "results/v6_2/phase2_5_data_onboarding/03_feature_or_layer_build/pseudobulk_addendum.GSE120575_rescue.tsv.gz",
        "orientation": "sample_by_gene",
        "role": "normalized_expression_rescue_support",
    },
    {
        "cohort_id": "GSE229772_rescue",
        "path": ROOT
        / "results/v6_2/phase2_5_data_onboarding/02_conversion_qc/pseudobulk_addendum.GSE229772_rescue.tsv.gz",
        "orientation": "gene_by_sample",
        "role": "support_only_expression_rescue_no_source_proven_response_label",
    },
    {
        "cohort_id": "GSE236581",
        "path": ROOT
        / "results/v6_2/phase2_5_data_onboarding/03_feature_or_layer_build/pseudobulk_addendum.GSE236581.tsv.gz",
        "orientation": "gene_by_sample",
        "role": "support_only_pathologic_responder_pole_crc_dmmr_confounded",
    },
    {
        "cohort_id": "GSE256326",
        "path": ROOT
        / "results/v6_2/phase2_5_data_onboarding/03_feature_or_layer_build/pseudobulk_addendum.GSE256326.tsv.gz",
        "orientation": "gene_by_sample",
        "role": "support_only_expression_addendum",
    },
    {
        "cohort_id": "GSE314072",
        "path": ROOT
        / "results/v6_2/phase2_5_data_onboarding/03_feature_or_layer_build/pseudobulk_addendum.GSE314072.tsv.gz",
        "orientation": "sample_by_gene",
        "role": "support_only_expression_addendum",
    },
]


def write_csv(df: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, index=False)


def write_yaml(obj: dict, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as fh:
        yaml.safe_dump(obj, fh, sort_keys=False, allow_unicode=True)


def read_yaml(path: Path) -> dict:
    if not path.exists():
        return {}
    with path.open() as fh:
        return yaml.safe_load(fh) or {}


def file_sha256(path: Path, nbytes: int = 4_194_304) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        h.update(fh.read(nbytes))
    return h.hexdigest()


def safe_token(value: object) -> str:
    token = re.sub(r"[^A-Za-z0-9]+", "_", str(value).strip()).strip("_").lower()
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
    for bad, repl in replacements.items():
        token = token.replace(bad, repl)
    return token or "unknown"


def normalize_gene_name(value: object) -> str:
    gene = str(value).strip()
    if not gene or gene.lower() in {"nan", "none", "null"}:
        return ""
    return gene.upper()


def collapse_gene_columns(frame: pd.DataFrame) -> pd.DataFrame:
    id_cols = ["cohort_id", "sample_key"]
    genes = [c for c in frame.columns if c not in id_cols]
    gene_part = frame[genes].copy()
    gene_part.columns = [normalize_gene_name(c) for c in gene_part.columns]
    keep = [bool(c) for c in gene_part.columns]
    gene_part = gene_part.loc[:, keep]
    gene_part = gene_part.apply(pd.to_numeric, errors="coerce").fillna(0.0)
    if gene_part.columns.duplicated().any():
        gene_part = gene_part.T.groupby(level=0).mean().T
    out = pd.concat([frame[id_cols].reset_index(drop=True), gene_part.reset_index(drop=True)], axis=1)
    return out


def load_source(src: dict) -> tuple[pd.DataFrame, dict]:
    path = Path(src["path"])
    if not path.exists():
        return pd.DataFrame(), {
            "cohort_id": src["cohort_id"],
            "source_path": str(path),
            "exists": "no",
            "status": "missing",
        }
    raw = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    if src["orientation"] == "sample_by_gene":
        sample_col = raw.columns[0]
        wide = raw.rename(columns={sample_col: "sample_id"}).copy()
    else:
        gene_col = "gene_symbol" if "gene_symbol" in raw.columns else raw.columns[0]
        transposed = raw.set_index(gene_col).T.reset_index()
        wide = transposed.rename(columns={"index": "sample_id"}).copy()

    wide["sample_id"] = wide["sample_id"].astype(str)
    wide.insert(0, "cohort_id", src["cohort_id"])
    wide.insert(1, "sample_key", wide["cohort_id"] + "::" + wide["sample_id"])
    wide = wide.drop(columns=["sample_id"])
    normalized = collapse_gene_columns(wide)
    genes = [c for c in normalized.columns if c not in {"cohort_id", "sample_key"}]
    audit = {
        "cohort_id": src["cohort_id"],
        "source_path": str(path),
        "exists": "yes",
        "status": "loaded",
        "orientation": src["orientation"],
        "role": src["role"],
        "n_samples": int(len(normalized)),
        "n_genes": int(len(genes)),
        "size_bytes": int(path.stat().st_size),
        "partial_sha256_4mb": file_sha256(path),
        "phase4b_use_status": "support_or_sensitivity_normalized_expression",
        "downgrade_reason": src["role"],
    }
    return normalized, audit


def combine_sources() -> tuple[pd.DataFrame, pd.DataFrame]:
    frames = []
    audit_rows = []
    for src in SOURCES:
        frame, audit = load_source(src)
        audit_rows.append(audit)
        if not frame.empty:
            frames.append(frame)
    if not frames:
        raise RuntimeError("No local rescue/addendum pseudobulk table loaded.")
    combined = pd.concat(frames, ignore_index=True, sort=False)
    id_cols = ["cohort_id", "sample_key"]
    gene_cols = [c for c in combined.columns if c not in id_cols]
    combined[id_cols] = combined[id_cols].astype(str)
    combined[gene_cols] = combined[gene_cols].apply(pd.to_numeric, errors="coerce").fillna(0.0)
    combined = combined.drop_duplicates(["cohort_id", "sample_key"], keep="first")
    return combined, pd.DataFrame(audit_rows)


def cohort_zscore(pb: pd.DataFrame, genes: list[str]) -> pd.DataFrame:
    out_parts = []
    for _, group in pb.groupby("cohort_id", sort=False):
        x = group[genes].astype(float)
        mu = x.mean(axis=0)
        sd = x.std(axis=0, ddof=0).replace(0, 1.0)
        z = (x - mu) / sd
        z.insert(0, "sample_key", group["sample_key"].to_numpy())
        z.insert(0, "cohort_id", group["cohort_id"].to_numpy())
        out_parts.append(z)
    return pd.concat(out_parts, ignore_index=True, sort=False)


def load_program_signatures() -> dict[str, dict[str, list[str]]]:
    path = ROOT / "mvp/outputs/program_signatures.json"
    if not path.exists():
        return {}
    raw = json.loads(path.read_text())
    signatures = {}
    for name, spec in raw.items():
        if not isinstance(spec, dict):
            continue
        signatures[name] = {
            "up": [normalize_gene_name(g) for g in spec.get("up_genes", []) if normalize_gene_name(g)],
            "down": [normalize_gene_name(g) for g in spec.get("down_genes", []) if normalize_gene_name(g)],
            "source": str(path),
        }
    return signatures


def score_signatures(pb: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    signatures = load_program_signatures()
    genes = [c for c in pb.columns if c not in {"cohort_id", "sample_key"}]
    z = cohort_zscore(pb, genes)
    matrix = pb[["sample_key", "cohort_id"]].copy()
    coverage_rows = []
    dict_rows = []
    for name, spec in signatures.items():
        up_present = [g for g in spec["up"] if g in z.columns]
        down_present = [g for g in spec["down"] if g in z.columns]
        total_genes = len(set(spec["up"] + spec["down"]))
        present = len(set(up_present + down_present))
        feature_id = f"signature_activity__{safe_token(name)}"
        coverage = present / total_genes if total_genes else 0.0
        coverage_rows.append(
            {
                "feature_id": feature_id,
                "resource": "mvp_program_signatures",
                "n_genes_total": total_genes,
                "n_genes_present": present,
                "coverage": coverage,
                "status": "scored" if present >= 2 else "blocked_low_gene_coverage",
            }
        )
        if present < 2:
            continue
        up_score = z[up_present].mean(axis=1) if up_present else 0.0
        down_score = z[down_present].mean(axis=1) if down_present else 0.0
        matrix[feature_id] = np.asarray(up_score) - np.asarray(down_score)
        dict_rows.append(
            feature_row(
                feature_id=feature_id,
                family="signature_activity",
                source_file="activity/signature_activity_matrix.sample_level.csv",
                source_layer="normalized_rescue_pseudobulk_cohort_zscore",
                method="cohort_zscore_signed_gene_set_mean",
                gene_set=name,
                allowed="support_or_sensitivity",
                notes="response-blind program signature from local rescue/addendum normalized pseudobulk",
            )
        )
    return matrix, pd.DataFrame(coverage_rows), pd.DataFrame(dict_rows)


def score_tf_activity(pb: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    path = (
        ROOT
        / "data/tf_regulon/consensus_tf_regulon_v1/saezlab_tf_regulon_consensus_v1.csv.gz"
    )
    if not path.exists():
        return (
            pb[["sample_key", "cohort_id"]].copy(),
            pd.DataFrame(
                [
                    {
                        "feature_id": "tf_activity",
                        "resource": str(path),
                        "n_genes_total": 0,
                        "n_genes_present": 0,
                        "coverage": 0,
                        "status": "missing_resource",
                    }
                ]
            ),
            pd.DataFrame(),
        )
    genes = [c for c in pb.columns if c not in {"cohort_id", "sample_key"}]
    gene_set = set(genes)
    z = cohort_zscore(pb, genes)
    regulon = pd.read_csv(path, dtype=str, keep_default_na=False)
    regulon = regulon[regulon["used_in_main"].astype(str).str.lower().isin({"true", "yes", "1"})].copy()
    regulon["gene_symbol"] = regulon["gene_symbol"].map(normalize_gene_name)
    regulon["signed_weight"] = pd.to_numeric(regulon["signed_weight"], errors="coerce").fillna(0.0)
    groups = []
    for feature_id, sub in regulon.groupby("feature_id"):
        present = sub[sub["gene_symbol"].isin(gene_set)].copy()
        if len(present) < 5:
            continue
        groups.append((feature_id, sub, present))
    groups = sorted(groups, key=lambda x: len(x[2]), reverse=True)[:MAX_TF_FEATURES]
    matrix = pb[["sample_key", "cohort_id"]].copy()
    coverage_rows = []
    dict_rows = []
    for feature_id, sub, present in groups:
        weights = present.set_index("gene_symbol")["signed_weight"]
        cols = list(weights.index)
        denom = float(np.abs(weights.to_numpy()).sum()) or 1.0
        matrix[feature_id] = (z[cols].to_numpy(float) @ weights.to_numpy(float)) / denom
        coverage_rows.append(
            {
                "feature_id": feature_id,
                "resource": "saezlab_tf_regulon_consensus_v1",
                "n_genes_total": int(sub["gene_symbol"].nunique()),
                "n_genes_present": int(present["gene_symbol"].nunique()),
                "coverage": float(present["gene_symbol"].nunique() / max(sub["gene_symbol"].nunique(), 1)),
                "status": "scored_sensitivity",
            }
        )
        dict_rows.append(
            feature_row(
                feature_id=feature_id,
                family="tf_activity",
                source_file="activity/tf_activity_matrix.sample_level.csv",
                source_layer="normalized_rescue_pseudobulk_cohort_zscore",
                method="cohort_zscore_signed_regulon_mean",
                gene_set=feature_id.replace("tf_activity__", ""),
                allowed="sensitivity_support_only",
                notes="limited TF sensitivity activity from consensus regulon and rescue/addendum expression",
            )
        )
    return matrix, pd.DataFrame(coverage_rows), pd.DataFrame(dict_rows)


def feature_row(
    feature_id: str,
    family: str,
    source_file: str,
    source_layer: str,
    method: str,
    gene_set: str,
    allowed: str,
    notes: str,
) -> dict:
    return {
        "feature_id": feature_id,
        "feature_name": feature_id,
        "feature_family": family,
        "source_file": source_file,
        "source_layer": source_layer,
        "cohort_scope": "support_rescue_addendum",
        "sample_scope": "sample_level",
        "resolution": "sample",
        "cell_state": "",
        "gene_set": gene_set,
        "method": method,
        "allowed_universe": "support_validation_pool",
        "allowed_downstream_use": allowed,
        "is_primary_feature": "no",
        "is_sensitivity_feature": "yes",
        "is_support_only": "yes",
        "is_qc_covariate": "no",
        "is_response_derived": "false",
        "missingness": "0.0",
        "cohort_dominance": "",
        "cancer_dominance": "",
        "platform_dependence": "",
        "cell_count_dependence": "",
        "notes": notes,
        "max_cohort_fraction": "",
        "max_cancer_fraction": "",
        "cell_count_correlation": "",
    }


def update_pseudobulk_registry(pb: pd.DataFrame) -> None:
    path = PSEUDOBULK / "cell_state_specific_pseudobulk_registry.csv"
    registry = pd.read_csv(path, dtype=str, keep_default_na=False) if path.exists() else pd.DataFrame()
    if not registry.empty:
        registry = registry[
            registry["pseudobulk_status"].ne("ready_support_normalized_rescue_addendum_expression")
        ].copy()
    sample_ids = pb["sample_key"].str.split("::", n=1).str[-1]
    add = pd.DataFrame(
        {
            "cohort_id": pb["cohort_id"],
            "sample_key": pb["sample_key"],
            "sample_id": sample_ids,
            "patient_key": pb["sample_key"],
            "timepoint": "unknown",
            "resolution": "sample_level",
            "cell_state": "bulk_sample_expression",
            "n_cells": "",
            "expression_layer_family": "normalized_expression_rescue_addendum",
            "raw_counts_available": "no",
            "normalized_expression_available": "yes",
            "pseudobulk_status": "ready_support_normalized_rescue_addendum_expression",
            "allowed_downstream_use": "support_or_sensitivity_only",
            "block_or_caveat": "not_raw_count;not_cell_state_specific;not_primary_supervised_feature",
        }
    )
    out = pd.concat([registry, add], ignore_index=True, sort=False)
    write_csv(out, path)


def update_support_matrix(activity: pd.DataFrame) -> None:
    path = MATRIX / "immune_state_feature_matrix.support.csv"
    support = pd.read_csv(path, dtype=str, keep_default_na=False) if path.exists() else pd.DataFrame()
    sample_id = activity["sample_key"].str.split("::", n=1).str[-1]
    add = activity.copy()
    add.insert(2, "sample_id", sample_id)
    add.insert(3, "patient_key", add["sample_key"])
    add.insert(4, "patient_id", sample_id)
    add.insert(5, "timepoint", "unknown")
    add.insert(6, "total_cells", "")
    if support.empty:
        merged = add
    else:
        merged = support.merge(add, on="sample_key", how="outer", suffixes=("", "_activity"))
        for col in ["cohort_id", "sample_id", "patient_key", "patient_id", "timepoint", "total_cells"]:
            alt = f"{col}_activity"
            if alt in merged.columns:
                merged[col] = merged[col].replace("", np.nan).fillna(merged[alt]).fillna("")
                merged = merged.drop(columns=[alt])
    write_csv(merged, path)
    by_family = MATRIX / "immune_state_feature_matrix_by_family"
    by_family.mkdir(parents=True, exist_ok=True)
    write_csv(activity, by_family / "support_activity_rescue_addendum.csv")


def update_feature_dictionaries(
    pseudobulk_dict: pd.DataFrame, signature_dict: pd.DataFrame, tf_dict: pd.DataFrame
) -> None:
    pb_rows = pd.DataFrame(
        [
            {
                "feature_id": "pseudobulk_normalized_expression__rescue_addendum_gene_matrix",
                "feature_name": "pseudobulk_normalized_expression__rescue_addendum_gene_matrix",
                "feature_family": "pseudobulk_normalized_expression",
                "source_file": "pseudobulk/pseudobulk_matrix_normalized_expression.parquet",
                "source_layer": "normalized_expression_rescue_addendum",
                "allowed_downstream_use": "support_or_sensitivity_expression_only",
                "is_response_derived": "false",
                "notes": "gene-level normalized/support pseudobulk table; not raw count and not cell-state-specific",
            },
            {
                "feature_id": "pseudobulk_raw_count_blocked_current_phase",
                "feature_name": "pseudobulk_raw_count_blocked_current_phase",
                "feature_family": "pseudobulk_raw_count",
                "source_file": "pseudobulk/pseudobulk_matrix_raw_count.parquet",
                "source_layer": "unavailable",
                "allowed_downstream_use": "blocked",
                "is_response_derived": "false",
                "notes": "raw-count cell-state-specific pseudobulk remains deferred",
            },
        ]
    )
    write_csv(pb_rows, PSEUDOBULK / "pseudobulk_feature_dictionary.csv")

    activity_rows = pd.concat([signature_dict, tf_dict], ignore_index=True, sort=False)
    if activity_rows.empty:
        activity_rows = pd.DataFrame(
            [
                {
                    "feature_id": "activity_blocked_no_scored_features",
                    "feature_name": "activity_blocked_no_scored_features",
                    "feature_family": "signature_pathway_tf_activity",
                    "source_file": "",
                    "source_layer": "",
                    "allowed_downstream_use": "blocked",
                    "is_response_derived": "false",
                    "notes": "No activity features scored.",
                }
            ]
        )
    write_csv(activity_rows, ACTIVITY / "signature_pathway_tf_feature_dictionary.csv")

    main_path = MATRIX / "feature_dictionary_v6_2.csv"
    main = pd.read_csv(main_path, dtype=str, keep_default_na=False) if main_path.exists() else pd.DataFrame()
    drop_families = {
        "pseudobulk",
        "pseudobulk_raw_count",
        "pseudobulk_normalized_expression",
        "signature_activity",
        "tf_activity",
        "signature_pathway_tf_activity",
    }
    if not main.empty:
        main = main[~main["feature_family"].isin(drop_families)].copy()
    full_new = pd.concat([pseudobulk_dict, signature_dict, tf_dict], ignore_index=True, sort=False)
    if main.empty:
        all_cols = list(full_new.columns)
    else:
        all_cols = list(main.columns)
    for col in all_cols:
        if col not in full_new.columns:
            full_new[col] = ""
    full_new = full_new[all_cols]
    out = pd.concat([main, full_new], ignore_index=True, sort=False)
    write_csv(out, main_path)


def support_feature_audit(activity: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    feature_cols = [c for c in activity.columns if c not in {"sample_key", "cohort_id"}]
    n = len(activity)
    cohort_counts = activity["cohort_id"].value_counts(normalize=True)
    max_cohort = float(cohort_counts.max()) if not cohort_counts.empty else 0.0
    missing_rows = []
    tag_rows = []
    for col in feature_cols:
        miss = float(activity[col].isna().mean()) if n else 1.0
        missing_rows.append(
            {
                "feature_id": col,
                "matrix": "support",
                "missingness": miss,
                "non_missing_sample_count": int(activity[col].notna().sum()),
                "non_missing_patient_count": "",
                "cohort_coverage": int(activity["cohort_id"].nunique()),
                "cancer_coverage": "",
                "platform_coverage": "",
                "max_cohort_fraction": max_cohort,
                "max_cancer_fraction": "",
                "platform_dependence": "",
                "cell_count_correlation": "",
                "cell_count_dependence": "not_assessed_support_expression",
                "high_missingness": "yes" if miss > 0.5 else "no",
                "high_cohort_dominance": "yes" if max_cohort > 0.5 else "no",
                "high_cancer_dominance": "not_assessed",
            }
        )
        tag_rows.append(
            {
                "feature_id": col,
                "high_missingness": "yes" if miss > 0.5 else "no",
                "high_cohort_dominance": "yes" if max_cohort > 0.5 else "no",
                "high_cancer_dominance": "not_assessed",
                "cell_count_dependence": "not_assessed_support_expression",
                "max_cohort_fraction": max_cohort,
                "max_cancer_fraction": "",
                "allowed_primary_after_audit": "no",
            }
        )
    missing = pd.DataFrame(missing_rows)
    tags = pd.DataFrame(tag_rows)
    family = (
        missing.assign(feature_family=missing["feature_id"].str.replace(r"__.*", "", regex=True))
        .groupby("feature_family", as_index=False)
        .agg(
            n_features=("feature_id", "count"),
            mean_missingness=("missingness", "mean"),
            high_cohort_dominance_features=("high_cohort_dominance", lambda x: int((x == "yes").sum())),
            high_cancer_dominance_features=("high_cancer_dominance", lambda x: int((x == "yes").sum())),
        )
    )
    return missing, tags, family


def append_audit_tables(activity: pd.DataFrame) -> None:
    missing, tags, family = support_feature_audit(activity)
    for path, add in [
        (AUDIT / "feature_missingness_report.csv", missing),
        (AUDIT / "feature_confounding_tags.csv", tags),
        (AUDIT / "feature_family_qc_summary.csv", family),
    ]:
        old = pd.read_csv(path, dtype=str, keep_default_na=False) if path.exists() else pd.DataFrame()
        if not old.empty and "feature_id" in old.columns:
            old = old[~old["feature_id"].isin(set(add.get("feature_id", [])))].copy()
        elif not old.empty and "feature_family" in old.columns:
            old = old[~old["feature_family"].isin(set(add.get("feature_family", [])))].copy()
        for col in old.columns:
            if col not in add.columns:
                add[col] = ""
        for col in add.columns:
            if col not in old.columns and not old.empty:
                old[col] = ""
        out = pd.concat([old, add[old.columns if not old.empty else add.columns]], ignore_index=True, sort=False)
        write_csv(out, path)


def update_manifests(pb: pd.DataFrame, sig: pd.DataFrame, tf: pd.DataFrame) -> None:
    tf_features = int(max(tf.shape[1] - 2, 0))
    feature_dict_path = MATRIX / "feature_dictionary_v6_2.csv"
    feature_dict_shape = {"rows": 0, "columns": 0}
    if feature_dict_path.exists():
        feature_dict = pd.read_csv(feature_dict_path, dtype=str, keep_default_na=False)
        feature_dict_shape = {"rows": int(len(feature_dict)), "columns": int(len(feature_dict.columns))}
    status_update = {
        "pseudobulk_raw_count": "blocked",
        "pseudobulk_normalized_expression": "support_or_sensitivity",
        "signature_activity": "support_or_sensitivity",
        "tf_activity": "sensitivity_support" if tf_features else "blocked",
        "pathway_activity": "blocked",
    }
    allowed_paths = [
        str(PSEUDOBULK / "pseudobulk_matrix_normalized_expression.parquet"),
        str(PSEUDOBULK / "pseudobulk_matrix_support_only.parquet"),
        str(ACTIVITY / "signature_activity_matrix.sample_level.csv"),
        str(MATRIX / "immune_state_feature_matrix.support.csv"),
    ]
    blocked_paths = [
        str(PSEUDOBULK / "pseudobulk_matrix_raw_count.parquet"),
        str(ACTIVITY / "pathway_activity_matrix.sample_level.csv"),
    ]
    if tf_features:
        allowed_paths.append(str(ACTIVITY / "tf_activity_matrix.sample_level.csv"))
    else:
        blocked_paths.append(str(ACTIVITY / "tf_activity_matrix.sample_level.csv"))
    for path in [HANDOFF / "phase4b_decision_manifest.yaml", HANDOFF / "phase4b_to_phase5_handoff.yaml"]:
        obj = read_yaml(path)
        obj.setdefault("feature_family_status", {}).update(status_update)
        obj["phase4b_blocker_patch"] = {
            "patch_id": "rescue_expression_support_patch",
            "created_at": TODAY,
            "n_support_expression_samples": int(len(pb)),
            "n_signature_features": int(max(sig.shape[1] - 2, 0)),
            "n_tf_features": int(max(tf.shape[1] - 2, 0)),
            "raw_count_pseudobulk_status": "blocked_cell_state_specific_raw_count_deferred",
            "supervised_feature_use": "not_allowed_for_support_expression_patch",
        }
        obj.setdefault("matrix_shapes", {})
        obj["matrix_shapes"]["feature_dictionary"] = feature_dict_shape
        obj["matrix_shapes"]["normalized_rescue_pseudobulk"] = {
            "rows": int(len(pb)),
            "columns": int(len(pb.columns)),
        }
        obj["matrix_shapes"]["signature_activity"] = {
            "rows": int(len(sig)),
            "columns": int(len(sig.columns)),
        }
        obj["matrix_shapes"]["tf_activity"] = {
            "rows": int(len(tf)),
            "columns": int(len(tf.columns)),
        }
        items = list(dict.fromkeys(obj.get("conditional_items", [])))
        items = [x for x in items if x not in {"pseudobulk_expression_aggregation_deferred", "signature_pathway_tf_activity_blocked_current_run"}]
        items.extend(
            [
                "raw_count_cell_state_specific_pseudobulk_deferred",
                "normalized_rescue_pseudobulk_support_available",
                "signature_activity_support_available",
                "tf_activity_deferred",
            ]
        )
        obj["conditional_items"] = list(dict.fromkeys(items))
        obj["phase5_allowed_inputs"] = list(dict.fromkeys(obj.get("phase5_allowed_inputs", []) + allowed_paths))
        obj["phase5_blocked_inputs"] = [x for x in obj.get("phase5_blocked_inputs", []) if x not in allowed_paths]
        obj["phase5_blocked_inputs"] = list(dict.fromkeys(obj.get("phase5_blocked_inputs", []) + blocked_paths))
        obj.setdefault("patch_queue", [])
        obj["patch_queue"] = [
            q
            for q in obj["patch_queue"]
            if q.get("issue") not in {"run_expression_pseudobulk_aggregation"}
        ]
        obj["patch_queue"].extend(
            [
                {
                    "issue": "run_cell_state_specific_raw_count_pseudobulk_aggregation",
                    "required_before": "raw-count module discovery or cell-state-specific expression claim",
                },
                {
                    "issue": "freeze_formal_pathway_resource",
                    "required_before": "pathway activity primary use",
                },
            ]
        )
        write_yaml(obj, path)


def update_report(pb: pd.DataFrame, sig: pd.DataFrame, tf: pd.DataFrame, audit: pd.DataFrame) -> None:
    report = PHASE / "PHASE4B_IMMUNE_STATE_FEATURE_CONSTRUCTION_REPORT.md"
    existing = report.read_text(encoding="utf-8") if report.exists() else ""
    marker = "<!-- phase4b_rescue_expression_patch -->"
    section = f"""{marker}

## Targeted Blocker Patch: Rescue/Addendum Expression Support

**Run date:** {TODAY}

- Standardized local rescue/addendum pseudobulk tables into `pseudobulk_matrix_normalized_expression.parquet`.
- Support expression samples: {len(pb)} across {pb['cohort_id'].nunique()} cohorts.
- Scored response-blind program signatures: {max(sig.shape[1] - 2, 0)} features.
- Scored limited TF sensitivity activities: {max(tf.shape[1] - 2, 0)} features.
- Raw-count cell-state-specific pseudobulk remains blocked/deferred; normalized rescue expression is support/sensitivity only.
- These features are not authorized for primary supervised response-direction claims.

Source cohorts:
{chr(10).join(f"- {r.cohort_id}: {r.phase4b_use_status} ({r.downgrade_reason})" for r in audit.itertuples())}
"""
    if marker in existing:
        existing = existing.split(marker)[0].rstrip() + "\n\n" + section
    else:
        existing = existing.rstrip() + "\n\n" + section
    report.write_text(existing.rstrip() + "\n", encoding="utf-8")


def refresh_output_index() -> None:
    rows = []
    for path in sorted(PHASE.rglob("*")):
        if path.is_file():
            rows.append(
                {
                    "path": str(path),
                    "relative_path": str(path.relative_to(PHASE)),
                    "size_bytes": path.stat().st_size,
                }
            )
    pd.DataFrame(rows).to_csv(HANDOFF / "phase4b_output_index.tsv", sep="\t", index=False)


def main() -> None:
    PSEUDOBULK.mkdir(parents=True, exist_ok=True)
    ACTIVITY.mkdir(parents=True, exist_ok=True)
    MATRIX.mkdir(parents=True, exist_ok=True)
    pb, source_audit = combine_sources()
    print(f"loaded support expression rows={len(pb)} genes={max(pb.shape[1] - 2, 0)}", flush=True)
    write_csv(source_audit, PSEUDOBULK / "pseudobulk_rescue_expression_source_audit.csv")
    normalized_path = PSEUDOBULK / "pseudobulk_matrix_normalized_expression.parquet"
    support_path = PSEUDOBULK / "pseudobulk_matrix_support_only.parquet"
    raw_path = PSEUDOBULK / "pseudobulk_matrix_raw_count.parquet"
    if not normalized_path.exists() or normalized_path.stat().st_size < 1024:
        pb.to_parquet(normalized_path, index=False)
    if not support_path.exists() or support_path.stat().st_size < 1024:
        shutil.copy2(normalized_path, support_path)
    if not raw_path.exists() or raw_path.stat().st_size < 1024:
        pd.DataFrame({"cohort_id": pd.Series(dtype=str), "sample_key": pd.Series(dtype=str)}).to_parquet(
            raw_path, index=False
        )
    print("pseudobulk parquet ready", flush=True)

    update_pseudobulk_registry(pb)
    print("pseudobulk registry ready", flush=True)
    pb_dict = pd.DataFrame(
        [
            feature_row(
                feature_id="pseudobulk_normalized_expression__rescue_addendum_gene_matrix",
                family="pseudobulk_normalized_expression",
                source_file="pseudobulk/pseudobulk_matrix_normalized_expression.parquet",
                source_layer="normalized_expression_rescue_addendum",
                method="local_addendum_pseudobulk_standardization",
                gene_set="",
                allowed="support_or_sensitivity_expression_only",
                notes="normalized/support pseudobulk only; not raw count and not cell-state-specific",
            )
        ]
    )

    sig_matrix, sig_cov, sig_dict = score_signatures(pb)
    print(f"signature scoring ready features={max(sig_matrix.shape[1] - 2, 0)}", flush=True)
    if RUN_TF_ACTIVITY:
        tf_matrix, tf_cov, tf_dict = score_tf_activity(pb)
    else:
        tf_matrix = pb[["sample_key", "cohort_id"]].copy()
        tf_cov = pd.DataFrame(
            [
                {
                    "feature_id": "tf_activity",
                    "resource": "saezlab_tf_regulon_consensus_v1",
                    "n_genes_total": 0,
                    "n_genes_present": 0,
                    "coverage": 0.0,
                    "status": "skipped_lightweight_patch",
                }
            ]
        )
        tf_dict = pd.DataFrame()
    write_csv(sig_matrix, ACTIVITY / "signature_activity_matrix.sample_level.csv")
    write_csv(tf_matrix, ACTIVITY / "tf_activity_matrix.sample_level.csv")
    write_csv(pd.DataFrame(columns=["sample_key", "cohort_id"]), ACTIVITY / "pathway_activity_matrix.sample_level.csv")
    coverage = pd.concat([sig_cov, tf_cov], ignore_index=True, sort=False)
    write_csv(coverage, ACTIVITY / "gene_set_coverage_report.csv")

    activity = sig_matrix.merge(tf_matrix, on=["sample_key", "cohort_id"], how="outer")
    update_support_matrix(activity)
    update_feature_dictionaries(pb_dict, sig_dict, tf_dict)
    append_audit_tables(activity)

    qc = f"""# Signature / Pathway / TF QC Report

**Run date:** {TODAY}

- Input expression family: normalized rescue/addendum pseudobulk.
- Signature features scored: {max(sig_matrix.shape[1] - 2, 0)}
- TF sensitivity features scored: {max(tf_matrix.shape[1] - 2, 0)}
- Pathway activity remains blocked because no formal pathway gene-set resource is frozen.
- Activity scores are response-blind cohort-wise z-score means.
- Use boundary: support/sensitivity only; not primary supervised response feature.
"""
    (ACTIVITY / "signature_pathway_tf_qc_report.md").write_text(qc, encoding="utf-8")

    update_manifests(pb, sig_matrix, tf_matrix)
    update_report(pb, sig_matrix, tf_matrix, source_audit)
    refresh_output_index()
    print(
        f"patched support expression: samples={len(pb)} cohorts={pb['cohort_id'].nunique()} "
        f"signatures={max(sig_matrix.shape[1] - 2, 0)} tf={max(tf_matrix.shape[1] - 2, 0)}"
    )


if __name__ == "__main__":
    main()
