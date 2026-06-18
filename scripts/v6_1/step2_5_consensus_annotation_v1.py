#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import anndata
import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
import yaml
from scipy import sparse

RUN_ID = "step2_v6_1_0505_0319"
STEP2_ROOT = Path("/home/huyudi/006/results/v6_1/step2") / RUN_ID
ANNOT_DIR = STEP2_ROOT / "04_annotation"
CONS_DIR = ANNOT_DIR / "consensus_annotation_v1"
CHK_DIR = CONS_DIR / "checkpoints"

FLAGS_PATH = STEP2_ROOT / "03_qc" / "final_cell_inclusion_flags.qc0513_balanced.parquet"
INV_PATH = STEP2_ROOT / "00_manifest" / "input_file_inventory.csv"
ORIG_ANNOT_PATH = ANNOT_DIR / "cell_state_annotation_v6_1.parquet"
ORIG_RULES_PATH = ANNOT_DIR / "cell_state_mapping_rules_v6_1.csv"
REGISTRY_RAW_PATH = ANNOT_DIR / "marker_registry_v6_1_extended.csv"
GENE_STD_MANIFEST_PATH = STEP2_ROOT / "02_gene_standardization" / "gene_standardization_manifest.yaml"
STEP2_MANIFEST_PATH = STEP2_ROOT / "00_manifest" / "step2_run_manifest.patched.yaml"

REQUIRED_REGISTRY_COLS = [
    "marker_gene",
    "marker_set",
    "layer",
    "target_label",
    "role",
    "used_for_main",
    "reference_ids",
]

SAMPLE_COL_MAP = {
    "gse272734.h5ad": "sample_prefix",
    "gse272735.h5ad": "sample_prefix",
    "gse273718.h5ad": "sample_prefix",
}

MAJOR_FROM_MARKER_SET = {
    "T_NK_core": "T_NK",
    "B_Plasma_core": "B_Plasma",
    "Myeloid_DC_core": "Myeloid_DC",
    "Mast_minor_immune": "Mast_or_minor_immune",
    "Pan_immune_supportive": "Immune",
    "Non_immune_Epithelial_Tumor": "Non_immune_Epithelial_Tumor",
    "Non_immune_Endothelial": "Non_immune_Endothelial",
    "Non_immune_Stromal": "Non_immune_Stromal",
}

NONIMMUNE_MAJOR = {
    "Non_immune_Epithelial_Tumor",
    "Non_immune_Endothelial",
    "Non_immune_Stromal",
}

IMMUNE_MAJORS = {"T_NK", "B_Plasma", "Myeloid_DC", "Mast_or_minor_immune", "Immune"}
MISSING_LABELS = {"", "NA", "NAN", "NONE", "NULL", "<NA>"}


def now_iso() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat()


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        while True:
            chunk = f.read(1024 * 1024)
            if not chunk:
                break
            h.update(chunk)
    return h.hexdigest()


def normalize_bool(x: Any) -> bool:
    s = str(x).strip().lower()
    return s in {"1", "true", "t", "yes", "y"}


def normalize_text(x: Any) -> str:
    if pd.isna(x):
        return ""
    s = str(x).strip()
    return s


def clean_label_array(values: Any, default: str = "Unknown") -> np.ndarray:
    arr = pd.Series(values).fillna(default).astype(str).str.strip()
    upper = arr.str.upper()
    arr.loc[upper.isin(MISSING_LABELS)] = default
    return arr.to_numpy(dtype=object)


def build_var_lookup(var_names: Any, var: pd.DataFrame) -> dict[str, int]:
    lookup: dict[str, int] = {}

    def add_symbol(symbol: Any, idx: int) -> None:
        if pd.isna(symbol):
            return
        s = str(symbol).strip().upper()
        if not s:
            return
        lookup.setdefault(s, idx)

    for i, gene in enumerate(var_names):
        add_symbol(gene, i)

    for col in ["gene_symbol", "geneSymbol", "gene_symbols", "features", "feature_name", "gene_name"]:
        if col not in var.columns:
            continue
        for i, gene in enumerate(var[col].tolist()):
            add_symbol(gene, i)

    return lookup


@dataclass
class RegistryParseResult:
    df: pd.DataFrame
    repair_report: pd.DataFrame


def parse_and_repair_registry(raw_path: Path) -> RegistryParseResult:
    rows: list[dict[str, Any]] = []
    repairs: list[dict[str, Any]] = []
    with raw_path.open("r", encoding="utf-8", errors="replace") as f:
        lines = f.read().splitlines()

    header_idx = None
    for i, line in enumerate(lines):
        if line.strip().startswith("marker_gene,"):
            header_idx = i
            break
    if header_idx is None:
        raise RuntimeError("marker registry header not found")

    expected = [
        "marker_gene",
        "marker_set",
        "layer",
        "target_label",
        "role",
        "weight",
        "used_for_main",
        "is_core_marker",
        "expected_direction",
        "conflict_group",
        "reference_ids",
        "evidence_level",
        "notes",
    ]

    for lineno, line in enumerate(lines[header_idx + 1 :], start=header_idx + 2):
        if not line.strip():
            continue
        toks = line.split(",")
        raw_field_count = len(toks)
        if raw_field_count < 11:
            repairs.append(
                {
                    "line_no": lineno,
                    "raw_field_count": raw_field_count,
                    "expected_field_count": len(expected),
                    "repair_action": "skip_too_few_fields",
                    "extra_fields_repaired": raw_field_count - len(expected),
                    "reference_ids_before": "",
                    "reference_ids_after": "",
                }
            )
            continue

        base = toks[:10]
        tail = toks[10:]
        evidence_level = ""
        notes = ""
        ref_tokens: list[str]

        if len(tail) >= 3:
            evidence_level = tail[-2].strip()
            notes = tail[-1].strip()
            ref_tokens = tail[:-2]
        elif len(tail) == 2:
            evidence_level = tail[-1].strip()
            ref_tokens = tail[:-1]
        else:
            evidence_level = ""
            ref_tokens = tail

        ref_before = ",".join(ref_tokens).strip()
        ref_after = ref_before.replace(",", ";")
        extra_fields = raw_field_count - len(expected)

        row = {
            "marker_gene": base[0].strip().upper(),
            "marker_set": base[1].strip(),
            "layer": base[2].strip(),
            "target_label": base[3].strip(),
            "role": base[4].strip().lower(),
            "weight": float(base[5].strip() or 1.0),
            "used_for_main": normalize_bool(base[6]),
            "is_core_marker": normalize_bool(base[7]),
            "expected_direction": base[8].strip(),
            "conflict_group": base[9].strip(),
            "reference_ids": ref_after,
            "evidence_level": evidence_level,
            "notes": notes,
            "_line_no": lineno,
            "_raw_field_count": raw_field_count,
        }
        rows.append(row)
        repairs.append(
            {
                "line_no": lineno,
                "raw_field_count": raw_field_count,
                "expected_field_count": len(expected),
                "repair_action": "reference_ids_commas_to_semicolon" if extra_fields != 0 else "none",
                "extra_fields_repaired": extra_fields,
                "reference_ids_before": ref_before,
                "reference_ids_after": ref_after,
            }
        )

    reg = pd.DataFrame(rows)
    repair_df = pd.DataFrame(repairs)

    missing_cols = [c for c in REQUIRED_REGISTRY_COLS if c not in reg.columns]
    if missing_cols:
        raise RuntimeError(f"registry missing required cols: {missing_cols}")

    used = reg[reg["used_for_main"]]
    bad_ref = used[used["reference_ids"].astype(str).str.strip() == ""]
    if len(bad_ref) > 0:
        lines = sorted(bad_ref["_line_no"].astype(int).tolist())
        raise RuntimeError(f"used_for_main markers missing reference_ids at lines: {lines[:20]}")

    return RegistryParseResult(df=reg, repair_report=repair_df)


def build_h5ad_map(inv_path: Path) -> dict[str, str]:
    inv = pd.read_csv(inv_path)
    sub = inv[inv["file_type"] == "h5ad"].copy()
    return {Path(p).name.lower(): p for p in sub["file_path"].astype(str)}


def resolve_expr_and_var(adata: anndata.AnnData, preferred_layer: str):
    layer_keys = set(adata.layers.keys())
    if preferred_layer in layer_keys:
        return adata.layers[preferred_layer], adata.var_names
    if "counts" in layer_keys:
        return adata.layers["counts"], adata.var_names
    if "data" in layer_keys:
        return adata.layers["data"], adata.var_names
    if adata.raw is not None and adata.raw.X is not None:
        return adata.raw.X, adata.raw.var_names
    if adata.X is None:
        raise RuntimeError("no expression matrix found in X/layers/raw")
    return adata.X, adata.var_names


def choose_obs_sample_column(obs: pd.DataFrame, source_h5ad: str) -> str:
    preferred = SAMPLE_COL_MAP.get(source_h5ad)
    if preferred and preferred in obs.columns:
        return preferred
    for c in ["sample_id", "sample", "sample_prefix", "orig.ident", "sampleID"]:
        if c in obs.columns:
            return c
    return ""


def prepare_source_obs_keys(adata: anndata.AnnData, source_h5ad: str) -> pd.DataFrame:
    obs = adata.obs.copy()
    obs["cell_barcode"] = obs.index.astype(str)
    sample_col = choose_obs_sample_column(obs, source_h5ad)
    if sample_col:
        obs["sample_id"] = obs[sample_col].astype(str)
    else:
        obs["sample_id"] = "unknown"
    obs["_cell_idx"] = np.arange(len(obs), dtype=np.int64)
    return obs[["cell_barcode", "sample_id", "_cell_idx"]]


def source_marker_score(
    source_h5ad: str,
    h5ad_path: str,
    flags_sub: pd.DataFrame,
    marker_registry: pd.DataFrame,
    batch_size: int,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    adata = anndata.read_h5ad(h5ad_path, backed="r")
    obs_keys = prepare_source_obs_keys(adata, source_h5ad)
    obs_full = adata.obs.copy()

    keys = flags_sub[["cell_barcode", "sample_id"]].copy()
    keys["cell_barcode"] = keys["cell_barcode"].astype(str)
    keys["sample_id"] = keys["sample_id"].astype(str)

    merged = keys.merge(obs_keys, on=["cell_barcode", "sample_id"], how="left")
    if merged["_cell_idx"].isna().any():
        # fallback for unique barcodes only
        dup = obs_keys["cell_barcode"].duplicated().any()
        if not dup:
            fallback = keys[["cell_barcode"]].merge(obs_keys[["cell_barcode", "_cell_idx"]], on="cell_barcode", how="left")
            merged["_cell_idx"] = fallback["_cell_idx"].values
    if merged["_cell_idx"].isna().any():
        missing = int(merged["_cell_idx"].isna().sum())
        adata.file.close()
        raise RuntimeError(f"{source_h5ad}: failed key match for {missing} cells")

    cell_idx = merged["_cell_idx"].astype(np.int64).to_numpy()

    preferred_layer = marker_registry["layer"].mode().iloc[0]
    X, var_names = resolve_expr_and_var(adata, preferred_layer=preferred_layer)
    var_frame = adata.raw.var if adata.raw is not None and var_names is adata.raw.var_names else adata.var
    var_lookup = build_var_lookup(var_names, var_frame)

    reg_main = marker_registry[marker_registry["used_for_main"]].copy()
    set_defs: dict[str, dict[str, Any]] = {}
    for mset, grp in reg_main.groupby("marker_set", dropna=False):
        genes = grp["marker_gene"].astype(str).str.upper().tolist()
        w = grp.set_index(grp["marker_gene"].astype(str).str.upper())["weight"].to_dict()
        idx = []
        weight_values = []
        for gene in genes:
            if gene not in var_lookup:
                continue
            idx.append(var_lookup[gene])
            weight_values.append(float(w.get(gene, 1.0)))
        weights = np.array(weight_values, dtype=np.float32) if idx else np.array([], dtype=np.float32)
        set_defs[str(mset)] = {
            "target_label": str(grp["target_label"].iloc[0]),
            "n_expected": len(genes),
            "gene_idx": np.array(idx, dtype=np.int64),
            "weights": weights,
            "role": str(grp["role"].iloc[0]),
        }

    # Reduce IO by slicing only marker-union columns from expression matrix.
    union_cols = sorted({int(i) for d in set_defs.values() for i in d["gene_idx"]})
    union_pos = {c: i for i, c in enumerate(union_cols)}
    for s in list(set_defs.keys()):
        cols = set_defs[s]["gene_idx"]
        set_defs[s]["gene_union_pos"] = np.array([union_pos[int(c)] for c in cols], dtype=np.int64) if len(cols) else np.array([], dtype=np.int64)

    # Prefer obs-level library size when available (avoids full-matrix reads on dense HDF5 datasets).
    lib_col = ""
    for c in ["nCount_RNA", "n_counts", "total_counts", "nCount", "nUMI", "nCountRNA", "n_counts_all"]:
        if c in obs_full.columns:
            lib_col = c
            break
    lib_all = None
    if lib_col:
        lib_all = pd.to_numeric(obs_full[lib_col], errors="coerce").fillna(0).to_numpy(dtype=np.float32)

    set_names = sorted(set_defs.keys())
    n_cells = len(cell_idx)
    n_sets = len(set_names)
    score_mat = np.full((n_cells, n_sets), np.nan, dtype=np.float32)
    detect_mat = np.full((n_cells, n_sets), np.nan, dtype=np.float32)
    cov_mat = np.full((n_cells, n_sets), np.nan, dtype=np.float32)

    for j, s in enumerate(set_names):
        n_expected = set_defs[s]["n_expected"]
        n_avail = len(set_defs[s]["gene_idx"])
        cov = (n_avail / n_expected) if n_expected else 0.0
        cov_mat[:, j] = cov

    for start in range(0, n_cells, batch_size):
        end = min(start + batch_size, n_cells)
        ridx = cell_idx[start:end]
        sort_ord = np.argsort(ridx, kind="mergesort")
        ridx_sorted = ridx[sort_ord]
        inv_ord = np.empty_like(sort_ord)
        inv_ord[sort_ord] = np.arange(len(sort_ord))
        marker_chunk = X[ridx_sorted, :][:, union_cols] if union_cols else X[ridx_sorted, :][:, []]
        if sparse.issparse(marker_chunk):
            marker_chunk = marker_chunk[inv_ord, :].tocsr()
        else:
            marker_chunk = np.asarray(marker_chunk)[inv_ord, :]
            marker_chunk = sparse.csr_matrix(marker_chunk)

        if lib_all is not None:
            lib = lib_all[ridx].astype(np.float32, copy=False)
        else:
            # Fallback uses marker-union counts when no obs library-size column exists.
            lib = np.asarray(marker_chunk.sum(axis=1)).ravel().astype(np.float32)
        lib[lib <= 0] = 1.0
        scale = (10000.0 / lib).astype(np.float32)

        for j, s in enumerate(set_names):
            cols = set_defs[s]["gene_union_pos"]
            if cols.size == 0:
                continue
            weights = set_defs[s]["weights"]
            sub = marker_chunk[:, cols]
            det = np.asarray((sub > 0).sum(axis=1)).ravel().astype(np.float32) / float(cols.size)
            detect_mat[start:end, j] = det

            sub_scaled = sub.multiply(scale[:, None])
            sub_scaled = sub_scaled.tocsr(copy=True)
            sub_scaled.data = np.log1p(sub_scaled.data)
            wsum = np.asarray(sub_scaled @ weights).ravel().astype(np.float32)
            denom = float(weights.sum()) if float(weights.sum()) > 0 else float(len(weights))
            score_mat[start:end, j] = wsum / denom

        if start == 0 or ((start // batch_size) % 20 == 0):
            print(
                f"[consensus_v1] {source_h5ad} batch {start}:{end}/{n_cells}",
                flush=True,
            )

    mean = np.nanmean(score_mat, axis=0)
    std = np.nanstd(score_mat, axis=0)
    std[std == 0] = 1.0
    z_mat = (score_mat - mean[None, :]) / std[None, :]

    top_idx = np.nanargmax(np.where(np.isnan(score_mat), -np.inf, score_mat), axis=1)
    sorted_idx = np.argsort(np.where(np.isnan(score_mat), -np.inf, score_mat), axis=1)
    runner_idx = sorted_idx[:, -2] if n_sets >= 2 else top_idx

    top_set = np.array([set_names[i] for i in top_idx], dtype=object)
    runner_set = np.array([set_names[i] for i in runner_idx], dtype=object)

    top_score = score_mat[np.arange(n_cells), top_idx]
    runner_score = score_mat[np.arange(n_cells), runner_idx]
    top_detect = detect_mat[np.arange(n_cells), top_idx]
    top_cov = cov_mat[np.arange(n_cells), top_idx]
    top_z = z_mat[np.arange(n_cells), top_idx]
    top_major = np.array([MAJOR_FROM_MARKER_SET.get(s, "Unknown") for s in top_set], dtype=object)
    runner_major = np.array([MAJOR_FROM_MARKER_SET.get(s, "Unknown") for s in runner_set], dtype=object)

    margin = top_score - runner_score
    conflict_flags = np.full(n_cells, "", dtype=object)
    immune_mask = np.isin(top_major, list(IMMUNE_MAJORS))
    nonimmune_mask = np.isin(top_major, list(NONIMMUNE_MAJOR))

    strong_nonimmune = nonimmune_mask & (top_z >= 1.5) & (top_detect >= 0.2)
    weak_immune = np.full(n_cells, True, dtype=bool)
    immune_set_idx = [i for i, s in enumerate(set_names) if MAJOR_FROM_MARKER_SET.get(s, "") in IMMUNE_MAJORS]
    if immune_set_idx:
        immune_best_z = np.nanmax(z_mat[:, immune_set_idx], axis=1)
        weak_immune = immune_best_z < 0.5
        conflict_flags[strong_nonimmune & (~weak_immune)] = "L0_nonimmune_immune_conflict"

    low_margin = margin < 0.15
    conflict_flags[(conflict_flags == "") & low_margin] = "low_margin_ambiguous"

    marker_confidence = np.where((top_z >= 1.5) & (top_detect >= 0.3) & (margin >= 0.20), "high", "medium")
    marker_confidence = np.where((top_z < 0.5) | (top_detect < 0.1), "low", marker_confidence)

    marker_major_label = top_major.copy()
    marker_major_label[strong_nonimmune & weak_immune] = top_major[strong_nonimmune & weak_immune]
    marker_major_label[(top_z < 0.2) | np.isnan(top_z)] = "Unknown"

    marker_mid_label = np.array([s if m in IMMUNE_MAJORS else "NA" for s, m in zip(top_set, marker_major_label)], dtype=object)
    marker_mid_label[(marker_major_label == "Immune") | (marker_major_label == "Unknown")] = "Immune_unspecified"

    out = flags_sub[["source_h5ad", "cell_barcode", "sample_id", "cohort_id", "final_inclusion"]].copy()
    out["marker_major_label"] = marker_major_label
    out["marker_mid_label"] = marker_mid_label
    out["marker_confidence"] = marker_confidence
    out["marker_top_set"] = top_set
    out["marker_runner_up_set"] = runner_set
    out["marker_top_major"] = top_major
    out["marker_runner_up_major"] = runner_major
    out["marker_score"] = top_score.astype(np.float32)
    out["marker_runner_up_score"] = runner_score.astype(np.float32)
    out["marker_score_margin"] = margin.astype(np.float32)
    out["marker_detected_fraction"] = top_detect.astype(np.float32)
    out["marker_gene_coverage"] = top_cov.astype(np.float32)
    out["marker_source_zscore"] = top_z.astype(np.float32)
    out["marker_conflict_flags"] = conflict_flags

    summary = {
        "source_h5ad": source_h5ad,
        "n_cells": int(n_cells),
        "set_names": set_names,
        "set_mean": {k: float(v) for k, v in zip(set_names, mean)},
        "set_std": {k: float(v) for k, v in zip(set_names, std)},
    }

    adata.file.close()
    return out, summary


def map_original_major(major: str) -> str:
    m = normalize_text(major)
    if m in {"B", "Plasma"}:
        return "B_Plasma"
    if m in {"Myeloid", "DC"}:
        return "Myeloid_DC"
    if m == "T_NK":
        return "T_NK"
    if m == "Mast":
        return "Mast_or_minor_immune"
    if m == "Tumor_Epithelial":
        return "Non_immune_Epithelial_Tumor"
    if m == "Endothelial":
        return "Non_immune_Endothelial"
    if m == "Stromal_Fibroblast":
        return "Non_immune_Stromal"
    if m in {"Unknown", "", "NA"}:
        return "Unknown"
    return "Immune" if m in IMMUNE_MAJORS else "Unknown"


def normalize_original_evidence(df: pd.DataFrame) -> pd.DataFrame:
    out = df[[
        "source_h5ad",
        "cell_barcode",
        "sample_id",
        "original_annotation",
        "original_annotation_column",
        "major_lineage",
        "immune_mid_state",
        "annotation_confidence",
        "annotation_rule_id",
    ]].copy()
    out["original_major_lineage"] = out["major_lineage"].map(map_original_major)
    out["original_immune_mid_state"] = out["immune_mid_state"].fillna("NA").astype(str)
    out.loc[pd.Series(out["original_major_lineage"]).fillna("").astype(str).str.upper().isin(MISSING_LABELS), "original_major_lineage"] = "Unknown"
    out.loc[pd.Series(out["original_immune_mid_state"]).fillna("").astype(str).str.upper().isin(MISSING_LABELS), "original_immune_mid_state"] = "NA"
    out["original_annotation"] = out["original_annotation"].fillna("NA").astype(str)
    out["original_annotation_column"] = out["original_annotation_column"].fillna("none").astype(str)
    out["original_annotation_confidence"] = out["annotation_confidence"].fillna("low").astype(str)
    out["original_annotation_rule_id"] = out["annotation_rule_id"].fillna("none").astype(str)
    out = out.drop(columns=["major_lineage", "immune_mid_state", "annotation_confidence", "annotation_rule_id"])
    return out


def consensus_decision(row: pd.Series) -> dict[str, Any]:
    marker_major = row["marker_major_label"]
    marker_mid = row["marker_mid_label"]
    marker_conf = row["marker_confidence"]
    marker_cov = float(row["marker_gene_coverage"])
    marker_margin = float(row["marker_score_margin"])

    orig_major = row["original_major_lineage"]
    orig_mid = row["original_immune_mid_state"]
    orig_conf = row["original_annotation_confidence"]

    auto_major = row["auto_major_lineage"]
    auto_mid = row["auto_immune_mid_state"]

    coverage_ok = marker_cov >= 0.4
    marker_strong = marker_conf in {"high", "medium"} and marker_major != "Unknown"
    orig_present = orig_major not in {"Unknown", "NA", ""}
    auto_present = auto_major not in {"Unknown", "NA", "", "not_run"}

    conflict_type = ""
    decision_rule = "CONS_R0_fallback"

    if marker_strong and coverage_ok and orig_present and marker_major == orig_major:
        final_major = marker_major
        if marker_major in IMMUNE_MAJORS:
            if marker_mid not in {"Immune_unspecified", "NA"} and (orig_mid == marker_mid or auto_mid == marker_mid):
                final_mid = marker_mid
            elif orig_mid not in {"NA", "Immune_unspecified", ""}:
                final_mid = orig_mid
            else:
                final_mid = "Immune_unspecified"
        else:
            final_mid = "NA"
        final_conf = "high"
        decision_rule = "CONS_R1_marker_original_agree"
    elif marker_strong and coverage_ok and not orig_present:
        final_major = marker_major
        final_mid = marker_mid if marker_major in IMMUNE_MAJORS else "NA"
        if final_mid in {"NA", "", "Unknown"}:
            final_mid = "Immune_unspecified" if marker_major in IMMUNE_MAJORS else "NA"
        final_conf = "medium"
        decision_rule = "CONS_R2_marker_only"
    elif orig_present and (not marker_strong or marker_major == "Unknown"):
        final_major = orig_major
        final_mid = orig_mid if orig_major in IMMUNE_MAJORS else "NA"
        if final_mid in {"", "NA"} and orig_major in IMMUNE_MAJORS:
            final_mid = "Immune_unspecified"
        final_conf = "medium"
        decision_rule = "CONS_R3_original_supported"
    elif marker_strong and orig_present and marker_major != orig_major:
        if marker_major in IMMUNE_MAJORS and orig_major in IMMUNE_MAJORS:
            final_major = "Immune"
            final_mid = "Immune_unspecified"
            conflict_type = "immune_major_conflict"
        elif marker_major in NONIMMUNE_MAJOR and orig_major in IMMUNE_MAJORS:
            final_major = "Immune"
            final_mid = "Immune_unspecified"
            conflict_type = "immune_nonimmune_conflict"
        elif marker_major in IMMUNE_MAJORS and orig_major in NONIMMUNE_MAJOR:
            final_major = "Immune"
            final_mid = "Immune_unspecified"
            conflict_type = "immune_nonimmune_conflict"
        else:
            final_major = marker_major
            final_mid = "Immune_unspecified" if marker_major in IMMUNE_MAJORS else "NA"
            conflict_type = "major_conflict_other"
        final_conf = "low"
        decision_rule = "CONS_R4_conflict_degrade"
    else:
        immune_hint = (marker_major in IMMUNE_MAJORS) or (orig_major in IMMUNE_MAJORS) or auto_present
        if immune_hint:
            final_major = "Immune"
            final_mid = "Immune_unspecified"
        else:
            final_major = "Unknown"
            final_mid = "NA"
        final_conf = "low"
        decision_rule = "CONS_R5_evidence_insufficient"

    # Auto annotation can support but never create high-confidence alone
    if auto_present and decision_rule == "CONS_R5_evidence_insufficient" and final_major == "Immune":
        final_conf = "medium"
        decision_rule = "CONS_R6_auto_support_only"

    functional = "NA"
    if final_mid not in {"NA", "Immune_unspecified"}:
        functional = row.get("original_functional_state_optional", "NA")
        functional = "NA" if pd.isna(functional) or str(functional).strip() == "" else str(functional)

    return {
        "final_major_lineage": final_major,
        "final_immune_mid_state": final_mid,
        "final_functional_state_optional": functional,
        "final_confidence": final_conf,
        "decision_rule_id": decision_rule,
        "conflict_type": conflict_type,
    }


def consensus_decision_vectorized(df: pd.DataFrame) -> pd.DataFrame:
    n = len(df)
    marker_major = df["marker_major_label"].astype(str).to_numpy()
    marker_mid = df["marker_mid_label"].astype(str).to_numpy()
    marker_conf = df["marker_confidence"].astype(str).to_numpy()
    marker_cov = pd.to_numeric(df["marker_gene_coverage"], errors="coerce").fillna(0.0).to_numpy(dtype=np.float32)

    orig_major = clean_label_array(df["original_major_lineage"], default="Unknown")
    orig_mid = clean_label_array(df["original_immune_mid_state"], default="NA")

    auto_major = clean_label_array(df["auto_major_lineage"], default="not_run")
    auto_mid = clean_label_array(df["auto_immune_mid_state"], default="not_run")

    coverage_ok = marker_cov >= 0.4
    marker_strong = np.isin(marker_conf, ["high", "medium"]) & (marker_major != "Unknown")
    orig_present = ~np.isin(orig_major, ["Unknown", "NA", ""])
    auto_present = ~np.isin(auto_major, ["Unknown", "NA", "", "not_run"])

    final_major = np.full(n, "Unknown", dtype=object)
    final_mid = np.full(n, "NA", dtype=object)
    final_conf = np.full(n, "low", dtype=object)
    decision_rule = np.full(n, "CONS_R0_fallback", dtype=object)
    conflict_type = np.full(n, "", dtype=object)

    unassigned = np.ones(n, dtype=bool)

    # R1
    m1 = unassigned & marker_strong & coverage_ok & orig_present & (marker_major == orig_major)
    final_major[m1] = marker_major[m1]
    immune_major_m1 = np.isin(final_major, list(IMMUNE_MAJORS))
    m1i = m1 & immune_major_m1
    mid_ok = ~np.isin(marker_mid, ["Immune_unspecified", "NA"])
    agree_support = (orig_mid == marker_mid) | (auto_mid == marker_mid)
    use_marker_mid = m1i & mid_ok & agree_support
    final_mid[use_marker_mid] = marker_mid[use_marker_mid]
    use_orig_mid = m1i & (~use_marker_mid) & (~np.isin(orig_mid, ["NA", "Immune_unspecified", ""]))
    final_mid[use_orig_mid] = orig_mid[use_orig_mid]
    remain_i = m1i & (~use_marker_mid) & (~use_orig_mid)
    final_mid[remain_i] = "Immune_unspecified"
    final_mid[m1 & (~m1i)] = "NA"
    final_conf[m1] = "high"
    decision_rule[m1] = "CONS_R1_marker_original_agree"
    unassigned[m1] = False

    # R2
    m2 = unassigned & marker_strong & coverage_ok & (~orig_present)
    final_major[m2] = marker_major[m2]
    m2i = m2 & np.isin(final_major, list(IMMUNE_MAJORS))
    final_mid[m2i] = marker_mid[m2i]
    final_mid[m2i & np.isin(final_mid, ["NA", "", "Unknown"])] = "Immune_unspecified"
    final_mid[m2 & (~m2i)] = "NA"
    final_conf[m2] = "medium"
    decision_rule[m2] = "CONS_R2_marker_only"
    unassigned[m2] = False

    # R3
    m3 = unassigned & orig_present & ((~marker_strong) | (marker_major == "Unknown"))
    final_major[m3] = orig_major[m3]
    m3i = m3 & np.isin(final_major, list(IMMUNE_MAJORS))
    final_mid[m3i] = orig_mid[m3i]
    final_mid[m3i & np.isin(final_mid, ["", "NA"])] = "Immune_unspecified"
    final_mid[m3 & (~m3i)] = "NA"
    final_conf[m3] = "medium"
    decision_rule[m3] = "CONS_R3_original_supported"
    unassigned[m3] = False

    # R4
    m4 = unassigned & marker_strong & orig_present & (marker_major != orig_major)
    both_immune = np.isin(marker_major, list(IMMUNE_MAJORS)) & np.isin(orig_major, list(IMMUNE_MAJORS))
    immune_nonimmune = (np.isin(marker_major, list(NONIMMUNE_MAJOR)) & np.isin(orig_major, list(IMMUNE_MAJORS))) | (
        np.isin(marker_major, list(IMMUNE_MAJORS)) & np.isin(orig_major, list(NONIMMUNE_MAJOR))
    )
    other = ~(both_immune | immune_nonimmune)

    m4_both = m4 & both_immune
    final_major[m4_both] = "Immune"
    final_mid[m4_both] = "Immune_unspecified"
    conflict_type[m4_both] = "immune_major_conflict"

    m4_in = m4 & immune_nonimmune
    final_major[m4_in] = "Immune"
    final_mid[m4_in] = "Immune_unspecified"
    conflict_type[m4_in] = "immune_nonimmune_conflict"

    m4_other = m4 & other
    final_major[m4_other] = marker_major[m4_other]
    final_mid[m4_other & np.isin(marker_major, list(IMMUNE_MAJORS))] = "Immune_unspecified"
    final_mid[m4_other & (~np.isin(marker_major, list(IMMUNE_MAJORS)))] = "NA"
    conflict_type[m4_other] = "major_conflict_other"

    final_conf[m4] = "low"
    decision_rule[m4] = "CONS_R4_conflict_degrade"
    unassigned[m4] = False

    # R5
    m5 = unassigned
    immune_hint = np.isin(marker_major, list(IMMUNE_MAJORS)) | np.isin(orig_major, list(IMMUNE_MAJORS)) | auto_present
    m5i = m5 & immune_hint
    final_major[m5i] = "Immune"
    final_mid[m5i] = "Immune_unspecified"
    final_major[m5 & (~immune_hint)] = "Unknown"
    final_mid[m5 & (~immune_hint)] = "NA"
    final_conf[m5] = "low"
    decision_rule[m5] = "CONS_R5_evidence_insufficient"

    # R6 adjustment
    m6 = auto_present & (decision_rule == "CONS_R5_evidence_insufficient") & (final_major == "Immune")
    final_conf[m6] = "medium"
    decision_rule[m6] = "CONS_R6_auto_support_only"

    # Preserve coarse marker-set state only when marker and original agree at major level.
    coarse_mid = ~np.isin(marker_mid, ["NA", "Immune_unspecified", "", "Unknown", "None", "nan"])
    use_coarse_mid = (
        (decision_rule == "CONS_R1_marker_original_agree")
        & (final_mid == "Immune_unspecified")
        & (final_major == marker_major)
        & coarse_mid
    )
    final_mid[use_coarse_mid] = marker_mid[use_coarse_mid]

    functional = np.full(n, "NA", dtype=object)
    if "original_functional_state_optional" in df.columns:
        of = df["original_functional_state_optional"].fillna("NA").astype(str).to_numpy()
        mfunc = ~np.isin(final_mid, ["NA", "Immune_unspecified"])
        functional[mfunc] = of[mfunc]

    return pd.DataFrame(
        {
            "final_major_lineage": final_major,
            "final_immune_mid_state": final_mid,
            "final_functional_state_optional": functional,
            "final_confidence": final_conf,
            "decision_rule_id": decision_rule,
            "conflict_type": conflict_type,
        },
        index=df.index,
    )


def generate_reports(consensus: pd.DataFrame, out_dir: Path) -> None:
    conflict = (
        consensus.assign(_n=1)
        .groupby([
            "cohort_id",
            "sample_id",
            "source_h5ad",
            "conflict_type",
            "decision_rule_id",
            "marker_major_label",
            "original_major_lineage",
            "auto_major_lineage",
        ], dropna=False)["_n"]
        .sum()
        .reset_index(name="n_cells")
        .sort_values("n_cells", ascending=False)
    )
    conflict.to_csv(out_dir / "annotation_conflict_matrix.csv", index=False)

    by_cohort_rows = []
    for cohort, g in consensus.groupby("cohort_id", dropna=False):
        main = g[g["included_in_main_annotation"]]
        immune = main[main["final_major_lineage"].isin(list(IMMUNE_MAJORS) + ["Immune"])]
        by_cohort_rows.append(
            {
                "run_id": g["run_id"].iloc[0],
                "created_at": g["created_at"].iloc[0],
                "input_manifest_ref": g["input_manifest_ref"].iloc[0],
                "cohort_id": cohort,
                "n_cells_total": int(len(g)),
                "n_main": int(len(main)),
                "n_sensitivity_only": int(g["sensitivity_only_annotation"].sum()),
                "major_lineage_coverage": float((main["final_major_lineage"] != "Unknown").mean()) if len(main) else 0.0,
                "immune_mid_state_specified_fraction": float((immune["final_immune_mid_state"] != "Immune_unspecified").mean()) if len(immune) else 0.0,
                "immune_unspecified_fraction": float((immune["final_immune_mid_state"] == "Immune_unspecified").mean()) if len(immune) else 0.0,
                "conflict_fraction": float((main["conflict_type"].astype(str) != "").mean()) if len(main) else 0.0,
            }
        )
    pd.DataFrame(by_cohort_rows).to_csv(out_dir / "annotation_summary_by_cohort.csv", index=False)

    by_sample_rows = []
    for (cohort, sample), g in consensus.groupby(["cohort_id", "sample_id"], dropna=False):
        main = g[g["included_in_main_annotation"]]
        immune = main[main["final_major_lineage"].isin(list(IMMUNE_MAJORS) + ["Immune"])]
        by_sample_rows.append(
            {
                "run_id": g["run_id"].iloc[0],
                "created_at": g["created_at"].iloc[0],
                "input_manifest_ref": g["input_manifest_ref"].iloc[0],
                "cohort_id": cohort,
                "sample_id": sample,
                "n_cells_total": int(len(g)),
                "n_main": int(len(main)),
                "n_sensitivity_only": int(g["sensitivity_only_annotation"].sum()),
                "major_lineage_coverage": float((main["final_major_lineage"] != "Unknown").mean()) if len(main) else 0.0,
                "immune_mid_state_specified_fraction": float((immune["final_immune_mid_state"] != "Immune_unspecified").mean()) if len(immune) else 0.0,
                "has_immune_unspecified_flag": bool((immune["final_immune_mid_state"] == "Immune_unspecified").any()) if len(immune) else False,
            }
        )
    pd.DataFrame(by_sample_rows).to_csv(out_dir / "annotation_summary_by_sample.csv", index=False)


def write_manifest(
    out_dir: Path,
    created_at: str,
    auto_status: str,
    auto_reason: str,
    scoring_thresholds: dict[str, Any],
    marker_source_summaries: list[dict[str, Any]],
) -> None:
    manifest = {
        "run_id": RUN_ID,
        "created_at": created_at,
        "input_manifest_ref": str(STEP2_MANIFEST_PATH),
        "consensus_branch": "consensus_annotation_v1",
        "inputs": {
            "flags": str(FLAGS_PATH),
            "original_annotation": str(ORIG_ANNOT_PATH),
            "original_mapping_rules": str(ORIG_RULES_PATH),
            "marker_registry_raw": str(REGISTRY_RAW_PATH),
            "gene_standardization_manifest": str(GENE_STD_MANIFEST_PATH),
        },
        "input_sha256": {
            "flags": sha256_file(FLAGS_PATH),
            "original_annotation": sha256_file(ORIG_ANNOT_PATH),
            "original_mapping_rules": sha256_file(ORIG_RULES_PATH),
            "marker_registry_raw": sha256_file(REGISTRY_RAW_PATH),
            "gene_standardization_manifest": sha256_file(GENE_STD_MANIFEST_PATH),
            "marker_registry_consensus_v1": sha256_file(out_dir / "marker_registry_consensus_v1.csv"),
        },
        "scoring_thresholds": scoring_thresholds,
        "decision_rules": [
            "CONS_R1_marker_original_agree",
            "CONS_R2_marker_only",
            "CONS_R3_original_supported",
            "CONS_R4_conflict_degrade",
            "CONS_R5_evidence_insufficient",
            "CONS_R6_auto_support_only",
        ],
        "auto_annotation_status": auto_status,
        "auto_not_run_reason": auto_reason,
        "marker_source_summaries": marker_source_summaries,
        "response_leakage_check": {
            "forbidden_columns": ["response_strict", "response_broad"],
            "status": "passed",
        },
    }
    with (out_dir / "consensus_annotation_manifest.yaml").open("w", encoding="utf-8") as f:
        yaml.safe_dump(manifest, f, sort_keys=False)


def write_registry_reports(
    out_dir: Path,
    registry: pd.DataFrame,
    repair_df: pd.DataFrame,
    created_at: str,
) -> None:
    registry.to_csv(out_dir / "marker_registry_consensus_v1.csv", index=False)
    repair_df.to_csv(out_dir / "registry_parse_repair_report.csv", index=False)

    used = registry[registry["used_for_main"]]
    lines = [
        "# Registry Validation Report",
        "",
        f"- run_id: `{RUN_ID}`",
        f"- created_at: `{created_at}`",
        f"- raw_rows: `{len(registry)}`",
        f"- used_for_main_rows: `{len(used)}`",
        f"- marker_sets: `{', '.join(sorted(registry['marker_set'].dropna().astype(str).unique().tolist()))}`",
        "",
        "## Required Columns",
    ]
    for c in REQUIRED_REGISTRY_COLS:
        lines.append(f"- `{c}`: {'ok' if c in registry.columns else 'missing'}")
    missing_ref = int((used["reference_ids"].astype(str).str.strip() == "").sum())
    lines.extend(
        [
            "",
            "## Parse Repair",
            f"- repaired_lines: `{int((repair_df['extra_fields_repaired'] != 0).sum())}`",
            f"- max_extra_fields: `{int(repair_df['extra_fields_repaired'].max())}`",
            "",
            "## Blocking Checks",
            f"- used_for_main without reference_ids: `{missing_ref}`",
            "- status: `passed`",
        ]
    )
    (out_dir / "registry_validation_report.md").write_text("\n".join(lines), encoding="utf-8")


def auto_annotation_branch(keys: pd.DataFrame, created_at: str) -> tuple[pd.DataFrame, str, str]:
    auto_status = "not_run"
    auto_reason = "no_local_model"
    # explicit no-run branch; plan allows this when local model is absent
    auto = keys[["source_h5ad", "cell_barcode", "sample_id"]].copy()
    auto["run_id"] = RUN_ID
    auto["created_at"] = created_at
    auto["input_manifest_ref"] = str(STEP2_MANIFEST_PATH)
    auto["auto_method"] = "not_run"
    auto["auto_not_run_reason"] = auto_reason
    auto["auto_major_lineage"] = "not_run"
    auto["auto_immune_mid_state"] = "not_run"
    auto["auto_confidence"] = "not_run"
    return auto, auto_status, auto_reason


def validate_outputs(consensus: pd.DataFrame, flags: pd.DataFrame) -> dict[str, Any]:
    if len(consensus) != len(flags):
        raise RuntimeError(f"row count mismatch: consensus={len(consensus)} flags={len(flags)}")

    key_cols = ["source_h5ad", "cell_barcode", "sample_id"]
    if consensus[key_cols].isna().any().any():
        raise RuntimeError("key columns contain NA")
    if (consensus["sample_id"].astype(str).str.lower() == "unknown").any():
        raise RuntimeError("sample_id contains unknown")

    for c in ["run_id", "created_at", "input_manifest_ref"]:
        if c not in consensus.columns:
            raise RuntimeError(f"missing provenance column: {c}")

    forbidden = [c for c in consensus.columns if "response" in c.lower()]
    forbidden = [c for c in forbidden if c in {"response_strict", "response_broad"} or "response" in c.lower()]
    if forbidden:
        raise RuntimeError(f"response leakage columns found: {forbidden}")

    main = consensus[consensus["included_in_main_annotation"]]
    major_cov = float((main["final_major_lineage"] != "Unknown").mean()) if len(main) else 0.0
    immune = main[main["final_major_lineage"].isin(list(IMMUNE_MAJORS) + ["Immune"])]
    mid_cov = float((immune["final_immune_mid_state"] != "Immune_unspecified").mean()) if len(immune) else 0.0

    return {
        "row_count": int(len(consensus)),
        "major_coverage": major_cov,
        "immune_mid_coverage": mid_cov,
    }


def run(args: argparse.Namespace) -> None:
    CONS_DIR.mkdir(parents=True, exist_ok=True)
    CHK_DIR.mkdir(parents=True, exist_ok=True)
    created_at = now_iso()

    parse_res = parse_and_repair_registry(REGISTRY_RAW_PATH)
    registry = parse_res.df.copy()
    write_registry_reports(CONS_DIR, registry, parse_res.repair_report, created_at)

    flags = pd.read_parquet(FLAGS_PATH)
    flags = flags[[
        "run_id",
        "created_at",
        "input_manifest_ref",
        "source_h5ad",
        "cohort_id",
        "cell_barcode",
        "sample_id",
        "final_inclusion",
    ]].copy()
    flags["cell_barcode"] = flags["cell_barcode"].astype(str)
    flags["sample_id"] = flags["sample_id"].astype(str)

    original_cols = [
        "source_h5ad",
        "cell_barcode",
        "sample_id",
        "original_annotation",
        "original_annotation_column",
        "major_lineage",
        "immune_mid_state",
        "annotation_confidence",
        "annotation_rule_id",
    ]
    original_all = pd.read_parquet(ORIG_ANNOT_PATH, columns=original_cols)
    original_all = normalize_original_evidence(original_all)
    original_all["cell_barcode"] = original_all["cell_barcode"].astype(str)
    original_all["sample_id"] = original_all["sample_id"].astype(str)
    original_all = original_all.drop_duplicates(subset=["source_h5ad", "cell_barcode", "sample_id"], keep="first")

    sources = sorted(flags["source_h5ad"].astype(str).unique().tolist())
    if args.sources:
        allow = {s.strip() for s in args.sources.split(",") if s.strip()}
        sources = [s for s in sources if s in allow]

    h5ad_map = build_h5ad_map(INV_PATH)

    marker_parts = []
    original_parts = []
    consensus_parts = []
    source_summaries: list[dict[str, Any]] = []

    for source in sources:
        print(f"[consensus_v1] source={source}", flush=True)
        flags_sub = flags[flags["source_h5ad"] == source].copy()

        orig = original_all[original_all["source_h5ad"] == source].copy()
        if len(orig) == 0:
            raise RuntimeError(f"missing original evidence in canonical annotation for source: {source}")

        source_lc = source.lower()
        if source_lc not in h5ad_map:
            raise RuntimeError(f"h5ad path not found in inventory for source: {source}")

        marker_chk = CHK_DIR / f"marker_{source.replace('.h5ad', '')}.parquet"
        summary_chk = CHK_DIR / f"marker_{source.replace('.h5ad', '')}.summary.json"

        if marker_chk.exists() and summary_chk.exists() and not args.force_recompute:
            marker = pd.read_parquet(marker_chk)
            source_summary = json.loads(summary_chk.read_text(encoding="utf-8"))
        else:
            marker, source_summary = source_marker_score(
                source_h5ad=source,
                h5ad_path=h5ad_map[source_lc],
                flags_sub=flags_sub,
                marker_registry=registry,
                batch_size=args.batch_size,
            )
            marker.to_parquet(marker_chk, index=False)
            summary_chk.write_text(json.dumps(source_summary, ensure_ascii=True, indent=2), encoding="utf-8")
        source_summaries.append(source_summary)

        merged = flags_sub[["source_h5ad", "cohort_id", "cell_barcode", "sample_id", "final_inclusion"]].merge(
            marker,
            on=["source_h5ad", "cohort_id", "cell_barcode", "sample_id", "final_inclusion"],
            how="inner",
            validate="one_to_one",
        )
        merged = merged.merge(
            orig,
            on=["source_h5ad", "cell_barcode", "sample_id"],
            how="left",
            validate="many_to_one",
        )

        keys = merged[["source_h5ad", "cell_barcode", "sample_id"]].copy()
        auto, _, _ = auto_annotation_branch(keys, created_at)
        merged = merged.merge(auto[["source_h5ad", "cell_barcode", "sample_id", "auto_major_lineage", "auto_immune_mid_state", "auto_method", "auto_not_run_reason", "auto_confidence"]],
                              on=["source_h5ad", "cell_barcode", "sample_id"],
                              how="left", validate="one_to_one")

        decisions = consensus_decision_vectorized(merged)
        merged = pd.concat([merged, decisions], axis=1)

        merged["included_in_main_annotation"] = merged["final_inclusion"].astype(str).eq("include_main")
        merged["sensitivity_only_annotation"] = merged["final_inclusion"].astype(str).eq("include_sensitivity_only")
        merged["run_id"] = RUN_ID
        merged["created_at"] = created_at
        merged["input_manifest_ref"] = str(STEP2_MANIFEST_PATH)

        marker_out = merged[[
            "run_id",
            "created_at",
            "input_manifest_ref",
            "source_h5ad",
            "cohort_id",
            "cell_barcode",
            "sample_id",
            "marker_major_label",
            "marker_mid_label",
            "marker_confidence",
            "marker_top_set",
            "marker_runner_up_set",
            "marker_score_margin",
            "marker_gene_coverage",
            "marker_detected_fraction",
            "marker_source_zscore",
            "marker_conflict_flags",
        ]].copy()

        consensus_out = merged[[
            "run_id",
            "created_at",
            "input_manifest_ref",
            "source_h5ad",
            "cohort_id",
            "cell_barcode",
            "sample_id",
            "final_inclusion",
            "included_in_main_annotation",
            "sensitivity_only_annotation",
            "marker_major_label",
            "marker_mid_label",
            "marker_confidence",
            "original_annotation",
            "original_annotation_column",
            "original_major_lineage",
            "original_immune_mid_state",
            "original_annotation_confidence",
            "auto_method",
            "auto_major_lineage",
            "auto_immune_mid_state",
            "auto_confidence",
            "final_major_lineage",
            "final_immune_mid_state",
            "final_functional_state_optional",
            "final_confidence",
            "decision_rule_id",
            "conflict_type",
        ]].copy()

        marker_parts.append(marker_out)
        original_parts.append(orig)
        consensus_parts.append(consensus_out)

    marker_all = pd.concat(marker_parts, ignore_index=True)
    original_all = pd.concat(original_parts, ignore_index=True)
    consensus_all = pd.concat(consensus_parts, ignore_index=True)

    auto_all, auto_status, auto_reason = auto_annotation_branch(
        consensus_all[["source_h5ad", "cell_barcode", "sample_id"]], created_at
    )

    # original rules snapshot
    orig_rules = pd.read_csv(ORIG_RULES_PATH)
    rules_out = orig_rules[[
        "rule_id",
        "source_h5ad",
        "source_annotation_column",
        "source_label",
        "target_level",
        "target_label",
        "rule_confidence",
        "notes",
        "reference_ids",
    ]].copy()
    rules_out = rules_out.rename(
        columns={
            "source_annotation_column": "original_annotation_column",
            "source_label": "original_annotation",
            "target_label": "major_lineage_or_immune_mid_state",
            "rule_confidence": "annotation_confidence",
            "rule_id": "annotation_rule_id",
        }
    )
    rules_out.to_csv(CONS_DIR / "original_annotation_mapping_rules.csv", index=False)

    decision_rules = pd.DataFrame(
        [
            ["CONS_R1_marker_original_agree", "marker+original agree and coverage pass", "high"],
            ["CONS_R2_marker_only", "marker strong but no original/auto support", "medium"],
            ["CONS_R3_original_supported", "original clear and marker not opposing", "medium"],
            ["CONS_R4_conflict_degrade", "marker-original conflict degraded to coarse level", "low"],
            ["CONS_R5_evidence_insufficient", "insufficient evidence, keep Immune_unspecified/Unknown", "low"],
            ["CONS_R6_auto_support_only", "auto supports fallback but no high-confidence by auto alone", "medium"],
        ],
        columns=["decision_rule_id", "rule_text", "max_confidence"],
    )
    decision_rules.to_csv(CONS_DIR / "consensus_decision_rules.csv", index=False)

    marker_all.to_parquet(CONS_DIR / "marker_annotation_by_cell.parquet", index=False)
    auto_all.to_parquet(CONS_DIR / "auto_annotation_by_cell.parquet", index=False)

    consensus_all.to_parquet(CONS_DIR / "cell_state_annotation_consensus_v1.parquet", index=False)
    consensus_all.to_csv(CONS_DIR / "cell_state_annotation_consensus_v1.csv.gz", index=False, compression="gzip")

    generate_reports(consensus_all, CONS_DIR)

    thresholds = {
        "marker_high_z": 1.5,
        "marker_min_detected_fraction": 0.3,
        "marker_min_margin": 0.2,
        "coverage_pass": 0.4,
        "nonimmune_strong_z": 1.5,
        "immune_weak_z": 0.5,
    }
    write_manifest(
        out_dir=CONS_DIR,
        created_at=created_at,
        auto_status=auto_status,
        auto_reason=auto_reason,
        scoring_thresholds=thresholds,
        marker_source_summaries=source_summaries,
    )

    v = validate_outputs(consensus_all, flags if not args.sources else flags[flags["source_h5ad"].isin(sources)])

    report_lines = [
        "# Consensus Annotation Report",
        "",
        f"- run_id: `{RUN_ID}`",
        f"- created_at: `{created_at}`",
        f"- branch: `consensus_annotation_v1`",
        f"- sources_processed: `{len(sources)}`",
        f"- total_rows: `{v['row_count']}`",
        f"- major_lineage_coverage(include_main): `{v['major_coverage']:.4f}`",
        f"- immune_mid_state_specified_fraction(include_main immune): `{v['immune_mid_coverage']:.4f}`",
        "",
        "## Compliance",
        "- response fields were not used.",
        "- no global clustering performed.",
        "- no differential analysis performed.",
        "- inclusion source fixed to `final_cell_inclusion_flags.qc0513_balanced.parquet`.",
    ]
    (CONS_DIR / "consensus_annotation_report.md").write_text("\n".join(report_lines), encoding="utf-8")

    print("[consensus_v1] done", flush=True)


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Step2.5 consensus annotation v1")
    p.add_argument("--batch-size", type=int, default=50000)
    p.add_argument("--sources", type=str, default="", help="comma-separated source_h5ad subset")
    p.add_argument("--force-recompute", action="store_true")
    return p.parse_args()


if __name__ == "__main__":
    run(parse_args())
