#!/usr/bin/env python3
from __future__ import annotations

import argparse
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
CONS_DIR = STEP2_ROOT / "04_annotation" / "consensus_annotation_v1"
OUT_DIR = STEP2_ROOT / "05_myeloid_qc" / "myeloid_adjudication_v1"
CHK_DIR = OUT_DIR / "checkpoints"

CONSENSUS_PATH = CONS_DIR / "cell_state_annotation_consensus_v1.parquet"
STEP2_5_MARKER_PATH = CONS_DIR / "marker_annotation_by_cell.parquet"
STEP2_5_MANIFEST_PATH = CONS_DIR / "consensus_annotation_manifest.yaml"
INV_PATH = STEP2_ROOT / "00_manifest" / "input_file_inventory.csv"
GENE_STD_MANIFEST_PATH = STEP2_ROOT / "02_gene_standardization" / "gene_standardization_manifest.yaml"
INPUT_MANIFEST_REF = STEP2_5_MANIFEST_PATH

KEY_COLS = ["source_h5ad", "cell_barcode", "sample_id"]
PROVENANCE_COLS = ["run_id", "created_at", "input_manifest_ref"]
RESPONSE_FORBIDDEN_SUBSTRINGS = ["response_strict", "response_broad"]

SAMPLE_COL_MAP = {
    "gse272734.h5ad": "sample_prefix",
    "gse272735.h5ad": "sample_prefix",
    "gse273718.h5ad": "sample_prefix",
}

MYELOID_MID_PATTERN = r"Myeloid|Mono|Macro|DC|pDC|cDC|Mast"
ALLOWED_STATUSES = {
    "accept_consensus",
    "refine_consensus",
    "downgrade_to_myeloid_unspecified",
    "conflict_review",
    "exclude_from_myeloid_main",
    "sensitivity_only",
}
ALLOWED_LABELS = {
    "Mono_FCN1",
    "Mono_inflammatory",
    "Macro_C1QC",
    "Macro_SPP1_TAM_like",
    "Macro_inflammatory",
    "cDC1",
    "cDC2",
    "pDC",
    "Mast_or_minor_immune",
    "Myeloid_DC_core",
    "Myeloid_unspecified",
    "Myeloid_conflict_review",
}

PANEL_DEFS: dict[str, list[str]] = {
    "pan_myeloid": ["LST1", "TYROBP", "FCER1G", "AIF1"],
    "mono_inflammatory": ["FCN1", "S100A8", "S100A9", "IL1B", "VCAN"],
    "macro": ["C1QA", "C1QB", "C1QC", "APOE", "MRC1"],
    "tam_like": ["SPP1", "MARCO", "TREM2", "CD163", "MSR1", "LGALS3"],
    "cdc1": ["CLEC9A", "XCR1", "BATF3"],
    "cdc2": ["CD1C", "FCER1A", "CLEC10A"],
    "pdc": ["LILRA4", "GZMB", "IRF7", "TCF4"],
    "apc": ["HLA-DRA", "HLA-DRB1", "CD74", "B2M"],
    "mast_minor_immune": ["TPSAB1", "TPSB2", "CPA3", "KIT", "MS4A2"],
    "nonimmune_exclusion": ["EPCAM", "KRT8", "KRT18", "KRT19", "PECAM1", "VWF", "COL1A1", "DCN", "LUM"],
}

PANEL_TO_LABEL = {
    "mono_inflammatory": "Mono_inflammatory",
    "macro": "Macro_C1QC",
    "tam_like": "Macro_SPP1_TAM_like",
    "cdc1": "cDC1",
    "cdc2": "cDC2",
    "pdc": "pDC",
    "mast_minor_immune": "Mast_or_minor_immune",
    "pan_myeloid": "Myeloid_DC_core",
    "apc": "Myeloid_DC_core",
}

MID_TO_LABEL = {
    "Mono_FCN1": "Mono_FCN1",
    "Macro_C1QC": "Macro_C1QC",
    "Macro_SPP1": "Macro_SPP1_TAM_like",
    "Macro_inflammatory": "Macro_inflammatory",
    "cDC1": "cDC1",
    "cDC2": "cDC2",
    "pDC": "pDC",
    "Myeloid_DC_core": "Myeloid_DC_core",
    "Mast_minor_immune": "Mast_or_minor_immune",
    "Mast_or_minor_immune": "Mast_or_minor_immune",
}

THRESHOLDS = {
    "coverage_min": 0.4,
    "support_z_min": 0.5,
    "support_detected_min": 0.15,
    "strong_z_min": 1.0,
    "strong_detected_min": 0.20,
    "margin_min": 0.10,
    "nonimmune_exclusion_z_min": 1.0,
    "nonimmune_exclusion_detected_min": 0.20,
    "tam_like_z_min": 0.75,
    "tam_like_detected_min": 0.20,
}

ENSEMBL_TO_SYMBOL_PATH = Path("/home/huyudi/006/data/ref/ensembl_to_symbol_gencode_v49.csv")


def now_iso() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat()


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def build_h5ad_map(inv_path: Path) -> dict[str, str]:
    inv = pd.read_csv(inv_path)
    sub = inv[inv["file_type"] == "h5ad"].copy()
    return {Path(p).name.lower(): p for p in sub["file_path"].astype(str)}


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
    obs["sample_id"] = obs[sample_col].astype(str) if sample_col else "unknown"
    obs["_cell_idx"] = np.arange(len(obs), dtype=np.int64)
    return obs[["cell_barcode", "sample_id", "_cell_idx"]]


def build_var_lookup(var_names: Any, var: pd.DataFrame) -> dict[str, int]:
    lookup: dict[str, int] = {}
    ensembl_map: dict[str, str] = {}

    var_names_list = list(var_names) if var_names is not None else []
    has_gene_id_col = "gene_id" in var.columns
    var_names_look_ensembl = (
        any(str(v).strip().startswith("ENSG") for v in var_names_list[:10])
        if len(var_names_list) > 0
        else False
    )

    if (has_gene_id_col or var_names_look_ensembl) and ENSEMBL_TO_SYMBOL_PATH.exists():
        df_map = pd.read_csv(ENSEMBL_TO_SYMBOL_PATH, header=None, names=["ensembl_id", "gene_symbol"])
        ensembl_map = dict(
            zip(
                df_map["ensembl_id"].astype(str).str.strip(),
                df_map["gene_symbol"].astype(str).str.strip().str.upper(),
            )
        )

    def add_symbol(symbol: Any, idx: int) -> None:
        if pd.isna(symbol):
            return
        s = str(symbol).strip().upper()
        if s:
            lookup.setdefault(s, idx)

    for i, gene in enumerate(var_names_list):
        add_symbol(gene, i)
        if ensembl_map:
            symbol = ensembl_map.get(str(gene).strip())
            if symbol:
                lookup.setdefault(symbol, i)
    for col in ["gene_symbol", "geneSymbol", "gene_symbols", "features", "feature_name", "gene_name"]:
        if col in var.columns:
            for i, gene in enumerate(var[col].tolist()):
                add_symbol(gene, i)
    if "gene_id" in var.columns and ensembl_map:
        for i, gid in enumerate(var["gene_id"].tolist()):
            gid_str = str(gid).strip()
            symbol = ensembl_map.get(gid_str)
            if symbol:
                lookup.setdefault(symbol, i)
    return lookup


def resolve_expr_and_var(adata: anndata.AnnData):
    for layer in ["counts", "data"]:
        if layer in set(adata.layers.keys()):
            return adata.layers[layer], adata.var_names, adata.var
    if adata.raw is not None and adata.raw.X is not None:
        return adata.raw.X, adata.raw.var_names, adata.raw.var
    if adata.X is None:
        raise RuntimeError("no expression matrix found in X/layers/raw")
    return adata.X, adata.var_names, adata.var


def select_myeloid_candidates(consensus: pd.DataFrame) -> pd.DataFrame:
    df = consensus.copy()
    included = df["included_in_main_annotation"].astype(bool)
    major = df["final_major_lineage"].fillna("Unknown").astype(str)
    mid = df["final_immune_mid_state"].fillna("NA").astype(str)
    myeloid_mid = mid.str.contains(MYELOID_MID_PATTERN, case=False, regex=True, na=False)
    candidate = included & ((major == "Myeloid_DC") | (major == "Immune") | myeloid_mid)
    out = df.loc[candidate].copy()
    out["candidate_selection_rule"] = np.select(
        [
            out["final_major_lineage"].astype(str).eq("Myeloid_DC"),
            out["final_major_lineage"].astype(str).eq("Immune"),
            out["final_immune_mid_state"].astype(str).str.contains(MYELOID_MID_PATTERN, case=False, regex=True, na=False),
        ],
        ["primary_final_major_myeloid_dc", "review_final_major_immune", "review_myeloid_mid_state"],
        default="review_other",
    )
    return out


def panel_support(z: pd.Series, det: pd.Series, cov: pd.Series, *, strong: bool = False) -> pd.Series:
    if strong:
        return (z >= THRESHOLDS["strong_z_min"]) & (det >= THRESHOLDS["strong_detected_min"]) & (cov >= THRESHOLDS["coverage_min"])
    return (z >= THRESHOLDS["support_z_min"]) & (det >= THRESHOLDS["support_detected_min"]) & (cov >= THRESHOLDS["coverage_min"])


def specific_panel_support(df: pd.DataFrame, panel: str) -> pd.Series:
    z = df[f"{panel}_zscore"]
    det = df[f"{panel}_detected_fraction"]
    cov = df[f"{panel}_gene_coverage"]
    if panel == "tam_like":
        return (
            (z >= THRESHOLDS["tam_like_z_min"])
            & (det >= THRESHOLDS["tam_like_detected_min"])
            & (cov >= THRESHOLDS["coverage_min"])
        )
    return panel_support(z, det, cov, strong=(panel in {"cdc1", "cdc2", "pdc"}))


def adjudicate_myeloid(candidates: pd.DataFrame, scores: pd.DataFrame) -> pd.DataFrame:
    df = candidates.merge(scores, on=KEY_COLS + ["cohort_id"], how="left", validate="one_to_one")
    n = len(df)
    status = np.full(n, "downgrade_to_myeloid_unspecified", dtype=object)
    reason = np.full(n, "insufficient_specific_marker_support", dtype=object)
    label = np.full(n, "Myeloid_unspecified", dtype=object)
    main = np.ones(n, dtype=bool)
    sensitivity = np.zeros(n, dtype=bool)

    final_mid = df["final_immune_mid_state"].fillna("Immune_unspecified").astype(str)
    final_major = df["final_major_lineage"].fillna("Unknown").astype(str)
    final_conf = df["final_confidence"].fillna("low").astype(str)
    conflict = df["conflict_type"].fillna("").astype(str)
    decision = df["decision_rule_id"].fillna("").astype(str)
    top_panel = df["myeloid_top_panel"].fillna("none").astype(str)

    nonimmune_support = (
        (df["nonimmune_exclusion_zscore"] >= THRESHOLDS["nonimmune_exclusion_z_min"])
        & (df["nonimmune_exclusion_detected_fraction"] >= THRESHOLDS["nonimmune_exclusion_detected_min"])
    )
    conflict_mask = conflict.isin(["immune_major_conflict", "immune_nonimmune_conflict"]) | decision.eq("CONS_R4_conflict_degrade")
    low_immune_review = final_major.eq("Immune") & final_conf.eq("low")

    status[nonimmune_support.to_numpy()] = "exclude_from_myeloid_main"
    reason[nonimmune_support.to_numpy()] = "nonimmune_exclusion_marker_support"
    label[nonimmune_support.to_numpy()] = "Myeloid_conflict_review"
    main[nonimmune_support.to_numpy()] = False
    sensitivity[nonimmune_support.to_numpy()] = True

    sens_mask = (conflict_mask | low_immune_review).to_numpy() & (~nonimmune_support.to_numpy())
    status[sens_mask] = "sensitivity_only"
    reason[sens_mask] = "conflict_or_low_confidence_immune_review"
    label[sens_mask] = "Myeloid_conflict_review"
    main[sens_mask] = False
    sensitivity[sens_mask] = True

    specific_panels = ["mono_inflammatory", "macro", "tam_like", "cdc1", "cdc2", "pdc", "mast_minor_immune"]
    panel_support_mask = np.zeros((n, len(specific_panels)), dtype=bool)

    for j, panel in enumerate(specific_panels):
        support = specific_panel_support(df, panel)
        panel_support_mask[:, j] = support.to_numpy() & (~nonimmune_support.to_numpy()) & (~sens_mask)

    # Tie-break multi-panel support by highest score (avoids fixed-order bias)
    score_cols = [f"{p}_score" for p in specific_panels]
    score_mat = df[score_cols].to_numpy().copy()
    score_mat[~panel_support_mask] = -np.inf
    score_mat = np.where(np.isnan(score_mat), -np.inf, score_mat)
    has_panel_support = panel_support_mask.any(axis=1)
    best_panel_idx = np.full(n, -1, dtype=np.int64)
    if has_panel_support.any():
        best_panel_idx[has_panel_support] = np.argmax(score_mat[has_panel_support], axis=1)

    for j, panel in enumerate(specific_panels):
        mask = has_panel_support & (best_panel_idx == j)
        if not mask.any():
            continue
        label[mask] = PANEL_TO_LABEL[panel]
        status[mask] = np.where(final_mid.map(MID_TO_LABEL).fillna("").to_numpy()[mask] == label[mask], "accept_consensus", "refine_consensus")
        reason[mask] = f"{panel}_marker_support"

    mono_fcn1 = final_mid.eq("Mono_FCN1") & specific_panel_support(df, "mono_inflammatory")
    mask = mono_fcn1.to_numpy() & np.isin(label, ["Myeloid_unspecified", "Mono_inflammatory"]) & (~nonimmune_support.to_numpy()) & (~sens_mask)
    label[mask] = "Mono_FCN1"
    status[mask] = "accept_consensus"
    reason[mask] = "consensus_mono_fcn1_with_marker_support"

    coarse_support = (
        panel_support(df["pan_myeloid_zscore"], df["pan_myeloid_detected_fraction"], df["pan_myeloid_gene_coverage"])
        | panel_support(df["apc_zscore"], df["apc_detected_fraction"], df["apc_gene_coverage"])
    )
    coarse_mask = coarse_support.to_numpy() & (label == "Myeloid_unspecified") & (~nonimmune_support.to_numpy()) & (~sens_mask)
    label[coarse_mask] = "Myeloid_DC_core"
    status[coarse_mask] = "accept_consensus"
    reason[coarse_mask] = "coarse_myeloid_or_apc_marker_support"

    no_support_mask = (label == "Myeloid_unspecified") & (~nonimmune_support.to_numpy()) & (~sens_mask)
    status[no_support_mask] = "downgrade_to_myeloid_unspecified"
    reason[no_support_mask] = "no_specific_panel_support"

    out = df.copy()
    out["myeloid_adjudication_status"] = status
    out["myeloid_adjudication_reason"] = reason
    out["final_myeloid_label"] = label
    out["included_in_myeloid_main"] = main
    out["myeloid_sensitivity_only"] = sensitivity
    return out


def score_source_panels(source_h5ad: str, h5ad_path: str, cand_sub: pd.DataFrame, batch_size: int) -> tuple[pd.DataFrame, dict[str, Any]]:
    adata = anndata.read_h5ad(h5ad_path, backed="r")
    obs_keys = prepare_source_obs_keys(adata, source_h5ad)
    keys = cand_sub[["cell_barcode", "sample_id"]].copy()
    keys["cell_barcode"] = keys["cell_barcode"].astype(str)
    keys["sample_id"] = keys["sample_id"].astype(str)
    merged = keys.merge(obs_keys, on=["cell_barcode", "sample_id"], how="left")
    if merged["_cell_idx"].isna().any() and not obs_keys["cell_barcode"].duplicated().any():
        fallback = keys[["cell_barcode"]].merge(obs_keys[["cell_barcode", "_cell_idx"]], on="cell_barcode", how="left")
        merged["_cell_idx"] = fallback["_cell_idx"].values
    if merged["_cell_idx"].isna().any():
        missing = int(merged["_cell_idx"].isna().sum())
        adata.file.close()
        raise RuntimeError(f"{source_h5ad}: failed key match for {missing} candidate cells")
    cell_idx = merged["_cell_idx"].astype(np.int64).to_numpy()

    X, var_names, var_frame = resolve_expr_and_var(adata)
    var_lookup = build_var_lookup(var_names, var_frame)
    panel_defs = {}
    for panel, genes in PANEL_DEFS.items():
        idx = [var_lookup[g] for g in genes if g in var_lookup]
        panel_defs[panel] = {"genes": genes, "idx": np.array(idx, dtype=np.int64), "n_expected": len(genes)}

    union_cols = sorted({int(i) for d in panel_defs.values() for i in d["idx"]})
    union_pos = {c: i for i, c in enumerate(union_cols)}
    for panel in panel_defs:
        panel_defs[panel]["union_pos"] = np.array([union_pos[int(c)] for c in panel_defs[panel]["idx"]], dtype=np.int64)

    obs_full = adata.obs.copy()
    lib_col = ""
    for c in ["nCount_RNA", "n_counts", "total_counts", "nCount", "nUMI", "nCountRNA", "n_counts_all"]:
        if c in obs_full.columns:
            lib_col = c
            break
    lib_all = pd.to_numeric(obs_full[lib_col], errors="coerce").fillna(0).to_numpy(dtype=np.float32) if lib_col else None

    panel_names = list(PANEL_DEFS.keys())
    n_cells = len(cell_idx)
    n_panels = len(panel_names)
    score_mat = np.full((n_cells, n_panels), np.nan, dtype=np.float32)
    detect_mat = np.full((n_cells, n_panels), np.nan, dtype=np.float32)
    cov_mat = np.full((n_cells, n_panels), np.nan, dtype=np.float32)
    for j, panel in enumerate(panel_names):
        n_expected = panel_defs[panel]["n_expected"]
        n_avail = len(panel_defs[panel]["idx"])
        cov_mat[:, j] = (n_avail / n_expected) if n_expected else 0.0

    for start in range(0, n_cells, batch_size):
        end = min(start + batch_size, n_cells)
        ridx = cell_idx[start:end]
        sort_ord = np.argsort(ridx, kind="mergesort")
        ridx_sorted = ridx[sort_ord]
        inv_ord = np.empty_like(sort_ord)
        inv_ord[sort_ord] = np.arange(len(sort_ord))
        chunk = X[ridx_sorted, :][:, union_cols] if union_cols else X[ridx_sorted, :][:, []]
        if sparse.issparse(chunk):
            chunk = chunk[inv_ord, :].tocsr()
        else:
            chunk = sparse.csr_matrix(np.asarray(chunk)[inv_ord, :])
        if lib_all is not None:
            lib = lib_all[ridx].astype(np.float32, copy=False)
        else:
            lib = np.asarray(chunk.sum(axis=1)).ravel().astype(np.float32)
        lib[lib <= 0] = 1.0
        scale = (10000.0 / lib).astype(np.float32)

        for j, panel in enumerate(panel_names):
            cols = panel_defs[panel]["union_pos"]
            if cols.size == 0:
                continue
            sub = chunk[:, cols]
            detect_mat[start:end, j] = np.asarray((sub > 0).sum(axis=1)).ravel().astype(np.float32) / float(cols.size)
            scaled = sub.multiply(scale[:, None]).tocsr(copy=True)
            scaled.data = np.log1p(scaled.data)
            score_mat[start:end, j] = np.asarray(scaled.mean(axis=1)).ravel().astype(np.float32)

        if start == 0 or ((start // batch_size) % 20 == 0):
            print(f"[step2.6] {source_h5ad} batch {start}:{end}/{n_cells}", flush=True)

    mean = np.nanmean(score_mat, axis=0)
    std = np.nanstd(score_mat, axis=0)
    std[std == 0] = 1.0
    z_mat = (score_mat - mean[None, :]) / std[None, :]
    top_idx = np.nanargmax(np.where(np.isnan(score_mat), -np.inf, score_mat), axis=1)
    sorted_idx = np.argsort(np.where(np.isnan(score_mat), -np.inf, score_mat), axis=1)
    runner_idx = sorted_idx[:, -2] if n_panels >= 2 else top_idx

    out = cand_sub[["source_h5ad", "cohort_id", "cell_barcode", "sample_id"]].copy()
    for j, panel in enumerate(panel_names):
        out[f"{panel}_score"] = score_mat[:, j].astype(np.float32)
        out[f"{panel}_zscore"] = z_mat[:, j].astype(np.float32)
        out[f"{panel}_detected_fraction"] = detect_mat[:, j].astype(np.float32)
        out[f"{panel}_gene_coverage"] = cov_mat[:, j].astype(np.float32)
    out["myeloid_top_panel"] = np.array([panel_names[i] for i in top_idx], dtype=object)
    out["myeloid_runner_up_panel"] = np.array([panel_names[i] for i in runner_idx], dtype=object)
    out["myeloid_score_margin"] = (score_mat[np.arange(n_cells), top_idx] - score_mat[np.arange(n_cells), runner_idx]).astype(np.float32)
    out["myeloid_ambiguity_flag"] = np.where(out["myeloid_score_margin"] < THRESHOLDS["margin_min"], "low_margin_ambiguous", "")

    summary = {
        "source_h5ad": source_h5ad,
        "n_cells": int(n_cells),
        "lib_size_column": lib_col or "marker_union_sum_fallback",
        "panel_gene_coverage": {p: float(cov_mat[0, j]) if n_cells else 0.0 for j, p in enumerate(panel_names)},
        "panel_mean": {p: float(mean[j]) for j, p in enumerate(panel_names)},
        "panel_std": {p: float(std[j]) for j, p in enumerate(panel_names)},
    }
    adata.file.close()
    return out, summary


def add_provenance(df: pd.DataFrame, created_at: str) -> pd.DataFrame:
    out = df.copy()
    out.insert(0, "input_manifest_ref", str(INPUT_MANIFEST_REF))
    out.insert(0, "created_at", created_at)
    out.insert(0, "run_id", RUN_ID)
    return out


def write_parquet_stream(path: Path, frames: list[pd.DataFrame]) -> None:
    writer = None
    try:
        for frame in frames:
            table = pa.Table.from_pandas(frame, preserve_index=False)
            if writer is None:
                writer = pq.ParquetWriter(path, table.schema, compression="zstd")
            writer.write_table(table)
    finally:
        if writer is not None:
            writer.close()


def validate_outputs(candidate: pd.DataFrame, scores: pd.DataFrame, adjud: pd.DataFrame) -> dict[str, Any]:
    forbidden_cols = [c for c in adjud.columns if any(x in c.lower() for x in RESPONSE_FORBIDDEN_SUBSTRINGS)]
    valid = {
        "candidate_rows": int(len(candidate)),
        "score_rows": int(len(scores)),
        "adjudication_rows": int(len(adjud)),
        "candidate_expected_rows": 1522971,
        "rows_match_candidate_score_adjudication": len(candidate) == len(scores) == len(adjud),
        "candidate_matches_locked_plan_count": len(candidate) == 1522971,
        "unique_keys": not adjud.duplicated(KEY_COLS).any(),
        "missing_key_any": bool(adjud[KEY_COLS].isna().any().any()),
        "unknown_sample_id": int((adjud["sample_id"].astype(str) == "unknown").sum()),
        "all_status_allowed": set(adjud["myeloid_adjudication_status"].astype(str).unique()).issubset(ALLOWED_STATUSES),
        "all_labels_allowed": set(adjud["final_myeloid_label"].astype(str).unique()).issubset(ALLOWED_LABELS),
        "forbidden_response_cols": forbidden_cols,
    }
    valid["passed"] = (
        valid["rows_match_candidate_score_adjudication"]
        and valid["candidate_matches_locked_plan_count"]
        and valid["unique_keys"]
        and not valid["missing_key_any"]
        and valid["unknown_sample_id"] == 0
        and valid["all_status_allowed"]
        and valid["all_labels_allowed"]
        and len(valid["forbidden_response_cols"]) == 0
    )
    return valid


def write_summaries(adjud: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    main = adjud[adjud["included_in_myeloid_main"].astype(bool)].copy()
    sample = (
        adjud.groupby(["cohort_id", "source_h5ad", "sample_id", "final_myeloid_label", "myeloid_adjudication_status"], dropna=False)
        .size()
        .reset_index(name="n_cells")
    )
    total = adjud.groupby(["cohort_id", "source_h5ad", "sample_id"], dropna=False).size().reset_index(name="candidate_cells")
    sample = sample.merge(total, on=["cohort_id", "source_h5ad", "sample_id"], how="left")
    sample["fraction_of_candidates"] = sample["n_cells"] / sample["candidate_cells"].replace(0, np.nan)

    cohort = (
        adjud.groupby(["cohort_id", "final_myeloid_label", "myeloid_adjudication_status"], dropna=False)
        .size()
        .reset_index(name="n_cells")
    )
    ctotal = adjud.groupby(["cohort_id"], dropna=False).size().reset_index(name="candidate_cells")
    cohort = cohort.merge(ctotal, on="cohort_id", how="left")
    cohort["fraction_of_candidates"] = cohort["n_cells"] / cohort["candidate_cells"].replace(0, np.nan)

    conflict = (
        adjud.groupby(["cohort_id", "source_h5ad", "conflict_type", "myeloid_adjudication_status", "final_myeloid_label"], dropna=False)
        .size()
        .reset_index(name="n_cells")
        .sort_values(["n_cells"], ascending=False)
    )
    return sample, cohort, conflict


def write_report(path: Path, validation: dict[str, Any], adjud: pd.DataFrame, source_summaries: list[dict[str, Any]], created_at: str) -> None:
    label_counts = adjud["final_myeloid_label"].value_counts(dropna=False).to_dict()
    status_counts = adjud["myeloid_adjudication_status"].value_counts(dropna=False).to_dict()

    zero_cov: list[str] = []
    low_cov: list[str] = []
    for s in source_summaries:
        covs = s.get("panel_gene_coverage", {})
        min_cov = min(covs.values()) if covs else 0.0
        if min_cov == 0.0:
            zero_cov.append(f"- `{s['source_h5ad']}`: {s['n_cells']} cells, all panels zero coverage")
        elif min_cov < 1.0:
            low_panels = ", ".join([f"{p}={v:.2f}" for p, v in covs.items() if v < 1.0])
            low_cov.append(f"- `{s['source_h5ad']}`: {s['n_cells']} cells, low panels: {low_panels}")

    lines = [
        "# Step2.6 Myeloid Focused Adjudication QC Report",
        "",
        f"- run_id: `{RUN_ID}`",
        f"- created_at: `{created_at}`",
        "- branch: `myeloid_adjudication_v1`",
        f"- candidate_rows: `{validation['candidate_rows']}`",
        f"- sources_processed: `{len(source_summaries)}`",
        f"- validation_passed: `{validation['passed']}`",
        "",
        "## Compliance",
        "- response fields were not used.",
        "- no global clustering performed.",
        "- no differential analysis performed.",
        "- input fixed to Step2.5 `consensus_annotation_v1`.",
        "",
        "## Status Counts",
        "",
    ]
    lines += [f"- `{k}`: `{v}`" for k, v in status_counts.items()]
    lines += ["", "## Label Counts", ""]
    lines += [f"- `{k}`: `{v}`" for k, v in label_counts.items()]
    lines += ["", "## Coverage Warnings", ""]
    if zero_cov:
        lines += ["### Zero Coverage Sources (marker genes unresolved)", ""]
        lines += zero_cov
    if low_cov:
        lines += ["### Low Coverage Sources (incomplete marker panels)", ""]
        lines += low_cov
    if not zero_cov and not low_cov:
        lines += ["- No coverage issues detected."]
    lines += ["", "## Validation", ""]
    lines += [f"- `{k}`: `{v}`" for k, v in validation.items()]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--batch-size", type=int, default=12000)
    parser.add_argument("--sources", type=str, default="", help="comma-separated source_h5ad subset for debugging")
    args = parser.parse_args()

    created_at = now_iso()
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    CHK_DIR.mkdir(parents=True, exist_ok=True)

    consensus_cols = [
        "source_h5ad", "cohort_id", "cell_barcode", "sample_id", "final_inclusion", "included_in_main_annotation",
        "sensitivity_only_annotation", "marker_major_label", "marker_mid_label", "marker_confidence",
        "original_annotation", "original_annotation_column", "original_major_lineage", "original_immune_mid_state",
        "original_annotation_confidence", "auto_method", "auto_major_lineage", "auto_immune_mid_state", "auto_confidence",
        "final_major_lineage", "final_immune_mid_state", "final_functional_state_optional", "final_confidence", "decision_rule_id", "conflict_type",
    ]
    consensus = pd.read_parquet(CONSENSUS_PATH, columns=consensus_cols)
    if args.sources:
        wanted = {s.strip() for s in args.sources.split(",") if s.strip()}
        consensus = consensus[consensus["source_h5ad"].isin(wanted)].copy()

    candidates = select_myeloid_candidates(consensus)
    candidates = add_provenance(candidates, created_at)
    candidates.to_parquet(OUT_DIR / "myeloid_candidate_cells.parquet", index=False)

    h5ad_map = build_h5ad_map(INV_PATH)
    score_frames: list[pd.DataFrame] = []
    source_summaries: list[dict[str, Any]] = []
    for source in sorted(candidates["source_h5ad"].astype(str).unique()):
        source_lc = source.lower()
        if source_lc not in h5ad_map:
            raise RuntimeError(f"h5ad path not found in inventory for source: {source}")
        sub = candidates[candidates["source_h5ad"] == source].copy()
        score_path = CHK_DIR / f"myeloid_scores_{source.replace('.h5ad', '')}.parquet"
        summary_path = CHK_DIR / f"myeloid_scores_{source.replace('.h5ad', '')}.summary.json"
        if score_path.exists() and summary_path.exists():
            scores = pd.read_parquet(score_path)
            summary = json.loads(summary_path.read_text(encoding="utf-8"))
            for col in PROVENANCE_COLS:
                if col in scores.columns:
                    scores = scores.drop(columns=[col])
            scores = add_provenance(scores, created_at)
        else:
            scores, summary = score_source_panels(source, h5ad_map[source_lc], sub, args.batch_size)
            scores = add_provenance(scores, created_at)
            scores.to_parquet(score_path, index=False)
            summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True), encoding="utf-8")
        score_frames.append(scores)
        source_summaries.append(summary)

    scores_all = pd.concat(score_frames, ignore_index=True) if score_frames else pd.DataFrame()
    scores_all.to_parquet(OUT_DIR / "myeloid_marker_scores_by_cell.parquet", index=False)

    # Remove duplicated provenance from scores before merge; candidates already carries provenance.
    score_payload_cols = [c for c in scores_all.columns if c not in PROVENANCE_COLS]
    adjud = adjudicate_myeloid(candidates, scores_all[score_payload_cols])
    adjud.to_parquet(OUT_DIR / "myeloid_annotation_adjudication.parquet", index=False)
    adjud.to_csv(OUT_DIR / "myeloid_annotation_adjudication.csv.gz", index=False, compression="gzip")

    sample, cohort, conflict = write_summaries(adjud)
    sample.to_csv(OUT_DIR / "myeloid_fraction_by_sample.csv", index=False)
    cohort.to_csv(OUT_DIR / "myeloid_fraction_by_cohort.csv", index=False)
    conflict.to_csv(OUT_DIR / "myeloid_conflict_audit.csv", index=False)

    validation = validate_outputs(candidates, scores_all, adjud)
    zero_cov_sources = []
    low_cov_sources = []
    for s in source_summaries:
        covs = s.get("panel_gene_coverage", {})
        min_cov = min(covs.values()) if covs else 0.0
        if min_cov == 0.0:
            zero_cov_sources.append(s["source_h5ad"])
        elif min_cov < 1.0:
            low_cov_sources.append({"source": s["source_h5ad"], "min_coverage": round(float(min_cov), 2)})
    validation["zero_coverage_sources"] = zero_cov_sources
    validation["low_coverage_sources"] = low_cov_sources
    validation["has_coverage_issues"] = bool(zero_cov_sources or low_cov_sources)
    manifest = {
        "run_id": RUN_ID,
        "created_at": created_at,
        "branch": "myeloid_adjudication_v1",
        "input_manifest_ref": str(INPUT_MANIFEST_REF),
        "input_paths": {
            "consensus_annotation": str(CONSENSUS_PATH),
            "step2_5_marker_annotation": str(STEP2_5_MARKER_PATH),
            "step2_5_manifest": str(STEP2_5_MANIFEST_PATH),
            "h5ad_inventory": str(INV_PATH),
            "gene_standardization_manifest": str(GENE_STD_MANIFEST_PATH),
        },
        "input_hashes": {p.name: sha256_file(p) for p in [CONSENSUS_PATH, STEP2_5_MARKER_PATH, STEP2_5_MANIFEST_PATH, INV_PATH, GENE_STD_MANIFEST_PATH]},
        "panel_definitions": PANEL_DEFS,
        "scoring_thresholds": THRESHOLDS,
        "candidate_selection": {
            "included_in_main_annotation": True,
            "final_major_lineage": ["Myeloid_DC", "Immune"],
            "myeloid_mid_state_pattern": MYELOID_MID_PATTERN,
        },
        "validation": validation,
        "source_summaries": source_summaries,
        "output_files": {
            "candidate_cells": str(OUT_DIR / "myeloid_candidate_cells.parquet"),
            "marker_scores": str(OUT_DIR / "myeloid_marker_scores_by_cell.parquet"),
            "adjudication_parquet": str(OUT_DIR / "myeloid_annotation_adjudication.parquet"),
            "adjudication_csv_gz": str(OUT_DIR / "myeloid_annotation_adjudication.csv.gz"),
            "fraction_by_sample": str(OUT_DIR / "myeloid_fraction_by_sample.csv"),
            "fraction_by_cohort": str(OUT_DIR / "myeloid_fraction_by_cohort.csv"),
            "conflict_audit": str(OUT_DIR / "myeloid_conflict_audit.csv"),
            "qc_report": str(OUT_DIR / "myeloid_qc_report.md"),
        },
        "response_used": False,
        "global_clustering_used": False,
        "differential_analysis_used": False,
    }
    with (OUT_DIR / "myeloid_decision_manifest.yaml").open("w", encoding="utf-8") as f:
        yaml.safe_dump(manifest, f, sort_keys=False, allow_unicode=True)
    write_report(OUT_DIR / "myeloid_qc_report.md", validation, adjud, source_summaries, created_at)

    if not validation["passed"] and not args.sources:
        raise RuntimeError(f"Step2.6 validation failed: {validation}")
    print(json.dumps(validation, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
