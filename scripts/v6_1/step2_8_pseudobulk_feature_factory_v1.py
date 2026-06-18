#!/usr/bin/env python3
from __future__ import annotations

import argparse
import hashlib
import json
import math
import re
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable

import anndata as ad
import h5py
import numpy as np
import pandas as pd
import yaml
from scipy import sparse

RUN_ID = "step2_v6_1_0505_0319"
REPO_ROOT = Path("/home/huyudi/006")
STEP2_ROOT = REPO_ROOT / "results" / "v6_1" / "step2" / RUN_ID
OUT_DIR = STEP2_ROOT / "07_pseudobulk" / "pseudobulk_feature_v1"
CONS_DIR = STEP2_ROOT / "04_annotation" / "consensus_annotation_v1"
MYELOID_DIR = STEP2_ROOT / "05_myeloid_qc" / "myeloid_adjudication_v1"
FRACTION_DIR = STEP2_ROOT / "06_fraction" / "fraction_feature_v1"
GENE_DIR = STEP2_ROOT / "02_gene_standardization"
MANIFEST_DIR = STEP2_ROOT / "00_manifest"

CONSENSUS_PATH = CONS_DIR / "cell_state_annotation_consensus_v1.parquet"
MYELOID_PATH = MYELOID_DIR / "myeloid_annotation_adjudication.parquet"
FRACTION_MANIFEST_PATH = FRACTION_DIR / "fraction_decision_manifest.yaml"
FLAGS_PATH = STEP2_ROOT / "03_qc" / "final_cell_inclusion_flags.qc0513_balanced.parquet"
INPUT_INVENTORY_PATH = MANIFEST_DIR / "input_file_inventory.csv"
GENE_MAPPING_PATH = GENE_DIR / "gene_id_mapping_by_cohort.csv"
CORE_GENES_PATH = GENE_DIR / "core_intersection_genes.txt"
IMMUNE_GENES_PATH = GENE_DIR / "immune_feature_genes.txt"
GENE_MANIFEST_PATH = GENE_DIR / "gene_standardization_manifest.yaml"
INPUT_MANIFEST_REF = FRACTION_MANIFEST_PATH

CELL_KEY_COLS = ["source_h5ad", "cohort_id", "cell_barcode", "sample_id"]
KEY_COLS = ["source_h5ad", "cell_barcode", "sample_id"]
ROW_META_COLS = [
    "run_id",
    "created_at",
    "input_manifest_ref",
    "row_id",
    "cohort_id",
    "sample_id",
    "source_h5ad",
    "aggregation_scope",
    "cell_state_label",
    "parent_lineage",
    "cell_count",
    "library_size",
    "detected_genes",
    "gene_universe",
    "raw_count_source",
    "raw_count_status",
]
FORBIDDEN_RESPONSE_TERMS = ["response_strict", "response_broad", "response_binary", "response_ordered", "response_raw", "response_"]
IMMUNE_MAJORS = {"T_NK", "B_Plasma", "Myeloid_DC", "Mast_or_minor_immune", "Immune"}
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
MISSING_LABELS = {"", "NA", "nan", "None", "Unknown"}
MIN_STATE_CELLS = 30
RAW_STATUSES_FOR_MAIN = {"confirmed_raw_counts", "raw_like_counts"}
MAX_DENSE_MAIN_ENTRIES = 500_000_000


def now_iso() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat()


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def sanitize_label(value: Any) -> str:
    text = "Unknown" if pd.isna(value) else str(value).strip()
    if text in MISSING_LABELS:
        text = "Unknown"
    text = re.sub(r"[^0-9A-Za-z]+", "_", text).strip("_")
    return text or "Unknown"


def read_gene_list(path: Path) -> list[str]:
    return [line.strip().upper() for line in path.read_text().splitlines() if line.strip()]


def matrix_to_csr(x: Any) -> sparse.csr_matrix:
    if sparse.issparse(x):
        return x.tocsr()
    return sparse.csr_matrix(np.asarray(x))


def get_matrix_slice(matrix: Any, obs_idx: np.ndarray | list[int], var_idx: np.ndarray | list[int] | slice) -> sparse.csr_matrix:
    # scipy sparse treats matrix[[rows], [cols]] as paired fancy indexing; slice in two steps
    # to preserve the rectangular row x gene block used by pseudobulk aggregation.
    if sparse.issparse(matrix):
        x = matrix[obs_idx, :][:, var_idx]
    else:
        try:
            x = matrix[obs_idx, var_idx]
            if len(getattr(x, "shape", ())) == 1 and not isinstance(var_idx, slice):
                x = matrix[obs_idx, :][:, var_idx]
        except Exception:
            x = matrix[obs_idx, :][:, var_idx]
    if hasattr(x, "to_memory"):
        x = x.to_memory()
    return matrix_to_csr(x)


def classify_raw_values(values: np.ndarray, source_name: str) -> dict[str, Any]:
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    if values.size == 0:
        return {
            "raw_count_status": "unavailable",
            "raw_count_accept_main": False,
            "raw_count_reason": "no finite values sampled",
            "sample_nonzero_n": 0,
            "sample_min": np.nan,
            "sample_max": np.nan,
            "sample_p95": np.nan,
            "integer_like_fraction": np.nan,
        }
    nonzero = values[values != 0]
    inspect = nonzero if nonzero.size else values
    min_v = float(np.min(inspect))
    max_v = float(np.max(inspect))
    p95 = float(np.percentile(inspect, 95))
    if min_v < -1e-8:
        status = "unclear_or_normalized"
        accept = False
        reason = "negative sampled values"
    else:
        integer_like = np.isclose(inspect, np.round(inspect), atol=1e-6)
        integer_like_fraction = float(integer_like.mean()) if inspect.size else np.nan
        if source_name == "layers:counts" and integer_like_fraction >= 0.95:
            status = "confirmed_raw_counts"
            accept = True
            reason = "counts layer present, nonnegative, and integer-like"
        elif source_name == "layers:counts" and max_v > 50 and p95 > 5:
            status = "raw_like_counts"
            accept = True
            reason = "counts layer present and large nonnegative raw-like values"
        elif source_name == "layers:counts":
            status = "unclear_or_normalized"
            accept = False
            reason = "counts layer present but sampled values are not integer-like raw counts"
        elif integer_like_fraction >= 0.95:
            status = "raw_like_counts"
            accept = True
            reason = "sampled values are integer-like and nonnegative"
        elif max_v > 50 and p95 > 5:
            status = "raw_like_counts"
            accept = True
            reason = "sampled values are large nonnegative raw-like counts"
        else:
            status = "unclear_or_normalized"
            accept = False
            reason = "sampled values look normalized or log transformed"
    return {
        "raw_count_status": status,
        "raw_count_accept_main": accept,
        "raw_count_reason": reason,
        "sample_nonzero_n": int(nonzero.size),
        "sample_min": min_v,
        "sample_max": max_v,
        "sample_p95": p95,
        "integer_like_fraction": float(np.isclose(inspect, np.round(inspect), atol=1e-6).mean()) if inspect.size else np.nan,
    }


@dataclass(frozen=True)
class RawMatrixChoice:
    matrix: Any
    raw_count_source: str
    raw_count_status: str
    raw_count_accept_main: bool


def choose_raw_matrix(adata: ad.AnnData, source_h5ad: str, audit: dict[str, Any]) -> RawMatrixChoice | None:
    if "counts" in adata.layers.keys():
        return RawMatrixChoice(adata.layers["counts"], "layers:counts", audit["raw_count_status"], bool(audit["raw_count_accept_main"]))
    if adata.raw is not None:
        return RawMatrixChoice(adata.raw.X, "raw:X", audit["raw_count_status"], bool(audit["raw_count_accept_main"]))
    return RawMatrixChoice(adata.X, "X", audit["raw_count_status"], bool(audit["raw_count_accept_main"]))


def h5_matrix_shape(group: Any) -> tuple[int, int]:
    if isinstance(group, h5py.Dataset):
        return tuple(group.shape)  # type: ignore[return-value]
    if "shape" in group.attrs:
        shape = tuple(int(x) for x in group.attrs["shape"])
        if len(shape) == 2:
            return shape  # type: ignore[return-value]
    if "data" in group and "indptr" in group:
        n_major = len(group["indptr"]) - 1
        return (int(n_major), 0)
    return (0, 0)


def h5_sample_values(group: Any, max_values: int = 20000) -> np.ndarray:
    if isinstance(group, h5py.Dataset):
        if len(group.shape) == 2:
            r = min(group.shape[0], 64)
            c = min(group.shape[1], 256)
            return np.asarray(group[:r, :c]).ravel()
        return np.asarray(group[: min(group.shape[0], max_values)]).ravel()
    if "data" in group:
        data = group["data"]
        return np.asarray(data[: min(data.shape[0], max_values)]).ravel()
    return np.array([])


def audit_one_h5ad(path: Path, source_h5ad: str, n_obs_in_universe: int = 0) -> dict[str, Any]:
    if not path.exists():
        return {
            "source_h5ad": source_h5ad,
            "h5ad_path": str(path),
            "file_exists": False,
            "n_obs": 0,
            "n_vars": 0,
            "matrix_storage_kind": "missing",
            "matrix_entries": 0,
            "n_obs_in_main_universe": int(n_obs_in_universe),
            "layers": "",
            "has_raw": False,
            "raw_count_source": "unavailable",
            "raw_count_status": "unavailable",
            "raw_count_accept_main": False,
            "raw_count_reason": "h5ad file missing",
        }
    with h5py.File(path, "r") as h5:
        layers = list(h5.get("layers", {}).keys()) if "layers" in h5 else []
        has_raw = "raw" in h5 and "X" in h5["raw"]
        if "counts" in layers:
            source_name = "layers:counts"
            matrix_group = h5["layers"]["counts"]
        elif has_raw:
            source_name = "raw:X"
            matrix_group = h5["raw"]["X"]
        elif "X" in h5:
            source_name = "X"
            matrix_group = h5["X"]
        else:
            source_name = "unavailable"
            matrix_group = None
        if matrix_group is None:
            shape = (0, 0)
            matrix_storage_kind = "missing"
            stats = classify_raw_values(np.array([]), source_name)
        else:
            shape = h5_matrix_shape(matrix_group)
            matrix_storage_kind = "dense_array" if isinstance(matrix_group, h5py.Dataset) else str(matrix_group.attrs.get("encoding-type", "group"))
            values = h5_sample_values(matrix_group)
            stats = classify_raw_values(values, source_name)
        matrix_entries = int(shape[0]) * int(shape[1])
        if matrix_storage_kind == "dense_array" and matrix_entries > MAX_DENSE_MAIN_ENTRIES and stats.get("raw_count_accept_main", False):
            stats["raw_count_status"] = "raw_like_dense_deferred"
            stats["raw_count_accept_main"] = False
            stats["raw_count_reason"] = f"dense raw-like matrix exceeds main pseudobulk scalability threshold ({matrix_entries} entries)"
        return {
            "source_h5ad": source_h5ad,
            "h5ad_path": str(path),
            "file_exists": True,
            "n_obs": int(shape[0]),
            "n_vars": int(shape[1]),
            "matrix_storage_kind": matrix_storage_kind,
            "matrix_entries": matrix_entries,
            "n_obs_in_main_universe": int(n_obs_in_universe),
            "layers": ";".join(layers),
            "has_raw": bool(has_raw),
            "raw_count_source": source_name,
            **stats,
        }

def validate_consensus_against_flags(consensus_all: pd.DataFrame) -> dict[str, Any]:
    flags = pd.read_parquet(FLAGS_PATH, columns=CELL_KEY_COLS + ["final_inclusion"])
    flags_main = flags[flags["final_inclusion"].astype(str).eq("include_main")][CELL_KEY_COLS].copy()
    consensus_main = consensus_all[consensus_all["included_in_main_annotation"].astype(bool)][CELL_KEY_COLS].copy()
    if len(flags_main) != len(consensus_main):
        raise RuntimeError(f"include_main row count mismatch: flags={len(flags_main)} consensus={len(consensus_main)}")
    if flags_main.duplicated(CELL_KEY_COLS).any() or consensus_main.duplicated(CELL_KEY_COLS).any():
        raise RuntimeError("include_main keys are not unique")
    merged = consensus_main.merge(flags_main, on=CELL_KEY_COLS, how="outer", indicator=True)
    mismatches = int((merged["_merge"] != "both").sum())
    if mismatches:
        raise RuntimeError(f"include_main key mismatch between consensus and flags: {mismatches}")
    return {
        "flags_path": str(FLAGS_PATH),
        "flags_include_main_rows": int(len(flags_main)),
        "consensus_include_main_rows": int(len(consensus_main)),
        "include_main_keys_match": True,
    }




def load_step2_7_flags_check() -> dict[str, Any]:
    with FRACTION_MANIFEST_PATH.open("r", encoding="utf-8") as f:
        manifest = yaml.safe_load(f) or {}
    check = manifest.get("flags_consensus_crosscheck", {}) if isinstance(manifest, dict) else {}
    if not check.get("include_main_keys_match", False):
        raise RuntimeError("Step2.7 manifest does not confirm flags/consensus include_main key match")
    return {**check, "source": str(FRACTION_MANIFEST_PATH), "reused_from_step2_7": True}

def load_source_counts_for_audit() -> tuple[dict[str, int], dict[str, Any]]:
    fraction_path = FRACTION_DIR / "cell_fraction_feature_matrix.parquet"
    cols = ["source_h5ad_set", "total_cells_used"]
    fraction = pd.read_parquet(fraction_path, columns=cols)
    flags_check = load_step2_7_flags_check()
    counts: dict[str, int] = {}
    for row in fraction.itertuples(index=False):
        sources = [s for s in str(row.source_h5ad_set).split(";") if s and s != "nan"]
        if len(sources) == 1:
            counts[sources[0]] = counts.get(sources[0], 0) + int(row.total_cells_used)
        else:
            share = int(row.total_cells_used) // max(len(sources), 1)
            for source in sources:
                counts[source] = counts.get(source, 0) + share
    return counts, {**flags_check, "audit_source_counts_from": str(fraction_path)}

def load_cell_universe() -> tuple[pd.DataFrame, dict[str, Any]]:
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
    flags_check = load_step2_7_flags_check()
    consensus = consensus_all[consensus_all["included_in_main_annotation"].astype(bool)].copy()
    myeloid = pd.read_parquet(MYELOID_PATH, columns=myeloid_cols)
    labels = apply_myeloid_overlay(consensus, myeloid)
    return labels, flags_check


def apply_myeloid_overlay(consensus: pd.DataFrame, myeloid: pd.DataFrame) -> pd.DataFrame:
    df = consensus.merge(
        myeloid[KEY_COLS + ["myeloid_adjudication_status", "final_myeloid_label", "included_in_myeloid_main", "myeloid_sensitivity_only"]],
        on=KEY_COLS,
        how="left",
        validate="one_to_one",
    )
    df["final_major_lineage"] = df["final_major_lineage"].fillna("Unknown").astype(str)
    df["final_immune_mid_state"] = df["final_immune_mid_state"].fillna("NA").astype(str)
    df["final_myeloid_label"] = df["final_myeloid_label"].fillna("NA").astype(str)
    df["myeloid_adjudication_status"] = df["myeloid_adjudication_status"].fillna("not_myeloid_candidate").astype(str)
    df["included_in_myeloid_main"] = df["included_in_myeloid_main"].fillna(False).astype(bool)
    state = df["final_major_lineage"].copy()
    parent = df["final_major_lineage"].copy()
    immune_mask = df["final_major_lineage"].isin(IMMUNE_MAJORS)
    mid = df["final_immune_mid_state"].where(~df["final_immune_mid_state"].isin(MISSING_LABELS), "Immune_unspecified")
    state.loc[immune_mask] = mid.loc[immune_mask]
    myeloid_main = df["included_in_myeloid_main"] & df["final_myeloid_label"].isin(MYELOID_MAIN_LABELS)
    state.loc[myeloid_main] = df.loc[myeloid_main, "final_myeloid_label"]
    parent.loc[myeloid_main] = "Myeloid_DC"
    myeloid_uncertain = df["myeloid_adjudication_status"].isin(["sensitivity_only", "exclude_from_myeloid_main", "downgrade_to_myeloid_unspecified"])
    myeloid_uncertain = myeloid_uncertain & df["final_major_lineage"].isin(["Myeloid_DC", "Immune"])
    state.loc[myeloid_uncertain & (~myeloid_main)] = "Myeloid_unspecified"
    parent.loc[myeloid_uncertain & (~myeloid_main)] = "Myeloid_DC"
    df["pseudobulk_state_label"] = state.map(sanitize_label)
    df["pseudobulk_parent_lineage"] = parent.map(sanitize_label)
    return df


def build_h5ad_inventory() -> dict[str, Path]:
    inv = pd.read_csv(INPUT_INVENTORY_PATH)
    h5ad = inv[inv["file_path"].astype(str).str.endswith(".h5ad", na=False)].copy()
    return {Path(p).name: Path(p) for p in h5ad["file_path"]}


def load_var_mapping() -> pd.DataFrame:
    cols = ["source_h5ad", "var_name", "standardized_gene_symbol", "symbol_status", "duplicate_symbol_flag"]
    mapping = pd.read_csv(GENE_MAPPING_PATH, usecols=cols)
    mapping["standardized_gene_symbol"] = mapping["standardized_gene_symbol"].fillna("").astype(str).str.upper().str.strip()
    mapping["var_name"] = mapping["var_name"].astype(str)
    mapping = mapping[(mapping["symbol_status"].astype(str).eq("ok")) & (mapping["standardized_gene_symbol"] != "")]
    mapping = mapping[~mapping["duplicate_symbol_flag"].astype(str).str.lower().eq("true")]
    return mapping


def source_var_symbols(source_h5ad: str, adata_var_names: Iterable[str], mapping: pd.DataFrame) -> list[str]:
    var_names = [str(v) for v in adata_var_names]
    sub = mapping[mapping["source_h5ad"].eq(source_h5ad)]
    if sub.empty:
        return [v.upper().strip() for v in var_names]
    lookup = dict(zip(sub["var_name"].astype(str), sub["standardized_gene_symbol"].astype(str)))
    return [lookup.get(v, v.upper().strip()) for v in var_names]


def build_gene_positions(var_symbols: list[str], universe: list[str] | None = None) -> tuple[list[str], list[int], list[int]]:
    if universe is None:
        genes = []
        positions = []
        seen = set()
        for i, gene in enumerate(var_symbols):
            if gene and gene not in seen:
                genes.append(gene)
                positions.append(i)
                seen.add(gene)
        return genes, positions, list(range(len(genes)))
    first_pos: dict[str, int] = {}
    for i, gene in enumerate(var_symbols):
        if gene and gene not in first_pos:
            first_pos[gene] = i
    genes = list(universe)
    present_positions = []
    output_positions = []
    for j, gene in enumerate(genes):
        pos = first_pos.get(gene)
        if pos is not None:
            present_positions.append(pos)
            output_positions.append(j)
    order = np.argsort(present_positions).tolist() if present_positions else []
    present_positions = [present_positions[i] for i in order]
    output_positions = [output_positions[i] for i in order]
    return genes, present_positions, output_positions


def make_groups(source_cells: pd.DataFrame, scope: str) -> tuple[pd.DataFrame, np.ndarray]:
    if scope == "sample":
        group_cols = ["cohort_id", "sample_id", "source_h5ad"]
        tmp = source_cells[group_cols].copy()
        tmp["aggregation_scope"] = "sample"
        tmp["cell_state_label"] = "all_cells"
        tmp["parent_lineage"] = "all"
    elif scope == "sample_cellstate":
        group_cols = ["cohort_id", "sample_id", "source_h5ad", "pseudobulk_state_label", "pseudobulk_parent_lineage"]
        tmp = source_cells[group_cols].copy()
        tmp = tmp.rename(columns={"pseudobulk_state_label": "cell_state_label", "pseudobulk_parent_lineage": "parent_lineage"})
        tmp["aggregation_scope"] = "sample_cellstate"
    else:
        raise ValueError(f"unknown scope: {scope}")
    key_cols = ["cohort_id", "sample_id", "source_h5ad", "aggregation_scope", "cell_state_label", "parent_lineage"]
    codes, uniques = pd.factorize(pd.MultiIndex.from_frame(tmp[key_cols]), sort=False)
    groups = pd.DataFrame(list(uniques), columns=key_cols)
    counts = pd.Series(codes).value_counts().sort_index().to_numpy()
    groups["cell_count"] = counts.astype(int)
    if scope == "sample_cellstate":
        keep = groups["cell_count"].ge(MIN_STATE_CELLS).to_numpy()
        remap = np.full(len(groups), -1, dtype=int)
        remap[np.where(keep)[0]] = np.arange(int(keep.sum()))
        groups = groups.loc[keep].reset_index(drop=True)
        codes = remap[codes]
    return groups, codes.astype(int)



def h5_raw_matrix_group(h5: h5py.File, raw_count_source: str) -> Any | None:
    if raw_count_source == "layers:counts" and "layers" in h5 and "counts" in h5["layers"]:
        return h5["layers"]["counts"]
    if raw_count_source == "raw:X" and "raw" in h5 and "X" in h5["raw"]:
        return h5["raw"]["X"]
    if raw_count_source == "X" and "X" in h5:
        return h5["X"]
    return None


def is_h5_csr(group: Any) -> bool:
    return hasattr(group, "attrs") and group.attrs.get("encoding-type", "") == "csr_matrix" and all(k in group for k in ["data", "indices", "indptr"])


def aggregate_h5_csr_universe(
    group: Any,
    obs_indices: np.ndarray,
    group_codes: np.ndarray,
    groups: pd.DataFrame,
    genes: list[str],
    present_positions: list[int],
    output_positions: list[int],
    batch_size: int,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    n_groups = len(groups)
    n_genes = len(genes)
    counts = np.full((n_groups, n_genes), np.nan, dtype=np.float32)
    if present_positions:
        counts[:, output_positions] = 0.0
    valid_mask = group_codes >= 0
    obs_indices = obs_indices[valid_mask].astype(np.int64)
    group_codes = group_codes[valid_mask].astype(np.int32)
    if len(obs_indices) == 0 or not present_positions:
        gene_df = pd.DataFrame(counts, columns=genes)
        return gene_df, {"library_size": np.nansum(counts, axis=1).astype(float), "detected_genes": np.sum(np.nan_to_num(counts, nan=0.0) > 0, axis=1).astype(int)}

    order = np.argsort(obs_indices, kind="mergesort")
    obs_indices = obs_indices[order]
    group_codes = group_codes[order]

    shape = tuple(int(x) for x in group.attrs["shape"])
    col_to_out = np.full(shape[1], -1, dtype=np.int32)
    col_to_out[np.asarray(present_positions, dtype=np.int64)] = np.asarray(output_positions, dtype=np.int32)
    indptr_ds = group["indptr"]
    indices_ds = group["indices"]
    data_ds = group["data"]

    for start in range(0, len(obs_indices), batch_size):
        end = min(start + batch_size, len(obs_indices))
        rows = obs_indices[start:end]
        codes = group_codes[start:end]
        ptr_start = int(indptr_ds[int(rows[0])])
        ptr_end = int(indptr_ds[int(rows[-1]) + 1])
        if ptr_end <= ptr_start:
            continue
        block_indices = indices_ds[ptr_start:ptr_end]
        block_data = data_ds[ptr_start:ptr_end]
        ptr_window = np.asarray(indptr_ds[int(rows[0]) : int(rows[-1]) + 2], dtype=np.int64)
        row_offsets = rows - int(rows[0])
        row_starts = ptr_window[row_offsets] - ptr_start
        row_ends = ptr_window[row_offsets + 1] - ptr_start
        for local_start, local_end, code in zip(row_starts, row_ends, codes):
            if local_end <= local_start:
                continue
            cols = block_indices[local_start:local_end]
            out_cols = col_to_out[cols]
            mask = out_cols >= 0
            if mask.any():
                np.add.at(counts[code], out_cols[mask], block_data[local_start:local_end][mask])
    gene_df = pd.DataFrame(counts, columns=genes)
    return gene_df, {"library_size": np.nansum(counts, axis=1).astype(float), "detected_genes": np.sum(np.nan_to_num(counts, nan=0.0) > 0, axis=1).astype(int)}

def aggregate_dispatch(
    matrix: Any,
    h5_group: Any | None,
    obs_indices: np.ndarray,
    group_codes: np.ndarray,
    groups: pd.DataFrame,
    genes: list[str],
    present_positions: list[int],
    output_positions: list[int],
    batch_size: int,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    if h5_group is not None and is_h5_csr(h5_group):
        return aggregate_h5_csr_universe(h5_group, obs_indices, group_codes, groups, genes, present_positions, output_positions, batch_size)
    effective_batch = min(batch_size, 512) if isinstance(h5_group, h5py.Dataset) and len(h5_group.shape) == 2 else batch_size
    return aggregate_universe(matrix, obs_indices, group_codes, groups, genes, present_positions, output_positions, effective_batch)

def aggregate_universe(
    matrix: Any,
    obs_indices: np.ndarray,
    group_codes: np.ndarray,
    groups: pd.DataFrame,
    genes: list[str],
    present_positions: list[int],
    output_positions: list[int],
    batch_size: int,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    n_groups = len(groups)
    n_genes = len(genes)
    counts = np.full((n_groups, n_genes), np.nan, dtype=np.float32)
    if present_positions:
        counts[:, output_positions] = 0.0
    valid_mask = group_codes >= 0
    obs_indices = obs_indices[valid_mask]
    group_codes = group_codes[valid_mask]
    for start in range(0, len(obs_indices), batch_size):
        end = min(start + batch_size, len(obs_indices))
        batch_obs = obs_indices[start:end]
        batch_groups = group_codes[start:end]
        if len(batch_obs) > 1:
            order = np.argsort(batch_obs, kind="mergesort")
            batch_obs = batch_obs[order]
            batch_groups = batch_groups[order]
        x = get_matrix_slice(matrix, batch_obs, np.asarray(present_positions, dtype=int)) if present_positions else sparse.csr_matrix((len(batch_obs), 0))
        indicator = sparse.csr_matrix((np.ones(len(batch_groups)), (batch_groups, np.arange(len(batch_groups)))), shape=(n_groups, len(batch_groups)))
        summed = indicator @ x
        if sparse.issparse(summed):
            summed = summed.toarray()
        counts[:, output_positions] += np.asarray(summed, dtype=np.float32)
    gene_df = pd.DataFrame(counts, columns=genes)
    library_size = np.nansum(counts, axis=1)
    detected_genes = np.sum(np.nan_to_num(counts, nan=0.0) > 0, axis=1)
    stats = {"library_size": library_size.astype(float), "detected_genes": detected_genes.astype(int)}
    return gene_df, stats


def make_count_matrix(
    groups: pd.DataFrame,
    gene_df: pd.DataFrame,
    stats: dict[str, Any],
    universe_name: str,
    raw_choice: RawMatrixChoice,
    created_at: str,
) -> pd.DataFrame:
    meta = groups.copy()
    meta.insert(0, "input_manifest_ref", str(INPUT_MANIFEST_REF))
    meta.insert(0, "created_at", created_at)
    meta.insert(0, "run_id", RUN_ID)
    meta["row_id"] = [
        f"{universe_name}|{r.source_h5ad}|{r.cohort_id}|{r.sample_id}|{r.aggregation_scope}|{r.cell_state_label}"
        for r in meta.itertuples(index=False)
    ]
    meta["library_size"] = stats["library_size"]
    meta["detected_genes"] = stats["detected_genes"]
    meta["gene_universe"] = universe_name
    meta["raw_count_source"] = raw_choice.raw_count_source
    meta["raw_count_status"] = raw_choice.raw_count_status
    meta = meta[ROW_META_COLS]
    return pd.concat([meta.reset_index(drop=True), gene_df.reset_index(drop=True)], axis=1)



def subset_count_matrix(counts: pd.DataFrame, genes: list[str], universe_name: str) -> pd.DataFrame:
    meta = counts[ROW_META_COLS].copy()
    available = [g for g in genes if g in counts.columns]
    gene_df = counts.reindex(columns=genes)
    vals = gene_df.to_numpy(dtype=float)
    meta["gene_universe"] = universe_name
    meta["library_size"] = np.nansum(vals, axis=1).astype(float)
    meta["detected_genes"] = np.sum(np.nan_to_num(vals, nan=0.0) > 0, axis=1).astype(int)
    meta["row_id"] = [
        f"{universe_name}|{r.source_h5ad}|{r.cohort_id}|{r.sample_id}|{r.aggregation_scope}|{r.cell_state_label}"
        for r in meta.itertuples(index=False)
    ]
    return pd.concat([meta.reset_index(drop=True), gene_df.reset_index(drop=True)], axis=1)


def ordered_union(left: list[str], right: list[str]) -> list[str]:
    out: list[str] = []
    seen: set[str] = set()
    for gene in list(left) + list(right):
        if gene not in seen:
            out.append(gene)
            seen.add(gene)
    return out

def counts_to_logcpm(counts: pd.DataFrame) -> pd.DataFrame:
    out = counts[ROW_META_COLS].copy()
    genes = [c for c in counts.columns if c not in ROW_META_COLS]
    vals = counts[genes].to_numpy(dtype=float)
    libs = counts["library_size"].to_numpy(dtype=float)
    with np.errstate(divide="ignore", invalid="ignore"):
        norm = np.log1p((vals / libs[:, None]) * 1_000_000.0)
    norm[~np.isfinite(norm)] = np.nan
    return pd.concat([out, pd.DataFrame(norm, columns=genes)], axis=1)


def forbidden_response_cols(df: pd.DataFrame) -> list[str]:
    return [c for c in df.columns if any(term in c.lower() for term in FORBIDDEN_RESPONSE_TERMS)]


def validate_pair(counts: pd.DataFrame, logcpm: pd.DataFrame) -> dict[str, Any]:
    gene_cols_counts = [c for c in counts.columns if c not in ROW_META_COLS]
    gene_cols_log = [c for c in logcpm.columns if c not in ROW_META_COLS]
    return {
        "rows": int(len(counts)),
        "gene_columns": int(len(gene_cols_counts)),
        "row_ids_match": counts["row_id"].tolist() == logcpm["row_id"].tolist(),
        "gene_columns_match": gene_cols_counts == gene_cols_log,
        "forbidden_response_cols": forbidden_response_cols(counts) + forbidden_response_cols(logcpm),
        "min_cell_count": int(counts["cell_count"].min()) if len(counts) else 0,
    }


def write_empty_matrix(path: Path, universe_name: str, created_at: str) -> pd.DataFrame:
    df = pd.DataFrame(columns=ROW_META_COLS)
    df.to_parquet(path, index=False)
    return df


def run(args: argparse.Namespace) -> None:
    global OUT_DIR
    if args.out_dir:
        OUT_DIR = Path(args.out_dir)
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    created_at = now_iso()
    if args.audit_only:
        labels = None
        source_counts, flags_check = load_source_counts_for_audit()
    else:
        labels, flags_check = load_cell_universe()
        source_counts = labels["source_h5ad"].value_counts().to_dict()
    h5ad_inventory = build_h5ad_inventory()
    mapping = load_var_mapping()
    core_genes = read_gene_list(CORE_GENES_PATH)
    immune_genes = read_gene_list(IMMUNE_GENES_PATH)

    audit_rows = []
    for source_h5ad, n in sorted(source_counts.items()):
        path = h5ad_inventory.get(source_h5ad)
        audit_rows.append(audit_one_h5ad(path or Path(source_h5ad), source_h5ad, int(n)))
    audit_df = pd.DataFrame(audit_rows)
    audit_df.to_csv(OUT_DIR / "raw_count_authority_by_source.csv", index=False)
    write_raw_authority_report(audit_df, OUT_DIR / "raw_count_authority_report.md", created_at)

    eligible_statuses = set(RAW_STATUSES_FOR_MAIN)
    if args.include_dense_deferred:
        eligible_statuses.add("raw_like_dense_deferred")
    eligible_sources = set(audit_df.loc[audit_df["raw_count_status"].isin(eligible_statuses), "source_h5ad"])
    if args.source_list:
        requested = {x.strip() for x in args.source_list.split(",") if x.strip()}
        eligible_sources = requested & eligible_sources
    if args.source:
        eligible_sources = {args.source} & eligible_sources
    blocked_sources = sorted(set(source_counts) - eligible_sources)
    if args.audit_only or not eligible_sources:
        manifest = build_manifest(created_at, flags_check, audit_df, {}, {}, audit_only=args.audit_only)
        write_manifest_and_report(manifest, {}, audit_df, pd.DataFrame(), created_at)
        return

    outputs: dict[str, list[pd.DataFrame]] = {
        "sample_native": [],
        "sample_core": [],
        "sample_cellstate_core": [],
        "sample_cellstate_immune": [],
    }
    qc_rows = []
    sample_meta_rows = []
    max_sources = args.max_sources if args.max_sources and args.max_sources > 0 else None
    processed = 0
    for source_h5ad in sorted(eligible_sources):
        if max_sources is not None and processed >= max_sources:
            break
        print(f"[Step2.8] processing {source_h5ad} ({processed + 1}/{len(eligible_sources)})", flush=True)
        source_cells = labels[labels["source_h5ad"].eq(source_h5ad)].copy()
        if source_cells.empty:
            continue
        path = h5ad_inventory[source_h5ad]
        audit_row = audit_df[audit_df["source_h5ad"].eq(source_h5ad)].iloc[0].to_dict()
        a = ad.read_h5ad(path, backed="r")
        h5 = h5py.File(path, "r")
        try:
            raw_choice = choose_raw_matrix(a, source_h5ad, audit_row)
            if raw_choice is not None and args.include_dense_deferred and audit_row.get("raw_count_status") == "raw_like_dense_deferred":
                raw_choice = RawMatrixChoice(raw_choice.matrix, raw_choice.raw_count_source, "raw_like_dense_rescued", True)
            if raw_choice is None or not raw_choice.raw_count_accept_main:
                continue
            h5_group = h5_raw_matrix_group(h5, raw_choice.raw_count_source)
            obs_lookup = pd.Series(np.arange(a.n_obs), index=a.obs_names.astype(str))
            if not obs_lookup.index.is_unique:
                obs_lookup = obs_lookup.groupby(level=0, sort=False).first()
            idx = obs_lookup.reindex(source_cells["cell_barcode"].astype(str)).to_numpy()
            missing = int(pd.isna(idx).sum())
            if missing:
                source_cells = source_cells.loc[~pd.isna(idx)].copy()
                idx = idx[~pd.isna(idx)]
            obs_indices = idx.astype(int)
            var_symbols = source_var_symbols(source_h5ad, a.var_names, mapping)

            sample_groups, sample_codes = make_groups(source_cells, "sample")
            if not sample_groups.empty:
                native_genes, native_positions, native_outpos = build_gene_positions(var_symbols, None)
                gene_df, stats = aggregate_dispatch(raw_choice.matrix, h5_group, obs_indices, sample_codes, sample_groups, native_genes, native_positions, native_outpos, args.batch_size)
                sample_native = make_count_matrix(sample_groups, gene_df, stats, "native", raw_choice, created_at)
                outputs["sample_native"].append(sample_native)
                outputs["sample_core"].append(subset_count_matrix(sample_native, core_genes, "core"))

            state_groups, state_codes = make_groups(source_cells, "sample_cellstate")
            if not state_groups.empty:
                state_union_genes = ordered_union(immune_genes, core_genes)
                state_genes, state_positions, state_outpos = build_gene_positions(var_symbols, state_union_genes)
                gene_df, stats = aggregate_dispatch(raw_choice.matrix, h5_group, obs_indices, state_codes, state_groups, state_genes, state_positions, state_outpos, args.batch_size)
                state_union = make_count_matrix(state_groups, gene_df, stats, "state_union", raw_choice, created_at)
                outputs["sample_cellstate_immune"].append(subset_count_matrix(state_union, immune_genes, "immune"))
                outputs["sample_cellstate_core"].append(subset_count_matrix(state_union, core_genes, "core"))
            sample_meta_rows.append(source_cells.groupby(["cohort_id", "sample_id", "source_h5ad"]).size().reset_index(name="cell_count"))
            qc_rows.append({
                "source_h5ad": source_h5ad,
                "h5ad_path": str(path),
                "cells_requested": int(len(labels[labels["source_h5ad"].eq(source_h5ad)])),
                "cells_matched_obs": int(len(source_cells)),
                "cells_missing_obs": missing,
                "raw_count_source": raw_choice.raw_count_source,
                "raw_count_status": raw_choice.raw_count_status,
            })
        finally:
            h5.close()
            a.file.close()
        processed += 1

    file_map = {
        "sample_native": ("pseudobulk_counts_by_sample.native.parquet", "pseudobulk_logcpm_by_sample.native.parquet"),
        "sample_core": ("pseudobulk_counts_by_sample_core.parquet", "pseudobulk_logcpm_by_sample_core.parquet"),
        "sample_cellstate_core": ("pseudobulk_counts_by_sample_cellstate_core.parquet", "pseudobulk_logcpm_by_sample_cellstate_core.parquet"),
        "sample_cellstate_immune": ("pseudobulk_counts_by_sample_cellstate_immune.parquet", "pseudobulk_logcpm_by_sample_cellstate_immune.parquet"),
    }
    validations = {}
    for key, (count_name, log_name) in file_map.items():
        if outputs[key]:
            counts = pd.concat(outputs[key], ignore_index=True, sort=False)
        else:
            counts = pd.DataFrame(columns=ROW_META_COLS)
        counts_path = OUT_DIR / count_name
        log_path = OUT_DIR / log_name
        counts.to_parquet(counts_path, index=False)
        logcpm = counts_to_logcpm(counts) if len(counts) else pd.DataFrame(columns=counts.columns)
        logcpm.to_parquet(log_path, index=False)
        validations[key] = validate_pair(counts, logcpm)

    sample_meta = pd.concat(sample_meta_rows, ignore_index=True) if sample_meta_rows else pd.DataFrame(columns=["cohort_id", "sample_id", "source_h5ad", "cell_count"])
    sample_meta.to_csv(OUT_DIR / "pseudobulk_sample_metadata.csv", index=False)
    qc = pd.DataFrame(qc_rows)
    if blocked_sources:
        blocked = audit_df[audit_df["source_h5ad"].isin(blocked_sources)][["source_h5ad", "raw_count_source", "raw_count_status", "raw_count_reason", "n_obs_in_main_universe"]]
        qc = pd.concat([qc, blocked.assign(cells_requested=blocked["n_obs_in_main_universe"], cells_matched_obs=0, cells_missing_obs=blocked["n_obs_in_main_universe"])], ignore_index=True, sort=False)
    qc.to_csv(OUT_DIR / "pseudobulk_qc_summary.csv", index=False)
    gene_manifest = build_gene_manifest(core_genes, immune_genes, audit_df, outputs)
    (OUT_DIR / "pseudobulk_gene_universe_manifest.yaml").write_text(yaml.safe_dump(gene_manifest, sort_keys=False), encoding="utf-8")
    manifest = build_manifest(created_at, flags_check, audit_df, validations, gene_manifest, audit_only=False)
    write_manifest_and_report(manifest, validations, audit_df, qc, created_at)


def build_gene_manifest(core_genes: list[str], immune_genes: list[str], audit_df: pd.DataFrame, outputs: dict[str, list[pd.DataFrame]]) -> dict[str, Any]:
    return {
        "native_gene_universe": "source-native standardized gene symbols; per-source columns are unioned with missing genes as NA",
        "core_intersection_gene_count": len(core_genes),
        "immune_feature_gene_count": len(immune_genes),
        "missing_gene_policy": "NA_not_zero",
        "measured_zero_policy": "0",
        "expression_imputation": "forbidden",
        "eligible_source_count": int(audit_df["raw_count_status"].isin(RAW_STATUSES_FOR_MAIN).sum()),
    }


def build_manifest(created_at: str, flags_check: dict[str, Any], audit_df: pd.DataFrame, validations: dict[str, Any], gene_manifest: dict[str, Any], audit_only: bool) -> dict[str, Any]:
    eligible = int(audit_df["raw_count_status"].isin(RAW_STATUSES_FOR_MAIN).sum()) if len(audit_df) else 0
    total = int(len(audit_df))
    if audit_only:
        gate = "raw_count_audit_only"
    elif eligible == 0:
        gate = "blocked_no_raw_count_sources"
    elif eligible < total:
        gate = "partial_ready_with_raw_count_exclusions"
    else:
        gate = "pass"
    input_paths = [CONSENSUS_PATH, MYELOID_PATH, FRACTION_MANIFEST_PATH, FLAGS_PATH, INPUT_INVENTORY_PATH, GENE_MAPPING_PATH, CORE_GENES_PATH, IMMUNE_GENES_PATH, GENE_MANIFEST_PATH]
    return {
        "run_id": RUN_ID,
        "step_id": "Step2.8_pseudobulk_feature_factory_v1",
        "created_at": created_at,
        "input_manifest_ref": str(INPUT_MANIFEST_REF),
        "gate_status": gate,
        "input_hashes": {str(p): sha256_file(p) for p in input_paths if p.exists()},
        "flags_consensus_crosscheck": flags_check,
        "raw_count_authority": {
            "total_sources": total,
            "eligible_sources": eligible,
            "blocked_sources": total - eligible,
            "eligible_statuses": sorted(RAW_STATUSES_FOR_MAIN),
            "dense_deferred_policy": "excluded_by_default; may be rescued with --include-dense-deferred",
        },
        "thresholds": {"sample_cellstate_min_cells": MIN_STATE_CELLS, "max_dense_main_entries": MAX_DENSE_MAIN_ENTRIES},
        "gene_universe": gene_manifest,
        "validation": validations,
        "forbidden_operations": ["response_leakage", "differential_expression", "model_training", "global_clustering"],
    }


def write_manifest_and_report(manifest: dict[str, Any], validations: dict[str, Any], audit_df: pd.DataFrame, qc: pd.DataFrame, created_at: str) -> None:
    (OUT_DIR / "pseudobulk_decision_manifest.yaml").write_text(yaml.safe_dump(manifest, sort_keys=False), encoding="utf-8")
    report = [
        "# Step2.8 Pseudobulk Feature Factory Report",
        "",
        f"- created_at: `{created_at}`",
        f"- gate_status: `{manifest['gate_status']}`",
        f"- total sources audited: `{manifest['raw_count_authority']['total_sources']}`",
        f"- eligible raw-count sources: `{manifest['raw_count_authority']['eligible_sources']}`",
        f"- blocked sources: `{manifest['raw_count_authority']['blocked_sources']}`",
        "- response usage: `forbidden_not_used`",
        "- differential analysis: `not_run`",
        "- global clustering: `not_run`",
        "- model training: `not_run`",
        "",
        "## Output Validation",
        "",
    ]
    if validations:
        for key, val in validations.items():
            report.append(f"- {key}: rows={val.get('rows')}, genes={val.get('gene_columns')}, row_ids_match={val.get('row_ids_match')}, gene_columns_match={val.get('gene_columns_match')}, forbidden_response_cols={val.get('forbidden_response_cols')}")
    else:
        report.append("- no pseudobulk matrix emitted in this run")
    if len(audit_df):
        report += ["", "## Raw Count Status", ""]
        for status, n in audit_df["raw_count_status"].value_counts().sort_index().items():
            report.append(f"- {status}: {int(n)}")
    (OUT_DIR / "pseudobulk_feature_factory_report.md").write_text("\n".join(report) + "\n", encoding="utf-8")


def write_raw_authority_report(audit_df: pd.DataFrame, path: Path, created_at: str) -> None:
    lines = [
        "# Step2.8 Raw Count Authority Report",
        "",
        f"- created_at: `{created_at}`",
        f"- sources audited: `{len(audit_df)}`",
        "- accepted for main: `confirmed_raw_counts` or `raw_like_counts`",
        "- blocked from main: `unclear_or_normalized` or `unavailable`",
        "",
        "## Status Counts",
        "",
    ]
    if len(audit_df):
        for status, n in audit_df["raw_count_status"].value_counts().sort_index().items():
            lines.append(f"- {status}: {int(n)}")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Step2.8 pseudobulk feature factory v1")
    parser.add_argument("--audit-only", action="store_true", help="run raw-count authority audit and manifests only")
    parser.add_argument("--batch-size", type=int, default=4096, help="cell batch size for source-wise aggregation")
    parser.add_argument("--max-sources", type=int, default=0, help="optional source limit for smoke runs")
    parser.add_argument("--out-dir", default="", help="optional output directory override for smoke runs")
    parser.add_argument("--source", default="", help="optional single source_h5ad filter for debugging")
    parser.add_argument("--source-list", default="", help="comma-separated source_h5ad filter for rescue runs")
    parser.add_argument("--include-dense-deferred", action="store_true", help="include raw_like_dense_deferred sources for dense rescue runs")
    return parser.parse_args()


if __name__ == "__main__":
    run(parse_args())
