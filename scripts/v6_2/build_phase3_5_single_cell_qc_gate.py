#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import os
import re
import sys
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable

import yaml

try:
    import numpy as np
except Exception as exc:  # pragma: no cover - checked at runtime
    np = None  # type: ignore
    NUMPY_IMPORT_ERROR = exc
else:
    NUMPY_IMPORT_ERROR = None

try:
    import anndata as ad
except Exception as exc:  # pragma: no cover - checked at runtime
    ad = None
    ANNDATA_IMPORT_ERROR = exc
else:
    ANNDATA_IMPORT_ERROR = None

try:
    import h5py
except Exception as exc:  # pragma: no cover - checked at runtime
    h5py = None
    H5PY_IMPORT_ERROR = exc
else:
    H5PY_IMPORT_ERROR = None


ROOT = Path(__file__).resolve().parents[2]
INPUT_VERSION = "frozen_v0"
OUT_DIR = ROOT / "results" / "v6_2" / "phase3_5_single_cell_processing_qc_gate"

FROZEN_SAMPLE = ROOT / "results" / "v6_2" / "sample_metadata_master.frozen_v0.csv"
FROZEN_PATIENT = ROOT / "results" / "v6_2" / "patient_metadata_master.frozen_v0.csv"
FROZEN_ROLE = ROOT / "results" / "v6_2" / "dataset_role_and_feature_eligibility.frozen_v0.csv"
FROZEN_RESPONSE = ROOT / "results" / "v6_2" / "response_label_environment.frozen_v0.csv"
FROZEN_LABEL_INV = ROOT / "results" / "v6_2" / "label_inventory_and_semantic_review.csv"
FROZEN_MANIFEST = ROOT / "results" / "v6_2" / "frozen_input_manifest.yaml"
H5AD_INVENTORY = ROOT / "results" / "v6_1" / "data_pool" / "snapshots" / "current_step2_input" / "h5ad_inventory.csv"
PHASE25 = ROOT / "results" / "v6_2" / "phase2_5_data_onboarding"
EXTERNAL = ROOT / "data" / "external_anchor_candidates"
GSE301741_DIR = EXTERNAL / "GSE301741_hnscc_pembrolizumab_scRNA"

CAUTION_PREFIXES = ("MT-", "RPL", "RPS", "HB")
CAUTION_GENES = {
    "MALAT1", "FOS", "JUN", "JUNB", "JUND", "HSPA1A", "HSPA1B", "HSP90AA1",
    "HSP90AB1", "MKI67", "TOP2A", "PCNA", "STMN1", "TUBB", "TUBA1B",
}
IMMUNE_CORE = {
    "CD3D", "CD3E", "CD3G", "CD4", "CD8A", "CD8B", "NKG7", "GNLY", "GZMB",
    "GZMA", "PRF1", "IFNG", "CXCL9", "CXCL10", "STAT1", "IRF1", "HLA-A",
    "HLA-B", "HLA-C", "B2M", "HLA-DRA", "HLA-DRB1", "CD74", "LST1",
    "LYZ", "S100A8", "S100A9", "FCGR3A", "MS4A1", "CD79A", "MZB1",
    "JCHAIN", "FOXP3", "IL2RA", "CTLA4", "PDCD1", "LAG3", "HAVCR2",
    "TIGIT", "CXCL13", "TCF7", "CCR7", "ITGAX", "FCER1A", "CLEC9A",
    "COL1A1", "COL1A2", "ACTA2", "PECAM1", "VWF", "KDR", "EPCAM", "KRT8",
    "KRT18", "KRT19", "ALB", "APOA1", "VEGFA", "TGFB1", "TGFB2",
}
CELL_STATE_MARKERS = {
    "T_NK": {"CD3D", "CD3E", "CD8A", "CD4", "NKG7", "GNLY", "GZMB", "PRF1"},
    "B_plasma": {"MS4A1", "CD79A", "CD79B", "MZB1", "JCHAIN"},
    "myeloid": {"LYZ", "LST1", "S100A8", "S100A9", "FCGR3A", "C1QA", "C1QB", "C1QC"},
    "DC_APC": {"HLA-DRA", "HLA-DRB1", "CD74", "ITGAX", "FCER1A", "CLEC9A"},
    "stromal_CAF": {"COL1A1", "COL1A2", "ACTA2", "DCN", "LUM"},
    "endothelial": {"PECAM1", "VWF", "KDR", "ENG"},
    "tumor_like": {"EPCAM", "KRT8", "KRT18", "KRT19", "ALB", "APOA1"},
}


def now_iso() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat()


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def write_csv(path: Path, rows: list[dict[str, Any]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({k: row.get(k, "") for k in fieldnames})


def write_tsv(path: Path, rows: list[dict[str, Any]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({k: row.get(k, "") for k in fieldnames})


def file_fingerprint(path: Path) -> tuple[str, str, int, str]:
    if not path.exists():
        return "missing", "", 0, ""
    st = path.stat()
    h = hashlib.sha256()
    with path.open("rb") as f:
        h.update(f.read(16 * 1024 * 1024))
    return "partial_sha256_16mb", h.hexdigest(), st.st_size, str(int(st.st_mtime))


def norm_path(path: str | Path) -> str:
    return str(path).replace(str(ROOT), "/home/huyudi/006")


def upper_gene(gene: str) -> str:
    return str(gene).strip().upper()


def is_caution_gene(gene: str) -> bool:
    g = upper_gene(gene)
    return g in CAUTION_GENES or g.startswith(CAUTION_PREFIXES)


# ---- v6.2 repair: robust gene-symbol resolution and count-like X detection ----

ENSEMBL_RE = re.compile(r"^ENSG\d{11}(\.\d+)?$", re.IGNORECASE)


def looks_like_ensembl(val: str) -> bool:
    return bool(ENSEMBL_RE.match(str(val).strip()))


def looks_like_numeric_ids(vals: list[str], n: int = 100) -> bool:
    sample = vals[:n]
    if not sample:
        return False
    return all(str(v).strip().isdigit() for v in sample)


def load_ensembl_to_symbol(gtf_path: Path | None = None) -> dict[str, str]:
    """Build Ensembl gene_id -> gene_name map from a Gencode GTF."""
    if gtf_path is None:
        candidates = [
            ROOT / "data" / "ref" / "gencode.v49.basic.annotation.gtf",
            ROOT / "data" / "ref" / "gencode.v49.annotation.gtf",
            ROOT / "data" / "db" / "tcga" / "gencode.v23.annotation.gtf",
        ]
        gtf_path = next((p for p in candidates if p.exists()), None)
    mapping: dict[str, str] = {}
    if gtf_path is None or not gtf_path.exists():
        return mapping
    opener = gzip.open if str(gtf_path).endswith(".gz") else open
    with opener(gtf_path, "rt", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9 or parts[2] != "gene":
                continue
            attrs = parts[8]
            gid = ""
            name = ""
            for token in attrs.split(";"):
                token = token.strip()
                if token.startswith("gene_id"):
                    gid = token.split("\"", 2)[1].split(".")[0] if "\"" in token else ""
                elif token.startswith("gene_name"):
                    name = token.split("\"", 2)[1] if "\"" in token else ""
                if gid and name:
                    mapping[gid] = name
                    break
    return mapping


# Lazy singleton so we only parse the GTF once per run.
_ENSEMBL_MAP: dict[str, str] | None = None


def ensembl_to_symbol_map() -> dict[str, str]:
    global _ENSEMBL_MAP
    if _ENSEMBL_MAP is None:
        _ENSEMBL_MAP = load_ensembl_to_symbol()
    return _ENSEMBL_MAP


def resolve_gene_symbols(path: Path, var_names: list[str]) -> tuple[list[str], str, str]:
    """Return (uppercased symbols, id_field, symbol_field).

    Fallback order:
    1. var_names if they already look like symbols.
    2. var columns named gene_symbol/symbol/Gene/etc.
    3. Ensembl -> symbol mapping from GTF.
    4. numeric ids -> var columns.
    """
    field_id = "var_names"
    field_symbol = "var_names"
    symbols: list[str] = []
    sample = [str(v).strip() for v in var_names[:200]]

    # Heuristic: symbols contain letters and are not Ensembl-like.
    likely_symbols = sample and not any(looks_like_ensembl(v) for v in sample) and not looks_like_numeric_ids(sample)
    if likely_symbols:
        return [upper_gene(v) for v in var_names], field_id, field_symbol

    # Need anndata-backed read for var columns.
    if ad is None:
        return [upper_gene(v) for v in var_names], field_id, "anndata_unavailable"

    try:
        a = ad.read_h5ad(path, backed="r")
    except Exception:
        return [upper_gene(v) for v in var_names], field_id, "read_failed"

    var_cols = list(a.var.columns)
    preferred = ["gene_symbol", "gene_symbols", "symbol", "Gene", "gene_name", "feature_name", "gene"]
    for col in preferred:
        if col in var_cols:
            try:
                symbols = [upper_gene(str(v)) for v in a.var[col].astype(str)]
                field_symbol = col
                field_id = "var_names" if likely_symbols else "ensembl_or_numeric"
                a.file.close()
                return symbols, field_id, field_symbol
            except Exception:
                continue

    # Ensembl mapping fallback.
    if sample and any(looks_like_ensembl(v) for v in sample):
        emap = ensembl_to_symbol_map()
        symbols = []
        for v in var_names:
            bare = str(v).strip().split(".")[0]
            symbols.append(upper_gene(emap.get(bare, bare)))
        field_id = "ensembl_id"
        field_symbol = "gencode_gtf_mapped"
        a.file.close()
        return symbols, field_id, field_symbol

    a.file.close()
    return [upper_gene(v) for v in var_names], field_id, "unresolved"


def is_matrix_countlike(mat, max_cells: int = 2000) -> bool:
    """Inspect a matrix (dense/sparse) and declare count-like if >99% non-zero values are integers and <1% negative."""
    if np is None:
        return False
    try:
        n_obs = mat.shape[0]
        n = min(max_cells, n_obs)
        if n <= 0:
            return False
        idx = np.sort(np.random.default_rng(42).choice(n_obs, size=n, replace=False))
        X = mat[idx]
        if hasattr(X, "toarray"):
            X = X.toarray()
        X = np.asarray(X)
        nz = X[X != 0]
        if nz.size == 0:
            return False
        int_frac = (nz == nz.astype(np.int64)).mean()
        neg_frac = (nz < 0).mean()
        return bool(int_frac >= 0.99 and neg_frac <= 0.01)
    except Exception:
        return False


def is_x_countlike(path: Path, max_cells: int = 2000) -> bool:
    """Sample X and declare count-like if >99% of non-zero values are integers and <1% negative."""
    if np is None or ad is None:
        return False
    try:
        a = ad.read_h5ad(path, backed="r")
        ok = is_matrix_countlike(a.X, max_cells)
        a.file.close()
        return ok
    except Exception:
        return False


def is_raw_countlike(path: Path, max_cells: int = 2000) -> bool:
    """Sample raw.X if a raw object exists and declare count-like."""
    if np is None or ad is None:
        return False
    try:
        a = ad.read_h5ad(path, backed="r")
        if a.raw is None:
            a.file.close()
            return False
        ok = is_matrix_countlike(a.raw.X, max_cells)
        a.file.close()
        return ok
    except Exception:
        return False


def infer_sample_column(columns: Iterable[str]) -> str:
    cols = list(columns)
    preferred = [
        "sample_key", "sample_id", "sample", "orig.ident", "orig_ident", "library_id",
        "patient_sample", "Sample", "sampleID", "sample_name", "donor", "patient_id",
    ]
    lower_map = {c.lower(): c for c in cols}
    for col in preferred:
        if col in cols:
            return col
        if col.lower() in lower_map:
            return lower_map[col.lower()]
    for c in cols:
        low = c.lower()
        if "sample" in low or "library" in low or "orig" in low:
            return c
    return ""


def infer_patient_column(columns: Iterable[str]) -> str:
    cols = list(columns)
    preferred = ["patient_key", "patient_id", "patient", "donor", "case", "Patient", "Subject"]
    lower_map = {c.lower(): c for c in cols}
    for col in preferred:
        if col in cols:
            return col
        if col.lower() in lower_map:
            return lower_map[col.lower()]
    for c in cols:
        low = c.lower()
        if "patient" in low or "donor" in low or "subject" in low:
            return c
    return ""


def find_annotation_columns(columns: Iterable[str]) -> list[str]:
    hits = []
    patterns = ("celltype", "cell_type", "annotation", "lineage", "cluster", "seurat_clusters", "major", "minor")
    for col in columns:
        low = col.lower()
        if any(p in low for p in patterns):
            hits.append(col)
    return hits[:12]


def detect_qc_columns(columns: Iterable[str]) -> dict[str, str]:
    cols = list(columns)
    lower = {c.lower(): c for c in cols}
    aliases = {
        "n_genes": ["n_genes_by_counts", "nfeature_rna", "n_genes", "detected_genes"],
        "total_counts": ["total_counts", "ncount_rna", "umi", "umis", "n_umi"],
        "pct_mito": ["pct_counts_mito", "percent.mt", "pct_mito", "mito_percent"],
        "pct_ribo": ["pct_counts_ribo", "percent.ribo", "pct_ribo"],
        "pct_hb": ["pct_counts_hb", "percent.hb", "pct_hb"],
        "doublet_score": ["doublet_score", "scrublet_score", "doublet_scores"],
    }
    result = {}
    for key, names in aliases.items():
        for name in names:
            if name in cols:
                result[key] = name
                break
            if name.lower() in lower:
                result[key] = lower[name.lower()]
                break
    return result


def decode_h5_value(value: Any) -> str:
    if isinstance(value, bytes):
        return value.decode("utf-8", errors="replace")
    return str(value)


def h5_dataframe_columns(group: Any) -> list[str]:
    try:
        order = group.attrs.get("column-order", [])
        return [decode_h5_value(x) for x in list(order)]
    except Exception:
        return [str(k) for k in group.keys() if str(k) != "_index"]


def h5_dataset_len(group: Any) -> int:
    if "_index" in group:
        return int(group["_index"].shape[0])
    for key in group.keys():
        obj = group[key]
        if hasattr(obj, "shape") and obj.shape:
            return int(obj.shape[0])
    return 0


def h5_read_string_dataset(group: Any, key: str, limit: int | None = None) -> list[str]:
    if key not in group:
        return []
    obj = group[key]
    try:
        data = obj[:limit] if limit else obj[:]
        return [decode_h5_value(x) for x in data]
    except Exception:
        return []


def h5_nunique(group: Any, key: str, max_values: int = 200000) -> int | str:
    vals = h5_read_string_dataset(group, key, max_values)
    if not vals:
        return ""
    return len(set(vals))


def safe_numeric_summary(series: Any) -> tuple[str, str, str]:
    try:
        import numpy as np
        import pandas as pd

        vals = pd.to_numeric(series, errors="coerce").dropna()
        if vals.empty:
            return "", "", ""
        return (
            f"{float(vals.median()):.4g}",
            f"{float(vals.quantile(0.05)):.4g}",
            f"{float(vals.quantile(0.95)):.4g}",
        )
    except Exception:
        return "", "", ""


def inspect_h5ad(path: Path, cohort_id: str) -> dict[str, Any]:
    row: dict[str, Any] = {
        "cohort_id": cohort_id,
        "object_path": norm_path(path),
        "object_type": "h5ad",
        "object_read_status": "not_read",
        "read_error": "",
        "n_cells": "",
        "n_genes": "",
        "obs_columns": "",
        "var_columns": "",
        "layers": "",
        "raw_available": "unknown",
        "selected_layer": "",
        "layer_decision": "support_only",
        "gene_id_field": "",
        "gene_symbol_field": "",
        "sample_column": "",
        "patient_column": "",
        "annotation_columns": "",
        "qc_columns": "",
        "n_samples_in_object": "",
        "n_patients_in_object": "",
        "n_gene_symbols": "",
        "n_duplicated_gene_symbols": "",
        "n_caution_genes": "",
        "immune_core_covered": "",
        "immune_core_total": len(IMMUNE_CORE),
        "major_cell_state_support": "",
        "hvg_source": "",
        "hvg_genes": "",
    }
    if h5py is None:
        row["object_read_status"] = "blocked_import_error"
        row["read_error"] = repr(H5PY_IMPORT_ERROR)
        return row
    try:
        f = h5py.File(path, "r")
        row["object_read_status"] = "readable"
        obs_group = f["obs"] if "obs" in f else None
        var_group = f["var"] if "var" in f else None
        obs_cols = h5_dataframe_columns(obs_group) if obs_group is not None else []
        var_cols = h5_dataframe_columns(var_group) if var_group is not None else []
        row["n_cells"] = h5_dataset_len(obs_group) if obs_group is not None else ""
        row["n_genes"] = h5_dataset_len(var_group) if var_group is not None else ""
        row["obs_columns"] = "|".join(obs_cols[:80])
        row["var_columns"] = "|".join(var_cols[:80])
        layer_keys = list(f["layers"].keys()) if "layers" in f else []
        row["layers"] = "|".join(map(str, layer_keys))
        row["raw_available"] = "yes" if "raw" in f else "no"
        sample_col = infer_sample_column(obs_cols)
        patient_col = infer_patient_column(obs_cols)
        anno_cols = find_annotation_columns(obs_cols)
        qc_cols = detect_qc_columns(obs_cols)
        row["sample_column"] = sample_col
        row["patient_column"] = patient_col
        row["annotation_columns"] = "|".join(anno_cols)
        row["qc_columns"] = json.dumps(qc_cols, sort_keys=True)
        if obs_group is not None and sample_col:
            row["n_samples_in_object"] = h5_nunique(obs_group, sample_col)
        if obs_group is not None and patient_col:
            row["n_patients_in_object"] = h5_nunique(obs_group, patient_col)
        layer_set = set(map(str, layer_keys))
        if "counts" in layer_set:
            row["selected_layer"] = "layers[counts]"
            row["layer_decision"] = "raw_counts"
        elif "raw_counts" in layer_set:
            row["selected_layer"] = "layers[raw_counts]"
            row["layer_decision"] = "raw_counts"
        elif "raw" in f:
            # v6.2 repair: raw.X is not always counts; inspect values.
            if is_raw_countlike(path):
                row["selected_layer"] = "raw.X"
                row["layer_decision"] = "raw_or_author_raw"
            elif is_x_countlike(path):
                row["selected_layer"] = "X"
                row["layer_decision"] = "raw_counts_in_X"
            else:
                row["selected_layer"] = "X"
                row["layer_decision"] = "processed_or_unknown_X"
        elif "X" in f:
            # v6.2 repair: inspect X values before declaring processed.
            if is_x_countlike(path):
                row["selected_layer"] = "X"
                row["layer_decision"] = "raw_counts_in_X"
            else:
                row["selected_layer"] = "X"
                row["layer_decision"] = "processed_or_unknown_X"
        else:
            row["selected_layer"] = ""
            row["layer_decision"] = "processed_or_unknown_X"
        var_names = h5_read_string_dataset(var_group, "_index") if var_group is not None else []
        # v6.2 repair: fallback to scanpy-backed read if h5py read returns empty
        if not var_names and ad is not None:
            try:
                a_tmp = ad.read_h5ad(path, backed="r")
                var_names = list(a_tmp.var_names)
                a_tmp.file.close()
            except Exception:
                pass
        upper, gene_id_field, gene_symbol_field = resolve_gene_symbols(path, var_names)
        row["gene_id_field"] = gene_id_field
        row["gene_symbol_field"] = gene_symbol_field
        row["n_gene_symbols"] = len(upper)
        row["n_duplicated_gene_symbols"] = len(upper) - len(set(upper))
        row["n_caution_genes"] = sum(1 for g in upper if is_caution_gene(g))
        immune_present = sorted(set(upper) & IMMUNE_CORE)
        row["immune_core_covered"] = len(immune_present)
        row["immune_core_genes_present"] = "|".join(immune_present)
        state_support = {}
        u = set(upper)
        for state, genes in CELL_STATE_MARKERS.items():
            state_support[state] = len(u & genes)
        row["major_cell_state_support"] = json.dumps(state_support, sort_keys=True)
        hvg_cols = [c for c in var_cols if "highly_variable" in c.lower() or c.lower() in {"hvg", "variable"}]
        row["hvg_source"] = "|".join(hvg_cols[:8])
        hvg_genes = []
        for col in hvg_cols[:2]:
            try:
                vals = f["var"][col][:]
                picked = [g for g, flag in zip(var_names, vals) if bool(flag)]
                hvg_genes.extend(picked[:2000])
            except Exception:
                continue
        row["hvg_genes"] = "|".join(dict.fromkeys(hvg_genes[:2000]))
        try:
            f.close()
        except Exception:
            pass
    except Exception as exc:
        row["object_read_status"] = "read_failed"
        row["read_error"] = f"{type(exc).__name__}: {exc}"
    return row


def inspect_pseudobulk(path: Path, cohort_id: str) -> dict[str, Any]:
    row = {
        "cohort_id": cohort_id,
        "object_path": norm_path(path),
        "object_type": "pseudobulk",
        "object_read_status": "not_read",
        "read_error": "",
        "n_cells": "",
        "n_genes": "",
        "obs_columns": "",
        "var_columns": "",
        "layers": "",
        "raw_available": "no",
        "selected_layer": "pseudobulk_table",
        "layer_decision": "pseudobulk_only",
        "gene_id_field": "",
        "gene_symbol_field": "",
        "sample_column": "",
        "patient_column": "",
        "annotation_columns": "",
        "qc_columns": "",
        "n_samples_in_object": "",
        "n_patients_in_object": "",
        "n_gene_symbols": "",
        "n_duplicated_gene_symbols": "",
        "n_caution_genes": "",
        "immune_core_covered": "",
        "immune_core_total": len(IMMUNE_CORE),
        "major_cell_state_support": "",
        "hvg_source": "",
        "hvg_genes": "",
    }
    try:
        opener = gzip.open if path.suffix == ".gz" else open
        with opener(path, "rt", encoding="utf-8", errors="replace") as f:
            header = f.readline().strip().split("\t")
            n = 0
            examples = []
            for line in f:
                if n < 10000:
                    examples.append(line.split("\t", 1)[0])
                n += 1
        row["object_read_status"] = "readable"
        row["obs_columns"] = "|".join(header[:80])
        row["n_genes"] = n if n else max(len(header) - 1, 0)
        genes = [upper_gene(x) for x in examples]
        row["n_caution_genes"] = sum(1 for g in genes if is_caution_gene(g))
        row["immune_core_covered"] = len(set(genes) & IMMUNE_CORE)
        row["n_gene_symbols"] = row["n_genes"]
        row["sample_column"] = infer_sample_column(header)
        row["patient_column"] = infer_patient_column(header)
        row["annotation_columns"] = "|".join(find_annotation_columns(header))
    except Exception as exc:
        row["object_read_status"] = "read_failed"
        row["read_error"] = f"{type(exc).__name__}: {exc}"
    return row


def inspect_seurat_rds(path: Path, cohort_id: str) -> dict[str, Any]:
    row = {
        "cohort_id": cohort_id,
        "object_path": norm_path(path),
        "object_type": "seurat_rds",
        "object_read_status": "readable_registered_not_loaded",
        "read_error": "",
        "n_cells": "",
        "n_genes": "",
        "obs_columns": "",
        "var_columns": "",
        "layers": "Seurat_RDS_assays_not_loaded",
        "raw_available": "unknown",
        "selected_layer": "Seurat_RDS_registered_object",
        "layer_decision": "seurat_rds_requires_conversion_or_R_loader",
        "gene_id_field": "",
        "gene_symbol_field": "",
        "sample_column": "",
        "patient_column": "",
        "annotation_columns": "",
        "qc_columns": "",
        "n_samples_in_object": "",
        "n_patients_in_object": "",
        "n_gene_symbols": "",
        "n_duplicated_gene_symbols": "",
        "n_caution_genes": "",
        "immune_core_covered": "",
        "immune_core_total": len(IMMUNE_CORE),
        "major_cell_state_support": "",
        "hvg_source": "",
        "hvg_genes": "",
    }
    # Avoid loading multi-GB Seurat objects in the manifest gate. Use local metadata
    # extracted during intake to make the object traceable without rewriting it.
    meta = path.with_name("GSE301741_metadata_key_fields.tsv")
    patient = path.with_name("GSE301741_patient_response_metadata.tsv")
    sample = path.with_name("GSE301741_sample_manifest.tsv")
    try:
        m = re.search(r"_(\d+)cells", path.name)
        if m:
            row["n_cells"] = int(m.group(1))
        if meta.exists():
            with meta.open(newline="", encoding="utf-8") as f:
                reader = csv.reader(f, delimiter="\t")
                header = next(reader)
                n_meta_rows = sum(1 for _ in reader)
            row["obs_columns"] = "|".join(header)
            row["sample_column"] = infer_sample_column(header)
            row["patient_column"] = infer_patient_column(header)
            row["annotation_columns"] = "|".join(find_annotation_columns(header))
            row["qc_columns"] = json.dumps(detect_qc_columns(header), sort_keys=True)
            row["n_samples_in_object"] = ""
            row["n_patients_in_object"] = ""
            row["read_error"] = f"rds_not_loaded_metadata_rows={n_meta_rows}"
        if patient.exists():
            with patient.open(newline="", encoding="utf-8") as f:
                n_patients = sum(1 for _ in csv.DictReader(f, delimiter="\t"))
            row["n_patients_in_object"] = n_patients
        if sample.exists():
            with sample.open(newline="", encoding="utf-8") as f:
                n_samples = sum(1 for _ in csv.DictReader(f, delimiter="\t"))
            row["n_samples_in_object"] = n_samples
    except Exception as exc:
        row["object_read_status"] = "metadata_read_failed"
        row["read_error"] = f"{type(exc).__name__}: {exc}"
    return row


def build_universe(sample_rows: list[dict[str, str]], role_rows: list[dict[str, str]]) -> set[str]:
    role_by = {r["cohort_id"]: r for r in role_rows}
    universe = set()
    for r in role_rows:
        ff = r.get("feature_family", "")
        if any(x in ff for x in ("scRNA", "tcr", "TCR", "pseudobulk", "normalized_expression")):
            if r.get("feature_available") == "yes" and r.get("dataset_role") != "excluded":
                universe.add(r["cohort_id"])
    for r in sample_rows:
        mod = r.get("modality", "")
        ff = r.get("feature_family_availability", "")
        if "scRNA" in mod or "scTCR" in mod or "TCR" in ff:
            rid = role_by.get(r["cohort_id"], {})
            if rid.get("dataset_role") != "excluded":
                universe.add(r["cohort_id"])
    return universe


def discover_objects(universe: set[str]) -> dict[str, list[Path]]:
    objects: dict[str, list[Path]] = defaultdict(list)
    if H5AD_INVENTORY.exists():
        for r in read_csv(H5AD_INVENTORY):
            cid = r.get("cohort_id", "")
            if cid in universe:
                objects[cid].append(Path(r["file_path"]))
    for path in sorted(PHASE25.glob("**/*.h5ad")):
        name = path.name
        m = re.search(r"(GSE\d+|TASK\d+|PRJCA\d+)", name)
        if m and m.group(1) in universe:
            objects[m.group(1)].append(path)
        elif "GSE291246" in str(path) and "GSE291246" in universe:
            objects["GSE291246"].append(path)
    for path in sorted(PHASE25.glob("**/pseudobulk_addendum.*.tsv.gz")):
        name = path.name
        m = re.search(r"pseudobulk_addendum\.(.+?)\.tsv\.gz", name)
        if m:
            cid = m.group(1).replace("_rescue", "")
            if cid in universe:
                objects[cid].append(path)
    if EXTERNAL.exists():
        for path in sorted(EXTERNAL.glob("**/*.h5ad")):
            for cid in universe:
                if cid.lower() in str(path).lower():
                    objects[cid].append(path)
                    break
        for path in sorted(EXTERNAL.glob("**/*.rds")):
            for cid in universe:
                if cid.lower() in str(path).lower():
                    objects[cid].append(path)
                    break
    for cid in list(objects):
        seen = []
        dedup = []
        for p in objects[cid]:
            s = str(p)
            if s not in seen:
                seen.append(s)
                dedup.append(p)
        objects[cid] = dedup
    return objects


def phase4a_status(inspected: list[dict[str, Any]], cohort_id: str, sample_count: int) -> tuple[str, str]:
    readable = [r for r in inspected if str(r.get("object_read_status", "")).startswith("readable")]
    if not readable:
        return "blocked_or_support_only", "no_readable_expression_object"
    decisions = {r.get("layer_decision") for r in readable}
    max_immune = max(int(r.get("immune_core_covered") or 0) for r in readable)
    max_cells = max(int(r.get("n_cells") or 0) for r in readable if str(r.get("n_cells") or "").isdigit()) if any(str(r.get("n_cells") or "").isdigit() for r in readable) else 0
    if decisions == {"pseudobulk_only"}:
        return "pseudobulk_only", "rescue_or_pseudobulk_input"
    if decisions == {"seurat_rds_requires_conversion_or_R_loader"}:
        return "reference_mapping_only", "seurat_rds_registered_requires_conversion_to_h5ad_or_R_loader"
    if max_immune < 10:
        return "reference_mapping_only", "low_immune_core_gene_coverage"
    if max_cells and max_cells < 500:
        return "coarse_integration", "low_cell_count"
    if sample_count < 2:
        return "reference_mapping_only", "single_sample_or_sample_mapping_limited"
    if any(d in decisions for d in ("raw_counts", "raw_or_author_raw", "raw_counts_in_X")):
        return "full_integration_candidate", "readable_count_like_layer_and_gene_coverage"
    return "coarse_integration", "processed_or_unknown_X_layer"


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out-dir", default=str(OUT_DIR))
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--limit-cohorts", default="")
    args = parser.parse_args()

    out_dir = Path(args.out_dir)
    sample_rows = read_csv(FROZEN_SAMPLE)
    patient_rows = read_csv(FROZEN_PATIENT)
    role_rows = read_csv(FROZEN_ROLE)
    response_rows = read_csv(FROZEN_RESPONSE)
    label_rows = read_csv(FROZEN_LABEL_INV)
    universe = build_universe(sample_rows, role_rows)
    if args.limit_cohorts:
        wanted = {x.strip() for x in args.limit_cohorts.split(",") if x.strip()}
        universe &= wanted
    objects = discover_objects(universe)

    samples_by_cohort = defaultdict(list)
    sample_keys_by_cohort = defaultdict(set)
    patient_keys_by_cohort = defaultdict(set)
    for row in sample_rows:
        cid = row["cohort_id"]
        if cid in universe:
            samples_by_cohort[cid].append(row)
            sample_keys_by_cohort[cid].add(row["sample_key"])
            patient_keys_by_cohort[cid].add(row["patient_key"])
    role_by = {r["cohort_id"]: r for r in role_rows}
    response_by = {r["cohort_id"]: r for r in response_rows}

    inspected_by_cohort: dict[str, list[dict[str, Any]]] = defaultdict(list)
    processing_rows = []
    object_registry_rows = []
    layer_rows = []
    gene_rows = []
    barcode_rows = []
    qc_rows = []
    filtering_rows = []
    anno_rows = []
    batch_rows = []
    tcr_rows = []
    all_hvg = Counter()
    gene_presence = Counter()
    caution_rows = []
    hard_blockers = []
    soft_warnings = []

    for cid in sorted(universe):
        obj_paths = objects.get(cid, [])
        if not obj_paths:
            obj_paths = []
            soft_warnings.append(f"{cid}: no local expression object discovered")
        for idx, path in enumerate(obj_paths or [Path("__missing__")]):
            exists = path.exists() if str(path) != "__missing__" else False
            hash_mode, hash_value, file_size, mtime = file_fingerprint(path) if exists else ("missing", "", 0, "")
            if not exists:
                info = {
                    "cohort_id": cid,
                    "object_path": "",
                    "object_type": "missing",
                    "object_read_status": "missing",
                    "read_error": "no object discovered",
                    "n_cells": "",
                    "n_genes": "",
                    "obs_columns": "",
                    "var_columns": "",
                    "layers": "",
                    "raw_available": "unknown",
                    "selected_layer": "",
                    "layer_decision": "support_only",
                    "gene_id_field": "",
                    "gene_symbol_field": "",
                    "sample_column": "",
                    "patient_column": "",
                    "annotation_columns": "",
                    "qc_columns": "",
                    "n_samples_in_object": "",
                    "n_patients_in_object": "",
                    "n_gene_symbols": "",
                    "n_duplicated_gene_symbols": "",
                    "n_caution_genes": "",
                    "immune_core_covered": "",
                    "immune_core_total": len(IMMUNE_CORE),
                    "major_cell_state_support": "",
                    "hvg_source": "",
                    "hvg_genes": "",
                }
            elif path.suffix.lower() == ".rds":
                info = inspect_seurat_rds(path, cid)
            elif path.suffix == ".gz" or "pseudobulk" in path.name:
                info = inspect_pseudobulk(path, cid)
            else:
                info = inspect_h5ad(path, cid)
            inspected_by_cohort[cid].append(info)
            role = role_by.get(cid, {})
            status, reason = phase4a_status([info], cid, len(samples_by_cohort[cid]))
            obj_id = f"{cid}::object_{idx + 1}"
            common = {
                "cohort_id": cid,
                "object_id": obj_id,
                "object_path": info.get("object_path", norm_path(path) if exists else ""),
                "object_type": info.get("object_type", ""),
                "object_read_status": info.get("object_read_status", ""),
                "file_hash_mode": hash_mode,
                "file_hash": hash_value,
                "file_size_bytes": file_size,
                "file_mtime_epoch": mtime,
                "input_version": INPUT_VERSION,
            }
            processing_rows.append({
                **common,
                "dataset_role": role.get("dataset_role", ""),
                "feature_family": role.get("feature_family", ""),
                "metadata_source": "frozen_v0",
                "assay_or_layer_summary": info.get("layers", ""),
                "sample_metadata_aligned": "yes" if samples_by_cohort[cid] else "no",
                "patient_metadata_aligned": "yes" if patient_keys_by_cohort[cid] else "no",
                "phase3_5_status": status,
                "downgrade_reason": reason,
            })
            object_registry_rows.append({
                **common,
                "n_cells": info.get("n_cells", ""),
                "n_genes": info.get("n_genes", ""),
                "selected_layer": info.get("selected_layer", ""),
                "layer_decision": info.get("layer_decision", ""),
                "raw_available": info.get("raw_available", ""),
                "sample_column": info.get("sample_column", ""),
                "patient_column": info.get("patient_column", ""),
                "annotation_columns": info.get("annotation_columns", ""),
                "gene_id_field": info.get("gene_id_field", ""),
                "gene_symbol_field": info.get("gene_symbol_field", ""),
                "phase4a_eligibility": status,
                "phase4b_eligibility": "yes" if status in {"full_integration_candidate", "coarse_integration"} and samples_by_cohort[cid] else "no",
                "downgrade_reason": reason,
            })
            layer_rows.append({
                "cohort_id": cid,
                "object_id": obj_id,
                "object_type": info.get("object_type", ""),
                "available_layers": info.get("layers", ""),
                "raw_available": info.get("raw_available", ""),
                "selected_layer": info.get("selected_layer", ""),
                "layer_decision": info.get("layer_decision", ""),
                "raw_counts_allowed": "yes" if info.get("layer_decision") in {"raw_counts", "raw_or_author_raw", "raw_counts_in_X"} else "no",
                "normalized_expression_allowed": "yes" if info.get("layer_decision") in {"processed_or_unknown_X", "pseudobulk_only"} else "no",
                "pseudobulk_allowed": "yes" if info.get("layer_decision") == "pseudobulk_only" else "no",
                "support_only_reason": "" if status != "blocked_or_support_only" else reason,
                "input_version": INPUT_VERSION,
            })
            gene_rows.append({
                "cohort_id": cid,
                "object_id": obj_id,
                "n_genes": info.get("n_genes", ""),
                "n_gene_symbols": info.get("n_gene_symbols", ""),
                "n_duplicated_gene_symbols": info.get("n_duplicated_gene_symbols", ""),
                "n_caution_genes": info.get("n_caution_genes", ""),
                "immune_core_covered": info.get("immune_core_covered", ""),
                "immune_core_total": len(IMMUNE_CORE),
                "immune_core_coverage_fraction": (
                    round(int(info.get("immune_core_covered") or 0) / len(IMMUNE_CORE), 4)
                    if str(info.get("immune_core_covered") or "").isdigit() else ""
                ),
                "major_cell_state_marker_support_json": info.get("major_cell_state_support", ""),
                "gene_id_field": info.get("gene_id_field", ""),
                "gene_symbol_field": info.get("gene_symbol_field", ""),
                "gene_id_harmonization_decision": (
                    "record_only_no_gene_rewrite" if info.get("gene_symbol_field") == "var_names"
                    else "repaired_symbol_mapping_in_memory" if info.get("gene_symbol_field")
                    else "record_only_no_gene_rewrite"
                ),
                "downgrade_reason": "low_immune_core_gene_coverage" if str(info.get("immune_core_covered") or "").isdigit() and int(info.get("immune_core_covered") or 0) < 10 else "",
                "input_version": INPUT_VERSION,
            })
            if info.get("hvg_genes"):
                for g in str(info["hvg_genes"]).split("|"):
                    if g:
                        all_hvg[g] += 1
            if info.get("immune_core_genes_present"):
                for g in str(info["immune_core_genes_present"]).split("|"):
                    if g:
                        gene_presence[g] += 1
            caution_rows.extend([
                {"gene": g, "reason": "mitochondrial_ribosomal_hb_stress_cell_cycle", "input_version": INPUT_VERSION}
                for g in list(CAUTION_GENES)
            ])
            if info.get("object_read_status") in {"read_failed", "metadata_read_failed"}:
                soft_warnings.append(f"{cid}: read_failed {info.get('read_error')}")

        # Cohort-level rows derived after all objects are inspected.
        status, reason = phase4a_status(inspected_by_cohort[cid], cid, len(samples_by_cohort[cid]))
        barcode_rows.append({
            "cohort_id": cid,
            "sample_key_count_frozen": len(sample_keys_by_cohort[cid]),
            "patient_key_count_frozen": len(patient_keys_by_cohort[cid]),
            "sample_column_candidates": "|".join(sorted({x.get("sample_column", "") for x in inspected_by_cohort[cid] if x.get("sample_column")})),
            "patient_column_candidates": "|".join(sorted({x.get("patient_column", "") for x in inspected_by_cohort[cid] if x.get("patient_column")})),
            "cell_barcode_context_required": "cohort_id+sample_id+cell_barcode",
            "join_to_frozen_v0_status": "aligned_by_frozen_sample_rows" if samples_by_cohort[cid] else "no_frozen_sample_rows",
            "phase4b_patient_feature_allowed": "yes" if status in {"full_integration_candidate", "coarse_integration"} and samples_by_cohort[cid] else "no",
            "downgrade_reason": reason if status not in {"full_integration_candidate", "coarse_integration"} else "",
            "input_version": INPUT_VERSION,
        })
        qc_rows.append({
            "cohort_id": cid,
            "n_frozen_samples": len(samples_by_cohort[cid]),
            "n_frozen_patients": len(patient_keys_by_cohort[cid]),
            "n_readable_objects": sum(1 for x in inspected_by_cohort[cid] if str(x.get("object_read_status", "")).startswith("readable")),
            "n_cells_total_object_reported": sum(int(x.get("n_cells") or 0) for x in inspected_by_cohort[cid] if str(x.get("n_cells") or "").isdigit()),
            "n_genes_max_object_reported": max([int(x.get("n_genes") or 0) for x in inspected_by_cohort[cid] if str(x.get("n_genes") or "").isdigit()] or [0]),
            "qc_metric_source": "existing_obs_metrics_or_object_level_metadata",
            "mitochondrial_fraction_status": "available" if any("pct_mito" in str(x.get("qc_columns")) for x in inspected_by_cohort[cid]) else "not_available",
            "ribosomal_fraction_status": "available" if any("pct_ribo" in str(x.get("qc_columns")) for x in inspected_by_cohort[cid]) else "not_available",
            "hemoglobin_fraction_status": "available" if any("pct_hb" in str(x.get("qc_columns")) for x in inspected_by_cohort[cid]) else "not_available",
            "doublet_status": "available" if any("doublet" in str(x.get("qc_columns")).lower() for x in inspected_by_cohort[cid]) else "not_available",
            "low_quality_fraction_status": "requires_phase4a_or_object_specific_filtering" if status != "blocked_or_support_only" else "not_estimated",
            "input_version": INPUT_VERSION,
        })
        filtering_rows.append({
            "cohort_id": cid,
            "filtering_action": "recommend_only_no_object_rewrite",
            "pre_filter_cells": qc_rows[-1]["n_cells_total_object_reported"],
            "post_filter_cells_estimate": "",
            "thresholds": "sample_adaptive_MAD_plus_hard_mito_gene_count_thresholds_recommended",
            "uses_author_qc": "yes_if_object_already_processed",
            "project_qc_added": "manifest_level_only",
            "downgrade_decision": status,
            "downgrade_reason": reason,
            "input_version": INPUT_VERSION,
        })
        anno_cols = sorted({c for x in inspected_by_cohort[cid] for c in str(x.get("annotation_columns", "")).split("|") if c})
        anno_rows.append({
            "cohort_id": cid,
            "annotation_columns": "|".join(anno_cols),
            "coarse_label_available": "yes" if anno_cols else "unknown",
            "mid_label_available": "unknown",
            "fine_label_available": "unknown",
            "annotation_source": "object_obs" if anno_cols else "not_found",
            "confidence": "medium" if anno_cols else "low",
            "phase4a_remap_required": "yes",
            "fine_grained_claim_allowed": "no",
            "input_version": INPUT_VERSION,
        })
        batch_rows.append({
            "cohort_id": cid,
            "allowed_integration_covariates": "cohort_id|platform|sample_key|patient_key|cancer_type|tissue_source|timepoint|treatment_context|object_id",
            "blocked_covariates": "response_raw|response_binary_harmonized|split|survival|outcome",
            "available_sample_covariates": "tissue_source|timepoint|treatment_context|dataset_role",
            "response_label_use": "evaluation_or_leakage_audit_only",
            "input_version": INPUT_VERSION,
        })
        tcr_available = role_by.get(cid, {}).get("tcr_available", "no")
        tcr_rows.append({
            "cohort_id": cid,
            "tcr_available_frozen": tcr_available,
            "gex_tcr_join_key": "cell_barcode+sample_context_required" if tcr_available == "yes" else "",
            "join_rate": "",
            "matched_cell_count": "",
            "unmatched_tcr_count": "",
            "unmatched_gex_count": "",
            "clonotype_fields_available": "unknown" if tcr_available == "yes" else "no",
            "tcr_feature_mainline_allowed": "no",
            "downgrade_reason": "join_rate_not_established_in_phase3_5" if tcr_available == "yes" else "no_tcr",
            "input_version": INPUT_VERSION,
        })

    hvg_global = []
    for gene, n in all_hvg.most_common(5000):
        if not is_caution_gene(gene):
            hvg_global.append({"gene": gene, "supporting_objects": n, "selection_rule": "existing_object_hvg_response_blind_reuse", "input_version": INPUT_VERSION})
    if not hvg_global:
        for gene, n in gene_presence.most_common():
            if not is_caution_gene(gene):
                hvg_global.append({"gene": gene, "supporting_objects": n, "selection_rule": "immune_core_presence_fallback_response_blind", "input_version": INPUT_VERSION})
    hvg_by_state = []
    for state, genes in CELL_STATE_MARKERS.items():
        for gene in sorted(genes):
            hvg_by_state.append({
                "major_cell_state": state,
                "gene": gene,
                "selection_rule": "curated_marker_fallback_response_blind_until_phase4a_annotation",
                "input_version": INPUT_VERSION,
            })
    caution_unique = sorted({r["gene"] for r in caution_rows} | {g for g in IMMUNE_CORE if is_caution_gene(g)})
    caution_rows = [{"gene": g, "reason": "mitochondrial_ribosomal_hb_stress_cell_cycle_or_platform_dominated_risk", "input_version": INPUT_VERSION} for g in caution_unique]

    readiness_rows = []
    for cid in sorted(universe):
        status, reason = phase4a_status(inspected_by_cohort[cid], cid, len(samples_by_cohort[cid]))
        readiness_rows.append({
            "cohort_id": cid,
            "phase4a_readiness": status,
            "phase4b_readiness": "patient_timepoint_features_allowed" if status in {"full_integration_candidate", "coarse_integration"} and samples_by_cohort[cid] else "not_allowed",
        "readable_objects": sum(1 for x in inspected_by_cohort[cid] if str(x.get("object_read_status", "")).startswith("readable")),
            "frozen_samples": len(samples_by_cohort[cid]),
            "frozen_patients": len(patient_keys_by_cohort[cid]),
            "downgrade_reason": reason,
            "input_version": INPUT_VERSION,
        })

    out_dir.mkdir(parents=True, exist_ok=True)
    write_csv(out_dir / "sc_processing_manifest.frozen_v0.csv", processing_rows, [
        "cohort_id", "object_id", "object_path", "object_type", "object_read_status", "file_hash_mode",
        "file_hash", "file_size_bytes", "file_mtime_epoch", "dataset_role", "feature_family",
        "metadata_source", "assay_or_layer_summary", "sample_metadata_aligned", "patient_metadata_aligned",
        "phase3_5_status", "downgrade_reason", "input_version",
    ])
    write_csv(out_dir / "analysis_object_registry.frozen_v0.csv", object_registry_rows, [
        "cohort_id", "object_id", "object_path", "object_type", "object_read_status", "file_hash_mode",
        "file_hash", "file_size_bytes", "file_mtime_epoch", "n_cells", "n_genes", "selected_layer",
        "layer_decision", "raw_available", "sample_column", "patient_column", "annotation_columns",
        "gene_id_field", "gene_symbol_field", "phase4a_eligibility", "phase4b_eligibility",
        "downgrade_reason", "input_version",
    ])
    write_csv(out_dir / "matrix_layer_decision_table.frozen_v0.csv", layer_rows, [
        "cohort_id", "object_id", "object_type", "available_layers", "raw_available", "selected_layer",
        "layer_decision", "raw_counts_allowed", "normalized_expression_allowed", "pseudobulk_allowed",
        "support_only_reason", "input_version",
    ])
    write_csv(out_dir / "gene_overlap_and_immune_core_coverage.csv", gene_rows, [
        "cohort_id", "object_id", "n_genes", "n_gene_symbols", "n_duplicated_gene_symbols",
        "n_caution_genes", "immune_core_covered", "immune_core_total", "immune_core_coverage_fraction",
        "major_cell_state_marker_support_json", "gene_id_field", "gene_symbol_field",
        "gene_id_harmonization_decision", "downgrade_reason", "input_version",
    ])
    write_csv(out_dir / "sample_cell_barcode_index.csv", barcode_rows, [
        "cohort_id", "sample_key_count_frozen", "patient_key_count_frozen", "sample_column_candidates",
        "patient_column_candidates", "cell_barcode_context_required", "join_to_frozen_v0_status",
        "phase4b_patient_feature_allowed", "downgrade_reason", "input_version",
    ])
    write_csv(out_dir / "cell_qc_summary_by_sample.csv", qc_rows, [
        "cohort_id", "n_frozen_samples", "n_frozen_patients", "n_readable_objects",
        "n_cells_total_object_reported", "n_genes_max_object_reported", "qc_metric_source",
        "mitochondrial_fraction_status", "ribosomal_fraction_status", "hemoglobin_fraction_status",
        "doublet_status", "low_quality_fraction_status", "input_version",
    ])
    write_csv(out_dir / "cell_filtering_decision_log.csv", filtering_rows, [
        "cohort_id", "filtering_action", "pre_filter_cells", "post_filter_cells_estimate", "thresholds",
        "uses_author_qc", "project_qc_added", "downgrade_decision", "downgrade_reason", "input_version",
    ])
    write_csv(out_dir / "gex_tcr_join_qc_report.csv", tcr_rows, [
        "cohort_id", "tcr_available_frozen", "gex_tcr_join_key", "join_rate", "matched_cell_count",
        "unmatched_tcr_count", "unmatched_gex_count", "clonotype_fields_available",
        "tcr_feature_mainline_allowed", "downgrade_reason", "input_version",
    ])
    write_csv(out_dir / "cell_annotation_source_inventory.csv", anno_rows, [
        "cohort_id", "annotation_columns", "coarse_label_available", "mid_label_available",
        "fine_label_available", "annotation_source", "confidence", "phase4a_remap_required",
        "fine_grained_claim_allowed", "input_version",
    ])
    write_tsv(out_dir / "hvg_global_pan_cancer.tsv", hvg_global[:5000], [
        "gene", "supporting_objects", "selection_rule", "input_version",
    ])
    write_tsv(out_dir / "hvg_by_major_cell_state.tsv", hvg_by_state, [
        "major_cell_state", "gene", "selection_rule", "input_version",
    ])
    write_tsv(out_dir / "hvg_caution_or_exclusion_list.tsv", caution_rows, [
        "gene", "reason", "input_version",
    ])
    write_csv(out_dir / "batch_covariate_registry.csv", batch_rows, [
        "cohort_id", "allowed_integration_covariates", "blocked_covariates", "available_sample_covariates",
        "response_label_use", "input_version",
    ])

    readable_cohorts = sum(1 for r in readiness_rows if int(r["readable_objects"]) > 0)
    full_candidates = sum(1 for r in readiness_rows if r["phase4a_readiness"] == "full_integration_candidate")
    coarse = sum(1 for r in readiness_rows if r["phase4a_readiness"] == "coarse_integration")
    pseudobulk = sum(1 for r in readiness_rows if r["phase4a_readiness"] == "pseudobulk_only")
    blocked = [r for r in readiness_rows if r["phase4a_readiness"] == "blocked_or_support_only"]
    verdict = "CONDITIONAL_GO_TO_PHASE4A"
    if readable_cohorts == 0 or readable_cohorts < max(1, len(universe) // 2):
        verdict = "BLOCKED"
        hard_blockers.append("less_than_half_single_cell_universe_has_readable_objects")
    if H5PY_IMPORT_ERROR is not None:
        verdict = "BLOCKED"
        hard_blockers.append(f"h5py_import_error::{H5PY_IMPORT_ERROR}")

    handoff = {
        "phase": "v6.2 Phase3.5 Single-cell Processing QC Gate",
        "verdict": verdict,
        "input_version": INPUT_VERSION,
        "created_at": now_iso(),
        "execution_mode": "manifest_first_no_large_object_rewrite",
        "environment": {
            "recommended_command_prefix": "MPLCONFIGDIR=/tmp/mplcfg NUMBA_CACHE_DIR=/tmp/numba_cache conda run -n sc_inte python",
            "h5py_available": H5PY_IMPORT_ERROR is None,
            "anndata_available": ANNDATA_IMPORT_ERROR is None,
        },
        "outputs": {
            "sc_processing_manifest": "sc_processing_manifest.frozen_v0.csv",
            "matrix_layer_decision": "matrix_layer_decision_table.frozen_v0.csv",
            "gene_id_harmonization_report": "gene_id_harmonization_report.md",
            "gene_overlap_and_immune_core_coverage": "gene_overlap_and_immune_core_coverage.csv",
            "cell_qc_summary_by_sample": "cell_qc_summary_by_sample.csv",
            "cell_filtering_decision_log": "cell_filtering_decision_log.csv",
            "sample_cell_barcode_index": "sample_cell_barcode_index.csv",
            "gex_tcr_join_qc_report": "gex_tcr_join_qc_report.csv",
            "cell_annotation_source_inventory": "cell_annotation_source_inventory.csv",
            "hvg_global_pan_cancer": "hvg_global_pan_cancer.tsv",
            "hvg_by_major_cell_state": "hvg_by_major_cell_state.tsv",
            "hvg_caution_or_exclusion_list": "hvg_caution_or_exclusion_list.tsv",
            "batch_covariate_registry": "batch_covariate_registry.csv",
            "integration_readiness_report": "integration_readiness_report.md",
            "analysis_object_registry": "analysis_object_registry.frozen_v0.csv",
            "handoff": "phase3_5_to_phase4_handoff.yaml",
        },
        "qc_summary": {
            "n_single_cell_or_rescue_cohorts": len(universe),
            "n_readable_object_cohorts": readable_cohorts,
            "n_full_integration_candidates": full_candidates,
            "n_coarse_integration": coarse,
            "n_pseudobulk_only": pseudobulk,
            "n_blocked_or_support_only": len(blocked),
            "response_or_split_used_for_hvg": False,
            "response_or_split_used_as_integration_covariate": False,
            "large_analysis_objects_rewritten": False,
        },
        "hard_blockers": hard_blockers,
        "soft_warnings": soft_warnings[:200],
        "phase4a_inputs": [
            {"cohort_id": r["cohort_id"], "readiness": r["phase4a_readiness"], "downgrade_reason": r["downgrade_reason"]}
            for r in readiness_rows
        ],
        "rules": {
            "direct_raw_download_to_phase4_allowed": False,
            "manifest_first": True,
            "response_blind_hvg": True,
            "response_blind_qc": True,
            "write_large_standard_objects": False,
        },
        "next_phase_readiness": {
            "enter_phase4a": verdict != "BLOCKED",
            "phase4a_allowed_to_read_only_phase3_5_registry": True,
            "required_patch_before_phase4a": [
                "establish_object_specific_cell_filtering_if_full_integration_requires_actual_filtered_h5ad",
                "curate_TCR_join_keys_before_mainline_TCR_features",
                "resolve_missing_objects_for_external_candidates_if_clean_anchor_full_integration_is_required",
            ],
        },
    }
    with (out_dir / "phase3_5_to_phase4_handoff.yaml").open("w", encoding="utf-8") as f:
        yaml.safe_dump(handoff, f, sort_keys=False, allow_unicode=True)

    report = f"""# v6.2 Phase3.5 Integration Readiness Report

Generated: {now_iso()}
Input version: `{INPUT_VERSION}`
Verdict: `{verdict}`
Execution mode: manifest-first, no large h5ad/zarr rewrite.

## Summary

- Single-cell / rescue cohorts in scope: {len(universe)}
- Cohorts with readable local objects: {readable_cohorts}
- Full integration candidates: {full_candidates}
- Coarse integration candidates: {coarse}
- Pseudobulk-only inputs: {pseudobulk}
- Blocked or support-only at object gate: {len(blocked)}

## Gate Rules

- Phase4A must read `analysis_object_registry.frozen_v0.csv` and `phase3_5_to_phase4_handoff.yaml`.
- Response, outcome, survival, and split are not allowed as QC, HVG, or integration covariates.
- Large analysis objects were not rewritten in Phase3.5.
- Cohorts without readable object, clear layer, gene coverage, or frozen sample/patient alignment are downgraded instead of silently dropped.

## Main Downgrades

"""
    for r in blocked[:80]:
        report += f"- `{r['cohort_id']}`: {r['downgrade_reason']}\n"
    if not blocked:
        report += "- None at cohort-level hard gate.\n"
    report += """
## Environment Note

Default base Python has NumPy ABI issues. Use:

```bash
MPLCONFIGDIR=/tmp/mplcfg NUMBA_CACHE_DIR=/tmp/numba_cache conda run -n sc_inte python scripts/v6_2/build_phase3_5_single_cell_qc_gate.py
```

## Phase4A Contract

Phase4A may use full/coarse/reference/pseudobulk eligible objects only through the registry and handoff. It must not directly scan raw download directories.
"""
    (out_dir / "integration_readiness_report.md").write_text(report, encoding="utf-8")

    gene_report = f"""# Gene ID Harmonization Report

Input version: `{INPUT_VERSION}`
Policy: record-only; Phase3.5 does not rewrite gene identifiers inside h5ad/pseudobulk objects.

## Decisions

- Gene symbols are read from `var_names` unless object-specific metadata exposes a clearer field.
- Duplicate gene symbols are counted and reported in `gene_overlap_and_immune_core_coverage.csv`.
- Mitochondrial, ribosomal, hemoglobin, stress, and cell-cycle genes are recorded in `hvg_caution_or_exclusion_list.tsv`.
- Immune-core coverage is reported per object and used for integration readiness downgrades.
- Any actual gene-name rewriting must be a later explicit object-building patch, not this manifest-first gate.

## Scope

- Cohorts in scope: {len(universe)}
- Object manifest rows: {len(processing_rows)}
- Gene coverage rows: {len(gene_rows)}
"""
    (out_dir / "gene_id_harmonization_report.md").write_text(gene_report, encoding="utf-8")

    print(json.dumps({
        "out_dir": str(out_dir),
        "verdict": verdict,
        "cohorts": len(universe),
        "readable_object_cohorts": readable_cohorts,
        "outputs": 16,
    }, indent=2))
    return 0 if verdict != "BLOCKED" or args.dry_run else 2


if __name__ == "__main__":
    raise SystemExit(main())
