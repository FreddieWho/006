#!/usr/bin/env python3
"""Shared utilities for Phase6 module algorithm benchmark."""
from __future__ import annotations

import gzip
import hashlib
import json
import os
import re
import subprocess
from datetime import date
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

ROOT = Path(__file__).resolve().parents[2]
PHASE6 = ROOT / "results/v6_2/phase6_response_blind_module_discovery"
_out_override = os.environ.get("PHASE6_ALGORITHM_BENCHMARK_OUT")
OUT = (ROOT / _out_override) if _out_override else (ROOT / "results/v6_2/phase6_module_algorithm_benchmark")
TODAY = str(date.today())
SEED = 1729
TOP_N = 40
EPS = 1e-9

META_COLUMNS = {
    "expression_unit_id",
    "cohort_id",
    "object_id",
    "sample_key",
    "patient_key",
    "cell_state_level",
    "cell_state",
    "layer_family",
    "n_cells_used",
    "library_size",
}
BANNED_TOKENS = [
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


def ensure_dirs() -> None:
    for sub in ["inputs", "methods", "score", "report"]:
        (OUT / sub).mkdir(parents=True, exist_ok=True)


def write_csv(df: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, index=False)


def write_yaml(obj: dict, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as fh:
        yaml.safe_dump(obj, fh, sort_keys=False, allow_unicode=True)


def write_md(text: str, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")


def safe_id(value: object) -> str:
    return re.sub(r"[^A-Za-z0-9]+", "_", str(value)).strip("_") or "unknown"


def sha256_path(path: Path, block: int = 4_194_304) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        while True:
            data = fh.read(block)
            if not data:
                break
            h.update(data)
    return h.hexdigest()


def gene_columns(df: pd.DataFrame) -> list[str]:
    return [c for c in df.columns if c not in META_COLUMNS]


def load_primary_tables() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, list[str]]:
    raw = pd.read_parquet(PHASE6 / "foundation/raw_count_cellstate_pseudobulk.v0.parquet")
    log = pd.read_parquet(PHASE6 / "foundation/logcpm_cellstate_pseudobulk.v0.parquet")
    registry = pd.read_csv(PHASE6 / "foundation/unified_expression_unit_registry.csv", dtype=str, keep_default_na=False)
    primary_ids = set(registry.loc[registry["inclusion_status"].eq("primary"), "expression_unit_id"].astype(str))
    raw = raw[raw["expression_unit_id"].astype(str).isin(primary_ids)].reset_index(drop=True)
    log = log[log["expression_unit_id"].astype(str).isin(primary_ids)].reset_index(drop=True)
    raw_genes = set(gene_columns(raw))
    log_genes = set(gene_columns(log))
    genes = sorted(raw_genes & log_genes)
    meta_cols = [c for c in raw.columns if c in META_COLUMNS]
    meta = raw[meta_cols].copy()
    raw = raw[meta_cols + genes].copy()
    log = log[[c for c in log.columns if c in META_COLUMNS] + genes].copy()
    return meta, raw, log, genes


def numeric_matrix(df: pd.DataFrame, genes: list[str]) -> np.ndarray:
    return df[genes].apply(pd.to_numeric, errors="coerce").fillna(0.0).to_numpy(dtype=float)


def balanced_indices(meta: pd.DataFrame, max_rows: int, seed: int = SEED) -> np.ndarray:
    rng = np.random.default_rng(seed)
    if len(meta) <= max_rows:
        return np.arange(len(meta))
    idx: list[int] = []
    for _, group in meta.groupby(["cohort_id", "cell_state"], sort=False):
        take = max(1, int(max_rows * len(group) / len(meta)))
        values = group.index.to_numpy()
        idx.extend(rng.choice(values, size=min(take, len(values)), replace=False).tolist())
    idx_arr = np.asarray(sorted(set(idx)), dtype=int)
    if len(idx_arr) > max_rows:
        idx_arr = rng.choice(idx_arr, size=max_rows, replace=False)
    return np.sort(idx_arr)


def top_gene_sets(components: np.ndarray, genes: list[str], top_n: int = TOP_N) -> list[set[str]]:
    out = []
    for row in components:
        weights = np.maximum(np.asarray(row, dtype=float), 0.0)
        if weights.max() <= 0:
            weights = np.abs(np.asarray(row, dtype=float))
        order = np.argsort(-weights)[: min(top_n, len(genes))]
        out.append({genes[i] for i in order if weights[i] > 0})
    return out


def max_jaccard(base: list[set[str]], other: list[set[str]]) -> list[float]:
    values = []
    for b in base:
        scores = []
        for o in other:
            union = len(b | o)
            scores.append(len(b & o) / union if union else 0.0)
        values.append(max(scores) if scores else 0.0)
    return values


def score_from_components(meta: pd.DataFrame, x: np.ndarray, components: np.ndarray) -> np.ndarray:
    weights = np.maximum(components.astype(float), 0.0)
    empty = weights.sum(axis=1) <= 0
    if empty.any():
        weights[empty] = np.abs(components[empty].astype(float))
    weights = weights / (weights.sum(axis=1, keepdims=True) + EPS)
    return x @ weights.T


def build_membership(method: str, run_id: str, components: np.ndarray, genes: list[str]) -> pd.DataFrame:
    rows = []
    for i, weights in enumerate(components):
        module_id = f"{method}__M{i + 1:02d}"
        pos = np.maximum(weights.astype(float), 0.0)
        if pos.max() <= 0:
            pos = np.abs(weights.astype(float))
        max_w = float(pos.max()) or 1.0
        order = set(np.argsort(-pos)[: min(TOP_N, len(genes))].tolist())
        for j, weight in enumerate(pos):
            if weight <= 0:
                continue
            rows.append(
                {
                    "method": method,
                    "run_id": run_id,
                    "module_id": module_id,
                    "gene": genes[j],
                    "membership_weight": float(weight),
                    "relative_weight": float(weight / max_w),
                    "is_top_gene": "yes" if j in order else "no",
                }
            )
    return pd.DataFrame(rows)


def build_score_matrix(method: str, run_id: str, meta: pd.DataFrame, scores: np.ndarray) -> pd.DataFrame:
    out = meta[
        ["expression_unit_id", "cohort_id", "object_id", "sample_key", "patient_key", "cell_state_level", "cell_state", "layer_family"]
    ].copy()
    for i in range(scores.shape[1]):
        out[f"{method}__M{i + 1:02d}"] = scores[:, i]
    out.insert(0, "run_id", run_id)
    out.insert(0, "method", method)
    return out


def leakage_columns(columns: list[str]) -> pd.DataFrame:
    rows = []
    for col in columns:
        lower = col.lower()
        hits = [t for t in BANNED_TOKENS if t in lower]
        rows.append({"column": col, "leakage_tokens": "|".join(hits), "leakage_risk": "yes" if hits else "no"})
    return pd.DataFrame(rows)


def write_method_outputs(
    method: str,
    run_id: str,
    membership: pd.DataFrame,
    scores: pd.DataFrame,
    run_qc: pd.DataFrame,
    manifest: dict,
) -> None:
    method_dir = OUT / "methods" / method
    method_dir.mkdir(parents=True, exist_ok=True)
    write_csv(membership, method_dir / f"module_membership.{method}.csv")
    scores.to_parquet(method_dir / f"module_score_matrix.{method}.parquet", index=False)
    write_csv(run_qc, method_dir / f"module_run_qc.{method}.csv")
    write_yaml(manifest | {"method": method, "run_id": run_id, "created_at": TODAY}, method_dir / f"module_method_manifest.{method}.yaml")


def write_blocked_method(method: str, reason: str, detail: str = "") -> None:
    method_dir = OUT / "methods" / method
    method_dir.mkdir(parents=True, exist_ok=True)
    run_qc = pd.DataFrame(
        [
            {
                "method": method,
                "run_id": f"{method}__blocked",
                "run_status": "blocked",
                "blocked_reason": reason,
                "detail": detail,
            }
        ]
    )
    write_csv(run_qc, method_dir / f"module_run_qc.{method}.csv")
    write_yaml(
        {
            "method": method,
            "run_id": f"{method}__blocked",
            "run_status": "blocked",
            "blocked_reason": reason,
            "detail": detail,
            "official_method_required": True,
        },
        method_dir / f"module_method_manifest.{method}.yaml",
    )


def run_cmd(cmd: list[str], cwd: Path | None = None, timeout: int | None = None) -> tuple[int, str, str]:
    proc = subprocess.run(cmd, cwd=str(cwd or ROOT), text=True, capture_output=True, timeout=timeout)
    return proc.returncode, proc.stdout, proc.stderr


def load_signatures() -> dict[str, set[str]]:
    path = ROOT / "mvp/outputs/program_signatures.json"
    if not path.exists():
        return {}
    data = json.loads(path.read_text())
    sigs: dict[str, set[str]] = {}
    if isinstance(data, dict):
        for key, value in data.items():
            if isinstance(value, list):
                sigs[str(key)] = {str(g).upper() for g in value}
            elif isinstance(value, dict):
                genes = value.get("genes") or value.get("signature") or []
                if not genes:
                    genes = list(value.get("up_genes") or []) + list(value.get("down_genes") or [])
                if isinstance(genes, list):
                    sigs[str(key)] = {str(g).upper() for g in genes}
    return sigs
