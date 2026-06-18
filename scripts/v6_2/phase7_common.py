#!/usr/bin/env python3
"""Shared, response-blind utilities for v6.2.1 Phase7."""

from __future__ import annotations

import hashlib
import json
import math
import os
import re
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd
import yaml


ROOT = Path(__file__).resolve().parents[2]
CONFIG_PATH = ROOT / "scripts/v6_2/phase7_v6_2_1_config.yaml"
MODULES = [f"FM{i:02d}" for i in range(1, 9)]


def load_config(path: Path | None = None) -> dict:
    if path is None:
        path = Path(os.environ.get("PHASE7_CONFIG_PATH", CONFIG_PATH))
    cfg = yaml.safe_load(path.read_text())
    cfg["root"] = ROOT
    cfg["out"] = ROOT / cfg["output_dir"]
    cfg["inputs"] = {k: ROOT / v for k, v in cfg["inputs"].items()}
    return cfg


def ensure_dirs(out: Path) -> None:
    for name in ("preflight", "measurement", "annotation", "sentinel", "identifiability", "audit", "handoff"):
        (out / name).mkdir(parents=True, exist_ok=True)


def sha256(path: Path, block: int = 1024 * 1024) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(block):
            h.update(chunk)
    return h.hexdigest()


def write_yaml(obj: dict, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(yaml.safe_dump(obj, sort_keys=False, allow_unicode=True))


def write_json(obj: dict, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(obj, indent=2, ensure_ascii=False) + "\n")


def file_record(path: Path, root: Path = ROOT) -> dict:
    return {
        "path": str(path.relative_to(root) if path.is_relative_to(root) else path),
        "exists": path.exists(),
        "bytes": path.stat().st_size if path.exists() else 0,
        "sha256": sha256(path) if path.exists() and path.is_file() else "",
    }


def normalize_timepoint(value: object, cohort_id: object | None = None) -> str:
    if pd.isna(value):
        return "unknown"
    s = str(value).strip().lower().replace("_", " ").replace("-", " ")
    if not s or s in {"nan", "na", "none", "unknown"}:
        return "unknown"
    if str(cohort_id).strip().upper() == "TASK01":
        task01 = {"t0": "baseline", "t1": "on_treatment", "t2": "post_treatment"}
        if s in task01:
            return task01[s]
    if any(x in s for x in ("baseline", "pre", "before", "week 0", "day 0")):
        return "baseline"
    if any(x in s for x in ("on treatment", "on tx", "during", "interim")):
        return "on_treatment"
    if any(x in s for x in ("post", "after", "resection", "surgery")):
        return "post_treatment"
    m = re.search(r"(?:week|wk|w)\s*(\d+)", s)
    if m:
        return f"week_{int(m.group(1))}"
    return re.sub(r"\s+", "_", s)


def robust_z(values: pd.Series, groups: pd.Series, min_n: int = 10) -> pd.Series:
    x = pd.to_numeric(values, errors="coerce")
    result = pd.Series(np.nan, index=x.index, dtype=float)
    for _, idx in groups.groupby(groups, dropna=False).groups.items():
        v = x.loc[idx]
        finite = v.dropna()
        if len(finite) < min_n:
            continue
        med = finite.median()
        scale = 1.4826 * np.median(np.abs(finite - med))
        if not np.isfinite(scale) or scale < 1e-10:
            scale = finite.std(ddof=0)
        if np.isfinite(scale) and scale >= 1e-10:
            result.loc[idx] = (v - med) / scale
    missing = result.isna() & x.notna()
    finite = x.dropna()
    if missing.any() and len(finite) > 1:
        med = finite.median()
        scale = 1.4826 * np.median(np.abs(finite - med))
        if not np.isfinite(scale) or scale < 1e-10:
            scale = finite.std(ddof=0)
        if np.isfinite(scale) and scale >= 1e-10:
            result.loc[missing] = (x.loc[missing] - med) / scale
    return result


def precision_weight(n_cells: pd.Series, reference: float = 100.0) -> pd.Series:
    n = pd.to_numeric(n_cells, errors="coerce").fillna(0).clip(lower=0)
    return np.minimum(1.0, np.sqrt(n / reference))


def collapse_categories(series: pd.Series, min_count: int = 20) -> pd.Series:
    s = series.fillna("unknown").astype(str)
    keep = s.value_counts()[lambda x: x >= min_count].index
    return s.where(s.isin(keep), "other")


def build_design(frame: pd.DataFrame, continuous: Iterable[str], categorical: Iterable[str]) -> tuple[np.ndarray, list[str]]:
    pieces = [pd.DataFrame({"intercept": np.ones(len(frame))}, index=frame.index)]
    for col in continuous:
        v = pd.to_numeric(frame.get(col), errors="coerce")
        med = v.median()
        v = v.fillna(med if np.isfinite(med) else 0.0)
        sd = v.std(ddof=0)
        pieces.append(pd.DataFrame({col: (v - v.mean()) / sd if sd > 1e-10 else 0.0}, index=frame.index))
    for col in categorical:
        if col not in frame:
            continue
        pieces.append(pd.get_dummies(collapse_categories(frame[col]), prefix=col, drop_first=True, dtype=float))
    design = pd.concat(pieces, axis=1)
    return design.to_numpy(dtype=float), design.columns.tolist()


def weighted_residual(y: np.ndarray, x: np.ndarray, weight: np.ndarray, ridge: float = 1e-6) -> tuple[np.ndarray, float]:
    valid = np.isfinite(y) & np.isfinite(weight) & (weight > 0) & np.isfinite(x).all(axis=1)
    residual = np.full(len(y), np.nan)
    if valid.sum() <= x.shape[1] + 2:
        return residual, np.nan
    xv, yv, wv = x[valid], y[valid], weight[valid]
    root_w = np.sqrt(wv)
    gram = (xv * root_w[:, None]).T @ (xv * root_w[:, None])
    penalty = np.eye(gram.shape[0]) * ridge
    penalty[0, 0] = 0.0
    beta = np.linalg.solve(gram + penalty, (xv * root_w[:, None]).T @ (yv * root_w))
    fitted = xv @ beta
    residual[valid] = yv - fitted
    mean = np.average(yv, weights=wv)
    denom = np.sum(wv * (yv - mean) ** 2)
    r2 = 1.0 - np.sum(wv * (yv - fitted) ** 2) / denom if denom > 0 else np.nan
    return residual, r2


def weighted_group_mean(frame: pd.DataFrame, value: str, weight: str, by: list[str]) -> pd.DataFrame:
    tmp = frame[by + [value, weight]].dropna(subset=[value]).copy()
    tmp["_wx"] = tmp[value] * tmp[weight]
    agg = tmp.groupby(by, dropna=False).agg(_wx=("_wx", "sum"), _w=(weight, "sum"), n_units=(value, "size")).reset_index()
    agg[value] = agg["_wx"] / agg["_w"].replace(0, np.nan)
    return agg.drop(columns="_wx")


def spearman_matrix(frame: pd.DataFrame, columns: list[str]) -> pd.DataFrame:
    return frame[columns].rank().corr(min_periods=20)


def spearman_pair(a: pd.Series, b: pd.Series) -> float:
    joined = pd.concat([pd.to_numeric(a, errors="coerce"), pd.to_numeric(b, errors="coerce")], axis=1).dropna()
    if len(joined) < 3:
        return np.nan
    return float(joined.iloc[:, 0].rank().corr(joined.iloc[:, 1].rank()))


def bh_fdr(pvalues: Iterable[float]) -> np.ndarray:
    p = np.asarray(list(pvalues), dtype=float)
    out = np.full(len(p), np.nan)
    valid = np.isfinite(p)
    if not valid.any():
        return out
    pv = p[valid]
    order = np.argsort(pv)
    ranked = pv[order] * len(pv) / np.arange(1, len(pv) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1].clip(max=1)
    restored = np.empty_like(ranked)
    restored[order] = ranked
    out[valid] = restored
    return out


def hypergeom_sf(k_minus_one: int, population: int, successes: int, draws: int) -> float:
    """P[X > k_minus_one] for a hypergeometric variate, without SciPy."""
    lower = max(k_minus_one + 1, 0, draws - (population - successes))
    upper = min(successes, draws)
    if lower > upper or population <= 0:
        return 0.0
    log_den = math.lgamma(population + 1) - math.lgamma(draws + 1) - math.lgamma(population - draws + 1)
    terms = []
    for i in range(lower, upper + 1):
        log_num = (
            math.lgamma(successes + 1) - math.lgamma(i + 1) - math.lgamma(successes - i + 1)
            + math.lgamma(population - successes + 1)
            - math.lgamma(draws - i + 1)
            - math.lgamma(population - successes - draws + i + 1)
        )
        terms.append(log_num - log_den)
    peak = max(terms)
    return float(min(1.0, math.exp(peak) * sum(math.exp(x - peak) for x in terms)))


def parse_gmt(path: Path) -> dict[str, set[str]]:
    sets: dict[str, set[str]] = {}
    with path.open() as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) >= 3:
                sets[fields[0]] = {g.upper() for g in fields[2:] if g}
    return sets


def forbidden_columns(columns: Iterable[str], tokens: Iterable[str], allow: Iterable[str] = ()) -> list[str]:
    allowed = {x.lower() for x in allow}
    bad = []
    for col in columns:
        low = str(col).lower()
        if low in allowed:
            continue
        if any(re.search(rf"(^|[^a-z0-9]){re.escape(t.lower())}([^a-z0-9]|$)", low) for t in tokens):
            bad.append(str(col))
    return bad


def rank_stability(reference: pd.Series, candidate: pd.Series) -> float:
    joined = pd.concat([reference, candidate], axis=1).dropna()
    if len(joined) < 3:
        return np.nan
    return float(np.corrcoef(joined.iloc[:, 0].rank(), joined.iloc[:, 1].rank())[0, 1])
