#!/usr/bin/env python3
"""Run higher-power bootstrap for Phase6 finalist methods."""
from __future__ import annotations

import argparse
import ctypes
import os
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

conda_lib = "/opt/anaconda3/lib"
if Path(conda_lib).exists() and conda_lib not in os.environ.get("LD_LIBRARY_PATH", "") and not os.environ.get("PHASE6_BOOT_REEXEC"):
    env = os.environ.copy()
    env["LD_LIBRARY_PATH"] = f"{conda_lib}:{env.get('LD_LIBRARY_PATH', '')}"
    env["PHASE6_BOOT_REEXEC"] = "1"
    os.execvpe(sys.executable, [sys.executable, *sys.argv], env)

conda_libstdcxx = Path("/opt/anaconda3/lib/libstdc++.so.6")
if conda_libstdcxx.exists():
    ctypes.CDLL(str(conda_libstdcxx), mode=ctypes.RTLD_GLOBAL)

from sklearn.decomposition import LatentDirichletAllocation, NMF

from phase6_algorithm_common import OUT, SEED, balanced_indices, load_primary_tables, max_jaccard, numeric_matrix, score_from_components, top_gene_sets, write_csv, write_md


def parse_nmf(run_id: str) -> dict:
    m = re.match(r"tuned_nmf__rank(\d+)__(.+?)__(cd|mu)__(.+)$", run_id)
    if not m:
        raise ValueError(f"cannot parse NMF run_id: {run_id}")
    return {"rank": int(m.group(1)), "init": m.group(2), "solver": m.group(3), "beta_loss": m.group(4)}


def parse_lda(run_id: str) -> dict:
    m = re.match(r"lda_topic__rank(\d+)__batch__(.+?)__tw(.+)$", run_id)
    if not m:
        m = re.match(r"lda_topic__rank(\d+)__batch__(.+)$", run_id)
        if not m:
            raise ValueError(f"cannot parse LDA run_id: {run_id}")
        return {"rank": int(m.group(1)), "prior": m.group(2), "topic_word_prior": 0.01}
    return {"rank": int(m.group(1)), "prior": m.group(2), "topic_word_prior": float(m.group(3))}


def lda_doc_prior(label: str, rank: int) -> float | None:
    if label in {"default", "None"}:
        return None
    if label == "inverse_k":
        return 1.0 / rank
    if label == "fixed_0_05":
        return 0.05
    return None


def fit_nmf(x: np.ndarray, params: dict, seed: int, max_iter: int) -> np.ndarray:
    model = NMF(
        n_components=params["rank"],
        init=params["init"],
        solver=params["solver"],
        beta_loss=params["beta_loss"],
        max_iter=max_iter,
        random_state=seed,
        l1_ratio=0.0,
        alpha_W=0.0,
        alpha_H=0.0,
    )
    model.fit(x)
    return model.components_.astype(float)


def fit_lda(x: np.ndarray, params: dict, seed: int, max_iter: int) -> np.ndarray:
    model = LatentDirichletAllocation(
        n_components=params["rank"],
        learning_method="batch",
        doc_topic_prior=lda_doc_prior(params["prior"], params["rank"]),
        topic_word_prior=float(params["topic_word_prior"]),
        max_iter=max_iter,
        random_state=seed,
        n_jobs=-1,
    )
    model.fit(x)
    return model.components_.astype(float)


def module_score_correlations(x: np.ndarray, base: np.ndarray, boot: np.ndarray) -> list[float]:
    base_scores = score_from_components(pd.DataFrame(index=range(x.shape[0])), x, base)
    boot_scores = score_from_components(pd.DataFrame(index=range(x.shape[0])), x, boot)
    base_df = pd.DataFrame(base_scores).rank()
    boot_df = pd.DataFrame(boot_scores).rank()
    corr = base_df.corrwith(boot_df, axis=0)
    if len(corr) >= min(base.shape[0], boot.shape[0]):
        return [float(x) if pd.notna(x) else 0.0 for x in corr.iloc[: base.shape[0]]]
    matrix = np.corrcoef(base_df.T, boot_df.T)[: base.shape[0], base.shape[0] :]
    return [float(np.nanmax(row)) if np.isfinite(row).any() else 0.0 for row in matrix]


def run(mode: str, bootstrap_n: int) -> None:
    score_dir = OUT / "score"
    summary_path = score_dir / "module_stability_rescore.csv"
    if not summary_path.exists():
        raise FileNotFoundError(f"missing score summary: {summary_path}")
    summary = pd.read_csv(summary_path)
    eligible = summary[summary["method"].isin(["tuned_nmf", "lda_topic"])].copy()
    eligible = eligible.sort_values(["selection_status", "composite_score", "median_bootstrap_jaccard"], ascending=[False, False, False]).head(3)
    meta, raw, log, genes = load_primary_tables()
    fit_rows = balanced_indices(meta, 1200)
    rng = np.random.default_rng(SEED + 9000)
    out_dir = OUT / "bootstrap"
    out_dir.mkdir(parents=True, exist_ok=True)
    summary_rows = []
    module_rows = []

    for _, row in eligible.iterrows():
        method = row["method"]
        run_id = row["selected_run"]
        if method == "tuned_nmf":
            x = numeric_matrix(log, genes)
            x_fit = x[fit_rows]
            params = parse_nmf(run_id)
            base = fit_nmf(x_fit, params, SEED + 1, 220 if mode == "strengthened" else 160)
            fit = lambda xb, seed: fit_nmf(xb, params, seed, 110 if mode == "strengthened" else 90)
        elif method == "lda_topic":
            x = numeric_matrix(raw, genes)
            x_fit = x[fit_rows]
            params = parse_lda(run_id)
            base = fit_lda(x_fit, params, SEED + 2, 40 if mode == "strengthened" else 25)
            fit = lambda xb, seed: fit_lda(xb, params, seed, 18 if mode == "strengthened" else 12)
        else:
            continue

        base40 = top_gene_sets(base, genes, top_n=40)
        base80 = top_gene_sets(base, genes, top_n=80)
        per_module_40 = [[] for _ in base40]
        per_module_80 = [[] for _ in base40]
        per_module_corr = [[] for _ in base40]
        for b in range(bootstrap_n):
            boot_rows = rng.integers(0, x_fit.shape[0], size=x_fit.shape[0])
            boot = fit(x_fit[boot_rows], SEED + b + 100)
            vals40 = max_jaccard(base40, top_gene_sets(boot, genes, top_n=40))
            vals80 = max_jaccard(base80, top_gene_sets(boot, genes, top_n=80))
            cors = module_score_correlations(x_fit, base, boot)
            for i, value in enumerate(vals40):
                per_module_40[i].append(value)
            for i, value in enumerate(vals80):
                per_module_80[i].append(value)
            for i, value in enumerate(cors[: len(per_module_corr)]):
                per_module_corr[i].append(value)

        all40 = [x for vals in per_module_40 for x in vals]
        all80 = [x for vals in per_module_80 for x in vals]
        allcorr = [x for vals in per_module_corr for x in vals]
        summary_rows.append(
            {
                "method": method,
                "run_id": run_id,
                "bootstrap_n": bootstrap_n,
                "median_top40_jaccard": float(np.median(all40)),
                "ci025_top40_jaccard": float(np.quantile(all40, 0.025)),
                "median_top80_jaccard": float(np.median(all80)),
                "median_score_spearman": float(np.nanmedian(allcorr)),
                "fraction_modules_top40_gt_0_40": float(np.mean([np.median(v) > 0.40 for v in per_module_40])),
            }
        )
        for i, vals in enumerate(per_module_40):
            module_rows.append(
                {
                    "method": method,
                    "run_id": run_id,
                    "module_index": i + 1,
                    "median_top40_jaccard": float(np.median(vals)),
                    "ci025_top40_jaccard": float(np.quantile(vals, 0.025)),
                    "median_top80_jaccard": float(np.median(per_module_80[i])),
                    "median_score_spearman": float(np.nanmedian(per_module_corr[i])),
                }
            )

    write_csv(pd.DataFrame(summary_rows), out_dir / "finalist_bootstrap_summary.csv")
    write_csv(pd.DataFrame(module_rows), out_dir / "finalist_bootstrap_module_stability.csv")
    write_md(
        "# Phase6 Finalist Bootstrap Report\n\n"
        f"- Mode: `{mode}`\n"
        f"- Bootstrap replicates per finalist: `{bootstrap_n}`\n"
        f"- Finalists: `{', '.join(eligible['method'].astype(str).tolist())}`\n",
        out_dir / "finalist_bootstrap_report.md",
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--mode", choices=["smoke", "full", "strengthened"], default="strengthened")
    parser.add_argument("--bootstrap-n", type=int, default=30)
    args = parser.parse_args()
    run(args.mode, args.bootstrap_n)
