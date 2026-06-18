#!/usr/bin/env python3
"""Run tuned NMF and LDA baselines for Phase6 module benchmark."""
from __future__ import annotations

import argparse
import ctypes
import itertools
import os
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

conda_lib = "/opt/anaconda3/lib"
if Path(conda_lib).exists() and conda_lib not in os.environ.get("LD_LIBRARY_PATH", "") and not os.environ.get("PHASE6_LIBSTDCPP_REEXEC"):
    env = os.environ.copy()
    env["LD_LIBRARY_PATH"] = f"{conda_lib}:{env.get('LD_LIBRARY_PATH', '')}"
    env["PHASE6_LIBSTDCPP_REEXEC"] = "1"
    os.execvpe(sys.executable, [sys.executable, *sys.argv], env)

conda_libstdcxx = Path("/opt/anaconda3/lib/libstdc++.so.6")
if conda_libstdcxx.exists():
    ctypes.CDLL(str(conda_libstdcxx), mode=ctypes.RTLD_GLOBAL)

from sklearn.decomposition import LatentDirichletAllocation, NMF

from phase6_algorithm_common import (
    SEED,
    balanced_indices,
    build_membership,
    build_score_matrix,
    ensure_dirs,
    load_primary_tables,
    max_jaccard,
    numeric_matrix,
    score_from_components,
    top_gene_sets,
    write_method_outputs,
)


def mode_config(mode: str, method: str) -> dict:
    if mode == "smoke":
        return {
            "fit_rows": 1200,
            "boot_n": 4,
            "nmf_ranks": [6, 8],
            "nmf_inits": ["nndsvda"],
            "nmf_specs": [("mu", "frobenius")],
            "nmf_iter": 300,
            "nmf_boot_iter": 120,
            "lda_ranks": [8, 12],
            "lda_priors": ["default"],
            "lda_topic_word_priors": [0.01],
            "lda_iter": 20,
            "lda_boot_iter": 10,
        }
    if mode == "strengthened":
        return {
            "fit_rows": 1200,
            "boot_n": 1,
            "nmf_ranks": [6, 8, 12, 16],
            "nmf_inits": ["nndsvda"],
            "nmf_specs": [("mu", "frobenius")],
            "nmf_iter": 100,
            "nmf_boot_iter": 40,
            "lda_ranks": [6, 8, 10, 12, 16, 20],
            "lda_priors": ["default", "inverse_k"],
            "lda_topic_word_priors": [0.01],
            "lda_iter": 25,
            "lda_boot_iter": 10,
        }
    return {
        "fit_rows": 1200,
        "boot_n": 2,
        "nmf_ranks": [6, 8, 10, 12],
        "nmf_inits": ["nndsvda"],
        "nmf_specs": [("mu", "frobenius")],
        "nmf_iter": 250,
        "nmf_boot_iter": 120,
        "lda_ranks": [8, 12, 16],
        "lda_priors": ["default", "inverse_k"],
        "lda_topic_word_priors": [0.01],
        "lda_iter": 50,
        "lda_boot_iter": 25,
    }


def bootstrap_jaccard(
    x_fit: np.ndarray,
    genes: list[str],
    base_components: np.ndarray,
    fit_func,
    n_boot: int,
    seed: int,
) -> tuple[float, float, float, int]:
    rng = np.random.default_rng(seed)
    base_sets_40 = top_gene_sets(base_components, genes, top_n=40)
    base_sets_80 = top_gene_sets(base_components, genes, top_n=80)
    vals_40: list[float] = []
    vals_80: list[float] = []
    for b in range(n_boot):
        rows = rng.integers(0, x_fit.shape[0], size=x_fit.shape[0])
        components = fit_func(x_fit[rows], seed + b + 1)
        vals_40.extend(max_jaccard(base_sets_40, top_gene_sets(components, genes, top_n=40)))
        vals_80.extend(max_jaccard(base_sets_80, top_gene_sets(components, genes, top_n=80)))
    if not vals_40:
        return 0.0, 0.0, 0.0, 0
    return (
        float(np.median(vals_40)),
        float(np.quantile(vals_40, 0.025)),
        float(np.median(vals_80)),
        len(vals_40),
    )


def lda_doc_prior(label: str, rank: int) -> float | None:
    if label == "default":
        return None
    if label == "inverse_k":
        return 1.0 / rank
    if label == "fixed_0_05":
        return 0.05
    raise ValueError(f"unknown LDA prior: {label}")


def run_nmf(mode: str) -> None:
    ensure_dirs()
    cfg = mode_config(mode, "tuned_nmf")
    meta, _, log, genes = load_primary_tables()
    x = numeric_matrix(log, genes)
    fit_rows = balanced_indices(meta, int(cfg["fit_rows"]))
    x_fit = x[fit_rows]
    rows = []
    best = None
    best_components = None

    for rank, init, (solver, beta_loss) in itertools.product(cfg["nmf_ranks"], cfg["nmf_inits"], cfg["nmf_specs"]):
        start = time.time()
        model = NMF(
            n_components=rank,
            init=init,
            solver=solver,
            beta_loss=beta_loss,
            max_iter=int(cfg["nmf_iter"]),
            random_state=SEED + rank + len(rows),
            l1_ratio=0.0,
            alpha_W=0.0,
            alpha_H=0.0,
        )
        model.fit(x_fit)
        components = model.components_.astype(float)

        def fit_func(xb: np.ndarray, seed: int) -> np.ndarray:
            m = NMF(
                n_components=rank,
                init=init,
                solver=solver,
                beta_loss=beta_loss,
                max_iter=int(cfg["nmf_boot_iter"]),
                random_state=seed,
                l1_ratio=0.0,
                alpha_W=0.0,
                alpha_H=0.0,
            )
            m.fit(xb)
            return m.components_.astype(float)

        median_j, ci_low, median_j80, n_vals = bootstrap_jaccard(
            x_fit, genes, components, fit_func, int(cfg["boot_n"]), SEED + rank * 1000 + len(rows)
        )
        rec = float(model.reconstruction_err_)
        row = {
            "method": "tuned_nmf",
            "run_id": f"tuned_nmf__rank{rank}__{init}__{solver}__{beta_loss}",
            "run_status": "complete",
            "rank": rank,
            "init": init,
            "solver": solver,
            "beta_loss": beta_loss,
            "reconstruction_error": rec,
            "median_bootstrap_jaccard": median_j,
            "bootstrap_jaccard_ci_low": ci_low,
            "median_bootstrap_jaccard_top80": median_j80,
            "n_bootstrap_values": n_vals,
            "fit_rows": len(fit_rows),
            "runtime_seconds": round(time.time() - start, 3),
        }
        rows.append(row)
        score_key = (median_j, median_j80, -rec, -rank)
        if best is None or score_key > best[0]:
            best = (score_key, row)
            best_components = components

    qc = pd.DataFrame(rows)
    best_row = best[1]
    run_id = best_row["run_id"]
    scores = score_from_components(meta, x, best_components)
    membership = build_membership("tuned_nmf", run_id, best_components, genes)
    score_df = build_score_matrix("tuned_nmf", run_id, meta, scores)
    write_method_outputs(
        "tuned_nmf",
        run_id,
        membership,
        score_df,
        qc,
        {"run_status": "complete", "selected_run": run_id, "mode": mode, "official_method_required": False},
    )


def run_lda(mode: str) -> None:
    ensure_dirs()
    cfg = mode_config(mode, "lda_topic")
    meta, raw, _, genes = load_primary_tables()
    x = numeric_matrix(raw, genes)
    fit_rows = balanced_indices(meta, int(cfg["fit_rows"]))
    x_fit = x[fit_rows]
    rows = []
    best = None
    best_components = None

    for rank, prior_label, topic_word_prior in itertools.product(
        cfg["lda_ranks"], cfg["lda_priors"], cfg["lda_topic_word_priors"]
    ):
        doc_topic_prior = lda_doc_prior(str(prior_label), int(rank))
        start = time.time()
        model = LatentDirichletAllocation(
            n_components=rank,
            learning_method="batch",
            doc_topic_prior=doc_topic_prior,
            topic_word_prior=float(topic_word_prior),
            max_iter=int(cfg["lda_iter"]),
            random_state=SEED + rank + len(rows),
            n_jobs=-1,
        )
        model.fit(x_fit)
        components = model.components_.astype(float)

        def fit_func(xb: np.ndarray, seed: int) -> np.ndarray:
            m = LatentDirichletAllocation(
                n_components=rank,
                learning_method="batch",
                doc_topic_prior=doc_topic_prior,
                topic_word_prior=float(topic_word_prior),
                max_iter=int(cfg["lda_boot_iter"]),
                random_state=seed,
                n_jobs=-1,
            )
            m.fit(xb)
            return m.components_.astype(float)

        median_j, ci_low, median_j80, n_vals = bootstrap_jaccard(
            x_fit, genes, components, fit_func, int(cfg["boot_n"]), SEED + rank * 2000 + len(rows)
        )
        perp = float(model.perplexity(x_fit))
        row = {
            "method": "lda_topic",
            "run_id": f"lda_topic__rank{rank}__batch__{prior_label}__tw{topic_word_prior}",
            "run_status": "complete",
            "rank": rank,
            "learning_method": "batch",
            "doc_topic_prior": prior_label,
            "topic_word_prior": topic_word_prior,
            "perplexity": perp,
            "median_bootstrap_jaccard": median_j,
            "bootstrap_jaccard_ci_low": ci_low,
            "median_bootstrap_jaccard_top80": median_j80,
            "n_bootstrap_values": n_vals,
            "fit_rows": len(fit_rows),
            "runtime_seconds": round(time.time() - start, 3),
        }
        rows.append(row)
        score_key = (median_j, median_j80, -perp, -rank)
        if best is None or score_key > best[0]:
            best = (score_key, row)
            best_components = components

    qc = pd.DataFrame(rows)
    best_row = best[1]
    run_id = best_row["run_id"]
    scores = score_from_components(meta, x, best_components)
    membership = build_membership("lda_topic", run_id, best_components, genes)
    score_df = build_score_matrix("lda_topic", run_id, meta, scores)
    write_method_outputs(
        "lda_topic",
        run_id,
        membership,
        score_df,
        qc,
        {"run_status": "complete", "selected_run": run_id, "mode": mode, "official_method_required": False},
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--mode", choices=["smoke", "full", "strengthened"], default="full")
    parser.add_argument("--method", choices=["all", "tuned_nmf", "lda_topic"], default="all")
    args = parser.parse_args()
    if args.method in {"all", "tuned_nmf"}:
        run_nmf(args.mode)
    if args.method in {"all", "lda_topic"}:
        run_lda(args.mode)


if __name__ == "__main__":
    main()
