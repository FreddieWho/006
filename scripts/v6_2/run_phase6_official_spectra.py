#!/usr/bin/env python3
"""Run official Spectra if available; otherwise record dependency block."""
from __future__ import annotations

import argparse
import importlib.util
import json
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd

from phase6_algorithm_common import (
    OUT,
    ROOT,
    build_membership,
    build_score_matrix,
    ensure_dirs,
    gene_columns,
    load_primary_tables,
    write_blocked_method,
    write_method_outputs,
)


def run(mode: str) -> None:
    ensure_dirs()
    method = "official_spectra"
    conda_lib = "/opt/anaconda3/lib"
    if Path(conda_lib).exists() and conda_lib not in os.environ.get("LD_LIBRARY_PATH", "") and not os.environ.get("PHASE6_SPECTRA_REEXEC"):
        env = os.environ.copy()
        env["LD_LIBRARY_PATH"] = f"{conda_lib}:{env.get('LD_LIBRARY_PATH', '')}"
        env["PHASE6_SPECTRA_REEXEC"] = "1"
        env["NUMBA_CACHE_DIR"] = "/tmp/phase6_numba_cache"
        env["MPLCONFIGDIR"] = "/tmp/phase6_mpl_cache"
        os.execvpe(sys.executable, [sys.executable, *sys.argv], env)
    os.environ.setdefault("NUMBA_CACHE_DIR", "/tmp/phase6_numba_cache")
    os.environ.setdefault("MPLCONFIGDIR", "/tmp/phase6_mpl_cache")
    Path(os.environ["NUMBA_CACHE_DIR"]).mkdir(parents=True, exist_ok=True)
    Path(os.environ["MPLCONFIGDIR"]).mkdir(parents=True, exist_ok=True)
    candidates = ["Spectra", "spectra", "scspectra"]
    found = [pkg for pkg in candidates if importlib.util.find_spec(pkg) is not None]
    if not found:
        write_blocked_method(
            method,
            "official_dependency_missing",
            "missing_python_package=spectra_or_scspectra; official Spectra required, no local replacement allowed",
        )
        return
    try:
        import anndata as ad
        import Spectra

        meta, _, log, all_genes = load_primary_tables()
        if mode == "smoke":
            log = log.head(1200).copy()
            meta = meta.head(1200).copy()
        elif mode == "strengthened":
            log = log.head(5000).copy()
            meta = meta.head(5000).copy()
        else:
            log = log.head(3000).copy()
            meta = meta.head(3000).copy()
        genes = gene_columns(log)
        data = json.loads((ROOT / "mvp/outputs/program_signatures.json").read_text())
        gene_sets = {}
        for key, value in data.items():
            gs = list(value.get("up_genes") or []) + list(value.get("down_genes") or [])
            gs = [g for g in gs if g in genes]
            if len(gs) >= 3:
                gene_sets[key] = gs
        adata = ad.AnnData(
            X=log[genes].apply(pd.to_numeric, errors="coerce").fillna(0.0).to_numpy(float),
            obs=meta.set_index(log["expression_unit_id"].astype(str)),
            var=pd.DataFrame(index=genes),
        )
        Spectra.est_spectra(
            adata,
            gene_sets,
            L=6 if mode == "smoke" else 16 if mode == "strengthened" else 12,
            cell_type_key="cell_state",
            use_cell_types=False,
            use_highly_variable=False,
            num_epochs=50 if mode == "smoke" else 700 if mode == "strengthened" else 500,
            lam=0.01,
            delta=0.001,
            rho=0.001,
            n_top_vals=50,
            overlap_threshold=0.2,
        )
        factors = adata.uns["SPECTRA_factors"]
        if isinstance(factors, pd.DataFrame):
            components = factors.to_numpy(dtype=float)
            factor_genes = [str(c) for c in factors.columns]
        else:
            components = np.asarray(factors, dtype=float)
            factor_genes = genes[: components.shape[1]]
        scores = np.asarray(adata.obsm["SPECTRA_cell_scores"], dtype=float)
        run_id = f"{method}__L{components.shape[0]}"
        membership = build_membership(method, run_id, components, factor_genes)
        score_df = build_score_matrix(method, run_id, meta.reset_index(drop=True), scores)
        run_qc = pd.DataFrame(
            [
                {
                    "method": method,
                    "run_id": run_id,
                    "run_status": "complete",
                    "rank": components.shape[0],
                    "n_gene_sets": len(gene_sets),
                    "num_epochs": 50 if mode == "smoke" else 700 if mode == "strengthened" else 500,
                }
            ]
        )
        write_method_outputs(
            method,
            run_id,
            membership,
            score_df,
            run_qc,
            {
                "run_status": "complete",
                "mode": mode,
                "official_method_required": True,
                "official_package": "Spectra",
                "prior_guided_module_candidate": True,
            },
        )
    except Exception as exc:  # noqa: BLE001
        write_blocked_method(
            method,
            "official_runtime_error",
            repr(exc)[:4000],
        )


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--mode", choices=["smoke", "full", "strengthened"], default="full")
    run(parser.parse_args().mode)
