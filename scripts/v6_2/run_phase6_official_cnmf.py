#!/usr/bin/env python3
"""Run official cNMF and export consensus modules."""
from __future__ import annotations

import argparse
import importlib.util
import os
import shutil
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd

from phase6_algorithm_common import (
    OUT,
    TOP_N,
    ensure_dirs,
    gene_columns,
    load_primary_tables,
    numeric_matrix,
    score_from_components,
    write_blocked_method,
    write_method_outputs,
)


META_COLS = ["expression_unit_id", "cohort_id", "object_id", "sample_key", "patient_key", "cell_state_level", "cell_state", "layer_family"]


def mode_config(mode: str) -> dict:
    if mode == "smoke":
        return {"max_units": 1200, "ks": ["6", "8"], "n_iter": "5", "numgenes": "500", "consensus_ks": ["6"]}
    if mode == "strengthened":
        return {
            "max_units": 1200,
            "ks": ["6", "8", "10"],
            "n_iter": "8",
            "numgenes": "1220",
            "consensus_ks": ["6", "8", "10"],
        }
    return {"max_units": 1500, "ks": ["6", "8"], "n_iter": "10", "numgenes": "1200", "consensus_ks": ["6"]}


def call(cmd: list[str], run_dir: Path, env: dict, timeout: int = 7200) -> tuple[int, str, str]:
    proc = subprocess.run(cmd, cwd=str(run_dir), env=env, text=True, capture_output=True, timeout=timeout)
    return proc.returncode, proc.stdout, proc.stderr


def build_membership_with_labels(method: str, run_id: str, components: np.ndarray, genes: list[str], labels: list[str]) -> pd.DataFrame:
    rows = []
    for i, weights in enumerate(components):
        module_id = labels[i]
        pos = np.maximum(weights.astype(float), 0.0)
        if pos.max() <= 0:
            pos = np.abs(weights.astype(float))
        max_w = float(pos.max()) or 1.0
        top = set(np.argsort(-pos)[: min(TOP_N, len(genes))].tolist())
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
                    "is_top_gene": "yes" if j in top else "no",
                }
            )
    return pd.DataFrame(rows)


def build_score_df(method: str, run_id: str, meta: pd.DataFrame, scores: np.ndarray, labels: list[str]) -> pd.DataFrame:
    out = meta[META_COLS].copy()
    for i, label in enumerate(labels):
        out[label] = scores[:, i]
    out.insert(0, "run_id", run_id)
    out.insert(0, "method", method)
    return out


def run(mode: str) -> None:
    ensure_dirs()
    method = "official_cnmf"
    cfg = mode_config(mode)
    missing = [pkg for pkg in ["cnmf", "anndata", "scanpy"] if importlib.util.find_spec(pkg) is None]
    cli = shutil.which("cnmf")
    if missing or cli is None:
        write_blocked_method(
            method,
            "official_dependency_missing",
            f"missing_python_packages={','.join(missing)}; cnmf_cli={cli or 'missing'}",
        )
        return

    run_dir = OUT / "methods" / method / "official_cnmf_run"
    if run_dir.exists():
        shutil.rmtree(run_dir)
    run_dir.mkdir(parents=True, exist_ok=True)

    raw = pd.read_parquet(OUT / "inputs/P0_pseudobulk_primary.raw_count.parquet").head(int(cfg["max_units"])).copy()
    genes = gene_columns(raw)
    input_path = run_dir / "P0_pseudobulk_primary.raw_count.cNMF.tsv"
    raw_matrix = raw.set_index("expression_unit_id")[genes].apply(pd.to_numeric, errors="coerce").fillna(0.0)
    raw_matrix = raw_matrix.loc[raw_matrix.sum(axis=1) > 0, raw_matrix.sum(axis=0) > 0]
    raw_matrix.to_csv(input_path, sep="\t")

    env = os.environ.copy()
    env["LD_LIBRARY_PATH"] = f"/opt/anaconda3/lib:{env.get('LD_LIBRARY_PATH', '')}"
    env["NUMBA_CACHE_DIR"] = "/tmp/phase6_numba_cache"
    env["MPLCONFIGDIR"] = "/tmp/phase6_mpl_cache"
    Path(env["NUMBA_CACHE_DIR"]).mkdir(parents=True, exist_ok=True)
    Path(env["MPLCONFIGDIR"]).mkdir(parents=True, exist_ok=True)

    prepare_cmd = [
        cli,
        "prepare",
        "--output-dir",
        str(run_dir),
        "--name",
        "phase6_cnmf",
        "--counts",
        str(input_path),
        "-k",
        *cfg["ks"],
        "--n-iter",
        str(cfg["n_iter"]),
        "--seed",
        "1729",
        "--total-workers",
        "1",
        "--numgenes",
        str(cfg["numgenes"]),
    ]
    code, stdout, stderr = call(prepare_cmd, run_dir, env)
    if code != 0:
        write_blocked_method(method, "official_cnmf_prepare_failed", stderr[-4000:])
        return

    logs = []
    for cmd in [
        [cli, "factorize", "--output-dir", str(run_dir), "--name", "phase6_cnmf", "--worker-index", "0"],
        [cli, "combine", "--output-dir", str(run_dir), "--name", "phase6_cnmf"],
    ]:
        code, so, se = call(cmd, run_dir, env)
        logs.append({"cmd": " ".join(cmd), "returncode": code, "stdout_tail": so[-1000:], "stderr_tail": se[-1000:]})
        if code != 0:
            write_blocked_method(method, "official_cnmf_runtime_failed", se[-4000:])
            return

    component_frames = []
    labels = []
    qc_rows = []
    for k in cfg["consensus_ks"]:
        cmd = [
            cli,
            "consensus",
            "--output-dir",
            str(run_dir),
            "--name",
            "phase6_cnmf",
            "--components",
            str(k),
            "--local-density-threshold",
            "0.01",
        ]
        code, so, se = call(cmd, run_dir, env)
        logs.append({"cmd": " ".join(cmd), "returncode": code, "stdout_tail": so[-1000:], "stderr_tail": se[-1000:]})
        if code != 0:
            qc_rows.append({"method": method, "run_id": f"{method}__k{k}", "run_status": "blocked", "rank": k, "blocked_reason": "consensus_failed"})
            continue
        consensus = sorted(run_dir.glob(f"**/phase6_cnmf.spectra.k_{k}.dt_0_01.consensus.txt"))
        if not consensus:
            qc_rows.append({"method": method, "run_id": f"{method}__k{k}", "run_status": "blocked", "rank": k, "blocked_reason": "no_consensus_spectra"})
            continue
        spectra_df = pd.read_csv(consensus[0], sep="\t", index_col=0)
        component_frames.append((k, spectra_df))
        labels.extend([f"{method}__k{k}_M{i + 1:02d}" for i in range(spectra_df.shape[0])])
        qc_rows.append(
            {
                "method": method,
                "run_id": f"{method}__k{k}",
                "run_status": "complete",
                "rank": k,
                "n_iter": cfg["n_iter"],
                "n_spectra_genes": spectra_df.shape[1],
                "spectra_file": str(consensus[0]),
            }
        )

    if not component_frames:
        write_blocked_method(method, "official_cnmf_no_consensus_spectra", "all requested k consensus runs failed")
        return

    meta, _, log, all_genes = load_primary_tables()
    use_genes = sorted(set.intersection(*(set(df.columns.astype(str)) for _, df in component_frames), set(all_genes)))
    components = np.vstack([df[use_genes].to_numpy(dtype=float) for _, df in component_frames])
    x = numeric_matrix(log, use_genes)
    run_id = f"{method}__strengthened_multi_k" if mode == "strengthened" else f"{method}__k{cfg['consensus_ks'][0]}"
    membership = build_membership_with_labels(method, run_id, components, use_genes, labels)
    score_df = build_score_df(method, run_id, meta, score_from_components(meta, x, components), labels)
    run_qc = pd.DataFrame(qc_rows)
    write_method_outputs(
        method,
        run_id,
        membership,
        score_df,
        run_qc,
        {
            "method": method,
            "run_status": "complete",
            "mode": mode,
            "official_method_required": True,
            "run_scope": "expanded_multi_k" if mode == "strengthened" else mode,
            "prepare_command": prepare_cmd,
            "commands": logs,
            "consensus_ks": cfg["consensus_ks"],
        },
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--mode", choices=["smoke", "full", "strengthened"], default="full")
    run(parser.parse_args().mode)
