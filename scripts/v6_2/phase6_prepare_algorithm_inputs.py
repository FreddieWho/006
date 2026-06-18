#!/usr/bin/env python3
"""Prepare fixed inputs for Phase6 module algorithm benchmark."""
from __future__ import annotations

import argparse
import gzip

import pandas as pd

from phase6_algorithm_common import OUT, ensure_dirs, gene_columns, load_primary_tables, sha256_path, write_csv, write_yaml


def compartment(cell_state: str) -> str:
    s = str(cell_state)
    if any(x in s for x in ["CD8", "CD4", "Treg", "NK", "T_NK"]):
        return "lymphoid"
    if any(x in s for x in ["Myeloid", "Macrophage", "DC", "Neutrophil", "APC"]):
        return "myeloid_apc"
    if any(x in s for x in ["Tumor", "CAF", "Stromal", "Endothelial"]):
        return "tumor_stromal"
    if any(x in s for x in ["B_cell", "B_Plasma", "Plasma"]):
        return "b_plasma"
    return "other"


def write_tsv_gz(df: pd.DataFrame, path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(path, "wt") as fh:
        df.to_csv(fh, sep="\t", index=False)


def run(mode: str) -> None:
    ensure_dirs()
    meta, raw, log, genes = load_primary_tables()
    meta = meta.copy()
    meta["compartment"] = meta["cell_state"].map(compartment)
    meta["benchmark_input_scope"] = "P0_pseudobulk_primary"

    input_dir = OUT / "inputs"
    raw_path = input_dir / "P0_pseudobulk_primary.raw_count.parquet"
    log_path = input_dir / "P0_pseudobulk_primary.logcpm.parquet"
    meta_path = input_dir / "P0_pseudobulk_primary.metadata.csv"
    raw_tsv = input_dir / "P0_pseudobulk_primary.raw_count.tsv.gz"
    log_tsv = input_dir / "P0_pseudobulk_primary.logcpm.tsv.gz"
    comp_path = input_dir / "P2_compartment_membership.csv"
    genes_path = input_dir / "gene_universe.txt"

    raw.to_parquet(raw_path, index=False)
    log.to_parquet(log_path, index=False)
    write_csv(meta, meta_path)
    write_csv(meta[["expression_unit_id", "cohort_id", "cell_state", "compartment"]], comp_path)
    genes_path.write_text("\n".join(genes) + "\n", encoding="utf-8")

    # R/offical CLI friendly wide matrices.
    keep = ["expression_unit_id"] + genes
    write_tsv_gz(raw[keep], raw_tsv)
    write_tsv_gz(log[keep], log_tsv)

    summary = pd.DataFrame(
        [
            {"field": "mode", "value": mode},
            {"field": "n_expression_units", "value": len(meta)},
            {"field": "n_genes", "value": len(genes)},
            {"field": "n_cohorts", "value": meta["cohort_id"].nunique()},
            {"field": "n_cell_states", "value": meta["cell_state"].nunique()},
            {"field": "n_compartments", "value": meta["compartment"].nunique()},
        ]
    )
    write_csv(summary, input_dir / "input_summary.csv")
    write_yaml(
        {
            "phase": "phase6_module_algorithm_benchmark_inputs",
            "mode": mode,
            "input_scope": "P0_pseudobulk_primary",
            "p1_cell_level_balanced_status": "not_prepared_in_default_run",
            "p1_reason": "official methods are first benchmarked on frozen Phase6 expression units; cell-level balanced extraction remains optional sensitivity",
            "paths": {
                "raw_count_parquet": str(raw_path),
                "logcpm_parquet": str(log_path),
                "metadata": str(meta_path),
                "raw_count_tsv_gz": str(raw_tsv),
                "logcpm_tsv_gz": str(log_tsv),
                "compartment_membership": str(comp_path),
                "gene_universe": str(genes_path),
            },
            "hashes": {
                "raw_count_parquet_sha256": sha256_path(raw_path),
                "logcpm_parquet_sha256": sha256_path(log_path),
                "metadata_sha256": sha256_path(meta_path),
            },
            "shape": {"n_expression_units": len(meta), "n_genes": len(genes)},
        },
        input_dir / "input_manifest.yaml",
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--mode", choices=["smoke", "full", "strengthened"], default="full")
    run(parser.parse_args().mode)
