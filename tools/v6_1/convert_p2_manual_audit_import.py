"""Assemble P2 manual-audit exports into h5ad."""
import sys
from pathlib import Path

import pandas as pd
import anndata
from scipy.io import mmread
from scipy import sparse

BASE = Path("/home/huyudi/006")
RAW_OUT = BASE / "data" / "processed" / "srt" / "raw"
RAW_OUT.mkdir(parents=True, exist_ok=True)

TMP_OUT = Path("/tmp/p2_manual_audit_export")


def write_h5ad(adata: anndata.AnnData, cohort_id: str):
    path = RAW_OUT / f"{cohort_id.lower()}.h5ad"
    for col in adata.obs.columns:
        if adata.obs[col].dtype == object:
            adata.obs[col] = adata.obs[col].astype(str)
    for col in adata.var.columns:
        if adata.var[col].dtype == object:
            adata.var[col] = adata.var[col].astype(str)
    adata.write_h5ad(path, compression="gzip")
    print(f"  written {path}  shape={adata.shape}")


def assemble(cid: str, subsets=None):
    if subsets is None:
        subsets = [""]
    adatas = []
    for subset in subsets:
        prefix = f"{cid}_{subset}" if subset else cid
        mtx_path = TMP_OUT / f"{prefix}_counts.mtx"
        obs_path = TMP_OUT / f"{prefix}_obs.csv"
        var_path = TMP_OUT / f"{prefix}_var.csv"
        if not mtx_path.exists():
            print(f"  WARN: {mtx_path} not found, skipping")
            continue
        print(f"[{cid}] reading {prefix} mtx …")
        with open(mtx_path, "rb") as fh:
            X = mmread(fh)
        X = sparse.csr_matrix(X.T)
        print(f"[{cid}] reading {prefix} obs …")
        obs = pd.read_csv(obs_path, index_col=0)
        print(f"[{cid}] reading {prefix} var …")
        var = pd.read_csv(var_path, index_col=0)
        ad = anndata.AnnData(X=X, obs=obs, var=var)
        adatas.append(ad)

    if len(adatas) == 0:
        raise ValueError(f"No adatas assembled for {cid}")
    if len(adatas) == 1:
        adata = adatas[0]
    else:
        print(f"[{cid}] concatenating {len(adatas)} subsets …")
        adata = anndata.concat(adatas, axis=0, join="outer", merge="same")
    write_h5ad(adata, cid)


if __name__ == "__main__":
    targets = sys.argv[1:] if sys.argv[1:] else ["krishna_2021_rcc", "lambrecht_brca"]
    for cid in targets:
        try:
            if cid == "lambrecht_brca":
                assemble(cid, subsets=["cells_cohort1", "cells_cohort2", "tcell_cohort1", "myeloid_cohort1", "DC_cohort1"])
            else:
                assemble(cid)
        except Exception as exc:
            print(f"ERROR assembling {cid}: {exc}", file=sys.stderr)
            import traceback
            traceback.print_exc()
    print("Done.")
