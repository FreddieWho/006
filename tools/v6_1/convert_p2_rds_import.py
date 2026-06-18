"""Assemble P2 RDS exports into h5ad."""
import gzip
from pathlib import Path

import pandas as pd
import scanpy as sc
import anndata
from scipy.io import mmread
from scipy import sparse

BASE = Path("/home/huyudi/006")
RAW_OUT = BASE / "data" / "processed" / "srt" / "raw"
RAW_OUT.mkdir(parents=True, exist_ok=True)

TMP_OUT = Path("/tmp/p2_rds_export")


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


def assemble(cid: str):
    print(f"[{cid}] reading mtx …")
    with open(TMP_OUT / f"{cid}_counts.mtx", "rb") as fh:
        X = mmread(fh)
    X = sparse.csr_matrix(X.T)  # cells x genes

    print(f"[{cid}] reading obs …")
    obs = pd.read_csv(TMP_OUT / f"{cid}_obs.csv", index_col=0)

    print(f"[{cid}] reading var …")
    var = pd.read_csv(TMP_OUT / f"{cid}_var.csv", index_col=0)

    print(f"[{cid}] building anndata …")
    adata = anndata.AnnData(X=X, obs=obs, var=var)
    write_h5ad(adata, cid)


if __name__ == "__main__":
    import sys
    targets = sys.argv[1:] if sys.argv[1:] else ["GSE225063", "GSE272993"]
    for cid in targets:
        try:
            assemble(cid)
        except Exception as exc:
            print(f"ERROR assembling {cid}: {exc}", file=sys.stderr)
            import traceback
            traceback.print_exc()
    print("Done.")
