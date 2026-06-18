"""Convert GSE176021 per-sample 10x tar.gz + annotations to h5ad."""
import gzip
import re
import sys
import tarfile
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import anndata

BASE = Path("/home/huyudi/006")
RAW_OUT = BASE / "data" / "processed" / "srt" / "raw"
RAW_OUT.mkdir(parents=True, exist_ok=True)

SRC_DIR = BASE / "data" / "imm" / "GSE176021"


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


def tar_to_prefix(tar_name: str) -> str:
    """Map tar inner dir name to annotation barcode prefix.
    MD01-024_tumor_1 -> MD01-024:tumor-1
    """
    parts = tar_name.split("_")
    return parts[0] + ":" + "-".join(parts[1:])


def convert_gse176021():
    cid = "GSE176021"
    print(f"[{cid}] discovering non-VDJ tar.gz files …")
    tar_paths = sorted(
        p for p in SRC_DIR.glob("GSM*.tar.gz")
        if ".vdj." not in p.name and ".vdt." not in p.name and ".tdj." not in p.name
    )
    print(f"  found {len(tar_paths)} samples")

    adatas = []
    for tar_path in tar_paths:
        # e.g. GSM5352886_MD01-024_tumor_1.tar.gz -> inner dir MD01-024_tumor_1
        sample_dir_name = tar_path.name.replace(".tar.gz", "").split("_", 1)[1]
        ann_prefix = tar_to_prefix(sample_dir_name)

        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            with tarfile.open(tar_path, "r:gz") as tf:
                for member in tf.getmembers():
                    # skip macOS resource forks
                    if member.name.startswith("._") or "/._" in member.name:
                        continue
                    tf.extract(member, path=tmpdir)
            mtx_files = list(tmpdir.rglob("matrix.mtx.gz"))
            if not mtx_files:
                print(f"  WARN: no matrix found in {tar_path.name}")
                continue
            mtx_dir = mtx_files[0].parent
            try:
                ad = sc.read_10x_mtx(mtx_dir, var_names="gene_symbols", make_unique=True)
            except Exception as exc:
                print(f"  WARN: read_10x_mtx failed for {tar_path.name}: {exc}")
                continue
            # construct full barcodes: prefix_rawbarcode-1
            ad.obs_names = [f"{ann_prefix}_{b}" for b in ad.obs_names]
            ad.obs["sample_prefix"] = sample_dir_name
            ad.obs["patient_id"] = sample_dir_name.split("_")[0]
            adatas.append(ad)
            print(f"  {sample_dir_name}: {ad.n_obs} cells")

    print(f"[{cid}] concatenating {len(adatas)} samples …")
    adata = anndata.concat(adatas, axis=0, join="outer", merge="same")
    print(f"  concat shape: {adata.shape}")

    # Convert RDS annotations to CSV via R, then load in Python
    print(f"[{cid}] converting RDS annotations to CSV …")
    import subprocess
    r_cmd = """
    library(Seurat)
    cd3 <- readRDS('/tmp/cd3_annotations.rds')
    write.csv(cd3, '/tmp/cd3_annotations.csv', row.names=FALSE)
    cd8 <- readRDS('/tmp/cd8_annotations.rds')
    write.csv(cd8, '/tmp/cd8_annotations.csv', row.names=FALSE)
    """
    subprocess.run(["Rscript", "-e", r_cmd], check=True, capture_output=True)

    print(f"[{cid}] loading CD3 annotations …")
    cd3 = pd.read_csv("/tmp/cd3_annotations.csv")
    cd3 = cd3.set_index("barcode")
    cd3 = cd3.add_prefix("cd3_")

    print(f"[{cid}] loading CD8 annotations …")
    cd8 = pd.read_csv("/tmp/cd8_annotations.csv")
    cd8 = cd8.set_index("barcode")
    cd8 = cd8.add_prefix("cd8_")

    adata.obs = adata.obs.join(cd3, how="left")
    adata.obs = adata.obs.join(cd8, how="left")

    adata.obs["cohort_id"] = cid

    def _tp_from_prefix(p):
        s = str(p)
        if "_W2_" in s or "_W4_" in s:
            return "on_treatment"
        if "_M3_" in s:
            return "post"
        return "post"

    adata.obs["timepoint_raw"] = adata.obs["sample_prefix"].apply(_tp_from_prefix)
    adata.obs["treatment_regimen_raw"] = "none"
    adata.obs["response_raw"] = "none"
    # celltype from CD3 if available, else CD8
    adata.obs["celltype_raw"] = adata.obs["cd3_CellType"].fillna(adata.obs["cd8_CellType"])
    adata.obs["tissue_source"] = adata.obs["sample_prefix"].str.split("_").str[1]
    write_h5ad(adata, cid)


if __name__ == "__main__":
    convert_gse176021()
