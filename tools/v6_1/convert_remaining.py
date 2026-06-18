"""Convert missing/broken P1 cohorts with robust dense reader."""
import gzip
import sys
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd
import scanpy as sc
import anndata
from scipy.io import mmread
from scipy import sparse

BASE = Path("/home/huyudi/006")
RAW_OUT = BASE / "data" / "processed" / "srt" / "raw"
RAW_OUT.mkdir(parents=True, exist_ok=True)


def write_h5ad(adata: anndata.AnnData, cohort_id: str):
    path = RAW_OUT / f"{cohort_id.lower()}.h5ad"
    adata.write_h5ad(path, compression="gzip")
    print(f"  written {path}  shape={adata.shape}")


def read_dense_txt_gz(path: Path, delimiter: str = "\t") -> anndata.AnnData:
    print(f"  reading dense matrix with pandas …")
    df = pd.read_csv(path, sep=delimiter, index_col=0, header=0, engine="c")
    print(f"  df shape: {df.shape}")
    values = df.values.astype(np.float32)
    adata = anndata.AnnData(
        X=values.T,
        obs=pd.DataFrame(index=df.columns),
        var=pd.DataFrame(index=df.index),
    )
    return adata


def convert_gse123813_bcc():
    cid = "GSE123813_bcc"
    counts_path = BASE / "data" / "imm" / "GSE123813" / "GSE123813_bcc_scRNA_counts.txt.gz"
    meta_path = BASE / "data" / "imm" / "GSE123813" / "GSE123813_bcc_all_metadata.txt.gz"
    print(f"[{cid}] reading counts …")
    adata = read_dense_txt_gz(counts_path, delimiter=" ")
    print(f"[{cid}] reading metadata …")
    meta = pd.read_csv(meta_path, sep="\t")
    meta = meta.set_index("cell.id")
    adata.obs = adata.obs.join(meta)
    adata.obs["cohort_id"] = cid
    adata.obs["patient_id"] = adata.obs["patient"].astype(str)
    adata.obs["treatment_regimen_raw"] = adata.obs["treatment"]
    adata.obs["celltype_raw"] = adata.obs["cluster"]
    adata.obs["response_raw"] = "none"
    adata.obs["timepoint_raw"] = adata.obs["treatment"].map(lambda x: "pre" if str(x).lower() == "pre" else "post")
    write_h5ad(adata, cid)


def convert_gse123813_scc():
    cid = "GSE123813_scc"
    counts_path = BASE / "data" / "imm" / "GSE123813" / "GSE123813_scc_scRNA_counts.txt.gz"
    meta_path = BASE / "data" / "imm" / "GSE123813" / "GSE123813_scc_metadata.txt.gz"
    print(f"[{cid}] reading counts …")
    adata = read_dense_txt_gz(counts_path, delimiter=" ")
    print(f"[{cid}] reading metadata …")
    meta = pd.read_csv(meta_path, sep="\t")
    meta = meta.set_index("cell.id")
    adata.obs = adata.obs.join(meta)
    adata.obs["cohort_id"] = cid
    adata.obs["patient_id"] = adata.obs["patient"].astype(str)
    adata.obs["treatment_regimen_raw"] = adata.obs["treatment"]
    adata.obs["celltype_raw"] = adata.obs["cluster"]
    adata.obs["response_raw"] = "none"
    adata.obs["timepoint_raw"] = adata.obs["treatment"].map(lambda x: "pre" if str(x).lower() == "pre" else "post")
    write_h5ad(adata, cid)


def convert_gse140228_droplet():
    cid = "GSE140228_droplet"
    mtx_path = BASE / "data" / "imm" / "GSE140228" / "GSE140228_UMI_counts_Droplet.mtx.gz"
    bc_path = BASE / "data" / "imm" / "GSE140228" / "GSE140228_UMI_counts_Droplet_barcodes.tsv.gz"
    gene_path = BASE / "data" / "imm" / "GSE140228" / "GSE140228_UMI_counts_Droplet_genes.tsv.gz"
    meta_path = BASE / "data" / "imm" / "GSE140228" / "GSE140228_UMI_counts_Droplet_cellinfo.tsv.gz"
    print(f"[{cid}] reading mtx …")
    with gzip.open(mtx_path, "rb") as fh:
        X = mmread(fh)
    X = sparse.csr_matrix(X.T)  # cells × genes
    barcodes = pd.read_csv(bc_path, header=None, names=["barcode"])
    genes = pd.read_csv(gene_path, sep="\t", header=0)
    gene_names = genes["SYMBOL"].astype(str).values
    adata = anndata.AnnData(
        X=X,
        obs=pd.DataFrame(index=barcodes["barcode"]),
        var=pd.DataFrame(index=gene_names),
    )
    # deduplicate var names
    adata.var_names_make_unique()
    print(f"[{cid}] reading metadata …")
    meta = pd.read_csv(meta_path, sep="\t")
    meta = meta.set_index("Barcode")
    adata.obs = adata.obs.join(meta, how="left")
    adata.obs["cohort_id"] = cid
    adata.obs["patient_id"] = adata.obs["Donor"].astype(str)
    adata.obs["tissue_source"] = adata.obs["Tissue"]
    adata.obs["treatment_regimen_raw"] = "none"
    adata.obs["response_raw"] = "none"
    adata.obs["timepoint_raw"] = "none"
    adata.obs["celltype_raw"] = adata.obs["celltype_global"]
    write_h5ad(adata, cid)


def convert_gse207422_sc():
    cid = "GSE207422_sc"
    counts_path = BASE / "data" / "combo" / "GSE207422" / "GSE207422_NSCLC_scRNAseq_UMI_matrix.txt.gz"
    meta_path = BASE / "data" / "combo" / "GSE207422" / "GSE207422_NSCLC_scRNAseq_metadata.xlsx"
    print(f"[{cid}] reading counts …")
    adata = read_dense_txt_gz(counts_path, delimiter="\t")
    print(f"[{cid}] reading metadata …")
    meta = pd.read_excel(meta_path, sheet_name="sheet1")
    meta = meta.set_index("Sample")
    adata.obs["sample_prefix"] = adata.obs.index.to_series().str.rsplit("_", n=1).str[0]
    adata.obs = adata.obs.join(meta, on="sample_prefix", how="left")
    adata.obs["cohort_id"] = cid
    adata.obs["patient_id"] = adata.obs["Patient"].astype(str)
    adata.obs["treatment_regimen_raw"] = adata.obs["PD1 Antibody"]
    adata.obs["response_raw"] = adata.obs["RECIST"]
    adata.obs["timepoint_raw"] = adata.obs["Resource"].map(
        lambda x: "pre" if "pre" in str(x).lower() else "post"
    )
    adata.obs["celltype_raw"] = np.nan
    adata.obs["tissue_source"] = "tumor"
    write_h5ad(adata, cid)


def convert_gse243013():
    cid = "GSE243013"
    mtx_path = BASE / "data" / "imm" / "GSE243013" / "GSE243013_NSCLC_immune_scRNA_counts.mtx.gz"
    bc_path = BASE / "data" / "imm" / "GSE243013" / "GSE243013_barcodes.csv.gz"
    gene_path = BASE / "data" / "imm" / "GSE243013" / "GSE243013_genes.csv.gz"
    meta_path = BASE / "data" / "imm" / "GSE243013" / "GSE243013_NSCLC_immune_scRNA_metadata.csv.gz"

    print(f"[{cid}] reading mtx …")
    with gzip.open(mtx_path, "rb") as fh:
        X = mmread(fh)
    print(f"  mmread shape: {X.shape}")
    X = sparse.csr_matrix(X)  # rows=cells, cols=genes (non-standard 10x orientation)
    print(f"  csr shape: {X.shape}")

    barcodes = pd.read_csv(bc_path, header=None, names=["barcode"])
    genes = pd.read_csv(gene_path, header=None, names=["geneSymbol"])

    adata = anndata.AnnData(X=X, obs=pd.DataFrame(index=barcodes["barcode"]), var=pd.DataFrame(index=genes["geneSymbol"]))

    print(f"[{cid}] reading metadata …")
    meta = pd.read_csv(meta_path)
    meta = meta.set_index("cellID")
    adata.obs = adata.obs.join(meta, how="left")
    adata.obs["cohort_id"] = cid
    adata.obs["patient_id"] = adata.obs["sampleID"].astype(str)
    adata.obs["treatment_regimen_raw"] = adata.obs["anti-PD1_therapy"]
    adata.obs["response_raw"] = adata.obs["pathological_response"]
    adata.obs["timepoint_raw"] = "pre"
    adata.obs["celltype_raw"] = adata.obs["major_cell_type"]
    adata.obs["tissue_source"] = "tumor"
    write_h5ad(adata, cid)


HANDLERS = {
    "GSE123813_bcc": convert_gse123813_bcc,
    "GSE123813_scc": convert_gse123813_scc,
    "GSE140228_droplet": convert_gse140228_droplet,
    "GSE207422_sc": convert_gse207422_sc,
    "GSE243013": convert_gse243013,
}

if __name__ == "__main__":
    targets = sys.argv[1:] if sys.argv[1:] else list(HANDLERS.keys())
    for cid in targets:
        handler = HANDLERS.get(cid)
        if handler is None:
            print(f"ERROR: unknown {cid}", file=sys.stderr)
            continue
        try:
            handler()
        except Exception as exc:
            print(f"ERROR converting {cid}: {exc}", file=sys.stderr)
            import traceback
            traceback.print_exc()
    print("Done.")
