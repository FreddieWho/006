"""Convert missing/broken P1 cohorts with fast line-by-line dense reader."""
import gzip
import sys
from pathlib import Path

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
    # h5py can't write mixed-type object columns; force all to string
    for col in adata.obs.columns:
        if adata.obs[col].dtype == object:
            adata.obs[col] = adata.obs[col].astype(str)
    for col in adata.var.columns:
        if adata.var[col].dtype == object:
            adata.var[col] = adata.var[col].astype(str)
    adata.write_h5ad(path, compression="gzip")
    print(f"  written {path}  shape={adata.shape}")


def read_dense_txt_gz_fast(path: Path, delimiter: str = "\t") -> anndata.AnnData:
    print(f"  counting rows …")
    with gzip.open(path, "rt") as fh:
        header = fh.readline().rstrip("\n").split(delimiter)
        first_data = fh.readline().rstrip("\n").split(delimiter)
        n_rows = 1 + sum(1 for _ in fh)
    # header may or may not have a gene-name placeholder column;
    # use the first data row to determine the number of value columns
    n_cols = len(first_data) - 1
    if len(header) == n_cols:
        obs_index = header
    elif len(header) == n_cols + 1:
        obs_index = header[1:]
    else:
        raise ValueError(
            f"header length ({len(header)}) inconsistent with data row length ({len(first_data)})"
        )
    print(f"  {n_rows} rows × {n_cols} cols")
    X = np.empty((n_rows, n_cols), dtype=np.float32)
    gene_names = [first_data[0]]
    X[0, :] = np.fromstring(" ".join(first_data[1:]), sep=" ", dtype=np.float32)
    print(f"  reading …")
    with gzip.open(path, "rt") as fh:
        next(fh)
        next(fh)
        for i, line in enumerate(fh, start=1):
            parts = line.rstrip("\n").split(delimiter)
            gene_names.append(parts[0])
            X[i, :] = np.fromstring(" ".join(parts[1:]), sep=" ", dtype=np.float32)
    print(f"  building anndata …")
    adata = anndata.AnnData(
        X=X.T,  # cells × genes
        obs=pd.DataFrame(index=obs_index),
        var=pd.DataFrame(index=gene_names),
    )
    return adata


def convert_gse123813_bcc():
    cid = "GSE123813_bcc"
    counts_path = BASE / "data" / "imm" / "GSE123813" / "GSE123813_bcc_scRNA_counts.txt.gz"
    meta_path = BASE / "data" / "imm" / "GSE123813" / "GSE123813_bcc_all_metadata.txt.gz"
    print(f"[{cid}] reading counts …")
    adata = read_dense_txt_gz_fast(counts_path, delimiter="\t")
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
    adata = read_dense_txt_gz_fast(counts_path, delimiter="\t")
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


def convert_gse140228_smartseq2():
    cid = "GSE140228_smartseq2"
    counts_path = BASE / "data" / "imm" / "GSE140228" / "GSE140228_read_counts_Smartseq2.csv.gz"
    gene_path = BASE / "data" / "imm" / "GSE140228" / "GSE140228_gene_info_Smartseq2.tsv.gz"

    print(f"[{cid}] reading dense CSV counts …")
    df = pd.read_csv(counts_path, index_col=0, header=0)
    print(f"  df shape: {df.shape}")
    values = df.values.astype(np.float32)
    adata = anndata.AnnData(
        X=values.T,
        obs=pd.DataFrame(index=df.columns),
        var=pd.DataFrame(index=df.index),
    )

    print(f"[{cid}] reading gene info …")
    gene_info = pd.read_csv(gene_path, sep="\t", index_col="SYMBOL")
    # join gene info onto var, keeping all genes from counts
    adata.var = adata.var.join(gene_info, how="left")

    adata.obs["cohort_id"] = cid
    adata.obs["patient_id"] = "unknown"
    adata.obs["treatment_regimen_raw"] = "none"
    adata.obs["response_raw"] = "none"
    adata.obs["timepoint_raw"] = "none"
    adata.obs["celltype_raw"] = np.nan
    adata.obs["tissue_source"] = "unknown"
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
    adata = read_dense_txt_gz_fast(counts_path, delimiter="\t")
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
    X = sparse.csr_matrix(X)  # rows=cells, cols=genes (non-standard orientation in this file)
    print(f"  csr shape: {X.shape}")

    barcodes = pd.read_csv(bc_path, header=0)
    genes = pd.read_csv(gene_path, header=0)

    adata = anndata.AnnData(X=X, obs=pd.DataFrame(index=barcodes.iloc[:, 0]), var=pd.DataFrame(index=genes.iloc[:, 0]))

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


def convert_gse235863_5pt():
    cid = "GSE235863_5pt_cd8"
    src = RAW_OUT / f"{cid.lower()}.h5ad"
    print(f"[{cid}] loading existing h5ad …")
    adata = anndata.read_h5ad(src)
    adata.obs["cohort_id"] = cid
    adata.obs["patient_id"] = adata.obs["patient"].astype(str)
    adata.obs["celltype_raw"] = adata.obs["sub_cluster"]
    adata.obs["tissue_source"] = adata.obs["tissue"]
    # derive timepoint from sample name: e.g. P53-post-P-CD8 → post
    adata.obs["timepoint_raw"] = adata.obs["sample"].str.split("-").str[1].map(
        lambda x: "pre" if str(x).lower() == "pre" else ("post" if str(x).lower() == "post" else str(x).lower())
    )
    adata.obs["treatment_regimen_raw"] = "none"
    adata.obs["response_raw"] = "none"
    write_h5ad(adata, cid)


def convert_gse235863_9pt():
    cid = "GSE235863_9pt_cd45"
    src = RAW_OUT / f"{cid.lower()}.h5ad"
    print(f"[{cid}] loading existing h5ad …")
    adata = anndata.read_h5ad(src)
    adata.obs["cohort_id"] = cid
    adata.obs["patient_id"] = adata.obs["patient"].astype(str)
    adata.obs["celltype_raw"] = adata.obs["major_cluster"]
    adata.obs["tissue_source"] = adata.obs["tissue"]
    adata.obs["timepoint_raw"] = adata.obs["sample"].str.split("-").str[1].map(
        lambda x: "pre" if str(x).lower() == "pre" else ("post" if str(x).lower() == "post" else str(x).lower())
    )
    adata.obs["treatment_regimen_raw"] = "none"
    adata.obs["response_raw"] = "none"
    write_h5ad(adata, cid)


HANDLERS = {
    "GSE123813_bcc": convert_gse123813_bcc,
    "GSE123813_scc": convert_gse123813_scc,
    "GSE140228_droplet": convert_gse140228_droplet,
    "GSE140228_smartseq2": convert_gse140228_smartseq2,
    "GSE207422_sc": convert_gse207422_sc,
    "GSE235863_5pt_cd8": convert_gse235863_5pt,
    "GSE235863_9pt_cd45": convert_gse235863_9pt,
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
