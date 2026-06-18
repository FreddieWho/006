"""
Convert P1 cohort local raw counts to canonical h5ad in data/processed/srt/raw/.

Each cohort handler reads local source files, builds anndata.AnnData, and writes
<cohort_id>.h5ad (lowercase with underscores).  Metadata columns are attached to
.obs so that step1_data_governance.py can derive sample_metadata_master without
re-reading the raw matrix.

Usage:
    conda run -n proj006 python tools/v6_1/convert_to_h5ad.py [cohort_id ...]

If no args, converts all P1 cohorts.
"""
import gzip
import sys
import re
import warnings
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd
import scanpy as sc
import anndata
from scipy.io import mmread
from scipy import sparse

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
BASE = Path("/home/huyudi/006")
RAW_OUT = BASE / "data" / "processed" / "srt" / "raw"
RAW_OUT.mkdir(parents=True, exist_ok=True)


# ---------------------------------------------------------------------------
# Utilities
# ---------------------------------------------------------------------------

def write_h5ad(adata: anndata.AnnData, cohort_id: str):
    """Write to canonical path; overwrite allowed."""
    path = RAW_OUT / f"{cohort_id.lower()}.h5ad"
    adata.write_h5ad(path, compression="gzip")
    print(f"  written {path}  shape={adata.shape}")


def read_gzip_tsv(path: Path, **kwargs) -> pd.DataFrame:
    with gzip.open(path, "rt") as fh:
        return pd.read_csv(fh, **kwargs)


def read_gzip_csv(path: Path, **kwargs) -> pd.DataFrame:
    with gzip.open(path, "rt") as fh:
        return pd.read_csv(fh, **kwargs)


# ---------------------------------------------------------------------------
# Handlers
# ---------------------------------------------------------------------------

def convert_gse123813_bcc():
    cid = "GSE123813_bcc"
    counts_path = BASE / "data" / "imm" / "GSE123813" / "GSE123813_bcc_scRNA_counts.txt.gz"
    meta_path = BASE / "data" / "imm" / "GSE123813" / "GSE123813_bcc_all_metadata.txt.gz"

    # scanpy.read_text expects rows=obs, cols=var.
    # The file is gene-major (genes as rows, cells as cols), so first read as-is then transpose.
    print(f"[{cid}] reading dense gene-major counts …")
    # header row contains cell IDs
    adata = sc.read_text(counts_path.as_posix(), delimiter=" ")
    adata = adata.T  # now cells × genes
    adata.var_names.name = "gene"

    print(f"[{cid}] reading metadata …")
    meta = read_gzip_tsv(meta_path, sep="\t")
    meta = meta.set_index("cell.id")
    adata.obs = adata.obs.join(meta)
    adata.obs["cohort_id"] = cid
    adata.obs["patient_id"] = adata.obs["patient"].astype(str)
    adata.obs["treatment_regimen_raw"] = adata.obs["treatment"]
    adata.obs["celltype_raw"] = adata.obs["cluster"]
    adata.obs["response_raw"] = "none"  # no response info
    adata.obs["timepoint_raw"] = adata.obs["treatment"].map(lambda x: "pre" if str(x).lower() == "pre" else "post")
    write_h5ad(adata, cid)


def convert_gse123813_scc():
    cid = "GSE123813_scc"
    counts_path = BASE / "data" / "imm" / "GSE123813" / "GSE123813_scc_scRNA_counts.txt.gz"
    meta_path = BASE / "data" / "imm" / "GSE123813" / "GSE123813_scc_metadata.txt.gz"

    print(f"[{cid}] reading dense gene-major counts …")
    adata = sc.read_text(counts_path.as_posix(), delimiter=" ")
    adata = adata.T
    adata.var_names.name = "gene"

    print(f"[{cid}] reading metadata …")
    meta = read_gzip_tsv(meta_path, sep="\t")
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
    mtx_dir = BASE / "data" / "imm" / "GSE140228"
    meta_path = BASE / "data" / "imm" / "GSE140228" / "GSE140228_UMI_counts_Droplet_cellinfo.tsv.gz"

    print(f"[{cid}] reading 10x mtx …")
    adata = sc.read_10x_mtx(
        mtx_dir.as_posix(),
        prefix="GSE140228_UMI_counts_Droplet_",
        gex_only=False,
        make_unique=True,
    )

    print(f"[{cid}] reading metadata …")
    meta = read_gzip_tsv(meta_path, sep="\t")
    meta = meta.set_index("Barcode")
    adata.obs = adata.obs.join(meta, how="left")
    adata.obs["cohort_id"] = cid
    adata.obs["patient_id"] = adata.obs["Donor"].astype(str)
    adata.obs["tissue_source"] = adata.obs["Tissue"]
    adata.obs["treatment_regimen_raw"] = "none"  # no treatment info
    adata.obs["response_raw"] = "none"
    adata.obs["timepoint_raw"] = "none"
    adata.obs["celltype_raw"] = adata.obs["celltype_global"]
    write_h5ad(adata, cid)


def convert_gse149614():
    cid = "GSE149614"
    counts_path = BASE / "data" / "imm" / "GSE149614" / "GSE149614_HCC.scRNAseq.S71915.count.txt.gz"
    meta_path = BASE / "data" / "imm" / "GSE149614" / "GSE149614_HCC.metadata.updated.txt.gz"

    print(f"[{cid}] reading dense gene-major counts …")
    adata = sc.read_text(counts_path.as_posix(), delimiter="\t")
    adata = adata.T
    adata.var_names.name = "gene"

    print(f"[{cid}] reading metadata …")
    meta = read_gzip_tsv(meta_path, sep="\t")
    meta = meta.set_index("Cell")
    adata.obs = adata.obs.join(meta)
    adata.obs["cohort_id"] = cid
    adata.obs["patient_id"] = adata.obs["patient"].astype(str)
    adata.obs["treatment_regimen_raw"] = "none"
    adata.obs["response_raw"] = "none"
    adata.obs["timepoint_raw"] = "pre"
    adata.obs["tissue_source"] = adata.obs["site"]
    adata.obs["celltype_raw"] = adata.obs["celltype"]
    write_h5ad(adata, cid)


def convert_gse207422_sc():
    cid = "GSE207422_sc"
    counts_path = BASE / "data" / "combo" / "GSE207422" / "GSE207422_NSCLC_scRNAseq_UMI_matrix.txt.gz"
    meta_path = BASE / "data" / "combo" / "GSE207422" / "GSE207422_NSCLC_scRNAseq_metadata.xlsx"

    print(f"[{cid}] reading dense gene-major counts …")
    adata = sc.read_text(counts_path.as_posix(), delimiter="\t")
    adata = adata.T
    adata.var_names.name = "gene"

    print(f"[{cid}] reading metadata …")
    meta = pd.read_excel(meta_path, sheet_name="sheet1")
    meta = meta.set_index("Sample")

    # cell IDs are e.g. BD_immune01_612637 → sample prefix = BD_immune01
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


def convert_gse235863_5pt():
    cid = "GSE235863_5pt_cd8"
    src = BASE / "data" / "combo" / "GSE235863" / "GSE235863_five_patients_scRNAseq_cd8t_raw_counts.h5ad.gz"
    out = RAW_OUT / f"{cid.lower()}.h5ad"
    print(f"[{cid}] decompressing h5ad.gz → {out}")
    import shutil
    with gzip.open(src, "rb") as fi, open(out, "wb") as fo:
        shutil.copyfileobj(fi, fo)
    adata = anndata.read_h5ad(out, backed="r")
    print(f"  verified shape={adata.shape}")


def convert_gse235863_9pt():
    cid = "GSE235863_9pt_cd45"
    src = BASE / "data" / "combo" / "GSE235863" / "GSE235863_nine_patients_scRNAseq_cd45_raw_counts.h5ad.gz"
    out = RAW_OUT / f"{cid.lower()}.h5ad"
    print(f"[{cid}] decompressing h5ad.gz → {out}")
    import shutil
    with gzip.open(src, "rb") as fi, open(out, "wb") as fo:
        shutil.copyfileobj(fi, fo)
    adata = anndata.read_h5ad(out, backed="r")
    print(f"  verified shape={adata.shape}")


def convert_gse243013():
    cid = "GSE243013"
    mtx_path = BASE / "data" / "imm" / "GSE243013" / "GSE243013_NSCLC_immune_scRNA_counts.mtx.gz"
    bc_path = BASE / "data" / "imm" / "GSE243013" / "GSE243013_barcodes.csv.gz"
    gene_path = BASE / "data" / "imm" / "GSE243013" / "GSE243013_genes.csv.gz"
    meta_path = BASE / "data" / "imm" / "GSE243013" / "GSE243013_NSCLC_immune_scRNA_metadata.csv.gz"

    print(f"[{cid}] reading mtx …")
    with gzip.open(mtx_path, "rb") as fh:
        X = mmread(fh)
    X = sparse.csr_matrix(X.T)  # mmread returns genes × cells; we need cells × genes

    print(f"[{cid}] reading barcodes/genes …")
    barcodes = read_gzip_csv(bc_path, header=None, names=["barcode"])
    genes = read_gzip_csv(gene_path, header=None, names=["geneSymbol"])

    adata = anndata.AnnData(X=X, obs=pd.DataFrame(index=barcodes["barcode"]), var=pd.DataFrame(index=genes["geneSymbol"]))

    print(f"[{cid}] reading metadata …")
    meta = read_gzip_csv(meta_path)
    meta = meta.set_index("cellID")
    adata.obs = adata.obs.join(meta, how="left")
    adata.obs["cohort_id"] = cid
    adata.obs["patient_id"] = adata.obs["sampleID"].astype(str)
    adata.obs["treatment_regimen_raw"] = adata.obs["anti-PD1_therapy"]
    adata.obs["response_raw"] = adata.obs["pathological_response"]
    adata.obs["timepoint_raw"] = "pre"  # all pre-treatment per metadata
    adata.obs["celltype_raw"] = adata.obs["major_cell_type"]
    adata.obs["tissue_source"] = "tumor"
    write_h5ad(adata, cid)


def convert_gse288199():
    cid = "GSE288199"
    data_dir = BASE / "data" / "combo" / "GSE288199"

    print(f"[{cid}] discovering sample-timepoint triples …")
    files = sorted(data_dir.glob("GSE288199_HN*_matrix.mtx.gz"))
    adatas = []
    for mtx_f in files:
        m = re.match(r"GSE288199_(HN\d+_(pre|post)_T)_matrix\.mtx\.gz", mtx_f.name)
        if not m:
            continue
        sample_time = m.group(1)
        prefix = f"GSE288199_{sample_time}_"
        print(f"  reading {sample_time} …")
        ad = sc.read_10x_mtx(
            data_dir.as_posix(),
            prefix=prefix,
            gex_only=False,
            make_unique=True,
        )
        ad.obs["sample_time"] = sample_time
        ad.obs["patient_id"] = sample_time.split("_")[0]
        ad.obs["timepoint_raw"] = sample_time.split("_")[1]
        adatas.append(ad)

    if not adatas:
        raise RuntimeError(f"No valid triples found in {data_dir}")

    print(f"[{cid}] concatenating {len(adatas)} samples …")
    adata = anndata.concat(adatas, join="outer", index_unique="-", merge="same")
    adata.obs["cohort_id"] = cid
    adata.obs["treatment_regimen_raw"] = "none"
    adata.obs["response_raw"] = "none"
    adata.obs["celltype_raw"] = np.nan
    adata.obs["tissue_source"] = "tumor"
    write_h5ad(adata, cid)


# ---------------------------------------------------------------------------
# Dispatch table
# ---------------------------------------------------------------------------
HANDLERS = {
    "GSE123813_bcc": convert_gse123813_bcc,
    "GSE123813_scc": convert_gse123813_scc,
    "GSE140228_droplet": convert_gse140228_droplet,
    "GSE149614": convert_gse149614,
    "GSE207422_sc": convert_gse207422_sc,
    "GSE235863_5pt_cd8": convert_gse235863_5pt,
    "GSE235863_9pt_cd45": convert_gse235863_9pt,
    "GSE243013": convert_gse243013,
    "GSE288199": convert_gse288199,
}


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main(argv=None):
    if argv is None:
        argv = sys.argv[1:]

    targets = argv if argv else list(HANDLERS.keys())
    errors = []
    for cid in targets:
        handler = HANDLERS.get(cid)
        if handler is None:
            print(f"ERROR: unknown cohort_id {cid}", file=sys.stderr)
            errors.append(cid)
            continue
        try:
            handler()
        except Exception as exc:
            print(f"ERROR converting {cid}: {exc}", file=sys.stderr)
            errors.append(cid)
    if errors:
        sys.exit(1)
    print("All conversions complete.")


if __name__ == "__main__":
    main()
