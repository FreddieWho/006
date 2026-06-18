"""Convert P2 RAW.tar cohorts to raw h5ad."""
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
    for col in adata.obs.columns:
        if adata.obs[col].dtype == object:
            adata.obs[col] = adata.obs[col].astype(str)
    for col in adata.var.columns:
        if adata.var[col].dtype == object:
            adata.var[col] = adata.var[col].astype(str)
    adata.write_h5ad(path, compression="gzip")
    print(f"  written {path}  shape={adata.shape}")


def read_dense_txt_gz_transpose(path: Path, delimiter: str = "\t") -> anndata.AnnData:
    """Read a genes x cells dense txt.gz and return cells x genes AnnData."""
    print(f"  counting …")
    with gzip.open(path, "rt") as fh:
        header = fh.readline().rstrip("\n").split(delimiter)
        first_data = fh.readline().rstrip("\n").split(delimiter)
        n_rows = 1 + sum(1 for _ in fh)
    if len(header) == len(first_data) - 1:
        # header is pure cell IDs, no gene-name placeholder
        obs_index = header
        n_cols = len(header)
        gene_names = [first_data[0]]
    elif len(header) == len(first_data):
        # header has a gene-name placeholder as first element
        obs_index = header[1:]
        n_cols = len(header) - 1
        gene_names = [header[0]]
    else:
        raise ValueError(
            f"header length ({len(header)}) inconsistent with data row length ({len(first_data)})"
        )
    print(f"  {n_rows} genes x {n_cols} cells")
    X = np.empty((n_rows, n_cols), dtype=np.float32)
    X[0, :] = np.fromstring(" ".join(first_data[1:]), sep=" ", dtype=np.float32)
    with gzip.open(path, "rt") as fh:
        next(fh)
        next(fh)
        for i, line in enumerate(fh, start=1):
            parts = line.rstrip("\n").split(delimiter)
            gene_names.append(parts[0])
            X[i, :] = np.fromstring(" ".join(parts[1:]), sep=" ", dtype=np.float32)
    print(f"  building anndata …")
    adata = anndata.AnnData(
        X=X.T,
        obs=pd.DataFrame(index=obs_index),
        var=pd.DataFrame(index=gene_names),
    )
    return adata


def convert_gse161801():
    cid = "GSE161801"
    src_dir = BASE / "data" / "imm" / cid
    meta_path = src_dir / "GSE161801_K43R_metadata_table.csv.gz"
    print(f"[{cid}] reading metadata …")
    meta = pd.read_csv(meta_path)
    meta = meta.set_index("Cell_barcode")

    print(f"[{cid}] discovering csv.gz count files …")
    csv_paths = sorted(src_dir.glob("GSM*.csv.gz"))
    print(f"  found {len(csv_paths)} files")

    adatas = []
    for p in csv_paths:
        # each csv: gene x cells
        df = pd.read_csv(p, index_col=0)
        ad = anndata.AnnData(
            X=df.values.T.astype(np.float32),
            obs=pd.DataFrame(index=df.columns),
            var=pd.DataFrame(index=df.index),
        )
        ad.obs["csv_file"] = p.name
        adatas.append(ad)

    print(f"[{cid}] concatenating …")
    adata = anndata.concat(adatas, axis=0, join="outer", merge="same")
    print(f"  concat shape: {adata.shape}")

    adata.obs = adata.obs.join(meta, how="left")
    adata.obs["cohort_id"] = cid
    adata.obs["patient_id"] = adata.obs["patient"].astype(str)
    adata.obs["sample_id"] = adata.obs["sample_id"].astype(str)
    adata.obs["timepoint_raw"] = adata.obs["timepoint"].astype(str)
    adata.obs["treatment_regimen_raw"] = "none"
    adata.obs["response_raw"] = "none"
    adata.obs["celltype_raw"] = adata.obs["sorting"].astype(str)
    adata.obs["tissue_source"] = "tumor"
    write_h5ad(adata, cid)


def convert_gse232240():
    cid = "GSE232240"
    src_dir = BASE / "data" / "combo" / cid
    counts_path = src_dir / "GSM7324294_Count_data_IMCISION.txt.gz"
    meta_path = src_dir / "GSM7324295_Meta_data_IMCISION.txt.gz"
    print(f"[{cid}] reading counts …")
    adata = read_dense_txt_gz_transpose(counts_path, delimiter="\t")
    print(f"[{cid}] reading metadata …")
    meta = pd.read_csv(meta_path, sep="\t")
    meta = meta.set_index("cell_id")
    adata.obs = adata.obs.join(meta, how="left")
    adata.obs["cohort_id"] = cid
    adata.obs["patient_id"] = adata.obs["patient"].astype(str)
    adata.obs["timepoint_raw"] = adata.obs["timepoint"].astype(str)
    adata.obs["treatment_regimen_raw"] = "none"
    adata.obs["response_raw"] = "none"
    adata.obs["celltype_raw"] = np.nan
    adata.obs["tissue_source"] = "tumor"
    write_h5ad(adata, cid)


def convert_gse272734():
    cid = "GSE272734"
    src_dir = BASE / "data" / "imm" / cid
    print(f"[{cid}] discovering 10x triples …")
    # group by sample prefix: GSM*_<sample>_barcodes.tsv.gz
    triples = {}
    for p in sorted(src_dir.glob("*_matrix.mtx.gz")):
        stem = p.name.replace("_matrix.mtx.gz", "")
        barcodes = src_dir / f"{stem}_barcodes.tsv.gz"
        features = src_dir / f"{stem}_features.tsv.gz"
        if barcodes.exists() and features.exists():
            triples[stem] = (barcodes, features, p)
    print(f"  found {len(triples)} samples")
    adatas = []
    for stem, (bc_path, feat_path, mtx_path) in triples.items():
        with gzip.open(mtx_path, "rb") as fh:
            X = mmread(fh)
        X = sparse.csr_matrix(X.T)
        barcodes = pd.read_csv(bc_path, header=None, names=["barcode"])
        features = pd.read_csv(feat_path, sep="\t", header=None)
        if features.shape[1] == 3:
            gene_ids = features.iloc[:, 0].astype(str).values
            gene_symbols = features.iloc[:, 1].astype(str).values
        else:
            gene_ids = features.iloc[:, 0].astype(str).values
            gene_symbols = gene_ids
        ad = anndata.AnnData(
            X=X,
            obs=pd.DataFrame(index=barcodes["barcode"]),
            var=pd.DataFrame({"gene_symbol": gene_symbols}, index=gene_ids),
        )
        ad.obs["sample_prefix"] = stem
        ad.var_names_make_unique()
        adatas.append(ad)
    adata = anndata.concat(adatas, axis=0, join="outer", merge="same")
    adata.obs["cohort_id"] = cid
    # sample_prefix format: GSM8409822_13-689-wk0 -> patient=13-689, week=0
    adata.obs["patient_id"] = adata.obs["sample_prefix"].str.extract(r"_([0-9]+-[0-9]+)-wk")[0]
    # infer timepoint from week suffix in sample_prefix
    def _tp_from_prefix(p):
        if "-wk0" in p:
            return "pre"
        if any(f"-wk{w}" in p for w in [3, 6, 9, 12]):
            return "on_treatment"
        return "none"
    adata.obs["timepoint_raw"] = adata.obs["sample_prefix"].apply(_tp_from_prefix)
    adata.obs["treatment_regimen_raw"] = "none"
    adata.obs["response_raw"] = "none"
    adata.obs["celltype_raw"] = np.nan
    adata.obs["tissue_source"] = "tumor"
    write_h5ad(adata, cid)


def convert_gse272735():
    cid = "GSE272735"
    src_dir = BASE / "data" / "imm" / cid
    print(f"[{cid}] discovering 10x triples …")
    triples = {}
    for p in sorted(src_dir.glob("*_matrix.mtx.gz")):
        stem = p.name.replace("_matrix.mtx.gz", "")
        barcodes = src_dir / f"{stem}_barcodes.tsv.gz"
        features = src_dir / f"{stem}_features.tsv.gz"
        if barcodes.exists() and features.exists():
            triples[stem] = (barcodes, features, p)
    print(f"  found {len(triples)} samples")
    adatas = []
    for stem, (bc_path, feat_path, mtx_path) in triples.items():
        with gzip.open(mtx_path, "rb") as fh:
            X = mmread(fh)
        X = sparse.csr_matrix(X.T)
        barcodes = pd.read_csv(bc_path, header=None, names=["barcode"])
        features = pd.read_csv(feat_path, sep="\t", header=None)
        if features.shape[1] == 3:
            gene_ids = features.iloc[:, 0].astype(str).values
            gene_symbols = features.iloc[:, 1].astype(str).values
        else:
            gene_ids = features.iloc[:, 0].astype(str).values
            gene_symbols = gene_ids
        ad = anndata.AnnData(
            X=X,
            obs=pd.DataFrame(index=barcodes["barcode"]),
            var=pd.DataFrame({"gene_symbol": gene_symbols}, index=gene_ids),
        )
        ad.obs["sample_prefix"] = stem
        ad.var_names_make_unique()
        adatas.append(ad)
    adata = anndata.concat(adatas, axis=0, join="outer", merge="same")
    adata.obs["cohort_id"] = cid
    # sample_prefix: GSM8409862_HD-P5-wk0 -> patient=HD-P5, tp from wk suffix
    adata.obs["patient_id"] = adata.obs["sample_prefix"].str.extract(r"GSM[0-9]+_([A-Z0-9\-]+)-wk")[0]

    def _tp_from_prefix(p):
        if "-wk0" in p:
            return "pre"
        if any(f"-wk{w}" in p for w in [3, 6, 9, 12]):
            return "on_treatment"
        return "none"

    adata.obs["timepoint_raw"] = adata.obs["sample_prefix"].apply(_tp_from_prefix)
    adata.obs["treatment_regimen_raw"] = "none"
    adata.obs["response_raw"] = "none"
    adata.obs["celltype_raw"] = np.nan
    adata.obs["tissue_source"] = "tumor"
    write_h5ad(adata, cid)


def convert_gse273718():
    cid = "GSE273718"
    src_dir = BASE / "data" / "imm" / cid
    print(f"[{cid}] discovering 10x triples …")
    triples = {}
    for p in sorted(src_dir.glob("*_matrix.mtx.gz")):
        stem = p.name.replace("_matrix.mtx.gz", "")
        barcodes = src_dir / f"{stem}_barcodes.tsv.gz"
        features = src_dir / f"{stem}_features.tsv.gz"
        if barcodes.exists() and features.exists():
            triples[stem] = (barcodes, features, p)
    print(f"  found {len(triples)} samples")
    adatas = []
    for stem, (bc_path, feat_path, mtx_path) in triples.items():
        with gzip.open(mtx_path, "rb") as fh:
            X = mmread(fh)
        X = sparse.csr_matrix(X.T)
        barcodes = pd.read_csv(bc_path, header=None, names=["barcode"])
        features = pd.read_csv(feat_path, sep="\t", header=None)
        if features.shape[1] == 3:
            gene_ids = features.iloc[:, 0].astype(str).values
            gene_symbols = features.iloc[:, 1].astype(str).values
        else:
            gene_ids = features.iloc[:, 0].astype(str).values
            gene_symbols = gene_ids
        ad = anndata.AnnData(
            X=X,
            obs=pd.DataFrame(index=barcodes["barcode"]),
            var=pd.DataFrame({"gene_symbol": gene_symbols}, index=gene_ids),
        )
        ad.obs["sample_prefix"] = stem
        ad.var_names_make_unique()
        adatas.append(ad)
    adata = anndata.concat(adatas, axis=0, join="outer", merge="same")
    adata.obs["cohort_id"] = cid
    # sample_prefix: GSM8436044_31_pre_PBMC -> patient=31, timepoint=pre, tissue=PBMC
    adata.obs["patient_id"] = adata.obs["sample_prefix"].str.split("_").str[1]
    adata.obs["timepoint_raw"] = adata.obs["sample_prefix"].str.split("_").str[2]
    adata.obs["tissue_source"] = adata.obs["sample_prefix"].str.split("_").str[3]
    adata.obs["treatment_regimen_raw"] = "none"
    adata.obs["response_raw"] = "none"
    adata.obs["celltype_raw"] = np.nan
    write_h5ad(adata, cid)


HANDLERS = {
    "GSE161801": convert_gse161801,
    "GSE232240": convert_gse232240,
    "GSE272734": convert_gse272734,
    "GSE272735": convert_gse272735,
    "GSE273718": convert_gse273718,
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
