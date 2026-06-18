"""Convert P2 manual-audit cohorts to raw h5ad."""
import sys
from pathlib import Path

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


def _read_sample_timepoint_xlsx(xlsx_path: Path) -> dict:
    """Read Sample -> timepoint mapping from xlsx."""
    import openpyxl
    wb = openpyxl.load_workbook(str(xlsx_path))
    ws = wb.active
    rows = list(ws.iter_rows(values_only=True))
    mapping = {}
    for r in rows[3:]:  # skip header rows
        if r[4] and r[3]:
            sample = str(r[4]).strip()
            tissue = str(r[3]).strip().lower()
            if "pre-treatment" in tissue:
                mapping[sample] = "pre"
            elif "on-treatment" in tissue:
                mapping[sample] = "on_treatment"
            elif "normal" in tissue or "tdln" in tissue:
                mapping[sample] = "pre"
    return mapping


def convert_lambrecht_hnscc():
    cid = "lambrecht_hnscc"
    src_dir = BASE / "data" / "combo" / cid
    adatas = []

    # Load timepoint mapping from xlsx
    xlsx_path = src_dir / "8243-Immunity_metadata.xlsx"
    tp_map = _read_sample_timepoint_xlsx(xlsx_path) if xlsx_path.exists() else {}

    # Tissue subset
    tissue_dir = src_dir / "tissue"
    print(f"[{cid}] reading tissue 10x …")
    ad_tissue = sc.read_10x_mtx(tissue_dir, var_names="gene_symbols", make_unique=True)
    meta_tissue = pd.read_csv(tissue_dir / "metadata.csv")
    meta_tissue = meta_tissue.set_index("Barcode")
    ad_tissue.obs = ad_tissue.obs.join(meta_tissue, how="left")
    ad_tissue.obs["tissue_source"] = "tumor"
    ad_tissue.obs["sample_subset"] = "tissue"
    # extract sample ID from barcode (e.g. 10X_HEN009_AAACCTG...) and map timepoint
    ad_tissue.obs["sample_id_from_barcode"] = ad_tissue.obs.index.str.split("_").str[:2].str.join("_")
    ad_tissue.obs["timepoint_raw"] = ad_tissue.obs["sample_id_from_barcode"].map(tp_map).fillna("none")
    adatas.append(ad_tissue)
    print(f"  tissue: {ad_tissue.n_obs} cells")

    # PBMC subset
    pbmc_dir = src_dir / "HNSCC_PBMC_counts"
    print(f"[{cid}] reading PBMC 10x …")
    ad_pbmc = sc.read_10x_mtx(pbmc_dir, var_names="gene_symbols", make_unique=True)
    meta_pbmc = pd.read_csv(pbmc_dir / "HNSCC_PBMC_metadata.csv")
    meta_pbmc = meta_pbmc.set_index("Barcode")
    ad_pbmc.obs = ad_pbmc.obs.join(meta_pbmc, how="left")
    ad_pbmc.obs["tissue_source"] = "pbmc"
    ad_pbmc.obs["sample_subset"] = "pbmc"
    ad_pbmc.obs["timepoint_raw"] = "none"
    adatas.append(ad_pbmc)
    print(f"  pbmc: {ad_pbmc.n_obs} cells")

    print(f"[{cid}] concatenating …")
    adata = anndata.concat(adatas, axis=0, join="outer", merge="same")
    adata.obs["cohort_id"] = cid
    adata.obs["patient_id"] = adata.obs["Patient"].astype(str)
    adata.obs["treatment_regimen_raw"] = adata.obs.get("Treatment", "none")
    adata.obs["response_raw"] = "none"
    adata.obs["celltype_raw"] = adata.obs["Annotation"]
    write_h5ad(adata, cid)


HANDLERS = {
    "lambrecht_hnscc": convert_lambrecht_hnscc,
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
