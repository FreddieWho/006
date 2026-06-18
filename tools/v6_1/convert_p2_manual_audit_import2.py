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


def assemble_krishna_2021_rcc():
    cid = "krishna_2021_rcc"
    print(f"[{cid}] reading mtx …")
    with open(TMP_OUT / f"{cid}_counts.mtx", "rb") as fh:
        X = mmread(fh)
    X = sparse.csr_matrix(X.T)

    print(f"[{cid}] reading cells …")
    cells = pd.read_csv(TMP_OUT / f"{cid}_cells.txt", header=None, names=["cell"])
    print(f"[{cid}] reading genes …")
    genes = pd.read_csv(TMP_OUT / f"{cid}_genes.txt", header=None, names=["gene_id"])

    print(f"[{cid}] reading annotations …")
    ann_path = BASE / "data" / "combo" / cid / "ccRCC_6pat_cell_annotations.txt"
    ann = pd.read_csv(ann_path, sep="\t")
    ann = ann.set_index("cell")

    obs = pd.DataFrame(index=cells["cell"])
    obs = obs.join(ann, how="left")
    obs["cohort_id"] = cid
    obs["patient_id"] = obs["Sample"].astype(str)
    # infer timepoint from Sample_name
    def _tp(name):
        s = str(name).lower()
        if "untreated" in s:
            return "pre"
        if "exposed" in s:
            return "on_treatment"
        if any(k in s for k in ["resist", "mixed", "cr"]):
            return "post"
        return "none"
    obs["timepoint_raw"] = obs["Sample_name"].apply(_tp)
    obs["treatment_regimen_raw"] = "none"
    obs["response_raw"] = "none"
    obs["celltype_raw"] = obs["cluster_name"].astype(str)
    obs["tissue_source"] = obs["type"].astype(str)

    var = pd.DataFrame(index=genes["gene_id"])
    adata = anndata.AnnData(X=X, obs=obs, var=var)
    write_h5ad(adata, cid)


def assemble_lambrecht_brca():
    src_dir = BASE / "data" / "combo" / "lambrecht_brca"

    subset_specs = [
        ("lambrecht_brca_cohort1_cells",   "1863-counts_cells_cohort1",   "1872-BIOKEY_metaData_cohort1_web.csv"),
        ("lambrecht_brca_cohort2_cells",   "1867-counts_cells_cohort2",   "1871-BIOKEY_metaData_cohort2_web.csv"),
        ("lambrecht_brca_cohort1_tcell",   "1864-counts_tcell_cohort1",   "1870-BIOKEY_metaData_tcells_cohort1_web.csv"),
        ("lambrecht_brca_cohort1_myeloid", "1865-counts_myeloid_cohort1", "1869-BIOKEY_metaData_myeloid_cohort1_web.csv"),
        ("lambrecht_brca_cohort1_dc",      "1866-counts_DC_cohort1",      "1868-BIOKEY_metaData_DC_cohort1_web.csv"),
    ]

    for subset_cid, mat_prefix, meta_file in subset_specs:
        mtx_path = TMP_OUT / f"lambrecht_brca_{mat_prefix}_counts.mtx"
        if not mtx_path.exists():
            print(f"  WARN: {mtx_path} missing, skipping {subset_cid}")
            continue
        print(f"[{subset_cid}] reading mtx …")
        with open(mtx_path, "rb") as fh:
            X = mmread(fh)
        X = sparse.csr_matrix(X.T)

        cells = pd.read_csv(TMP_OUT / f"lambrecht_brca_{mat_prefix}_cells.txt", header=None, names=["cell"])
        genes = pd.read_csv(TMP_OUT / f"lambrecht_brca_{mat_prefix}_genes.txt", header=None, names=["gene_symbol"])

        meta_path = src_dir / meta_file
        if meta_path.exists():
            meta = pd.read_csv(meta_path)
            meta = meta.set_index("Cell")
        else:
            meta = pd.DataFrame(index=cells["cell"])

        obs = pd.DataFrame(index=cells["cell"])
        obs = obs.join(meta, how="left")
        obs["cohort_id"] = subset_cid
        obs["patient_id"] = obs.get("patient_id", "unknown").astype(str)
        obs["timepoint_raw"] = obs.get("timepoint", "none").astype(str)
        obs["treatment_regimen_raw"] = "none"
        obs["response_raw"] = "none"
        obs["celltype_raw"] = obs.get("cellType", obs.get("cellSubType", "unknown")).astype(str)
        obs["tissue_source"] = "tumor"

        var = pd.DataFrame(index=genes["gene_symbol"])
        ad = anndata.AnnData(X=X, obs=obs, var=var)
        print(f"  {subset_cid}: {ad.n_obs} cells")
        write_h5ad(ad, subset_cid)


if __name__ == "__main__":
    targets = sys.argv[1:] if sys.argv[1:] else ["krishna_2021_rcc", "lambrecht_brca"]
    for cid in targets:
        try:
            if cid == "krishna_2021_rcc":
                assemble_krishna_2021_rcc()
            elif cid == "lambrecht_brca":
                assemble_lambrecht_brca()
            else:
                print(f"ERROR: unknown {cid}", file=sys.stderr)
        except Exception as exc:
            print(f"ERROR assembling {cid}: {exc}", file=sys.stderr)
            import traceback
            traceback.print_exc()
    print("Done.")
