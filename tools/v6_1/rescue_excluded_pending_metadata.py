"""Rescue GSE169246, GSE115978, GSE229772 from excluded_pending_metadata.

Patches h5ad obs metadata directly, then step1 freeze scripts must be re-run.
"""
import gzip
from pathlib import Path

import pandas as pd
import anndata

BASE = Path("/home/huyudi/006")
RAW_OUT = BASE / "data" / "processed" / "srt" / "raw"


def write_h5ad(adata: anndata.AnnData, path: Path):
    for col in adata.obs.columns:
        if adata.obs[col].dtype == object:
            adata.obs[col] = adata.obs[col].astype(str)
    for col in adata.var.columns:
        if adata.var[col].dtype == object:
            adata.var[col] = adata.var[col].astype(str)
    adata.write_h5ad(path, compression="gzip")
    print(f"  written {path}  shape={adata.shape}")


# ---------------------------------------------------------------------------
# GSE169246
# ---------------------------------------------------------------------------
def rescue_gse169246():
    cid = "GSE169246"
    path = RAW_OUT / f"{cid.lower()}.h5ad"
    print(f"[{cid}] loading …")
    ad = anndata.read_h5ad(path)

    # Parse timepoint from cell_id: gse169246_AAACCTGAGGTTACCT_Pre_P007_b
    def _parse_tp(cell_id):
        s = str(cell_id)
        if "_Pre_" in s:
            return "pre"
        if "_Post_" in s:
            return "post"
        if "_Prog_" in s:
            return "on_treatment"
        return "none"

    ad.obs["timepoint_raw"] = ad.obs["cell_id"].apply(_parse_tp)

    # Build proper sample_id from donor_id + timepoint
    def _build_sid(row):
        donor = str(row["donor_id"])
        tp = str(row["timepoint_raw"])
        return f"{donor}_{tp}"

    ad.obs["sample_id"] = ad.obs.apply(_build_sid, axis=1)
    ad.obs["patient_id"] = ad.obs["donor_id"].astype(str)

    # Tissue source from organ if available, else tissue
    organ = ad.obs.get("organ", pd.Series("", index=ad.obs.index))
    tissue = ad.obs.get("tissue", pd.Series("", index=ad.obs.index))
    ad.obs["tissue_source"] = organ.where(organ != "", tissue).astype(str)

    print(f"  timepoint distribution:")
    print(ad.obs["timepoint_raw"].value_counts())
    print(f"  unique samples: {ad.obs['sample_id'].nunique()}")
    print(f"  unique patients: {ad.obs['patient_id'].nunique()}")

    write_h5ad(ad, path)


# ---------------------------------------------------------------------------
# GSE115978
# ---------------------------------------------------------------------------
def rescue_gse115978():
    cid = "GSE115978"
    path = RAW_OUT / f"{cid.lower()}.h5ad"
    print(f"[{cid}] loading …")
    ad = anndata.read_h5ad(path)

    # Load metadata CSV
    csv_path = BASE / "data" / "combo" / "GSE115978" / "gse115978.csv"
    meta = pd.read_csv(csv_path)
    meta.columns = [c.lstrip("﻿") for c in meta.columns]

    # Build timepoint mapping from CSV
    tp_map = {"Pre": "pre", "Post": "post"}
    sample_to_tp = {
        row["sample_id"]: tp_map.get(row["timepoint"], "none")
        for _, row in meta.iterrows()
    }

    # Add timepoint_raw to h5ad
    ad.obs["timepoint_raw"] = ad.obs["sample_id"].map(sample_to_tp).fillna("none")

    # Ensure patient_id is set (use donor_id if present, else sample_id)
    if "patient_id" not in ad.obs.columns:
        ad.obs["patient_id"] = ad.obs.get("donor_id", ad.obs["sample_id"]).astype(str)

    # Ensure tissue_source is set (prefer organ, fall back to tissue)
    if "tissue_source" not in ad.obs.columns:
        organ = ad.obs.get("organ", pd.Series("", index=ad.obs.index))
        tissue = ad.obs.get("tissue", pd.Series("", index=ad.obs.index))
        ad.obs["tissue_source"] = organ.where(organ != "", tissue).astype(str)

    print(f"  timepoint distribution:")
    print(ad.obs["timepoint_raw"].value_counts())
    print(f"  unique samples: {ad.obs['sample_id'].nunique()}")
    print(f"  unique patients: {ad.obs['patient_id'].nunique()}")

    write_h5ad(ad, path)


# ---------------------------------------------------------------------------
# GSE229772
# ---------------------------------------------------------------------------
def rescue_gse229772():
    cid = "GSE229772"
    path = RAW_OUT / f"{cid.lower()}.h5ad"
    print(f"[{cid}] loading …")
    ad = anndata.read_h5ad(path)

    # Load id -> sample mapping
    map_path = BASE / "data" / "imm" / "GSE229772" / "GSE229772_scRNA_id_name_match.txt.gz"
    with gzip.open(map_path, "rt") as fh:
        id_sample = pd.read_csv(fh, sep="\t")
    id_sample = id_sample.set_index("id")

    # Extract suffix number from obs_name: gse229772_AAACCTGAGCTCTCGG-1  -> 1
    def _extract_suffix(name):
        parts = str(name).split("-")
        return int(parts[-1]) if parts[-1].isdigit() else -1

    ad.obs["_suffix_id"] = ad.obs_names.map(_extract_suffix)

    # Map suffix_id -> sample name
    ad.obs["sample_id"] = ad.obs["_suffix_id"].map(id_sample["sample"]).astype(str)

    # Parse patient and timepoint from sample name: C26-t0 -> patient=C26, tp=t0
    def _parse_sample(name):
        s = str(name)
        if "-" in s:
            patient, tp = s.split("-", 1)
            return patient, tp
        return s, "none"

    parsed = ad.obs["sample_id"].apply(_parse_sample)
    ad.obs["patient_id"] = parsed.apply(lambda x: x[0])
    tp_raw = parsed.apply(lambda x: x[1])

    def _tp_map(tp):
        if tp == "t0":
            return "pre"
        if tp.startswith("t"):
            return "on_treatment"
        return "none"

    ad.obs["timepoint_raw"] = tp_raw.apply(_tp_map)
    ad.obs["tissue_source"] = "tumor"

    ad.obs.drop(columns=["_suffix_id"], inplace=True, errors="ignore")

    print(f"  timepoint distribution:")
    print(ad.obs["timepoint_raw"].value_counts())
    print(f"  unique samples: {ad.obs['sample_id'].nunique()}")
    print(f"  unique patients: {ad.obs['patient_id'].nunique()}")
    print(f"  sample examples: {sorted(ad.obs['sample_id'].unique())[:10]}")

    write_h5ad(ad, path)


if __name__ == "__main__":
    rescue_gse169246()
    rescue_gse115978()
    rescue_gse229772()
    print("Done.")
