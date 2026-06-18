#!/usr/bin/env python3
"""Build a cell-level h5ad for GSE301741 from RAW 10x matrices.

The deposited Seurat RDS is 13 GB and too large to load comfortably, so we
rebuild from the RAW tar `filtered_feature_bc_matrix.h5` files and merge them
with the frozen sample/patient metadata.
"""

import re
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc

warnings.filterwarnings("ignore", category=FutureWarning)

ROOT = Path(__file__).resolve().parents[1]
EXTRACT_DIR = ROOT / "data/external_anchor_candidates/GSE301741_hnscc_pembrolizumab_scRNA/_extracted"
SAMPLE_MANIFEST = ROOT / "data/external_anchor_candidates/GSE301741_hnscc_pembrolizumab_scRNA/GSE301741_sample_manifest.tsv"
SAMPLE_META = ROOT / "results/v6_2/sample_metadata_master.frozen_v0.csv"
PATIENT_META = ROOT / "results/v6_2/patient_metadata_master.frozen_v0.csv"
OUT_H5AD = ROOT / "data/processed/srt/raw/gse301741.h5ad"


def main():
    sample_manifest = pd.read_csv(SAMPLE_MANIFEST, sep="\t")
    sample_meta = pd.read_csv(SAMPLE_META)
    patient_meta = pd.read_csv(PATIENT_META)

    # Build GSM -> (patient_id, timepoint, library_type, cell_type)
    gsm_info = {}
    for _, r in sample_manifest.iterrows():
        gsm = r["geo_accession"]
        if not isinstance(gsm, str):
            continue
        gsm_info[gsm] = {
            "patient_id": str(r.get("patient_id", "")),
            "timepoint": str(r.get("timepoint", "")),
            "library_type": str(r.get("library_type", "")),
            "cell_type": str(r.get("cell_type", "")),
        }

    # Frozen sample metadata for GSE301741
    frozen_samples = sample_meta[sample_meta["cohort_id"] == "GSE301741"].copy()
    frozen_samples["patient_timepoint"] = (
        frozen_samples["patient_id"] + "_" + frozen_samples["timepoint"]
    )

    adatas = []
    skipped = []
    h5_files = sorted(EXTRACT_DIR.glob("*_filtered_feature_bc_matrix.h5"))
    if not h5_files:
        print("No filtered_feature_bc_matrix.h5 files found", file=sys.stderr)
        sys.exit(1)

    for h5_path in h5_files:
        m = re.match(r"(GSM\d+)_.*", h5_path.name)
        if not m:
            skipped.append((h5_path.name, "no_GSM"))
            continue
        gsm = m.group(1)
        info = gsm_info.get(gsm)
        if info is None:
            skipped.append((h5_path.name, "no_manifest"))
            continue
        patient_id = info["patient_id"]
        timepoint = str(info["timepoint"]).strip().lower()
        sample_id = f"gse301741_{patient_id}_{timepoint}"

        ad = sc.read_10x_h5(str(h5_path), gex_only=False)
        ad.var_names_make_unique()
        # obs_names: keep original barcode in a column, then make unique per sample
        ad.obs["cell_barcode_original"] = ad.obs_names.astype(str)
        ad.obs_names = sample_id + "_" + ad.obs_names.astype(str)
        ad.obs["sample_id"] = sample_id
        ad.obs["patient_id"] = patient_id
        ad.obs["patient_key"] = f"GSE301741::{patient_id}"
        ad.obs["timepoint"] = timepoint
        ad.obs["geo_accession"] = gsm
        ad.obs["library_type"] = info["library_type"]
        ad.obs["cell_type_enrichment"] = info["cell_type"]

        # Merge frozen sample-level metadata
        sm = frozen_samples[frozen_samples["sample_id"] == sample_id]
        if not sm.empty:
            for col in [
                "tissue_source", "treatment_context", "treatment_raw",
                "response_raw", "response_endpoint_type", "response_binary_harmonized",
                "response_harmonization_confidence", "split_eligibility",
                "supervised_use_allowed", "support_use_allowed",
            ]:
                if col in sm.columns:
                    ad.obs[col] = sm[col].iloc[0]
        else:
            skipped.append((h5_path.name, "no_frozen_sample_match"))
            continue

        adatas.append(ad)
        print(f"Loaded {h5_path.name}: {ad.n_obs} cells -> {sample_id}")

    if skipped:
        print("\nSkipped files:")
        for name, reason in skipped:
            print(f"  {name}: {reason}")

    if not adatas:
        print("No adatas loaded", file=sys.stderr)
        sys.exit(1)

    print("\nConcatenating...")
    adata = sc.concat(adatas, join="inner")
    # obs_names are already unique because we prefixed them with sample_id.
    adata.obs["cohort_id"] = "GSE301741"

    # Ensure X is raw counts (integer sparse)
    if not hasattr(adata.X, "toarray"):
        # dense case
        pass

    # Store counts layer explicitly for downstream
    adata.layers["counts"] = adata.X.copy()

    # Add minimal patient metadata
    patient_cols = ["patient_key", "patient_id", "sex", "age", "cancer_type", "treatment_setting"]
    avail_patient_cols = [c for c in patient_cols if c in patient_meta.columns]
    pat = patient_meta[patient_meta["cohort_id"] == "GSE301741"][avail_patient_cols].copy()
    if not pat.empty:
        adata.obs = adata.obs.merge(
            pat, on="patient_key", how="left", suffixes=("", "_patient")
        )

    OUT_H5AD.parent.mkdir(parents=True, exist_ok=True)
    print(f"\nWriting {OUT_H5AD} ({adata.n_obs} cells x {adata.n_vars} genes)")
    adata.write(OUT_H5AD, compression="gzip")
    print("Done.")


if __name__ == "__main__":
    main()
