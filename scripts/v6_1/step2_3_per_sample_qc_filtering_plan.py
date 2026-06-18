from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import anndata as ad
import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
import scanpy as sc

ROOT = Path(__file__).resolve().parents[2]


def now_iso() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat()


def csv_write(path: Path, fieldnames: List[str], rows: List[Dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def load_patched_manifest(run_root: Path) -> Dict[str, str]:
    patched = run_root / "00_manifest" / "step2_run_manifest.patched.yaml"
    orig = run_root / "00_manifest" / "step2_run_manifest.yaml"
    path = patched if patched.exists() else orig
    result: Dict[str, str] = {}
    with path.open("r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            if ":" in line:
                key, val = line.split(":", 1)
                result[key.strip()] = val.strip()
    return result


def check_gate(run_root: Path) -> bool:
    gate_path = run_root / "00_manifest" / "step_specific_gate_status.csv"
    if not gate_path.exists():
        return False
    with gate_path.open("r", encoding="utf-8") as f:
        for row in csv.DictReader(f):
            if row.get("step_name") == "Step2.3_qc":
                return row.get("status") == "ready"
    return False


def compute_mad(vals: np.ndarray) -> float:
    """Median absolute deviation."""
    median = float(np.median(vals))
    return float(np.median(np.abs(vals - median)))


def subset_to_memory(adata: ad.AnnData, obs_idx: np.ndarray) -> ad.AnnData:
    """Subset backed AnnData and return in-memory copy, handling missing X layer."""
    sub = adata[obs_idx]
    try:
        # Try standard to_memory first
        return sub.to_memory()
    except Exception:
        # X layer missing — reconstruct from layers
        layer_key = None
        for key in ("counts", "data"):
            if key in sub.layers:
                layer_key = key
                break
        if layer_key is None:
            raise RuntimeError(f"No X layer and no usable layer (counts/data) found in {adata.filename}")
        X = sub.layers[layer_key]
        if hasattr(X, "to_memory"):
            X = X.to_memory()
        return ad.AnnData(
            X=X,
            obs=sub.obs.copy(),
            var=sub.var.copy(),
        )


def select_counts_anndata(adata: ad.AnnData) -> Tuple[ad.AnnData, str]:
    """Return an AnnData with counts-like matrix in X."""
    if "counts" in adata.layers:
        return ad.AnnData(X=adata.layers["counts"], obs=adata.obs.copy(), var=adata.var.copy()), "layers[counts]"
    if adata.raw is not None and getattr(adata.raw, "X", None) is not None:
        return ad.AnnData(X=adata.raw.X, obs=adata.obs.copy(), var=adata.raw.var.copy()), "raw.X"
    return adata, "X"


def harmonize_effective_sample_ids(mapped_df: pd.DataFrame) -> Tuple[pd.DataFrame, List[Dict[str, object]]]:
    """Collapse pathological cell-level sample proxies to a pooled source sample."""
    df = mapped_df.copy()
    rules: List[Dict[str, object]] = []
    df["effective_sample_id"] = df["mapped_sample_id"].astype(str)
    df["sample_granularity_rule"] = "as_mapped"

    for source_name, g in df.groupby("source_h5ad"):
        n_cells = len(g)
        if n_cells == 0:
            continue
        n_samples = g["mapped_sample_id"].astype(str).nunique(dropna=False)
        unknown_patient_frac = float((g["mapped_patient_id"].astype(str) == "unknown").mean())
        sample_to_cell_ratio = n_samples / n_cells
        # Heuristic: likely cell-level pseudo-sample labels, not true biological samples.
        if sample_to_cell_ratio > 0.8 and unknown_patient_frac > 0.95:
            cohort_id = str(g["mapped_cohort_id"].iloc[0])
            pooled_id = f"{cohort_id}__source_pooled"
            mask = df["source_h5ad"] == source_name
            df.loc[mask, "effective_sample_id"] = pooled_id
            df.loc[mask, "sample_granularity_rule"] = "pooled_source_due_cell_level_sample_proxy"
            rules.append({
                "source_h5ad": source_name,
                "cohort_id": cohort_id,
                "n_cells": n_cells,
                "n_samples_before": n_samples,
                "n_samples_after": 1,
                "sample_to_cell_ratio": round(sample_to_cell_ratio, 6),
                "unknown_patient_frac": round(unknown_patient_frac, 6),
                "rule": "pooled_source_due_cell_level_sample_proxy",
            })
    return df, rules


def compute_cell_qc(adata: ad.AnnData) -> pd.DataFrame:
    """Compute QC metrics for all cells in an AnnData object using scanpy."""
    var_names = [str(v).upper() for v in adata.var_names]

    # Set up qc_vars for scanpy
    adata.var["mito"] = [g.startswith("MT-") for g in var_names]
    adata.var["ribo"] = [g.startswith(("RPL", "RPS")) for g in var_names]
    adata.var["hb"] = [g.startswith("HB") for g in var_names]
    adata.var["malat1"] = [g == "MALAT1" for g in var_names]

    sc.pp.calculate_qc_metrics(
        adata,
        qc_vars=["mito", "ribo", "hb", "malat1"],
        inplace=True,
    )

    df = pd.DataFrame({
        "cell_barcode": adata.obs_names.astype(str),
        "total_counts": adata.obs["total_counts"].values,
        "n_genes_by_counts": adata.obs["n_genes_by_counts"].values,
        "pct_mito": adata.obs["pct_counts_mito"].values,
        "pct_ribo": adata.obs["pct_counts_ribo"].values,
        "pct_hb": adata.obs["pct_counts_hb"].values,
        "malat1_fraction": adata.obs["pct_counts_malat1"].values,
        "top20_gene_fraction": np.nan,
    })
    return df


def run_step2_3_qc(
    root_dir: Path | str = ROOT,
    run_id: str | None = None,
) -> Dict[str, object]:
    root = Path(root_dir)
    if run_id is None:
        now = datetime.now(timezone.utc)
        run_id = f"step2_v6_1_{now.strftime('%m%d_%H%M')}"

    run_root = root / "results" / "v6_1" / "step2" / run_id
    out_dir = run_root / "03_qc"
    out_dir.mkdir(parents=True, exist_ok=True)
    checkpoint_dir = out_dir / "checkpoints"
    checkpoint_dir.mkdir(parents=True, exist_ok=True)
    fig_dir = out_dir / "qc_figures"
    fig_dir.mkdir(parents=True, exist_ok=True)

    manifest = load_patched_manifest(run_root)
    manifest_ref = str(run_root / "00_manifest" / "step2_run_manifest.patched.yaml")
    if not Path(manifest_ref).exists():
        manifest_ref = str(run_root / "00_manifest" / "step2_run_manifest.yaml")

    gate_ready = check_gate(run_root)
    if not gate_ready:
        print("Gate not ready for Step2.3; writing degraded note.")
        (out_dir / "qc_report.md").write_text("# QC Report\n\nGate not ready.\n", encoding="utf-8")
        return {"run_id": run_id, "status": "gate_not_ready"}

    # Read mapping parquet
    mapping_path = run_root / "01_mapping_audit" / "cell_to_sample_mapping_draft.parquet"
    if not mapping_path.exists():
        print("Mapping parquet not found.")
        return {"run_id": run_id, "status": "missing_mapping"}

    print("Reading mapping parquet ...")
    mapping_table = pq.read_table(
        mapping_path,
        columns=[
            "cell_barcode",
            "source_h5ad",
            "mapped_cohort_id",
            "mapped_sample_id",
            "mapped_patient_id",
            "mapping_status",
            "mapping_reason",
        ],
    )
    mapping_df = mapping_table.to_pandas()
    mapped_df = mapping_df[mapping_df["mapping_status"] == "mapped"].copy()
    # Deduplicate: some h5ad files have duplicate obs indices causing duplicate mapping rows
    before_dedup = len(mapped_df)
    mapped_df = mapped_df.drop_duplicates(subset=["cell_barcode", "source_h5ad", "mapped_sample_id"]).copy()
    after_dedup = len(mapped_df)
    if before_dedup != after_dedup:
        print(f"Deduplicated: {before_dedup} -> {after_dedup} rows")
    print(f"Total mapped cells: {len(mapped_df)}")
    mapped_df, granularity_rules = harmonize_effective_sample_ids(mapped_df)
    if granularity_rules:
        print(f"Applied sample granularity harmonization rules to {len(granularity_rules)} sources")

    # Group by source
    source_groups = mapped_df.groupby("source_h5ad")

    # Process each source into chunk files for resumable/recoverable execution.
    total_qc_rows = 0
    chunk_dir = out_dir / "chunks"
    chunk_dir.mkdir(parents=True, exist_ok=True)
    chunk_paths = sorted(chunk_dir.glob("cell_qc_*.parquet"))
    reuse_existing_chunks = len(chunk_paths) > 0
    cell_qc_path = out_dir / "cell_qc_metrics.parquet"
    if cell_qc_path.exists():
        cell_qc_path.unlink()

    if reuse_existing_chunks:
        print(f"Reusing existing QC chunks: n_files={len(chunk_paths)}")
        for chunk_path in chunk_paths:
            total_qc_rows += pq.ParquetFile(chunk_path).metadata.num_rows
    else:
        h5ad_dir = root / "data" / "processed" / "srt" / "raw"
        for idx, (source_name, group) in enumerate(source_groups):
            path = h5ad_dir / source_name
            if not path.exists():
                print(f"  [{idx+1}/{len(source_groups)}] {source_name} not found, skipping.")
                continue

            print(f"  [{idx+1}/{len(source_groups)}] Processing {source_name} ({len(group)} cells) ...")

            try:
                # Read h5ad
                adata = ad.read_h5ad(path, backed="r")

                # Subset to mapped cells only
                cell_barcodes = group["cell_barcode"].tolist()
                obs_idx = adata.obs_names.isin(cell_barcodes)
                if not obs_idx.any():
                    print(f"    No matching cells found.")
                    adata.file.close()
                    continue

                # Read subset into memory
                adata_sub = subset_to_memory(adata, obs_idx)
                adata.file.close()

                # Compute QC on counts-like matrix (not normalized X when avoidable).
                adata_qc, matrix_layer_used = select_counts_anndata(adata_sub)
                qc_df = compute_cell_qc(adata_qc)
                qc_df["source_h5ad"] = source_name
                qc_df["cohort_id"] = group["mapped_cohort_id"].iloc[0]
                qc_df["sample_id"] = group["effective_sample_id"].values
                qc_df["sample_granularity_rule"] = group["sample_granularity_rule"].values
                qc_df["original_mapped_sample_id"] = group["mapped_sample_id"].values

                src_df = pd.DataFrame({
                    "run_id": run_id,
                    "created_at": now_iso(),
                    "input_manifest_ref": manifest_ref,
                    "source_h5ad": source_name,
                    "cohort_id": qc_df["cohort_id"].astype(str).values,
                    "cell_barcode": qc_df["cell_barcode"].astype(str).values,
                    "sample_id": qc_df["sample_id"].astype(str).values,
                    "total_counts": qc_df["total_counts"].astype(float).values,
                    "n_genes_by_counts": qc_df["n_genes_by_counts"].astype(int).values,
                    "pct_mito": qc_df["pct_mito"].astype(float).values,
                    "pct_ribo": qc_df["pct_ribo"].astype(float).values,
                    "pct_hb": qc_df["pct_hb"].astype(float).values,
                    "malat1_fraction": qc_df["malat1_fraction"].astype(float).values,
                    "top20_gene_fraction": qc_df["top20_gene_fraction"].astype(float).values,
                    "matrix_layer_used": matrix_layer_used,
                    "matrix_scale_uncertain": False,
                    "sample_granularity_rule": qc_df["sample_granularity_rule"].astype(str).values,
                    "original_mapped_sample_id": qc_df["original_mapped_sample_id"].astype(str).values,
                })
                table = pa.Table.from_pandas(src_df, preserve_index=False)
                chunk_path = chunk_dir / f"cell_qc_{idx+1:03d}_{source_name.replace('.h5ad','')}.parquet"
                pq.write_table(table, chunk_path)
                chunk_paths.append(chunk_path)
                total_qc_rows += len(src_df)
                print(f"    -> QC rows: {len(qc_df)}")

            except Exception as e:
                print(f"    ERROR: {e}")
                continue

    print(f"Total QC rows: {total_qc_rows}")

    if total_qc_rows == 0:
        print("No QC rows generated.")
        return {"run_id": run_id, "status": "no_qc_rows"}

    # Materialize a single canonical cell_qc_metrics.parquet from chunks.
    final_writer = None
    for chunk_path in sorted(chunk_paths):
        table = pq.read_table(chunk_path)
        if final_writer is None:
            final_writer = pq.ParquetWriter(cell_qc_path, table.schema)
        final_writer.write_table(table)
    if final_writer is not None:
        final_writer.close()

    print("Computing per-sample thresholds from cell_qc_metrics.parquet in chunks ...")
    sample_acc: Dict[str, Dict[str, object]] = {}
    for chunk_path in sorted(chunk_paths):
        pf = pq.ParquetFile(chunk_path)
        for batch in pf.iter_batches(batch_size=500_000):
            bdf = batch.to_pandas()
            for sample_id, sdf in bdf.groupby("sample_id"):
                acc = sample_acc.setdefault(sample_id, {
                    "cohort_id": str(sdf["cohort_id"].iloc[0]),
                    "source_h5ad": str(sdf["source_h5ad"].iloc[0]),
                    "sample_granularity_rule": str(sdf["sample_granularity_rule"].iloc[0]),
                    "counts": [],
                    "genes": [],
                    "mito": [],
                    "ribo": [],
                    "hb": [],
                })
                acc["counts"].append(sdf["total_counts"].to_numpy(dtype=float))
                acc["genes"].append(sdf["n_genes_by_counts"].to_numpy(dtype=float))
                acc["mito"].append(sdf["pct_mito"].to_numpy(dtype=float))
                acc["ribo"].append(sdf["pct_ribo"].to_numpy(dtype=float))
                acc["hb"].append(sdf["pct_hb"].to_numpy(dtype=float))

    sample_metrics = []
    threshold_rows = []
    thresh_lookup: Dict[str, Dict[str, float]] = {}

    for sample_id, acc in sample_acc.items():
        counts = np.concatenate(acc["counts"])
        genes = np.concatenate(acc["genes"])
        mito = np.concatenate(acc["mito"])
        ribo = np.concatenate(acc["ribo"])
        hb = np.concatenate(acc["hb"])
        n_cells = len(counts)

        median_counts = float(np.median(counts))
        mad_counts = compute_mad(counts)
        median_genes = float(np.median(genes))
        mad_genes = compute_mad(genes)
        median_mito = float(np.median(mito))
        mad_mito = compute_mad(mito)

        min_genes = 200
        min_counts = 500
        max_genes = median_genes + 5 * mad_genes
        max_counts = median_counts + 5 * mad_counts
        max_pct_mito = min(25.0, median_mito + 5 * mad_mito)
        if max_genes < min_genes * 2:
            max_genes = min_genes * 2
        if max_counts < min_counts * 2:
            max_counts = min_counts * 2

        keep_mask = (
            (genes >= min_genes) &
            (counts >= min_counts) &
            (genes <= max_genes) &
            (counts <= max_counts) &
            (mito <= max_pct_mito)
        )
        n_keep = int(keep_mask.sum())

        flags = []
        if n_cells < 200:
            flags.append("sample_too_small_blocking")
        if n_keep < 200:
            flags.append("post_filter_too_small_blocking")

        sample_metrics.append({
            "run_id": run_id,
            "created_at": now_iso(),
            "input_manifest_ref": manifest_ref,
            "cohort_id": acc["cohort_id"],
            "sample_id": sample_id,
            "source_h5ad": acc["source_h5ad"],
            "n_cells_raw": n_cells,
            "median_counts": median_counts,
            "median_genes": median_genes,
            "median_pct_mito": median_mito,
            "median_pct_ribo": float(np.median(ribo)),
            "median_pct_hb": float(np.median(hb)),
            "library_size": float(counts.sum()),
            "detected_genes": int(np.unique(genes).size),
            "n_cells_keep_expected": n_keep,
            "immune_cells_keep_expected": "",
            "immune_support_method": "",
            "sample_flags": ";".join(flags) if flags else "none",
            "matrix_scale_uncertain": False,
            "sample_granularity_rule": acc["sample_granularity_rule"],
        })

        thresh_lookup[sample_id] = {
            "min_genes": float(min_genes),
            "min_counts": float(min_counts),
            "max_genes": float(max_genes),
            "max_counts": float(max_counts),
            "max_pct_mito": float(max_pct_mito),
        }
        threshold_rows.append({
            "run_id": run_id,
            "created_at": now_iso(),
            "input_manifest_ref": manifest_ref,
            "sample_id": sample_id,
            "cohort_id": acc["cohort_id"],
            "min_genes": min_genes,
            "min_counts": min_counts,
            "max_genes": round(max_genes, 1),
            "max_counts": round(max_counts, 1),
            "max_pct_mito": round(max_pct_mito, 2),
            "min_genes_method": "fixed_default_200",
            "min_counts_method": "fixed_default_500",
            "max_genes_method": "sample_median_plus_5mad",
            "max_counts_method": "sample_median_plus_5mad",
            "max_pct_mito_method": "min_25pct_or_sample_median_plus_5mad",
            "sample_median_counts": round(median_counts, 1),
            "sample_mad_counts": round(mad_counts, 1),
            "sample_median_genes": round(median_genes, 1),
            "sample_mad_genes": round(mad_genes, 1),
            "sample_median_pct_mito": round(median_mito, 4),
            "source_h5ad": acc["source_h5ad"],
            "n_cells_raw": n_cells,
            "threshold_method_summary": "per_sample_fixed_mins_plus_median_5mad_plus_mito_cap",
            "sample_granularity_rule": acc["sample_granularity_rule"],
        })

    # Write sample_qc_summary.csv
    csv_write(
        out_dir / "sample_qc_summary.csv",
        [
            "run_id", "created_at", "input_manifest_ref", "cohort_id", "sample_id",
            "source_h5ad", "n_cells_raw", "median_counts", "median_genes",
            "median_pct_mito", "median_pct_ribo", "median_pct_hb", "library_size",
            "detected_genes", "n_cells_keep_expected", "immune_cells_keep_expected",
            "immune_support_method", "sample_flags", "matrix_scale_uncertain",
            "sample_granularity_rule",
        ],
        sample_metrics,
    )

    # Write per_sample_qc_thresholds.csv
    csv_write(
        out_dir / "per_sample_qc_thresholds.csv",
        [
            "run_id", "created_at", "input_manifest_ref", "sample_id", "cohort_id",
            "min_genes", "min_counts", "max_genes", "max_counts", "max_pct_mito",
            "min_genes_method", "min_counts_method", "max_genes_method",
            "max_counts_method", "max_pct_mito_method", "sample_median_counts",
            "sample_mad_counts", "sample_median_genes", "sample_mad_genes",
            "sample_median_pct_mito", "source_h5ad", "n_cells_raw",
            "threshold_method_summary", "sample_granularity_rule",
        ],
        threshold_rows,
    )

    # Build filtering plan in chunks
    filtering_path = out_dir / "cell_filtering_plan.parquet"
    if filtering_path.exists():
        filtering_path.unlink()
    filter_writer = None
    n_keep = 0
    n_remove = 0
    n_review = 0
    for chunk_path in sorted(chunk_paths):
        pf = pq.ParquetFile(chunk_path)
        for batch in pf.iter_batches(batch_size=500_000):
            bdf = batch.to_pandas()
            out = pd.DataFrame({
                "run_id": run_id,
                "created_at": now_iso(),
                "input_manifest_ref": manifest_ref,
                "source_h5ad": bdf["source_h5ad"].astype(str),
                "cohort_id": bdf["cohort_id"].astype(str),
                "cell_barcode": bdf["cell_barcode"].astype(str),
                "sample_id": bdf["sample_id"].astype(str),
                "filter_plan": "keep_high_confidence",
                "plan_reason": "within_all_thresholds",
            })
            sample_ids = bdf["sample_id"].astype(str)
            min_genes = sample_ids.map(lambda s: thresh_lookup.get(s, {}).get("min_genes", np.nan)).to_numpy()
            min_counts = sample_ids.map(lambda s: thresh_lookup.get(s, {}).get("min_counts", np.nan)).to_numpy()
            max_genes = sample_ids.map(lambda s: thresh_lookup.get(s, {}).get("max_genes", np.nan)).to_numpy()
            max_counts = sample_ids.map(lambda s: thresh_lookup.get(s, {}).get("max_counts", np.nan)).to_numpy()
            max_mito = sample_ids.map(lambda s: thresh_lookup.get(s, {}).get("max_pct_mito", np.nan)).to_numpy()

            genes = bdf["n_genes_by_counts"].to_numpy(dtype=float)
            counts = bdf["total_counts"].to_numpy(dtype=float)
            mito = bdf["pct_mito"].to_numpy(dtype=float)

            missing_mask = np.isnan(min_genes)
            remove_mask = (
                (genes < min_genes) |
                (counts < min_counts) |
                (genes > max_genes) |
                (counts > max_counts) |
                (mito > max_mito)
            ) & (~missing_mask)
            out.loc[missing_mask, "filter_plan"] = "review_borderline"
            out.loc[missing_mask, "plan_reason"] = "missing_threshold"
            out.loc[remove_mask, "filter_plan"] = "remove_low_quality"
            out.loc[remove_mask, "plan_reason"] = "outside_thresholds"

            remaining = ~(missing_mask | remove_mask)
            borderline_mask = (
                (((genes < (min_genes * 1.1)) & (genes >= min_genes)) |
                 ((counts < (min_counts * 1.1)) & (counts >= min_counts)) |
                 ((genes > (max_genes * 0.9)) & (genes <= max_genes)) |
                 ((counts > (max_counts * 0.9)) & (counts <= max_counts)) |
                 ((mito > (max_mito * 0.9)) & (mito <= max_mito)))
                & remaining
            )
            out.loc[borderline_mask, "filter_plan"] = "review_borderline"
            out.loc[borderline_mask, "plan_reason"] = "near_threshold_edge"

            n_keep += int((out["filter_plan"] == "keep_high_confidence").sum())
            n_remove += int((out["filter_plan"] == "remove_low_quality").sum())
            n_review += int((out["filter_plan"] == "review_borderline").sum())

            table = pa.Table.from_pandas(out, preserve_index=False)
            if filter_writer is None:
                filter_writer = pq.ParquetWriter(filtering_path, table.schema)
            filter_writer.write_table(table)
    if filter_writer is not None:
        filter_writer.close()
    print("Wrote cell_filtering_plan.parquet")

    # Write report
    n_samples = len(sample_metrics)
    n_cells_total = total_qc_rows
    n_blocking = sum(1 for s in sample_metrics if "sample_too_small_blocking" in s["sample_flags"] or "post_filter_too_small_blocking" in s["sample_flags"])

    granularity_block = ""
    if granularity_rules:
        granularity_block = "\n## Sample Granularity Harmonization\n\n"
        for rule in granularity_rules:
            granularity_block += (
                f"- {rule['source_h5ad']}: {rule['n_samples_before']} -> {rule['n_samples_after']} "
                f"({rule['rule']}, sample_to_cell_ratio={rule['sample_to_cell_ratio']}, "
                f"unknown_patient_frac={rule['unknown_patient_frac']})\n"
            )

    report = f"""# QC Filtering Plan Report

- **run_id**: {run_id}
- **gate_ready**: {gate_ready}
- **n_sources**: {len(source_groups)}
- **n_samples**: {n_samples}
- **n_cells_total**: {n_cells_total}
- **n_keep_high_confidence**: {n_keep}
- **n_remove_low_quality**: {n_remove}
- **n_review_borderline**: {n_review}
- **n_blocking_samples**: {n_blocking}

## Threshold Policy

- min_genes = 200 (fixed)
- min_counts = 500 (fixed)
- max_genes = sample_median + 5 * MAD
- max_counts = sample_median + 5 * MAD
- max_pct_mito = min(25%, sample_median + 5 * MAD)

## Cell Filtering Summary

| Plan | Count | Fraction |
|------|-------|----------|
| keep_high_confidence | {n_keep} | {n_keep/n_cells_total:.4f} |
| remove_low_quality | {n_remove} | {n_remove/n_cells_total:.4f} |
| review_borderline | {n_review} | {n_review/n_cells_total:.4f} |

## Sample Blocking Summary

- Samples with blocking flags: {n_blocking}

## Key Parameter References

- Per-cell QC should be done before downstream interpretation and should use library complexity and mitochondrial proportion as core indicators:
  Luecken & Theis 2019 (Mol Syst Biol, doi:10.15252/msb.20188746)
- Multifeature low-quality-cell identification rationale:
  Ilicic et al. 2016 (Genome Biol, doi:10.1186/s13059-016-0888-1)

{granularity_block}
"""
    (out_dir / "qc_report.md").write_text(report, encoding="utf-8")

    # Checkpoint
    (checkpoint_dir / "step2_3.checkpoint.yaml").write_text(
        f"run_id: {run_id}\ncreated_at: {now_iso()}\nsource_input_ref: {manifest_ref}\n"
        f"rows_written: {total_qc_rows}\nstatus: complete\n",
        encoding="utf-8",
    )

    print(f"Done. samples={n_samples}, cells={n_cells_total}, keep={n_keep}, remove={n_remove}, review={n_review}")
    return {
        "run_id": run_id,
        "status": "complete",
        "n_samples": n_samples,
        "n_cells": n_cells_total,
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Step2.3 per-sample QC filtering plan")
    parser.add_argument("--run-id", required=True, help="run_id from Step2.0")
    parser.add_argument("--root-dir", default=str(ROOT), help="project root")
    args = parser.parse_args()
    result = run_step2_3_qc(root_dir=args.root_dir, run_id=args.run_id)
    print(result)
