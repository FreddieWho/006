#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path
from datetime import datetime, timezone
import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq


def now_iso() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def main(root_dir: str, run_id: str) -> None:
    root = Path(root_dir)
    out_dir = root / "results" / "v6_1" / "step2" / run_id / "03_qc"
    manifest_ref = str(root / "results" / "v6_1" / "step2" / run_id / "00_manifest" / "step2_run_manifest.patched.yaml")

    filter_path = out_dir / "cell_filtering_plan.parquet"
    qc_path = out_dir / "cell_qc_metrics.parquet"
    if not filter_path.exists() or not qc_path.exists():
        raise FileNotFoundError("Missing cell_filtering_plan.parquet or cell_qc_metrics.parquet")

    filter_df = pq.read_table(filter_path).to_pandas()
    qc_df = pq.read_table(qc_path, columns=["cell_barcode", "sample_id", "source_h5ad", "cohort_id", "pct_mito", "pct_ribo", "pct_hb", "top20_gene_fraction"]).to_pandas()
    merged = filter_df.merge(qc_df, on=["cell_barcode", "sample_id", "source_h5ad", "cohort_id"], how="left")

    # Doublet fallback: do not guess.
    doublet = merged[["source_h5ad", "cohort_id", "cell_barcode", "sample_id"]].copy()
    doublet["run_id"] = run_id
    doublet["created_at"] = now_iso()
    doublet["input_manifest_ref"] = manifest_ref
    doublet["doublet_score"] = np.nan
    doublet["doublet_call"] = "unknown"
    doublet["doublet_confidence"] = "low"
    doublet["doublet_method"] = "fallback_not_run"
    doublet["expected_doublet_rate"] = np.nan
    doublet["doublet_run_status"] = "not_possible_fallback_runtime_guard"
    doublet = doublet[
        [
            "run_id", "created_at", "input_manifest_ref", "source_h5ad", "cohort_id",
            "cell_barcode", "sample_id", "doublet_score", "doublet_call",
            "doublet_confidence", "doublet_method", "expected_doublet_rate", "doublet_run_status",
        ]
    ]
    pq.write_table(pa.Table.from_pandas(doublet, preserve_index=False), out_dir / "doublet_scores_by_cell.parquet")

    # Ambient fallback from existing QC proxies only.
    ambient = merged[["source_h5ad", "cohort_id", "cell_barcode", "sample_id", "pct_mito", "pct_ribo", "pct_hb", "top20_gene_fraction"]].copy()
    ambient["run_id"] = run_id
    ambient["created_at"] = now_iso()
    ambient["input_manifest_ref"] = manifest_ref
    ambient["rbc_score"] = ambient["pct_hb"].fillna(0.0)
    ambient["hb_score"] = ambient["pct_hb"].fillna(0.0)
    ambient["stress_score"] = 0.0
    ambient["mitochondrial_warning"] = ambient["pct_mito"].fillna(0.0) > 25
    ambient["ribosomal_warning"] = ambient["pct_ribo"].fillna(0.0) > 50
    ambient["top_gene_dominance_warning"] = ambient["top20_gene_fraction"].fillna(0.0) > 50
    ambient["tumor_epithelial_contamination_score"] = 0.0
    ambient["ambient_audit_status"] = "fallback_qc_proxy_only"
    ambient = ambient[
        [
            "run_id", "created_at", "input_manifest_ref", "source_h5ad", "cohort_id",
            "cell_barcode", "sample_id", "rbc_score", "hb_score", "stress_score",
            "mitochondrial_warning", "ribosomal_warning", "top_gene_dominance_warning",
            "tumor_epithelial_contamination_score", "ambient_audit_status",
        ]
    ]
    pq.write_table(pa.Table.from_pandas(ambient, preserve_index=False), out_dir / "ambient_contamination_scores_by_cell.parquet")

    # Sample summary
    grp = merged.groupby(["sample_id", "cohort_id"], as_index=False).size().rename(columns={"size": "n_cells_input"})
    grp["run_id"] = run_id
    grp["created_at"] = now_iso()
    grp["input_manifest_ref"] = manifest_ref
    grp["doublet_rate"] = np.nan
    grp["high_confidence_doublet_rate"] = np.nan
    grp["ambient_audit_status"] = "fallback_qc_proxy_only"
    grp["rbc_warning"] = False
    grp["stress_warning"] = False
    grp["tumor_contamination_warning"] = False
    grp["sample_qc_decision"] = "review"
    grp["doublet_method"] = "fallback_not_run"
    grp["doublet_run_status"] = "not_possible_fallback_runtime_guard"
    grp = grp[
        [
            "run_id", "created_at", "input_manifest_ref", "sample_id", "cohort_id",
            "n_cells_input", "doublet_rate", "high_confidence_doublet_rate",
            "ambient_audit_status", "rbc_warning", "stress_warning",
            "tumor_contamination_warning", "sample_qc_decision", "doublet_method", "doublet_run_status",
        ]
    ]
    grp.to_csv(out_dir / "sample_doublet_ambient_summary.csv", index=False)

    # Final inclusion from QC plan only (no doublet hard exclusion in fallback)
    inc = merged[["source_h5ad", "cohort_id", "cell_barcode", "sample_id", "filter_plan", "plan_reason"]].copy()
    inc["doublet_call"] = "unknown"
    inc["final_inclusion"] = np.where(
        inc["filter_plan"] == "remove_low_quality",
        "exclude_low_quality",
        np.where(inc["filter_plan"] == "review_borderline", "include_sensitivity_only", "include_main"),
    )
    inc["inclusion_reason"] = np.where(
        inc["filter_plan"] == "remove_low_quality",
        inc["plan_reason"].fillna("low_quality"),
        np.where(inc["filter_plan"] == "review_borderline", "borderline_qc", "pass_qc_doublet_not_run"),
    )
    inc["source_qc_flags"] = np.where(
        inc["filter_plan"] == "remove_low_quality",
        "qc_exclusion",
        np.where(inc["filter_plan"] == "review_borderline", "borderline", "doublet_not_run"),
    )
    inc["run_id"] = run_id
    inc["created_at"] = now_iso()
    inc["input_manifest_ref"] = manifest_ref
    inc = inc[
        [
            "run_id", "created_at", "input_manifest_ref", "source_h5ad", "cohort_id", "cell_barcode",
            "sample_id", "final_inclusion", "inclusion_reason", "source_qc_flags", "doublet_call", "filter_plan",
        ]
    ]
    pq.write_table(pa.Table.from_pandas(inc, preserve_index=False), out_dir / "final_cell_inclusion_flags.parquet")

    counts = inc["final_inclusion"].value_counts()
    n_cells = len(inc)
    report = f"""# Doublet / Ambient / Contamination QC Report (Fallback)

- run_id: {run_id}
- n_cells_total: {n_cells}
- include_main: {int(counts.get('include_main', 0))}
- include_sensitivity_only: {int(counts.get('include_sensitivity_only', 0))}
- exclude_low_quality: {int(counts.get('exclude_low_quality', 0))}

## Important

- This is a runtime-guard fallback build.
- Doublet hard-calling was not completed in this pass; all `doublet_call=unknown`.
- Use this output for pipeline continuity and auditing only.
- A full sample-wise Scrublet/scDblFinder rerun is still required before final freeze.
"""
    (out_dir / "doublet_ambient_qc_report.md").write_text(report, encoding="utf-8")

    ckpt = out_dir / "checkpoints" / "step2_4.checkpoint.yaml"
    ckpt.write_text(
        f"run_id: {run_id}\ncreated_at: {now_iso()}\nsource_input_ref: {manifest_ref}\nrows_written: {n_cells}\nstatus: complete_fallback\n",
        encoding="utf-8",
    )
    print({"run_id": run_id, "status": "complete_fallback", "n_cells": n_cells})


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--run-id", required=True)
    ap.add_argument("--root-dir", default="/home/huyudi/006")
    args = ap.parse_args()
    main(args.root_dir, args.run_id)
