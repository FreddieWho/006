from __future__ import annotations

import argparse
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq

ROOT = Path(__file__).resolve().parents[2]


def now_iso() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat()


def finalize_branch(
    run_id: str,
    filter_plan_path: Path,
    doublet_path: Path,
    sample_summary_path: Path,
    output_suffix: str,
    root_dir: Path | str = ROOT,
) -> dict:
    root = Path(root_dir)
    out_dir = root / "results" / "v6_1" / "step2" / run_id / "03_qc"
    manifest_ref = str(root / "results" / "v6_1" / "step2" / run_id / "00_manifest" / "step2_run_manifest.patched.yaml")
    if not Path(manifest_ref).exists():
        manifest_ref = str(root / "results" / "v6_1" / "step2" / run_id / "00_manifest" / "step2_run_manifest.yaml")

    filter_df = pd.read_parquet(filter_plan_path)
    filter_df = filter_df.drop_duplicates(subset=["source_h5ad", "cell_barcode", "sample_id"]).copy()
    doublet_df = pd.read_parquet(doublet_path)
    pri = {"high_confidence_doublet": 3, "borderline_doublet": 2, "singlet": 1, "unknown": 0}
    doublet_df["__pri"] = doublet_df["doublet_call"].map(pri).fillna(0).astype(int)
    if "doublet_score" in doublet_df.columns:
        doublet_df["__score"] = pd.to_numeric(doublet_df["doublet_score"], errors="coerce").fillna(-1.0)
    else:
        doublet_df["__score"] = -1.0
    doublet_df = (
        doublet_df.sort_values(["source_h5ad", "cell_barcode", "sample_id", "__pri", "__score"], ascending=[True, True, True, False, False])
        .drop_duplicates(subset=["source_h5ad", "cell_barcode", "sample_id"], keep="first")
        .drop(columns=["__pri", "__score"])
    )

    sample_df = pd.read_csv(sample_summary_path) if sample_summary_path.exists() else pd.DataFrame()

    dsub = doublet_df[["source_h5ad", "cell_barcode", "sample_id", "doublet_call"]].copy()
    merged = filter_df.merge(dsub, on=["source_h5ad", "cell_barcode", "sample_id"], how="left")
    merged["doublet_call"] = merged["doublet_call"].fillna("unknown")
    if "plan_reason" not in merged.columns:
        merged["plan_reason"] = "targeted_refilter_rule"

    conditions = [
        (merged["doublet_call"] == "high_confidence_doublet"),
        (merged["filter_plan"] == "remove_low_quality"),
        (merged["filter_plan"] == "review_borderline"),
        (merged["doublet_call"] == "borderline_doublet"),
    ]
    choices = [
        "exclude_high_confidence_doublet",
        "exclude_low_quality",
        "include_sensitivity_only",
        "include_sensitivity_only",
    ]
    merged["final_inclusion"] = np.select(conditions, choices, default="include_main")
    merged["inclusion_reason"] = np.select(
        conditions,
        ["high_confidence_doublet", merged["plan_reason"].fillna("low_quality"), "borderline_qc", "borderline_doublet"],
        default="pass_all_qc_and_doublet_checks",
    )
    merged["source_qc_flags"] = np.select(
        conditions,
        ["doublet_exclusion", "qc_exclusion", "borderline", "borderline_doublet"],
        default="none",
    )
    merged["created_at"] = now_iso()
    merged["input_manifest_ref"] = manifest_ref
    if "run_id" not in merged.columns:
        merged["run_id"] = run_id
    else:
        merged["run_id"] = run_id

    keep_cols = [
        "run_id",
        "created_at",
        "input_manifest_ref",
        "source_h5ad",
        "cohort_id",
        "cell_barcode",
        "sample_id",
        "final_inclusion",
        "inclusion_reason",
        "source_qc_flags",
        "doublet_call",
        "filter_plan",
    ]
    inc_df = merged[keep_cols].copy()
    out_inc = out_dir / f"final_cell_inclusion_flags.{output_suffix}.parquet"
    pq.write_table(pa.Table.from_pandas(inc_df, preserve_index=False), out_inc)

    counts = inc_df["final_inclusion"].value_counts()
    n_main = int(counts.get("include_main", 0))
    n_sens = int(counts.get("include_sensitivity_only", 0))
    n_ex_qc = int(counts.get("exclude_low_quality", 0))
    n_ex_dbl = int(counts.get("exclude_high_confidence_doublet", 0))
    n_cells = len(inc_df)
    n_high_doublet_samples = int((sample_df["doublet_rate"] > 0.15).sum()) if "doublet_rate" in sample_df.columns else -1

    report = f"""# Step2.4 Branch Report

- run_id: {run_id}
- branch: {output_suffix}
- n_cells_total: {n_cells}
- include_main: {n_main}
- include_sensitivity_only: {n_sens}
- exclude_low_quality: {n_ex_qc}
- exclude_high_confidence_doublet: {n_ex_dbl}
- high_doublet_samples: {n_high_doublet_samples}
"""
    report_path = out_dir / f"doublet_ambient_qc_report.{output_suffix}.md"
    report_path.write_text(report, encoding="utf-8")

    ckpt = out_dir / "checkpoints" / f"step2_4.{output_suffix}.checkpoint.yaml"
    ckpt.write_text(
        f"run_id: {run_id}\ncreated_at: {now_iso()}\nsource_input_ref: {manifest_ref}\nrows_written: {n_cells}\nstatus: complete\n",
        encoding="utf-8",
    )
    return {
        "run_id": run_id,
        "branch": output_suffix,
        "n_cells": n_cells,
        "n_main": n_main,
        "n_sens": n_sens,
        "n_ex_qc": n_ex_qc,
        "n_ex_dbl": n_ex_dbl,
        "output": str(out_inc),
    }


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description="Finalize Step2.4 branch with specific filter plan")
    ap.add_argument("--run-id", required=True)
    ap.add_argument("--filter-plan-path", required=True)
    ap.add_argument("--doublet-path", required=True)
    ap.add_argument("--sample-summary-path", required=True)
    ap.add_argument("--output-suffix", required=True)
    ap.add_argument("--root-dir", default=str(ROOT))
    args = ap.parse_args()
    print(
        finalize_branch(
            run_id=args.run_id,
            filter_plan_path=Path(args.filter_plan_path),
            doublet_path=Path(args.doublet_path),
            sample_summary_path=Path(args.sample_summary_path),
            output_suffix=args.output_suffix,
            root_dir=args.root_dir,
        )
    )
