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


def finalize_inclusion_flags(run_id: str, root_dir: Path | str = ROOT) -> dict:
    root = Path(root_dir)
    out_dir = root / "results" / "v6_1" / "step2" / run_id / "03_qc"

    manifest_ref = str(root / "results" / "v6_1" / "step2" / run_id / "00_manifest" / "step2_run_manifest.patched.yaml")
    if not Path(manifest_ref).exists():
        manifest_ref = str(root / "results" / "v6_1" / "step2" / run_id / "00_manifest" / "step2_run_manifest.yaml")

    print("Reading cell_filtering_plan.parquet ...")
    filter_df = pd.read_parquet(out_dir / "cell_filtering_plan.parquet")
    filter_df = filter_df.drop_duplicates(subset=["cell_barcode", "sample_id"]).copy()
    print(f"  -> {len(filter_df)} rows")

    print("Reading doublet_scores_by_cell.parquet ...")
    doublet_df = pd.read_parquet(out_dir / "doublet_scores_by_cell.parquet")
    # Deduplicate potential shard overlap rows before merge.
    # Priority: high_confidence_doublet > borderline_doublet > singlet > unknown.
    pri = {"high_confidence_doublet": 3, "borderline_doublet": 2, "singlet": 1, "unknown": 0}
    doublet_df["__pri"] = doublet_df["doublet_call"].map(pri).fillna(0).astype(int)
    if "doublet_score" in doublet_df.columns:
        doublet_df["__score"] = pd.to_numeric(doublet_df["doublet_score"], errors="coerce").fillna(-1.0)
    else:
        doublet_df["__score"] = -1.0
    doublet_df = (
        doublet_df.sort_values(["cell_barcode", "sample_id", "__pri", "__score"], ascending=[True, True, False, False])
        .drop_duplicates(subset=["cell_barcode", "sample_id"], keep="first")
        .drop(columns=["__pri", "__score"])
        .copy()
    )
    print(f"  -> {len(doublet_df)} rows")

    print("Reading sample_doublet_ambient_summary.csv ...")
    sample_df = pd.read_csv(out_dir / "sample_doublet_ambient_summary.csv")
    print(f"  -> {len(sample_df)} rows")

    # Merge doublet calls into filter plan
    print("Merging ...")
    doublet_subset = doublet_df[["cell_barcode", "sample_id", "doublet_call"]].copy()
    merged = filter_df.merge(
        doublet_subset,
        on=["cell_barcode", "sample_id"],
        how="left",
    )
    merged["doublet_call"] = merged["doublet_call"].fillna("unknown")

    # Vectorized inclusion logic
    print("Computing inclusion flags ...")
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

    reason_conditions = [
        (merged["doublet_call"] == "high_confidence_doublet"),
        (merged["filter_plan"] == "remove_low_quality"),
        (merged["filter_plan"] == "review_borderline"),
        (merged["doublet_call"] == "borderline_doublet"),
    ]
    reason_choices = [
        "high_confidence_doublet",
        merged["plan_reason"],
        "borderline_qc",
        "borderline_doublet",
    ]
    merged["inclusion_reason"] = np.select(reason_conditions, reason_choices, default="pass_all_qc_and_doublet_checks")

    # Flags
    flag_conditions = [
        (merged["doublet_call"] == "high_confidence_doublet"),
        (merged["filter_plan"] == "remove_low_quality"),
        (merged["filter_plan"] == "review_borderline"),
        (merged["doublet_call"] == "borderline_doublet"),
    ]
    flag_choices = [
        "doublet_exclusion",
        "qc_exclusion",
        "borderline",
        "borderline_doublet",
    ]
    merged["source_qc_flags"] = np.select(flag_conditions, flag_choices, default="none")

    # Build output
    inc_df = merged[[
        "run_id", "created_at", "input_manifest_ref", "source_h5ad", "cohort_id",
        "cell_barcode", "sample_id", "final_inclusion", "inclusion_reason",
        "source_qc_flags", "doublet_call", "filter_plan",
    ]].copy()
    inc_df["created_at"] = now_iso()
    inc_df["input_manifest_ref"] = manifest_ref

    # Write parquet
    print("Writing final_cell_inclusion_flags.parquet ...")
    table = pa.Table.from_pandas(inc_df, preserve_index=False)
    pq.write_table(table, out_dir / "final_cell_inclusion_flags.parquet")
    print("  -> done")

    # Report stats
    n_main = int((inc_df["final_inclusion"] == "include_main").sum())
    n_sens = int((inc_df["final_inclusion"] == "include_sensitivity_only").sum())
    n_ex_qc = int((inc_df["final_inclusion"] == "exclude_low_quality").sum())
    n_ex_dbl = int((inc_df["final_inclusion"] == "exclude_high_confidence_doublet").sum())
    n_high_doublet_samples = int((sample_df["doublet_rate"] > 0.15).sum())

    report = f"""# Doublet / Ambient / Contamination QC Report

- **run_id**: {run_id}
- **n_samples**: {len(sample_df)}
- **n_cells_total**: {len(inc_df)}
- **n_include_main**: {n_main}
- **n_include_sensitivity_only**: {n_sens}
- **n_exclude_low_quality**: {n_ex_qc}
- **n_exclude_high_confidence_doublet**: {n_ex_dbl}
- **n_high_doublet_samples**: {n_high_doublet_samples}

## Cell Inclusion Summary

| Inclusion | Count | Fraction |
|-----------|-------|----------|
| include_main | {n_main} | {n_main/len(inc_df):.4f} |
| include_sensitivity_only | {n_sens} | {n_sens/len(inc_df):.4f} |
| exclude_low_quality | {n_ex_qc} | {n_ex_qc/len(inc_df):.4f} |
| exclude_high_confidence_doublet | {n_ex_dbl} | {n_ex_dbl/len(inc_df):.4f} |

## Sample Doublet Summary

- Samples with doublet rate > 15%: {n_high_doublet_samples}
"""
    (out_dir / "doublet_ambient_qc_report.md").write_text(report, encoding="utf-8")
    print("Wrote doublet_ambient_qc_report.md")

    # Checkpoint
    checkpoint_dir = out_dir / "checkpoints"
    checkpoint_dir.mkdir(parents=True, exist_ok=True)
    (checkpoint_dir / "step2_4.checkpoint.yaml").write_text(
        f"run_id: {run_id}\ncreated_at: {now_iso()}\nsource_input_ref: {manifest_ref}\n"
        f"rows_written: {len(inc_df)}\nstatus: complete\n",
        encoding="utf-8",
    )
    print("Wrote checkpoint")

    return {
        "run_id": run_id,
        "status": "complete",
        "n_cells": len(inc_df),
        "n_samples": len(sample_df),
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Finalize step2.4 inclusion flags")
    parser.add_argument("--run-id", required=True, help="run_id from Step2.0")
    parser.add_argument("--root-dir", default=str(ROOT), help="project root")
    args = parser.parse_args()
    result = finalize_inclusion_flags(run_id=args.run_id, root_dir=args.root_dir)
    print(result)
