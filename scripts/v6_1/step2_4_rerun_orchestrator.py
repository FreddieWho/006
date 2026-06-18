from __future__ import annotations

import argparse
from pathlib import Path
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq

from step2_4_doublet_ambient_contamination_audit import run_step2_4
from step2_4_finalize_branch import finalize_branch


def merge_parquet(parts: list[Path], out_path: Path) -> int:
    writer = None
    rows = 0
    for p in parts:
        t = pq.read_table(p)
        rows += t.num_rows
        if writer is None:
            writer = pq.ParquetWriter(out_path, t.schema)
        writer.write_table(t)
    if writer is not None:
        writer.close()
    return rows


def main(run_id: str, root_dir: str, suffix: str, n_shards: int) -> dict:
    root = Path(root_dir)
    out = root / "results" / "v6_1" / "step2" / run_id / "03_qc"
    bdir = out / "doublet_batches_rerun"
    bdir.mkdir(parents=True, exist_ok=True)

    # clean old shard files for this suffix
    for p in bdir.glob(f"*{suffix}.shard*.parquet"):
        p.unlink()
    for p in bdir.glob(f"*{suffix}.shard*.csv"):
        p.unlink()

    results = []
    for i in range(n_shards):
        res = run_step2_4(
            root_dir=root,
            run_id=run_id,
            batch_tag=f"shard{i}",
            sample_shard_index=i,
            sample_shard_count=n_shards,
            output_suffix=suffix,
            batch_dir_name="doublet_batches_rerun",
        )
        results.append(res)

    d_parts = sorted(bdir.glob(f"doublet_scores_by_cell.{suffix}.shard*.parquet"))
    a_parts = sorted(bdir.glob(f"ambient_contamination_scores_by_cell.{suffix}.shard*.parquet"))
    s_parts = sorted(bdir.glob(f"sample_doublet_ambient_summary.{suffix}.shard*.csv"))
    if not d_parts or not a_parts:
        raise RuntimeError("No shard outputs produced")

    d_out = out / f"doublet_scores_by_cell.{suffix}.parquet"
    a_out = out / f"ambient_contamination_scores_by_cell.{suffix}.parquet"
    if d_out.exists():
        d_out.unlink()
    if a_out.exists():
        a_out.unlink()
    d_rows = merge_parquet(d_parts, d_out)
    a_rows = merge_parquet(a_parts, a_out)

    s_out = out / f"sample_doublet_ambient_summary.{suffix}.csv"
    if s_parts:
        s = pd.concat([pd.read_csv(p) for p in s_parts], ignore_index=True)
        if "created_at" in s.columns:
            s = s.sort_values("created_at")
        s = s.drop_duplicates(subset=["sample_id", "cohort_id"], keep="last")
        s.to_csv(s_out, index=False)
    else:
        s_out.write_text("run_id,created_at,input_manifest_ref,sample_id,cohort_id,n_cells_input,doublet_rate,high_confidence_doublet_rate,ambient_audit_status,rbc_warning,stress_warning,tumor_contamination_warning,sample_qc_decision,doublet_method,doublet_run_status\n")

    # finalize both branches with same doublet results
    baseline = finalize_branch(
        run_id=run_id,
        filter_plan_path=out / "cell_filtering_plan.parquet",
        doublet_path=d_out,
        sample_summary_path=s_out,
        output_suffix=f"baseline_{suffix}",
        root_dir=root,
    )
    rescue = finalize_branch(
        run_id=run_id,
        filter_plan_path=out / "cell_filtering_plan.targeted_refilter_v2.parquet",
        doublet_path=d_out,
        sample_summary_path=s_out,
        output_suffix=f"rescue_v2_{suffix}",
        root_dir=root,
    )
    return {
        "run_id": run_id,
        "suffix": suffix,
        "n_shards": n_shards,
        "doublet_rows": d_rows,
        "ambient_rows": a_rows,
        "baseline": baseline,
        "rescue": rescue,
    }


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--run-id", required=True)
    ap.add_argument("--root-dir", default="/home/huyudi/006")
    ap.add_argument("--suffix", default="rerun_0509c")
    ap.add_argument("--n-shards", type=int, default=6)
    args = ap.parse_args()
    print(main(args.run_id, args.root_dir, args.suffix, args.n_shards))

