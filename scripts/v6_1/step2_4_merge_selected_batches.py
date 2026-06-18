from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq

from step2_4_finalize_branch import finalize_branch

ROOT = Path(__file__).resolve().parents[2]


DEDUP_KEYS = ["source_h5ad", "cell_barcode", "sample_id"]


def merge_parquet(parts: list[Path], out_path: Path) -> tuple[int, int]:
    writer = None
    rows = 0
    deduped_rows = 0
    for p in parts:
        t = pq.read_table(p)
        rows += t.num_rows
        df = t.to_pandas()
        df = df.drop_duplicates(subset=DEDUP_KEYS, keep="first")
        deduped_rows += len(df)
        t = pa.Table.from_pandas(df, preserve_index=False)
        if writer is None:
            writer = pq.ParquetWriter(out_path, t.schema)
        writer.write_table(t)
    if writer is not None:
        writer.close()
    return rows, deduped_rows


def parse_parts(values: list[str]) -> list[Path]:
    return [Path(v) for v in values]


def main(
    run_id: str,
    output_suffix: str,
    doublet_parts: list[Path],
    ambient_parts: list[Path],
    summary_parts: list[Path],
    root_dir: Path | str = ROOT,
) -> dict:
    root = Path(root_dir)
    out = root / "results" / "v6_1" / "step2" / run_id / "03_qc"
    out.mkdir(parents=True, exist_ok=True)

    d_out = out / f"doublet_scores_by_cell.{output_suffix}.parquet"
    a_out = out / f"ambient_contamination_scores_by_cell.{output_suffix}.parquet"
    s_out = out / f"sample_doublet_ambient_summary.{output_suffix}.csv"

    for p in [d_out, a_out, s_out]:
        if p.exists():
            p.unlink()

    d_rows_raw, d_rows = merge_parquet(doublet_parts, d_out)
    a_rows_raw, a_rows = merge_parquet(ambient_parts, a_out)

    s = pd.concat([pd.read_csv(p) for p in summary_parts], ignore_index=True)
    if "created_at" in s.columns:
        s = s.sort_values("created_at")
    s = s.drop_duplicates(subset=["sample_id", "cohort_id"], keep="last")
    s.to_csv(s_out, index=False)

    baseline = finalize_branch(
        run_id=run_id,
        filter_plan_path=out / "cell_filtering_plan.parquet",
        doublet_path=d_out,
        sample_summary_path=s_out,
        output_suffix=f"baseline_{output_suffix}",
        root_dir=root,
    )
    rescue = finalize_branch(
        run_id=run_id,
        filter_plan_path=out / "cell_filtering_plan.targeted_refilter_v2.parquet",
        doublet_path=d_out,
        sample_summary_path=s_out,
        output_suffix=f"rescue_v2_{output_suffix}",
        root_dir=root,
    )
    return {
        "run_id": run_id,
        "output_suffix": output_suffix,
        "doublet_rows_raw": d_rows_raw,
        "doublet_rows": d_rows,
        "ambient_rows_raw": a_rows_raw,
        "ambient_rows": a_rows,
        "baseline": baseline,
        "rescue": rescue,
    }


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--run-id", required=True)
    ap.add_argument("--output-suffix", required=True)
    ap.add_argument("--root-dir", default=str(ROOT))
    ap.add_argument("--doublet-part", action="append", required=True)
    ap.add_argument("--ambient-part", action="append", required=True)
    ap.add_argument("--summary-part", action="append", required=True)
    args = ap.parse_args()
    print(
        main(
            run_id=args.run_id,
            output_suffix=args.output_suffix,
            doublet_parts=parse_parts(args.doublet_part),
            ambient_parts=parse_parts(args.ambient_part),
            summary_parts=parse_parts(args.summary_part),
            root_dir=args.root_dir,
        )
    )
