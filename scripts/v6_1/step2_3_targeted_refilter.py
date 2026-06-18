from __future__ import annotations

import argparse
from datetime import datetime, timezone
from pathlib import Path

import duckdb
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq

ROOT = Path(__file__).resolve().parents[2]


def now_iso() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat()


def main(run_id: str, root_dir: Path | str = ROOT) -> dict:
    root = Path(root_dir)
    run_root = root / "results" / "v6_1" / "step2" / run_id
    qdir = run_root / "03_qc"
    rdir = run_root / "10_reports"
    rdir.mkdir(parents=True, exist_ok=True)

    manifest_ref = str(run_root / "00_manifest" / "step2_run_manifest.patched.yaml")
    if not Path(manifest_ref).exists():
        manifest_ref = str(run_root / "00_manifest" / "step2_run_manifest.yaml")

    inc = pd.read_parquet(qdir / "final_cell_inclusion_flags.parquet", columns=["cohort_id", "final_inclusion"])
    base = (
        inc.groupby("cohort_id", as_index=False)
        .agg(n_cells=("final_inclusion", "size"), n_exclude=("final_inclusion", lambda s: (s == "exclude_low_quality").sum()))
    )
    base["exclude_frac_before"] = base["n_exclude"] / base["n_cells"]

    cohorts_gt_02 = set(base.loc[base["exclude_frac_before"] > 0.2, "cohort_id"].astype(str))
    top3_gt_06 = (
        base.loc[base["exclude_frac_before"] > 0.6]
        .sort_values("exclude_frac_before", ascending=False)
        .head(3)["cohort_id"]
        .astype(str)
        .tolist()
    )

    thr = pd.read_csv(qdir / "per_sample_qc_thresholds.csv")
    ss = pd.read_csv(qdir / "sample_qc_summary.csv", usecols=["sample_id", "cohort_id"])
    if "cohort_id" not in thr.columns:
        thr = thr.merge(ss, on="sample_id", how="left")
    else:
        thr = thr.merge(ss, on="sample_id", how="left", suffixes=("", "_from_summary"))
        if "cohort_id_from_summary" in thr.columns:
            thr["cohort_id"] = thr["cohort_id"].fillna(thr["cohort_id_from_summary"])
            thr = thr.drop(columns=["cohort_id_from_summary"])
    thr["cohort_id"] = thr["cohort_id"].astype(str)

    # phase-1 adjustment for exclude_frac > 0.2 cohorts
    m1 = thr["cohort_id"].isin(cohorts_gt_02)
    thr["min_counts_adj"] = thr["min_counts"]
    thr["min_genes_adj"] = thr["min_genes"]
    thr["max_counts_adj"] = thr["max_counts"]
    thr["max_genes_adj"] = thr["max_genes"]
    thr["max_pct_mito_adj"] = thr["max_pct_mito"]
    thr.loc[m1, "min_counts_adj"] = (thr.loc[m1, "min_counts"] * 0.7).clip(lower=200)
    thr.loc[m1, "min_genes_adj"] = (thr.loc[m1, "min_genes"] * 0.75).clip(lower=100)
    thr.loc[m1, "max_counts_adj"] = thr.loc[m1, "max_counts"] * 1.2
    thr.loc[m1, "max_genes_adj"] = thr.loc[m1, "max_genes"] * 1.2
    thr.loc[m1, "max_pct_mito_adj"] = (thr.loc[m1, "max_pct_mito"] + 5).clip(upper=35)
    thr["adjustment_phase"] = "none"
    thr.loc[m1, "adjustment_phase"] = "phase1"

    con = duckdb.connect()
    con.execute("PRAGMA threads=8")
    con.register("thr_adj", thr[[
        "sample_id", "cohort_id", "min_counts_adj", "min_genes_adj", "max_counts_adj", "max_genes_adj", "max_pct_mito_adj", "adjustment_phase"
    ]])

    con.execute(f"""
    CREATE OR REPLACE TABLE cell_d AS
    SELECT * EXCLUDE (rn) FROM (
      SELECT *, row_number() OVER (PARTITION BY cell_barcode, sample_id ORDER BY created_at DESC) AS rn
      FROM read_parquet('{(qdir / "cell_qc_metrics.parquet").as_posix()}')
    ) t WHERE rn=1
    """)

    con.execute("""
    CREATE OR REPLACE TABLE pass1 AS
    SELECT
      c.run_id,
      c.source_h5ad,
      c.cohort_id,
      c.cell_barcode,
      c.sample_id,
      CASE
        WHEN c.n_genes_by_counts < t.min_genes_adj THEN 'remove_low_quality'
        WHEN c.total_counts < t.min_counts_adj THEN 'remove_low_quality'
        WHEN c.n_genes_by_counts > t.max_genes_adj THEN 'remove_low_quality'
        WHEN c.total_counts > t.max_counts_adj THEN 'remove_low_quality'
        WHEN c.pct_mito > t.max_pct_mito_adj THEN 'remove_low_quality'
        WHEN (
          (c.n_genes_by_counts < t.min_genes_adj * 1.1 AND c.n_genes_by_counts >= t.min_genes_adj) OR
          (c.total_counts < t.min_counts_adj * 1.1 AND c.total_counts >= t.min_counts_adj) OR
          (c.n_genes_by_counts > t.max_genes_adj * 0.9 AND c.n_genes_by_counts <= t.max_genes_adj) OR
          (c.total_counts > t.max_counts_adj * 0.9 AND c.total_counts <= t.max_counts_adj) OR
          (c.pct_mito > t.max_pct_mito_adj * 0.9 AND c.pct_mito <= t.max_pct_mito_adj)
        ) THEN 'review_borderline'
        ELSE 'keep_high_confidence'
      END AS filter_plan_targeted,
      t.adjustment_phase
    FROM cell_d c
    JOIN thr_adj t USING(sample_id)
    """)

    p1 = con.execute("""
    SELECT cohort_id,
           count(*) AS n_cells,
           sum(CASE WHEN filter_plan_targeted='remove_low_quality' THEN 1 ELSE 0 END) AS n_exclude
    FROM pass1 GROUP BY cohort_id
    """).df()
    p1["exclude_frac_pass1"] = p1["n_exclude"] / p1["n_cells"]

    need_phase2 = set(
        p1[p1["cohort_id"].astype(str).isin(top3_gt_06) & (p1["exclude_frac_pass1"] > 0.4)]["cohort_id"].astype(str)
    )

    # phase-2 stronger relaxation on top3 still > 0.4
    if need_phase2:
        m2 = thr["cohort_id"].isin(need_phase2)
        thr.loc[m2, "min_counts_adj"] = (thr.loc[m2, "min_counts"] * 0.5).clip(lower=100)
        thr.loc[m2, "min_genes_adj"] = (thr.loc[m2, "min_genes"] * 0.5).clip(lower=80)
        thr.loc[m2, "max_counts_adj"] = thr.loc[m2, "max_counts"] * 1.4
        thr.loc[m2, "max_genes_adj"] = thr.loc[m2, "max_genes"] * 1.4
        thr.loc[m2, "max_pct_mito_adj"] = (thr.loc[m2, "max_pct_mito"] + 10).clip(upper=40)
        thr.loc[m2, "adjustment_phase"] = "phase2"

        con.unregister("thr_adj")
        con.register("thr_adj", thr[[
            "sample_id", "cohort_id", "min_counts_adj", "min_genes_adj", "max_counts_adj", "max_genes_adj", "max_pct_mito_adj", "adjustment_phase"
        ]])

        con.execute("""
        CREATE OR REPLACE TABLE pass1 AS
        SELECT
          c.run_id,
          c.source_h5ad,
          c.cohort_id,
          c.cell_barcode,
          c.sample_id,
          CASE
            WHEN c.n_genes_by_counts < t.min_genes_adj THEN 'remove_low_quality'
            WHEN c.total_counts < t.min_counts_adj THEN 'remove_low_quality'
            WHEN c.n_genes_by_counts > t.max_genes_adj THEN 'remove_low_quality'
            WHEN c.total_counts > t.max_counts_adj THEN 'remove_low_quality'
            WHEN c.pct_mito > t.max_pct_mito_adj THEN 'remove_low_quality'
            WHEN (
              (c.n_genes_by_counts < t.min_genes_adj * 1.1 AND c.n_genes_by_counts >= t.min_genes_adj) OR
              (c.total_counts < t.min_counts_adj * 1.1 AND c.total_counts >= t.min_counts_adj) OR
              (c.n_genes_by_counts > t.max_genes_adj * 0.9 AND c.n_genes_by_counts <= t.max_genes_adj) OR
              (c.total_counts > t.max_counts_adj * 0.9 AND c.total_counts <= t.max_counts_adj) OR
              (c.pct_mito > t.max_pct_mito_adj * 0.9 AND c.pct_mito <= t.max_pct_mito_adj)
            ) THEN 'review_borderline'
            ELSE 'keep_high_confidence'
          END AS filter_plan_targeted,
          t.adjustment_phase
        FROM cell_d c
        JOIN thr_adj t USING(sample_id)
        """)

    # write targeted filtering plan
    targeted_plan = con.execute("""
    SELECT
      run_id,
      ? AS created_at,
      ? AS input_manifest_ref,
      source_h5ad,
      cohort_id,
      cell_barcode,
      sample_id,
      filter_plan_targeted AS filter_plan,
      adjustment_phase
    FROM pass1
    """, [now_iso(), manifest_ref]).fetch_arrow_table()
    pq.write_table(targeted_plan, qdir / "cell_filtering_plan.targeted_refilter.parquet")

    # build report
    after = con.execute("""
    SELECT cohort_id,
           count(*) AS n_cells,
           sum(CASE WHEN filter_plan_targeted='remove_low_quality' THEN 1 ELSE 0 END) AS n_exclude_after
    FROM pass1 GROUP BY cohort_id
    """).df()
    after["exclude_frac_after"] = after["n_exclude_after"] / after["n_cells"]

    rep = base[["cohort_id", "n_cells", "exclude_frac_before"]].merge(after[["cohort_id", "exclude_frac_after"]], on="cohort_id", how="left")
    rep["delta"] = rep["exclude_frac_after"] - rep["exclude_frac_before"]
    rep = rep.sort_values("exclude_frac_before", ascending=False)
    rep.insert(0, "run_id", run_id)
    rep.insert(1, "created_at", now_iso())
    rep.insert(2, "input_manifest_ref", manifest_ref)
    rep.to_csv(rdir / "step2_3_targeted_refilter_cohort_before_after.csv", index=False)

    thr_out = thr[[
        "sample_id",
        "cohort_id",
        "min_genes",
        "min_counts",
        "max_genes",
        "max_counts",
        "max_pct_mito",
        "min_genes_adj",
        "min_counts_adj",
        "max_genes_adj",
        "max_counts_adj",
        "max_pct_mito_adj",
        "adjustment_phase",
    ]].copy()
    thr_out.insert(0, "run_id", run_id)
    thr_out.insert(1, "created_at", now_iso())
    thr_out.insert(2, "input_manifest_ref", manifest_ref)
    thr_out.to_csv(rdir / "step2_3_targeted_refilter_thresholds_by_sample.csv", index=False)

    md = [
        "# Step2.3 Targeted Refilter Revision Draft",
        "",
        f"- run_id: `{run_id}`",
        f"- created_at: `{now_iso()}`",
        "",
        f"- cohorts with exclude_frac_before > 0.2: {len(cohorts_gt_02)}",
        f"- top3 cohorts with exclude_frac_before > 0.6: {', '.join(top3_gt_06) if top3_gt_06 else 'none'}",
        f"- phase2 applied cohorts (still >0.4 after pass1 among top3): {', '.join(sorted(need_phase2)) if need_phase2 else 'none'}",
        "",
        "## Outputs",
        "",
        "- 03_qc/cell_filtering_plan.targeted_refilter.parquet",
        "- 10_reports/step2_3_targeted_refilter_cohort_before_after.csv",
        "- 10_reports/step2_3_targeted_refilter_thresholds_by_sample.csv",
        "",
        "## Notes",
        "",
        "- This is a targeted sensitivity refilter plan and does not overwrite the baseline filtering files.",
        "- No response-derived fields were used.",
    ]
    (rdir / "step2_3_targeted_refilter_revision.md").write_text("\n".join(md), encoding="utf-8")

    return {
        "run_id": run_id,
        "n_cohorts_gt_02": len(cohorts_gt_02),
        "top3_gt_06": top3_gt_06,
        "phase2_applied": sorted(list(need_phase2)),
    }


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description="Step2.3 targeted refilter")
    ap.add_argument("--run-id", required=True)
    ap.add_argument("--root-dir", default=str(ROOT))
    args = ap.parse_args()
    print(main(run_id=args.run_id, root_dir=args.root_dir))
