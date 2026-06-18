from __future__ import annotations

import argparse
from datetime import datetime, timezone
from pathlib import Path

import duckdb
import pandas as pd
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

    # baseline cohorts by exclude frac
    inc = pd.read_parquet(qdir / "final_cell_inclusion_flags.parquet", columns=["cohort_id", "sample_id", "cell_barcode", "final_inclusion", "doublet_call"])
    base = (
        inc.groupby("cohort_id", as_index=False)
        .agg(n_cells=("cell_barcode", "size"), n_exclude=("final_inclusion", lambda s: (s == "exclude_low_quality").sum()))
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

    # thresholds
    thr = pd.read_csv(qdir / "per_sample_qc_thresholds.csv")
    ss = pd.read_csv(qdir / "sample_qc_summary.csv", usecols=["sample_id", "cohort_id"])
    if "cohort_id" in thr.columns:
        thr = thr.merge(ss, on="sample_id", how="left", suffixes=("", "_s"))
        if "cohort_id_s" in thr.columns:
            thr["cohort_id"] = thr["cohort_id"].fillna(thr["cohort_id_s"])
            thr = thr.drop(columns=["cohort_id_s"])
    else:
        thr = thr.merge(ss, on="sample_id", how="left")
    thr["cohort_id"] = thr["cohort_id"].astype(str)

    thr["min_counts_adj"] = thr["min_counts"]
    thr["min_genes_adj"] = thr["min_genes"]
    thr["max_counts_adj"] = thr["max_counts"]
    thr["max_genes_adj"] = thr["max_genes"]
    thr["max_pct_mito_adj"] = thr["max_pct_mito"]
    thr["adjustment_phase"] = "none"

    # phase1 mild for all exclude_frac>0.2
    m = thr["cohort_id"].isin(cohorts_gt_02)
    thr.loc[m, "min_counts_adj"] = (thr.loc[m, "min_counts"] * 0.8).clip(lower=250)
    thr.loc[m, "min_genes_adj"] = (thr.loc[m, "min_genes"] * 0.8).clip(lower=120)
    thr.loc[m, "max_counts_adj"] = thr.loc[m, "max_counts"] * 1.2
    thr.loc[m, "max_genes_adj"] = thr.loc[m, "max_genes"] * 1.2
    thr.loc[m, "max_pct_mito_adj"] = (thr.loc[m, "max_pct_mito"] + 5).clip(upper=35)
    thr.loc[m, "adjustment_phase"] = "phase1"

    # phase1 strong for top3
    top3 = set(top3_gt_06)
    m3 = thr["cohort_id"].isin(top3)
    # GSE205506 mito-driven: mostly mito cap; others low-count-driven
    m3_mito = m3 & (thr["cohort_id"] == "GSE205506")
    m3_low = m3 & (~m3_mito)
    thr.loc[m3_low, "min_counts_adj"] = 1
    thr.loc[m3_low, "min_genes_adj"] = 1
    thr.loc[m3_low, "max_counts_adj"] = thr.loc[m3_low, "max_counts"] * 1.5
    thr.loc[m3_low, "max_genes_adj"] = thr.loc[m3_low, "max_genes"] * 1.5
    thr.loc[m3_low, "max_pct_mito_adj"] = (thr.loc[m3_low, "max_pct_mito"] + 10).clip(upper=45)
    thr.loc[m3_mito, "min_counts_adj"] = (thr.loc[m3_mito, "min_counts"] * 0.7).clip(lower=200)
    thr.loc[m3_mito, "min_genes_adj"] = (thr.loc[m3_mito, "min_genes"] * 0.7).clip(lower=100)
    thr.loc[m3_mito, "max_pct_mito_adj"] = (thr.loc[m3_mito, "max_pct_mito"] + 15).clip(upper=45)
    thr.loc[m3, "adjustment_phase"] = "phase1_top3"

    con = duckdb.connect()
    con.execute("PRAGMA threads=8")
    con.register("thr_adj", thr[["sample_id", "cohort_id", "min_counts_adj", "min_genes_adj", "max_counts_adj", "max_genes_adj", "max_pct_mito_adj", "adjustment_phase"]])
    con.execute(f"""
    CREATE OR REPLACE TABLE cell_d AS
    SELECT * EXCLUDE (rn) FROM (
      SELECT *, row_number() OVER (PARTITION BY cell_barcode, sample_id ORDER BY created_at DESC) AS rn
      FROM read_parquet('{(qdir / "cell_qc_metrics.parquet").as_posix()}')
    ) t WHERE rn=1
    """)
    con.execute(f"""
    CREATE OR REPLACE TABLE base_inc AS
    SELECT * EXCLUDE (rn) FROM (
      SELECT *, row_number() OVER (PARTITION BY cell_barcode, sample_id ORDER BY created_at DESC) AS rn
      FROM read_parquet('{(qdir / "final_cell_inclusion_flags.parquet").as_posix()}')
    ) t WHERE rn=1
    """)
    con.execute("""
    CREATE OR REPLACE TABLE pass1 AS
    SELECT
      c.run_id, c.source_h5ad, c.cohort_id, c.cell_barcode, c.sample_id,
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
      END AS filter_plan,
      t.adjustment_phase
    FROM cell_d c
    JOIN thr_adj t USING(sample_id)
    """)
    p1 = con.execute("""
    SELECT cohort_id, count(*) n_cells, sum(CASE WHEN filter_plan='remove_low_quality' THEN 1 ELSE 0 END) n_ex
    FROM pass1 GROUP BY cohort_id
    """).df()
    p1["exclude_frac_pass1"] = p1["n_ex"] / p1["n_cells"]
    need_phase2 = set(p1[p1["cohort_id"].astype(str).isin(top3) & (p1["exclude_frac_pass1"] > 0.4)]["cohort_id"].astype(str))

    if need_phase2:
        m2 = thr["cohort_id"].isin(need_phase2)
        # stronger salvage branch
        thr.loc[m2 & (thr["cohort_id"] != "GSE205506"), "min_counts_adj"] = 0
        thr.loc[m2 & (thr["cohort_id"] != "GSE205506"), "min_genes_adj"] = 0
        thr.loc[m2 & (thr["cohort_id"] != "GSE205506"), "max_pct_mito_adj"] = 50
        thr.loc[m2 & (thr["cohort_id"] != "GSE205506"), "max_counts_adj"] = thr.loc[m2 & (thr["cohort_id"] != "GSE205506"), "max_counts"] * 2.0
        thr.loc[m2 & (thr["cohort_id"] != "GSE205506"), "max_genes_adj"] = thr.loc[m2 & (thr["cohort_id"] != "GSE205506"), "max_genes"] * 2.0
        thr.loc[m2 & (thr["cohort_id"] == "GSE205506"), "max_pct_mito_adj"] = 55
        thr.loc[m2, "adjustment_phase"] = "phase2_top3"

        con.unregister("thr_adj")
        con.register("thr_adj", thr[["sample_id", "cohort_id", "min_counts_adj", "min_genes_adj", "max_counts_adj", "max_genes_adj", "max_pct_mito_adj", "adjustment_phase"]])
        con.execute("""
        CREATE OR REPLACE TABLE pass1 AS
        SELECT
          c.run_id, c.source_h5ad, c.cohort_id, c.cell_barcode, c.sample_id,
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
          END AS filter_plan,
          t.adjustment_phase
        FROM cell_d c
        JOIN thr_adj t USING(sample_id)
        """)

    # save targeted filter plan
    plan_tbl = con.execute("""
    SELECT run_id, ? AS created_at, ? AS input_manifest_ref,
           source_h5ad, cohort_id, cell_barcode, sample_id, filter_plan, adjustment_phase
    FROM pass1
    """, [now_iso(), manifest_ref]).fetch_arrow_table()
    pq.write_table(plan_tbl, qdir / "cell_filtering_plan.targeted_refilter_v2.parquet")

    # combine with doublet call from baseline inclusion
    con.execute("""
    CREATE OR REPLACE TABLE final_v2 AS
    SELECT
      p.run_id,
      ? AS created_at,
      ? AS input_manifest_ref,
      p.source_h5ad, p.cohort_id, p.cell_barcode, p.sample_id,
      CASE
        WHEN b.doublet_call='high_confidence_doublet' THEN 'exclude_high_confidence_doublet'
        WHEN p.filter_plan='remove_low_quality' THEN 'exclude_low_quality'
        WHEN p.filter_plan='review_borderline' THEN 'include_sensitivity_only'
        WHEN b.doublet_call='borderline_doublet' THEN 'include_sensitivity_only'
        ELSE 'include_main'
      END AS final_inclusion,
      COALESCE(b.doublet_call,'unknown') AS doublet_call,
      p.filter_plan
    FROM pass1 p
    LEFT JOIN base_inc b USING(cell_barcode, sample_id)
    """, [now_iso(), manifest_ref])
    final_tbl = con.execute("SELECT * FROM final_v2").fetch_arrow_table()
    pq.write_table(final_tbl, qdir / "final_cell_inclusion_flags.targeted_refilter_v2.parquet")

    # rescue stats
    bf = con.execute("""
    SELECT cohort_id,
           count(*) AS n_cells,
           sum(CASE WHEN final_inclusion='exclude_low_quality' THEN 1 ELSE 0 END) AS n_ex_before
    FROM base_inc GROUP BY cohort_id
    """).df()
    af = con.execute("""
    SELECT cohort_id,
           count(*) AS n_cells_after,
           sum(CASE WHEN final_inclusion='exclude_low_quality' THEN 1 ELSE 0 END) AS n_ex_after
    FROM final_v2 GROUP BY cohort_id
    """).df()
    rep = bf.merge(af, on="cohort_id", how="inner")
    rep["exclude_frac_before"] = rep["n_ex_before"] / rep["n_cells"]
    rep["exclude_frac_after"] = rep["n_ex_after"] / rep["n_cells_after"]
    rep["rescued_cells"] = rep["n_ex_before"] - rep["n_ex_after"]
    rep["rescued_frac_of_total"] = rep["rescued_cells"] / rep["n_cells"]
    rep.insert(0, "run_id", run_id)
    rep.insert(1, "created_at", now_iso())
    rep.insert(2, "input_manifest_ref", manifest_ref)
    rep = rep.sort_values("rescued_cells", ascending=False)
    rep.to_csv(rdir / "step2_3_targeted_refilter_v2_rescue_by_cohort.csv", index=False)

    total_rescued = int(rep["rescued_cells"].sum())
    top3_after = rep[rep["cohort_id"].isin(top3_gt_06)][["cohort_id", "exclude_frac_after"]]
    top3_after.to_csv(rdir / "step2_3_targeted_refilter_v2_top3_after.csv", index=False)

    md = [
        "# Step2.3 Targeted Refilter V2 (Top3 Intensive Reaudit)",
        "",
        f"- run_id: `{run_id}`",
        f"- created_at: `{now_iso()}`",
        f"- top3 initial (>0.6): {', '.join(top3_gt_06) if top3_gt_06 else 'none'}",
        f"- phase2 applied cohorts: {', '.join(sorted(need_phase2)) if need_phase2 else 'none'}",
        "",
        f"- total rescued cells (exclude_low_quality -> non-exclude_low_quality): **{total_rescued}**",
        "",
        "## Outputs",
        "- 03_qc/cell_filtering_plan.targeted_refilter_v2.parquet",
        "- 03_qc/final_cell_inclusion_flags.targeted_refilter_v2.parquet",
        "- 10_reports/step2_3_targeted_refilter_v2_rescue_by_cohort.csv",
        "- 10_reports/step2_3_targeted_refilter_v2_top3_after.csv",
        "",
        "## Note",
        "- This is a sensitivity rescue branch and does not overwrite baseline files.",
    ]
    (rdir / "step2_3_targeted_refilter_v2_report.md").write_text("\n".join(md), encoding="utf-8")

    return {
        "run_id": run_id,
        "top3": top3_gt_06,
        "phase2_applied": sorted(list(need_phase2)),
        "total_rescued_cells": total_rescued,
    }


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--run-id", required=True)
    ap.add_argument("--root-dir", default=str(ROOT))
    args = ap.parse_args()
    print(main(run_id=args.run_id, root_dir=args.root_dir))

