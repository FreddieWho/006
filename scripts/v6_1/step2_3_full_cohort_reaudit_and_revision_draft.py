from __future__ import annotations

import argparse
from datetime import datetime, timezone
from pathlib import Path

import duckdb
import matplotlib.pyplot as plt
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]


def now_iso() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat()


def main(run_id: str, root_dir: Path | str = ROOT) -> dict:
    root = Path(root_dir)
    run_root = root / "results" / "v6_1" / "step2" / run_id
    qc_dir = run_root / "03_qc"
    out_dir = run_root / "10_reports" / "step2_3_full_cohort_reaudit"
    fig_dir = out_dir / "figures"
    out_dir.mkdir(parents=True, exist_ok=True)
    fig_dir.mkdir(parents=True, exist_ok=True)

    manifest_ref = str(run_root / "00_manifest" / "step2_run_manifest.patched.yaml")
    if not Path(manifest_ref).exists():
        manifest_ref = str(run_root / "00_manifest" / "step2_run_manifest.yaml")

    cell_qc = qc_dir / "cell_qc_metrics.parquet"
    thr_csv = qc_dir / "per_sample_qc_thresholds.csv"
    inc_pq = qc_dir / "final_cell_inclusion_flags.parquet"
    sample_sum = qc_dir / "sample_qc_summary.csv"

    con = duckdb.connect()
    con.execute("PRAGMA threads=8")

    # Dedup joins on (cell_barcode, sample_id) to align with finalization logic.
    sql = f"""
    WITH cell_d AS (
      SELECT * EXCLUDE (rn) FROM (
        SELECT *,
               row_number() OVER (PARTITION BY cell_barcode, sample_id ORDER BY created_at DESC) AS rn
        FROM read_parquet('{cell_qc.as_posix()}')
      ) t WHERE rn = 1
    ),
    thr AS (
      SELECT sample_id,
             min_genes::DOUBLE AS min_genes,
             min_counts::DOUBLE AS min_counts,
             max_genes::DOUBLE AS max_genes,
             max_counts::DOUBLE AS max_counts,
             max_pct_mito::DOUBLE AS max_pct_mito
      FROM read_csv_auto('{thr_csv.as_posix()}', header=true)
    ),
    inc_d AS (
      SELECT * EXCLUDE (rn) FROM (
        SELECT *,
               row_number() OVER (PARTITION BY cell_barcode, sample_id ORDER BY created_at DESC) AS rn
        FROM read_parquet('{inc_pq.as_posix()}')
      ) t WHERE rn = 1
    ),
    m AS (
      SELECT
        c.cohort_id,
        c.sample_id,
        c.cell_barcode,
        c.total_counts::DOUBLE AS total_counts,
        c.n_genes_by_counts::DOUBLE AS n_genes_by_counts,
        c.pct_mito::DOUBLE AS pct_mito,
        t.min_genes,
        t.min_counts,
        t.max_genes,
        t.max_counts,
        t.max_pct_mito,
        i.final_inclusion
      FROM cell_d c
      JOIN thr t USING(sample_id)
      JOIN inc_d i USING(cell_barcode, sample_id)
    )
    SELECT
      cohort_id,
      count(*) AS n_cells,
      sum(CASE WHEN final_inclusion='exclude_low_quality' THEN 1 ELSE 0 END) AS n_exclude_low_quality,
      sum(CASE WHEN final_inclusion='include_main' THEN 1 ELSE 0 END) AS n_include_main,
      sum(CASE WHEN final_inclusion='include_sensitivity_only' THEN 1 ELSE 0 END) AS n_include_sensitivity_only,
      sum(CASE WHEN n_genes_by_counts < min_genes THEN 1 ELSE 0 END) AS fail_min_genes,
      sum(CASE WHEN total_counts < min_counts THEN 1 ELSE 0 END) AS fail_min_counts,
      sum(CASE WHEN n_genes_by_counts > max_genes THEN 1 ELSE 0 END) AS fail_max_genes,
      sum(CASE WHEN total_counts > max_counts THEN 1 ELSE 0 END) AS fail_max_counts,
      sum(CASE WHEN pct_mito > max_pct_mito THEN 1 ELSE 0 END) AS fail_mito,
      quantile_cont(total_counts, 0.5) AS q50_counts,
      quantile_cont(total_counts, 0.1) AS q10_counts,
      quantile_cont(total_counts, 0.9) AS q90_counts,
      quantile_cont(n_genes_by_counts, 0.5) AS q50_genes,
      quantile_cont(n_genes_by_counts, 0.1) AS q10_genes,
      quantile_cont(n_genes_by_counts, 0.9) AS q90_genes,
      quantile_cont(pct_mito, 0.5) AS q50_mito,
      quantile_cont(pct_mito, 0.9) AS q90_mito
    FROM m
    GROUP BY cohort_id
    ORDER BY n_exclude_low_quality DESC
    """
    cohort = con.execute(sql).df()
    cohort["exclude_low_quality_frac"] = cohort["n_exclude_low_quality"] / cohort["n_cells"]
    for c in ["fail_min_genes", "fail_min_counts", "fail_max_genes", "fail_max_counts", "fail_mito"]:
        cohort[f"{c}_frac"] = cohort[c] / cohort["n_cells"]

    # Excluded-only trigger contribution
    ex_sql = f"""
    WITH cell_d AS (
      SELECT * EXCLUDE (rn) FROM (
        SELECT *, row_number() OVER (PARTITION BY cell_barcode, sample_id ORDER BY created_at DESC) AS rn
        FROM read_parquet('{cell_qc.as_posix()}')
      ) t WHERE rn = 1
    ),
    thr AS (
      SELECT sample_id,
             min_genes::DOUBLE AS min_genes,
             min_counts::DOUBLE AS min_counts,
             max_genes::DOUBLE AS max_genes,
             max_counts::DOUBLE AS max_counts,
             max_pct_mito::DOUBLE AS max_pct_mito
      FROM read_csv_auto('{thr_csv.as_posix()}', header=true)
    ),
    inc_d AS (
      SELECT * EXCLUDE (rn) FROM (
        SELECT *, row_number() OVER (PARTITION BY cell_barcode, sample_id ORDER BY created_at DESC) AS rn
        FROM read_parquet('{inc_pq.as_posix()}')
      ) t WHERE rn = 1
    ),
    m AS (
      SELECT c.cohort_id, c.sample_id, c.cell_barcode,
             c.total_counts::DOUBLE AS total_counts,
             c.n_genes_by_counts::DOUBLE AS n_genes_by_counts,
             c.pct_mito::DOUBLE AS pct_mito,
             t.min_genes, t.min_counts, t.max_genes, t.max_counts, t.max_pct_mito,
             i.final_inclusion
      FROM cell_d c
      JOIN thr t USING(sample_id)
      JOIN inc_d i USING(cell_barcode, sample_id)
    )
    SELECT
      cohort_id,
      count(*) AS ex_n,
      sum(CASE WHEN n_genes_by_counts < min_genes THEN 1 ELSE 0 END) AS ex_fail_min_genes,
      sum(CASE WHEN total_counts < min_counts THEN 1 ELSE 0 END) AS ex_fail_min_counts,
      sum(CASE WHEN n_genes_by_counts > max_genes THEN 1 ELSE 0 END) AS ex_fail_max_genes,
      sum(CASE WHEN total_counts > max_counts THEN 1 ELSE 0 END) AS ex_fail_max_counts,
      sum(CASE WHEN pct_mito > max_pct_mito THEN 1 ELSE 0 END) AS ex_fail_mito
    FROM m
    WHERE final_inclusion='exclude_low_quality'
    GROUP BY cohort_id
    """
    ex = con.execute(ex_sql).df()
    for c in ["ex_fail_min_genes", "ex_fail_min_counts", "ex_fail_max_genes", "ex_fail_max_counts", "ex_fail_mito"]:
        ex[f"{c}_share_in_excluded"] = ex[c] / ex["ex_n"]

    # Sample-level distribution summary
    ss = pd.read_csv(sample_sum)
    dist = (
        ss.groupby("cohort_id", as_index=False)
        .agg(
            n_samples=("sample_id", "nunique"),
            n_cells_raw_sum=("n_cells_raw", "sum"),
            median_counts_p50=("median_counts", "median"),
            median_counts_p10=("median_counts", lambda x: x.quantile(0.1)),
            median_counts_p90=("median_counts", lambda x: x.quantile(0.9)),
            median_genes_p50=("median_genes", "median"),
            median_genes_p10=("median_genes", lambda x: x.quantile(0.1)),
            median_genes_p90=("median_genes", lambda x: x.quantile(0.9)),
            median_pct_mito_p50=("median_pct_mito", "median"),
            median_pct_mito_p90=("median_pct_mito", lambda x: x.quantile(0.9)),
        )
    )

    out = cohort.merge(ex, on="cohort_id", how="left").merge(dist, on="cohort_id", how="left")
    out.insert(0, "run_id", run_id)
    out.insert(1, "created_at", now_iso())
    out.insert(2, "input_manifest_ref", manifest_ref)

    # Revision draft
    draft = out[["run_id", "created_at", "input_manifest_ref", "cohort_id", "n_cells", "n_exclude_low_quality", "exclude_low_quality_frac"]].copy()
    draft["dominant_trigger"] = (
        out[
            [
                "ex_fail_min_counts_share_in_excluded",
                "ex_fail_min_genes_share_in_excluded",
                "ex_fail_mito_share_in_excluded",
                "ex_fail_max_counts_share_in_excluded",
                "ex_fail_max_genes_share_in_excluded",
            ]
        ]
        .idxmax(axis=1)
        .fillna("unknown")
        .str.replace("_share_in_excluded", "", regex=False)
    )
    draft["suggested_revision_class"] = "keep_main_thresholds"
    draft.loc[(draft["exclude_low_quality_frac"] > 0.70) & (draft["dominant_trigger"].str.contains("min_counts|min_genes", regex=True)), "suggested_revision_class"] = "mapping_or_counts_layer_audit_first"
    draft.loc[(draft["exclude_low_quality_frac"] > 0.60) & (draft["dominant_trigger"].str.contains("mito", regex=True)), "suggested_revision_class"] = "mito_cap_sensitivity_review"
    draft.loc[(draft["exclude_low_quality_frac"] > 0.25) & (draft["dominant_trigger"].str.contains("max_counts|max_genes", regex=True)), "suggested_revision_class"] = "upper_tail_doublet_or_outlier_review"
    draft["suggested_action"] = "no_change_main"
    draft.loc[draft["suggested_revision_class"] == "mapping_or_counts_layer_audit_first", "suggested_action"] = "verify sample mapping and raw counts source before threshold relaxation"
    draft.loc[draft["suggested_revision_class"] == "mito_cap_sensitivity_review", "suggested_action"] = "add sensitivity branch with max_pct_mito up to 30 for this cohort only"
    draft.loc[draft["suggested_revision_class"] == "upper_tail_doublet_or_outlier_review", "suggested_action"] = "tighten doublet adjudication for high-count tails before relaxing QC"

    # Write tables
    out.to_csv(out_dir / "cohort_qc_reaudit_summary.csv", index=False)
    out[
        [
            "run_id",
            "created_at",
            "input_manifest_ref",
            "cohort_id",
            "ex_n",
            "ex_fail_min_counts_share_in_excluded",
            "ex_fail_min_genes_share_in_excluded",
            "ex_fail_mito_share_in_excluded",
            "ex_fail_max_counts_share_in_excluded",
            "ex_fail_max_genes_share_in_excluded",
        ]
    ].to_csv(out_dir / "cohort_trigger_contribution.csv", index=False)
    out[
        [
            "run_id",
            "created_at",
            "input_manifest_ref",
            "cohort_id",
            "n_samples",
            "n_cells_raw_sum",
            "median_counts_p10",
            "median_counts_p50",
            "median_counts_p90",
            "median_genes_p10",
            "median_genes_p50",
            "median_genes_p90",
            "median_pct_mito_p50",
            "median_pct_mito_p90",
        ]
    ].to_csv(out_dir / "cohort_metric_distribution_summary.csv", index=False)
    draft.to_csv(out_dir / "cohort_revision_draft.csv", index=False)

    # Visuals
    vis = out.sort_values("exclude_low_quality_frac", ascending=False).copy()
    plt.figure(figsize=(14, 6))
    plt.bar(vis["cohort_id"], vis["exclude_low_quality_frac"])
    plt.xticks(rotation=90, fontsize=7)
    plt.ylabel("exclude_low_quality_frac")
    plt.title("Exclude Low-Quality Fraction by Cohort")
    plt.tight_layout()
    plt.savefig(fig_dir / "exclude_low_quality_frac_by_cohort.png", dpi=180)
    plt.close()

    top = vis.head(20).copy()
    comp_cols = [
        "ex_fail_min_counts_share_in_excluded",
        "ex_fail_min_genes_share_in_excluded",
        "ex_fail_mito_share_in_excluded",
        "ex_fail_max_counts_share_in_excluded",
        "ex_fail_max_genes_share_in_excluded",
    ]
    labels = ["min_counts", "min_genes", "mito", "max_counts", "max_genes"]
    plt.figure(figsize=(14, 6))
    bottom = pd.Series([0.0] * len(top), index=top.index)
    for c, lab in zip(comp_cols, labels):
        vals = top[c].fillna(0.0).reindex(top.index)
        plt.bar(top["cohort_id"], vals.values, bottom=bottom.values, label=lab)
        bottom = bottom + vals
    plt.xticks(rotation=90, fontsize=7)
    plt.ylabel("share_in_excluded_cells")
    plt.title("Trigger Contribution in Excluded Cells (Top 20 Cohorts)")
    plt.legend(loc="upper right", ncol=3, fontsize=8)
    plt.tight_layout()
    plt.savefig(fig_dir / "trigger_contribution_stacked_top20.png", dpi=180)
    plt.close()

    # Heatmap-like matrix (simple)
    hm = vis[["cohort_id", "exclude_low_quality_frac", "q50_counts", "q50_genes", "q50_mito"]].copy()
    hm = hm.set_index("cohort_id")
    hm_norm = (hm - hm.min()) / (hm.max() - hm.min() + 1e-9)
    plt.figure(figsize=(12, 10))
    plt.imshow(hm_norm.values, aspect="auto")
    plt.yticks(range(len(hm_norm.index)), hm_norm.index, fontsize=6)
    plt.xticks(range(len(hm_norm.columns)), hm_norm.columns, rotation=45, ha="right")
    plt.title("Cohort Metric Distribution (Normalized)")
    plt.colorbar(label="normalized value")
    plt.tight_layout()
    plt.savefig(fig_dir / "metric_distribution_heatmap.png", dpi=180)
    plt.close()

    md = [
        "# Step2.3 Full-Cohort Reaudit And Revision Draft",
        "",
        f"- run_id: `{run_id}`",
        f"- created_at: `{now_iso()}`",
        f"- input_manifest_ref: `{manifest_ref}`",
        "",
        "## Summary",
        "",
        f"- total cohorts reviewed: {len(out)}",
        f"- highest exclude_low_quality cohort: {vis.iloc[0]['cohort_id']} ({vis.iloc[0]['exclude_low_quality_frac']:.4f})",
        "",
        "## Outputs",
        "",
        "- cohort_qc_reaudit_summary.csv",
        "- cohort_trigger_contribution.csv",
        "- cohort_metric_distribution_summary.csv",
        "- cohort_revision_draft.csv",
        "- figures/exclude_low_quality_frac_by_cohort.png",
        "- figures/trigger_contribution_stacked_top20.png",
        "- figures/metric_distribution_heatmap.png",
        "",
        "## Notes",
        "",
        "- This is a revision draft; no main-threshold overwrite is executed here.",
        "- No response-derived columns were used.",
    ]
    (out_dir / "revision_rationale.md").write_text("\n".join(md), encoding="utf-8")

    return {"run_id": run_id, "cohorts": int(len(out)), "out_dir": str(out_dir)}


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description="Step2.3 full-cohort reaudit and revision draft")
    ap.add_argument("--run-id", required=True)
    ap.add_argument("--root-dir", default=str(ROOT))
    args = ap.parse_args()
    res = main(run_id=args.run_id, root_dir=args.root_dir)
    print(res)
