from __future__ import annotations

from datetime import datetime, timezone
from pathlib import Path

import duckdb
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq

ROOT = Path(__file__).resolve().parents[2]
RUN_ID = "step2_v6_1_0505_0319"


def now_iso() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat()


BRANCHES = {
    "conservative": {
        "min_counts_floor": 250,
        "min_genes_floor": 100,
        "min_counts_multiplier": 0.70,
        "min_genes_multiplier": 0.70,
        "upper_multiplier": 1.50,
        "mito_add": 10.0,
        "mito_cap": 40.0,
        "borderline_margin": 0.10,
        "description": "Moderate rescue: lower hard minima and relax upper-tail/mitochondrial thresholds.",
    },
    "balanced": {
        "min_counts_floor": 150,
        "min_genes_floor": 80,
        "min_counts_multiplier": 0.50,
        "min_genes_multiplier": 0.50,
        "upper_multiplier": 2.00,
        "mito_add": 15.0,
        "mito_cap": 45.0,
        "borderline_margin": 0.10,
        "description": "Recommended branch: rescues low-depth tissue cells while retaining extreme-tail review.",
    },
    "max_rescue_review": {
        "min_counts_floor": 100,
        "min_genes_floor": 50,
        "min_counts_multiplier": 0.35,
        "min_genes_multiplier": 0.35,
        "upper_multiplier": 3.00,
        "mito_add": 25.0,
        "mito_cap": 55.0,
        "borderline_margin": 0.15,
        "description": "Sensitivity branch: broad rescue, intended for manual review and not default main analysis.",
    },
}


def manifest_ref(run_root: Path) -> str:
    patched = run_root / "00_manifest" / "step2_run_manifest.patched.yaml"
    original = run_root / "00_manifest" / "step2_run_manifest.yaml"
    return str(patched if patched.exists() else original)


def reclassify_doublet_calls(qdir: Path, created_at: str, mref: str) -> Path:
    src = qdir / "doublet_scores_by_cell.parquet"
    dst = qdir / "doublet_scores_by_cell.qc0513_reclassified.parquet"
    df = pd.read_parquet(src)
    status = df["doublet_run_status"].astype(str)
    call = df["doublet_call"].astype(str)
    new_call = call.copy()
    new_call[(call == "unknown") & status.str.startswith("failed_with_reason")] = "not_called_failed"
    new_call[(call == "unknown") & status.str.contains("not_possible_too_few_cells", regex=False)] = "not_called_too_few_cells"
    new_call[(call == "unknown") & status.str.contains("implausible_predicted_doublet_rate", regex=False)] = "not_called_implausible_scrublet_rate"
    new_call[(call == "unknown") & (new_call == "unknown")] = "not_called_unresolved"
    df["doublet_call_original"] = call
    df["doublet_call"] = new_call
    df["doublet_reclassification_reason"] = "unchanged"
    df.loc[call != new_call, "doublet_reclassification_reason"] = "scrublet_not_called_not_filtered"
    df["created_at"] = created_at
    df["input_manifest_ref"] = mref
    pq.write_table(pa.Table.from_pandas(df, preserve_index=False), dst)
    return dst


def finalize_from_tables(qdir: Path, branch: str, filter_path: Path, doublet_path: Path, created_at: str, mref: str) -> Path:
    con = duckdb.connect()
    con.execute("PRAGMA threads=8")
    con.execute(f"""
    CREATE OR REPLACE TABLE filt AS
    SELECT * FROM read_parquet('{filter_path.as_posix()}')
    """)
    con.execute(f"""
    CREATE OR REPLACE TABLE dbl AS
    SELECT source_h5ad, cell_barcode, sample_id, doublet_call
    FROM read_parquet('{doublet_path.as_posix()}')
    """)
    out = qdir / f"final_cell_inclusion_flags.qc0513_{branch}.parquet"
    con.execute("""
    CREATE OR REPLACE TABLE final AS
    SELECT
      f.run_id,
      ? AS created_at,
      ? AS input_manifest_ref,
      f.source_h5ad,
      f.cohort_id,
      f.cell_barcode,
      f.sample_id,
      CASE
        WHEN d.doublet_call='high_confidence_doublet' THEN 'exclude_high_confidence_doublet'
        WHEN f.filter_plan='remove_low_quality' THEN 'exclude_low_quality'
        WHEN f.filter_plan='review_borderline' THEN 'include_sensitivity_only'
        WHEN d.doublet_call='borderline_doublet' THEN 'include_sensitivity_only'
        ELSE 'include_main'
      END AS final_inclusion,
      CASE
        WHEN d.doublet_call='high_confidence_doublet' THEN 'high_confidence_doublet'
        WHEN f.filter_plan='remove_low_quality' THEN f.plan_reason
        WHEN f.filter_plan='review_borderline' THEN 'borderline_qc'
        WHEN d.doublet_call='borderline_doublet' THEN 'borderline_doublet'
        WHEN starts_with(COALESCE(d.doublet_call,'not_called_missing'), 'not_called') THEN 'pass_qc_doublet_not_called_degraded'
        ELSE 'pass_all_qc_and_doublet_checks'
      END AS inclusion_reason,
      CASE
        WHEN d.doublet_call='high_confidence_doublet' THEN 'doublet_exclusion'
        WHEN f.filter_plan='remove_low_quality' THEN 'qc_exclusion'
        WHEN f.filter_plan='review_borderline' THEN 'borderline_qc'
        WHEN d.doublet_call='borderline_doublet' THEN 'borderline_doublet'
        WHEN starts_with(COALESCE(d.doublet_call,'not_called_missing'), 'not_called') THEN 'doublet_not_called_degraded'
        ELSE 'none'
      END AS source_qc_flags,
      COALESCE(d.doublet_call, 'not_called_missing') AS doublet_call,
      f.filter_plan
    FROM filt f
    LEFT JOIN dbl d USING(source_h5ad, cell_barcode, sample_id)
    """, [created_at, mref])
    con.execute(f"COPY final TO '{out.as_posix()}' (FORMAT PARQUET)")
    con.close()
    return out


def write_branch_filter(qdir: Path, branch: str, cfg: dict, created_at: str, mref: str) -> Path:
    thresholds = pd.read_csv(qdir / "per_sample_qc_thresholds.csv")
    thresholds["min_counts_adj"] = (thresholds["min_counts"] * cfg["min_counts_multiplier"]).clip(lower=cfg["min_counts_floor"])
    thresholds["min_genes_adj"] = (thresholds["min_genes"] * cfg["min_genes_multiplier"]).clip(lower=cfg["min_genes_floor"])
    thresholds["max_counts_adj"] = thresholds["max_counts"] * cfg["upper_multiplier"]
    thresholds["max_genes_adj"] = thresholds["max_genes"] * cfg["upper_multiplier"]
    thresholds["max_pct_mito_adj"] = (thresholds["max_pct_mito"] + cfg["mito_add"]).clip(upper=cfg["mito_cap"])
    thresholds["branch"] = branch
    threshold_out = qdir / f"per_sample_qc_thresholds.qc0513_{branch}.csv"
    thresholds.to_csv(threshold_out, index=False)

    con = duckdb.connect()
    con.execute("PRAGMA threads=8")
    con.register(
        "thr",
        thresholds[[
            "sample_id",
            "min_counts_adj",
            "min_genes_adj",
            "max_counts_adj",
            "max_genes_adj",
            "max_pct_mito_adj",
        ]],
    )
    qc_path = qdir / "cell_qc_metrics.parquet"
    out = qdir / f"cell_filtering_plan.qc0513_{branch}.parquet"
    m = float(cfg["borderline_margin"])
    con.execute(f"""
    CREATE OR REPLACE TABLE plan AS
    SELECT
      c.run_id,
      ? AS created_at,
      ? AS input_manifest_ref,
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
          (c.n_genes_by_counts < t.min_genes_adj * (1 + {m}) AND c.n_genes_by_counts >= t.min_genes_adj) OR
          (c.total_counts < t.min_counts_adj * (1 + {m}) AND c.total_counts >= t.min_counts_adj) OR
          (c.n_genes_by_counts > t.max_genes_adj * (1 - {m}) AND c.n_genes_by_counts <= t.max_genes_adj) OR
          (c.total_counts > t.max_counts_adj * (1 - {m}) AND c.total_counts <= t.max_counts_adj) OR
          (c.pct_mito > t.max_pct_mito_adj * (1 - {m}) AND c.pct_mito <= t.max_pct_mito_adj)
        ) THEN 'review_borderline'
        ELSE 'keep_high_confidence'
      END AS filter_plan,
      CASE
        WHEN c.n_genes_by_counts < t.min_genes_adj THEN 'low_genes'
        WHEN c.total_counts < t.min_counts_adj THEN 'low_counts'
        WHEN c.n_genes_by_counts > t.max_genes_adj THEN 'high_genes'
        WHEN c.total_counts > t.max_counts_adj THEN 'high_counts'
        WHEN c.pct_mito > t.max_pct_mito_adj THEN 'high_mito'
        WHEN (
          (c.n_genes_by_counts < t.min_genes_adj * (1 + {m}) AND c.n_genes_by_counts >= t.min_genes_adj) OR
          (c.total_counts < t.min_counts_adj * (1 + {m}) AND c.total_counts >= t.min_counts_adj) OR
          (c.n_genes_by_counts > t.max_genes_adj * (1 - {m}) AND c.n_genes_by_counts <= t.max_genes_adj) OR
          (c.total_counts > t.max_counts_adj * (1 - {m}) AND c.total_counts <= t.max_counts_adj) OR
          (c.pct_mito > t.max_pct_mito_adj * (1 - {m}) AND c.pct_mito <= t.max_pct_mito_adj)
        ) THEN 'near_threshold_edge'
        ELSE 'within_all_thresholds'
      END AS plan_reason
    FROM read_parquet('{qc_path.as_posix()}') c
    JOIN thr t USING(sample_id)
    """, [created_at, mref])
    con.execute(f"COPY plan TO '{out.as_posix()}' (FORMAT PARQUET)")
    con.close()
    return out


def count_parquet(path: Path, col: str) -> pd.Series:
    con = duckdb.connect()
    df = con.execute(f"""
    SELECT {col}, count(*) AS n
    FROM read_parquet('{path.as_posix()}')
    GROUP BY {col}
    ORDER BY n DESC
    """).df()
    con.close()
    return df


def md_table(df: pd.DataFrame) -> str:
    cols = [str(c) for c in df.columns]
    lines = ["| " + " | ".join(cols) + " |", "| " + " | ".join(["---"] * len(cols)) + " |"]
    for _, row in df.iterrows():
        vals = [str(row[c]) for c in df.columns]
        lines.append("| " + " | ".join(vals) + " |")
    return "\n".join(lines)


def main() -> None:
    run_root = ROOT / "results" / "v6_1" / "step2" / RUN_ID
    qdir = run_root / "03_qc"
    rdir = run_root / "10_reports"
    rdir.mkdir(parents=True, exist_ok=True)
    created_at = now_iso()
    mref = manifest_ref(run_root)

    doublet_reclassified = reclassify_doublet_calls(qdir, created_at, mref)
    rows = []
    reports = [
        "# Step2.3/2.4 QC Doublet Refilter 0513",
        "",
        f"- run_id: `{RUN_ID}`",
        f"- created_at: `{created_at}`",
        "- response fields were not read or used.",
        "- canonical Step2.3/2.4 files were not overwritten.",
        "",
        "## Branches",
    ]

    baseline_final = qdir / "final_cell_inclusion_flags.parquet"
    baseline_counts = count_parquet(baseline_final, "final_inclusion")
    reports.extend(["", "### Baseline Final Inclusion", "", md_table(baseline_counts), ""])

    for branch, cfg in BRANCHES.items():
        filter_path = write_branch_filter(qdir, branch, cfg, created_at, mref)
        final_path = finalize_from_tables(qdir, branch, filter_path, doublet_reclassified, created_at, mref)
        f_counts = count_parquet(filter_path, "filter_plan")
        i_counts = count_parquet(final_path, "final_inclusion")
        d_counts = count_parquet(final_path, "doublet_call")
        rows.append({
            "run_id": RUN_ID,
            "created_at": created_at,
            "input_manifest_ref": mref,
            "branch": branch,
            "description": cfg["description"],
            "filter_plan_path": str(filter_path),
            "final_inclusion_path": str(final_path),
            "threshold_path": str(qdir / f"per_sample_qc_thresholds.qc0513_{branch}.csv"),
            "include_main": int(i_counts.loc[i_counts["final_inclusion"] == "include_main", "n"].sum()),
            "include_sensitivity_only": int(i_counts.loc[i_counts["final_inclusion"] == "include_sensitivity_only", "n"].sum()),
            "exclude_low_quality": int(i_counts.loc[i_counts["final_inclusion"] == "exclude_low_quality", "n"].sum()),
            "exclude_high_confidence_doublet": int(i_counts.loc[i_counts["final_inclusion"] == "exclude_high_confidence_doublet", "n"].sum()),
        })
        reports.extend([
            f"### {branch}",
            "",
            cfg["description"],
            "",
            "Filter plan:",
            "",
            md_table(f_counts),
            "",
            "Final inclusion:",
            "",
            md_table(i_counts),
            "",
            "Doublet call after reclassification:",
            "",
            md_table(d_counts),
            "",
        ])

    opt = pd.DataFrame(rows)
    opt.to_csv(qdir / "qc_refilter_options_0513.csv", index=False)
    reports.extend([
        "## Recommendation",
        "",
        "`balanced` is the recommended next Step2.5 input candidate. It rescues many basic-QC excluded cells while preserving hard exclusion for high-confidence doublets and keeping borderline cells out of the main matrix.",
        "",
        "Use `max_rescue_review` only for manual sensitivity review.",
    ])
    (rdir / "step2_3_4_qc_doublet_refilter_0513_report.md").write_text("\n".join(reports), encoding="utf-8")


if __name__ == "__main__":
    main()
