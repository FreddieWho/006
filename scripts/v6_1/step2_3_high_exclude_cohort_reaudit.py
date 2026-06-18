from __future__ import annotations

import argparse
from datetime import datetime, timezone
from pathlib import Path

import duckdb
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]


def now_iso() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat()


def main(run_id: str, root_dir: Path | str = ROOT) -> dict:
    root = Path(root_dir)
    run_root = root / "results" / "v6_1" / "step2" / run_id
    qdir = run_root / "03_qc"
    out_dir = run_root / "10_reports" / "step2_3_high_exclude_cohort_reaudit"
    out_dir.mkdir(parents=True, exist_ok=True)
    per_dir = out_dir / "per_cohort"
    per_dir.mkdir(parents=True, exist_ok=True)

    manifest_ref = str(run_root / "00_manifest" / "step2_run_manifest.patched.yaml")
    if not Path(manifest_ref).exists():
        manifest_ref = str(run_root / "00_manifest" / "step2_run_manifest.yaml")

    inc = pd.read_parquet(qdir / "final_cell_inclusion_flags.parquet", columns=["cohort_id", "sample_id", "cell_barcode", "final_inclusion"])
    base = (
        inc.groupby("cohort_id", as_index=False)
        .agg(n_cells=("cell_barcode", "size"), n_exclude=("final_inclusion", lambda s: (s == "exclude_low_quality").sum()))
    )
    base["exclude_frac"] = base["n_exclude"] / base["n_cells"]
    target = sorted(base.loc[base["exclude_frac"] > 0.2, "cohort_id"].astype(str).tolist())

    con = duckdb.connect()
    con.execute("PRAGMA threads=8")
    con.execute(f"""
    CREATE OR REPLACE TABLE cell_d AS
    SELECT * EXCLUDE (rn) FROM (
      SELECT *, row_number() OVER (PARTITION BY cell_barcode, sample_id ORDER BY created_at DESC) AS rn
      FROM read_parquet('{(qdir / "cell_qc_metrics.parquet").as_posix()}')
    ) t WHERE rn=1
    """)
    con.execute(f"""
    CREATE OR REPLACE TABLE thr AS
    SELECT sample_id,
           min_genes::DOUBLE AS min_genes,
           min_counts::DOUBLE AS min_counts,
           max_genes::DOUBLE AS max_genes,
           max_counts::DOUBLE AS max_counts,
           max_pct_mito::DOUBLE AS max_pct_mito
    FROM read_csv_auto('{(qdir / "per_sample_qc_thresholds.csv").as_posix()}', header=true)
    """)
    con.execute(f"""
    CREATE OR REPLACE TABLE inc_d AS
    SELECT * EXCLUDE (rn) FROM (
      SELECT *, row_number() OVER (PARTITION BY cell_barcode, sample_id ORDER BY created_at DESC) AS rn
      FROM read_parquet('{(qdir / "final_cell_inclusion_flags.parquet").as_posix()}')
    ) t WHERE rn=1
    """)
    con.execute("""
    CREATE OR REPLACE TABLE m AS
    SELECT c.cohort_id, c.sample_id, c.cell_barcode,
           c.total_counts::DOUBLE AS total_counts,
           c.n_genes_by_counts::DOUBLE AS n_genes_by_counts,
           c.pct_mito::DOUBLE AS pct_mito,
           t.min_genes, t.min_counts, t.max_genes, t.max_counts, t.max_pct_mito,
           i.final_inclusion
    FROM cell_d c
    JOIN thr t USING(sample_id)
    JOIN inc_d i USING(cell_barcode, sample_id)
    """)

    summary_rows = []
    for cohort in target:
        c = cohort.replace("'", "''")
        agg = con.execute(f"""
        SELECT
          '{cohort}' AS cohort_id,
          count(*) AS n_cells,
          sum(CASE WHEN final_inclusion='exclude_low_quality' THEN 1 ELSE 0 END) AS n_ex_low,
          sum(CASE WHEN n_genes_by_counts < min_genes THEN 1 ELSE 0 END) AS fail_min_genes,
          sum(CASE WHEN total_counts < min_counts THEN 1 ELSE 0 END) AS fail_min_counts,
          sum(CASE WHEN n_genes_by_counts > max_genes THEN 1 ELSE 0 END) AS fail_max_genes,
          sum(CASE WHEN total_counts > max_counts THEN 1 ELSE 0 END) AS fail_max_counts,
          sum(CASE WHEN pct_mito > max_pct_mito THEN 1 ELSE 0 END) AS fail_mito,
          quantile_cont(total_counts, 0.5) AS q50_counts,
          quantile_cont(total_counts, 0.1) AS q10_counts,
          quantile_cont(n_genes_by_counts, 0.5) AS q50_genes,
          quantile_cont(n_genes_by_counts, 0.1) AS q10_genes,
          quantile_cont(pct_mito, 0.5) AS q50_mito,
          quantile_cont(pct_mito, 0.9) AS q90_mito
        FROM m
        WHERE cohort_id='{c}'
        """).df().iloc[0].to_dict()
        ex = con.execute(f"""
        SELECT
          sum(CASE WHEN n_genes_by_counts < min_genes THEN 1 ELSE 0 END) AS ex_fail_min_genes,
          sum(CASE WHEN total_counts < min_counts THEN 1 ELSE 0 END) AS ex_fail_min_counts,
          sum(CASE WHEN n_genes_by_counts > max_genes THEN 1 ELSE 0 END) AS ex_fail_max_genes,
          sum(CASE WHEN total_counts > max_counts THEN 1 ELSE 0 END) AS ex_fail_max_counts,
          sum(CASE WHEN pct_mito > max_pct_mito THEN 1 ELSE 0 END) AS ex_fail_mito
        FROM m
        WHERE cohort_id='{c}' AND final_inclusion='exclude_low_quality'
        """).df().iloc[0].to_dict()
        agg.update(ex)
        n_ex = max(int(agg["n_ex_low"]), 1)
        agg["exclude_frac"] = float(agg["n_ex_low"]) / float(agg["n_cells"])
        for k in ["ex_fail_min_genes", "ex_fail_min_counts", "ex_fail_max_genes", "ex_fail_max_counts", "ex_fail_mito"]:
            agg[k + "_share"] = float(agg[k]) / n_ex
        # sample-level anomalies
        samp = con.execute(f"""
        SELECT sample_id,
               count(*) AS n_cells,
               quantile_cont(total_counts,0.5) AS med_counts,
               quantile_cont(n_genes_by_counts,0.5) AS med_genes,
               quantile_cont(pct_mito,0.5) AS med_mito
        FROM m WHERE cohort_id='{c}'
        GROUP BY sample_id
        ORDER BY n_cells DESC
        """).df()
        anomaly = samp[(samp["n_cells"] >= 10000) & ((samp["med_counts"] <= 2) | (samp["med_genes"] <= 2))]
        agg["n_large_samples_low_complexity"] = int(len(anomaly))
        summary_rows.append(agg)
        # per-cohort files
        samp.to_csv(per_dir / f"{cohort}.sample_distribution.csv", index=False)
        anomaly.to_csv(per_dir / f"{cohort}.anomaly_samples.csv", index=False)
        dom = max(
            [
                ("min_counts", agg["ex_fail_min_counts_share"]),
                ("min_genes", agg["ex_fail_min_genes_share"]),
                ("mito", agg["ex_fail_mito_share"]),
                ("max_counts", agg["ex_fail_max_counts_share"]),
                ("max_genes", agg["ex_fail_max_genes_share"]),
            ],
            key=lambda x: x[1],
        )[0]
        rec = "mapping_or_counts_layer_audit_first" if dom in {"min_counts", "min_genes"} and agg["exclude_frac"] > 0.6 else (
            "mito_threshold_sensitivity_branch" if dom == "mito" else "upper_tail_and_doublet_review"
        )
        md = [
            f"# Cohort Reaudit: {cohort}",
            "",
            f"- run_id: `{run_id}`",
            f"- created_at: `{now_iso()}`",
            f"- input_manifest_ref: `{manifest_ref}`",
            "",
            "## Core Stats",
            f"- n_cells: {int(agg['n_cells'])}",
            f"- exclude_low_quality: {int(agg['n_ex_low'])}",
            f"- exclude_frac: {agg['exclude_frac']:.4f}",
            "",
            "## Excluded Trigger Shares",
            f"- min_counts: {agg['ex_fail_min_counts_share']:.4f}",
            f"- min_genes: {agg['ex_fail_min_genes_share']:.4f}",
            f"- mito: {agg['ex_fail_mito_share']:.4f}",
            f"- max_counts: {agg['ex_fail_max_counts_share']:.4f}",
            f"- max_genes: {agg['ex_fail_max_genes_share']:.4f}",
            "",
            "## Distribution Snapshot",
            f"- q10/q50 counts: {agg['q10_counts']:.2f} / {agg['q50_counts']:.2f}",
            f"- q10/q50 genes: {agg['q10_genes']:.2f} / {agg['q50_genes']:.2f}",
            f"- q50/q90 mito: {agg['q50_mito']:.2f} / {agg['q90_mito']:.2f}",
            "",
            "## Anomaly Check",
            f"- n_large_samples_low_complexity: {agg['n_large_samples_low_complexity']}",
            "",
            "## Suggested Action",
            f"- dominant_trigger: `{dom}`",
            f"- recommendation: `{rec}`",
        ]
        (per_dir / f"{cohort}.reaudit.md").write_text("\n".join(md), encoding="utf-8")

    out = pd.DataFrame(summary_rows).sort_values("exclude_frac", ascending=False)
    out.insert(0, "run_id", run_id)
    out.insert(1, "created_at", now_iso())
    out.insert(2, "input_manifest_ref", manifest_ref)
    out.to_csv(out_dir / "high_exclude_cohort_reaudit_summary.csv", index=False)

    idx = [
        "# Step2.3 High-Exclude Cohort Reaudit Index",
        "",
        f"- run_id: `{run_id}`",
        f"- created_at: `{now_iso()}`",
        f"- cohorts reviewed (`exclude_frac > 0.2`): {len(target)}",
        "",
        "## Summary Table",
        "",
        f"- [high_exclude_cohort_reaudit_summary.csv]({(out_dir / 'high_exclude_cohort_reaudit_summary.csv').as_posix()})",
        "",
        "## Per-Cohort Reports",
    ]
    for c in target:
        idx.append(f"- [{c}.reaudit.md]({(per_dir / f'{c}.reaudit.md').as_posix()})")
    (out_dir / "index.md").write_text("\n".join(idx), encoding="utf-8")
    return {"run_id": run_id, "n_cohorts": len(target), "out_dir": str(out_dir)}


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description="Reaudit cohorts with exclude_frac > 0.2")
    ap.add_argument("--run-id", required=True)
    ap.add_argument("--root-dir", default=str(ROOT))
    args = ap.parse_args()
    print(main(run_id=args.run_id, root_dir=args.root_dir))

