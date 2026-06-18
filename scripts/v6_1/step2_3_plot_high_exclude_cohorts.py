from __future__ import annotations

import argparse
from pathlib import Path

import duckdb
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

ROOT = Path(__file__).resolve().parents[2]


def main(run_id: str, root_dir: Path | str = ROOT) -> dict:
    root = Path(root_dir)
    run_root = root / "results" / "v6_1" / "step2" / run_id
    qdir = run_root / "03_qc"
    ra_dir = run_root / "10_reports" / "step2_3_high_exclude_cohort_reaudit"
    fig_dir = ra_dir / "figures"
    fig_dir.mkdir(parents=True, exist_ok=True)

    s = pd.read_csv(ra_dir / "high_exclude_cohort_reaudit_summary.csv")
    cohorts = s["cohort_id"].astype(str).tolist()

    con = duckdb.connect()
    con.execute("PRAGMA threads=8")
    cohort_list_sql = ",".join(["'" + c.replace("'", "''") + "'" for c in cohorts])
    df = con.execute(
        f"""
        WITH cell_d AS (
          SELECT * EXCLUDE (rn) FROM (
            SELECT *, row_number() OVER (PARTITION BY cell_barcode, sample_id ORDER BY created_at DESC) AS rn
            FROM read_parquet('{(qdir / "cell_qc_metrics.parquet").as_posix()}')
          ) t WHERE rn=1
        ),
        inc_d AS (
          SELECT * EXCLUDE (rn) FROM (
            SELECT *, row_number() OVER (PARTITION BY cell_barcode, sample_id ORDER BY created_at DESC) AS rn
            FROM read_parquet('{(qdir / "final_cell_inclusion_flags.parquet").as_posix()}')
          ) t WHERE rn=1
        )
        SELECT
          c.cohort_id,
          c.total_counts::DOUBLE AS total_counts,
          c.n_genes_by_counts::DOUBLE AS n_genes_by_counts,
          c.pct_mito::DOUBLE AS pct_mito,
          i.final_inclusion
        FROM cell_d c
        JOIN inc_d i USING(cell_barcode, sample_id)
        WHERE c.cohort_id IN ({cohort_list_sql})
        """
    ).df()

    # Sample to keep figures light
    sampled = (
        df.groupby(["cohort_id", "final_inclusion"], group_keys=False)
        .apply(lambda x: x.sample(min(len(x), 30000), random_state=42))
        .reset_index(drop=True)
    )

    sns.set_theme(style="whitegrid")

    # Boxplots
    for metric, ylim in [
        ("total_counts", None),
        ("n_genes_by_counts", None),
        ("pct_mito", (0, min(100, sampled["pct_mito"].quantile(0.99) * 1.1))),
    ]:
        plt.figure(figsize=(16, 6))
        sns.boxplot(
            data=sampled,
            x="cohort_id",
            y=metric,
            hue="final_inclusion",
            showfliers=False,
        )
        plt.xticks(rotation=75, ha="right", fontsize=8)
        if ylim is not None:
            plt.ylim(*ylim)
        plt.title(f"{metric} Boxplot by Cohort and Inclusion")
        plt.tight_layout()
        plt.savefig(fig_dir / f"boxplot_{metric}_by_cohort_inclusion.png", dpi=180)
        plt.close()

    # Violin (trimmed for readability)
    for metric in ["total_counts", "n_genes_by_counts", "pct_mito"]:
        plt.figure(figsize=(16, 6))
        sns.violinplot(
            data=sampled,
            x="cohort_id",
            y=metric,
            hue="final_inclusion",
            cut=0,
            scale="width",
            inner="quartile",
        )
        plt.xticks(rotation=75, ha="right", fontsize=8)
        plt.title(f"{metric} Violin by Cohort and Inclusion")
        plt.tight_layout()
        plt.savefig(fig_dir / f"violin_{metric}_by_cohort_inclusion.png", dpi=180)
        plt.close()

    # Exclude fraction bar
    sb = s.sort_values("exclude_frac", ascending=False)
    plt.figure(figsize=(12, 5))
    sns.barplot(data=sb, x="cohort_id", y="exclude_frac", color="#c0392b")
    plt.xticks(rotation=75, ha="right")
    plt.ylim(0, 1)
    plt.title("Exclude Low-Quality Fraction (High-Exclude Cohorts)")
    plt.tight_layout()
    plt.savefig(fig_dir / "bar_exclude_frac_high_exclude_cohorts.png", dpi=180)
    plt.close()

    return {"run_id": run_id, "cohorts": len(cohorts), "fig_dir": str(fig_dir)}


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description="Plot high-exclude cohorts")
    ap.add_argument("--run-id", required=True)
    ap.add_argument("--root-dir", default=str(ROOT))
    args = ap.parse_args()
    print(main(args.run_id, args.root_dir))
