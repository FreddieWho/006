from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns


RUN_ID = "step2_v6_1_0505_0319"
RUN_ROOT = Path("/home/huyudi/006/results/v6_1/step2") / RUN_ID
CONS_DIR = RUN_ROOT / "04_annotation" / "consensus_annotation_v1"
FIG_DIR = CONS_DIR / "figures"
TABLE_DIR = CONS_DIR / "figure_tables"

LINEAGE_ORDER = [
    "T_NK",
    "Myeloid_DC",
    "B_Plasma",
    "Mast_or_minor_immune",
    "Immune",
    "Non_immune_Epithelial_Tumor",
    "Non_immune_Stromal",
    "Non_immune_Endothelial",
    "Unknown",
]

MID_ORDER = [
    "CD8_T",
    "CD4_T",
    "Treg",
    "NK",
    "T_NK_core",
    "B_naive_memory",
    "B_Plasma_core",
    "Plasma",
    "Mono_FCN1",
    "Macro_C1QC",
    "Macro_SPP1",
    "Myeloid_DC_core",
    "cDC1",
    "cDC2",
    "pDC",
    "Mast",
    "Mast_minor_immune",
    "Immune_unspecified",
]

DECISION_ORDER = [
    "CONS_R1_marker_original_agree",
    "CONS_R2_marker_only",
    "CONS_R3_original_supported",
    "CONS_R4_conflict_degrade",
    "CONS_R5_evidence_insufficient",
    "CONS_R6_auto_support_only",
]

LINEAGE_COLORS = {
    "T_NK": "#2C7FB8",
    "Myeloid_DC": "#41AB5D",
    "B_Plasma": "#756BB1",
    "Mast_or_minor_immune": "#C51B8A",
    "Immune": "#9E9AC8",
    "Non_immune_Epithelial_Tumor": "#D95F0E",
    "Non_immune_Stromal": "#8C6D31",
    "Non_immune_Endothelial": "#1B9E77",
    "Unknown": "#BDBDBD",
}

DECISION_COLORS = {
    "CONS_R1_marker_original_agree": "#1B9E77",
    "CONS_R2_marker_only": "#66A61E",
    "CONS_R3_original_supported": "#7570B3",
    "CONS_R4_conflict_degrade": "#D95F02",
    "CONS_R5_evidence_insufficient": "#BDBDBD",
    "CONS_R6_auto_support_only": "#E6AB02",
}


def save_figure(fig: plt.Figure, name: str) -> None:
    fig.tight_layout()
    for ext in ["png", "pdf"]:
        fig.savefig(FIG_DIR / f"{name}.{ext}", dpi=220, bbox_inches="tight")
    plt.close(fig)


def fraction_table(counts: pd.DataFrame, index_cols: list[str], value_col: str) -> pd.DataFrame:
    total = counts.groupby(index_cols)["n_cells"].transform("sum")
    out = counts.copy()
    out[value_col] = out["n_cells"] / total
    return out


def plot_quality_overview(summary: pd.DataFrame) -> None:
    df = summary.copy()
    df["log10_n_main"] = np.log10(df["n_main"].clip(lower=1))
    df.to_csv(TABLE_DIR / "quality_overview_by_cohort.csv", index=False)

    fig, ax = plt.subplots(figsize=(9.5, 6.2))
    sc = ax.scatter(
        df["major_lineage_coverage"],
        df["immune_mid_state_specified_fraction"],
        s=20 + 65 * df["log10_n_main"],
        c=df["conflict_fraction"],
        cmap="mako_r",
        alpha=0.82,
        edgecolor="#2B2B2B",
        linewidth=0.4,
    )
    ax.axvline(0.90, color="#444444", linestyle="--", linewidth=1)
    ax.axhline(0.80, color="#444444", linestyle="--", linewidth=1)
    ax.set_xlim(0, 1.02)
    ax.set_ylim(0, 1.02)
    ax.set_xlabel("Major lineage coverage")
    ax.set_ylabel("Immune mid-state specified fraction")
    ax.set_title("Consensus Annotation Coverage by Cohort")
    cb = fig.colorbar(sc, ax=ax)
    cb.set_label("Conflict fraction")
    low = df.sort_values(["major_lineage_coverage", "immune_mid_state_specified_fraction"]).head(8)
    for r in low.itertuples():
        ax.text(r.major_lineage_coverage + 0.01, r.immune_mid_state_specified_fraction + 0.01, r.cohort_id, fontsize=7)
    save_figure(fig, "01_quality_overview_by_cohort")


def plot_lineage_composition(lineage_counts: pd.DataFrame, summary: pd.DataFrame) -> None:
    main_counts = lineage_counts[lineage_counts["final_inclusion"] == "include_main"].copy()
    frac = fraction_table(main_counts, ["cohort_id"], "fraction")
    wide = (
        frac.pivot_table(index="cohort_id", columns="final_major_lineage", values="fraction", fill_value=0)
        .reindex(columns=LINEAGE_ORDER, fill_value=0)
        .join(summary.set_index("cohort_id")["n_main"])
        .sort_values("n_main", ascending=False)
        .drop(columns=["n_main"])
    )
    wide.to_csv(TABLE_DIR / "major_lineage_fraction_by_cohort.csv")

    fig, ax = plt.subplots(figsize=(12.5, 10))
    left = np.zeros(len(wide))
    y = np.arange(len(wide))
    for lineage in LINEAGE_ORDER:
        vals = wide[lineage].to_numpy()
        ax.barh(y, vals, left=left, color=LINEAGE_COLORS[lineage], label=lineage, height=0.82)
        left += vals
    ax.set_yticks(y)
    ax.set_yticklabels(wide.index, fontsize=7)
    ax.invert_yaxis()
    ax.set_xlim(0, 1)
    ax.set_xlabel("Fraction of include_main cells")
    ax.set_title("Major Lineage Composition by Cohort")
    ax.legend(ncol=3, fontsize=8, loc="lower center", bbox_to_anchor=(0.5, -0.18), frameon=False)
    save_figure(fig, "02_major_lineage_composition_by_cohort")


def plot_midstate_heatmap(mid_counts: pd.DataFrame, summary: pd.DataFrame) -> None:
    main = mid_counts[mid_counts["final_inclusion"] == "include_main"].copy()
    immune = main[main["final_immune_mid_state"].isin(MID_ORDER)].copy()
    frac = fraction_table(immune, ["cohort_id"], "fraction")
    wide = (
        frac.pivot_table(index="cohort_id", columns="final_immune_mid_state", values="fraction", fill_value=0)
        .reindex(columns=MID_ORDER, fill_value=0)
        .join(summary.set_index("cohort_id")["n_main"])
        .sort_values("n_main", ascending=False)
        .drop(columns=["n_main"])
    )
    wide.to_csv(TABLE_DIR / "immune_mid_state_fraction_by_cohort.csv")

    fig, ax = plt.subplots(figsize=(12.5, 9))
    sns.heatmap(wide, ax=ax, cmap="viridis", vmin=0, vmax=max(0.05, float(wide.max().max())), linewidths=0.0)
    ax.set_xlabel("Final immune mid-state")
    ax.set_ylabel("Cohort")
    ax.set_title("Immune Mid-State Composition by Cohort")
    ax.tick_params(axis="x", labelrotation=45, labelsize=8)
    ax.tick_params(axis="y", labelsize=7)
    save_figure(fig, "03_immune_mid_state_heatmap_by_cohort")


def plot_decision_rules(decision_counts: pd.DataFrame, summary: pd.DataFrame) -> None:
    main = decision_counts[decision_counts["final_inclusion"] == "include_main"].copy()
    frac = fraction_table(main, ["cohort_id"], "fraction")
    wide = (
        frac.pivot_table(index="cohort_id", columns="decision_rule_id", values="fraction", fill_value=0)
        .reindex(columns=DECISION_ORDER, fill_value=0)
        .join(summary.set_index("cohort_id")["n_main"])
        .sort_values("n_main", ascending=False)
        .drop(columns=["n_main"])
    )
    wide.to_csv(TABLE_DIR / "decision_rule_fraction_by_cohort.csv")

    fig, ax = plt.subplots(figsize=(12.5, 10))
    left = np.zeros(len(wide))
    y = np.arange(len(wide))
    for rule in DECISION_ORDER:
        vals = wide[rule].to_numpy()
        ax.barh(y, vals, left=left, color=DECISION_COLORS[rule], label=rule, height=0.82)
        left += vals
    ax.set_yticks(y)
    ax.set_yticklabels(wide.index, fontsize=7)
    ax.invert_yaxis()
    ax.set_xlim(0, 1)
    ax.set_xlabel("Fraction of include_main cells")
    ax.set_title("Consensus Decision Rule Composition by Cohort")
    ax.legend(ncol=2, fontsize=7, loc="lower center", bbox_to_anchor=(0.5, -0.17), frameon=False)
    save_figure(fig, "04_decision_rule_composition_by_cohort")


def plot_conflict_and_samples(summary: pd.DataFrame, sample_summary: pd.DataFrame) -> None:
    top = summary.sort_values("conflict_fraction", ascending=False).head(20)
    top.to_csv(TABLE_DIR / "top20_conflict_fraction_cohorts.csv", index=False)

    fig, ax = plt.subplots(figsize=(9.5, 6.8))
    ax.barh(top["cohort_id"], top["conflict_fraction"], color="#D95F02")
    ax.invert_yaxis()
    ax.set_xlabel("Conflict fraction among include_main cells")
    ax.set_title("Top Cohorts by Annotation Conflict Fraction")
    save_figure(fig, "05_top20_conflict_fraction_cohorts")

    sample_summary.to_csv(TABLE_DIR / "sample_coverage_distribution.csv", index=False)
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.8), sharey=False)
    sns.histplot(sample_summary["major_lineage_coverage"], bins=30, color="#2C7FB8", ax=axes[0])
    axes[0].axvline(0.90, color="#444444", linestyle="--", linewidth=1)
    axes[0].set_xlabel("Major lineage coverage")
    axes[0].set_title("Sample-Level Major Coverage")
    sns.histplot(sample_summary["immune_mid_state_specified_fraction"], bins=30, color="#41AB5D", ax=axes[1])
    axes[1].axvline(0.80, color="#444444", linestyle="--", linewidth=1)
    axes[1].set_xlabel("Immune mid-state specified fraction")
    axes[1].set_title("Sample-Level Immune Mid-State Coverage")
    save_figure(fig, "06_sample_level_coverage_distributions")


def write_index(figures: list[str]) -> None:
    lines = [
        "# Step2.5 Consensus Annotation Figures",
        "",
        f"- run_id: `{RUN_ID}`",
        "- branch: `consensus_annotation_v1`",
        "- source table: `cell_state_annotation_consensus_v1.parquet`",
        "- scope: `include_main` unless stated otherwise",
        "- response labels were not used.",
        "",
        "## Figures",
    ]
    lines.extend(f"- `{name}.png` / `{name}.pdf`" for name in figures)
    lines.extend(
        [
            "",
            "## Source Tables",
            "- `figure_tables/quality_overview_by_cohort.csv`",
            "- `figure_tables/major_lineage_fraction_by_cohort.csv`",
            "- `figure_tables/immune_mid_state_fraction_by_cohort.csv`",
            "- `figure_tables/decision_rule_fraction_by_cohort.csv`",
            "- `figure_tables/top20_conflict_fraction_cohorts.csv`",
            "- `figure_tables/sample_coverage_distribution.csv`",
        ]
    )
    (FIG_DIR / "index.md").write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    TABLE_DIR.mkdir(parents=True, exist_ok=True)
    sns.set_theme(style="whitegrid", context="notebook")

    summary = pd.read_csv(CONS_DIR / "annotation_summary_by_cohort.csv")
    sample_summary = pd.read_csv(CONS_DIR / "annotation_summary_by_sample.csv")
    cols = [
        "cohort_id",
        "final_inclusion",
        "final_major_lineage",
        "final_immune_mid_state",
        "decision_rule_id",
    ]
    annot = pd.read_parquet(CONS_DIR / "cell_state_annotation_consensus_v1.parquet", columns=cols)

    lineage_counts = (
        annot.groupby(["cohort_id", "final_inclusion", "final_major_lineage"], dropna=False)
        .size()
        .reset_index(name="n_cells")
    )
    mid_counts = (
        annot.groupby(["cohort_id", "final_inclusion", "final_immune_mid_state"], dropna=False)
        .size()
        .reset_index(name="n_cells")
    )
    decision_counts = (
        annot.groupby(["cohort_id", "final_inclusion", "decision_rule_id"], dropna=False)
        .size()
        .reset_index(name="n_cells")
    )
    lineage_counts.to_csv(TABLE_DIR / "major_lineage_counts_by_cohort.csv", index=False)
    mid_counts.to_csv(TABLE_DIR / "immune_mid_state_counts_by_cohort.csv", index=False)
    decision_counts.to_csv(TABLE_DIR / "decision_rule_counts_by_cohort.csv", index=False)

    plot_quality_overview(summary)
    plot_lineage_composition(lineage_counts, summary)
    plot_midstate_heatmap(mid_counts, summary)
    plot_decision_rules(decision_counts, summary)
    plot_conflict_and_samples(summary, sample_summary)

    figures = [
        "01_quality_overview_by_cohort",
        "02_major_lineage_composition_by_cohort",
        "03_immune_mid_state_heatmap_by_cohort",
        "04_decision_rule_composition_by_cohort",
        "05_top20_conflict_fraction_cohorts",
        "06_sample_level_coverage_distributions",
    ]
    write_index(figures)
    print(f"Wrote {len(figures)} figure sets to {FIG_DIR}")


if __name__ == "__main__":
    main()
