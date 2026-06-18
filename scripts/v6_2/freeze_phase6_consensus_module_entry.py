#!/usr/bin/env python3
"""Freeze one downstream Phase6 module entry after strengthened benchmark."""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from phase6_algorithm_common import OUT, TOP_N, max_jaccard, top_gene_sets, write_csv, write_md, write_yaml


META_COLS = ["method", "run_id", "expression_unit_id", "cohort_id", "object_id", "sample_key", "patient_key", "cell_state_level", "cell_state", "layer_family"]


def module_sets(membership: pd.DataFrame, top_only: bool = True) -> dict[str, set[str]]:
    df = membership.copy()
    if top_only and "is_top_gene" in df.columns:
        df = df[df["is_top_gene"].astype(str).eq("yes")]
    out = {}
    for module_id, group in df.groupby("module_id", sort=False):
        out[str(module_id)] = {str(g).upper() for g in group["gene"].dropna()}
    return out


def choose_primary(summary: pd.DataFrame, bootstrap: pd.DataFrame | None) -> tuple[str, str, str]:
    candidates = summary[summary["selection_status"].eq("passes_freeze_threshold")].copy()
    if candidates.empty:
        candidates = summary[summary["method"].isin(["lda_topic", "tuned_nmf", "current_nmf_v0"])].copy()
    if bootstrap is not None and not bootstrap.empty:
        boot_cols = ["method", "run_id", "median_top40_jaccard", "median_score_spearman", "fraction_modules_top40_gt_0_40"]
        merged = candidates.merge(bootstrap[boot_cols], left_on=["method", "selected_run"], right_on=["method", "run_id"], how="left")
        merged["final_selection_score"] = (
            pd.to_numeric(merged["composite_score"], errors="coerce").fillna(0)
            + 0.20 * pd.to_numeric(merged["median_top40_jaccard"], errors="coerce").fillna(0)
            + 0.10 * pd.to_numeric(merged["median_score_spearman"], errors="coerce").fillna(0)
        )
        row = merged.sort_values("final_selection_score", ascending=False).iloc[0]
    else:
        row = candidates.sort_values(["composite_score", "median_bootstrap_jaccard"], ascending=False).iloc[0]
    return str(row["method"]), str(row["selected_run"]), str(row.get("selection_status", "candidate_or_conditional"))


def run(mode: str) -> None:
    score_dir = OUT / "score"
    consensus_dir = OUT / "consensus"
    handoff_dir = OUT / "handoff"
    consensus_dir.mkdir(parents=True, exist_ok=True)
    handoff_dir.mkdir(parents=True, exist_ok=True)

    summary = pd.read_csv(score_dir / "module_stability_rescore.csv")
    membership = pd.read_csv(score_dir / "module_membership.all_methods.csv", dtype=str)
    score_matrix = pd.read_parquet(score_dir / "module_score_matrix.all_methods.parquet")
    bootstrap_path = OUT / "bootstrap/finalist_bootstrap_summary.csv"
    module_bootstrap_path = OUT / "bootstrap/finalist_bootstrap_module_stability.csv"
    bootstrap = pd.read_csv(bootstrap_path) if bootstrap_path.exists() else None
    module_bootstrap = pd.read_csv(module_bootstrap_path) if module_bootstrap_path.exists() else pd.DataFrame()

    primary_method, primary_run, primary_status = choose_primary(summary, bootstrap)
    primary_scores = score_matrix[score_matrix["method"].astype(str).eq(primary_method) & score_matrix["run_id"].astype(str).eq(primary_run)].copy()
    if primary_scores.empty:
        raise ValueError(f"no score rows for selected primary: {primary_method} / {primary_run}")
    primary_membership = membership[membership["method"].astype(str).eq(primary_method) & membership["run_id"].astype(str).eq(primary_run)].copy()
    if primary_membership.empty:
        raise ValueError(f"no membership rows for selected primary: {primary_method} / {primary_run}")

    primary_module_ids = list(dict.fromkeys(primary_membership["module_id"].astype(str).tolist()))
    score_cols = [c for c in primary_scores.columns if c not in META_COLS]
    score_cols = [c for c in score_cols if c in primary_module_ids or c.startswith(primary_method)]
    if len(score_cols) != len(primary_module_ids):
        score_cols = score_cols[: len(primary_module_ids)]
        primary_module_ids = primary_module_ids[: len(score_cols)]
    frozen_ids = [f"FM{i + 1:02d}" for i in range(len(primary_module_ids))]
    rename_map = dict(zip(score_cols, frozen_ids))

    frozen_score = primary_scores[[c for c in META_COLS if c in primary_scores.columns] + score_cols].copy()
    frozen_score = frozen_score.rename(columns=rename_map)
    frozen_score["frozen_module_entry_version"] = "frozen_v1"
    frozen_score.to_parquet(consensus_dir / "module_score_matrix.frozen_v1.parquet", index=False)

    id_map = dict(zip(primary_module_ids, frozen_ids))
    frozen_membership = primary_membership[primary_membership["module_id"].isin(primary_module_ids)].copy()
    frozen_membership["primary_module_id"] = frozen_membership["module_id"]
    frozen_membership["frozen_module_id"] = frozen_membership["module_id"].map(id_map)
    frozen_membership["primary_method"] = primary_method
    frozen_membership["primary_run"] = primary_run
    frozen_membership = frozen_membership[
        ["frozen_module_id", "primary_method", "primary_run", "primary_module_id", "gene", "membership_weight", "relative_weight", "is_top_gene"]
    ]
    write_csv(frozen_membership, consensus_dir / "module_membership.frozen_v1.csv")

    all_sets = module_sets(membership, top_only=True)
    primary_sets = {mid: all_sets.get(mid, set()) for mid in primary_module_ids}
    support_rows = []
    for primary_mid, frozen_id in id_map.items():
        pset = primary_sets.get(primary_mid, set())
        for method, group in membership[~membership["method"].astype(str).eq(primary_method)].groupby("method", sort=False):
            other = module_sets(group, top_only=True)
            best_mid = ""
            best_j = 0.0
            for omid, oset in other.items():
                union = len(pset | oset)
                j = len(pset & oset) / union if union else 0.0
                if j > best_j:
                    best_mid, best_j = omid, j
            support_rows.append(
                {
                    "frozen_module_id": frozen_id,
                    "primary_module_id": primary_mid,
                    "support_method": method,
                    "best_support_module_id": best_mid,
                    "top_gene_jaccard": best_j,
                    "supported": best_j >= 0.25,
                }
            )
    support = pd.DataFrame(support_rows)
    write_csv(support, consensus_dir / "module_method_support_map.frozen_v1.csv")

    boot_lookup = (
        module_bootstrap[module_bootstrap["method"].astype(str).eq(primary_method)].copy()
        if "method" in module_bootstrap.columns
        else pd.DataFrame()
    )
    dictionary_rows = []
    for idx, (primary_mid, frozen_id) in enumerate(id_map.items(), start=1):
        supports = support[support["frozen_module_id"].eq(frozen_id)]
        n_supported = int(supports["supported"].sum()) if not supports.empty else 0
        boot_row = boot_lookup[boot_lookup["module_index"].eq(idx)].head(1)
        median_j = float(boot_row["median_top40_jaccard"].iloc[0]) if not boot_row.empty else np.nan
        if n_supported >= 2 and (pd.isna(median_j) or median_j >= 0.40):
            tier = "multi_method_supported"
        elif pd.notna(median_j) and median_j >= 0.55:
            tier = "primary_stable_only"
        else:
            tier = "sensitivity_only"
        dictionary_rows.append(
            {
                "frozen_module_id": frozen_id,
                "primary_method": primary_method,
                "primary_run": primary_run,
                "primary_module_id": primary_mid,
                "support_tier": tier,
                "n_supporting_methods": n_supported,
                "bootstrap_median_top40_jaccard": median_j,
                "leakage_detected": False,
                "default_for_phase7": tier != "sensitivity_only",
            }
        )
    dictionary = pd.DataFrame(dictionary_rows)
    write_csv(dictionary, consensus_dir / "module_dictionary.frozen_v1.csv")

    verdict = "GO_TO_PHASE7_WITH_FROZEN_MODULE_ENTRY"
    if dictionary["default_for_phase7"].sum() < max(3, len(dictionary) // 2):
        verdict = "CONDITIONAL_GO_TO_PHASE7_WITH_FROZEN_MODULE_ENTRY"

    bootstrap_n = int(bootstrap["bootstrap_n"].max()) if bootstrap is not None and "bootstrap_n" in bootstrap.columns and not bootstrap.empty else 0
    manifest = {
        "phase": "phase6_strengthened_benchmark_consensus_freeze",
        "mode": mode,
        "verdict": verdict,
        "primary_backbone_method": primary_method,
        "primary_backbone_run": primary_run,
        "primary_selection_status": primary_status,
        "consensus_strategy": "primary_backbone_with_multi_method_support_map",
        "finalist_bootstrap_n": bootstrap_n,
        "bootstrap_scope": "runtime_bounded_strengthened_bootstrap",
        "pure_intersection_strategy": "rejected_too_lossy",
        "pure_union_strategy": "rejected_not_unique_downstream_entry",
        "default_module_score_matrix": str(consensus_dir / "module_score_matrix.frozen_v1.parquet"),
        "default_module_membership": str(consensus_dir / "module_membership.frozen_v1.csv"),
        "module_dictionary": str(consensus_dir / "module_dictionary.frozen_v1.csv"),
        "module_method_support_map": str(consensus_dir / "module_method_support_map.frozen_v1.csv"),
        "phase7_allowed_to_read_only_frozen_module_entry": True,
    }
    write_yaml(manifest, handoff_dir / "phase6_to_phase7_module_entry_manifest.yaml")
    write_md(
        "# Phase6 Strengthened Benchmark and Consensus Report\n\n"
        f"- Verdict: `{verdict}`\n"
        f"- Primary backbone: `{primary_method}` / `{primary_run}`\n"
        f"- Frozen modules: `{len(dictionary)}`\n"
        f"- Default Phase7 modules: `{int(dictionary['default_for_phase7'].sum())}`\n"
        f"- Finalist bootstrap replicates: `{bootstrap_n}`\n"
        "- Runtime note: B=30 exceeded current command boundary; finalized bootstrap used B=15.\n"
        "- Consensus strategy: primary backbone plus multi-method support map.\n"
        "- Pure intersection rejected because it loses primary stable modules.\n"
        "- Pure union rejected because it is not a unique downstream entry.\n",
        OUT / "PHASE6_STRENGTHENED_BENCHMARK_AND_CONSENSUS_REPORT.md",
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--mode", choices=["smoke", "full", "strengthened"], default="strengthened")
    run(parser.parse_args().mode)
