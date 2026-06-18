#!/usr/bin/env python3
"""Score Phase6 module benchmark outputs across methods."""
from __future__ import annotations

import argparse
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

from phase6_algorithm_common import (
    OUT,
    PHASE6,
    TOP_N,
    load_primary_tables,
    load_signatures,
    max_jaccard,
    top_gene_sets,
    write_csv,
    write_md,
    write_yaml,
)


def read_yaml(path: Path) -> dict:
    if not path.exists():
        return {}
    with path.open() as fh:
        return yaml.safe_load(fh) or {}


def current_v0() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    membership = pd.read_csv(PHASE6 / "modules/module_membership.v0.csv", dtype=str)
    membership["method"] = "current_nmf_v0"
    membership["run_id"] = "current_nmf_v0__rank8"
    membership["module_id"] = membership["module_id"].map(lambda x: f"current_nmf_v0__{x}")
    score = pd.read_csv(PHASE6 / "modules/module_score_matrix.v0.csv", dtype=str)
    for col in list(score.columns):
        if col.startswith("module_M"):
            score = score.rename(columns={col: f"current_nmf_v0__{col}"})
    score.insert(0, "run_id", "current_nmf_v0__rank8")
    score.insert(0, "method", "current_nmf_v0")
    stab = pd.read_csv(PHASE6 / "modules/module_stability_metrics.csv", dtype=str)
    selected = stab[stab["selected_for_freeze"].eq("yes")].copy()
    qc = pd.DataFrame(
        [
            {
                "method": "current_nmf_v0",
                "run_id": "current_nmf_v0__rank8",
                "run_status": "complete",
                "rank": selected["rank"].iloc[0] if not selected.empty else "",
                "median_bootstrap_jaccard": selected["median_bootstrap_jaccard"].iloc[0] if not selected.empty else "",
                "reconstruction_error": selected["reconstruction_error"].iloc[0] if not selected.empty else "",
            }
        ]
    )
    return membership, score, qc


def collect_method_outputs() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    memberships = []
    scores = []
    qcs = []
    m0, s0, q0 = current_v0()
    memberships.append(m0)
    scores.append(s0)
    qcs.append(q0)
    methods_dir = OUT / "methods"
    if methods_dir.exists():
        for method_dir in sorted(methods_dir.iterdir()):
            if not method_dir.is_dir():
                continue
            method = method_dir.name
            m_path = method_dir / f"module_membership.{method}.csv"
            q_path = method_dir / f"module_run_qc.{method}.csv"
            if q_path.exists():
                qcs.append(pd.read_csv(q_path, dtype=str))
            if m_path.exists():
                memberships.append(pd.read_csv(m_path, dtype=str))
            for s_path in [method_dir / f"module_score_matrix.{method}.parquet", method_dir / f"module_score_matrix.{method}.csv"]:
                if s_path.exists():
                    scores.append(pd.read_parquet(s_path) if s_path.suffix == ".parquet" else pd.read_csv(s_path, dtype=str))
                    break
    return (
        pd.concat(memberships, ignore_index=True, sort=False) if memberships else pd.DataFrame(),
        pd.concat(scores, ignore_index=True, sort=False) if scores else pd.DataFrame(),
        pd.concat(qcs, ignore_index=True, sort=False) if qcs else pd.DataFrame(),
    )


def signature_scores(membership: pd.DataFrame) -> pd.DataFrame:
    sigs = load_signatures()
    rows = []
    if membership.empty:
        return pd.DataFrame()
    top = membership[membership["is_top_gene"].astype(str).eq("yes")].copy()
    for (method, module_id), group in top.groupby(["method", "module_id"], sort=False):
        genes = {str(g).upper() for g in group["gene"]}
        best_name = ""
        best_j = 0.0
        for name, sig_genes in sigs.items():
            union = len(genes | sig_genes)
            j = len(genes & sig_genes) / union if union else 0.0
            if j > best_j:
                best_name, best_j = name, j
        rows.append({"method": method, "module_id": module_id, "best_signature": best_name, "best_signature_jaccard": best_j})
    return pd.DataFrame(rows)


def redundancy_scores(membership: pd.DataFrame) -> pd.DataFrame:
    rows = []
    top = membership[membership["is_top_gene"].astype(str).eq("yes")].copy()
    for method, group in top.groupby("method", sort=False):
        module_sets = []
        for module_id, mod in group.groupby("module_id", sort=False):
            module_sets.append((module_id, {str(g) for g in mod["gene"]}))
        vals = []
        for (_, a), (_, b) in combinations(module_sets, 2):
            union = len(a | b)
            vals.append(len(a & b) / union if union else 0.0)
        rows.append(
            {
                "method": method,
                "n_modules": len(module_sets),
                "median_inter_module_jaccard": float(np.median(vals)) if vals else 0.0,
                "max_inter_module_jaccard": float(np.max(vals)) if vals else 0.0,
            }
        )
    return pd.DataFrame(rows)


def cohort_robustness(scores: pd.DataFrame) -> pd.DataFrame:
    rows = []
    if scores.empty:
        return pd.DataFrame()
    _, _, _, genes = load_primary_tables()
    meta_cols = {"method", "run_id", "expression_unit_id", "cohort_id", "object_id", "sample_key", "patient_key", "cell_state_level", "cell_state", "layer_family"}
    for method, group in scores.groupby("method", sort=False):
        if "cohort_id" not in group.columns:
            rows.append({"method": method, "max_cohort_fraction": 1.0, "cohort_dominance_flag": "yes", "n_score_rows": len(group)})
            continue
        n = len(group)
        max_frac = group["cohort_id"].value_counts(normalize=True).max() if n else 1.0
        module_cols = [c for c in group.columns if c not in meta_cols and str(c).startswith(method)]
        rows.append(
            {
                "method": method,
                "max_cohort_fraction": float(max_frac),
                "cohort_dominance_flag": "yes" if max_frac > 0.5 else "no",
                "n_score_rows": n,
                "n_score_modules": len(module_cols),
            }
        )
    return pd.DataFrame(rows)


def leakage_audit(scores: pd.DataFrame) -> pd.DataFrame:
    banned = ["response", "responder", "outcome", "survival", "progression", "recist", "split", "train", "test"]
    rows = []
    for col in scores.columns:
        lower = str(col).lower()
        hits = [b for b in banned if b in lower]
        rows.append({"column": col, "leakage_tokens": "|".join(hits), "leakage_risk": "yes" if hits else "no"})
    return pd.DataFrame(rows)


def method_summary(qc: pd.DataFrame, bio: pd.DataFrame, red: pd.DataFrame, rob: pd.DataFrame, leak: pd.DataFrame) -> pd.DataFrame:
    rows = []
    methods = sorted(set(qc.get("method", pd.Series(dtype=str)).dropna().astype(str)))
    for method in methods:
        q = qc[qc["method"].eq(method)].copy()
        complete = q[q.get("run_status", "").astype(str).eq("complete")]
        if complete.empty:
            reason = q.get("blocked_reason", pd.Series(["unknown"])).iloc[0] if not q.empty else "no_output"
            rows.append({"method": method, "selection_status": "blocked", "blocked_reason": reason})
            continue
        q["median_bootstrap_jaccard_num"] = pd.to_numeric(q.get("median_bootstrap_jaccard", ""), errors="coerce")
        q["_selection_stability"] = q["median_bootstrap_jaccard_num"].fillna(-1.0)
        best_q = q.sort_values("_selection_stability", ascending=False).iloc[0]
        b = bio[bio["method"].eq(method)]
        r = red[red["method"].eq(method)]
        c = rob[rob["method"].eq(method)]
        bio_frac = float((pd.to_numeric(b["best_signature_jaccard"], errors="coerce").fillna(0) > 0).mean()) if not b.empty else 0.0
        red_j = float(r["median_inter_module_jaccard"].iloc[0]) if not r.empty else 1.0
        max_cohort = float(c["max_cohort_fraction"].iloc[0]) if not c.empty else 1.0
        stability_raw = pd.to_numeric(pd.Series([best_q.get("median_bootstrap_jaccard_num", np.nan)]), errors="coerce").iloc[0]
        stability = 0.0 if pd.isna(stability_raw) else float(stability_raw)
        leak_fail = leak["leakage_risk"].eq("yes").any()
        redundancy_score = max(0.0, 1.0 - min(red_j, 1.0))
        cohort_score = max(0.0, 1.0 - min(max_cohort, 1.0))
        plausibility_score = 1.0
        rank = pd.to_numeric(pd.Series([best_q.get("rank", np.nan)]), errors="coerce").iloc[0]
        if pd.notna(rank) and (rank < 6 or rank > 24):
            plausibility_score = 0.5
        composite_score = (
            0.45 * stability
            + 0.15 * bio_frac
            + 0.15 * redundancy_score
            + 0.10 * cohort_score
            + 0.05 * plausibility_score
        )
        passes = stability >= 0.45 and red_j < 0.35 and max_cohort < 0.5 and bio_frac >= 0.7 and not leak_fail
        rows.append(
            {
                "method": method,
                "selection_status": "passes_freeze_threshold" if passes else "candidate_or_conditional",
                "selected_run": best_q.get("run_id", ""),
                "median_bootstrap_jaccard": stability,
                "bootstrap_jaccard_ci_low": best_q.get("bootstrap_jaccard_ci_low", ""),
                "median_bootstrap_jaccard_top80": best_q.get("median_bootstrap_jaccard_top80", ""),
                "biological_mappability_fraction": bio_frac,
                "median_inter_module_jaccard": red_j,
                "max_cohort_fraction": max_cohort,
                "composite_score": composite_score,
                "leakage_detected": str(leak_fail).lower(),
            }
        )
    return pd.DataFrame(rows)


def run(mode: str) -> None:
    (OUT / "score").mkdir(parents=True, exist_ok=True)
    membership, scores, qc = collect_method_outputs()
    write_csv(qc, OUT / "score/algorithm_run_registry.csv")
    write_csv(membership, OUT / "score/module_membership.all_methods.csv")
    if not scores.empty:
        scores.to_parquet(OUT / "score/module_score_matrix.all_methods.parquet", index=False)
    bio = signature_scores(membership)
    red = redundancy_scores(membership)
    rob = cohort_robustness(scores)
    leak = leakage_audit(scores)
    summary = method_summary(qc, bio, red, rob, leak)
    write_csv(bio, OUT / "score/module_biological_mappability_rescore.csv")
    write_csv(red, OUT / "score/module_redundancy_rescore.csv")
    write_csv(rob, OUT / "score/module_cross_cohort_robustness.csv")
    write_csv(leak, OUT / "score/module_algorithm_leakage_audit.csv")
    write_csv(summary, OUT / "score/module_stability_rescore.csv")
    selected = summary[summary["selection_status"].eq("passes_freeze_threshold")]
    verdict = "GO_WITH_SELECTED_V1" if not selected.empty else "KEEP_V0_CONDITIONAL_WITH_BENCHMARKED_ALTERNATIVES"
    selected_method = selected.sort_values(["composite_score", "median_bootstrap_jaccard"], ascending=False)["method"].iloc[0] if not selected.empty else "current_nmf_v0"
    run_scope = {
        "scope": "strengthened" if mode == "strengthened" else "full_constrained" if mode == "full" else "smoke",
        "note": (
            "Official packages were used where required. Strengthened mode expands "
            "rank/prior grids, adds bootstrap evidence, optimizes GeneNMF meta-programs, "
            "and prepares one downstream module entry."
        ),
    }
    write_yaml(
        {
            "phase": "phase6_module_algorithm_benchmark",
            "mode": mode,
            "run_scope": run_scope,
            "verdict": verdict,
            "selected_method": selected_method,
            "official_methods_required_not_reimplemented": True,
            "outputs": {
                "algorithm_run_registry": str(OUT / "score/algorithm_run_registry.csv"),
                "module_stability_rescore": str(OUT / "score/module_stability_rescore.csv"),
                "module_membership_all_methods": str(OUT / "score/module_membership.all_methods.csv"),
                "module_score_matrix_all_methods": str(OUT / "score/module_score_matrix.all_methods.parquet"),
                "downstream_entry_manifest": str(OUT / "handoff/phase6_to_phase7_module_entry_manifest.yaml"),
            },
        },
        OUT / "score/module_algorithm_selection_manifest.yaml",
    )
    complete_runs = int(qc["run_status"].astype(str).eq("complete").sum()) if "run_status" in qc.columns else 0
    leakage_detected = bool(leak["leakage_risk"].eq("yes").any()) if not leak.empty else False
    blocked = qc[qc.get("run_status", pd.Series(dtype=str)).astype(str).eq("blocked")].copy()
    method_lines = []
    for row in summary.sort_values(["selection_status", "method"]).to_dict("records"):
        method_lines.append(
            "- `{method}`: `{status}`, run `{run}`, composite `{comp}`, stability `{stab}`, bio-map `{bio}`, redundancy `{red}`, max cohort `{cohort}`".format(
                method=row.get("method", ""),
                status=row.get("selection_status", ""),
                run=row.get("selected_run", ""),
                comp="" if pd.isna(row.get("composite_score")) else round(float(row.get("composite_score")), 3),
                stab="" if pd.isna(row.get("median_bootstrap_jaccard")) else round(float(row.get("median_bootstrap_jaccard")), 3),
                bio="" if pd.isna(row.get("biological_mappability_fraction")) else round(float(row.get("biological_mappability_fraction")), 3),
                red="" if pd.isna(row.get("median_inter_module_jaccard")) else round(float(row.get("median_inter_module_jaccard")), 3),
                cohort="" if pd.isna(row.get("max_cohort_fraction")) else round(float(row.get("max_cohort_fraction")), 3),
            )
        )
    blocked_lines = []
    for row in blocked.to_dict("records"):
        blocked_lines.append(f"- `{row.get('method', '')}`: `{row.get('blocked_reason', row.get('run_status', ''))}`")
    if not blocked_lines:
        blocked_lines.append("- None")
    write_md(
        "# Phase6 Module Algorithm Benchmark Report\n\n"
        f"- Mode: `{mode}`\n"
        f"- Run scope: `{run_scope['scope']}`\n"
        f"- Verdict: `{verdict}`\n"
        f"- Selected method: `{selected_method}`\n"
        f"- Methods scored: `{summary['method'].nunique() if not summary.empty else 0}`\n"
        f"- Complete runs: `{complete_runs}`\n"
        f"- Leakage detected: `{str(leakage_detected).lower()}`\n"
        "- Official methods with missing dependencies are recorded as blocked, not replaced by local implementations.\n\n"
        "## Method Summary\n\n"
        + "\n".join(method_lines)
        + "\n\n## Blocked Runs\n\n"
        + "\n".join(blocked_lines)
        + "\n\n## Scope Note\n\n"
        + run_scope["note"]
        + "\n",
        OUT / "PHASE6_MODULE_ALGORITHM_BENCHMARK_REPORT.md",
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--mode", choices=["smoke", "full", "strengthened"], default="full")
    run(parser.parse_args().mode)
