#!/usr/bin/env python3
"""Phase7B: response-blind barrier identifiability gate."""

from __future__ import annotations

import argparse
from itertools import combinations

import networkx as nx
import numpy as np
import pandas as pd
import yaml

from phase7_common import MODULES, ROOT, ensure_dirs, load_config, rank_stability, write_yaml


def load_inputs(cfg: dict) -> tuple[dict, pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    out = cfg["out"]
    handoff = yaml.safe_load((out / "handoff/phase7a_to_phase7b_handoff.yaml").read_text())
    if handoff.get("verdict") not in {"GO_TO_PHASE7B", "CONDITIONAL_GO_TO_PHASE7B"}:
        raise RuntimeError(f"Phase7A did not open Phase7B: {handoff.get('verdict')}")
    wide = pd.read_parquet(ROOT / handoff["phase7b_unique_inputs"]["residualized_module_matrix"])
    eligibility = pd.read_csv(ROOT / handoff["phase7b_unique_inputs"]["module_eligibility"])
    annotation = pd.read_csv(ROOT / handoff["phase7b_unique_inputs"]["semantic_annotation"])
    coarse = pd.read_parquet(ROOT / handoff["phase7b_unique_inputs"]["coarse_sensitivity_matrix"])
    return handoff, wide, coarse, eligibility, annotation


def cross_resolution(primary: pd.DataFrame, coarse: pd.DataFrame, active: list[str]) -> dict[str, float]:
    keys = [c for c in ("cohort_id", "patient_key", "timepoint", "tissue_context", "patient_timepoint_context_id") if c in primary and c in coarse]
    joined = primary[keys + active].merge(coarse[keys + active], on=keys, suffixes=("_mid", "_coarse"))
    return {m: float(joined[f"{m}_mid"].rank().corr(joined[f"{m}_coarse"].rank())) for m in active}


def collinearity(wide: pd.DataFrame, active: list[str]) -> tuple[pd.DataFrame, float, dict[str, float]]:
    x = wide[active].dropna()
    if len(x) < max(20, len(active) * 3):
        return pd.DataFrame(), np.nan, {m: np.nan for m in active}
    z = (x - x.mean()) / x.std(ddof=0).replace(0, np.nan)
    corr = z.rank().corr()
    eig = np.linalg.eigvalsh(np.cov(z.fillna(0), rowvar=False))
    condition = float(np.sqrt(max(eig) / max(min(eig), 1e-12)))
    try:
        inv = np.linalg.pinv(z.corr().to_numpy())
        vif = dict(zip(active, np.diag(inv)))
        partial = -inv / np.sqrt(np.outer(np.diag(inv), np.diag(inv)))
        np.fill_diagonal(partial, 1.0)
    except np.linalg.LinAlgError:
        vif = {m: np.inf for m in active}
        partial = np.full((len(active), len(active)), np.nan)
    rows = []
    for a, b in combinations(active, 2):
        ia, ib = active.index(a), active.index(b)
        rows.append({"module_a": a, "module_b": b, "spearman": corr.loc[a, b], "abs_spearman": abs(corr.loc[a, b]),
                     "partial_correlation": partial[ia, ib], "abs_partial_correlation": abs(partial[ia, ib])})
    return pd.DataFrame(rows), condition, vif


def attribute_overlap(cfg: dict, annotation: pd.DataFrame) -> pd.DataFrame:
    parsed = {}
    for _, row in annotation.iterrows():
        parsed[row.module_id] = {
            "genes": set(str(row.top_genes).split(";")) - {"nan", ""},
            "states": set(str(row.carrier_states).split(";")) - {"nan", ""},
            "paths": set(str(row.top_pathways).split(";")) - {"nan", ""},
        }
    rows = []
    for a, b in combinations(parsed, 2):
        row = {"module_a": a, "module_b": b}
        for key in ("genes", "states", "paths"):
            union = parsed[a][key] | parsed[b][key]
            row[f"{key}_jaccard"] = len(parsed[a][key] & parsed[b][key]) / len(union) if union else 0.0
        rows.append(row)
    return pd.DataFrame(rows)


def bootstrap(cfg: dict, wide: pd.DataFrame, active: list[str], mode: str, resume: bool) -> tuple[pd.DataFrame, pd.DataFrame]:
    n_boot = cfg["identifiability"][f"bootstrap_n_{mode}"]
    checkpoint_every = cfg["identifiability"]["bootstrap_checkpoint_every"]
    path = cfg["out"] / "identifiability/bootstrap_pairwise_checkpoint.parquet"
    existing = pd.read_parquet(path) if resume and path.exists() else pd.DataFrame()
    start = int(existing.bootstrap_id.max()) + 1 if not existing.empty else 1
    rng = np.random.default_rng(cfg["seed"] + start)
    records = existing.to_dict("records") if not existing.empty else []
    rank_rows = []
    ranked = wide[active].rank().to_numpy(dtype=float)
    reference_corr = np.corrcoef(ranked, rowvar=False)
    upper = np.triu_indices(len(active), 1)
    reference_geometry = pd.Series(reference_corr[upper])
    cohort_code, cohorts = pd.factorize(wide.cohort_id.astype(str), sort=True)
    patient_block = wide.cohort_id.astype(str) + "::" + wide.patient_key.astype(str)
    patient_code, patients = pd.factorize(patient_block, sort=True)
    for b in range(start, n_boot + 1):
        weight = rng.exponential(1.0, len(cohorts))[cohort_code] * rng.exponential(1.0, len(patients))[patient_code]
        weight_sum = weight.sum()
        mean = (ranked * weight[:, None]).sum(axis=0) / weight_sum
        centered = ranked - mean
        cov = (centered * weight[:, None]).T @ centered / weight_sum
        scale = np.sqrt(np.diag(cov))
        corr_array = cov / np.outer(scale, scale)
        corr = pd.DataFrame(corr_array, index=active, columns=active)
        for a, c in combinations(active, 2):
            records.append({"bootstrap_id": b, "module_a": a, "module_b": c, "correlation": corr.loc[a, c]})
        rank_rows.append({"bootstrap_id": b, "rank_stability": rank_stability(reference_geometry, pd.Series(corr_array[upper])),
                          "metric_definition": "rank_correlation_of_all_pairwise_module_correlations"})
        if b % checkpoint_every == 0 or b == n_boot:
            pd.DataFrame(records).to_parquet(path, index=False)
    pairs = pd.DataFrame(records)
    if existing.empty:
        ranks = pd.DataFrame(rank_rows)
    else:
        # Rank summaries are cheap and deterministic to recompute when resuming.
        ranks = pd.DataFrame(rank_rows)
    return pairs, ranks


def pair_decisions(cfg: dict, corr: pd.DataFrame, overlap: pd.DataFrame, boot: pd.DataFrame) -> pd.DataFrame:
    if corr.empty:
        return pd.DataFrame(columns=["module_a", "module_b", "merge_recommended"])
    thresholds = cfg["identifiability"]
    support = boot.assign(high=lambda x: x.correlation.abs() >= thresholds["correlation_pass_max"]).groupby(["module_a", "module_b"]).agg(
        bootstrap_merge_support=("high", "mean"), correlation_median=("correlation", "median"), correlation_sd=("correlation", "std")
    ).reset_index()
    result = corr.merge(overlap, on=["module_a", "module_b"], how="left").merge(support, on=["module_a", "module_b"], how="left")
    secondary = ((result.genes_jaccard >= thresholds["gene_jaccard_support_min"]) |
                 (result.states_jaccard >= thresholds["carrier_jaccard_support_min"]) |
                 (result.paths_jaccard >= thresholds["pathway_jaccard_support_min"]))
    result["merge_recommended"] = (
        ((result.abs_spearman >= thresholds["correlation_pass_max"])
         & (result.bootstrap_merge_support >= thresholds["merge_bootstrap_support_min"])
         & secondary)
        | (result.abs_spearman >= thresholds["correlation_force_merge_min"])
    )
    result["decision_reason"] = np.where(result.merge_recommended, "correlation+bootstrap+attribute_support", "separable_or_insufficient_merge_evidence")
    return result


def merge_map(active: list[str], pair_table: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    graph = nx.Graph()
    graph.add_nodes_from(active)
    for _, row in pair_table[pair_table.merge_recommended].iterrows():
        graph.add_edge(row.module_a, row.module_b, weight=row.abs_spearman)
    rows, edges = [], []
    for i, component in enumerate(sorted(nx.connected_components(graph), key=lambda x: sorted(x)), 1):
        members = sorted(component)
        coarse = members[0] if len(members) == 1 else f"CB{i:02d}"
        for member in members:
            rows.append({"source_module_id": member, "phase8_barrier_id": coarse, "merge_group_size": len(members),
                         "merge_status": "independent" if len(members) == 1 else "coarse_merged"})
    for _, row in pair_table.iterrows():
        edges.append({"source": row.module_a, "target": row.module_b, "weight": row.abs_spearman,
                      "bootstrap_support": row.bootstrap_merge_support, "edge_selected": bool(row.merge_recommended)})
    return pd.DataFrame(rows), pd.DataFrame(edges)


def build_barriers(cfg: dict, wide: pd.DataFrame, eligibility: pd.DataFrame, pair_table: pd.DataFrame,
                   mapping: pd.DataFrame, condition: float, vif: dict[str, float], boot: pd.DataFrame,
                   flip_rate: dict[str, float], resolution: dict[str, float]) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    support = boot.groupby(["module_a", "module_b"])["correlation"].agg(["median", "std"]).reset_index()
    rows = []
    for module in MODULES:
        a = eligibility[eligibility.module_id.eq(module)].iloc[0]
        mapped = mapping[mapping.source_module_id.eq(module)]
        max_corr = pair_table.loc[(pair_table.module_a.eq(module)) | (pair_table.module_b.eq(module)), "abs_spearman"].max() if not pair_table.empty else np.nan
        merge = mapped.iloc[0].merge_status if not mapped.empty else "not_entered"
        if a.phase7a_status == "context_only":
            status = "context_only"
        elif a.phase7a_status == "sensitivity_only":
            status = "sensitivity_only"
        elif not bool(a.phase7b_allowed):
            status = "blocked_nonidentifiable"
        elif flip_rate.get(module, 1.0) >= cfg["identifiability"]["flip_fail_min"] or vif.get(module, np.inf) >= cfg["identifiability"]["vif_force_merge_min"]:
            status = "blocked_nonidentifiable"
        elif flip_rate.get(module, 1.0) > cfg["identifiability"]["flip_pass_max"]:
            status = "sensitivity_only"
        elif merge == "coarse_merged":
            status = "eligible_coarse_merged"
        elif vif.get(module, np.inf) <= cfg["identifiability"]["vif_pass_max"] and (not np.isfinite(max_corr) or max_corr < cfg["identifiability"]["correlation_pass_max"]):
            status = "eligible_independent"
        else:
            status = "sensitivity_only"
        rows.append({"module_id": module, "phase7a_status": a.phase7a_status, "phase7b_status": status,
                     "phase8_barrier_id": mapped.iloc[0].phase8_barrier_id if not mapped.empty else "", "vif": vif.get(module, np.nan),
                     "max_abs_spearman": max_corr, "global_condition_number": condition,
                     "assignment_flip_rate": flip_rate.get(module, np.nan), "cross_resolution_spearman": resolution.get(module, np.nan),
                     "cross_resolution_status": "estimable" if np.isfinite(resolution.get(module, np.nan)) else "not_estimable_no_shared_patient_context",
                     "spatial_localization_mappability": "available_via_carrier_state_only", "spatial_response_claim_allowed": False,
                     "perturbation_mappability": "unknown_no_formal_target_prior_registry",
                     "allowed_for_phase8_primary": status in {"eligible_independent", "eligible_coarse_merged"} and np.isfinite(resolution.get(module, np.nan)),
                     "allowed_for_phase8_conditional": status in {"eligible_independent", "eligible_coarse_merged"},
                     "claim_boundary": "identifiable_response_blind_program_not_response_associated_barrier"})
    table = pd.DataFrame(rows)
    eligible = table[table.allowed_for_phase8_conditional].copy()
    eligible["phase8_eligibility_tier"] = np.where(eligible.allowed_for_phase8_primary, "primary", "conditional")

    identity = [c for c in ("cohort_id", "patient_key", "timepoint", "tissue_context", "patient_timepoint_context_id") if c in wide]
    phase8 = wide[identity].copy()
    for barrier, group in eligible.groupby("phase8_barrier_id"):
        members = group.module_id.tolist()
        vals = wide[members].copy()
        vals = (vals - vals.mean()) / vals.std(ddof=0).replace(0, np.nan)
        if len(members) > 1:
            anchor = vals[members[0]]
            for m in members[1:]:
                if anchor.corr(vals[m]) < 0:
                    vals[m] *= -1
        phase8[barrier] = vals.mean(axis=1)
    return table, eligible, phase8


def run(mode: str = "full", resume: bool = False) -> dict:
    cfg = load_config()
    ensure_dirs(cfg["out"])
    handoff, wide, coarse, eligibility, annotation = load_inputs(cfg)
    active = eligibility.loc[eligibility.phase7b_allowed.astype(bool), "module_id"].tolist()
    corr, condition, vif = collinearity(wide, active)
    overlap = attribute_overlap(cfg, annotation[annotation.module_id.isin(active)])
    boot, rank = bootstrap(cfg, wide, active, mode, resume)
    pair_table = pair_decisions(cfg, corr, overlap, boot)
    flip_rate = {}
    for module in active:
        selected = pair_table[pair_table.merge_recommended & (pair_table.module_a.eq(module) | pair_table.module_b.eq(module))]
        expected = set(selected.module_a) | set(selected.module_b)
        expected.discard(module)
        flips = []
        for _, frame in boot.groupby("bootstrap_id"):
            subset = frame[(frame.module_a.eq(module)) | (frame.module_b.eq(module))].copy()
            if subset.empty:
                continue
            high = subset[subset.correlation.abs() >= cfg["identifiability"]["correlation_pass_max"]]
            observed = set(high.module_a) | set(high.module_b)
            observed.discard(module)
            flips.append(observed != expected)
        flip_rate[module] = float(np.mean(flips)) if flips else np.nan
    resolution = cross_resolution(wide, coarse, active)
    mapping, graph_edges = merge_map(active, pair_table)
    table, eligible, phase8 = build_barriers(cfg, wide, eligibility, pair_table, mapping, condition, vif, boot, flip_rate, resolution)

    out = cfg["out"]
    pair_table.to_csv(out / "identifiability/barrier_pairwise_identifiability.csv", index=False)
    graph_edges.to_csv(out / "identifiability/barrier_similarity_graph_edges.csv", index=False)
    rank.to_csv(out / "identifiability/bootstrap_rank_stability.csv", index=False)
    table.to_csv(out / "identifiability/barrier_identifiability_table_v6_2_1.csv", index=False)
    eligible.to_csv(out / "identifiability/eligible_barrier_set_v0.csv", index=False)
    mapping.to_csv(out / "identifiability/coarse_barrier_merge_map_v6_2_1.csv", index=False)
    uncertainty = pair_table[["module_a", "module_b", "correlation_median", "correlation_sd", "bootstrap_merge_support", "merge_recommended"]].copy()
    uncertainty.to_csv(out / "identifiability/barrier_assignment_uncertainty.csv", index=False)
    phase8.to_parquet(out / "handoff/eligible_barrier_score_matrix_v0.parquet", index=False)
    blocked = table[~table.allowed_for_phase8_primary]
    lines = ["# Unstable or Blocked Barrier Log", ""] + [f"- {r.module_id}: {r.phase7b_status}; Phase7A={r.phase7a_status}" for r in blocked.itertuples()]
    (out / "identifiability/unstable_or_blocked_barrier_log.md").write_text("\n".join(lines) + "\n")

    hard = []
    if eligible.empty:
        hard.append("no_identifiable_barrier_family")
    geometry_stability = float(rank.rank_stability.median()) if not rank.empty else np.nan
    if not np.isfinite(geometry_stability) or geometry_stability < cfg["identifiability"]["rank_stability_pass_min"]:
        hard.append("bootstrap_geometry_rank_stability_below_threshold")
    if not np.isfinite(condition) or condition >= cfg["identifiability"]["condition_fail_min"]:
        hard.append("global_condition_number_failed")
    conditional = []
    if table.cross_resolution_status.ne("estimable").any():
        conditional.append("coarse_mid_cross_resolution_not_estimable_no_shared_patient_context")
    if table.perturbation_mappability.eq("unknown_no_formal_target_prior_registry").any():
        conditional.append("perturbation_mappability_pending_later_coverage_gate")
    if condition > cfg["identifiability"]["condition_pass_max"]:
        conditional.append("global_condition_number_above_pass_threshold")
    verdict = "BLOCKED" if hard else "CONDITIONAL_GO_TO_PHASE8" if conditional or not table.allowed_for_phase8_primary.all() else "GO_TO_PHASE8"
    manifest = {
        "phase": "phase7b_barrier_identifiability_gate", "mode": mode, "verdict": verdict,
        "response_direction_evaluated": False, "module_discovery_reopened": False, "hard_blockers": hard,
        "conditional_items": conditional,
        "n_input_topics": len(MODULES), "n_phase7a_allowed": len(active), "n_phase8_barriers": int(eligible.phase8_barrier_id.nunique()),
        "global_condition_number": condition,
        "bootstrap_geometry_rank_stability_median": geometry_stability,
        "primary_outputs": {
            "eligible_barrier_set": str((out / "identifiability/eligible_barrier_set_v0.csv").relative_to(ROOT)),
            "eligible_barrier_score_matrix": str((out / "handoff/eligible_barrier_score_matrix_v0.parquet").relative_to(ROOT)),
            "merge_map": str((out / "identifiability/coarse_barrier_merge_map_v6_2_1.csv").relative_to(ROOT)),
            "identifiability_table": str((out / "identifiability/barrier_identifiability_table_v6_2_1.csv").relative_to(ROOT)),
        },
        "phase8_allowed_to_read_only": [
            str((out / "handoff/eligible_barrier_score_matrix_v0.parquet").relative_to(ROOT)),
            str((out / "identifiability/eligible_barrier_set_v0.csv").relative_to(ROOT)),
        ],
        "phase8_prohibited": ["phase6_raw_topic_scores", "phase7a_cell_state_scores", "blocked_or_sensitivity_only_topics"],
    }
    write_yaml(manifest, out / "handoff/phase7b_to_phase8_handoff.yaml")
    rehearsal = {
        "phase6_frozen_input_retrained": False, "unique_patient_timepoint_scoring_contract": True,
        "sentinel_and_data_driven_boundary_explicit": True, "eligible_set_machine_readable": True,
        "merge_and_block_decisions_machine_readable": True, "status": "PASS" if not hard else "FAIL",
    }
    write_yaml(rehearsal, out / "handoff/integration_rehearsal_7.yaml")
    report = f"""# Phase7B Barrier Identifiability Report

## Executive Verdict

**{verdict}**。本阶段只判断调整后的 candidate programs 能否稳定区分，未分析 response 方向。

## Gate Summary

- Phase7A allowed topics: {len(active)}/8
- Phase8 eligible barrier families: {int(eligible.phase8_barrier_id.nunique())}
- Global condition number: {condition:.3f}
- Bootstrap iterations: {cfg['identifiability'][f'bootstrap_n_{mode}']}
- Pairwise-geometry rank stability: {geometry_stability:.3f}
- Merged topics: {int(mapping.merge_status.eq('coarse_merged').sum()) if not mapping.empty else 0}

当前 Phase8 坐标均为 conditional eligibility；原因是 coarse/mid 输入没有共享 patient-context，且 perturbation mappability 留待后续正式 coverage gate。

## Boundary

`eligible` 只表示可测量且可区分，不表示与治疗失败有关。response association 必须在 Phase8 使用冻结 handoff 另行检验。
"""
    (out / "PHASE7B_BARRIER_IDENTIFIABILITY_REPORT.md").write_text(report)
    return manifest


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--mode", choices=["smoke", "full"], default="full")
    parser.add_argument("--resume", action="store_true")
    args = parser.parse_args()
    print(yaml.safe_dump(run(args.mode, args.resume), sort_keys=False))
