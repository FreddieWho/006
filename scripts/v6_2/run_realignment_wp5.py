#!/usr/bin/env python3
"""Response-blind WP5 representation necessity evaluation.

This script deliberately accepts only the frozen Phase7 wide module matrix and
non-response cohort metadata. It never loads response/outcome columns.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
from typing import Dict, Iterable, List, Tuple

import numpy as np
import pandas as pd
import yaml


MODULES = [f"FM{i:02d}" for i in range(1, 9)]
META = ["cohort_id", "patient_key", "timepoint", "tissue_context", "patient_timepoint_context_id"]
FORBIDDEN_TOKENS = ("response", "outcome", "label", "endpoint", "responder", "failure")
GROUPS = {
    "FM01": ["FM01"],
    "CB02": ["FM02", "FM03"],
    "FM04": ["FM04"],
    "CB04": ["FM05", "FM08"],
    "FM06": ["FM06"],
    "FM07": ["FM07"],
}


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def rank_array(values: np.ndarray) -> np.ndarray:
    order = np.argsort(values, kind="mergesort")
    ranks = np.empty(len(values), dtype=float)
    ranks[order] = np.arange(len(values), dtype=float)
    return ranks


def safe_corr(a: np.ndarray, b: np.ndarray) -> float:
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    if a.size < 3 or np.std(a) == 0 or np.std(b) == 0:
        return 1.0 if np.allclose(a, b) else 0.0
    return float(np.corrcoef(a, b)[0, 1])


def spearman(a: np.ndarray, b: np.ndarray) -> float:
    return safe_corr(rank_array(a), rank_array(b))


def eta_squared(values: np.ndarray, groups: Iterable[str]) -> float:
    values = np.asarray(values, dtype=float)
    groups = np.asarray(list(groups), dtype=object)
    if values.ndim == 1:
        values = values[:, None]
    total = np.sum((values - np.mean(values, axis=0)) ** 2, axis=0)
    between = np.zeros(values.shape[1], dtype=float)
    for group in pd.unique(groups):
        mask = groups == group
        if mask.sum() == 0:
            continue
        between += mask.sum() * (np.mean(values[mask], axis=0) - np.mean(values, axis=0)) ** 2
    ratios = np.divide(between, total, out=np.zeros_like(between), where=total > 0)
    return float(np.mean(ratios))


def stable_fold(value: str, n: int = 5) -> int:
    return int(hashlib.sha256(str(value).encode("utf-8")).hexdigest()[:12], 16) % n


def load_inputs(wide_path: Path, registry_path: Path) -> Tuple[pd.DataFrame, pd.DataFrame]:
    columns = META + MODULES
    df = pd.read_parquet(wide_path, columns=columns)
    bad_columns = [c for c in df.columns if any(token in c.lower() for token in FORBIDDEN_TOKENS)]
    if bad_columns:
        raise RuntimeError(f"response-lock violation in matrix columns: {bad_columns}")
    missing = [c for c in columns if c not in df.columns]
    if missing:
        raise RuntimeError(f"missing frozen matrix columns: {missing}")
    registry = pd.read_csv(
        registry_path,
        usecols=["cohort_id", "disease", "cancer_group", "data_modality"],
        dtype={"cohort_id": str, "disease": str, "cancer_group": str, "data_modality": str},
    ).drop_duplicates("cohort_id")
    # These are context/availability fields only; no response labels are loaded.
    df = df.merge(registry, on="cohort_id", how="left", validate="many_to_one")
    for column in ["disease", "cancer_group", "data_modality"]:
        df[column] = df[column].fillna("unknown").astype(str)
    return df, registry


class Route:
    def __init__(self, route_id: str):
        self.route_id = route_id
        self.mean = None
        self.scale = None
        self.components = None
        self.blocks = None

    def fit(self, x: np.ndarray) -> "Route":
        self.mean = np.nanmean(x, axis=0)
        self.scale = np.nanstd(x, axis=0)
        self.scale[self.scale == 0] = 1.0
        z = np.where(np.isfinite(x), x, self.mean)
        z = (z - self.mean) / self.scale
        if self.route_id == "linear_latent":
            _, _, vt = np.linalg.svd(z, full_matrices=False)
            self.components = vt[: min(4, z.shape[1])]
        elif self.route_id == "srb_stage_a":
            self.blocks = []
            for block in GROUPS.values():
                indices = [MODULES.index(m) for m in block]
                if len(indices) == 1:
                    self.blocks.append((indices, np.ones((1, 1))))
                    continue
                _, _, vt = np.linalg.svd(z[:, indices], full_matrices=False)
                component = vt[:1]
                # Canonical sign prevents seed/fold sign flips in audit files.
                if component[0, 0] < 0:
                    component = -component
                self.blocks.append((indices, component))
        return self

    def _standardize(self, x: np.ndarray) -> np.ndarray:
        z = np.where(np.isfinite(x), x, self.mean)
        return (z - self.mean) / self.scale

    def transform_reconstruct(self, x: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
        z = self._standardize(x)
        if self.route_id == "direct_module":
            return z, z.copy()
        if self.route_id == "linear_latent":
            latent = z @ self.components.T
            return latent, latent @ self.components
        if self.route_id == "hierarchical_module":
            latent_parts = []
            recon = np.zeros_like(z)
            for block in GROUPS.values():
                indices = [MODULES.index(m) for m in block]
                value = np.mean(z[:, indices], axis=1, keepdims=True)
                latent_parts.append(value)
                recon[:, indices] = value
            return np.concatenate(latent_parts, axis=1), recon
        if self.route_id == "srb_stage_a":
            latent_parts = []
            recon = np.zeros_like(z)
            for indices, component in self.blocks:
                block_z = z[:, indices]
                latent = block_z @ component.T
                latent_parts.append(latent)
                recon[:, indices] = latent @ component
            return np.concatenate(latent_parts, axis=1), recon
        raise ValueError(self.route_id)

    def interpretability(self) -> float:
        if self.route_id in {"direct_module", "hierarchical_module"}:
            return 1.0
        if self.route_id == "srb_stage_a":
            return 1.0  # every latent is confined to one frozen module/block
        loadings = np.abs(self.components)
        return float(np.mean(np.max(loadings, axis=1)))


def metric_row(route: Route, x_test: np.ndarray, meta: pd.DataFrame, split: str, fold: int, train_n: int) -> Dict[str, object]:
    latent, recon = route.transform_reconstruct(x_test)
    z_test = route._standardize(x_test)
    mse = float(np.mean((z_test - recon) ** 2))
    rank_scores = [spearman(z_test[:, j], recon[:, j]) for j in range(z_test.shape[1])]
    block_scores = []
    for block in GROUPS.values():
        indices = [MODULES.index(m) for m in block]
        block_scores.append(spearman(np.mean(z_test[:, indices], axis=1), np.mean(recon[:, indices], axis=1)))
    rng = np.random.default_rng(20260821 + fold)
    masked = x_test.copy()
    mask = rng.random(masked.shape) < 0.10
    masked[mask] = np.nan
    _, masked_recon = route.transform_reconstruct(masked)
    mask_mse = float(np.mean((z_test - masked_recon) ** 2))
    return {
        "route_id": route.route_id,
        "split": split,
        "fold": fold,
        "n_train": train_n,
        "n_test": len(x_test),
        "effective_units": latent.shape[1],
        "reconstruction_mse": mse,
        "reconstruction_fidelity": float(1.0 / (1.0 + mse)),
        "module_rank_preservation": float(np.mean(rank_scores)),
        "coarse_block_preservation": float(np.mean(block_scores)),
        "mask_sensitivity_mse": mask_mse,
        "mask_sensitivity_delta": float(mask_mse - mse),
        "cohort_leakage_eta2": eta_squared(latent, meta["cohort_id"]),
        "cancer_leakage_eta2": eta_squared(latent, meta["cancer_group"]),
        "platform_leakage_eta2": np.nan,
        "interpretability_score": route.interpretability(),
        "coverage_fraction": float(np.isfinite(x_test).mean()),
    }


def bootstrap_stability(x: np.ndarray, route_id: str, n_bootstrap: int = 12) -> float:
    reconstructions = []
    rng = np.random.default_rng(20260821)
    for _ in range(n_bootstrap):
        indices = rng.integers(0, len(x), size=len(x))
        route = Route(route_id).fit(x[indices])
        _, recon = route.transform_reconstruct(x)
        reconstructions.append(recon.ravel())
    correlations = []
    for i in range(len(reconstructions)):
        for j in range(i + 1, len(reconstructions)):
            correlations.append(safe_corr(reconstructions[i], reconstructions[j]))
    return float(np.median(correlations))


def run(args: argparse.Namespace) -> None:
    wide_path = Path(args.wide_matrix)
    registry_path = Path(args.registry)
    outdir = Path(args.output_dir)
    outdir.mkdir(parents=True, exist_ok=True)
    df, _ = load_inputs(wide_path, registry_path)
    x = df[MODULES].to_numpy(dtype=float)
    route_ids = ["direct_module", "linear_latent", "hierarchical_module", "srb_stage_a"]
    rows: List[Dict[str, object]] = []
    split_specs = []
    for split, column in [("patient_grouped", "patient_key"), ("cohort_held_out", "cohort_id"), ("cancer_held_out", "cancer_group")]:
        folds = df[column].map(stable_fold).to_numpy()
        split_specs.append((split, folds))
    for split, folds in split_specs:
        for fold in range(5):
            test = folds == fold
            if test.sum() < 3 or (~test).sum() < 10:
                continue
            for route_id in route_ids:
                route = Route(route_id).fit(x[~test])
                rows.append(metric_row(route, x[test], df.loc[test], split, fold, int((~test).sum())))
    metrics = pd.DataFrame(rows)
    if metrics.empty:
        raise RuntimeError("no valid response-blind folds")
    metrics.to_csv(outdir / "representation_baseline_comparison.csv", index=False)
    leakage_cols = ["route_id", "split", "fold", "cohort_leakage_eta2", "cancer_leakage_eta2", "platform_leakage_eta2"]
    metrics[leakage_cols].to_csv(outdir / "representation_leakage_audit.csv", index=False)
    interp_cols = ["route_id", "split", "fold", "module_rank_preservation", "coarse_block_preservation", "interpretability_score", "effective_units"]
    metrics[interp_cols].to_csv(outdir / "representation_interpretability_audit.csv", index=False)
    stability_rows = []
    for route_id in route_ids:
        stability_rows.append({"route_id": route_id, "bootstrap_n": 12, "bootstrap_reconstruction_stability": bootstrap_stability(x, route_id)})
    stability = pd.DataFrame(stability_rows)
    stability.to_csv(outdir / "representation_stability_audit.csv", index=False)
    aggregate = metrics.groupby("route_id", as_index=False).agg(
        reconstruction_mse=("reconstruction_mse", "mean"),
        reconstruction_fidelity=("reconstruction_fidelity", "mean"),
        module_rank_preservation=("module_rank_preservation", "mean"),
        coarse_block_preservation=("coarse_block_preservation", "mean"),
        mask_sensitivity_delta=("mask_sensitivity_delta", "mean"),
        cohort_leakage_eta2=("cohort_leakage_eta2", "mean"),
        cancer_leakage_eta2=("cancer_leakage_eta2", "mean"),
        interpretability_score=("interpretability_score", "mean"),
        effective_units=("effective_units", "mean"),
    ).merge(stability, on="route_id", how="left")
    direct = aggregate.loc[aggregate.route_id == "direct_module"].iloc[0]
    verdicts = []
    for _, row in aggregate.iterrows():
        fidelity_noninferior = bool(row.reconstruction_fidelity >= direct.reconstruction_fidelity - 0.05)
        leakage_advantage = bool(
            row.cohort_leakage_eta2 + row.cancer_leakage_eta2
            <= direct.cohort_leakage_eta2 + direct.cancer_leakage_eta2 - 0.05
        )
        interpretability_noninferior = bool(row.interpretability_score >= direct.interpretability_score - 0.10)
        stability_noninferior = bool(row.bootstrap_reconstruction_stability >= direct.bootstrap_reconstruction_stability - 0.05)
        verdicts.append({
            "route_id": row.route_id,
            "fidelity_noninferior_to_direct": fidelity_noninferior,
            "leakage_advantage_vs_direct": leakage_advantage,
            "interpretability_noninferior": interpretability_noninferior,
            "stability_noninferior": stability_noninferior,
            "passes_srb_adoption_gate": bool(fidelity_noninferior and leakage_advantage and interpretability_noninferior and stability_noninferior),
        })
    verdict_df = pd.DataFrame(verdicts)
    aggregate.to_csv(outdir / "representation_route_aggregate.csv", index=False)
    srb_pass = bool(verdict_df.loc[verdict_df.route_id == "srb_stage_a", "passes_srb_adoption_gate"].iloc[0])
    selected = "srb_stage_a" if srb_pass else "direct_module"
    route_verdict = {
        "schema_version": "v1",
        "status": "COMPLETE_RESPONSE_BLIND",
        "run_id": args.run_id,
        "response_lock": "LOCKED",
        "response_values_read": False,
        "outcome_values_read": False,
        "routes_compared": route_ids,
        "split_schemes": [x[0] for x in split_specs],
        "selection_rule": {
            "fidelity_noninferiority_tolerance": 0.05,
            "leakage_advantage_minimum_eta2_sum": 0.05,
            "interpretability_noninferiority_tolerance": 0.10,
            "stability_noninferiority_tolerance": 0.05,
        },
        "srb_adoption": "ADOPT" if srb_pass else "DO_NOT_ADOPT",
        "selected_backbone": selected,
        "verdict_reason": "SRB must improve leakage while preserving fidelity, interpretability and stability; otherwise direct module remains the transparent backbone.",
        "aggregate_metrics": json.loads(aggregate.round(6).to_json(orient="records")),
        "route_verdicts": json.loads(verdict_df.to_json(orient="records")),
        "input_hashes": {"wide_matrix": sha256(wide_path), "cohort_registry": sha256(registry_path)},
    }
    with (outdir / "representation_route_verdict.yaml").open("w", encoding="utf-8") as handle:
        yaml.safe_dump(route_verdict, handle, sort_keys=False, allow_unicode=True)
    manifest = {
        "schema_version": "v1",
        "status": "FROZEN_RESPONSE_BLIND_REPRESENTATION",
        "run_id": args.run_id,
        "selected_backbone": selected,
        "representation_route_verdict": "representation_route_verdict.yaml",
        "module_ids": MODULES,
        "n_patient_timepoint_contexts": int(len(df)),
        "n_unique_patients": int(df.patient_key.nunique()),
        "n_cohorts": int(df.cohort_id.nunique()),
        "response_lock": True,
        "response_values_read": False,
        "outcome_values_read": False,
        "membership_hash": sha256(Path(args.membership)),
        "input_hashes": {"wide_matrix": sha256(wide_path), "cohort_registry": sha256(registry_path)},
        "route_outputs": [
            "representation_baseline_comparison.csv",
            "representation_leakage_audit.csv",
            "representation_interpretability_audit.csv",
            "representation_stability_audit.csv",
            "representation_route_aggregate.csv",
            "representation_route_verdict.yaml",
        ],
        "limitations": [
            "SRB Stage-A here is a response-blind module-factorized linear candidate, not the later nonlinear supervised SRB head.",
            "Platform metadata was unavailable in the selected non-response registry and is marked not evaluated.",
            "No response/outcome values were loaded or used.",
        ],
    }
    with (outdir / "frozen_representation_manifest.yaml").open("w", encoding="utf-8") as handle:
        yaml.safe_dump(manifest, handle, sort_keys=False, allow_unicode=True)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--wide-matrix", required=True)
    parser.add_argument("--registry", required=True)
    parser.add_argument("--membership", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--run-id", required=True)
    run(parser.parse_args())


if __name__ == "__main__":
    main()
