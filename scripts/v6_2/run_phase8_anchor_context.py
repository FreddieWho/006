#!/usr/bin/env python3
"""v6.2.1 Phase8: anchor preflight, bounded projections, and transparent association."""

from __future__ import annotations

import argparse
import hashlib
import math
import os
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

from phase7_common import ROOT, bh_fdr, build_design, normalize_timepoint, rank_stability, robust_z, sha256, write_yaml


CONFIG = ROOT / "scripts/v6_2/phase8_v6_2_1_config.yaml"


def load_config() -> dict:
    config_path = Path(os.environ.get("PHASE8_CONFIG_PATH", CONFIG))
    cfg = yaml.safe_load(config_path.read_text())
    cfg["out"] = ROOT / cfg["output_dir"]
    cfg["inputs"] = {k: ROOT / v for k, v in cfg["inputs"].items()}
    for d in ("preflight", "projection", "hcc_sensitivity", "association", "meta", "robustness", "negative_controls", "audit", "handoff"):
        (cfg["out"] / d).mkdir(parents=True, exist_ok=True)
    return cfg


def auc_rank(y: np.ndarray, score: np.ndarray) -> float:
    valid = np.isfinite(y) & np.isfinite(score)
    y, score = y[valid], score[valid]
    n1, n0 = int((y == 1).sum()), int((y == 0).sum())
    if n1 == 0 or n0 == 0:
        return np.nan
    ranks = pd.Series(score).rank(method="average").to_numpy()
    return float((ranks[y == 1].sum() - n1 * (n1 + 1) / 2) / (n1 * n0))


def spearman(a: pd.Series, b: pd.Series) -> float:
    x = pd.concat([pd.to_numeric(a, errors="coerce"), pd.to_numeric(b, errors="coerce")], axis=1).dropna()
    return float(x.iloc[:, 0].rank().corr(x.iloc[:, 1].rank())) if len(x) >= 3 else np.nan


def hedges_g(data: pd.DataFrame, feature: str) -> dict:
    r = pd.to_numeric(data.loc[data.response_binary == 1, feature], errors="coerce").dropna()
    nr = pd.to_numeric(data.loc[data.response_binary == 0, feature], errors="coerce").dropna()
    n1, n0 = len(r), len(nr)
    if min(n1, n0) < 2:
        return {"effect": np.nan, "se": np.nan, "n_R": n1, "n_NR": n0}
    df = n1 + n0 - 2
    pooled = math.sqrt(((n1 - 1) * r.var(ddof=1) + (n0 - 1) * nr.var(ddof=1)) / df)
    if not np.isfinite(pooled) or pooled <= 1e-12:
        return {"effect": np.nan, "se": np.nan, "n_R": n1, "n_NR": n0}
    d = (r.mean() - nr.mean()) / pooled
    correction = 1 - 3 / (4 * df - 1) if df > 1 else 1.0
    g = correction * d
    variance = (n1 + n0) / (n1 * n0) + g * g / (2 * df)
    return {"effect": float(g), "se": float(math.sqrt(variance)), "n_R": n1, "n_NR": n0}


def random_effects(effects: pd.DataFrame) -> dict:
    x = effects.dropna(subset=["effect", "se"]).copy()
    k = len(x)
    if k == 0:
        return {"k": 0}
    y, v = x.effect.to_numpy(), x.se.to_numpy() ** 2
    fixed_w = 1 / v
    fixed = np.sum(fixed_w * y) / np.sum(fixed_w)
    q = np.sum(fixed_w * (y - fixed) ** 2)
    c = np.sum(fixed_w) - np.sum(fixed_w ** 2) / np.sum(fixed_w)
    tau2 = max(0.0, (q - (k - 1)) / c) if k > 1 and c > 0 else 0.0
    w = 1 / (v + tau2)
    mu = np.sum(w * y) / np.sum(w)
    se = math.sqrt(1 / np.sum(w))
    i2 = max(0.0, (q - (k - 1)) / q) if q > 0 and k > 1 else 0.0
    pred = math.sqrt(se * se + tau2)
    return {"k": k, "pooled_effect": mu, "pooled_se": se, "ci_low": mu - 1.96 * se, "ci_high": mu + 1.96 * se,
            "tau2": tau2, "Q": q, "I2": i2, "prediction_low": mu - 1.96 * pred, "prediction_high": mu + 1.96 * pred}


def bayesian_partial_pool(effects: pd.DataFrame, cfg: dict, mu_prior_sd: float = 2.0, tau_prior_sd: float = 1.0) -> dict:
    x = effects.dropna(subset=["effect", "se"])
    if len(x) == 0:
        return {"k": 0}
    y, v = x.effect.to_numpy(), x.se.to_numpy() ** 2
    taus = np.linspace(0, cfg["analysis"]["meta_tau_grid_max"], cfg["analysis"]["meta_tau_grid_n"])
    logp, means, variances = [], [], []
    for tau in taus:
        w = 1 / (v + tau * tau)
        prior_precision = 1 / (mu_prior_sd ** 2)
        post_var = 1 / (w.sum() + prior_precision)
        post_mean = post_var * np.sum(w * y)
        means.append(post_mean); variances.append(post_var)
        log_det = np.log(v + tau * tau).sum()
        quad = np.sum(w * (y - post_mean) ** 2) + post_mean * post_mean * prior_precision
        log_integral = 0.5 * math.log(post_var / (mu_prior_sd ** 2))
        log_prior = -0.5 * (tau / tau_prior_sd) ** 2 - math.log(tau_prior_sd)
        logp.append(-0.5 * (log_det + quad) + log_integral + log_prior)
    p = np.exp(np.asarray(logp) - np.max(logp)); p /= p.sum()
    means, variances = np.asarray(means), np.asarray(variances)
    mu = float(np.sum(p * means))
    var = float(np.sum(p * (variances + means ** 2)) - mu ** 2)
    sd = math.sqrt(max(var, 0))
    prob_pos = 0.5 * (1 + math.erf(mu / (sd * math.sqrt(2)))) if sd > 0 else float(mu > 0)
    return {"k": len(x), "posterior_mean": mu, "posterior_sd": sd, "credible_low": mu - 1.96 * sd,
            "credible_high": mu + 1.96 * sd, "probability_positive": prob_pos,
            "probability_negative": 1 - prob_pos, "posterior_tau_mean": float(np.sum(p * taus)),
            "mu_prior_sd": mu_prior_sd, "tau_prior_sd": tau_prior_sd}


def build_anchor_environment(cfg: dict) -> pd.DataFrame:
    barriers = pd.read_parquet(cfg["inputs"]["barrier_matrix"])
    barriers["timepoint"] = [
        normalize_timepoint(value, cohort)
        for value, cohort in zip(barriers["timepoint"], barriers["cohort_id"])
    ]
    patient = pd.read_csv(cfg["inputs"]["patient_metadata"])
    sample = pd.read_csv(cfg["inputs"]["sample_metadata"])
    cohorts = pd.read_csv(cfg["inputs"]["cohort_registry"])
    roles = pd.read_csv(cfg["inputs"]["dataset_roles"])

    pcols = ["patient_key", "cancer_type", "response_raw_summary", "response_endpoint_type", "response_binary_harmonized",
             "response_harmonization_confidence", "supervised_use_allowed", "support_use_allowed", "split"]
    env = barriers.merge(patient[pcols], on="patient_key", how="left", validate="many_to_one")
    sample["timepoint_normalized"] = [
        normalize_timepoint(value, cohort)
        for value, cohort in zip(sample["timepoint"], sample["cohort_id"])
    ]
    summary = sample.groupby(["patient_key", "timepoint_normalized"], dropna=False).agg(
        treatment_raw=("treatment_raw", lambda x: "|".join(sorted(set(x.dropna().astype(str))))),
        treatment_context=("treatment_context", lambda x: "|".join(sorted(set(x.dropna().astype(str))))),
        n_registered_samples=("sample_key", "nunique"),
    ).reset_index()
    env = env.merge(summary, left_on=["patient_key", "timepoint"], right_on=["patient_key", "timepoint_normalized"], how="left")
    env = env.merge(cohorts[["cohort_id", "anchor_lane_status", "dataset_role"]], on="cohort_id", how="left", validate="many_to_one")
    env = env.merge(roles[["cohort_id", "clean_anchor_lane", "support_lane"]], on="cohort_id", how="left", validate="many_to_one")
    response_env = pd.read_csv(cfg["inputs"]["response_environment"])[["cohort_id", "endpoint_type", "confidence", "evidence_source", "use_boundary"]]
    response_env = response_env.drop_duplicates("cohort_id").rename(columns={"endpoint_type": "environment_endpoint_type", "confidence": "environment_confidence",
        "evidence_source": "environment_evidence_source", "use_boundary": "environment_use_boundary"})
    env = env.merge(response_env, on="cohort_id", how="left", validate="many_to_one")

    g286 = pd.read_csv(cfg["inputs"]["gse286827_patient_metadata"], sep="\t")
    g286["g286_patient_id"] = g286.Patient.astype(str).str.extract(r"(P\d+)$", expand=False)
    env["g286_patient_id"] = env.patient_key.astype(str).str.extract(r"(P\d+)$", expand=False)
    env = env.merge(g286[["g286_patient_id", "Neoadj_type"]], on="g286_patient_id", how="left")
    text = (env.treatment_raw.fillna("") + "|" + env.treatment_context.fillna("")).str.lower()
    env["treatment_arm_class"] = np.where(text.str.contains(r"d\+t|dual|combo|ctla4|extension", regex=True), "combination",
                                           np.where(text.str.contains(r"pd1|durvalumab|pembrolizumab|nivolumab|monotherapy|anchor", regex=True), "monotherapy", "unknown"))
    mask286 = env.cohort_id.eq("GSE286827")
    env.loc[mask286, "treatment_arm_class"] = np.where(env.loc[mask286, "Neoadj_type"].eq("D"), "monotherapy_D_only",
                                                        np.where(env.loc[mask286, "Neoadj_type"].eq("D+T"), "combination_D_plus_T", "unknown"))
    env["response_known"] = env.response_binary_harmonized.isin(["responder", "non_responder"])
    env["response_binary"] = env.response_binary_harmonized.map({"responder": 1, "non_responder": 0})
    env["score_resolution"] = "mid"
    env["barrier_score_available"] = env[cfg["analysis"]["barriers"]].notna().all(axis=1)
    env["baseline_primary_eligible"] = env.timepoint.eq(cfg["analysis"]["primary_timepoint"])
    clean = env.clean_anchor_lane.astype(str).str.lower().isin({"yes", "conditional"})
    supervised = env.supervised_use_allowed.astype(str).str.lower().eq("yes")
    mono = env.treatment_arm_class.isin(["monotherapy", "monotherapy_D_only"])
    endpoint_ok = ~env.response_endpoint_type.fillna("unknown_or_mixed").isin(["unknown_or_mixed", "unknown"])
    endpoint_consistent = (env.environment_endpoint_type.isna() | env.response_endpoint_type.eq(env.environment_endpoint_type))
    environment_supported = env.environment_confidence.fillna("low").isin(["high", "medium"])
    env["phase8_analysis_role"] = "registry_or_support_only"
    env.loc[env.response_known & env.barrier_score_available, "phase8_analysis_role"] = "sensitivity"
    env.loc[env.response_known & clean & supervised & mono & endpoint_ok & endpoint_consistent & environment_supported & env.baseline_primary_eligible, "phase8_analysis_role"] = "conditional_anchor_primary"
    env["eligibility_reason"] = np.select(
        [~env.barrier_score_available, ~env.response_known, ~clean, ~supervised, ~mono, ~endpoint_ok, ~endpoint_consistent, ~environment_supported, ~env.baseline_primary_eligible],
        ["missing_barrier_score", "response_unknown", "not_clean_anchor_lane", "supervised_use_not_allowed", "non_monotherapy_or_unknown_arm", "unknown_or_mixed_endpoint",
         "endpoint_environment_mismatch", "response_environment_low_confidence", "nonbaseline_timepoint"],
        default="eligible_conditional_anchor",
    )
    env["input_version"] = "phase8_anchor_environment_v6_2_1"
    return env


def project_gse301741(cfg: dict, env: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    import anndata as ad
    membership = pd.read_csv(cfg["inputs"]["phase6_membership"])
    adata = ad.read_h5ad(cfg["inputs"]["gse301741_h5ad"], backed="r")
    genes = pd.Index(adata.var_names.astype(str))
    sample_ids = adata.obs.sample_id.astype(str)
    module_weights = {}
    selected = set()
    for module in cfg["analysis"]["barriers"]:
        m = membership[membership.frozen_module_id.eq(module)].set_index("gene").membership_weight
        available = [g for g in m.index if g in genes]
        pos = genes.get_indexer(available)
        w = m.loc[available].to_numpy(dtype=float); w /= w.sum()
        module_weights[module] = (pos, w); selected.update(pos.tolist())
    rows = []
    layer = adata.layers["counts"]
    for sid in sorted(sample_ids.unique()):
        idx = np.flatnonzero(sample_ids.to_numpy() == sid)
        library = 0.0
        sums = {m: 0.0 for m in cfg["analysis"]["barriers"]}
        for start in range(0, len(idx), 5000):
            block = layer[idx[start:start + 5000], :]
            library += float(block.sum())
            for module, (pos, w) in module_weights.items():
                sums[module] += float(np.asarray(block[:, pos].sum(axis=0)).ravel() @ w)
        obs = adata.obs.iloc[idx[0]]
        row = {"sample_key": f"GSE301741::{sid}", "sample_id": sid, "patient_key": str(obs.patient_key),
               "timepoint": normalize_timepoint(obs.timepoint), "tissue_context": str(obs.tissue_source), "n_cells": len(idx), "library_size": library}
        for module in cfg["analysis"]["barriers"]:
            row[module] = math.log1p(1e6 * sums[module] / library) if library > 0 else np.nan
        rows.append(row)
    adata.file.close()
    projection = pd.DataFrame(rows)
    for module in cfg["analysis"]["barriers"]:
        projection[f"{module}_projection_z"] = robust_z(projection[module], pd.Series("GSE301741", index=projection.index), min_n=10)
    projection["patient_timepoint_context_id"] = projection.patient_key + "::" + projection.timepoint + "::" + projection.tissue_context

    mid = env[env.cohort_id.eq("GSE301741")][["patient_timepoint_context_id"] + cfg["analysis"]["barriers"]]
    bridge = projection.merge(mid, on="patient_timepoint_context_id", suffixes=("_sample", "_mid"))
    metrics = []
    for module in cfg["analysis"]["barriers"]:
        metrics.append({"module_id": module, "n_shared_contexts": len(bridge),
                        "sample_mid_spearman": spearman(bridge[f"{module}_projection_z"], bridge[f"{module}_mid"]) if len(bridge) else np.nan})
    audit = pd.DataFrame(metrics)
    passed = len(bridge) >= cfg["analysis"]["projection_bridge_min_contexts"] and (audit.sample_mid_spearman >= cfg["analysis"]["projection_bridge_spearman_min"]).all()
    audit["bridge_status"] = "projection_support" if passed else "support_only_insufficient_or_unstable_bridge"
    projection["final_role"] = audit.bridge_status.iloc[0]
    return projection, audit


def patient_level(data: pd.DataFrame, barriers: list[str], baseline_only: bool) -> pd.DataFrame:
    x = data.copy()
    if baseline_only:
        x = x[x.timepoint.eq("baseline")]
    keys = ["cohort_id", "patient_key", "cancer_type", "response_endpoint_type", "treatment_arm_class", "response_binary"]
    return x.groupby(keys, dropna=False)[barriers].mean().reset_index()


def cross_fitted_residualization(cfg: dict, env: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    from run_phase7a_measurement import make_activity, prepare_units, read_phase6
    from phase7_common import load_config as load_phase7_config

    cfg7 = load_phase7_config()
    _, score, _, _ = read_phase6(cfg7)
    units = prepare_units(cfg7, score)
    activity = make_activity(cfg7, units)
    continuous = ["state_fraction", "log_library", "log_cells", "annotation_coverage", "low_quality_fraction", "mitochondrial_fraction"]
    categorical = ["cell_state", "cancer_type", "platform", "tissue_source", "timepoint", "treatment_context"]
    x_raw, _ = build_design(units, continuous, categorical)
    _, r = np.linalg.qr(x_raw, mode="reduced")
    keep = np.abs(np.diag(r)) > max(x_raw.shape) * np.finfo(float).eps * np.abs(np.diag(r)).max()
    x = x_raw[:, keep]
    targets = ["TASK01", "GSE286827", "GSE301741"]
    records = []
    for target in targets:
        train = units.cohort_id.ne(target).to_numpy()
        held = units.cohort_id.eq(target).to_numpy()
        weights = units.precision_weight.to_numpy(dtype=float)
        for module in cfg["analysis"]["barriers"]:
            y = activity.loc[activity.module_id.eq(module), "activity_global_robust_z"].to_numpy(dtype=float)
            valid = train & np.isfinite(y) & (weights > 0) & np.isfinite(x).all(axis=1)
            root_w = np.sqrt(weights[valid]); xv, yv = x[valid], y[valid]
            gram = (xv * root_w[:, None]).T @ (xv * root_w[:, None])
            penalty = np.eye(gram.shape[0]) * 1e-6; penalty[0, 0] = 0
            beta = np.linalg.solve(gram + penalty, (xv * root_w[:, None]).T @ (yv * root_w))
            idx = np.flatnonzero(held)
            tmp = units.iloc[idx][["cohort_id", "sample_key", "patient_key", "timepoint", "tissue_context", "patient_timepoint_context_id", "binding_status",
                                          "cell_state_level", "cell_state", "state_fraction", "eligible_fraction_denominator", "precision_weight"]].copy()
            tmp["barrier_id"] = module
            tmp["crossfit_residual"] = y[idx] - x[idx] @ beta
            records.append(tmp)
    long = pd.concat(records, ignore_index=True)
    long = long[long.binding_status.eq("bound_for_patient_context") & long.cell_state_level.eq("mid")]
    sample_keys = ["cohort_id", "sample_key", "patient_key", "timepoint", "tissue_context", "patient_timepoint_context_id", "barrier_id"]
    coverage = long.groupby(sample_keys, dropna=False).agg(represented=("state_fraction", "sum"), denominator=("eligible_fraction_denominator", "max"),
                                                           n_states=("cell_state", "nunique")).reset_index()
    coverage["coverage"] = coverage.represented / coverage.denominator.replace(0, np.nan)
    valid = coverage[(coverage.coverage >= 0.5) & (coverage.n_states >= 2)][sample_keys]
    long = long.merge(valid.assign(_valid=True), on=sample_keys, how="inner")
    long["_wx"] = long.crossfit_residual * long.precision_weight
    sample = long.groupby(sample_keys, dropna=False).agg(_wx=("_wx", "sum"), _w=("precision_weight", "sum")).reset_index()
    sample["score"] = sample._wx / sample._w
    context_keys = ["cohort_id", "patient_key", "timepoint", "tissue_context", "patient_timepoint_context_id", "barrier_id"]
    context = sample.groupby(context_keys, dropna=False).score.mean().reset_index()
    wide = context.pivot_table(index=context_keys[:-1], columns="barrier_id", values="score").reset_index()

    labels = env[["patient_key", "cancer_type", "response_endpoint_type", "treatment_arm_class", "response_binary", "phase8_analysis_role"]].drop_duplicates("patient_key")
    joined = wide.merge(labels, on="patient_key", how="left")
    primary = patient_level(joined[joined.phase8_analysis_role.eq("conditional_anchor_primary")], cfg["analysis"]["barriers"], baseline_only=True)
    assoc = association_table(primary, cfg, "outer_cohort_cross_fitted_residualization")
    assoc["crossfit_status"] = "complete_response_blind_outer_cohort_fit"
    return wide, assoc


def association_table(data: pd.DataFrame, cfg: dict, analysis_scope: str) -> pd.DataFrame:
    barriers = cfg["analysis"]["barriers"]
    rows = []
    group = ["cohort_id", "cancer_type", "response_endpoint_type", "treatment_arm_class"]
    for keys, frame in data.groupby(group, dropna=False):
        for barrier in barriers:
            est = hedges_g(frame, barrier)
            rows.append(dict(zip(group, keys)) | {"barrier_id": barrier, "analysis_scope": analysis_scope, **est,
                                                  "direction": "responder_higher" if est["effect"] > 0 else "nonresponder_higher" if est["effect"] < 0 else "undetermined"})
    return pd.DataFrame(rows)


def meta_tables(assoc: pd.DataFrame, cfg: dict) -> tuple[pd.DataFrame, pd.DataFrame]:
    usable = assoc[(assoc.n_R >= cfg["analysis"]["minimum_patients_per_class"]) & (assoc.n_NR >= cfg["analysis"]["minimum_patients_per_class"])]
    frequent, bayes = [], []
    for scope_cols, scope_name in [(["barrier_id", "response_endpoint_type"], "endpoint_specific"), (["barrier_id"], "cross_endpoint_sensitivity")]:
        for keys, frame in usable.groupby(scope_cols, dropna=False):
            keys = keys if isinstance(keys, tuple) else (keys,)
            base = dict(zip(scope_cols, keys)) | {"meta_scope": scope_name}
            estimate = random_effects(frame)
            z = abs(estimate.get("pooled_effect", np.nan) / estimate.get("pooled_se", np.nan))
            estimate["p_value"] = math.erfc(z / math.sqrt(2)) if np.isfinite(z) else np.nan
            frequent.append(base | estimate)
            bayes.append(base | bayesian_partial_pool(frame, cfg))
    frequent = pd.DataFrame(frequent)
    frequent["fdr_within_meta_scope"] = frequent.groupby("meta_scope")["p_value"].transform(lambda x: bh_fdr(x))
    return frequent, pd.DataFrame(bayes)


def hierarchical_prior_sensitivity(assoc: pd.DataFrame, cfg: dict) -> tuple[pd.DataFrame, pd.DataFrame]:
    usable = assoc[(assoc.n_R >= cfg["analysis"]["minimum_patients_per_class"]) & (assoc.n_NR >= cfg["analysis"]["minimum_patients_per_class"])].dropna(subset=["effect", "se"])
    rows, ident = [], []
    for barrier, frame in usable.groupby("barrier_id"):
        environment = frame[["cohort_id", "cancer_type", "response_endpoint_type"]].drop_duplicates()
        cancer_endpoint_confounded = environment.cancer_type.nunique() == environment.response_endpoint_type.nunique() == len(environment)
        status = "not_estimable_cancer_endpoint_cohort_confounded" if len(environment) < 3 or cancer_endpoint_confounded else "estimable"
        ident.append({"barrier_id": barrier, "n_environments": len(environment), "n_cancers": environment.cancer_type.nunique(),
                      "n_endpoints": environment.response_endpoint_type.nunique(), "shared_context_separation_status": status})
        for mu_sd in (0.5, 1.0, 2.0, 4.0):
            for tau_sd in (0.5, 1.0, 2.0):
                rows.append({"barrier_id": barrier, "model_status": status, **bayesian_partial_pool(frame, cfg, mu_sd, tau_sd)})
    return pd.DataFrame(rows), pd.DataFrame(ident)


def leave_out(assoc: pd.DataFrame, cfg: dict) -> pd.DataFrame:
    usable = assoc.dropna(subset=["effect", "se"])
    rows = []
    for barrier, data in usable.groupby("barrier_id"):
        for field in ("cohort_id", "cancer_type", "response_endpoint_type"):
            for value in data[field].dropna().unique():
                est = random_effects(data[data[field] != value])
                rows.append({"barrier_id": barrier, "leave_out_type": field, "left_out_value": value, **est})
    return pd.DataFrame(rows)


def stratified_meta_permutation(primary: pd.DataFrame, observed_assoc: pd.DataFrame, cfg: dict) -> pd.DataFrame:
    observed = {}
    usable = observed_assoc[(observed_assoc.n_R >= cfg["analysis"]["minimum_patients_per_class"]) & (observed_assoc.n_NR >= cfg["analysis"]["minimum_patients_per_class"])]
    for barrier in cfg["analysis"]["barriers"]:
        observed[barrier] = random_effects(usable[usable.barrier_id.eq(barrier)]).get("pooled_effect", np.nan)
    rng = np.random.default_rng(cfg["seed"] + 811)
    null = {b: [] for b in cfg["analysis"]["barriers"]}
    strata = ["cohort_id", "response_endpoint_type", "treatment_arm_class"]
    for _ in range(cfg["analysis"]["permutation_n"]):
        perm = primary.copy()
        perm["response_binary"] = perm.groupby(strata, dropna=False).response_binary.transform(lambda x: rng.permutation(x.to_numpy()))
        pa = association_table(perm, cfg, "stratified_permutation")
        pu = pa[(pa.n_R >= cfg["analysis"]["minimum_patients_per_class"]) & (pa.n_NR >= cfg["analysis"]["minimum_patients_per_class"])]
        for barrier in cfg["analysis"]["barriers"]:
            null[barrier].append(random_effects(pu[pu.barrier_id.eq(barrier)]).get("pooled_effect", np.nan))
    rows = []
    for barrier in cfg["analysis"]["barriers"]:
        values = np.asarray(null[barrier], dtype=float); values = values[np.isfinite(values)]
        p = (1 + np.sum(np.abs(values) >= abs(observed[barrier]))) / (len(values) + 1) if len(values) and np.isfinite(observed[barrier]) else np.nan
        rows.append({"barrier_id": barrier, "observed_random_effects_meta": observed[barrier], "permutation_n": len(values),
                     "empirical_two_sided_p": p, "exchangeability_block": "cohort_endpoint_treatment_arm_patient"})
    return pd.DataFrame(rows)


def negative_controls(primary: pd.DataFrame, cfg: dict) -> tuple[pd.DataFrame, pd.DataFrame]:
    rows = []
    y = primary.response_binary.to_numpy(dtype=float)
    for field in ("cohort_id", "response_endpoint_type", "treatment_arm_class"):
        group_sum = primary.groupby(field).response_binary.transform("sum") - primary.response_binary
        group_n = primary.groupby(field).response_binary.transform("count") - 1
        pred = group_sum / group_n.replace(0, np.nan)
        rows.append({"control_model": field + "_only", "n_patients": len(primary), "auc": auc_rank(y, pred.to_numpy())})
    for barrier in cfg["analysis"]["barriers"]:
        rows.append({"control_model": barrier + "_only", "n_patients": len(primary), "auc": auc_rank(y, primary[barrier].to_numpy())})
    fold_rows = []
    for cohort, frame in primary.groupby("cohort_id"):
        for barrier in cfg["analysis"]["barriers"]:
            raw = hedges_g(frame, barrier)
            ranked = frame.copy(); ranked[barrier] = ranked[barrier].rank(pct=True)
            rank_est = hedges_g(ranked, barrier)
            fold_rows.append({"cohort_id": cohort, "barrier_id": barrier, "raw_effect": raw["effect"], "fold_aware_rank_effect": rank_est["effect"],
                              "sign_agreement": np.sign(raw["effect"]) == np.sign(rank_est["effect"]) if np.isfinite(raw["effect"]) and np.isfinite(rank_est["effect"]) else False,
                              "method": "response_blind_within_cohort_rank_recalibration"})
    return pd.DataFrame(rows), pd.DataFrame(fold_rows)


def gse206325_sensitivity(cfg: dict) -> pd.DataFrame:
    coarse = pd.read_parquet(cfg["inputs"]["coarse_matrix"])
    binding = pd.read_csv(cfg["inputs"]["response_binding"])
    binding = binding[binding.cohort_id.eq("GSE206325")].groupby("patient_key", dropna=False).agg(
        response_binary_harmonized=("response_binary_harmonized", lambda x: x.dropna().iloc[0] if len(x.dropna()) else "unknown"),
        response_endpoint_type=("response_endpoint_type", lambda x: x.dropna().iloc[0] if len(x.dropna()) else "unknown_or_mixed"),
    ).reset_index()
    x = coarse[coarse.cohort_id.eq("GSE206325")].merge(binding, on="patient_key", how="left")
    x["response_binary"] = x.response_binary_harmonized.map({"responder": 1, "non_responder": 0})
    rows = []
    for scope, frame in [("baseline_only_unavailable", x[x.timepoint.eq("baseline")]), ("post_treatment_patient_mean", x[x.timepoint.eq("post_treatment")].groupby(["patient_key", "response_binary", "response_endpoint_type"], dropna=False)[cfg["analysis"]["barriers"]].mean().reset_index())]:
        for barrier in cfg["analysis"]["barriers"]:
            est = hedges_g(frame, barrier)
            rows.append({"cohort_id": "GSE206325", "barrier_id": barrier, "analysis_scope": scope, "score_resolution": "coarse",
                         "endpoint_type": "unknown_or_mixed", **est, "use_boundary": "post_treatment_coarse_HCC_sensitivity_only_not_predictive_not_pooled_with_mid"})
    return pd.DataFrame(rows)


def hcc_context_adjudication(assoc: pd.DataFrame, hcc: pd.DataFrame, cfg: dict) -> pd.DataFrame:
    task = assoc[(assoc.cohort_id.eq("TASK01")) & assoc.analysis_scope.eq("baseline_conditional_anchor")]
    coarse = hcc[hcc.analysis_scope.eq("post_treatment_patient_mean")]
    rows = []
    for barrier in cfg["analysis"]["barriers"]:
        t = task[task.barrier_id.eq(barrier)]; c = coarse[coarse.barrier_id.eq(barrier)]
        te = float(t.effect.iloc[0]) if len(t) else np.nan; ce = float(c.effect.iloc[0]) if len(c) else np.nan
        rows.append({"barrier_id": barrier, "task01_mid_baseline_effect": te, "gse206325_coarse_post_effect": ce,
                     "direction_concordant": bool(np.sign(te) == np.sign(ce)) if np.isfinite(te) and np.isfinite(ce) else False,
                     "hcc_residual_estimable": False,
                     "estimability_status": "not_estimable_resolution_timepoint_endpoint_confounded",
                     "r2_use_boundary": "directional_sensitivity_only_cannot_define_HCC_residual"})
    return pd.DataFrame(rows)


def route_flags(statistical_pass: bool, shared_estimable: bool, hcc_estimable: bool, fold_pass: bool,
                same_endpoint_replicated: bool, endpoint_direction_consistent: bool) -> dict:
    return {
        "R1_pass": bool(statistical_pass and shared_estimable and fold_pass and same_endpoint_replicated),
        "R1_lite_pass": bool(statistical_pass and shared_estimable and fold_pass and endpoint_direction_consistent and not same_endpoint_replicated),
        "R2_pass": bool(statistical_pass and hcc_estimable and fold_pass),
        "route_adjudication_complete": bool(shared_estimable and hcc_estimable and fold_pass),
    }


def same_endpoint_same_direction_replicated(effects: pd.DataFrame, minimum_cohorts: int = 2, minimum_patients_per_class: int = 3) -> bool:
    valid = effects.dropna(subset=["effect"])
    if {"n_R", "n_NR"}.issubset(valid.columns):
        valid = valid[(valid.n_R >= minimum_patients_per_class) & (valid.n_NR >= minimum_patients_per_class)]
    for _, frame in valid.groupby("response_endpoint_type", dropna=False):
        cohort_effects = frame.groupby("cohort_id").effect.mean()
        signs = np.sign(cohort_effects[cohort_effects != 0])
        if len(signs) >= minimum_cohorts and signs.nunique() == 1:
            return True
    return False


def phase9_gate(meta: pd.DataFrame, prior: pd.DataFrame, ident: pd.DataFrame, hcc_context: pd.DataFrame, loco: pd.DataFrame,
                assoc: pd.DataFrame, controls: pd.DataFrame, perms: pd.DataFrame, fold: pd.DataFrame, cfg: dict) -> tuple[pd.DataFrame, dict]:
    rows = []
    usable = assoc.dropna(subset=["effect"])
    usable = usable[(usable.n_R >= cfg["analysis"]["minimum_patients_per_class"]) & (usable.n_NR >= cfg["analysis"]["minimum_patients_per_class"])]
    metadata_values = controls[controls.control_model.str.endswith("_only") & ~controls.control_model.str.startswith(tuple(cfg["analysis"]["barriers"]))].auc
    metadata_strength = float(np.max(np.abs(metadata_values - 0.5))) if len(metadata_values) else np.nan
    for barrier in cfg["analysis"]["barriers"]:
        a = usable[usable.barrier_id.eq(barrier)]
        signs = np.sign(a.effect)
        direction_count = max(int((signs > 0).sum()), int((signs < 0).sum())) if len(a) else 0
        b = prior[prior.barrier_id.eq(barrier)]
        probability = float(np.min(np.maximum(b.probability_positive, b.probability_negative))) if len(b) else np.nan
        l = loco[loco.barrier_id.eq(barrier)].dropna(subset=["pooled_effect"])
        overall = meta[(meta.barrier_id.eq(barrier)) & meta.meta_scope.eq("cross_endpoint_sensitivity")]
        overall_sign = np.sign(overall.pooled_effect.iloc[0]) if len(overall) else 0
        loo_stability = float((np.sign(l.pooled_effect) == overall_sign).mean()) if len(l) and overall_sign != 0 else 0.0
        barrier_auc_row = controls[controls.control_model.eq(barrier + "_only")]
        barrier_auc = float(barrier_auc_row.auc.iloc[0]) if len(barrier_auc_row) else np.nan
        auc_advantage = abs(barrier_auc - 0.5) - metadata_strength if np.isfinite(barrier_auc) and np.isfinite(metadata_strength) else np.nan
        perm = perms[perms.barrier_id.eq(barrier)]
        permutation_p = float(perm.empirical_two_sided_p.iloc[0]) if len(perm) else np.nan
        shared_status = ident.loc[ident.barrier_id.eq(barrier), "shared_context_separation_status"]
        shared_estimable = len(shared_status) and shared_status.iloc[0] == "estimable"
        hcc_status = hcc_context.loc[hcc_context.barrier_id.eq(barrier), "hcc_residual_estimable"]
        hcc_estimable = bool(hcc_status.iloc[0]) if len(hcc_status) else False
        f = fold[(fold.barrier_id.eq(barrier)) & (fold.n_R >= cfg["analysis"]["minimum_patients_per_class"]) & (fold.n_NR >= cfg["analysis"]["minimum_patients_per_class"])].dropna(subset=["effect"])
        original_by_cohort = a.set_index("cohort_id").effect
        comparisons = [np.sign(row.effect) == np.sign(original_by_cohort.get(row.cohort_id, np.nan)) for row in f.itertuples() if row.cohort_id in original_by_cohort]
        fold_sign_agreement = float(np.mean(comparisons)) if comparisons else np.nan
        fold_complete = f.cohort_id.nunique() >= 2 and np.isfinite(fold_sign_agreement) and fold_sign_agreement >= cfg["gate"]["leave_one_out_sign_stability_min"]
        statistical_pass = (direction_count >= cfg["gate"]["minimum_independent_cohorts_same_direction"]
                  and probability >= cfg["gate"]["direction_probability_min"]
                  and loo_stability >= cfg["gate"]["leave_one_out_sign_stability_min"]
                  and auc_advantage >= cfg["gate"]["metadata_negative_control_margin_min"]
                  and permutation_p <= 0.05)
        same_endpoint_replicated = same_endpoint_same_direction_replicated(a, cfg["gate"]["minimum_independent_cohorts_same_direction"], cfg["analysis"]["minimum_patients_per_class"])
        endpoint_effects = a.groupby("response_endpoint_type").effect.mean().dropna()
        endpoint_direction_consistent = len(endpoint_effects) >= 2 and np.abs(np.sign(endpoint_effects)).sum() == len(endpoint_effects) and np.sign(endpoint_effects).nunique() == 1
        routes = route_flags(statistical_pass, shared_estimable, hcc_estimable, fold_complete, same_endpoint_replicated, endpoint_direction_consistent)
        rows.append({"barrier_id": barrier, "n_estimable_strata": len(a), "max_same_direction_strata": direction_count,
                     "posterior_direction_probability": probability, "leave_out_sign_stability": loo_stability,
                     "barrier_auc": barrier_auc, "maximum_metadata_auc_deviation": metadata_strength, "auc_advantage_over_metadata": auc_advantage,
                     "shuffled_response_empirical_p": permutation_p,
                     "shared_context_separation_estimable": shared_estimable, "hcc_residual_estimable": hcc_estimable,
                     "fold_aware_residualization_complete": fold_complete, "fold_aware_sign_agreement": fold_sign_agreement,
                     "same_endpoint_replicated": same_endpoint_replicated, "endpoint_direction_consistent": endpoint_direction_consistent,
                     "statistical_screen_pass": statistical_pass, **routes,
                     "phase9_gate_pass": routes["R1_pass"] or routes["R1_lite_pass"] or routes["R2_pass"],
                     "failure_reason": "" if (routes["R1_pass"] or routes["R1_lite_pass"] or routes["R2_pass"]) else "phase8_incomplete_context_or_fold_evidence"})
    table = pd.DataFrame(rows)
    passed = table[table.phase9_gate_pass].barrier_id.tolist()
    phase8_complete = table.route_adjudication_complete.all()
    blockers = []
    if not table.shared_context_separation_estimable.all():
        blockers.append("shared_context_effect_not_identifiable")
    if not table.hcc_residual_estimable.all():
        blockers.append("HCC_R2_route_not_identifiable")
    if not table.fold_aware_residualization_complete.all():
        blockers.append("fold_aware_residualization_not_completed")
    verdict = "GO_TO_PHASE9" if passed else "STOP_AT_RESULT_LAYER_A" if phase8_complete else "PHASE8_INCOMPLETE_PHASE9_BLOCKED"
    manifest = {"phase": "phase8_anchor_context_adjudication", "verdict": verdict, "eligible_barriers_for_phase9": passed,
                "hard_blockers": blockers,
                "phase9_allowed": bool(passed), "srb_training_allowed": bool(passed), "counterfactual_repair_allowed": False,
                "route_verdict": "R1_or_R2_candidate" if passed else "Layer_A_only" if phase8_complete else "ROUTE_NOT_ADJUDICATED",
                "fallback_if_no_barrier_passes": "response_blind_pan_cancer_state_structure_and_measurement_framework"}
    return table, manifest


def run() -> dict:
    cfg = load_config(); out = cfg["out"]
    env = build_anchor_environment(cfg)
    env.to_csv(out / "preflight/anchor_environment_v6_2_1.csv", index=False)
    summary = env.groupby(["cohort_id", "response_endpoint_type", "treatment_arm_class", "score_resolution", "phase8_analysis_role"], dropna=False).agg(
        contexts=("patient_timepoint_context_id", "nunique"), patients=("patient_key", "nunique"),
        n_R=("response_binary", lambda x: int((x == 1).sum())), n_NR=("response_binary", lambda x: int((x == 0).sum())), n_unknown=("response_known", lambda x: int((~x).sum()))).reset_index()
    summary.to_csv(out / "preflight/anchor_surface_summary.csv", index=False)

    projection, bridge = project_gse301741(cfg, env)
    projection.to_csv(out / "projection/gse301741_sample_level_frozen_module_projection.csv", index=False)
    bridge.to_csv(out / "projection/gse301741_projection_bridge_audit.csv", index=False)
    hcc = gse206325_sensitivity(cfg); hcc.to_csv(out / "hcc_sensitivity/gse206325_coarse_response_sensitivity.csv", index=False)

    primary_context = env[env.phase8_analysis_role.eq("conditional_anchor_primary")].copy()
    primary = patient_level(primary_context, cfg["analysis"]["barriers"], baseline_only=True)
    all_known = patient_level(env[env.response_known & env.treatment_arm_class.isin(["monotherapy", "monotherapy_D_only"])], cfg["analysis"]["barriers"], baseline_only=False)
    primary.to_csv(out / "association/phase8_primary_patient_table.csv", index=False)
    assoc_primary = association_table(primary, cfg, "baseline_conditional_anchor")
    assoc_sens = association_table(all_known, cfg, "all_timepoint_monotherapy_sensitivity")
    assoc = pd.concat([assoc_primary, assoc_sens], ignore_index=True)
    assoc.to_csv(out / "association/barrier_response_association_by_cohort.csv", index=False)

    meta, bayes = meta_tables(assoc_primary, cfg)
    meta.to_csv(out / "meta/barrier_response_meta_analysis.csv", index=False)
    bayes.to_csv(out / "meta/barrier_response_bayesian_partial_pooling.csv", index=False)
    prior, ident = hierarchical_prior_sensitivity(assoc_primary, cfg)
    prior.to_csv(out / "meta/hierarchical_prior_sensitivity.csv", index=False)
    ident.to_csv(out / "meta/shared_context_identifiability.csv", index=False)
    loco = leave_out(assoc_primary, cfg); loco.to_csv(out / "robustness/leave_one_environment_analysis.csv", index=False)
    controls, rank_recalibration = negative_controls(primary, cfg)
    perms = stratified_meta_permutation(primary, assoc_primary, cfg)
    controls.to_csv(out / "negative_controls/metadata_only_negative_controls.csv", index=False)
    perms.to_csv(out / "negative_controls/stratified_shuffled_response_meta_audit.csv", index=False)
    rank_recalibration.to_csv(out / "robustness/response_blind_rank_recalibration_sensitivity.csv", index=False)
    crossfit_matrix, fold = cross_fitted_residualization(cfg, env)
    crossfit_matrix.to_parquet(out / "robustness/outer_cohort_crossfit_barrier_matrix.parquet", index=False)
    fold.to_csv(out / "robustness/fold_aware_residualization_sensitivity.csv", index=False)
    coverage = env.groupby("cohort_id", dropna=False).agg(total_contexts=("patient_timepoint_context_id", "nunique"), total_patients=("patient_key", "nunique"),
        response_known_contexts=("response_known", "sum"), primary_contexts=("phase8_analysis_role", lambda x: int((x == "conditional_anchor_primary").sum())),
        primary_patients=("patient_key", lambda x: x[env.loc[x.index, "phase8_analysis_role"].eq("conditional_anchor_primary")].nunique())).reset_index()
    coverage["response_known_fraction"] = coverage.response_known_contexts / coverage.total_contexts
    coverage["primary_context_fraction"] = coverage.primary_contexts / coverage.total_contexts
    coverage.to_csv(out / "negative_controls/coverage_selection_audit.csv", index=False)

    hcc_context = hcc_context_adjudication(assoc_primary, hcc, cfg)
    hcc_context.to_csv(out / "hcc_sensitivity/hcc_index_context_residuals_v6_2_1.csv", index=False)
    gate, manifest = phase9_gate(meta, prior, ident, hcc_context, loco, assoc_primary, controls, perms, fold, cfg)
    gate.to_csv(out / "handoff/phase9_entry_gate_by_barrier.csv", index=False)
    shared = gate[["barrier_id", "statistical_screen_pass", "shared_context_separation_estimable", "R1_pass", "R1_lite_pass"]].copy()
    shared.to_csv(out / "handoff/shared_antipd1_barrier_candidates.csv", index=False)
    context = gate[["barrier_id", "hcc_residual_estimable", "R2_pass"]].merge(hcc_context, on="barrier_id", how="left")
    context.to_csv(out / "handoff/context_modulated_barrier_candidates.csv", index=False)
    direction = meta[meta.meta_scope.eq("cross_endpoint_sensitivity")][["barrier_id", "pooled_effect", "ci_low", "ci_high"]].copy()
    direction["working_direction"] = np.where(direction.pooled_effect > 0, "responder_compatible", "failure_compatible")
    direction["use_boundary"] = "descriptive_cross_endpoint_sensitivity_not_route_qualified"
    direction.to_csv(out / "handoff/responder_compatible_direction_v0.csv", index=False)
    manifest["inputs"] = {k: str(v.relative_to(ROOT)) for k, v in cfg["inputs"].items()}
    manifest["outputs"] = {"anchor_environment": str((out / "preflight/anchor_environment_v6_2_1.csv").relative_to(ROOT)),
                           "association": str((out / "association/barrier_response_association_by_cohort.csv").relative_to(ROOT)),
                           "meta_analysis": str((out / "meta/barrier_response_meta_analysis.csv").relative_to(ROOT)),
                           "phase9_gate": str((out / "handoff/phase9_entry_gate_by_barrier.csv").relative_to(ROOT))}
    manifest["gse301741_projection_status"] = bridge.bridge_status.iloc[0]
    manifest["gse206325_role"] = "coarse_HCC_sensitivity_only"
    write_yaml(manifest, out / "handoff/phase8_to_phase9_handoff.yaml")
    write_yaml({"phase": "phase8_route_adjudication", "verdict": manifest["route_verdict"], "phase8_complete": manifest["verdict"] != "PHASE8_INCOMPLETE_PHASE9_BLOCKED",
                "R1_candidates": gate.loc[gate.R1_pass, "barrier_id"].tolist(), "R1_lite_candidates": gate.loc[gate.R1_lite_pass, "barrier_id"].tolist(),
                "R2_candidates": gate.loc[gate.R2_pass, "barrier_id"].tolist(), "blocking_reasons": manifest["hard_blockers"]}, out / "handoff/route_verdict_v6_2_1.yaml")

    gate_disk = pd.read_csv(out / "handoff/phase9_entry_gate_by_barrier.csv")
    handoff_disk = yaml.safe_load((out / "handoff/phase8_to_phase9_handoff.yaml").read_text())
    route_disk = yaml.safe_load((out / "handoff/route_verdict_v6_2_1.yaml").read_text())

    feature = pd.read_parquet(cfg["inputs"]["barrier_matrix"])
    g301_patients = projection.patient_key.nunique()
    hcc_all = hcc[hcc.analysis_scope.eq("post_treatment_patient_mean")]
    primary_g286 = env[(env.cohort_id == "GSE286827") & (env.phase8_analysis_role == "conditional_anchor_primary")]
    gate_recomputed = gate.phase9_gate_pass.eq(gate.R1_pass | gate.R1_lite_pass | gate.R2_pass).all()
    expected_verdict = "GO_TO_PHASE9" if gate_disk.phase9_gate_pass.any() else "STOP_AT_RESULT_LAYER_A" if gate_disk.route_adjudication_complete.all() else "PHASE8_INCOMPLETE_PHASE9_BLOCKED"
    expected_blockers = []
    if not gate_disk.shared_context_separation_estimable.all(): expected_blockers.append("shared_context_effect_not_identifiable")
    if not gate_disk.hcc_residual_estimable.all(): expected_blockers.append("HCC_R2_route_not_identifiable")
    if not gate_disk.fold_aware_residualization_complete.all(): expected_blockers.append("fold_aware_residualization_not_completed")
    estimable_primary = assoc_primary[(assoc_primary.n_R >= cfg["analysis"]["minimum_patients_per_class"]) & (assoc_primary.n_NR >= cfg["analysis"]["minimum_patients_per_class"])].cohort_id.nunique()
    reliability = f"""# Anchor Reliability Report v6.2.1

- Scored baseline patients: {len(primary)}
- Estimable primary cohorts: {estimable_primary}
- GSE301741 bridge: {manifest['gse301741_projection_status']}
- GSE206325 boundary: post-treatment coarse HCC sensitivity only
- Shared effect separation: not estimable because cohort, cancer and endpoint are confounded
- HCC residual: not estimable because TASK01 is baseline-mid while GSE206325 is post-treatment-coarse
- Outer-cohort response-blind crossfit: completed
- Route verdict: {manifest['route_verdict']}
"""
    (out / "anchor_reliability_report_v6_2_1.md").write_text(reliability)
    audit_rows = [
        {"check": "phase7_barriers_only", "passed": set(cfg["analysis"]["barriers"]) == {"FM01", "FM04", "FM07"}},
        {"check": "environment_left_join_complete", "passed": len(env) == len(feature) == 1036},
        {"check": "environment_context_key_unique", "passed": env.patient_timepoint_context_id.is_unique},
        {"check": "no_response_in_barrier_matrix", "passed": not any("response" in c.lower() for c in feature.columns)},
        {"check": "gse286827_D_only_primary", "passed": len(primary_g286) == 13 and not primary_g286.Neoadj_type.ne("D").any()},
        {"check": "gse243013_excluded_from_primary", "passed": not env[(env.cohort_id == "GSE243013") & (env.phase8_analysis_role == "conditional_anchor_primary")].shape[0]},
        {"check": "role_reason_lane_consistent", "passed": env.loc[env.phase8_analysis_role.eq("conditional_anchor_primary"), "eligibility_reason"].eq("eligible_conditional_anchor").all() and env.loc[env.clean_anchor_lane.astype(str).str.lower().eq("no"), "eligibility_reason"].isin(["not_clean_anchor_lane", "response_unknown", "missing_barrier_score"]).all()},
        {"check": "response_environment_consumed", "passed": env.loc[env.phase8_analysis_role.eq("conditional_anchor_primary"), "environment_confidence"].isin(["high", "medium"]).all() and env.loc[env.phase8_analysis_role.eq("conditional_anchor_primary"), "response_endpoint_type"].eq(env.loc[env.phase8_analysis_role.eq("conditional_anchor_primary"), "environment_endpoint_type"]).all()},
        {"check": "gse206325_not_pooled", "passed": "GSE206325" not in set(primary.cohort_id)},
        {"check": "gse206325_expected_response_surface", "passed": len(hcc_all) == 3 and hcc_all.n_R.eq(7).all() and hcc_all.n_NR.eq(17).all()},
        {"check": "gse301741_projection_complete", "passed": len(projection) == 27 and g301_patients == 16 and projection[cfg["analysis"]["barriers"]].notna().all().all()},
        {"check": "projection_not_primary", "passed": len(bridge) == 3 and bridge.n_shared_contexts.eq(4).all() and not bridge.bridge_status.eq("projection_support").any()},
        {"check": "primary_surface_expected", "passed": len(primary) == 39 and int((primary.response_binary == 1).sum()) == 17 and int((primary.response_binary == 0).sum()) == 22},
        {"check": "meta_and_bayesian_outputs_complete", "passed": set(meta.barrier_id) == set(cfg["analysis"]["barriers"]) and set(bayes.barrier_id) == set(cfg["analysis"]["barriers"])},
        {"check": "hierarchical_prior_sensitivity_complete", "passed": len(prior) == 36 and set(ident.barrier_id) == set(cfg["analysis"]["barriers"])},
        {"check": "shared_context_identifiability_explicit", "passed": ident.shared_context_separation_status.notna().all()},
        {"check": "hcc_context_branch_completed", "passed": len(hcc_context) == 3 and not hcc_context.hcc_residual_estimable.any()},
        {"check": "stratified_permutation_complete", "passed": len(perms) == 3 and perms.permutation_n.eq(cfg["analysis"]["permutation_n"]).all()},
        {"check": "fold_aware_crossfit_completed", "passed": gate.fold_aware_residualization_complete.all() and len(fold) >= 6},
        {"check": "required_route_outputs_exist", "passed": all((out / f).exists() for f in ["hcc_sensitivity/hcc_index_context_residuals_v6_2_1.csv", "handoff/shared_antipd1_barrier_candidates.csv", "handoff/context_modulated_barrier_candidates.csv", "handoff/responder_compatible_direction_v0.csv", "handoff/route_verdict_v6_2_1.yaml", "robustness/outer_cohort_crossfit_barrier_matrix.parquet", "anchor_reliability_report_v6_2_1.md"])},
        {"check": "negative_controls_complete", "passed": set(cfg["analysis"]["barriers"]).issubset(set(controls.control_model.str.replace("_only", "", regex=False))) and len(perms) == 3},
        {"check": "phase9_gate_recomputed_consistent", "passed": gate_recomputed},
        {"check": "handoff_verdict_consistent", "passed": handoff_disk["verdict"] == expected_verdict and handoff_disk["hard_blockers"] == expected_blockers and route_disk["blocking_reasons"] == expected_blockers},
        {"check": "phase9_gate_machine_readable", "passed": len(gate) == 3},
    ]
    audit = pd.DataFrame(audit_rows); audit.to_csv(out / "audit/phase8_strict_audit.csv", index=False)
    audit_verdict = "PASS" if audit.passed.all() else "FAIL"
    write_yaml({"phase": "phase8_strict_audit", "verdict": audit_verdict, "n_checks": len(audit), "n_failures": int((~audit.passed).sum())}, out / "audit/phase8_strict_audit_manifest.yaml")
    if audit_verdict != "PASS":
        raise RuntimeError("Phase8 strict audit failed")

    report = f"""# Phase8 Anchor and Context Adjudication Report

## Verdict

**{manifest['verdict']}**

Phase8 先冻结真实 anchor surface，再对 FM01、FM04、FM07 进行 cohort/endpoint/treatment-stratified association。未训练 SRB，未执行 counterfactual repair。

## Actual Surface

- Barrier contexts: {len(env):,}
- Primary baseline patients: {len(primary):,}
- Registered/scored primary cohorts: {primary.cohort_id.nunique()}
- Estimable primary cohorts: {estimable_primary}
- GSE301741 projection: {manifest['gse301741_projection_status']}
- GSE206325: coarse HCC sensitivity only

## Gate

进入 Phase9 的 barriers: {', '.join(manifest['eligible_barriers_for_phase9']) if manifest['eligible_barriers_for_phase9'] else 'none'}。

Outer-cohort response-blind crossfit sensitivity 已完成。当前 shared/context separation 与 HCC residual 仍不可辨识，因此 Phase8 未完成、Phase9 blocked。不得训练 supervised SRB，也不得进入 X-class repair；此状态尚不能最终排除 R2 或裁定 Layer-A-only。
"""
    (out / "PHASE8_ANCHOR_CONTEXT_ADJUDICATION_REPORT.md").write_text(report)
    return manifest


if __name__ == "__main__":
    parser = argparse.ArgumentParser(); parser.parse_args()
    print(yaml.safe_dump(run(), sort_keys=False, allow_unicode=True))
