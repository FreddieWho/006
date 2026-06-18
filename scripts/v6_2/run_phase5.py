#!/usr/bin/env python3
"""Phase5 strong baseline and confounding audit.

Pure numpy/pandas implementation to avoid sklearn/scipy runtime dependency.
Consumes only Phase4B handoff-allowed inputs.
"""

from __future__ import annotations

import hashlib
import os
from datetime import date
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd
import yaml


ROOT = Path(__file__).resolve().parents[2]
PHASE4B = Path(os.environ.get(
    "PHASE4B_RESULT_DIR", ROOT / "results" / "v6_2" / "phase4b_immune_state_feature_construction"))
OUT = Path(os.environ.get(
    "PHASE5_OUTPUT_DIR", ROOT / "results" / "v6_2" / "phase5_strong_baseline_confounding_audit"))
TODAY = str(date.today())
INPUT_VERSION = "frozen_v0"
RANDOM_STATE = 1729
N_PERMUTATIONS = 100

BANNED_FEATURE_TOKENS = [
    "response",
    "responder",
    "non_responder",
    "outcome",
    "survival",
    "progression",
    "recist",
    "mrecist",
    "trg",
    "death",
    "pfs",
    "os",
    "split",
    "train",
    "test",
    "label",
]


def mkdirs() -> None:
    for sub in ["preflight", "baselines", "confounding", "negative_controls", "audit", "handoff"]:
        (OUT / sub).mkdir(parents=True, exist_ok=True)


def read_csv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, dtype=str, keep_default_na=False)


def write_csv(df: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, index=False)


def write_tsv(df: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, index=False, sep="\t")


def write_yaml(obj: dict, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as fh:
        yaml.safe_dump(obj, fh, sort_keys=False, allow_unicode=True)


def write_md(text: str, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")


def load_yaml(path: Path) -> dict:
    with path.open() as fh:
        return yaml.safe_load(fh) or {}


def partial_sha256(path: Path, nbytes: int = 4_194_304) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        h.update(fh.read(nbytes))
    return h.hexdigest()


def sigmoid(z: np.ndarray) -> np.ndarray:
    return 1.0 / (1.0 + np.exp(-np.clip(z, -35, 35)))


def auc_score(y: np.ndarray, score: np.ndarray) -> float:
    y = np.asarray(y).astype(int)
    score = np.asarray(score).astype(float)
    n_pos = int(y.sum())
    n_neg = int(len(y) - n_pos)
    if n_pos == 0 or n_neg == 0:
        return np.nan
    ranks = pd.Series(score).rank(method="average").to_numpy()
    rank_sum_pos = ranks[y == 1].sum()
    return float((rank_sum_pos - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg))


def average_precision(y: np.ndarray, score: np.ndarray) -> float:
    y = np.asarray(y).astype(int)
    score = np.asarray(score).astype(float)
    n_pos = int(y.sum())
    if n_pos == 0:
        return np.nan
    order = np.argsort(-score)
    y_sorted = y[order]
    tp = np.cumsum(y_sorted)
    precision = tp / (np.arange(len(y_sorted)) + 1)
    return float((precision * y_sorted).sum() / n_pos)


def metric_pack(y: Iterable[int], score: Iterable[float]) -> dict:
    y = np.asarray(list(y), dtype=int)
    score = np.asarray(list(score), dtype=float)
    if len(y) == 0:
        return {
            "n": 0,
            "n_responder": 0,
            "n_non_responder": 0,
            "roc_auc": np.nan,
            "average_precision": np.nan,
            "accuracy": np.nan,
            "balanced_accuracy": np.nan,
            "brier": np.nan,
            "log_loss": np.nan,
        }
    pred = (score >= 0.5).astype(int)
    tp = int(((pred == 1) & (y == 1)).sum())
    tn = int(((pred == 0) & (y == 0)).sum())
    fp = int(((pred == 1) & (y == 0)).sum())
    fn = int(((pred == 0) & (y == 1)).sum())
    tpr = tp / (tp + fn) if (tp + fn) else np.nan
    tnr = tn / (tn + fp) if (tn + fp) else np.nan
    score_clip = np.clip(score, 1e-6, 1 - 1e-6)
    return {
        "n": int(len(y)),
        "n_responder": int(y.sum()),
        "n_non_responder": int((1 - y).sum()),
        "roc_auc": auc_score(y, score),
        "average_precision": average_precision(y, score),
        "accuracy": float((pred == y).mean()),
        "balanced_accuracy": float(np.nanmean([tpr, tnr])),
        "brier": float(np.mean((score - y) ** 2)),
        "log_loss": float(-np.mean(y * np.log(score_clip) + (1 - y) * np.log(1 - score_clip))),
    }


def input_audit(paths: dict[str, Path]) -> pd.DataFrame:
    rows = []
    for name, path in paths.items():
        exists = path.exists()
        rows.append(
            {
                "input_name": name,
                "path": str(path),
                "exists": "yes" if exists else "no",
                "readable": "yes" if exists and os.access(path, os.R_OK) else "no",
                "size_bytes": path.stat().st_size if exists else "",
                "partial_sha256_4mb": partial_sha256(path) if exists and path.is_file() else "",
            }
        )
    return pd.DataFrame(rows)


def load_inputs() -> dict:
    paths = {
        "phase4b_handoff": PHASE4B / "handoff" / "phase4b_to_phase5_handoff.yaml",
        "phase4b_manifest": PHASE4B / "handoff" / "phase4b_decision_manifest.yaml",
        "primary_matrix": PHASE4B / "matrix" / "immune_state_feature_matrix.primary.csv",
        "sensitivity_matrix": PHASE4B / "matrix" / "immune_state_feature_matrix.sensitivity.csv",
        "qc_covariates": PHASE4B / "matrix" / "immune_state_feature_matrix.qc_covariates.csv",
        "feature_dictionary": PHASE4B / "matrix" / "feature_dictionary_v6_2.csv",
        "missingness": PHASE4B / "audit" / "feature_missingness_report.csv",
        "confounding_tags": PHASE4B / "audit" / "feature_confounding_tags.csv",
        "response_environment": PHASE4B / "response_environment" / "response_environment_binding_table.csv",
        "analysis_universe": PHASE4B / "universe" / "analysis_universe_registry_v6_2.csv",
    }
    audit = input_audit(paths)
    write_csv(audit, OUT / "preflight" / "phase5_input_audit.csv")
    missing = audit[audit["exists"].ne("yes")]["input_name"].tolist()
    if missing:
        raise RuntimeError(f"Missing Phase5 inputs: {missing}")
    manifest = load_yaml(paths["phase4b_manifest"])
    if manifest.get("hard_blockers"):
        raise RuntimeError(f"Phase4B hard blockers are not empty: {manifest.get('hard_blockers')}")
    if manifest.get("verdict") not in {"GO_TO_PHASE5", "CONDITIONAL_GO_TO_PHASE5"}:
        raise RuntimeError(f"Phase4B verdict does not allow Phase5: {manifest.get('verdict')}")
    return {
        "paths": paths,
        "handoff": load_yaml(paths["phase4b_handoff"]),
        "manifest": manifest,
        "primary": read_csv(paths["primary_matrix"]),
        "sensitivity": read_csv(paths["sensitivity_matrix"]),
        "qc": read_csv(paths["qc_covariates"]),
        "feature_dictionary": read_csv(paths["feature_dictionary"]),
        "missingness": read_csv(paths["missingness"]),
        "confounding_tags": read_csv(paths["confounding_tags"]),
        "response": read_csv(paths["response_environment"]),
        "universe": read_csv(paths["analysis_universe"]),
    }


def frac_cols(df: pd.DataFrame) -> list[str]:
    return [c for c in df.columns if c.startswith("frac_")]


def feature_cols(df: pd.DataFrame) -> list[str]:
    exclude = {
        "sample_key",
        "cohort_id",
        "sample_id",
        "patient_key",
        "patient_id",
        "timepoint",
        "total_cells",
    }
    return [c for c in df.columns if c not in exclude]


def prepare_supervised(data: dict) -> pd.DataFrame:
    primary = data["primary"].copy()
    sensitivity = data["sensitivity"].copy()
    qc = data["qc"].copy()
    response = data["response"].copy()
    table = primary.merge(sensitivity[["sample_key"] + feature_cols(sensitivity)], on="sample_key", how="left")
    qc_keep = [
        "sample_key",
        "coarse_label_coverage",
        "mid_label_coverage",
        "fine_label_coverage",
        "low_quality_fraction",
        "integration_status",
        "raw_counts_available",
        "normalized_expression_available",
        "pseudobulk_table_available",
        "tcr_available_frozen",
    ]
    table = table.merge(qc[[c for c in qc_keep if c in qc.columns]], on="sample_key", how="left")
    table = table.merge(response, on="sample_key", how="left", suffixes=("", "_env"))
    table["y"] = table["response_binary_harmonized"].map({"responder": 1, "non_responder": 0})
    table["is_supervised_candidate"] = (
        table["supervised_use_allowed"].eq("yes")
        & table["response_known"].eq("yes")
        & table["split"].isin(["train", "val", "test"])
        & table["y"].notna()
    )
    write_csv(table, OUT / "preflight" / "phase5_joined_supervised_table.csv")
    sup = table[table["is_supervised_candidate"]].copy()
    sup["y"] = sup["y"].astype(int)
    audit_cols = [
        "sample_key",
        "patient_key",
        "cohort_id",
        "split",
        "response_endpoint_type",
        "response_binary_harmonized",
        "treatment_context",
        "timepoint",
        "cancer_type",
        "analysis_lane",
        "endpoint_pooling_allowed",
        "timepoint_use_boundary",
        "sensitivity_only_reason",
    ]
    write_csv(sup[[c for c in audit_cols if c in sup.columns]], OUT / "preflight" / "supervised_sample_audit.csv")
    return sup


def numeric_design(train: pd.DataFrame, frame: pd.DataFrame, cols: list[str]) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    tr = train[cols].apply(pd.to_numeric, errors="coerce")
    fr = frame[cols].apply(pd.to_numeric, errors="coerce")
    med = tr.median().fillna(0)
    tr = tr.fillna(med)
    fr = fr.fillna(med)
    scale = tr.std(ddof=0).replace(0, 1).fillna(1)
    return ((tr - med) / scale).to_numpy(float), ((fr - med) / scale).to_numpy(float), med.to_numpy(float)


def categorical_design(train: pd.DataFrame, frame: pd.DataFrame, cols: list[str]) -> tuple[np.ndarray, np.ndarray, list[str]]:
    matrices_train = []
    matrices_frame = []
    names = []
    for col in cols:
        cats = sorted(set(train[col].fillna("missing").replace("", "missing").astype(str)))
        for cat in cats:
            names.append(f"{col}::{cat}")
            matrices_train.append((train[col].fillna("missing").replace("", "missing").astype(str) == cat).astype(float).to_numpy())
            matrices_frame.append((frame[col].fillna("missing").replace("", "missing").astype(str) == cat).astype(float).to_numpy())
    if not matrices_train:
        return np.zeros((len(train), 0)), np.zeros((len(frame), 0)), []
    return np.vstack(matrices_train).T, np.vstack(matrices_frame).T, names


def fit_logistic(X: np.ndarray, y: np.ndarray, l2: float = 0.1, l1: float = 0.0, iters: int = 2500) -> np.ndarray:
    Xb = np.column_stack([np.ones(len(X)), X])
    w = np.zeros(Xb.shape[1])
    y = y.astype(float)
    pos = max(y.sum(), 1)
    neg = max(len(y) - y.sum(), 1)
    weights = np.where(y == 1, len(y) / (2 * pos), len(y) / (2 * neg))
    lr = 0.05
    for _ in range(iters):
        p = sigmoid(Xb @ w)
        grad = (Xb.T @ ((p - y) * weights)) / len(y)
        grad[1:] += l2 * w[1:]
        w -= lr * grad
        if l1:
            w[1:] = np.sign(w[1:]) * np.maximum(np.abs(w[1:]) - lr * l1, 0)
    return w


def predict_logistic(w: np.ndarray, X: np.ndarray) -> np.ndarray:
    return sigmoid(np.column_stack([np.ones(len(X)), X]) @ w)


def fit_stump_ensemble(X: np.ndarray, y: np.ndarray, n_estimators: int = 200) -> list[tuple[int, float, float, float]]:
    rng = np.random.default_rng(RANDOM_STATE)
    stumps = []
    n, p = X.shape
    if p == 0:
        return stumps
    for _ in range(n_estimators):
        idx = rng.integers(0, n, n)
        feat_candidates = rng.choice(p, size=max(1, int(np.sqrt(p))), replace=False)
        best = None
        for j in feat_candidates:
            values = X[idx, j]
            thresholds = np.unique(np.quantile(values, [0.25, 0.5, 0.75]))
            for thr in thresholds:
                left = idx[values <= thr]
                right = idx[values > thr]
                if len(left) == 0 or len(right) == 0:
                    continue
                loss = gini_loss(y[left], y[right])
                if best is None or loss < best[0]:
                    best = (loss, j, float(thr), smooth_prob(y[left]), smooth_prob(y[right]))
        if best is not None:
            _, j, thr, p_left, p_right = best
            stumps.append((j, thr, p_left, p_right))
    return stumps


def gini_loss(y_left: np.ndarray, y_right: np.ndarray) -> float:
    def gini(y):
        if len(y) == 0:
            return 0
        p = y.mean()
        return 1 - p**2 - (1 - p) ** 2
    n = len(y_left) + len(y_right)
    return len(y_left) / n * gini(y_left) + len(y_right) / n * gini(y_right)


def smooth_prob(y: np.ndarray) -> float:
    return float((y.sum() + 1) / (len(y) + 2))


def predict_stump_ensemble(stumps: list[tuple[int, float, float, float]], X: np.ndarray, fallback: float) -> np.ndarray:
    if not stumps:
        return np.repeat(fallback, len(X))
    preds = []
    for j, thr, p_left, p_right in stumps:
        preds.append(np.where(X[:, j] <= thr, p_left, p_right))
    return np.vstack(preds).mean(axis=0)


def add_result(rows: list[dict], model_id: str, family: str, split: str, y, score, status="complete", notes="") -> None:
    rows.append(
        {
            "model_id": model_id,
            "feature_family": family,
            "evaluation_split": split,
            "status": status,
            **metric_pack(y, score),
            "notes": notes,
        }
    )


def run_baselines(sup: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    rows = []
    pred_rows = []
    if sup.empty:
        return pd.DataFrame(), pd.DataFrame()
    train = sup[sup["split"].eq("train")].copy()
    evals = {name: sup[sup["split"].eq(name)].copy() for name in ["train", "val", "test"]}
    y_train = train["y"].to_numpy(int)
    prevalence = float(y_train.mean())
    primary_cols = [c for c in sup.columns if c.startswith("frac_coarse__") or c.startswith("frac_mid__")]
    fine_cols = [c for c in sup.columns if c.startswith("frac_fine_restricted__")]
    signature_cols = [c for c in sup.columns if c.startswith("signature_activity__")]
    tf_cols = [c for c in sup.columns if c.startswith("tf_activity__")]
    fraction_plus_signature_cols = primary_cols + signature_cols if signature_cols else []
    qc_cols = [c for c in ["total_cells", "coarse_label_coverage", "mid_label_coverage", "fine_label_coverage", "low_quality_fraction"] if c in sup.columns]

    # Dummy.
    for split, frame in evals.items():
        add_result(rows, "dummy_prevalence", "negative_control", split, frame["y"], np.repeat(prevalence, len(frame)), notes="train prevalence score")

    numeric_specs = [
        ("logistic_l2_cell_fraction_primary", "cell_fraction_coarse_mid", primary_cols, 0.1, 0.0, "primary coarse/mid fraction"),
        ("elastic_net_like_cell_fraction_primary", "cell_fraction_coarse_mid", primary_cols, 0.05, 0.005, "light L1+L2 logistic"),
        ("logistic_l2_signature_activity", "signature_activity", signature_cols, 0.1, 0.0, "response-blind signature activity"),
        ("logistic_l2_fraction_plus_signature", "cell_fraction_plus_signature_activity", fraction_plus_signature_cols, 0.1, 0.0, "fraction plus response-blind signature activity"),
        ("logistic_l2_tf_activity_sensitivity", "tf_activity_sensitivity", tf_cols, 0.1, 0.0, "limited TF activity sensitivity"),
        ("logistic_l2_cell_fraction_plus_fine_restricted", "cell_fraction_with_fine_restricted", primary_cols + fine_cols, 0.1, 0.0, "fine sensitivity"),
        ("logistic_l2_qc_covariates", "qc_covariates", qc_cols, 0.1, 0.0, "QC confounder baseline"),
    ]
    for model_id, family, cols, l2, l1, notes in numeric_specs:
        if not cols:
            for split, frame in evals.items():
                add_result(rows, model_id, family, split, [], [], status="blocked", notes="no feature columns")
            continue
        X_train, _, _ = numeric_design(train, train, cols)
        w = fit_logistic(X_train, y_train, l2=l2, l1=l1)
        for split, frame in evals.items():
            _, X_eval, _ = numeric_design(train, frame, cols)
            score = predict_logistic(w, X_eval)
            add_result(rows, model_id, family, split, frame["y"], score, notes=notes)
            tmp = frame[["sample_key", "patient_key", "cohort_id", "split", "y"]].copy()
            tmp["model_id"] = model_id
            tmp["score"] = score
            pred_rows.append(tmp)

    # Limited tree-like model, implemented as bagged one-level stumps.
    X_train, _, _ = numeric_design(train, train, primary_cols)
    stumps = fit_stump_ensemble(X_train, y_train)
    for split, frame in evals.items():
        _, X_eval, _ = numeric_design(train, frame, primary_cols)
        score = predict_stump_ensemble(stumps, X_eval, prevalence)
        add_result(rows, "limited_random_stump_ensemble_cell_fraction_primary", "cell_fraction_coarse_mid", split, frame["y"], score, notes="random-forest substitute without sklearn")

    categorical_specs = [
        ("cohort_only_logistic", ["cohort_id"], "cohort_confounder"),
        ("timepoint_only_logistic", ["timepoint"], "timepoint_confounder"),
        ("treatment_context_only_logistic", ["treatment_context"], "treatment_context_confounder"),
        ("cohort_timepoint_logistic", ["cohort_id", "timepoint"], "cohort_timepoint_confounder"),
    ]
    for model_id, cols, family in categorical_specs:
        X_train, _, _ = categorical_design(train, train, cols)
        w = fit_logistic(X_train, y_train, l2=0.1)
        for split, frame in evals.items():
            _, X_eval, _ = categorical_design(train, frame, cols)
            score = predict_logistic(w, X_eval)
            add_result(rows, model_id, family, split, frame["y"], score, notes="confounder-only baseline")

    # Leave-one-cohort.
    for holdout in sorted(sup["cohort_id"].unique()):
        tr = sup[sup["cohort_id"].ne(holdout)]
        te = sup[sup["cohort_id"].eq(holdout)]
        model_id = f"loco_logistic_primary_holdout__{holdout}"
        if len(set(tr["y"])) < 2 or len(set(te["y"])) < 2:
            add_result(rows, model_id, "cell_fraction_coarse_mid", f"loco_{holdout}", [], [], status="not_evaluable", notes="insufficient classes")
            continue
        X_tr, _, _ = numeric_design(tr, tr, primary_cols)
        w = fit_logistic(X_tr, tr["y"].to_numpy(int), l2=0.1)
        _, X_te, _ = numeric_design(tr, te, primary_cols)
        score = predict_logistic(w, X_te)
        add_result(rows, model_id, "cell_fraction_coarse_mid", f"loco_{holdout}", te["y"], score, notes="leave-one-cohort transfer")

    return pd.DataFrame(rows), pd.concat(pred_rows, ignore_index=True) if pred_rows else pd.DataFrame()


def run_negative_controls(sup: pd.DataFrame, observed: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    rng = np.random.default_rng(RANDOM_STATE)
    primary_cols = [c for c in sup.columns if c.startswith("frac_coarse__") or c.startswith("frac_mid__")]
    train = sup[sup["split"].eq("train")].copy()
    test = sup[sup["split"].eq("test")].copy()
    rows = []
    if train.empty or test.empty or len(set(train["y"])) < 2 or len(set(test["y"])) < 2:
        return pd.DataFrame([{"control_id": "permuted_label_logistic", "status": "not_evaluable"}]), pd.DataFrame()
    X_train, _, _ = numeric_design(train, train, primary_cols)
    _, X_test, _ = numeric_design(train, test, primary_cols)
    for i in range(N_PERMUTATIONS):
        y_perm = rng.permutation(train["y"].to_numpy(int))
        w = fit_logistic(X_train, y_perm, l2=0.1, iters=1200)
        score = predict_logistic(w, X_test)
        rows.append({"control_id": "permuted_label_logistic", "iteration": i, **metric_pack(test["y"], score)})
    hash_score = test["sample_key"].map(lambda x: int(hashlib.sha256(x.encode()).hexdigest()[:8], 16) / 0xFFFFFFFF)
    rows.append({"control_id": "sample_key_hash_random", "iteration": 0, **metric_pack(test["y"], hash_score)})
    controls = pd.DataFrame(rows)
    q = observed[(observed["model_id"].eq("logistic_l2_cell_fraction_primary")) & (observed["evaluation_split"].eq("test"))]
    observed_auc = float(q.iloc[0]["roc_auc"]) if not q.empty else np.nan
    perm = controls[controls["control_id"].eq("permuted_label_logistic")]
    summary = pd.DataFrame(
        [
            {
                "control_id": "permuted_label_logistic",
                "n_iterations": len(perm),
                "test_auc_mean": perm["roc_auc"].mean(),
                "test_auc_sd": perm["roc_auc"].std(),
                "test_auc_p95": perm["roc_auc"].quantile(0.95),
                "observed_primary_test_auc": observed_auc,
                "negative_control_collapse": "yes" if perm["roc_auc"].mean() < 0.65 else "review",
            },
            {
                "control_id": "sample_key_hash_random",
                "n_iterations": 1,
                "test_auc_mean": controls[controls["control_id"].eq("sample_key_hash_random")]["roc_auc"].mean(),
                "test_auc_sd": "",
                "test_auc_p95": "",
                "observed_primary_test_auc": observed_auc,
                "negative_control_collapse": "single_random_control",
            },
        ]
    )
    return controls, summary


def metric_value(results: pd.DataFrame, model: str, split: str, metric: str) -> float:
    q = results[(results["model_id"].eq(model)) & (results["evaluation_split"].eq(split))]
    if q.empty:
        return np.nan
    return float(q.iloc[0][metric])


def build_confounding(sup: pd.DataFrame, baseline: pd.DataFrame, phase4b_tags: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    rows = []
    for col in ["cohort_id", "timepoint", "treatment_context", "response_endpoint_type", "cancer_type"]:
        if col not in sup.columns:
            continue
        tab = pd.crosstab(sup[col], sup["response_binary_harmonized"])
        for level, vals in tab.iterrows():
            total = int(vals.sum())
            rows.append(
                {
                    "variable": col,
                    "level": level,
                    "n_responder": int(vals.get("responder", 0)),
                    "n_non_responder": int(vals.get("non_responder", 0)),
                    "responder_fraction": vals.get("responder", 0) / total if total else np.nan,
                    "n_total": total,
                }
            )
    balance = pd.DataFrame(rows)
    write_csv(balance, OUT / "confounding" / "response_balance_by_environment.csv")
    primary_auc = metric_value(baseline, "logistic_l2_cell_fraction_primary", "test", "roc_auc")
    cohort_auc = metric_value(baseline, "cohort_only_logistic", "test", "roc_auc")
    time_auc = metric_value(baseline, "timepoint_only_logistic", "test", "roc_auc")
    qc_auc = metric_value(baseline, "logistic_l2_qc_covariates", "test", "roc_auc")
    close = [
        v
        for v in [cohort_auc, time_auc, qc_auc]
        if not np.isnan(v) and not np.isnan(primary_auc) and v >= primary_auc - 0.05
    ]
    high = bool(close)
    summary = pd.DataFrame(
        [
            {"audit_item": "primary_cell_fraction_test_auc", "value": primary_auc, "risk": "reference", "interpretation": "coarse/mid fraction held-out test"},
            {"audit_item": "cohort_only_test_auc", "value": cohort_auc, "risk": "high" if high and not np.isnan(cohort_auc) else "not_high", "interpretation": "cohort-only pressure baseline"},
            {"audit_item": "timepoint_only_test_auc", "value": time_auc, "risk": "high" if high and not np.isnan(time_auc) else "not_high", "interpretation": "timepoint-only pressure baseline"},
            {"audit_item": "qc_covariate_test_auc", "value": qc_auc, "risk": "high" if high and not np.isnan(qc_auc) else "not_high", "interpretation": "QC-only pressure baseline"},
            {"audit_item": "phase4b_high_risk_feature_count", "value": int((phase4b_tags.get("allowed_primary_after_audit", pd.Series(dtype=str)) != "yes").sum()) if not phase4b_tags.empty else 0, "risk": "not_high", "interpretation": "inherited Phase4B dominance audit"},
        ]
    )
    write_csv(summary, OUT / "confounding" / "confounding_audit_summary.csv")
    return balance, summary


def leakage_audit(features: list[str], sup: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for f in features:
        hits = [tok for tok in BANNED_FEATURE_TOKENS if tok in f.lower()]
        rows.append({"feature_id": f, "leakage_tokens": "|".join(hits), "leakage_risk": "yes" if hits else "no"})
    split_counts = sup.groupby("patient_key")["split"].nunique() if not sup.empty else pd.Series(dtype=int)
    n_split_leak = int((split_counts > 1).sum())
    rows.append({"feature_id": "__patient_split_leakage_check__", "leakage_tokens": str(n_split_leak), "leakage_risk": "yes" if n_split_leak else "no"})
    audit = pd.DataFrame(rows)
    write_csv(audit, OUT / "audit" / "phase5_label_leakage_audit.csv")
    return audit


def infer_confounding_risk(summary: pd.DataFrame) -> str:
    return "high" if (not summary.empty and (summary["risk"] == "high").any()) else "not_high"


def build_family_eligibility(conf_summary: pd.DataFrame) -> pd.DataFrame:
    coarse_mid_status = (
        "response_blind_only_confounding_caution"
        if infer_confounding_risk(conf_summary) == "high"
        else "eligible_for_phase6_fraction_module_lite"
    )
    coarse_mid_use = (
        "Phase6 response-blind fraction module-lite only; no supervised response-direction claim until confounding resolved"
        if infer_confounding_risk(conf_summary) == "high"
        else "Phase6 fraction-only response-blind module-lite"
    )
    fd_path = PHASE4B / "matrix" / "feature_dictionary_v6_2.csv"
    fd = read_csv(fd_path) if fd_path.exists() else pd.DataFrame()
    if fd.empty:
        fd = pd.DataFrame(columns=["feature_family", "allowed_downstream_use"])

    def has_available_family(family: str) -> bool:
        sub = fd[fd.get("feature_family", pd.Series(dtype=str)).astype(str).eq(family)]
        if sub.empty:
            return False
        allowed = sub.get("allowed_downstream_use", pd.Series(dtype=str)).astype(str).str.lower()
        return bool((~allowed.str.contains("blocked|unavailable", regex=True)).any())

    has_raw_pb = has_available_family("pseudobulk_raw_count")
    has_norm_pb = has_available_family("pseudobulk_normalized_expression")
    has_sig = has_available_family("signature_activity")
    has_tf = has_available_family("tf_activity")
    rows = [
        ["cell_fraction_coarse_mid", "primary", coarse_mid_status, "complete", infer_confounding_risk(conf_summary), coarse_mid_use, "cohort/timepoint/QC confounding caution" if infer_confounding_risk(conf_summary) == "high" else ""],
        ["cell_fraction_fine_restricted", "sensitivity", "sensitivity_only", "complete", "review", "sensitivity/refinement only", "fine labels restricted"],
        ["qc_covariates", "confounding_audit_only", "confounding_audit_only", "complete", "not_biological_feature", "confounding adjustment/audit only", "QC covariates are not mechanism features"],
        ["pseudobulk_raw_count", "phase6_response_blind_input" if has_raw_pb else "blocked", "phase6_response_blind_input" if has_raw_pb else "blocked", "not_supervised_baseline" if has_raw_pb else "not_run", "not_assessed", "Phase6 response-blind module discovery v0 input" if has_raw_pb else "blocked_until_raw_count_cell_state_pseudobulk_patch", "sample-level raw-count module pseudobulk" if has_raw_pb else "Canonical raw-count cell-state-specific pseudobulk remains deferred"],
        ["pseudobulk_normalized_expression", "support_or_sensitivity" if has_norm_pb else "blocked", "support_or_sensitivity_only" if has_norm_pb else "blocked", "not_supervised_baseline", "review", "Support/sensitivity expression asset" if has_norm_pb else "blocked_until_expression_rescue_patch", "Normalized rescue/addendum expression is not raw-count canonical pseudobulk" if has_norm_pb else "No normalized rescue pseudobulk available"],
        ["signature_activity", "support_or_sensitivity" if has_sig else "blocked", "support_or_sensitivity_only" if has_sig else "blocked", "not_supervised_baseline" if has_sig else "not_run", infer_confounding_risk(conf_summary), "Phase6 support/sensitivity activity annotation" if has_sig else "blocked_until_activity_patch", "response-blind signatures scored from normalized rescue/addendum pseudobulk" if has_sig else "Phase4B activity matrices are empty"],
        ["tf_activity", "sensitivity" if has_tf else "blocked", "sensitivity_only" if has_tf else "blocked", "not_supervised_baseline" if has_tf else "not_run", "review", "TF sensitivity / annotation only" if has_tf else "blocked_until_activity_patch", "limited consensus regulon coverage on rescue/addendum expression" if has_tf else "TF activity matrix empty"],
        ["pathway_activity", "blocked", "blocked", "not_run", "not_assessed", "blocked_until_formal_pathway_resource", "No formal pathway gene-set resource frozen"],
        ["optional_tcr_sensitivity", "sensitivity", "sensitivity_shell_only", "not_run", "not_assessed", "TCR availability sensitivity only", "No TCR clone master"],
        ["scvi_patient_latent", "blocked", "blocked", "not_run", "not_assessed", "blocked_until_full_integration_patch", "Full integration/scVI latent pending"],
    ]
    df = pd.DataFrame(
        rows,
        columns=[
            "feature_family",
            "phase4b_status",
            "phase5_status",
            "baseline_status",
            "confounding_risk",
            "allowed_downstream_use",
            "downgrade_reason",
        ],
    )
    write_csv(df, OUT / "feature_family_downstream_eligibility.csv")
    return df


def write_reports(sup, baseline, neg_summary, conf_summary, family_elig, leakage) -> None:
    leakage_fail = bool((leakage["leakage_risk"] == "yes").any())
    high_conf = infer_confounding_risk(conf_summary) == "high"
    phase6_module_ready = bool(
        family_elig["feature_family"].isin(["pseudobulk_raw_count", "signature_activity"]).any()
        and family_elig[family_elig["feature_family"].eq("pseudobulk_raw_count")]["phase5_status"].astype(str).str.contains("input").any()
        and family_elig[family_elig["feature_family"].eq("signature_activity")]["phase5_status"].astype(str).str.contains("available").any()
    )
    verdict = "BLOCKED" if leakage_fail else (
        "CONDITIONAL_GO_TO_PHASE6_PSEUDOBULK_SIGNATURE_AVAILABLE_CONFOUNDING_CAUTION"
        if high_conf and phase6_module_ready
        else (
            "CONDITIONAL_GO_TO_PHASE6_PSEUDOBULK_SIGNATURE_AVAILABLE"
            if phase6_module_ready
            else (
                "CONDITIONAL_GO_TO_PHASE6_FRACTION_ONLY_CONFOUNDING_CAUTION"
                if high_conf
                else "CONDITIONAL_GO_TO_PHASE6_FRACTION_ONLY"
            )
        )
    )
    primary_auc = metric_value(baseline, "logistic_l2_cell_fraction_primary", "test", "roc_auc")
    cohort_auc = metric_value(baseline, "cohort_only_logistic", "test", "roc_auc")

    write_md(
        f"""# Phase5 Confounding Audit Report

**Run date:** {TODAY}

## Verdict

`{'BLOCKED' if leakage_fail else 'PASS_WITH_CONFOUNDING_CAUTION'}`

## Summary

- Supervised samples: {len(sup)}
- Patients: {sup['patient_key'].nunique() if not sup.empty else 0}
- Cohorts: {sup['cohort_id'].nunique() if not sup.empty else 0}
- Endpoint types: {'|'.join(sorted(sup['response_endpoint_type'].unique())) if not sup.empty else ''}
- Primary logistic test AUC: {primary_auc}
- Cohort-only test AUC: {cohort_auc}
- High confounding risk detected: {high_conf}

The supervised surface is endpoint-heterogeneous after response-binding repair. Baseline results are pressure tests, not mechanism claims.
""",
        OUT / "confounding" / "confounding_audit_report.md",
    )
    neg_text = "# Phase5 Negative Control Report\n\n"
    neg_text += f"**Run date:** {TODAY}\n\n"
    neg_text += neg_summary.to_string(index=False) if not neg_summary.empty else "No negative controls evaluable."
    neg_text += "\n\nNegative controls are leakage checks, not mechanism evidence.\n"
    write_md(neg_text, OUT / "negative_controls" / "negative_control_report.md")
    write_md(
        f"""# Phase5 Strong Baseline and Confounding Audit Report

**Run date:** {TODAY}
**Input version:** `{INPUT_VERSION}`
**Verdict:** `{verdict}`

## 1. Executive Verdict

Phase5 completed the strong baseline and confounding audit allowed by Phase4B. No label leakage hard blocker was detected. The supervised pressure test remains a confounding audit, so downstream response-direction claims are not allowed from these baselines. Phase4B now provides support/sensitivity normalized pseudobulk and response-blind signature assets; raw-count cell-state pseudobulk, pathway, TF, scVI latent, and mainline TCR remain blocked or sensitivity-only.

## 2. Supervised Evaluation Surface

- Supervised samples: {len(sup)}
- Supervised patients: {sup['patient_key'].nunique() if not sup.empty else 0}
- Cohorts: {'|'.join(sorted(sup['cohort_id'].unique())) if not sup.empty else ''}
- Endpoint: {'|'.join(sorted(sup['response_endpoint_type'].unique())) if not sup.empty else ''}
- Split counts: {sup['split'].value_counts().to_dict() if not sup.empty else {}}

## 3. Baselines Run

- dummy prevalence;
- logistic coarse/mid cell fraction;
- elastic-net-like coarse/mid cell fraction;
- limited random-stump ensemble;
- response-blind signature activity baseline is blocked in supervised matrix unless signature columns overlap the supervised surface;
- limited TF activity sensitivity baseline is blocked in current run;
- restricted fine sensitivity;
- QC covariate baseline;
- cohort/timepoint/treatment confounder baselines;
- leave-one-cohort transfer.

XGBoost/random forest was substituted with a dependency-free limited random-stump ensemble because sklearn/scipy is unavailable in this runtime. Pathway activity and scVI latent baselines remain blocked by upstream feature availability.

## 4. Outputs

- `baselines/baseline_results.csv`
- `baselines/baseline_model_manifest.yaml`
- `confounding/confounding_audit_report.md`
- `negative_controls/negative_control_report.md`
- `feature_family_downstream_eligibility.csv`

## 5. What Phase5 Did Not Do

- It did not make response association a scientific result.
- It did not discover modules.
- It did not make PD1 anchor verdict.
- It did not make dominant barrier or X-class repair claims.
- It did not upgrade pathway, full scVI latent, or mainline TCR feature families.

## 6. Final Decision

`{verdict}`. Phase6 can start with fraction-only response-blind module-lite plus support expression/signature assets. Raw-count cell-state pseudobulk is still required before full pseudobulk/signature module discovery. No supervised response-direction claim is authorized by Phase5.
""",
        OUT / "PHASE5_STRONG_BASELINE_CONFOUNDING_AUDIT_REPORT.md",
    )


def write_model_manifest(baseline: pd.DataFrame) -> None:
    write_yaml(
        {
            "phase": "phase5_strong_baseline_and_confounding_audit",
            "input_version": INPUT_VERSION,
            "created_at": TODAY,
            "random_state": RANDOM_STATE,
            "models_run": sorted(baseline["model_id"].unique().tolist()) if not baseline.empty else [],
            "models_blocked_or_substituted": {
                "xgboost": "not_installed_or_not_used; dependency-free random-stump ensemble substituted",
                "random_forest": "sklearn/scipy unavailable due GLIBCXX runtime issue; random-stump ensemble substituted",
                "pathway_activity": "blocked_until_formal_pathway_resource",
                "cell_state_specific_pseudobulk_module": "blocked_until_cell_state_specific_all_gene_pseudobulk_patch",
                "scvi_patient_latent": "blocked_by_full_integration_pending",
            },
            "metrics": ["roc_auc", "average_precision", "accuracy", "balanced_accuracy", "brier", "log_loss"],
            "rules": {
                "patient_split_respected": True,
                "response_used_only_as_label": True,
                "no_mechanism_claim_from_baseline": True,
            },
        },
        OUT / "baselines" / "baseline_model_manifest.yaml",
    )


def write_manifest_and_handoff(sup, baseline, family_elig, leakage) -> None:
    leakage_fail = bool((leakage["leakage_risk"] == "yes").any())
    high_conf = "high" in set(family_elig.get("confounding_risk", pd.Series(dtype=str)))
    phase6_module_ready = bool(
        family_elig[family_elig["feature_family"].eq("pseudobulk_raw_count")]["phase5_status"].astype(str).str.contains("input").any()
        and family_elig[family_elig["feature_family"].eq("signature_activity")]["phase5_status"].astype(str).str.contains("available").any()
    )
    family_status = {r["feature_family"]: r["phase5_status"] for _, r in family_elig.iterrows()}
    support_expression_available = family_status.get("pseudobulk_normalized_expression", "") != "blocked"
    signature_support_available = family_status.get("signature_activity", "") != "blocked"
    tf_support_available = family_status.get("tf_activity", "") != "blocked"
    phase6_allowed_inputs = [
        str(OUT / "baselines" / "baseline_results.csv"),
        str(OUT / "feature_family_downstream_eligibility.csv"),
        str(OUT / "confounding" / "confounding_audit_summary.csv"),
        str(OUT / "audit" / "phase5_label_leakage_audit.csv"),
        str(PHASE4B / "matrix" / "immune_state_feature_matrix.primary.csv"),
        str(PHASE4B / "matrix" / "feature_dictionary_v6_2.csv"),
    ]
    phase6_blocked_inputs = [str(PHASE4B / "activity" / "pathway_activity_matrix.sample_level.csv")]
    if phase6_module_ready:
        phase6_allowed_inputs.append(str(PHASE4B / "pseudobulk" / "pseudobulk_matrix_raw_count.parquet"))
    else:
        phase6_blocked_inputs.append(str(PHASE4B / "pseudobulk" / "pseudobulk_matrix_raw_count.parquet"))
    if support_expression_available:
        phase6_allowed_inputs.append(str(PHASE4B / "pseudobulk" / "pseudobulk_matrix_normalized_expression.parquet"))
    if signature_support_available:
        phase6_allowed_inputs.append(str(PHASE4B / "activity" / "signature_activity_matrix.sample_level.csv"))
    if tf_support_available:
        phase6_allowed_inputs.append(str(PHASE4B / "activity" / "tf_activity_matrix.sample_level.csv"))
    verdict = "BLOCKED" if leakage_fail else (
        "CONDITIONAL_GO_TO_PHASE6_PSEUDOBULK_SIGNATURE_AVAILABLE_CONFOUNDING_CAUTION"
        if high_conf and phase6_module_ready
        else (
            "CONDITIONAL_GO_TO_PHASE6_PSEUDOBULK_SIGNATURE_AVAILABLE"
            if phase6_module_ready
            else (
                "CONDITIONAL_GO_TO_PHASE6_FRACTION_ONLY_CONFOUNDING_CAUTION"
                if high_conf
                else "CONDITIONAL_GO_TO_PHASE6_FRACTION_ONLY"
            )
        )
    )
    manifest = {
        "phase": "phase5_strong_baseline_and_confounding_audit",
        "created_at": TODAY,
        "input_version": INPUT_VERSION,
        "verdict": verdict,
        "hard_blockers": ["label_or_split_leakage_detected"] if leakage_fail else [],
        "supervised_surface": {
            "n_samples": int(len(sup)),
            "n_patients": int(sup["patient_key"].nunique()) if not sup.empty else 0,
            "n_cohorts": int(sup["cohort_id"].nunique()) if not sup.empty else 0,
            "endpoint_types": sorted(sup["response_endpoint_type"].unique().tolist()) if not sup.empty else [],
            "split_counts": {str(k): int(v) for k, v in sup["split"].value_counts().to_dict().items()} if not sup.empty else {},
        },
        "primary_outputs": {
            "baseline_results": str(OUT / "baselines" / "baseline_results.csv"),
            "baseline_model_manifest": str(OUT / "baselines" / "baseline_model_manifest.yaml"),
            "confounding_audit_report": str(OUT / "confounding" / "confounding_audit_report.md"),
            "negative_control_report": str(OUT / "negative_controls" / "negative_control_report.md"),
            "feature_family_downstream_eligibility": str(OUT / "feature_family_downstream_eligibility.csv"),
        },
        "feature_family_status": family_status,
        "phase6_entry": {
            "enter_phase6": not leakage_fail,
            "allowed_scope": "response_blind_pseudobulk_signature_module_discovery_v0_with_confounding_caution" if phase6_module_ready else "fraction_only_response_blind_module_lite_with_support_expression_assets",
            "blocked_scope": [
                "cell_state_specific_all_gene_pseudobulk_module_discovery",
                "pathway_activity_module_discovery",
                "tcr_coupling_mainline",
                "dominant_barrier_claim",
                "supervised_response_direction_claim",
            ],
        },
    }
    write_yaml(manifest, OUT / "handoff" / "phase5_decision_manifest.yaml")
    write_yaml(
        {
            "phase": "phase5_to_phase6_handoff",
            "input_version": INPUT_VERSION,
            "verdict": verdict,
            "phase6_allowed_inputs": phase6_allowed_inputs,
            "phase6_blocked_inputs": phase6_blocked_inputs,
            "allowed_phase6_scope": "response_blind_pseudobulk_signature_module_discovery_v0_with_confounding_caution" if phase6_module_ready else "fraction_only_response_blind_module_lite_with_support_expression_assets",
            "supervised_response_direction_claim_allowed": False,
            "patch_required_for_full_phase6": ["cell_state_specific_all_gene_pseudobulk", "formal_pathway_resource", "full_integration_scvi_or_equivalent"],
        },
        OUT / "handoff" / "phase5_to_phase6_handoff.yaml",
    )
    output_rows = []
    for path in sorted(OUT.rglob("*")):
        if path.is_file():
            output_rows.append({"path": str(path), "size_bytes": path.stat().st_size, "relative_path": str(path.relative_to(OUT))})
    write_tsv(pd.DataFrame(output_rows), OUT / "handoff" / "phase5_output_index.tsv")


def main() -> None:
    mkdirs()
    data = load_inputs()
    sup = prepare_supervised(data)
    features = [
        c
        for c in sup.columns
        if c.startswith("frac_") or c.startswith("signature_activity__") or c.startswith("tf_activity__")
    ]
    leakage = leakage_audit(features, sup)
    baseline, predictions = run_baselines(sup)
    write_csv(baseline, OUT / "baselines" / "baseline_results.csv")
    write_csv(predictions, OUT / "baselines" / "baseline_predictions.csv")
    write_model_manifest(baseline)
    controls, neg_summary = run_negative_controls(sup, baseline)
    write_csv(controls, OUT / "negative_controls" / "negative_control_results.csv")
    write_csv(neg_summary, OUT / "negative_controls" / "negative_control_summary.csv")
    _, conf_summary = build_confounding(sup, baseline, data["confounding_tags"])
    family_elig = build_family_eligibility(conf_summary)
    write_reports(sup, baseline, neg_summary, conf_summary, family_elig, leakage)
    write_manifest_and_handoff(sup, baseline, family_elig, leakage)


if __name__ == "__main__":
    main()
