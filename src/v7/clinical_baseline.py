"""Patient-level Stage 5 clinical baseline for v7.

This module reads the audited treatment/endpoint table only for the clinical
anchor task.  It never pools cells as patients, selects features from response,
or treats a model score as a causal treatment effect.  The primary model uses
the D2 ``measurable`` feature subset; all 39 features are a pre-declared
sensitivity analysis.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Iterable

import numpy as np
import pandas as pd
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import average_precision_score, brier_score_loss, roc_auc_score
from sklearn.model_selection import StratifiedKFold, cross_val_predict
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler

from .spatial_control import sha256_file, stable_hash


PRIMARY_DATASETS = (
    "TASK01",
    "TASK02",
    "LAMBRECHT_HCC",
    "GSE301741",
    "GSE286827",
)
PRE_TIMEPOINTS = frozenset({"t0", "baseline", "pre", "pretreatment", "pre_treatment"})
RESPONSE_MAP = {"responder": 1, "non_responder": 0, "non-responder": 0}


def _normal(value: object) -> str:
    return str(value).strip().lower()


def _match_dataset(value: object) -> str | None:
    key = _normal(value)
    aliases = {
        "task01": "TASK01",
        "task02": "TASK02",
        "lambrecht_hcc": "LAMBRECHT_HCC",
        "gse301741": "GSE301741",
        "gse286827": "GSE286827",
    }
    return aliases.get(key)


def _pre_flag(value: object) -> str:
    key = _normal(value)
    if key in PRE_TIMEPOINTS:
        return "primary_pre"
    if "pre" in key and "post" not in key and "on_treatment" not in key:
        return "primary_pre"
    if "pre" in key:
        return "ambiguous_pre"
    return "not_pre"


def _read_response_table(path: str | Path) -> pd.DataFrame:
    columns = [
        "dataset_id",
        "patient_id",
        "treatment_values",
        "timepoint_values",
        "endpoint_values",
        "response_values",
        "n_samples",
        "n_timepoints",
        "conflict_flag",
        "complete_for_response",
    ]
    frame = pd.read_csv(path, sep="\t", usecols=columns, dtype=str, keep_default_na=False)
    frame["analysis_dataset"] = frame.dataset_id.map(_match_dataset)
    frame = frame[frame.analysis_dataset.notna()].copy()
    frame["response_binary"] = frame.response_values.map(RESPONSE_MAP)
    # Patient identifiers in the legacy registry differ in prefix case from
    # the Stage2 score export; case-folding is an identity normalization, not
    # a merge across distinct patients.
    frame["patient_key"] = frame.patient_id.astype(str).str.casefold()
    return frame


def _read_abundance_scores(path: str | Path) -> pd.DataFrame:
    columns = [
        "cohort_id",
        "patient_key",
        "timepoint",
        "cell_state_level",
        "feature_id",
        "coverage",
        "score_native",
        "score_standardized",
        "n_cells_used",
        "technical_status",
        "aggregation_method",
        "cell_state",
    ]
    frame = pd.read_parquet(path, columns=columns)
    frame["analysis_dataset"] = frame.cohort_id.map(_match_dataset)
    frame["patient_key"] = frame.patient_key.astype(str).str.casefold()
    frame = frame[
        frame.analysis_dataset.notna()
        & frame.cell_state_level.eq("coarse")
        & frame.aggregation_method.eq("abundance_weighted")
        & frame.cell_state.eq("ALL_STATES")
        & frame.technical_status.eq("valid")
    ].copy()
    frame["pre_status"] = frame.timepoint.map(_pre_flag)
    # There should be one row per patient/timepoint/feature.  Mean aggregation
    # is only a defensive collapse for duplicated exported rows.
    keys = ["analysis_dataset", "patient_key", "timepoint", "feature_id"]
    value_cols = ["coverage", "score_native", "score_standardized", "n_cells_used"]
    return frame.groupby(keys, as_index=False, dropna=False)[value_cols].mean().merge(
        frame[keys + ["pre_status"]].drop_duplicates(keys), on=keys, how="left"
    )


def _read_cell_composition(path: str | Path, allowed_patients: set[str]) -> pd.DataFrame:
    columns = ["cohort_id", "patient_key", "timepoint", "cell_state_level", "cell_state", "n_cells_used"]
    frame = pd.read_parquet(path, columns=columns)
    frame["analysis_dataset"] = frame.cohort_id.map(_match_dataset)
    frame["patient_key"] = frame.patient_key.astype(str).str.casefold()
    frame = frame[
        frame.analysis_dataset.notna()
        & frame.cell_state_level.eq("coarse")
        & frame.patient_key.isin(allowed_patients)
    ].copy()
    # n_cells_used repeats once per feature in the upstream score export.
    frame = frame.drop_duplicates(
        ["analysis_dataset", "patient_key", "timepoint", "cell_state"]
    )
    frame["n_cells_used"] = pd.to_numeric(frame.n_cells_used, errors="coerce").fillna(0.0)
    total = frame.groupby(["analysis_dataset", "patient_key", "timepoint"], dropna=False).n_cells_used.transform("sum")
    frame["fraction"] = np.divide(
        frame.n_cells_used,
        total,
        out=np.full(len(frame), np.nan, dtype=float),
        where=total.to_numpy(float) > 0,
    )
    # Cell-state choice is prevalence based and response-blind.  Retain the
    # largest 12 coarse states to keep the clinical models low capacity.
    keep = (
        frame.groupby("cell_state", dropna=False).n_cells_used.sum().sort_values(ascending=False).head(12).index
    )
    frame = frame[frame.cell_state.isin(keep)].copy()
    frame["composition_feature"] = "composition__" + frame.cell_state.astype(str)
    pivot = frame.pivot_table(
        index=["analysis_dataset", "patient_key", "timepoint"],
        columns="composition_feature",
        values="fraction",
        aggfunc="mean",
    ).reset_index()
    pivot.columns.name = None
    return pivot


def _make_patient_table(
    response: pd.DataFrame,
    scores: pd.DataFrame,
    composition: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    pre_scores = scores[scores.pre_status.eq("primary_pre")].copy()
    pre_scores = pre_scores.merge(
        response[["analysis_dataset", "patient_key", "response_binary", "treatment_values", "endpoint_values", "complete_for_response", "conflict_flag"]],
        on=["analysis_dataset", "patient_key"],
        how="inner",
    )
    pre_scores = pre_scores[pre_scores.response_binary.notna()].copy()
    # If a patient has more than one primary pre record, average within patient
    # and retain the number of source rows as a QC field.
    abundance = pre_scores.pivot_table(
        index=["analysis_dataset", "patient_key"],
        columns="feature_id",
        values="score_standardized",
        aggfunc="mean",
    ).reset_index()
    abundance.columns = [str(c) if c in ("analysis_dataset", "patient_key") else f"abundance__{c}" for c in abundance.columns]
    qc = pre_scores.groupby(["analysis_dataset", "patient_key"], as_index=False).agg(
        qc__n_features=("feature_id", "nunique"),
        qc__coverage_median=("coverage", "median"),
        qc__n_cells=("n_cells_used", "median"),
        qc__n_pre_rows=("timepoint", "nunique"),
        treatment_values=("treatment_values", "first"),
        endpoint_values=("endpoint_values", "first"),
        complete_for_response=("complete_for_response", "first"),
        conflict_flag=("conflict_flag", "first"),
        response_binary=("response_binary", "first"),
    )
    table = qc.merge(abundance, on=["analysis_dataset", "patient_key"], how="left")
    pre_composition = composition[composition.timepoint.map(_pre_flag).eq("primary_pre")].copy()
    pre_composition = pre_composition.drop(columns=["timepoint"], errors="ignore")
    pre_composition = pre_composition.groupby(["analysis_dataset", "patient_key"], as_index=False).mean(numeric_only=True)
    table = table.merge(pre_composition, on=["analysis_dataset", "patient_key"], how="left")
    table["patient_row_id"] = table.analysis_dataset.astype(str) + "::" + table.patient_key.astype(str)
    return table.sort_values(["analysis_dataset", "patient_key"], kind="mergesort").reset_index(drop=True), pre_scores


def _pre_selection_audit(response: pd.DataFrame, scores: pd.DataFrame) -> pd.DataFrame:
    """Record why each audited patient/timepoint did or did not enter baseline."""
    joined = scores[["analysis_dataset", "patient_key", "timepoint", "pre_status"]].drop_duplicates()
    joined = joined.merge(
        response[["analysis_dataset", "patient_key", "response_binary", "complete_for_response", "conflict_flag"]],
        on=["analysis_dataset", "patient_key"],
        how="right",
    )
    return (
        joined.groupby(["analysis_dataset", "patient_key"], as_index=False, dropna=False)
        .agg(
            n_score_timepoints=("timepoint", "nunique"),
            n_primary_pre=("pre_status", lambda x: int((x == "primary_pre").sum())),
            n_ambiguous_pre=("pre_status", lambda x: int((x == "ambiguous_pre").sum())),
            response_binary=("response_binary", "first"),
            complete_for_response=("complete_for_response", "first"),
            conflict_flag=("conflict_flag", "first"),
        )
        .assign(
            baseline_status=lambda x: np.select(
                [
                    x.n_primary_pre.ge(1) & x.response_binary.notna(),
                    x.n_primary_pre.eq(0) & x.n_ambiguous_pre.ge(1),
                    x.response_binary.isna(),
                ],
                ["INCLUDED_PRIMARY_PRE", "EXCLUDED_AMBIGUOUS_PRE", "EXCLUDED_RESPONSE_UNAVAILABLE"],
                default="EXCLUDED_PRE_SCORE_UNAVAILABLE",
            )
        )
        .sort_values(["analysis_dataset", "patient_key"], kind="mergesort")
        .reset_index(drop=True)
    )


def _feature_sets(table: pd.DataFrame, reliability: pd.DataFrame) -> dict[str, list[str]]:
    abundance_all = sorted(c for c in table.columns if c.startswith("abundance__"))
    measurable = set(reliability.loc[reliability.D2_class.eq("measurable"), "feature_id"].astype(str))
    abundance_primary = sorted(f"abundance__{feature}" for feature in measurable if f"abundance__{feature}" in table)
    composition = sorted(c for c in table.columns if c.startswith("composition__"))
    qc = sorted(c for c in table.columns if c.startswith("qc__"))
    return {
        "qc": qc,
        "composition": composition,
        "abundance_primary": abundance_primary,
        "abundance_all": abundance_all,
        "joint_primary": sorted(set(qc + composition + abundance_primary)),
        "joint_all": sorted(set(qc + composition + abundance_all)),
    }


def _model_pipeline() -> Pipeline:
    return Pipeline(
        [
            ("impute", SimpleImputer(strategy="median", add_indicator=True)),
            ("scale", StandardScaler(with_mean=False)),
            (
                "model",
                LogisticRegression(
                    C=1.0,
                    penalty="l2",
                    solver="liblinear",
                    class_weight="balanced",
                    max_iter=1000,
                    random_state=20260911,
                ),
            ),
        ]
    )


def _bootstrap_prediction_intervals(
    y: np.ndarray,
    prediction: np.ndarray,
    *,
    seed: int,
    draws: int = 1000,
) -> dict[str, float | int | str]:
    """Patient-resampling intervals around fixed out-of-fold predictions."""
    if len(y) < 6 or len(np.unique(y)) < 2:
        return {
            "bootstrap_status": "NOT_ESTIMABLE",
            "bootstrap_draws": draws,
            "auroc_ci_low": np.nan,
            "auroc_ci_high": np.nan,
            "average_precision_ci_low": np.nan,
            "average_precision_ci_high": np.nan,
        }
    rng = np.random.default_rng(seed)
    aurocs: list[float] = []
    aps: list[float] = []
    for _ in range(draws):
        indices = rng.integers(0, len(y), size=len(y))
        sampled_y = y[indices]
        if len(np.unique(sampled_y)) < 2:
            continue
        sampled_prediction = prediction[indices]
        aurocs.append(float(roc_auc_score(sampled_y, sampled_prediction)))
        aps.append(float(average_precision_score(sampled_y, sampled_prediction)))
    if not aurocs:
        return {
            "bootstrap_status": "NOT_ESTIMABLE",
            "bootstrap_draws": draws,
            "auroc_ci_low": np.nan,
            "auroc_ci_high": np.nan,
            "average_precision_ci_low": np.nan,
            "average_precision_ci_high": np.nan,
        }
    return {
        "bootstrap_status": "ESTIMABLE",
        "bootstrap_draws": draws,
        "auroc_ci_low": float(np.quantile(aurocs, 0.025)),
        "auroc_ci_high": float(np.quantile(aurocs, 0.975)),
        "average_precision_ci_low": float(np.quantile(aps, 0.025)),
        "average_precision_ci_high": float(np.quantile(aps, 0.975)),
    }


def _evaluate_model(
    frame: pd.DataFrame,
    features: list[str],
    *,
    environment: str,
    model_name: str,
    seed: int,
    bootstrap_draws: int = 1000,
    cv_splits: list[tuple[np.ndarray, np.ndarray]] | None = None,
) -> tuple[dict[str, Any], pd.DataFrame]:
    subset = frame.dropna(subset=["response_binary"]).copy()
    y = subset.response_binary.astype(int).to_numpy()
    result: dict[str, Any] = {
        "environment": environment,
        "model": model_name,
        "n_patients": int(len(subset)),
        "n_features": int(len(features)),
        "n_responder": int(y.sum()),
        "n_non_responder": int((1 - y).sum()),
        "status": "NOT_TESTABLE",
        "auroc": np.nan,
        "average_precision": np.nan,
        "brier": np.nan,
        "folds": 0,
    }
    if len(features) == 0 or len(subset) < 6 or min(int(y.sum()), int((1 - y).sum())) < 3:
        return result, pd.DataFrame()
    n_splits = min(5, int(y.sum()), int((1 - y).sum()))
    splitter = cv_splits if cv_splits is not None else list(
        StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=seed).split(subset, y)
    )
    try:
        prediction = cross_val_predict(
            _model_pipeline(),
            subset[features],
            y,
            cv=splitter,
            method="predict_proba",
        )[:, 1]
    except (ValueError, RuntimeError) as exc:
        result["status"] = f"NOT_TESTABLE_{type(exc).__name__}"
        return result, pd.DataFrame()
    result.update(
        {
            "status": "ESTIMABLE",
            "auroc": float(roc_auc_score(y, prediction)),
            "average_precision": float(average_precision_score(y, prediction)),
            "brier": float(brier_score_loss(y, prediction)),
            "folds": n_splits,
        }
    )
    result.update(
        _bootstrap_prediction_intervals(
            y,
            prediction,
            seed=seed ^ 0x5A17,
            draws=bootstrap_draws,
        )
    )
    predictions = pd.DataFrame(
        {
            "patient_row_id": subset.patient_row_id.to_numpy(),
            "environment": environment,
            "model": model_name,
            "response_binary": y,
            "prediction": prediction,
        }
    )
    fold_ids = np.full(len(subset), -1, dtype=int)
    for fold, (_, test) in enumerate(splitter):
        fold_ids[test] = fold
    predictions["fold"] = fold_ids
    return result, predictions


def _null_metric(
    frame: pd.DataFrame, environment: str,
    cv_splits: list[tuple[np.ndarray, np.ndarray]] | None = None,
) -> dict[str, Any]:
    y = frame.response_binary.astype(int).to_numpy()
    prevalence = float(np.mean(y)) if len(y) else np.nan
    prediction = np.full(len(y), prevalence)
    if cv_splits is not None:
        for train, test in cv_splits:
            prediction[test] = float(np.mean(y[train]))
    return {
        "environment": environment,
        "model": "null_prevalence",
        "n_patients": int(len(frame)),
        "n_features": 0,
        "n_responder": int(y.sum()),
        "n_non_responder": int((1 - y).sum()),
        "status": "ESTIMABLE" if len(y) else "NOT_TESTABLE",
        "auroc": float(roc_auc_score(y, prediction)) if len(np.unique(y)) == 2 else np.nan,
        "average_precision": float(average_precision_score(y, prediction)) if len(y) else np.nan,
        "brier": float(brier_score_loss(y, prediction)) if len(y) else np.nan,
        "folds": len(cv_splits) if cv_splits is not None else 0,
    }


def _environment_groups(table: pd.DataFrame):
    """Do not pool treatment arms or endpoints inside an accession."""
    for keys, frame in table.groupby(
        ["analysis_dataset", "treatment_values", "endpoint_values"], sort=True, dropna=False
    ):
        yield "|".join(str(key) for key in keys), frame


def _response_effects(table: pd.DataFrame, feature_sets: dict[str, list[str]]) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    for environment, frame in _environment_groups(table):
        valid = frame[frame.response_binary.notna()].copy()
        for set_name in ("abundance_primary", "composition"):
            for feature in feature_sets[set_name]:
                x = pd.to_numeric(valid[feature], errors="coerce")
                r = x[valid.response_binary.eq(1)].dropna().to_numpy(float)
                n = x[valid.response_binary.eq(0)].dropna().to_numpy(float)
                if len(r) < 2 or len(n) < 2:
                    rows.append({"environment": environment, "feature": feature, "feature_set": set_name, "status": "NOT_ESTIMABLE_TOO_FEW_PER_GROUP"})
                    continue
                pooled = np.sqrt(((len(r) - 1) * np.var(r, ddof=1) + (len(n) - 1) * np.var(n, ddof=1)) / max(len(r) + len(n) - 2, 1))
                rows.append(
                    {
                        "environment": environment,
                        "feature": feature,
                        "feature_set": set_name,
                        "n_responder": len(r),
                        "n_non_responder": len(n),
                        "responder_median": float(np.median(r)),
                        "non_responder_median": float(np.median(n)),
                        "median_difference": float(np.median(r) - np.median(n)),
                        "standardized_difference": float((np.mean(r) - np.mean(n)) / pooled) if pooled > 0 else np.nan,
                        "status": "DESCRIPTIVE_ENVIRONMENT_ASSOCIATION",
                    }
                )
    return pd.DataFrame(rows)


def _paired_displacement(scores: pd.DataFrame, response: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    score = scores.copy()
    score = score.merge(response[["analysis_dataset", "patient_key", "response_binary"]], on=["analysis_dataset", "patient_key"], how="inner")
    for environment in ("TASK01", "TASK02"):
        group = score[score.analysis_dataset.eq(environment)]
        for target in sorted(set(group.timepoint) - {"T0"}):
            pre = group[group.timepoint.eq("T0")]
            post = group[group.timepoint.eq(target)]
            if pre.empty or post.empty:
                continue
            merged = pre.merge(post, on=["analysis_dataset", "patient_key", "feature_id"], suffixes=("_pre", "_post"))
            merged["delta"] = merged.score_standardized_post - merged.score_standardized_pre
            for feature, values in merged.groupby("feature_id", sort=True):
                for response_value, sub in values.groupby("response_binary_pre", dropna=False):
                    rows.append(
                        {
                            "environment": environment,
                            "target_timepoint": target,
                            "feature_id": feature,
                            "response_binary": int(response_value),
                            "n_patients": int(sub.patient_key.nunique()),
                            "median_delta": float(sub.delta.median()) if not sub.empty else np.nan,
                            "status": "DESCRIPTIVE_PAIRED_DISPLACEMENT" if sub.patient_key.nunique() >= 3 else "NOT_ESTIMABLE_TOO_FEW_PATIENTS",
                        }
                    )
    return pd.DataFrame(rows)


def run_clinical_baseline(
    *,
    output_root: str | Path = "results/v7/clinical_anchor",
    response_table_path: str | Path = "results/v7/registry/treatment_response_completeness.tsv",
    patient_scores_path: str | Path = "results/v7/ontology/scores/patient_timepoint_scores.parquet",
    expression_scores_path: str | Path = "results/v7/ontology/scores/expression_unit_scores.parquet",
    reliability_path: str | Path = "results/v7/ontology/measurement_reliability.tsv",
    seed: int = 20260911,
    response_permutations: int = 50,
    bootstrap_draws: int = 1000,
) -> dict[str, Any]:
    """Run environment-specific patient-level clinical baselines."""

    output = Path(output_root)
    output.mkdir(parents=True, exist_ok=True)
    response = _read_response_table(response_table_path)
    scores = _read_abundance_scores(patient_scores_path)
    allowed = set(response.patient_key)
    composition = _read_cell_composition(expression_scores_path, allowed)
    table, pre_scores = _make_patient_table(response, scores, composition)
    pre_audit = _pre_selection_audit(response, scores)
    reliability = pd.read_csv(reliability_path, sep="\t")
    feature_sets = _feature_sets(table, reliability)

    table.to_csv(output / "patient_level_table.tsv", sep="\t", index=False)
    pre_audit.to_csv(output / "pre_selection_audit.tsv", sep="\t", index=False)
    pre_scores.to_parquet(output / "pre_score_long.parquet", index=False)
    composition.to_csv(output / "composition_features.tsv", sep="\t", index=False)
    pd.DataFrame(
        [{"feature_set": key, "feature": feature} for key, values in feature_sets.items() for feature in values]
    ).to_csv(output / "feature_selection.tsv", sep="\t", index=False)

    metric_rows: list[dict[str, Any]] = []
    prediction_rows: list[pd.DataFrame] = []
    for environment, frame in _environment_groups(table):
        y = frame.response_binary.astype(int).to_numpy()
        fold_seed = seed ^ int(stable_hash({"environment": environment})[:8], 16)
        splits = list(StratifiedKFold(
            n_splits=min(5, int(y.sum()), int((1-y).sum())),
            shuffle=True, random_state=fold_seed,
        ).split(frame, y)) if min(int(y.sum()), int((1-y).sum())) >= 3 else None
        metric_rows.append(_null_metric(frame, environment, splits))
        for name, features in feature_sets.items():
            result, predictions = _evaluate_model(
                frame,
                features,
                environment=environment,
                model_name=name,
                seed=fold_seed,
                bootstrap_draws=bootstrap_draws,
                cv_splits=splits,
            )
            metric_rows.append(result)
            if not predictions.empty:
                prediction_rows.append(predictions)
    metrics = pd.DataFrame(metric_rows)
    metrics.to_csv(output / "environment_effects.tsv", sep="\t", index=False)
    predictions = pd.concat(prediction_rows, ignore_index=True) if prediction_rows else pd.DataFrame()
    predictions.to_parquet(output / "out_of_fold_predictions.parquet", index=False)

    effects = _response_effects(table, feature_sets)
    effects.to_csv(output / "responder_compatible_ecology.tsv", sep="\t", index=False)
    paired = _paired_displacement(scores, response)
    paired.to_csv(output / "paired_displacement.tsv", sep="\t", index=False)

    # Response permutation is a negative control for the primary abundance
    # model.  It is performed inside each environment and never used to tune C.
    permutation_rows: list[dict[str, Any]] = []
    for environment, frame in _environment_groups(table):
        features = feature_sets["abundance_primary"]
        valid = frame.dropna(subset=["response_binary"]).copy()
        if len(features) == 0 or len(valid) < 6 or valid.response_binary.nunique() < 2:
            permutation_rows.append({"environment": environment, "status": "NOT_TESTABLE", "n_permutations": response_permutations})
            continue
        rng = np.random.default_rng(seed ^ int(stable_hash({"environment": environment})[:8], 16))
        values: list[float] = []
        for _ in range(response_permutations):
            shuffled = valid.copy()
            shuffled["response_binary"] = rng.permutation(shuffled.response_binary.to_numpy())
            result, _ = _evaluate_model(
                shuffled,
                features,
                environment=environment,
                model_name="abundance_primary_response_permuted",
                seed=seed + 17,
                bootstrap_draws=0,
            )
            if result["status"] == "ESTIMABLE":
                values.append(float(result["auroc"]))
        permutation_rows.append(
            {
                "environment": environment,
                "status": "ESTIMABLE" if values else "NOT_TESTABLE",
                "n_permutations": response_permutations,
                "n_successful": len(values),
                "null_auroc_median": float(np.median(values)) if values else np.nan,
                "null_auroc_q025": float(np.quantile(values, 0.025)) if values else np.nan,
                "null_auroc_q975": float(np.quantile(values, 0.975)) if values else np.nan,
            }
        )
    pd.DataFrame(permutation_rows).to_csv(output / "confounding_baselines.tsv", sep="\t", index=False)

    manifest = {
        "schema": "v7.stage5.clinical_baseline.run.v1",
        "status": "S5_BASELINE_COMPLETE_WITH_LIMITATIONS",
        "datasets": sorted(table.analysis_dataset.unique()),
        "n_patients": int(table.patient_row_id.nunique()),
        "n_patients_by_dataset": table.groupby("analysis_dataset").patient_row_id.nunique().to_dict(),
        "response_read": True,
        "environment_policy": "dataset_x_treatment_x_endpoint",
        "fold_policy": "shared_patient_folds_across_models",
        "baseline_policy": "training_fold_prevalence",
        "spatial_architecture_read": False,
        "feature_sets": {key: len(value) for key, value in feature_sets.items()},
        "response_permutations": response_permutations,
        "bootstrap_draws": bootstrap_draws,
        "response_table_sha256": sha256_file(response_table_path),
        "patient_scores_sha256": sha256_file(patient_scores_path),
        "expression_scores_sha256": sha256_file(expression_scores_path),
        "reliability_sha256": sha256_file(reliability_path),
        "paired_spatial_rewiring_status": "NOT_IDENTIFIABLE_WITH_CURRENT_AUDITED_ASSETS",
    }
    (output / "STAGE5_RUN_MANIFEST.json").write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    report = f"""# v7 Stage 5 患者级临床基线报告

## 状态

**S5_BASELINE_COMPLETE_WITH_LIMITATIONS**。本轮将疗效字段仅用于临床锚定，按治疗环境和癌种分开分析；没有读取空间结构候选，也没有把细胞当作独立患者。

## 数据与方法

- 纳入 {manifest['n_patients']} 位患者，队列为 {"、".join(manifest['datasets'])}。
- 预治疗样本按 T0/baseline/pre 口径选择；预治疗缺失、重复或冲突在患者表中保留。
- 主要分子模型使用 D2 `measurable` 特征；39 特征模型作为预先定义的敏感性分析。
- 细胞组成使用 scRNA pseudobulk 的 coarse-state 细胞数比例，按无 response 的总体丰度选择最多12个状态；这是患者级组成特征，不是疗效标签。
- 每个环境内部做患者分层交叉验证，报告 AUROC、PR-AUC 和 Brier；同时运行 response permutation 负对照。

## 结论边界

输出可以回答哪些分子/细胞状态在各治疗环境内与疗效共同变化，以及空间结果是否值得继续桥接。它不能单独证明 PD-1 因果机制，也不能证明 X 的增量疗效。TASK01/TASK02 的配对位移仅作描述；直接纵向空间重排仍是 `NOT_IDENTIFIABLE_WITH_CURRENT_AUDITED_ASSETS`。

## 复现

```bash
env LD_LIBRARY_PATH=/opt/anaconda3/lib PYTHONPATH=src \\
  python scripts/v7/clinical/run_stage5.py
```
"""
    (output / "CLINICAL_ANCHOR_REPORT.md").write_text(report, encoding="utf-8")
    return {**manifest, "output_root": str(output.resolve())}


__all__ = ["PRIMARY_DATASETS", "run_clinical_baseline"]
