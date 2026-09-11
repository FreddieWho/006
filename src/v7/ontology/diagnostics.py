"""Response-blind reliability diagnostics for Stage 2 measurements."""

from __future__ import annotations

from collections.abc import Iterable
from itertools import combinations

import numpy as np
import pandas as pd


def spearman_correlation(left: Iterable[float], right: Iterable[float]) -> float:
    frame = pd.DataFrame(
        {
            "left": pd.to_numeric(pd.Series(left), errors="coerce"),
            "right": pd.to_numeric(pd.Series(right), errors="coerce"),
        }
    )
    frame = frame.loc[np.isfinite(frame.left) & np.isfinite(frame.right)]
    if len(frame) < 3 or frame.left.nunique() < 2 or frame.right.nunique() < 2:
        return float("nan")
    ranks = frame.rank(method="average")
    return float(np.corrcoef(ranks.left.to_numpy(), ranks.right.to_numpy())[0, 1])


def bootstrap_median_interval(
    values: Iterable[float],
    *,
    replicates: int = 1000,
    seed: int = 20260831,
    alpha: float = 0.05,
) -> tuple[float, float, float]:
    array = np.asarray(list(values), dtype=float)
    array = array[np.isfinite(array)]
    if array.size == 0:
        return float("nan"), float("nan"), float("nan")
    if array.size == 1:
        value = float(array[0])
        return value, value, value
    rng = np.random.default_rng(seed)
    samples = rng.choice(array, size=(replicates, array.size), replace=True)
    estimates = np.median(samples, axis=1)
    return (
        float(np.median(array)),
        float(np.quantile(estimates, alpha / 2)),
        float(np.quantile(estimates, 1 - alpha / 2)),
    )


def eta_squared(values: Iterable[float], groups: Iterable[str]) -> float:
    frame = pd.DataFrame({"value": values, "group": groups}).dropna()
    if len(frame) < 3 or frame.value.nunique() < 2 or frame.group.nunique() < 2:
        return float("nan")
    overall = float(frame.value.mean())
    between = sum(len(group) * (float(group.value.mean()) - overall) ** 2 for _, group in frame.groupby("group"))
    total = float(((frame.value - overall) ** 2).sum())
    return float(between / total) if total > 0 else float("nan")


def robust_standardize_scores(
    scores: pd.DataFrame,
    *,
    value_column: str = "score_native",
    group_columns: tuple[str, ...] = (
        "cohort_id",
        "modality",
        "aggregation_level",
        "cell_state_level",
        "cell_state",
        "feature_id",
    ),
) -> pd.DataFrame:
    result = scores.copy()
    groups = [column for column in group_columns if column in result]
    values = pd.to_numeric(result[value_column], errors="coerce")
    grouped_values = values.groupby(
        [result[column] for column in groups], dropna=False, sort=False
    )
    median = grouped_values.transform("median")
    absolute_deviation = (values - median).abs()
    mad = absolute_deviation.groupby(
        [result[column] for column in groups], dropna=False, sort=False
    ).transform("median")
    q75 = grouped_values.transform("quantile", q=0.75)
    q25 = grouped_values.transform("quantile", q=0.25)
    std = grouped_values.transform("std", ddof=0)
    scale = 1.4826 * mad
    invalid = ~np.isfinite(scale) | scale.eq(0)
    scale = scale.where(~invalid, (q75 - q25) / 1.349)
    invalid = ~np.isfinite(scale) | scale.eq(0)
    scale = scale.where(~invalid, std)
    invalid = ~np.isfinite(scale) | scale.eq(0)
    scale = scale.where(~invalid, 1.0)
    result["score_standardized"] = (values - median) / scale
    return result


def patient_timepoint_scores(expression_scores: pd.DataFrame) -> pd.DataFrame:
    required = {
        "cohort_id",
        "study_family",
        "patient_key",
        "timepoint",
        "cell_state_level",
        "cell_state",
        "feature_family",
        "feature_id",
        "score_native",
        "n_cells_used",
    }
    missing = sorted(required - set(expression_scores.columns))
    if missing:
        raise ValueError(f"patient aggregation missing columns: {missing}")
    source = expression_scores.loc[
        expression_scores.patient_key.fillna("").astype(str).str.strip().ne("")
    ].copy()
    keys = [
        "cohort_id",
        "study_family",
        "patient_key",
        "timepoint",
        "cell_state_level",
        "feature_family",
        "feature_id",
    ]
    source["_score"] = pd.to_numeric(source.score_native, errors="coerce")
    source["_weight"] = pd.to_numeric(source.n_cells_used, errors="coerce").fillna(0.0)
    source["_valid"] = np.isfinite(source._score)
    source["_valid_score"] = source._score.where(source._valid, 0.0)
    source["_valid_weight"] = source._weight.where(source._valid, 0.0)
    source["_weighted_score"] = source._valid_score * source._valid_weight
    source["_valid_state"] = source.cell_state.where(source._valid)
    source["_coverage"] = pd.to_numeric(source.get("coverage", np.nan), errors="coerce")
    grouped = (
        source.groupby(keys, dropna=False, sort=True)
        .agg(
            score_sum=("_valid_score", "sum"),
            weighted_score_sum=("_weighted_score", "sum"),
            valid_count=("_valid", "sum"),
            valid_weight=("_valid_weight", "sum"),
            n_states=("_valid_state", "nunique"),
            coverage=("_coverage", "median"),
        )
        .reset_index()
    )
    equal_score = np.divide(
        grouped.score_sum,
        grouped.valid_count,
        out=np.full(len(grouped), np.nan, dtype=float),
        where=grouped.valid_count.to_numpy() > 0,
    )
    weighted_score = np.divide(
        grouped.weighted_score_sum,
        grouped.valid_weight,
        out=np.asarray(equal_score, dtype=float).copy(),
        where=grouped.valid_weight.to_numpy() > 0,
    )
    frames: list[pd.DataFrame] = []
    for method, values in (
        ("abundance_weighted", weighted_score),
        ("equal_state", equal_score),
    ):
        frame = grouped[keys + ["n_states", "coverage"]].copy()
        frame["patient_timepoint_id"] = (
            frame.cohort_id.astype(str)
            + "::"
            + frame.patient_key.astype(str)
            + "::"
            + frame.timepoint.astype(str)
        )
        frame["aggregation_level"] = "patient_timepoint"
        frame["aggregation_method"] = method
        frame["modality"] = "scRNA"
        frame["cell_state"] = "ALL_STATES"
        frame["score_native"] = values
        frame["n_cells_used"] = grouped.valid_weight.astype(int)
        frame["technical_status"] = np.where(
            np.isfinite(values), "valid", "reject"
        )
        frames.append(frame)
    return pd.concat(frames, ignore_index=True).sort_values(
        keys + ["aggregation_method"]
    ).reset_index(drop=True)


def leave_family_out_state_profiles(scores: pd.DataFrame) -> pd.DataFrame:
    source = scores.loc[
        scores.technical_status.eq("valid")
        & scores.cell_state.notna()
        & ~scores.cell_state.isin(["Unknown", "Low_quality_or_ambient"])
    ].copy()
    profile = (
        source.groupby(["feature_id", "study_family", "cell_state"], as_index=False, dropna=False)
        .score_native.median()
    )
    profile = robust_standardize_scores(
        profile,
        value_column="score_native",
        group_columns=("feature_id", "study_family"),
    ).rename(columns={"score_standardized": "state_profile_score"})
    rows: list[dict[str, object]] = []
    for feature_id, feature in profile.groupby("feature_id", sort=True):
        for family, heldout in feature.groupby("study_family", sort=True):
            reference = (
                feature.loc[~feature.study_family.eq(family)]
                .groupby("cell_state", as_index=False)
                .state_profile_score.median()
                .rename(columns={"state_profile_score": "reference_score"})
            )
            paired = heldout.merge(reference, on="cell_state", how="inner")
            correlation = spearman_correlation(
                paired.state_profile_score, paired.reference_score
            )
            rows.append(
                {
                    "feature_id": feature_id,
                    "heldout_study_family": family,
                    "n_shared_states": len(paired),
                    "state_profile_spearman": correlation,
                    "diagnostic_status": (
                        "PASS"
                        if len(paired) >= 3 and np.isfinite(correlation)
                        else "NOT_TESTABLE"
                    ),
                }
            )
    return pd.DataFrame(rows)


def carrier_diagnostics(patient_scores: pd.DataFrame) -> pd.DataFrame:
    index = ["cohort_id", "patient_key", "timepoint", "cell_state_level", "feature_id"]
    pivot = patient_scores.pivot_table(index=index, columns="aggregation_method", values="score_native", aggfunc="first").reset_index()
    rows = []
    for feature_id, frame in pivot.groupby("feature_id", sort=True):
        correlation = spearman_correlation(frame.get("abundance_weighted"), frame.get("equal_state"))
        rows.append(
            {
                "feature_id": feature_id,
                "n_patient_timepoints": len(frame),
                "abundance_equal_state_spearman": correlation,
                "carrier_dependency_status": (
                    "NOT_TESTABLE"
                    if not np.isfinite(correlation)
                    else "DEPENDENT_OR_MIXED"
                    if correlation <= 0
                    else "CONCORDANT_WITH_RESIDUAL_COMPOSITION_RISK"
                ),
            }
        )
    return pd.DataFrame(rows)


def ambient_diagnostics(expression_scores: pd.DataFrame) -> pd.DataFrame:
    """Contrast the observed low-quality state with other states in matched contexts.

    This is a native-score sensitivity diagnostic.  It treats the annotated
    ``Low_quality_or_ambient`` state as an observed proxy; it does not estimate
    ambient contamination or a causal contamination effect.
    """

    required = {"feature_id", "cell_state", "score_native"}
    missing = sorted(required - set(expression_scores.columns))
    if missing:
        raise ValueError(f"ambient diagnostic missing columns: {missing}")
    context_candidates = (
        "cohort_id",
        "study_family",
        "object_id",
        "sample_key",
        "patient_key",
        "timepoint",
        "tissue_context",
        "platform",
        "modality",
        "aggregation_level",
        "cell_state_level",
    )
    context_columns = [column for column in context_candidates if column in expression_scores]
    source = expression_scores.copy()
    source["_native"] = pd.to_numeric(source.score_native, errors="coerce")
    source = source.loc[np.isfinite(source._native)].copy()
    if not context_columns:
        source["_matched_context"] = "ALL"
        context_columns = ["_matched_context"]
    source["_ambient"] = source._native.where(
        source.cell_state.eq("Low_quality_or_ambient")
    )
    source["_reference"] = source._native.where(
        ~source.cell_state.isin(["Low_quality_or_ambient", "Unknown"])
    )
    context = (
        source.groupby(["feature_id", *context_columns], dropna=False, sort=False)
        .agg(
            ambient_median=("_ambient", "median"),
            reference_median=("_reference", "median"),
            n_ambient=("_ambient", "count"),
            n_reference=("_reference", "count"),
        )
        .reset_index()
    )
    matched = context.loc[
        np.isfinite(context.ambient_median)
        & np.isfinite(context.reference_median)
    ].copy()
    matched["_difference"] = matched.ambient_median - matched.reference_median
    summary = (
        matched.groupby("feature_id", as_index=False, sort=True)
        .agg(
            n_matched_contexts=("_difference", "size"),
            n_ambient_units=("n_ambient", "sum"),
            n_reference_units=("n_reference", "sum"),
            ambient_minus_other_median_native=("_difference", "median"),
        )
    )
    result = pd.DataFrame(
        {"feature_id": sorted(expression_scores.feature_id.astype(str).unique())}
    ).merge(summary, on="feature_id", how="left", validate="one_to_one")
    for column in ("n_matched_contexts", "n_ambient_units", "n_reference_units"):
        result[column] = result[column].fillna(0).astype(int)
    # Compatibility field retained but deliberately not populated: the revised
    # diagnostic is not on a z-score scale.
    result["ambient_minus_other_median_z"] = np.nan
    result["ambient_metric_scale"] = "native_score"
    result["ambient_interpretation"] = (
        "observed_low_quality_state_proxy_not_contamination_causality"
    )
    result["ambient_status"] = np.where(
        result.n_matched_contexts.gt(0),
        "OBSERVED_LOW_QUALITY_STATE_PROXY_AVAILABLE",
        "NOT_TESTABLE_NO_MATCHED_CONTEXT",
    )
    return result


def confounding_diagnostics(expression_scores: pd.DataFrame) -> pd.DataFrame:
    """Estimate residual cohort/platform association without cohort centering.

    Native scores are first centered only by the global cell-state/level median.
    Residuals are then averaged once per patient-timepoint so each clinical unit
    receives equal weight in the cohort eta-squared calculation.
    """

    required = {
        "feature_id",
        "cohort_id",
        "cell_state_level",
        "cell_state",
        "score_native",
    }
    missing = sorted(required - set(expression_scores.columns))
    if missing:
        raise ValueError(f"confounding diagnostic missing columns: {missing}")

    rows: list[dict[str, object]] = []
    for feature_id, frame in expression_scores.groupby("feature_id", sort=True):
        source = frame.copy()
        source["_native"] = pd.to_numeric(source.score_native, errors="coerce")
        source = source.loc[np.isfinite(source._native)].copy()
        state_keys = ["cell_state_level", "cell_state"]
        state_median = source.groupby(state_keys, dropna=False)["_native"].transform("median")
        source["_state_residual"] = source._native - state_median

        if "patient_timepoint_id" in source:
            supplied = source.patient_timepoint_id.fillna("").astype(str).str.strip()
        else:
            supplied = pd.Series("", index=source.index, dtype=str)
        patient = (
            source.get("patient_key", pd.Series("", index=source.index))
            .fillna("")
            .astype(str)
            .str.strip()
        )
        timepoint = (
            source.get("timepoint", pd.Series("unknown", index=source.index))
            .fillna("unknown")
            .astype(str)
        )
        fallback_sample = None
        for candidate in ("sample_key", "object_id", "expression_unit_id"):
            if candidate in source:
                fallback_sample = source[candidate].fillna("").astype(str)
                break
        if fallback_sample is None:
            fallback_sample = pd.Series(source.index.astype(str), index=source.index)
        clinical_id = np.where(
            patient.ne(""),
            source.cohort_id.astype(str) + "::" + patient + "::" + timepoint,
            source.cohort_id.astype(str) + "::sample::" + fallback_sample,
        )
        source["_clinical_unit"] = supplied.where(supplied.ne(""), clinical_id)

        clinical_groups = source.groupby(
            ["_clinical_unit", "cohort_id"], dropna=False, sort=False
        )
        patient_timepoints = clinical_groups["_state_residual"].mean().reset_index()
        for column in ("platform", "tissue_context"):
            if column not in source:
                continue
            context = clinical_groups[column].agg(
                lambda values: (
                    str(values.dropna().iloc[0])
                    if values.dropna().astype(str).nunique() == 1
                    else "MULTIPLE_OR_MISSING"
                )
            )
            patient_timepoints = patient_timepoints.merge(
                context.rename(column).reset_index(),
                on=["_clinical_unit", "cohort_id"],
                how="left",
                validate="one_to_one",
            )
        n_cohorts = int(patient_timepoints.cohort_id.nunique(dropna=True))
        cohort_eta = eta_squared(
            patient_timepoints._state_residual, patient_timepoints.cohort_id
        )
        cohort_status = (
            "COMPUTED_PATIENT_TIMEPOINT_EQUAL_WEIGHT"
            if n_cohorts >= 2 and np.isfinite(cohort_eta)
            else "NOT_TESTABLE_SINGLE_COHORT"
            if n_cohorts < 2
            else "NOT_TESTABLE_NO_RESIDUAL_VARIATION"
        )

        if "platform" not in patient_timepoints:
            platform_eta = float("nan")
            platform_status = "NOT_TESTABLE_MISSING_PLATFORM"
            n_platforms = 0
        else:
            n_platforms = int(patient_timepoints.platform.nunique(dropna=True))
            platform_eta = eta_squared(
                patient_timepoints._state_residual, patient_timepoints.platform
            )
            platform_status = (
                "NOT_TESTABLE_SINGLE_PLATFORM"
                if n_platforms < 2
                else "COMPUTED_PATIENT_TIMEPOINT_EQUAL_WEIGHT"
                if np.isfinite(platform_eta)
                else "NOT_TESTABLE_NO_RESIDUAL_VARIATION"
            )
            if n_platforms < 2:
                platform_eta = float("nan")

        if "tissue_context" in patient_timepoints:
            tissue_eta = eta_squared(
                patient_timepoints._state_residual,
                patient_timepoints.tissue_context,
            )
        else:
            tissue_eta = float("nan")
        library = (
            source.library_size
            if "library_size" in source
            else pd.Series(np.nan, index=source.index)
        )
        n_cells = (
            source.n_cells_used
            if "n_cells_used" in source
            else pd.Series(np.nan, index=source.index)
        )
        rows.append(
            {
                "feature_id": feature_id,
                "n_patient_timepoints": int(patient_timepoints._clinical_unit.nunique()),
                "n_cohorts": n_cohorts,
                "n_platforms": n_platforms,
                "cohort_eta_squared": cohort_eta,
                "cohort_status": cohort_status,
                "platform_eta_squared": platform_eta,
                "platform_status": platform_status,
                "tissue_eta_squared_context_only": tissue_eta,
                "library_size_spearman": spearman_correlation(
                    source._state_residual, library
                ),
                "n_cells_spearman": spearman_correlation(
                    source._state_residual, n_cells
                ),
                "residualization": "global_cell_state_level_median_on_native_score",
                "status": "DIAGNOSTIC_NOT_AUTOMATICALLY_REGRESSED",
            }
        )
    return pd.DataFrame(rows)


def method_agreement_diagnostics(
    expression_scores: pd.DataFrame,
    *,
    seed: int = 20260831,
    bootstrap_replicates: int = 1000,
) -> pd.DataFrame:
    rows = []
    for feature_id, frame in expression_scores.groupby("feature_id", sort=True):
        by_family = [
            spearman_correlation(group.score_native, group.score_sensitivity)
            for _, group in frame.groupby("study_family", sort=True)
        ]
        estimate, lower, upper = bootstrap_median_interval(
            by_family,
            replicates=bootstrap_replicates,
            seed=seed,
        )
        rows.append(
            {
                "feature_id": feature_id,
                "method_agreement_median_spearman": estimate,
                "method_agreement_ci_low": lower,
                "method_agreement_ci_high": upper,
                "method_agreement_status": "PASS" if np.isfinite(estimate) else "NOT_TESTABLE",
            }
        )
    return pd.DataFrame(rows)


def technical_panel_reliability(
    panel: pd.DataFrame,
    *,
    seed: int = 20260831,
    bootstrap_replicates: int = 1000,
) -> pd.DataFrame:
    required = {"feature_id", "cell_mean_vs_pseudobulk_spearman"}
    missing = sorted(required - set(panel.columns))
    if missing:
        raise ValueError(f"technical panel missing columns: {missing}")
    rows: list[dict[str, object]] = []
    for feature_id, frame in panel.groupby("feature_id", sort=True):
        values = pd.to_numeric(
            frame.cell_mean_vs_pseudobulk_spearman, errors="coerce"
        )
        estimate, lower, upper = bootstrap_median_interval(
            values,
            replicates=bootstrap_replicates,
            seed=seed,
        )
        finite = values[np.isfinite(values)]
        rows.append(
            {
                "feature_id": feature_id,
                "technical_panel_median_spearman": estimate,
                "technical_panel_ci_low": lower,
                "technical_panel_ci_high": upper,
                "n_technical_panel_comparisons": int(len(finite)),
                "n_negative_technical_panel_comparisons": int((finite < 0).sum()),
                "technical_panel_status": (
                    "PASS_DIRECTIONAL_STABILITY"
                    if np.isfinite(lower) and lower > 0
                    else "MIXED_OR_NOT_TESTABLE"
                ),
            }
        )
    return pd.DataFrame(rows)


def summarize_reliability(
    expression_scores: pd.DataFrame,
    leave_family_out: pd.DataFrame,
    carrier: pd.DataFrame,
    ambient: pd.DataFrame,
    confounding: pd.DataFrame,
    method_agreement: pd.DataFrame,
    *,
    null_intervals: pd.DataFrame | None = None,
    technical_panel: pd.DataFrame | None = None,
    seed: int = 20260831,
    bootstrap_replicates: int = 1000,
) -> pd.DataFrame:
    rows = []
    null_lookup = (
        null_intervals.set_index("feature_id").to_dict("index")
        if null_intervals is not None and not null_intervals.empty
        else {}
    )
    technical_lookup = (
        technical_panel.set_index("feature_id").to_dict("index")
        if technical_panel is not None and not technical_panel.empty
        else {}
    )
    for feature_id, frame in expression_scores.groupby("feature_id", sort=True):
        lfo_values = leave_family_out.loc[
            leave_family_out.feature_id.eq(feature_id), "state_profile_spearman"
        ]
        estimate, lower, upper = bootstrap_median_interval(
            lfo_values,
            replicates=bootstrap_replicates,
            seed=seed,
        )
        null = null_lookup.get(feature_id, {})
        technical_valid = frame.technical_status.eq("valid").any() and frame.coverage.fillna(0).gt(0).any()
        null_upper = float(
            null.get(
                "matched_null_permutation_high",
                null.get("matched_null_ci_high", np.nan),
            )
        )
        coherence_lower = float(
            null.get(
                "observed_coherence_stability_low",
                null.get("observed_coherence_ci_low", np.nan),
            )
        )
        method_row = method_agreement.loc[method_agreement.feature_id.eq(feature_id)]
        method_lower = (
            float(method_row.method_agreement_ci_low.iloc[0]) if not method_row.empty else np.nan
        )
        panel = technical_lookup.get(feature_id, {})
        panel_lower = float(panel.get("technical_panel_ci_low", np.nan))
        elif_condition = (
            np.isfinite(lower)
            and lower > 0
            and np.isfinite(coherence_lower)
            and np.isfinite(null_upper)
            and coherence_lower > null_upper
            and np.isfinite(method_lower)
            and method_lower > 0
            and np.isfinite(panel_lower)
            and panel_lower > 0
        )
        if not technical_valid:
            technical_status = "reject"
            d2 = "annotation_only_or_reject"
            reason = "technical_invalid_or_zero_coverage"
        elif elif_condition:
            technical_status = "valid"
            d2 = "measurable"
            reason = "positive_study_family_outer_stability_with_technical_route_concordance_and_null_separation"
        else:
            technical_status = "valid"
            d2 = "joint_representation_with_uncertainty"
            reason = "scoreable_but_stability_or_null_evidence_is_incomplete_or_overlapping"
        row = {
            "feature_id": feature_id,
            "feature_family": frame.feature_family.iloc[0],
            "n_expression_units": int(frame.expression_unit_id.nunique()),
            "n_study_families": int(frame.study_family.nunique()),
            "median_coverage": float(pd.to_numeric(frame.coverage, errors="coerce").median()),
            "leave_family_out_median_spearman": estimate,
            "leave_family_out_ci_low": lower,
            "leave_family_out_ci_high": upper,
            "matched_null_ci_high": null_upper,
            "observed_coherence_ci_low": coherence_lower,
            "matched_null_permutation_high": null_upper,
            "observed_coherence_stability_low": coherence_lower,
            "technical_status": technical_status,
            "D2_class": d2,
            "D2_reason": reason,
        }
        rows.append(row)
    result = pd.DataFrame(rows)
    for diagnostic in (
        carrier,
        ambient,
        confounding,
        method_agreement,
        technical_panel,
    ):
        if diagnostic is None or diagnostic.empty:
            continue
        result = result.merge(diagnostic, on="feature_id", how="left", validate="one_to_one")
    return result


def _profile_score(
    matrix: np.ndarray,
    columns: list[int],
    directions: list[int],
    weights: list[float],
) -> np.ndarray:
    if not columns:
        return np.full(matrix.shape[0], np.nan)
    direction_array = np.asarray(directions, dtype=int)
    weight_array = np.asarray(weights, dtype=float)
    values = np.asarray(matrix[:, columns], dtype=float)

    def component(mask: np.ndarray) -> np.ndarray | None:
        if not mask.any():
            return None
        selected_weights = weight_array[mask]
        if not np.isfinite(selected_weights).all() or selected_weights.sum() <= 0:
            selected_weights = np.ones(mask.sum(), dtype=float)
        selected = values[:, mask]
        observed = np.isfinite(selected)
        numerator = np.where(observed, selected, 0.0) @ selected_weights
        denominator = observed @ selected_weights
        return np.divide(
            numerator,
            denominator,
            out=np.full(values.shape[0], np.nan, dtype=float),
            where=denominator > 0,
        )

    positive = component(direction_array > 0)
    negative = component(direction_array < 0)
    if positive is not None and negative is not None:
        return positive - negative
    if positive is not None:
        return positive
    if negative is not None:
        return -negative
    return np.full(matrix.shape[0], np.nan)


def matched_null_coherence(
    profile_matrix: np.ndarray,
    profile_symbols: list[str],
    membership: pd.DataFrame,
    dictionary: pd.DataFrame,
    *,
    replicates: int = 100,
    seed: int = 20260831,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Compare within-program split-half coherence with expression-matched nulls."""

    if profile_matrix.shape[1] != len(profile_symbols):
        raise ValueError("profile symbol count does not match matrix")
    finite = np.asarray(profile_matrix, dtype=float)
    observation_mask = np.isfinite(finite)
    observation_count = observation_mask.sum(axis=0)
    gene_mean = np.divide(
        np.where(observation_mask, finite, 0.0).sum(axis=0),
        observation_count,
        out=np.full(finite.shape[1], np.nan, dtype=float),
        where=observation_count > 0,
    )
    gene_detection = np.divide(
        ((finite > 0) & observation_mask).sum(axis=0),
        observation_count,
        out=np.full(finite.shape[1], np.nan, dtype=float),
        where=observation_count > 0,
    )
    stats = pd.DataFrame(
        {
            "gene_symbol": profile_symbols,
            "mean_logcpm": gene_mean,
            "detection_fraction": gene_detection,
            "n_finite_observations": observation_count,
            "observation_mask_key": [
                np.packbits(observation_mask[:, index]).tobytes()
                for index in range(finite.shape[1])
            ],
        }
    )
    for column, output in (("mean_logcpm", "mean_bin"), ("detection_fraction", "detection_bin")):
        ranks = stats[column].rank(method="first")
        stats[output] = pd.qcut(ranks, q=min(5, len(stats)), labels=False, duplicates="drop").fillna(0).astype(int)
    stats["match_bin"] = stats.mean_bin.astype(str) + "::" + stats.detection_bin.astype(str)
    position = {gene: index for index, gene in enumerate(profile_symbols)}
    by_mask = {
        key: group.index.to_numpy(dtype=np.int64)
        for key, group in stats.loc[stats.n_finite_observations.gt(0)].groupby(
            "observation_mask_key", sort=False
        )
    }
    scoreable = set(
        dictionary.loc[
            dictionary.scoreable.astype(str).str.lower().isin({"true", "1", "yes"}), "feature_id"
        ].astype(str)
    )
    rng = np.random.default_rng(seed)
    diagnostics: list[dict[str, object]] = []
    null_records: list[dict[str, object]] = []
    for feature_id, raw in membership.loc[membership.feature_id.astype(str).isin(scoreable)].groupby("feature_id", sort=True):
        feature_family = str(raw.feature_family.iloc[0])
        use = raw.copy()
        top_gene = (
            use.is_top_gene.astype(str)
            .str.lower()
            .isin({"yes", "true", "1"})
        )
        if feature_family == "legacy_fm" and top_gene.any():
            use = use.loc[top_gene]
        use = use.loc[use.gene_symbol.astype(str).isin(position)].drop_duplicates("gene_symbol")
        if len(use) < 4:
            diagnostics.append(
                {
                    "feature_id": feature_id,
                    "n_profile_genes": len(use),
                    "observed_coherence_median": np.nan,
                    "observed_coherence_ci_low": np.nan,
                    "observed_coherence_ci_high": np.nan,
                    "matched_null_median": np.nan,
                    "matched_null_ci_low": np.nan,
                    "matched_null_ci_high": np.nan,
                    "empirical_null_exceedance": np.nan,
                    "diagnostic_status": "NOT_TESTABLE_TOO_FEW_GENES",
                }
            )
            continue
        genes = use.gene_symbol.astype(str).tolist()
        directions = np.where(pd.to_numeric(use.direction, errors="coerce").fillna(1).to_numpy() < 0, -1, 1).tolist()
        raw_weights = pd.to_numeric(use.weight, errors="coerce").fillna(1.0).clip(lower=0).to_numpy(float)
        if not np.any(raw_weights > 0):
            raw_weights = np.ones(len(use), dtype=float)
        targets = np.asarray([position[gene] for gene in genes], dtype=np.int64)
        target_set = set(map(int, targets))
        target_masks = stats.loc[targets, "observation_mask_key"].tolist()
        insufficient_mask_pool = False
        for mask_key in sorted(set(target_masks)):
            n_targets_for_mask = sum(mask == mask_key for mask in target_masks)
            available = [
                int(index)
                for index in by_mask.get(mask_key, np.asarray([], dtype=np.int64))
                if int(index) not in target_set
            ]
            if len(available) < n_targets_for_mask:
                insufficient_mask_pool = True
                break
        if insufficient_mask_pool:
            diagnostics.append(
                {
                    "feature_id": feature_id,
                    "n_profile_genes": len(use),
                    "observed_coherence_median": np.nan,
                    "observed_coherence_stability_low": np.nan,
                    "observed_coherence_stability_high": np.nan,
                    "observed_coherence_ci_low": np.nan,
                    "observed_coherence_ci_high": np.nan,
                    "matched_null_median": np.nan,
                    "matched_null_permutation_low": np.nan,
                    "matched_null_permutation_high": np.nan,
                    "matched_null_ci_low": np.nan,
                    "matched_null_ci_high": np.nan,
                    "empirical_null_exceedance": np.nan,
                    "null_matching_method": "exact_observation_mask_nearest_expression_detection",
                    "interval_semantics": "algorithmic_stability_and_permutation_envelopes_not_biological_confidence_intervals",
                    "diagnostic_status": "NOT_TESTABLE_INSUFFICIENT_EXACT_MASK_NULLS",
                }
            )
            continue
        observed_values: list[float] = []
        null_values: list[float] = []
        for replicate in range(replicates):
            left_indices: list[int] = []
            right_indices: list[int] = []
            for direction in (-1, 1):
                members = np.flatnonzero(np.asarray(directions) == direction)
                if len(members) == 1:
                    left_indices.extend(members.tolist())
                    right_indices.extend(members.tolist())
                elif len(members) > 1:
                    shuffled = rng.permutation(members)
                    cut = max(1, len(shuffled) // 2)
                    left_indices.extend(shuffled[:cut].tolist())
                    right_indices.extend(shuffled[cut:].tolist())
            left = _profile_score(
                finite,
                [position[genes[index]] for index in left_indices],
                [directions[index] for index in left_indices],
                [raw_weights[index] for index in left_indices],
            )
            right = _profile_score(
                finite,
                [position[genes[index]] for index in right_indices],
                [directions[index] for index in right_indices],
                [raw_weights[index] for index in right_indices],
            )
            observed_values.append(spearman_correlation(left, right))

            matched_array = np.full(len(targets), -1, dtype=np.int64)
            used: set[int] = set(target_set)
            for mask_key in sorted(set(target_masks)):
                member_indices = np.asarray(
                    [index for index, value in enumerate(target_masks) if value == mask_key],
                    dtype=np.int64,
                )
                # Candidate eligibility is never widened: every null gene must
                # have exactly the same finite-row mask as its target gene.
                pool = np.asarray(by_mask[mask_key], dtype=np.int64)
                available = np.asarray(
                    [index for index in pool if int(index) not in used],
                    dtype=np.int64,
                )
                # Randomize target order to avoid systematic greedy priority;
                # within each target choose the nearest expression/detection
                # profile among the still-unused exact-mask candidates.
                for member_index in rng.permutation(member_indices):
                    target = int(targets[member_index])
                    candidates = np.asarray(
                        [index for index in available if int(index) not in used],
                        dtype=np.int64,
                    )
                    if not len(candidates):
                        # Static pool sufficiency is checked above.  This is a
                        # defensive invariant rather than a relaxed fallback.
                        raise RuntimeError(
                            f"exact-mask null assignment invariant failed for {feature_id}"
                        )
                    mean_delta = (
                        stats.loc[candidates, "mean_logcpm"].to_numpy(float)
                        - float(stats.loc[target, "mean_logcpm"])
                    )
                    detection_delta = (
                        stats.loc[candidates, "detection_fraction"].to_numpy(float)
                        - float(stats.loc[target, "detection_fraction"])
                    )
                    distance = mean_delta**2 + detection_delta**2
                    finite_distance = np.isfinite(distance)
                    if not finite_distance.any():
                        raise RuntimeError(
                            f"exact-mask null statistics invariant failed for {feature_id}"
                        )
                    best = candidates[finite_distance][
                        np.flatnonzero(
                            np.isclose(
                                distance[finite_distance],
                                np.nanmin(distance[finite_distance]),
                                rtol=0,
                                atol=1e-12,
                            )
                        )
                    ]
                    choice = int(rng.choice(best))
                    used.add(choice)
                    matched_array[member_index] = choice
            if np.any(matched_array < 0):
                raise RuntimeError(f"incomplete exact-mask null assignment for {feature_id}")
            matched = matched_array.tolist()
            null_left = _profile_score(
                finite,
                [matched[index] for index in left_indices],
                [directions[index] for index in left_indices],
                [raw_weights[index] for index in left_indices],
            )
            null_right = _profile_score(
                finite,
                [matched[index] for index in right_indices],
                [directions[index] for index in right_indices],
                [raw_weights[index] for index in right_indices],
            )
            value = spearman_correlation(null_left, null_right)
            null_values.append(value)
            null_records.append(
                {
                    "feature_id": feature_id,
                    "null_id": replicate,
                    "split_half_spearman": value,
                    "n_target_genes": len(targets),
                    "n_null_genes": len(matched),
                    "target_gene_symbols": "|".join(genes),
                    "null_gene_symbols": "|".join(
                        str(profile_symbols[index]) for index in matched
                    ),
                    "observation_mask_match": bool(
                        all(
                            np.array_equal(
                                observation_mask[:, target],
                                observation_mask[:, null],
                            )
                            for target, null in zip(targets, matched)
                        )
                    ),
                    "null_matching_method": "exact_observation_mask_nearest_expression_detection",
                }
            )
        observed_array = np.asarray(observed_values, dtype=float)
        observed_finite = observed_array[np.isfinite(observed_array)]
        null_array = np.asarray(null_values, dtype=float)
        null_finite = null_array[np.isfinite(null_array)]
        observed = (
            (
                float(np.median(observed_finite)),
                float(np.quantile(observed_finite, 0.025)),
                float(np.quantile(observed_finite, 0.975)),
            )
            if observed_finite.size
            else (np.nan, np.nan, np.nan)
        )
        null = (
            (
                float(np.median(null_finite)),
                float(np.quantile(null_finite, 0.025)),
                float(np.quantile(null_finite, 0.975)),
            )
            if null_finite.size
            else (np.nan, np.nan, np.nan)
        )
        diagnostic_status = (
            "PASS_COMPUTED"
            if observed_finite.size and null_finite.size
            else "NOT_TESTABLE_NO_PAIRWISE_COMPLETE_PROFILE"
        )
        diagnostics.append(
            {
                "feature_id": feature_id,
                "n_profile_genes": len(use),
                "observed_coherence_median": observed[0],
                "observed_coherence_stability_low": observed[1],
                "observed_coherence_stability_high": observed[2],
                # Compatibility aliases.  These are algorithm-repeat envelopes,
                # not confidence intervals for biological effects.
                "observed_coherence_ci_low": observed[1],
                "observed_coherence_ci_high": observed[2],
                "matched_null_median": null[0],
                "matched_null_permutation_low": null[1],
                "matched_null_permutation_high": null[2],
                "matched_null_ci_low": null[1],
                "matched_null_ci_high": null[2],
                "empirical_null_exceedance": (
                    float(np.mean(null_finite >= np.median(observed_finite)))
                    if observed_finite.size and null_finite.size
                    else np.nan
                ),
                "null_matching_method": "exact_observation_mask_nearest_expression_detection",
                "interval_semantics": "algorithmic_stability_and_permutation_envelopes_not_biological_confidence_intervals",
                "diagnostic_status": diagnostic_status,
            }
        )
    return pd.DataFrame(diagnostics), pd.DataFrame(null_records)


def gse207422_paired_descriptive_diagnostics(
    paired_scores: pd.DataFrame,
    *,
    value_column: str = "score_native",
) -> pd.DataFrame:
    """Summarize response-blind scRNA/bulk pairs already joined to the crosswalk.

    Only ``exact_matched``/``paired`` patient-timepoints enter the primary
    comparison.  Timepoint-mismatched rows remain counted as support-only and
    cannot influence the paired correlation.
    """

    required = {"feature_id", "modality", "match_status", value_column}
    missing = sorted(required - set(paired_scores.columns))
    if missing:
        raise ValueError(f"GSE207422 paired diagnostic missing columns: {missing}")
    pair_column = next(
        (
            column
            for column in ("pair_id", "patient_id", "patient_key")
            if column in paired_scores
        ),
        None,
    )
    if pair_column is None:
        raise ValueError(
            "GSE207422 paired diagnostic requires pair_id, patient_id, or patient_key"
        )
    source = paired_scores.copy()
    source[pair_column] = source[pair_column].astype(str)
    source["_value"] = pd.to_numeric(source[value_column], errors="coerce")
    role = (
        source.bridge_role.astype(str)
        if "bridge_role" in source
        else pd.Series(
            np.where(source.match_status.eq("exact_matched"), "paired", "support_only"),
            index=source.index,
        )
    )
    source["_primary"] = source.match_status.eq("exact_matched") & role.eq("paired")
    source["_support"] = ~source._primary

    rows: list[dict[str, object]] = []
    for feature_id, frame in source.groupby("feature_id", sort=True):
        primary = frame.loc[frame._primary]
        primary_pivot = primary.pivot_table(
            index=pair_column,
            columns="modality",
            values="_value",
            aggfunc="median",
        )
        if {"scRNA", "bulkRNA"}.issubset(primary_pivot.columns):
            primary_pivot = primary_pivot.loc[
                np.isfinite(primary_pivot.scRNA) & np.isfinite(primary_pivot.bulkRNA)
            ]
        else:
            primary_pivot = primary_pivot.iloc[0:0]
        support = frame.loc[frame._support]
        support_pivot = support.pivot_table(
            index=pair_column,
            columns="modality",
            values="_value",
            aggfunc="median",
        )
        if {"scRNA", "bulkRNA"}.issubset(support_pivot.columns):
            support_pivot = support_pivot.loc[
                np.isfinite(support_pivot.scRNA) & np.isfinite(support_pivot.bulkRNA)
            ]
        else:
            support_pivot = support_pivot.iloc[0:0]
        n_exact = int(len(primary_pivot))
        correlation = (
            spearman_correlation(primary_pivot.scRNA, primary_pivot.bulkRNA)
            if n_exact >= 3
            else float("nan")
        )
        status = (
            "paired_descriptive_available_non_directional"
            if n_exact >= 3 and np.isfinite(correlation)
            else "paired_descriptive_insufficient_n"
        )
        rows.append(
            {
                "feature_id": feature_id,
                "n_exact_pairs": n_exact,
                "n_support_only_pairs": int(len(support_pivot)),
                "primary_pair_ids": "|".join(
                    sorted(primary_pivot.index.astype(str).tolist())
                ),
                "support_only_pair_ids": "|".join(
                    sorted(support_pivot.index.astype(str).tolist())
                ),
                "paired_spearman": correlation,
                "paired_diagnostic_status": status,
                "paired_interpretation": "response_blind_descriptive_only_no_response_direction",
            }
        )
    return pd.DataFrame(rows)


def gse207422_pair_vector_diagnostics(
    paired_scores: pd.DataFrame,
    *,
    value_column: str = "score_native",
) -> pd.DataFrame:
    """Describe route-level feature-vector agreement without validating features."""

    required = {"feature_id", "modality", "match_status", value_column}
    missing = sorted(required - set(paired_scores.columns))
    if missing:
        raise ValueError(f"GSE207422 pair-vector diagnostic missing columns: {missing}")
    pair_column = next(
        (
            column
            for column in ("pair_id", "patient_id", "patient_key")
            if column in paired_scores
        ),
        None,
    )
    if pair_column is None:
        raise ValueError(
            "GSE207422 pair-vector diagnostic requires pair_id, patient_id, or patient_key"
        )
    source = paired_scores.copy()
    source[pair_column] = source[pair_column].astype(str)
    source["_value"] = pd.to_numeric(source[value_column], errors="coerce")
    source["_primary"] = source.match_status.eq("exact_matched")
    if "bridge_role" in source:
        source["_primary"] &= source.bridge_role.astype(str).eq("paired")
    rows: list[dict[str, object]] = []
    for pair_id, frame in source.groupby(pair_column, sort=True):
        pivot = frame.pivot_table(
            index="feature_id", columns="modality", values="_value", aggfunc="median"
        )
        if {"scRNA", "bulkRNA"}.issubset(pivot.columns):
            complete = np.isfinite(pivot.scRNA) & np.isfinite(pivot.bulkRNA)
            correlation = spearman_correlation(
                pivot.loc[complete, "scRNA"], pivot.loc[complete, "bulkRNA"]
            )
            n_features = int(complete.sum())
        else:
            correlation, n_features = float("nan"), 0
        primary = bool(frame._primary.all())
        rows.append(
            {
                "comparison_id": str(pair_id),
                "comparison_type": "exact_pair_feature_vector"
                if primary
                else "support_only_pair_feature_vector",
                "n_features": n_features,
                "feature_vector_spearman": correlation,
                "comparison_status": "PAIRED_DESCRIPTIVE"
                if primary and np.isfinite(correlation)
                else "SUPPORT_ONLY"
                if not primary
                else "NOT_TESTABLE",
                "interpretation": "route_level_descriptive_only_not_feature_validation",
            }
        )

    exact = source.loc[source._primary]
    exact_pivot = exact.pivot_table(
        index=[pair_column, "feature_id"],
        columns="modality",
        values="_value",
        aggfunc="median",
    )
    pair_ids = sorted(exact[pair_column].astype(str).unique())
    if len(pair_ids) == 2 and {"scRNA", "bulkRNA"}.issubset(exact_pivot.columns):
        left = exact_pivot.xs(pair_ids[0], level=pair_column)
        right = exact_pivot.xs(pair_ids[1], level=pair_column)
        shared = left.index.intersection(right.index)
        sc_delta = right.loc[shared, "scRNA"] - left.loc[shared, "scRNA"]
        bulk_delta = right.loc[shared, "bulkRNA"] - left.loc[shared, "bulkRNA"]
        complete = np.isfinite(sc_delta) & np.isfinite(bulk_delta)
        rows.append(
            {
                "comparison_id": f"{pair_ids[1]}-minus-{pair_ids[0]}",
                "comparison_type": "two_pair_feature_contrast",
                "n_features": int(complete.sum()),
                "feature_vector_spearman": spearman_correlation(
                    sc_delta.loc[complete], bulk_delta.loc[complete]
                ),
                "comparison_status": "PAIRED_DESCRIPTIVE_INSUFFICIENT_PATIENT_N",
                "interpretation": "route_level_descriptive_only_not_feature_validation",
            }
        )
    return pd.DataFrame(rows)


def gse193736_repeatability_diagnostics(
    scores: pd.DataFrame,
    design: pd.DataFrame | None = None,
    *,
    value_column: str = "score_native",
) -> pd.DataFrame:
    """Measure repeatability of condition profiles across GSE193736 replicates.

    The metric is the median pairwise Spearman correlation between replicate
    profiles over lineage x perturbation x rest/stim conditions.  It evaluates
    repeatability only and makes no perturbation or response-direction claim.
    """

    required_scores = {"feature_id", value_column}
    missing = sorted(required_scores - set(scores.columns))
    if missing:
        raise ValueError(f"GSE193736 repeatability missing score columns: {missing}")
    source = scores.copy()
    design_columns = {"lineage", "perturbation", "rest_stim", "replicate"}
    if not design_columns.issubset(source.columns):
        if design is None:
            missing_design = sorted(design_columns - set(source.columns))
            raise ValueError(
                f"GSE193736 repeatability missing design columns: {missing_design}"
            )
        score_sample = next(
            (column for column in ("sample_key", "sample_id") if column in source),
            None,
        )
        design_sample = next(
            (column for column in ("sample_id", "sample_key") if column in design),
            None,
        )
        if score_sample is None or design_sample is None:
            raise ValueError("GSE193736 repeatability requires sample_key or sample_id")
        design_required = design_columns | {design_sample}
        missing_design = sorted(design_required - set(design.columns))
        if missing_design:
            raise ValueError(
                f"GSE193736 repeatability missing design columns: {missing_design}"
            )
        source = source.merge(
            design[list(design_required)].rename(columns={design_sample: score_sample}),
            on=score_sample,
            how="left",
            validate="many_to_one",
        )
    source["_value"] = pd.to_numeric(source[value_column], errors="coerce")

    rows: list[dict[str, object]] = []
    condition_columns = ["lineage", "perturbation", "rest_stim"]
    for feature_id, frame in source.groupby("feature_id", sort=True):
        pivot = frame.pivot_table(
            index=condition_columns,
            columns="replicate",
            values="_value",
            aggfunc="median",
        ).sort_index()
        correlations: list[float] = []
        complete_counts: list[int] = []
        for left, right in combinations(pivot.columns.tolist(), 2):
            complete = np.isfinite(pivot[left]) & np.isfinite(pivot[right])
            complete_counts.append(int(complete.sum()))
            correlations.append(
                spearman_correlation(pivot.loc[complete, left], pivot.loc[complete, right])
            )
        finite = np.asarray(correlations, dtype=float)
        finite = finite[np.isfinite(finite)]
        repeatability = float(np.median(finite)) if finite.size else float("nan")
        rows.append(
            {
                "feature_id": feature_id,
                "n_conditions": int(len(pivot)),
                "n_replicates": int(len(pivot.columns)),
                "n_replicate_pairs": int(len(finite)),
                "minimum_pairwise_complete_conditions": (
                    min(complete_counts) if complete_counts else 0
                ),
                "replicate_profile_median_spearman": repeatability,
                "repeatability_metric": "median_pairwise_spearman_across_condition_profiles",
                "repeatability_status": (
                    "REPEATABILITY_COMPUTED_NO_DIRECTIONAL_CLAIM"
                    if finite.size
                    else "NOT_TESTABLE_REPEATABILITY"
                ),
            }
        )
    return pd.DataFrame(rows)


def cross_modal_mapping(
    scores: pd.DataFrame,
    reliability: pd.DataFrame,
    *,
    paired_bulk: pd.DataFrame | None = None,
    perturb_repeatability: pd.DataFrame | None = None,
) -> pd.DataFrame:
    grouped = (
        scores.groupby(
            ["feature_family", "feature_id", "modality", "aggregation_level"],
            as_index=False,
            dropna=False,
        )
        .agg(
            n_observations=("score_native", "size"),
            n_scoreable=("score_native", lambda values: int(pd.to_numeric(values, errors="coerce").notna().sum())),
            n_datasets=("cohort_id", "nunique"),
            n_patients=("patient_key", "nunique"),
            median_coverage=("coverage", "median"),
        )
    )
    mapped = grouped.merge(
        reliability[
            ["feature_id", "technical_status", "D2_class", "D2_reason"]
        ].rename(
            columns={
                "technical_status": "_scrna_technical_status",
                "D2_class": "_scrna_D2_class",
                "D2_reason": "_scrna_D2_reason",
            }
        ),
        on="feature_id",
        how="left",
        validate="many_to_one",
    )
    if paired_bulk is not None and not paired_bulk.empty:
        paired_columns = [
            column
            for column in (
                "feature_id",
                "n_exact_pairs",
                "n_support_only_pairs",
                "paired_spearman",
                "paired_diagnostic_status",
            )
            if column in paired_bulk
        ]
        mapped = mapped.merge(
            paired_bulk[paired_columns],
            on="feature_id",
            how="left",
            validate="many_to_one",
        )
    if perturb_repeatability is not None and not perturb_repeatability.empty:
        repeatability_columns = [
            column
            for column in (
                "feature_id",
                "replicate_profile_median_spearman",
                "n_replicate_pairs",
                "repeatability_status",
            )
            if column in perturb_repeatability
        ]
        mapped = mapped.merge(
            perturb_repeatability[repeatability_columns],
            on="feature_id",
            how="left",
            validate="many_to_one",
        )

    for column, default in (
        ("n_exact_pairs", 0),
        ("n_support_only_pairs", 0),
        ("paired_spearman", np.nan),
        ("paired_diagnostic_status", "paired_descriptive_insufficient_n"),
        ("replicate_profile_median_spearman", np.nan),
        ("n_replicate_pairs", 0),
        ("repeatability_status", "NOT_TESTABLE_REPEATABILITY"),
    ):
        if column not in mapped:
            mapped[column] = default

    scrna = mapped.modality.eq("scRNA")
    bulk = mapped.modality.eq("bulkRNA")
    perturb = mapped.modality.eq("bulk_perturbation")
    mapped["technical_status"] = "modality_specific_unclassified"
    mapped["D2_class"] = "modality_specific_evidence_unclassified"
    mapped["D2_reason"] = "no_modality_specific_evidence_contract"
    mapped["modality_evidence_status"] = "UNCLASSIFIED"

    mapped.loc[scrna, "technical_status"] = mapped.loc[
        scrna, "_scrna_technical_status"
    ].fillna("reject")
    mapped.loc[scrna, "D2_class"] = mapped.loc[scrna, "_scrna_D2_class"].fillna(
        "annotation_only_or_reject"
    )
    mapped.loc[scrna, "D2_reason"] = mapped.loc[scrna, "_scrna_D2_reason"].fillna(
        "missing_scRNA_reliability"
    )
    mapped.loc[scrna, "modality_evidence_status"] = "scRNA_D2_inherited"

    mapped.loc[bulk, "technical_status"] = "descriptive_only"
    mapped.loc[bulk, "D2_class"] = "paired_descriptive_insufficient_n"
    mapped.loc[bulk, "D2_reason"] = (
        "GSE207422_exact_pairs_are_response_blind_and_insufficient_for_directional_validation"
    )
    mapped.loc[bulk, "modality_evidence_status"] = "paired_descriptive_insufficient_n"

    mapped.loc[perturb, "technical_status"] = "interface_only"
    mapped.loc[perturb, "D2_class"] = "interface_only_no_directional_validation"
    mapped.loc[perturb, "D2_reason"] = (
        "GSE193736_repeatability_does_not_supply_perturbation_or_response_direction"
    )
    mapped.loc[perturb, "modality_evidence_status"] = (
        "interface_only_no_directional_validation"
    )

    mapped["n_paired_bulk_evidence"] = np.where(
        bulk, pd.to_numeric(mapped.n_exact_pairs, errors="coerce").fillna(0), 0
    ).astype(int)
    mapped["perturbation_repeatability"] = np.where(
        perturb,
        pd.to_numeric(
            mapped.replicate_profile_median_spearman, errors="coerce"
        ),
        np.nan,
    )
    return mapped.drop(
        columns=[
            "_scrna_technical_status",
            "_scrna_D2_class",
            "_scrna_D2_reason",
        ]
    )
