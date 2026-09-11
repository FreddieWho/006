"""Patient-first aggregation for nested spatial observations."""

from __future__ import annotations

from collections.abc import Sequence

import numpy as np
import pandas as pd


def patient_nested_summary(
    section_statistics: pd.DataFrame,
    *,
    value_column: str = "estimate",
    status_column: str = "status",
    patient_column: str = "patient_id",
    section_column: str = "section_id",
    cohort_column: str = "cohort_id",
    group_columns: Sequence[str] = (),
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Average sections within patient, then patients within cohort equally."""

    required = {
        value_column,
        status_column,
        patient_column,
        section_column,
        cohort_column,
        *group_columns,
    }
    missing = sorted(required - set(section_statistics.columns))
    if missing:
        raise ValueError(f"section statistics missing columns: {missing}")
    source = section_statistics.copy()
    source["_value"] = pd.to_numeric(source[value_column], errors="coerce")
    source["_valid"] = source[status_column].eq("ESTIMABLE") & np.isfinite(source._value)
    patient_keys = [cohort_column, patient_column, *group_columns]
    patient_rows: list[dict[str, object]] = []
    for key, frame in source.groupby(patient_keys, dropna=False, sort=True):
        if not isinstance(key, tuple):
            key = (key,)
        row = dict(zip(patient_keys, key))
        valid = frame.loc[frame._valid, "_value"]
        row.update(
            {
                "estimate": float(valid.mean()) if len(valid) else np.nan,
                "section_sd": float(valid.std(ddof=1)) if len(valid) > 1 else np.nan,
                "n_sections": int(frame[section_column].nunique()),
                "n_estimable_sections": int(len(valid)),
                "status": (
                    "ESTIMABLE" if len(valid) else "NOT_ESTIMABLE_NO_VALID_SECTIONS"
                ),
            }
        )
        patient_rows.append(row)
    patient = pd.DataFrame(patient_rows)

    cohort_keys = [cohort_column, *group_columns]
    cohort_rows: list[dict[str, object]] = []
    for key, frame in patient.groupby(cohort_keys, dropna=False, sort=True):
        if not isinstance(key, tuple):
            key = (key,)
        row = dict(zip(cohort_keys, key))
        valid = pd.to_numeric(
            frame.loc[frame.status.eq("ESTIMABLE"), "estimate"], errors="coerce"
        ).dropna()
        standard_error = (
            float(valid.std(ddof=1) / np.sqrt(len(valid))) if len(valid) > 1 else np.nan
        )
        row.update(
            {
                "estimate": float(valid.mean()) if len(valid) else np.nan,
                "patient_sd": float(valid.std(ddof=1)) if len(valid) > 1 else np.nan,
                "patient_se": standard_error,
                "n_patients": int(frame[patient_column].nunique()),
                "n_estimable_patients": int(len(valid)),
                "status": (
                    "ESTIMABLE" if len(valid) else "NOT_ESTIMABLE_NO_VALID_PATIENTS"
                ),
            }
        )
        cohort_rows.append(row)
    return patient, pd.DataFrame(cohort_rows)
