"""Patient-blocked folds for leakage-safe spatial model evaluation."""

from __future__ import annotations

from collections.abc import Sequence

import numpy as np
import pandas as pd


def patient_block_folds(
    patient_ids: np.ndarray | Sequence[object],
    *,
    n_splits: int = 5,
    seed: int = 20260904,
    strata: np.ndarray | Sequence[object] | None = None,
) -> np.ndarray:
    """Assign every observation from one patient to the same deterministic fold."""

    patients = np.asarray(patient_ids)
    if patients.ndim != 1 or len(patients) == 0:
        raise ValueError("patient_ids must be a non-empty one-dimensional array")
    if not isinstance(n_splits, (int, np.integer)) or n_splits < 2:
        raise ValueError("n_splits must be an integer >= 2")
    if any(value is None or str(value).strip() == "" for value in patients):
        raise ValueError("patient_ids contain missing or empty values")
    if strata is None:
        stratum_values = np.repeat("ALL", len(patients))
    else:
        stratum_values = np.asarray(strata)
        if stratum_values.ndim == 1:
            if len(stratum_values) != len(patients):
                raise ValueError("strata and patient_ids must be aligned")
            stratum_values = np.asarray([str(value) for value in stratum_values])
        elif stratum_values.ndim == 2 and len(stratum_values) == len(patients):
            stratum_values = np.asarray(
                ["\x1f".join(map(str, row)) for row in stratum_values]
            )
        else:
            raise ValueError("strata must be an aligned one- or two-dimensional array")

    patient_table = pd.DataFrame({"patient": patients.astype(str), "stratum": stratum_values})
    consistency = patient_table.groupby("patient").stratum.nunique(dropna=False)
    if consistency.gt(1).any():
        bad = consistency.index[consistency.gt(1)].tolist()
        raise ValueError(f"patients map to multiple strata: {bad[:5]}")
    unique = patient_table.drop_duplicates("patient").reset_index(drop=True)
    rng = np.random.default_rng(seed)
    fold_counts = np.zeros(n_splits, dtype=int)
    assignment: dict[str, int] = {}
    for _, group in unique.groupby("stratum", dropna=False, sort=True):
        patient_group = group.patient.to_numpy(copy=True)
        rng.shuffle(patient_group)
        for patient in patient_group:
            candidates = np.flatnonzero(fold_counts == fold_counts.min())
            fold = int(rng.choice(candidates))
            assignment[str(patient)] = fold
            fold_counts[fold] += 1
    return np.asarray([assignment[str(patient)] for patient in patients], dtype=int)


def patient_block_split(
    frame: pd.DataFrame,
    *,
    patient_column: str = "patient_id",
    n_splits: int = 5,
    seed: int = 20260904,
    stratify_columns: Sequence[str] = (),
    fold_column: str = "fold",
) -> pd.DataFrame:
    """Return a copy with patient-blocked fold assignments."""

    required = {patient_column, *stratify_columns}
    missing = sorted(required - set(frame.columns))
    if missing:
        raise ValueError(f"patient split missing columns: {missing}")
    strata = (
        frame[list(stratify_columns)].astype(str).to_numpy()
        if stratify_columns
        else None
    )
    result = frame.copy()
    result[fold_column] = patient_block_folds(
        result[patient_column].to_numpy(),
        n_splits=n_splits,
        seed=seed,
        strata=strata,
    )
    return result
