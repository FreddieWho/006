from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from src.v7.spatial_stats import patient_block_folds, patient_block_split


def test_patient_block_folds_never_split_one_patient_across_sections() -> None:
    patient_ids = np.array(["P1", "P1", "P2", "P2", "P2", "P3", "P4"])

    folds = patient_block_folds(patient_ids, n_splits=3, seed=12)

    for patient in np.unique(patient_ids):
        assert len(np.unique(folds[patient_ids == patient])) == 1
    assert np.array_equal(folds, patient_block_folds(patient_ids, n_splits=3, seed=12))


def test_patient_block_split_preserves_rows_and_stratifies_at_patient_level() -> None:
    frame = pd.DataFrame(
        {
            "patient_id": np.repeat(["P1", "P2", "P3", "P4", "P5", "P6"], 2),
            "section_id": [f"S{i}" for i in range(12)],
            "cohort_id": np.repeat(["A", "A", "A", "B", "B", "B"], 2),
        }
    )

    split = patient_block_split(
        frame,
        patient_column="patient_id",
        n_splits=3,
        seed=5,
        stratify_columns=("cohort_id",),
    )

    assert len(split) == len(frame)
    assert split.groupby("patient_id").fold.nunique().eq(1).all()
    assert set(split.fold) == {0, 1, 2}


def test_patient_block_split_rejects_inconsistent_patient_strata() -> None:
    frame = pd.DataFrame(
        {"patient_id": ["P1", "P1"], "cohort_id": ["A", "B"]}
    )

    with pytest.raises(ValueError, match="multiple strata"):
        patient_block_split(
            frame,
            patient_column="patient_id",
            stratify_columns=("cohort_id",),
        )
