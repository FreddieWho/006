from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
from scipy import sparse

from pathlib import Path

from src.v7.ontology.aggregate import (
    ObjectSpec,
    _aggregate_level,
    _annotation_partition,
    shard_input_fingerprint,
)
from src.v7.ontology.contracts import Stage2Error


def _identity() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "_row": [0, 1, 2, 3],
            "_resolved": [True, True, True, False],
            "analysis_unit_key": ["S2", "S1", "S1", "S3"],
            "patient_key": ["P2", "P1", "P1", "P3"],
            "timepoint": ["pre", "pre", "pre", "pre"],
            "tissue_context": ["tumor"] * 4,
            "coarse": ["T_NK", "Myeloid", "Myeloid", "Tumor_like"],
        }
    )


def test_raw_count_aggregation_is_sorted_and_conserves_eligible_counts() -> None:
    matrix = sparse.csr_matrix(
        np.asarray(
            [
                [1, 0, 2],
                [0, 3, 0],
                [4, 0, 5],
                [99, 99, 99],
            ],
            dtype=float,
        )
    )

    aggregated, units, audit = _aggregate_level(
        matrix, _identity(), "coarse", pd.Index(["A", "B", "C"]), chunk_size=2
    )

    assert units.analysis_unit_key.tolist() == ["S1", "S2"]
    assert units.n_cells_used.tolist() == [2, 1]
    assert np.array_equal(
        aggregated.toarray(), np.asarray([[4, 3, 5], [1, 0, 2]], dtype=float)
    )
    assert audit["n_eligible_cells"] == 3
    assert audit["eligible_cell_total"] == 15
    assert audit["pseudobulk_total"] == 15
    assert audit["counts_conserved"] is True


def test_raw_count_aggregation_rejects_negative_declared_counts() -> None:
    matrix = sparse.csr_matrix(np.asarray([[1, -1], [2, 3], [1, 1], [1, 1]]))

    with pytest.raises(Stage2Error, match="INVALID_COUNTS"):
        _aggregate_level(
            matrix, _identity(), "coarse", pd.Index(["A", "B"]), chunk_size=4
        )


def test_annotation_partition_allows_unique_case_only_cohort_alias(tmp_path: Path) -> None:
    annotation_root = tmp_path / "annotations"
    annotation_root.mkdir()
    frame = pd.DataFrame(
        {
            "source_object_id": ["lambrecht_hcc::object_1"],
            "source_cell_id": ["cell_1"],
            "analysis_unit_key": ["S1"],
            "study_subject_key": ["P1"],
            "patient_key": ["P1"],
            "patient_id": ["P1"],
            "normalized_timepoint": ["pre"],
            "timepoint": ["pre"],
            "tissue_context": ["tumor"],
            "assignment_status": ["resolved"],
            "harmonized_coarse_label": ["T_NK"],
            "harmonized_mid_label": ["CD8_T"],
            "harmonized_fine_label": ["Cytotoxic"],
            "original_annotation": ["CD8"],
            "mapping_confidence": ["high"],
            "low_quality_flag": [0],
            "ambient_rna_risk_flag": [0],
        }
    )
    frame.to_parquet(annotation_root / "cohort_id=lambrecht_hcc.parquet", index=False)
    spec = ObjectSpec(
        cohort_id="LAMBRECHT_HCC",
        object_id="lambrecht_hcc::object_1",
        h5ad=tmp_path / "unused.h5ad",
        count_layer="counts",
        identity_sidecar=tmp_path / "unused.parquet",
        annotation_path=annotation_root,
    )

    observed = _annotation_partition(spec)

    assert observed.source_cell_id.tolist() == ["cell_1"]


def test_shard_fingerprint_detects_same_size_source_and_semantic_changes(
    tmp_path: Path,
) -> None:
    h5ad = tmp_path / "object.h5ad"
    h5ad.write_bytes(b"AAAA")
    identity = tmp_path / "identity"
    identity.mkdir()
    (identity / "part.parquet").write_bytes(b"BBBB")
    annotation = tmp_path / "annotation"
    annotation.mkdir()
    (annotation / "cohort_id=C.parquet").write_bytes(b"CCCC")
    spec = ObjectSpec(
        cohort_id="C",
        object_id="C::1",
        h5ad=h5ad,
        count_layer="counts",
        identity_sidecar=identity,
        annotation_path=annotation,
    )
    first = shard_input_fingerprint(
        spec, study_family="F", semantic_context={"gene_map": "one"}
    )
    h5ad.write_bytes(b"DDDD")
    same_size_changed = shard_input_fingerprint(
        spec, study_family="F", semantic_context={"gene_map": "one"}
    )
    assert first["fingerprint_sha256"] != same_size_changed["fingerprint_sha256"]

    semantic_changed = shard_input_fingerprint(
        spec, study_family="F", semantic_context={"gene_map": "two"}
    )
    assert (
        same_size_changed["fingerprint_sha256"]
        != semantic_changed["fingerprint_sha256"]
    )
