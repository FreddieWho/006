"""Contract tests for response-blind spatial inputs."""

import unittest

import numpy as np
import pandas as pd
from scipy import sparse

from src.v7.spatial_io import (
    RAW_INTEGER_COUNTS,
    SpatialContractError,
    SpatialIdentity,
    SpatialUnitData,
    assert_model_safe_fields,
    audit_raw_integer_counts,
    validate_spatial_unit,
)


def identity() -> SpatialIdentity:
    return SpatialIdentity(
        dataset_id="synthetic",
        opaque_patient_id="patient::001",
        opaque_block_id="block::001",
        opaque_section_id="section::001",
        capture_id="capture::001",
        leakage_group_id="patient::001",
        relationship_provenance="synthetic_fixture",
        relationship_confidence="high",
    )


def valid_unit() -> SpatialUnitData:
    ident = identity()
    return SpatialUnitData(
        counts=sparse.csr_matrix([[1, 0], [0, 2]]),
        observations=pd.DataFrame(
            {
                "observation_id": ["capture::001::AAAC-1", "capture::001::TTTG-1"],
                "native_observation_id": ["AAAC-1", "TTTG-1"],
                "capture_id": ident.capture_id,
                "section_id": ident.opaque_section_id,
                "native_x": [1.0, 2.0],
                "native_y": [3.0, 4.0],
                "native_coordinate_unit": "pixel",
                "analysis_x": [1.0, 2.0],
                "analysis_y": [3.0, 4.0],
                "analysis_coordinate_unit": "pixel",
                "in_tissue": [True, True],
                "segmentation_id": pd.Series([pd.NA, pd.NA], dtype="string"),
                "observation_qc": "PASS",
            }
        ),
        features=pd.DataFrame(
            {
                "feature_id": ["ENSG1", "ENSG2"],
                "native_feature_id": ["ENSG1", "ENSG2"],
                "gene_symbol": ["A", "B"],
                "feature_type": "Gene Expression",
                "is_measured": True,
            }
        ),
        identity=ident,
        platform="synthetic",
        modality="spatial_transcriptomics",
        resolution="spot",
        counts_semantics=RAW_INTEGER_COUNTS,
        coordinate_system="pixel",
        source_assets=("/tmp/synthetic.h5",),
        audit={"adapter_status": "PASS"},
    )


class SpatialContractTest(unittest.TestCase):
    def test_valid_contract_preserves_count_total(self) -> None:
        unit = validate_spatial_unit(valid_unit())
        self.assertTrue(sparse.isspmatrix_csr(unit.counts))
        self.assertEqual(unit.count_total, 3)
        self.assertEqual(unit.n_observations, 2)
        self.assertEqual(unit.n_features, 2)

    def test_raw_count_audit_rejects_normalized_negative_and_nonfinite(self) -> None:
        invalid = (
            np.array([[0.5, 1.0]]),
            np.array([[-1.0, 1.0]]),
            np.array([[np.nan, 1.0]]),
        )
        for matrix in invalid:
            with self.subTest(matrix=matrix):
                with self.assertRaises(SpatialContractError):
                    audit_raw_integer_counts(matrix, "fixture")

    def test_contract_rejects_alignment_and_coordinate_failures(self) -> None:
        wrong_shape = valid_unit()
        wrong_shape.counts = sparse.csr_matrix([[1, 0]])
        with self.assertRaisesRegex(SpatialContractError, "shape"):
            validate_spatial_unit(wrong_shape)

        degenerate = valid_unit()
        degenerate.observations.loc[:, ["analysis_x", "analysis_y"]] = [1.0, 1.0]
        with self.assertRaisesRegex(SpatialContractError, "degenerate"):
            validate_spatial_unit(degenerate)

    def test_forbidden_response_and_ground_truth_fields_fail_closed(self) -> None:
        for field in ("response", "clinical_outcome", "GT", "structure_gt", "tls_distance"):
            with self.subTest(field=field):
                with self.assertRaises(SpatialContractError):
                    assert_model_safe_fields(["barcode", field], "fixture")

        unit = valid_unit()
        unit.observations["response_label"] = ["x", "y"]
        with self.assertRaises(SpatialContractError):
            validate_spatial_unit(unit)

    def test_unmeasured_panel_genes_must_stay_absent(self) -> None:
        unit = valid_unit()
        unit.features.loc[1, "is_measured"] = False
        with self.assertRaisesRegex(SpatialContractError, "stay absent"):
            validate_spatial_unit(unit)


if __name__ == "__main__":
    unittest.main()
