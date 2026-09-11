"""Patient/capture identity and native observation ID invariants."""

import unittest

from src.v7.spatial_io import (
    SpatialContractError,
    SpatialIdentity,
    namespace_observation_ids,
    validate_identity_collection,
)


def identity(
    *,
    patient: str,
    section: str,
    capture: str,
    leakage_group: str | None = None,
) -> SpatialIdentity:
    return SpatialIdentity(
        dataset_id="dataset",
        opaque_patient_id=patient,
        opaque_block_id=f"block::{patient}",
        opaque_section_id=section,
        capture_id=capture,
        leakage_group_id=leakage_group or patient,
        relationship_provenance="synthetic_fixture",
        relationship_confidence="high",
    )


class SpatialIdentityTest(unittest.TestCase):
    def test_namespace_preserves_native_barcode_suffix(self) -> None:
        table = namespace_observation_ids(["AAAC-1", "AAAC-2"], "capture::A")
        self.assertEqual(table["native_observation_id"].tolist(), ["AAAC-1", "AAAC-2"])
        self.assertEqual(
            table["observation_id"].tolist(),
            ["capture::A::AAAC-1", "capture::A::AAAC-2"],
        )

    def test_duplicate_native_observations_fail(self) -> None:
        with self.assertRaisesRegex(SpatialContractError, "duplicated"):
            namespace_observation_ids(["AAAC-1", "AAAC-1"], "capture::A")

    def test_multiple_captures_for_one_patient_share_outer_group(self) -> None:
        identities = [
            identity(patient="patient::1", section="section::1", capture="capture::1"),
            identity(patient="patient::1", section="section::2", capture="capture::2"),
            identity(patient="patient::2", section="section::3", capture="capture::3"),
        ]
        validate_identity_collection(identities)

    def test_section_or_leakage_group_cannot_span_patients(self) -> None:
        same_section = [
            identity(patient="patient::1", section="section::shared", capture="capture::1"),
            identity(patient="patient::2", section="section::shared", capture="capture::2"),
        ]
        with self.assertRaisesRegex(SpatialContractError, "section identity conflict"):
            validate_identity_collection(same_section)

        same_leakage = [
            identity(
                patient="patient::1",
                section="section::1",
                capture="capture::1",
                leakage_group="outer::shared",
            ),
            identity(
                patient="patient::2",
                section="section::2",
                capture="capture::2",
                leakage_group="outer::shared",
            ),
        ]
        with self.assertRaisesRegex(SpatialContractError, "leakage group identity conflict"):
            validate_identity_collection(same_leakage)

    def test_duplicate_capture_fails(self) -> None:
        duplicate = [
            identity(patient="patient::1", section="section::1", capture="capture::1"),
            identity(patient="patient::1", section="section::1", capture="capture::1"),
        ]
        with self.assertRaisesRegex(SpatialContractError, "duplicate capture"):
            validate_identity_collection(duplicate)


if __name__ == "__main__":
    unittest.main()
