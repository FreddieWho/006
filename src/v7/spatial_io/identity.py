"""Identity helpers that prevent cross-capture and cross-patient leakage."""

from __future__ import annotations

from collections import defaultdict
from typing import Iterable, Sequence

import pandas as pd

from .contracts import SpatialContractError, SpatialIdentity


def namespace_observation_ids(
    native_ids: Sequence[object] | pd.Index,
    capture_id: str,
) -> pd.DataFrame:
    """Namespace IDs while retaining every native barcode/cell suffix verbatim."""

    if not isinstance(capture_id, str) or not capture_id.strip():
        raise SpatialContractError("capture_id must be a non-empty string")
    native = pd.Series(list(native_ids), dtype="string")
    if native.isna().any() or native.eq("").any():
        raise SpatialContractError("native observation IDs contain missing or empty values")
    if native.duplicated().any():
        duplicated = native.loc[native.duplicated(False)].astype(str).tolist()
        raise SpatialContractError(f"native observation IDs are duplicated: {duplicated[:5]}")
    native_text = native.astype(str)
    return pd.DataFrame(
        {
            "observation_id": capture_id + "::" + native_text,
            "native_observation_id": native_text,
        }
    )


def validate_identity_collection(identities: Iterable[SpatialIdentity]) -> None:
    """Enforce one patient envelope for each section, capture and leakage group."""

    identities = list(identities)
    if not identities:
        raise SpatialContractError("identity collection is empty")
    captures: set[str] = set()
    section_patients: dict[tuple[str, str], set[str]] = defaultdict(set)
    leakage_patients: dict[tuple[str, str], set[str]] = defaultdict(set)
    patient_leakage: dict[tuple[str, str], set[str]] = defaultdict(set)
    for identity in identities:
        capture_key = f"{identity.dataset_id}::{identity.capture_id}"
        if capture_key in captures:
            raise SpatialContractError(f"duplicate capture identity: {capture_key}")
        captures.add(capture_key)
        section_patients[(identity.dataset_id, identity.opaque_section_id)].add(
            identity.opaque_patient_id
        )
        leakage_patients[(identity.dataset_id, identity.leakage_group_id)].add(
            identity.opaque_patient_id
        )
        patient_leakage[(identity.dataset_id, identity.opaque_patient_id)].add(
            identity.leakage_group_id
        )

    for label, mapping in (
        ("section", section_patients),
        ("leakage group", leakage_patients),
        ("patient", patient_leakage),
    ):
        conflicts = {key: sorted(values) for key, values in mapping.items() if len(values) != 1}
        if conflicts:
            raise SpatialContractError(f"{label} identity conflict: {conflicts}")
