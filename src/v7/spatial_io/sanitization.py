"""Create an auditable, response-blind h5ad derivative for discovery.

Some upstream conversion files carry outcome annotations in ``uns`` even when
the matrix and coordinates are otherwise suitable.  The discovery adapters
must continue to fail closed on such files.  This helper therefore makes an
explicit HDF5-level copy that removes only forbidden ``uns`` branches, without
materialising their values, and records the source/destination hashes.  Any
forbidden field outside ``uns`` remains a hard block.
"""

from __future__ import annotations

import hashlib
from pathlib import Path
from typing import Any

import h5py

from .contracts import SpatialContractError, forbidden_field_names


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _copy_group_filtered(
    source: h5py.Group,
    destination: h5py.Group,
    prefix: str,
    removed: list[str],
) -> None:
    """Copy a group while omitting forbidden paths below ``uns``."""

    for key, value in source.attrs.items():
        destination.attrs[key] = value
    for name in source:
        path = f"{prefix}/{name}" if prefix else name
        if forbidden_field_names([path]):
            removed.append(path)
            continue
        source.copy(name, destination, name=name)


def sanitize_h5ad_for_discovery(
    source: str | Path,
    destination: str | Path,
) -> dict[str, Any]:
    """Copy an h5ad after removing only forbidden ``uns`` branches.

    The source is never modified.  Field names are inspected before any
    AnnData values are read; endpoint-bearing fields in ``obs``, ``var``,
    ``obsm``, ``layers`` or other non-``uns`` namespaces are not silently
    dropped and instead raise ``SpatialContractError``.
    """

    source_path = Path(source).resolve()
    destination_path = Path(destination).resolve()
    if not source_path.is_file():
        raise SpatialContractError(f"h5ad source is not a readable file: {source_path}")
    destination_path.parent.mkdir(parents=True, exist_ok=True)

    with h5py.File(source_path, "r") as source_handle:
        names: list[str] = []
        source_handle.visit(names.append)
        forbidden = sorted(forbidden_field_names(names))
        outside_uns = [name for name in forbidden if not name.startswith("uns/")]
        if outside_uns:
            raise SpatialContractError(
                "cannot sanitize forbidden fields outside h5ad uns: "
                + ",".join(outside_uns)
            )
        if not forbidden:
            # Keep a normal copied artifact so downstream provenance is always
            # explicit, even when a future input no longer carries ``uns``
            # annotations.
            with h5py.File(destination_path, "w") as destination_handle:
                for key, value in source_handle.attrs.items():
                    destination_handle.attrs[key] = value
                for name in source_handle:
                    source_handle.copy(name, destination_handle, name=name)
            removed: list[str] = []
        else:
            with h5py.File(destination_path, "w") as destination_handle:
                for key, value in source_handle.attrs.items():
                    destination_handle.attrs[key] = value
                for name in source_handle:
                    if name == "uns":
                        destination_handle.create_group("uns")
                        _copy_group_filtered(
                            source_handle[name], destination_handle["uns"], "uns", []
                        )
                    else:
                        source_handle.copy(name, destination_handle, name=name)
            # Re-read names only; no values from removed branches are loaded.
            with h5py.File(destination_path, "r") as destination_handle:
                kept_names: list[str] = []
                destination_handle.visit(kept_names.append)
            removed = [name for name in forbidden if name not in kept_names]

    return {
        "source_locator": str(source_path),
        "source_sha256": _sha256(source_path),
        "destination_locator": str(destination_path),
        "destination_sha256": _sha256(destination_path),
        "removed_fields": removed,
        "status": "SANITIZED_UNS_ONLY" if removed else "COPIED_NO_FORBIDDEN_FIELDS",
    }


__all__ = ["sanitize_h5ad_for_discovery"]
