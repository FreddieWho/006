from __future__ import annotations

import h5py

from v7.spatial_io.sanitization import sanitize_h5ad_for_discovery


def test_sanitizer_removes_only_forbidden_uns_branches(tmp_path) -> None:
    source = tmp_path / "source.h5ad"
    destination = tmp_path / "safe.h5ad"
    with h5py.File(source, "w") as handle:
        handle.create_dataset("X", data=[1, 2, 3])
        uns = handle.create_group("uns")
        uns.create_dataset("response", data=[1])
        uns.create_dataset("response_provenance", data=[2])
        uns.create_dataset("technical_note", data=[3])

    audit = sanitize_h5ad_for_discovery(source, destination)
    assert audit["status"] == "SANITIZED_UNS_ONLY"
    assert audit["removed_fields"] == ["uns/response", "uns/response_provenance"]
    with h5py.File(destination, "r") as handle:
        names: list[str] = []
        handle.visit(names.append)
    assert "uns/technical_note" in names
    assert "uns/response" not in names
    assert "uns/response_provenance" not in names
