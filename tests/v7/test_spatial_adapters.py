"""Synthetic h5ad, Visium and Xenium adapter tests."""

import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import anndata as ad
import h5py
import numpy as np
import pandas as pd
from scipy import sparse

from src.v7.spatial_io import (
    SpatialContractError,
    SpatialIdentity,
    load_h5ad_counts,
    load_visium_10x,
    load_xenium_h5ad,
)


def identity(capture: str = "capture::001") -> SpatialIdentity:
    return SpatialIdentity(
        dataset_id="synthetic",
        opaque_patient_id="patient::001",
        opaque_block_id="block::001",
        opaque_section_id="section::001",
        capture_id=capture,
        leakage_group_id="patient::001",
        relationship_provenance="synthetic_fixture",
        relationship_confidence="high",
    )


def write_h5ad(path: Path, *, fractional: bool = False, forbidden: bool = False) -> None:
    values = np.array([[1.5 if fractional else 1, 0], [2, 3]], dtype=float)
    obs = pd.DataFrame(index=pd.Index(["AAAC-1", "TTTG-1"], name="cell_id"))
    obs["in_tissue"] = [1, 1]
    if forbidden:
        obs["response"] = ["sealed-a", "sealed-b"]
    var = pd.DataFrame(
        {"gene_symbol": ["CXCL13", "TCF7"]},
        index=pd.Index(["ENSG1", "ENSG2"], name="gene_id"),
    )
    adata = ad.AnnData(X=sparse.csr_matrix(values), obs=obs, var=var)
    adata.obsm["spatial"] = np.array([[10.0, 20.0], [30.0, 40.0]])
    adata.write_h5ad(path)


def write_10x_h5(path: Path) -> None:
    # 2 genes x 2 observations in the native 10x CSC orientation.
    matrix = sparse.csc_matrix(np.array([[1, 0], [2, 3]], dtype=np.int32))
    with h5py.File(path, "w") as handle:
        group = handle.create_group("matrix")
        group.create_dataset("data", data=matrix.data)
        group.create_dataset("indices", data=matrix.indices)
        group.create_dataset("indptr", data=matrix.indptr)
        group.create_dataset("shape", data=np.asarray(matrix.shape, dtype=np.int64))
        group.create_dataset("barcodes", data=np.asarray([b"AAAC-1", b"TTTG-1"]))
        features = group.create_group("features")
        features.create_dataset("id", data=np.asarray([b"ENSG1", b"ENSG2"]))
        features.create_dataset("name", data=np.asarray([b"CXCL13", b"TCF7"]))
        features.create_dataset(
            "feature_type", data=np.asarray([b"Gene Expression", b"Gene Expression"])
        )


class SpatialAdapterTest(unittest.TestCase):
    def test_h5ad_explicit_x_preserves_counts_ids_and_coordinates(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            path = Path(temp) / "spatial.h5ad"
            write_h5ad(path)
            unit = load_h5ad_counts(
                path,
                identity(),
                counts_layer="X",
                coordinate_unit="pixel",
                gene_symbol_column="gene_symbol",
            )
        self.assertEqual(unit.counts.toarray().tolist(), [[1, 0], [2, 3]])
        self.assertEqual(unit.observations["native_observation_id"].tolist(), ["AAAC-1", "TTTG-1"])
        self.assertEqual(unit.observations["observation_id"].tolist()[0], "capture::001::AAAC-1")
        self.assertEqual(unit.observations[["native_x", "native_y"]].values.tolist(), [[10.0, 20.0], [30.0, 40.0]])
        self.assertEqual(unit.audit["counts_source"], "X")
        self.assertEqual(unit.audit["count_total"], 6)

    def test_h5ad_rejects_normalized_counts_and_forbidden_fields(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            normalized = Path(temp) / "normalized.h5ad"
            write_h5ad(normalized, fractional=True)
            with self.assertRaisesRegex(SpatialContractError, "not raw integer"):
                load_h5ad_counts(
                    normalized,
                    identity(),
                    counts_layer="X",
                    coordinate_unit="pixel",
                )

            forbidden = Path(temp) / "forbidden.h5ad"
            write_h5ad(forbidden, forbidden=True)
            with patch("src.v7.spatial_io.adapters.h5ad_counts.ad.read_h5ad") as reader:
                with self.assertRaisesRegex(SpatialContractError, "forbidden"):
                    load_h5ad_counts(
                        forbidden,
                        identity(),
                        counts_layer="X",
                        coordinate_unit="pixel",
                    )
                reader.assert_not_called()

            valid = Path(temp) / "valid.h5ad"
            write_h5ad(valid)
            with self.assertRaisesRegex(SpatialContractError, "missing counts layer"):
                load_h5ad_counts(
                    valid,
                    identity(),
                    counts_layer="counts",
                    coordinate_unit="pixel",
                )

    def test_visium_h5_positions_join_preserves_all_native_coordinates(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            h5 = root / "filtered_feature_bc_matrix.h5"
            positions = root / "tissue_positions.csv"
            scales = root / "scalefactors_json.json"
            write_10x_h5(h5)
            pd.DataFrame(
                {
                    "barcode": ["TTTG-1", "EXTRA-1", "AAAC-1"],
                    "in_tissue": [1, 0, 1],
                    "array_row": [2, 9, 1],
                    "array_col": [4, 9, 3],
                    "pxl_row_in_fullres": [200, 900, 100],
                    "pxl_col_in_fullres": [400, 900, 300],
                }
            ).to_csv(positions, index=False)
            scales.write_text(json.dumps({"tissue_hires_scalef": 0.25}), encoding="utf-8")
            unit = load_visium_10x(h5, positions, identity(), scalefactors_path=scales)
        self.assertEqual(unit.counts.toarray().tolist(), [[1, 2], [0, 3]])
        self.assertEqual(unit.observations["native_observation_id"].tolist(), ["AAAC-1", "TTTG-1"])
        self.assertEqual(unit.observations["array_row"].tolist(), [1, 2])
        self.assertEqual(unit.observations["native_x"].tolist(), [300.0, 400.0])
        self.assertEqual(unit.observations["native_y"].tolist(), [100.0, 200.0])
        self.assertEqual(unit.audit["n_positions_extra"], 1)
        self.assertEqual(unit.audit["scale_status"], "pixel_scalefactors_available_no_micron_claim")

    def test_visium_missing_or_duplicate_barcode_fails_closed(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            h5 = root / "matrix.h5"
            positions = root / "positions.csv"
            write_10x_h5(h5)
            pd.DataFrame(
                {
                    "barcode": ["AAAC-1", "AAAC-1"],
                    "in_tissue": [1, 1],
                    "array_row": [1, 2],
                    "array_col": [1, 2],
                    "pxl_row_in_fullres": [10, 20],
                    "pxl_col_in_fullres": [10, 20],
                }
            ).to_csv(positions, index=False)
            with self.assertRaisesRegex(SpatialContractError, "duplicate barcodes"):
                load_visium_10x(h5, positions, identity())

    def test_xenium_uses_native_cell_id_as_segmentation_and_keeps_panel_only(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            path = Path(temp) / "xenium.h5ad"
            write_h5ad(path)
            unit = load_xenium_h5ad(
                path,
                identity("xenium::001"),
                counts_layer="X",
                gene_symbol_column="gene_symbol",
            )
        self.assertEqual(unit.platform, "Xenium")
        self.assertEqual(unit.resolution, "cell")
        self.assertEqual(unit.coordinate_system, "micrometer")
        self.assertEqual(unit.observations["segmentation_id"].tolist(), ["AAAC-1", "TTTG-1"])
        self.assertEqual(unit.features["gene_symbol"].tolist(), ["CXCL13", "TCF7"])
        self.assertEqual(unit.audit["unmeasured_genes"], "absent_not_zero")


if __name__ == "__main__":
    unittest.main()
