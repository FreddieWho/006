"""Fail-closed reader for a Visium 10x H5 matrix and tissue positions."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import h5py
import numpy as np
import pandas as pd
from scipy import sparse

from ..contracts import SpatialContractError, SpatialIdentity, SpatialUnitData, assert_model_safe_fields
from .base import build_feature_table, build_observation_table, make_spatial_unit


_POSITION_COLUMNS = (
    "barcode",
    "in_tissue",
    "array_row",
    "array_col",
    "pxl_row_in_fullres",
    "pxl_col_in_fullres",
)


def _decode(values: Any) -> list[str]:
    array = np.asarray(values)
    return [value.decode("utf-8") if isinstance(value, bytes) else str(value) for value in array]


def _read_10x_h5(path: Path) -> tuple[sparse.csr_matrix, list[str], list[str], list[str], list[str]]:
    with h5py.File(path, "r") as handle:
        names: list[str] = []
        handle.visit(names.append)
        assert_model_safe_fields(names, "10x H5 fields")
        if "matrix" not in handle:
            raise SpatialContractError("10x H5 is missing the matrix group")
        group = handle["matrix"]
        required = ("data", "indices", "indptr", "shape", "barcodes")
        missing = [name for name in required if name not in group]
        if missing:
            raise SpatialContractError(f"10x H5 matrix group is missing: {missing}")
        shape = tuple(int(value) for value in np.asarray(group["shape"][:]))
        if len(shape) != 2:
            raise SpatialContractError("10x H5 matrix shape must have two dimensions")
        gene_by_observation = sparse.csc_matrix(
            (
                np.asarray(group["data"][:]),
                np.asarray(group["indices"][:]),
                np.asarray(group["indptr"][:]),
            ),
            shape=shape,
        )
        counts = gene_by_observation.T.tocsr()
        barcodes = _decode(group["barcodes"][:])
        if len(barcodes) != counts.shape[0]:
            raise SpatialContractError("10x H5 barcode count does not match matrix shape")

        features = group.get("features")
        if features is not None:
            if "id" not in features or "name" not in features:
                raise SpatialContractError("10x H5 features group requires id and name")
            feature_ids = _decode(features["id"][:])
            gene_symbols = _decode(features["name"][:])
            feature_types = (
                _decode(features["feature_type"][:])
                if "feature_type" in features
                else ["Gene Expression"] * len(feature_ids)
            )
        elif "genes" in group and "gene_names" in group:
            feature_ids = _decode(group["genes"][:])
            gene_symbols = _decode(group["gene_names"][:])
            feature_types = ["Gene Expression"] * len(feature_ids)
        else:
            raise SpatialContractError("10x H5 is missing feature identifiers and names")
        if not (len(feature_ids) == len(gene_symbols) == len(feature_types) == counts.shape[1]):
            raise SpatialContractError("10x H5 feature annotations do not match matrix shape")
    return counts, barcodes, feature_ids, gene_symbols, feature_types


def _read_positions(path: Path) -> pd.DataFrame:
    header = pd.read_csv(path, nrows=0)
    if "barcode" in header.columns:
        assert_model_safe_fields(header.columns, "Visium positions")
        missing = [column for column in _POSITION_COLUMNS if column not in header.columns]
        if missing:
            raise SpatialContractError(f"Visium positions are missing columns: {missing}")
        positions = pd.read_csv(path, usecols=list(_POSITION_COLUMNS))
    else:
        first_row = pd.read_csv(path, header=None, nrows=1)
        if first_row.shape[1] != len(_POSITION_COLUMNS):
            raise SpatialContractError(
                "headerless Visium positions must contain exactly the six Space Ranger columns"
            )
        positions = pd.read_csv(path, header=None, names=list(_POSITION_COLUMNS))
    positions["barcode"] = positions["barcode"].astype("string").str.strip()
    if positions["barcode"].isna().any() or positions["barcode"].eq("").any():
        raise SpatialContractError("Visium positions contain missing barcodes")
    if positions["barcode"].duplicated().any():
        raise SpatialContractError("Visium positions contain duplicate barcodes")
    for column in _POSITION_COLUMNS[1:]:
        positions[column] = pd.to_numeric(positions[column], errors="coerce")
    if positions.loc[:, list(_POSITION_COLUMNS[1:])].isna().any().any():
        raise SpatialContractError("Visium positions contain non-numeric coordinates or masks")
    return positions


def load_visium_10x(
    h5_path: str | Path,
    positions_path: str | Path,
    identity: SpatialIdentity,
    *,
    scalefactors_path: str | Path | None = None,
) -> SpatialUnitData:
    """Load one Visium capture, preserving barcode, array and pixel coordinates."""

    matrix_path = Path(h5_path)
    position_path = Path(positions_path)
    if not matrix_path.is_file() or not position_path.is_file():
        raise SpatialContractError("Visium H5 and positions files must both be readable")
    counts, barcodes, feature_ids, gene_symbols, feature_types = _read_10x_h5(matrix_path)
    if len(barcodes) != len(set(barcodes)):
        raise SpatialContractError("10x H5 contains duplicate barcodes")
    positions = _read_positions(position_path).set_index("barcode", drop=False)
    missing = sorted(set(barcodes) - set(positions.index))
    if missing:
        raise SpatialContractError(
            f"Visium positions are missing {len(missing)} matrix barcodes; examples={missing[:5]}"
        )
    aligned = positions.loc[barcodes].reset_index(drop=True)
    coordinates = aligned.loc[:, ["pxl_col_in_fullres", "pxl_row_in_fullres"]].to_numpy()
    observations = build_observation_table(
        native_ids=barcodes,
        coordinates=coordinates,
        identity=identity,
        native_coordinate_unit="fullres_pixel",
        in_tissue=aligned["in_tissue"],
    )
    observations["array_row"] = aligned["array_row"].to_numpy()
    observations["array_col"] = aligned["array_col"].to_numpy()
    observations["pixel_row_fullres"] = aligned["pxl_row_in_fullres"].to_numpy()
    observations["pixel_col_fullres"] = aligned["pxl_col_in_fullres"].to_numpy()
    features = build_feature_table(
        native_feature_ids=feature_ids,
        gene_symbols=gene_symbols,
        feature_types=feature_types,
    )

    source_assets: list[str | Path] = [matrix_path, position_path]
    scale_audit: dict[str, Any] = {"scale_status": "not_provided"}
    if scalefactors_path is not None:
        scale_path = Path(scalefactors_path)
        if not scale_path.is_file():
            raise SpatialContractError(f"Visium scalefactors file is not readable: {scale_path}")
        with scale_path.open(encoding="utf-8") as handle:
            scale_factors = json.load(handle)
        if not isinstance(scale_factors, dict):
            raise SpatialContractError("Visium scalefactors must be a JSON object")
        assert_model_safe_fields(scale_factors.keys(), "Visium scalefactors")
        invalid = {
            key: value
            for key, value in scale_factors.items()
            if not isinstance(value, (int, float)) or not np.isfinite(value)
        }
        if invalid:
            raise SpatialContractError(
                f"Visium scalefactors contain non-finite/non-numeric entries: {sorted(invalid)}"
            )
        source_assets.append(scale_path)
        scale_audit = {
            "scale_status": "pixel_scalefactors_available_no_micron_claim",
            "scalefactors": dict(sorted(scale_factors.items())),
        }

    return make_spatial_unit(
        counts=counts,
        observations=observations,
        features=features,
        identity=identity,
        platform="Visium",
        resolution="spot",
        coordinate_system="fullres_pixel",
        source_assets=source_assets,
        audit={
            "counts_source": "10x_h5_raw_integer_matrix",
            "coordinate_source": position_path.name,
            "position_join_status": "PASS_ONE_TO_ONE_FOR_MATRIX_BARCODES",
            "n_positions_extra": int(len(positions) - len(barcodes)),
            "storage_dtype": str(counts.dtype),
            "image_status": "not_loaded",
            "segmentation_status": "not_applicable",
            **scale_audit,
        },
    )
