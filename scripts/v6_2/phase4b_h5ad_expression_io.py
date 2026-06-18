#!/usr/bin/env python3
"""Small h5ad expression reader for Phase4B pseudobulk patch.

This module intentionally avoids anndata. It reads the minimal AnnData/H5AD
layout needed by the v6.2 objects: dense datasets or CSR-like groups under X or
layers[counts].
"""
from __future__ import annotations

from collections import defaultdict
from pathlib import Path
from typing import Iterable

import h5py
import numpy as np
import pandas as pd
from scipy import sparse


def _decode(values: np.ndarray) -> np.ndarray:
    if values.dtype.kind in {"S", "O"}:
        return np.asarray([x.decode() if isinstance(x, bytes) else str(x) for x in values])
    return values.astype(str)


def read_axis_strings(handle: h5py.File, axis: str, column: str) -> np.ndarray:
    node = handle[axis][column]
    if isinstance(node, h5py.Dataset):
        return _decode(node[()])
    if isinstance(node, h5py.Group):
        if "codes" in node and "categories" in node:
            codes = node["codes"][()]
            cats = _decode(node["categories"][()])
            out = np.asarray(["" if c < 0 else cats[int(c)] for c in codes], dtype=object)
            return out
        if "values" in node:
            return _decode(node["values"][()])
    raise ValueError(f"Unsupported h5ad axis column layout: {axis}/{column}")


def read_var_names(path: Path) -> list[str]:
    with h5py.File(path, "r") as handle:
        return read_axis_strings(handle, "var", "_index").astype(str).tolist()


def matrix_node(handle: h5py.File, selected_layer: str):
    layer = str(selected_layer)
    if layer == "X":
        return handle["X"]
    if layer.startswith("layers[") and layer.endswith("]"):
        name = layer[len("layers[") : -1]
        return handle["layers"][name]
    if layer == "raw.X":
        return handle["raw"]["X"]
    raise ValueError(f"Unsupported selected_layer: {selected_layer}")


def matrix_shape(node, n_obs: int, n_vars: int) -> tuple[int, int]:
    if isinstance(node, h5py.Dataset):
        return tuple(node.shape)
    shape = node.attrs.get("shape")
    if shape is not None:
        return int(shape[0]), int(shape[1])
    return n_obs, n_vars


def read_matrix_chunk(node, row_start: int, row_end: int, selected_cols: np.ndarray, n_vars: int):
    if isinstance(node, h5py.Dataset):
        return np.asarray(node[row_start:row_end, selected_cols])
    required = {"data", "indices", "indptr"}
    if isinstance(node, h5py.Group) and required.issubset(set(node.keys())):
        indptr = node["indptr"][row_start : row_end + 1]
        start = int(indptr[0])
        stop = int(indptr[-1])
        data = node["data"][start:stop]
        indices = node["indices"][start:stop]
        rel_indptr = indptr - start
        mat = sparse.csr_matrix((data, indices, rel_indptr), shape=(row_end - row_start, n_vars))
        return mat[:, selected_cols]
    raise ValueError("Unsupported h5ad matrix layout")


def aggregate_sample_pseudobulk(
    object_path: Path,
    selected_layer: str,
    sample_column: str,
    sample_to_key: dict[str, str],
    gene_symbols: Iterable[str],
    chunk_rows: int = 8192,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Aggregate sample-level pseudobulk for selected genes.

    Returns:
      wide_counts: sample_key + one column per gene.
      sample_counts: sample_key, n_cells_used.
    """
    genes = list(dict.fromkeys(str(g) for g in gene_symbols))
    with h5py.File(object_path, "r") as handle:
        var_names = read_axis_strings(handle, "var", "_index").astype(str)
        gene_to_idx = {g: i for i, g in enumerate(var_names)}
        selected = [(g, gene_to_idx[g]) for g in genes if g in gene_to_idx]
        if not selected:
            return pd.DataFrame(columns=["sample_key"]), pd.DataFrame(columns=["sample_key", "n_cells_used"])
        selected_genes = [g for g, _ in selected]
        selected_cols = np.asarray([i for _, i in selected], dtype=int)
        sample_values = read_axis_strings(handle, "obs", sample_column).astype(str)
        sample_keys = np.asarray([sample_to_key.get(v, "") for v in sample_values], dtype=object)
        node = matrix_node(handle, selected_layer)
        n_obs, n_vars = matrix_shape(node, len(sample_values), len(var_names))
        sums: dict[str, np.ndarray] = defaultdict(lambda: np.zeros(len(selected_genes), dtype=np.float64))
        counts: dict[str, int] = defaultdict(int)
        for start in range(0, n_obs, chunk_rows):
            end = min(start + chunk_rows, n_obs)
            keys = sample_keys[start:end]
            valid = keys != ""
            if not valid.any():
                continue
            mat = read_matrix_chunk(node, start, end, selected_cols, n_vars)
            unique_keys = np.unique(keys[valid])
            for sample_key in unique_keys:
                mask = keys == sample_key
                block = mat[mask]
                summed = np.asarray(block.sum(axis=0)).ravel()
                sums[str(sample_key)] += summed
                counts[str(sample_key)] += int(mask.sum())
    wide = pd.DataFrame.from_dict(sums, orient="index", columns=selected_genes).reset_index()
    wide = wide.rename(columns={"index": "sample_key"})
    n_cells = pd.DataFrame(
        [{"sample_key": key, "n_cells_used": value} for key, value in counts.items()]
    )
    return wide, n_cells
