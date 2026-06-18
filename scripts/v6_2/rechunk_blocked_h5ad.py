#!/usr/bin/env python3
"""Rechunk Phase6-blocked h5ad CSR matrices without requiring anndata.

This script creates sidecar h5ad files. Original objects are read-only. The
selected raw-count layer is copied into X with small CSR chunks so downstream
chunked reads do not repeatedly decompress huge HDF5 chunks.
"""
from __future__ import annotations

import argparse
import hashlib
import os
from datetime import date
from pathlib import Path

import h5py
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
PHASE35 = ROOT / "results/v6_2/phase3_5_single_cell_processing_qc_gate"
PHASE6 = ROOT / "results/v6_2/phase6_response_blind_module_discovery"
OUT = ROOT / "results/v6_2/phase6_h5ad_rechunk_patch"
RECHUNK_DIR = ROOT / "data/processed/srt/rechunked_h5ad"

DATA_CHUNK = 1_000_000
INDPTR_CHUNK = 250_000
COPY_BLOCK = 50_000_000
TODAY = str(date.today())

METADATA_GROUPS = ("obs", "var", "uns", "obsm", "varm", "obsp", "varp")


def read_csv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, dtype=str, keep_default_na=False)


def write_csv(df: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, index=False)


def partial_sha256(path: Path, nbytes: int = 4_194_304) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        h.update(fh.read(nbytes))
    return h.hexdigest()


def safe_file_id(value: object) -> str:
    token = "".join(ch if ch.isalnum() else "_" for ch in str(value)).strip("_")
    return token or "unknown"


def matrix_node(handle: h5py.File, selected_layer: str):
    layer = str(selected_layer)
    if layer == "X":
        return handle["X"]
    if layer.startswith("layers[") and layer.endswith("]"):
        return handle["layers"][layer[len("layers[") : -1]]
    if layer == "raw.X":
        return handle["raw"]["X"]
    raise ValueError(f"Unsupported selected_layer: {selected_layer}")


def matrix_shape(node) -> tuple[int, int]:
    if isinstance(node, h5py.Dataset):
        return int(node.shape[0]), int(node.shape[1])
    shape = node.attrs.get("shape")
    if shape is None:
        raise ValueError("CSR matrix group missing shape attr")
    return int(shape[0]), int(shape[1])


def require_csr(node) -> None:
    if not isinstance(node, h5py.Group) or not {"data", "indices", "indptr"}.issubset(node.keys()):
        raise ValueError("Selected matrix is not CSR-like h5ad group")


def copy_root_attrs(src: h5py.File, dst: h5py.File) -> None:
    for key, value in src.attrs.items():
        dst.attrs[key] = value


def copy_metadata_groups(src: h5py.File, dst: h5py.File) -> list[str]:
    copied = []
    for name in METADATA_GROUPS:
        if name in src:
            src.copy(name, dst, name=name)
            copied.append(name)
    return copied


def create_csr_dataset(
    dst_x: h5py.Group,
    name: str,
    src_ds: h5py.Dataset,
    chunk_size: int,
    compression: str | None,
) -> h5py.Dataset:
    shape = src_ds.shape
    chunks = (max(1, min(int(shape[0]), int(chunk_size))),)
    kwargs = {"shape": shape, "dtype": src_ds.dtype, "chunks": chunks}
    if compression:
        kwargs["compression"] = compression
    return dst_x.create_dataset(name, **kwargs)


def copy_vector(src_ds: h5py.Dataset, dst_ds: h5py.Dataset, block_size: int) -> None:
    n = int(src_ds.shape[0])
    for start in range(0, n, block_size):
        end = min(start + block_size, n)
        dst_ds[start:end] = src_ds[start:end]


def rechunk_one(
    source_path: Path,
    selected_layer: str,
    target_path: Path,
    copy_block_size: int,
    compression: str | None,
) -> dict[str, object]:
    tmp_path = target_path.with_suffix(target_path.suffix + ".tmp")
    if tmp_path.exists():
        tmp_path.unlink()
    target_path.parent.mkdir(parents=True, exist_ok=True)

    with h5py.File(source_path, "r") as src, h5py.File(tmp_path, "w") as dst:
        node = matrix_node(src, selected_layer)
        require_csr(node)
        n_obs, n_vars = matrix_shape(node)
        copied_groups = copy_metadata_groups(src, dst)
        copy_root_attrs(src, dst)
        dst_x = dst.create_group("X")
        for key, value in node.attrs.items():
            dst_x.attrs[key] = value
        dst_x.attrs["encoding-type"] = node.attrs.get("encoding-type", "csr_matrix")
        dst_x.attrs["encoding-version"] = node.attrs.get("encoding-version", "0.1.0")
        dst_x.attrs["shape"] = np.asarray([n_obs, n_vars], dtype=np.int64)

        dst_data = create_csr_dataset(dst_x, "data", node["data"], DATA_CHUNK, compression)
        dst_indices = create_csr_dataset(dst_x, "indices", node["indices"], DATA_CHUNK, compression)
        dst_indptr = create_csr_dataset(dst_x, "indptr", node["indptr"], INDPTR_CHUNK, compression)
        copy_vector(node["data"], dst_data, copy_block_size)
        copy_vector(node["indices"], dst_indices, copy_block_size)
        copy_vector(node["indptr"], dst_indptr, copy_block_size)

    tmp_path.replace(target_path)
    return validate_rechunk(source_path, selected_layer, target_path) | {
        "metadata_groups_copied": ";".join(copied_groups),
        "compression": compression or "none",
    }


def validate_rechunk(source_path: Path, selected_layer: str, target_path: Path) -> dict[str, object]:
    with h5py.File(source_path, "r") as src, h5py.File(target_path, "r") as dst:
        src_node = matrix_node(src, selected_layer)
        dst_node = dst["X"]
        require_csr(src_node)
        require_csr(dst_node)
        src_shape = matrix_shape(src_node)
        dst_shape = matrix_shape(dst_node)
        src_nnz = int(src_node["data"].shape[0])
        dst_nnz = int(dst_node["data"].shape[0])
        src_indptr_n = int(src_node["indptr"].shape[0])
        dst_indptr_n = int(dst_node["indptr"].shape[0])
        data_sum_equal = bool(np.isclose(src_node["data"][:].sum(), dst_node["data"][:].sum()))
        first_last_equal = bool(
            src_node["indptr"][0] == dst_node["indptr"][0]
            and src_node["indptr"][-1] == dst_node["indptr"][-1]
        )
        data_chunk = int(dst_node["data"].chunks[0]) if dst_node["data"].chunks else 0
        indices_chunk = int(dst_node["indices"].chunks[0]) if dst_node["indices"].chunks else 0
        indptr_chunk = int(dst_node["indptr"].chunks[0]) if dst_node["indptr"].chunks else 0
    valid = (
        src_shape == dst_shape
        and src_nnz == dst_nnz
        and src_indptr_n == dst_indptr_n
        and data_sum_equal
        and first_last_equal
        and 0 < data_chunk <= DATA_CHUNK
        and 0 < indices_chunk <= DATA_CHUNK
        and 0 < indptr_chunk <= INDPTR_CHUNK
    )
    return {
        "validation_status": "pass" if valid else "fail",
        "source_shape": f"{src_shape[0]}x{src_shape[1]}",
        "target_shape": f"{dst_shape[0]}x{dst_shape[1]}",
        "source_nnz": src_nnz,
        "target_nnz": dst_nnz,
        "data_sum_equal": str(data_sum_equal).lower(),
        "indptr_boundary_equal": str(first_last_equal).lower(),
        "target_data_chunk": data_chunk,
        "target_indices_chunk": indices_chunk,
        "target_indptr_chunk": indptr_chunk,
    }


def load_blocked_plan() -> pd.DataFrame:
    qc_path = PHASE6 / "foundation/object_join_qc.csv"
    layer_path = PHASE35 / "matrix_layer_decision_table.frozen_v0.csv"
    registry_path = PHASE35 / "analysis_object_registry.frozen_v0.csv"
    qc_tables = []
    if qc_path.exists():
        qc_tables.append(read_csv(qc_path))
    parts_dir = PHASE6 / "foundation/object_parts"
    for part_qc in parts_dir.glob("*.qc.csv"):
        try:
            qc_tables.append(read_csv(part_qc))
        except Exception:
            pass
    qc = pd.concat(qc_tables, ignore_index=True, sort=False) if qc_tables else pd.DataFrame()
    layer = read_csv(layer_path)
    registry = read_csv(registry_path)
    if qc.empty:
        return pd.DataFrame(columns=["cohort_id", "object_id", "status"])
    blocked = qc[qc["status"].str.startswith("blocked_hdf5", na=False)][["cohort_id", "object_id", "status"]].drop_duplicates()
    keep = ["cohort_id", "object_id", "selected_layer", "raw_counts_allowed"]
    blocked = blocked.merge(layer[[c for c in keep if c in layer.columns]], on=["cohort_id", "object_id"], how="left")
    reg_keep = ["cohort_id", "object_id", "object_path", "n_cells", "n_genes"]
    blocked = blocked.merge(registry[[c for c in reg_keep if c in registry.columns]], on=["cohort_id", "object_id"], how="left")
    return blocked


def run(args: argparse.Namespace) -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    RECHUNK_DIR.mkdir(parents=True, exist_ok=True)
    plan = load_blocked_plan()
    if args.object_id:
        wanted = set(args.object_id)
        plan = plan[plan["object_id"].isin(wanted)].copy()
    plan["n_cells_num"] = pd.to_numeric(plan.get("n_cells", ""), errors="coerce").fillna(0).astype(int)
    plan = plan.sort_values(["n_cells_num", "cohort_id", "object_id"]).reset_index(drop=True)
    if args.limit:
        plan = plan.head(args.limit).copy()
    write_csv(plan, OUT / "rechunk_plan.csv")

    existing = read_csv(OUT / "rechunk_manifest.csv") if (OUT / "rechunk_manifest.csv").exists() else pd.DataFrame()
    rows = []
    if not existing.empty:
        rows.extend(existing.to_dict("records"))
    seen = {(r.get("cohort_id"), r.get("object_id")) for r in rows}

    for _, row in plan.iterrows():
        cohort_id = row["cohort_id"]
        object_id = row["object_id"]
        source_path = Path(row["object_path"])
        selected_layer = row["selected_layer"]
        target_path = RECHUNK_DIR / f"{safe_file_id(cohort_id)}__{safe_file_id(object_id)}.rechunked_X.h5ad"
        key = (cohort_id, object_id)
        if key in seen and not args.force:
            current = [r for r in rows if (r.get("cohort_id"), r.get("object_id")) == key][-1]
            if current.get("copy_status") == "complete" and Path(current.get("rechunked_h5ad", "")).exists():
                continue
        base = {
            "patch_id": "phase6_h5ad_rechunk_patch",
            "created_at": TODAY,
            "cohort_id": cohort_id,
            "object_id": object_id,
            "original_h5ad": str(source_path),
            "original_selected_layer": selected_layer,
            "rechunked_h5ad": str(target_path),
            "phase6_selected_layer_after_rechunk": "X",
            "copy_status": "",
            "source_block_reason": row.get("status", ""),
            "source_partial_sha256_4mb": partial_sha256(source_path) if source_path.exists() else "",
            "target_size_bytes": "",
            "error": "",
        }
        try:
            if target_path.exists() and not args.force:
                qc = validate_rechunk(source_path, selected_layer, target_path)
            else:
                qc = rechunk_one(
                    source_path=source_path,
                    selected_layer=selected_layer,
                    target_path=target_path,
                    copy_block_size=args.copy_block_size,
                    compression=args.compression,
                )
            base.update(qc)
            base["copy_status"] = "complete" if base.get("validation_status") == "pass" else "failed_validation"
            base["target_size_bytes"] = target_path.stat().st_size if target_path.exists() else ""
        except Exception as exc:  # noqa: BLE001 - manifest must capture per-object failure.
            base["copy_status"] = "failed"
            base["error"] = repr(exc)
        rows = [r for r in rows if (r.get("cohort_id"), r.get("object_id")) != key]
        rows.append(base)
        write_csv(pd.DataFrame(rows), OUT / "rechunk_manifest.csv")

    manifest = pd.DataFrame(rows)
    write_csv(manifest, OUT / "rechunk_manifest.csv")
    summary = (
        manifest.groupby("copy_status", dropna=False)
        .size()
        .reset_index(name="n_objects")
        .sort_values("copy_status")
        if not manifest.empty
        else pd.DataFrame(columns=["copy_status", "n_objects"])
    )
    write_csv(summary, OUT / "rechunk_status_summary.csv")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--object-id", action="append", default=[], help="Object id to rechunk; can be repeated.")
    parser.add_argument("--limit", type=int, default=0)
    parser.add_argument("--force", action="store_true")
    parser.add_argument("--copy-block-size", type=int, default=COPY_BLOCK)
    parser.add_argument("--compression", choices=["gzip", "lzf"], default=None)
    return parser.parse_args()


if __name__ == "__main__":
    run(parse_args())
