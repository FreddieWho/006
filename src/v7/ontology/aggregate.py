"""Streaming raw-count pseudobulk aggregation for v7 Stage 2.

The aggregation layer is deliberately response blind.  It consumes only the
repaired identity and cell-state fields required to bind a cell to an
expression unit.  Outcome columns are neither requested nor propagated.
"""

from __future__ import annotations

import json
import hashlib
import os
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable

import numpy as np
import pandas as pd
import yaml
from scipy import sparse

from .contracts import Stage2Error, sha256, write_json, write_tsv


SAFE_IDENTITY_COLUMNS = (
    "source_object_id",
    "source_cell_id",
    "analysis_unit_key",
    "study_subject_key",
    "patient_key",
    "patient_id",
    "normalized_timepoint",
    "timepoint",
    "tissue_context",
    "tissue_source",
    "assignment_status",
)


@dataclass(frozen=True)
class ObjectSpec:
    cohort_id: str
    object_id: str
    h5ad: Path
    count_layer: str
    identity_sidecar: Path
    annotation_path: Path
    h5ad_cell_id_column: str | None = None
    h5ad_gene_column: str | None = None


@dataclass(frozen=True)
class AggregatedShard:
    cohort_id: str
    object_id: str
    level: str
    matrix_path: Path
    units_path: Path
    genes_path: Path
    audit_path: Path


def load_object_specs(projection_config: Path) -> list[ObjectSpec]:
    config = yaml.safe_load(projection_config.read_text(encoding="utf-8")) or {}
    inputs = config.get("inputs") or {}
    annotation_path = Path(inputs.get("phase4a_annotations", ""))
    if not annotation_path.exists():
        raise Stage2Error("STAGE2_BLOCKED_MISSING_ANNOTATION_MASTER")
    objects = inputs.get("objects") or []
    specs: list[ObjectSpec] = []
    for raw in objects:
        required = ("cohort_id", "h5ad", "count_layer", "identity_sidecar")
        if any(not raw.get(key) for key in required):
            raise Stage2Error("STAGE2_BLOCKED_OBJECT_CONFIG: incomplete object entry")
        spec = ObjectSpec(
            cohort_id=str(raw["cohort_id"]),
            object_id=str(raw.get("object_id", raw["cohort_id"])),
            h5ad=Path(raw["h5ad"]),
            count_layer=str(raw["count_layer"]),
            identity_sidecar=Path(raw["identity_sidecar"]),
            annotation_path=annotation_path,
            h5ad_cell_id_column=raw.get("h5ad_cell_id_column"),
            h5ad_gene_column=raw.get("h5ad_gene_column"),
        )
        if spec.count_layer not in {"X", "counts"}:
            raise Stage2Error(f"STAGE2_BLOCKED_LAYER: {spec.object_id}={spec.count_layer}")
        for path in (spec.h5ad, spec.identity_sidecar):
            if not path.exists():
                raise Stage2Error(f"STAGE2_BLOCKED_MISSING_OBJECT_INPUT: {path}")
        specs.append(spec)
    if not specs or len({item.object_id for item in specs}) != len(specs):
        raise Stage2Error("STAGE2_BLOCKED_OBJECT_CONFIG: object IDs are empty or duplicated")
    return sorted(specs, key=lambda item: (item.cohort_id, item.object_id))


def _read_table(path: Path, columns: Iterable[str] | None = None) -> pd.DataFrame:
    requested = list(columns) if columns is not None else None
    if path.is_dir() or path.suffix == ".parquet":
        return pd.read_parquet(path, columns=requested)
    return pd.read_csv(path, sep="\t" if ".tsv" in path.suffixes else ",", usecols=requested)


def _annotation_partition_path(spec: ObjectSpec) -> Path:
    direct = spec.annotation_path / f"cohort_id={spec.cohort_id}.parquet"
    if not direct.exists():
        matches = [
            path
            for path in spec.annotation_path.glob("cohort_id=*.parquet")
            if path.stem.removeprefix("cohort_id=").casefold()
            == spec.cohort_id.casefold()
        ]
        if len(matches) != 1:
            reason = "AMBIGUOUS" if matches else "MISSING"
            raise Stage2Error(
                f"STAGE2_BLOCKED_{reason}_ANNOTATION_PARTITION: {spec.cohort_id}"
            )
        direct = matches[0]
    return direct


def _annotation_partition(spec: ObjectSpec) -> pd.DataFrame:
    direct = _annotation_partition_path(spec)
    columns = [
        "source_object_id",
        "source_cell_id",
        "analysis_unit_key",
        "study_subject_key",
        "patient_key",
        "patient_id",
        "normalized_timepoint",
        "tissue_context",
        "assignment_status",
        "harmonized_coarse_label",
        "harmonized_mid_label",
        "harmonized_fine_label",
        "original_annotation",
        "mapping_confidence",
        "low_quality_flag",
        "ambient_rna_risk_flag",
    ]
    frame = pd.read_parquet(direct, columns=columns)
    if "source_object_id" in frame:
        frame = frame.loc[frame.source_object_id.astype(str).eq(spec.object_id)].copy()
    if frame.empty:
        raise Stage2Error(f"STAGE2_BLOCKED_EMPTY_ANNOTATION_PARTITION: {spec.object_id}")
    return frame


def _identity_table(spec: ObjectSpec, obs: pd.DataFrame) -> tuple[pd.DataFrame, dict[str, Any]]:
    sidecar = _read_table(spec.identity_sidecar)
    if "source_object_id" in sidecar:
        sidecar = sidecar.loc[sidecar.source_object_id.astype(str).eq(spec.object_id)].copy()
    annotation = _annotation_partition(spec)
    for name, frame in (("sidecar", sidecar), ("annotation", annotation)):
        if "source_cell_id" not in frame or frame.source_cell_id.astype(str).duplicated().any():
            raise Stage2Error(f"STAGE2_BLOCKED_DUPLICATE_CELL_ID: {spec.object_id}:{name}")

    if spec.h5ad_cell_id_column == "__row_index__":
        source_cells = pd.Series([f"row::{i}" for i in range(len(obs))])
    elif spec.h5ad_cell_id_column:
        if spec.h5ad_cell_id_column not in obs:
            raise Stage2Error(f"STAGE2_BLOCKED_MISSING_H5AD_CELL_ID: {spec.object_id}")
        source_cells = obs[spec.h5ad_cell_id_column].astype(str).reset_index(drop=True)
    else:
        source_cells = pd.Series(obs.index.astype(str))

    table = pd.DataFrame({"_row": np.arange(len(obs)), "source_cell_id": source_cells.to_numpy()})
    side_keep = [column for column in SAFE_IDENTITY_COLUMNS if column in sidecar]
    table = table.merge(sidecar[side_keep], on="source_cell_id", how="left", validate="one_to_one")
    annotation_keep = [
        "source_cell_id",
        "analysis_unit_key",
        "study_subject_key",
        "patient_key",
        "patient_id",
        "normalized_timepoint",
        "tissue_context",
        "assignment_status",
        "harmonized_coarse_label",
        "harmonized_mid_label",
        "harmonized_fine_label",
        "original_annotation",
        "mapping_confidence",
        "low_quality_flag",
        "ambient_rna_risk_flag",
    ]
    annotation = annotation[annotation_keep].rename(
        columns={column: f"ann__{column}" for column in annotation_keep if column != "source_cell_id"}
    )
    table = table.merge(annotation, on="source_cell_id", how="left", validate="one_to_one")

    ann_resolved = table.get("ann__assignment_status", pd.Series("", index=table.index)).fillna("").astype(str).str.lower().str.startswith("resolved")
    for column in (
        "analysis_unit_key",
        "study_subject_key",
        "patient_key",
        "patient_id",
        "normalized_timepoint",
        "tissue_context",
        "assignment_status",
    ):
        incoming = f"ann__{column}"
        if incoming in table:
            if column not in table:
                table[column] = ""
            table.loc[ann_resolved, column] = table.loc[ann_resolved, incoming]

    table["patient_key"] = table.get("patient_key", pd.Series("", index=table.index)).fillna("").astype(str)
    subject = table.get("study_subject_key", pd.Series("", index=table.index)).fillna("").astype(str)
    patient_id = table.get("patient_id", pd.Series("", index=table.index)).fillna("").astype(str)
    table["patient_key"] = table.patient_key.where(table.patient_key.str.strip().ne(""), subject)
    table["patient_key"] = table.patient_key.where(table.patient_key.str.strip().ne(""), patient_id)
    table["timepoint"] = table.get("normalized_timepoint", pd.Series("unknown", index=table.index)).fillna("unknown").astype(str)
    table["tissue_context"] = table.get("tissue_context", pd.Series("unknown", index=table.index)).fillna("unknown").astype(str)
    table["coarse"] = table.get("ann__harmonized_coarse_label", pd.Series("Unknown", index=table.index)).fillna("Unknown").astype(str)
    table["mid"] = table.get("ann__harmonized_mid_label", pd.Series("Unknown", index=table.index)).fillna("Unknown").astype(str)
    table["fine"] = table.get("ann__harmonized_fine_label", pd.Series("Unknown", index=table.index)).fillna("Unknown").astype(str)
    table["original_label"] = table.get("ann__original_annotation", pd.Series("", index=table.index)).fillna("").astype(str)
    table["mapping_confidence"] = table.get("ann__mapping_confidence", pd.Series("unmapped", index=table.index)).fillna("unmapped").astype(str)
    table["low_quality_flag"] = pd.to_numeric(
        table.get("ann__low_quality_flag", pd.Series(0, index=table.index)), errors="coerce"
    ).fillna(0).astype(int)
    table["ambient_rna_risk_flag"] = pd.to_numeric(
        table.get("ann__ambient_rna_risk_flag", pd.Series(0, index=table.index)), errors="coerce"
    ).fillna(0).astype(int)
    status = table.get("assignment_status", pd.Series("", index=table.index)).fillna("").astype(str).str.lower()
    table["_resolved"] = status.str.startswith(("resolved", "eligible"))
    audit = {
        "total_cells": int(len(table)),
        "resolved_cells": int(table._resolved.sum()),
        "unresolved_cells": int((~table._resolved).sum()),
        "missing_analysis_unit": int(table.analysis_unit_key.fillna("").astype(str).str.strip().eq("").sum()),
        "missing_patient": int(table.patient_key.str.strip().eq("").sum()),
    }
    return table, audit


def _raw_gene_ids(spec: ObjectSpec, adata: Any) -> pd.Index:
    if spec.h5ad_gene_column:
        if spec.h5ad_gene_column not in adata.var:
            raise Stage2Error(f"STAGE2_BLOCKED_MISSING_GENE_COLUMN: {spec.object_id}")
        genes = pd.Index(adata.var[spec.h5ad_gene_column].astype(str))
    else:
        genes = pd.Index(adata.var_names.astype(str))
    if len(genes) != adata.n_vars or (genes.str.strip() == "").any():
        raise Stage2Error(f"STAGE2_BLOCKED_INVALID_GENE_NAMESPACE: {spec.object_id}")
    return genes


def _atomic_sparse_write(matrix: sparse.spmatrix, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, name = tempfile.mkstemp(prefix=f".{path.stem}.", suffix=".npz", dir=path.parent)
    os.close(fd)
    temporary = Path(name)
    try:
        sparse.save_npz(temporary, matrix.tocsr(), compressed=True)
        temporary.replace(path)
    finally:
        temporary.unlink(missing_ok=True)


def _aggregate_level(
    source: Any,
    identity: pd.DataFrame,
    level: str,
    genes: pd.Index,
    chunk_size: int,
) -> tuple[sparse.csr_matrix, pd.DataFrame, dict[str, Any]]:
    eligible = (
        identity._resolved
        & identity.analysis_unit_key.notna()
        & identity.patient_key.astype(str).str.strip().ne("")
        & identity.analysis_unit_key.astype(str).str.strip().ne("")
        & identity[level].astype(str).str.strip().ne("")
    )
    selected = identity.loc[eligible].copy()
    if selected.empty:
        raise Stage2Error(f"STAGE2_BLOCKED_NO_ELIGIBLE_CELLS: level={level}")
    key_columns = ["analysis_unit_key", "patient_key", "timepoint", "tissue_context", level]
    keys = selected[key_columns].astype(str).drop_duplicates().sort_values(key_columns).reset_index(drop=True)
    key_to_code = {tuple(row): index for index, row in keys.iterrows()}
    selected["_code"] = np.fromiter(
        (key_to_code[tuple(row)] for _, row in selected[key_columns].astype(str).iterrows()),
        dtype=np.int64,
        count=len(selected),
    )
    total = sparse.csr_matrix((len(keys), len(genes)), dtype=np.float64)
    source_total = 0.0
    for start in range(0, len(selected), chunk_size):
        chunk = selected.iloc[start : start + chunk_size]
        matrix = source[chunk._row.to_numpy()]
        matrix = matrix.tocsr() if sparse.issparse(matrix) else sparse.csr_matrix(matrix)
        matrix = matrix.astype(np.float64, copy=False)
        if matrix.data.size:
            if np.any(matrix.data < 0) or not np.all(np.isfinite(matrix.data)):
                raise Stage2Error("STAGE2_BLOCKED_INVALID_COUNTS: negative or non-finite value")
            if not np.allclose(matrix.data, np.rint(matrix.data), rtol=0, atol=1e-6):
                raise Stage2Error("STAGE2_BLOCKED_LAYER_SEMANTICS: declared counts are non-integer")
        source_total += float(matrix.sum())
        selector = sparse.csr_matrix(
            (np.ones(len(chunk)), (chunk._code.to_numpy(), np.arange(len(chunk)))),
            shape=(len(keys), len(chunk)),
        )
        total += selector @ matrix

    units = keys.rename(columns={level: "cell_state"})
    units["cell_state_level"] = level
    units["n_cells_used"] = np.bincount(selected._code.to_numpy(), minlength=len(keys))
    units["library_size"] = np.asarray(total.sum(axis=1)).ravel()
    units["n_genes_detected"] = np.asarray((total > 0).sum(axis=1)).ravel()
    units["low_quality_or_ambient_state"] = units.cell_state.eq("Low_quality_or_ambient")
    pseudobulk_total = float(total.sum())
    audit = {
        "cell_state_level": level,
        "n_eligible_cells": int(len(selected)),
        "n_expression_units": int(len(units)),
        "eligible_cell_total": source_total,
        "pseudobulk_total": pseudobulk_total,
        "counts_conserved": bool(np.isclose(source_total, pseudobulk_total, rtol=0, atol=0)),
    }
    if not audit["counts_conserved"]:
        raise Stage2Error("STAGE2_BLOCKED_COUNTS_CONSERVATION")
    return total, units, audit


def aggregate_object(
    spec: ObjectSpec,
    output_root: Path,
    chunk_size: int = 50_000,
    study_family: str | None = None,
    input_fingerprint: dict[str, Any] | None = None,
) -> tuple[list[AggregatedShard], list[dict[str, Any]]]:
    import anndata as ad

    adata = ad.read_h5ad(spec.h5ad, backed="r")
    try:
        if spec.count_layer != "X" and spec.count_layer not in adata.layers:
            raise Stage2Error(f"STAGE2_BLOCKED_MISSING_COUNT_LAYER: {spec.object_id}")
        source = adata.X if spec.count_layer == "X" else adata.layers[spec.count_layer]
        identity, identity_audit = _identity_table(spec, adata.obs)
        genes = _raw_gene_ids(spec, adata)
        shards: list[AggregatedShard] = []
        audits: list[dict[str, Any]] = []
        safe_id = spec.object_id.replace("::", "__").replace("/", "_")
        for level in ("coarse", "mid"):
            matrix, units, audit = _aggregate_level(source, identity, level, genes, chunk_size)
            units.insert(0, "expression_unit_id", [
                f"{spec.cohort_id}::{spec.object_id}::{row.analysis_unit_key}::{level}::{row.cell_state}::count_pseudobulk"
                for row in units.itertuples()
            ])
            units.insert(1, "cohort_id", spec.cohort_id)
            units.insert(2, "study_family", study_family or spec.cohort_id)
            units.insert(3, "object_id", spec.object_id)
            units["sample_key"] = units.analysis_unit_key
            units["modality"] = "scRNA"
            units["platform"] = "scRNA"
            units["expression_layer"] = "counts"
            base = output_root / "cache" / "pseudobulk" / f"{safe_id}__{level}"
            matrix_path = base.with_suffix(".counts.npz")
            units_path = base.with_suffix(".units.parquet")
            genes_path = base.with_suffix(".genes.json")
            audit_path = base.with_suffix(".audit.json")
            _atomic_sparse_write(matrix, matrix_path)
            units_path.parent.mkdir(parents=True, exist_ok=True)
            units.to_parquet(units_path, index=False)
            write_json({"genes": genes.astype(str).tolist()}, genes_path)
            record = {
                **identity_audit,
                **audit,
                "cohort_id": spec.cohort_id,
                "study_family": study_family or spec.cohort_id,
                "object_id": spec.object_id,
                "h5ad_path": str(spec.h5ad),
                "h5ad_bytes": spec.h5ad.stat().st_size,
                "count_layer": spec.count_layer,
                "n_raw_features": len(genes),
                "input_fingerprint": input_fingerprint
                or shard_input_fingerprint(spec, study_family=study_family),
                "matrix_path": str(matrix_path),
                "units_path": str(units_path),
                "genes_path": str(genes_path),
                "status": "PASS",
            }
            write_json(record, audit_path)
            audits.append(record)
            shards.append(AggregatedShard(spec.cohort_id, spec.object_id, level, matrix_path, units_path, genes_path, audit_path))
        return shards, audits
    finally:
        if getattr(adata, "file", None) is not None:
            adata.file.close()


def load_shards(output_root: Path) -> list[AggregatedShard]:
    shards: list[AggregatedShard] = []
    for audit_path in sorted((output_root / "cache" / "pseudobulk").glob("*.audit.json")):
        record = json.loads(audit_path.read_text(encoding="utf-8"))
        shard = AggregatedShard(
            cohort_id=record["cohort_id"],
            object_id=record["object_id"],
            level=record["cell_state_level"],
            matrix_path=Path(record["matrix_path"]),
            units_path=Path(record["units_path"]),
            genes_path=Path(record["genes_path"]),
            audit_path=audit_path,
        )
        if not all(path.exists() for path in (shard.matrix_path, shard.units_path, shard.genes_path)):
            raise Stage2Error(f"STAGE2_BLOCKED_INCOMPLETE_SHARD: {audit_path}")
        shards.append(shard)
    return shards


def _path_state(path: Path) -> object:
    if path.is_file():
        stat = path.stat()
        return {
            "kind": "file",
            "size": stat.st_size,
            "mtime_ns": stat.st_mtime_ns,
        }
    if not path.is_dir():
        raise Stage2Error(f"STAGE2_BLOCKED_FINGERPRINT_PATH: {path}")
    return {
        "kind": "directory",
        "files": [
            {
                "path": str(item.relative_to(path)),
                "size": item.stat().st_size,
                "mtime_ns": item.stat().st_mtime_ns,
            }
            for item in sorted(candidate for candidate in path.rglob("*") if candidate.is_file())
        ],
    }


def _path_sha256(path: Path) -> str:
    if path.is_file():
        return sha256(path)
    if not path.is_dir():
        raise Stage2Error(f"STAGE2_BLOCKED_FINGERPRINT_PATH: {path}")
    records = [
        (str(item.relative_to(path)), sha256(item))
        for item in sorted(candidate for candidate in path.rglob("*") if candidate.is_file())
    ]
    return hashlib.sha256(
        json.dumps(records, separators=(",", ":"), ensure_ascii=False).encode("utf-8")
    ).hexdigest()


def shard_input_fingerprint(
    spec: ObjectSpec,
    *,
    study_family: str | None = None,
    semantic_context: dict[str, Any] | None = None,
    previous: dict[str, Any] | None = None,
) -> dict[str, Any]:
    """Hash every source and semantic choice that can change a cached shard."""

    annotation = _annotation_partition_path(spec)
    previous = previous or {}

    def hashed_path(path: Path, key: str) -> tuple[str, object, str]:
        state = _path_state(path)
        previous_state = previous.get(f"{key}_state")
        previous_hash = str(previous.get(f"{key}_sha256", ""))
        same_state = previous_state is not None and previous_state == state
        same_path = str(previous.get(key, "")) == str(path)
        if (
            previous_hash
            and previous_hash != "partitioned_dataset"
            and same_path
            and (same_state or previous_state is None)
        ):
            provenance = "reused_previous_hash_stat_verified" if same_state else "reused_legacy_hash_manifest"
            return previous_hash, state, provenance
        return _path_sha256(path), state, "computed_current"

    h5ad_hash, h5ad_state, h5ad_provenance = hashed_path(spec.h5ad, "h5ad")
    identity_hash, identity_state, identity_provenance = hashed_path(
        spec.identity_sidecar, "identity_sidecar"
    )
    annotation_hash, annotation_state, annotation_provenance = hashed_path(
        annotation, "annotation_partition"
    )
    payload = {
        "cohort_id": spec.cohort_id,
        "object_id": spec.object_id,
        "h5ad": str(spec.h5ad),
        "h5ad_sha256": h5ad_hash,
        "h5ad_state": h5ad_state,
        "identity_sidecar": str(spec.identity_sidecar),
        "identity_sidecar_sha256": identity_hash,
        "identity_sidecar_state": identity_state,
        "annotation_partition": str(annotation),
        "annotation_partition_sha256": annotation_hash,
        "annotation_partition_state": annotation_state,
        "count_layer": spec.count_layer,
        "h5ad_cell_id_column": spec.h5ad_cell_id_column or "",
        "h5ad_gene_column": spec.h5ad_gene_column or "",
        "study_family": study_family or spec.cohort_id,
        "semantic_context": semantic_context or {},
        "hash_provenance": {
            "h5ad": h5ad_provenance,
            "identity_sidecar": identity_provenance,
            "annotation_partition": annotation_provenance,
        },
    }
    fingerprint_payload = {
        key: value for key, value in payload.items() if key != "hash_provenance"
    }
    payload["fingerprint_sha256"] = hashlib.sha256(
        json.dumps(
            fingerprint_payload,
            sort_keys=True,
            separators=(",", ":"),
            ensure_ascii=False,
        ).encode("utf-8")
    ).hexdigest()
    return payload
