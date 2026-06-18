#!/usr/bin/env python3
"""Build obs-first context_repair_v1 sidecars without mutating source h5ad."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import yaml

from identity import IDENTITY_VERSION, apply_conflict_quarantine, build_cell_identity, cardinality_audit, cohort_from_obs


ROOT = Path(__file__).resolve().parents[3]


def _decode(values: np.ndarray) -> list[str]:
    return [value.decode() if isinstance(value, bytes) else str(value) for value in values]


def read_obs(path: Path) -> tuple[str, pd.DataFrame]:
    """Read only h5ad /obs, including AnnData categorical encodings."""
    with h5py.File(path, "r") as handle:
        obs_group = handle["obs"]
        index_key = obs_group.attrs.get("_index", "_index")
        index_key = index_key.decode() if isinstance(index_key, bytes) else str(index_key)
        index = _decode(obs_group[index_key][()])
        columns: dict[str, list[object]] = {}
        for name, node in obs_group.items():
            if name == index_key:
                continue
            if isinstance(node, h5py.Group) and {"categories", "codes"}.issubset(node.keys()):
                categories = _decode(node["categories"][()])
                codes = node["codes"][()]
                columns[name] = [categories[code] if code >= 0 else None for code in codes]
            elif isinstance(node, h5py.Dataset) and len(node.shape) == 1:
                columns[name] = _decode(node[()])
    obs = pd.DataFrame(columns, index=pd.Index(index, name="obs_name"))
    return cohort_from_obs(path.stem, obs), obs


def write_cell_partition(frame: pd.DataFrame, path: Path, cohort_id: str, source_object_id: str) -> None:
    """Write one object at a time; never concatenate all cell sidecars in RAM."""
    safe_object = source_object_id.replace("/", "_").replace(":", "_")
    partition = path / f"cohort_id={cohort_id}"
    partition.mkdir(parents=True, exist_ok=True)
    # pyarrow materializes string columns during conversion; chunking keeps a
    # large object's metadata-only sidecar below the process memory ceiling.
    for start in range(0, len(frame), 100_000):
        frame.iloc[start:start + 100_000].to_parquet(partition / f"{safe_object}.part-{start // 100_000:04d}.parquet", index=False)


def build(config: dict, output_dir: Path, cohorts: set[str] | None = None) -> dict:
    raw_dir = ROOT / config["inputs"]["raw_h5ad_dir"]
    registry_path = ROOT / config["inputs"].get(
        "object_registry", "results/v6_2/phase3_5_single_cell_processing_qc_gate/analysis_object_registry.frozen_v0.csv")
    registry = pd.read_csv(registry_path) if registry_path.exists() else pd.DataFrame()
    object_by_path = {
        str(Path(row.object_path).resolve()): str(row.object_id)
        for row in registry.itertuples() if getattr(row, "object_path", "")
    }
    objects, audits, samples, contexts = [], [], [], []
    n_cells = resolved_cells = quarantined_cells = 0
    cell_root = output_dir / "cell_identity.parquet"
    for path in sorted(raw_dir.glob("*.h5ad")):
        if cohorts is not None and path.stem.lower() not in {name.lower() for name in cohorts}:
            continue
        try:
            cohort_id, obs = read_obs(path)
            if obs.index.duplicated().any():
                obs["_source_barcode"] = obs.index.astype(str)
                obs.index = pd.Index([f"row::{i}" for i in range(len(obs))], name="source_row_id")
        except Exception as exc:
            objects.append({"cohort_id": path.stem, "source_object_id": f"unreadable::{path.stem}", "h5ad_path": str(path.relative_to(ROOT)), "n_obs": 0, "metadata_only": True, "assignment_status": "quarantined", "quarantine_reason": f"obs_read_failed:{type(exc).__name__}"})
            continue
            continue
        source_object_id = object_by_path.get(str(path.resolve()), f"{cohort_id}::object_1")
        cells = apply_conflict_quarantine(build_cell_identity(obs, cohort_id, source_object_id))
        write_cell_partition(cells, cell_root, cohort_id, source_object_id)
        audits.append(cardinality_audit(cells))
        sample_rows = cells.loc[:, ["source_sample_id", "biological_sample_key", "cohort_id", "study_subject_key", "library_key", "demux_id", "normalized_timepoint", "tissue_context", "lesion_context", "treatment_arm", "analysis_unit_key", "assignment_status", "assignment_provenance", "quarantine_reason"]].copy()
        sample_rows = sample_rows.rename(columns={"source_sample_id": "sample_id"})
        cardinality = sample_rows.groupby(["cohort_id", "sample_id"], dropna=False).study_subject_key.transform("nunique")
        samples.append(sample_rows[cardinality.eq(1)].drop_duplicates(["cohort_id", "sample_id"]))
        contexts.append(cells.drop_duplicates("analysis_unit_key").loc[:, ["analysis_unit_key", "cohort_id", "study_subject_key", "normalized_timepoint", "tissue_context", "lesion_context", "treatment_arm", "assignment_status", "assignment_provenance", "quarantine_reason"]])
        n_cells += len(cells)
        resolved_cells += int(cells.assignment_status.eq("resolved").sum())
        quarantined_cells += int(cells.assignment_status.eq("quarantined").sum())
        objects.append({"cohort_id": cohort_id, "source_object_id": source_object_id, "h5ad_path": str(path.relative_to(ROOT)), "n_obs": len(obs), "metadata_only": True, "assignment_status": "resolved", "quarantine_reason": ""})
    if not audits:
        raise RuntimeError(f"no h5ad files in {raw_dir}")
    output_dir.mkdir(parents=True, exist_ok=True)
    pd.concat(samples, ignore_index=True).drop_duplicates(["cohort_id", "sample_id"]).to_csv(output_dir / "sample_identity_sidecar.csv", index=False)
    pd.concat(contexts, ignore_index=True).drop_duplicates("analysis_unit_key").to_csv(output_dir / "analysis_context_sidecar.csv", index=False)
    audit = pd.concat(audits, ignore_index=True).sort_values("cohort_id")
    audit.to_csv(output_dir / "cohort_cardinality_audit.csv", index=False)
    pd.DataFrame(objects).sort_values("cohort_id").to_csv(output_dir / "object_metadata_audit.csv", index=False)
    manifest = {"identity_version": IDENTITY_VERSION, "metadata_only": True, "response_fields_used": False, "raw_h5ad_mutated": False, "n_objects": len(objects), "n_cells": n_cells, "resolved_cells": resolved_cells, "quarantined_cells": quarantined_cells}
    (output_dir / "manifest.yaml").write_text(yaml.safe_dump(manifest, sort_keys=False))
    (output_dir / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    return manifest


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", type=Path, default=ROOT / "scripts/v6_2/context_repair_v1.yaml")
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--cohorts", help="optional comma-separated smoke subset; omitted audits every raw h5ad")
    args = parser.parse_args()
    config = yaml.safe_load(args.config.read_text())
    try:
        cohorts = set(args.cohorts.split(",")) if args.cohorts else None
        print(json.dumps(build(config, args.output_dir, cohorts), sort_keys=True))
    except Exception as exc:
        print(f"context_repair_v1 failed: {exc}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
