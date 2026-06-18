#!/usr/bin/env python3
"""Merge bounded cohort-level identity repairs into the canonical sidecar root."""

from __future__ import annotations

import argparse
import json
import shutil
from pathlib import Path

import pandas as pd
import yaml


ROOT = Path(__file__).resolve().parents[3]
DEFAULT_TARGET = ROOT / "results/v6_2/data_interface_repair_v1"
TABLES = (
    "cohort_cardinality_audit.csv",
    "object_metadata_audit.csv",
    "sample_identity_sidecar.csv",
    "analysis_context_sidecar.csv",
)


def parse_source(value: str) -> tuple[Path, set[str]]:
    root, separator, cohorts = value.partition(":")
    if not separator or not cohorts:
        raise argparse.ArgumentTypeError("source must be ROOT:COHORT[,COHORT...]")
    source = Path(root).resolve()
    if cohorts.strip() == "*":
        audit = source / "cohort_cardinality_audit.csv"
        if not audit.exists():
            raise argparse.ArgumentTypeError(f"cannot expand * without {audit}")
        selected = set(pd.read_csv(audit)["cohort_id"].astype(str))
    else:
        selected = {item.strip() for item in cohorts.split(",") if item.strip()}
    return source, selected


def merge_table(
    target: Path,
    sources: list[tuple[Path, set[str]]],
    allow_empty_replacement: bool = False,
) -> None:
    current = pd.read_csv(target)
    replacements = []
    replaced: set[str] = set()
    for root, cohorts in sources:
        source = pd.read_csv(root / target.name)
        requested = {item.casefold(): item for item in cohorts}
        selected = source[source["cohort_id"].astype(str).str.casefold().isin(requested)].copy()
        missing = set(requested) - set(selected["cohort_id"].astype(str).str.casefold())
        if missing and not allow_empty_replacement:
            raise RuntimeError(f"{root / target.name}: missing cohorts {sorted(missing)}")
        replacements.append(selected)
        replaced.update(requested)
    merged = pd.concat(
        [current[~current["cohort_id"].astype(str).str.casefold().isin(replaced)], *replacements],
        ignore_index=True,
        sort=False,
    )
    merged = merged.sort_values(["cohort_id"] + (["source_object_id"] if "source_object_id" in merged else []))
    staged = target.with_suffix(target.suffix + ".new")
    merged.to_csv(staged, index=False)
    staged.replace(target)


def replace_partition(target_root: Path, source_root: Path, cohort: str) -> None:
    source = source_root / "cell_identity.parquet" / f"cohort_id={cohort}"
    target = target_root / "cell_identity.parquet" / f"cohort_id={cohort}"
    if not source.is_dir() or not any(source.glob("*.parquet")):
        raise RuntimeError(f"missing repaired cell partition: {source}")
    staged = target.with_name(target.name + ".new")
    if staged.exists():
        shutil.rmtree(staged)
    shutil.copytree(source, staged)
    if target.exists():
        shutil.rmtree(target)
    staged.replace(target)


def validate_and_manifest(target: Path) -> dict:
    objects = pd.read_csv(target / "object_metadata_audit.csv")
    cohorts = pd.read_csv(target / "cohort_cardinality_audit.csv")
    if objects["source_object_id"].astype(str).duplicated().any():
        raise RuntimeError("duplicate source_object_id after merge")
    if cohorts["cohort_id"].astype(str).duplicated().any():
        raise RuntimeError("duplicate cohort_id in cardinality audit after merge")
    if not (cohorts["n_cells"].astype(int) == cohorts["resolved_cells"].astype(int) + cohorts["quarantined_cells"].astype(int)).all():
        raise RuntimeError("cell conservation failed in cohort cardinality audit")
    manifest = {
        "identity_version": "context_repair_v1",
        "metadata_only": True,
        "response_fields_used": False,
        "raw_h5ad_mutated": False,
        "n_objects": int(len(objects)),
        "n_cells": int(cohorts["n_cells"].sum()),
        "resolved_cells": int(cohorts["resolved_cells"].sum()),
        "quarantined_cells": int(cohorts["quarantined_cells"].sum()),
        "n_cohorts": int(cohorts["cohort_id"].nunique()),
    }
    (target / "manifest.yaml").write_text(yaml.safe_dump(manifest, sort_keys=False))
    (target / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    return manifest


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--target", type=Path, default=DEFAULT_TARGET)
    parser.add_argument("--source", action="append", type=parse_source, required=True)
    args = parser.parse_args()
    target = args.target.resolve()
    seen: set[str] = set()
    for source, cohorts in args.source:
        overlap = seen & cohorts
        if overlap:
            raise SystemExit(f"cohort appears in multiple sources: {sorted(overlap)}")
        seen.update(cohorts)
        for cohort in sorted(cohorts):
            replace_partition(target, source, cohort)
    for name in TABLES:
        # A cohort can have no unique source-sample row while still having
        # resolved patient-timepoint analysis units in the context sidecar.
        merge_table(
            target / name,
            args.source,
            allow_empty_replacement=name == "sample_identity_sidecar.csv",
        )
    print(json.dumps(validate_and_manifest(target), indent=2))


if __name__ == "__main__":
    main()
