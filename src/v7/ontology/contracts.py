"""Shared contracts and deterministic I/O for v7 Stage 2."""

from __future__ import annotations

import hashlib
import json
import os
import re
import tempfile
from pathlib import Path
from typing import Any, Iterable

import pandas as pd
import yaml


class Stage2Error(RuntimeError):
    """Fail-closed Stage 2 contract violation."""


def resolve_path(value: str | Path, base: Path) -> Path:
    path = Path(value)
    return path if path.is_absolute() else (base / path).resolve()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def load_config(path: Path) -> dict[str, Any]:
    config_path = path.resolve()
    raw = yaml.safe_load(config_path.read_text(encoding="utf-8")) or {}
    if raw.get("stage") != "v7_stage2_ontology_and_measurement":
        raise Stage2Error("STAGE2_BLOCKED_CONFIG: unexpected stage identifier")
    base = config_path.parent
    inputs = raw.get("inputs") or {}
    if not inputs:
        raise Stage2Error("STAGE2_BLOCKED_CONFIG: inputs are missing")
    raw["inputs"] = {key: resolve_path(value, base) for key, value in inputs.items()}
    raw["output_root"] = resolve_path(raw.get("output_root", "../../results/v7/ontology"), base)
    return raw


def prohibited_columns(columns: Iterable[str], tokens: Iterable[str]) -> list[str]:
    hits: list[str] = []
    for column in columns:
        text = str(column).lower()
        for token in tokens:
            pattern = rf"(^|[^a-z0-9]){re.escape(str(token).lower())}([^a-z0-9]|$)"
            if re.search(pattern, text):
                hits.append(str(column))
                break
    return sorted(set(hits))


def require_inputs(config: dict[str, Any]) -> pd.DataFrame:
    records: list[dict[str, Any]] = []
    missing: list[str] = []
    for role, path in sorted(config["inputs"].items()):
        exists = path.exists()
        if not exists:
            missing.append(role)
        records.append(
            {
                "input_role": role,
                "path": str(path),
                "exists": exists,
                "bytes": path.stat().st_size if exists and path.is_file() else "",
                "sha256": sha256(path) if exists and path.is_file() else "",
                "status": "PASS" if exists else "MISSING",
            }
        )
    if missing:
        raise Stage2Error(f"STAGE2_BLOCKED_MISSING_INPUT: {','.join(missing)}")
    return pd.DataFrame(records)


def _atomic_replace(path: Path, writer: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, name = tempfile.mkstemp(prefix=f".{path.name}.", suffix=".tmp", dir=path.parent)
    os.close(fd)
    temporary = Path(name)
    try:
        writer(temporary)
        temporary.replace(path)
    finally:
        temporary.unlink(missing_ok=True)


def write_tsv(frame: pd.DataFrame, path: Path, sort_by: list[str] | None = None) -> None:
    material = frame.copy()
    order = sort_by or list(material.columns)
    order = [column for column in order if column in material]
    if order and not material.empty:
        material = material.sort_values(order, kind="mergesort", na_position="last").reset_index(drop=True)
    _atomic_replace(path, lambda target: material.to_csv(target, sep="\t", index=False, lineterminator="\n"))


def write_json(payload: Any, path: Path) -> None:
    text = json.dumps(payload, indent=2, ensure_ascii=False, sort_keys=True) + "\n"
    _atomic_replace(path, lambda target: target.write_text(text, encoding="utf-8"))


def write_yaml(payload: Any, path: Path) -> None:
    text = yaml.safe_dump(payload, sort_keys=False, allow_unicode=True)
    _atomic_replace(path, lambda target: target.write_text(text, encoding="utf-8"))


def content_manifest(paths: Iterable[Path]) -> list[dict[str, Any]]:
    records = []
    for path in sorted(paths, key=lambda item: str(item)):
        records.append(
            {
                "path": str(path),
                "bytes": path.stat().st_size if path.exists() and path.is_file() else "",
                "sha256": sha256(path) if path.exists() and path.is_file() else "",
                "exists": path.exists(),
            }
        )
    return records


def architecture_surrogate_contract(version: str) -> dict[str, Any]:
    return {
        "contract": "v7_architecture_surrogate",
        "version": version,
        "stage2_status": "SCHEMA_ONLY_NOT_MODEL_SELECTION",
        "target": {
            "status": "UNSELECTED_UNTIL_STAGE4",
            "required_semantics": "response_blind_patient_level_topology",
        },
        "allowed_predictors": [
            "D2_measurable_or_uncertain_module_scores",
            "canonical_cell_state_fractions",
            "state_specific_module_scores",
            "pretreatment_covariates_available_before_outcome",
        ],
        "forbidden_predictors": [
            "response",
            "RECIST",
            "survival",
            "post_outcome_variable",
            "dataset_id_as_predictor",
            "response_proxy",
            "direct_topology_derived_predictor",
            "ground_truth_or_GT_derived_feature",
        ],
        "evaluation": {
            "outer_unit": "patient",
            "same_patient_or_leakage_group_across_folds": False,
            "dataset_holdout_required": True,
            "cross_fitting_required": True,
            "uncertainty_or_abstention_required": True,
            "abundance_only_baseline_required": True,
        },
    }
