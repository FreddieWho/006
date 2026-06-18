#!/usr/bin/env python3
"""Freeze context_repair_v1 handoffs and Phase7/8 configuration files."""
from __future__ import annotations

import argparse
import hashlib
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd
import yaml


ROOT = Path(__file__).resolve().parents[3]
ORIGINAL_P6 = ROOT / "results/v6_2/phase6_module_algorithm_benchmark_strengthened"
ORIGINAL_P7_CFG = ROOT / "scripts/v6_2/phase7_v6_2_1_config.yaml"
ORIGINAL_P8_CFG = ROOT / "scripts/v6_2/phase8_v6_2_1_config.yaml"


def now_iso() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat()


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(1024 * 1024):
            h.update(chunk)
    return h.hexdigest()


def one(root: Path, names: tuple[str, ...]) -> Path:
    for name in names:
        path = root / name
        if path.exists():
            return path.resolve()
    found = []
    for name in names:
        found.extend(root.glob(f"**/{Path(name).name}"))
    if len(found) == 1:
        return found[0].resolve()
    raise FileNotFoundError(f"Expected exactly one of {names} under {root}; found={found}")


def metadata_path(root: Path, stem: str) -> Path:
    candidates = [
        root / f"{stem}.effective_v1.csv",
        root / f"{stem}.frozen_v1.csv",
    ]
    for path in candidates:
        if path.exists():
            return path.resolve()
    found = list(root.glob(f"{stem}*effective*v1*.csv")) + list(root.glob(f"{stem}*frozen*v1*.csv"))
    if found:
        return found[0].resolve()
    raise FileNotFoundError(f"Missing v1 metadata table: {stem}")


def build(args: argparse.Namespace) -> dict:
    score = one(args.phase6_root, (
        "consensus/module_score_matrix.frozen_v1.parquet",
        "frozen_module_score_matrix.parquet",
        "frozen_module_score_matrix.csv",
    ))
    if score.suffix == ".csv":
        converted = args.phase6_root / "consensus/module_score_matrix.frozen_v1.parquet"
        converted.parent.mkdir(parents=True, exist_ok=True)
        pd.read_csv(score).to_parquet(converted, index=False)
        score = converted.resolve()
    registry = one(args.phase6_root, ("foundation/unified_expression_unit_registry.csv", "unified_expression_unit_registry.csv", "expression_unit_registry.csv"))
    logcpm = one(args.phase6_root, (
        "foundation/logcpm_cellstate_pseudobulk.v1.parquet",
        "logcpm_cellstate_pseudobulk.parquet",
        "repaired_cell_state_pseudobulk.logcpm.parquet",
    ))
    membership = ORIGINAL_P6 / "consensus/module_membership.frozen_v1.csv"
    dictionary = ORIGINAL_P6 / "consensus/module_dictionary.frozen_v1.csv"
    support = ORIGINAL_P6 / "consensus/module_method_support_map.frozen_v1.csv"
    original_hash = sha256(membership)
    projection_manifest = one(args.phase6_root, ("projection_manifest.yaml", "handoff/projection_manifest.yaml"))
    projection = yaml.safe_load(projection_manifest.read_text()) or {}
    if projection.get("verdict") != "PASS":
        raise RuntimeError(f"Projection manifest is not PASS: {projection.get('verdict')}")
    recorded_hash = projection.get("membership_sha256") or projection.get("membership_sha256_after")
    if recorded_hash and recorded_hash != original_hash:
        raise RuntimeError("Frozen membership hash differs from projection manifest")
    handoff_dir = args.phase6_root / "handoff"
    handoff_dir.mkdir(parents=True, exist_ok=True)
    handoff_path = handoff_dir / "phase6_to_phase7_module_entry_manifest.yaml"
    handoff = {
        "phase": "phase6_frozen_measurement_context_repair_v1",
        "mode": "frozen_projection_only",
        "verdict": "GO_TO_PHASE7_WITH_FROZEN_MODULE_ENTRY",
        "primary_backbone_method": "lda_topic",
        "primary_backbone_run": "lda_topic__rank8__batch__inverse_k__tw0.01",
        "primary_selection_status": "membership_unchanged_remeasurement_only",
        "consensus_strategy": "frozen_primary_backbone_no_retraining",
        "default_module_score_matrix": str(score),
        "default_module_membership": str(membership.resolve()),
        "module_dictionary": str(dictionary.resolve()),
        "module_method_support_map": str(support.resolve()),
        "phase7_allowed_to_read_only_frozen_module_entry": True,
        "membership_sha256": original_hash,
        "module_retraining_performed": False,
    }
    handoff_path.write_text(yaml.safe_dump(handoff, sort_keys=False, allow_unicode=True))

    cohort = metadata_path(args.metadata_root, "cohort_registry")
    sample = args.phase4b_root / "response_environment/sample_metadata_analysis_unit_v1.csv"
    patient = args.phase4b_root / "response_environment/patient_metadata_analysis_unit_v1.csv"
    if not sample.exists() or not patient.exists():
        raise FileNotFoundError("Phase4B analysis-unit sample/patient metadata must be built before downstream configs")
    roles = metadata_path(args.metadata_root, "dataset_role_and_feature_eligibility")
    response = metadata_path(args.metadata_root, "response_label_environment")
    p7 = yaml.safe_load(ORIGINAL_P7_CFG.read_text())
    p7["output_dir"] = str(args.phase7_root.relative_to(ROOT))
    p7["inputs"].update({
        "phase6_handoff": str(handoff_path.relative_to(ROOT)),
        "phase6_unit_registry": str(registry.relative_to(ROOT)),
        "phase6_logcpm": str(logcpm.relative_to(ROOT)),
        "sample_metadata": str(sample.relative_to(ROOT)),
        "patient_metadata": str(patient.relative_to(ROOT)),
        "cohort_registry": str(cohort.relative_to(ROOT)),
        "sample_fractions": str((args.phase4b_root / "fractions/cell_state_fraction_matrix.sample_level.csv").relative_to(ROOT)),
        "patient_fractions": str((args.phase4b_root / "fractions/cell_state_fraction_matrix.patient_timepoint_level.csv").relative_to(ROOT)),
        "qc_covariates": str((args.phase4b_root / "qc_covariates/sample_patient_timepoint_qc_covariates.csv").relative_to(ROOT)),
        "cell_state_reliability": str((ROOT / "results/v6_2/phase4a_cell_state_harmonization/qc/cell_state_reliability_scores.csv").relative_to(ROOT)),
    })
    args.config_root.mkdir(parents=True, exist_ok=True)
    p7_path = args.config_root / "phase7_context_repair_v1.yaml"
    p7_path.write_text(yaml.safe_dump(p7, sort_keys=False, allow_unicode=True))

    p8 = yaml.safe_load(ORIGINAL_P8_CFG.read_text())
    p8["output_dir"] = str(args.phase8_root.relative_to(ROOT))
    p8["inputs"].update({
        "phase7_handoff": str((args.phase7_root / "handoff/phase7b_to_phase8_handoff.yaml").relative_to(ROOT)),
        "barrier_matrix": str((args.phase7_root / "handoff/eligible_barrier_score_matrix_v0.parquet").relative_to(ROOT)),
        "barrier_table": str((args.phase7_root / "identifiability/eligible_barrier_set_v0.csv").relative_to(ROOT)),
        "coarse_matrix": str((args.phase7_root / "measurement/module_activity_by_patient_timepoint_residualized.coarse_sensitivity.parquet").relative_to(ROOT)),
        "phase6_membership": str(membership.relative_to(ROOT)),
        "sample_metadata": str(sample.relative_to(ROOT)),
        "patient_metadata": str(patient.relative_to(ROOT)),
        "cohort_registry": str(cohort.relative_to(ROOT)),
        "dataset_roles": str(roles.relative_to(ROOT)),
        "response_environment": str(response.relative_to(ROOT)),
        "response_binding": str((args.phase4b_root / "response_environment/response_environment_binding_table.csv").relative_to(ROOT)),
        "sample_metadata": str(sample.relative_to(ROOT)),
        "patient_metadata": str(patient.relative_to(ROOT)),
    })
    p8_path = args.config_root / "phase8_context_repair_v1.yaml"
    p8_path.write_text(yaml.safe_dump(p8, sort_keys=False, allow_unicode=True))
    manifest = {
        "release": "context_repair_v1", "created_at": now_iso(),
        "module_membership_sha256": original_hash, "module_retraining_performed": False,
        "phase6_handoff": str(handoff_path), "phase7_config": str(p7_path), "phase8_config": str(p8_path),
        "phase7_output": str(args.phase7_root), "phase8_output": str(args.phase8_root),
    }
    (args.config_root / "downstream_config_manifest.yaml").write_text(
        yaml.safe_dump(manifest, sort_keys=False, allow_unicode=True))
    return manifest


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--phase6-root", type=Path, default=ROOT / "results/v6_2/phase6_measurement_foundation_repair_v1")
    parser.add_argument("--phase4b-root", type=Path, default=ROOT / "results/v6_2/phase4b_immune_state_feature_construction_repair_v1")
    parser.add_argument("--metadata-root", type=Path, default=ROOT / "results/v6_2/data_interface_repair_v1/metadata")
    parser.add_argument("--phase7-root", type=Path, default=ROOT / "results/v6_2/phase7_module_measurement_and_barrier_identifiability_repair_v1")
    parser.add_argument("--phase8-root", type=Path, default=ROOT / "results/v6_2/phase8_anchor_repair_and_statistical_strengthening_repair_v1")
    parser.add_argument("--config-root", type=Path, default=ROOT / "results/v6_2/data_interface_repair_v1/configs")
    args = parser.parse_args()
    for name in ("phase6_root", "phase4b_root", "metadata_root", "phase7_root", "phase8_root", "config_root"):
        setattr(args, name, getattr(args, name).resolve())
    print(yaml.safe_dump(build(args), sort_keys=False, allow_unicode=True))


if __name__ == "__main__":
    main()
