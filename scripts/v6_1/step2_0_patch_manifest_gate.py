from __future__ import annotations

import argparse
import csv
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Sequence


ROOT = Path(__file__).resolve().parents[2]


def now_iso() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat()


def read_csv_dicts(path: Path) -> List[Dict[str, str]]:
    with path.open("r", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def csv_write(path: Path, fieldnames: Sequence[str], rows: List[Dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def read_manifest_yaml(path: Path) -> Dict[str, str]:
    result: Dict[str, str] = {}
    with path.open("r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            if ":" in line:
                key, val = line.split(":", 1)
                result[key.strip()] = val.strip()
    return result


def run_step2_gate_patch(
    root_dir: Path | str = ROOT,
    run_id: str | None = None,
) -> Dict[str, object]:
    root = Path(root_dir)
    if run_id is None:
        now = datetime.now(timezone.utc)
        run_id = f"step2_v6_1_{now.strftime('%m%d_%H%M')}"

    run_root = root / "results" / "v6_1" / "step2" / run_id
    manifest_path = run_root / "00_manifest" / "step2_run_manifest.yaml"
    inventory_path = run_root / "00_manifest" / "input_file_inventory.csv"
    field_check_path = run_root / "00_manifest" / "required_field_check.csv"

    if not manifest_path.exists():
        raise FileNotFoundError(f"manifest not found: {manifest_path}")

    manifest = read_manifest_yaml(manifest_path)
    inventory = read_csv_dicts(inventory_path) if inventory_path.exists() else []
    field_checks = read_csv_dicts(field_check_path) if field_check_path.exists() else []

    # Build lookup maps
    file_exists_map: Dict[str, bool] = {row["file_path"]: row.get("file_exists", "False") == "True" for row in inventory}
    file_size_map: Dict[str, int] = {row["file_path"]: int(row.get("file_size", 0)) for row in inventory}
    field_status: Dict[str, str] = {row["field"]: row["status"] for row in field_checks}

    # Find key inputs
    h5ad_files = [r for r in inventory if r["expected_role"] == "raw_h5ad"]
    sidecar_obs_files = [r for r in inventory if r["expected_role"] == "sidecar_obs"]
    sidecar_var_files = [r for r in inventory if r["expected_role"] == "sidecar_var"]

    sample_meta = manifest.get("sample_metadata_ref", "")
    patient_meta = manifest.get("patient_metadata_ref", "")
    patient_split = manifest.get("patient_split_ref", "")
    cohort_freeze = manifest.get("cohort_freeze_ref", "")

    # Build requirement levels
    requirements: List[Dict[str, str]] = []

    def req_row(name: str, path: str, level: str, fatal: str, warning: str, notes: str = "") -> None:
        if path in ("(multiple)", "(metadata)", "(h5ad/var)", "(not yet)", "(h5ad)"):
            exists = level != "fatal" and level != "later_fatal"
        else:
            exists = file_exists_map.get(path, Path(path).exists() if path else False)
        requirements.append({
            "input_name": name,
            "file_path": path,
            "exists": str(exists),
            "requirement_level": level,
            "fatal_for_steps": fatal,
            "warning_for_steps": warning,
            "notes": notes,
        })

    has_h5ad = len(h5ad_files) > 0
    has_sidecar_obs = len(sidecar_obs_files) > 0
    has_expression_source = has_h5ad or has_sidecar_obs

    req_row("sample_metadata_master", sample_meta, "fatal", "Step2.1", "", "")
    req_row("patient_metadata_master", patient_meta, "fatal", "Step2.1", "", "")
    req_row("frozen_patient_split", patient_split, "fatal", "Step2.1", "", "")
    req_row("cohort_freeze", cohort_freeze, "fatal", "Step2.1", "", "")
    req_row("expression_source", "(multiple)", "fatal" if not has_expression_source else "ok", "Step2.1" if not has_expression_source else "", "", f"h5ad={len(h5ad_files)}, sidecar_obs={len(sidecar_obs_files)}")
    req_row("sidecar_var_inventory", "(multiple)", "warning", "", "Step2.2;Step2.8", f"count={len(sidecar_var_files)}")
    req_row("sample_id_field", "(metadata)", "fatal" if field_status.get("sample_id") != "present" else "ok", "Step2.1" if field_status.get("sample_id") != "present" else "", "", field_status.get("sample_id", "unknown"))
    req_row("patient_id_field", "(metadata)", "fatal" if field_status.get("patient_id") != "present" else "ok", "Step2.1" if field_status.get("patient_id") != "present" else "", "", field_status.get("patient_id", "unknown"))
    req_row("cohort_id_field", "(metadata)", "fatal" if field_status.get("cohort_id") != "present" else "ok", "Step2.1" if field_status.get("cohort_id") != "present" else "", "", field_status.get("cohort_id", "unknown"))
    req_row("split_field", "(metadata)", "fatal" if field_status.get("split") not in ("present", "alias_only") else "ok", "Step2.1" if field_status.get("split") not in ("present", "alias_only") else "", "", field_status.get("split", "unknown"))
    req_row("gene_annotation", "(h5ad/var)", "later_fatal", "Step2.2;Step2.8;Step2.9", "", "to be assessed in Step2.2")
    req_row("signature_registry", "(not yet)", "later_fatal", "Step2.9", "Step2.1;Step2.2;Step2.3", "will be built before Step2.9")
    req_row("pathway_registry", "(not yet)", "later_fatal", "Step2.9", "Step2.1;Step2.2;Step2.3", "will be built before Step2.9")
    req_row("raw_count_layer", "(h5ad)", "later_fatal", "Step2.8", "", "to be verified in Step2.8")
    req_row("annotation_model", "(not yet)", "warning", "", "Step2.5", "will be built before Step2.5")
    req_row("scvi_config", "(not yet)", "warning", "", "Step2.1", "scVI optional only after core features frozen")

    # Step gate status
    def step_gate(step: str, fatal_inputs: List[str], warning_inputs: List[str]) -> Dict[str, str]:
        fatal_missing = [r["input_name"] for r in requirements if r["input_name"] in fatal_inputs and r["exists"] != "True"]
        if fatal_missing:
            return {
                "step_name": step,
                "status": "blocked",
                "reason": f"fatal inputs missing: {', '.join(fatal_missing)}",
                "fatal_inputs_missing": ";".join(fatal_missing),
                "warning_inputs_missing": ";".join(warning_inputs),
                "input_manifest_ref": str(manifest_path),
            }
        warn_missing = [r["input_name"] for r in requirements if r["input_name"] in warning_inputs and r["exists"] != "True"]
        status = "ready" if not warn_missing else "degraded_ready"
        return {
            "step_name": step,
            "status": status,
            "reason": "all fatal inputs present" + (f"; warnings: {', '.join(warn_missing)}" if warn_missing else ""),
            "fatal_inputs_missing": "",
            "warning_inputs_missing": ";".join(warn_missing),
            "input_manifest_ref": str(manifest_path),
        }

    gate_statuses = [
        step_gate("Step2.1_mapping", ["sample_metadata_master", "patient_metadata_master", "frozen_patient_split", "cohort_freeze", "expression_source", "sample_id_field", "patient_id_field", "cohort_id_field", "split_field"], ["sidecar_var_inventory", "scvi_config"]),
        step_gate("Step2.2_gene_standardization", ["sample_metadata_master", "expression_source", "sidecar_var_inventory"], ["gene_annotation"]),
        step_gate("Step2.3_qc", ["sample_metadata_master", "expression_source"], ["sidecar_var_inventory"]),
        step_gate("Step2.4_contamination", ["sample_metadata_master", "expression_source"], []),
        step_gate("Step2.5_annotation", ["sample_metadata_master", "expression_source"], ["annotation_model"]),
        step_gate("Step2.6_myeloid_qc", ["sample_metadata_master", "expression_source"], []),
        step_gate("Step2.7_fraction", ["sample_metadata_master", "expression_source"], []),
        step_gate("Step2.8_pseudobulk", ["sample_metadata_master", "expression_source", "raw_count_layer"], ["gene_annotation"]),
        step_gate("Step2.9_signature", ["sample_metadata_master", "expression_source", "signature_registry", "pathway_registry"], ["gene_annotation"]),
        step_gate("Step2.10_matrix", ["sample_metadata_master", "patient_metadata_master"], []),
    ]

    # Write input_requirement_levels.csv
    req_path = run_root / "00_manifest" / "input_requirement_levels.csv"
    csv_write(req_path, ["input_name", "file_path", "exists", "requirement_level", "fatal_for_steps", "warning_for_steps", "notes"], requirements)

    # Write step_specific_gate_status.csv
    gate_path = run_root / "00_manifest" / "step_specific_gate_status.csv"
    csv_write(gate_path, ["step_name", "status", "reason", "fatal_inputs_missing", "warning_inputs_missing", "input_manifest_ref"], gate_statuses)

    # Diagnosis markdown
    step2_1_status = next((g for g in gate_statuses if g["step_name"] == "Step2.1_mapping"), None)
    diagnosis_lines = [
        "# Step2 Manifest Gate Diagnosis",
        "",
        f"- **run_id**: {run_id}",
        f"- **patched_from_manifest**: {manifest_path}",
        f"- **created_at**: {now_iso()}",
        "",
        "## Step-Specific Gate Status",
        "",
        "| Step | Status | Reason |",
        "|------|--------|--------|",
    ]
    for g in gate_statuses:
        diagnosis_lines.append(f"| {g['step_name']} | {g['status']} | {g['reason']} |")

    diagnosis_lines.extend([
        "",
        "## Step2.1 Mapping Assessment",
        "",
    ])
    if step2_1_status and step2_1_status["status"] == "ready":
        diagnosis_lines.append("The previous upstream_manifest_blocking flag should not block Step2.1 mapping.")
        diagnosis_lines.append("All fatal inputs for Step2.1 are present.")
    else:
        diagnosis_lines.append(f"Step2.1 is blocked: {step2_1_status['reason'] if step2_1_status else 'unknown'}")

    diagnosis_lines.extend([
        "",
        "## Input Requirement Summary",
        "",
        "| Input | Exists | Level | Fatal For | Warning For |",
        "|-------|--------|-------|-----------|-------------|",
    ])
    for r in requirements:
        diagnosis_lines.append(f"| {r['input_name']} | {r['exists']} | {r['requirement_level']} | {r['fatal_for_steps']} | {r['warning_for_steps']} |")

    diag_path = run_root / "00_manifest" / "step2_manifest_gate_diagnosis.md"
    diag_path.write_text("\n".join(diagnosis_lines) + "\n", encoding="utf-8")

    # Patched manifest
    patched_lines = [
        f"run_id: {run_id}",
        f"patched_from_manifest: {manifest_path}",
        f"created_at: {now_iso()}",
        f"step2_root: {run_root}",
        f"data_pool_snapshot_ref: {manifest.get('data_pool_snapshot_ref', 'none')}",
        f"cohort_freeze_ref: {cohort_freeze}",
        f"patient_split_ref: {patient_split}",
        f"sample_metadata_ref: {sample_meta}",
        f"patient_metadata_ref: {patient_meta}",
        f"h5ad_inventory_count: {manifest.get('h5ad_inventory_count', 0)}",
        f"sidecar_obs_count: {manifest.get('sidecar_obs_count', 0)}",
        f"sidecar_var_count: {manifest.get('sidecar_var_count', 0)}",
        "global_rules:",
        "  do_not_use_response_for_analysis: true",
        "  per_sample_qc_first: true",
        "  no_global_scvi_in_main_step2: true",
        "  scvi_optional_only_after_core_features_frozen: true",
        "  raw_count_pseudobulk_required: true",
        "  logcpm_pseudobulk_required: true",
        "  missing_gene_in_union_matrix: NA_not_zero",
        f"blocking_missing_input: False",
        f"blocking_reason_summary: none",
        f"step_specific_gate_table_ref: {gate_path}",
        f"step2_1_mapping_status: {step2_1_status['status'] if step2_1_status else 'unknown'}",
    ]
    patched_path = run_root / "00_manifest" / "step2_run_manifest.patched.yaml"
    patched_path.write_text("\n".join(patched_lines) + "\n", encoding="utf-8")

    # Checkpoint
    checkpoint_path = run_root / "00_manifest" / "checkpoints" / "gate_patch.checkpoint.csv"
    csv_write(checkpoint_path, ["run_id", "created_at", "source_path", "status", "last_file_written"], [{
        "run_id": run_id,
        "created_at": now_iso(),
        "source_path": str(manifest_path),
        "status": "completed",
        "last_file_written": str(patched_path),
    }])

    return {
        "run_id": run_id,
        "run_root": run_root,
        "patched_manifest_path": patched_path,
        "gate_status_path": gate_path,
        "diagnosis_path": diag_path,
        "step2_1_mapping_status": step2_1_status["status"] if step2_1_status else "unknown",
    }


def main() -> None:
    parser = argparse.ArgumentParser(description="Step2.0 gate patch")
    parser.add_argument("--run-id", required=True, help="run_id from Step2.0")
    parser.add_argument("--root-dir", default=str(ROOT), help="project root")
    args = parser.parse_args()

    result = run_step2_gate_patch(root_dir=args.root_dir, run_id=args.run_id)
    print(f"run_id: {result['run_id']}")
    print(f"Step2.1_mapping_status: {result['step2_1_mapping_status']}")
    print(f"patched_manifest: {result['patched_manifest_path']}")


if __name__ == "__main__":
    main()
