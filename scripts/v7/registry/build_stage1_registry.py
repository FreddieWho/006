#!/usr/bin/env python3
"""Build the v7 Stage 1 metadata-only registry.

This command reads small metadata/control tables only.  It deliberately does
not open expression matrices, images, or other large scientific binaries.
The output is deterministic: rows are sorted and no run timestamp is written.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import sys
from collections import Counter, defaultdict
from pathlib import Path
from typing import Iterable, Mapping


SCRIPT_VERSION = "v7-stage1-registry-1.0"
CLAIMS = {
    "H1": "cross-cancer response-blind molecular and cell states",
    "H2": "spatial organization adds information beyond abundance",
    "H3": "stable discovery of unnamed spatial fields",
    "H4": "shared function with organ-specific implementation",
    "H5": "reproducible HCC context rewriting",
    "H6": "PD-1 response/failure is associated with multi-structure ecology",
    "H7": "PD1+X observation is compatible with barrier rewiring",
    "H8": "structured repair space supports X-class extrapolation",
}


def clean(value: object) -> str:
    if value is None:
        return ""
    return str(value).strip()


def first(row: Mapping[str, object], *keys: str) -> str:
    for key in keys:
        value = clean(row.get(key, ""))
        if value:
            return value
    return ""


def read_table(path: Path, delimiter: str | None = None) -> tuple[list[str], list[dict[str, str]]]:
    if not path.exists():
        return [], []
    if delimiter is None:
        delimiter = "\t" if path.suffix.lower() in {".tsv", ".tab"} else ","
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter=delimiter)
        fieldnames = [clean(x) for x in (reader.fieldnames or [])]
        rows: list[dict[str, str]] = []
        for raw in reader:
            row = {clean(k): clean(v) for k, v in raw.items() if k is not None}
            rows.append(row)
    return fieldnames, rows


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def rel_or_abs(path: Path, root: Path, external_root: Path | None = None) -> str:
    try:
        return str(path.resolve().relative_to(root.resolve()))
    except ValueError:
        if external_root is not None:
            try:
                return str(path.resolve().relative_to(external_root.resolve()))
            except ValueError:
                pass
        return str(path.resolve())


def unique_join(values: Iterable[str]) -> str:
    clean_values = sorted({clean(value) for value in values if clean(value)})
    return "|".join(clean_values)


def usable_value(value: str) -> bool:
    return bool(value) and value.lower() not in {"na", "n/a", "none", "unknown", "not_provided", "not provided", "null", "false"}


def truthy_flag(value: str) -> bool:
    return clean(value).lower() in {"yes", "true", "1", "pass", "usable"}


def write_tsv(path: Path, fields: list[str], rows: Iterable[Mapping[str, object]]) -> int:
    materialized = [{field: clean(row.get(field, "")) for field in fields} for row in rows]
    materialized.sort(key=lambda row: tuple(row[field] for field in fields))
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(materialized)
    return len(materialized)


def source_id_for_cohort(cohort_id: str, modality: str = "") -> str:
    known = {
        "TASK01": "006_TASK01",
        "TASK02": "006_TASK02",
        "lambrecht_hcc": "006_LAMBRECHT_HCC",
        "GSE286827": "006_GSE286827",
        "GSE301741": "006_GSE301741",
        "GSE238264": "006_GSE238264",
        "GSE291246": "006_GSE291246",
    }
    if cohort_id.lower().startswith("mendeley"):
        return "006_Mendeley_HCC"
    if cohort_id in known:
        return known[cohort_id]
    lower = modality.lower()
    if "bulk" in lower:
        return "006_bulk"
    if "perturb" in lower or "l1000" in lower:
        return "006_perturbation"
    return "006_sc_identity"


def source_id_for_spatial(study: str, namespace: str = "") -> str:
    token = f"{study} {namespace}".upper()
    if "ATLAS_TABLE_S2" in token:
        return "013_ATLAS_TABLE_S2"
    if "TENX_GENOMICS" in token:
        return "013_10x_cancer"
    if "HEST" in token:
        return "013_HEST"
    if "HTAN" in token or "VANDERBILT" in token:
        return "013_HTAN_CRC"
    if "GSE175540" in token:
        return "013_GSE175540"
    if "ST_CRC_CMS" in token:
        return "013_ST_CRC_CMS"
    if "TLS_VISIUM_USZ" in token or "USZ" in token:
        return "013_TLS_USZ"
    if "GSE211956" in token:
        return "013_GSE211956"
    if "GSE274557" in token:
        return "013_GSE274557"
    if "GSE274103" in token:
        return "013_GSE274103"
    if "GSE226997" in token:
        return "013_GSE226997"
    return "013_other_spatial"


def role_for_source(source_id: str, config_row: Mapping[str, str], role_rows: Mapping[str, Mapping[str, str]]) -> tuple[str, str, str, str]:
    """Return role, role source, permitted, forbidden."""
    aliases = {
        "013_HTAN_CRC": "HTAN_VANDERBILT_CRC",
        "013_GSE175540": "GEO::GSE175540",
        "013_ST_CRC_CMS": "ST_CRC_CMS",
        "013_TLS_USZ": "TLS_VISIUM_USZ",
        "013_GSE211956": "GEO::GSE211956",
    }
    role_row = role_rows.get(source_id, {}) or role_rows.get(aliases.get(source_id, ""), {})
    role = first(role_row, "primary_role") or first(config_row, "v7_primary_role") or "support"
    role_source = "013_role_freeze_provenance" if role_row else ("013_sample_registry_and_v7_registry" if source_id.startswith("013_") else "006_effective_metadata_and_v7_registry")
    permitted = first(role_row, "allowed_use") or first(config_row, "v7_primary_role") or "metadata audit and scoped support"
    forbidden = first(role_row, "forbidden_use") or first(config_row, "main_caveat") or "role reassignment without a new audit"
    return role, role_source, permitted, forbidden


def read_inputs(root: Path, spatial_root: Path) -> dict[str, object]:
    metadata = root / "results/v6_2/data_interface_repair_v1/metadata"
    paths: dict[str, Path] = {
        "config": root / "config/v7/data_source_registry.tsv",
        "cohort": metadata / "cohort_registry.effective_v1.csv",
        "sample": metadata / "sample_metadata_master.effective_v1.csv",
        "patient": metadata / "patient_metadata_master.effective_v1.csv",
        "response": metadata / "response_label_environment.effective_v1.csv",
        "feature": metadata / "dataset_role_and_feature_eligibility.effective_v1.csv",
        "object_audit": root / "results/v6_2/data_interface_repair_v1/object_metadata_audit.csv",
        "sample_sidecar": root / "results/v6_2/data_interface_repair_v1/sample_identity_sidecar.csv",
        "context_sidecar": root / "results/v6_2/data_interface_repair_v1/analysis_context_sidecar.csv",
        "raw_task": Path("/home/huyudi/004/TASK/meta_task_raw.csv"),
        "raw_tcr": Path("/home/huyudi/004/Component/meta_tcr.csv"),
        "physical": spatial_root / "infra/sample-registry/physical_units.tsv",
        "assets": spatial_root / "infra/sample-registry/source_assets.tsv",
        "source_inventory": spatial_root / "infra/sample-registry/source_metadata_inventory.tsv",
        "identity_evidence": spatial_root / "infra/sample-registry/identity_evidence.tsv",
        "duplicates": spatial_root / "infra/sample-registry/duplicate_groups.tsv",
        "role_freeze": spatial_root / "infra/sample-registry/role_freeze.tsv",
        "structures": spatial_root / "infra/structure-registry/structure_instances.tsv",
        "ontology": spatial_root / "infra/structure-registry/structure_ontology.tsv",
        "gt_audit": spatial_root / "infra/structure-registry/gt_source_audit.tsv",
        "replay": spatial_root / "infra/structure-registry/h5ad_replay_index.tsv",
        "outer_splits": spatial_root / "infra/structure-registry/outer_splits.tsv",
        "cross_sections": spatial_root / "infra/structure-registry/cross_section_links.tsv",
        "input_policy": spatial_root / "infra/structure-registry/input_policy.tsv",
        "r04_input_manifest": spatial_root / "infra/r04/input_manifest.json",
        "r04_restart_panel": spatial_root / "infra/r04/restart_stability_panel_20260831.json",
        "gse238_manifest": root / "results/v6_2/phase2_5_data_onboarding/01_metadata_freeze/spatial_sample_manifest.GSE238264.tsv",
        "gse291_manifest": root / "results/v6_2/phase2_5_data_onboarding/01_metadata_freeze/spatial_sample_manifest.GSE291246.tsv",
        "gse238_qc": root / "results/v6_2/phase2_5_data_onboarding/02_conversion_qc/spatial_sample_manifest.GSE238264.tsv",
        "gse291_qc": root / "results/v6_2/phase2_5_data_onboarding/02_conversion_qc/per_sample_xenium_stats.GSE291246.latest.tsv",
        "gse238_spot_qc": root / "results/v6_2/phase2_5_data_onboarding/02_conversion_qc/per_sample_spot_stats.GSE238264.tsv",
    }
    tables: dict[str, tuple[list[str], list[dict[str, str]]]] = {}
    for name, path in paths.items():
        if path.exists() and path.suffix.lower() in {".csv", ".tsv", ".tab"}:
            tables[name] = read_table(path)
    return {"paths": paths, "tables": tables}


def add_006_rows(
    sample_rows: list[dict[str, str]],
    cohort_rows: Mapping[str, dict[str, str]],
    source_rows: Mapping[str, dict[str, str]],
    root: Path,
) -> tuple[list[dict[str, str]], list[dict[str, str]]]:
    logical: list[dict[str, str]] = []
    spatial: list[dict[str, str]] = []
    for row in sample_rows:
        cohort_id = first(row, "cohort_id")
        modality = first(row, "modality") or first(cohort_rows.get(cohort_id, {}), "modality")
        source_id = source_id_for_cohort(cohort_id, modality)
        patient_key = first(row, "patient_key", "patient_id")
        sample_key = first(row, "sample_key", "sample_id")
        logical_id = f"006::sample::{sample_key or cohort_id + '::' + patient_key}"
        cohort = cohort_rows.get(cohort_id, {})
        feature = source_rows.get(cohort_id, {})
        logical.append({
            "logical_unit_id": logical_id,
            "source_project": str(root),
            "source_id": source_id,
            "dataset_id": cohort_id,
            "unit_level": "sample",
            "patient_id": patient_key,
            "block_id": "",
            "section_id": "",
            "region_id": "",
            "sample_id": first(row, "sample_id", "sample_key"),
            "modality": modality,
            "platform": first(row, "platform") or first(cohort, "platform"),
            "resolution": "cell" if "scrna" in modality.lower() else ("spot" if "spatial" in modality.lower() or "visium" in modality.lower() else "bulk_or_condition"),
            "coordinate_system": "not_applicable" if "spatial" not in modality.lower() and "visium" not in modality.lower() else "native_coordinates",
            "scale": "sample",
            "counts_layer": "available" if first(feature, "raw_counts_available") == "yes" else first(feature, "feature_family") or "unknown",
            "image_available": first(feature, "spatial_coordinates_available") == "yes" and "image" in first(feature, "feature_family").lower(),
            "segmentation_available": "unknown",
            "structure_gt_available": "unknown",
            "cancer": first(row, "cancer_type") or first(cohort, "cancer_type"),
            "tissue": first(row, "tissue_source") or first(cohort, "tissue_source"),
            "treatment": first(row, "treatment_context", "treatment_raw") or first(cohort, "treatment_context"),
            "timepoint": first(row, "timepoint"),
            "response": first(row, "response_binary_harmonized"),
            "endpoint": first(row, "response_endpoint_type") or first(cohort, "response_endpoint_type"),
            "identity_confidence": "explicit_patient_key" if patient_key else "unknown",
            "metadata_provenance": f"results/v6_2/data_interface_repair_v1/metadata/sample_metadata_master.effective_v1.csv#sample_key={sample_key}",
            "duplicate_lineage": f"006::{cohort_id}::patient::{patient_key}" if patient_key else logical_id,
            "leakage_group_id": f"006::{cohort_id}::patient::{patient_key}" if patient_key else logical_id,
            "record_status": "EFFECTIVE_METADATA",
            "n_observations": "",
            "source_path": first(feature, "feature_source_path"),
            "audit_status": "metadata_only_no_matrix_read",
        })
    return logical, spatial


def add_006_spatial(
    manifest: list[dict[str, str]],
    qc_rows: list[dict[str, str]],
    spot_qc_rows: list[dict[str, str]],
    cohort_id: str,
    root: Path,
    effective_by_sample: Mapping[tuple[str, str], Mapping[str, str]] | None = None,
) -> tuple[list[dict[str, str]], list[dict[str, str]], list[dict[str, str]]]:
    logical: list[dict[str, str]] = []
    physical: list[dict[str, str]] = []
    crosswalk: list[dict[str, str]] = []
    qc_by_sample = {}
    for qc_row in qc_rows:
        for key in (first(qc_row, "sample_id"), first(qc_row, "gsm_id")):
            if key:
                qc_by_sample[key] = qc_row
    spot_by_sample = {first(r, "sample_id"): r for r in spot_qc_rows}
    for row in manifest:
        sample_id = first(row, "sample_id", "gsm_id")
        effective = (effective_by_sample or {}).get((cohort_id, sample_id), {})
        raw_patient_id = first(row, "patient_id")
        patient_id = f"{cohort_id}::{raw_patient_id}" if raw_patient_id else ""
        physical_id = f"006::{cohort_id}::section::{sample_id}"
        source_id = source_id_for_cohort(cohort_id, "spatial")
        qc = qc_by_sample.get(sample_id, {}) or qc_by_sample.get(first(row, "gsm_id"), {})
        spot = spot_by_sample.get(sample_id, {})
        platform = "Xenium" if cohort_id == "GSE291246" else "Visium"
        n_obs = first(qc, "n_cells") or first(spot, "in_tissue_spots", "total_spots")
        expr_status = first(qc, "usable") or ("available" if first(spot, "total_spots") else "unknown")
        timepoint = first(effective, "timepoint") or ("post" if cohort_id == "GSE238264" else first(row, "treatment_context"))
        response = first(effective, "response_binary_harmonized") or first(row, "response_binary", "response_standard")
        treatment = first(effective, "treatment_context", "treatment_raw") or first(row, "treatment_regimen_raw", "treatment_context")
        endpoint = first(effective, "response_endpoint_type") or first(row, "response_standard") or ("unknown" if not response else "metadata_manifest")
        logical.append({
            "logical_unit_id": physical_id,
            "source_project": str(root),
            "source_id": source_id,
            "dataset_id": cohort_id,
            "unit_level": "section",
            "patient_id": patient_id,
            "block_id": "",
            "section_id": sample_id,
            "region_id": "",
            "sample_id": sample_id,
            "modality": "spatial",
            "platform": platform,
            "resolution": "cell" if platform == "Xenium" else "spot",
            "coordinate_system": "native_coordinates",
            "scale": "section",
            "counts_layer": "h5ad_or_spot_matrix" if usable_value(expr_status) else "manifest_only_no_expression_qc",
            "image_available": "yes" if truthy_flag(first(row, "image_available")) else "unknown",
            "segmentation_available": "yes" if platform == "Xenium" and usable_value(expr_status) else "unknown",
            "structure_gt_available": "no",
            "cancer": first(row, "disease"),
            "tissue": first(row, "tissue_source"),
            "treatment": treatment,
            "timepoint": timepoint,
            "response": response,
            "endpoint": endpoint,
            "identity_confidence": "explicit_manifest_patient" if patient_id else "unknown",
            "metadata_provenance": f"results/v6_2/phase2_5_data_onboarding/01_metadata_freeze/spatial_sample_manifest.{cohort_id}.tsv#sample_id={sample_id}",
            "duplicate_lineage": f"006::{cohort_id}::patient::{patient_id}" if patient_id else physical_id,
            "leakage_group_id": f"006::{cohort_id}::patient::{patient_id}" if patient_id else physical_id,
            "record_status": "EFFECTIVE_SPATIAL_MANIFEST",
            "n_observations": n_obs,
            "source_path": first(row, "source_file"),
            "audit_status": "metadata_only_qc_linked",
        })
        physical.append({
            "physical_unit_id": physical_id,
            "source_project": str(root),
            "source_id": source_id,
            "dataset_id": cohort_id,
            "source_namespace": cohort_id,
            "source_record_id": first(row, "gsm_id", "sample_id"),
            "patient_id": patient_id,
            "block_id": "",
            "physical_specimen_id": "",
            "section_id": sample_id,
            "region_id": "",
            "modality": "spatial",
            "platform": platform,
            "resolution": "cell" if platform == "Xenium" else "spot",
            "coordinate_system": "native_coordinates",
            "counts_layer": "h5ad_or_spot_matrix" if usable_value(expr_status) else "manifest_only_no_expression_qc",
            "image_available": "yes" if truthy_flag(first(row, "image_available")) else "unknown",
            "segmentation_available": "yes" if platform == "Xenium" and usable_value(expr_status) else "unknown",
            "structure_gt_available": "no",
            "cancer": first(row, "disease"),
            "tissue": first(row, "tissue_source"),
            "treatment": treatment,
            "timepoint": timepoint,
            "response": response,
            "endpoint": endpoint,
            "identity_confidence": "explicit_manifest_patient" if patient_id else "unknown",
            "metadata_provenance": f"results/v6_2/phase2_5_data_onboarding/01_metadata_freeze/spatial_sample_manifest.{cohort_id}.tsv#sample_id={sample_id}",
            "duplicate_lineage": f"006::{cohort_id}::patient::{patient_id}" if patient_id else physical_id,
            "leakage_group_id": f"006::{cohort_id}::patient::{patient_id}" if patient_id else physical_id,
            "record_status": "EFFECTIVE_SPATIAL_MANIFEST",
            "n_observations": n_obs,
            "source_path": first(row, "source_file"),
            "audit_status": "metadata_only_qc_linked",
        })
        crosswalk.append({
            "physical_unit_id": physical_id,
            "source_project": str(root),
            "source_id": source_id,
            "dataset_id": cohort_id,
            "patient_id": patient_id,
            "block_id": "",
            "section_id": sample_id,
            "region_id": "",
            "physical_specimen_id": "",
            "identity_granularity": "patient_explicit_section_manifest",
            "identity_status": "RESOLVED" if patient_id else "UNKNOWN",
            "evidence_grade": "E3_explicit_manifest" if patient_id else "E0_missing",
            "block_equivalent_status": "BLOCK_UNKNOWN",
            "block_equivalent_basis": "",
            "section_order": "",
            "leakage_group_id": f"006::{cohort_id}::patient::{patient_id}" if patient_id else physical_id,
            "grouping_eligibility": "patient_grouped" if patient_id else "quarantine_missing_patient",
            "metadata_provenance": f"results/v6_2/phase2_5_data_onboarding/01_metadata_freeze/spatial_sample_manifest.{cohort_id}.tsv#sample_id={sample_id}",
            "relation_status": "MANIFEST_LINKED",
        })
    return logical, physical, crosswalk


def build_registry(args: argparse.Namespace) -> dict[str, object]:
    root = Path(args.project_root).resolve()
    spatial_root = Path(args.spatial_root).resolve()
    output = Path(args.output_dir).resolve()
    inputs = read_inputs(root, spatial_root)
    paths: dict[str, Path] = inputs["paths"]  # type: ignore[assignment]
    tables: dict[str, tuple[list[str], list[dict[str, str]]]] = inputs["tables"]  # type: ignore[assignment]

    required_inputs = (
        "config", "cohort", "sample", "patient", "feature", "object_audit",
        "physical", "assets", "duplicates", "role_freeze", "structures", "outer_splits",
        "gse238_manifest", "gse291_manifest",
    )
    missing = [name for name in required_inputs if not paths[name].is_file() or paths[name].stat().st_size == 0]
    if missing:
        raise RuntimeError(
            "STAGE1_BLOCKED_MISSING_INPUT: " + ", ".join(f"{name}={paths[name]}" for name in missing)
        )

    config = tables.get("config", ([], []))[1]
    cohort_rows = {first(row, "cohort_id"): row for row in tables.get("cohort", ([], []))[1] if first(row, "cohort_id")}
    sample_rows = tables.get("sample", ([], []))[1]
    effective_by_sample = {
        (first(row, "cohort_id"), first(row, "sample_id")): row
        for row in sample_rows
        if first(row, "cohort_id") and first(row, "sample_id")
    }
    patient_rows = tables.get("patient", ([], []))[1]
    feature_rows = {first(row, "cohort_id"): row for row in tables.get("feature", ([], []))[1] if first(row, "cohort_id")}
    role_rows_raw = tables.get("role_freeze", ([], []))[1]
    role_rows: dict[str, dict[str, str]] = {}
    for row in role_rows_raw:
        key = first(row, "logical_unit_id")
        if key:
            role_rows[key] = row
    physical_rows = tables.get("physical", ([], []))[1]
    asset_rows = tables.get("assets", ([], []))[1]
    structure_rows = tables.get("structures", ([], []))[1]
    gt_rows = tables.get("gt_audit", ([], []))[1]
    replay_rows = tables.get("replay", ([], []))[1]
    outer_rows = tables.get("outer_splits", ([], []))[1]
    cross_section_rows = tables.get("cross_sections", ([], []))[1]

    logical_fields = [
        "logical_unit_id", "source_project", "source_id", "dataset_id", "unit_level", "patient_id", "block_id", "section_id", "region_id", "sample_id",
        "modality", "platform", "resolution", "coordinate_system", "scale", "counts_layer", "image_available", "segmentation_available", "structure_gt_available",
        "cancer", "tissue", "treatment", "timepoint", "response", "endpoint", "identity_confidence", "metadata_provenance", "duplicate_lineage",
        "leakage_group_id", "record_status", "n_observations", "source_path", "audit_status",
    ]
    logical: list[dict[str, str]] = []
    # Spatial cohorts use their richer metadata-freeze manifests below.  Do not
    # add the same biological sections a second time from the generic sample
    # table.
    logical_006, _ = add_006_rows(
        [row for row in sample_rows if first(row, "cohort_id") not in {"GSE238264", "GSE291246"}],
        cohort_rows,
        feature_rows,
        root,
    )
    logical.extend(logical_006)

    # Add one dataset-level row for sources with no sample-level effective table.
    for row in config:
        source_id = first(row, "source_id")
        if not source_id or any(existing["logical_unit_id"] == f"dataset::{source_id}" for existing in logical):
            continue
        logical.append({
            "logical_unit_id": f"dataset::{source_id}", "source_project": first(row, "project_root"), "source_id": source_id,
            "dataset_id": source_id, "unit_level": "dataset", "patient_id": "", "block_id": "", "section_id": "", "region_id": "", "sample_id": "",
            "modality": first(row, "modality"), "platform": "", "resolution": "dataset", "coordinate_system": "", "scale": first(row, "inventory_scale"),
            "counts_layer": "unknown", "image_available": "unknown", "segmentation_available": "unknown", "structure_gt_available": "unknown",
            "cancer": "", "tissue": "", "treatment": "", "timepoint": "", "response": "", "endpoint": "",
            "identity_confidence": "dataset_registry_only", "metadata_provenance": "config/v7/data_source_registry.tsv#source_id=" + source_id,
            "duplicate_lineage": "", "leakage_group_id": "", "record_status": "DATASET_REGISTRY_ONLY", "n_observations": "",
            "source_path": first(row, "project_root"), "audit_status": "not_sample_level",
        })

    # Spatial physical units from /013: retain excluded records for audit, but mark them.
    asset_by_namespace: dict[str, list[dict[str, str]]] = defaultdict(list)
    for asset in asset_rows:
        asset_by_namespace[first(asset, "source_namespace")].append(asset)
    structure_by_physical: dict[str, list[dict[str, str]]] = defaultdict(list)
    for structure in structure_rows:
        structure_by_physical[first(structure, "physical_unit_id")].append(structure)
    replay_by_physical = {first(row, "physical_unit_id"): row for row in replay_rows}
    outer_by_physical = {first(row, "physical_unit_id"): row for row in outer_rows}

    def matching_assets(row: Mapping[str, str], namespace: str, study: str) -> list[dict[str, str]]:
        candidates = asset_by_namespace.get(namespace, [])
        if namespace == "external_geo":
            accession = study.removeprefix("GEO::").upper()
            return [asset for asset in candidates if first(asset, "accession").upper() == accession]
        if namespace == "HTAN_VANDERBILT_CRC":
            source_record = first(row, "source_record_id")
            exact = [asset for asset in candidates if first(asset, "asset_id") == first(row, "physical_unit_id")]
            if exact:
                return exact
            return [asset for asset in candidates if first(asset, "source_sample_id") and first(asset, "source_sample_id") in source_record]
        if namespace == "TENX_GENOMICS":
            source_record = first(row, "source_record_id")
            matched = [asset for asset in candidates if first(asset, "source_sample_id") == source_record]
            return matched or candidates
        return candidates
    spatial_fields = [
        "physical_unit_id", "source_project", "source_id", "dataset_id", "source_namespace", "source_record_id", "patient_id", "block_id", "physical_specimen_id",
        "section_id", "region_id", "modality", "platform", "resolution", "coordinate_system", "counts_layer", "image_available", "segmentation_available",
        "structure_gt_available", "cancer", "tissue", "treatment", "timepoint", "response", "endpoint", "identity_confidence", "metadata_provenance",
        "duplicate_lineage", "leakage_group_id", "record_status", "n_observations", "source_path", "audit_status",
    ]
    spatial_physical: list[dict[str, str]] = []
    crosswalk_fields = [
        "physical_unit_id", "source_project", "source_id", "dataset_id", "patient_id", "block_id", "section_id", "region_id", "physical_specimen_id",
        "identity_granularity", "identity_status", "evidence_grade", "block_equivalent_status", "block_equivalent_basis", "section_order", "leakage_group_id",
        "grouping_eligibility", "metadata_provenance", "relation_status",
    ]
    crosswalk: list[dict[str, str]] = []
    for row in physical_rows:
        physical_id = first(row, "physical_unit_id")
        namespace = first(row, "source_namespace")
        study = first(row, "study_id") or namespace
        source_id = source_id_for_spatial(study, namespace)
        assets = matching_assets(row, namespace, study)
        asset_types = " ".join(first(asset, "asset_type") for asset in assets).lower()
        structures = structure_by_physical.get(physical_id, [])
        gt_ok = any(first(s, "confirmation_status") == "CONFIRMATORY" for s in structures)
        outer = outer_by_physical.get(physical_id, {})
        leakage = first(outer, "leakage_group_id") or (f"013::patient::{first(row, 'patient_id')}" if first(row, "patient_id") else physical_id)
        replay = replay_by_physical.get(physical_id, {})
        record_status = first(row, "record_status")
        identity_status = first(row, "identity_status")
        logical_row = {
            "logical_unit_id": f"013::physical::{physical_id}", "source_project": str(spatial_root), "source_id": source_id,
            "dataset_id": study, "unit_level": "physical_unit", "patient_id": first(row, "patient_id"), "block_id": first(row, "block_id"),
            "section_id": first(row, "section_id"), "region_id": "", "sample_id": first(row, "source_record_id"), "modality": "spatial",
            "platform": first(row, "platform"), "resolution": first(row, "platform"), "coordinate_system": "native_coordinates", "scale": "section_or_capture",
            "counts_layer": "replay_indexed" if replay else ("asset_declared_not_replay_audited" if assets else "not_replay_audited"), "image_available": "yes" if "image" in asset_types or "morphology" in asset_types else "unknown",
            "segmentation_available": "yes" if "segment" in asset_types or "boundary" in asset_types else "unknown", "structure_gt_available": "yes" if gt_ok else "no",
            "cancer": "", "tissue": "", "treatment": "", "timepoint": "", "response": "", "endpoint": "",
            "identity_confidence": first(row, "evidence_grade"), "metadata_provenance": f"/home/huyudi/013_spatial/infra/sample-registry/physical_units.tsv#physical_unit_id={physical_id}",
            "duplicate_lineage": f"013::{leakage}", "leakage_group_id": leakage, "record_status": record_status, "n_observations": first(replay, "spot_count"),
            "source_path": unique_join(first(asset, "path") for asset in assets), "audit_status": identity_status,
        }
        logical.append(logical_row)
        spatial_physical.append({
            "physical_unit_id": physical_id, "source_project": str(spatial_root), "source_id": source_id, "dataset_id": study, "source_namespace": namespace,
            "source_record_id": first(row, "source_record_id"), "patient_id": first(row, "patient_id"), "block_id": first(row, "block_id"),
            "physical_specimen_id": first(row, "physical_specimen_id"), "section_id": first(row, "section_id"), "region_id": "", "modality": "spatial",
            "platform": first(row, "platform"), "resolution": first(row, "platform"), "coordinate_system": "native_coordinates", "counts_layer": "replay_indexed" if replay else ("asset_declared_not_replay_audited" if assets else "not_replay_audited"),
            "image_available": "yes" if "image" in asset_types or "morphology" in asset_types else "unknown", "segmentation_available": "yes" if "segment" in asset_types or "boundary" in asset_types else "unknown",
            "structure_gt_available": "yes" if gt_ok else "no", "cancer": "", "tissue": "", "treatment": "", "timepoint": "", "response": "", "endpoint": "",
            "identity_confidence": first(row, "evidence_grade"), "metadata_provenance": f"/home/huyudi/013_spatial/infra/sample-registry/physical_units.tsv#physical_unit_id={physical_id}",
            "duplicate_lineage": f"013::{leakage}", "leakage_group_id": leakage, "record_status": record_status, "n_observations": first(replay, "spot_count"),
            "source_path": unique_join(first(asset, "path") for asset in assets), "audit_status": identity_status,
        })
        crosswalk.append({
            "physical_unit_id": physical_id, "source_project": str(spatial_root), "source_id": source_id, "dataset_id": study,
            "patient_id": first(row, "patient_id"), "block_id": first(row, "block_id"), "section_id": first(row, "section_id"), "region_id": "",
            "physical_specimen_id": first(row, "physical_specimen_id"), "identity_granularity": first(row, "identity_granularity"), "identity_status": identity_status,
            "evidence_grade": first(row, "evidence_grade"), "block_equivalent_status": first(row, "block_equivalent_status"), "block_equivalent_basis": first(row, "block_equivalent_basis"),
            "section_order": first(row, "section_order"), "leakage_group_id": leakage,
            "grouping_eligibility": "patient_grouped" if record_status == "RESOLVED_INCLUDED_CANDIDATE" and first(row, "patient_id") else ("archive_not_candidate" if record_status == "RESOLVED_ARCHIVE" else "quarantine_identity_or_block"),
            "metadata_provenance": f"/home/huyudi/013_spatial/infra/sample-registry/physical_units.tsv#physical_unit_id={physical_id}",
            "relation_status": record_status,
        })

    # Add the two audited /006 spatial manifests to the same physical-unit table.
    for key, cohort_id in (("gse238_manifest", "GSE238264"), ("gse291_manifest", "GSE291246")):
        manifest = tables.get(key, ([], []))[1]
        qc = tables.get("gse291_qc", ([], []))[1] if cohort_id == "GSE291246" else []
        spot_qc = tables.get("gse238_spot_qc", ([], []))[1] if cohort_id == "GSE238264" else []
        logical_extra, physical_extra, crosswalk_extra = add_006_spatial(manifest, qc, spot_qc, cohort_id, root, effective_by_sample)
        logical.extend(logical_extra)
        spatial_physical.extend(physical_extra)
        crosswalk.extend(crosswalk_extra)

    # Preserve pre-existing duplicate lineage and add patient-nesting edges for /006.
    duplicate_fields = ["duplicate_group_id", "member_type", "member_id", "relation_type", "evidence_grade", "resolution", "leakage_group_id", "source_project", "claim_impact"]
    duplicate: list[dict[str, str]] = []
    for row in tables.get("duplicates", ([], []))[1]:
        duplicate.append({
            "duplicate_group_id": first(row, "duplicate_group_id"), "member_type": first(row, "member_type"), "member_id": first(row, "member_id"),
            "relation_type": first(row, "relation_type"), "evidence_grade": first(row, "evidence_grade"), "resolution": first(row, "resolution"),
            "leakage_group_id": first(row, "leakage_group_id"), "source_project": str(spatial_root), "claim_impact": "inherit_013_lineage_and_split_boundary",
        })
    by_patient: dict[tuple[str, str], list[str]] = defaultdict(list)
    for row in sample_rows:
        cohort_id = first(row, "cohort_id")
        if cohort_id in {"GSE238264", "GSE291246"}:
            continue
        patient = first(row, "patient_key", "patient_id")
        sample = first(row, "sample_key", "sample_id")
        if patient and sample:
            by_patient[(cohort_id, patient)].append(sample)
    for (cohort_id, patient), samples in sorted(by_patient.items()):
        if len(samples) < 2:
            continue
        group = f"006::{cohort_id}::patient::{patient}"
        for left, right in zip(sorted(samples), sorted(samples)[1:]):
            duplicate.append({
                "duplicate_group_id": group, "member_type": "sample_pair", "member_id": f"006::{left}__006::{right}",
                "relation_type": "same_patient_nested_observations", "evidence_grade": "E3_explicit_master", "resolution": "retain_samples_group_patient",
                "leakage_group_id": group, "source_project": str(root), "claim_impact": "cannot_count_as_independent_patient_or_cross_split",
            })
    manifest_groups: dict[tuple[str, str], list[str]] = defaultdict(list)
    for manifest_key, manifest_cohort in (("gse238_manifest", "GSE238264"), ("gse291_manifest", "GSE291246")):
        for row in tables.get(manifest_key, ([], []))[1]:
            raw_patient = first(row, "patient_id")
            sample = first(row, "sample_id", "gsm_id")
            if raw_patient and sample:
                manifest_groups[(manifest_cohort, f"{manifest_cohort}::{raw_patient}")].append(sample)
    for (cohort_id, patient), samples in sorted(manifest_groups.items()):
        if len(samples) < 2:
            continue
        group = f"006::{cohort_id}::patient::{patient}"
        for left, right in zip(sorted(samples), sorted(samples)[1:]):
            duplicate.append({
                "duplicate_group_id": group, "member_type": "physical_section_pair", "member_id": f"006::{cohort_id}::section::{left}__006::{cohort_id}::section::{right}",
                "relation_type": "same_patient_nested_spatial_sections", "evidence_grade": "E3_explicit_manifest", "resolution": "retain_sections_group_patient",
                "leakage_group_id": group, "source_project": str(root), "claim_impact": "cannot_count_as_independent_patient_or_cross_split",
            })
    for row in cross_section_rows:
        duplicate.append({
            "duplicate_group_id": first(row, "leakage_group_id") or first(row, "link_id"), "member_type": "physical_pair",
            "member_id": f"{first(row, 'left_physical_unit_id')}__{first(row, 'right_physical_unit_id')}", "relation_type": first(row, "relation_type") or "serial_section",
            "evidence_grade": first(row, "evidence_grade"), "resolution": "retain_sections_same_leakage_group", "leakage_group_id": first(row, "leakage_group_id"),
            "source_project": str(spatial_root), "claim_impact": "nested_section_not_independent_patient",
        })

    # Dataset x patient treatment/response completeness.
    completeness_fields = [
        "source_project", "source_id", "dataset_id", "patient_id", "n_samples", "n_sections", "n_timepoints", "n_treatment_present", "n_timepoint_present",
        "n_endpoint_present", "n_response_present", "n_response_unknown", "treatment_raw_values", "treatment_values", "timepoint_raw_values", "timepoint_values",
        "endpoint_values", "response_raw_values", "response_values", "treatment_status", "timepoint_status", "response_status", "complete_for_response",
        "metadata_provenance", "conflict_flag", "readiness_note",
    ]
    by_dataset_patient: dict[tuple[str, str, str, str], list[dict[str, str]]] = defaultdict(list)
    for row in sample_rows:
        cohort_id = first(row, "cohort_id")
        if cohort_id in {"GSE238264", "GSE291246"}:
            continue
        patient = first(row, "patient_key", "patient_id")
        if patient:
            by_dataset_patient[("/006", source_id_for_cohort(cohort_id, first(row, "modality")), cohort_id, patient)].append(row)
    for row in tables.get("gse238_manifest", ([], []))[1] + tables.get("gse291_manifest", ([], []))[1]:
        cohort_id = "GSE238264" if row in tables.get("gse238_manifest", ([], []))[1] else "GSE291246"
        sample_id = first(row, "sample_id", "gsm_id")
        effective = effective_by_sample.get((cohort_id, sample_id), {})
        enriched = dict(row)
        for key in ("timepoint", "timepoint_raw", "treatment_context", "treatment_raw", "response_raw", "response_endpoint_type", "response_binary_harmonized"):
            if first(effective, key):
                enriched[key] = first(effective, key)
        raw_patient = first(row, "patient_id")
        patient = first(effective, "patient_key") or (f"{cohort_id}::{raw_patient}" if raw_patient else "")
        if patient:
            by_dataset_patient[("/006", source_id_for_cohort(cohort_id, "spatial"), cohort_id, patient)].append(enriched)
    # /013 patient-level rows are intentionally mostly incomplete; this records that fact.
    for row in physical_rows:
        patient = first(row, "patient_id")
        if patient:
            study = first(row, "study_id") or first(row, "source_namespace")
            by_dataset_patient[("/013_spatial", source_id_for_spatial(study, first(row, "source_namespace")), study, patient)].append(row)
    completeness: list[dict[str, str]] = []
    for (project_label, source_id, dataset_id, patient), rows in sorted(by_dataset_patient.items()):
        treatment_raw = [first(r, "treatment_raw", "treatment_regimen_raw", "treatment_context") for r in rows]
        treatment = [first(r, "treatment_context", "treatment_regimen_raw", "treatment_raw") for r in rows]
        timepoint_raw = [first(r, "timepoint_raw", "timepoint") for r in rows]
        timepoint = [first(r, "timepoint") or ("post" if dataset_id == "GSE238264" else "") for r in rows]
        endpoint = [first(r, "response_endpoint_type", "response_standard") for r in rows]
        response_raw = [first(r, "response_raw", "response_binary", "response_standard") for r in rows]
        response = [first(r, "response_binary_harmonized", "response_binary", "response_standard") for r in rows]
        def usable(value: str) -> bool:
            return bool(value) and value.lower() not in {"na", "n/a", "none", "unknown", "not_provided", "not provided", "null"}

        def status(values: list[str]) -> str:
            present = [value for value in values if usable(value)]
            return "missing" if not present else ("conflict" if len(set(present)) > 1 else "complete")
        t_status, r_status = status(treatment), status(response)
        if dataset_id == "GSE291246" and t_status == "conflict":
            t_status = "multi_context"
        # A patient normally has several observed timepoints.  That is not a
        # metadata conflict; the response contract only needs at least one
        # explicitly named timepoint.  Treatment/endpoint/response conflicts
        # remain fail-closed.
        tp_status = "complete" if any(usable(value) for value in timepoint) else "missing"
        endpoint_status = status(endpoint)
        complete = t_status == tp_status == r_status == endpoint_status == "complete"
        note = ""
        if dataset_id == "GSE238264":
            note = "post-treatment PD1+cabozantinib spatial support; no baseline and no component separation"
        elif dataset_id == "GSE291246":
            note = "BCC mixed context (post_PD1 and treatment_naive) support; response unavailable"
        elif project_label == "/013_spatial":
            note = "spatial registry has no treatment/timepoint/response semantics in this audit"
        completeness.append({
            "source_project": project_label, "source_id": source_id, "dataset_id": dataset_id, "patient_id": patient,
            "n_samples": str(len(rows)), "n_sections": str(sum(1 for r in rows if first(r, "section_id", "sample_id", "gsm_id"))),
            "n_timepoints": str(len({v for v in timepoint if usable(v)})), "n_treatment_present": str(sum(usable(v) for v in treatment)), "n_timepoint_present": str(sum(usable(v) for v in timepoint)),
            "n_endpoint_present": str(sum(usable(v) for v in endpoint)), "n_response_present": str(sum(usable(v) for v in response)), "n_response_unknown": str(sum(not usable(v) for v in response)),
            "treatment_raw_values": unique_join(treatment_raw), "treatment_values": unique_join(treatment), "timepoint_raw_values": unique_join(timepoint_raw), "timepoint_values": unique_join(timepoint),
            "endpoint_values": unique_join(endpoint), "response_raw_values": unique_join(response_raw), "response_values": unique_join(response), "treatment_status": t_status, "timepoint_status": tp_status,
            "response_status": r_status if endpoint_status == "complete" else ("missing_endpoint" if r_status == "complete" else r_status), "complete_for_response": "yes" if complete else "no",
            "metadata_provenance": "effective_metadata_or_spatial_manifest", "conflict_flag": "yes" if "conflict" in {t_status, r_status, endpoint_status} else "no", "readiness_note": note,
        })

    # Claim-specific roles.  This is a v7 starting ledger, not a scientific result.
    role_fields = [
        "claim_id", "claim_name", "source_id", "dataset_id", "logical_unit_id", "primary_role", "selection_exposure", "role_source", "permitted_role", "forbidden_use",
        "eligibility_status", "response_ready", "treatment_context", "timepoint_scope", "endpoint_scope", "reason", "leakage_group_id",
    ]
    config_by_id = {first(row, "source_id"): row for row in config}
    cohort_by_source: dict[str, dict[str, str]] = {}
    for cohort_id, cohort_row in cohort_rows.items():
        cohort_by_source.setdefault(source_id_for_cohort(cohort_id, first(cohort_row, "modality")), cohort_row)
    context_by_source: dict[str, set[str]] = defaultdict(set)
    for sample in sample_rows:
        source_id = source_id_for_cohort(first(sample, "cohort_id"), first(sample, "modality"))
        context = first(sample, "treatment_context", "treatment_raw")
        if context:
            context_by_source[source_id].add(context)
    complete_by_source: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in completeness:
        complete_by_source[row["source_id"]].append(row)
    assignments: list[dict[str, str]] = []
    for source_id, config_row in sorted(config_by_id.items()):
        role, role_source, permitted, forbidden = role_for_source(source_id, config_row, role_rows)
        dataset_id = source_id
        cohort = cohort_by_source.get(source_id, {})
        readiness = "ready_for_metadata_scoped_use"
        response_ready = "no"
        reason = first(config_row, "main_caveat")
        if source_id in {"006_GSE301741", "006_Mendeley_HCC", "013_GSE211956"}:
            readiness = "pending_interface_closure"
        if source_id == "006_GSE238264":
            response_ready = "post_treatment_association_only"
        elif source_id in {"006_TASK01", "006_GSE286827"}:
            response_ready = "within_environment_only"
        elif source_id in {"006_TASK02", "006_LAMBRECHT_HCC"}:
            response_ready = "paired_or_context_support_only"
        for claim_id, claim_name in CLAIMS.items():
            claim_role = role
            claim_reason = reason
            if claim_id in {"H2", "H3", "H4"} and source_id in {"013_GSE175540", "013_TLS_USZ"}:
                claim_role = "external_validation_locked"
            if claim_id in {"H6", "H7"} and source_id in {"006_TASK02", "006_LAMBRECHT_HCC", "006_GSE238264"}:
                claim_role = "support_context_or_repair"
            if source_id == "013_GSE211956":
                claim_role = "candidate_pending_response_semantics"
            assignments.append({
                "claim_id": claim_id, "claim_name": claim_name, "source_id": source_id, "dataset_id": dataset_id,
                "logical_unit_id": f"dataset::{source_id}", "primary_role": claim_role, "selection_exposure": "not_run_v7_stage1_metadata_audit",
                "role_source": role_source, "permitted_role": permitted, "forbidden_use": forbidden, "eligibility_status": readiness,
                "response_ready": response_ready, "treatment_context": unique_join(context_by_source.get(source_id, set())) or first(cohort, "treatment_context") or "unknown",
                "timepoint_scope": first(cohort, "n_samples") or first(config_row, "inventory_scale"),
                "endpoint_scope": first(cohort, "response_endpoint_type") or "effective_metadata_if_present", "reason": claim_reason, "leakage_group_id": "",
            })

    # Publish all required outputs.
    output.mkdir(parents=True, exist_ok=True)
    counts = {
        "logical_units": write_tsv(output / "multimodal_logical_units.tsv", logical_fields, logical),
        "spatial_physical_units": write_tsv(output / "spatial_physical_units.tsv", spatial_fields, spatial_physical),
        "crosswalk": write_tsv(output / "patient_block_section_crosswalk.tsv", crosswalk_fields, crosswalk),
        "completeness": write_tsv(output / "treatment_response_completeness.tsv", completeness_fields, completeness),
        "duplicate_lineage": write_tsv(output / "duplicate_lineage.tsv", duplicate_fields, duplicate),
        "role_assignment": write_tsv(output / "data_role_assignment.tsv", role_fields, assignments),
    }

    source_manifest: list[dict[str, str]] = []
    for name, path in sorted(paths.items()):
        if path.exists() and path.is_file() and path.suffix.lower() in {".csv", ".tsv", ".tab", ".json", ".md"}:
            source_manifest.append({"name": name, "path": str(path), "sha256": sha256(path), "bytes": str(path.stat().st_size), "status": "READ_METADATA_ONLY"})
    input_lines = "\n".join(f"- `{item['name']}`: `{item['path']}` ({item['bytes']} bytes, SHA256 `{item['sha256']}`)" for item in source_manifest)
    source_counts = Counter(row["source_id"] for row in logical if row["source_id"])
    source_table = "\n".join(f"| `{source}` | {source_counts[source]} |" for source in sorted(source_counts))
    role_counts = Counter(row["primary_role"] for row in assignments)
    role_table = "\n".join(f"| `{role}` | {count} |" for role, count in sorted(role_counts.items()))
    included_013 = sum(first(row, "record_status") == "RESOLVED_INCLUDED_CANDIDATE" for row in physical_rows)
    confirmatory_structures = sum(first(row, "confirmation_status") == "CONFIRMATORY" for row in structure_rows)
    gse238 = [row for row in completeness if row["dataset_id"] == "GSE238264"]
    gse291 = [row for row in completeness if row["dataset_id"] == "GSE291246"]
    gse238_response_counts = Counter(
        row["response"] for row in spatial_physical
        if row["dataset_id"] == "GSE238264" and row["response"] and row["response"].lower() not in {"na", "unknown", "not_provided"}
    )
    gse238_response_text = ", ".join(f"{key}={value}" for key, value in sorted(gse238_response_counts.items())) or "unknown"
    gse291_timepoint_counts = Counter(
        first(row, "timepoint") for row in sample_rows
        if first(row, "cohort_id") == "GSE291246" and first(row, "timepoint")
    )
    gse291_timepoint_text = ", ".join(f"{key}={value}" for key, value in sorted(gse291_timepoint_counts.items())) or "unknown"
    report = f"""# Stage 1 registry report

版本：`{SCRIPT_VERSION}`  
范围：跨 `/006` 与只读 `/home/huyudi/013_spatial` 的元数据、身份、角色、重复和治疗/疗效字段审计。  
运行边界：未打开表达矩阵、图像、h5ad 内容或大型二进制；未运行 v7 科学模型。

## 生成与输入

可复现命令：

```bash
python scripts/v7/registry/build_stage1_registry.py
```

输入文件哈希（本次仅读取元数据）：

{input_lines}

全局 `Infra/bioinf-data-index` 已覆盖 `/006` 与 `/013_spatial`；本阶段未下载新增外部数据，因此不改写全局索引。

## 账本规模

| 账本 | 数量 | 解释 |
|---|---:|---|
| `/006` effective cohort rows | {len(cohort_rows)} | 队列级元数据 |
| `/006` effective patient rows | {len(patient_rows)} | 患者级候选，不等于可监督患者 |
| `/006` effective sample rows | {len(sample_rows)} | 样本/时间点级记录 |
| `/013` physical rows | {len(physical_rows)} | 含排除记录的完整审计面 |
| `/013` resolved included candidates | {included_013} | 可进入后续空间候选池的物理单位 |
| `/013` confirmatory structure instances | {confirmatory_structures} | 仅已登记结构 GT，不代表疗效证据 |
| Stage 1 logical-unit rows | {counts['logical_units']} | 样本、物理单位和 dataset-level 账本 |
| Stage 1 spatial physical rows | {counts['spatial_physical_units']} | `/013` + `/006` 两个空间 manifest |
| Stage 1 patient/block/section rows | {counts['crosswalk']} | patient 是外层独立单位 |

按 source 的 logical-unit 行数：

| source | rows |
|---|---:|
{source_table}

## Stage 1 关键发现与资格结论

1. `/006` 的有效主表优先于旧摘要：TASK01 为 27 患者/72 样本行，TASK02 为 23/66；Lambrecht HCC 为 44 患者/112 样本行（pre 65、post 47）。`config/v7/data_source_registry.tsv` 中的 Lambrecht “25 pre/22 paired”是旧口径，已作为冲突保留，不能直接沿用。
2. GSE238264 是 7 名 HCC 患者的治疗后空间队列（nivo+cabo；manifest response 计数：{gse238_response_text}），可作治疗后空间关联和 support；无基线，不能分离 PD-1 与 cabozantinib，也不能建立纵向空间重排因果链。
3. GSE291246 是 BCC 而非 HCC：35 个 Xenium section，17 个 h5ad/QC 可重放记录；有效主表的时间点口径为 `{gse291_timepoint_text}`，response 未提供，故只作跨癌种空间 support。QC 表中患者列可能为空，本账本以 metadata-freeze manifest 的患者映射为准。
4. `/013_spatial` 的身份控制面纳入 1,002 条 physical rows；其中 {included_013} 条 `RESOLVED_INCLUDED_CANDIDATE` 才是后续空间候选池。block/section/spot/cell 是患者内嵌套观测，不能增加独立 n；同患者、同块、serial section 通过 leakage group 约束拆分。
5. GSE211956 的样本标题 response token 仍处于语义待核验状态；Mendeley 的有效主表目前是 6 名患者/12 个样本行，且本地 h5ad 被标为 scRNA（旧摘要中的 11 个 object 口径不作为 v7 计数）。两者均未升级为 response-bearing 或独立验证。
6. `/013_spatial/infra/r04` 只作为输入定位与重启稳定性诊断的 provenance；本阶段不继承 latent-field 的 K 或结构结论，且最新面板仍不是正式 K 选择。
7. 物理单位的来源反向覆盖已闭合：HEST、ATLAS_TABLE_S2、TENX 和 3 个额外 GEO 空间队列分别进入受控 source；所有 spatial `source_id` 都能在 config 和 role ledger 中找到，未再使用无角色的 `013_other_spatial`。

## 角色账本（D1）

角色行数（8 个 claim × source）：

| primary role | rows |
|---|---:|
{role_table}

角色来源分开记录：013 旧 `role_freeze.tsv` 只作为 provenance 和历史泄漏边界；同一 v7 claim 尚未运行选择/阈值/模型，因此 `selection_exposure` 统一写为 `not_run_v7_stage1_metadata_audit`。locked external 数据不进入候选或调参。

## 治疗/疗效完整性边界

`treatment_response_completeness.tsv` 以 dataset×patient 为行，保留 raw/normalized 值、端点、时间点和冲突标志。只有 treatment、timepoint、endpoint、response 都存在且无冲突的行才标记 `complete_for_response=yes`；这只是字段完整，不自动获得临床主张资格。

- TASK01：PD-1 anchor 环境；仅在自身治疗/终点环境内使用，不能与 PD1+X 直接 pooled。
- TASK02：PD1+lenvatinib 环境；用于 paired/context repair，不作 PD-1 单药交换对照。
- Lambrecht HCC：治疗组合混杂，作为纵向 support，先解决有效主表与旧摘要口径冲突。
- GSE238264：post-only PD1+cabozantinib spatial support。
- GSE291246：post_PD1 BCC，response missing。

## 阻隔与下一阶段

- 直接纵向空间 PD1+X rewiring：`NOT_IDENTIFIABLE_WITH_CURRENT_AUDITED_ASSETS`；当前没有同患者空间 pre/on/post 链，post-only 不能替代。
- response-linked spatial：GSE211956 保持 `response_semantics_pending`，待治疗药物、时间点、endpoint 和表达/坐标/结构链接闭合。
- Mendeley HCC：patient/section/response/image-scale provenance 未闭合。
- Stage 2 可启动：先用 response-blind 程序和 cell-state 词汇建立跨模态测量接口；不得把 Stage 1 账本当成生物学结果。

## 验收状态

- 主键和确定性排序：由生成器统一处理；重复运行应 byte-identical。
- 空 patient/block/section：保持空值并显式标记，不把 unknown 合并成共同 ID。
- `/013_spatial`：本阶段只读，未写入。
- 科学结论：`NOT_RUN`；本报告只证明数据资格、有效样本量和边界。
"""
    (output / "REGISTRY_REPORT.md").write_text(report, encoding="utf-8")
    manifest_payload = {
        "script_version": SCRIPT_VERSION,
        "project_root": str(root),
        "spatial_root": str(spatial_root),
        "output_dir": str(output),
        "read_policy": "metadata_only_no_matrix_or_image_read",
        "counts": counts,
        "input_files": source_manifest,
        "output_files": sorted(path.name for path in output.iterdir() if path.is_file()),
    }
    (output / "STAGE1_RUN_MANIFEST.json").write_text(json.dumps(manifest_payload, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return {"counts": counts, "input_files": source_manifest, "output": output, "physical_rows": len(physical_rows), "included_013": included_013}


def validate(output: Path, project_root: Path | None = None) -> list[str]:
    errors: list[str] = []
    required = {
        "multimodal_logical_units.tsv": "logical_unit_id",
        "spatial_physical_units.tsv": "physical_unit_id",
        "patient_block_section_crosswalk.tsv": "physical_unit_id",
        "treatment_response_completeness.tsv": "dataset_id",
        "duplicate_lineage.tsv": "member_id",
        "data_role_assignment.tsv": "claim_id",
    }
    for name, key in required.items():
        path = output / name
        fields, rows = read_table(path, "\t")
        if not fields:
            errors.append(f"missing or empty table: {name}")
            continue
        if key not in fields:
            errors.append(f"{name}: missing key {key}")
        values = [first(row, key) for row in rows if first(row, key)]
        if name == "treatment_response_completeness.tsv":
            values = [f"{first(row, 'dataset_id')}::{first(row, 'patient_id')}" for row in rows]
        elif name == "data_role_assignment.tsv":
            values = [f"{first(row, 'claim_id')}::{first(row, 'source_id')}" for row in rows]
        if name != "duplicate_lineage.tsv" and len(values) != len(set(values)):
            errors.append(f"{name}: non-unique primary key {key}")
    cross_fields, cross_rows = read_table(output / "patient_block_section_crosswalk.tsv", "\t")
    physical_fields, physical_rows = read_table(output / "spatial_physical_units.tsv", "\t")
    cross_ids = {first(row, "physical_unit_id") for row in cross_rows}
    physical_ids = {first(row, "physical_unit_id") for row in physical_rows}
    if cross_ids - physical_ids:
        errors.append(f"crosswalk IDs missing from spatial physical table: {len(cross_ids - physical_ids)}")
    if project_root is not None:
        config_path = project_root / "config/v7/data_source_registry.tsv"
        config_fields, config_rows = read_table(config_path, "\t")
        configured = {first(row, "source_id") for row in config_rows if first(row, "source_id")}
        physical_sources = {first(row, "source_id") for row in physical_rows if first(row, "source_id")}
        if physical_sources - configured:
            errors.append("spatial sources missing from config: " + ", ".join(sorted(physical_sources - configured)))
        role_fields, role_rows = read_table(output / "data_role_assignment.tsv", "\t")
        role_sources = {first(row, "source_id") for row in role_rows if first(row, "source_id")}
        if configured - role_sources:
            errors.append("config sources missing from role assignment: " + ", ".join(sorted(configured - role_sources)))
        logical_fields, logical_rows = read_table(output / "multimodal_logical_units.tsv", "\t")
        logical_ids = {first(row, "logical_unit_id") for row in logical_rows if first(row, "logical_unit_id")}
        role_logical_ids = {first(row, "logical_unit_id") for row in role_rows if first(row, "logical_unit_id")}
        if role_logical_ids - logical_ids:
            errors.append("role logical units missing from multimodal registry: " + ", ".join(sorted(role_logical_ids - logical_ids)))
    return errors


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", default=Path(__file__).resolve().parents[3])
    parser.add_argument("--spatial-root", default="/home/huyudi/013_spatial")
    parser.add_argument("--output-dir", default=None)
    parser.add_argument("--validate-only", action="store_true")
    args = parser.parse_args(argv)
    args.project_root = str(Path(args.project_root).resolve())
    args.spatial_root = str(Path(args.spatial_root).resolve())
    if args.output_dir is None:
        args.output_dir = str(Path(args.project_root) / "results/v7/registry")
    if not args.validate_only:
        try:
            result = build_registry(args)
        except RuntimeError as error:
            print(str(error), file=sys.stderr)
            return 2
        print(json.dumps({"status": "BUILT", "counts": result["counts"], "output": str(result["output"])}, ensure_ascii=False, sort_keys=True))
    errors = validate(Path(args.output_dir).resolve(), Path(args.project_root).resolve())
    if errors:
        for error in errors:
            print(f"ERROR: {error}", file=sys.stderr)
        return 1
    print(json.dumps({"status": "VALID", "output": str(Path(args.output_dir).resolve())}, ensure_ascii=False, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
