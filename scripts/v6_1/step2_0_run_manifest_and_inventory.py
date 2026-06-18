from __future__ import annotations

import argparse
import csv
import hashlib
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, Iterable, List, Sequence, Tuple


ROOT = Path(__file__).resolve().parents[2]
STEP2_ROOT = ROOT / "results" / "v6_1" / "step2"
SNAPSHOT_DIR = ROOT / "results" / "v6_1" / "data_pool" / "snapshots" / "current_step2_input"

REQUIRED_METADATA_FIELDS = [
    "cohort_id",
    "patient_id",
    "sample_id",
    "split",
    "disease",
    "tissue_source",
    "timepoint",
    "treatment_context",
]

FIELD_ALIASES: dict[str, list[str]] = {
    "split": ["split_label", "split_group"],
    "timepoint": ["timepoint_normalized", "timepoint_raw", "timepoint_standardized", "response_source_timepoint"],
    "treatment_context": ["treatment_context_primary", "treatment_context_from_flags"],
}

RESPONSE_COLUMNS = ["response_raw", "response_binary", "response_strict", "response_broad"]

SUBDIRS = [
    "00_manifest",
    "01_mapping_audit",
    "02_gene_standardization",
    "03_qc",
    "04_annotation",
    "05_myeloid_qc",
    "06_fraction",
    "07_pseudobulk",
    "08_signature_pathway_tf",
    "09_feature_matrix",
    "10_reports",
]


def now_iso() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat()


def _sha256_file(path_str: str) -> Tuple[str, str]:
    path = Path(path_str)
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(262144), b""):
            h.update(chunk)
    return (path_str, h.hexdigest())


def sha256_file(path: Path) -> str:
    return _sha256_file(str(path))[1]


def sha256_files_parallel(paths: List[Path], workers: int | None = None) -> Dict[str, str]:
    if workers is None:
        import os
        workers = max(1, os.cpu_count() or 1)
    with ThreadPoolExecutor(max_workers=workers) as executor:
        results = list(executor.map(_sha256_file, [str(p) for p in paths]))
    return {path: h for path, h in results}


def read_csv_dicts(path: Path, delimiter: str | None = ",") -> List[Dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        if delimiter is None:
            sample = handle.read(4096)
            handle.seek(0)
            try:
                dialect = csv.Sniffer().sniff(sample, delimiters=",\t")
                delimiter = dialect.delimiter
            except Exception:
                delimiter = "\t" if path.suffix.lower() == ".tsv" else ","
        return list(csv.DictReader(handle, delimiter=delimiter))


def csv_write(path: Path, fieldnames: Sequence[str], rows: Iterable[Dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def write_yaml(path: Path, lines: List[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def discover_files(root: Path, workers: int | None = None) -> List[Dict[str, str]]:
    files_to_hash: List[Path] = []
    file_roles: List[Tuple[Path, str]] = []

    # snapshot / legacy step1
    snapshot = root / "results" / "v6_1" / "data_pool" / "snapshots" / "current_step2_input"
    legacy = root / "results" / "v6_1" / "step1"

    refs = [
        (snapshot / "step1_exports" / "frozen_cohort_inclusion_v6_1.csv", "frozen_cohort_inclusion"),
        (legacy / "frozen_cohort_inclusion_v6_1.csv", "frozen_cohort_inclusion"),
        (snapshot / "step1_exports" / "frozen_patient_split_v6_1.csv", "frozen_patient_split"),
        (legacy / "frozen_patient_split_v6_1.csv", "frozen_patient_split"),
        (snapshot / "step1_exports" / "sample_metadata_master_v6_1.broad_response.csv", "sample_metadata_broad_response"),
        (legacy / "sample_metadata_master_v6_1.broad_response.csv", "sample_metadata_broad_response"),
        (snapshot / "step1_exports" / "patient_metadata_master_v6_1.broad_response.csv", "patient_metadata_broad_response"),
        (legacy / "patient_metadata_master_v6_1.broad_response.csv", "patient_metadata_broad_response"),
        (snapshot / "snapshot.yaml", "snapshot_manifest"),
        (snapshot / "h5ad_inventory.csv", "h5ad_inventory"),
    ]
    for p, role in refs:
        if p.exists():
            files_to_hash.append(p)
            file_roles.append((p, role))

    # raw h5ad
    h5ad_dir = root / "data" / "processed" / "srt" / "raw"
    if h5ad_dir.exists():
        for f in sorted(h5ad_dir.glob("*.h5ad")):
            files_to_hash.append(f)
            file_roles.append((f, "raw_h5ad"))

    # sidecars
    sidecar_dir = root / "mvp" / "sidecars"
    if sidecar_dir.exists():
        for cohort_dir in sorted(sidecar_dir.iterdir()):
            if not cohort_dir.is_dir():
                continue
            for f in sorted(cohort_dir.iterdir()):
                if f.name == "obs.parquet":
                    files_to_hash.append(f)
                    file_roles.append((f, "sidecar_obs"))
                elif f.name == "var.parquet":
                    files_to_hash.append(f)
                    file_roles.append((f, "sidecar_var"))
                elif f.name == "sample_summary.csv":
                    files_to_hash.append(f)
                    file_roles.append((f, "sidecar_sample_summary"))

    # configs
    config_dir = root / "config"
    if config_dir.exists():
        for f in sorted(config_dir.rglob("*")):
            if f.is_file():
                files_to_hash.append(f)
                file_roles.append((f, "project_config"))

    # requirements
    for req in sorted(root.glob("requirements*.txt")):
        files_to_hash.append(req)
        file_roles.append((req, "project_config"))
    pyright = root / "mvp" / "pyrightconfig.json"
    if pyright.exists():
        files_to_hash.append(pyright)
        file_roles.append((pyright, "project_config"))

    # parallel hash
    hash_map = sha256_files_parallel(files_to_hash, workers=workers)

    entries: List[Dict[str, str]] = []
    for path, role in file_roles:
        stat = path.stat()
        entries.append({
            "file_path": str(path),
            "file_exists": "True",
            "file_size": str(stat.st_size),
            "modified_time": datetime.fromtimestamp(stat.st_mtime, tz=timezone.utc).isoformat(),
            "sha256": hash_map.get(str(path), ""),
            "file_type": path.suffix.lstrip(".") or "unknown",
            "expected_role": role,
        })

    return entries


def check_metadata_fields(sample_path: Path, patient_path: Path) -> List[Dict[str, str]]:
    sample_header = set()
    patient_header = set()

    if sample_path.exists():
        with sample_path.open("r", encoding="utf-8-sig", newline="") as f:
            sample_header = set(next(csv.reader(f)))
    if patient_path.exists():
        with patient_path.open("r", encoding="utf-8-sig", newline="") as f:
            patient_header = set(next(csv.reader(f)))

    results: List[Dict[str, str]] = []
    for field in REQUIRED_METADATA_FIELDS:
        alias_list = FIELD_ALIASES.get(field, [])
        in_sample = field in sample_header or any(a in sample_header for a in alias_list)
        in_patient = field in patient_header or any(a in patient_header for a in alias_list)
        if in_sample or in_patient:
            status = "present"
            if (field not in sample_header and field not in patient_header):
                status = "alias_only"
        else:
            status = "missing"
        results.append({
            "field": field,
            "status": status,
            "in_sample": str(in_sample),
            "in_patient": str(in_patient),
            "aliases_checked": ";".join(alias_list) if alias_list else "",
        })

    # response columns presence check (only presence, no values)
    for rc in RESPONSE_COLUMNS:
        in_s = rc in sample_header
        in_p = rc in patient_header
        results.append({
            "field": rc,
            "status": "present" if (in_s or in_p) else "missing",
            "in_sample": str(in_s),
            "in_patient": str(in_p),
            "aliases_checked": "",
        })

    return results


def run_step2_manifest_inventory(
    root_dir: Path | str = ROOT,
    run_id: str | None = None,
    snapshot_path: Path | str | None = None,
    workers: int | None = None,
) -> Dict[str, object]:
    root = Path(root_dir)
    if run_id is None:
        now = datetime.now(timezone.utc)
        run_id = f"step2_v6_1_{now.strftime('%m%d_%H%M')}"

    run_root = root / "results" / "v6_1" / "step2" / run_id
    run_root.mkdir(parents=True, exist_ok=True)

    # create subdirs
    for sub in SUBDIRS:
        (run_root / sub).mkdir(parents=True, exist_ok=True)
        (run_root / sub / "checkpoints").mkdir(parents=True, exist_ok=True)
        (run_root / sub / "logs").mkdir(parents=True, exist_ok=True)
    (run_root / "03_qc" / "qc_figures").mkdir(parents=True, exist_ok=True)

    # resolve refs
    snapshot = root / "results" / "v6_1" / "data_pool" / "snapshots" / "current_step2_input"
    legacy = root / "results" / "v6_1" / "step1"

    if snapshot_path is not None:
        sp = Path(snapshot_path)
    elif (snapshot / "snapshot.yaml").exists():
        sp = snapshot / "snapshot.yaml"
    else:
        sp = None

    def pick(pref_snap: Path, pref_legacy: Path) -> Path:
        if pref_snap.exists():
            return pref_snap
        return pref_legacy

    cohort_freeze = pick(snapshot / "step1_exports" / "frozen_cohort_inclusion_v6_1.csv", legacy / "frozen_cohort_inclusion_v6_1.csv")
    patient_split = pick(snapshot / "step1_exports" / "frozen_patient_split_v6_1.csv", legacy / "frozen_patient_split_v6_1.csv")
    sample_meta = pick(snapshot / "step1_exports" / "sample_metadata_master_v6_1.broad_response.csv", legacy / "sample_metadata_master_v6_1.broad_response.csv")
    patient_meta = pick(snapshot / "step1_exports" / "patient_metadata_master_v6_1.broad_response.csv", legacy / "patient_metadata_master_v6_1.broad_response.csv")

    # file inventory
    inventory = discover_files(root, workers=workers)
    h5ad_count = sum(1 for e in inventory if e["expected_role"] == "raw_h5ad")
    sidecar_obs_count = sum(1 for e in inventory if e["expected_role"] == "sidecar_obs")
    sidecar_var_count = sum(1 for e in inventory if e["expected_role"] == "sidecar_var")

    inventory_path = run_root / "00_manifest" / "input_file_inventory.csv"
    csv_write(inventory_path, ["file_path", "file_exists", "file_size", "modified_time", "sha256", "file_type", "expected_role"], inventory)

    # field check
    field_checks = check_metadata_fields(sample_meta, patient_meta)
    field_check_path = run_root / "00_manifest" / "required_field_check.csv"
    csv_write(field_check_path, ["field", "status", "in_sample", "in_patient", "aliases_checked"], field_checks)

    missing_required = [f["field"] for f in field_checks if f["status"] == "missing" and f["field"] in REQUIRED_METADATA_FIELDS]

    # blocking logic
    blocking = False
    blocking_reasons: List[str] = []
    if not sample_meta.exists():
        blocking = True
        blocking_reasons.append("sample_metadata_master missing")
    if not patient_meta.exists():
        blocking = True
        blocking_reasons.append("patient_metadata_master missing")
    if not patient_split.exists():
        blocking = True
        blocking_reasons.append("patient_split missing")
    if not cohort_freeze.exists():
        blocking = True
        blocking_reasons.append("cohort_freeze missing")
    if h5ad_count == 0 and sidecar_obs_count == 0:
        blocking = True
        blocking_reasons.append("no raw_h5ad and no sidecar_obs")

    # scope lock markdown
    scope_lines = [
        "# Step2 Scope Lock",
        "",
        f"- **run_id**: {run_id}",
        f"- **created_at**: {now_iso()}",
        "",
        "## Frozen Rules",
        "",
        "- do_not_use_response_for_analysis: true",
        "- per_sample_qc_first: true",
        "- no_global_scvi_in_main_step2: true",
        "- scvi_optional_only_after_core_features_frozen: true",
        "- raw_count_pseudobulk_required: true",
        "- logcpm_pseudobulk_required: true",
        "- missing_gene_in_union_matrix: NA_not_zero",
        "",
        "## What Step2 Does NOT Do",
        "",
        "- No responder vs non-responder comparison",
        "- No model training",
        "- No differential analysis",
        "- No module discovery",
        "- No mechanistic conclusion",
        "",
        "## Blocking Status",
        "",
        f"- blocking_missing_input: {blocking}",
        f"- blocking_reasons: {', '.join(blocking_reasons) if blocking_reasons else 'none'}",
    ]
    scope_path = run_root / "00_manifest" / "step2_scope_lock.md"
    scope_path.write_text("\n".join(scope_lines) + "\n", encoding="utf-8")

    # manifest yaml
    manifest_lines = [
        f"run_id: {run_id}",
        f"created_at: {now_iso()}",
        f"step2_root: {run_root}",
        f"data_pool_snapshot_ref: {sp if sp else 'none'}",
        f"cohort_freeze_ref: {cohort_freeze}",
        f"patient_split_ref: {patient_split}",
        f"sample_metadata_ref: {sample_meta}",
        f"patient_metadata_ref: {patient_meta}",
        f"h5ad_inventory_count: {h5ad_count}",
        f"sidecar_obs_count: {sidecar_obs_count}",
        f"sidecar_var_count: {sidecar_var_count}",
        "global_rules:",
        "  do_not_use_response_for_analysis: true",
        "  per_sample_qc_first: true",
        "  no_global_scvi_in_main_step2: true",
        "  scvi_optional_only_after_core_features_frozen: true",
        "  raw_count_pseudobulk_required: true",
        "  logcpm_pseudobulk_required: true",
        "  missing_gene_in_union_matrix: NA_not_zero",
        f"blocking_missing_input: {blocking}",
        f"blocking_reason_summary: {'; '.join(blocking_reasons) if blocking_reasons else 'none'}",
    ]
    manifest_path = run_root / "00_manifest" / "step2_run_manifest.yaml"
    write_yaml(manifest_path, manifest_lines)

    # checkpoint
    checkpoint_path = run_root / "00_manifest" / "checkpoints" / "input_file_inventory.checkpoint.csv"
    csv_write(checkpoint_path, ["run_id", "created_at", "source_path", "status", "last_file_written"], [{
        "run_id": run_id,
        "created_at": now_iso(),
        "source_path": str(inventory_path),
        "status": "completed" if not blocking else "blocking",
        "last_file_written": str(manifest_path),
    }])

    return {
        "run_id": run_id,
        "run_root": run_root,
        "blocking_missing_input": blocking,
        "blocking_reasons": blocking_reasons,
        "h5ad_count": h5ad_count,
        "sidecar_obs_count": sidecar_obs_count,
        "sidecar_var_count": sidecar_var_count,
        "manifest_path": manifest_path,
        "inventory_path": inventory_path,
        "field_check_path": field_check_path,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description="Step2.0 run manifest and inventory")
    parser.add_argument("--run-id", default=None, help="run_id (default auto)")
    parser.add_argument("--snapshot-path", default=None, help="path to snapshot.yaml")
    parser.add_argument("--root-dir", default=str(ROOT), help="project root")
    parser.add_argument("--workers", type=int, default=None, help="parallel sha256 workers (default cpu_count-1)")
    args = parser.parse_args()

    result = run_step2_manifest_inventory(
        root_dir=args.root_dir,
        run_id=args.run_id,
        snapshot_path=args.snapshot_path,
        workers=args.workers,
    )
    print(f"run_id: {result['run_id']}")
    print(f"run_root: {result['run_root']}")
    print(f"blocking: {result['blocking_missing_input']}")
    if result['blocking_reasons']:
        print(f"reasons: {', '.join(result['blocking_reasons'])}")


if __name__ == "__main__":
    main()
