from __future__ import annotations

import argparse
import csv
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, Iterable, List, Sequence


ROOT = Path(__file__).resolve().parents[2]
STEP1_DIR = ROOT / "results" / "v6_1" / "step1"
DATA_POOL_DIR = ROOT / "results" / "v6_1" / "data_pool"


def now_iso() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat()


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


def filter_rows(rows: List[Dict[str, str]], cohort_ids: set[str]) -> List[Dict[str, str]]:
    return [row for row in rows if row.get("cohort_id", "") in cohort_ids]


def select_snapshot_cohorts(cohort_rows: List[Dict[str, str]], allowed_statuses: set[str]) -> List[Dict[str, str]]:
    return [
        row
        for row in cohort_rows
        if row.get("cohort_id")
        and row.get("inclusion_status", "") in allowed_statuses
        and str(row.get("sample_ingested", "")).lower() == "true"
    ]


def build_project_data_pool(
    root_dir: Path | str = ROOT,
    snapshot_id: str = "current_step2_input",
    allowed_statuses: Sequence[str] = ("main", "sensitivity"),
) -> Dict[str, object]:
    root = Path(root_dir)
    step1_dir = root / "results" / "v6_1" / "step1"
    data_pool_dir = root / "results" / "v6_1" / "data_pool"
    registry_dir = data_pool_dir / "registry"
    packages_dir = data_pool_dir / "packages"
    snapshot_dir = data_pool_dir / "snapshots" / snapshot_id
    step1_exports_dir = snapshot_dir / "step1_exports"

    cohort_freeze_rows = read_csv_dicts(step1_dir / "frozen_cohort_inclusion_v6_1.csv")
    patient_split_rows = read_csv_dicts(step1_dir / "frozen_patient_split_v6_1.csv")
    sample_broad_rows = read_csv_dicts(step1_dir / "sample_metadata_master_v6_1.broad_response.csv")
    patient_broad_rows = read_csv_dicts(step1_dir / "patient_metadata_master_v6_1.broad_response.csv")
    source_rows = read_csv_dicts(step1_dir / "source_records_local.tsv", delimiter=None)

    snapshot_cohort_rows = select_snapshot_cohorts(cohort_freeze_rows, set(allowed_statuses))
    snapshot_cohort_ids = {row["cohort_id"] for row in snapshot_cohort_rows}

    sample_broad_by_cohort: Dict[str, List[Dict[str, str]]] = {}
    patient_broad_by_cohort: Dict[str, List[Dict[str, str]]] = {}
    split_by_cohort: Dict[str, List[Dict[str, str]]] = {}
    for cohort_id in snapshot_cohort_ids:
        sample_broad_by_cohort[cohort_id] = filter_rows(sample_broad_rows, {cohort_id})
        patient_broad_by_cohort[cohort_id] = filter_rows(patient_broad_rows, {cohort_id})
        split_by_cohort[cohort_id] = filter_rows(patient_split_rows, {cohort_id})

    package_rows: List[Dict[str, object]] = []
    asset_registry_rows: List[Dict[str, object]] = []
    h5ad_inventory_rows: List[Dict[str, object]] = []

    for cohort_id in sorted(snapshot_cohort_ids):
        cohort_package_dir = packages_dir / cohort_id
        cohort_package_dir.mkdir(parents=True, exist_ok=True)

        sample_rows = sample_broad_by_cohort[cohort_id]
        patient_rows = patient_broad_by_cohort[cohort_id]
        split_rows = split_by_cohort[cohort_id]
        local_rows = [row for row in source_rows if row.get("cohort_id") == cohort_id]

        sample_path = cohort_package_dir / "sample_metadata_master_v6_1.broad_response.csv"
        patient_path = cohort_package_dir / "patient_metadata_master_v6_1.broad_response.csv"
        split_path = cohort_package_dir / "frozen_patient_split_v6_1.csv"

        if sample_rows:
            csv_write(sample_path, list(sample_rows[0].keys()), sample_rows)
        if patient_rows:
            csv_write(patient_path, list(patient_rows[0].keys()), patient_rows)
        if split_rows:
            csv_write(split_path, list(split_rows[0].keys()), split_rows)

        authoritative_h5ad = ""
        for row in local_rows:
            source_path = row.get("source_path", "")
            asset_registry_rows.append(
                {
                    "cohort_id": cohort_id,
                    "source_type": row.get("source_type", ""),
                    "source_path": source_path,
                    "status": row.get("status", ""),
                    "purpose": row.get("purpose", ""),
                    "file_exists": Path(source_path).exists() if source_path else False,
                }
            )
            if row.get("source_type") == "local_h5ad" and row.get("status") == "used" and source_path:
                authoritative_h5ad = source_path

        if authoritative_h5ad:
            h5ad_inventory_rows.append({"cohort_id": cohort_id, "file_path": authoritative_h5ad})

        package_manifest = cohort_package_dir / "package_manifest.yaml"
        write_yaml(
            package_manifest,
            [
                'manifest_version: "v6.1"',
                f'cohort_id: "{cohort_id}"',
                f'created_at: "{now_iso()}"',
                f'authoritative_h5ad: "{authoritative_h5ad}"',
                f'sample_metadata_ref: "{sample_path.as_posix()}"',
                f'patient_metadata_ref: "{patient_path.as_posix()}"',
                f'patient_split_ref: "{split_path.as_posix()}"',
            ],
        )

        package_rows.append(
            {
                "cohort_id": cohort_id,
                "authoritative_h5ad": authoritative_h5ad,
                "sample_metadata_ref": sample_path.as_posix(),
                "patient_metadata_ref": patient_path.as_posix(),
                "patient_split_ref": split_path.as_posix(),
                "package_manifest_ref": package_manifest.as_posix(),
            }
        )

    cohort_freeze_snapshot = step1_exports_dir / "frozen_cohort_inclusion_v6_1.csv"
    patient_split_snapshot = step1_exports_dir / "frozen_patient_split_v6_1.csv"
    sample_broad_snapshot = step1_exports_dir / "sample_metadata_master_v6_1.broad_response.csv"
    patient_broad_snapshot = step1_exports_dir / "patient_metadata_master_v6_1.broad_response.csv"

    if snapshot_cohort_rows:
        csv_write(cohort_freeze_snapshot, list(snapshot_cohort_rows[0].keys()), snapshot_cohort_rows)
    if filter_rows(patient_split_rows, snapshot_cohort_ids):
        split_rows = filter_rows(patient_split_rows, snapshot_cohort_ids)
        csv_write(patient_split_snapshot, list(split_rows[0].keys()), split_rows)
    if filter_rows(sample_broad_rows, snapshot_cohort_ids):
        rows = filter_rows(sample_broad_rows, snapshot_cohort_ids)
        csv_write(sample_broad_snapshot, list(rows[0].keys()), rows)
    if filter_rows(patient_broad_rows, snapshot_cohort_ids):
        rows = filter_rows(patient_broad_rows, snapshot_cohort_ids)
        csv_write(patient_broad_snapshot, list(rows[0].keys()), rows)

    h5ad_inventory_path = snapshot_dir / "h5ad_inventory.csv"
    if h5ad_inventory_rows:
        csv_write(h5ad_inventory_path, list(h5ad_inventory_rows[0].keys()), h5ad_inventory_rows)

    cohort_registry_path = registry_dir / "cohort_registry.csv"
    asset_registry_path = registry_dir / "asset_registry.csv"
    package_registry_path = registry_dir / "package_registry.csv"
    if cohort_freeze_rows:
        csv_write(cohort_registry_path, list(cohort_freeze_rows[0].keys()), cohort_freeze_rows)
    if asset_registry_rows:
        csv_write(asset_registry_path, list(asset_registry_rows[0].keys()), asset_registry_rows)
    if package_rows:
        csv_write(package_registry_path, list(package_rows[0].keys()), package_rows)

    snapshot_manifest = snapshot_dir / "snapshot.yaml"
    write_yaml(
        snapshot_manifest,
        [
            'manifest_version: "v6.1"',
            f'snapshot_id: "{snapshot_id}"',
            f'created_at: "{now_iso()}"',
            f'cohort_freeze_ref: "{cohort_freeze_snapshot.as_posix()}"',
            f'patient_split_ref: "{patient_split_snapshot.as_posix()}"',
            f'sample_metadata_ref: "{sample_broad_snapshot.as_posix()}"',
            f'patient_metadata_ref: "{patient_broad_snapshot.as_posix()}"',
            f'h5ad_inventory_ref: "{h5ad_inventory_path.as_posix()}"',
        ],
    )

    return {
        "data_pool_dir": data_pool_dir,
        "snapshot_dir": snapshot_dir,
        "snapshot_manifest": snapshot_manifest,
        "cohort_count": len(snapshot_cohort_ids),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description="Build a project-local data pool snapshot from current Step1 frozen outputs.")
    parser.add_argument("--root-dir", default=str(ROOT))
    parser.add_argument("--snapshot-id", default="current_step2_input")
    args = parser.parse_args()
    result = build_project_data_pool(root_dir=args.root_dir, snapshot_id=args.snapshot_id)
    print(f"snapshot_dir={result['snapshot_dir']}")
    print(f"cohort_count={result['cohort_count']}")


if __name__ == "__main__":
    main()
