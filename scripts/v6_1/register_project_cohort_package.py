from __future__ import annotations

import argparse
import csv
import shutil
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, Iterable, List, Sequence


ROOT = Path(__file__).resolve().parents[2]


def now_iso() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat()


def read_csv_dicts(path: Path) -> List[Dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def csv_write(path: Path, fieldnames: Sequence[str], rows: Iterable[Dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def upsert_rows(existing: List[Dict[str, str]], incoming: List[Dict[str, str]], key_fields: Sequence[str]) -> List[Dict[str, str]]:
    by_key = {tuple(row.get(k, "") for k in key_fields): row for row in existing}
    for row in incoming:
        by_key[tuple(row.get(k, "") for k in key_fields)] = row
    return [by_key[key] for key in sorted(by_key)]


def copy_csv(src: Path, dst: Path) -> List[Dict[str, str]]:
    rows = read_csv_dicts(src)
    if rows:
        csv_write(dst, list(rows[0].keys()), rows)
    else:
        dst.write_text("", encoding="utf-8")
    return rows


def write_yaml(path: Path, lines: List[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def register_project_cohort_package(
    root_dir: Path | str,
    cohort_id: str,
    h5ad_path: Path,
    sample_metadata_csv: Path,
    patient_metadata_csv: Path,
    patient_split_csv: Path,
    cohort_freeze_csv: Path,
    snapshot_id: str = "current_step2_input",
) -> Dict[str, object]:
    root = Path(root_dir)
    data_pool_dir = root / "results" / "v6_1" / "data_pool"
    registry_dir = data_pool_dir / "registry"
    packages_dir = data_pool_dir / "packages"
    snapshot_dir = data_pool_dir / "snapshots" / snapshot_id
    step1_exports = snapshot_dir / "step1_exports"

    package_dir = packages_dir / cohort_id
    package_dir.mkdir(parents=True, exist_ok=True)

    dst_sample = package_dir / "sample_metadata_master_v6_1.broad_response.csv"
    dst_patient = package_dir / "patient_metadata_master_v6_1.broad_response.csv"
    dst_split = package_dir / "frozen_patient_split_v6_1.csv"
    dst_freeze = package_dir / "frozen_cohort_inclusion_v6_1.csv"

    sample_rows = copy_csv(sample_metadata_csv, dst_sample)
    patient_rows = copy_csv(patient_metadata_csv, dst_patient)
    split_rows = copy_csv(patient_split_csv, dst_split)
    freeze_rows = copy_csv(cohort_freeze_csv, dst_freeze)

    package_manifest = package_dir / "package_manifest.yaml"
    write_yaml(
        package_manifest,
        [
            'manifest_version: "v6.1"',
            f'cohort_id: "{cohort_id}"',
            f'created_at: "{now_iso()}"',
            f'authoritative_h5ad: "{h5ad_path.as_posix()}"',
            f'sample_metadata_ref: "{dst_sample.as_posix()}"',
            f'patient_metadata_ref: "{dst_patient.as_posix()}"',
            f'patient_split_ref: "{dst_split.as_posix()}"',
            f'cohort_freeze_ref: "{dst_freeze.as_posix()}"',
        ],
    )

    asset_registry_path = registry_dir / "asset_registry.csv"
    package_registry_path = registry_dir / "package_registry.csv"
    cohort_registry_path = registry_dir / "cohort_registry.csv"

    asset_rows = read_csv_dicts(asset_registry_path) if asset_registry_path.exists() else []
    asset_rows = upsert_rows(
        asset_rows,
        [
            {
                "cohort_id": cohort_id,
                "source_type": "local_h5ad",
                "source_path": h5ad_path.as_posix(),
                "status": "used",
                "purpose": "manual_project_data_pool_registration",
                "file_exists": str(h5ad_path.exists()),
            }
        ],
        ["cohort_id", "source_type", "source_path"],
    )
    csv_write(asset_registry_path, list(asset_rows[0].keys()), asset_rows)

    package_rows = read_csv_dicts(package_registry_path) if package_registry_path.exists() else []
    package_rows = upsert_rows(
        package_rows,
        [
            {
                "cohort_id": cohort_id,
                "authoritative_h5ad": h5ad_path.as_posix(),
                "sample_metadata_ref": dst_sample.as_posix(),
                "patient_metadata_ref": dst_patient.as_posix(),
                "patient_split_ref": dst_split.as_posix(),
                "package_manifest_ref": package_manifest.as_posix(),
            }
        ],
        ["cohort_id"],
    )
    csv_write(package_registry_path, list(package_rows[0].keys()), package_rows)

    cohort_registry_rows = read_csv_dicts(cohort_registry_path) if cohort_registry_path.exists() else []
    cohort_registry_rows = upsert_rows(cohort_registry_rows, freeze_rows, ["cohort_id"])
    csv_write(cohort_registry_path, list(cohort_registry_rows[0].keys()), cohort_registry_rows)

    snapshot_freeze_path = step1_exports / "frozen_cohort_inclusion_v6_1.csv"
    snapshot_split_path = step1_exports / "frozen_patient_split_v6_1.csv"
    snapshot_sample_path = step1_exports / "sample_metadata_master_v6_1.broad_response.csv"
    snapshot_patient_path = step1_exports / "patient_metadata_master_v6_1.broad_response.csv"
    snapshot_h5ad_inventory = snapshot_dir / "h5ad_inventory.csv"

    snapshot_freeze_rows = read_csv_dicts(snapshot_freeze_path) if snapshot_freeze_path.exists() else []
    snapshot_split_rows = read_csv_dicts(snapshot_split_path) if snapshot_split_path.exists() else []
    snapshot_sample_rows = read_csv_dicts(snapshot_sample_path) if snapshot_sample_path.exists() else []
    snapshot_patient_rows = read_csv_dicts(snapshot_patient_path) if snapshot_patient_path.exists() else []
    snapshot_h5ad_rows = read_csv_dicts(snapshot_h5ad_inventory) if snapshot_h5ad_inventory.exists() else []

    snapshot_freeze_rows = upsert_rows(snapshot_freeze_rows, freeze_rows, ["cohort_id"])
    snapshot_split_rows = upsert_rows(snapshot_split_rows, split_rows, ["cohort_id", "patient_id"])
    snapshot_sample_rows = upsert_rows(snapshot_sample_rows, sample_rows, ["cohort_id", "patient_id", "sample_id"])
    snapshot_patient_rows = upsert_rows(snapshot_patient_rows, patient_rows, ["cohort_id", "patient_id"])
    snapshot_h5ad_rows = upsert_rows(snapshot_h5ad_rows, [{"cohort_id": cohort_id, "file_path": h5ad_path.as_posix()}], ["cohort_id"])

    csv_write(snapshot_freeze_path, list(snapshot_freeze_rows[0].keys()), snapshot_freeze_rows)
    csv_write(snapshot_split_path, list(snapshot_split_rows[0].keys()), snapshot_split_rows)
    csv_write(snapshot_sample_path, list(snapshot_sample_rows[0].keys()), snapshot_sample_rows)
    csv_write(snapshot_patient_path, list(snapshot_patient_rows[0].keys()), snapshot_patient_rows)
    csv_write(snapshot_h5ad_inventory, list(snapshot_h5ad_rows[0].keys()), snapshot_h5ad_rows)

    snapshot_manifest = snapshot_dir / "snapshot.yaml"
    write_yaml(
        snapshot_manifest,
        [
            'manifest_version: "v6.1"',
            f'snapshot_id: "{snapshot_id}"',
            f'created_at: "{now_iso()}"',
            f'cohort_freeze_ref: "{snapshot_freeze_path.as_posix()}"',
            f'patient_split_ref: "{snapshot_split_path.as_posix()}"',
            f'sample_metadata_ref: "{snapshot_sample_path.as_posix()}"',
            f'patient_metadata_ref: "{snapshot_patient_path.as_posix()}"',
            f'h5ad_inventory_ref: "{snapshot_h5ad_inventory.as_posix()}"',
        ],
    )

    return {
        "cohort_id": cohort_id,
        "package_dir": package_dir,
        "snapshot_dir": snapshot_dir,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description="Register one cohort package into the project-local data pool without rerunning global Step1.")
    parser.add_argument("--root-dir", default=str(ROOT))
    parser.add_argument("--cohort-id", required=True)
    parser.add_argument("--h5ad-path", required=True)
    parser.add_argument("--sample-metadata-csv", required=True)
    parser.add_argument("--patient-metadata-csv", required=True)
    parser.add_argument("--patient-split-csv", required=True)
    parser.add_argument("--cohort-freeze-csv", required=True)
    parser.add_argument("--snapshot-id", default="current_step2_input")
    args = parser.parse_args()

    result = register_project_cohort_package(
        root_dir=args.root_dir,
        cohort_id=args.cohort_id,
        h5ad_path=Path(args.h5ad_path),
        sample_metadata_csv=Path(args.sample_metadata_csv),
        patient_metadata_csv=Path(args.patient_metadata_csv),
        patient_split_csv=Path(args.patient_split_csv),
        cohort_freeze_csv=Path(args.cohort_freeze_csv),
        snapshot_id=args.snapshot_id,
    )
    print(f"cohort_id={result['cohort_id']}")
    print(f"package_dir={result['package_dir']}")
    print(f"snapshot_dir={result['snapshot_dir']}")


if __name__ == "__main__":
    main()
