"""Deterministic unpack stage for GEO RAW.tar cohorts.

Extracts RAW.tar into a deterministic working directory structure:
  data/imm/<cohort_id>/

Usage:
    python unpack_raw_tar.py <cohort_id> [cohort_id ...]
"""
import sys
import tarfile
from pathlib import Path

BASE = Path("/home/huyudi/006")


def unpack_cohort(cohort_id: str):
    src_dir = BASE / "data" / "imm" / cohort_id
    tar_path = src_dir / f"{cohort_id}_RAW.tar"
    if not tar_path.exists():
        print(f"[{cohort_id}] ERROR: {tar_path} not found", file=sys.stderr)
        return False

    print(f"[{cohort_id}] unpacking {tar_path} …")
    with tarfile.open(tar_path, "r") as tf:
        for member in tf.getmembers():
            # skip macOS resource forks
            if member.name.startswith("._") or "/._" in member.name:
                continue
            tf.extract(member, path=src_dir)
    print(f"[{cohort_id}] done")
    return True


if __name__ == "__main__":
    targets = sys.argv[1:]
    if not targets:
        print("Usage: python unpack_raw_tar.py <cohort_id> [cohort_id ...]", file=sys.stderr)
        sys.exit(1)
    ok = 0
    for cid in targets:
        if unpack_cohort(cid):
            ok += 1
    print(f"Unpacked {ok}/{len(targets)} cohorts.")
