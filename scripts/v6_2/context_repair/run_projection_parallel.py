#!/usr/bin/env python3
"""Run immutable projection by cohort shards, then merge audited outputs."""
from __future__ import annotations

import argparse
import concurrent.futures
import hashlib
import subprocess
import sys
from pathlib import Path

import pandas as pd
import yaml

from frozen_projection import load_config, sha256


def safe(value: str) -> str:
    return hashlib.sha256(value.encode()).hexdigest()[:12]


def run_one(cli: Path, config: Path, shard_root: Path, cohort: str, resume: bool) -> tuple[str, int, str]:
    out = shard_root / f"{safe(cohort)}__{cohort.replace('/', '_')}"
    command = [sys.executable, str(cli), "--config", str(config), "--output-root", str(out), "--cohort", cohort]
    if resume:
        command.append("--resume")
    result = subprocess.run(command, text=True, capture_output=True)
    return cohort, result.returncode, result.stdout + result.stderr


def merge_shards(config: Path, output: Path, cohorts: list[str]) -> dict:
    _, inputs = load_config(config)
    shard_root = output / "shards"
    dirs = [shard_root / f"{safe(c)}__{c.replace('/', '_')}" for c in cohorts]
    manifests = [yaml.safe_load((path / "projection_manifest.yaml").read_text()) or {} for path in dirs]
    if any(item.get("verdict") != "PASS" for item in manifests):
        raise RuntimeError("At least one projection shard is not PASS")
    expected_hash = sha256(inputs.membership)
    if any(item.get("membership_sha256_after") != expected_hash for item in manifests):
        raise RuntimeError("Shard membership hash mismatch")
    registry = pd.concat([pd.read_csv(path / "unified_expression_unit_registry.csv") for path in dirs], ignore_index=True)
    score = pd.concat([pd.read_parquet(path / "frozen_module_score_matrix.parquet") for path in dirs], ignore_index=True)
    expression = pd.concat([pd.read_parquet(path / "repaired_cell_state_pseudobulk.logcpm.parquet") for path in dirs], ignore_index=True, sort=False)
    coverage = pd.concat([pd.read_csv(path / "gene_overlap_coverage_audit.csv") for path in dirs], ignore_index=True)
    audit = pd.concat([pd.read_csv(path / "counts_conservation_and_block_audit.csv") for path in dirs], ignore_index=True)
    if audit.status.isin({"BLOCK", "HARD_FAIL"}).any():
        raise RuntimeError("Merged projection audit contains blocking rows")
    if registry.expression_unit_id.duplicated().any() or score.expression_unit_id.duplicated().any():
        raise RuntimeError("Duplicate expression_unit_id across projection shards")
    registry.to_csv(output / "unified_expression_unit_registry.csv", index=False)
    score.to_parquet(output / "frozen_module_score_matrix.parquet", index=False)
    score.to_csv(output / "frozen_module_score_matrix.csv", index=False)
    expression.to_parquet(output / "repaired_cell_state_pseudobulk.logcpm.parquet", index=False)
    coverage.to_csv(output / "gene_overlap_coverage_audit.csv", index=False)
    audit.to_csv(output / "counts_conservation_and_block_audit.csv", index=False)
    manifest = {
        "verdict": "PASS", "execution": "parallel_cohort_shards",
        "membership_sha256_before": expected_hash, "membership_sha256_after": expected_hash,
        "membership_path": str(inputs.membership), "cohorts": cohorts,
        "n_shards": len(dirs), "n_expression_units": len(registry), "n_score_rows": len(score),
        "phase7_logcpm_parquet": str(output / "repaired_cell_state_pseudobulk.logcpm.parquet"),
        "prohibitions": ["module_fit", "module_retraining", "gene_selection", "membership_mutation"],
    }
    (output / "projection_manifest.yaml").write_text(yaml.safe_dump(manifest, sort_keys=False))
    return manifest


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--output-root", required=True, type=Path)
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--resume", action="store_true")
    args = parser.parse_args()
    config = args.config.resolve()
    output = args.output_root.resolve()
    output.mkdir(parents=True, exist_ok=True)
    (output / "shards").mkdir(exist_ok=True)
    _, inputs = load_config(config)
    cohorts = sorted({item.cohort_id for item in inputs.objects})
    cli = Path(__file__).with_name("run_frozen_projection.py")
    failures = []
    with concurrent.futures.ThreadPoolExecutor(max_workers=max(1, args.workers)) as pool:
        futures = [pool.submit(run_one, cli, config, output / "shards", cohort, args.resume) for cohort in cohorts]
        for future in concurrent.futures.as_completed(futures):
            cohort, code, log = future.result()
            (output / "shards" / f"{safe(cohort)}__{cohort.replace('/', '_')}.log").write_text(log)
            if code:
                failures.append(cohort)
    if failures:
        (output / "projection_manifest.yaml").write_text(yaml.safe_dump({
            "verdict": "BLOCKED", "failed_cohorts": failures, "formal_outputs_published": False,
        }, sort_keys=False))
        raise SystemExit(f"Projection shards failed: {failures}")
    print(yaml.safe_dump(merge_shards(config, output, cohorts), sort_keys=False))


if __name__ == "__main__":
    main()
