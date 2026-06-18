#!/usr/bin/env python3
"""CLI for immutable Phase6 frozen-module projection."""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import yaml

from frozen_projection import ProjectionError, load_config, run_projection, sha256


def main() -> int:
    parser = argparse.ArgumentParser(description="Rebuild count pseudobulks and project immutable frozen modules.")
    parser.add_argument("--config", required=True, type=Path, help="YAML with h5ad, identity sidecar, Phase4A annotations.")
    parser.add_argument("--output-root", required=True, type=Path)
    parser.add_argument("--cohort", help="Optional cohort_id filter; requires sidecar cohort_id column.")
    parser.add_argument("--resume", action="store_true", help="Reuse a completed output only if membership SHA256 still matches.")
    args = parser.parse_args()
    try:
        config, inputs = load_config(args.config.resolve())
        out = args.output_root.resolve()
        manifest_path = out / "projection_manifest.yaml"
        current_hash = sha256(inputs.membership)
        if args.resume and manifest_path.exists():
            prior = yaml.safe_load(manifest_path.read_text()) or {}
            if prior.get("membership_sha256_before") != current_hash or prior.get("membership_sha256_after") != current_hash:
                raise ProjectionError("--resume refused: frozen membership SHA256 differs from existing manifest")
            required = ["unified_expression_unit_registry.csv", "frozen_module_score_matrix.parquet",
                        "repaired_cell_state_pseudobulk.logcpm.parquet", "gene_overlap_coverage_audit.csv",
                        "counts_conservation_and_block_audit.csv"]
            if all((out / name).exists() for name in required):
                print(f"resume: verified immutable projection at {out}")
                return 0
        manifest = run_projection(config, inputs, out, args.cohort)
        print(yaml.safe_dump(manifest, sort_keys=False))
        return 0
    except ProjectionError as error:
        print(f"frozen projection failed: {error}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
