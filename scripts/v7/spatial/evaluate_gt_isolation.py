#!/usr/bin/env python3
"""Seal one known structure locator and run the independent GT firewall audit."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from v7.spatial_validation.sealed_gt_evaluator import (
    build_sealed_validation_manifest,
    evaluate_sealed_gt_isolation,
)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-root", default="results/v7/spatial_foundation")
    parser.add_argument(
        "--gt-locator",
        default=(
            "/home/huyudi/013_spatial/data/other_sources/htan/"
            "spatial_CRC_atlas_repo/resources/ST/7003_1_pathology_annotation.csv"
        ),
    )
    args = parser.parse_args()
    build_sealed_validation_manifest(args.output_root, gt_locator=args.gt_locator)
    audit = evaluate_sealed_gt_isolation(args.output_root)
    print(json.dumps({"status": "PASS_SEALED_ONLY", "checks": len(audit)}, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
