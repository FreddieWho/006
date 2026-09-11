#!/usr/bin/env python3
"""Run the v7 Stage 3 response-blind spatial foundation replay."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from v7.spatial_pipeline import run_spatial_foundation


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-root", default="results/v7/spatial_foundation")
    parser.add_argument(
        "--permutations",
        type=int,
        default=99,
        help="section-preserving Moran null permutations per feature (default: 99)",
    )
    parser.add_argument("--bootstrap-draws", type=int, default=1000)
    parser.add_argument("--seed", type=int, default=20260904)
    args = parser.parse_args()
    result = run_spatial_foundation(
        output_root=Path(args.output_root),
        permutations=args.permutations,
        bootstrap_draws=args.bootstrap_draws,
        seed=args.seed,
    )
    summary = {
        key: result[key]
        for key in (
            "status",
            "n_logical_pilot_rows",
            "n_physical_captures",
            "n_patients",
            "platforms",
            "datasets",
            "discovery_output_hash",
        )
    }
    print(json.dumps(summary, ensure_ascii=False, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
