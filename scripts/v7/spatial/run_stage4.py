#!/usr/bin/env python3
"""Run the v7 response-blind spatial discovery pass."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from v7.spatial_discovery import run_spatial_discovery


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-root", default="results/v7/spatial_discovery")
    parser.add_argument("--permutations", type=int, default=19)
    parser.add_argument("--seed", type=int, default=20260911)
    parser.add_argument("--nmf-sample-limit", type=int, default=1500)
    parser.add_argument("--config", default="config/v7/stage4.yaml")
    args = parser.parse_args()
    result = run_spatial_discovery(
        output_root=Path(args.output_root),
        permutations=args.permutations,
        seed=args.seed,
        nmf_sample_limit=args.nmf_sample_limit,
        config_path=Path(args.config),
    )
    print(
        json.dumps(
            {
                key: result[key]
                for key in (
                    "status",
                    "n_logical_units",
                    "n_physical_captures",
                    "n_patients",
                    "datasets",
                    "platforms",
                    "surrogate_status",
                )
            },
            ensure_ascii=False,
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
