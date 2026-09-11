#!/usr/bin/env python3
"""Run v7 Stage 7 molecular repair baseline."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from v7.repair_analysis import run_repair_analysis


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-root", default="results/v7/repair")
    parser.add_argument("--response-permutations", type=int, default=50)
    parser.add_argument("--config", default="config/v7/stage7.yaml")
    parser.add_argument("--seed", type=int, default=20260911)
    args = parser.parse_args()
    result = run_repair_analysis(
        output_root=Path(args.output_root),
        config_path=Path(args.config),
        permutations=args.response_permutations,
        seed=args.seed,
    )
    print(
        json.dumps(
            {
                key: result[key]
                for key in (
                    "status",
                    "molecular_branch",
                    "spatial_branch",
                    "n_paired_rows",
                    "n_paired_patients",
                    "n_paired_patients_by_environment",
                )
            },
            ensure_ascii=False,
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
