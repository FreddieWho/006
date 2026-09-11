#!/usr/bin/env python3
"""Run the v7 patient-level clinical baseline."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from v7.clinical_baseline import run_clinical_baseline


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-root", default="results/v7/clinical_anchor")
    parser.add_argument("--response-permutations", type=int, default=50)
    parser.add_argument("--bootstrap-draws", type=int, default=1000)
    parser.add_argument("--seed", type=int, default=20260911)
    args = parser.parse_args()
    result = run_clinical_baseline(
        output_root=Path(args.output_root),
        response_permutations=args.response_permutations,
        bootstrap_draws=args.bootstrap_draws,
        seed=args.seed,
    )
    print(
        json.dumps(
            {
                key: result[key]
                for key in (
                    "status",
                    "datasets",
                    "n_patients",
                    "n_patients_by_dataset",
                    "feature_sets",
                    "paired_spatial_rewiring_status",
                )
            },
            ensure_ascii=False,
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
