#!/usr/bin/env python3
"""Run v7 Stage 6 context decomposition and bridge audit."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from v7.context_analysis import run_context_analysis


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-root", default="results/v7/context")
    parser.add_argument("--stage4-root", default="results/v7/spatial_discovery")
    parser.add_argument("--stage5-root", default="results/v7/clinical_anchor")
    args = parser.parse_args()
    result = run_context_analysis(
        output_root=Path(args.output_root),
        stage4_root=Path(args.stage4_root),
        stage5_root=Path(args.stage5_root),
    )
    print(
        json.dumps(
            {
                key: result[key]
                for key in (
                    "status",
                    "n_locked_candidates",
                    "n_context_implementation_rows",
                    "n_hcc_rows",
                    "n_clinical_annotation_rows",
                    "spatial_response_bridge_status",
                )
            },
            ensure_ascii=False,
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
