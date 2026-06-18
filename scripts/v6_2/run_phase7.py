#!/usr/bin/env python3
"""Single controller entry point for v6.2.1 Phase7."""

from __future__ import annotations

import argparse
import os
from pathlib import Path

from audit_phase7 import run as run_audit
from run_phase7a_measurement import run as run_7a
from run_phase7b_identifiability import run as run_7b


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--stage", choices=["7a", "7b", "audit", "all"], default="all")
    parser.add_argument("--mode", choices=["smoke", "full"], default="full")
    parser.add_argument("--resume", action="store_true")
    parser.add_argument("--config", help="Override the Phase7 YAML configuration path.")
    args = parser.parse_args()
    if args.config:
        os.environ["PHASE7_CONFIG_PATH"] = str(Path(args.config).resolve())
    if args.stage in {"7a", "all"}:
        run_7a(args.mode)
    if args.stage in {"7b", "all"}:
        run_7b(args.mode, args.resume)
    if args.stage in {"audit", "all"}:
        result = run_audit()
        if result["verdict"] != "PASS":
            raise SystemExit(1)


if __name__ == "__main__":
    main()
