#!/usr/bin/env python3
"""Compatibility entrypoint: rescore outputs and write selection manifest."""
from __future__ import annotations

import argparse

from score_phase6_algorithm_modules import run


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--mode", choices=["smoke", "full", "strengthened"], default="full")
    run(parser.parse_args().mode)
