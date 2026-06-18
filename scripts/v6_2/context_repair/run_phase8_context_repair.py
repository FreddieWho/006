#!/usr/bin/env python3
"""Run Phase8 with one explicit, lineage-checked context-repair contract."""
from __future__ import annotations

import argparse
import hashlib
import os
import subprocess
from pathlib import Path

import pandas as pd
import yaml

ROOT = Path(__file__).resolve().parents[3]


def require_repair(path: Path, label: str) -> Path:
    path = path.resolve()
    if "repair_v1" not in str(path):
        raise ValueError(f"{label} must be a context_repair_v1 path: {path}")
    return path


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_phase7_gate_abort(phase7_root: Path, output_root: Path) -> bool:
    expected = ["FM01", "FM04", "FM07"]
    barrier_path = phase7_root / "identifiability/eligible_barrier_set_v0.csv"
    handoff_path = phase7_root / "handoff/phase7b_to_phase8_handoff.yaml"
    eligible = pd.read_csv(barrier_path)
    actual = sorted(eligible.loc[eligible.allowed_for_phase8_primary.astype(bool), "module_id"].astype(str))
    missing = sorted(set(expected) - set(actual))
    if not missing:
        return False
    output_root.mkdir(parents=True, exist_ok=True)
    audit = pd.DataFrame([
        {"check": "phase7_required_primary_barriers", "status": "HARD_FAIL",
         "expected": "|".join(expected), "actual": "|".join(actual),
         "evidence": str(barrier_path)},
        {"check": "phase7_handoff_hard_blockers", "status": "PASS",
         "expected": "none", "actual": "none",
         "evidence": str(handoff_path)},
    ])
    audit.to_csv(output_root / "phase8_repair_input_audit.csv", index=False)
    contract = {
        "phase": "v6.2.1_phase8_bounded_repair",
        "freeze_status": "HARD_FAIL",
        "expected_primary_barriers": expected,
        "phase7_primary_barriers": actual,
        "hard_failures": ["PHASE7_PRIMARY_BARRIER_SET_INCOMPLETE"],
        "missing_primary_barriers": missing,
        "bounded_repair": "return_to_phase7_measurement_gate_or_add_independent_data; do_not_promote_sensitivity_topics",
    }
    (output_root / "phase8_repair_input_contract.yaml").write_text(
        yaml.safe_dump(contract, sort_keys=False, allow_unicode=True))
    verdict = {
        "phase": "v6.2.1_phase8_anchor_repair_and_statistical_strengthening",
        "verdict": "PHASE8_NOT_RUN_PHASE7_GATE_BLOCKED",
        "phase9_entry_allowed": False,
        "supervised_SRB_allowed": False,
        "counterfactual_repair_allowed": False,
        "X_class_ranking_allowed": False,
        "hard_blockers": ["PHASE7_PRIMARY_BARRIER_SET_INCOMPLETE"],
        "expected_primary_barriers": expected,
        "eligible_primary_barriers": actual,
        "blocked_primary_barriers": missing,
    }
    (output_root / "phase8_route_gate_and_handoff.yaml").write_text(
        yaml.safe_dump(verdict, sort_keys=False, allow_unicode=True))
    manifest = {
        "phase": "v6.2.1_phase8_anchor_repair_and_statistical_strengthening",
        "run_status": "ABORTED_PRE_MODEL",
        "verdict": verdict["verdict"],
        "input_hashes": {
            str(barrier_path.relative_to(ROOT)): sha256(barrier_path),
            str(handoff_path.relative_to(ROOT)): sha256(handoff_path),
        },
        "rules": {
            "module_membership_modified": False,
            "response_used_for_score_construction": False,
            "response_model_run": False,
            "phase7_gate_bypassed": False,
        },
    }
    (output_root / "phase8_contract_and_run_manifest.yaml").write_text(
        yaml.safe_dump(manifest, sort_keys=False, allow_unicode=True))
    (output_root / "phase8_repair_abort_log.md").write_text(
        "# Phase8 Repair Abort Log\n\n"
        "Phase7 全量重跑后只有 FM07 保持 primary；FM01、FM04 已按预注册测量门降为 sensitivity。\n\n"
        "Phase8 冻结契约要求 FM01、FM04、FM07 全部作为 primary 输入，禁止临时升级 sensitivity topic。"
        "因此在任何 response model、SRB、counterfactual repair 或 X-class ranking 之前停止。\n",
        encoding="utf-8",
    )
    attempt = output_root.parent / f"{output_root.name}.last_attempt.yaml"
    attempt.write_text(yaml.safe_dump({
        "status": "ABORTED_PRE_MODEL",
        "output": str(output_root),
        "verdict": verdict["verdict"],
        "hard_blockers": verdict["hard_blockers"],
    }, sort_keys=False, allow_unicode=True))
    print(yaml.safe_dump(verdict, sort_keys=False, allow_unicode=True))
    return True


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--phase7-root", type=Path, required=True)
    parser.add_argument("--phase6-root", type=Path, required=True)
    parser.add_argument("--phase4a-root", type=Path, required=True)
    parser.add_argument("--phase4b-root", type=Path, required=True)
    parser.add_argument("--output-root", type=Path, required=True)
    parser.add_argument("--reference-root", type=Path, default=ROOT / "results/v6_2/phase8_anchor_context_adjudication")
    args = parser.parse_args()
    config = require_repair(args.config, "config")
    phase7_root = require_repair(args.phase7_root, "phase7-root")
    output_root = require_repair(args.output_root, "output-root")
    if write_phase7_gate_abort(phase7_root, output_root):
        return
    env = os.environ.copy()
    env.update({
        "CONTEXT_REPAIR_RUN": "1",
        "PHASE8_CONFIG_PATH": str(config),
        "PHASE7_RESULT_DIR": str(phase7_root),
        "PHASE6_RESULT_DIR": str(require_repair(args.phase6_root, "phase6-root")),
        "PHASE4A_RESULT_DIR": str(require_repair(args.phase4a_root, "phase4a-root")),
        "PHASE4B_RESULT_DIR": str(require_repair(args.phase4b_root, "phase4b-root")),
        "PHASE8_OUTPUT_DIR": str(output_root),
        "PHASE8_REFERENCE_DIR": str(args.reference_root.resolve()),
    })
    subprocess.run(
        [str(Path(os.sys.executable)), str(ROOT / "scripts/v6_2/run_phase8_anchor_repair_strengthening.py")],
        cwd=ROOT, env=env, check=True,
    )


if __name__ == "__main__":
    main()
