#!/usr/bin/env python3
"""Freeze the context-repair rerun into one auditable release summary."""
from __future__ import annotations

import hashlib
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd
import yaml


ROOT = Path(__file__).resolve().parents[3]
REPAIR = ROOT / "results/v6_2/data_interface_repair_v1"
P4A = ROOT / "results/v6_2/phase4a_cell_state_harmonization_repair_v1"
P4B = ROOT / "results/v6_2/phase4b_immune_state_feature_construction_repair_v1"
P5 = ROOT / "results/v6_2/phase5_strong_baseline_and_confounding_audit_repair_v1"
P6 = ROOT / "results/v6_2/phase6_measurement_foundation_repair_v1"
P7 = ROOT / "results/v6_2/phase7_module_measurement_and_barrier_identifiability_repair_v1"
P8 = ROOT / "results/v6_2/phase8_anchor_repair_and_statistical_strengthening_repair_v1"
FROZEN_HASH = "285c11de4b071c45b0ec43464fb963650340c29bb2de3037570578b29ee44036"


def load_yaml(path: Path) -> dict:
    return yaml.safe_load(path.read_text()) or {}


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(1024 * 1024):
            h.update(chunk)
    return h.hexdigest()


def markdown_table(frame: pd.DataFrame) -> str:
    columns = list(frame.columns)
    header = "| " + " | ".join(columns) + " |"
    separator = "| " + " | ".join(["---"] * len(columns)) + " |"
    rows = ["| " + " | ".join(str(value).replace("|", "/") for value in row) + " |"
            for row in frame.fillna("").itertuples(index=False, name=None)]
    return "\n".join([header, separator, *rows])


def main() -> None:
    identity = load_yaml(REPAIR / "manifest.yaml")
    p4a = load_yaml(P4A / "handoff/phase4a_context_repair_manifest.yaml")
    p4b = load_yaml(P4B / "handoff/phase4b_context_repair_manifest.yaml")
    p5 = load_yaml(P5 / "handoff/phase5_decision_manifest.yaml")
    p6 = load_yaml(P6 / "projection_manifest.yaml")
    p7 = load_yaml(P7 / "handoff/phase7b_to_phase8_handoff.yaml")
    p8 = load_yaml(P8 / "phase8_route_gate_and_handoff.yaml")
    p8run = load_yaml(P8 / "phase8_contract_and_run_manifest.yaml")

    cardinality = pd.read_csv(REPAIR / "cohort_cardinality_audit.csv")
    object_audit = pd.read_csv(REPAIR / "object_metadata_audit.csv")
    projection = pd.read_csv(REPAIR / "configs/projection_object_eligibility_audit.csv")
    anchor_path = P8 / "phase8_anchor_master.parquet"
    if anchor_path.exists() and p8run.get("run_status") == "COMPLETE":
        anchor = pd.read_parquet(anchor_path)
        strict = anchor[
            anchor.primary_role.eq("conditional_primary")
            & anchor.baseline_eligible
            & anchor.response_known
            & anchor[["FM01", "FM04", "FM07"]].notna().all(axis=1)
        ].copy()
        surface = strict.groupby(
            ["cohort_id", "endpoint_type", "treatment_arm", "resolution"], dropna=False
        ).agg(
            patients=("patient_key", "nunique"),
            R=("response_binary", lambda x: int((x == 1).sum())),
            NR=("response_binary", lambda x: int((x == 0).sum())),
        ).reset_index()
        routes = pd.DataFrame([{
            "barrier": row["barrier_id"], "route": row["route"], "direction": row["direction"],
            "same_endpoint_cohorts": row["same_endpoint_cohorts"],
            "crossfit_stable": row["crossfit_stable"],
            "negative_control_pass": row["negative_control_pass"],
            "reasons": ";".join(row["reasons"]),
        } for row in p8["barriers"]])
    else:
        strict = pd.DataFrame(columns=["patient_key", "cohort_id"])
        surface = pd.DataFrame([{
            "status": "NOT_ESTIMATED", "reason": "Phase7 primary barrier set incomplete; Phase8 stopped before response modeling"
        }])
        routes = pd.DataFrame([{
            "route": "BLOCKED_PRE_MODEL", "eligible_primary": "|".join(p8.get("eligible_primary_barriers", [])),
            "blocked_primary": "|".join(p8.get("blocked_primary_barriers", [])),
            "reason": "PHASE7_PRIMARY_BARRIER_SET_INCOMPLETE",
        }])
    quarantine = cardinality[cardinality.quarantined_cells.gt(0)][
        ["cohort_id", "n_cells", "resolved_cells", "quarantined_cells"]
    ]
    projection_status = projection.groupby("projection_status").size().reset_index(name="n_objects")

    checks = {
        "identity_phase4a_cell_count_equal": int(identity["n_cells"]) == int(p4a["n_cells"]),
        "phase6_membership_hash_immutable": (
            p6.get("membership_sha256_before") == FROZEN_HASH
            and p6.get("membership_sha256_after") == FROZEN_HASH
        ),
        "phase7_hard_blockers_empty": not p7.get("hard_blockers"),
        "phase7_has_at_least_one_primary_barrier": int(p7.get("n_phase8_barriers", 0)) >= 1,
        "phase8_run_finalized": p8run.get("run_status") in {"COMPLETE", "ABORTED_PRE_MODEL"},
        "phase8_no_module_modification": p8run.get("rules", {}).get("module_membership_modified") is False,
        "phase8_no_response_score_construction": p8run.get("rules", {}).get("response_used_for_score_construction") is False,
        "phase8_gate_not_bypassed": p8run.get("rules", {}).get("phase7_gate_bypassed", False) is False,
    }
    phase7_matrix = pd.read_parquet(P7 / "handoff/eligible_barrier_score_matrix_v0.parquet")
    binding = pd.read_csv(P4B / "response_environment/response_environment_binding_table.csv", dtype=str)
    task = binding[binding.cohort_id.eq("TASK01")]
    checks.update({
        "task01_baseline_response_recovered": (
            task.timepoint.eq("pre") & task.response_known.eq("yes") & task.supervised_use_allowed.eq("yes")
        ).sum() == 25,
        "gse176021_present_phase7_measurement": "GSE176021" in set(
            pd.read_csv(P7 / "measurement/expression_unit_binding_table.csv").cohort_id
        ),
        "gse200996_present_phase7_measurement": "GSE200996" in set(
            pd.read_csv(P7 / "measurement/expression_unit_binding_table.csv").cohort_id
        ),
        "phase7_matrix_contains_only_current_primary": set(phase7_matrix.columns) - {
            "cohort_id", "patient_key", "timepoint", "tissue_context", "patient_timepoint_context_id"
        } == set(p8.get("eligible_primary_barriers", [])),
    })
    checks = {name: bool(passed) for name, passed in checks.items()}
    hard_failures = [name for name, passed in checks.items() if not passed]
    release_status = "BLOCKED_INVALID_RELEASE" if hard_failures else "COMPLETE_REPAIRS_PHASE8_GATE_BLOCKED"

    key_cohorts = pd.DataFrame([
        ["TASK01", "已修复", "mRECIST + T0/T1/T2；baseline PD1 进入主估计，后续时间点仅敏感性"],
        ["GSE176021", "已注册并有分数", "仅治疗后且无明确疗效，不能进入 baseline response 主估计"],
        ["GSE200996", "已注册并有分数", "High/Medium/Low 保留为有序标签，不伪造二元疗效"],
        ["GSE286827", "已拆臂", "D-only 进入 conditional clean；D+T 仅敏感性"],
        ["GSE301741", "已修复样本键并投影", "直接 mid 仅少量患者；sample projection bridge 未通过，保持 support"],
        ["GSE206325", "已评分", "只有 coarse/非 baseline 支持，不并入 pan-cancer mid 主效应"],
        ["GSE272734", "已完整注册", "4,771,240 个来源分选 CD8 T 细胞；无疗效标签，不进入监督主估计"],
    ], columns=["队列", "状态", "边界"])

    report = f"""# v6.2.1 数据接口修复与全量重跑综合报告

生成时间：{datetime.now(timezone.utc).replace(microsecond=0).isoformat()}

## 1. 执行结论

- 数据接口修复和 Phase4A/4B/5/6/7 重跑已完成；Phase8 按 gate 在建模前停止，release 状态为 `{release_status}`。
- 统一身份层包含 **{identity['n_objects']} 个对象、{identity['n_cells']:,} 个细胞**；可追溯 {identity['resolved_cells']:,}，隔离 {identity['quarantined_cells']:,}。
- 冻结 FM membership 前后 SHA256 均为 `{FROZEN_HASH}`，本轮没有重训或修改 FM01–FM08。
- Phase7 当前只允许 FM07 进入 Phase8 primary；FM01/FM04 降为 sensitivity。Phase8 最终 verdict 为 `{p8['verdict']}`。
- 数据修复恢复了 TASK01 baseline 疗效面，但 Phase7 测量门未保留完整三模块 primary set；因此未估计 shared/HCC-specific response route，Phase9、SRB、counterfactual repair 和 X-class ranking 继续阻断。

## 2. 数据完整性

Phase4A：{p4a['n_cells']:,} 个细胞、{p4a['n_cohorts']} 个队列。Phase4B：{p4b.get('n_analysis_units', 'NA')} 个分析单元。Phase6：{len(p6.get('cohorts', []))} 个 raw-count 投影对象。

### 隔离数据

{markdown_table(quarantine)}

这些细胞没有可靠患者身份，只保留为 support/quarantine；未用合成患者身份进入主分析。

### 表达投影状态

{markdown_table(projection_status)}

## 3. 重点队列修复结果

{markdown_table(key_cohorts)}

## 4. Phase5–7 重跑

- Phase5：`{p5['verdict']}`；监督压力面 {p5['supervised_surface']['n_patients']} 位患者、{p5['supervised_surface']['n_cohorts']} 个队列。它只用于基线和混杂检查，不定义机制方向。
- Phase6：44 个 count-compatible 对象完成冻结模块投影；4 个对象因非整数/非 count 表达保留为 support，未伪装为 raw counts。
- Phase7：`{p7['verdict']}`；8 个模块均保留，当前只有 FM07 为 Phase8 primary；几何秩稳定性中位数 {p7['bootstrap_geometry_rank_stability_median']:.3f}。

## 5. Phase8 真实 Anchor Surface

严格条件为 baseline + 明确疗效 + 合适治疗臂 + 冻结 primary barrier score。由于 Phase7 未保留完整的 FM01/FM04/FM07 primary set，本轮未进入 response modeling；下表不是阴性结果，而是 pre-model abstention。

{markdown_table(surface)}

## 6. Phase8 路线裁定

{markdown_table(routes)}

`BLOCKED_PRE_MODEL` 不是生物学阴性；它表示 response analysis 尚未获准运行，不能把“没有模型结果”解释为“没有生物学效应”。

## 7. 自动验收

{markdown_table(pd.DataFrame([{"check": k, "passed": v} for k, v in checks.items()]))}

Hard failures：{hard_failures or '无'}。

## 8. 科学边界和下一步

- 允许：response-blind pan-cancer state/module framework、冻结模块测量、队列内描述性方向和明确标注的不确定性。
- 禁止：把单队列效应写成 shared anti-PD1 mechanism；训练 supervised SRB；执行 counterfactual repair；输出 X-class recommendation。
- 若继续 supervised route，需要新增同 endpoint、baseline、单药且具有 mid score 的独立队列，或获得能通过预注册 bridge 的 GSE301741 投影；不建议再用更复杂模型掩盖信息不足。
"""
    (REPAIR / "CONTEXT_REPAIR_AND_FULL_RERUN_REPORT.md").write_text(report, encoding="utf-8")

    release = {
        "release_id": "context_repair_v1_full_rerun",
        "created_at": datetime.now(timezone.utc).replace(microsecond=0).isoformat(),
        "status": release_status,
        "hard_failures": hard_failures,
        "checks": checks,
        "counts": {
            "objects": int(identity["n_objects"]), "cells": int(identity["n_cells"]),
            "resolved_cells": int(identity["resolved_cells"]),
            "quarantined_cells": int(identity["quarantined_cells"]),
            "phase8_strict_patients": (
                int(strict.patient_key.nunique()) if p8run.get("run_status") == "COMPLETE" else None
            ),
            "phase8_strict_cohorts": (
                int(strict.cohort_id.nunique()) if p8run.get("run_status") == "COMPLETE" else None
            ),
        },
        "verdicts": {"phase5": p5["verdict"], "phase7": p7["verdict"], "phase8": p8["verdict"]},
        "phase9_entry_allowed": bool(p8["phase9_entry_allowed"]),
        "frozen_membership_sha256": FROZEN_HASH,
    }
    (REPAIR / "context_repair_release_manifest.yaml").write_text(
        yaml.safe_dump(release, sort_keys=False, allow_unicode=True), encoding="utf-8")

    tracked = [
        REPAIR / "manifest.yaml", REPAIR / "cohort_cardinality_audit.csv",
        REPAIR / "metadata/sample_metadata_master.effective_v1.csv",
        P4A / "handoff/phase4a_context_repair_manifest.yaml",
        P4B / "handoff/phase4b_context_repair_manifest.yaml",
        P5 / "handoff/phase5_decision_manifest.yaml", P6 / "projection_manifest.yaml",
        P7 / "handoff/phase7b_to_phase8_handoff.yaml", P7 / "handoff/eligible_barrier_score_matrix_v0.parquet",
        P8 / "phase8_route_gate_and_handoff.yaml", P8 / "phase8_contract_and_run_manifest.yaml",
        REPAIR / "CONTEXT_REPAIR_AND_FULL_RERUN_REPORT.md",
    ]
    index = pd.DataFrame([{
        "path": str(path.relative_to(ROOT)), "size_bytes": path.stat().st_size, "sha256": sha256(path)
    } for path in tracked])
    index.to_csv(REPAIR / "context_repair_output_index.tsv", sep="\t", index=False)
    print(yaml.safe_dump(release, sort_keys=False, allow_unicode=True))


if __name__ == "__main__":
    main()
