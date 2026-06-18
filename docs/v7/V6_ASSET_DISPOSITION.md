# v6 → v7 资产处置

日期：2026-08-31

## 总原则

v7 继承经过验证的数据和测量能力，不继承尚未成立的科学结论。历史资产不物理删除；`reuse` 表示可以作为 v7 输入，`reinterpret` 表示只能在新语义下使用，`archive` 表示仅用于溯源或风险审计，`reject` 表示禁止进入模型或结论。

## 科学资产

| 资产 | v7 状态 | 用途或原因 |
|---|---|---|
| `results/v6_2/data_interface_repair_v1/` | reuse | 当前最可靠的 patient/sample/timepoint/treatment/expression 接口 |
| `phase6_measurement_foundation_repair_v1/` | reuse | count-compatible 对象上的统一投影和表达单位 |
| `phase6_module_algorithm_benchmark_strengthened/` | reuse + reinterpret | 模块算法比较与稳定性可复用；成员和命名允许 v7 重审 |
| `phase7_*_repair_v1/` | reuse + reinterpret | 测量、missingness、route disagreement、共线和 identifiability 数值可复用 |
| v6 的 8 个 FM | reinterpret | 经验分子程序基线，不是 PD-1 failure barrier |
| FM07 | reinterpret | 旧条件下的独立归因候选，不是唯一有信息的程序 |
| `phase5_*confounding*` | archive/risk evidence | cohort-only 信号强，作为 v7 混杂基线和反例 |
| `phase8_anchor_context_adjudication/` | archive | 旧路线结果，不能作为 v7 response 证据 |
| `phase8_bounded_repair_and_readjudication/` | archive | 只保留 hypothesis seeds 与风险记录 |
| `phase8_anchor_repair_and_statistical_strengthening*_v1/` | archive | 正式状态为 pre-model blocked/not run，不是阴性生物学结果 |
| `scientific_engineering_realignment_v1/` | reuse decisions + archive gates | 架构漂移审计可复用；预注册、unlock、固定 gate 不进入 v7 |
| TASK01/TASK02 patient-timepoint 资产 | reuse after v7 audit | PD-1 reference 与 PD1+X longitudinal；保持治疗与队列语义 |
| GSE301741 repair 包 | reuse | 身份链和 coverage 根因；未修复测量不能升级为正式 anchor |
| v6.1 hash-derived perturbation scores | reject | 数值不是实验效应；只可保留人工 target/pathway annotation |
| 真实 perturbation/L1000/XAtlas | reuse after audit | 重算 intervention→module/cell-state direction |

## 代码资产

| 路径 | v7 状态 | 说明 |
|---|---|---|
| `scripts/v6_2/context_repair/` | reuse | 身份、metadata、projection 和接口修复工具 |
| `scripts/v6_2/run_phase4b.py` | reuse components | cell-state/pseudobulk/feature 工厂可拆分复用 |
| `scripts/v6_2/run_phase6.py` | reuse components | response-blind foundation 与 balanced fitting 思路 |
| `scripts/v6_2/run_phase7*.py` | reuse components | measurement/identifiability 审计，不继承 hard filter 语义 |
| `scripts/v6_2/run_phase8*.py` | archive | 旧 anchor/contract 路线，不作为 v7 入口 |
| `scripts/v6_2/*patch*`, `*rescue*`, `*official*`, `*tuned*` | archive candidate | 保留故障和算法比较溯源，不作为主入口 |
| `scripts/archive/v6_2_diagnostics/` | archive | 根目录临时检查和调试脚本 |
| `tests/v6_2/` | archive + regression | 保留 v6 接口回归；v7 新增独立测试 |

## 不继承为事实的陈述

- shared anti-PD1 failure direction 已建立；
- FM07 是唯一 barrier；
- dominant barrier 已确定；
- SRB 已实现或已被完整否定；
- X repair 已证明；
- spatial niche 已正式验证；
- zero-shot X-class 已可训练。

这些问题在 v7 中重新建立证据。
