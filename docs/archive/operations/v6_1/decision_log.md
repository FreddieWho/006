# decision_log (v6.1)

## Log rule
- 状态字段：`locked` / `active` / `deprecated`
- 仅允许在 scope review meeting 中修改 `locked` 决策，并必须记录触发证据。

## Decisions
| ID | Date | Status | Decision | Rationale | Execution impact | Revisit trigger |
|---|---|---|---|---|---|---|
| D-001 | 2026-04-27 | locked | 采用 **PD-1-centered continuum**，不拆成 PD-1 与 PD-1+X 两个独立课题 | 生物学上 X 的价值必须定义为对 PD-1 残余障碍的修复 | 所有主结论必须包含 `PD-1 anchor -> X repair` 链条 | 若 PD-1 anchor 长期不稳定 |
| D-002 | 2026-04-27 | locked | Context C（高混杂联合治疗）只做 sensitivity/external support | 标签混杂过高，进入主训练将污染因果叙事 | Context C 结果默认进附录或敏感性章节 | 若出现高质量去混杂标注可重审 |
| D-003 | 2026-04-27 | locked | 证据等级仅保留 L1/L2/L3 三档 | 简化沟通并提高可执行性 | 证据卡强制输出三档字段 | 若审稿要求更细等级 |
| D-004 | 2026-04-27 | locked | 强基线前置为硬门槛 | 防止复杂模型“无对手”叙事 | Step2 不完成，禁止发布高级模型主结果 | 若任务转为纯描述性研究 |
| D-005 | 2026-04-27 | locked | 空间/组织数据定位为 biological adjudicator | TME 机制必须在组织生态层可解释 | 无空间裁判，不得升格 L2/L3 主结论 | 若缺失可用空间数据，仅可降级到 L1 |
| D-006 | 2026-04-27 | locked | MVP 采取审慎继承：底座可继承，机制演示层降级 | MVP 完成工程原型但机制闭环不足 | MVP 文件默认 proxy/sanity，非主证据 | 若某 MVP 输出被独立复现并通过新链路验证 |
| D-007 | 2026-04-27 | locked | 禁止 LINCS reversal 单独提名药物 | 单源方向证据不足以支持转化提名 | LINCS 只能作为 L2 辅助证据 | 若出现独立 patient+spatial 同向支持 |
| D-008 | 2026-04-27 | locked | 禁止 marker 均值/简单组合直接输出机制结论 | 过度简化，易产生伪机制 | marker 仅可用于 QC 或补充可解释性 | 若与模块/空间/扰动证据完全一致且可复现 |
| D-009 | 2026-04-27 | locked | 输出定位为 mechanism nomination，不做临床推荐系统叙事 | 当前证据链不支持患者级处方推荐 | 主文案使用 nomination 语义 | 若进入前瞻性临床验证阶段 |
| D-010 | 2026-04-27 | locked | mIMS/MRI/mICS 为主读数；ESR/LRG-spatial 为辅读数 | 与 v6.1 简化方案一致，便于工程落地 | `M_signals_v6_1.csv` 强制包含 5 项读数 | 若 pre/post 动态证据充分，可增开 FG/RT 等 |
| D-011 | 2026-04-27 | active | 英文用于术语与工程键；中文为文档主语言 | 兼顾科研表达准确性与团队协作效率 | 合同与范围文档以中文主导，术语双语 | 若投稿阶段切换为英文主稿 |
| D-012 | 2026-04-27 | active | 统一输出清单与 manifest 字段，所有主结论可追溯 | 降低后期审计与复现实验成本 | 每步必须填 `output_manifest` | 若 pipeline 重构导致字段不足 |
| D-013 | 2026-04-27 | locked | 任何主结论不得绕过中间表链路 | 防止“结果正确但不可复现” | analysis contract 中定义 no-bypass list | 若 contract 本身被正式修订 |
| D-014 | 2026-04-27 | locked | 关闭逻辑优先于叙事完整性 | 科学有效性优先于故事完整性 | 触发关闭条件时自动降级或回滚 | 若新增独立证据解除关闭条件 |

## Open items
| ID | Priority | Item | Required evidence to close |
|---|---|---|---|
| O-001 | high | PD-1 anchor 的最小 cohort 清单与样本量阈值尚未最终锁定 | Step0 cohort audit + label completeness report |
| O-002 | high | Spatial adjudication 的最小可接受数据密度阈值尚未数值化 | 空间数据 QC 基线与失败样本审计 |
| O-003 | medium | mICS 权重策略（expert prior vs data-driven）未最终定版 | 消融实验 + 外部锚点一致性比较 |

## Change protocol
1. 修改任一 `locked` 决策前，必须提交对照证据（旧链路 vs 新链路）。
2. 修改后必须同步更新：`scope_lock_v6_1.md`、`analysis_contract.md`、`terminology_dictionary.md`。
3. 所有变更必须记录日期、责任人、影响步骤、回滚条件。
