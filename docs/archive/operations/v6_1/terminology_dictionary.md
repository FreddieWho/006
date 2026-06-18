# terminology_dictionary (v6.1)

## Usage rule
- 本词典用于统一项目内方法描述、表头命名、图注与讨论语义。
- 当术语用于潜在论文正文（Methods/Figure legend）时，优先采用本词典的中英文对应定义。

## Core terms

| Term | 中文定义 | English definition | Operational boundary |
|---|---|---|---|
| **PD-1 anchor** | 指在 PD-1/PD-L1 单药或较干净 ICI 场景中学习到的基础免疫敏感性与原发耐药参照系。 | Reference immune sensitivity/resistance axis learned from PD-1/PD-L1 monotherapy or clean ICI cohorts. | 必须来自 Context A；不能由高混杂联合治疗直接定义。 |
| **PD-1+X repair logic** | 基于 PD-1 anchor 识别残余障碍，并解释 X 修复了哪一类障碍。 | A mechanism chain that maps residual barriers from the PD-1 anchor to the repair function of X. | 每条 X nomination 必须可追溯到具体 anchor barrier。 |
| **shared immune module** | 跨癌种、跨中心可复现并与 ICI 敏感性/耐药相关的免疫状态模块。 | Cross-cancer and cross-center reproducible immune module associated with ICI sensitivity/resistance. | 至少满足 L1 稳健性标准。 |
| **HCC-specific barrier module** | 在 HCC 组织生态中表现出特异性且与 PD-1 耐药相关的障碍模块。 | HCC-enriched barrier module linked to PD-1 resistance within liver-specific immune ecology. | 必须报告与 pan-cancer 共享模块的差异证据。 |
| **mIMS** | module immune mechanism score；患者在某机制模块上的激活/抑制强度读数。 | Module Immune Mechanism Score; per-patient activation/suppression intensity for a mechanism module. | 属于基础主读数，不等同于干预可行性。 |
| **MRI** | mechanism restoration index；某 X 机制将患者状态推向 responder-like 的恢复程度。 | Mechanism Restoration Index; extent to which an X mechanism shifts a patient toward a responder-like state. | 必须相对 PD-1 anchor 定义；非绝对疗效概率。 |
| **mICS** | module intervention candidate score；整合中心性、方向一致性、druggability、空间支持与外部锚点的一体化提名分数。 | Module Intervention Candidate Score; integrated nomination score combining centrality, directionality, druggability, spatial support, and external anchors. | 仅用于 nomination 排序，不直接代表临床获益。 |
| **ESR** | edge/evidence support rate；模块边级证据覆盖度与一致性比例。 | Edge/Evidence Support Rate; proportion and consistency of supported edges within a module. | 用于机制可信度补充，不单独定义候选优先级。 |
| **LRG-spatial** | ligand-receptor gain with spatial support；具备空间共定位基础的 LR 增益证据。 | Ligand-Receptor Gain with Spatial Support; LR evidence conditioned on spatial co-localization. | 无空间支撑时不得标记为 LRG-spatial positive。 |
| **biological adjudication** | 由空间/组织层证据对机制假设进行“裁判式”确认或降级。 | A referee-like biological decision layer where spatial/tissue evidence confirms or downgrades mechanism claims. | 结果可提升也可否决机制等级。 |
| **translational nomination** | 结合 L1/L2 证据、外部锚点与最小实验可行性后形成的候选机制/靶点提名。 | Candidate mechanism/target nomination supported by L1/L2 evidence, external anchors, and minimal experimental feasibility. | 属于研究提名，不等价于临床处方推荐。 |

## Synonym control
- `PD1 anchor`, `PD-1 anchor`, `monotherapy anchor` -> 统一写作 **PD-1 anchor**。
- `repair module`, `X repair`, `combination repair logic` -> 统一写作 **PD-1+X repair logic**。
- `spatial validation`（若语义为裁判）-> 建议改写为 **biological adjudication**。
- `candidate ranking`（机制提名语境）-> 建议改写为 **translational nomination ranking**。

## Not-recommended phrasing
- “clinical recommendation model” -> 在 v6.1 阶段替换为 “mechanism nomination framework”。
- “causal proof” -> 若仅有统计与先验证据，应写作 “causal-candidate evidence” 或 “mechanistically supported association”。
