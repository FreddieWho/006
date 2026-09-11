# Stage 3 阶段摘要：空间基础设施

- 科学问题：空间输入、坐标、图和患者层级能否可靠回放？（节点计划 Stage 3）
- 范围版本：plan_v7.0；8 逻辑 pilot、11 物理捕获、7 患者身份包（Visium/Xenium；HTAN CRC、GSE238264 HCC、GSE291246 BCC）。
- 现状：`SPATIAL_FOUNDATION_REPORT.md`＋D3_GATE（`D3_PASS_WITH_LIMITATIONS`）＋run manifest 已存在；本次补齐五件套，不重跑全量回放。
- 方法：同一 adapter/身份/分割/坐标合同；matrix/counts 读取、spot/cell/bin identity、坐标系与单位、tissue mask、gene mapping、region/GT 隔离、patient/block/section 嵌套；每切片独立构图（section-local kNN），禁跨 section/patient 连边；新增或改变路径才做技术回放。
- 主要结果：11/11 捕获通过 identity/coordinate/raw-count/GT isolation；1 共同特征跨 Visium/Xenium 完整可测并生成非退化统计；患者 split 同患者不跨 fold。
- 负结果/限制：Xenium image/composition 未提供、部分特征面板覆盖不足，保留限制不插值；kNN 邻接≠真实细胞接触；Meylan/USZ/ST_CRC_CMS 锁定为后续验证，本阶段不读 response。
- 判定：合格输入可重放，计数与坐标语义闭合；不在此判定屏障/疗效/潜在场成立。D3 `PASS_WITH_LIMITATIONS`，可进 Stage 4。
- 下游用途：给 Stage 4 精确到 unit 的可用输入、特征覆盖与不允许的结构解释（见 NEXT）。
