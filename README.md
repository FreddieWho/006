# 006 v7 — PD-1 空间屏障结构与 X-class 修复

本项目研究：跨癌种免疫失败程序如何形成可重复的空间屏障结构，器官环境如何改变这些结构，以及 PD1+X 是否与其修复性重排相容。

v7 保留泛癌共性、HCC context、PD1+X repair 和 X-class 外推的宏观目标，并把空间数据从末端验证提升为主线。空间主线同时承担开放潜在场发现、可解释屏障结构、跨癌/跨平台迁移、HCC 组织实现和 repair 方向约束。

## 当前入口

- [科学方案](docs/plan/v7/plan.md)
- [执行路线图](docs/plan/v7/roadmap.md)
- [当前状态](results/v7/STATUS.md)
- [2026-09-11修复审核与实际重算](docs/v7/STAGE_CLOSURE_REVIEW_2026-09-11.md)
- [v6 资产处置](docs/v7/V6_ASSET_DISPOSITION.md)
- [空间资产接入](docs/v7/SPATIAL_ASSET_INTAKE.md)
- [仓库迁移说明](docs/v7/REPOSITORY_MIGRATION.md)
- [科学一致性复核](docs/v7/V7_PLAN_SCIENTIFIC_REVIEW.md)
- [项目索引](INDEX.md)
- [目录规范](STRUCTURE.md)

## 当前事实

- `/006` 已有可追溯的多癌种单细胞、bulk、perturbation、TCR 和少量空间资产；
- v6 已建立数据接口、response-blind 分子程序与测量审计，但没有正式建立 shared PD-1 failure、dominant barrier 或 X repair；
- `~/013_spatial` 约有 1.9 TB 多平台空间资产，其中 HTAN CRC、GSE175540、ST_CRC_CMS 和 USZ TLS 具有较清晰的患者/切片/结构锚点；
- 当前 R-04 已完成多 restart 稳定性诊断，但状态仍为 `RESTART_STABILITY_DIAGNOSTIC_COMPLETE_NOT_FORMAL_K_SELECTION`；只复用接口和审计，不继承潜在场或有效 K 结论；
- v7基础层已有产物；本轮补了空间置换/低秩留出、临床共同折、配对/组成/方向和真实扰动等重算，当前状态为 `V7_REVIEW_REMEDIATION_PARTIAL_SCIENTIFIC_OPEN`。完整空间屏障、独立验证与repair证据链仍未闭合。
- 当前已审计资产没有闭合的同患者纵向 spatial PD1+X 链；直接空间 repair 暂为 `NOT_IDENTIFIABLE`，不能用纵向 scRNA 与 post-only spatial 拼接替代。

## 数据与 Git

原始数据、大型 h5ad/parquet、图像、checkpoint 和历史结果保留在本地，通过 manifest 和 SHA-256 引用，不进入普通 GitHub 提交。Git 仓库只维护代码、测试、方案、小型配置和可审计清单。

## 结论边界

项目输出是可复现的空间屏障结构、context modulation 和可实验检验的 repair hypothesis，不是临床用药推荐。没有可交换 treatment arms 时，PD1+X 结果只能称为 combination-associated 或 repair-compatible rewiring。

历史 v5/v6 文档位于 [`docs/archive/`](docs/archive/README.md)。
