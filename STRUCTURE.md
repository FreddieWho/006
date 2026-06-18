# 006 v7 目录规范

版本：3.0
日期：2026-08-31

## 顶层结构

```text
006/
├── README.md
├── INDEX.md
├── STRUCTURE.md
├── AGENTS.md
├── config/
│   └── v7/                  # 小型配置与数据角色
├── data/                    # 本地大型数据；Git 忽略
├── docs/
│   ├── plan/v7/            # 唯一活跃科学方案与路线图
│   ├── v7/                 # 当前资产、状态和迁移说明
│   └── archive/            # v5/v6 历史文档
├── scripts/
│   ├── v7/                 # v7 命令入口
│   ├── v6_1/               # 历史可重放代码
│   ├── v6_2/               # 历史可重放代码
│   ├── preprocessing/      # 数据源预处理
│   └── archive/            # 临时诊断/废弃入口
├── src/
│   └── v7/                 # v7 可复用实现
├── tests/
│   ├── v7/                 # v7 测试
│   └── v6_2/               # 历史回归测试
└── results/
│   ├── v7/                 # v7 状态、manifest 和结果
│   ├── v6_2/               # 本地历史结果；Git 忽略
│   └── v6_1/               # 本地历史结果；Git 忽略
```

## 文件归属

### 科学目标或假设改变

修改 `docs/plan/v7/plan.md`，并在提交说明中写明原因与结论影响。

### 阶段、依赖或执行顺序改变

修改 `docs/plan/v7/roadmap.md`。

### 当前实际状态改变

更新 `results/v7/STATUS.md`。状态必须区分 `NOT_RUN`、`RUNNING`、`COMPLETE`、`BLOCKED` 和 `NOT_IDENTIFIABLE`。

### 数据或外部来源改变

更新 v7 registry/manifest；新增外部生信数据时，交付前同步 `/home/huyudi/Infra/bioinf-data-index/`。

### 代码

- 可重放命令放 `scripts/v7/<stage>/`；
- 可复用函数放 `src/v7/<module>/`；
- 临时检查不放根目录，完成后进入 `scripts/archive/` 或删除；
- v6 代码不就地改造成 v7，先复制最小组件并注明来源。

### 结果

每个 v7 阶段目录至少包含：

- `PHASE_SUMMARY.md`
- `OUTPUT_MANIFEST.yaml`
- `DECISION_LOG.md`
- `EVIDENCE_AND_CONFLICTS.md`
- `NEXT_PHASE_READINESS.yaml`

大型矩阵、图像和 checkpoint 可以保留本地，但 manifest 必须记录路径、hash、参数、环境和可重算入口。

## 空间数据结构

空间 logical unit 必须形成：

```text
patient → block → section → region → spot/cell/bin
```

每层保留原始 ID、来源、坐标系统、比例尺、平台、counts 语义、图像/分割可用性、治疗/时间点/response 和 duplicate lineage。跨患者不得直接连图；相邻切片不得作为独立患者验证。

## 归档原则

- 文档和小型代码做物理归档；
- TB/GB 级数据和结果做逻辑归档，不移动，以免破坏绝对路径和可重放性；
- archive 表示不再作为当前入口，不表示结论有效；
- 禁止从历史结果复制结论而不重审其数据、模型和 claim boundary。

## Git 边界

普通 Git 只维护代码、测试、Markdown、配置和小型 manifest。下列内容不提交：raw/processed 大数据、h5ad/rds/qs/RData、大型 parquet/CSV、图像、checkpoint、scratch、缓存、日志和本地 IDE/agent 索引。
