# 006 v7 工作规则

## 当前依据

- 科学目标：`docs/plan/v7/plan.md`
- 执行顺序：`docs/plan/v7/roadmap.md`
- 当前状态：`results/v7/STATUS.md`
- 目录规范：`STRUCTURE.md`

v7 是探索性生信研究，不使用预注册、response unlock 或旧 frozen gate 作为新工作的硬约束。历史约束只作 provenance。

## 科学纪律

- 优先回答“免疫程序如何形成空间屏障、context 如何改写、X 是否与 repair-compatible rewiring 相容”。
- patient 是临床独立单位；block/section/spot/cell 是嵌套观测。
- 无 response 的空间数据可用于发现和结构验证，不能冒充疗效证据。
- 图注意力、邻近、表达相关和模型 transport 不等于因果关系。
- 先比较透明基线；复杂模型只有提供稳定增量才采用。
- 核心计算不能静默跳过、替代或缩减；无法运行时保持阶段未完成并说明结论影响。

## 数据与复用

- 先查 manifest、身份链、counts 语义和重复来源，再复用 v6/013 资产。
- `/home/huyudi/013_spatial` 是只读外部项目来源；不修改其数据、代码或运行任务。
- v6 的数据接口和测量资产可复用；旧 Phase8、SRB/FM07-only 与 hash perturbation 分数不继承为 v7 事实。
- 新增外部生信数据时，交付前更新 `/home/huyudi/Infra/bioinf-data-index/`。

## 工程与交付

- v7 命令入口放 `scripts/v7/`，实现放 `src/v7/`，测试放 `tests/v7/`。
- 阶段结束提供 summary、manifest、decision、evidence/conflicts 和 readiness。
- 大型数据、结果和 checkpoint 不进普通 Git；提交路径、hash 和可重算入口。
- 历史资产可归档或废弃，但不破坏 provenance，不把归档结果包装成现行结论。
