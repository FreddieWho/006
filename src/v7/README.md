# v7 reusable implementation

未来模块按职责拆分：

- `registry`：patient/block/section 和 duplicate lineage；
- `spatial_io`：平台专属读取、坐标、尺度和图像接口；
- `spatial_stats`：距离、邻接、边界、permutation 和患者级汇总；
- `ontology`：module/cell-state 跨模态映射；
- `clinical`：环境内模型和 meta-analysis；
- `repair`：paired displacement 与 X-class compatibility；
- `evidence`：Evidence Card 与 conflict register。

当前为空结构说明，不代表已经实现模型。
# v7 source modules

`ontology/` contains Stage 2 measurement and cross-modal vocabulary code.
`spatial_io/`, `spatial_stats/` and `spatial_pipeline.py` implement the Stage 3
input contract. `spatial_discovery.py` is the response-blind Stage 4 discovery
pass, `clinical_baseline.py` is the patient-level Stage 5 clinical anchor, and
`context_analysis.py` performs the Stage 6 context decomposition and bridge
audit, and `repair_analysis.py` performs the Stage 7 patient-level molecular
repair baseline while fail-closing spatial rewiring when paired assets are
missing.

Stage 4 marker scores are explicitly composition proxies, not inferred cell
fractions. Stage 5 models read audited response fields only in the clinical
anchor path; neither module converts spatial associations into causal claims.
