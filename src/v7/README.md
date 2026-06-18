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
