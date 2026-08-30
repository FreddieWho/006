# Stage 1 registry builder

`build_stage1_registry.py` 是 v7 的第一阶段入口。它只读取小型 CSV/TSV 元数据和控制面文件，建立跨 `/006` 与只读 `/home/huyudi/013_spatial` 的身份、空间物理单位、患者—区块—切片、治疗/疗效完整性、重复血缘和 claim-specific role 账本。

默认运行：

```bash
python scripts/v7/registry/build_stage1_registry.py
```

检查已有输出：

```bash
python scripts/v7/registry/build_stage1_registry.py --validate-only
```

大矩阵、图像、h5ad 内容和模型均不会在此阶段读取。缺失或冲突字段保留为空/显式状态；脚本不根据文件名补患者、治疗或响应语义。输出排序固定，可在相同输入下复现。
