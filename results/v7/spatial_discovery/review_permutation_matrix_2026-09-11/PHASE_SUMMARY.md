# Stage 4｜999次空间置换

执行：COMPUTATION_COMPLETE。71捕获、43患者、2769特征×捕获项；2752项各999次，17项观测不足，不补零或缩减次数。原19次结果保留为历史探索。

FDR族为每捕获39特征；结论仅限捕获内自相关。真实topology增量、跨患者推断及独立结构GT仍未完成，科学阶段OPEN。

12个完成捕获经输入hash、种子、39特征与次数检查后从中断运行复用；其余捕获采用共享缺失掩码的向量置换。每特征边际null保持一致，联合随机数耦合不同，已在null_kernel列记录。两个中断目录保留源代码快照及状态。

重算入口：`scripts/v7/review/run_closure_remediation.py spatial_permutation --output <new_directory>`；可不使用resume从头运行全量。原启动命令、输入和输出hash见RUN_RECEIPT.json。
