# v7 tests

Stage 1–3 首批测试应覆盖：

- patient/block/section identity conservation；
- 同患者不跨 split；
- coordinate scale 与单位；
- counts/normalized layer 语义；
- ground-truth 和派生变量隔离；
- duplicate lineage；
- treatment/timepoint signal 不被预处理误删；
- spot/cell 不被当作独立患者。

当前尚无 v7 科学实现，因此没有伪造空测试结果。
