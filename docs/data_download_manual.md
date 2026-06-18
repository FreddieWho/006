# PD-1/PD-L1 免疫治疗单细胞/空间组学数据下载手册

**项目**: 006  
**更新时间**: 2026-05-05  
**数据来源**: 合并精简后的补全报告

---

## 目录

1. [数据概览](#数据概览)
2. [单细胞 RNA/TCR 数据集](#单细胞-rnatchr-数据集)
3. [空间转录组数据集](#空间转录组数据集)
4. [下载工具与命令](#下载工具与命令)
5. [数据预处理建议](#数据预处理建议)
6. [注意事项](#注意事项)

---

## 数据概览

| 类型 | 公开可下载 | 需申请 | 跟踪线索 | 已存在 |
|------|-----------|--------|----------|--------|
| 单细胞 RNA/TCR | 4 | 2 | 1 | 0 |
| 空间转录组 | 7 | 1 | 0 | 1 |
| **总计** | **11** | **3** | **1** | **1** |

**新增数据集**: 17个（不含已存在的 mendeley_skrx2fz79n）

---

## 单细胞 RNA/TCR 数据集

### 1. GSE236581 — 结直肠癌新辅助抗 PD-1

**基本信息**:
- 癌种: 结直肠癌 (CRC)
- 患者数: 22 例 (169 个单细胞样本)
- 数据类型: scRNA-seq + 配对 scTCR-seq
- 治疗方案: 新辅助抗 PD-1 单药 (Pembrolizumab/Sintilimab)
- 时间点: 治疗前、中、后多时间点
- 来源: GEO, Cancer Cell 2024

**下载方式**:
```bash
# 方法1: GEO FTP 直接下载
wget -r -np -nH --cut-dirs=5 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE236nnn/GSE236581/

# 方法2: 使用 GEOquery (R)
R -e 'library(GEOquery); getGEO("GSE236581", destdir="./data")'

# 方法3: 使用 pysradb (Python)
pip install pysradb
pysradb gse-to-srp GSE236581
pysradb sra-to-srr SRPXXXXXX
```

**文件结构**:
```
GSE236581/
├── GSMXXXXXXX_*.h5          # 10x HDF5 格式
├── GSMXXXXXXX_*.csv.gz      # 计数矩阵
├── GSMXXXXXXX_*_TCR.csv.gz  # TCR 注释
└── metadata.txt             # 样本元数据
```

**注意事项**:
- 数据包含原发肿瘤、邻近正常、外周血多时间点
- TCR 注释文件需单独下载
- 建议使用 CellRanger 或 STAR-solo 重新比对以获得 BAM 文件

---

### 2. HRA003591 — 食管鳞癌 PD-1 + 化疗

**基本信息**:
- 癌种: 食管鳞癌 (ESCC)
- 患者数: 10 例
- 数据类型: scRNA + scTCR
- 治疗: PD-1 抑制剂联合化疗
- 标签: TRG (肿瘤退缩分级), 治疗前后配对
- 来源: GSA-Human, Nat Commun 2024

**下载方式**:
```bash
# 需要机构邮箱申请 Controlled Access
# 步骤:
# 1. 访问 https://ngdc.cncb.ac.cn/gsa-human/browse/HRA003591
# 2. 点击 "Apply for Access"
# 3. 使用机构邮箱注册并提交申请
# 4. 等待审核 (通常 1-3 个工作日)
# 5. 获批后使用 Aspera 或 FTP 下载

# 下载命令 (获批后):
ascp -QT -l 300m -P 33001 -i ~/.aspera/connect/etc/asperaweb_id_dsa.openssh \
    era-fasp@fasp.cncb.ac.cn:/path/to/HRA003591 ./data/
```

**注意事项**:
- Controlled Access 数据需签署数据使用协议
- 仅限非商业研究用途
- 数据不可再分发

---

### 3. OMIX005710 — 食管鳞癌 PD-1 + 化疗

**基本信息**:
- 癌种: 食管鳞癌
- 患者数: 22 例
- 数据类型: scRNA-seq
- 治疗: PD-1 + 化疗
- 标签: pCR/MPR/IPR, 部分治疗前后配对
- 来源: OMIX/GSA, 2024

**下载方式**:
```bash
# 访问 OMIX 数据库
# https://ngdc.cncb.ac.cn/omix/browse/OMIX005710

# 根据页面提示下载
# 可能需要注册账号
```

**注意事项**:
- 报告标记为可获取
- 具体访问权限需查看页面说明

---

### 4. GSE221561 — 食管鳞癌 PD-1 + 化疗/放疗

**基本信息**:
- 癌种: 食管鳞癌
- 患者数: 7 例
- 数据类型: scRNA-seq
- 治疗: PD-1 + 化疗/放疗
- 标签: TRG, 部分配对
- 来源: GEO, 2025

**下载方式**:
```bash
# GEO FTP 直接下载
wget -r -np -nH --cut-dirs=5 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE221nnn/GSE221561/

# 或使用 GEOquery
R -e 'library(GEOquery); getGEO("GSE221561", destdir="./data")'
```

---

### 5. HRA007492 — 宫颈癌 PD-1 + 化疗

**基本信息**:
- 癌种: 宫颈癌
- 患者数: 5 例 (10 个样本)
- 数据类型: scRNA + scTCR
- 治疗: PD-1 + 化疗
- 标签: MPR/NMPR, 治疗前后配对
- 来源: GSA-Human, Cancer Res 2024

**下载方式**:
```bash
# 需要机构邮箱申请 Controlled Access
# 步骤同 HRA003591

# 访问: https://ngdc.cncb.ac.cn/gsa-human/browse/HRA007492
```

**注意事项**:
- 与空间转录组 HRA007492_ST 来自同一研究
- 可同时申请 scRNA 和空间数据

---

### 6. PRJNA932556 — 结直肠癌 MSI-H 抗 PD-1 单药

**基本信息**:
- 癌种: 结直肠癌 (MSI-H)
- 患者数: 6 例
- 数据类型: scRNA-seq
- 治疗: PD-1 单药
- 标签: 敏感/耐药 (S/R), 非配对
- 来源: SRA/BioProject, BMC Med 2023

**下载方式**:
```bash
# 方法1: SRA Toolkit
prefetch PRJNA932556
fasterq-dump --split-files SRRXXXXXXX

# 方法2: pysradb
pysradb metadata PRJNA932556

# 方法3: 直接访问 BioProject 页面
# https://www.ncbi.nlm.nih.gov/bioproject/PRJNA932556
```

---

### 7. GSE256326 — 肾髓质癌 nivolumab + ipilimumab

**基本信息**:
- 癌种: 肾髓质癌 (RMC, 极罕见)
- 样本: 15 个样本 (9 个 scRNA)
- 治疗方案: nivolumab + ipilimumab (PD-1 + CTLA-4)
- 标签: 治疗前 vs 治疗后高进展
- 来源: GEO

**下载方式**:
```bash
wget -r -np -nH --cut-dirs=5 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE256nnn/GSE256326/
```

**注意事项**:
- 包含 7 例免疫治疗 naïve 和 2 例 nivolumab+ipilimumab 后高进展病灶
- 可用于极端耐药对照研究

---

### 8. GSE314072 — 肾细胞癌 (仅作参考)

**基本信息**:
- 癌种: 肾细胞癌 (RCC)
- 数据类型: scRNA + scTCR, 24 样本
- **无 PD-1/PD-L1 治疗信息**

**下载方式**:
```bash
wget -r -np -nH --cut-dirs=5 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE314nnn/GSE314072/
```

**注意事项**:
- ⚠️ 不能直接作为 PD-1 轴治疗队列
- 仅作免疫图谱参考

---

## 空间转录组数据集

### 9. mendeley_skrx2fz79n — 肝癌 Visium (已存在)

**基本信息**:
- 平台: Visium v1
- 癌种: 肝癌 (HCC)
- 患者数: 8
- 治疗: PD-1 单药
- 标签: 响应/不应答, 非配对
- 来源: Mendeley

**状态**: 已在 data_collection.csv 中 (第12-13行)

---

### 10. GSE289745 — 皮肤鳞癌 Visium

**基本信息**:
- 癌种: 皮肤鳞癌 (cSCC)
- 样本数: 11 (Visium)
- 治疗关联: PD-1/PD-L1 阻断相关微环境研究
- 标签: 需从原文手动对应样本治疗信息
- 来源: GEO, Nature 2025

**下载方式**:
```bash
wget -r -np -nH --cut-dirs=5 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE289nnn/GSE289745/
```

**注意事项**:
- GEO 元数据缺少逐样本的前/后或应答标注
- 需要从原文手动对应样本治疗信息

---

### 11. GSE291246 — 基底细胞癌 Xenium

**基本信息**:
- 平台: 10x Xenium (原位空间)
- 癌种: 基底细胞癌
- 患者数: 35
- 治疗: PD-1, 有临床响应标签, 配对
- 来源: GEO, PNAS 2025

**下载方式**:
```bash
wget -r -np -nH --cut-dirs=5 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE291nnn/GSE291246/
```

---

### 12. GSE235672 — 胶质母细胞瘤 Visium

**基本信息**:
- 癌种: 胶质母细胞瘤
- 患者数: 20
- 治疗: PD-1, 有响应/不应答标签, 部分配对
- 来源: GEO, Nat Cancer 2023

**下载方式**:
```bash
wget -r -np -nH --cut-dirs=5 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE235nnn/GSE235672/
```

---

### 13. GSE238264 — 肝癌 Visium

**基本信息**:
- 癌种: 肝癌 (HCC)
- 患者数: 7
- 治疗: PD-1 + TKI, 有病理响应标签, 非配对
- 来源: GEO, Genome Med 2023

**下载方式**:
```bash
wget -r -np -nH --cut-dirs=5 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE238nnn/GSE238264/
```

---

### 14. PRJCA039752 — 食管鳞癌 Visium HD

**基本信息**:
- 平台: Visium HD
- 癌种: 食管鳞癌
- 患者数: 4
- 治疗: PD-1 + 化疗, 有响应/不应答标签, 配对
- 来源: CNCB, 2024

**下载方式**:
```bash
# 访问 CNCB
# https://ngdc.cncb.ac.cn/gsa/browse/PRJCA039752

# 使用 Aspera 下载
ascp -QT -l 300m -P 33001 -i ~/.aspera/connect/etc/asperaweb_id_dsa.openssh \
    era-fasp@fasp.cncb.ac.cn:/path/to/PRJCA039752 ./data/
```

---

### 15. HRA007492_ST — 宫颈癌 Visium

**基本信息**:
- 癌种: 宫颈癌
- 患者数: 5
- 治疗: PD-1 + 化疗, MPR/NMPR 标签, 配对
- 来源: GSA-Human, Cancer Res 2024

**下载方式**:
```bash
# 需要机构邮箱申请 Controlled Access
# 与 scRNA 数据 HRA007492 来自同一研究
# 访问: https://ngdc.cncb.ac.cn/gsa-human/browse/HRA007492
```

---

### 16. NCT02451982 — 胰腺癌 Visium FFPE

**基本信息**:
- 癌种: 胰腺癌
- 患者数: 12
- 治疗: PD-1 + 疫苗, 有病理响应标签, 配对
- 来源: Cancer Cell 2024

**下载方式**:
```bash
# 需要从原文获取具体 GEO/SRA 编号
# 临床试验编号: NCT02451982
# 建议搜索: https://clinicaltrials.gov/ct2/show/NCT02451982
```

---

### 17. GSE177043 — 三阴性乳腺癌 Visium

**基本信息**:
- 癌种: TNBC
- 样本数: >20
- 治疗: 抗 PD-1, 有明确的应答/不应答分组
- 平台: 10x Visium
- 来源: GEO

**下载方式**:
```bash
wget -r -np -nH --cut-dirs=5 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE177nnn/GSE177043/
```

**注意事项**:
- 经典旧队列, 空间数据质量高且应答标签完备
- 可作为重要的补充队列

---

### 18. HRA011002 — NSCLC 脑转移 (跟踪线索)

**基本信息**:
- 癌种: NSCLC 脑转移
- 数据类型: scRNA
- 治疗: 免疫治疗
- 来源: GSA-Human, Nat Commun 2026

**下载方式**:
```bash
# 目前尚未完全公开
# 访问: https://ngdc.cncb.ac.cn/gsa-human/browse/HRA011002
# 持续关注是否转为公开或 Controlled Access
```

**注意事项**:
- 跟踪线索，暂不纳入正式数据集
- 整合了部分已有的 NSCLC ICB 队列

---

## 下载工具与命令

### 1. GEO 数据下载

```bash
# 使用 wget 递归下载
wget -r -np -nH --cut-dirs=5 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE236nnn/GSE236581/

# 使用 GEOquery (R)
R -e 'library(GEOquery); getGEO("GSE236581", destdir="./data")'

# 使用 pysradb (Python)
pip install pysradb
pysradb gse-to-srp GSE236581
pysradb sra-to-srr SRPXXXXXX
```

### 2. SRA 数据下载

```bash
# 安装 SRA Toolkit
conda install -c bioconda sra-tools

# 下载 SRA 数据
prefetch PRJNA932556
fasterq-dump --split-files SRRXXXXXXX

# 转换为 FASTQ
fastq-dump --split-files SRRXXXXXXX
```

### 3. GSA-Human 数据下载 (Controlled Access)

```bash
# 1. 申请访问权限
# 访问 https://ngdc.cncb.ac.cn/gsa-human/
# 使用机构邮箱注册并提交申请

# 2. 下载 (获批后)
ascp -QT -l 300m -P 33001 -i ~/.aspera/connect/etc/asperaweb_id_dsa.openssh \
    era-fasp@fasp.cncb.ac.cn:/path/to/HRA003591 ./data/

# 或使用 FTP
wget -r -np -nH ftp://ftp.cncb.ac.cn/gsa-human/HRA003591/
```

### 4. Mendeley 数据下载

```bash
# 访问 Mendeley Data 页面
# https://data.mendeley.com/datasets/skrx2fz79n

# 或使用 curl
curl -L "https://data.mendeley.com/datasets/skrx2fz79n/draft" -o skrx2fz79n.zip
```

---

## 数据预处理建议

### 1. 单细胞数据

```bash
# 使用 CellRanger 比对
cellranger count --id=sample1 \
    --transcriptome=/path/to/refdata-gex-GRCh38-2020-A \
    --fastqs=/path/to/fastqs \
    --sample=sample1

# 或使用 STAR-solo
STAR --runThreadN 8 \
    --genomeDir /path/to/genomeDir \
    --readFilesIn sample1_R2.fastq.gz sample1_R1.fastq.gz \
    --readFilesCommand zcat \
    --soloType CB_UMI_Simple \
    --soloCBwhitelist /path/to/737K-august-2016.txt
```

### 2. 空间转录组数据

```bash
# 使用 SpaceRanger (10x Visium)
spaceranger count --id=sample1 \
    --transcriptome=/path/to/refdata-gex-GRCh38-2020-A \
    --fastqs=/path/to/fastqs \
    --sample=sample1 \
    --slide=V19J01-123 \
    --area=A1

# 使用 Xenium Ranger (10x Xenium)
xeniumranger import --id=sample1 \
    --xenium-bundle=/path/to/xenium-bundle
```

### 3. 数据整合

```bash
# 使用 Seurat (R)
library(Seurat)
obj <- Read10X(data.dir = "./data/sample1/outs/filtered_feature_bc_matrix")
seurat_obj <- CreateSeuratObject(counts = obj, project = "sample1")

# 使用 Scanpy (Python)
import scanpy as sc
adata = sc.read_10x_mtx("./data/sample1/outs/filtered_feature_bc_matrix")
```

---

## 注意事项

### 1. 数据访问权限

- **公开数据** (GEO/SRA/CNCB): 可直接下载，无需申请
- **Controlled Access** (GSA-Human): 需通过机构邮箱申请，签署数据使用协议
- **跟踪线索**: 数据尚未完全公开，需持续关注

### 2. 数据使用协议

- Controlled Access 数据仅限非商业研究用途
- 数据不可再分发
- 需在发表时引用原始数据来源

### 3. 存储空间估算

| 数据类型 | 单样本大小 | 总计 (17 个新数据集) |
|----------|-----------|-------------------|
| scRNA (FASTQ) | 10-50 GB | 170-850 GB |
| scRNA (BAM) | 10-30 GB | 170-510 GB |
| ST-Visium (FASTQ) | 5-20 GB | 85-340 GB |
| ST-Visium (BAM) | 5-15 GB | 85-255 GB |

**建议**: 准备至少 2 TB 存储空间用于完整下载和预处理。

### 4. 计算资源

- **CellRanger/SpaceRanger**: 需要 8+ CPU 核心, 64+ GB 内存
- **STAR-solo**: 需要 8+ CPU 核心, 32+ GB 内存
- **数据整合**: 建议使用 GPU 加速 (如 NVIDIA A100)

### 5. 数据质量检查

下载后建议进行以下检查:
```bash
# 检查文件完整性
md5sum -c checksums.md5

# 检查 FASTQ 格式
fastqc sample1_R1.fastq.gz sample1_R2.fastq.gz

# 检查 BAM 文件
samtools flagstat sample1.bam
```

---

## 附录: 空间转录组平台对比

| 平台 | 分辨率 | 基因检测 | 成本 | 适用场景 |
|------|--------|----------|------|----------|
| Visium | 55 μm | 全转录组 | 中 | 组织结构分析 |
| Visium HD | 2 μm | 全转录组 | 高 | 单细胞级空间分析 |
| Xenium | 亚细胞 | 靶向 panel | 高 | 原位单细胞分析 |
| CosMx | 亚细胞 | 靶向 panel | 高 | 多重蛋白+RNA |

---

## 参考文献

1. Cancer Cell 2024 - GSE236581 结直肠癌新辅助抗 PD-1
2. Nat Commun 2024 - HRA003591 食管鳞癌 PD-1 + 化疗
3. Cancer Res 2024 - HRA007492 宫颈癌 PD-1 + 化疗
4. BMC Med 2023 - PRJNA932556 结直肠癌 MSI-H 抗 PD-1
5. Nature 2025 - GSE289745 皮肤鳞癌 Visium
6. PNAS 2025 - GSE291246 基底细胞癌 Xenium
7. Nat Cancer 2023 - GSE235672 胶质母细胞瘤 Visium
8. Genome Med 2023 - GSE238264 肝癌 Visium
9. Cancer Cell 2024 - NCT02451982 胰腺癌 Visium FFPE

---

**手册维护者**: Hermes Agent  
**最后更新**: 2026-05-05  
**版本**: 1.0  
**项目**: 006
