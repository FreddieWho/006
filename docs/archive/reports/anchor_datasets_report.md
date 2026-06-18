# Anti-PD-1/PD-L1 Clean Anchor 数据集检索报告

**检索日期**: 2026-06-26
**检索范围**: GEO, Zenodo, PubMed, TISCH2, TIGER, CellxGene, Broad SCP, IMMUcan
**检索代理**: 8个并行搜索代理 + 6个深度验证代理
**总搜索量**: 15+ 数据库, 50+ 关键词组合, 200+ 候选初步筛选

---

## EXECUTIVE SUMMARY

本次检索在已有 A/B/C 类清单基础上，通过系统性多数据库并行搜索，发现 **1个全新 A 类 clean anchor 数据集** 和 **4 个 B 类 near-miss 候选**。

| 分类 | 数量 | 关键发现 |
|------|------|---------|
| **A 类 (Clean Anchor)** | 1 new | GSE301741 (HNSCC, pembro mono, 16 pts, 137k cells, pre+post) |
| **B 类 (Near-Miss)** | 4 new | GSE286827, GSE288199, GSE274975, SCP3121 |
| **C 类 (Excluded)** | 2 new | GSE272347, SCP3341 |

---

## A 类：可直接补主干的 Clean Anchor

### A-NEW-1: GSE301741 (HIGHEST PRIORITY)

```
dataset_id: GSE301741
repository: NCBI GEO
cancer_type: Head and neck squamous cell carcinoma (HNSCC), HPV-negative, Stage III/IV
modality: scRNA-seq (10x Genomics 3'/5') + TCR-seq (VDJ)
tissue_type: TUMOR - Fresh primary tumor (Oral Cavity, Oropharynx, Larynx)
treatment_regimen: Pembrolizumab (anti-PD-1) monotherapy, 200mg q3wk x2 cycles neoadjuvant
monotherapy_or_clean_ici: YES - Confirmed clean anti-PD-1 monotherapy
response_definition: Pathologic Tumor Response (pTR): pTR-2 (>50% regression) or pTR-1 (10-50%) = R; <10% = NR
response_groups_counts: 5 responders + 5 non-responders (matched pre/post cohort of 10); all 16 patients have labels in Seurat object
timepoint: Pre (baseline, n=11) + Post (after neoadjuvant, n=16), 11 matched pre/post
patient_count: 16
sample_count: 58 GEO libraries (20 pre + 38 post)
cells_profiled: 137,020 QC-passing cells
processed_files_available: YES
  - Seurat RDS object (12.4 Gb) with FULL metadata
  - Individual H5 matrices (586.5 Mb tar)
  - Series matrix, SOFT, MINiML
direct_download_urls:
  - Seurat RDS: https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE301741&format=file&file=GSE301741%5FSeurat%5FObject%5FQCpass%5F137020cells%5FwithMetaData%2Erds
  - FTP mirror: ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE301nnn/GSE301741/suppl/GSE301741_Seurat_Object_QCpass_137020cells_withMetaData.rds
  - RAW tar: https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE301741&format=file
metadata_file_with_response: YES - GSE301741_Seurat_Object_QCpass_137020cells_withMetaData.rds contains cell-level response labels
metadata_file_with_patient_sample_mapping: YES - Same Seurat RDS contains Patient ID and Timepoint for each cell; Series Matrix maps GSM to patient IDs
request_only: NO
raw_fastq_required: NO - Processed matrices and Seurat object directly available
trunk_eligible: YES - ALL hard criteria satisfied
if_not_trunk_eligible_why: N/A
confounder_risks:
  - Chemotherapy: NONE in neoadjuvant window (pembro only pre-surgery)
  - TKI: NONE
  - Combination ICI: NONE
  - Minor: Mix of 3' and 5' chemistry; different cell sorting strategies (CD45/CD3/unsorted)
why_it_is_better_than_GSE236581: Confirmed pembro monotherapy vs unknown regimen; 137k cells with malignant+immune; pre+post paired design; response labels in downloadable Seurat object; practice-changing KEYNOTE-689 Phase 2 cohort
why_it_is_better_than_bulk-only_candidates: scRNA-seq with cellular resolution + TCR-seq; fresh tumor not FFPE; includes malignant cells; pathologic response from surgical specimen
publication: Oliveira et al., Cell Reports Medicine 2026 Apr 21;7(4):102715. PMID: 41923630




















































clinical_trial: Phase 2 multicenter trial (Cohort 2) -> basis for KEYNOTE-689 Phase 3
SRA: PRJNA1283925
BioProject: PRJNA1283925
```

**A 类验证清单**:
- [x] Public direct download: YES
- [x] Non-FASTQ-only: YES (Seurat RDS + H5 matrices)
- [x] Anti-PD-1 monotherapy: YES (pembrolizumab only)
- [x] Bilateral response: YES (5R + 5NR matched)
- [x] Machine-readable response labels: YES (in Seurat metadata)
- [x] Tumor-derived: YES (fresh primary tumor)

---

## B 类：高价值 Near-Miss 候选

### B-NEW-1: GSE286827 (CONDITIONAL - monotherapy subset)

```
dataset_id: GSE286827
repository: NCBI GEO
cancer_type: HNSCC (oral cavity, oropharynx, hypopharynx, larynx)
modality: scRNA-seq (10x 5' v2) + scTCR-seq
tissue_type: TUMOR - Fresh tumor biopsies (NOT PBMC)
treatment_regimen: Durvalumab (anti-PD-L1) monotherapy OR Durvalumab + Tremelimumab (anti-CTLA-4)
monotherapy_or_clean_ici: CONDITIONAL - 13/29 (45%) monotherapy, separable via Neoadj_type column
response_definition: Pathological tumor regression (>=50% = R, <50% = NR)
response_groups_counts: All 29 pts: 8R + 21NR; Monotherapy subset: 4R + 9NR (or 5R+8NR per Table S1)
timepoint: Pre + Post, 27 patients paired
patient_count: 29 total; 13 monotherapy subset
sample_count: 87 total (57 GEX + 30 TCR)
processed_files_available: YES
  - 7 cell-type-specific countdata RDS files (Bcells, CD4T, CD8T, Malignant, TAMs, TumorSpecific_CD8T, immunecells_broad)
  - 7 metadata RDS files (all contain response_final column: R/NR)
  - RAW tar (2.5 Gb MTX/TSV)
direct_download_urls:
  - Count data: https://ftp.ncbi.nlm.nih.gov/geo/series/GSE286nnn/GSE286827/suppl/GSE286827_countdata_*.rds.gz (7 files)
  - Metadata: https://ftp.ncbi.nlm.nih.gov/geo/series/GSE286nnn/GSE286827/suppl/GSE286827_metadata_*.rds.gz (7 files)
  - RAW: https://ftp.ncbi.nlm.nih.gov/geo/series/GSE286nnn/GSE286827/suppl/GSE286827_RAW.tar
  - Table S1: https://europepmc.org/api/fulltextRepo?pmcId=PMC12629806



















        &type=FILE&fileName=mmc2.xlsx
metadata_file_with_response: YES - All 7 metadata RDS files contain response_final column (R/NR)
metadata_file_with_patient_sample_mapping: YES - All metadata contain Patient and Time columns
request_only: NO
raw_fastq_required: NO
trunk_eligible: YES (IF using monotherapy subset only, Neoadj_type == "D")
if_not_trunk_eligible_why: 55% of patients received dual ICI (D+T); must exclude dual ICI patients for clean monotherapy analysis
confounder_risks: Dual ICI ratio concern (55% combo); small monotherapy subset (n=13); response imbalance (4R vs 9NR); response discrepancy P1 between Table S1 and metadata
why_it_is_better_than_GSE236581: Larger than typical excluded datasets; clean separable monotherapy subset; tumor tissue; pre+post paired; pathological response; cell-type-specific matrices with metadata
publication: Cha J et al., Cell Reports Medicine 2025 Oct 21;6(10):102408. PMID: 41045934


















clinical_trial: NCT03737968
```

**进入 A 类条件**: 仅使用 Neoadj_type == "D" 的 13 个 monotherapy patients (4R+9NR)。需要在分析时排除所有 D+T 患者。

---

### B-NEW-2: GSE288199 (3-arm combination trial)

```
dataset_id: GSE288199
repository: NCBI GEO
cancer_type: HNSCC (Stage III/IVa)
modality: scRNA-seq + CITE-seq + scTCR-seq
tissue_type: TUMOR - FACS-sorted CD45+CD3+ TIL
treatment_regimen: 3-arm trial - Nivolumab mono (n=14) / Nivolumab+Relatlimab (n=15) / Nivolumab+Ipilimumab (n=12)
monotherapy_or_clean_ici: NO - 3-arm combination therapy trial; monotherapy arm too small
response_definition: pTR-0 (91-100% viable), pTR-1 (51-90%), pTR-2 (<=50% viable)
response_groups_counts: 42 randomized; Mono arm: ~12 evaluable, pTR-2 rate only 8.4%
timepoint: Pre + Post, 20 matched pre/post for scRNA-seq
patient_count: 42 randomized; 35 with scRNA-seq
sample_count: 372,914 CD45+CD3+ TIL
cells_profiled: 372,914 TIL
processed_files_available: YES - MTX format matrices (barcodes/features/matrix per sample)
direct_download_urls: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE288199 (180 individual sample files)
metadata_file_with_response: YES - pTR labels in sample names and characteristics
metadata_file_with_patient_sample_mapping: YES - Sample naming: GSE288199_HN{patient}_{pre/post}_T_...
request_only: NO
raw_fastq_required: NO
trunk_eligible: NO
if_not_trunk_eligible_why:
  1. 3-arm trial design - 2/3 arms are combination therapies
  2. Monotherapy arm too small (n=14 randomized, ~12 evaluable)
  3. HPV+ patients all on mono arm (confounding)
  4. Trial designed to compare combos; extracting only mono loses primary context
confounder_risks: Multi-arm combination design; anti-LAG-3 and anti-CTLA-4 combination arms; small mono cohort
why_it_is_better_than_GSE236581: High-quality Cancer Cell 2025 paper; scRNA+CITE+TCR multi-modal; 372k TIL; pathologic response; could be used IF trunk criteria expanded to include combo trials
publication: Li H et al., Cancer Cell 2025 Apr 14;43(4):757-775.e8. PMID: 40086437


















clinical_trial: NCT04080804
SRA: PRJNA1216557
```

**为何是 B 类而非 C 类**: 虽然是组合试验，但数据质量极高，如果未来主干标准放宽（如允许组合试验作为 support 数据集），可以重新评估。Mono arm (n=14) 理论上可提取但太小。

---

### B-NEW-3: GSE274975 (BULK RNA-seq anchor)

```
dataset_id: GSE274975
repository: NCBI GEO
cancer_type: NSCLC (adenocarcinoma 24, squamous 27, mixed 9, neuroendocrine 1)
modality: BULK RNA-seq (FFPE tissue, NOT scRNA-seq)
tissue_type: FFPE tumor tissue (>50% tumor cell content)
treatment_regimen: Anti-PD-1/PD-L1 ICI - 30 monotherapy / 31 combination
monotherapy_or_clean_ici: YES for monotherapy subset (30 pts: pembro 11, nivo 9, tislelizumab 9, atezo 1)
response_definition: RECIST v1.1 (CR/PR/SD/PD)
response_groups_counts: All 61: CR 3, PR 12, SD 25, PD 17, missing 4; Monotherapy subset: ~5 R, ~22 NR (22/30 have response data)
timepoint: Pretreatment (prior to ICI)
patient_count: 61 total; 30 monotherapy subset
sample_count: 61
processed_files_available: YES - Gene-level raw counts (2.7 Mb)
direct_download_urls:
  - Raw counts: https://ftp.ncbi.nlm.nih.gov/geo/series/GSE274nnn/GSE274975/suppl/GSE274975_raw_counts.tsv.gz
  - Clinical Table 1: https://public-pages-files-2025.frontiersin.org/articles/1493877/file/Table_1.xlsx/1493877_table_1/1
metadata_file_with_response: PARTIALLY - Response labels in paper's Supplementary Table 1 (NOT embedded in GEO metadata); must merge by Sample_ID
metadata_file_with_patient_sample_mapping: YES - GEO sample IDs (OB_pat_LuC_XXX) perfectly match Table 1 Sample_IDs
request_only: NO
raw_fastq_required: NO - Raw counts directly available
trunk_eligible: NO - BULK RNA-seq, not scRNA-seq
if_not_trunk_eligible_why: BULK RNA-seq modality (FFPE tissue, gene-level counts); not single-cell resolution
confounder_risks: FFPE tissue (lower RNA quality); 8/30 monotherapy patients lack response labels; mix of histotypes; mix of primary (49) and metastatic (12) sites; includes non-FDA-approved drug (tislelizumab)
why_it_is_better_than_GSE236581: Monotherapy subset cleanly separable; pretreatment samples; RECIST labels; raw counts available; Russian cohort adds ethnic diversity; 30 monotherapy patients usable
publication: Poddubskaya EV et al., Frontiers in Immunology 2024 Dec 12;15:1493877. PMID: 39723204


















SRA: PRJNA1148665
```

**进入 A 类条件**: 不可能（bulk RNA-seq 不可能成为 scRNA anchor）。但作为 BULK anchor 可用：使用 30 个 monotherapy patients，合并 Table 1 临床数据与表达矩阵。

---

### B-NEW-4: SCP3121 / IMvigor210 (Conditional on response linkage)

```
dataset_id: SCP3121
repository: Broad Single Cell Portal (+ GEO GSE313827 for CAF cell lines only)
cancer_type: Pan-cancer (breast, renal, colorectal, liver, lung, ovarian, pancreatic, prostate, bladder + normal)
modality: scRNA-seq
tissue_type: TUMOR (cancer/epithelial + CAF cells confirmed)
treatment_regimen: Study focuses on atezolizumab (anti-PD-L1) monotherapy (IMvigor210/211)
monotherapy_or_clean_ici: PARTIALLY - sc data is pan-cancer meta-analysis, not single trial
response_definition: RECIST v1.1 (at BULK level only, via IMvigor210CoreBiologies R package)
response_groups_counts: UNKNOWN at single-cell level; BULK: 2,800 patients with RECIST data
timepoint: UNKNOWN (not confirmed pre-treatment)
patient_count: 200 samples (mix of patients + healthy donors); BULK: ~2,800
sample_count: 200 scRNA libraries
cells_profiled: 884,573 cells (15,007 genes)
processed_files_available: YES on Broad SCP (requires login); GEO has only 18 CAF cell line bulk samples
direct_download_urls:
  - Broad SCP: https://singlecell.broadinstitute.org/single_cell/study/SCP3121 (requires Google auth)
  - GEO (CAF only): https://ftp.ncbi.nlm.nih.gov/geo/series/GSE313nnn/GSE313827/suppl/GSE313827_cpm_counts.tsv.gz
metadata_file_with_response: NO - Response labels NOT linked to single-cell samples; bulk data has response via R package
metadata_file_with_patient_sample_mapping: PARTIALLY - donor_id and biosample_id exist but not verified
request_only: NO - Public but requires Broad SCP login
raw_fastq_required: NO
trunk_eligible: NO
if_not_trunk_eligible_why:
  1. Response labels (RECIST) NOT linked to individual single-cell samples (CRITICAL)
  2. Pan-cancer meta-analysis, not single focused trial
  3. Pre-treatment status unconfirmed
  4. Download requires authentication
confounder_risks: Multi-study heterogeneity; healthy donor samples mixed; unknown treatment timing across compiled studies
why_it_is_better_than_GSE236581: 884k cells from 200 samples; tumor tissue; clean anti-PD-L1 focus; associated with well-characterized IMvigor210 clinical trial
publication: Roels J et al., Cancer Research 2026 Jun 1;86(11):2592-2605. PMID: 41784658




















GEO: GSE313827 (CAF cell lines only, NOT scRNA-seq data)
```

**进入 A 类条件**: 需要联系作者 (Juliette Roels, roels.juliette@gene.com) 获取 donor-level metadata 文件，将 RECIST response 关联到 single-cell samples。如果确认response labels可用且为pre-treatment样本，可升级为 A 类。

---

## C 类：明确不应纳入主干的数据

### C-NEW-1: GSE272347

| 排除原因 | 详情 |
|---------|------|
| **TKI+ICI 组合治疗** | Cabozantinib (TKI) + nivolumab; GitHub代码显示混合regimens ('anti-PD1+TKI' 和 'anti-PD1') |
| **多数为PBMC** | 仅1例患者(OT6)有肿瘤scRNA-seq; 其余为外周血 |
| **仅post-treatment** | 所有样本为新辅助治疗后手术切除标本，无pre-treatment |
| **样本量小** | 仅14例患者 |
| **TLS-density为主要分析轴** | 主要比较TLS-high vs TLS-low，而非responder vs non-responder |
| **分类** | C类 - 排除 |

- Paper: Shu DH et al., Nature Immunology 2024 Nov;25(11):2110-2123. PMID: 39455893




















- 注意: 与GSE206325完全不同（不同机构、不同药物、不同科学问题）

---

### C-NEW-2: SCP3341

| 排除原因 | 详情 |
|---------|------|
| **未确认单药治疗** | 使用通用"ICB-treated"标签，可能包含多种regimens (nivolumab, pembrolizumab, vaccine combos) |
| **无二元响应标签** | 仅survival-based (MES subtype预测OS)，无R/NR分类 |
| **GBM ICB响应极差** | 已知对ICB单药响应极低的癌种，不适合构建通用anchor |
| **snRNA-seq非scRNA-seq** | 单核RNA-seq (nuclei而非whole cells) |
| **分类** | C类 - 排除 |

- Paper: Ghannam JY et al., PMID: 42237038



















        , Dana-Farber/Harvard/Broad

---

## 已确认 A 类（原有清单）更新状态

| 数据集 | 状态 | 备注 |
|--------|------|------|
| GSE206325 (HCC) | **仍有效** | 本次检索未找到替代/更新版本 |
| GSE123813 (BCC) | **仍有效** | TISCH2确认可下载; Zenodo 10407125有RDS |
| GSE176021 (NSCLC) | **仍有效** | CheckMate 159 MPR/NMPR |

---

## 已确认 B 类（原有清单）更新状态

| 数据集 | 状态 | 备注 |
|--------|------|------|
| GSE120575 (Mel) | 仍有效 | 25% combo率 |
| GSE115978 (Mel) | 仍有效 | combo |
| GSE202069 (HCC bulk) | 仍有效 | labels需从PDF解析 |
| GSE145996 (Mel bulk) | 仍有效 | 21R+31NR |
| GSE205335 (NSCLC) | 仍有效 | ICI +/- combo |
| GSE235863 (HCC) | 仍有效 | anti-PD-1+lenvatinib |
| GSE243013 (NSCLC) | 仍有效 | anti-PD-1+chemo, 234 pts |
| GSE232240 (HNSCC) | 仍有效 | anti-PD-1+CTLA-4 |
| GSE279750 (HCC bulk) | 仍有效 | anti-PD-L1 combo |
| GSE151530 (HCC/iCCA) | 仍有效 | iCCA subset待确认 |

---

## 数据完整性总结

### A 类 Clean Anchor 汇总（现有+新增）

| # | Dataset | Cancer | Modality | Treatment | Patients | Cells | Response | Pre/Post |
|---|---------|--------|----------|-----------|----------|-------|----------|----------|
| 1 | GSE206325 | HCC | scRNA | Cemiplimab mono | ~50 | ~510K | Yes | Baseline+On |
| 2 | GSE123813 | BCC | scRNA | Anti-PD-1 mono | 11+4 | ~79K | Yes | Pre+Post |
| 3 | GSE176021 | NSCLC | scRNA | Nivolumab mono | ~16 | ~816K | MPR/NMPR | Surgical |
| **4** | **GSE301741** | **HNSCC** | **scRNA+TCR** | **Pembro mono** | **16** | **137K** | **pTR** | **Pre+Post** |

### 推荐下载优先级

1. **GSE301741 Seurat RDS** (12.4 Gb): 包含全部137k细胞+metadata+response labels
   - URL: `https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE301741&format=file&file=GSE301741%5FSeurat%5FObject%5FQCpass%5F137020cells%5FwithMetaData%2Erds`

2. **GSE286827 monotherapy subset** (条件性): 使用 Neoadj_type=="D" 筛选13例monotherapy patients
   - 下载7个countdata RDS + 7个metadata RDS
   - 或下载RAW tar (2.5 Gb)

3. **GSE274975 bulk counts** (如需要bulk anchor):
   - URL: `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE274nnn/GSE274975/suppl/GSE274975_raw_counts.tsv.gz`
   - 需配合Table 1.xlsx合并临床数据

---

## 关键经验教训验证

| # | 经验教训 | 本次验证 |
|---|---------|---------|
| 1 | Response labels可能隐藏在sample names中 | **CONFIRMED** - GSE288199使用 HN{patient}_{pre/post}_T 命名 |
| 2 | Zenodo可能是processed data的唯一来源 | CONFIRMED - Zenodo 10407125/14511579有RDS文件 |
| 3 | 论文Supplementary PDF可能是唯一response label来源 | CONFIRMED - GSE274975需Table 1.xlsx |
| 4 | cancer-cell-specific不等于完整TME | CONFIRMED - GSE301741明确包含malignant cells |
| 5 | TISCH2的ICB数据集列表是重要索引 | CONFIRMED - 命名规律{CancerType}_{GSE}_aPD1 |
| 6 | TIGER提供expression + phenotype一站式下载 | **PARTIALLY** - TIGER不提供scRNA expression直接下载，仅bulk |
| 7 | CellxGene可在线预览确认labels | CONFIRMED - Collection 61e422dd可用 |
| 8 | 不要遗漏pathological response关键词 | **CRITICAL** - GSE301741靠pTR捕获，是本次最大发现 |
| 9 | 新发表论文可能是最新数据源 | CONFIRMED - GSE301741来自2026 Cell Rep Med |
| 10 | 检查combination therapy的精确比例 | CONFIRMED - GSE286827 55% dual ICI需排除 |
