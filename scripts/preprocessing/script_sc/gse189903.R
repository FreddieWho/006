id = 'gse189903'

library(Seurat)
library(anndataR) 
source('~/006/scripts/utils/gene_name_align.R')
# Step 1: read cts and data -----------------------------------------------

cts <- Read10X('~/006/data/combo/GSE189903/')
meta <- read.csv('~/006/data/combo/GSE189903/GSE189903_Info.txt.gz',sep = '\t')

# gene name align ---------------------------------------------------------
colnames(cts) <- paste0(id,'_',colnames(cts))
genes <- rownames(cts)
res<- gene_symbol_align(genes,gene_anno_list,alias_map)
res_type_flt <- subset(res,gene_type %in% c('protein_coding','IG_C_gene','TR_C_gene'))

dup_gene <- res_type_flt[duplicated(res_type_flt$aligned_symbol),'aligned_symbol'] 

res_undup <- subset(res_type_flt,!(aligned_symbol %in% dup_gene))

res_dup <- subset(res_type_flt,aligned_symbol %in% dup_gene)
res_dup$meanCts <- cts[match(res_dup$original_id,genes),] %>% rowMeans()
res_dup %>% 
  group_by(aligned_symbol) %>% 
  slice_max(meanCts) %>% 
  select(-meanCts)-> res_dup
res_dup <- subset(res_dup,!(original_id %in% c('SPANXB1','XAGE2')))

res_final <- rbind(res_dup,res_undup)

cts_alignGene <- cts[match(res_final$original_id,genes),]
rownames(cts_alignGene) <- res_final$aligned_symbol
dim(cts_alignGene)

# clean meta --------------------------------------------------------------
# fields: cell_id, sample_id, donor_id,cell_sorting, sex, age,disease,organ,tissue,
#         lib,target,timepoint, surv,time_type,surv_type,surv_status, 
#         response,resp_standard,anno_orig
meta$cell_id <- paste0(id,'_',meta$Cell)
meta$sample_id <- paste0(id,'_',meta$Sample)
meta$donor_id <- paste0(id,'_',substr(meta$Sample,0,1))

meta$cell_sorting <- 'whole'

meta$disease <- ifelse(substr(meta$Sample,2,2) == 'C','icc','hcc')
meta$organ <- 'liver'
meta$tissue <- case_when(
  substr(meta$Sample,3,3) == 'B' ~ 'border',
  substr(meta$Sample,3,3) == 'T' ~ 'tumor',
  substr(meta$Sample,3,3) == 'N' ~ 'adj'
)

meta$lib <- '10x_5_v1'
meta$timepoint <- 'Pre'
meta$anno_orig <- meta$Type


selInfo <- c('cell_id', 'sample_id', 'donor_id','cell_sorting','disease','organ','tissue','lib','timepoint',
             'anno_orig')
selInfo[!(selInfo %in% colnames(meta))]
meta_final <- meta[,selInfo]
rownames(meta_final) <- meta_final$cell_id

# final -------------------------------------------------------------------

srt <-  CreateSeuratObject(cts_alignGene,meta.data = meta_final,min.cells = 0,min.features = 0,project = id)
qs::qsave(srt,'~/006/data/processed/srt/raw/gse189903.qs')
write_h5ad(srt, '~/006/data/processed/srt/raw/gse189903.h5ad')
