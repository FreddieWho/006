id = 'gse151530'

library(Seurat)
library(anndataR) 
source('~/006/scripts/utils/gene_name_align.R')
# Step 1: read cts and data -----------------------------------------------

cts <- Read10X('~/006/data/combo/GSE151530/')
meta <- read.csv('~/006/data/combo/GSE151530/GSE151530_Info.txt',
                 sep = '\t')
meta_cli <- read.csv('~/006/data/combo/GSE151530/gse151530.csv')
# clean meta --------------------------------------------------------------
# fields: cell_id, sample_id, donor_id,sex,age,disease,organ,tissue,lib,target,timepoint,
#         surv,time_type,surv_type,surv_status, response,resp_standard,anno_orig
meta$cell_id <- meta$Cell
meta$sample_id <- meta$Sample
meta$donor_id <- gsub('[a-z]$','',meta$Sample)
meta <- left_join(meta,meta_cli)

meta$cell_id <- paste0(id,'_',meta$cell_id)
meta$sample_id <- paste0(id,'_',meta$sample_id)
meta$donor_id <- paste0(id,'_',meta$donor_id)

meta$organ <- 'liver'
meta$tissue <- 'tumor'

meta$lib <- '10x_3_v2'

rownames(meta) <- meta$cell_id

selInfo  <- c('cell_id','sample_id','donor_id','sex','age','disease','disease',
              'organ','tissue','lib','target','timepoint','stage')
meta_final <- meta[,selInfo]

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

srt <-  CreateSeuratObject(cts_alignGene,meta.data = meta_final,min.cells = 0,min.features = 0,project = id)
qs::qsave(srt,'~/006/data/processed/srt/raw/gse151530.qs')
write_h5ad(srt, '~/006/data/processed/srt/raw/gse151530.h5ad')
