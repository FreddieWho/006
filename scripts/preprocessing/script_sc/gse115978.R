library(Seurat)
library(anndataR) 
source('~/006/scripts/utils/gene_name_align.R')
# Step 1: read cts and data -----------------------------------------------
cts <- data.table::fread('~/006/data/combo/GSE115978/GSE115978_counts.csv') %>% as.data.frame()

genes <- cts$V1
cts <- cts[,-1]
res<- gene_symbol_align(genes,gene_anno_list,alias_map)
res_type_flt <- subset(res,gene_type %in% c('protein_coding','IG_C_gene','TR_C_gene'))

dup_gene <- res_type_flt[duplicated(res_type_flt$aligned_symbol),'aligned_symbol'] 

res_undup <- subset(res_type_flt,!(aligned_symbol %in% dup_gene))

res_dup <- subset(res_type_flt,aligned_symbol %in% dup_gene)
res_dup$meanCts <- cts[match(res_dup$original_id,genes),] %>% rowMeans()
res_dup %>% 
  group_by(aligned_symbol) %>% 
  slice_max(meanCts) %>% 
  select(-meanCts) -> res_dup
res_dup <- subset(res_dup,!(original_id %in% c('SPANXB2','SPANXF1')))

res_final <- rbind(res_dup,res_undup)

cts_alignGene <- cts[match(res_final$original_id,genes),]
rownames(cts_alignGene) <- res_final$aligned_symbol


# clean meta -------------------------------------------------------------
# fields: cell_id, sample_id, doron_id,disease,organ,tissue,lib,target,timepoint,response,resp_standard,anno_orig
colnames(cts_alignGene) <- paste0('gse115978_',gsub('\\.','_',colnames(cts_alignGene)))

meta <- data.table::fread('~/006/data/combo/GSE115978/GSE115978_cell.annotations.csv') %>% as.data.frame()

meta_cli <- read.csv('~/006/data/combo/GSE115978/gse115978.csv')

meta$sample_id <- meta$samples

rownames(meta) <- paste0('gse115978_',gsub('\\.','_',meta$cells))
meta$cell_id <- rownames(meta)

meta <- merge(meta,meta_cli,by = c('sample_id','Cohort'),all.x = T)

meta$disease <- 'mela'
meta$lib <- 'ss2'
meta$resp_standard <- 'RECIST'
meta$anno_orig <- meta$cell.types

selInfo <- c('cell_id','sample_id','donor_id','disease','organ','tissue','lib','target','response','resp_standard','anno_orig')
selInfo[!(selInfo %in% colnames(meta))]
meta_final <- meta[,selInfo]

meta$sample_id <- paste0('gse115978_',meta$sample_id)
meta$donor_id <- paste0('gse115978_',meta$donor_id)

srt <-  CreateSeuratObject(cts_alignGene,meta.data = meta_final,min.cells = 0,min.features = 0,project = 'gse115978')
qs::qsave(srt,'~/006/data/processed/srt/raw/gse115978.qs')
write_h5ad(srt, '~/006/data/processed/srt/raw/gse115978.h5ad')
