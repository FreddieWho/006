id = 'gse169246'
library(Seurat)
library(anndataR) 
source('~/006/scripts/utils/gene_name_align.R')
# Step 1: read cts and data -----------------------------------------------

cts <- Matrix::readMM('~/006/data/combo/GSE169246/matrix.mtx')
barcodes <- read.csv('~/006/data/combo/GSE169246/barcodes.tsv',header = F)[,1,drop = T]
features <- read.csv('~/006/data/combo/GSE169246/features.tsv',header = F)[,1,drop = T]
rownames(cts) <- features
colnames(cts) <- barcodes

# gene name align ---------------------------------------------------------
colnames(cts) <- paste0(id,'_',gsub('\\.','_',colnames(cts)))
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

# clean meta --------------------------------------------------------------
# fields: cell_id, sample_id, donor_id,cell_sorting, sex, age,disease,organ,tissue,
#         lib,target,timepoint, surv,time_type,surv_type,surv_status, 
#         response,resp_standard,anno_orig
meta_cli <- read.csv('~/006/data/combo/GSE169246/gse169246.csv')
meta_cli$donor_id <- paste0(id,'_',meta_cli$donor_id)

meta <- data.frame(cell_id = colnames(cts),
                   row.names = colnames(cts))
meta$sample_id <- paste0(id,'_',gsub('.*\\.','',meta$cell_id))
meta$donor_id <- paste0(id,'_',gsub('.*_','',gsub('_[tb]','',meta$sample_id)))

meta <- left_join(meta,meta_cli)

meta$cell_sorting <- 'CD45p'
meta$sex <- 'F'
meta$disease <- 'tnbc'
meta$tissue <- ifelse(stringr::str_detect(meta$sample_id,'_b^'),'blood','tumor')

meta$lib <- '10x_5_v2'
meta$target <- meta$targe

meta$timepoint <- gsub('_.*','',gsub('.*\\.','',meta$cell_id))

meta$resp_standard <- 'RECIST'

selInfo <- c('cell_id', 'sample_id', 'donor_id','cell_sorting','sex', 'age','disease','organ','tissue','lib','target','timepoint',
             'surv','time_type','surv_type','surv_status', 'response','resp_standard')
selInfo[!(selInfo %in% colnames(meta))]
meta_final <- meta[,selInfo]
rownames(meta_final) <- meta_final$cell_id

# final -------------------------------------------------------------------

srt <-  CreateSeuratObject(cts_alignGene,meta.data = meta_final,min.cells = 0,min.features = 0,project = id)
qs::qsave(srt,'~/006/data/processed/srt/raw/gse169246.qs')
write_h5ad(srt, '~/006/data/processed/srt/raw/gse169246.h5ad')
