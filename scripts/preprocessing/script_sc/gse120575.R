library(Seurat)
library(anndataR) 
source('~/006/scripts/utils/gene_name_align.R')
library(data.table)
# Step 1: read cts and data -----------------------------------------------

header <- fread('~/006/data/combo/GSE120575/GSE120575_Sade_Feldman_melanoma_single_cells_TPM_GEO.txt', 
                data.table = FALSE,sep = '\t',nrow = 2)[,-1] %>% as.data.frame()
cts <-fread('~/006/data/combo/GSE120575/GSE120575_Sade_Feldman_melanoma_single_cells_TPM_GEO.txt',
            skip = 2,header = FALSE) %>% as.data.frame()

genes <- cts$V1
cts <- cts[,-1]
cts <- cts[,-16292]


# gene name align ---------------------------------------------------------

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

colnames(cts_alignGene) <- unlist(header[1,])

colnames(cts_alignGene) <- paste0('gse120575_',gsub('\\.','_',colnames(cts_alignGene)))

# clean meta --------------------------------------------------------------
# fields: cell_id, sample_id, donor_id,disease,organ,tissue,lib,target,timepoint,
#         surv,time_type,surv_type,surv_status, response,resp_standard,anno_orig
meta <- fread('~/006/data/combo/GSE120575/GSE120575_patient_ID_single_cells.txt',skip = 19,header = T) %>% 
  as.data.frame() %>% subset(`characteristics: therapy` != '')
meta_cli <- read.csv('~/006/data/combo/GSE120575/gse120575.csv')

meta$cell_id <- paste0('gse120575_',gsub('\\.','_',meta$title))
meta$sample_id <- meta$`characteristics: patinet ID (Pre=baseline; Post= on treatment)`
meta$donor_id <- gsub('_.*','',gsub('.*_','',meta$`characteristics: patinet ID (Pre=baseline; Post= on treatment)`))
meta <- left_join(meta,meta_cli)
meta$disease <- 'mela'
meta$organ <- 'unassign'
meta$tissue <- 'tumor'
meta$lib <- 'ss2'
meta$timepoint <- gsub('_.*','',meta$`characteristics: patinet ID (Pre=baseline; Post= on treatment)`)
meta$response <- meta$`characteristics: response`
meta$resp_standard <- 'RECIST'

meta$cell_sorting <- 'CD45p'

selInfo <- c('cell_id', 'sample_id', 'donor_id','cell_sorting','disease','organ','tissue','lib','target','timepoint',
             'surv','time_type','surv_type','surv_status', 'response','resp_standard')
selInfo[!(selInfo %in% colnames(meta))]
meta_final <- meta[,selInfo]
rownames(meta_final) <- meta_final$cell_id

meta_final$sample_id <- paste0('gse120575_',meta_final$sample_id)
meta_final$donor_id <- paste0('gse120575_',meta_final$donor_id)

srt <-  CreateSeuratObject(cts_alignGene,meta.data = meta_final,min.cells = 0,min.features = 0,project = 'gse120575')
qs::qsave(srt,'~/006/data/processed/srt/raw/gse120575.qs')
write_h5ad(srt, '~/006/data/processed/srt/raw/gse120575.h5ad')
