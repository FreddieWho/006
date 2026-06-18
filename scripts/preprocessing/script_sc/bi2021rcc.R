library(Seurat)
library(anndataR) 
source('~/006/scripts/utils/gene_name_align.R')
# Step 1: read cts and data -----------------------------------------------
cts <- Read10X('~/006/data/combo/bi_2021_rcc/',strip.suffix = T)

# Step 2: align gene symbol -----------------------------------------------
res<- gene_symbol_align(rownames(cts),gene_anno_list,alias_map)
res_type_flt <- subset(res,gene_type %in% c('protein_coding','IG_C_gene','TR_C_gene'))

dup_gene <- res_type_flt[duplicated(res_type_flt$aligned_symbol),'aligned_symbol'] 

res_undup <- subset(res_type_flt,!(aligned_symbol %in% dup_gene))
res_dup <- subset(res_type_flt,aligned_symbol %in% dup_gene)
res_dup$meanCts <- cts[res_dup$original_id,] %>% rowMeans()

res_dup %>% 
  group_by(aligned_symbol) %>% 
  slice_max(meanCts) %>% 
  select(-meanCts) -> res_dup

res_final <- rbind(res_dup,res_undup)

cts_alignGene <- cts[res_final$original_id,]


# Step 3: refine cell id --------------------------------------------------
colnames(cts_alignGene) <- paste0('bi2021rcc_',gsub('\\.','_',colnames(cts_alignGene)))


# Step 4: patient info ----------------------------------------------------
# fields: cell_id, sample_id, donor_id,disease,organ,tissue,lib,target,timepoint,response,resp_standard,anno_orig

meta <- read.csv('~/006/data/combo/bi_2021_rcc/Final_SCP_Metadata.txt',
                 sep = '\t')[-1,]
meta$cell_id_raw <- meta$NAME
meta$cell_id <- paste0('bi2021rcc_',gsub('\\.','_',meta$NAME))
meta$sample_id <- meta$donor_id
meta$disease <- 'ccrcc'
meta$organ <-  meta$organ__ontology_label
meta$tissue <- 'tumor'
meta$lib <- '10x_3_v2'
meta$target <- case_when(
  meta$sample_id %in% c('P55','P906') ~ 'PD1',
  meta$sample_id %in% c('P912','P913') ~ 'PD1;TKI',
  meta$sample_id == 'P915' ~ 'PD1;CTLA4',
  T ~ 'none'
)
meta$timepoint <- 'Post'
meta$anno_orig <- meta$FinalCellType

meta$response <- case_when(
  meta$ICB_Response == 'ICB_PR' ~ 'PR',
  meta$ICB_Response == 'ICB_PD' ~ 'PD',
  meta$ICB_Response == 'ICB_SD' ~ 'SD',
  meta$ICB_Response == 'ICB_NE'  ~ 'NE',
  T ~ 'none'
)
meta$resp_standard <- 'RECIST_LIKE'

selInfo <- c('cell_id_raw','cell_id','sample_id','donor_id','disease','organ','tissue','lib','target','timepoint','response','resp_standard','anno_orig')
meta <- meta[,selInfo]
rownames(meta) <- meta$cell_id

meta$sample_id <- paste0('bi2021rcc_',meta$sample_id)
meta$donor_id <- paste0('bi2021rcc_',meta$donor_id)
srt <-  CreateSeuratObject(cts_alignGene,meta.data = meta,min.cells = 0,min.features = 0,project = 'bi2021rcc')
qs::qsave(srt,'~/006/data/processed/srt/raw/bi2021rcc.qs')
write_h5ad(srt, '~/006/data/processed/srt/raw/bi2021rcc.h5ad')

