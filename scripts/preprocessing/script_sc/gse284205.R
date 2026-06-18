id <- 'gse284205'

library(Seurat)
library(anndataR)
library(Matrix)
library(dplyr)
library(stringr)
library(qs)

gene_anno_list <- readRDS('~/006/data/ref/gene_anno_list.rds')
alias_map <- readRDS('~/006/data/ref/alias_map.rds')

detect_id_type <- function(ids) {
  sapply(ids, function(x) {
    if (grepl('^ENS[A-Z]*G[0-9]+$', x)) return('ensembl')
    if (grepl('^[0-9]+$', x)) return('entrez')
    'symbol'
  }, USE.NAMES = FALSE)
}

convert_to_symbol <- function(ids, id_type, gene_anno_list, alias_map) {
  gene_v19 <- gene_anno_list$v19
  gene_v38 <- gene_anno_list$v38
  out <- character(length(ids))

  ens_idx <- which(id_type == 'ensembl')
  ent_idx <- which(id_type == 'entrez')
  sym_idx <- which(id_type == 'symbol')

  if (length(ens_idx) > 0) {
    ens_ids <- ids[ens_idx]
    map38 <- gene_v38$gene_name[match(ens_ids, gene_v38$gene_id)]
    miss <- is.na(map38)
    if (any(miss)) {
      map19 <- gene_v19$gene_name[match(ens_ids[miss], gene_v19$gene_id)]
      map38[miss] <- map19
    }
    out[ens_idx] <- map38
  }

  if (length(ent_idx) > 0) out[ent_idx] <- NA_character_

  if (length(sym_idx) > 0) {
    sym_ids <- ids[sym_idx]
    in_v38 <- sym_ids %in% gene_v38$gene_name
    in_v19 <- sym_ids %in% gene_v19$gene_name
    alias_in <- sym_ids %in% names(alias_map)
    sym_res <- character(length(sym_ids))
    sym_res[in_v38] <- sym_ids[in_v38]
    sym_res[!in_v38 & in_v19] <- sym_ids[!in_v38 & in_v19]
    sym_res[!in_v38 & !in_v19 & alias_in] <- alias_map[sym_ids[!in_v38 & !in_v19 & alias_in]]
    sym_res[!(in_v38 | in_v19 | alias_in)] <- sym_ids[!(in_v38 | in_v19 | alias_in)]
    out[sym_idx] <- sym_res
  }

  out
}

get_gene_type <- function(symbols, gene_anno_list) {
  gene_v19 <- gene_anno_list$v19
  gene_v38 <- gene_anno_list$v38
  idx38 <- match(symbols, gene_v38$gene_name)
  type38 <- gene_v38$gene_type[idx38]
  miss <- is.na(type38)
  idx19 <- match(symbols[miss], gene_v19$gene_name)
  type38[miss] <- gene_v19$gene_type[idx19]
  type38
}

gene_symbol_align <- function(input_ids, gene_anno_list, alias_map) {
  id_types <- detect_id_type(input_ids)
  aligned_symbols <- convert_to_symbol(input_ids, id_types, gene_anno_list, alias_map)
  gene_types <- get_gene_type(aligned_symbols, gene_anno_list)
  data.frame(
    original_id = input_ids,
    aligned_symbol = aligned_symbols,
    gene_type = gene_types,
    conversion_status = ifelse(is.na(aligned_symbols), 'failed', 'success'),
    stringsAsFactors = FALSE
  )
}

input_dir <- '~/006/data/combo/GSE284205'
out_dir <- '~/006/data/processed/srt/raw'
tmp_dir <- '~/006/tmp'

h5_files <- sort(Sys.glob(file.path(input_dir, 'GSM*_filtered_feature_bc_matrix.h5')))
if (length(h5_files) == 0) stop('No .h5 files found in input_dir')

sample_tokens <- gsub('_filtered_feature_bc_matrix.h5$', '', basename(h5_files))

sample_meta <- data.frame(sample_token = sample_tokens, stringsAsFactors = FALSE) %>%
  mutate(
    gsm_id = str_extract(sample_token, '^GSM[0-9]+'),
    patient_id = str_match(sample_token, '^GSM[0-9]+_liver_biopsy_(P[0-9]+)')[, 2],
    replicate = str_match(sample_token, '_rep([0-9]+)$')[, 2],
    timepoint = ifelse(str_detect(sample_token, '_rep1$') & !str_detect(sample_token, '_F1_'), 'Pre', 'Post')
  ) %>%
  mutate(
    sample_id = paste0(id, '_', gsm_id),
    donor_id = paste0(id, '_', patient_id),
    disease = 'icc',
    organ = 'liver',
    tissue = 'tumor',
    lib = '10x_5_v2',
    target = 'VEGFA;CTLA4;PDL1',
    response = 'none',
    resp_standard = NA_character_,
    treatment = ifelse(timepoint == 'Pre', 'Untreated', 'Immunotherapy'),
    sex = NA_character_,
    age = NA_real_,
    stage = NA_character_,
    anno_orig = NA_character_
  )

if (any(is.na(sample_meta$gsm_id) | is.na(sample_meta$patient_id) | is.na(sample_meta$timepoint))) {
  stop('Failed to parse sample metadata from filenames')
}

mat_list <- vector('list', length(h5_files))
names(mat_list) <- sample_tokens
feature_ref <- NULL

for (i in seq_along(h5_files)) {
  tok <- sample_tokens[i]
  x <- Read10X_h5(h5_files[i])
  cts <- if (is.list(x)) {
    if ('Gene Expression' %in% names(x)) x[['Gene Expression']] else x[[1]]
  } else {
    x
  }

  if (is.null(rownames(cts)) || is.null(colnames(cts))) stop('Read10X_h5 returned matrix without dimnames: ', tok)

  if (is.null(feature_ref)) {
    feature_ref <- rownames(cts)
  } else if (!identical(feature_ref, rownames(cts))) {
    stop('Feature set/order differs across samples at: ', tok)
  }

  colnames(cts) <- paste0(id, '__', tok, '__', gsub('\\.', '_', colnames(cts)))
  mat_list[[tok]] <- cts
}

cts_all <- do.call(cbind, mat_list)
all_cells <- unlist(lapply(mat_list, colnames), use.names = FALSE)
dimnames(cts_all) <- list(feature_ref, all_cells)

if (is.null(rownames(cts_all)) || is.null(colnames(cts_all))) stop('Combined matrix dimnames lost')

genes <- rownames(cts_all)
res <- gene_symbol_align(genes, gene_anno_list, alias_map)
res_type_flt <- subset(res, gene_type %in% c('protein_coding', 'IG_C_gene', 'TR_C_gene'))
if (nrow(res_type_flt) == 0) stop('No genes left after biotype filter')

dup_gene <- res_type_flt[duplicated(res_type_flt$aligned_symbol), 'aligned_symbol']
res_undup <- subset(res_type_flt, !(aligned_symbol %in% dup_gene))
res_dup <- subset(res_type_flt, aligned_symbol %in% dup_gene)

if (nrow(res_dup) > 0) {
  res_dup$meanCts <- cts_all[match(res_dup$original_id, genes), ] %>% rowMeans()
  res_dup <- res_dup %>%
    group_by(aligned_symbol) %>%
    slice_max(meanCts, n = 1, with_ties = FALSE) %>%
    select(-meanCts) %>%
    as.data.frame()
  res_dup <- subset(res_dup, !(original_id %in% c('SPANXB1', 'XAGE2')))
}

res_final <- rbind(res_dup, res_undup)
sel_idx <- match(res_final$original_id, genes)
cts_alignGene <- cts_all[sel_idx, , drop = FALSE]
cts_alignGene <- Matrix::Matrix(cts_alignGene, sparse = TRUE)
dimnames(cts_alignGene) <- list(res_final$aligned_symbol, colnames(cts_all))

if (any(is.na(rownames(cts_alignGene))) || any(rownames(cts_alignGene) == '')) stop('Aligned gene names contain NA/empty')
if (any(duplicated(rownames(cts_alignGene)))) stop('Aligned gene names duplicated after dedup')

cell_meta <- data.frame(cell_id = colnames(cts_alignGene), stringsAsFactors = FALSE) %>%
  mutate(sample_token = str_match(cell_id, paste0('^', id, '__(.+?)__'))[, 2])

meta <- cell_meta %>%
  left_join(sample_meta, by = 'sample_token')

if (any(is.na(meta$sample_id))) stop('Failed metadata join for some cells')

selInfo <- c(
  'cell_id', 'sample_id', 'donor_id', 'disease', 'organ', 'tissue', 'lib',
  'target', 'timepoint', 'response', 'resp_standard', 'sex', 'age', 'stage', 'treatment', 'anno_orig'
)
meta_final <- meta[, selInfo]
rownames(meta_final) <- meta_final$cell_id

if (!identical(colnames(cts_alignGene), rownames(meta_final))) stop('Metadata rownames mismatch count columns')

srt <- CreateSeuratObject(
  counts = cts_alignGene,
  meta.data = meta_final,
  min.cells = 0,
  min.features = 0,
  project = id
)

qs_path <- file.path(out_dir, paste0(id, '.qs'))
h5ad_path <- file.path(out_dir, paste0(id, '.h5ad'))

qs::qsave(srt, qs_path)
write_h5ad(srt, h5ad_path)

sample_counts <- as.data.frame(table(meta_final$sample_id), stringsAsFactors = FALSE)
colnames(sample_counts) <- c('sample_id', 'n_cells')
sample_counts <- sample_counts[order(sample_counts$sample_id), ]

report_path <- file.path(tmp_dir, 'gse284205_alignment_report.md')
lines <- c(
  '# gse284205 Alignment Report',
  '',
  paste0('- Generated at: ', format(Sys.time(), '%Y-%m-%d %H:%M:%S')),
  paste0('- Input dir: ', input_dir),
  paste0('- Output qs: ', qs_path),
  paste0('- Output h5ad: ', h5ad_path),
  '',
  '## 1) Dataset scale',
  paste0('- Samples detected: ', length(sample_tokens)),
  paste0('- Cells (raw): ', ncol(cts_all)),
  paste0('- Genes (raw): ', nrow(cts_all)),
  paste0('- Genes (aligned+filtered): ', nrow(cts_alignGene)),
  paste0('- Cells (final): ', ncol(cts_alignGene)),
  '',
  '## 2) Mapping policy used (user confirmed)',
  '- disease = icc',
  '- organ = liver, tissue = tumor',
  '- target = VEGFA;CTLA4;PDL1',
  '- timepoint = Pre/Post from filename / rep pattern',
  '- treatment = Untreated (Pre), Immunotherapy (Post)',
  '- response = none',
  '- resp_standard = NA',
  "- lib = 10x_5_v2 (from GEO soft: '10X Genomics 5’ v2 Single Cell User Guide' and 5' mRNA library text)",
  '',
  '## 3) Metadata completeness checks',
  paste0('- Missing sample_id: ', sum(is.na(meta_final$sample_id))),
  paste0('- Missing donor_id: ', sum(is.na(meta_final$donor_id))),
  paste0('- Missing disease: ', sum(is.na(meta_final$disease))),
  paste0('- Missing organ: ', sum(is.na(meta_final$organ))),
  paste0('- Missing tissue: ', sum(is.na(meta_final$tissue))),
  paste0('- Missing lib: ', sum(is.na(meta_final$lib))),
  paste0('- Missing target: ', sum(is.na(meta_final$target))),
  paste0('- Missing timepoint: ', sum(is.na(meta_final$timepoint))),
  paste0('- Missing response: ', sum(is.na(meta_final$response))),
  paste0('- Missing resp_standard: ', sum(is.na(meta_final$resp_standard))),
  '',
  '## 4) Controlled-value snapshots',
  paste0('- organ: ', paste(sort(unique(meta_final$organ)), collapse = ', ')),
  paste0('- tissue: ', paste(sort(unique(meta_final$tissue)), collapse = ', ')),
  paste0('- lib: ', paste(sort(unique(meta_final$lib)), collapse = ', ')),
  paste0('- target: ', paste(sort(unique(meta_final$target)), collapse = ', ')),
  paste0('- timepoint: ', paste(sort(unique(meta_final$timepoint)), collapse = ', ')),
  paste0('- response: ', paste(sort(unique(meta_final$response)), collapse = ', ')),
  '',
  '## 5) Evidence source',
  '- GEO soft: https://ftp.ncbi.nlm.nih.gov/geo/series/GSE284nnn/GSE284205/soft/GSE284205_family.soft.gz',
  '- Series summary: liver cancer tumor biopsy samples (pre- and post-treatment) from three patients',
  '- Sample-level titles: P1699/P1700/P1779 with Pre/Post and replicate labels',
  '',
  '## 6) Per-sample cell counts',
  '',
  '| sample_id | n_cells |',
  '|---|---:|'
)

tab_lines <- apply(sample_counts, 1, function(x) paste0('| ', x[[1]], ' | ', x[[2]], ' |'))
lines <- c(lines, tab_lines)

writeLines(lines, report_path)

cat('Done\n')
cat('Report: ', report_path, '\n')
