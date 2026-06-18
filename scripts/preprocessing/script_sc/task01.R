id <- 'task01'

library(Seurat)
library(anndataR)
library(Matrix)
library(dplyr)
library(stringr)
library(qs)

# Internal gene alignment utilities (embedded to avoid modifying shared script)
gene_anno_list <- readRDS('~/006/data/ref/gene_anno_list.rds')

read_hgnc_alias_local <- function(hgnc_file) {
  hgnc <- read.delim(hgnc_file, stringsAsFactors = FALSE, quote = "", check.names = FALSE)
  split_alias <- function(symbols, raw_alias) {
    keep <- !is.na(raw_alias) & raw_alias != ''
    if (!any(keep)) {
      return(data.frame(symbol = character(), alias = character(), stringsAsFactors = FALSE))
    }
    raw_alias <- gsub('"', '', raw_alias[keep])
    symbol <- symbols[keep]
    pieces <- strsplit(raw_alias, '\\|', fixed = FALSE)
    alias <- unlist(pieces, use.names = FALSE)
    symbol_rep <- rep(symbol, lengths(pieces))
    data.frame(symbol = symbol_rep, alias = trimws(alias), stringsAsFactors = FALSE)
  }
  alias1 <- split_alias(hgnc$symbol, hgnc$alias_symbol)
  alias2 <- split_alias(hgnc$symbol, hgnc$prev_symbol)
  alias_df <- unique(rbind(alias1, alias2))
  alias_df <- alias_df[alias_df$alias != '' & !is.na(alias_df$alias), , drop = FALSE]
  alias_map <- setNames(alias_df$symbol, alias_df$alias)
  return(alias_map)
}

alias_map <- read_hgnc_alias_local('~/006/data/ref/hgnc_complete_set.txt')

detect_id_type <- function(ids) {
  sapply(ids, function(id) {
    if (grepl("^ENS[A-Z]*G[0-9]+$", id)) return("ensembl")
    if (grepl("^[0-9]+$", id)) return("entrez")
    return("symbol")
  }, USE.NAMES = FALSE)
}

convert_to_symbol <- function(ids, id_type, gene_anno_list, alias_map) {
  gene_v19 <- gene_anno_list$v19
  gene_v38 <- gene_anno_list$v38
  symbol_out <- character(length(ids))
  ens_idx <- which(id_type == "ensembl")
  ent_idx <- which(id_type == "entrez")
  sym_idx <- which(id_type == "symbol")
  if(length(ens_idx) > 0) {
    ens_ids <- ids[ens_idx]
    map38 <- gene_v38$gene_name[match(ens_ids, gene_v38$gene_id)]
    miss38 <- is.na(map38)
    if(any(miss38)) {
      map19 <- gene_v19$gene_name[match(ens_ids[miss38], gene_v19$gene_id)]
      map38[miss38] <- map19
    }
    symbol_out[ens_idx] <- map38
  }
  if(length(ent_idx) > 0) symbol_out[ent_idx] <- NA_character_
  if(length(sym_idx) > 0) {
    sym_ids <- ids[sym_idx]
    in_v38 <- sym_ids %in% gene_v38$gene_name
    in_v19 <- sym_ids %in% gene_v19$gene_name
    alias_in <- sym_ids %in% names(alias_map)
    sym_res <- character(length(sym_ids))
    sym_res[in_v38] <- sym_ids[in_v38]
    sym_res[!in_v38 & in_v19] <- sym_ids[!in_v38 & in_v19]
    sym_res[!in_v38 & !in_v19 & alias_in] <- alias_map[sym_ids[!in_v38 & !in_v19 & alias_in]]
    unmatched <- !(in_v38 | in_v19 | alias_in)
    sym_res[unmatched] <- sym_ids[unmatched]
    symbol_out[sym_idx] <- sym_res
  }
  return(symbol_out)
}

get_gene_type <- function(symbols, gene_anno_list) {
  gene_v19 <- gene_anno_list$v19
  gene_v38 <- gene_anno_list$v38
  idx38 <- match(symbols, gene_v38$gene_name)
  type38 <- gene_v38$gene_type[idx38]
  missing <- is.na(type38)
  idx19 <- match(symbols[missing], gene_v19$gene_name)
  type19 <- gene_v19$gene_type[idx19]
  type38[missing] <- type19
  return(type38)
}

gene_symbol_align <- function(input_ids, gene_anno_list, alias_map) {
  message("Detecting input ID types...")
  id_types <- detect_id_type(input_ids)
  message("Mapping IDs to approved gene symbols...")
  aligned_symbols <- convert_to_symbol(input_ids, id_types, gene_anno_list, alias_map)
  message("Getting gene biotype annotations...")
  gene_types <- get_gene_type(aligned_symbols, gene_anno_list)
  data.frame(
    original_id = input_ids,
    aligned_symbol = aligned_symbols,
    gene_type = gene_types,
    conversion_status = ifelse(is.na(aligned_symbols), "failed", "success"),
    stringsAsFactors = FALSE
  )
}

input_dir <- '~/004.1/zz.CELLRANGER/TASK01'
meta_csv <- '~/004.1/zz.META/task01_raw.csv'
response_csv <- '~/004.1/zz.META/task_response.csv'
out_dir <- '~/006/data/processed/srt/raw'
tmp_dir <- '~/006/tmp'

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)

read_10x_counts <- function(path) {
  x <- Read10X(path)
  if (is.list(x)) {
    if ('Gene Expression' %in% names(x)) return(x[['Gene Expression']])
    return(x[[1]])
  }
  x
}

standardize_task_meta <- function(meta_csv, response_csv, proj_id, target_value, treatment_value) {
  meta <- read.csv(meta_csv, stringsAsFactors = FALSE)
  meta <- subset(meta, rawID != '')
  meta$pt_label <- sub('_T[0-9]+$', '', meta$rawID)
  
  resp <- read.csv(response_csv, stringsAsFactors = FALSE)
  resp <- subset(resp, Proj == 'TASK01')
  resp <- resp[, c('PtName', 'Response')]
  names(resp) <- c('pt_label', 'response')
  
  meta <- dplyr::left_join(meta, resp, by = 'pt_label')
  if (any(is.na(meta$response))) stop('Missing response mapping for task01')
  
  meta %>%
    mutate(
      sample_token = dirname,
      sample_id = paste0(proj_id, '_', rawID),
      donor_id = paste0(proj_id, '_', pt_label),
      disease = 'hcc',
      organ = 'liver',
      tissue = 'tumor',
      lib = '10x_5',
      target = target_value,
      timepoint = str_extract(rawID, 'T[0-9]+$'),
      resp_standard = 'mRECIST',
      treatment = treatment_value,
      sex = NA_character_,
      age = NA_real_,
      stage = NA_character_,
      anno_orig = NA_character_
    )
}

sample_meta <- standardize_task_meta(meta_csv, response_csv, id, 'PD1', 'PD1')
sample_meta$timepoint[is.na(sample_meta$timepoint)] <- 'NA'

if (any(is.na(sample_meta$sample_id) | is.na(sample_meta$donor_id) | is.na(sample_meta$response))) {
  stop('Task01 metadata is incomplete after standardization')
}

sample_tokens <- sample_meta$sample_token
sample_meta$sample_token <- sample_tokens

# Read all samples with union of genes
mat_list <- vector('list', length(sample_tokens))
names(mat_list) <- sample_tokens
all_genes_list <- list()

for (i in seq_along(sample_tokens)) {
  tok <- sample_tokens[i]
  dir <- sprintf('%s/%s/', input_dir, tok)
  if (dir.exists(file.path(dir, 'filtered_feature_bc_matrix'))) {
    cts <- read_10x_counts(file.path(dir, 'filtered_feature_bc_matrix'))
  } else {
    cts <- read_10x_counts(file.path(dir, 'outs', 'filtered_feature_bc_matrix'))
  }
  if (is.null(rownames(cts)) || is.null(colnames(cts))) stop('Missing dimnames for sample: ', tok)
  all_genes_list[[tok]] <- rownames(cts)
  colnames(cts) <- paste0(id, '__', tok, '__', gsub('\\.', '_', colnames(cts)))
  mat_list[[tok]] <- cts
}

# Union alignment: take union of all genes, fill missing with 0
all_genes <- unique(unlist(all_genes_list))
message('Union gene count: ', length(all_genes))

# Create union matrix with all genes
union_mat_list <- lapply(names(mat_list), function(tok) {
  cts <- mat_list[[tok]]
  missing_genes <- setdiff(all_genes, rownames(cts))
  if (length(missing_genes) > 0) {
    zero_mat <- Matrix::Matrix(0, nrow = length(missing_genes), ncol = ncol(cts), sparse = TRUE)
    rownames(zero_mat) <- missing_genes
    colnames(zero_mat) <- colnames(cts)
    cts <- rbind(cts, zero_mat)
  }
  cts <- cts[all_genes, , drop = FALSE]
  cts
})

cts_all <- do.call(cbind, union_mat_list)
message('Combined matrix dimensions: ', nrow(cts_all), ' x ', ncol(cts_all))

# Gene symbol alignment
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

# Build cell metadata
cell_meta <- data.frame(cell_id = colnames(cts_alignGene), stringsAsFactors = FALSE) %>%
  mutate(sample_token = str_match(cell_id, paste0('^', id, '__(.+?)__'))[, 2])

meta <- dplyr::left_join(cell_meta, sample_meta, by = 'sample_token')
if (any(is.na(meta$sample_id))) stop('Failed metadata join for task01')

selInfo <- c('cell_id', 'sample_id', 'donor_id', 'disease', 'organ', 'tissue', 'lib',
             'target', 'timepoint', 'response', 'resp_standard', 'sex', 'age', 'stage', 'treatment', 'anno_orig')
meta_final <- meta[, selInfo]
rownames(meta_final) <- meta_final$cell_id

if (!identical(colnames(cts_alignGene), rownames(meta_final))) stop('Metadata rownames mismatch count columns')

srt <- CreateSeuratObject(counts = cts_alignGene, meta.data = meta_final, min.cells = 0, min.features = 0, project = id)

qs_path <- file.path(out_dir, paste0(id, '.qs'))
h5ad_path <- file.path(out_dir, paste0(id, '.h5ad'))
qs::qsave(srt, qs_path, nthreads = 20)
write_h5ad(srt, h5ad_path)

sample_counts <- as.data.frame(table(meta_final$sample_id), stringsAsFactors = FALSE)
names(sample_counts) <- c('sample_id', 'n_cells')
sample_counts <- sample_counts[order(sample_counts$sample_id), ]

report_path <- file.path(tmp_dir, 'task01_alignment_report.md')
lines <- c(
  '# task01 Alignment Report',
  '',
  paste0('- Generated at: ', format(Sys.time(), '%Y-%m-%d %H:%M:%S')),
  paste0('- Input dir: ', input_dir),
  paste0('- Output qs: ', qs_path),
  paste0('- Output h5ad: ', h5ad_path),
  paste0('- Raw dimensions: ', length(all_genes), ' genes x ', ncol(cts_all), ' cells'),
  paste0('- Aligned dimensions: ', nrow(cts_alignGene), ' genes x ', ncol(cts_alignGene), ' cells'),
  '',
  '## Required-field completeness',
  paste0('- sample_id missing: ', sum(is.na(meta_final$sample_id))),
  paste0('- donor_id missing: ', sum(is.na(meta_final$donor_id))),
  paste0('- disease missing: ', sum(is.na(meta_final$disease))),
  paste0('- organ missing: ', sum(is.na(meta_final$organ))),
  paste0('- tissue missing: ', sum(is.na(meta_final$tissue))),
  paste0('- lib missing: ', sum(is.na(meta_final$lib))),
  '',
  '## Controlled vocabulary snapshot',
  paste0('- organ: ', paste(sort(unique(meta_final$organ)), collapse = ', ')),
  paste0('- tissue: ', paste(sort(unique(meta_final$tissue)), collapse = ', ')),
  paste0('- lib: ', paste(sort(unique(meta_final$lib)), collapse = ', ')),
  paste0('- target: ', paste(sort(unique(meta_final$target)), collapse = ', ')),
  paste0('- timepoint: ', paste(sort(unique(meta_final$timepoint)), collapse = ', ')),
  paste0('- response: ', paste(sort(unique(meta_final$response)), collapse = ', ')),
  paste0('- resp_standard: ', paste(sort(unique(meta_final$resp_standard)), collapse = ', ')),
  '',
  '## Per-sample cell counts',
  paste(capture.output(print(sample_counts, row.names = FALSE)), collapse = '\n'),
  '',
  '## ID integrity',
  paste0('- cell_id unique: ', !anyDuplicated(meta_final$cell_id)),
  paste0('- rownames(meta_final) aligned to count cols: ', identical(colnames(cts_alignGene), rownames(meta_final))),
  paste0('- final gene names unique: ', !anyDuplicated(rownames(cts_alignGene)))
)
writeLines(lines, report_path)

message('Saved task01 QS: ', qs_path)
message('Saved task01 H5AD: ', h5ad_path)
message('Saved task01 report: ', report_path)
