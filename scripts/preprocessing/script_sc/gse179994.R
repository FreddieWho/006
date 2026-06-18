id <- 'gse179994'

library(Seurat)
library(anndataR)
library(Matrix)
library(dplyr)
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

align_matrix <- function(cts, prefix) {
  genes <- rownames(cts)
  res <- gene_symbol_align(genes, gene_anno_list, alias_map)
  res_type_flt <- subset(res, gene_type %in% c('protein_coding', 'IG_C_gene', 'TR_C_gene'))
  if (nrow(res_type_flt) == 0) stop('No genes left after biotype filter for ', prefix)

  dup_gene <- res_type_flt[duplicated(res_type_flt$aligned_symbol), 'aligned_symbol']
  res_undup <- subset(res_type_flt, !(aligned_symbol %in% dup_gene))
  res_dup <- subset(res_type_flt, aligned_symbol %in% dup_gene)

  if (nrow(res_dup) > 0) {
    res_dup$meanCts <- cts[match(res_dup$original_id, genes), ] %>% rowMeans()
    res_dup <- res_dup %>%
      group_by(aligned_symbol) %>%
      slice_max(meanCts, n = 1, with_ties = FALSE) %>%
      select(-meanCts) %>%
      as.data.frame()
  }

  res_final <- rbind(res_dup, res_undup)
  sel_idx <- match(res_final$original_id, genes)
  cts_align <- cts[sel_idx, , drop = FALSE]
  cts_align <- Matrix::Matrix(cts_align, sparse = TRUE)
  dimnames(cts_align) <- list(res_final$aligned_symbol, colnames(cts))
  if (any(duplicated(rownames(cts_align)))) stop('Duplicated aligned genes in ', prefix)
  cts_align
}

input_counts <- '~/006/data/imm/GSE179994/GSE179994_all.Tcell.rawCounts.rds.gz'
input_meta <- '~/006/data/imm/GSE179994/GSE179994_Tcell.metadata.tsv.gz'
input_tcr <- '~/006/data/imm/GSE179994/GSE179994_all.scTCR.tsv.gz'
out_dir <- '~/006/data/processed/srt/raw'
tmp_dir <- '~/006/tmp'

if (!file.exists(path.expand(input_counts))) stop('Missing counts: ', input_counts)
if (!file.exists(path.expand(input_meta))) stop('Missing metadata: ', input_meta)
if (!file.exists(path.expand(input_tcr))) stop('Missing TCR file: ', input_tcr)

message('Loading counts and metadata...')
counts_con <- gzfile(path.expand(input_counts), 'rb')
counts <- readRDS(counts_con)
close(counts_con)
meta <- read.delim(path.expand(input_meta), sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)

if (!inherits(counts, 'dgCMatrix') && !inherits(counts, 'dgTMatrix') && !inherits(counts, 'matrix')) {
  counts <- as(counts, 'CsparseMatrix')
}

if (!identical(colnames(counts), meta$cellid)) {
  stop('Counts columns do not match metadata cellid order')
}

message('Aligning genes...')
counts <- align_matrix(counts, id)

meta$cell_id_raw <- meta$cellid
meta$cell_id <- paste0(id, '_', meta$cellid)
meta$sample_id_raw <- meta$sample
meta$sample_id <- paste0(id, '_', meta$sample)
meta$donor_id <- paste0(id, '_', meta$patient)
meta$disease <- 'nsclc'
meta$organ <- 'lung'
meta$tissue <- 'tumor'
meta$lib <- '10x_5'
meta$target <- 'PD1'
meta$timepoint <- ifelse(grepl('post', meta$sample), 'post', 'pre')
meta$response <- NA_character_
meta$resp_standard <- 'RECIST'
meta$cell_sorting <- 'Tcell'
meta$source_tag <- 'rawCounts_rds'
meta$anno_orig <- ifelse(!is.na(meta$cluster) & meta$cluster != 'NA' & nzchar(meta$cluster), meta$cluster, meta$celltype)

selInfo <- c(
  'cell_id_raw', 'cell_id', 'sample_id_raw', 'sample_id', 'donor_id', 'cell_sorting',
  'disease', 'organ', 'tissue', 'lib', 'target', 'timepoint', 'response', 'resp_standard',
  'anno_orig', 'patient', 'sample', 'celltype', 'cluster', 'source_tag'
)
selInfo <- selInfo[selInfo %in% colnames(meta)]
meta_final <- meta[, selInfo, drop = FALSE]
rownames(meta_final) <- meta_final$cell_id

counts_col_ids <- paste0(id, '_', colnames(counts))
colnames(counts) <- counts_col_ids

if (!identical(colnames(counts), rownames(meta_final))) stop('Metadata rownames mismatch count columns')

message('Creating Seurat object...')
srt <- CreateSeuratObject(
  counts = counts,
  meta.data = meta_final,
  min.cells = 0,
  min.features = 0,
  project = id
)

message('Exporting...')
qs_path <- file.path(out_dir, paste0(id, '.qs'))
h5ad_path <- file.path(out_dir, paste0(id, '.h5ad'))
qs::qsave(srt, qs_path)
if (file.exists(h5ad_path)) file.remove(h5ad_path)
write_h5ad(srt, h5ad_path)

sample_counts <- as.data.frame(table(meta_final$sample_id), stringsAsFactors = FALSE)
colnames(sample_counts) <- c('sample_id', 'n_cells')
sample_counts <- sample_counts[order(sample_counts$sample_id), ]

report_path <- file.path(tmp_dir, 'gse179994_alignment_report.md')
lines <- c(
  '# gse179994 Alignment Report',
  '',
  paste0('- Generated at: ', format(Sys.time(), '%Y-%m-%d %H:%M:%S')),
  paste0('- Counts input: ', path.expand(input_counts)),
  paste0('- Metadata input: ', path.expand(input_meta)),
  paste0('- TCR input (kept separate): ', path.expand(input_tcr)),
  '',
  '## 1) Evidence and mapping',
  '- GEO accession: GSE179994',
  '- human-only lung cancer T cells',
  '- scRNA-seq and scTCR-seq; TCR is exported separately',
  '- raw expression input is rawCounts.rds.gz (dgCMatrix)',
  '- cell_id and metadata order are already aligned in the source files',
  '- disease = nsclc',
  '- organ = lung',
  '- tissue = tumor',
  '- lib = 10x_5',
  '- target = PD1',
  '- timepoint = pre/post from sample names',
  '- resp_standard = RECIST',
  '- response = NA (not encoded in local metadata)',
  '- anno_orig = cluster if available else celltype',
  '',
  '## 2) Dataset scale',
  paste0('- Unique cells final: ', ncol(counts)),
  paste0('- Genes final (union): ', nrow(counts)),
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
  paste0('- disease: ', paste(sort(unique(meta_final$disease)), collapse = ', ')),
  paste0('- organ: ', paste(sort(unique(meta_final$organ)), collapse = ', ')),
  paste0('- tissue: ', paste(sort(unique(meta_final$tissue)), collapse = ', ')),
  paste0('- lib: ', paste(sort(unique(meta_final$lib)), collapse = ', ')),
  paste0('- target: ', paste(sort(unique(meta_final$target)), collapse = ', ')),
  paste0('- timepoint: ', paste(sort(unique(meta_final$timepoint)), collapse = ', ')),
  paste0('- response: ', paste(sort(unique(na.omit(meta_final$response))), collapse = ', ')),
  paste0('- resp_standard: ', paste(sort(unique(na.omit(meta_final$resp_standard))), collapse = ', ')),
  '',
  '## 5) Evidence source',
  '- GSE179994_all.Tcell.rawCounts.rds.gz',
  '- GSE179994_Tcell.metadata.tsv.gz',
  '- GSE179994_all.scTCR.tsv.gz (kept separate)',
  '',
  '## 6) Per-sample cell counts',
  '',
  '| sample_id | n_cells |',
  '|---|---:|'
)
lines <- c(lines, apply(sample_counts, 1, function(x) paste0('| ', x[[1]], ' | ', x[[2]], ' |')))
writeLines(lines, report_path)

cat('Done\n')
cat('Report: ', report_path, '\n')
