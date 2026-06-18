id <- 'gse229772'

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
  list(cts = cts_align, genes = rownames(cts_align), stats = list(raw_genes = nrow(cts), aligned_genes = nrow(cts_align)))
}

read_sparse_norm_matrix <- function(path) {
  con <- gzfile(path, 'rt')
  on.exit(close(con), add = TRUE)

  header <- strsplit(readLines(con, n = 1), '\t', fixed = TRUE)[[1]]
  cell_ids <- sub('\\.([0-9]+)$', '-\\1', header)
  n_cells <- length(cell_ids)

  gene_ids <- character()
  row_idx <- 0L
  block_size <- 500L
  i_chunks <- list()
  j_chunks <- list()
  x_chunks <- list()
  tmp_i <- integer()
  tmp_j <- integer()
  tmp_x <- numeric()

  repeat {
    lines <- readLines(con, n = block_size)
    if (length(lines) == 0) break
    for (ln in lines) {
      row_idx <- row_idx + 1L
      parts <- strsplit(ln, '\t', fixed = TRUE)[[1]]
      gene_ids[row_idx] <- parts[1]
      vals <- suppressWarnings(as.numeric(parts[-1]))
      nz <- which(vals != 0)
      if (length(nz) > 0) {
        tmp_i <- c(tmp_i, rep.int(row_idx, length(nz)))
        tmp_j <- c(tmp_j, nz)
        tmp_x <- c(tmp_x, vals[nz])
      }
    }
    if (length(tmp_x) > 0) {
      i_chunks[[length(i_chunks) + 1L]] <- tmp_i
      j_chunks[[length(j_chunks) + 1L]] <- tmp_j
      x_chunks[[length(x_chunks) + 1L]] <- tmp_x
      tmp_i <- integer()
      tmp_j <- integer()
      tmp_x <- numeric()
    }
  }

  if (length(tmp_x) > 0) {
    i_chunks[[length(i_chunks) + 1L]] <- tmp_i
    j_chunks[[length(j_chunks) + 1L]] <- tmp_j
    x_chunks[[length(x_chunks) + 1L]] <- tmp_x
  }

  i <- if (length(i_chunks) > 0) unlist(i_chunks, use.names = FALSE) else integer()
  j <- if (length(j_chunks) > 0) unlist(j_chunks, use.names = FALSE) else integer()
  x <- if (length(x_chunks) > 0) unlist(x_chunks, use.names = FALSE) else numeric()

  mat <- Matrix::sparseMatrix(
    i = i,
    j = j,
    x = x,
    dims = c(row_idx, n_cells),
    dimnames = list(gene_ids, cell_ids)
  )
  mat
}

input_dir <- '~/006/data/imm/GSE229772'
out_dir <- '~/006/data/processed/srt/raw'
tmp_dir <- '~/006/tmp'

counts_path <- path.expand(file.path(input_dir, 'GSE229772_norm_counts.txt.gz'))
subtype_path <- path.expand(file.path(input_dir, 'GSE229772_cell_subtypes.txt.gz'))
sample_map_path <- path.expand(file.path(input_dir, 'GSE229772_scRNA_id_name_match.txt.gz'))

if (!file.exists(counts_path)) stop('Missing normalized counts: ', counts_path)
if (!file.exists(subtype_path)) stop('Missing cell subtype file: ', subtype_path)
if (!file.exists(sample_map_path)) stop('Missing sample map file: ', sample_map_path)

message('Reading normalized count matrix...')
cts_raw <- read_sparse_norm_matrix(counts_path)
message('Reading cell subtype annotations...')
cell_subtypes <- read.csv(subtype_path, sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)
message('Reading sample name map...')
sample_map <- read.csv(sample_map_path, sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)

if (!all(c('sample', 'subtype') %in% colnames(cell_subtypes))) {
  stop('Cell subtype file missing expected columns: sample, subtype')
}
if (!all(c('id', 'sample') %in% colnames(sample_map))) {
  stop('Sample map file missing expected columns: id, sample')
}

if (!identical(colnames(cts_raw), sub('\\.([0-9]+)$', '-\\1', cell_subtypes$sample))) {
  message('Cell subtype order differs from matrix columns; realigning by barcode')
}

subtype_idx <- match(colnames(cts_raw), sub('\\.([0-9]+)$', '-\\1', cell_subtypes$sample))
if (any(is.na(subtype_idx))) {
  stop('Some matrix cells are missing from the subtype file: ', paste(head(colnames(cts_raw)[is.na(subtype_idx)], 10), collapse = ', '))
}

cell_subtypes <- cell_subtypes[subtype_idx, , drop = FALSE]
rownames(cell_subtypes) <- sub('\\.([0-9]+)$', '-\\1', cell_subtypes$sample)

message('Aligning genes...')
aligned <- align_matrix(cts_raw, id)
counts <- aligned$cts

cell_meta <- data.frame(
  cell_id_raw = colnames(counts),
  cell_id = paste0(id, '_', colnames(counts)),
  sample_id_raw = NA_character_,
  sample_id = NA_character_,
  donor_id = NA_character_,
  disease = 'liver_cancer',
  organ = 'liver',
  tissue = 'tumor',
  lib = '10x_3p',
  target = 'immunotherapy',
  timepoint = NA_character_,
  response = NA_character_,
  resp_standard = NA_character_,
  treatment = 'immunotherapy',
  anno_orig = cell_subtypes$subtype,
  source_tag = 'normalized_counts',
  stringsAsFactors = FALSE
)
rownames(cell_meta) <- cell_meta$cell_id
colnames(counts) <- cell_meta$cell_id

sample_map$sample_name <- as.character(sample_map$sample)

selInfo <- c(
  'cell_id_raw', 'cell_id', 'sample_id_raw', 'sample_id', 'donor_id', 'disease', 'organ', 'tissue', 'lib',
  'target', 'timepoint', 'response', 'resp_standard', 'treatment', 'anno_orig', 'source_tag'
)
meta_final <- cell_meta[, selInfo, drop = FALSE]

if (!identical(colnames(counts), rownames(cell_meta))) stop('Metadata rownames mismatch count columns')

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

report_path <- file.path(tmp_dir, 'gse229772_alignment_report.md')
subtype_counts <- as.data.frame(table(cell_subtypes$subtype), stringsAsFactors = FALSE)
colnames(subtype_counts) <- c('subtype', 'n_cells')
subtype_counts <- subtype_counts[order(subtype_counts$subtype), ]

lines <- c(
  '# gse229772 Alignment Report',
  '',
  paste0('- Generated at: ', format(Sys.time(), '%Y-%m-%d %H:%M:%S')),
  paste0('- Input counts: ', counts_path),
  paste0('- Cell subtypes: ', subtype_path),
  paste0('- Sample map: ', sample_map_path),
  '',
  '## 1) Evidence and mapping',
  '- GEO title: Lineage and ecology define liver tumor evolution in response to treatment',
  '- human-only scRNA-seq only',
  '- accessible GEO matrix is normalized counts (not raw counts)',
  '- per-cell sample linkage is not preserved in the local files; sample_id fields are left NA',
  '- tissue = tumor',
  '- organ = liver',
  '- disease = liver_cancer',
  '- lib = 10x_3p',
  '- target = immunotherapy',
  '- anno_orig = cell_subtypes subtype',
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
  paste0('- response: ', paste(sort(unique(na.omit(meta_final$response))), collapse = ', ')),
  paste0('- resp_standard: ', paste(sort(unique(na.omit(meta_final$resp_standard))), collapse = ', ')),
  '',
  '## 5) Evidence source',
  '- GEO accession record and publication context indicate human scRNA-seq only',
  '- local accessible files: normalized matrix, cell subtypes, sample-name map',
  '',
  '## 6) Per-subtype cell counts',
  '',
  '| subtype | n_cells |',
  '|---|---:|'
)

lines <- c(lines, apply(subtype_counts, 1, function(x) paste0('| ', x[[1]], ' | ', x[[2]], ' |')))
writeLines(lines, report_path)

cat('Done\n')
cat('Report: ', report_path, '\n')
