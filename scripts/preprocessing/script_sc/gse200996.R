id <- 'gse200996'

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
  symbol_out <- character(length(ids))

  ens_idx <- which(id_type == 'ensembl')
  ent_idx <- which(id_type == 'entrez')
  sym_idx <- which(id_type == 'symbol')

  if (length(ens_idx) > 0) {
    ens_ids <- ids[ens_idx]
    map38 <- gene_v38$gene_name[match(ens_ids, gene_v38$gene_id)]
    miss38 <- is.na(map38)
    if (any(miss38)) {
      map19 <- gene_v19$gene_name[match(ens_ids[miss38], gene_v19$gene_id)]
      map38[miss38] <- map19
    }
    symbol_out[ens_idx] <- map38
  }

  if (length(ent_idx) > 0) symbol_out[ent_idx] <- NA_character_

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
    symbol_out[sym_idx] <- sym_res
  }

  symbol_out
}

get_gene_type <- function(symbols, gene_anno_list) {
  gene_v19 <- gene_anno_list$v19
  gene_v38 <- gene_anno_list$v38
  idx38 <- match(symbols, gene_v38$gene_name)
  type38 <- gene_v38$gene_type[idx38]
  missing <- is.na(type38)
  idx19 <- match(symbols[missing], gene_v19$gene_name)
  type38[missing] <- gene_v19$gene_type[idx19]
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

input_dir <- '~/006/data/combo/GSE200996'
out_dir <- '~/006/data/processed/srt/raw'
tmp_dir <- '~/006/tmp'

if (!dir.exists(path.expand(input_dir))) stop('Missing input dir: ', input_dir)

merge_label <- function(x) {
  vals <- unique(na.omit(as.character(x)))
  if (length(vals) == 0) return(NA_character_)
  paste(sort(vals), collapse = ';')
}

first_nonempty <- function(x) {
  x <- as.character(x)
  x <- x[!is.na(x) & nzchar(x)]
  if (length(x) == 0) NA_character_ else x[[1]]
}

parse_subset_label <- function(path) {
  nm <- basename(path)
  nm <- sub('^GSE200996_', '', nm)
  nm <- sub('\\.single\\.cell\\.meta\\.data\\.txt\\.gz$', '', nm)
  nm <- sub('\\.meta\\.data\\.txt\\.gz$', '', nm)
  nm <- gsub('\\.', '_', nm)
  nm
}

normalize_token <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x <- gsub('[[:space:]]+', '', x)
  tolower(x)
}

parse_matrix_label <- function(path) {
  nm <- basename(path)
  nm <- sub('^GSM[0-9]+_raw_feature_bc_matrix_', '', nm)
  nm <- sub('_GEX_sc_.*$', '', nm)
  nm
}

build_pbmc_pool_map <- function(path) {
  if (!file.exists(path)) stop('Missing PBMC hashtag mapping file: ', path)
  df <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  fixed_cols <- c('id', 'name', 'read', 'pattern', 'sequence', 'feature_type')
  pool_cols <- setdiff(colnames(df), fixed_cols)
  if (length(pool_cols) == 0) stop('No PBMC pool columns found in hashtag mapping file')

  mapped <- do.call(rbind, lapply(pool_cols, function(pool) {
    vals <- trimws(as.character(df[[pool]]))
    vals <- unique(vals[!is.na(vals) & nzchar(vals)])
    if (length(vals) == 0) return(NULL)
    data.frame(
      sample_id_raw = gsub('-', '_', vals),
      matrix_label = pool,
      stringsAsFactors = FALSE
    )
  }))

  if (is.null(mapped) || nrow(mapped) == 0) stop('No PBMC sample-to-pool mappings extracted')

  mapped$sample_token_norm <- normalize_token(mapped$sample_id_raw)
  mapped$matrix_label_norm <- normalize_token(mapped$matrix_label)
  mapped <- unique(mapped[, c('sample_id_raw', 'matrix_label', 'sample_token_norm', 'matrix_label_norm')])

  amb <- mapped %>%
    group_by(sample_token_norm) %>%
    summarise(n_pool = n_distinct(matrix_label_norm), pools = merge_label(matrix_label), .groups = 'drop') %>%
    filter(n_pool > 1)
  if (nrow(amb) > 0) {
    stop('Ambiguous PBMC sample-to-pool mapping for sample token(s): ', paste(head(amb$sample_token_norm, 10), collapse = ', '))
  }

  mapped
}

infer_matrix_label <- function(sample_id_raw, source_tissue, pbmc_pool_map) {
  sample_norm <- normalize_token(sample_id_raw)
  out <- rep(NA_character_, length(sample_norm))
  is_pbmc <- tolower(as.character(source_tissue)) == 'pbmc'

  if (any(is_pbmc)) {
    idx_pbmc <- match(sample_norm[is_pbmc], pbmc_pool_map$sample_token_norm)
    out[is_pbmc] <- pbmc_pool_map$matrix_label_norm[idx_pbmc]
  }

  if (any(!is_pbmc)) {
    out[!is_pbmc] <- sample_norm[!is_pbmc]
  }

  out
}

read_meta_file <- function(path) {
  meta <- read.delim(path, row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)
  rn <- rownames(meta)
  meta <- as.data.frame(lapply(meta, as.character), stringsAsFactors = FALSE)
  rownames(meta) <- rn
  meta$cell_id_raw <- rownames(meta)
  meta$barcode_key <- sub('_.*$', '', rownames(meta))
  meta$sample_id_raw <- sub('^[^_]+_', '', rownames(meta))
  meta$source_meta_file <- basename(path)
  meta$source_subset <- parse_subset_label(path)
  meta$source_tissue <- ifelse(grepl('PBMC', meta$source_subset, ignore.case = TRUE), 'pbmc', 'tumor')
  meta$source_organ <- ifelse(meta$source_tissue == 'pbmc', 'pbmc', 'liver')
  meta$anno_orig <- ifelse(is.na(meta$CellType_ID) | !nzchar(as.character(meta$CellType_ID)), meta$source_subset, paste0(meta$source_subset, ':', meta$CellType_ID))
  meta
}

read_gex_h5 <- function(path) {
  x <- Read10X_h5(path)
  if (is.list(x)) {
    if ('Gene Expression' %in% names(x)) return(x[['Gene Expression']])
    return(x[[1]])
  }
  x
}

align_counts <- function(cts) {
  genes <- rownames(cts)
  res <- gene_symbol_align(genes, gene_anno_list, alias_map)
  res_type_flt <- subset(res, gene_type %in% c('protein_coding', 'IG_C_gene', 'TR_C_gene'))
  if (nrow(res_type_flt) == 0) stop('No genes left after biotype filter')

  dup_gene <- unique(res_type_flt[duplicated(res_type_flt$aligned_symbol), 'aligned_symbol'])
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

  res_final <- rbind(as.data.frame(res_dup), as.data.frame(res_undup))
  idx <- match(res_final$original_id, genes)
  cts2 <- cts[idx, , drop = FALSE]
  cts2 <- Matrix::Matrix(cts2, sparse = TRUE)
  rownames(cts2) <- res_final$aligned_symbol
  cts2
}

build_union_matrix <- function(mat_list, union_genes) {
  lapply(mat_list, function(x) {
    miss <- setdiff(union_genes, rownames(x))
    if (length(miss) > 0) {
      z <- Matrix(0, nrow = length(miss), ncol = ncol(x), sparse = TRUE)
      rownames(z) <- miss
      colnames(z) <- colnames(x)
      x <- rbind(x, z)
    }
    x[union_genes, , drop = FALSE]
  })
}

meta_files <- list.files(input_dir, pattern = 'single\\.cell\\.meta\\.data\\.txt\\.gz$', full.names = TRUE)
h5_files <- list.files(input_dir, pattern = 'raw_feature_bc_matrix.*\\.h5$', full.names = TRUE)
pbmc_hashtag_path <- file.path(path.expand(input_dir), 'GSE200996_hashtag_info_sc_GEX_PBMC.csv')

if (length(meta_files) == 0) stop('No metadata files found in ', input_dir)
if (length(h5_files) == 0) stop('No H5 expression files found in ', input_dir)

message('Reading metadata files...')
meta_all <- bind_rows(lapply(meta_files, read_meta_file))

pbmc_pool_map <- build_pbmc_pool_map(pbmc_hashtag_path)
meta_all$matrix_label_norm <- infer_matrix_label(meta_all$sample_id_raw, meta_all$source_tissue, pbmc_pool_map)
if (any(is.na(meta_all$matrix_label_norm) | !nzchar(meta_all$matrix_label_norm))) {
  bad <- unique(meta_all$sample_id_raw[is.na(meta_all$matrix_label_norm) | !nzchar(meta_all$matrix_label_norm)])
  stop('Failed to infer matrix context label for sample_id_raw values: ', paste(head(bad, 12), collapse = ', '))
}
meta_all$join_key <- paste(meta_all$matrix_label_norm, meta_all$barcode_key, sep = '::')

raw_barcode_multi_context <- meta_all %>%
  group_by(barcode_key) %>%
  summarise(n_context = n_distinct(matrix_label_norm), .groups = 'drop') %>%
  filter(n_context > 1)

meta_index <- meta_all %>%
  group_by(join_key, barcode_key, matrix_label_norm) %>%
  summarise(
    cell_id_raw = first_nonempty(cell_id_raw),
    sample_id_raw = first_nonempty(sample_id_raw),
    source_meta_file = merge_label(source_meta_file),
    source_subset = merge_label(source_subset),
    source_tissue = first_nonempty(source_tissue),
    source_organ = first_nonempty(source_organ),
    Patient_ID = first_nonempty(Patient_ID),
    Stage = first_nonempty(Stage),
    Cohort = first_nonempty(Cohort),
    Path_response = first_nonempty(Path_response),
    CellType_ID = first_nonempty(CellType_ID),
    UMAP_1 = first_nonempty(UMAP_1),
    UMAP_2 = first_nonempty(UMAP_2),
    anno_orig = merge_label(anno_orig),
    .groups = 'drop'
  )

meta_index <- as.data.frame(meta_index, stringsAsFactors = FALSE)

rownames(meta_index) <- meta_index$cell_id_raw

message('Reading expression matrices...')
counts_list <- list()
meta_list <- list()
drop_stats <- data.frame(file = character(), n_total = integer(), n_matched = integer(), stringsAsFactors = FALSE)
drop_stats$matrix_label <- character(0)
drop_stats <- drop_stats[, c('file', 'matrix_label', 'n_total', 'n_matched')]

for (h5 in sort(h5_files)) {
  message('  ', basename(h5))
  cts <- read_gex_h5(h5)
  if (!inherits(cts, 'dgCMatrix') && !inherits(cts, 'dgTMatrix') && !inherits(cts, 'matrix')) {
    cts <- as(cts, 'CsparseMatrix')
  }

  if (nrow(cts) == 0 || ncol(cts) == 0) next

  matrix_label_norm <- normalize_token(parse_matrix_label(h5))
  col_keys <- sub('-1$', '', colnames(cts))
  join_keys <- paste(matrix_label_norm, col_keys, sep = '::')
  idx <- match(join_keys, meta_index$join_key)
  keep <- !is.na(idx)
  drop_stats <- rbind(drop_stats, data.frame(file = basename(h5), matrix_label = matrix_label_norm, n_total = length(col_keys), n_matched = sum(keep), stringsAsFactors = FALSE))
  if (!any(keep)) next

  cts <- cts[, keep, drop = FALSE]
  matched_meta <- meta_index[idx[keep], , drop = FALSE]
  matched_meta <- as.data.frame(matched_meta, stringsAsFactors = FALSE)
  matched_meta$cell_id <- paste0(id, '_', matched_meta$cell_id_raw)
  matched_meta$sample_id <- paste0(id, '_', matched_meta$sample_id_raw)
  matched_meta$donor_id <- paste0(id, '_', matched_meta$Patient_ID)
  matched_meta$disease <- 'hcc'
  matched_meta$organ <- matched_meta$source_organ
  matched_meta$tissue <- matched_meta$source_tissue
  matched_meta$lib <- '10x'
  matched_meta$target <- matched_meta$Cohort
  matched_meta$timepoint <- matched_meta$Stage
  matched_meta$response <- matched_meta$Path_response
  matched_meta$resp_standard <- 'PATHOLOGICAL_RESPONSE'
  matched_meta$source_matrix <- basename(h5)
  matched_meta$celltype_raw <- matched_meta$CellType_ID
  matched_meta$stage_raw <- matched_meta$Stage
  matched_meta$cohort_raw <- matched_meta$Cohort
  matched_meta$path_response_raw <- matched_meta$Path_response
  rownames(matched_meta) <- matched_meta$cell_id

  colnames(cts) <- matched_meta$cell_id

  cts <- align_counts(cts)
  counts_list[[basename(h5)]] <- cts
  meta_list[[basename(h5)]] <- matched_meta[, c('cell_id', 'cell_id_raw', 'sample_id', 'sample_id_raw', 'donor_id', 'disease', 'organ', 'tissue', 'lib', 'target', 'timepoint', 'response', 'resp_standard', 'source_matrix', 'source_meta_file', 'source_subset', 'anno_orig', 'celltype_raw', 'stage_raw', 'cohort_raw', 'path_response_raw', 'Patient_ID', 'Stage', 'Cohort', 'Path_response', 'UMAP_1', 'UMAP_2'), drop = FALSE]
}

if (length(counts_list) == 0) stop('No matrices were matched to metadata')

union_genes <- sort(unique(unlist(lapply(counts_list, rownames), use.names = FALSE)))
counts_list <- build_union_matrix(counts_list, union_genes)

combined_counts <- do.call(cbind, counts_list)
combined_meta <- bind_rows(meta_list)

combined_meta <- combined_meta %>%
  group_by(cell_id) %>%
  summarise(
    cell_id_raw = first_nonempty(cell_id_raw),
    sample_id = first_nonempty(sample_id),
    sample_id_raw = first_nonempty(sample_id_raw),
    donor_id = first_nonempty(donor_id),
    disease = first_nonempty(disease),
    organ = first_nonempty(organ),
    tissue = first_nonempty(tissue),
    lib = first_nonempty(lib),
    target = first_nonempty(target),
    timepoint = first_nonempty(timepoint),
    response = first_nonempty(response),
    resp_standard = first_nonempty(resp_standard),
    source_matrix = merge_label(source_matrix),
    source_meta_file = merge_label(source_meta_file),
    source_subset = merge_label(source_subset),
    anno_orig = merge_label(anno_orig),
    celltype_raw = merge_label(celltype_raw),
    stage_raw = first_nonempty(stage_raw),
    cohort_raw = first_nonempty(cohort_raw),
    path_response_raw = first_nonempty(path_response_raw),
    Patient_ID = first_nonempty(Patient_ID),
    Stage = first_nonempty(Stage),
    Cohort = first_nonempty(Cohort),
    Path_response = first_nonempty(Path_response),
    UMAP_1 = first_nonempty(UMAP_1),
    UMAP_2 = first_nonempty(UMAP_2),
    .groups = 'drop'
  )

combined_meta <- as.data.frame(combined_meta, stringsAsFactors = FALSE)
combined_meta$cell_id <- paste0(id, '_', combined_meta$cell_id_raw)
rownames(combined_meta) <- combined_meta$cell_id

keep_cols <- !duplicated(colnames(combined_counts))
if (any(!keep_cols)) {
  combined_counts <- combined_counts[, keep_cols, drop = FALSE]
}

meta_idx <- match(colnames(combined_counts), rownames(combined_meta))
keep_final <- !is.na(meta_idx)
combined_counts <- combined_counts[, keep_final, drop = FALSE]
combined_meta <- combined_meta[meta_idx[keep_final], , drop = FALSE]
stopifnot(all(colnames(combined_counts) == rownames(combined_meta)))

srt <- CreateSeuratObject(
  counts = combined_counts,
  meta.data = combined_meta,
  min.cells = 0,
  min.features = 0,
  project = id
)

qs_path <- file.path(out_dir, paste0(id, '.qs'))
h5ad_path <- file.path(out_dir, paste0(id, '.h5ad'))
report_path <- file.path(tmp_dir, 'gse200996_alignment_report.md')

qs::qsave(srt, qs_path)
if (file.exists(h5ad_path)) file.remove(h5ad_path)
write_h5ad(srt, h5ad_path)

sample_counts <- as.data.frame(table(combined_meta$sample_id), stringsAsFactors = FALSE)
colnames(sample_counts) <- c('sample_id', 'n_cells')
sample_counts <- sample_counts[order(sample_counts$sample_id), ]

report_lines <- c(
  '# gse200996 Alignment Report',
  '',
  paste0('- Generated at: ', format(Sys.time(), '%Y-%m-%d %H:%M:%S')),
  paste0('- Metadata files: ', length(meta_files)),
  paste0('- H5 files: ', length(h5_files)),
  paste0('- Matched matrices: ', length(counts_list)),
  '',
  '## Evidence and mapping',
  '- disease = hcc',
  '- organ = liver / pbmc (source-specific)',
  '- tissue = tumor / pbmc (source-specific)',
  '- lib = 10x',
  '- target = Cohort raw values (Combo / Mono)',
  '- timepoint = Stage raw values',
   '- response = Path_response raw values',
   '- resp_standard = PATHOLOGICAL_RESPONSE',
   '- anno_orig = source subset + CellType_ID',
   '- barcode mapping uses matrix/library context + barcode key (not bare global 16nt barcode)',
   '- PBMC context mapping source: GSE200996_hashtag_info_sc_GEX_PBMC.csv',
   paste0('- Raw barcodes seen in multiple contexts: ', nrow(raw_barcode_multi_context)),
   '',
  '## Dataset scale',
  paste0('- Final cells: ', ncol(combined_counts)),
  paste0('- Final genes: ', nrow(combined_counts)),
  paste0('- Cells after final cell-id match: ', ncol(combined_counts)),
  paste0('- Dropped unmatched cells: ', sum(drop_stats$n_total - drop_stats$n_matched)),
  '',
  '## Required field missingness',
  paste0('- cell_id: ', sum(is.na(combined_meta$cell_id))),
  paste0('- sample_id: ', sum(is.na(combined_meta$sample_id))),
  paste0('- donor_id: ', sum(is.na(combined_meta$donor_id))),
  paste0('- disease: ', sum(is.na(combined_meta$disease))),
  paste0('- organ: ', sum(is.na(combined_meta$organ))),
  paste0('- tissue: ', sum(is.na(combined_meta$tissue))),
  paste0('- lib: ', sum(is.na(combined_meta$lib))),
  '',
  '## Controlled vocabulary snapshots',
  paste0('- organ: ', paste(sort(unique(na.omit(combined_meta$organ))), collapse = ', ')),
  paste0('- tissue: ', paste(sort(unique(na.omit(combined_meta$tissue))), collapse = ', ')),
  paste0('- lib: ', paste(sort(unique(na.omit(combined_meta$lib))), collapse = ', ')),
  paste0('- target: ', paste(sort(unique(na.omit(combined_meta$target))), collapse = ', ')),
  paste0('- timepoint: ', paste(sort(unique(na.omit(combined_meta$timepoint))), collapse = ', ')),
  paste0('- response: ', paste(sort(unique(na.omit(combined_meta$response))), collapse = ', ')),
  '',
  '## Per-sample cell counts',
  '',
  '| sample_id | n_cells |',
  '|---|---:|'
)
report_lines <- c(report_lines, apply(sample_counts, 1, function(x) paste0('| ', x[[1]], ' | ', x[[2]], ' |')))
writeLines(report_lines, report_path)

cat('Done\n')
cat('Report: ', report_path, '\n')
