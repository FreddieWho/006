id <- 'gse246613'

library(Seurat)
library(anndataR)
library(Matrix)
library(dplyr)
library(stringr)
library(qs)
library(rhdf5)

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

read_h5_index <- function(file, group) {
  file <- path.expand(file)
  ls <- h5ls(file, recursive = TRUE)
  if (any(ls$group == group & ls$name == '_index')) return(as.character(h5read(file, paste0(group, '/_index'))))
  if (any(ls$group == group & ls$name == 'index')) return(as.character(h5read(file, paste0(group, '/index'))))
  stop('No index dataset found under ', group, ' in ', file)
}

read_obs_col <- function(file, col) {
  file <- path.expand(file)
  ls <- h5ls(file, recursive = TRUE)
  val <- h5read(file, paste0('/obs/', col))
  if (is.list(val) && all(c('categories', 'codes') %in% names(val))) {
    return(as.character(val$categories[val$codes + 1]))
  }
  if ((is.integer(val) || is.numeric(val)) && any(ls$group == '/obs/__categories' & ls$name == col)) {
    cats <- h5read(file, paste0('/obs/__categories/', col))
    return(as.character(cats[val + 1]))
  }
  as.character(val)
}

read_sparse_group <- function(file, x_group, var_group, obs_group = '/obs') {
  file <- path.expand(file)
  data <- h5read(file, paste0(x_group, '/data'))
  indices <- h5read(file, paste0(x_group, '/indices'))
  indptr <- h5read(file, paste0(x_group, '/indptr'))
  gene_names <- read_h5_index(file, var_group)
  cell_names <- read_h5_index(file, obs_group)
  mat <- Matrix::sparseMatrix(
    i = as.integer(indices) + 1L,
    p = as.integer(indptr),
    x = as.numeric(data),
    dims = c(length(gene_names), length(cell_names)),
    dimnames = list(gene_names, cell_names)
  )
  mat
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

build_meta <- function(md, subset_label, source_tag, response_col = NULL) {
  md <- as.data.frame(md, stringsAsFactors = FALSE)
  md$cell_id_raw <- rownames(md)
  md$cell_id <- paste0(id, '_', md$cell_id_raw)
  md$sample_id <- paste0(id, '_', md$batch)
  md$donor_id <- paste0(id, '_', md$cohort)
  md$disease <- 'tnbc'
  md$organ <- 'breast'
  md$tissue <- 'tumor'
  md$lib <- '10x'
  md$target <- 'PD1+chemo'
  md$timepoint <- ifelse(md$treatment == 'Base', 'pre', 'on')
  md$response <- if (!is.null(response_col) && response_col %in% colnames(md)) md[[response_col]] else NA_character_
  md$resp_standard <- 'pCR'
  md$treatment <- md$treatment
  md$anno_orig <- subset_label
  md$source_tag <- source_tag
  rownames(md) <- md$cell_id
  md
}

merge_label <- function(existing, new_label) {
  x <- unique(c(unlist(strsplit(existing, ';', fixed = TRUE)), new_label))
  paste(sort(x[nzchar(x)]), collapse = ';')
}

check_same_sample_meta <- function(old_row, new_row, fields, cell_id) {
  for (f in fields) {
    if (!identical(old_row[[f]], new_row[[f]])) {
      stop('Metadata conflict for duplicated cell ', cell_id, ' on field ', f)
    }
  }
}

input_dir <- '~/006/data/combo/GSE246613'
out_dir <- '~/006/data/processed/srt/raw'
tmp_dir <- '~/006/tmp'

immune_path <- path.expand(file.path(input_dir, 'GSE246613_PembroRT_immune_R100_final.h5ad'))
nonimmune_path <- path.expand(file.path(input_dir, 'GSE246613_PembroRT_non_immune_cells.h5ad'))

if (!file.exists(immune_path)) stop('Missing immune h5ad: ', immune_path)
if (!file.exists(nonimmune_path)) stop('Missing non-immune h5ad: ', nonimmune_path)

subset_specs <- list(
  immune = list(path = immune_path, x_group = '/X', var_group = '/var', label = 'immune', response_col = 'pCR'),
  non_immune = list(path = nonimmune_path, x_group = '/raw/X', var_group = '/raw/var', label = 'non_immune', response_col = 'response_group')
)

message('Pass 1: Computing union gene set...')
pass1 <- list()
all_genes <- character(0)
align_stats <- list()

for (nm in names(subset_specs)) {
  sp <- subset_specs[[nm]]
  message('  Aligning ', nm)
  cts <- read_sparse_group(sp$path, sp$x_group, sp$var_group)
  meta_cols <- c('cohort', 'batch', 'celltype', 'subcluster', 'pCR', 'patient_treatment', 'hormone_receptor', 'CD45_enrich', 'treatment', 'combined_tcr', 'celltype_new', 'subtype_new', 'response_group', 'majority_voting', 'predicted_labels', 'conf_score')
  meta_cols <- meta_cols[meta_cols %in% unique(h5ls(sp$path, recursive = TRUE)$name)]
  md <- NULL
  for (col in unique(c('cohort', 'batch', 'treatment', 'patient_treatment', 'pCR', 'response_group', 'celltype', 'subcluster', 'celltype_new', 'subtype_new', 'majority_voting', 'predicted_labels', 'hormone_receptor', 'CD45_enrich', 'combined_tcr'))) {
    if (any(h5ls(sp$path, recursive = TRUE)$group == '/obs' & h5ls(sp$path, recursive = TRUE)$name == col)) {
      md[[col]] <- read_obs_col(sp$path, col)
    }
  }
  md <- as.data.frame(md, stringsAsFactors = FALSE)
  rownames(md) <- read_h5_index(sp$path, '/obs')
  aligned <- align_matrix(cts, nm)
  all_genes <- union(all_genes, aligned$genes)
  pass1[[nm]] <- list(cts = aligned$cts, genes = aligned$genes, meta = md, stats = aligned$stats)
  align_stats[[nm]] <- c(aligned$stats, raw_cells = ncol(cts), n_new = ncol(cts))
  rm(cts, md, aligned); gc()
}

union_genes <- sort(all_genes)
message('  Union genes: ', length(union_genes))

message('Pass 2: Combining subsets...')
combined_counts <- NULL
combined_meta <- NULL
for (nm in names(pass1)) {
  ar <- pass1[[nm]]
  cts <- ar$cts[intersect(union_genes, ar$genes), , drop = FALSE]
  missing_genes <- setdiff(union_genes, rownames(cts))
  if (length(missing_genes) > 0) {
    zero_mat <- Matrix::sparseMatrix(
      i = integer(0), j = integer(0), x = numeric(0),
      dims = c(length(missing_genes), ncol(cts)),
      dimnames = list(missing_genes, colnames(cts))
    )
    cts <- rbind(cts, zero_mat)
  }
  cts <- cts[union_genes, , drop = FALSE]

  md <- ar$meta
  md <- md[colnames(cts), , drop = FALSE]
  md <- build_meta(md, subset_label = nm, source_tag = nm, response_col = if (nm == 'immune') 'pCR' else 'response_group')
  colnames(cts) <- md$cell_id

  if (is.null(combined_counts)) {
    combined_counts <- cts
    combined_meta <- md
  } else {
    dup_ids <- intersect(colnames(cts), colnames(combined_counts))
    new_ids <- setdiff(colnames(cts), dup_ids)
    if (length(dup_ids) > 0) {
      for (cid in dup_ids) {
        old_row <- combined_meta[combined_meta$cell_id == cid, , drop = FALSE]
        new_row <- md[md$cell_id == cid, , drop = FALSE]
        check_same_sample_meta(old_row, new_row, c('cell_id_raw', 'sample_id', 'donor_id', 'disease', 'organ', 'tissue', 'lib', 'target', 'timepoint', 'treatment'), cid)
        combined_meta$anno_orig[combined_meta$cell_id == cid] <- merge_label(combined_meta$anno_orig[combined_meta$cell_id == cid], nm)
      }
    }
    if (length(new_ids) > 0) {
      combined_counts <- cbind(combined_counts, cts[, new_ids, drop = FALSE])
      combined_meta <- bind_rows(combined_meta, md[md$cell_id %in% new_ids, , drop = FALSE])
    }
  }
  rm(cts, md); gc()
}

rownames(combined_meta) <- combined_meta$cell_id
combined_meta <- combined_meta[colnames(combined_counts), , drop = FALSE]

selInfo <- c(
  'cell_id_raw', 'cell_id', 'sample_id', 'donor_id', 'disease', 'organ', 'tissue', 'lib',
  'target', 'timepoint', 'response', 'resp_standard', 'treatment', 'anno_orig',
  'batch', 'cohort', 'pCR', 'response_group', 'patient_treatment', 'hormone_receptor', 'CD45_enrich',
  'celltype', 'subcluster', 'celltype_new', 'subtype_new', 'majority_voting', 'predicted_labels', 'combined_tcr'
)
selInfo <- selInfo[selInfo %in% colnames(combined_meta)]
meta_final <- combined_meta[, selInfo, drop = FALSE]

if (!identical(colnames(combined_counts), rownames(meta_final))) stop('Metadata rownames mismatch count columns')
if (!identical(rownames(meta_final), meta_final$cell_id)) stop('Metadata rownames do not match cell_id')

message('Creating Seurat object...')
srt <- CreateSeuratObject(
  counts = combined_counts,
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

report_path <- file.path(tmp_dir, 'gse246613_alignment_report.md')
lines <- c(
  '# gse246613 Alignment Report',
  '',
  paste0('- Generated at: ', format(Sys.time(), '%Y-%m-%d %H:%M:%S')),
  paste0('- Input immune: ', immune_path),
  paste0('- Input non-immune: ', nonimmune_path),
  '',
  '## 1) Evidence and mapping',
  '- Registry: TNBC / breast / tumor / PD1+chemo / pre-on / pCR',
  '- Merged mode = immune + non-immune human cohort, mouse excluded',
  '- disease = tnbc',
  '- organ = breast',
  '- tissue = tumor',
  '- lib = 10x',
  '- target = PD1+chemo',
  '- timepoint = pre/on from treatment (Base -> pre, PD1/RTPD1 -> on)',
  '- response = pCR / response_group from source obs (preserved)',
  '- resp_standard = pCR',
  '- anno_orig = immune / non_immune',
  '',
  '## 2) Dataset scale',
  paste0('- Unique cells final: ', ncol(combined_counts)),
  paste0('- Genes final (union): ', nrow(combined_counts)),
  '',
  '## 3) Per-subset stats',
  ''
)

for (nm in names(pass1)) {
  st <- align_stats[[nm]]
  lines <- c(lines,
    paste0('### ', nm),
    paste0('- Cells raw: ', st$raw_cells),
    paste0('- Genes raw: ', st$raw_genes),
    paste0('- Genes aligned+filtered: ', st$aligned_genes),
    paste0('- New unique cells added: ', st$n_new),
    ''
  )
}

lines <- c(lines,
  '## 4) Metadata completeness checks',
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
  '## 5) Controlled-value snapshots',
  paste0('- disease: ', paste(sort(unique(meta_final$disease)), collapse = ', ')),
  paste0('- organ: ', paste(sort(unique(meta_final$organ)), collapse = ', ')),
  paste0('- tissue: ', paste(sort(unique(meta_final$tissue)), collapse = ', ')),
  paste0('- lib: ', paste(sort(unique(meta_final$lib)), collapse = ', ')),
  paste0('- target: ', paste(sort(unique(meta_final$target)), collapse = ', ')),
  paste0('- timepoint: ', paste(sort(unique(meta_final$timepoint)), collapse = ', ')),
  paste0('- response: ', paste(sort(unique(na.omit(meta_final$response))), collapse = ', ')),
  paste0('- resp_standard: ', paste(sort(unique(na.omit(meta_final$resp_standard))), collapse = ', ')),
  '',
  '## 6) Evidence source',
  '- Registry table: ~/006/docs/data_collection.csv',
  '- Mendeley Data: https://data.mendeley.com/datasets/skrx2fz79n/1',
  '- Supplementary workbook: ~/006/data/combo/GSE246613/mmc2 (3).xlsx',
  '- TCR table: ~/006/data/combo/GSE246613/GSE246613_TCRtable_annotation_expansion_detection.csv.gz',
  '',
  '## 7) Per-sample cell counts',
  '',
  '| sample_id | n_cells |',
  '|---|---:|'
)

lines <- c(lines, apply(sample_counts, 1, function(x) paste0('| ', x[[1]], ' | ', x[[2]], ' |')))
writeLines(lines, report_path)

cat('Done\n')
cat('Report: ', report_path, '\n')
