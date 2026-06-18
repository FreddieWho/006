id <- 'prjca001063'

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
  if (any(ls$group == group & ls$name == '_index')) {
    return(as.character(h5read(file, paste0(group, '/_index'))))
  }
  if (any(ls$group == group & ls$name == 'index')) {
    return(as.character(h5read(file, paste0(group, '/index'))))
  }
  stop('No index dataset found under ', group, ' in ', file)
}

read_obs_col <- function(file, col) {
  file <- path.expand(file)
  val <- h5read(file, paste0('/obs/', col))
  if ((is.integer(val) || is.numeric(val)) && any(h5ls(file, recursive = TRUE)$group == '/obs/__categories' & h5ls(file, recursive = TRUE)$name == col)) {
    cats <- h5read(file, paste0('/obs/__categories/', col))
    return(as.character(cats[val + 1]))
  }
  as.character(val)
}

read_h5ad_matrix <- function(file) {
  file <- path.expand(file)
  data <- h5read(file, '/X/data')
  indices <- h5read(file, '/X/indices')
  indptr <- h5read(file, '/X/indptr')
  gene_names <- read_h5_index(file, '/var')
  cell_names <- read_h5_index(file, '/obs')
  mat <- Matrix::sparseMatrix(
    i = as.integer(indices) + 1L,
    p = as.integer(indptr),
    x = as.numeric(data),
    dims = c(length(gene_names), length(cell_names)),
    dimnames = list(gene_names, cell_names)
  )
  mat
}

input_dir <- '~/006/data/imm/PRJCA001063'
out_dir <- '~/006/data/processed/srt/raw'
tmp_dir <- '~/006/tmp'
raw_path <- path.expand(file.path(input_dir, 'StdWf1_PRJCA001063_CRC_besca2.raw.h5ad'))
ann_path <- path.expand(file.path(input_dir, 'StdWf1_PRJCA001063_CRC_besca2.annotated.h5ad'))

if (!file.exists(raw_path)) stop('Missing raw h5ad: ', raw_path)

cts_all <- read_h5ad_matrix(raw_path)

raw_cell_cols <- c('CELL', 'CONDITION', 'Patient', 'Type', 'Cell_type')
raw_meta <- data.frame(
  CELL = read_obs_col(raw_path, 'CELL'),
  CONDITION = read_obs_col(raw_path, 'CONDITION'),
  Patient = read_obs_col(raw_path, 'Patient'),
  Type = read_obs_col(raw_path, 'Type'),
  Cell_type = read_obs_col(raw_path, 'Cell_type'),
  stringsAsFactors = FALSE
)

if (any(duplicated(raw_meta$CELL))) stop('Duplicated CELL IDs in raw metadata')
if (!identical(raw_meta$CELL, colnames(cts_all))) stop('CELL order mismatch between obs and matrix')

ann_cols <- c('percent_mito', 'n_counts', 'n_genes', 'leiden', 'celltype0', 'celltype1', 'celltype2', 'celltype3', 'dblabel')
ann_meta <- NULL
if (file.exists(ann_path)) {
  ann_cells <- read_h5_index(ann_path, '/obs')
  ann_meta <- data.frame(
    CELL = ann_cells,
    percent_mito = suppressWarnings(as.numeric(read_obs_col(ann_path, 'percent_mito'))),
    n_counts = suppressWarnings(as.numeric(read_obs_col(ann_path, 'n_counts'))),
    n_genes = suppressWarnings(as.numeric(read_obs_col(ann_path, 'n_genes'))),
    leiden = read_obs_col(ann_path, 'leiden'),
    celltype0 = read_obs_col(ann_path, 'celltype0'),
    celltype1 = read_obs_col(ann_path, 'celltype1'),
    celltype2 = read_obs_col(ann_path, 'celltype2'),
    celltype3 = read_obs_col(ann_path, 'celltype3'),
    dblabel = read_obs_col(ann_path, 'dblabel'),
    stringsAsFactors = FALSE
  )
}

meta <- raw_meta %>%
  mutate(
    cell_id_raw = CELL,
    cell_id = paste0(id, '_', CELL),
    sample_id = paste0(id, '_', Patient),
    donor_id = paste0(id, '_', sub('^[TN]', '', Patient)),
    disease = 'pdac',
    organ = 'pancreas',
    tissue = ifelse(Type == 'T', 'tumor', 'normal'),
    lib = '10x',
    target = NA_character_,
    timepoint = NA_character_,
    response = 'none',
    resp_standard = NA_character_,
    treatment = NA_character_,
    anno_orig = Cell_type
  )

if (!is.null(ann_meta)) {
  meta <- meta %>% left_join(ann_meta, by = 'CELL')
}

meta$cell_id <- paste0(id, '_', meta$CELL)
rownames(meta) <- meta$cell_id

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
dimnames(cts_alignGene) <- list(res_final$aligned_symbol, meta$cell_id)

if (any(duplicated(rownames(cts_alignGene)))) stop('Aligned gene names duplicated after dedup')
if (!identical(colnames(cts_alignGene), meta$cell_id)) stop('Metadata rownames mismatch count columns')

meta_final <- meta[, c(
  'cell_id_raw', 'cell_id', 'sample_id', 'donor_id', 'disease', 'organ', 'tissue', 'lib',
  'target', 'timepoint', 'response', 'resp_standard', 'treatment', 'CELL', 'CONDITION', 'Patient',
  'Type', 'Cell_type', 'anno_orig', 'percent_mito', 'n_counts', 'n_genes', 'leiden',
  'celltype0', 'celltype1', 'celltype2', 'celltype3', 'dblabel'
), drop = FALSE]

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

report_path <- file.path(tmp_dir, 'prjca001063_alignment_report.md')
lines <- c(
  '# prjca001063 Alignment Report',
  '',
  paste0('- Generated at: ', format(Sys.time(), '%Y-%m-%d %H:%M:%S')),
  paste0('- Input raw: ', raw_path),
  paste0('- Input annotated: ', ann_path),
  '',
  '## 1) Evidence and mapping',
  '- BioProject/PMC evidence supports PDAC / pancreatic cancer',
  '- organ = pancreas',
  '- tissue = tumor/normal from Type',
  '- lib = 10x (coarse; exact chemistry not explicit in available evidence)',
  '- target = NA (no treatment target)',
  '- timepoint = NA (not a longitudinal cohort)',
  '- response = none',
  '- resp_standard = NA (missing-by-design)',
  '- anno_orig = Cell_type',
  '',
  '## 2) Dataset scale',
  paste0('- Raw cells: ', ncol(cts_all)),
  paste0('- Raw genes: ', nrow(cts_all)),
  paste0('- Aligned+filtered genes: ', nrow(cts_alignGene)),
  paste0('- Final cells: ', ncol(cts_alignGene)),
  paste0('- Raw/annotated cell overlap: ', length(intersect(raw_meta$CELL, if (!is.null(ann_meta)) ann_meta$CELL else character(0)))),
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
  paste0('- response: ', paste(sort(unique(meta_final$response)), collapse = ', ')),
  '',
  '## 5) Per-sample cell counts',
  '',
  '| sample_id | n_cells |',
  '|---|---:|'
)

lines <- c(lines, apply(sample_counts, 1, function(x) paste0('| ', x[[1]], ' | ', x[[2]], ' |')))
writeLines(lines, report_path)

cat('Done\n')
cat('Report: ', report_path, '\n')
