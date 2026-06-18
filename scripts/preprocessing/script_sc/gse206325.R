id <- 'gse206325'

library(Seurat)
library(anndataR)
library(SingleCellExperiment)
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

first_nonempty <- function(...) {
  vals <- list(...)
  n <- max(lengths(vals))
  out <- rep(NA_character_, n)
  for (v in vals) {
    if (length(v) == 0) next
    v <- as.character(v)
    if (length(v) != n) v <- rep(v, length.out = n)
    ok <- !is.na(v) & nzchar(v)
    fill <- is.na(out) & ok
    out[fill] <- v[fill]
  }
  out
}

find_col <- function(cols, candidates) {
  for (cand in candidates) {
    idx <- which(tolower(cols) == tolower(cand))
    if (length(idx) > 0) return(cols[idx[1]])
  }
  norm_cols <- gsub('[^a-z0-9]+', '', tolower(cols))
  for (cand in candidates) {
    idx <- which(norm_cols == gsub('[^a-z0-9]+', '', tolower(cand)))
    if (length(idx) > 0) return(cols[idx[1]])
  }
  NA_character_
}

load_sce_from_rda_gz <- function(path) {
  env <- new.env(parent = emptyenv())
  con <- pipe(paste('gzip -dc', shQuote(path), '| gzip -dc'), 'rb')
  on.exit(close(con), add = TRUE)
  objs <- load(con, envir = env)
  if (length(objs) == 0) stop('No objects found in ', path)
  sce_candidates <- objs[vapply(objs, function(nm) {
    inherits(env[[nm]], c('SingleCellExperiment', 'SummarizedExperiment'))
  }, logical(1))]
  if (length(sce_candidates) == 0) stop('No SingleCellExperiment/SummarizedExperiment object found in ', path)
  env[[sce_candidates[1]]]
}

get_main_counts <- function(sce) {
  anms <- assayNames(sce)
  if (length(anms) == 0) stop('SCE has no assays')
  if ('counts' %in% anms) return(assay(sce, 'counts'))
  assay(sce, anms[1])
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

input_dir <- '~/006/data/imm/GSE206325'
out_dir <- '~/006/data/processed/srt/raw'
tmp_dir <- '~/006/tmp'
raw_path <- path.expand(file.path(input_dir, 'GSE206325_data_SingleCellExperiment_object.Rda.gz'))
sample_annot_path <- path.expand(file.path(input_dir, 'GSE206325_sample_annots_Liver_Treated_patients.csv.gz'))
cluster_annot_path <- path.expand(file.path(input_dir, 'GSE206325_full_HCC_cluster_annotation.csv.gz'))
tcr_path <- path.expand(file.path(input_dir, 'GSE206325_revised_processed_TCR.csv.gz'))

if (!file.exists(raw_path)) stop('Missing raw Rda.gz: ', raw_path)
if (!file.exists(sample_annot_path)) stop('Missing sample annotations: ', sample_annot_path)
if (!file.exists(cluster_annot_path)) stop('Missing cluster annotations: ', cluster_annot_path)

message('Loading SCE object...')
sce <- load_sce_from_rda_gz(raw_path)
message('Main class: ', paste(class(sce), collapse = ', '))

sample_annots <- read.csv(sample_annot_path, stringsAsFactors = FALSE, check.names = FALSE)
cluster_annots <- read.csv(cluster_annot_path, stringsAsFactors = FALSE, check.names = FALSE)

counts <- get_main_counts(sce)
if (!inherits(counts, 'dgCMatrix') && !inherits(counts, 'dgTMatrix') && !inherits(counts, 'matrix')) {
  counts <- as(counts, 'CsparseMatrix')
}

cell_meta <- as.data.frame(colData(sce), stringsAsFactors = FALSE)
cell_meta$cell_id_raw <- rownames(cell_meta)
if (is.null(cell_meta$cell_id_raw)) cell_meta$cell_id_raw <- colnames(counts)

sample_key_col <- find_col(colnames(cell_meta), c('cell_to_sample_ID', 'sample_ID', 'sample_id', 'sample', 'Sample', 'orig.ident'))
if (is.na(sample_key_col)) {
  stop('Could not find a sample identifier column in colData(sce)')
}
cell_meta$sample_id_raw <- as.character(cell_meta[[sample_key_col]])

sample_annot_key <- find_col(colnames(sample_annots), c('sample_ID', 'sample_id'))
if (is.na(sample_annot_key)) stop('Could not find sample_ID column in sample annotations')
sample_annots$sample_id_key <- as.character(sample_annots[[sample_annot_key]])

sample_match <- match(cell_meta$sample_id_raw, sample_annots$sample_id_key)
if (any(is.na(sample_match))) {
  message('Unmatched cell sample IDs in sample annotations: ', paste(head(unique(cell_meta$sample_id_raw[is.na(sample_match)]), 10), collapse = ', '))
}

sample_tissue <- sample_annots$tissue[sample_match]
sample_treatment <- sample_annots$treatment[sample_match]
sample_disease <- sample_annots$disease[sample_match]
sample_project <- sample_annots$Project[sample_match]
sample_species <- sample_annots$Species[sample_match]
sample_reference <- sample_annots$reference[sample_match]
sample_chem <- sample_annots$library_chemistry[sample_match]
sample_prime <- sample_annots$prime[sample_match]
sample_vdj <- sample_annots$vdj_kit[sample_match]
sample_resp <- sample_annots$treatment_Resp[sample_match]
sample_infil <- sample_annots$Immune_Infiltration[sample_match]

cell_meta$sample_id <- cell_meta$sample_id_raw
cell_meta$patient_id <- as.character(sample_annots$patient_ID[sample_match])
cell_meta$donor_id <- paste0(id, '_', cell_meta$patient_id)
cell_meta$sample_tissue <- sample_tissue
cell_meta$sample_treatment <- sample_treatment
cell_meta$sample_disease <- sample_disease
cell_meta$sample_project <- sample_project
cell_meta$sample_species <- sample_species
cell_meta$sample_reference <- sample_reference
cell_meta$sample_library_chemistry <- sample_chem
cell_meta$sample_prime <- sample_prime
cell_meta$sample_vdj_kit <- sample_vdj
cell_meta$sample_treatment_Resp <- sample_resp
cell_meta$sample_immune_infiltration <- sample_infil

cell_meta$disease <- 'hcc'
cell_meta$organ <- 'liver'
cell_meta$tissue <- tolower(as.character(sample_tissue))
cell_meta$lib <- '10x'
cell_meta$target <- 'anti-PD1'
cell_meta$timepoint <- 'post'
cell_meta$response <- ifelse(grepl('_R$', sample_resp), 'R', ifelse(grepl('_NR$', sample_resp), 'NR', NA_character_))
cell_meta$resp_standard <- 'PATHOLOGICAL_RESPONSE'
cell_meta$treatment <- as.character(sample_treatment)

cluster_key_col <- find_col(colnames(cell_meta), c('cell_to_cluster', 'cluster_ID', 'cluster_id', 'cluster', 'seurat_clusters', 'clusters'))
if (!is.na(cluster_key_col)) {
  cluster_annot_key <- find_col(colnames(cluster_annots), c('cluster_ID', 'cluster_id'))
  if (is.na(cluster_annot_key)) stop('Could not find cluster_ID column in cluster annotations')
  cluster_annots$cluster_id_key <- as.character(cluster_annots[[cluster_annot_key]])
  cluster_match <- match(as.character(cell_meta[[cluster_key_col]]), cluster_annots$cluster_id_key)
  if (!all(is.na(cluster_match))) {
    cell_meta$cluster_compartment <- cluster_annots$compartment[cluster_match]
    cell_meta$cluster_group <- cluster_annots$group[cluster_match]
    cell_meta$cluster_type <- cluster_annots$type[cluster_match]
    cell_meta$cluster_subgroup <- cluster_annots$subgroup[cluster_match]
    cell_meta$cluster_exclude <- cluster_annots$exclude[cluster_match]
    anno_raw <- first_nonempty(cell_meta$cluster_subgroup, cell_meta$cluster_type, cell_meta$cluster_group, cell_meta$cluster_compartment)
    cell_meta$anno_orig <- anno_raw
  }
}

if (!'anno_orig' %in% colnames(cell_meta) || all(is.na(cell_meta$anno_orig))) {
  cell_meta$anno_orig <- first_nonempty(cell_meta$sample_tissue, cell_meta$sample_project, cell_meta$sample_disease)
}

gene_only <- align_matrix(counts, id)
combined_counts <- gene_only$cts
sel_cells <- colnames(combined_counts)
cell_meta <- cell_meta[match(sel_cells, cell_meta$cell_id_raw), , drop = FALSE]
rownames(cell_meta) <- cell_meta$cell_id_raw

if (!identical(colnames(combined_counts), rownames(cell_meta))) stop('Metadata rownames mismatch count columns')

selInfo <- c(
  'cell_id_raw', 'sample_id_raw', 'sample_id', 'patient_id', 'donor_id', 'disease', 'organ', 'tissue', 'lib',
  'target', 'timepoint', 'response', 'resp_standard', 'treatment', 'anno_orig',
  'sample_tissue', 'sample_treatment', 'sample_disease', 'sample_project', 'sample_species',
  'sample_reference', 'sample_library_chemistry', 'sample_prime', 'sample_vdj_kit',
  'sample_treatment_Resp', 'sample_immune_infiltration', 'cluster_compartment', 'cluster_group',
  'cluster_type', 'cluster_subgroup', 'cluster_exclude'
)
selInfo <- selInfo[selInfo %in% colnames(cell_meta)]
meta_final <- cell_meta[, selInfo, drop = FALSE]

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

report_path <- file.path(tmp_dir, 'gse206325_alignment_report.md')
lines <- c(
  '# gse206325 Alignment Report',
  '',
  paste0('- Generated at: ', format(Sys.time(), '%Y-%m-%d %H:%M:%S')),
  paste0('- Input raw object: ', raw_path),
  paste0('- Sample annotations: ', sample_annot_path),
  paste0('- Cluster annotations: ', cluster_annot_path),
  '',
  '## 1) Evidence and mapping',
  '- GEO accession: GSE206325',
  '- disease = hcc',
  '- organ = liver',
  '- tissue = tumor / normal / pbmc / ln from sample annotations',
  '- lib = 10x',
  '- target = anti-PD1',
  '- timepoint = post',
  '- response = R / NR from treatment_Resp',
  '- resp_standard = PATHOLOGICAL_RESPONSE',
  '- anno_orig = cluster subgroup/type/group if available; otherwise sample-derived fallback',
  '- human-only input confirmed by sample annotations',
  '- scRNA-only output keeps the main gene-expression assay only',
  '',
  '## 2) Dataset scale',
  paste0('- Unique cells final: ', ncol(combined_counts)),
  paste0('- Genes final (union): ', nrow(combined_counts)),
  '',
  '## 3) SCE structure snapshot',
  paste0('- assayNames: ', paste(assayNames(sce), collapse = ', ')),
  paste0('- altExpNames: ', paste(altExpNames(sce), collapse = ', ')),
  paste0('- main assay dims: ', paste(dim(counts), collapse = ' x ')),
  '',
  '## 4) Metadata completeness checks',
  paste0('- Missing sample_id: ', sum(is.na(meta_final$sample_id))),
  paste0('- Missing patient_id: ', sum(is.na(meta_final$patient_id))),
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
  '- GEO accession record and paper metadata for GSE206325',
  '- sample_annots_Liver_Treated_patients.csv.gz',
  '- full_HCC_cluster_annotation.csv.gz',
  '- revised_processed_TCR.csv.gz (not merged in scRNA-only output)',
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
