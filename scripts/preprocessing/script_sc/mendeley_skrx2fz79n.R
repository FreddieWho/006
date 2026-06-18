id <- 'mendeley_skrx2fz79n'

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

input_dir <- '~/006/data/imm/mendeley_skrx2fz79n'
out_dir <- '~/006/data/processed/srt/raw'
tmp_dir <- '~/006/tmp'
raw_path <- path.expand(file.path(input_dir, 'WholeTissueList.rds'))

if (!file.exists(raw_path)) stop('Missing raw RDS: ', raw_path)

seurat_list <- readRDS(raw_path)
if (!is.list(seurat_list) || length(seurat_list) == 0) stop('Raw RDS is not a non-empty list')

celltype_names <- names(seurat_list)
if (is.null(celltype_names) || any(!nzchar(celltype_names))) stop('Expected named Seurat list for cell types')

pass1 <- list()
all_aligned <- character(0)

message('Pass 1: Computing union gene set...')
for (nm in celltype_names) {
  message('  Aligning ', nm)
  obj <- seurat_list[[nm]]
  cts <- GetAssayData(obj, assay = 'RNA', slot = 'counts')
  if (!inherits(cts, 'dgCMatrix') && !inherits(cts, 'dgTMatrix') && !inherits(cts, 'matrix')) {
    cts <- as(cts, 'CsparseMatrix')
  }
  genes <- rownames(cts)
  res <- gene_symbol_align(genes, gene_anno_list, alias_map)
  res_type_flt <- subset(res, gene_type %in% c('protein_coding', 'IG_C_gene', 'TR_C_gene'))
  if (nrow(res_type_flt) == 0) stop('No genes left after biotype filter in ', nm)

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
  all_aligned <- union(all_aligned, res_final$aligned_symbol)
  pass1[[nm]] <- list(cts = cts, res_final = res_final, meta = obj@meta.data)
  rm(obj, cts, res, res_type_flt, res_undup, res_dup, res_final); gc()
}

union_genes <- sort(all_aligned)
message('  Union genes: ', length(union_genes))

message('Pass 2: Combining subsets...')
combined_counts <- NULL
combined_meta <- NULL
align_stats <- list()

for (nm in celltype_names) {
  ar <- pass1[[nm]]
  cts <- ar$cts
  genes <- rownames(cts)
  res_final <- ar$res_final
  md <- ar$meta

  present <- intersect(union_genes, res_final$aligned_symbol)
  missing_genes <- setdiff(union_genes, present)

  cts_sub <- cts[match(res_final$original_id, genes), , drop = FALSE]
  rownames(cts_sub) <- res_final$aligned_symbol
  cts_sub <- cts_sub[present, , drop = FALSE]

  if (length(missing_genes) > 0) {
    zero_mat <- Matrix::sparseMatrix(
      i = integer(0), j = integer(0), x = numeric(0),
      dims = c(length(missing_genes), ncol(cts_sub)),
      dimnames = list(missing_genes, colnames(cts_sub))
    )
    cts_full <- rbind(cts_sub, zero_mat)
  } else {
    cts_full <- cts_sub
  }
  cts_full <- cts_full[union_genes, , drop = FALSE]

  sample_meta <- as.data.frame(md, stringsAsFactors = FALSE)
  sample_meta$cell_id_raw <- rownames(sample_meta)
  sample_meta$cell_id <- paste0(id, '_', sample_meta$cell_id_raw)
  sample_meta$sample_id <- paste0(id, '_', as.character(sample_meta$orig.ident), '_', as.character(sample_meta$Tissues))
  sample_meta$donor_id <- paste0(id, '_', as.character(sample_meta$orig.ident))
  sample_meta$disease <- 'hcc'
  sample_meta$organ <- 'liver'
  sample_meta$tissue <- ifelse(as.character(sample_meta$Tissues) == 'T', 'tumor', 'adj')
  sample_meta$lib <- '10x'
  sample_meta$target <- 'PD1'
  sample_meta$timepoint <- 'post'
  sample_meta$response <- NA_character_
  sample_meta$resp_standard <- 'RECIST'
  sample_meta$treatment <- 'PD1'
  sample_meta$anno_orig <- nm

  rownames(sample_meta) <- sample_meta$cell_id
  colnames(cts_full) <- sample_meta$cell_id

  if (is.null(combined_counts)) {
    combined_counts <- cts_full
    combined_meta <- sample_meta
    dup_ids <- character(0)
    new_ids <- colnames(cts_full)
  } else {
    dup_ids <- intersect(colnames(cts_full), colnames(combined_counts))
    new_ids <- setdiff(colnames(cts_full), dup_ids)

    if (length(dup_ids) > 0) {
      for (cid in dup_ids) {
        old_idx <- which(combined_meta$cell_id == cid)
        new_row <- sample_meta[sample_meta$cell_id == cid, , drop = FALSE]
        old_row <- combined_meta[old_idx, , drop = FALSE]
        check_same_sample_meta(old_row, new_row, c('cell_id_raw', 'sample_id', 'donor_id', 'disease', 'organ', 'tissue', 'lib', 'target', 'timepoint', 'treatment'), cid)
        combined_meta$anno_orig[old_idx] <- merge_label(combined_meta$anno_orig[old_idx], nm)
      }
    }

    if (length(new_ids) > 0) {
      combined_counts <- cbind(combined_counts, cts_full[, new_ids, drop = FALSE])
      combined_meta <- bind_rows(combined_meta, sample_meta[sample_meta$cell_id %in% new_ids, , drop = FALSE])
    }
  }

  align_stats[[nm]] <- list(
    raw_cells = ncol(ar$cts),
    raw_genes = nrow(ar$cts),
    aligned_genes = nrow(res_final),
    union_genes = length(union_genes),
    present_genes = length(present),
    n_dups = length(dup_ids),
    n_new = length(new_ids)
  )

  rm(cts_full, sample_meta); gc()
}

rownames(combined_meta) <- combined_meta$cell_id
combined_meta <- combined_meta[colnames(combined_counts), , drop = FALSE]

selInfo <- c(
  'cell_id_raw', 'cell_id', 'sample_id', 'donor_id', 'disease', 'organ', 'tissue', 'lib',
  'target', 'timepoint', 'response', 'resp_standard', 'treatment', 'anno_orig',
  'orig.ident', 'nCount_RNA', 'nFeature_RNA', 'Tissues', 'Pat_Tissues', 'percent.mito', 'percent.RP', 'DefineTypes'
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

report_path <- file.path(tmp_dir, 'mendeley_skrx2fz79n_alignment_report.md')
lines <- c(
  '# mendeley_skrx2fz79n Alignment Report',
  '',
  paste0('- Generated at: ', format(Sys.time(), '%Y-%m-%d %H:%M:%S')),
  paste0('- Input raw: ', raw_path),
  '',
  '## 1) Evidence and mapping',
  '- Mendeley Data title: Integrating single-cell and spatial transcriptomics to elucidate the cross-talk of SPP1+ macrophage and cancer associated fibroblast in HCC with immune excluded microenvironment',
  '- disease = hcc',
  '- organ = liver',
  '- tissue = tumor/adj from Tissues (T -> tumor, N -> adj)',
  '- lib = 10x',
  '- target = PD1',
  '- timepoint = post (registry-level context)',
  '- response = NA (not cell-mapped in local Seurat objects)',
  '- resp_standard = RECIST (registry-level context)',
  '- anno_orig = source Seurat list name',
  '- merged mode = duplicate cell IDs collapsed across subsets',
  '',
  '## 2) Dataset scale',
  paste0('- Unique cells final: ', ncol(combined_counts)),
  paste0('- Genes final (union): ', nrow(combined_counts)),
  '',
  '## 3) Per-celltype subset stats',
  ''
)

for (nm in names(align_stats)) {
  st <- align_stats[[nm]]
  lines <- c(lines,
    paste0('### ', nm),
    paste0('- Cells raw: ', st$raw_cells),
    paste0('- Genes raw: ', st$raw_genes),
    paste0('- Genes aligned+filtered: ', st$aligned_genes),
    paste0('- Genes present in union: ', st$present_genes),
    paste0('- Duplicate cells collapsed: ', st$n_dups),
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
  paste0('- response: ', if (all(is.na(meta_final$response))) 'NA' else paste(sort(unique(na.omit(meta_final$response))), collapse = ', ')),
  paste0('- resp_standard: ', paste(sort(unique(na.omit(meta_final$resp_standard))), collapse = ', ')),
  '',
  '## 6) Evidence source',
  '- Local evidence: ~/006/data/imm/mendeley_skrx2fz79n/WholeTissueList.rds',
  '- Registry evidence: ~/006/docs/data_collection.csv',
  '- Mendeley Data: https://data.mendeley.com/datasets/skrx2fz79n/1',
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
