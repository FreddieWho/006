id <- 'gse266919'

library(Seurat)
library(anndataR)
library(SingleCellExperiment)
library(Matrix)
library(dplyr)
library(stringr)
library(qs)
library(readxl)

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

read_sce_file <- function(path) {
  if (grepl('\\.(gz|GZ)$', path)) {
    con <- pipe(paste('gzip -dc', shQuote(path), '| gzip -dc'), 'rb')
    on.exit(close(con), add = TRUE)
    return(readRDS(con))
  }
  readRDS(path)
}

align_sce_genes <- function(cts, gene_anno_list, alias_map) {
  genes <- rownames(cts)
  res <- gene_symbol_align(genes, gene_anno_list, alias_map)
  res_type_flt <- subset(res, gene_type %in% c('protein_coding', 'IG_C_gene', 'TR_C_gene'))
  if (nrow(res_type_flt) == 0) stop('No genes left after biotype filter')

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
    res_dup <- subset(res_dup, !(original_id %in% c('SPANXB1', 'XAGE2')))
  }

  res_final <- rbind(res_dup, res_undup)
  list(res_final = res_final, genes = genes, cts = cts)
}

input_dir <- '~/006/data/combo/GSE266919'
out_dir <- '~/006/data/processed/srt/raw'
tmp_dir <- '~/006/tmp'

subset_specs <- list(
  Bcell = 'GSE266919_Bcell.rds',
  CD4Tcell = 'GSE266919_CD4Tcell.rds.gz',
  CD8Tcell = 'GSE266919_CD8Tcell.rds',
  Myeloid = 'GSE266919_Myeloid.rds.gz',
  NKcell = 'GSE266919_NKcell.rds.gz'
)

# Pass 1: compute union of aligned genes across all subsets
message('Pass 1: Computing union gene set...')
all_aligned <- character(0)
align_results <- list()

for (celltype in names(subset_specs)) {
  file <- path.expand(file.path(input_dir, subset_specs[[celltype]]))
  message('  Aligning ', celltype)
  sce <- read_sce_file(file)
  cts <- assay(sce, 'counts')
  if (!inherits(cts, 'dgCMatrix') && !inherits(cts, 'dgTMatrix') && !inherits(cts, 'matrix')) {
    cts <- as(cts, 'CsparseMatrix')
  }
  ar <- align_sce_genes(cts, gene_anno_list, alias_map)
  all_aligned <- union(all_aligned, ar$res_final$aligned_symbol)
  align_results[[celltype]] <- ar
  rm(sce, cts); gc()
}

union_genes <- sort(all_aligned)
message('  Union genes: ', length(union_genes))

# Pass 2: build per-subset matrices with union genes, then combine
message('Pass 2: Combining subsets...')
combined_counts <- NULL
combined_meta <- NULL
align_stats <- list()

for (celltype in names(subset_specs)) {
  ar <- align_results[[celltype]]
  cts <- ar$cts
  genes <- ar$genes
  res_final <- ar$res_final

  # Map aligned genes to union_genes
  present <- intersect(union_genes, res_final$aligned_symbol)
  missing_genes <- setdiff(union_genes, present)

  # Build matrix with union_genes rows
  cts_sub <- cts[match(res_final$original_id, genes), , drop = FALSE]
  rownames(cts_sub) <- res_final$aligned_symbol
  cts_sub <- cts_sub[present, , drop = FALSE]

  if (length(missing_genes) > 0) {
    zero_mat <- Matrix::sparseMatrix(i = integer(0), j = integer(0), x = numeric(0),
                                     dims = c(length(missing_genes), ncol(cts_sub)),
                                     dimnames = list(missing_genes, colnames(cts_sub)))
    cts_full <- rbind(cts_sub, zero_mat)
  } else {
    cts_full <- cts_sub
  }
  cts_full <- cts_full[union_genes, , drop = FALSE]

  # Metadata
  coldata <- as.data.frame(colData(sce <- read_sce_file(path.expand(file.path(input_dir, subset_specs[[celltype]])))))
  sample_meta <- data.frame(cell_id_raw = colnames(cts_full), stringsAsFactors = FALSE)
  for (c in colnames(coldata)) sample_meta[[c]] <- coldata[[c]]
  sample_meta <- sample_meta %>%
    mutate(
      cell_id = paste0(id, '_', cell_id_raw),
      sample_id = paste0(id, '_', as.character(Sample)),
      donor_id = paste0(id, '_', as.character(Patient)),
      disease = 'tnbc',
      organ = 'breast',
      tissue = 'tumor',
      lib = '10x',
      target = ifelse(str_detect(as.character(Treatment), 'Anti-PD-L1'), 'PDL1', 'none'),
      timepoint = ifelse(as.character(Group) == 'Pre-treatment', 'Pre', 'Post'),
      response = ifelse(as.character(Response) == 'R', 'PR', 'NE'),
      resp_standard = ifelse(as.character(Response) == 'R', 'TUMOR_SIZE_CHANGE',
                             ifelse(as.character(Response) == 'NR', 'TUMOR_SIZE_CHANGE', NA_character_)),
      treatment = as.character(Treatment),
      anno_orig = celltype
    )

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
        old_idx <- which(combined_meta$cell_id_raw == cid)
        new_row <- sample_meta[sample_meta$cell_id_raw == cid, , drop = FALSE]
        old_row <- combined_meta[old_idx, , drop = FALSE]
        check_same_sample_meta(old_row, new_row, c('Sample', 'Patient', 'Group', 'Tissue', 'Treatment', 'Efficacy', 'Response'), cid)
        combined_meta$anno_orig[old_idx] <- merge_label(combined_meta$anno_orig[old_idx], celltype)
      }
    }

    if (length(new_ids) > 0) {
      combined_counts <- cbind(combined_counts, cts_full[, new_ids, drop = FALSE])
      combined_meta <- bind_rows(combined_meta, sample_meta[sample_meta$cell_id_raw %in% new_ids, , drop = FALSE])
    }
  }

  align_stats[[celltype]] <- list(
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

rownames(combined_meta) <- combined_meta$cell_id_raw

# Debug: check for mismatches
missing_in_meta <- setdiff(colnames(combined_counts), rownames(combined_meta))
missing_in_counts <- setdiff(rownames(combined_meta), colnames(combined_counts))
if (length(missing_in_meta) > 0) {
  message('ERROR: Cells in counts but not in meta (first 10): ', paste(head(missing_in_meta, 10), collapse = ', '))
}
if (length(missing_in_counts) > 0) {
  message('ERROR: Cells in meta but not in counts (first 10): ', paste(head(missing_in_counts, 10), collapse = ', '))
}

combined_meta <- combined_meta[colnames(combined_counts), , drop = FALSE]

extra_cols <- setdiff(colnames(combined_meta), c(
  'cell_id_raw', 'cell_id', 'sample_id', 'donor_id', 'disease', 'organ', 'tissue', 'lib',
  'target', 'timepoint', 'response', 'resp_standard', 'sex', 'age', 'stage', 'treatment', 'anno_orig'
))

selInfo <- c('cell_id_raw', 'cell_id', 'sample_id', 'donor_id', 'disease', 'organ', 'tissue', 'lib',
             'target', 'timepoint', 'response', 'resp_standard', 'treatment', 'anno_orig', extra_cols)
selInfo <- selInfo[selInfo %in% colnames(combined_meta)]
meta_final <- combined_meta[, selInfo, drop = FALSE]

if (!identical(colnames(combined_counts), rownames(meta_final))) stop('Metadata rownames mismatch count columns')

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

report_path <- file.path(tmp_dir, 'gse266919_alignment_report.md')
lines <- c(
  '# gse266919 Alignment Report',
  '',
  paste0('- Generated at: ', format(Sys.time(), '%Y-%m-%d %H:%M:%S')),
  paste0('- Input dir: ', input_dir),
  '',
  '## 1) Evidence and mapping',
  '- GEO title: Distinct Cellular Mechanisms Underlie Chemotherapies and Their Combinations with PD-L1 Checkpoint Inhibitor in TNBC',
  '- disease = tnbc',
  '- organ = breast',
  '- tissue = tumor',
  '- lib = 10x',
  '- target = PDL1 (Anti-PD-L1 treated) or none',
  '- timepoint = Pre/Post from Group',
  '- response = R/NR from tumor-size change',
  '- anno_orig = cell-type subset name (merged if cell in multiple subsets)',
  '- merged mode = duplicate cells collapsed across subsets',
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
  paste0('- response: ', paste(sort(unique(meta_final$response)), collapse = ', ')),
  paste0('- resp_standard: ', paste(sort(unique(na.omit(meta_final$resp_standard))), collapse = ', ')),
  '',
  '## 6) Evidence source',
  '- GEO: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE266919',
  '- OmicsDI: https://www.omicsdi.org/dataset/geo/GSE266919',
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
