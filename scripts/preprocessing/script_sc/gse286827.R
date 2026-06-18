id <- 'gse286827'

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

read_double_gzip_rds <- function(path) {
  path <- path.expand(path)
  tmp <- tempfile(fileext = '.rds')
  on.exit(unlink(tmp), add = TRUE)
  status <- system2('gzip', c('-dc', path), stdout = tmp, stderr = FALSE)
  if (!is.numeric(status) || status != 0) stop('Failed to decompress file: ', path)
  readRDS(tmp)
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

parse_geo_soft <- function(url) {
  con <- gzcon(url(url, open = 'rb'))
  on.exit(close(con), add = TRUE)
  lines <- readLines(con, warn = FALSE)

  samples <- list()
  cur <- NULL
  flush_cur <- function() {
    if (!is.null(cur)) samples[[length(samples) + 1]] <<- cur
  }

  for (line in lines) {
    if (startsWith(line, '^SAMPLE = ')) {
      flush_cur()
      cur <- list(chars = character())
    } else if (startsWith(line, '!Sample_geo_accession = ')) {
      cur$gsm <- sub('^!Sample_geo_accession = ', '', line)
    } else if (startsWith(line, '!Sample_title = ')) {
      cur$title <- sub('^!Sample_title = ', '', line)
    } else if (startsWith(line, '!Sample_source_name_ch1 = ')) {
      cur$source_name <- sub('^!Sample_source_name_ch1 = ', '', line)
    } else if (startsWith(line, '!Sample_description = ')) {
      cur$description <- sub('^!Sample_description = ', '', line)
    } else if (startsWith(line, '!Sample_characteristics_ch1 = ')) {
      cur$chars <- c(cur$chars, sub('^!Sample_characteristics_ch1 = ', '', line))
    } else if (startsWith(line, '!Sample_extract_protocol_ch1 = ')) {
      cur$extract_protocol <- sub('^!Sample_extract_protocol_ch1 = ', '', line)
    } else if (startsWith(line, '!Sample_data_processing = ')) {
      cur$data_processing <- sub('^!Sample_data_processing = ', '', line)
    }
  }
  flush_cur()

  out <- lapply(samples, function(s) {
    chars <- s$chars
    get_char <- function(key) {
      hits <- chars[str_detect(chars, paste0('^', key, ': '))]
      if (length(hits) == 0) return(NA_character_)
      sub(paste0('^', key, ': '), '', hits[1])
    }
    data.frame(
      gsm = s$gsm,
      title = s$title,
      source_name = s$source_name,
      description = s$description,
      treatment = get_char('treatment'),
      tissue = get_char('tissue'),
      hpv_p16 = get_char('hpv p16'),
      ici_duration_days = get_char('ici duration_days'),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, out)
}

geo_url <- 'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE286nnn/GSE286827/soft/GSE286827_family.soft.gz'
geo_all <- parse_geo_soft(geo_url)
geo_gex <- subset(geo_all, str_detect(title, '_GEX$'))
geo_gex$sample_token <- sub('_GEX$', '', geo_gex$title)
geo_gex$patient_id <- str_match(geo_gex$sample_token, '^(P[0-9]+)')[, 2]
geo_gex$timepoint <- str_match(geo_gex$sample_token, '_(pre|post2?|post)$')[, 2]
geo_gex$timepoint <- ifelse(geo_gex$timepoint == 'pre', 'Pre', ifelse(geo_gex$timepoint == 'post2', 'Post2', 'Post'))
geo_gex$disease <- 'hnscc'
geo_gex$organ <- tolower(geo_gex$source_name)
geo_gex$tissue <- 'tumor'
geo_gex$lib <- '10x_5_v2'
geo_gex$response <- 'none'
geo_gex$resp_standard <- NA_character_
geo_gex$target <- ifelse(geo_gex$treatment == 'D', 'PDL1', 'PDL1;CTLA4')
geo_gex$treatment_label <- ifelse(geo_gex$treatment == 'D', 'Durvalumab', 'Durvalumab+Tremelimumab')

input_dir <- '~/006/data/combo/GSE286827'
out_dir <- '~/006/data/processed/srt/raw'
tmp_dir <- '~/006/tmp'

subset_specs <- list(
  immunecells_broad = 'GSE286827_countdata_immunecells_broad.rds.gz',
  Malignant = 'GSE286827_countdata_Malignant.rds.gz',
  TAMs = 'GSE286827_countdata_TAMs.rds.gz',
  TumorSpecific_CD8T = 'GSE286827_countdata_TumorSpecific_CD8T.rds.gz',
  Bcells = 'GSE286827_countdata_Bcells.rds.gz',
  CD4T = 'GSE286827_countdata_CD4T.rds.gz',
  CD8T = 'GSE286827_countdata_CD8T.rds.gz'
)

combined_counts <- NULL
combined_meta <- NULL
source_map <- data.frame()
gene_ref <- NULL
align_stats <- list()

for (celltype in names(subset_specs)) {
  file <- file.path(input_dir, subset_specs[[celltype]])
  message('Processing ', celltype, ' from ', basename(file))
  cts <- read_double_gzip_rds(file)
  if (!(inherits(cts, 'dgCMatrix') || inherits(cts, 'dgTMatrix') || inherits(cts, 'matrix'))) {
    stop('Unsupported object type in ', file)
  }
  if (is.null(rownames(cts)) || is.null(colnames(cts))) stop('Matrix missing dimnames: ', file)

  sample_token <- str_match(colnames(cts), '^(P[0-9]+_(?:pre|post2?|post))_')[, 2]
  if (any(is.na(sample_token))) stop('Failed to parse sample_token from cell IDs in ', file)

  sample_meta <- data.frame(cell_id = colnames(cts), sample_token = sample_token, stringsAsFactors = FALSE) %>%
    left_join(geo_gex, by = 'sample_token')
  if (any(is.na(sample_meta$gsm))) stop('Unmapped cells to GEO sample metadata in ', file)

  sample_meta <- sample_meta %>%
    mutate(
      sample_id = paste0(id, '_', sample_token),
      donor_id = paste0(id, '_', patient_id),
      disease = disease,
      organ = tolower(source_name),
      tissue = 'tumor',
      target = target,
      timepoint = timepoint,
      response = response,
      resp_standard = resp_standard,
      lib = lib,
      treatment = treatment_label,
      sex = NA_character_,
      age = NA_real_,
      stage = NA_character_,
      anno_orig = celltype
    )

  genes <- rownames(cts)
  res <- gene_symbol_align(genes, gene_anno_list, alias_map)
  res_type_flt <- subset(res, gene_type %in% c('protein_coding', 'IG_C_gene', 'TR_C_gene'))
  if (nrow(res_type_flt) == 0) stop('No genes left after biotype filter: ', file)

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
  sel_idx <- match(res_final$original_id, genes)
  cts_align <- cts[sel_idx, , drop = FALSE]
  cts_align <- Matrix::Matrix(cts_align, sparse = TRUE)
  dimnames(cts_align) <- list(res_final$aligned_symbol, colnames(cts))
  if (any(duplicated(rownames(cts_align)))) stop('Duplicated aligned genes in ', file)

  if (is.null(gene_ref)) {
    gene_ref <- rownames(cts_align)
  } else if (!identical(gene_ref, rownames(cts_align))) {
    stop('Gene reference mismatch across cell-type subsets: ', celltype)
  }

  # merge duplicate cells across subsets by exact cell_id and annotate provenance in anno_orig
  if (is.null(combined_counts)) {
    combined_counts <- cts_align
    combined_meta <- sample_meta
    dup_ids <- character(0)
    new_ids <- colnames(cts_align)
  } else {
    dup_ids <- intersect(colnames(cts_align), colnames(combined_counts))
    new_ids <- setdiff(colnames(cts_align), dup_ids)

    if (length(dup_ids) > 0) {
      for (cid in dup_ids) {
        old_idx <- which(combined_meta$cell_id == cid)
        new_row <- sample_meta[sample_meta$cell_id == cid, , drop = FALSE]
        old_row <- combined_meta[old_idx, , drop = FALSE]
        check_same_sample_meta(old_row, new_row, c('sample_token', 'sample_id', 'donor_id', 'disease', 'organ', 'tissue', 'lib', 'target', 'timepoint', 'response', 'resp_standard', 'treatment'), cid)
        combined_meta$anno_orig[old_idx] <- merge_label(combined_meta$anno_orig[old_idx], celltype)
      }
    }

    if (length(new_ids) > 0) {
      combined_counts <- cbind(combined_counts, cts_align[, new_ids, drop = FALSE])
      combined_meta <- bind_rows(combined_meta, sample_meta[sample_meta$cell_id %in% new_ids, , drop = FALSE])
    }
  }

  align_stats[[celltype]] <- list(
    raw_cells = ncol(cts),
    raw_genes = nrow(cts),
    aligned_genes = nrow(cts_align),
    n_dups = length(dup_ids),
    n_new = length(new_ids)
  )
}

rownames(combined_meta) <- combined_meta$cell_id
combined_meta <- combined_meta[colnames(combined_counts), , drop = FALSE]

selInfo <- c('cell_id', 'sample_id', 'donor_id', 'disease', 'organ', 'tissue', 'lib', 'target', 'timepoint', 'response', 'resp_standard', 'sex', 'age', 'stage', 'treatment', 'anno_orig', 'sample_token')
meta_final <- combined_meta[, selInfo]
rownames(meta_final) <- meta_final$cell_id

if (!identical(colnames(combined_counts), rownames(meta_final))) stop('Metadata rownames mismatch count columns')

srt <- CreateSeuratObject(
  counts = combined_counts,
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

report_path <- file.path(tmp_dir, 'gse286827_alignment_report.md')
lines <- c(
  '# gse286827 Alignment Report',
  '',
  paste0('- Generated at: ', format(Sys.time(), '%Y-%m-%d %H:%M:%S')),
  paste0('- Input dir: ', input_dir),
  paste0('- GEO evidence: ', geo_url),
  '',
  '## 1) Series-level evidence',
  '- disease = hnscc',
  '- tissue = tumor',
  '- organ = site-specific head & neck anatomy (from source_name_ch1)',
  '- treatment groups = D / D+T',
  '- timepoint = Pre/Post (with P8_post2 preserved as Post2)',
  '- lib = 10x_5_v2 (from explicit 5′ v2 protocol)',
  '- response = none (no patient-level response label in GEO soft)',
  '- resp_standard = NA (missing-by-design)',
  '- merged mode = immunecells_broad-first; duplicate cell_ids from later subsets collapse into anno_orig only',
  '- orig provenance field = anno_orig',
  '',
  '## 2) Dataset scale',
  paste0('- Unique cells final: ', ncol(combined_counts)),
  paste0('- Genes final: ', nrow(combined_counts)),
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
  paste0('- organ: ', paste(sort(unique(meta_final$organ)), collapse = ', ')),
  paste0('- tissue: ', paste(sort(unique(meta_final$tissue)), collapse = ', ')),
  paste0('- lib: ', paste(sort(unique(meta_final$lib)), collapse = ', ')),
  paste0('- target: ', paste(sort(unique(meta_final$target)), collapse = ', ')),
  paste0('- timepoint: ', paste(sort(unique(meta_final$timepoint)), collapse = ', ')),
  paste0('- response: ', paste(sort(unique(meta_final$response)), collapse = ', ')),
  '',
  '## 6) Evidence source',
  '- GEO soft: https://ftp.ncbi.nlm.nih.gov/geo/series/GSE286nnn/GSE286827/soft/GSE286827_family.soft.gz',
  '- Sample titles encode P#/pre/post/GEX, and extract protocol states Chromium Next GEM Single Cell 5′ v2',
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
