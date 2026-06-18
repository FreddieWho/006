id <- 'gse205506'

library(Seurat)
library(anndataR)
library(Matrix)
library(data.table)
library(dplyr)
library(stringr)
library(readxl)
library(tidyr)
library(qs)

gene_anno_list <- readRDS('~/006/data/ref/gene_anno_list.rds')
alias_map <- readRDS('~/006/data/ref/alias_map.rds')

detect_id_type <- function(ids) {
  sapply(ids, function(id) {
    if (grepl('^ENS[A-Z]*G[0-9]+$', id)) {
      'ensembl'
    } else if (grepl('^[0-9]+$', id)) {
      'entrez'
    } else {
      'symbol'
    }
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

  if (length(ent_idx) > 0) {
    symbol_out[ent_idx] <- NA_character_
  }

  if (length(sym_idx) > 0) {
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

  symbol_out
}

get_gene_type <- function(symbols, gene_anno_list) {
  gene_v19 <- gene_anno_list$v19
  gene_v38 <- gene_anno_list$v38
  idx38 <- match(symbols, gene_v38$gene_name)
  type38 <- gene_v38$gene_type[idx38]
  miss <- is.na(type38)
  idx19 <- match(symbols[miss], gene_v19$gene_name)
  type19 <- gene_v19$gene_type[idx19]
  type38[miss] <- type19
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

input_dir <- '~/006/data/combo/GSE205506'
out_dir <- '~/006/data/processed/srt/raw'
tmp_dir <- '~/006/tmp'

matrix_files <- Sys.glob(file.path(input_dir, '*_matrix.mtx.gz'))
if (length(matrix_files) == 0) stop('No *_matrix.mtx.gz files found in input_dir')

sample_tokens <- gsub('_matrix.mtx.gz$', '', basename(matrix_files))
sample_tokens <- sort(sample_tokens)

get_row <- function(lines, key) {
  idx <- which(startsWith(lines, key))
  if (length(idx) == 0) return(NULL)
  strsplit(lines[idx[1]], '\t')[[1]][-1]
}

parse_geo_matrix <- function(url) {
  con <- gzcon(url(url, open = 'rb'))
  on.exit(close(con), add = TRUE)
  lines <- readLines(con, warn = FALSE)

  gsm <- get_row(lines, '!Sample_geo_accession')
  title <- get_row(lines, '!Sample_title')
  if (is.null(gsm) || is.null(title)) stop('Failed to parse GEO matrix metadata: ', url)

  gsm <- gsub('"', '', gsm)
  title <- gsub('"', '', title)

  m <- stringr::str_match(title, '^([0-9]+):\\s*(P[0-9]+),\\s*([^,]+),\\s*([^,]+),\\s*([^,]+),\\s*scRNA-seq$')
  if (any(is.na(m[, 1]))) {
    bad <- title[is.na(m[, 1])]
    stop('Unparsed sample titles in GEO metadata: ', paste(bad, collapse = ' | '))
  }

  data.frame(
    gsm_id = gsm,
    sample_code = m[, 2],
    patient_id = m[, 3],
    site_raw = m[, 4],
    genotype_raw = tolower(m[, 5]),
    treatment_raw = m[, 6],
    stringsAsFactors = FALSE
  )
}

geo_1 <- parse_geo_matrix('https://ftp.ncbi.nlm.nih.gov/geo/series/GSE205nnn/GSE205506/matrix/GSE205506-GPL24676_series_matrix.txt.gz')
geo_2 <- parse_geo_matrix('https://ftp.ncbi.nlm.nih.gov/geo/series/GSE205nnn/GSE205506/matrix/GSE205506-GPL29480_series_matrix.txt.gz')
geo_map <- bind_rows(geo_1, geo_2)

if (nrow(geo_map) != 40) stop('Expected 40 GEO samples, got: ', nrow(geo_map))
if (any(duplicated(geo_map$gsm_id))) stop('Duplicated GSM IDs in GEO metadata')

meta_xlsx <- read_xlsx(file.path(input_dir, 'GSE205506.xlsx'))
colnames(meta_xlsx)[9:12] <- c('Tumor_pre', 'Tumor_post', 'Normal_post', 'Treatment_response')
meta_xlsx <- meta_xlsx[!is.na(meta_xlsx$`Patient ID`), ]
meta_xlsx <- tidyr::fill(meta_xlsx, Treatment, .direction = 'down')

patient_meta <- meta_xlsx %>%
  transmute(
    patient_id = `Patient ID`,
    treatment_group = Treatment,
    location = `Tumor anatomical location`,
    mismatch_repair = `Mismatch repair defective protein`,
    microsatellite_status = `Microsatellite status`,
    stage = `Stage(cTMN)`,
    sex = Sex,
    age = Age,
    treatment_response_raw = Treatment_response
  )

if (any(duplicated(patient_meta$patient_id))) stop('Duplicated patient IDs in xlsx metadata')

token_df <- data.frame(sample_token = sample_tokens, stringsAsFactors = FALSE) %>%
  mutate(
    gsm_id = stringr::str_extract(sample_token, '^GSM[0-9]+'),
    version = stringr::str_extract(sample_token, 'v[0-9]+$')
  )

sample_meta <- token_df %>%
  left_join(geo_map, by = 'gsm_id') %>%
  left_join(patient_meta, by = 'patient_id') %>%
  mutate(
    sample_id = paste0(id, '_', sample_token),
    donor_id = paste0(id, '_', patient_id),
    disease = 'crc',
    organ = ifelse(str_detect(tolower(site_raw), 'rectum'), 'rectum', 'colon'),
    tissue = ifelse(genotype_raw == 'normal', 'adj', 'tumor'),
    lib = case_when(
      version == 'v2' ~ '10x_5_v2',
      version == 'v3' ~ '10x_5_v3',
      TRUE ~ '10x_5'
    ),
    target = case_when(
      str_detect(tolower(treatment_raw), 'untreated') ~ 'none',
      str_detect(tolower(treatment_raw), 'celecoxib') ~ 'PD1;COX2',
      str_detect(tolower(treatment_raw), 'anti-pd-1') ~ 'PD1',
      TRUE ~ 'none'
    ),
    timepoint = ifelse(str_detect(tolower(treatment_raw), 'untreated'), 'Pre', 'Post'),
    response = ifelse(treatment_response_raw == 'pCR', 'PR', 'NE'),
    resp_standard = paste0('PATHOLOGICAL_RESPONSE:', treatment_response_raw),
    treatment = treatment_raw,
    anno_orig = NA_character_
  )

required_map <- c('sample_token', 'sample_id', 'donor_id', 'patient_id', 'disease', 'organ', 'tissue', 'lib',
                  'target', 'timepoint', 'response', 'resp_standard', 'sex', 'age', 'stage', 'treatment')
missing_required <- setdiff(required_map, colnames(sample_meta))
if (length(missing_required) > 0) stop('Missing required metadata columns: ', paste(missing_required, collapse = ','))

if (any(is.na(sample_meta$patient_id))) {
  bad <- sample_meta$sample_token[is.na(sample_meta$patient_id)]
  stop('Unmapped GSM to patient metadata for samples: ', paste(bad, collapse = ', '))
}

mat_list <- vector('list', length(sample_tokens))
names(mat_list) <- sample_tokens
feature_ref <- NULL

for (tok in sample_tokens) {
  mtx_f <- file.path(input_dir, paste0(tok, '_matrix.mtx.gz'))
  bar_f <- file.path(input_dir, paste0(tok, '_barcodes.tsv.gz'))
  fea_f <- file.path(input_dir, paste0(tok, '_features.tsv.gz'))

  cts <- Matrix::readMM(mtx_f)
  barcodes <- fread(bar_f, header = FALSE)[[1]]
  fea_dt <- fread(fea_f, header = FALSE)
  features <- if (ncol(fea_dt) >= 2) fea_dt[[2]] else fea_dt[[1]]

  if (length(features) != nrow(cts)) stop('Feature length mismatch in sample: ', tok)
  if (length(barcodes) != ncol(cts)) stop('Barcode length mismatch in sample: ', tok)

  rownames(cts) <- features
  colnames(cts) <- paste0(id, '__', tok, '__', gsub('\\.', '_', barcodes))

  if (is.null(feature_ref)) {
    feature_ref <- rownames(cts)
  } else {
    if (!identical(feature_ref, rownames(cts))) stop('Feature set/order differs across samples at: ', tok)
  }

  mat_list[[tok]] <- cts
}

cts_all <- do.call(cbind, mat_list)
all_cell_ids <- unlist(lapply(mat_list, colnames), use.names = FALSE)
dimnames(cts_all) <- list(feature_ref, all_cell_ids)

if (is.null(rownames(cts_all)) || is.null(colnames(cts_all))) {
  stop('Combined raw matrix lost dimnames during cbind')
}

genes <- rownames(cts_all)
res <- gene_symbol_align(genes, gene_anno_list, alias_map)
res_type_flt <- subset(res, gene_type %in% c('protein_coding', 'IG_C_gene', 'TR_C_gene'))

dup_gene <- res_type_flt[duplicated(res_type_flt$aligned_symbol), 'aligned_symbol']
res_undup <- subset(res_type_flt, !(aligned_symbol %in% dup_gene))
res_dup <- subset(res_type_flt, aligned_symbol %in% dup_gene)
res_dup$meanCts <- cts_all[match(res_dup$original_id, genes), ] %>% rowMeans()
res_dup <- res_dup %>%
  group_by(aligned_symbol) %>%
  slice_max(meanCts, n = 1, with_ties = FALSE) %>%
  select(-meanCts) %>%
  as.data.frame()
res_dup <- subset(res_dup, !(original_id %in% c('SPANXB1', 'XAGE2')))

res_final <- rbind(res_dup, res_undup)
sel_idx <- match(res_final$original_id, genes)
aligned_genes <- res_final$aligned_symbol
aligned_cells <- colnames(cts_all)
cts_alignGene <- cts_all[sel_idx, , drop = FALSE]
cts_alignGene <- Matrix::Matrix(cts_alignGene, sparse = TRUE)
dimnames(cts_alignGene) <- list(aligned_genes, aligned_cells)

if (is.null(rownames(cts_alignGene)) || is.null(colnames(cts_alignGene))) {
  cat('DEBUG dimnames lost after alignment\n')
  cat('class:', class(cts_alignGene), '\n')
  cat('rownames NULL:', is.null(rownames(cts_alignGene)), ' colnames NULL:', is.null(colnames(cts_alignGene)), '\n')
  cat('aligned_genes length:', length(aligned_genes), ' aligned_cells length:', length(aligned_cells), '\n')
  dimnames(cts_alignGene) <- list(aligned_genes, aligned_cells)
  if (is.null(rownames(cts_alignGene)) || is.null(colnames(cts_alignGene))) {
    stop('Count matrix dimnames are NULL after alignment')
  }
}
if (any(is.na(rownames(cts_alignGene))) || any(rownames(cts_alignGene) == '')) {
  stop('Aligned gene names contain NA/empty values')
}
if (any(duplicated(rownames(cts_alignGene)))) {
  stop('Aligned gene names still duplicated after deduplication')
}

cell_meta <- data.frame(cell_id = colnames(cts_alignGene), stringsAsFactors = FALSE) %>%
  mutate(
    sample_token = stringr::str_match(cell_id, paste0('^', id, '__([^_]+(?:_[^_]+)*)__'))[, 2]
  )

if (any(is.na(cell_meta$sample_token))) stop('Failed to parse sample_token from cell_id for some cells')

meta <- cell_meta %>%
  left_join(sample_meta, by = 'sample_token')

if (any(is.na(meta$sample_id))) stop('Cell metadata join failed for some cells')

selInfo <- c(
  'cell_id', 'sample_id', 'donor_id', 'disease', 'organ', 'tissue', 'lib',
  'target', 'timepoint', 'response', 'resp_standard', 'sex', 'age', 'stage', 'treatment', 'anno_orig'
)
meta_final <- meta[, selInfo]
rownames(meta_final) <- meta_final$cell_id

if (!identical(colnames(cts_alignGene), rownames(meta_final))) {
  stop('Metadata rownames do not match count matrix columns')
}

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

sample_cell_counts <- as.data.frame(table(meta_final$sample_id), stringsAsFactors = FALSE)
colnames(sample_cell_counts) <- c('sample_id', 'n_cells')
sample_cell_counts <- sample_cell_counts[order(sample_cell_counts$sample_id), ]

report_path <- file.path(tmp_dir, 'gse205506_alignment_report.md')
lines <- c(
  '# gse205506 Alignment Report',
  '',
  paste0('- Generated at: ', format(Sys.time(), '%Y-%m-%d %H:%M:%S')),
  paste0('- Input dir: ', input_dir),
  paste0('- Output qs: ', qs_path),
  paste0('- Output h5ad: ', h5ad_path),
  '',
  '## 1) Dataset scale',
  paste0('- Samples detected: ', length(sample_tokens)),
  paste0('- Cells (raw): ', ncol(cts_all)),
  paste0('- Genes (raw): ', nrow(cts_all)),
  paste0('- Genes (aligned+filtered): ', nrow(cts_alignGene)),
  paste0('- Cells (final): ', ncol(cts_alignGene)),
  '',
  '## 2) Mapping policy used (user confirmed)',
  '- disease = crc',
  '- timepoint: untreated -> Pre; treated -> Post',
  '- target: anti-PD-1 -> PD1; anti-PD-1+celecoxib -> PD1;COX2; untreated -> none',
  '- tissue: genotype normal -> adj; tumor -> tumor',
  '- organ: keep original site (colon/rectum)',
  '- response: pCR -> PR; non-pCR -> NE',
  '- resp_standard: PATHOLOGICAL_RESPONSE:<raw pCR/non-pCR>',
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
  paste0('- organ: ', paste(sort(unique(meta_final$organ)), collapse = ', ')),
  paste0('- tissue: ', paste(sort(unique(meta_final$tissue)), collapse = ', ')),
  paste0('- lib: ', paste(sort(unique(meta_final$lib)), collapse = ', ')),
  paste0('- target: ', paste(sort(unique(meta_final$target)), collapse = ', ')),
  paste0('- timepoint: ', paste(sort(unique(meta_final$timepoint)), collapse = ', ')),
  paste0('- response: ', paste(sort(unique(meta_final$response)), collapse = ', ')),
  '',
  '## 5) Per-sample cell counts',
  '',
  '| sample_id | n_cells |',
  '|---|---:|'
)

tab_lines <- apply(sample_cell_counts, 1, function(x) paste0('| ', x[[1]], ' | ', x[[2]], ' |'))
lines <- c(lines, tab_lines)

writeLines(lines, report_path)

cat('Done\n')
cat('Report: ', report_path, '\n')
