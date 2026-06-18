id <- 'gse200996'

library(dplyr)
library(qs)

input_dir <- '~/006/data/combo/GSE200996'
out_dir <- '~/006/data/processed/tcr/raw'
tmp_dir <- '~/006/tmp'

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

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

parse_sample_label <- function(path) {
  nm <- basename(path)
  nm <- sub('^GSM[0-9]+_', '', nm)
  nm <- sub('^filtered_contig_annotations_', '', nm)
  nm <- sub('_TCR_sc_.*$', '', nm)
  nm
}

read_tcr_file <- function(path) {
  df <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  df$sample_label <- parse_sample_label(path)
  df$source_file <- basename(path)
  df
}

files <- list.files(input_dir, pattern = 'filtered_contig_annotations_.*\\.csv\\.gz$', full.names = TRUE)
if (length(files) == 0) stop('No contig annotation files found in ', input_dir)

all_tcr <- bind_rows(lapply(sort(files), read_tcr_file))

all_tcr$barcode <- as.character(all_tcr$barcode)
all_tcr$cell_id <- paste0(id, '_', all_tcr$sample_label, '_', all_tcr$barcode)

summary_df <- all_tcr %>%
  mutate(
    chain = as.character(chain),
    v_gene = as.character(v_gene),
    j_gene = as.character(j_gene),
    cdr3 = as.character(cdr3),
    cdr3_nt = as.character(cdr3_nt),
    raw_clonotype_id = as.character(raw_clonotype_id),
    full_length = as.character(full_length),
    productive = as.character(productive),
    reads = as.numeric(reads),
    umis = as.numeric(umis)
  ) %>%
  group_by(cell_id) %>%
  summarise(
    barcode = first_nonempty(barcode),
    sample_id = first_nonempty(sample_label),
    clone_id = merge_label(raw_clonotype_id),
    tra_cdr3 = merge_label(cdr3[chain == 'TRA']),
    trb_cdr3 = merge_label(cdr3[chain == 'TRB']),
    tra_v_gene = merge_label(v_gene[chain == 'TRA']),
    tra_j_gene = merge_label(j_gene[chain == 'TRA']),
    trb_v_gene = merge_label(v_gene[chain == 'TRB']),
    trb_j_gene = merge_label(j_gene[chain == 'TRB']),
    tra_full = any(full_length[chain == 'TRA'] == 'True', na.rm = TRUE),
    trb_full = any(full_length[chain == 'TRB'] == 'True', na.rm = TRUE),
    tra_present = any(chain == 'TRA', na.rm = TRUE),
    trb_present = any(chain == 'TRB', na.rm = TRUE),
    paired = tra_present & trb_present,
    full_length = paired & tra_full & trb_full,
    productive_status = dplyr::case_when(
      paired & full_length ~ 'paired_full_length',
      paired ~ 'paired_partial',
      tra_present | trb_present ~ 'single_chain',
      TRUE ~ 'none'
    ),
    source = merge_label(source_file),
    .groups = 'drop'
  )

summary_df <- as.data.frame(summary_df, stringsAsFactors = FALSE)
rownames(summary_df) <- summary_df$cell_id

out_path <- file.path(out_dir, 'gse200996_tcr_cellid_barcode.tsv.gz')
out_con <- gzfile(out_path, open = 'wt')
write.table(
  summary_df[, c('cell_id', 'barcode', 'sample_id', 'clone_id', 'tra_cdr3', 'trb_cdr3', 'tra_v_gene', 'tra_j_gene', 'trb_v_gene', 'trb_j_gene', 'productive_status', 'full_length', 'paired', 'source')],
  out_con,
  sep = '\t', quote = FALSE, row.names = FALSE
)
close(out_con)

sample_counts <- as.data.frame(table(summary_df$sample_id), stringsAsFactors = FALSE)
colnames(sample_counts) <- c('sample_id', 'n_cells')
sample_counts <- sample_counts[order(sample_counts$sample_id), ]

report_path <- file.path(tmp_dir, 'gse200996_tcr_report.md')
lines <- c(
  '# gse200996 TCR Export Report',
  '',
  paste0('- Input files: ', length(files)),
  paste0('- Rows: ', nrow(summary_df)),
  paste0('- Columns: ', ncol(summary_df)),
  '',
  '## Controlled values',
  paste0('- paired: ', paste(sort(unique(summary_df$paired)), collapse = ', ')),
  paste0('- full_length: ', paste(sort(unique(summary_df$full_length)), collapse = ', ')),
  paste0('- productive_status: ', paste(sort(unique(summary_df$productive_status)), collapse = ', ')),
  '',
  '## Per-sample cell counts',
  '',
  '| sample_id | n_cells |',
  '|---|---:|'
)
lines <- c(lines, apply(sample_counts, 1, function(x) paste0('| ', x[[1]], ' | ', x[[2]], ' |')))
writeLines(lines, report_path)

cat('Done\n')
cat('Rows: ', nrow(summary_df), '\n')
cat('Report: ', report_path, '\n')
