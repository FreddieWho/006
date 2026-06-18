id <- 'gse179994'

library(qs)

input_path <- path.expand('~/006/data/imm/GSE179994/GSE179994_all.scTCR.tsv.gz')
out_dir <- path.expand('~/006/data/processed/tcr/raw')
report_path <- path.expand('~/006/tmp/gse179994_tcr_report.md')
out_path <- file.path(out_dir, 'gse179994_tcr_cellid_barcode.tsv.gz')

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

read_tcr <- function(path) {
  read.delim(path, sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)
}

norm_na <- function(x) {
  x <- as.character(x)
  x[x %in% c('', 'NA', 'NaN')] <- NA_character_
  x
}

collapse_unique <- function(...) {
  vals <- unique(na.omit(unlist(list(...), use.names = FALSE)))
  if (length(vals) == 0) return(NA_character_)
  paste(vals, collapse = ';')
}

any_true <- function(x) {
  x <- toupper(as.character(x))
  any(x == 'TRUE', na.rm = TRUE)
}

extract_tcr_cell_map <- function(path) {
  tcr <- read_tcr(path)
  required_cols <- c('CellName', 'sample', 'clone.id')
  missing <- setdiff(required_cols, colnames(tcr))
  if (length(missing) > 0) stop('Missing required TCR columns: ', paste(missing, collapse = ', '))

  summarize_row <- function(a1, a2, b1, b2, a1id, a2id, b1id, b2id,
                            a1full, a2full, b1full, b2full,
                            av1, av2, aj1, aj2, bv1, bv2, bj1, bj2) {
    alpha_present <- any(c(a1id, a2id) != 'NA' & !is.na(c(a1id, a2id)))
    beta_present <- any(c(b1id, b2id) != 'NA' & !is.na(c(b1id, b2id)))
    alpha_full <- any(c(a1full, a2full) == 'TRUE', na.rm = TRUE)
    beta_full <- any(c(b1full, b2full) == 'TRUE', na.rm = TRUE)
    paired <- alpha_present && beta_present
    full_length <- paired && alpha_full && beta_full
    productive_status <- if (paired && full_length) {
      'paired_full_length'
    } else if (paired) {
      'paired_partial'
    } else if (alpha_present || beta_present) {
      'single_chain'
    } else {
      'none'
    }
    c(
      tra_cdr3 = collapse_unique(a1, a2),
      trb_cdr3 = collapse_unique(b1, b2),
      tra_v_gene = collapse_unique(av1, av2),
      tra_j_gene = collapse_unique(aj1, aj2),
      trb_v_gene = collapse_unique(bv1, bv2),
      trb_j_gene = collapse_unique(bj1, bj2),
      productive_status = productive_status,
      full_length = if (full_length) 'TRUE' else 'FALSE',
      paired = if (paired) 'TRUE' else 'FALSE'
    )
  }

  a1 <- norm_na(tcr[['Identifier(Alpha1)']]); a2 <- norm_na(tcr[['Identifier(Alpha2)']])
  b1 <- norm_na(tcr[['Identifier(Beta1)']]);  b2 <- norm_na(tcr[['Identifier(Beta2)']])
  a1f <- norm_na(tcr[['Full_length(Alpha1)']]); a2f <- norm_na(tcr[['Full_length(Alpha2)']])
  b1f <- norm_na(tcr[['Full_length(Beta1)']]);  b2f <- norm_na(tcr[['Full_length(Beta2)']])

  row_summaries <- mapply(
    summarize_row,
    tcr[['CDR3(Alpha1)']], tcr[['CDR3(Alpha2)']],
    tcr[['CDR3(Beta1)']], tcr[['CDR3(Beta2)']],
    a1, a2, b1, b2,
    a1f, a2f, b1f, b2f,
    tcr[['V_gene(Alpha1)']], tcr[['V_gene(Alpha2)']],
    tcr[['J_gene(Alpha1)']], tcr[['J_gene(Alpha2)']],
    tcr[['V_gene(Beta1)']], tcr[['V_gene(Beta2)']],
    tcr[['J_gene(Beta1)']], tcr[['J_gene(Beta2)']],
    SIMPLIFY = FALSE
  )

  summary_df <- as.data.frame(do.call(rbind, row_summaries), stringsAsFactors = FALSE)

  tcr_map <- data.frame(
    cell_id = tcr[['CellName']],
    barcode = sub('^.*\\.', '', tcr[['CellName']]),
    sample_id = tcr[['sample']],
    clone_id = tcr[['clone.id']],
    summary_df,
    source = basename(path),
    stringsAsFactors = FALSE
  )

  tcr_map <- tcr_map[!duplicated(tcr_map$cell_id), , drop = FALSE]
  tcr_map <- tcr_map[, c('cell_id', 'barcode', 'sample_id', 'clone_id', 'tra_cdr3', 'trb_cdr3',
                         'tra_v_gene', 'tra_j_gene', 'trb_v_gene', 'trb_j_gene',
                         'productive_status', 'full_length', 'paired', 'source')]

  tcr_map
}

tcr_map <- extract_tcr_cell_map(input_path)
out_con <- gzfile(out_path, open = 'wt')
write.table(tcr_map, out_con, sep = '\t', quote = FALSE, row.names = FALSE)
close(out_con)

sample_counts <- as.data.frame(table(tcr_map$sample_id), stringsAsFactors = FALSE)
colnames(sample_counts) <- c('sample_id', 'n_cells')
sample_counts <- sample_counts[order(sample_counts$sample_id), ]

report_lines <- c(
  '# gse179994 TCR Export Report',
  '',
  paste0('- Input: ', input_path),
  paste0('- Output: ', out_path),
  paste0('- Rows: ', nrow(tcr_map)),
  paste0('- Columns: ', ncol(tcr_map)),
  '',
  '## Controlled values',
  paste0('- paired: ', paste(sort(unique(tcr_map$paired)), collapse = ', ')),
  paste0('- full_length: ', paste(sort(unique(tcr_map$full_length)), collapse = ', ')),
  paste0('- productive_status: ', paste(sort(unique(tcr_map$productive_status)), collapse = ', ')),
  '',
  '## Per-sample cell counts',
  '',
  '| sample_id | n_cells |',
  '|---|---:|'
)
report_lines <- c(report_lines, apply(sample_counts, 1, function(x) paste0('| ', x[[1]], ' | ', x[[2]], ' |')))
writeLines(report_lines, report_path)

cat('Done\n')
cat('Rows: ', nrow(tcr_map), '\n')
cat('Report: ', report_path, '\n')
