id <- 'lambrecht_hcc'

library(Seurat)
library(anndataR)
source('~/006/scripts/utils/gene_name_align.R')
# dplyr and Matrix are loaded by gene_name_align.R / Seurat

# --- Report helper --------------------------------------------------------
report_lines <- c()
rpt <- function(...) {

  msg <- paste0(...)
  report_lines <<- c(report_lines, msg)
  message(msg)
}

rpt("=== lambrecht_hcc raw preprocessing ===")
rpt(paste0("Date: ", Sys.time()))

# Step 0: Extract zips to temp location ------------------------------------
tmp_dir <- file.path(tempdir(), 'lambrecht_hcc')
dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)

hcc_zip  <- '~/006/data/combo/lambrecht_hcc/7634-HCC_TME.zip'
pbmc_zip <- '~/006/data/combo/lambrecht_hcc/4667-PBMC (1).zip'

unzip(hcc_zip,  exdir = tmp_dir, overwrite = TRUE)
unzip(pbmc_zip, exdir = tmp_dir, overwrite = TRUE)

# Step 1: Read HCC 10x data -----------------------------------------------
# features.tsv has only 1 column (gene symbol); use feature.column = 1 (G6)
hcc_mtx_dir <- file.path(tmp_dir,
  'HCC_data_VIB_CCB_portal/processed_feature_bc_matrix')

hcc_cts <- ReadMtx(
  mtx      = file.path(hcc_mtx_dir, 'matrix.mtx.gz'),
  cells    = file.path(hcc_mtx_dir, 'barcodes.tsv.gz'),
  features = file.path(hcc_mtx_dir, 'features.tsv.gz'),
  feature.column = 1
)

rpt(paste0("HCC raw dims: ", nrow(hcc_cts), " genes x ", ncol(hcc_cts), " cells"))

# Step 2: Read PBMC 10x data ----------------------------------------------
pbmc_mtx_dir <- file.path(tmp_dir, 'PBMC/PBMC_scRNAseq_NKT')

pbmc_cts <- ReadMtx(
  mtx      = file.path(pbmc_mtx_dir, 'matrix.mtx'),
  cells    = file.path(pbmc_mtx_dir, 'barcodes.tsv'),
  features = file.path(pbmc_mtx_dir, 'features.tsv'),
  feature.column = 1
)

rpt(paste0("PBMC raw dims: ", nrow(pbmc_cts), " genes x ", ncol(pbmc_cts), " cells"))

# Step 3: Read metadata ----------------------------------------------------
hcc_meta <- read.delim(
  file.path(tmp_dir, 'HCC_data_VIB_CCB_portal/meta_data.txt'),
  row.names = 1, check.names = FALSE, stringsAsFactors = FALSE
)

pbmc_meta <- read.delim(
  file.path(tmp_dir, 'PBMC/PBMC_meta_data.txt'),
  row.names = 1, check.names = FALSE, stringsAsFactors = FALSE
)

rpt(paste0("HCC meta rows: ", nrow(hcc_meta)))
rpt(paste0("PBMC meta rows: ", nrow(pbmc_meta)))

# G1 input gate
stopifnot(nrow(hcc_cts) > 0, ncol(hcc_cts) > 0)
stopifnot(nrow(pbmc_cts) > 0, ncol(pbmc_cts) > 0)

# Step 4: Align barcodes to metadata --------------------------------------
# HCC: exact match expected
hcc_meta <- hcc_meta[colnames(hcc_cts), , drop = FALSE]
stopifnot(!any(is.na(rownames(hcc_meta))))
stopifnot(all(colnames(hcc_cts) == rownames(hcc_meta)))

# PBMC: 250 barcodes without metadata → keep only cells with metadata
pbmc_common  <- intersect(colnames(pbmc_cts), rownames(pbmc_meta))
pbmc_dropped <- ncol(pbmc_cts) - length(pbmc_common)
rpt(paste0("PBMC cells dropped (no metadata): ", pbmc_dropped))
pbmc_cts  <- pbmc_cts[, pbmc_common]
pbmc_meta <- pbmc_meta[pbmc_common, , drop = FALSE]
stopifnot(all(colnames(pbmc_cts) == rownames(pbmc_meta)))

# Step 5: Verify same gene set and merge matrices --------------------------
stopifnot(identical(rownames(hcc_cts), rownames(pbmc_cts)))
cts <- cbind(hcc_cts, pbmc_cts)
rpt(paste0("Combined raw dims: ", nrow(cts), " genes x ", ncol(cts), " cells"))

# G7 dimnames integrity after merge
stopifnot(!is.null(rownames(cts)), !is.null(colnames(cts)))
stopifnot(length(rownames(cts)) == nrow(cts))
stopifnot(length(colnames(cts)) == ncol(cts))
stopifnot(!any(duplicated(colnames(cts))))

# Step 6: Build unified metadata -------------------------------------------
# Vocabulary maps (registry-aligned)
target_map <- c(
  'atezo+bev'  = 'PD1+VEGF',
  'atezo+cabo' = 'PD1+TKI',
  'ICP'        = 'PD1',
  'TKI'        = 'TKI',
  'None'       = 'none'
)

response_map <- c(
  'Responder'           = 'R',
  'NonResponder'        = 'NR',
  'DeathBeforeImaging'  = 'DBI'
)

tp_map <- c('W0' = 'pre', 'W3' = 'post', 'W6' = 'post')

# --- HCC ---
hcc_m <- data.frame(
  cell_id        = paste0(id, '_', rownames(hcc_meta)),
  sample_id      = paste0(id, '_', hcc_meta$Tumour_ID),
  donor_id       = paste0(id, '_', hcc_meta$Patient_ID),
  disease        = 'hcc',
  organ          = 'liver',
  tissue         = 'tumor',
  lib            = '10x',
  target         = unname(target_map[hcc_meta$Treatment]),
  timepoint      = unname(tp_map[hcc_meta$Week]),
  response       = unname(response_map[hcc_meta$Response]),
  resp_standard  = 'RECIST',
  week_raw       = hcc_meta$Week,
  treatment_raw  = hcc_meta$Treatment,
  stringsAsFactors = FALSE
)
rownames(hcc_m) <- hcc_m$cell_id

# --- PBMC ---
pbmc_m <- data.frame(
  cell_id        = paste0(id, '_', rownames(pbmc_meta)),
  sample_id      = paste0(id, '_', pbmc_meta$Sample_ID),
  donor_id       = paste0(id, '_', pbmc_meta$Patient_ID),
  disease        = 'hcc',
  organ          = 'pbmc',
  tissue         = 'pbmc',
  lib            = '10x',
  target         = unname(target_map[pbmc_meta$Treatment]),
  timepoint      = unname(tp_map[pbmc_meta$Week]),
  response       = unname(response_map[pbmc_meta$Response]),
  resp_standard  = 'RECIST',
  week_raw       = pbmc_meta$Week,
  treatment_raw  = pbmc_meta$Treatment,
  stringsAsFactors = FALSE
)
rownames(pbmc_m) <- pbmc_m$cell_id

# Combine metadata
meta <- rbind(hcc_m, pbmc_m)

# Prefix cell IDs on count matrix
colnames(cts) <- paste0(id, '_', colnames(cts))

# G4 ID integrity
stopifnot(all(colnames(cts) == meta$cell_id))
stopifnot(all(colnames(cts) == rownames(meta)))
stopifnot(!any(duplicated(meta$cell_id)))
stopifnot(all(grepl(paste0('^', id, '_'), meta$cell_id)))

rpt(paste0("Unified meta rows: ", nrow(meta)))

# Step 7: Select final metadata columns ------------------------------------
selInfo <- c('cell_id', 'sample_id', 'donor_id', 'disease', 'organ', 'tissue',
             'lib', 'target', 'timepoint', 'response', 'resp_standard',
             'week_raw', 'treatment_raw')
meta_final <- meta[, selInfo]

# G3 metadata gate
stopifnot(nrow(meta_final) == ncol(cts))
req_fields <- c('cell_id', 'sample_id', 'donor_id', 'disease', 'organ', 'tissue', 'lib')
for (f in req_fields) {
  stopifnot(sum(is.na(meta_final[[f]])) == 0)
}

# G8 required-vs-conditional missingness
rpt("--- G8 required field NA counts ---")
for (f in req_fields) {
  rpt(paste0("  ", f, ": ", sum(is.na(meta_final[[f]])), " NA"))
}
cond_fields <- c('target', 'timepoint', 'response', 'resp_standard',
                 'week_raw', 'treatment_raw')
for (f in cond_fields) {
  rpt(paste0("  ", f, ": ", sum(is.na(meta_final[[f]])), " NA"))
}

# Step 8: Gene symbol alignment --------------------------------------------
genes <- rownames(cts)
res <- gene_symbol_align(genes, gene_anno_list, alias_map)

n_success <- sum(res$conversion_status == 'success')
rpt(paste0("Gene alignment: ", n_success, "/", nrow(res), " success"))

# Filter biotypes: protein_coding, IG_C_gene, TR_C_gene
res_type_flt <- subset(res, gene_type %in% c('protein_coding', 'IG_C_gene', 'TR_C_gene'))
rpt(paste0("After biotype filter: ", nrow(res_type_flt), " genes"))

# G6 post-filter check
stopifnot(nrow(res_type_flt) > 0)

# Resolve duplicated aligned symbols with max rowMeans
dup_gene <- unique(res_type_flt$aligned_symbol[duplicated(res_type_flt$aligned_symbol)])

res_undup <- subset(res_type_flt, !(aligned_symbol %in% dup_gene))

res_dup <- subset(res_type_flt, aligned_symbol %in% dup_gene)
if (nrow(res_dup) > 0) {
  res_dup$meanCts <- cts[match(res_dup$original_id, genes), ] %>% rowMeans()
  res_dup %>%
    group_by(aligned_symbol) %>%
    slice_max(meanCts, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    select(-meanCts) %>%
    as.data.frame() -> res_dup
}

res_final <- rbind(res_dup, res_undup)
rpt(paste0("After dedup: ", nrow(res_final), " unique genes"))

# G2 gene gate
stopifnot(!any(duplicated(res_final$aligned_symbol)))

# Align count matrix
cts_alignGene <- cts[match(res_final$original_id, genes), ]
rownames(cts_alignGene) <- res_final$aligned_symbol

# G7 post-alignment integrity
stopifnot(!is.null(rownames(cts_alignGene)), !is.null(colnames(cts_alignGene)))
stopifnot(length(rownames(cts_alignGene)) == nrow(cts_alignGene))
stopifnot(length(colnames(cts_alignGene)) == ncol(cts_alignGene))
stopifnot(!any(duplicated(rownames(cts_alignGene))))
stopifnot(!any(duplicated(colnames(cts_alignGene))))
stopifnot(all(colnames(cts_alignGene) == rownames(meta_final)))

rpt(paste0("Final aligned dims: ", nrow(cts_alignGene), " genes x ",
           ncol(cts_alignGene), " cells"))

# Step 9: Create Seurat object ---------------------------------------------
srt <- CreateSeuratObject(
  cts_alignGene,
  meta.data  = meta_final,
  min.cells  = 0,
  min.features = 0,
  project    = id
)

rpt(paste0("Seurat object: ", ncol(srt), " cells x ", nrow(srt), " genes"))

# Step 10: Save outputs ----------------------------------------------------
out_qs   <- paste0('~/006/data/processed/srt/raw/', id, '.qs')
out_h5ad <- paste0('~/006/data/processed/srt/raw/', id, '.h5ad')

qs::qsave(srt, out_qs)
rpt(paste0("Saved: ", out_qs))

write_h5ad(srt, out_h5ad)
rpt(paste0("Saved: ", out_h5ad))

# G5 output gate
stopifnot(file.exists(normalizePath(out_qs, mustWork = FALSE)))
stopifnot(file.info(normalizePath(out_qs, mustWork = FALSE))$size > 0)
stopifnot(file.exists(normalizePath(out_h5ad, mustWork = FALSE)))
stopifnot(file.info(normalizePath(out_h5ad, mustWork = FALSE))$size > 0)
rpt("G5 output gate: PASS")

# Step 11: Reload verification ---------------------------------------------
srt_reload <- qs::qread(out_qs)
stopifnot(ncol(srt_reload) == ncol(srt))
stopifnot(nrow(srt_reload) == nrow(srt))
rpt("QS reload check: PASS")

# Step 12: Full audit report -----------------------------------------------
rpt("")
rpt("=== AUDIT REPORT ===")
rpt(paste0("dataset_id: ", id))
rpt(paste0("raw_dims: ", nrow(cts), " genes x ", ncol(cts), " cells"))
rpt(paste0("aligned_dims: ", nrow(cts_alignGene), " genes x ",
           ncol(cts_alignGene), " cells"))

rpt("")
rpt("--- Per-source cell counts ---")
rpt(paste0("  HCC TME: ", ncol(hcc_cts), " cells"))
rpt(paste0("  PBMC:    ", ncol(pbmc_cts), " cells"))

rpt("")
rpt("--- Per-sample cell counts ---")
sample_counts <- sort(table(meta_final$sample_id))
for (s in names(sample_counts)) {
  rpt(paste0("  ", s, ": ", sample_counts[s]))
}

rpt("")
rpt("--- Controlled vocabulary snapshots ---")
rpt(paste0("  organ:          ", paste(sort(unique(meta_final$organ)), collapse = ", ")))
rpt(paste0("  tissue:         ", paste(sort(unique(meta_final$tissue)), collapse = ", ")))
rpt(paste0("  lib:            ", paste(sort(unique(meta_final$lib)), collapse = ", ")))
rpt(paste0("  target:         ", paste(sort(unique(meta_final$target)), collapse = ", ")))
rpt(paste0("  timepoint:      ", paste(sort(unique(meta_final$timepoint)), collapse = ", ")))
rpt(paste0("  response:       ", paste(sort(unique(na.omit(meta_final$response))), collapse = ", ")))
rpt(paste0("  resp_standard:  ", paste(sort(unique(meta_final$resp_standard)), collapse = ", ")))

rpt("")
rpt("--- ID integrity ---")
rpt(paste0("  cell_id unique:          ", !any(duplicated(meta_final$cell_id))))
rpt(paste0("  prefix consistent:       ", all(grepl(paste0('^', id, '_'), meta_final$cell_id))))
rpt(paste0("  rownames == cell_id:     ", all(rownames(meta_final) == meta_final$cell_id)))
rpt(paste0("  colnames == rownames:    ", all(colnames(cts_alignGene) == rownames(meta_final))))

rpt("")
rpt("--- Output files ---")
rpt(paste0("  ", out_qs,   " exists: ", file.exists(normalizePath(out_qs, mustWork = FALSE))))
rpt(paste0("  ", out_h5ad, " exists: ", file.exists(normalizePath(out_h5ad, mustWork = FALSE))))

rpt("")
rpt("--- Mapping policy (user-confirmed defaults) ---")
rpt("  target_map:   atezo+bev->PD1+VEGF, atezo+cabo->PD1+TKI, ICP->PD1, TKI->TKI, None->none")
rpt("  response_map: Responder->R, NonResponder->NR, DeathBeforeImaging->DBI")
rpt("  timepoint_map: W0->pre, W3->post, W6->post")
rpt("  tissue_map:   HCC->tumor, PBMC->pbmc")
rpt("  organ_map:    HCC->liver, PBMC->pbmc")
rpt("  lib:          10x (exact version unknown; coarse label per policy)")
rpt("  resp_standard: RECIST (from registry)")
rpt("  features.tsv: 1-column (gene symbol only); used column 1 per G6")
rpt("  PBMC barcode-meta mismatch: 250 cells dropped (no metadata)")
rpt("  Evidence: local zip metadata + docs/data_collection.csv registry")
rpt("  Paper: https://www.nature.com/articles/s41467-023-43381-1")

# Write report to tmp/
report_dir <- normalizePath('~/006/tmp', mustWork = FALSE)
dir.create(report_dir, recursive = TRUE, showWarnings = FALSE)
report_path <- file.path(report_dir, paste0(id, '_alignment_report.txt'))
writeLines(report_lines, report_path)
rpt(paste0("Report written to: ", report_path))

message("=== Done! ===")
