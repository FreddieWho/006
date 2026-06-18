# Export manual-audit matrices to mtx + cell lists
library(Matrix)

BASE <- "/home/huyudi/006"
TMP_OUT <- "/tmp/p2_manual_audit_export"
dir.create(TMP_OUT, showWarnings = FALSE, recursive = TRUE)

export_krishna_2021_rcc <- function() {
  cid <- "krishna_2021_rcc"
  rds_path <- file.path(BASE, "data", "combo", cid, "ccRCC_6pat_Seurat")
  message("[", cid, "] reading RDS …")
  obj <- readRDS(rds_path)
  message("  converting counts …")
  counts <- obj@raw.data
  if (!inherits(counts, "sparseMatrix")) {
    counts <- as(counts, "sparseMatrix")
  }
  message("  exporting counts …")
  writeMM(counts, file.path(TMP_OUT, paste0(cid, "_counts.mtx")))
  message("  exporting cells …")
  writeLines(colnames(counts), file.path(TMP_OUT, paste0(cid, "_cells.txt")))
  message("  exporting genes …")
  writeLines(rownames(counts), file.path(TMP_OUT, paste0(cid, "_genes.txt")))
  message("  done")
}

export_lambrecht_brca <- function() {
  cid <- "lambrecht_brca"
  src_dir <- file.path(BASE, "data", "combo", cid)
  mats <- c(
    "1863-counts_cells_cohort1.rds",
    "1867-counts_cells_cohort2.rds",
    "1864-counts_tcell_cohort1.rds",
    "1865-counts_myeloid_cohort1.rds",
    "1866-counts_DC_cohort1.rds"
  )
  for (m in mats) {
    mat_path <- file.path(src_dir, m)
    if (!file.exists(mat_path)) next
    subset <- sub("\\.rds$", "", m)
    message("[", cid, "] reading ", subset, " …")
    counts <- readRDS(mat_path)
    message("  dims: ", nrow(counts), " x ", ncol(counts))
    writeMM(counts, file.path(TMP_OUT, paste0(cid, "_", subset, "_counts.mtx")))
    writeLines(colnames(counts), file.path(TMP_OUT, paste0(cid, "_", subset, "_cells.txt")))
    writeLines(rownames(counts), file.path(TMP_OUT, paste0(cid, "_", subset, "_genes.txt")))
    message("  exported")
  }
  message("  done")
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) args <- c("krishna_2021_rcc", "lambrecht_brca")

for (cid in args) {
  if (cid == "krishna_2021_rcc") {
    tryCatch(export_krishna_2021_rcc(), error = function(e) message("ERROR ", cid, ": ", e$message))
  } else if (cid == "lambrecht_brca") {
    tryCatch(export_lambrecht_brca(), error = function(e) message("ERROR ", cid, ": ", e$message))
  } else {
    message("ERROR: unknown ", cid)
  }
}

message("Done.")
