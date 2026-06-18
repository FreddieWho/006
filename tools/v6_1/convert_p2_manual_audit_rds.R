# Export manual-audit RDS objects to mtx + csv for Python assembly
library(Matrix)

BASE <- "/home/huyudi/006"
TMP_OUT <- "/tmp/p2_manual_audit_export"
dir.create(TMP_OUT, showWarnings = FALSE, recursive = TRUE)

export_krishna_2021_rcc <- function() {
  cid <- "krishna_2021_rcc"
  rds_path <- file.path(BASE, "data", "combo", cid, "ccRCC_6pat_Seurat")
  message("[", cid, "] reading RDS …")
  obj <- readRDS(rds_path)
  message("  exporting counts …")
  counts <- obj@raw.data
  writeMM(counts, file.path(TMP_OUT, paste0(cid, "_counts.mtx")))
  message("  exporting obs …")
  obs <- data.frame(
    cell = colnames(counts),
    stringsAsFactors = FALSE
  )
  rownames(obs) <- obs$cell
  # bind annotations if present
  ann_path <- file.path(BASE, "data", "combo", cid, "ccRCC_6pat_cell_annotations.txt")
  if (file.exists(ann_path)) {
    ann <- read.delim(ann_path, stringsAsFactors = FALSE)
    rownames(ann) <- ann$cell
    obs <- cbind(obs, ann[rownames(obs), setdiff(colnames(ann), "cell")])
  }
  obs$cohort_id <- cid
  obs$patient_id <- as.character(obs$Sample)
  obs$timepoint_raw <- "none"
  obs$treatment_regimen_raw <- "none"
  obs$response_raw <- "none"
  obs$celltype_raw <- as.character(obs$cluster_name)
  obs$tissue_source <- as.character(obs$type)
  write.csv(obs, file.path(TMP_OUT, paste0(cid, "_obs.csv")))
  message("  exporting var …")
  var <- data.frame(gene_id = rownames(counts))
  write.csv(var, file.path(TMP_OUT, paste0(cid, "_var.csv")))
  message("  done")
}

export_lambrecht_brca <- function() {
  cid <- "lambrecht_brca"
  src_dir <- file.path(BASE, "data", "combo", cid)

  # map matrix files to metadata files
  pairs <- list(
    list(mat = "1863-counts_cells_cohort1.rds", meta = "1872-BIOKEY_metaData_cohort1_web.csv", subset = "cells_cohort1"),
    list(mat = "1867-counts_cells_cohort2.rds", meta = "1871-BIOKEY_metaData_cohort2_web.csv", subset = "cells_cohort2"),
    list(mat = "1864-counts_tcell_cohort1.rds", meta = "1870-BIOKEY_metaData_tcells_cohort1_web.csv", subset = "tcell_cohort1"),
    list(mat = "1865-counts_myeloid_cohort1.rds", meta = "1869-BIOKEY_metaData_myeloid_cohort1_web.csv", subset = "myeloid_cohort1"),
    list(mat = "1866-counts_DC_cohort1.rds", meta = "1868-BIOKEY_metaData_DC_cohort1_web.csv", subset = "DC_cohort1")
  )

  adatas <- list()
  for (p in pairs) {
    mat_path <- file.path(src_dir, p$mat)
    meta_path <- file.path(src_dir, p$meta)
    if (!file.exists(mat_path)) next
    message("[", cid, "] reading ", p$subset, " …")
    counts <- readRDS(mat_path)
    message("  dims: ", nrow(counts), " x ", ncol(counts))

    obs <- NULL
    if (file.exists(meta_path)) {
      meta <- read.csv(meta_path, stringsAsFactors = FALSE)
      rownames(meta) <- meta$Cell
      obs <- meta
    } else {
      obs <- data.frame(cell = colnames(counts), stringsAsFactors = FALSE)
      rownames(obs) <- obs$cell
    }
    obs$cohort_id <- cid
    obs$subset <- p$subset
    obs$patient_id <- as.character(obs$patient_id)
    obs$timepoint_raw <- as.character(obs$timepoint)
    obs$treatment_regimen_raw <- "none"
    obs$response_raw <- "none"
    obs$celltype_raw <- as.character(obs$cellType)
    obs$tissue_source <- "tumor"

    # save per-subset files
    out_prefix <- file.path(TMP_OUT, paste0(cid, "_", p$subset))
    writeMM(counts, paste0(out_prefix, "_counts.mtx"))
    write.csv(obs, paste0(out_prefix, "_obs.csv"))
    var <- data.frame(gene_symbol = rownames(counts))
    write.csv(var, paste0(out_prefix, "_var.csv"))
    message("  exported ", p$subset)
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
