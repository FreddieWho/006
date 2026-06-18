# Convert P2 RDS Seurat objects to h5ad
library(Seurat)
library(zellkonverter)
library(SingleCellExperiment)

BASE <- "/home/huyudi/006"
RAW_OUT <- file.path(BASE, "data", "processed", "srt", "raw")

safe_char <- function(vec) {
  if (is.null(vec)) return(rep(NA_character_, length(vec)))
  as.character(vec)
}

convert_gse225063 <- function() {
  cid <- "GSE225063"
  rds_path <- file.path(BASE, "data", "combo", cid, "GSE225063_18071_cd8_analysis_object_231127.rds.gz")
  out_path <- file.path(RAW_OUT, "gse225063.h5ad")
  message("[", cid, "] decompressing RDS …")
  tmp_rds <- tempfile(fileext = ".rds")
  system(paste("zcat", shQuote(rds_path), ">", shQuote(tmp_rds)))
  message("[", cid, "] reading RDS …")
  obj <- readRDS(tmp_rds)
  unlink(tmp_rds)
  message("  dims: ", nrow(obj), " x ", ncol(obj))
  message("[", cid, "] converting to SCE …")
  sce <- as.SingleCellExperiment(obj, assay = "RNA")
  meta <- as.data.frame(colData(sce))
  meta$cohort_id <- cid
  meta$patient_id <- safe_char(meta$patient_id)
  meta$timepoint_raw <- safe_char(meta$timepoint)
  meta$treatment_regimen_raw <- safe_char(meta$treatment)
  meta$response_raw <- if ("response" %in% colnames(meta)) safe_char(meta$response) else "none"
  meta$celltype_raw <- safe_char(meta$final_clustering_results)
  meta$tissue_source <- safe_char(meta$sample_type)
  colData(sce) <- DataFrame(meta)
  message("[", cid, "] writing h5ad …")
  writeH5AD(sce, out_path, compression = "gzip")
  message("  written ", out_path)
}

convert_gse272993 <- function() {
  cid <- "GSE272993"
  rds_path <- file.path(BASE, "data", "combo", cid, "GSE272993_cd8_nn_labeled_FINAL.RDS")
  out_path <- file.path(RAW_OUT, "gse272993.h5ad")
  message("[", cid, "] reading RDS …")
  obj <- readRDS(rds_path)
  message("  dims: ", nrow(obj), " x ", ncol(obj))
  message("[", cid, "] converting to SCE …")
  sce <- as.SingleCellExperiment(obj, assay = "RNA")
  meta <- as.data.frame(colData(sce))
  meta$cohort_id <- cid
  meta$patient_id <- safe_char(meta$patient_alias)
  meta$timepoint_raw <- safe_char(meta$timepoint)
  meta$treatment_regimen_raw <- safe_char(meta$treatment)
  meta$response_raw <- safe_char(meta$response)
  meta$celltype_raw <- safe_char(meta$gene_cluster_name)
  meta$tissue_source <- "tumor"
  colData(sce) <- DataFrame(meta)
  message("[", cid, "] writing h5ad …")
  writeH5AD(sce, out_path, compression = "gzip")
  message("  written ", out_path)
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) args <- c("GSE225063", "GSE272993")

for (cid in args) {
  if (cid == "GSE225063") {
    tryCatch(convert_gse225063(), error = function(e) message("ERROR ", cid, ": ", e$message))
  } else if (cid == "GSE272993") {
    tryCatch(convert_gse272993(), error = function(e) message("ERROR ", cid, ": ", e$message))
  } else {
    message("ERROR: unknown ", cid)
  }
}

message("Done.")
