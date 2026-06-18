# Export P2 RDS Seurat objects to mtx + csv for Python h5ad assembly
library(Seurat)
library(Matrix)

BASE <- "/home/huyudi/006"
TMP_OUT <- "/tmp/p2_rds_export"
dir.create(TMP_OUT, showWarnings = FALSE, recursive = TRUE)

export_gse225063 <- function() {
  cid <- "GSE225063"
  rds_path <- file.path(BASE, "data", "combo", cid, "GSE225063_18071_cd8_analysis_object_231127.rds.gz")
  message("[", cid, "] decompressing RDS …")
  tmp_rds <- tempfile(fileext = ".rds")
  system(paste("zcat", shQuote(rds_path), ">", shQuote(tmp_rds)))
  message("[", cid, "] reading RDS …")
  obj <- readRDS(tmp_rds)
  unlink(tmp_rds)
  message("  dims: ", nrow(obj), " x ", ncol(obj))

  message("[", cid, "] exporting counts …")
  counts <- GetAssayData(obj, assay = "RNA", layer = "counts")
  writeMM(counts, file.path(TMP_OUT, paste0(cid, "_counts.mtx")))

  message("[", cid, "] exporting obs …")
  obs <- obj@meta.data
  obs$cohort_id <- cid
  obs$patient_id <- as.character(obs$patient_id)
  obs$timepoint_raw <- as.character(obs$timepoint)
  obs$treatment_regimen_raw <- as.character(obs$treatment)
  obs$response_raw <- if ("response" %in% colnames(obs)) as.character(obs$response) else "none"
  obs$celltype_raw <- as.character(obs$final_clustering_results)
  obs$tissue_source <- as.character(obs$sample_type)
  write.csv(obs, file.path(TMP_OUT, paste0(cid, "_obs.csv.gz")))

  message("[", cid, "] exporting var …")
  var <- data.frame(gene_symbol = rownames(obj))
  write.csv(var, file.path(TMP_OUT, paste0(cid, "_var.csv.gz")))
  message("  done")
}

export_gse272993 <- function() {
  cid <- "GSE272993"
  rds_path <- file.path(BASE, "data", "combo", cid, "GSE272993_cd8_nn_labeled_FINAL.RDS")
  message("[", cid, "] reading RDS …")
  obj <- readRDS(rds_path)
  message("  dims: ", nrow(obj), " x ", ncol(obj))

  message("[", cid, "] exporting counts …")
  counts <- GetAssayData(obj, assay = "RNA", layer = "counts")
  writeMM(counts, file.path(TMP_OUT, paste0(cid, "_counts.mtx")))

  message("[", cid, "] exporting obs …")
  obs <- obj@meta.data
  obs$cohort_id <- cid
  obs$patient_id <- as.character(obs$patient_alias)
  obs$timepoint_raw <- as.character(obs$timepoint)
  obs$treatment_regimen_raw <- as.character(obs$treatment)
  obs$response_raw <- as.character(obs$response)
  obs$celltype_raw <- as.character(obs$gene_cluster_name)
  obs$tissue_source <- "tumor"
  write.csv(obs, file.path(TMP_OUT, paste0(cid, "_obs.csv.gz")))

  message("[", cid, "] exporting var …")
  var <- data.frame(gene_symbol = rownames(obj))
  write.csv(var, file.path(TMP_OUT, paste0(cid, "_var.csv.gz")))
  message("  done")
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) args <- c("GSE225063", "GSE272993")

for (cid in args) {
  if (cid == "GSE225063") {
    tryCatch(export_gse225063(), error = function(e) message("ERROR ", cid, ": ", e$message))
  } else if (cid == "GSE272993") {
    tryCatch(export_gse272993(), error = function(e) message("ERROR ", cid, ": ", e$message))
  } else {
    message("ERROR: unknown ", cid)
  }
}

message("Done.")
