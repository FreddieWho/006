# Reproducible conversion path for Seurat / RDS-backed cohorts into AnnData h5ad.
#
# Usage:
#   Rscript rds_to_h5ad.R <cohort_id> <rds_path> <out_h5ad_path> [assay_name]
#
# Arguments:
#   cohort_id     - cohort identifier (written into obs$cohort_id)
#   rds_path      - path to the Seurat / RDS object
#   out_h5ad_path - output h5ad path
#   assay_name    - assay to extract (default: "RNA")
#
# Dependencies: Seurat, zellkonverter, SingleCellExperiment

library(Seurat)
library(zellkonverter)
library(SingleCellExperiment)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: Rscript rds_to_h5ad.R <cohort_id> <rds_path> <out_h5ad_path> [assay_name]")
}

cid        <- args[1]
rds_path   <- args[2]
out_path   <- args[3]
assay_name <- if (length(args) >= 4) args[4] else "RNA"

if (!file.exists(rds_path)) {
  stop(paste("RDS file not found:", rds_path))
}

message("[", cid, "] reading RDS …")
obj <- readRDS(rds_path)

message("[", cid, "] converting to SingleCellExperiment …")
if (inherits(obj, "Seurat")) {
  sce <- as.SingleCellExperiment(obj, assay = assay_name)
} else if (inherits(obj, "SingleCellExperiment")) {
  sce <- obj
} else {
  stop("RDS object is neither Seurat nor SingleCellExperiment")
}

meta <- as.data.frame(colData(sce))
meta$cohort_id <- cid

# Preserve common metadata columns if present
safe_char <- function(vec) {
  if (is.null(vec)) return(rep(NA_character_, length(vec)))
  as.character(vec)
}

if ("patient_id" %in% colnames(meta)) {
  meta$patient_id <- safe_char(meta$patient_id)
} else if ("patient" %in% colnames(meta)) {
  meta$patient_id <- safe_char(meta$patient_id)
}
if ("timepoint" %in% colnames(meta)) {
  meta$timepoint_raw <- safe_char(meta$timepoint)
}
if ("treatment" %in% colnames(meta)) {
  meta$treatment_regimen_raw <- safe_char(meta$treatment)
}
if ("response" %in% colnames(meta)) {
  meta$response_raw <- safe_char(meta$response)
}

colData(sce) <- DataFrame(meta)

message("[", cid, "] writing h5ad …")
writeH5AD(sce, out_path, compression = "gzip")
message("  written ", out_path, "  dims: ", nrow(sce), " x ", ncol(sce))
message("Done.")
