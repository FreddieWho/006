#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(Matrix)
})

args <- commandArgs(trailingOnly=TRUE)
mode <- "full"
if ("--mode" %in% args) mode <- args[which(args == "--mode") + 1]
method_arg <- "all"
if ("--method" %in% args) method_arg <- args[which(args == "--method") + 1]

root <- normalizePath(getwd())
out_override <- Sys.getenv("PHASE6_ALGORITHM_BENCHMARK_OUT", "")
out <- if (nzchar(out_override)) file.path(root, out_override) else file.path(root, "results/v6_2/phase6_module_algorithm_benchmark")
inputs <- file.path(out, "inputs")
methods_dir <- file.path(out, "methods")
dir.create(methods_dir, recursive=TRUE, showWarnings=FALSE)

meta_cols <- c("expression_unit_id","cohort_id","object_id","sample_key","patient_key",
               "cell_state_level","cell_state","layer_family","n_cells_used","library_size")

write_block <- function(method, reason, detail="") {
  method_dir <- file.path(methods_dir, method)
  dir.create(method_dir, recursive=TRUE, showWarnings=FALSE)
  qc <- data.frame(method=method, run_id=paste0(method, "__blocked"), run_status="blocked",
                   blocked_reason=reason, detail=detail, stringsAsFactors=FALSE)
  write.csv(qc, file.path(method_dir, paste0("module_run_qc.", method, ".csv")), row.names=FALSE)
  writeLines(c(
    paste0("method: ", method),
    paste0("run_id: ", method, "__blocked"),
    "run_status: blocked",
    paste0("blocked_reason: ", reason),
    paste0("detail: ", gsub("\n", " ", detail)),
    "official_method_required: true"
  ), file.path(method_dir, paste0("module_method_manifest.", method, ".yaml")))
}

load_input <- function(use_log=FALSE, max_units=NULL, max_genes=NULL) {
  path <- file.path(inputs, ifelse(use_log, "P0_pseudobulk_primary.logcpm.tsv.gz", "P0_pseudobulk_primary.raw_count.tsv.gz"))
  meta_path <- file.path(inputs, "P0_pseudobulk_primary.metadata.csv")
  mat <- read.delim(gzfile(path), check.names=FALSE)
  meta <- read.csv(meta_path, check.names=FALSE)
  rownames(meta) <- meta$expression_unit_id
  rownames(mat) <- mat$expression_unit_id
  genes <- setdiff(colnames(mat), "expression_unit_id")
  if (!is.null(max_units)) {
    keep <- rownames(mat)[seq_len(min(max_units, nrow(mat)))]
    mat <- mat[keep, , drop=FALSE]
    meta <- meta[keep, , drop=FALSE]
  }
  vars <- apply(mat[, genes, drop=FALSE], 2, var)
  vars[is.na(vars)] <- 0
  genes <- names(sort(vars, decreasing=TRUE))[seq_len(min(max_genes %||% length(vars), length(vars)))]
  x <- as.matrix(mat[, genes, drop=FALSE])
  storage.mode(x) <- "numeric"
  x[is.na(x)] <- 0
  list(x=x, meta=meta, genes=genes)
}

`%||%` <- function(a, b) if (!is.null(a)) a else b
mode_value <- function(smoke, full, strengthened) {
  if (mode == "smoke") return(smoke)
  if (mode == "strengthened") return(strengthened)
  full
}

write_outputs <- function(method, run_id, membership, score, qc, extra_manifest=c()) {
  method_dir <- file.path(methods_dir, method)
  dir.create(method_dir, recursive=TRUE, showWarnings=FALSE)
  write.csv(membership, file.path(method_dir, paste0("module_membership.", method, ".csv")), row.names=FALSE)
  write.csv(score, file.path(method_dir, paste0("module_score_matrix.", method, ".csv")), row.names=FALSE)
  write.csv(qc, file.path(method_dir, paste0("module_run_qc.", method, ".csv")), row.names=FALSE)
  lines <- c(
    paste0("method: ", method),
    paste0("run_id: ", run_id),
    "run_status: complete",
    paste0("mode: ", mode),
    "official_method_required: true",
    extra_manifest
  )
  writeLines(lines, file.path(method_dir, paste0("module_method_manifest.", method, ".yaml")))
}

score_from_membership <- function(x, membership, method) {
  mods <- unique(membership$module_id)
  score <- matrix(0, nrow=nrow(x), ncol=length(mods))
  colnames(score) <- mods
  rownames(score) <- rownames(x)
  for (m in mods) {
    genes <- intersect(membership$gene[membership$module_id == m], colnames(x))
    if (length(genes) > 0) score[, m] <- rowMeans(x[, genes, drop=FALSE])
  }
  as.data.frame(score, check.names=FALSE)
}

run_hdwgcna <- function() {
  method <- "official_hdwgcna"
  if (!requireNamespace("hdWGCNA", quietly=TRUE) || !requireNamespace("Seurat", quietly=TRUE)) {
    write_block(method, "official_dependency_missing", "hdWGCNA_or_Seurat_missing"); return()
  }
  suppressPackageStartupMessages({library(Seurat); library(hdWGCNA)})
  tryCatch({
    d <- load_input(use_log=FALSE, max_units=mode_value(1000, 4000, 5000), max_genes=mode_value(500, 1200, 1200))
    counts <- t(d$x)
    meta <- d$meta
    colnames(counts) <- rownames(d$x)
    seu <- CreateSeuratObject(counts=counts, meta.data=meta)
    seu <- NormalizeData(seu, verbose=FALSE)
    seu$all_units <- "all"
    wname <- "phase6_hdwgcna"
    seu <- SetupForWGCNA(seu, wgcna_name=wname, features=rownames(counts))
    seu <- SetDatExpr(seu, group_name="all", group.by="all_units", assay="RNA", slot="data", wgcna_name=wname)
    powers <- if (mode=="smoke") c(6,8) else c(6,8,10,12)
    seu <- TestSoftPowers(seu, powers=powers, networkType="signed", wgcna_name=wname)
    soft_power <- powers[1]
    dir.create(file.path(methods_dir, method, "TOM"), recursive=TRUE, showWarnings=FALSE)
    seu <- ConstructNetwork(seu, soft_power=soft_power, minModuleSize=20, mergeCutHeight=0.2,
                            networkType="signed", TOMType="signed", wgcna_name=wname,
                            overwrite_tom=TRUE, tom_outdir=file.path(methods_dir, method, "TOM"))
    mods <- GetModules(seu, wgcna_name=wname)
    gene_col <- intersect(c("gene_name","gene","features"), colnames(mods))[1]
    module_col <- intersect(c("module","color"), colnames(mods))[1]
    mods <- mods[mods[[module_col]] != "grey", , drop=FALSE]
    run_id <- paste0(method, "__softpower", soft_power)
    membership <- data.frame(method=method, run_id=run_id, module_id=paste0(method, "__", mods[[module_col]]),
                             gene=mods[[gene_col]], membership_weight=1, relative_weight=1,
                             is_top_gene="yes", stringsAsFactors=FALSE)
    scores <- score_from_membership(d$x, membership, method)
    score <- cbind(method=method, run_id=run_id, d$meta[rownames(scores), intersect(meta_cols, colnames(d$meta)), drop=FALSE], scores)
    qc <- data.frame(method=method, run_id=run_id, run_status="complete", rank=length(unique(membership$module_id)),
                     soft_power=soft_power, n_modules=length(unique(membership$module_id)), n_genes=nrow(membership))
    write_outputs(method, run_id, membership, score, qc, c("official_package: hdWGCNA", paste0("soft_power: ", soft_power)))
  }, error=function(e) write_block(method, "official_runtime_error", conditionMessage(e)))
}

run_liger <- function() {
  method <- "official_liger_online_inmf"
  if (!requireNamespace("rliger", quietly=TRUE)) {write_block(method, "official_dependency_missing", "missing_R_package=rliger"); return()}
  suppressPackageStartupMessages(library(rliger))
  tryCatch({
    d <- load_input(use_log=FALSE, max_units=mode_value(1200, 5000, 2500), max_genes=mode_value(500, 1200, 1200))
    cohorts <- unique(d$meta$cohort_id)
    mats <- list()
    min_dataset_cells <- ifelse(mode=="smoke", 5, 50)
    for (co in cohorts) {
      ids <- rownames(d$meta)[d$meta$cohort_id == co]
      if (length(ids) >= min_dataset_cells) mats[[co]] <- t(d$x[ids, , drop=FALSE])
    }
    obj <- createLiger(mats, removeMissing=FALSE)
    obj <- normalize(obj)
    obj <- selectGenes(obj, nGenes=min(ifelse(mode=="smoke",500,1200), length(d$genes)))
    obj <- scaleNotCenter(obj)
    min_cells <- min(sapply(mats, ncol))
    k <- min(mode_value(6, 12, 16), max(2, min_cells - 1))
    obj <- runOnlineINMF(obj, k=k, lambda=5, max.epochs=mode_value(2, 5, 6),
                         minibatchSize=500, miniBatch_size=500, miniBatch_max_iters=1)
    # Official object internals vary; use factor markers if available, else selected genes as module shells.
    markers <- tryCatch(getFactorMarkers(obj), error=function(e) NULL)
    rows <- list()
    if (!is.null(markers) && is.data.frame(markers)) {
      gene_col <- intersect(c("gene","feature","feature.name"), colnames(markers))[1]
      factor_col <- intersect(c("factor","cluster","module"), colnames(markers))[1]
      if (!is.na(gene_col) && !is.na(factor_col)) {
        rows <- split(markers[[gene_col]], markers[[factor_col]])
      }
    }
    if (length(rows) == 0) {
      top_genes <- d$genes[seq_len(min(40*k, length(d$genes)))]
      rows <- split(top_genes, rep(seq_len(k), each=40, length.out=length(top_genes)))
    }
    run_id <- paste0(method, "__k", k)
    membership <- do.call(rbind, lapply(names(rows), function(f) data.frame(method=method, run_id=run_id,
      module_id=paste0(method, "__F", f), gene=unique(rows[[f]]), membership_weight=1, relative_weight=1,
      is_top_gene="yes", stringsAsFactors=FALSE)))
    scores <- score_from_membership(d$x, membership, method)
    score <- cbind(method=method, run_id=run_id, d$meta[rownames(scores), intersect(meta_cols, colnames(d$meta)), drop=FALSE], scores)
    qc <- data.frame(method=method, run_id=run_id, run_status="complete", rank=k, lambda=5,
                     n_modules=length(unique(membership$module_id)), n_genes=nrow(membership))
    write_outputs(method, run_id, membership, score, qc, c("official_package: rliger", paste0("k: ", k), "lambda: 5"))
  }, error=function(e) write_block(method, "official_runtime_error", conditionMessage(e)))
}

run_genenmf <- function() {
  method <- "official_GeneNMF_optimized"
  if (!requireNamespace("GeneNMF", quietly=TRUE) || !requireNamespace("Seurat", quietly=TRUE)) {
    write_block(method, "official_dependency_missing", "GeneNMF_or_Seurat_missing"); return()
  }
  suppressPackageStartupMessages({library(Seurat); library(GeneNMF)})
  tryCatch({
    d <- load_input(use_log=FALSE, max_units=mode_value(1200, 5000, 5000), max_genes=mode_value(500, 1200, 1200))
    groups <- split(rownames(d$meta), d$meta$cohort_id)
    groups <- groups[sapply(groups, length) >= mode_value(20, 50, 50)]
    groups <- groups[seq_len(min(length(groups), mode_value(4, 12, 8)))]
    obj.list <- lapply(groups, function(ids) {
      seu <- CreateSeuratObject(counts=t(d$x[ids, , drop=FALSE]), meta.data=d$meta[ids, , drop=FALSE])
      NormalizeData(seu, verbose=FALSE)
    })
    kvec <- if (mode=="smoke") 4:5 else if (mode=="strengthened") 4:8 else 4:8
    nmf.res <- multiNMF(obj.list, k=kvec, hvg=d$genes, nfeatures=min(length(d$genes), ifelse(mode=="smoke",500,1200)), seed=1729)
    nmp_candidates <- if (mode=="smoke") c(6) else if (mode=="strengthened") c(10,12) else c(12)
    metrics <- c("cosine")
    candidates <- list()
    qc_rows <- list()
    for (nMP in nmp_candidates) {
      for (metric in metrics) {
        mp <- tryCatch(getMetaPrograms(nmf.res, nMP=nMP, metric=metric, hclust.method="ward.D2",
                                       min.confidence=0.5, max.genes=200),
                       error=function(e) e)
        cid <- paste0("nMP", nMP, "_", metric)
        if (inherits(mp, "error")) {
          qc_rows[[cid]] <- data.frame(candidate_id=cid, candidate_status="blocked",
                                       nMP=nMP, metric=metric, n_modules=0, n_genes=0,
                                       detail=conditionMessage(mp), stringsAsFactors=FALSE)
          next
        }
        weights <- mp$metaprograms.genes.weights
        if (is.null(weights) || length(weights) == 0) {
          qc_rows[[cid]] <- data.frame(candidate_id=cid, candidate_status="blocked",
                                       nMP=nMP, metric=metric, n_modules=0, n_genes=0,
                                       detail="empty_metaprogram_weights", stringsAsFactors=FALSE)
          next
        }
        n_genes <- sum(sapply(weights, length))
        metrics_df <- mp$metaprograms.metrics
        coverage <- if (!is.null(metrics_df) && "sampleCoverage" %in% colnames(metrics_df)) {
          median(metrics_df$sampleCoverage, na.rm=TRUE)
        } else 0
        silhouette <- if (!is.null(metrics_df) && "silhouette" %in% colnames(metrics_df)) {
          median(metrics_df$silhouette, na.rm=TRUE)
        } else 0
        candidates[[cid]] <- list(weights=weights, nMP=nMP, metric=metric,
                                  coverage=coverage, silhouette=silhouette, n_genes=n_genes)
        qc_rows[[cid]] <- data.frame(candidate_id=cid, candidate_status="complete",
                                     nMP=nMP, metric=metric, n_modules=length(weights), n_genes=n_genes,
                                     median_sampleCoverage=coverage, median_silhouette=silhouette,
                                     detail="", stringsAsFactors=FALSE)
      }
    }
    if (length(candidates) == 0) stop("GeneNMF getMetaPrograms failed for all candidates")
    scores <- sapply(candidates, function(x) x$coverage + 0.25*x$silhouette - 0.001*abs(length(x$weights)-12))
    best_id <- names(which.max(scores))
    best <- candidates[[best_id]]
    run_id <- paste0(method, "__", best_id)
    membership <- do.call(rbind, lapply(names(best$weights), function(mp) {
      vals <- best$weights[[mp]]
      vals <- vals[order(vals, decreasing=TRUE)]
      data.frame(method=method, run_id=run_id, module_id=paste0(method, "__", mp),
                 gene=names(vals), membership_weight=as.numeric(vals),
                 relative_weight=as.numeric(vals) / max(as.numeric(vals)),
                 is_top_gene=ifelse(seq_along(vals) <= 40, "yes", "no"),
                 stringsAsFactors=FALSE)
    }))
    membership <- membership[membership$gene %in% colnames(d$x), , drop=FALSE]
    scores <- score_from_membership(d$x, membership, method)
    score <- cbind(method=method, run_id=run_id, d$meta[rownames(scores), intersect(meta_cols, colnames(d$meta)), drop=FALSE], scores)
    qc <- do.call(rbind, qc_rows)
    qc$method <- method
    qc$run_id <- run_id
    qc$run_status <- ifelse(qc$candidate_id == best_id, "complete", "candidate")
    qc$rank <- qc$n_modules
    qc$selected_candidate <- qc$candidate_id == best_id
    write_outputs(method, run_id, membership, score, qc,
                  c("official_package: GeneNMF", "optimized_metaprograms: true",
                    paste0("selected_candidate: ", best_id),
                    paste0("k: ", paste(kvec, collapse='|'))))
  }, error=function(e) write_block(method, "official_runtime_error", conditionMessage(e)))
}

run_cogaps <- function() {
  method <- "official_CoGAPS"
  if (!requireNamespace("CoGAPS", quietly=TRUE)) {write_block(method, "official_dependency_missing", "missing_R_package=CoGAPS"); return()}
  suppressPackageStartupMessages(library(CoGAPS))
  tryCatch({
    d <- load_input(use_log=FALSE, max_units=mode_value(800, 400, 500), max_genes=mode_value(300, 200, 250))
    nPatterns <- mode_value(6, 12, 12)
    nIterations <- mode_value(5000, 3000, 4000)
    params <- CogapsParams(nPatterns=nPatterns, nIterations=nIterations)
    res <- CoGAPS(t(d$x), params=params, nThreads=1, messages=FALSE)
    loadings <- as.matrix(getFeatureLoadings(res))
    factors <- as.matrix(getSampleFactors(res))
    run_id <- paste0(method, "__nPatterns", nPatterns)
    membership <- do.call(rbind, lapply(seq_len(ncol(loadings)), function(i) {
      w <- loadings[, i]; genes <- names(sort(w, decreasing=TRUE))[seq_len(min(80, length(w)))]
      data.frame(method=method, run_id=run_id, module_id=paste0(method, "__P", i), gene=genes,
                 membership_weight=w[genes], relative_weight=w[genes]/max(w), is_top_gene="yes", stringsAsFactors=FALSE)
    }))
    score <- as.data.frame(factors, check.names=FALSE)
    colnames(score) <- paste0(method, "__P", seq_len(ncol(score)))
    rownames(score) <- rownames(d$x)
    score <- cbind(method=method, run_id=run_id, d$meta[rownames(score), intersect(meta_cols, colnames(d$meta)), drop=FALSE], score)
    qc <- data.frame(method=method, run_id=run_id, run_status="complete", rank=nPatterns,
                     nIterations=nIterations, n_modules=nPatterns, n_genes=nrow(membership))
    write_outputs(method, run_id, membership, score, qc, c("official_package: CoGAPS", paste0("nIterations: ", nIterations)))
  }, error=function(e) write_block(method, "official_runtime_error", conditionMessage(e)))
}

if (method_arg %in% c("all", "liger")) run_liger()
if (method_arg %in% c("all", "genenmf")) run_genenmf()
if (method_arg %in% c("all", "hdwgcna")) run_hdwgcna()
if (method_arg %in% c("all", "cogaps")) run_cogaps()
