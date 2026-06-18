#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(AnnotationDbi)
  library(org.Hs.eg.db)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("usage: build_ensembl_symbol_map.R OUTPUT.csv")
output <- args[[1L]]
dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
keys <- keys(org.Hs.eg.db, keytype = "ENSEMBL")
mapped <- select(org.Hs.eg.db, keys = keys, keytype = "ENSEMBL", columns = "SYMBOL")
mapped <- mapped[!is.na(mapped$SYMBOL) & nzchar(mapped$SYMBOL), c("ENSEMBL", "SYMBOL")]
mapped <- unique(mapped[order(mapped$ENSEMBL, mapped$SYMBOL), ])
mapped$mapping_rank <- ave(mapped$SYMBOL, mapped$ENSEMBL, FUN = seq_along)
write.csv(mapped, output, row.names = FALSE, quote = TRUE)
