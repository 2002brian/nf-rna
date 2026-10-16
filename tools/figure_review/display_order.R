#!/usr/bin/env Rscript
# Recompute, from a run's source-data tables, the row/column order the legacy code path draws for the
# sample correlation heatmap (pheatmap's own clustering) and for every DEG heatmap (nf_rna_cluster on
# the unclipped z-scores).  Prints JSON {figure path: {rows, columns}}; review_matrix.py compares it
# with the order recorded in the standard run's figure manifest.
#   display_order.R RUN_DIR
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) stop("usage: display_order.R RUN_DIR")
script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[[1]])
repo <- normalizePath(file.path(dirname(normalizePath(script_file)), "..", ".."))
suppressPackageStartupMessages({ library(jsonlite); library(pheatmap) })
source(file.path(repo, "src", "rnaseq", "r", "figures.R"))
root <- normalizePath(args[[1]]); out <- list()
for (path in list.files(root, "^sample_correlation\\.tsv$", recursive = TRUE, full.names = TRUE)) {
  table <- read.delim(path, check.names = FALSE); m <- as.matrix(table[, -1]); rownames(m) <- table$sample_id
  heat <- pheatmap(m, silent = TRUE)
  out[[sub(paste0("^", root, "/"), "", sub("\\.tsv$", "", path))]] <- list(rows = rownames(m)[heat$tree_row$order], columns = colnames(m)[heat$tree_col$order])
}
for (path in list.files(root, "^heatmap_zscores\\.tsv$", recursive = TRUE, full.names = TRUE)) {
  table <- read.delim(path, check.names = FALSE); z <- as.matrix(table[, -1]); rownames(z) <- table$gene_id
  z <- z[rowSums(is.finite(z)) >= 2, , drop = FALSE]
  rows <- nf_rna_cluster(z); columns <- nf_rna_cluster(t(z))
  out[[sub(paste0("^", root, "/"), "", sub("heatmap_zscores\\.tsv$", "heatmap", path))]] <- list(
    rows = if (is.null(rows)) rownames(z) else rownames(z)[rows$order], columns = if (is.null(columns)) colnames(z) else colnames(z)[columns$order])
}
cat(toJSON(out, auto_unbox = FALSE))
