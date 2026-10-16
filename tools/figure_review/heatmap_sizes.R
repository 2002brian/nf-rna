#!/usr/bin/env Rscript
# Review aid: draw the standard-preset DEG and correlation heatmaps for synthetic matrices of
# 20/50/100 genes and 6/12/26/60 samples with short and long labels, using figures.R's own
# functions, and write a size table (heatmap_sizes.tsv) plus the PNG/PDF files.
#   heatmap_sizes.R OUTPUT_DIR [profile]
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("usage: heatmap_sizes.R OUTPUT_DIR [profile]")
out <- args[[1]]; profile <- if (length(args) >= 2) args[[2]] else "general"
dir.create(out, recursive = TRUE, showWarnings = FALSE)
script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[[1]])
repo <- normalizePath(file.path(dirname(normalizePath(script_file)), "..", ".."))
suppressPackageStartupMessages({ library(jsonlite); library(ggplot2); library(pheatmap) })
source(file.path(repo, "src", "rnaseq", "r", "figures.R"))
base <- nf_rna_figure_context(list(figure_style = list(preset = "standard", profile = profile)))
hfig <- nf_rna_figure_scope(base, "heatmap"); cfig <- nf_rna_figure_scope(base, "correlation")
rows <- list()
set.seed(11)
for (genes in c(20, 50, 100)) for (samples in c(6, 12, 26)) for (label in c("short", "long")) {
  ids <- if (label == "short") sprintf("Gene%05d", seq_len(genes)) else sprintf("ENSMUSG%011d_Long_Symbol", seq_len(genes))
  cols <- if (label == "short") paste0("S", seq_len(samples)) else sprintf("Organoid_Patient_Derived_Line_%02d", seq_len(samples))
  group <- factor(rep(c("Control", "Treatment"), length.out = samples)); names(group) <- cols
  z <- matrix(rnorm(genes * samples), genes, samples, dimnames = list(ids, cols)); z[, group == "Treatment"] <- z[, group == "Treatment"] + rep(c(1, -1), length.out = genes)
  z <- pmax(pmin(z, 2), -2)
  annotation <- data.frame(condition = group, row.names = cols)
  colours <- list(condition = nf_rna_group_colours(levels(group), context = hfig))
  fit <- nf_rna_deg_heatmap(hfig, z, annotation, colours, nf_rna_cluster(z), nf_rna_cluster(t(z)), "Treatment vs Control\nTop DEGs · row z-score of VST", 2)
  stem <- sprintf("deg_%dgenes_%dsamples_%s", genes, samples, label)
  nf_rna_save_figure(hfig, fit$heat, file.path(out, stem), fit$width_mm, fit$height_mm, out)
  rows[[length(rows) + 1]] <- data.frame(figure = stem, width_mm = round(fit$width_mm, 1), height_mm = round(fit$height_mm, 1), cell_w_mm = round(fit$cw, 2), cell_h_mm = round(fit$ch, 2),
    gene_labels = fit$show_rows, sample_labels = fit$show_cols, row_pt = fit$row_pt, max_w = hfig$profile$max_width_mm, max_h = hfig$profile$max_height_mm)
}
for (samples in c(6, 12, 26, 60)) for (label in c("short", "long")) {
  cols <- if (label == "short") paste0("S", seq_len(samples)) else sprintf("Organoid_Patient_Derived_Line_%02d", seq_len(samples))
  m <- matrix(rnorm(2000 * samples, 8, 2), 2000, samples, dimnames = list(NULL, cols)); m <- m + rep(rnorm(samples, 0, 0.3), each = 2000)
  corr <- cor(m); group <- factor(rep(c("Control", "Treatment"), length.out = samples)); annotation <- data.frame(condition = group, row.names = cols)
  colours <- list(condition = nf_rna_group_colours(levels(group), context = cfig))
  breaks <- seq(floor(min(corr) * 100) / 100, 1, length.out = 101)
  fit <- nf_rna_correlation_heatmap(cfig, corr, nf_rna_sequential_palette(100, cfig), breaks, annotation, colours, "Pearson correlation (blind VST)")
  stem <- sprintf("corr_%dsamples_%s", samples, label)
  nf_rna_save_figure(cfig, fit$heat, file.path(out, stem), fit$width_mm, fit$height_mm, out)
  rows[[length(rows) + 1]] <- data.frame(figure = stem, width_mm = round(fit$width_mm, 1), height_mm = round(fit$height_mm, 1), cell_w_mm = round(fit$cw, 2), cell_h_mm = round(fit$ch, 2),
    gene_labels = fit$show_rows, sample_labels = fit$show_cols, row_pt = fit$row_pt, max_w = cfig$profile$max_width_mm, max_h = cfig$profile$max_height_mm)
}
table <- do.call(rbind, rows)
write.table(table, file.path(out, "heatmap_sizes.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
print(table, row.names = FALSE)
