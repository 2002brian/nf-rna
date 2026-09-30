args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2 || args[[1]] != "--config") stop("usage: l1_analysis.R --config CONFIG.json")
script_arg <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", script_arg[grep("^--file=", script_arg)][[1]])
source(file.path(dirname(normalizePath(script_file)), "provenance.R"))
source(file.path(dirname(normalizePath(script_file)), "design_metadata.R"))
source(file.path(dirname(normalizePath(script_file)), "sample_inputs.R"))
source(file.path(dirname(normalizePath(script_file)), "figures.R"))
suppressPackageStartupMessages({ library(jsonlite); library(yaml); library(DESeq2); library(ggplot2); library(pheatmap); library(ggrepel) })
cfg <- fromJSON(args[[2]], simplifyVector = FALSE)
dir.create(cfg$output_dir, recursive = TRUE, showWarnings = FALSE)
samples <- unlist(cfg$samples, use.names = FALSE)
metadata <- nf_rna_design_metadata(cfg, samples)
formula <- as.formula(cfg$formula)
if (cfg$source_type == "raw_counts" || cfg$source_type == "featurecounts_raw_counts") {
  table <- read.csv(cfg$counts, check.names = FALSE, stringsAsFactors = FALSE)
  ids <- table[[1]]
  matrix_counts <- as.matrix(table[, samples, drop = FALSE])
  rownames(matrix_counts) <- ids
  storage.mode(matrix_counts) <- "numeric"
  dds <- DESeqDataSetFromMatrix(countData = matrix_counts, colData = metadata, design = formula)
  source_counts <- matrix_counts
} else if (cfg$source_type == "salmon_tximport") {
  suppressPackageStartupMessages(library(tximport))
  files <- nf_rna_quant_sf_files(cfg$quant_sf, samples)
  tx2gene <- read.delim(cfg$tx2gene, check.names = FALSE, stringsAsFactors = FALSE)[, 1:2]
  txi <- tximport(files, type = "salmon", tx2gene = tx2gene)
  source_counts <- nf_rna_tximport_counts(txi, samples)
  dds <- DESeqDataSetFromTximport(txi, colData = metadata, design = formula)
} else stop("unsupported source_type")
all_zero <- rowSums(source_counts) == 0
low_total <- !all_zero & rowSums(source_counts) < cfg$filter$minimum_total_count
keep <- !(all_zero | low_total)
if (sum(keep) < 2) stop("Filtering left fewer than two genes; L1 QC cannot proceed.")
dds <- dds[keep, ]
dds <- estimateSizeFactors(dds)
normalized <- counts(dds, normalized = TRUE)
# DESeq2's fast vst() expects a large enough feature set for its default
# subsampling.  The exact VST remains the appropriate deterministic fallback
# for compact, legitimate smoke fixtures.
vsd <- if (nrow(dds) < 1000) varianceStabilizingTransformation(dds, blind = TRUE) else vst(dds, blind = TRUE)
vst_matrix <- assay(vsd)
write_matrix <- function(x, path) write.table(data.frame(gene_id = rownames(x), x, check.names = FALSE), path, sep = "\t", quote = FALSE, row.names = FALSE)
write_matrix(normalized, file.path(cfg$output_dir, "normalized_counts.tsv"))
write_matrix(vst_matrix, file.path(cfg$output_dir, "vst.tsv"))
write.csv(data.frame(gene_id = rownames(source_counts), source_counts, check.names = FALSE), file.path(cfg$output_dir, "source_counts.csv"), row.names = FALSE, quote = FALSE)
write.csv(data.frame(gene_id = rownames(vst_matrix), vst_matrix, check.names = FALSE), file.path(cfg$output_dir, "vst.csv"), row.names = FALSE, quote = FALSE)
size_factor_values <- sizeFactors(dds)
# DESeqDataSetFromTximport may use per-gene normalization factors (offsets)
# rather than scalar size factors; retain that distinction in the QC table.
if (is.null(size_factor_values)) size_factor_values <- rep(NA_real_, length(samples))
library_size <- data.frame(sample_id = samples, input_total = as.numeric(colSums(source_counts)), retained_total = as.numeric(colSums(source_counts[keep, , drop = FALSE])), size_factor = as.numeric(size_factor_values), normalized_total = as.numeric(colSums(normalized)), check.names = FALSE)
write.table(library_size, file.path(cfg$output_dir, "library_size_qc.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
pca <- prcomp(t(vst_matrix), center = TRUE, scale. = FALSE)
variance <- pca$sdev^2 / sum(pca$sdev^2)
scores <- data.frame(sample_id = samples, PC1 = pca$x[, 1], PC2 = if (ncol(pca$x) >= 2) pca$x[, 2] else 0, check.names = FALSE)
if (!is.null(cfg$design_variable_types)) for (variable in names(cfg$design_variable_types)) scores[[variable]] <- metadata[scores$sample_id, variable]
write.table(scores, file.path(cfg$output_dir, "pca_scores.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(data.frame(component = paste0("PC", seq_along(variance)), proportion_variance = variance), file.path(cfg$output_dir, "pca_variance.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
color_name <- tail(all.vars(formula), 1)
scores$group <- as.factor(metadata[scores$sample_id, color_name])
fig <- nf_rna_figure_context(cfg)
pc_labels <- labs(x = sprintf("PC1 (%.1f%%)", 100 * variance[[1]]), y = sprintf("PC2 (%.1f%%)", 100 * ifelse(length(variance) >= 2, variance[[2]], 0)))
# PC scores share one unit, so both axes use the same scale (coord_fixed).
NF_RNA_PCA_LABEL_SEED <- 20260930L
pca_plot <- function(colour_by, colour_name, shape_by, shape_name) {
  colours <- nf_rna_categorical_colors(levels(scores[[colour_by]]))
  shapes <- nf_rna_categorical_shapes(levels(scores[[shape_by]]))
  mapping <- if (is.null(shapes)) aes(x = PC1, y = PC2, color = .data[[colour_by]]) else aes(x = PC1, y = PC2, color = .data[[colour_by]], shape = .data[[shape_by]])
  # Repelled labels keep every sample ID legible in tight clusters.  The seed and
  # an iteration-only stopping rule (no wall-clock limit) make the layout
  # deterministic; labels stay inside the panel and never drop a sample.
  p <- ggplot(scores, mapping) + geom_point(size = 2) +
    geom_text_repel(aes(label = sample_id), size = nf_rna_pt_to_mm(fig$profile$min_pt), family = fig$font$family, show.legend = FALSE,
                    seed = NF_RNA_PCA_LABEL_SEED, max.time = Inf, max.iter = 10000, max.overlaps = Inf,
                    box.padding = 0.3, point.padding = 0.3, min.segment.length = 0.4, segment.size = nf_rna_linewidth(0.3), segment.colour = "grey50") +
    scale_color_manual(values = colours, name = colour_name) +
    scale_x_continuous(expand = expansion(mult = 0.05)) + scale_y_continuous(expand = expansion(mult = 0.15)) +
    coord_fixed(clip = "off") + pc_labels + nf_rna_theme(fig)
  if (!is.null(shapes)) p <- p + scale_shape_manual(values = shapes, name = shape_name)
  list(plot = p, colours = colours, shapes = shapes)
}
# With equal scaling the panel's aspect ratio is the data's; size the figure
# to it (instead of padding a fixed height with white space).
pca_height_mm <- function() {
  x_range <- max(diff(range(scores$PC1)), 1e-9); y_range <- max(diff(range(scores$PC2)), x_range * 0.02)
  panel_width <- fig$profile$double_col_mm - 45
  min(fig$profile$max_height_mm, max(50, 25 + panel_width * (y_range * 1.3) / (x_range * 1.1)))
}
pca_meta <- function(encoding) list(
  kind = "PCA", data = "blind VST, all retained genes, prcomp(center = TRUE, scale. = FALSE)", axes = "PC1/PC2 with proportion of variance explained; equal axis scaling",
  encoding = encoding, ellipses = "none", paired_lines = "none",
  labels = list(text = "every sample_id, unabbreviated", placement = sprintf("ggrepel::geom_text_repel, seed %d, iteration-limited (max.iter 10000, no time limit), kept inside the panel; leader lines only for displaced labels", NF_RNA_PCA_LABEL_SEED)),
  source_data = list(nf_rna_source_data(file.path(cfg$output_dir, "pca_scores.tsv"), cfg$output_dir, "PC coordinates"), nf_rna_source_data(file.path(cfg$output_dir, "pca_variance.tsv"), cfg$output_dir, "proportion of variance"))
)
main_pca <- pca_plot("group", color_name, "group", color_name)
nf_rna_save_figure(fig, main_pca$plot, file.path(cfg$output_dir, "pca"), fig$profile$double_col_mm, pca_height_mm(), cfg$output_dir,
  pca_meta(list(colour = color_name, shape = if (is.null(main_pca$shapes)) "none (too many levels)" else color_name, colours = as.list(main_pca$colours), palette = nf_rna_categorical_palette_name(length(main_pca$colours)))))
if (!is.null(cfg$design_variable_types) && identical(cfg$design_variable_types[["batch"]], "categorical")) {
  scores$batch <- as.factor(metadata[scores$sample_id, "batch"])
  # Batch is the colour; the biological group keeps its meaning as the shape.
  batch_pca <- pca_plot("batch", "batch", "group", color_name)
  nf_rna_save_figure(fig, batch_pca$plot, file.path(cfg$output_dir, "pca_by_batch"), fig$profile$double_col_mm, pca_height_mm(), cfg$output_dir,
    pca_meta(list(colour = "batch", shape = if (is.null(batch_pca$shapes)) "none (too many levels)" else color_name, colours = as.list(batch_pca$colours), palette = nf_rna_categorical_palette_name(length(batch_pca$colours)))))
}
corr <- cor(vst_matrix, method = "pearson")
write.table(data.frame(sample_id = rownames(corr), corr, check.names = FALSE), file.path(cfg$output_dir, "sample_correlation.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
# Pearson r has a meaningful zero: a diverging scale centred on 0 is used only
# when negative correlations occur; otherwise a sequential scale spans the
# observed range up to 1 so that QC differences stay visible.
if (all(corr >= 0, na.rm = TRUE)) {
  corr_scale <- nf_rna_sequential_limits(floor(min(corr, na.rm = TRUE) * 100) / 100, upper = 1)
  corr_breaks <- seq(corr_scale$limits[[1]], corr_scale$limits[[2]], length.out = 101); corr_colours <- nf_rna_sequential_palette(100); corr_kind <- "sequential"
} else {
  corr_breaks <- seq(-1, 1, length.out = 101); corr_colours <- nf_rna_diverging_palette(100); corr_kind <- "diverging, midpoint 0"
}
corr_plot <- pheatmap(corr, color = corr_colours, breaks = corr_breaks, border_color = NA, main = "Pearson correlation (blind VST)",
  fontsize = fig$profile$base_pt, fontsize_row = fig$profile$small_pt, fontsize_col = fig$profile$small_pt, silent = TRUE)
nf_rna_save_figure(fig, corr_plot, file.path(cfg$output_dir, "sample_correlation"), fig$profile$double_col_mm, 150, cfg$output_dir,
  list(kind = "sample correlation heatmap", data = "Pearson correlation of blind VST", clustering = "pheatmap default: complete linkage of Euclidean distances between correlation profiles",
       colour_scale = list(kind = corr_kind, limits = range(corr_breaks), na_colour = "pheatmap default"),
       source_data = list(nf_rna_source_data(file.path(cfg$output_dir, "sample_correlation.tsv"), cfg$output_dir, "correlation matrix"))))
nf_rna_write_figure_manifest(fig, cfg$output_dir)
flags <- character()
median_library <- median(library_size$input_total)
for (i in seq_len(nrow(library_size))) if (library_size$input_total[[i]] < median_library * 0.25) flags <- c(flags, paste0("LOW_LIBRARY_SIZE: ", library_size$sample_id[[i]], " is below 25% of the median input total."))
summary <- list(source_type = cfg$source_type, samples = samples, genes_input = nrow(source_counts), genes_removed_all_zero = sum(all_zero), genes_removed_low_total = sum(low_total), genes_retained = sum(keep), filter = cfg$filter, normalization = list(method = if (cfg$source_type == "salmon_tximport") "DESeq2 with tximport average-transcript-length normalization offsets" else "DESeq2 median-ratio size factors", scalar_size_factors_available = !all(is.na(size_factor_values))), vst = list(method = if (nrow(dds) < 1000) "DESeq2 varianceStabilizingTransformation (exact small-matrix path)" else "DESeq2 vst"), package_versions = list(R = R.version.string, DESeq2 = as.character(packageVersion("DESeq2")), tximport = if (cfg$source_type == "salmon_tximport") as.character(packageVersion("tximport")) else NULL), pca_proportion_variance = as.list(variance), flags = as.list(flags))
write(toJSON(summary, auto_unbox = TRUE, pretty = TRUE, null = "null"), file.path(cfg$output_dir, "backend_summary.json"))
nf_rna_write_provenance(cfg, cfg$output_dir, "L1", "SUCCESS", c("DESeq2", "tximport", "ggplot2", "ggrepel", "pheatmap", "scales", "systemfonts", "yaml", "jsonlite"), summary)
