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
pfig <- nf_rna_figure_scope(fig, "pca")  # the context PCA figures are drawn with (standard preset when requested)
pc_labels <- labs(x = sprintf("PC1 (%.1f%%)", 100 * variance[[1]]), y = sprintf("PC2 (%.1f%%)", 100 * ifelse(length(variance) >= 2, variance[[2]], 0)))
# PC scores share one unit, so both axes use the same scale (coord_fixed).
NF_RNA_PCA_LABEL_SEED <- 20260930L
pca_y <- nf_rna_pca_y_axis(scores$PC1, scores$PC2)
pca_plot <- function(colour_by, colour_name, shape_by, shape_name) {
  colours <- if (identical(colour_by, "batch")) nf_rna_batch_colours(levels(scores[[colour_by]]), context = pfig) else nf_rna_group_colours(levels(scores[[colour_by]]), context = pfig)
  shapes <- nf_rna_categorical_shapes(levels(scores[[shape_by]]), pfig)
  mapping <- if (is.null(shapes)) aes(x = PC1, y = PC2, color = .data[[colour_by]]) else aes(x = PC1, y = PC2, color = .data[[colour_by]], shape = .data[[shape_by]])
  # Repelled labels keep every sample ID legible in tight clusters.  The seed and
  # an iteration-only stopping rule (no wall-clock limit) make the layout
  # deterministic; labels stay inside the panel and never drop a sample.
  p <- ggplot(scores, mapping) + geom_point(size = 2) +
    geom_text_repel(aes(label = sample_id), size = nf_rna_pt_to_mm(nf_rna_pca_label_pt(pfig)), family = pfig$font$family, show.legend = FALSE,
                    seed = NF_RNA_PCA_LABEL_SEED, max.time = Inf, max.iter = 10000, max.overlaps = Inf,
                    box.padding = 0.3, point.padding = 0.3, min.segment.length = 0.4, segment.size = nf_rna_linewidth(0.3), segment.colour = "grey50") +
    scale_color_manual(values = colours, name = colour_name) +
    scale_x_continuous(expand = expansion(mult = 0.05)) + scale_y_continuous(breaks = pca_y$breaks, expand = expansion(mult = 0.15)) +
    # A constant PC2 (for example two samples) gets a readable fixed extent and one break; see nf_rna_pca_y_axis().
    expand_limits(y = pca_y$limits) +
    coord_fixed(clip = "off") + pc_labels + nf_rna_theme(pfig) + pca_legend
  if (!is.null(shapes)) p <- p + scale_shape_manual(values = shapes, name = shape_name)
  list(plot = p, colours = colours, shapes = shapes)
}
# With equal scaling the panel's aspect ratio is the data's; size the figure
# to it (instead of padding a fixed height with white space).
# The legacy layout is the double column with the legend on the right; the standard preset
# uses the single column (its width comes from the preset's layout) with the legend below.
pca_width_mm <- nf_rna_layout_width_mm(pfig, "pca", fig$profile$double_col_mm)
pca_legend <- if (identical(pfig$preset, "standard")) theme(legend.position = "bottom", legend.box = "vertical") else NULL
pca_height_mm <- function() {
  x_range <- max(diff(range(scores$PC1)), 1e-9); y_range <- pca_y$extent
  if (identical(pfig$preset, "standard")) {
    panel_width <- pca_width_mm - 18; overhead <- 30
    return(min(pfig$profile$max_height_mm, max(nf_rna_layout_value(pfig, "pca", "min_height_mm", 50), overhead + panel_width * (y_range * 1.3) / (x_range * 1.1))))
  }
  panel_width <- fig$profile$double_col_mm - 45
  min(fig$profile$max_height_mm, max(50, 25 + panel_width * (y_range * 1.3) / (x_range * 1.1)))
}
pca_meta <- function(encoding) list(
  kind = "PCA", data = "blind VST, all retained genes, prcomp(center = TRUE, scale. = FALSE)", axes = paste0("PC1/PC2 with proportion of variance explained; equal axis scaling",
    if (pca_y$constant) sprintf("; PC2 is constant (no variance): y range fixed at +/-20%% of the PC1 span with a single break at %g", pca_y$breaks) else ""),
  encoding = encoding, ellipses = "none", paired_lines = "none",
  labels = list(text = "every sample_id, unabbreviated", placement = sprintf("ggrepel::geom_text_repel, seed %d, iteration-limited (max.iter 10000, no time limit), kept inside the panel; leader lines only for displaced labels", NF_RNA_PCA_LABEL_SEED)),
  source_data = list(nf_rna_source_data(file.path(cfg$output_dir, "pca_scores.tsv"), cfg$output_dir, "PC coordinates"), nf_rna_source_data(file.path(cfg$output_dir, "pca_variance.tsv"), cfg$output_dir, "proportion of variance"))
)
main_pca <- pca_plot("group", color_name, "group", color_name)
nf_rna_save_figure(pfig, main_pca$plot, file.path(cfg$output_dir, "pca"), pca_width_mm, pca_height_mm(), cfg$output_dir,
  pca_meta(list(colour = color_name, shape = if (is.null(main_pca$shapes)) "none (too many levels)" else color_name, colours = as.list(main_pca$colours), palette = nf_rna_categorical_palette_name(length(main_pca$colours), pfig))))
if (!is.null(cfg$design_variable_types) && identical(cfg$design_variable_types[["batch"]], "categorical")) {
  scores$batch <- as.factor(metadata[scores$sample_id, "batch"])
  # Batch is the colour; the biological group keeps its meaning as the shape.
  batch_pca <- pca_plot("batch", "batch", "group", color_name)
  nf_rna_save_figure(pfig, batch_pca$plot, file.path(cfg$output_dir, "pca_by_batch"), pca_width_mm, pca_height_mm(), cfg$output_dir,
    pca_meta(list(colour = "batch", shape = if (is.null(batch_pca$shapes)) "none (too many levels)" else color_name, colours = as.list(batch_pca$colours), palette = nf_rna_categorical_palette_name(length(batch_pca$colours), pfig, batch = TRUE))))
}
corr <- cor(vst_matrix, method = "pearson")
write.table(data.frame(sample_id = rownames(corr), corr, check.names = FALSE), file.path(cfg$output_dir, "sample_correlation.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
cfig <- nf_rna_figure_scope(fig, "correlation")  # the context the correlation heatmap is drawn with
# Pearson r has a meaningful zero: a diverging scale centred on 0 is used only
# when negative correlations occur; otherwise a sequential scale spans the
# observed range up to 1 so that QC differences stay visible.
if (all(corr >= 0, na.rm = TRUE)) {
  corr_scale <- nf_rna_sequential_limits(floor(min(corr, na.rm = TRUE) * 100) / 100, upper = 1)
  corr_breaks <- seq(corr_scale$limits[[1]], corr_scale$limits[[2]], length.out = 101); corr_colours <- nf_rna_sequential_palette(100, cfig); corr_kind <- "sequential"
} else {
  corr_breaks <- seq(-1, 1, length.out = 101); corr_colours <- nf_rna_diverging_palette(100, cfig); corr_kind <- "diverging, midpoint 0"
}
box_group <- nf_rna_sample_annotation(metadata, samples, color_name)
# The batch annotation is QC-only and independent of the DE design: any metadata column named
# "batch" with a value for every plotted sample is shown, whether or not the formula uses it.
box_batch <- NULL
box_batch_note <- "none (no batch column in metadata)"
if ("batch" %in% colnames(metadata)) {
  batch_values <- trimws(as.character(metadata[samples, "batch"]))
  if (anyNA(batch_values) || any(batch_values == "")) box_batch_note <- "none (batch column has missing values for some samples)"
  else box_batch <- nf_rna_sample_annotation(metadata, samples, "batch")
}
corr_meta <- function() list(kind = "sample correlation heatmap", data = "Pearson correlation of blind VST", clustering = "pheatmap default: complete linkage of Euclidean distances between correlation profiles",
       colour_scale = list(kind = corr_kind, limits = range(corr_breaks), na_colour = "pheatmap default"),
       source_data = list(nf_rna_source_data(file.path(cfg$output_dir, "sample_correlation.tsv"), cfg$output_dir, "correlation matrix")))
if (!identical(cfig$preset, "standard")) {
  corr_plot <- pheatmap(corr, color = corr_colours, breaks = corr_breaks, border_color = NA, main = "Pearson correlation (blind VST)",
    fontsize = fig$profile$base_pt, fontsize_row = fig$profile$small_pt, fontsize_col = fig$profile$small_pt, silent = TRUE)
  nf_rna_save_figure(fig, corr_plot, file.path(cfg$output_dir, "sample_correlation"), fig$profile$double_col_mm, 150, cfg$output_dir, corr_meta())
} else {
  # Same matrix, same pheatmap clustering; the standard style adds the group (and batch) annotation in the
  # PCA colours and sizes the figure from its content.
  corr_universe <- levels(box_group)
  corr_annotation <- data.frame(metadata[colnames(corr), color_name], row.names = colnames(corr), check.names = FALSE); names(corr_annotation) <- color_name
  corr_annotation_colours <- setNames(list(nf_rna_group_colours(corr_universe, unique(as.character(corr_annotation[[color_name]])), cfig)), color_name)
  if (!is.null(box_batch)) {
    corr_batch <- box_batch[match(colnames(corr), samples)]
    corr_annotation$batch <- corr_batch
    corr_annotation_colours$batch <- nf_rna_batch_colours(levels(box_batch), as.character(unique(corr_batch)), cfig)
  }
  pt <- cfig$profile
  fit <- nf_rna_correlation_heatmap(cfig, corr, corr_colours, corr_breaks, corr_annotation, corr_annotation_colours, "Pearson correlation (blind VST)")
  corr_heat <- fit$heat
  nf_rna_save_figure(cfig, corr_heat, file.path(cfg$output_dir, "sample_correlation"), fit$width_mm, fit$height_mm, cfg$output_dir,
    c(corr_meta(), list(
      annotation = list(group = color_name, colours = as.list(corr_annotation_colours[[color_name]]), batch = if (is.null(corr_annotation_colours$batch)) "none" else as.list(corr_annotation_colours$batch)),
      layout = list(cell_mm = fit$cw, sample_labels = if (fit$show_rows) sprintf("%.1f pt sample IDs", pt$small_pt) else "hidden: too many samples to label legibly within the profile maximum"),
      display_order = list(rows = rownames(corr)[fit$heat$tree_row$order], columns = colnames(corr)[fit$heat$tree_col$order]))))
}
# Sample-level expression distributions (QC only: they say nothing about batch
# effects).  The same retained genes and the same sample order (cfg$samples) are
# used before and after normalization.
raw_log2 <- log2(source_counts[keep, samples, drop = FALSE] + 1)
distribution_stats <- rbind(nf_rna_sample_box_stats(raw_log2, samples, "log2_raw_count_plus_1"), nf_rna_sample_box_stats(vst_matrix, samples, "blind_vst"))
write.table(distribution_stats, file.path(cfg$output_dir, "sample_distribution_qc.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
box_meta <- function(data, y_label, matrix_label, box) list(
  kind = "sample expression distribution boxplot", data = data, y_axis = y_label,
  box = "one box per sample: Tukey hinges (quartiles), median line, whiskers to the most extreme value within 1.5 x IQR; outliers not drawn",
  order = "cfg sample order, identical to the PCA/correlation tables", genes = "all genes retained by the L1 filter",
  encoding = list(fill = color_name, colours = as.list(box$colours), batch_marker = if (is.null(box_batch)) box_batch_note else if (is.null(box$shapes)) "none (too many batch levels)" else "black marker shape below each box", batch_in_de_design = identical(cfg$design_variable_types[["batch"]], "categorical") || "batch" %in% all.vars(formula), batch_values = "plotted values are not batch-corrected"),
  labels = if (box$labels_shown) sprintf("every sample_id, rotated %d degrees", box$rotation) else "hidden: more than 100 samples",
  interpretation = "distribution-level QC only; boxplots neither show nor exclude a batch effect",
  source_data = list(nf_rna_source_data(file.path(cfg$output_dir, "sample_distribution_qc.tsv"), cfg$output_dir, sprintf("box statistics (matrix = %s)", matrix_label)))
)
bfig <- nf_rna_figure_scope(fig, "boxplot")  # the context the sample boxplots are drawn with
for (spec in list(
  list(stem = "boxplot_log2_raw_counts", values = "raw", y = "log2(raw count + 1)", data = "log2(raw count + 1) of the L1-retained genes, before size-factor normalization", matrix = "log2_raw_count_plus_1"),
  list(stem = "boxplot_vst", values = "vst", y = "VST expression (blind, log2 scale)", data = "blind VST of the L1-retained genes (the matrix used for PCA and correlation)", matrix = "blind_vst"))) {
  box <- nf_rna_sample_boxplot(bfig, distribution_stats[distribution_stats$matrix == spec$matrix, , drop = FALSE], box_group, color_name, spec$y, box_batch, "batch")
  nf_rna_save_figure(bfig, box$plot, file.path(cfg$output_dir, spec$stem), box$width_mm %||% bfig$profile$double_col_mm, box$height_mm, cfg$output_dir, box_meta(spec$data, spec$y, spec$matrix, box))
}
nf_rna_write_figure_manifest(fig, cfg$output_dir)
flags <- character()
median_library <- median(library_size$input_total)
for (i in seq_len(nrow(library_size))) if (library_size$input_total[[i]] < median_library * 0.25) flags <- c(flags, paste0("LOW_LIBRARY_SIZE: ", library_size$sample_id[[i]], " is below 25% of the median input total."))
summary <- list(source_type = cfg$source_type, samples = samples, genes_input = nrow(source_counts), genes_removed_all_zero = sum(all_zero), genes_removed_low_total = sum(low_total), genes_retained = sum(keep), filter = cfg$filter, normalization = list(method = if (cfg$source_type == "salmon_tximport") "DESeq2 with tximport average-transcript-length normalization offsets" else "DESeq2 median-ratio size factors", scalar_size_factors_available = !all(is.na(size_factor_values))), vst = list(method = if (nrow(dds) < 1000) "DESeq2 varianceStabilizingTransformation (exact small-matrix path)" else "DESeq2 vst"), package_versions = list(R = R.version.string, DESeq2 = as.character(packageVersion("DESeq2")), tximport = if (cfg$source_type == "salmon_tximport") as.character(packageVersion("tximport")) else NULL), pca_proportion_variance = as.list(variance), flags = as.list(flags))
write(toJSON(summary, auto_unbox = TRUE, pretty = TRUE, null = "null"), file.path(cfg$output_dir, "backend_summary.json"))
nf_rna_write_provenance(cfg, cfg$output_dir, "L1", "SUCCESS", c("DESeq2", "tximport", "ggplot2", "ggrepel", "pheatmap", "scales", "systemfonts", "yaml", "jsonlite"), summary)
