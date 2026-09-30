args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2 || args[[1]] != "--config") stop("usage: l2_analysis.R --config CONFIG.json")
script_arg <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", script_arg[grep("^--file=", script_arg)][[1]])
source(file.path(dirname(normalizePath(script_file)), "provenance.R"))
source(file.path(dirname(normalizePath(script_file)), "design_metadata.R"))
source(file.path(dirname(normalizePath(script_file)), "sample_inputs.R"))
source(file.path(dirname(normalizePath(script_file)), "figures.R"))
suppressPackageStartupMessages({ library(jsonlite); library(DESeq2); library(ggplot2); library(pheatmap) })
cfg <- fromJSON(args[[2]], simplifyVector = FALSE)
dir.create(cfg$output_dir, recursive = TRUE, showWarnings = FALSE)
samples <- unlist(cfg$samples, use.names = FALSE)
metadata <- nf_rna_design_metadata(cfg, samples)
for (factor_name in unique(vapply(cfg$contrasts, function(item) item$factor, character(1)))) {
  if (!(factor_name %in% colnames(metadata))) stop(paste("contrast factor is absent from metadata:", factor_name))
  if (!is.null(cfg$design_variable_types) && identical(cfg$design_variable_types[[factor_name]], "continuous")) stop(paste("continuous variable cannot be used as contrast factor:", factor_name))
  metadata[[factor_name]] <- factor(metadata[[factor_name]])
}
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
if (sum(keep) < 2) stop("Filtering left fewer than two genes; L2 inference cannot proceed.")
dds <- dds[keep, ]
# The fit is deliberately performed once; every supported contrast is then extracted from it.
# Start with DESeq2 defaults.  The documented gene-wise fallback is only used
# for the exceptional case where DESeq2 says a dispersion trend cannot be fit.
fit_method <- "DESeq2::DESeq() default dispersion fit"
fit <- tryCatch(DESeq(dds), error = function(e) e)
if (inherits(fit, "error")) {
  if (!grepl("all gene-wise dispersion estimates", conditionMessage(fit), fixed = TRUE)) stop(fit)
  dds <- estimateSizeFactors(dds)
  dds <- estimateDispersionsGeneEst(dds)
  dispersions(dds) <- mcols(dds)$dispGeneEst
  dds <- nbinomWaldTest(dds)
  fit_method <- "DESeq2 gene-wise dispersion fallback after default trend-fit failure"
} else dds <- fit

# This is intentionally derived after DESeq2 has constructed and fitted the
# object.  It therefore records the actual typed colData and model matrix, not
# a parallel reconstruction in the control plane.
model_metadata <- as.data.frame(colData(dds))
model_variables <- all.vars(design(dds))
variable_details <- list()
for (variable in model_variables) {
  values <- model_metadata[[variable]]
  if (is.factor(values)) {
    variable_details[[variable]] <- list(type = "categorical", levels = as.character(levels(values)))
  } else if (is.numeric(values)) {
    numeric_values <- as.numeric(values)
    variable_details[[variable]] <- list(
      type = "continuous",
      n = length(numeric_values),
      min = min(numeric_values),
      max = max(numeric_values),
      mean = mean(numeric_values)
    )
  } else {
    stop(paste("fitted design variable has unsupported R representation:", variable))
  }
}
model_matrix <- model.matrix(design(dds), model_metadata)
design_details <- list(
  declared_variable_types = cfg$design_variable_types %||% list()
)
if (!is.null(cfg$pair_id)) {
  design_details$pair_id <- cfg$pair_id
}
design_details <- c(design_details, list(
  variables = variable_details,
  model_matrix = list(
    rank = qr(model_matrix)$rank,
    column_count = ncol(model_matrix),
    full_rank = qr(model_matrix)$rank == ncol(model_matrix),
    column_names = colnames(model_matrix)
  )
))
vst_table <- read.delim(cfg$l1_vst, check.names = FALSE, stringsAsFactors = FALSE)
rownames(vst_table) <- vst_table[[1]]
vst_matrix <- as.matrix(vst_table[, samples, drop = FALSE])
write_results <- function(df, path) write.table(df, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
fig <- nf_rna_figure_context(cfg)
for (item in cfg$contrasts) {
  contrast_id <- item$contrast_id; factor <- item$factor; numerator <- item$numerator; denominator <- item$denominator
  if (!(factor %in% names(metadata))) stop(paste("contrast factor absent from metadata:", factor))
  levels <- unique(metadata[[factor]])
  if (!(numerator %in% levels) || !(denominator %in% levels) || numerator == denominator) stop(paste("invalid contrast:", contrast_id))
  directory <- file.path(cfg$output_dir, "contrasts", contrast_id); dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  result <- results(dds, contrast = c(factor, numerator, denominator), independentFiltering = TRUE)
  result_df <- data.frame(gene_id = rownames(result), as.data.frame(result), check.names = FALSE)
  result_df <- result_df[, c("gene_id", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")]
  write_results(result_df, file.path(directory, "all_genes.tsv"))
  significant <- !is.na(result_df$padj) & !is.na(result_df$log2FoldChange) & result_df$padj < cfg$thresholds$padj & abs(result_df$log2FoldChange) >= cfg$thresholds$abs_log2fc
  significant_df <- result_df[significant, , drop = FALSE]
  up_df <- significant_df[significant_df$log2FoldChange > 0, , drop = FALSE]
  down_df <- significant_df[significant_df$log2FoldChange < 0, , drop = FALSE]
  write_results(significant_df, file.path(directory, "significant.tsv")); write_results(up_df, file.path(directory, "upregulated.tsv")); write_results(down_df, file.path(directory, "downregulated.tsv"))
  # Display policy only: padj == 0 (underflow) is drawn at the smallest positive
  # double and marked as capped; the tables keep the reported padj.
  padj_floor <- .Machine$double.xmin
  direction <- ifelse(significant & result_df$log2FoldChange > 0, "Up", ifelse(significant & result_df$log2FoldChange < 0, "Down", "Not significant"))
  plot_df <- data.frame(lfc = result_df$log2FoldChange, plot_padj = ifelse(is.na(result_df$padj), NA_real_, pmax(result_df$padj, padj_floor)), group = direction,
                        capped = !is.na(result_df$padj) & result_df$padj < padj_floor)
  plot_df <- plot_df[order(plot_df$group != "Not significant"), , drop = FALSE]  # significant points drawn on top, no transparency
  drawable <- is.finite(result_df$log2FoldChange) & !is.na(result_df$padj)
  counts <- table(factor(direction[drawable], levels = names(NF_RNA_DIRECTION_COLOURS)))
  direction_labels <- c(Up = sprintf("Higher in %s (%d)", numerator, counts[["Up"]]), Down = sprintf("Higher in %s (%d)", denominator, counts[["Down"]]), "Not significant" = sprintf("Not significant (%d)", counts[["Not significant"]]))
  x_extent <- max(c(abs(plot_df$lfc[is.finite(plot_df$lfc)]), cfg$thresholds$abs_log2fc, 1), na.rm = TRUE) * 1.05
  line <- nf_rna_linewidth(0.5)
  p <- ggplot(plot_df, aes(x = lfc, y = -log10(plot_padj), color = group, shape = capped)) +
    geom_hline(yintercept = -log10(cfg$thresholds$padj), linetype = "dashed", linewidth = line, colour = "grey40") +
    (if (cfg$thresholds$abs_log2fc > 0) geom_vline(xintercept = c(-1, 1) * cfg$thresholds$abs_log2fc, linetype = "dashed", linewidth = line, colour = "grey40")) +
    geom_point(na.rm = TRUE, size = 0.8) +
    scale_color_manual(values = NF_RNA_DIRECTION_COLOURS, labels = direction_labels, breaks = names(NF_RNA_DIRECTION_COLOURS), name = NULL) +
    scale_shape_manual(values = c("FALSE" = 16, "TRUE" = 17), breaks = "TRUE", labels = "padj = 0, drawn at display floor", name = NULL, guide = if (any(plot_df$capped)) "legend" else "none") +
    scale_x_continuous(limits = c(-x_extent, x_extent)) +
    labs(x = paste0("log2 fold change (", numerator, " / ", denominator, ")"), y = "-log10(adjusted P, BH)") + nf_rna_theme(fig) + theme(legend.position = "bottom")
  nf_rna_save_figure(fig, p, file.path(directory, "volcano"), fig$profile$double_col_mm, 120, cfg$output_dir, list(
    kind = "volcano", contrast = list(contrast_id = contrast_id, factor = factor, numerator = numerator, denominator = denominator, direction = "log2FC > 0: higher in numerator"),
    statistics = list(fold_change = "DESeq2 unshrunken log2FoldChange", adjusted_p = "DESeq2 padj (Benjamini-Hochberg, independent filtering)"),
    thresholds = list(padj_below = cfg$thresholds$padj, abs_log2fc_at_least = cfg$thresholds$abs_log2fc, lines = "dashed lines at the significance thresholds"),
    display = list(padj_zero_floor = padj_floor, capped_points = sum(plot_df$capped), x_axis = "symmetric about 0", genes_without_padj = "not drawn", labels = "none; no genes are labelled", legend_counts = "drawn genes only (finite log2FC and non-NA padj)", genes_not_drawn = sum(!drawable)),
    counts = as.list(counts), colours = as.list(NF_RNA_DIRECTION_COLOURS),
    source_data = list(nf_rna_source_data(file.path(directory, "all_genes.tsv"), cfg$output_dir, "all DESeq2 results for this contrast"))
  ))
  heatmap_status <- "NOT_APPLICABLE"; heatmap_reason <- "no significant genes"
  if (nrow(significant_df) > 0) {
    ordered <- significant_df[order(significant_df$padj, significant_df$gene_id), , drop = FALSE]
    selected <- head(ordered, cfg$heatmap_top_n)
    write.table(data.frame(gene_id = selected$gene_id), file.path(directory, "heatmap_genes.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
    heatmap_samples <- samples[metadata[samples, factor] %in% c(numerator, denominator)]
    heatmap_values <- vst_matrix[selected$gene_id, heatmap_samples, drop = FALSE]
    # Undefined z-scores (missing values, constant rows) stay NA; they are never
    # drawn as 0.  Rows with no defined z-score are left out of the drawing and
    # listed in heatmap_zscores.tsv/figure_manifest.json.
    z <- nf_rna_row_zscore(heatmap_values)
    write_results(data.frame(gene_id = rownames(z), z, check.names = FALSE), file.path(directory, "heatmap_zscores.tsv"))
    undefined <- attr(z, "undefined_rows")
    drawn <- z[!(rownames(z) %in% undefined), , drop = FALSE]
    if (nrow(drawn) > 0) {
      z_limit <- 2
      display <- nf_rna_saturate(drawn, z_limit)
      row_tree <- nf_rna_cluster(drawn); column_tree <- nf_rna_cluster(t(drawn))
      group_colours <- nf_rna_categorical_colors(levels(metadata[[factor]]))  # same levels and order as the L1 PCA
      annotation <- data.frame(metadata[heatmap_samples, factor], row.names = heatmap_samples, check.names = FALSE); names(annotation) <- factor
      annotation_colours <- setNames(list(group_colours[intersect(names(group_colours), unique(annotation[[factor]]))]), factor)
      profile <- fig$profile
      overhead_mm <- 40 + nf_rna_pt_to_mm(profile$small_pt) * 0.6 * max(nchar(heatmap_samples))
      row_pt <- profile$small_pt
      fits <- function(pt) overhead_mm + nrow(display) * nf_rna_pt_to_mm(pt) * 1.25 <= profile$max_height_mm
      while (!fits(row_pt) && row_pt > profile$min_pt) row_pt <- row_pt - 0.5
      show_rows <- fits(row_pt)
      height_mm <- if (show_rows) max(80, overhead_mm + nrow(display) * nf_rna_pt_to_mm(row_pt) * 1.25) else profile$max_height_mm
      legend_at <- seq(-z_limit, z_limit, by = 1)
      heat <- pheatmap(display, color = nf_rna_diverging_palette(100), breaks = seq(-z_limit, z_limit, length.out = 101), na_col = NF_RNA_NA_COLOUR,
        cluster_rows = if (is.null(row_tree)) FALSE else row_tree, cluster_cols = if (is.null(column_tree)) FALSE else column_tree,
        annotation_col = annotation, annotation_colors = annotation_colours, border_color = NA, show_rownames = show_rows,
        legend_breaks = legend_at, legend_labels = c(paste0("\u2264", -z_limit), legend_at[-c(1, length(legend_at))], paste0("\u2265", z_limit)),
        main = if (nrow(display) == nrow(selected)) sprintf("%s: row z-score of VST, top %d genes by padj", contrast_id, nrow(selected)) else sprintf("%s: row z-score of VST, top %d genes by padj (%d drawn)", contrast_id, nrow(selected), nrow(display)),
        fontsize = profile$base_pt, fontsize_row = row_pt, fontsize_col = profile$small_pt, silent = TRUE)
      nf_rna_save_figure(fig, heat, file.path(directory, "heatmap"), profile$double_col_mm, height_mm, cfg$output_dir, list(
        kind = "DEG heatmap", contrast = list(contrast_id = contrast_id, factor = factor, numerator = numerator, denominator = denominator),
        display_subset = list(rule = sprintf("top %d significant genes by padj, ties by gene_id", cfg$heatmap_top_n), selected = nrow(selected), drawn = nrow(display), significant_total = nrow(significant_df),
                              samples = "samples in the numerator and denominator groups"),
        transformation = list(values = "blind VST from L1", row_zscore = "per gene across the drawn samples: (x - mean) / sd (n - 1)",
                              undefined_policy = "missing values and rows with <2 finite values or zero variance stay NA; rows with no defined z-score are not drawn",
                              undefined_rows = as.list(undefined), missing_cells_drawn = sum(is.na(drawn)), na_colour = NF_RNA_NA_COLOUR,
                              colour_saturation = sprintf("visual only: colours saturate at |z| >= %d (legend \u2264-%d / \u2265%d); heatmap_zscores.tsv keeps the unclipped values", z_limit, z_limit, z_limit)),
        clustering = list(rows = if (is.null(row_tree)) "off (undefined distances or <2 rows)" else "complete linkage, Euclidean distance on unclipped z-scores",
                          columns = if (is.null(column_tree)) "off (undefined distances or <2 samples)" else "complete linkage, Euclidean distance on unclipped z-scores"),
        colour_scale = list(kind = "diverging, midpoint 0", limits = c(-z_limit, z_limit), colours = as.list(NF_RNA_DIVERGING_ENDS)), group_colours = as.list(group_colours),
        row_labels = if (show_rows) sprintf("%.1f pt gene IDs", row_pt) else "hidden: too many rows to label legibly at this height",
        source_data = list(nf_rna_source_data(file.path(directory, "heatmap_zscores.tsv"), cfg$output_dir, "z-scores (unclipped, NA preserved)"),
                           nf_rna_source_data(file.path(directory, "heatmap_genes.tsv"), cfg$output_dir, "selected genes"),
                           nf_rna_source_data(file.path(directory, "significant.tsv"), cfg$output_dir, "complete significant result set"))
      ))
      heatmap_status <- "AVAILABLE"; heatmap_reason <- NULL
    } else heatmap_reason <- "no selected gene has a defined row z-score"
  }
  tested <- sum(!is.na(result_df$pvalue))
  summary <- list(input_genes = nrow(source_counts), filtered_genes = sum(!keep), retained_genes = sum(keep), tested_genes = tested, pvalue_na = sum(is.na(result_df$pvalue)), padj_na = sum(is.na(result_df$padj)), independent_filtering = TRUE, multiple_testing_method = "Benjamini-Hochberg (DESeq2 default)", fit_method = fit_method, heatmap_status = heatmap_status, heatmap_top_n = cfg$heatmap_top_n)
  if (!is.null(heatmap_reason) && nrow(significant_df) > 0) summary$heatmap_reason <- heatmap_reason
  write(toJSON(summary, auto_unbox = TRUE, pretty = TRUE), file.path(directory, "backend_summary.json"))
}
nf_rna_write_figure_manifest(fig, cfg$output_dir)
nf_rna_write_provenance(
  cfg, cfg$output_dir, "L2", "SUCCESS",
  c("DESeq2", "tximport", "ggplot2", "pheatmap", "scales", "systemfonts", "jsonlite"),
  list(contrast_count = length(cfg$contrasts), independent_filtering = TRUE, fit_method = fit_method),
  design_details
)
