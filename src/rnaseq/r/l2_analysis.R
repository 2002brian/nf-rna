args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2 || args[[1]] != "--config") stop("usage: l2_analysis.R --config CONFIG.json")
script_arg <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", script_arg[grep("^--file=", script_arg)][[1]])
source(file.path(dirname(normalizePath(script_file)), "provenance.R"))
source(file.path(dirname(normalizePath(script_file)), "design_metadata.R"))
source(file.path(dirname(normalizePath(script_file)), "sample_inputs.R"))
source(file.path(dirname(normalizePath(script_file)), "figures.R"))
source(file.path(dirname(normalizePath(script_file)), "gene_annotation.R"))
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
# Descriptive outputs of the fitted object.  None of them enters the DEG definition.
normalized_counts <- counts(dds, normalized = TRUE)[, samples, drop = FALSE]  # size factors of the complete fitted dds
annotation_result <- nf_rna_gene_annotation(rownames(dds), cfg$gene_annotation %||% list())
gene_annotation <- annotation_result$table
write_results(data.frame(gene_annotation[, c("gene_id", "symbol")], normalized_counts, check.names = FALSE), file.path(cfg$output_dir, "normalized_counts.tsv"))
gene_label <- function(ids) { symbol <- gene_annotation$symbol[match(ids, gene_annotation$gene_id)]; ifelse(is.na(symbol) | symbol == "", ids, symbol) }
same_values <- function(a, b) identical(is.na(a), is.na(b)) && isTRUE(all.equal(a[!is.na(a)], b[!is.na(b)], tolerance = 1e-4))  # IRLS convergence noise of a reparametrized refit is ~1e-6
# Display-only LFC shrinkage.  apeglm on an existing coefficient; otherwise
# apeglm on a releveled copy whose raw estimates must reproduce the original
# fit; otherwise ashr on the contrast.  A failure leaves NA, never raw LFC.
shrink_lfc <- function(result, factor, numerator, denominator) {
  coef <- paste0(factor, "_", numerator, "_vs_", denominator)
  info <- list(method = NULL, coef = NULL, contrast = NULL, releveled = FALSE, relevel_verification = NULL, status = "failed", error = NULL)
  shrunk <- tryCatch({
    if (coef %in% resultsNames(dds)) {
      info$method <- "apeglm"; info$coef <- coef
      lfcShrink(dds, coef = coef, type = "apeglm", quiet = TRUE)
    } else {
      copy <- dds
      copy[[factor]] <- relevel(copy[[factor]], ref = denominator)
      copy <- nbinomWaldTest(copy, quiet = TRUE)
      check <- if (coef %in% resultsNames(copy)) results(copy, name = coef) else NULL
      verified <- !is.null(check) && same_values(check$log2FoldChange, result$log2FoldChange) && same_values(check$pvalue, result$pvalue)
      info$relevel_verification <- if (is.null(check)) "coefficient absent after relevel; copy not used" else if (verified) "raw log2FoldChange and pvalue of the releveled copy match the original fit (all.equal mean relative difference <= 1e-4)" else "raw log2FoldChange or pvalue of the releveled copy differ from the original fit; copy not used"
      if (verified) {
        info$method <- "apeglm"; info$coef <- coef; info$releveled <- TRUE
        lfcShrink(copy, coef = coef, type = "apeglm", quiet = TRUE)
      } else {
        info$method <- "ashr"; info$contrast <- c(factor, numerator, denominator)
        lfcShrink(dds, contrast = c(factor, numerator, denominator), type = "ashr", quiet = TRUE)
      }
    }
  }, error = function(e) { info$error <<- conditionMessage(e); NULL })
  lfc <- rep(NA_real_, nrow(result))
  if (!is.null(shrunk)) { lfc <- as.numeric(shrunk[rownames(result), "log2FoldChange"]); info$status <- "success" }
  list(lfc = lfc, info = info)
}
shrinkage_versions <- nf_rna_package_versions(c("DESeq2", "apeglm", "ashr"))
fig <- nf_rna_figure_context(cfg)
hfig <- nf_rna_figure_scope(fig, "heatmap")  # the context DEG heatmaps are drawn with (standard preset when requested)
vfig <- nf_rna_figure_scope(fig, "volcano")  # the context volcano and MA plots are drawn with (standard preset when requested)
deg_figure_formats <- c("png", "png_raster", "tiff", "pdf")  # adds a raster-dpi PNG beside the report preview
volcano_top_labels <- 10
for (item in cfg$contrasts) {
  contrast_id <- item$contrast_id; factor <- item$factor; numerator <- item$numerator; denominator <- item$denominator
  if (!(factor %in% names(metadata))) stop(paste("contrast factor absent from metadata:", factor))
  levels <- unique(metadata[[factor]])
  if (!(numerator %in% levels) || !(denominator %in% levels) || numerator == denominator) stop(paste("invalid contrast:", contrast_id))
  directory <- file.path(cfg$output_dir, "contrasts", contrast_id); dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  # Independent filtering maximises rejections at padj < alpha, so alpha is the
  # declared DEG FDR threshold (DESeq2 vignette); the DESeq2 default 0.1 is not used.
  result <- results(dds, contrast = c(factor, numerator, denominator), independentFiltering = TRUE, alpha = cfg$thresholds$padj)
  filter_threshold <- unname(S4Vectors::metadata(result)$filterThreshold)
  result_df <- data.frame(gene_id = rownames(result), as.data.frame(result), check.names = FALSE)
  result_df <- result_df[, c("gene_id", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")]
  stopifnot(identical(result_df$gene_id, gene_annotation$gene_id))
  shrink <- shrink_lfc(result, factor, numerator, denominator)
  level_of <- as.character(colData(dds)[[factor]])
  # Descriptive group means of normalized counts, not model estimates.
  deg_df <- data.frame(gene_annotation[, c("gene_id", "symbol", "description")], result_df[, -1, drop = FALSE],
                       mean_normcount_numerator = unname(rowMeans(normalized_counts[, level_of == numerator, drop = FALSE])),
                       mean_normcount_denominator = unname(rowMeans(normalized_counts[, level_of == denominator, drop = FALSE])),
                       log2FoldChange_shrunk = shrink$lfc, check.names = FALSE)
  write_results(deg_df, file.path(directory, "all_genes.tsv"))
  # The DEG definition is unchanged: padj and the raw (unshrunken) log2FoldChange.
  significant <- !is.na(result_df$padj) & !is.na(result_df$log2FoldChange) & result_df$padj < cfg$thresholds$padj & abs(result_df$log2FoldChange) >= cfg$thresholds$abs_log2fc
  significant_df <- deg_df[significant, , drop = FALSE]
  up_df <- significant_df[significant_df$log2FoldChange > 0, , drop = FALSE]
  down_df <- significant_df[significant_df$log2FoldChange < 0, , drop = FALSE]
  write_results(significant_df, file.path(directory, "significant.tsv")); write_results(up_df, file.path(directory, "upregulated.tsv")); write_results(down_df, file.path(directory, "downregulated.tsv"))
  # Display policy only: padj == 0 (underflow) is drawn at the smallest positive
  # double and marked as capped; the tables keep the reported padj.
  padj_floor <- .Machine$double.xmin
  direction <- ifelse(significant & result_df$log2FoldChange > 0, "Up", ifelse(significant & result_df$log2FoldChange < 0, "Down", "Not significant"))
  # Figures show the shrunken LFC; if shrinkage failed they fall back to the raw
  # LFC and say so.  Colour (significance) always follows the raw DEG definition.
  shrunk_shown <- identical(shrink$info$status, "success")
  shown_lfc <- if (shrunk_shown) deg_df$log2FoldChange_shrunk else result_df$log2FoldChange
  lfc_name <- if (shrunk_shown) sprintf("shrunken log2 fold change (%s)", shrink$info$method) else "raw log2 fold change (shrinkage failed)"
  deg_rule <- sprintf("DEG: padj < %s and |raw log2FC| \u2265 %s; LFC shrinkage: %s", cfg$thresholds$padj, cfg$thresholds$abs_log2fc, if (shrunk_shown) shrink$info$method else "failed, raw LFC shown")
  drawable <- is.finite(shown_lfc) & !is.na(result_df$padj)
  by_padj <- order(result_df$padj, result_df$gene_id)
  label_ids <- head(result_df$gene_id[by_padj][(significant & drawable)[by_padj]], volcano_top_labels)
  plot_df <- data.frame(lfc = shown_lfc, plot_padj = ifelse(is.na(result_df$padj), NA_real_, pmax(result_df$padj, padj_floor)), group = direction,
                        capped = !is.na(result_df$padj) & result_df$padj < padj_floor,
                        label = ifelse(result_df$gene_id %in% label_ids, gene_label(result_df$gene_id), NA_character_))
  plot_df <- plot_df[order(plot_df$group != "Not significant"), , drop = FALSE]  # significant points drawn on top, no transparency
  counts <- table(factor(direction[drawable], levels = names(nf_rna_direction_colours(vfig))))
  direction_labels <- c(Up = sprintf("Higher in %s (%d)", numerator, counts[["Up"]]), Down = sprintf("Higher in %s (%d)", denominator, counts[["Down"]]), "Not significant" = sprintf("Not significant (%d)", counts[["Not significant"]]))
  x_extent <- max(c(abs(plot_df$lfc[is.finite(plot_df$lfc)]), if (!shrunk_shown) cfg$thresholds$abs_log2fc, 1), na.rm = TRUE) * 1.05
  line <- nf_rna_linewidth(0.5)
  volcano_width_mm <- nf_rna_layout_width_mm(vfig, "volcano", fig$profile$double_col_mm)
  volcano_height_mm <- if (is.na(nf_rna_layout_aspect(vfig, "volcano"))) 120 else min(vfig$profile$max_height_mm, volcano_width_mm * nf_rna_layout_aspect(vfig, "volcano"))
  p <- ggplot(plot_df, aes(x = lfc, y = -log10(plot_padj), color = group, shape = capped)) +
    geom_hline(yintercept = -log10(cfg$thresholds$padj), linetype = "dashed", linewidth = line, colour = "grey40") +
    # |log2FC| threshold lines belong to the raw scale only; on the shrunken axis they would misstate the rule.
    (if (!shrunk_shown && cfg$thresholds$abs_log2fc > 0) geom_vline(xintercept = c(-1, 1) * cfg$thresholds$abs_log2fc, linetype = "dashed", linewidth = line, colour = "grey40")) +
    geom_point(na.rm = TRUE, size = 0.8, alpha = nf_rna_point_alpha(vfig)) +
    ggrepel::geom_text_repel(data = plot_df[!is.na(plot_df$label), , drop = FALSE], aes(label = label), na.rm = TRUE, show.legend = FALSE, seed = 1,
                             size = nf_rna_pt_to_mm(nf_rna_pca_label_pt(vfig)), family = vfig$font$family, max.overlaps = Inf, min.segment.length = 0) +
    scale_color_manual(values = nf_rna_direction_colours(vfig), labels = direction_labels, breaks = names(nf_rna_direction_colours(vfig)), name = NULL) +
    scale_shape_manual(values = c("FALSE" = 16, "TRUE" = 17), breaks = "TRUE", labels = "padj = 0, drawn at display floor", name = NULL, guide = if (any(plot_df$capped)) "legend" else "none") +
    scale_x_continuous(limits = c(-x_extent, x_extent)) +
    labs(x = paste0(lfc_name, " (", numerator, " / ", denominator, ")"), y = "-log10(adjusted P, BH)", caption = deg_rule) + nf_rna_theme(vfig) + theme(legend.position = "bottom") +
    (if (!is.null(nf_rna_layout_value(vfig, "volcano", "legend_ncol"))) theme(legend.key.height = grid::unit(3, "mm"), legend.key.spacing.y = grid::unit(0.5, "mm"))) +
    (if (!is.null(nf_rna_layout_value(vfig, "volcano", "legend_ncol"))) guides(color = guide_legend(ncol = nf_rna_layout_value(vfig, "volcano", "legend_ncol"), order = 1), shape = guide_legend(order = 2)))
  nf_rna_save_figure(vfig, p, file.path(directory, "volcano"), volcano_width_mm, volcano_height_mm, cfg$output_dir, list(
    kind = "volcano", contrast = list(contrast_id = contrast_id, factor = factor, numerator = numerator, denominator = denominator, direction = "log2FC > 0: higher in numerator"),
    statistics = list(fold_change = if (shrunk_shown) sprintf("log2FoldChange_shrunk (%s); display only", shrink$info$method) else "DESeq2 unshrunken log2FoldChange (fallback: shrinkage failed)",
                      significance_colour = "raw DESeq2 log2FoldChange and padj (unchanged DEG definition)", adjusted_p = "DESeq2 padj (Benjamini-Hochberg, independent filtering)"),
    thresholds = list(padj_below = cfg$thresholds$padj, abs_log2fc_at_least = cfg$thresholds$abs_log2fc,
                      lines = if (shrunk_shown) "dashed padj line only; the |log2FC| threshold applies to the raw LFC and is stated in the caption" else "dashed lines at the significance thresholds"),
    display = list(padj_zero_floor = padj_floor, capped_points = sum(plot_df$capped), x_axis = "symmetric about 0", genes_without_padj = "not drawn",
                   labels = sprintf("top %d significant drawn genes by padj (ties by gene_id), labelled by symbol or gene_id when symbol is NA", volcano_top_labels),
                   labelled_genes = as.list(label_ids), legend_counts = "drawn genes only (finite shown LFC and non-NA padj)", genes_not_drawn = sum(!drawable)),
    counts = as.list(counts), colours = as.list(nf_rna_direction_colours(vfig)),
    source_data = list(nf_rna_source_data(file.path(directory, "all_genes.tsv"), cfg$output_dir, "all DESeq2 results for this contrast"))
  ), formats = deg_figure_formats)
  ma_df <- data.frame(x = log10(result_df$baseMean), y = shown_lfc, group = direction)
  ma_drawable <- is.finite(ma_df$x) & is.finite(ma_df$y)
  ma_df <- ma_df[ma_drawable, , drop = FALSE]
  ma_df <- ma_df[order(ma_df$group != "Not significant"), , drop = FALSE]
  ma_labels <- c(Up = sprintf("Higher in %s", numerator), Down = sprintf("Higher in %s", denominator), "Not significant" = "Not significant")
  p_ma <- ggplot(ma_df, aes(x = x, y = y, color = group)) +
    geom_hline(yintercept = 0, linewidth = line, colour = "grey40") +
    geom_point(size = 0.8, alpha = nf_rna_point_alpha(vfig)) +
    scale_color_manual(values = nf_rna_direction_colours(vfig), labels = ma_labels, breaks = names(nf_rna_direction_colours(vfig)), name = NULL) +
    labs(x = "log10(baseMean)", y = paste0(lfc_name, "\n(", numerator, " / ", denominator, ")"), caption = deg_rule) + nf_rna_theme(vfig) + theme(legend.position = "bottom")
  nf_rna_save_figure(vfig, p_ma, file.path(directory, "ma_plot"), volcano_width_mm, volcano_height_mm, cfg$output_dir, list(
    kind = "MA plot", contrast = list(contrast_id = contrast_id, factor = factor, numerator = numerator, denominator = denominator, direction = "log2FC > 0: higher in numerator"),
    statistics = list(x = "log10(DESeq2 baseMean)", y = if (shrunk_shown) sprintf("log2FoldChange_shrunk (%s); display only", shrink$info$method) else "DESeq2 unshrunken log2FoldChange (fallback: shrinkage failed)",
                      significance_colour = "raw DESeq2 log2FoldChange and padj (unchanged DEG definition)"),
    display = list(genes_drawn = nrow(ma_df), genes_not_drawn = sum(!ma_drawable), genes_not_drawn_rule = "non-finite log10(baseMean) or shown LFC"),
    colours = as.list(nf_rna_direction_colours(vfig)),
    source_data = list(nf_rna_source_data(file.path(directory, "all_genes.tsv"), cfg$output_dir, "all DESeq2 results for this contrast"))
  ), formats = deg_figure_formats)
  heatmap_status <- "NOT_APPLICABLE"; heatmap_reason <- "no significant genes"
  if (nrow(significant_df) > 0) {
    ordered <- significant_df[order(significant_df$padj, significant_df$gene_id), , drop = FALSE]
    selected <- head(ordered, cfg$heatmap_top_n)
    write_results(data.frame(gene_id = selected$gene_id, symbol = selected$symbol), file.path(directory, "heatmap_genes.tsv"))
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
      # Drawn row labels: symbol, else gene_id; a repeated symbol is disambiguated with its gene_id.
      row_labels <- gene_label(rownames(display))
      repeated <- duplicated(row_labels) | duplicated(row_labels, fromLast = TRUE)
      row_labels[repeated] <- paste0(row_labels[repeated], " (", rownames(display)[repeated], ")")
      labelled <- display; rownames(labelled) <- row_labels
      row_tree <- nf_rna_cluster(drawn); column_tree <- nf_rna_cluster(t(drawn))
      group_colours <- nf_rna_group_colours(levels(metadata[[factor]]), context = hfig)  # the same colour per group as the L1 PCA and boxplots
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
      standard_heatmap <- identical(hfig$preset, "standard")
      heat_title <- nf_rna_heatmap_title(hfig, numerator, denominator, width_mm = if (standard_heatmap) nf_rna_single_column_mm(hfig) else profile$double_col_mm,
        detail = if (nrow(display) == nrow(selected)) sprintf("Top %d DEGs \u00b7 row z-score of VST", nrow(selected)) else sprintf("Top %d DEGs (%d drawn) \u00b7 row z-score of VST", nrow(selected), nrow(display)))
      heat_args <- list(mat = labelled, color = nf_rna_diverging_palette(100, hfig), breaks = seq(-z_limit, z_limit, length.out = 101), na_col = nf_rna_na_colour(hfig),
        cluster_rows = if (is.null(row_tree)) FALSE else row_tree, cluster_cols = if (is.null(column_tree)) FALSE else column_tree,
        annotation_col = annotation, annotation_colors = annotation_colours, border_color = NA,
        legend_breaks = legend_at, legend_labels = c(paste0("\u2264", -z_limit), legend_at[-c(1, length(legend_at))], paste0("\u2265", z_limit)),
        main = heat_title, fontsize = profile$base_pt, fontsize_col = profile$small_pt, silent = TRUE)
      if (!standard_heatmap) {
        heat <- do.call(pheatmap, c(heat_args, list(show_rownames = show_rows, fontsize_row = row_pt)))
      } else {
        # Content-aware size: absolute cells, measured labels/trees/legend, never above the profile maximum.
        fit <- nf_rna_deg_heatmap(hfig, labelled, annotation, annotation_colours, row_tree, column_tree, heat_title, z_limit)
        heat <- fit$heat; show_rows <- fit$show_rows; row_pt <- fit$row_pt; height_mm <- fit$height_mm
      }
      if (!standard_heatmap) heat <- nf_rna_span_pheatmap_title(heat)
      heat_width_mm <- if (standard_heatmap) fit$width_mm else profile$double_col_mm
      nf_rna_save_figure(hfig, heat, file.path(directory, "heatmap"), heat_width_mm, height_mm, cfg$output_dir, c(list(
        kind = "DEG heatmap", title = heat_title, title_layout = "contrast line(s) over a detail line, centred across the full figure width",
        contrast = list(contrast_id = contrast_id, factor = factor, numerator = numerator, denominator = denominator),
        display_subset = list(rule = sprintf("top %d significant genes by padj, ties by gene_id", cfg$heatmap_top_n), selected = nrow(selected), drawn = nrow(display), significant_total = nrow(significant_df),
                              samples = "samples in the numerator and denominator groups"),
        transformation = list(values = "blind VST from L1", row_zscore = "per gene across the drawn samples: (x - mean) / sd (n - 1)",
                              undefined_policy = "missing values and rows with <2 finite values or zero variance stay NA; rows with no defined z-score are not drawn",
                              undefined_rows = as.list(undefined), missing_cells_drawn = sum(is.na(drawn)), na_colour = nf_rna_na_colour(fig),
                              colour_saturation = sprintf("visual only: colours saturate at |z| >= %d (legend \u2264-%d / \u2265%d); heatmap_zscores.tsv keeps the unclipped values", z_limit, z_limit, z_limit)),
        clustering = list(rows = if (is.null(row_tree)) "off (undefined distances or <2 rows)" else "complete linkage, Euclidean distance on unclipped z-scores",
                          columns = if (is.null(column_tree)) "off (undefined distances or <2 samples)" else "complete linkage, Euclidean distance on unclipped z-scores"),
        colour_scale = list(kind = "diverging, midpoint 0", limits = c(-z_limit, z_limit), colours = as.list(nf_rna_diverging_ends(fig))), group_colours = as.list(group_colours),
        row_labels = if (show_rows) sprintf("%.1f pt gene symbols (gene_id when the symbol is NA; gene_id appended to a repeated symbol)", row_pt) else "hidden: too many rows to label legibly at this height",
        source_data = list(nf_rna_source_data(file.path(directory, "heatmap_zscores.tsv"), cfg$output_dir, "z-scores (unclipped, NA preserved)"),
                           nf_rna_source_data(file.path(directory, "heatmap_genes.tsv"), cfg$output_dir, "selected genes"),
                           nf_rna_source_data(file.path(directory, "significant.tsv"), cfg$output_dir, "complete significant result set"))
      ), if (standard_heatmap) list(
        layout = list(cell_width_mm = fit$cw, cell_height_mm = fit$ch, sample_labels = if (fit$show_cols) "shown" else "hidden: too many samples to label legibly within the profile maximum"),
        display_order = list(rows = if (is.null(row_tree)) rownames(display) else rownames(display)[heat$tree_row$order],
                             columns = if (is.null(column_tree)) colnames(display) else colnames(display)[heat$tree_col$order]))), formats = deg_figure_formats)
      heatmap_status <- "AVAILABLE"; heatmap_reason <- NULL
    } else heatmap_reason <- "no selected gene has a defined row z-score"
  }
  tested <- sum(!is.na(result_df$pvalue))
  summary <- list(input_genes = nrow(source_counts), filtered_genes = sum(!keep), retained_genes = sum(keep), tested_genes = tested, pvalue_na = sum(is.na(result_df$pvalue)), padj_na = sum(is.na(result_df$padj)), independent_filtering = TRUE, independent_filtering_alpha = cfg$thresholds$padj, independent_filtering_threshold = filter_threshold, multiple_testing_method = "Benjamini-Hochberg (DESeq2 default)", fit_method = fit_method, heatmap_status = heatmap_status, heatmap_top_n = cfg$heatmap_top_n)
  if (!is.null(heatmap_reason) && nrow(significant_df) > 0) summary$heatmap_reason <- heatmap_reason
  write(toJSON(summary, auto_unbox = TRUE, pretty = TRUE), file.path(directory, "backend_summary.json"))
  contrast_metadata <- list(
    contrast_id = contrast_id, factor = factor,
    numerator = list(level = numerator, n_samples = sum(level_of == numerator)),
    denominator = list(level = denominator, n_samples = sum(level_of == denominator)),
    design_formula = paste(deparse(design(dds)), collapse = ""), fit_method = fit_method,
    deg_definition = list(padj_below = cfg$thresholds$padj, abs_log2fc_at_least = cfg$thresholds$abs_log2fc,
                          rule = "padj < threshold AND abs(raw log2FoldChange) >= threshold; log2FoldChange_shrunk is not used",
                          testing = "Wald test of LFC = 0, then effect-size filtering (no lfcThreshold)",
                          independent_filtering = list(alpha = cfg$thresholds$padj, filter_statistic = "baseMean", filter_threshold = filter_threshold,
                                                       rule = "DESeq2 results(alpha = padj threshold): the baseMean quantile maximising rejections at padj < alpha")),
    shrinkage = c(shrink$info, list(column = "log2FoldChange_shrunk", use = "display and ranking only")),
    group_means = list(columns = c("mean_normcount_numerator", "mean_normcount_denominator"),
                       definition = "arithmetic mean of counts(dds, normalized = TRUE) over the analysed samples of each level; descriptive, not a model estimate"),
    annotation = annotation_result$info,
    software = shrinkage_versions
  )
  write(toJSON(contrast_metadata, auto_unbox = TRUE, pretty = TRUE, null = "null", na = "null"), file.path(directory, "contrast_metadata.json"))
}
nf_rna_write_figure_manifest(fig, cfg$output_dir)
nf_rna_write_provenance(
  cfg, cfg$output_dir, "L2", "SUCCESS",
  c("DESeq2", "apeglm", "ashr", "AnnotationDbi", cfg$gene_annotation$orgdb_package, "tximport", "ggplot2", "ggrepel", "pheatmap", "scales", "systemfonts", "jsonlite"),
  list(contrast_count = length(cfg$contrasts), independent_filtering = TRUE, fit_method = fit_method, gene_annotation = annotation_result$info),
  design_details
)
