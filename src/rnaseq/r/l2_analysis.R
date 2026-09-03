args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2 || args[[1]] != "--config") stop("usage: l2_analysis.R --config CONFIG.json")
suppressPackageStartupMessages({ library(jsonlite); library(DESeq2); library(ggplot2); library(pheatmap) })
cfg <- fromJSON(args[[2]], simplifyVector = FALSE)
dir.create(cfg$output_dir, recursive = TRUE, showWarnings = FALSE)
metadata <- read.csv(cfg$metadata, check.names = FALSE, stringsAsFactors = FALSE)
rownames(metadata) <- metadata$sample_id
samples <- unlist(cfg$samples, use.names = FALSE)
metadata <- metadata[samples, , drop = FALSE]
formula <- as.formula(cfg$formula)
if (cfg$source_type == "raw_counts") {
  table <- read.csv(cfg$counts, check.names = FALSE, stringsAsFactors = FALSE)
  ids <- table[[1]]
  matrix_counts <- as.matrix(table[, samples, drop = FALSE])
  rownames(matrix_counts) <- ids
  storage.mode(matrix_counts) <- "numeric"
  dds <- DESeqDataSetFromMatrix(countData = matrix_counts, colData = metadata, design = formula)
  source_counts <- matrix_counts
} else if (cfg$source_type == "salmon_tximport") {
  suppressPackageStartupMessages(library(tximport))
  files <- unlist(cfg$quant_sf, use.names = FALSE); names(files) <- samples
  tx2gene <- read.delim(cfg$tx2gene, check.names = FALSE, stringsAsFactors = FALSE)[, 1:2]
  txi <- tximport(files, type = "salmon", tx2gene = tx2gene)
  source_counts <- txi$counts; colnames(source_counts) <- samples
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
vst_table <- read.delim(cfg$l1_vst, check.names = FALSE, stringsAsFactors = FALSE)
rownames(vst_table) <- vst_table[[1]]
vst_matrix <- as.matrix(vst_table[, samples, drop = FALSE])
write_results <- function(df, path) write.table(df, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
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
  plot_df <- data.frame(lfc = result_df$log2FoldChange, plot_padj = ifelse(is.na(result_df$padj), NA_real_, pmax(result_df$padj, .Machine$double.xmin)), group = ifelse(significant & result_df$log2FoldChange > 0, "Up", ifelse(significant & result_df$log2FoldChange < 0, "Down", "Not significant")))
  p <- ggplot(plot_df, aes(x = lfc, y = -log10(plot_padj), color = group)) + geom_point(na.rm = TRUE, alpha = 0.75, size = 1.5) + scale_color_manual(values = c("Up" = "#b2182b", "Down" = "#2166ac", "Not significant" = "#8c8c8c")) + labs(x = paste0("log2 fold change (", numerator, " / ", denominator, ")"), y = "-log10(padj)", color = NULL) + theme_minimal()
  ggsave(file.path(directory, "volcano.png"), p, width = 7, height = 5, dpi = 150)
  ggsave(file.path(directory, "volcano.tiff"), p, width = 7, height = 5, dpi = 300, compression = "lzw")
  heatmap_status <- "NOT_APPLICABLE"
  if (nrow(significant_df) > 0) {
    ordered <- significant_df[order(significant_df$padj, significant_df$gene_id), , drop = FALSE]
    selected <- head(ordered, cfg$heatmap_top_n)
    write.table(data.frame(gene_id = selected$gene_id), file.path(directory, "heatmap_genes.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
    heatmap_samples <- samples[metadata[samples, factor] %in% c(numerator, denominator)]
    heatmap_values <- vst_matrix[selected$gene_id, heatmap_samples, drop = FALSE]
    z <- t(scale(t(heatmap_values))); z[is.na(z)] <- 0
    png(file.path(directory, "heatmap.png"), width = 1050, height = max(450, 18 * nrow(z) + 250), res = 150)
    pheatmap(z, cluster_rows = nrow(z) > 1, cluster_cols = length(heatmap_samples) > 1, main = paste0("Significant genes: ", contrast_id), fontsize_row = ifelse(nrow(z) > 30, 6, 9))
    dev.off(); heatmap_status <- "AVAILABLE"
    tiff(file.path(directory, "heatmap.tiff"), width = 7, height = max(3, 0.12 * nrow(z) + 1.7), units = "in", res = 300, compression = "lzw")
    pheatmap(z, cluster_rows = nrow(z) > 1, cluster_cols = length(heatmap_samples) > 1, main = paste0("Significant genes: ", contrast_id), fontsize_row = ifelse(nrow(z) > 30, 6, 9))
    dev.off()
  }
  tested <- sum(!is.na(result_df$pvalue))
  summary <- list(input_genes = nrow(source_counts), filtered_genes = sum(!keep), retained_genes = sum(keep), tested_genes = tested, pvalue_na = sum(is.na(result_df$pvalue)), padj_na = sum(is.na(result_df$padj)), independent_filtering = TRUE, multiple_testing_method = "Benjamini-Hochberg (DESeq2 default)", fit_method = fit_method, heatmap_status = heatmap_status, heatmap_top_n = cfg$heatmap_top_n)
  write(toJSON(summary, auto_unbox = TRUE, pretty = TRUE), file.path(directory, "backend_summary.json"))
}
