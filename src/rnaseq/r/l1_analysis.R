args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2 || args[[1]] != "--config") stop("usage: l1_analysis.R --config CONFIG.json")
suppressPackageStartupMessages({ library(jsonlite); library(yaml); library(DESeq2); library(ggplot2); library(pheatmap) })
cfg <- fromJSON(args[[2]], simplifyVector = FALSE)
dir.create(cfg$output_dir, recursive = TRUE, showWarnings = FALSE)
metadata <- read.csv(cfg$metadata, check.names = FALSE, stringsAsFactors = FALSE)
rownames(metadata) <- metadata$sample_id
samples <- unlist(cfg$samples, use.names = FALSE)
metadata <- metadata[samples, , drop = FALSE]
if (!is.null(cfg$pairing_column)) {
  if (!(cfg$pairing_column %in% colnames(metadata))) stop("configured pairing_column is absent from metadata")
  metadata[[cfg$pairing_column]] <- factor(metadata[[cfg$pairing_column]])
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
  files <- unlist(cfg$quant_sf, use.names = FALSE)
  names(files) <- samples
  tx2gene <- read.delim(cfg$tx2gene, check.names = FALSE, stringsAsFactors = FALSE)[, 1:2]
  txi <- tximport(files, type = "salmon", tx2gene = tx2gene)
  source_counts <- txi$counts
  colnames(source_counts) <- samples
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
size_factor_values <- sizeFactors(dds)
# DESeqDataSetFromTximport may use per-gene normalization factors (offsets)
# rather than scalar size factors; retain that distinction in the QC table.
if (is.null(size_factor_values)) size_factor_values <- rep(NA_real_, length(samples))
library_size <- data.frame(sample_id = samples, input_total = as.numeric(colSums(source_counts)), retained_total = as.numeric(colSums(source_counts[keep, , drop = FALSE])), size_factor = as.numeric(size_factor_values), normalized_total = as.numeric(colSums(normalized)), check.names = FALSE)
write.table(library_size, file.path(cfg$output_dir, "library_size_qc.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
pca <- prcomp(t(vst_matrix), center = TRUE, scale. = FALSE)
variance <- pca$sdev^2 / sum(pca$sdev^2)
scores <- data.frame(sample_id = samples, PC1 = pca$x[, 1], PC2 = if (ncol(pca$x) >= 2) pca$x[, 2] else 0, check.names = FALSE)
write.table(scores, file.path(cfg$output_dir, "pca_scores.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(data.frame(component = paste0("PC", seq_along(variance)), proportion_variance = variance), file.path(cfg$output_dir, "pca_variance.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
color_name <- tail(all.vars(formula), 1)
scores$group <- as.factor(metadata[scores$sample_id, color_name])
p <- ggplot(scores, aes(x = PC1, y = PC2, color = group, label = sample_id)) + geom_point(size = 3) + geom_text(vjust = -0.8, show.legend = FALSE) + labs(color = color_name, x = sprintf("PC1 (%.1f%%)", 100 * variance[[1]]), y = sprintf("PC2 (%.1f%%)", 100 * ifelse(length(variance) >= 2, variance[[2]], 0))) + theme_minimal()
ggsave(file.path(cfg$output_dir, "pca.png"), p, width = 7, height = 5, dpi = 150)
ggsave(file.path(cfg$output_dir, "pca.tiff"), p, width = 7, height = 5, dpi = 300, compression = "lzw")
corr <- cor(vst_matrix, method = "pearson")
write.table(data.frame(sample_id = rownames(corr), corr, check.names = FALSE), file.path(cfg$output_dir, "sample_correlation.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
png(file.path(cfg$output_dir, "sample_correlation.png"), width = 1050, height = 900, res = 150)
pheatmap(corr, main = "Sample correlation (blind VST)")
dev.off()
tiff(file.path(cfg$output_dir, "sample_correlation.tiff"), width = 7, height = 6, units = "in", res = 300, compression = "lzw")
pheatmap(corr, main = "Sample correlation (blind VST)")
dev.off()
flags <- character()
median_library <- median(library_size$input_total)
for (i in seq_len(nrow(library_size))) if (library_size$input_total[[i]] < median_library * 0.25) flags <- c(flags, paste0("LOW_LIBRARY_SIZE: ", library_size$sample_id[[i]], " is below 25% of the median input total."))
summary <- list(source_type = cfg$source_type, samples = samples, genes_input = nrow(source_counts), genes_removed_all_zero = sum(all_zero), genes_removed_low_total = sum(low_total), genes_retained = sum(keep), filter = cfg$filter, normalization = list(method = if (cfg$source_type == "salmon_tximport") "DESeq2 with tximport average-transcript-length normalization offsets" else "DESeq2 median-ratio size factors", scalar_size_factors_available = !all(is.na(size_factor_values))), vst = list(method = if (nrow(dds) < 1000) "DESeq2 varianceStabilizingTransformation (exact small-matrix path)" else "DESeq2 vst"), package_versions = list(R = R.version.string, DESeq2 = as.character(packageVersion("DESeq2")), tximport = if (cfg$source_type == "salmon_tximport") as.character(packageVersion("tximport")) else NULL), pca_proportion_variance = as.list(variance), flags = as.list(flags))
write(toJSON(summary, auto_unbox = TRUE, pretty = TRUE, null = "null"), file.path(cfg$output_dir, "backend_summary.json"))
