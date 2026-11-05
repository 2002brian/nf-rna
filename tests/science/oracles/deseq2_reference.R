# Independent DESeq2 reference for the scientific regression gate.
#
# Deliberately does NOT source anything under src/rnaseq/r: it re-implements
# nf-rna's documented L2 contract from the scientific contract alone
# (docs/scientific_contract.md) and calls DESeq2 directly:
#   * remove genes with total count < 10 (all-zero genes included);
#   * one DESeq() fit of the declared formula, categorical variables as factors;
#   * results(contrast = c(factor, numerator, denominator), alpha = padj threshold);
#   * apeglm shrinkage of the condition coefficient for display.
# Usage: Rscript deseq2_reference.R COUNTS.csv METADATA.csv FORMULA CATEGORICAL(comma) ALPHA OUT_PREFIX
suppressPackageStartupMessages(library(DESeq2))
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 6)
counts <- read.csv(args[[1]], check.names = FALSE, row.names = 1)
metadata <- read.csv(args[[2]], check.names = FALSE, colClasses = "character", row.names = 1)
formula <- as.formula(args[[3]])
categorical <- strsplit(args[[4]], ",", fixed = TRUE)[[1]]
alpha <- as.numeric(args[[5]])
prefix <- args[[6]]
metadata <- metadata[colnames(counts), , drop = FALSE]
for (name in categorical) metadata[[name]] <- factor(metadata[[name]], levels = sort(unique(metadata[[name]]), method = "radix"))
metadata$condition <- relevel(metadata$condition, ref = "Control")
counts <- as.matrix(counts)
keep <- rowSums(counts) >= 10
dds <- DESeqDataSetFromMatrix(counts[keep, , drop = FALSE], metadata, formula)
dds <- DESeq(dds, quiet = TRUE)
res <- results(dds, contrast = c("condition", "Treatment", "Control"), alpha = alpha)
shrunk <- lfcShrink(dds, coef = "condition_Treatment_vs_Control", type = "apeglm", quiet = TRUE)
out <- data.frame(gene_id = rownames(res), baseMean = res$baseMean, log2FoldChange = res$log2FoldChange, lfcSE = res$lfcSE,
                  stat = res$stat, pvalue = res$pvalue, padj = res$padj,
                  log2FoldChange_shrunk = shrunk[rownames(res), "log2FoldChange"], check.names = FALSE)
write.table(out, paste0(prefix, "_results.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
info <- list(retained = sum(keep), filter_threshold = unname(metadata(res)$filterThreshold),
             cooks_applied = any(!is.na(mcols(dds)$maxCooks)), size_factors = as.list(sizeFactors(dds)))
writeLines(jsonlite::toJSON(info, auto_unbox = TRUE, digits = NA), paste0(prefix, "_info.json"))
