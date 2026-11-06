# Independent tximport + DESeq2 reference (does NOT source src/rnaseq/r).
#
# Writes two result tables for the same quant.sf files:
#   *_offset.tsv    DESeqDataSetFromTximport (average-transcript-length offsets), the documented nf-rna path;
#   *_no_offset.tsv DESeqDataSetFromMatrix(round(counts)), which ignores the lengths.
# Usage: Rscript tximport_reference.R FILES.tsv(sample<TAB>path<TAB>condition) TX2GENE.tsv OUT_PREFIX
suppressPackageStartupMessages({ library(DESeq2); library(tximport) })
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 3)
table <- read.delim(args[[1]], header = FALSE, col.names = c("sample", "path", "condition"), colClasses = "character")
files <- setNames(table$path, table$sample)
tx2gene <- read.delim(args[[2]], check.names = FALSE, stringsAsFactors = FALSE)[, 1:2]
txi <- tximport(files, type = "salmon", tx2gene = tx2gene)
metadata <- data.frame(condition = relevel(factor(table$condition), ref = "Control"), row.names = table$sample)
keep <- rowSums(txi$counts) >= 10
write_result <- function(dds, path) {
  dds <- DESeq(dds, quiet = TRUE)
  res <- results(dds, contrast = c("condition", "Treatment", "Control"), alpha = 0.05)
  out <- data.frame(gene_id = rownames(res), as.data.frame(res))[, c("gene_id", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")]
  write.table(out, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
}
write_result(DESeqDataSetFromTximport(txi, metadata, ~ condition)[keep, ], paste0(args[[3]], "_offset.tsv"))
write_result(DESeqDataSetFromMatrix(round(txi$counts[keep, ]), metadata, ~ condition), paste0(args[[3]], "_no_offset.tsv"))
write.table(data.frame(gene_id = rownames(txi$counts), txi$counts, check.names = FALSE), paste0(args[[3]], "_counts.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
