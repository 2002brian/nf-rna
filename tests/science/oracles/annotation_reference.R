# Independent AnnotationDbi reference for ORA/GSEA identifier handling.
#
# Does NOT source src/rnaseq/r.  For the given Ensembl IDs (version suffix
# removed for the query only) it writes every Ensembl -> Entrez pair from
# org.Hs.eg.db, the Entrez IDs annotated (with ancestors) in each GO ontology,
# and the propagated (GOALL) Entrez members of the requested GO terms.
# Usage: Rscript annotation_reference.R IDS.txt TERMS.txt OUT_PREFIX
suppressPackageStartupMessages({ library(AnnotationDbi); library(org.Hs.eg.db) })
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 3)
ids <- readLines(args[[1]])
terms <- readLines(args[[2]])
prefix <- args[[3]]
keys <- intersect(unique(sub("\\.[0-9]+$", "", ids)), keys(org.Hs.eg.db, "ENSEMBL"))
pairs <- suppressMessages(select(org.Hs.eg.db, keys = keys, keytype = "ENSEMBL", columns = "ENTREZID"))
pairs <- unique(pairs[!is.na(pairs$ENTREZID), c("ENSEMBL", "ENTREZID")])
write.table(pairs, paste0(prefix, "_mapping.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
for (ontology in c("BP", "MF", "CC")) {
  annotated <- unique(suppressMessages(select(org.Hs.eg.db, keys = ontology, keytype = "ONTOLOGYALL", columns = "ENTREZID"))$ENTREZID)
  writeLines(sort(annotated), paste0(prefix, "_annotated_", ontology, ".txt"))
}
if (length(terms)) {
  members <- suppressMessages(select(org.Hs.eg.db, keys = terms, keytype = "GOALL", columns = "ENTREZID"))
  members <- unique(members[!is.na(members$ENTREZID), c("GOALL", "ENTREZID")])
  write.table(members, paste0(prefix, "_term_members.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
}
