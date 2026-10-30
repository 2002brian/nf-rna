# Shared ORA universe selection.  Both GO and KEGG must use this exact
# definition (v1.5): a retained gene enters the universe only when its DESeq2
# padj is not NA, i.e. it passed independent filtering and outlier handling.
# The Entrez mapping and the intersection with the gene-set database follow in
# nf_rna_ora_universe().
nf_rna_ora_gene_sets <- function(table) {
  if (!all(c("gene_id", "pvalue", "padj") %in% names(table))) {
    stop("L2 all_genes.tsv must contain gene_id, pvalue and padj for ORA")
  }
  gene_id <- as.character(table$gene_id)
  retained_mask <- !is.na(gene_id) & gene_id != ""
  retained <- unique(gene_id[retained_mask])
  pvalue <- suppressWarnings(as.numeric(table$pvalue))
  padj <- suppressWarnings(as.numeric(table$padj))
  tested_mask <- retained_mask & !is.na(padj)
  tested <- unique(gene_id[tested_mask])

  count_unique <- function(mask) length(unique(gene_id[retained_mask & mask]))
  list(
    retained = retained,
    tested = tested,
    counts = list(
      post_l1_retained_genes = length(retained),
      finite_pvalue_genes = count_unique(is.finite(pvalue)),
      padj_non_na_genes = length(tested),
      padj_na_excluded_genes = length(setdiff(retained, tested))
    )
  )
}
