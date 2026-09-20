# Shared ORA universe selection.  Both GO and KEGG must use this exact
# definition: a retained gene is statistically tested only when its L2 p-value
# is finite.  This excludes NA, NaN, Inf, and -Inf deterministically.
nf_rna_ora_gene_sets <- function(table) {
  if (!all(c("gene_id", "pvalue") %in% names(table))) {
    stop("L2 all_genes.tsv must contain gene_id and pvalue for ORA")
  }
  gene_id <- as.character(table$gene_id)
  retained_mask <- !is.na(gene_id) & gene_id != ""
  retained <- unique(gene_id[retained_mask])
  pvalue <- suppressWarnings(as.numeric(table$pvalue))
  tested_mask <- retained_mask & is.finite(pvalue)
  tested <- unique(gene_id[tested_mask])

  count_unique <- function(mask) length(unique(gene_id[retained_mask & mask]))
  list(
    retained = retained,
    tested = tested,
    counts = list(
      post_l1_retained_genes = length(retained),
      statistically_tested_genes = length(tested),
      pvalue_na_excluded_genes = count_unique(is.na(pvalue) & !is.nan(pvalue)),
      pvalue_nan_excluded_genes = count_unique(is.nan(pvalue)),
      pvalue_infinite_excluded_genes = count_unique(is.infinite(pvalue)),
      pvalue_nonfinite_excluded_genes = length(setdiff(retained, tested))
    )
  )
}
