# Shared post-calculation filtering for GO and KEGG preranked GSEA.
# clusterProfiler performs the enrichment calculation; nf-rna alone applies
# the configured reporting/significance predicate to its returned terms.
gsea_term_tables <- function(raw, empty_terms, pvalue_cutoff, padj_cutoff) {
  wanted <- names(empty_terms)
  if (is.null(raw) || nrow(raw) == 0) {
    terms <- empty_terms
  } else {
    for (name in setdiff(wanted, names(raw))) raw[[name]] <- NA
    terms <- raw[, wanted, drop=FALSE]
    terms <- terms[order(terms$p.adjust, terms$ID), , drop=FALSE]
  }
  # This remains nf-rna's established public reporting predicate. The
  # calculation call itself always uses pvalueCutoff=1 and does not prefilter.
  significant <- terms[
    !is.na(terms$p.adjust) & terms$p.adjust <= as.numeric(padj_cutoff) &
      !is.na(terms$pvalue) & terms$pvalue <= as.numeric(pvalue_cutoff),
    , drop=FALSE
  ]
  positive <- significant[!is.na(significant$NES) & significant$NES > 0, , drop=FALSE]
  negative <- significant[!is.na(significant$NES) & significant$NES < 0, , drop=FALSE]
  list(terms=terms, significant=significant, positive=positive, negative=negative)
}
