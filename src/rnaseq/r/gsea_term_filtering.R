# Shared post-calculation filtering for GO and KEGG preranked GSEA.
# clusterProfiler performs the enrichment calculation with pvalueCutoff = 1;
# nf-rna alone marks significance (v1.5: p.adjust <= 0.05, enrichment_common.R).
gsea_term_tables <- function(raw, empty_terms) {
  wanted <- setdiff(names(empty_terms), "significant")
  if (is.null(raw) || nrow(raw) == 0) {
    terms <- empty_terms
  } else {
    for (name in setdiff(wanted, names(raw))) raw[[name]] <- NA
    terms <- nf_rna_flag_significant(raw[, wanted, drop=FALSE])
    terms <- terms[order(terms$p.adjust, terms$ID), , drop=FALSE]
  }
  significant <- terms[terms$significant, , drop=FALSE]
  positive <- significant[!is.na(significant$NES) & significant$NES > 0, , drop=FALSE]
  negative <- significant[!is.na(significant$NES) & significant$NES < 0, , drop=FALSE]
  list(terms=terms, significant=significant, positive=positive, negative=negative)
}
