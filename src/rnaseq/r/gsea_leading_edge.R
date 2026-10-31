# Leading edge of preranked GSEA results by the standard running-sum definition
# (Subramanian et al. 2005, PNAS 102:15545), replacing DOSE's extraction.
#
# DOSE 4.4.0 leading_edge() keeps only the running score after each hit and, for
# ES < 0, takes which.min() over those post-step values.  The running-sum minimum
# lies on a miss immediately before a hit, so DOSE can select an earlier hit and
# report a different core_enrichment, tags and signal.  ES, NES and p-values come
# from fgsea and are not touched here.
#
# Definition (weight p = 1, the gseaParam/exponent nf-rna uses):
#   hits h_1 < ... < h_k are the positions of the set's genes in the ranked list;
#   the running sum after hit j is top_j = sum_{i<=j} |r_i| / N_R - (h_j - j) / (N - k),
#   just before hit j it is bottom_j = top_j - |r_j| / N_R.
#   ES > 0: the maximum is top_i; the leading edge is hits 1..i.
#   ES < 0: the minimum is bottom_i; the leading edge is hits i..k.
# A tied extremum resolves to its first occurrence along the list (as fgsea and
# GSEA desktop).  Gene positions are those of the gseaResult geneList, the exact
# ordered vector passed to fgsea (DESeq2 stat descending, ties by Entrez ID).

NF_RNA_LEADING_EDGE_DEFINITION <- paste(
  "standard running-sum leading edge (Subramanian et al. 2005), recomputed by nf-rna from the ranked list and gene set:",
  "ES > 0, hits up to the running-sum maximum; ES < 0, hits from the running-sum minimum to the end;",
  "first extremum on ties; the recomputed extremum must equal the fgsea ES")
NF_RNA_LEADING_EDGE_ES_TOLERANCE <- 1e-9

nf_rna_leading_edge <- function(gene_list, members, es) {
  n <- length(gene_list)
  hits <- sort(unique(match(as.character(members), names(gene_list))))
  hits <- hits[!is.na(hits)]
  k <- length(hits)
  if (k == 0L || k >= n) stop("leading edge: gene set must contain at least one and fewer than all ranked genes")
  weight <- abs(unname(gene_list[hits])); total <- sum(weight)
  if (!is.finite(total) || total <= 0) stop("leading edge: gene set has no non-zero ranking statistic")
  top <- cumsum(weight) / total - (hits - seq_len(k)) / (n - k)
  bottom <- top - weight / total
  if (es >= 0) {
    i <- which.max(top); core <- hits[seq_len(i)]; extremum <- top[[i]]; rank <- hits[[i]]
  } else {
    # The minimum sits at position h_i - 1; DOSE counts its rank from there to the end.
    i <- which.min(bottom); core <- hits[i:k]; extremum <- bottom[[i]]; rank <- n - hits[[i]] + 2L
  }
  list(core=names(gene_list)[core], extremum=extremum, rank=as.integer(rank), set_size=k, ranked_genes=n)
}

# Replaces rank, leading_edge and core_enrichment of a gseaResult in place and
# returns it with an audit of how many terms differ from DOSE's extraction.
nf_rna_correct_leading_edge <- function(result) {
  audit <- list(definition=NF_RNA_LEADING_EDGE_DEFINITION, terms=0L, core_enrichment_changed=0L,
                core_enrichment_changed_positive_es=0L, core_enrichment_changed_negative_es=0L, max_abs_es_difference=0)
  if (is.null(result) || !nrow(result@result)) return(list(result=result, audit=audit))
  exponent <- result@params$exponent
  if (!is.null(exponent) && !isTRUE(all.equal(as.numeric(exponent), 1))) stop("leading edge: only the GSEA weight exponent 1 is supported")
  terms <- result@result
  fixed <- lapply(seq_len(nrow(terms)), function(i) nf_rna_leading_edge(result@geneList, result@geneSets[[terms$ID[[i]]]], terms$enrichmentScore[[i]]))
  es_difference <- abs(vapply(fixed, `[[`, numeric(1), "extremum") - terms$enrichmentScore)
  if (any(!(es_difference <= NF_RNA_LEADING_EDGE_ES_TOLERANCE))) stop(sprintf("leading edge: recomputed running-sum extremum differs from the fgsea ES for %s", terms$ID[[which(!(es_difference <= NF_RNA_LEADING_EDGE_ES_TOLERANCE))[[1]]]]))
  set_size <- vapply(fixed, `[[`, integer(1), "set_size")
  if (!identical(as.integer(terms$setSize), set_size)) stop("leading edge: gene-set membership in the ranked list disagrees with the fgsea setSize")
  n <- length(result@geneList)
  core <- vapply(fixed, function(x) paste(x$core, collapse="/"), character(1))
  rank <- vapply(fixed, `[[`, integer(1), "rank")
  tags <- lengths(lapply(fixed, `[[`, "core")) / set_size
  list_fraction <- rank / n
  signal <- tags * (1 - list_fraction) * (n / (n - set_size))
  percent <- function(x) paste0(round(x * 100), "%")
  changed <- core != as.character(terms$core_enrichment)
  terms$rank <- rank
  terms$leading_edge <- paste0("tags=", percent(tags), ", list=", percent(list_fraction), ", signal=", percent(signal))
  terms$core_enrichment <- core
  result@result <- terms
  audit$terms <- nrow(terms); audit$core_enrichment_changed <- sum(changed)
  audit$core_enrichment_changed_positive_es <- sum(changed & terms$enrichmentScore >= 0)
  audit$core_enrichment_changed_negative_es <- sum(changed & terms$enrichmentScore < 0)
  audit$max_abs_es_difference <- max(es_difference)
  list(result=result, audit=audit)
}
