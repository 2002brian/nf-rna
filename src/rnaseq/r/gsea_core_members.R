# Bounded, deterministic GO GSEA core-member table writer.

write_gsea_core_members <- function(terms, ranked_mapping, path, max_gs_size) {
  header <- data.frame(
    GO_ID=character(), GO_term=character(), entrez_id=character(),
    original_gene_id=character(), symbol=character(), stat=numeric(),
    stringsAsFactors=FALSE
  )
  write.table(header, path, sep="\t", quote=FALSE, row.names=FALSE, na="NA")
  if (nrow(terms) == 0) return(list(rows=0L, pathways=0L))
  required <- c("mapped_entrez_id", "original_gene_id", "mapped_symbol", "stat")
  if (!all(required %in% names(ranked_mapping))) stop("Core-member lookup is missing required ranked-mapping columns.")
  if (anyNA(ranked_mapping$mapped_entrez_id) || anyDuplicated(ranked_mapping$mapped_entrez_id)) stop("Core-member lookup must contain one nonblank retained mapping per Entrez ID.")
  if (any(!is.finite(ranked_mapping$stat)) || any(is.na(ranked_mapping$original_gene_id) | ranked_mapping$original_gene_id == "")) stop("Core-member lookup contains missing source gene IDs or non-finite statistics.")
  lookup <- match
  total_rows <- 0L
  for (i in order(as.character(terms$ID), na.last=TRUE)) {
    core <- terms$core_enrichment[[i]]
    if (is.na(core) || !nzchar(core)) next
    tokens <- strsplit(as.character(core), "/", fixed=TRUE)[[1]]
    ids <- sort(unique(tokens[!is.na(tokens) & nzchar(tokens)]))
    set_size <- suppressWarnings(as.integer(terms$setSize[[i]]))
    if (length(ids) > max_gs_size) stop(sprintf("Core-member invariant failed for %s: %d members exceeds max_gs_size %d.", terms$ID[[i]], length(ids), max_gs_size))
    if (!is.na(set_size) && length(ids) > set_size) stop(sprintf("Core-member invariant failed for %s: %d core members exceeds setSize %d.", terms$ID[[i]], length(ids), set_size))
    positions <- lookup(ids, ranked_mapping$mapped_entrez_id)
    if (anyNA(positions)) {
      missing <- paste(head(ids[is.na(positions)], 3), collapse=", ")
      stop(sprintf("Core-member invariant failed for %s: ranked mapping is missing Entrez ID(s): %s.", terms$ID[[i]], missing))
    }
    related <- ranked_mapping[positions, c("mapped_entrez_id", "original_gene_id", "mapped_symbol", "stat"), drop=FALSE]
    if (any(!is.finite(related$stat)) || any(is.na(related$original_gene_id) | related$original_gene_id == "")) stop(sprintf("Core-member invariant failed for %s: emitted rows require finite stat and original_gene_id.", terms$ID[[i]]))
    block <- data.frame(
      GO_ID=rep(terms$ID[[i]], length(ids)), GO_term=rep(terms$Description[[i]], length(ids)),
      entrez_id=related$mapped_entrez_id, original_gene_id=related$original_gene_id,
      symbol=related$mapped_symbol, stat=related$stat, stringsAsFactors=FALSE
    )
    block <- block[order(block$GO_ID, block$entrez_id, block$original_gene_id, block$symbol, block$stat, na.last=TRUE), , drop=FALSE]
    if (anyDuplicated(paste(block$GO_ID, block$entrez_id, sep="\r"))) stop(sprintf("Core-member invariant failed for %s: duplicate pathway/Entrez rows.", terms$ID[[i]]))
    total_rows <- total_rows + nrow(block)
    if (total_rows > nrow(terms) * max_gs_size) stop("Core-member invariant failed: total rows exceed returned pathways * max_gs_size.")
    write.table(block, path, sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE, append=TRUE, na="NA")
    rm(tokens, related, block)
  }
  list(rows=total_rows, pathways=nrow(terms))
}
