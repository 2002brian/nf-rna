# Bounded, deterministic KEGG member-table writer for GSEA and legacy callers.

write_kegg_member_table <- function(terms, mapping, member_column, path, strict=FALSE, max_members=Inf) {
  header <- data.frame(
    pathway_id=character(), pathway_description=character(), entrez_id=character(),
    symbol=character(), original_gene_id=character(), stat=numeric(),
    stringsAsFactors=FALSE
  )
  write.table(header, path, sep="\t", quote=FALSE, row.names=FALSE, na="NA")
  if (nrow(terms) == 0) return(list(rows=0L, pathways=0L, pathway_entrez_pairs=0L))
  required <- c("mapped_entrez_id", "mapped_symbol", "original_gene_id", "stat")
  if (!all(required %in% names(mapping))) stop("KEGG member lookup is missing required mapping columns.")
  if (!(member_column %in% names(terms))) stop("KEGG member table is missing its member column.")

  mapped_ids <- as.character(mapping$mapped_entrez_id)
  valid <- !is.na(mapped_ids) & nzchar(mapped_ids)
  lookup <- split(which(valid), mapped_ids[valid], drop=TRUE)
  total_rows <- 0L
  total_pairs <- 0L
  for (i in order(as.character(terms$ID), na.last=TRUE)) {
    raw_members <- terms[[member_column]][[i]]
    tokens <- if (is.na(raw_members) || !nzchar(as.character(raw_members))) character() else strsplit(as.character(raw_members), "/", fixed=TRUE)[[1]]
    ids <- sort(unique(tokens[!is.na(tokens) & nzchar(tokens)]))
    if (strict) {
      set_size <- suppressWarnings(as.integer(terms$setSize[[i]]))
      if (length(ids) > max_members) stop(sprintf("KEGG core-member invariant failed for %s: %d members exceeds max_gs_size %d.", terms$ID[[i]], length(ids), max_members))
      if (!is.na(set_size) && length(ids) > set_size) stop(sprintf("KEGG core-member invariant failed for %s: %d core members exceeds setSize %d.", terms$ID[[i]], length(ids), set_size))
    }
    emitted_ids <- character()
    for (id in ids) {
      positions <- lookup[[id]]
      if (is.null(positions)) {
        if (strict) stop(sprintf("KEGG core-member invariant failed for %s: mapped universe is missing Entrez ID %s.", terms$ID[[i]], id))
        related <- data.frame(mapped_entrez_id=id, mapped_symbol=NA_character_, original_gene_id=NA_character_, stat=NA_real_, stringsAsFactors=FALSE)
      } else {
        related <- mapping[positions, c("mapped_entrez_id", "mapped_symbol", "original_gene_id", "stat"), drop=FALSE]
        if (strict && (any(!is.finite(related$stat)) || any(is.na(related$original_gene_id) | related$original_gene_id == ""))) stop(sprintf("KEGG core-member invariant failed for %s: emitted rows require finite stat and original_gene_id.", terms$ID[[i]]))
        related <- related[order(as.character(related$original_gene_id), as.character(related$mapped_symbol), related$stat, na.last=TRUE), , drop=FALSE]
      }
      block <- data.frame(
        pathway_id=rep(as.character(terms$ID[[i]]), nrow(related)),
        pathway_description=rep(as.character(terms$Description[[i]]), nrow(related)),
        entrez_id=as.character(related$mapped_entrez_id), symbol=as.character(related$mapped_symbol),
        original_gene_id=as.character(related$original_gene_id), stat=as.numeric(related$stat),
        stringsAsFactors=FALSE
      )
      if (strict && (nrow(block) > length(positions) || any(block$entrez_id != id))) stop(sprintf("KEGG core-member invariant failed for %s: member lookup produced unrelated rows.", terms$ID[[i]]))
      write.table(block, path, sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE, append=TRUE, na="NA")
      total_rows <- total_rows + nrow(block)
      total_pairs <- total_pairs + 1L
      emitted_ids <- c(emitted_ids, id)
      rm(related, block, positions)
    }
    if (strict && length(unique(emitted_ids)) != length(ids)) stop(sprintf("KEGG core-member invariant failed for %s: emitted Entrez IDs disagree with core_enrichment.", terms$ID[[i]]))
    if (is.finite(max_members) && total_pairs > nrow(terms) * max_members) stop("KEGG core-member invariant failed: pathway/Entrez pairs exceed the configured cardinality bound.")
    rm(raw_members, tokens, ids, emitted_ids)
    gc(verbose=FALSE)
  }
  list(rows=total_rows, pathways=nrow(terms), pathway_entrez_pairs=total_pairs)
}
