# Shared v1.5 enrichment semantics for GO/KEGG ORA and GO/KEGG preranked GSEA.
# clusterProfiler is always called with pvalueCutoff = 1 and qvalueCutoff = 1
# so it returns every evaluated term; nf-rna alone decides significance.

NF_RNA_SIGNIFICANCE_SEMANTICS <- "v1.5"
NF_RNA_ENRICHMENT_PADJ <- 0.05
NF_RNA_ENRICHMENT_ADJUST <- "BH"
NF_RNA_GSEA_NPERM_SIMPLE <- 10000L
NF_RNA_RUNNING_SCORE_TOP <- 5L
NF_RNA_SIMPLIFY_CUTOFF <- 0.7
NF_RNA_DOTPLOT_TOP <- 15L
NF_RNA_ENRICHMENT_FORMATS <- c("png", "png_raster", "tiff", "pdf")

nf_rna_significance_record <- function() list(
  significance_semantics=NF_RNA_SIGNIFICANCE_SEMANTICS,
  rule=sprintf("significant = p.adjust (%s) <= %s", NF_RNA_ENRICHMENT_ADJUST, NF_RNA_ENRICHMENT_PADJ),
  clusterprofiler_call=list(pAdjustMethod=NF_RNA_ENRICHMENT_ADJUST, pvalueCutoff=1, qvalueCutoff=1),
  returned_terms="every term clusterProfiler evaluated (gene set within minGSSize..maxGSSize and, for ORA, at least one foreground gene)",
  legacy_config_cutoffs="annotation.enrichment.*.pvalue_cutoff/qvalue_cutoff/padj_cutoff are not used for significance under v1.5 semantics"
)

nf_rna_flag_significant <- function(terms) {
  terms$significant <- if (nrow(terms)) !is.na(terms$p.adjust) & terms$p.adjust <= NF_RNA_ENRICHMENT_PADJ else logical()
  terms
}

nf_rna_gsea_settings <- function(seed, min_gs_size, max_gs_size) list(
  ranking_metric="DESeq2 Wald statistic (stat) of the contrast, unmodified; log2FoldChange_shrunk is not used",
  na_handling="genes with a non-finite stat are excluded",
  duplicate_entrez="largest |stat| retained; ties by gene_id",
  seed=as.integer(seed), rng_kind=RNGkind(),
  seed_call="set.seed(seed) immediately before the call, seed = FALSE (DOSE 4.4.0 seed = TRUE re-seeds with set.seed(.Random.seed), a fixed state that ignores the configured seed)",
  leading_edge=NF_RNA_LEADING_EDGE_DEFINITION,
  nPermSimple=NF_RNA_GSEA_NPERM_SIMPLE, eps=0, minGSSize=as.integer(min_gs_size), maxGSSize=as.integer(max_gs_size)
)

# gene_id -> symbol / log2FoldChange / padj from the Phase 1 all_genes.tsv.
nf_rna_gene_lookup <- function(table) {
  column <- function(name, cast) if (name %in% names(table)) suppressWarnings(cast(table[[name]])) else cast(rep(NA, nrow(table)))
  data.frame(gene_id=as.character(table$gene_id), symbol=column("symbol", as.character),
             log2FoldChange=column("log2FoldChange", as.numeric), padj=column("padj", as.numeric), stringsAsFactors=FALSE)
}

# ORA universe and its record.  `source` are the padj-non-NA gene_ids, `mapping`
# their Entrez mapping, `database_genes` every Entrez ID of the gene-set source.
nf_rna_ora_universe <- function(source, mapping, database_genes) {
  mapped <- unique(mapping$original_gene_id[!is.na(mapping$mapped_entrez_id)])
  entrez <- unique(as.character(mapping$mapped_entrez_id[!is.na(mapping$mapped_entrez_id)]))
  universe <- sort(intersect(entrez, database_genes))
  list(genes=universe, record=list(
    definition="gene_id with non-NA DESeq2 padj, mapped to Entrez, intersected with the gene-set database genes",
    padj_non_na_genes=length(source), entrez_mapped_genes=length(mapped),
    entrez_mapping_rate=if (length(source)) length(mapped) / length(source) else 0,
    unique_entrez_ids=length(entrez), database_genes=length(unique(database_genes)), universe_genes=length(universe)))
}

# BgRatio "M/N": N is the universe clusterProfiler actually used.
nf_rna_bgratio_n <- function(terms) {
  if (!nrow(terms) || !"BgRatio" %in% names(terms)) return(NA_integer_)
  unique(as.integer(sub("^[0-9]+/", "", terms$BgRatio)))
}

# Long-format member table: one row per pathway x gene_id.  `pairs` holds
# pathway_id, entrez_id, gene_id, fallback_symbol and (GSEA) leading_edge.
nf_rna_write_members <- function(terms, pairs, lookup, path, gsea=FALSE) {
  columns <- c("pathway_id", "pathway_description", "gene_id", "symbol", "entrez_id", "log2FoldChange", "padj", "pathway_p.adjust", if (gsea) c("NES", "leading_edge"))
  if (!nrow(pairs)) {
    write.table(as.data.frame(setNames(replicate(length(columns), character(), simplify=FALSE), columns)), path, sep="\t", quote=FALSE, row.names=FALSE)
    return(list(rows=0L, pathways=0L, genes=0L, symbol_unmapped=0L, symbol_unmapped_rate=NA_real_))
  }
  term_index <- match(pairs$pathway_id, terms$ID); gene_index <- match(pairs$gene_id, lookup$gene_id)
  symbol <- lookup$symbol[gene_index]; symbol[is.na(symbol) | symbol == ""] <- pairs$fallback_symbol[is.na(symbol) | symbol == ""]
  out <- data.frame(pathway_id=pairs$pathway_id, pathway_description=terms$Description[term_index], gene_id=pairs$gene_id, symbol=symbol,
                    entrez_id=pairs$entrez_id, log2FoldChange=lookup$log2FoldChange[gene_index], padj=lookup$padj[gene_index],
                    pathway_p.adjust=terms$p.adjust[term_index], stringsAsFactors=FALSE, check.names=FALSE)
  if (gsea) { out$NES <- terms$NES[term_index]; out$leading_edge <- pairs$leading_edge }
  out <- out[order(term_index, out$gene_id, out$entrez_id), , drop=FALSE]
  write.table(out, path, sep="\t", quote=FALSE, row.names=FALSE, na="NA")
  unmapped <- sum(is.na(out$symbol) | out$symbol == "")
  list(rows=nrow(out), pathways=length(unique(out$pathway_id)), genes=length(unique(out$gene_id)), symbol_unmapped=unmapped, symbol_unmapped_rate=unmapped / nrow(out))
}

# ORA pairs: foreground genes (geneID, Entrez) of each listed term, expanded to
# every foreground gene_id that maps to that Entrez ID.
nf_rna_ora_member_pairs <- function(terms, mapping) {
  if (!nrow(terms)) return(data.frame(pathway_id=character(), entrez_id=character(), gene_id=character(), fallback_symbol=character()))
  pairs <- do.call(rbind, lapply(seq_len(nrow(terms)), function(i) data.frame(pathway_id=terms$ID[[i]], entrez_id=strsplit(terms$geneID[[i]], "/", fixed=TRUE)[[1]], stringsAsFactors=FALSE)))
  mapped <- mapping[!is.na(mapping$mapped_entrez_id), c("mapped_entrez_id", "original_gene_id", "mapped_symbol"), drop=FALSE]
  out <- merge(pairs, mapped, by.x="entrez_id", by.y="mapped_entrez_id", sort=FALSE)
  unique(data.frame(pathway_id=out$pathway_id, entrez_id=out$entrez_id, gene_id=out$original_gene_id, fallback_symbol=out$mapped_symbol, stringsAsFactors=FALSE))
}

# GSEA pairs: every ranked gene of each listed gene set, with its retained
# gene_id; leading_edge marks membership in core_enrichment.
nf_rna_gsea_member_pairs <- function(result, terms, ranked_mapping) {
  if (!nrow(terms)) return(data.frame(pathway_id=character(), entrez_id=character(), gene_id=character(), fallback_symbol=character(), leading_edge=logical()))
  ranked <- names(result@geneList)
  do.call(rbind, lapply(seq_len(nrow(terms)), function(i) {
    ids <- sort(intersect(as.character(result@geneSets[[terms$ID[[i]]]]), ranked))
    core <- strsplit(as.character(terms$core_enrichment[[i]]), "/", fixed=TRUE)[[1]]
    at <- match(ids, ranked_mapping$mapped_entrez_id)
    if (anyNA(at)) stop(sprintf("ranked mapping is missing Entrez IDs of %s", terms$ID[[i]]))
    data.frame(pathway_id=terms$ID[[i]], entrez_id=ids, gene_id=ranked_mapping$original_gene_id[at], fallback_symbol=ranked_mapping$mapped_symbol[at], leading_edge=ids %in% core, stringsAsFactors=FALSE)
  }))
}

# Up to NF_RNA_RUNNING_SCORE_TOP significant sets per NES sign, by p.adjust
# (ties by ID).  Returns the figure records' reason when nothing is drawn.
nf_rna_running_score_plots <- function(fig, result, significant, contrast, directory, output_root, meta=list()) {
  if (!nrow(significant)) return(list(drawn=0L, reason="no significant gene sets (p.adjust <= 0.05)"))
  dir.create(directory, recursive=TRUE, showWarnings=FALSE)
  ordered <- significant[order(significant$p.adjust, significant$ID), , drop=FALSE]
  picks <- rbind(head(ordered[ordered$NES > 0, , drop=FALSE], NF_RNA_RUNNING_SCORE_TOP), head(ordered[ordered$NES < 0, , drop=FALSE], NF_RNA_RUNNING_SCORE_TOP))
  for (i in seq_len(nrow(picks))) {
    term <- picks[i, , drop=FALSE]; sign <- if (term$NES > 0) "positive" else "negative"
    rank <- sum((picks$NES[seq_len(i)] > 0) == (term$NES > 0))
    stem <- file.path(directory, sprintf("%s_%d_%s", sign, rank, gsub("[^A-Za-z0-9._-]", "_", term$ID)))
    nf_rna_gsea_running_score(fig, result, term, contrast, stem, output_root, c(list(selection=sprintf("top %d significant gene sets per NES sign by p.adjust, ties by ID", NF_RNA_RUNNING_SCORE_TOP)), meta))
  }
  list(drawn=nrow(picks), positive=sum(picks$NES > 0), negative=sum(picks$NES < 0), ids=as.list(picks$ID))
}

# GO redundancy reduction on the significant terms only; the full result is untouched.
nf_rna_simplify_go <- function(result, significant) {
  if (nrow(significant) < 2) return(list(table=significant, reason=if (nrow(significant)) "single significant term; nothing to reduce" else "no significant terms"))
  reduced <- result
  reduced@result <- reduced@result[reduced@result$ID %in% significant$ID, , drop=FALSE]
  reduced <- clusterProfiler::simplify(reduced, cutoff=NF_RNA_SIMPLIFY_CUTOFF, by="p.adjust", select_fun=min, measure="Wang")
  table <- significant[significant$ID %in% reduced@result$ID, , drop=FALSE]
  list(table=table, reason=NULL)
}

nf_rna_simplify_record <- function(before, simplified) list(
  method=sprintf("clusterProfiler::simplify(cutoff = %s, by = \"p.adjust\", select_fun = min, measure = \"Wang\") on significant terms only", NF_RNA_SIMPLIFY_CUTOFF),
  significant_terms=nrow(before), simplified_terms=nrow(simplified$table), reason=simplified$reason)
