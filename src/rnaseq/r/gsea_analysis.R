args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2 || args[[1]] != "--config") stop("usage: gsea_analysis.R --config CONFIG.json")
script_arg <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", script_arg[grep("^--file=", script_arg)][[1]])
source(file.path(dirname(normalizePath(script_file)), "gsea_core_members.R"))
suppressPackageStartupMessages({ library(jsonlite); library(AnnotationDbi); library(clusterProfiler); library(ggplot2) })

cfg <- fromJSON(args[[2]], simplifyVector = FALSE)
dir.create(cfg$output_dir, recursive = TRUE, showWarnings = FALSE)
orgdb <- get(cfg$orgdb_package, envir = asNamespace(cfg$orgdb_package))
input_type <- cfg$annotation$input_id_type
gsea_cfg <- cfg$annotation$enrichment$gsea
normalize_id <- function(x) if (input_type == "ENSEMBL") sub("\\.[0-9]+$", "", x) else x

empty_terms <- function() data.frame(
  ID=character(), Description=character(), setSize=integer(), enrichmentScore=numeric(), NES=numeric(),
  pvalue=numeric(), p.adjust=numeric(), qvalue=numeric(), rank=integer(), leading_edge=character(), core_enrichment=character(),
  stringsAsFactors=FALSE
)
write_table <- function(x, path) write.table(x, path, sep="\t", quote=FALSE, row.names=FALSE, na="NA")

map_ranked_sources <- function(source, stat) {
  normalized <- normalize_id(source)
  valid_keys <- AnnotationDbi::keys(orgdb, keytype=input_type)
  requested <- intersect(unique(normalized), valid_keys)
  selected <- if (length(requested) == 0) data.frame() else suppressMessages(AnnotationDbi::select(orgdb, keys=requested, keytype=input_type, columns=c("ENTREZID", "SYMBOL")))
  if (nrow(selected) == 0) selected <- data.frame(key=character(), ENTREZID=character(), SYMBOL=character(), stringsAsFactors=FALSE)
  names(selected)[[1]] <- input_type
  selected <- selected[!is.na(selected$ENTREZID) & selected$ENTREZID != "", , drop=FALSE]
  rows <- lapply(seq_along(source), function(i) {
    hit <- selected[selected[[input_type]] == normalized[[i]], c("ENTREZID", "SYMBOL"), drop=FALSE]
    targets <- unique(hit$ENTREZID)
    status <- if (length(targets) == 0) "UNMAPPED" else if (length(targets) == 1) "MAPPED_UNIQUE" else "ONE_TO_MANY"
    if (nrow(hit) == 0) data.frame(original_gene_id=source[[i]], normalized_gene_id=normalized[[i]], mapped_entrez_id=NA_character_, mapped_symbol=NA_character_, stat=stat[[i]], mapping_status=status, stringsAsFactors=FALSE) else data.frame(original_gene_id=source[[i]], normalized_gene_id=normalized[[i]], mapped_entrez_id=as.character(hit$ENTREZID), mapped_symbol=as.character(hit$SYMBOL), stat=stat[[i]], mapping_status=status, stringsAsFactors=FALSE)
  })
  mapping <- do.call(rbind, rows)
  mapping$target_duplicate <- FALSE
  mapping$rank_status <- ifelse(is.na(mapping$mapped_entrez_id), "UNMAPPED", "CANDIDATE")
  mapped <- which(!is.na(mapping$mapped_entrez_id))
  if (length(mapped) > 0) {
    duplicate_targets <- names(table(mapping$mapped_entrez_id[mapped]))[table(mapping$mapped_entrez_id[mapped]) > 1]
    mapping$target_duplicate[mapped] <- mapping$mapped_entrez_id[mapped] %in% duplicate_targets
    for (target in unique(mapping$mapped_entrez_id[mapped])) {
      idx <- mapped[mapping$mapped_entrez_id[mapped] == target]
      candidates <- mapping[idx, , drop=FALSE]
      order_idx <- order(-abs(candidates$stat), candidates$original_gene_id, candidates$normalized_gene_id)
      winner <- idx[[order_idx[[1]]]]
      mapping$rank_status[idx] <- "COLLAPSED_DUPLICATE_TARGET"
      mapping$rank_status[winner] <- "RETAINED"
    }
  }
  mapping
}

run_ontology <- function(gene_list, ranked_mapping, ontology, root) {
  directory <- file.path(root, ontology); dir.create(directory, recursive=TRUE, showWarnings=FALSE)
  result <- raw <- terms <- significant <- positive <- negative <- top <- p <- core_audit <- NULL
  on.exit({ rm(result, raw, terms, significant, positive, negative, top, p, core_audit); gc(verbose=FALSE) }, add=TRUE)
  tryCatch({
    set.seed(as.integer(gsea_cfg$seed))
    result <- suppressMessages(gseGO(
      geneList=gene_list, OrgDb=orgdb, keyType="ENTREZID", ont=ontology,
      minGSSize=as.integer(gsea_cfg$min_gs_size), maxGSSize=as.integer(gsea_cfg$max_gs_size),
      pvalueCutoff=as.numeric(gsea_cfg$pvalue_cutoff), pAdjustMethod=gsea_cfg$p_adjust_method,
      eps=0, verbose=FALSE, seed=TRUE
    ))
    raw <- if (is.null(result)) data.frame() else as.data.frame(result)
    wanted <- names(empty_terms())
    if (nrow(raw) == 0) {
      terms <- empty_terms()
    } else {
      for (name in setdiff(wanted, names(raw))) raw[[name]] <- NA
      terms <- raw[, wanted, drop=FALSE]
      terms <- terms[order(terms$p.adjust, terms$ID), , drop=FALSE]
    }
    significant <- terms[!is.na(terms$p.adjust) & terms$p.adjust <= as.numeric(gsea_cfg$padj_cutoff) & !is.na(terms$pvalue) & terms$pvalue <= as.numeric(gsea_cfg$pvalue_cutoff), , drop=FALSE]
    positive <- significant[!is.na(significant$NES) & significant$NES > 0, , drop=FALSE]
    negative <- significant[!is.na(significant$NES) & significant$NES < 0, , drop=FALSE]
    write_table(terms, file.path(directory, "all_terms.tsv"))
    write_table(significant, file.path(directory, "significant.tsv"))
    write_table(positive, file.path(directory, "positive_enrichment.tsv"))
    write_table(negative, file.path(directory, "negative_enrichment.tsv"))
    core_audit <- write_gsea_core_members(terms, ranked_mapping, file.path(directory, "gsea_core_members.tsv"), as.integer(gsea_cfg$max_gs_size))
    plot_terms <- significant[!is.na(significant$NES), , drop=FALSE]
    if (nrow(plot_terms) > 0) {
      top <- head(plot_terms[order(-abs(plot_terms$NES), plot_terms$ID), , drop=FALSE], 15)
      top$Description <- factor(top$Description, levels=rev(top$Description))
      p <- ggplot(top, aes(x=Description, y=NES, size=setSize, color=NES)) + geom_point() + coord_flip() + scale_color_gradient2(low="#b2182b", mid="grey90", high="#2166ac") + labs(x=NULL, y="Normalized enrichment score") + theme_minimal()
    ggsave(file.path(directory, "dotplot.png"), p, width=7, height=5, dpi=150)
    ggsave(file.path(directory, "dotplot.tiff"), p, width=7, height=5, dpi=300, compression="lzw")
      rm(plot_terms)
    }
    summary <- list(status=if (nrow(significant) == 0) "NO_SIGNIFICANT_TERMS" else "SUCCESS", all_terms=nrow(terms), significant_terms=nrow(significant), positive_terms=nrow(positive), negative_terms=nrow(negative), na_pathways=sum(is.na(terms$pvalue) | is.na(terms$p.adjust) | is.na(terms$NES)), core_member_rows=core_audit$rows)
    rm(result, raw, terms, significant, positive, negative, top, p, core_audit)
    gc(verbose=FALSE)
    summary
  }, error=function(error) list(status="FAILED", reason=conditionMessage(error), all_terms=0, significant_terms=0, positive_terms=0, negative_terms=0, na_pathways=0))
}

contrast_summaries <- list()
blocked_reasons <- character()
for (contrast in cfg$contrasts) {
  root <- file.path(cfg$output_dir, contrast$contrast_id); dir.create(root, recursive=TRUE, showWarnings=FALSE)
  source_table <- read.delim(contrast$all_genes, check.names=FALSE, stringsAsFactors=FALSE)
  if (!all(c("gene_id", "stat") %in% names(source_table))) stop("M4A all_genes.tsv must contain gene_id and stat columns")
  raw_stat <- suppressWarnings(as.numeric(source_table$stat))
  finite <- is.finite(raw_stat) & !is.na(source_table$gene_id) & source_table$gene_id != ""
  mapping <- map_ranked_sources(as.character(source_table$gene_id[finite]), raw_stat[finite])
  write_table(mapping, file.path(root, "gsea_rank_mapping.tsv"))
  retained <- mapping[mapping$rank_status == "RETAINED", , drop=FALSE]
  retained <- retained[order(-retained$stat, retained$mapped_entrez_id), , drop=FALSE]
  ranked <- data.frame(rank=seq_len(nrow(retained)), entrez_id=retained$mapped_entrez_id, stat=retained$stat, original_gene_id=retained$original_gene_id, symbol=retained$mapped_symbol, stringsAsFactors=FALSE)
  write_table(ranked, file.path(root, "ranked_gene_list.tsv"))
  mapped_sources <- unique(mapping$original_gene_id[!is.na(mapping$mapped_entrez_id)])
  source_count <- sum(finite)
  mapping_rate <- if (source_count == 0) 0 else length(mapped_sources) / source_count
  ranking <- list(
    all_genes_rows=nrow(source_table), finite_stat_source_genes=source_count,
    nonfinite_or_missing_stat_source_genes=nrow(source_table)-source_count,
    mapped_source_genes=length(mapped_sources), unmapped_source_genes=source_count-length(mapped_sources),
    mapping_rate=mapping_rate, one_to_many_source_ids=length(unique(mapping$original_gene_id[mapping$mapping_status == "ONE_TO_MANY"])),
    duplicate_target_ids=sum(table(mapping$mapped_entrez_id[!is.na(mapping$mapped_entrez_id)]) > 1),
    duplicate_target_rows_collapsed=sum(mapping$rank_status == "COLLAPSED_DUPLICATE_TARGET"),
    final_ranked_genes=nrow(ranked), positive_stats=sum(ranked$stat > 0), negative_stats=sum(ranked$stat < 0), zero_stats=sum(ranked$stat == 0), stat_ties=sum(duplicated(ranked$stat)), tie_handling="DESeq2 stat is unmodified; ties are ordered by ascending Entrez ID before fgsea"
  )
  rank_reason <- NULL
  if (mapping_rate < cfg$annotation$minimum_mapping_rate) rank_reason <- paste0("GO preranked GSEA blocked: mapped ", sprintf("%.1f%%", 100*mapping_rate), " of finite tested genes; required minimum is ", sprintf("%.1f%%", 100*cfg$annotation$minimum_mapping_rate), ".")
  if (is.null(rank_reason) && nrow(ranked) < as.integer(gsea_cfg$minimum_ranked_genes)) rank_reason <- paste0("GO preranked GSEA blocked: ", nrow(ranked), " unique ranked Entrez genes; minimum is ", as.integer(gsea_cfg$minimum_ranked_genes), ".")
  if (!is.null(rank_reason)) {
    write(toJSON(c(ranking, list(status="BLOCKED", reason=rank_reason)), auto_unbox=TRUE, pretty=TRUE), file.path(root, "gsea_ranking_summary.json"))
    contrast_summaries[[length(contrast_summaries)+1]] <- list(contrast_id=contrast$contrast_id, ranking=ranking, status="BLOCKED", reason=rank_reason, ontologies=list())
    blocked_reasons <- c(blocked_reasons, rank_reason)
    next
  }
  ranked_mapping <- mapping[mapping$rank_status == "RETAINED" & !is.na(mapping$mapped_entrez_id), , drop=FALSE]
  if (anyDuplicated(ranked_mapping$mapped_entrez_id) || any(!is.finite(ranked_mapping$stat))) stop("Deterministic retained rank mapping is invalid for GSEA core-member provenance.")
  gene_list <- ranked$stat; names(gene_list) <- ranked$entrez_id
  ontologies <- list(); completed_ontologies <- character(); failed_ontology <- NULL
  for (ontology in c("BP", "MF", "CC")) {
    outcome <- run_ontology(gene_list, ranked_mapping, ontology, root)
    ontologies[[ontology]] <- outcome
    if (outcome$status == "FAILED") { failed_ontology <- ontology; break }
    completed_ontologies <- c(completed_ontologies, ontology)
    rm(outcome); gc(verbose=FALSE)
  }
  contrast_status <- if (!is.null(failed_ontology)) "FAILED" else if (all(vapply(ontologies, function(x) x$status == "NO_SIGNIFICANT_TERMS", logical(1)))) "NO_SIGNIFICANT_TERMS" else "SUCCESS"
  write(toJSON(c(ranking, list(status=contrast_status, completed_ontologies=as.list(completed_ontologies), failed_ontology=failed_ontology)), auto_unbox=TRUE, pretty=TRUE, null="null"), file.path(root, "gsea_ranking_summary.json"))
  contrast_summaries[[length(contrast_summaries)+1]] <- list(contrast_id=contrast$contrast_id, ranking=ranking, status=contrast_status, completed_ontologies=as.list(completed_ontologies), failed_ontology=failed_ontology, ontologies=ontologies)
  rm(ranked_mapping)
}
overall <- if (any(vapply(contrast_summaries, function(x) x$status == "FAILED", logical(1)))) "FAILED" else if (length(blocked_reasons) > 0) "BLOCKED" else if (all(vapply(contrast_summaries, function(x) x$status == "NO_SIGNIFICANT_TERMS", logical(1)))) "NO_SIGNIFICANT_TERMS" else "SUCCESS"
summary <- list(status=overall, reason=if (length(blocked_reasons)>0) paste(blocked_reasons, collapse=" ") else NULL, annotation_database=cfg$orgdb_package, annotation_database_version=as.character(packageVersion(cfg$orgdb_package)), clusterProfiler_version=as.character(packageVersion("clusterProfiler")), contrasts=contrast_summaries)
write(toJSON(summary, auto_unbox=TRUE, pretty=TRUE, null="null"), file.path(cfg$output_dir, "gsea_backend_summary.json"))
