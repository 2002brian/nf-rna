args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2 || args[[1]] != "--config") stop("usage: gsea_analysis.R --config CONFIG.json")
script_arg <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", script_arg[grep("^--file=", script_arg)][[1]])
source(file.path(dirname(normalizePath(script_file)), "provenance.R"))
source(file.path(dirname(normalizePath(script_file)), "gene_annotation.R"))
source(file.path(dirname(normalizePath(script_file)), "gsea_core_members.R"))
source(file.path(dirname(normalizePath(script_file)), "annotation_mapping_qc.R"))
source(file.path(dirname(normalizePath(script_file)), "gsea_term_filtering.R"))
source(file.path(dirname(normalizePath(script_file)), "figures.R"))
source(file.path(dirname(normalizePath(script_file)), "enrichment_common.R"))
suppressPackageStartupMessages({ library(jsonlite); library(AnnotationDbi); library(clusterProfiler); library(ggplot2) })

cfg <- fromJSON(args[[2]], simplifyVector = FALSE)
dir.create(cfg$output_dir, recursive = TRUE, showWarnings = FALSE)
fig <- nf_rna_figure_context(cfg)
orgdb <- get(cfg$orgdb_package, envir = asNamespace(cfg$orgdb_package))
input_type <- cfg$annotation$input_id_type
gsea_cfg <- cfg$annotation$enrichment$gsea
normalize_id <- function(x) nf_rna_annotation_key(x, input_type)

empty_terms <- function() data.frame(
  ID=character(), Description=character(), setSize=integer(), enrichmentScore=numeric(), NES=numeric(),
  pvalue=numeric(), p.adjust=numeric(), qvalue=numeric(), rank=integer(), leading_edge=character(), core_enrichment=character(),
  significant=logical(), stringsAsFactors=FALSE
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

run_ontology <- function(gene_list, ranked_mapping, ontology, root, contrast, lookup) {
  directory <- file.path(root, ontology); dir.create(directory, recursive=TRUE, showWarnings=FALSE)
  result <- raw <- terms <- significant <- positive <- negative <- top <- p <- core_audit <- filtered <- NULL
  # The success path already removes these; on.exit releases whatever an error left behind.
  on.exit({ rm(list=intersect(c("result", "raw", "terms", "significant", "positive", "negative", "top", "p", "core_audit", "filtered"), ls())); gc(verbose=FALSE) }, add=TRUE)
  tryCatch({
    set.seed(as.integer(gsea_cfg$seed))
    result <- suppressMessages(gseGO(
      geneList=gene_list, OrgDb=orgdb, keyType="ENTREZID", ont=ontology,
      minGSSize=as.integer(gsea_cfg$min_gs_size), maxGSSize=as.integer(gsea_cfg$max_gs_size),
      pvalueCutoff=1, pAdjustMethod=NF_RNA_ENRICHMENT_ADJUST,
      eps=0, verbose=FALSE, seed=TRUE, nPermSimple=NF_RNA_GSEA_NPERM_SIMPLE
    ))
    raw <- if (is.null(result)) data.frame() else as.data.frame(result)
    filtered <- gsea_term_tables(raw, empty_terms())
    terms <- filtered$terms; significant <- filtered$significant
    positive <- filtered$positive; negative <- filtered$negative
    write_table(terms, file.path(directory, "all_terms.tsv"))
    write_table(significant, file.path(directory, "significant.tsv"))
    write_table(positive, file.path(directory, "positive_enrichment.tsv"))
    write_table(negative, file.path(directory, "negative_enrichment.tsv"))
    core_audit <- write_gsea_core_members(terms, ranked_mapping, file.path(directory, "gsea_core_members.tsv"), as.integer(gsea_cfg$max_gs_size))
    member_record <- nf_rna_write_members(significant, nf_rna_gsea_member_pairs(result, significant, ranked_mapping), lookup, file.path(directory, "members.tsv"), gsea=TRUE)
    plot_terms <- significant[!is.na(significant$NES), , drop=FALSE]
    figures <- list(dotplot="not drawn: no significant terms (p.adjust <= 0.05)", dotplot_simplified="not drawn: no significant terms (p.adjust <= 0.05)")
    simplified_record <- list(significant_terms=0, simplified_terms=0, reason="no significant terms")
    dotplot_meta <- function(n, path, what) list(ontology=ontology, seed=as.integer(gsea_cfg$seed), significance_semantics=NF_RNA_SIGNIFICANCE_SEMANTICS,
      display_subset=sprintf("top %d %s by |NES|, ties by ID, of %d (p.adjust <= %s)", NF_RNA_DOTPLOT_TOP, what, n, NF_RNA_ENRICHMENT_PADJ),
      source_tables=list(nf_rna_source_data(path, cfg$output_dir, what), nf_rna_source_data(file.path(root, "ranked_gene_list.tsv"), cfg$output_dir, "ranked gene list")))
    if (nrow(plot_terms) > 0) {
      top <- head(plot_terms[order(-abs(plot_terms$NES), plot_terms$ID), , drop=FALSE], NF_RNA_DOTPLOT_TOP)
      nf_rna_gsea_dotplot(fig, top, contrast, NF_RNA_ENRICHMENT_PADJ, NF_RNA_ENRICHMENT_ADJUST, file.path(directory, "dotplot"), cfg$output_dir,
        dotplot_meta(nrow(plot_terms), file.path(directory, "significant.tsv"), "significant terms"), formats=NF_RNA_ENRICHMENT_FORMATS)
      simplified <- nf_rna_simplify_go(result, significant)
      write_table(simplified$table, file.path(directory, "significant_simplified.tsv"))
      simple_terms <- simplified$table[!is.na(simplified$table$NES), , drop=FALSE]
      nf_rna_gsea_dotplot(fig, head(simple_terms[order(-abs(simple_terms$NES), simple_terms$ID), , drop=FALSE], NF_RNA_DOTPLOT_TOP), contrast, NF_RNA_ENRICHMENT_PADJ, NF_RNA_ENRICHMENT_ADJUST,
        file.path(directory, "dotplot_simplified"), cfg$output_dir, dotplot_meta(nrow(simple_terms), file.path(directory, "significant_simplified.tsv"), "simplified significant terms"), formats=NF_RNA_ENRICHMENT_FORMATS)
      simplified_record <- nf_rna_simplify_record(significant, simplified)
      figures <- list(dotplot="drawn", dotplot_simplified="drawn")
      rm(plot_terms)
    }
    running <- nf_rna_running_score_plots(fig, result, significant, contrast, file.path(directory, "running_score"), cfg$output_dir, list(ontology=ontology))
    figures$running_score <- running
    summary <- list(status=if (nrow(significant) == 0) "NO_SIGNIFICANT_TERMS" else "SUCCESS", evaluated_terms=nrow(terms), all_terms=nrow(terms), significant_terms=nrow(significant), positive_terms=nrow(positive), negative_terms=nrow(negative), p_adjust_method=NF_RNA_ENRICHMENT_ADJUST, significance_rule=nf_rna_significance_record()$rule, calculation_pvalue_cutoff=1, min_gs_size=as.integer(gsea_cfg$min_gs_size), max_gs_size=as.integer(gsea_cfg$max_gs_size), na_pathways=sum(is.na(terms$pvalue) | is.na(terms$p.adjust) | is.na(terms$NES)), core_member_rows=core_audit$rows, members=member_record, simplified=simplified_record, figures=figures)
    rm(result, raw, terms, significant, positive, negative, top, p, core_audit, filtered)
    gc(verbose=FALSE)
    summary
  }, error=function(error) list(status="FAILED", reason=conditionMessage(error), evaluated_terms=0, all_terms=0, significant_terms=0, positive_terms=0, negative_terms=0, p_adjust_method=NF_RNA_ENRICHMENT_ADJUST, calculation_pvalue_cutoff=1, min_gs_size=as.integer(gsea_cfg$min_gs_size), max_gs_size=as.integer(gsea_cfg$max_gs_size), na_pathways=0))
}

contrast_summaries <- list()
blocked_reasons <- character()
for (contrast in cfg$contrasts) {
  root <- file.path(cfg$output_dir, contrast$contrast_id); dir.create(root, recursive=TRUE, showWarnings=FALSE)
  source_table <- read.delim(contrast$all_genes, check.names=FALSE, stringsAsFactors=FALSE)
  if (!all(c("gene_id", "stat") %in% names(source_table))) stop("M4A all_genes.tsv must contain gene_id and stat columns")
  lookup <- nf_rna_gene_lookup(source_table)
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
  mapping_qc <- annotation_mapping_qc(mapping_rate, cfg$annotation)
  ranking <- list(
    all_genes_rows=nrow(source_table), finite_stat_source_genes=source_count,
    nonfinite_or_missing_stat_source_genes=nrow(source_table)-source_count,
    mapped_source_genes=length(mapped_sources), unmapped_source_genes=source_count-length(mapped_sources),
    mapping_rate=mapping_rate, one_to_many_source_ids=length(unique(mapping$original_gene_id[mapping$mapping_status == "ONE_TO_MANY"])),
    duplicate_target_ids=sum(table(mapping$mapped_entrez_id[!is.na(mapping$mapped_entrez_id)]) > 1),
    duplicate_target_rows_collapsed=sum(mapping$rank_status == "COLLAPSED_DUPLICATE_TARGET"),
    final_ranked_genes=nrow(ranked), positive_stats=sum(ranked$stat > 0), negative_stats=sum(ranked$stat < 0), zero_stats=sum(ranked$stat == 0), stat_ties=sum(duplicated(ranked$stat)), tie_handling="DESeq2 stat is unmodified; ties are ordered by ascending Entrez ID before fgsea", annotation_qc=mapping_qc
  )
  rank_reason <- NULL
  if (mapping_qc$status == "BLOCKED") rank_reason <- paste0("GO preranked GSEA blocked: ", mapping_qc$reason)
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
    outcome <- run_ontology(gene_list, ranked_mapping, ontology, root, contrast, lookup)
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
annotation_contract <- list(input_id_type=cfg$annotation$input_id_type, target_id_type=cfg$annotation$target_id_type, warning_threshold=as.numeric(cfg$annotation$mapping_warning_rate), blocking_threshold=as.numeric(cfg$annotation$minimum_mapping_rate))
summary <- list(status=overall, reason=if (length(blocked_reasons)>0) paste(blocked_reasons, collapse=" ") else NULL,
                significance_semantics=NF_RNA_SIGNIFICANCE_SEMANTICS, significance=nf_rna_significance_record(),
                gsea=nf_rna_gsea_settings(gsea_cfg$seed, gsea_cfg$min_gs_size, gsea_cfg$max_gs_size),
                gene_set_source=list(kind="GO via OrgDb", orgdb=cfg$orgdb_package, orgdb_version=as.character(packageVersion(cfg$orgdb_package)), go_db_version=as.character(packageVersion("GO.db")),
                                     clusterProfiler_version=as.character(packageVersion("clusterProfiler")), snapshot="none; GO results depend on these package versions"),
                annotation_qc_status=annotation_qc_status(contrast_summaries), annotation=annotation_contract, annotation_database=cfg$orgdb_package, annotation_database_version=as.character(packageVersion(cfg$orgdb_package)), clusterProfiler_version=as.character(packageVersion("clusterProfiler")), contrasts=contrast_summaries)
write(toJSON(summary, auto_unbox=TRUE, pretty=TRUE, null="null"), file.path(cfg$output_dir, "gsea_backend_summary.json"))
nf_rna_write_figure_manifest(fig, cfg$output_dir)
nf_rna_write_provenance(cfg, cfg$output_dir, "GO_GSEA", overall, c("AnnotationDbi", "clusterProfiler", "GO.db", "GOSemSim", "DOSE", "fgsea", "enrichplot", "ggplot2", "scales", "systemfonts", "jsonlite", cfg$orgdb_package), summary)
