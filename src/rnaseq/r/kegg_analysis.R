args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2 || args[[1]] != "--config") stop("usage: kegg_analysis.R --config CONFIG.json")
script_arg <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", script_arg[grep("^--file=", script_arg)][[1]])
source(file.path(dirname(normalizePath(script_file)), "provenance.R"))
source(file.path(dirname(normalizePath(script_file)), "gene_annotation.R"))
source(file.path(dirname(normalizePath(script_file)), "ora_helpers.R"))
source(file.path(dirname(normalizePath(script_file)), "kegg_core_members.R"))
source(file.path(dirname(normalizePath(script_file)), "annotation_mapping_qc.R"))
source(file.path(dirname(normalizePath(script_file)), "gsea_term_filtering.R"))
source(file.path(dirname(normalizePath(script_file)), "figures.R"))
source(file.path(dirname(normalizePath(script_file)), "enrichment_common.R"))
source(file.path(dirname(normalizePath(script_file)), "kegg_snapshot.R"))
suppressPackageStartupMessages({ library(jsonlite); library(AnnotationDbi); library(clusterProfiler); library(ggplot2) })
cfg <- fromJSON(args[[2]], simplifyVector = FALSE)
dir.create(cfg$output_dir, recursive=TRUE, showWarnings=FALSE)
fig <- nf_rna_figure_context(cfg)
orgdb <- get(cfg$orgdb_package, envir=asNamespace(cfg$orgdb_package))
input_type <- cfg$annotation$input_id_type
normalize_id <- function(x) nf_rna_annotation_key(x, input_type)
write_table <- function(x, path) write.table(x, path, sep="\t", quote=FALSE, row.names=FALSE, na="NA")
packages <- c("AnnotationDbi", "clusterProfiler", "DOSE", "fgsea", "enrichplot", "KEGGREST", "openssl", "ggplot2", "scales", "systemfonts", "jsonlite", cfg$orgdb_package)
module <- if (cfg$mode == "ora") "KEGG_ORA" else "KEGG_GSEA"
summary_name <- if (cfg$mode == "ora") "kegg_backend_summary.json" else "gsea_kegg_backend_summary.json"
ora_cfg <- cfg$annotation$enrichment$kegg$ora; gsea_cfg <- cfg$annotation$enrichment$kegg$gsea

map_sources <- function(source, stat=NULL) {
  normalized <- normalize_id(source); valid <- AnnotationDbi::keys(orgdb, keytype=input_type)
  requested <- intersect(unique(normalized), valid)
  selected <- if (length(requested) == 0) data.frame() else suppressMessages(AnnotationDbi::select(orgdb, keys=requested, keytype=input_type, columns=c("ENTREZID", "SYMBOL")))
  if (nrow(selected) == 0) selected <- data.frame(key=character(), ENTREZID=character(), SYMBOL=character(), stringsAsFactors=FALSE)
  names(selected)[[1]] <- input_type; selected <- selected[!is.na(selected$ENTREZID) & selected$ENTREZID != "", , drop=FALSE]
  rows <- lapply(seq_along(source), function(i) {
    hits <- selected[selected[[input_type]] == normalized[[i]], c("ENTREZID", "SYMBOL"), drop=FALSE]
    targets <- unique(hits$ENTREZID); status <- if (length(targets)==0) "UNMAPPED" else if (length(targets)==1) "MAPPED_UNIQUE" else "ONE_TO_MANY"
    value <- if (is.null(stat)) NA_real_ else stat[[i]]
    if (nrow(hits)==0) data.frame(original_gene_id=source[[i]], normalized_gene_id=normalized[[i]], mapped_entrez_id=NA_character_, mapped_symbol=NA_character_, stat=value, mapping_status=status, stringsAsFactors=FALSE) else data.frame(original_gene_id=source[[i]], normalized_gene_id=normalized[[i]], mapped_entrez_id=as.character(hits$ENTREZID), mapped_symbol=as.character(hits$SYMBOL), stat=value, mapping_status=status, stringsAsFactors=FALSE)
  })
  do.call(rbind, rows)
}

mapping_stats <- function(mapping, source) {
  mapped_sources <- unique(mapping$original_gene_id[!is.na(mapping$mapped_entrez_id)])
  targets <- unique(mapping$mapped_entrez_id[!is.na(mapping$mapped_entrez_id)])
  list(source_genes=length(source), mapped_source_genes=length(mapped_sources), unmapped_source_genes=length(source)-length(mapped_sources), unique_target_genes=length(targets), one_to_many_source_ids=length(unique(mapping$original_gene_id[mapping$mapping_status=="ONE_TO_MANY"])), duplicate_target_ids=sum(table(mapping$mapped_entrez_id[!is.na(mapping$mapped_entrez_id)]) > 1), mapping_rate=if(length(source)==0) 0 else length(mapped_sources)/length(source), targets=as.character(targets))
}
without_targets <- function(x) { x$targets <- NULL; x }
empty_ora <- function() data.frame(ID=character(), Description=character(), GeneRatio=character(), BgRatio=character(), pvalue=numeric(), p.adjust=numeric(), qvalue=numeric(), geneID=character(), Count=integer(), significant=logical(), stringsAsFactors=FALSE)
empty_gsea <- function() data.frame(ID=character(), Description=character(), setSize=integer(), enrichmentScore=numeric(), NES=numeric(), pvalue=numeric(), p.adjust=numeric(), qvalue=numeric(), rank=integer(), leading_edge=character(), core_enrichment=character(), significant=logical(), stringsAsFactors=FALSE)
common_summary <- function(snapshot_record) list(significance_semantics=NF_RNA_SIGNIFICANCE_SEMANTICS, significance=nf_rna_significance_record(), gene_set_source=snapshot_record,
                                                 clusterProfiler_version=as.character(packageVersion("clusterProfiler")), r_version=R.version.string)

# One snapshot for this analysis; a failure is reported, never replaced by a live query.
snapshot <- tryCatch(nf_rna_resolve_kegg_snapshot(cfg$kegg$snapshot_root, cfg$kegg$organism_code, cfg$kegg$snapshot), error=function(e) e)
if (inherits(snapshot, "error")) {
  reason <- paste0("KEGG gene-set snapshot unavailable: ", conditionMessage(snapshot))
  message(reason)
  resource <- list(provider="kegg_rest_snapshot", status="NETWORK_UNAVAILABLE", requested_snapshot=cfg$kegg$snapshot, error=conditionMessage(snapshot),
                   reproducible=FALSE, reproducibility="NOT_FULLY_REPRODUCIBLE: no verified KEGG snapshot was available; no KEGG result was produced")
  summary <- c(list(status="NETWORK_UNAVAILABLE", reason=reason, resource=resource), common_summary(resource), list(contrasts=list()))
  write(toJSON(summary, auto_unbox=TRUE, pretty=TRUE, null="null"), file.path(cfg$output_dir, summary_name))
  nf_rna_write_figure_manifest(fig, cfg$output_dir); nf_rna_write_provenance(cfg, cfg$output_dir, module, "NETWORK_UNAVAILABLE", packages, summary); quit(status=0)
}
snapshot_record <- nf_rna_snapshot_record(snapshot)
file.copy(file.path(snapshot$directory, "manifest.json"), file.path(cfg$output_dir, "kegg_snapshot_manifest.json"), overwrite=TRUE)
term2gene <- snapshot$term2gene; term2name <- snapshot$term2name
kegg_genes <- unique(term2gene$gene)

if (cfg$mode == "ora") {
  min_gs <- as.integer(ora_cfg$min_gs_size); max_gs <- as.integer(ora_cfg$max_gs_size)
  run_fg <- function(ids, name, root, mapping, tested, universe, contrast, lookup) {
    ids <- unique(as.character(ids[!is.na(ids) & ids!=""])); if (!all(ids %in% tested)) stop(paste0("KEGG ORA foreground ", name, " is not a subset of the padj-non-NA genes"))
    submap <- mapping[mapping$original_gene_id %in% ids,,drop=FALSE]; stats <- mapping_stats(submap, ids); targets <- stats$targets; stats <- without_targets(stats); dir <- file.path(root,name); dir.create(dir,recursive=TRUE,showWarnings=FALSE)
    head_fields <- list(foreground_significant_genes=length(ids), successfully_mapped_foreground_genes=length(targets))
    if (length(targets) < cfg$annotation$minimum_mapped_foreground) {
      write_table(empty_ora(),file.path(dir,"kegg_all_terms.tsv")); write_table(empty_ora(),file.path(dir,"kegg.tsv")); write_kegg_member_table(empty_ora(), mapping, "geneID", file.path(dir,"kegg_gene_members.tsv"))
      nf_rna_write_members(empty_ora(), nf_rna_ora_member_pairs(empty_ora(), submap), NULL, file.path(dir, "kegg_members.tsv"))
      return(c(head_fields, stats, list(status="NOT_APPLICABLE", reason=paste0("only ",length(targets)," mapped foreground genes"), tested_terms=0, significant_terms=0, figures=list(dotplot="not drawn: too few mapped foreground genes"))))
    }
    result <- suppressMessages(enricher(gene=targets, universe=universe$genes, TERM2GENE=term2gene, TERM2NAME=term2name, minGSSize=min_gs, maxGSSize=max_gs, pAdjustMethod=NF_RNA_ENRICHMENT_ADJUST, pvalueCutoff=1, qvalueCutoff=1))
    terms <- if (is.null(result)) empty_ora() else nf_rna_flag_significant(as.data.frame(result))
    if (nrow(terms)) terms <- terms[order(terms$p.adjust, terms$ID),,drop=FALSE]
    significant <- terms[terms$significant,,drop=FALSE]
    write_table(terms,file.path(dir,"kegg_all_terms.tsv")); write_table(significant,file.path(dir,"kegg.tsv")); write_kegg_member_table(significant,mapping,"geneID",file.path(dir,"kegg_gene_members.tsv"))
    member_record <- nf_rna_write_members(significant, nf_rna_ora_member_pairs(significant, submap), lookup, file.path(dir, "kegg_members.tsv"))
    figures <- list(dotplot="not drawn: no significant pathways (p.adjust <= 0.05)")
    if (nrow(significant)) {
      nf_rna_ora_dotplot(fig, head(significant, NF_RNA_DOTPLOT_TOP), sprintf("-log10(KEGG adjusted P, %s)", NF_RNA_ENRICHMENT_ADJUST), file.path(dir,"dotplot"), cfg$output_dir,
        list(contrast=list(contrast_id=contrast$contrast_id, numerator=contrast$numerator, denominator=contrast$denominator), foreground=name, significance_semantics=NF_RNA_SIGNIFICANCE_SEMANTICS, kegg_snapshot=snapshot$id,
             display_subset=sprintf("top %d significant pathways by adjusted P, ties by ID, of %d significant (p.adjust <= %s)", NF_RNA_DOTPLOT_TOP, nrow(significant), NF_RNA_ENRICHMENT_PADJ),
             source_tables=list(nf_rna_source_data(file.path(dir,"kegg.tsv"), cfg$output_dir, "all significant pathways"))),
        subtitle=nf_rna_ora_subtitle(fig, contrast, name), formats=NF_RNA_ENRICHMENT_FORMATS)
      figures$dotplot <- "drawn"
    }
    universe_record <- c(universe$record, list(foreground_entrez_ids=length(targets), foreground_in_universe=length(intersect(targets, universe$genes)), clusterprofiler_bgratio_n=nf_rna_bgratio_n(terms)))
    c(head_fields, stats, list(status=if(nrow(significant))"SUCCESS"else"NO_SIGNIFICANT_TERMS", tested_terms=nrow(terms), significant_terms=nrow(significant), universe=universe_record, members=member_record, figures=figures))
  }
  contrast_summaries <- list()
  for (contrast in cfg$contrasts) {
    root <- file.path(cfg$output_dir,contrast$contrast_id); dir.create(root,recursive=TRUE,showWarnings=FALSE)
    tab <- read.delim(contrast$all_genes,check.names=FALSE,stringsAsFactors=FALSE); lookup <- nf_rna_gene_lookup(tab)
    gene_sets <- nf_rna_ora_gene_sets(tab); retained <- gene_sets$retained; tested <- gene_sets$tested
    mapping <- map_sources(retained); write_table(mapping,file.path(root,"kegg_mapping.tsv"))
    tested_mapping <- mapping[mapping$original_gene_id%in%tested,,drop=FALSE]; stats <- mapping_stats(tested_mapping,tested)
    universe <- nf_rna_ora_universe(tested, tested_mapping, kegg_genes)
    mapping_qc <- annotation_mapping_qc(stats$mapping_rate,cfg$annotation)
    universe_summary <- c(gene_sets$counts, list(successfully_mapped_tested_genes=length(stats$targets)), without_targets(stats), list(annotation_qc=mapping_qc, kegg=universe$record))
    fgs <- list(significant=contrast$significant,up=contrast$up,down=contrast$down)
    values <- lapply(names(fgs), function(n) {
      if (mapping_qc$status=="BLOCKED") list(status="BLOCKED", reason=paste0("KEGG ORA blocked: ",mapping_qc$reason))
      else run_fg(read.delim(fgs[[n]],check.names=FALSE,stringsAsFactors=FALSE)$gene_id, n, root, tested_mapping, tested, universe, contrast, lookup)
    })
    names(values) <- names(fgs)
    contrast_summaries[[length(contrast_summaries)+1]] <- list(contrast_id=contrast$contrast_id, universe=universe_summary, universe_size=length(universe$genes), foregrounds=values)
  }
  statuses <- unlist(lapply(contrast_summaries,function(x)vapply(x$foregrounds,function(y)y$status,character(1))))
  overall <- if(any(statuses=="BLOCKED"))"BLOCKED"else if(all(statuses=="NOT_APPLICABLE"))"NOT_APPLICABLE"else if(all(statuses!="SUCCESS"))"NO_SIGNIFICANT_TERMS"else"SUCCESS"
  annotation_contract <- list(input_id_type=cfg$annotation$input_id_type,target_id_type=cfg$annotation$target_id_type,warning_threshold=as.numeric(cfg$annotation$mapping_warning_rate),blocking_threshold=as.numeric(cfg$annotation$minimum_mapping_rate))
  summary <- c(list(status=overall), common_summary(snapshot_record),
               list(ora=list(foregrounds="significant = significant.tsv, up = upregulated.tsv (log2FoldChange > 0), down = downregulated.tsv (log2FoldChange < 0); the unchanged L2 DEG definition",
                             universe="padj-non-NA genes mapped to Entrez, intersected with the genes of the KEGG snapshot", min_gs_size=min_gs, max_gs_size=max_gs),
                    annotation_qc_status=annotation_qc_status(contrast_summaries), annotation=annotation_contract, annotation_database=cfg$orgdb_package, annotation_database_version=as.character(packageVersion(cfg$orgdb_package)),
                    resource=snapshot_record, contrasts=contrast_summaries))
} else {
  min_gs <- as.integer(gsea_cfg$min_gs_size); max_gs <- as.integer(gsea_cfg$max_gs_size)
  contrast_summaries <- list(); blocked <- character()
  for (contrast in cfg$contrasts) {
    root <- file.path(cfg$output_dir,contrast$contrast_id); dir.create(root,recursive=TRUE,showWarnings=FALSE)
    source_table <- read.delim(contrast$all_genes,check.names=FALSE,stringsAsFactors=FALSE); if(!all(c("gene_id","stat") %in% names(source_table)))stop("M4A all_genes.tsv must contain gene_id and stat")
    lookup <- nf_rna_gene_lookup(source_table)
    stat <- suppressWarnings(as.numeric(source_table$stat)); finite <- is.finite(stat)&!is.na(source_table$gene_id)&source_table$gene_id!=""
    mapping <- map_sources(as.character(source_table$gene_id[finite]),stat[finite]); mapping$target_duplicate <- FALSE; mapping$rank_status <- ifelse(is.na(mapping$mapped_entrez_id),"UNMAPPED","CANDIDATE")
    mapped <- which(!is.na(mapping$mapped_entrez_id))
    if (length(mapped)>0) for (target in unique(mapping$mapped_entrez_id[mapped])) {
      idx <- mapped[mapping$mapped_entrez_id[mapped]==target]; cand <- mapping[idx,,drop=FALSE]
      winner <- idx[[order(-abs(cand$stat),cand$original_gene_id,cand$normalized_gene_id)[[1]]]]
      mapping$rank_status[idx] <- "COLLAPSED_DUPLICATE_TARGET"; mapping$rank_status[winner] <- "RETAINED"
    }
    write_table(mapping,file.path(root,"gsea_kegg_rank_mapping.tsv"))
    retained <- mapping[mapping$rank_status=="RETAINED",,drop=FALSE]; retained <- retained[order(-retained$stat,retained$mapped_entrez_id),,drop=FALSE]
    ranked <- data.frame(rank=seq_len(nrow(retained)),entrez_id=retained$mapped_entrez_id,stat=retained$stat,original_gene_id=retained$original_gene_id,symbol=retained$mapped_symbol,stringsAsFactors=FALSE); write_table(ranked,file.path(root,"ranked_gene_list.tsv"))
    mapped_sources <- unique(mapping$original_gene_id[!is.na(mapping$mapped_entrez_id)]); rate <- if(sum(finite)==0)0 else length(mapped_sources)/sum(finite); mapping_qc <- annotation_mapping_qc(rate,cfg$annotation)
    ranking <- list(all_genes_rows=nrow(source_table),finite_stat_source_genes=sum(finite),nonfinite_or_missing_stat_source_genes=nrow(source_table)-sum(finite),mapped_source_genes=length(mapped_sources),unmapped_source_genes=sum(finite)-length(mapped_sources),mapping_rate=rate,one_to_many_source_ids=length(unique(mapping$original_gene_id[mapping$mapping_status=="ONE_TO_MANY"])),duplicate_target_ids=sum(table(mapping$mapped_entrez_id[!is.na(mapping$mapped_entrez_id)])>1),duplicate_target_rows_collapsed=sum(mapping$rank_status=="COLLAPSED_DUPLICATE_TARGET"),final_ranked_genes=nrow(ranked),ranked_genes_in_kegg_snapshot=sum(ranked$entrez_id %in% kegg_genes),positive_stats=sum(ranked$stat>0),negative_stats=sum(ranked$stat<0),zero_stats=sum(ranked$stat==0),stat_ties=sum(duplicated(ranked$stat)),tie_handling="DESeq2 stat is unmodified; ties are ordered by ascending Entrez ID before fgsea",annotation_qc=mapping_qc)
    reason <- NULL
    if (mapping_qc$status=="BLOCKED") reason <- paste0("KEGG GSEA blocked: ",mapping_qc$reason)
    if (is.null(reason) && nrow(ranked) < as.integer(gsea_cfg$minimum_ranked_genes)) reason <- paste0("KEGG GSEA blocked: ",nrow(ranked)," unique ranked Entrez genes; minimum is ",as.integer(gsea_cfg$minimum_ranked_genes),".")
    if (!is.null(reason)) { contrast_summaries[[length(contrast_summaries)+1]] <- list(contrast_id=contrast$contrast_id,ranking=ranking,status="BLOCKED",reason=reason); blocked <- c(blocked,reason); next }
    ranked_mapping <- mapping[mapping$rank_status=="RETAINED" & !is.na(mapping$mapped_entrez_id),,drop=FALSE]
    gene_list <- ranked$stat; names(gene_list) <- ranked$entrez_id
    set.seed(as.integer(gsea_cfg$seed))
    result <- suppressMessages(GSEA(geneList=gene_list, TERM2GENE=term2gene, TERM2NAME=term2name, minGSSize=min_gs, maxGSSize=max_gs, eps=0, pvalueCutoff=1, pAdjustMethod=NF_RNA_ENRICHMENT_ADJUST, verbose=FALSE, seed=TRUE, nPermSimple=NF_RNA_GSEA_NPERM_SIMPLE))
    raw <- if (is.null(result)) data.frame() else as.data.frame(result)
    filtered <- gsea_term_tables(raw, empty_gsea()); terms <- filtered$terms; significant <- filtered$significant; positive <- filtered$positive; negative <- filtered$negative
    write_table(terms,file.path(root,"all_terms.tsv")); write_table(significant,file.path(root,"significant.tsv")); write_table(positive,file.path(root,"positive_enrichment.tsv")); write_table(negative,file.path(root,"negative_enrichment.tsv"))
    core_audit <- write_kegg_member_table(terms,mapping,"core_enrichment",file.path(root,"core_members.tsv"),strict=TRUE,max_members=max_gs)
    member_record <- nf_rna_write_members(significant, nf_rna_gsea_member_pairs(result, significant, ranked_mapping), lookup, file.path(root, "members.tsv"), gsea=TRUE)
    figures <- list(dotplot="not drawn: no significant pathways (p.adjust <= 0.05)")
    if (nrow(significant)>0) {
      top <- head(significant[order(-abs(significant$NES),significant$ID),,drop=FALSE], NF_RNA_DOTPLOT_TOP)
      nf_rna_gsea_dotplot(fig, top, contrast, NF_RNA_ENRICHMENT_PADJ, NF_RNA_ENRICHMENT_ADJUST, file.path(root,"dotplot"), cfg$output_dir,
        list(seed=as.integer(gsea_cfg$seed), significance_semantics=NF_RNA_SIGNIFICANCE_SEMANTICS, kegg_snapshot=snapshot$id,
             display_subset=sprintf("top %d significant pathways by |NES|, ties by ID, of %d significant (p.adjust <= %s)", NF_RNA_DOTPLOT_TOP, nrow(significant), NF_RNA_ENRICHMENT_PADJ),
             source_tables=list(nf_rna_source_data(file.path(root,"significant.tsv"), cfg$output_dir, "all significant pathways"), nf_rna_source_data(file.path(root,"ranked_gene_list.tsv"), cfg$output_dir, "ranked gene list"))),
        formats=NF_RNA_ENRICHMENT_FORMATS)
      figures$dotplot <- "drawn"
    }
    figures$running_score <- nf_rna_running_score_plots(fig, result, significant, contrast, file.path(root, "running_score"), cfg$output_dir, list(kegg_snapshot=snapshot$id))
    contrast_summaries[[length(contrast_summaries)+1]] <- list(contrast_id=contrast$contrast_id,ranking=ranking,status=if(nrow(significant)==0)"NO_SIGNIFICANT_TERMS"else"SUCCESS",evaluated_terms=nrow(terms),all_terms=nrow(terms),significant_terms=nrow(significant),positive_terms=nrow(positive),negative_terms=nrow(negative),p_adjust_method=NF_RNA_ENRICHMENT_ADJUST,significance_rule=nf_rna_significance_record()$rule,calculation_pvalue_cutoff=1,min_gs_size=min_gs,max_gs_size=max_gs,core_member_rows=core_audit$rows,members=member_record,figures=figures)
  }
  overall <- if(length(blocked)>0)"BLOCKED"else if(all(vapply(contrast_summaries,function(x)x$status=="NO_SIGNIFICANT_TERMS",logical(1))))"NO_SIGNIFICANT_TERMS"else"SUCCESS"
  annotation_contract <- list(input_id_type=cfg$annotation$input_id_type,target_id_type=cfg$annotation$target_id_type,warning_threshold=as.numeric(cfg$annotation$mapping_warning_rate),blocking_threshold=as.numeric(cfg$annotation$minimum_mapping_rate))
  summary <- c(list(status=overall, reason=if(length(blocked)>0)paste(blocked,collapse=" ")else NULL), common_summary(snapshot_record),
               list(gsea=nf_rna_gsea_settings(gsea_cfg$seed, min_gs, max_gs), annotation_qc_status=annotation_qc_status(contrast_summaries), annotation=annotation_contract, annotation_database=cfg$orgdb_package,
                    annotation_database_version=as.character(packageVersion(cfg$orgdb_package)), resource=snapshot_record, contrasts=contrast_summaries))
}
write(toJSON(summary,auto_unbox=TRUE,pretty=TRUE,null="null"),file.path(cfg$output_dir,summary_name))
nf_rna_write_figure_manifest(fig,cfg$output_dir)
nf_rna_write_provenance(cfg,cfg$output_dir,module,summary$status,packages,summary)
