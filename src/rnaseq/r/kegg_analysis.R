args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2 || args[[1]] != "--config") stop("usage: kegg_analysis.R --config CONFIG.json")
script_arg <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", script_arg[grep("^--file=", script_arg)][[1]])
source(file.path(dirname(normalizePath(script_file)), "kegg_core_members.R"))
suppressPackageStartupMessages({ library(jsonlite); library(AnnotationDbi); library(clusterProfiler); library(ggplot2) })
cfg <- fromJSON(args[[2]], simplifyVector = FALSE)
dir.create(cfg$output_dir, recursive=TRUE, showWarnings=FALSE)
orgdb <- get(cfg$orgdb_package, envir=asNamespace(cfg$orgdb_package))
input_type <- cfg$annotation$input_id_type
normalize_id <- function(x) if (input_type == "ENSEMBL") sub("\\.[0-9]+$", "", x) else x
write_table <- function(x, path) write.table(x, path, sep="\t", quote=FALSE, row.names=FALSE, na="NA")
utc_now <- function() format(Sys.time(), tz="UTC", usetz=TRUE)

probe_resource <- function() {
  tryCatch({
    con <- url(cfg$kegg$probe_endpoint, open="rt", encoding="UTF-8"); on.exit(close(con), add=TRUE)
    first <- readLines(con, n=1, warn=FALSE)
    if (length(first) == 0) stop("KEGG REST returned no pathway records")
    list(status="AVAILABLE", provider=cfg$kegg$provider, retrieval_timestamp=utc_now(), probe_endpoint=cfg$kegg$probe_endpoint, metadata=list(first_pathway_record=first[[1]]))
  }, error=function(e) list(status="NETWORK_UNAVAILABLE", provider=cfg$kegg$provider, retrieval_timestamp=utc_now(), probe_endpoint=cfg$kegg$probe_endpoint, error=conditionMessage(e)))
}

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
empty_ora <- function() data.frame(ID=character(), Description=character(), GeneRatio=character(), BgRatio=character(), pvalue=numeric(), p.adjust=numeric(), qvalue=numeric(), geneID=character(), Count=integer(), stringsAsFactors=FALSE)
empty_gsea <- function() data.frame(ID=character(), Description=character(), setSize=integer(), enrichmentScore=numeric(), NES=numeric(), pvalue=numeric(), p.adjust=numeric(), qvalue=numeric(), rank=integer(), leading_edge=character(), core_enrichment=character(), stringsAsFactors=FALSE)

network_summary <- function(resource) list(status="NETWORK_UNAVAILABLE", reason=paste0("KEGG resource unavailable: ", resource$error), resource=resource, clusterProfiler_version=as.character(packageVersion("clusterProfiler")), r_version=R.version.string, contrasts=list())
resource <- probe_resource()
if (resource$status != "AVAILABLE") {
  message("KEGG resource unavailable: ", resource$error)
  summary <- network_summary(resource); name <- if(cfg$mode=="ora") "kegg_backend_summary.json" else "gsea_kegg_backend_summary.json"; write(toJSON(summary, auto_unbox=TRUE, pretty=TRUE, null="null"), file.path(cfg$output_dir, name)); quit(status=0)
}

if (cfg$mode == "ora") {
  all <- read.delim(cfg$contrasts[[1]]$all_genes, check.names=FALSE, stringsAsFactors=FALSE)$gene_id
  source <- unique(as.character(all)); mapping <- map_sources(source); write_table(mapping, file.path(cfg$output_dir, "kegg_mapping.tsv"))
  tested <- mapping_stats(mapping, source); universe <- unique(tested$targets); mapping_summary <- without_targets(tested)
  if (tested$mapping_rate < cfg$annotation$minimum_mapping_rate) {
    summary <- list(status="BLOCKED", reason=paste0("KEGG ORA blocked: mapped ", sprintf("%.1f%%",100*tested$mapping_rate), " of tested genes; required minimum is ", sprintf("%.1f%%",100*cfg$annotation$minimum_mapping_rate), "."), resource=resource, clusterProfiler_version=as.character(packageVersion("clusterProfiler")), r_version=R.version.string, mapping=mapping_summary, contrasts=list()); write(toJSON(summary,auto_unbox=TRUE,pretty=TRUE),file.path(cfg$output_dir,"kegg_backend_summary.json")); quit(status=0)
  }
  run_fg <- function(ids, name, root) {
    submap <- mapping[mapping$original_gene_id %in% ids,,drop=FALSE]; targets <- unique(submap$mapped_entrez_id[!is.na(submap$mapped_entrez_id)]); dir <- file.path(root,name); dir.create(dir,recursive=TRUE,showWarnings=FALSE)
    stats <- mapping_stats(submap, unique(as.character(ids))); stats$targets <- NULL
    if (length(targets) < cfg$annotation$minimum_mapped_foreground) { write_table(empty_ora(),file.path(dir,"kegg.tsv")); write_kegg_member_table(empty_ora(), mapping, "geneID", file.path(dir,"kegg_gene_members.tsv")); return(c(stats,list(status="NOT_APPLICABLE",reason=paste0("only ",length(targets)," mapped foreground genes"),pathway_count=0))) }
    if (!all(targets %in% universe)) stop("KEGG ORA foreground is not a subset of the mapped tested-gene universe")
    result <- tryCatch(suppressMessages(enrichKEGG(gene=targets, organism=cfg$kegg$organism_code, keyType="ncbi-geneid", universe=universe, pvalueCutoff=cfg$annotation$enrichment$kegg$ora$pvalue_cutoff, pAdjustMethod=cfg$annotation$enrichment$kegg$ora$p_adjust_method, qvalueCutoff=cfg$annotation$enrichment$kegg$ora$qvalue_cutoff, minGSSize=as.integer(cfg$annotation$enrichment$kegg$ora$min_gs_size), maxGSSize=as.integer(cfg$annotation$enrichment$kegg$ora$max_gs_size), use_internal_data=FALSE)), error=function(e)e)
    if (inherits(result,"error")) return(c(stats,list(status="NETWORK_UNAVAILABLE",reason=paste0("KEGG service error after successful probe: ",conditionMessage(result)),pathway_count=0)))
    terms <- if(is.null(result)) empty_ora() else as.data.frame(result); if(nrow(terms)>0) terms <- terms[order(terms$p.adjust,terms$ID),,drop=FALSE]
    write_table(terms,file.path(dir,"kegg.tsv")); write_kegg_member_table(terms,mapping,"geneID",file.path(dir,"kegg_gene_members.tsv"))
    if(nrow(terms)>0) { top<-head(terms,15); top$Description<-factor(top$Description,levels=rev(top$Description)); p<-ggplot(top,aes(x=Description,y=-log10(p.adjust),size=Count))+geom_point(color="#2166ac")+coord_flip()+theme_minimal()+labs(x=NULL,y="-log10(KEGG adjusted p-value)"); ggsave(file.path(dir,"dotplot.png"),p,width=7,height=5,dpi=150); ggsave(file.path(dir,"dotplot.tiff"),p,width=7,height=5,dpi=300,compression="lzw") }
    c(stats,list(status=if(nrow(terms)==0) "NO_SIGNIFICANT_TERMS" else "SUCCESS", pathway_count=nrow(terms)))
  }
  contrast_summaries<-list(); for(contrast in cfg$contrasts) { root<-file.path(cfg$output_dir,contrast$contrast_id); dir.create(root,recursive=TRUE,showWarnings=FALSE); fgs<-list(significant=contrast$significant,up=contrast$up,down=contrast$down); values<-lapply(names(fgs),function(n)run_fg(read.delim(fgs[[n]],check.names=FALSE,stringsAsFactors=FALSE)$gene_id,n,root)); names(values)<-names(fgs); contrast_summaries[[length(contrast_summaries)+1]]<-list(contrast_id=contrast$contrast_id,universe_size=length(universe),foregrounds=values) }
  overall<-if(any(vapply(contrast_summaries,function(x)any(vapply(x$foregrounds,function(y)y$status=="NETWORK_UNAVAILABLE",logical(1))),logical(1)))) "NETWORK_UNAVAILABLE" else if(all(vapply(contrast_summaries,function(x)all(vapply(x$foregrounds,function(y)y$status!="SUCCESS",logical(1))),logical(1)))) "NO_SIGNIFICANT_TERMS" else "SUCCESS"; summary<-list(status=overall,resource=resource,clusterProfiler_version=as.character(packageVersion("clusterProfiler")),r_version=R.version.string,mapping=mapping_summary,contrasts=contrast_summaries); write(toJSON(summary,auto_unbox=TRUE,pretty=TRUE,null="null"),file.path(cfg$output_dir,"kegg_backend_summary.json"))
} else {
  contrast_summaries<-list(); blocked<-character()
  for(contrast in cfg$contrasts) {
    root<-file.path(cfg$output_dir,contrast$contrast_id);dir.create(root,recursive=TRUE,showWarnings=FALSE); source_table<-read.delim(contrast$all_genes,check.names=FALSE,stringsAsFactors=FALSE); if(!all(c("gene_id","stat") %in% names(source_table)))stop("M4A all_genes.tsv must contain gene_id and stat")
    stat<-suppressWarnings(as.numeric(source_table$stat)); finite<-is.finite(stat)&!is.na(source_table$gene_id)&source_table$gene_id!=""; mapping<-map_sources(as.character(source_table$gene_id[finite]),stat[finite]); mapping$target_duplicate<-FALSE;mapping$rank_status<-ifelse(is.na(mapping$mapped_entrez_id),"UNMAPPED","CANDIDATE"); mapped<-which(!is.na(mapping$mapped_entrez_id));if(length(mapped)>0)for(target in unique(mapping$mapped_entrez_id[mapped])){idx<-mapped[mapping$mapped_entrez_id[mapped]==target];cand<-mapping[idx,,drop=FALSE];winner<-idx[[order(-abs(cand$stat),cand$original_gene_id,cand$normalized_gene_id)[[1]]]];mapping$rank_status[idx]<-"COLLAPSED_DUPLICATE_TARGET";mapping$rank_status[winner]<-"RETAINED"}
    write_table(mapping,file.path(root,"gsea_kegg_rank_mapping.tsv")); retained<-mapping[mapping$rank_status=="RETAINED",,drop=FALSE]; retained<-retained[order(-retained$stat,retained$mapped_entrez_id),,drop=FALSE]; ranked<-data.frame(rank=seq_len(nrow(retained)),entrez_id=retained$mapped_entrez_id,stat=retained$stat,original_gene_id=retained$original_gene_id,symbol=retained$mapped_symbol,stringsAsFactors=FALSE);write_table(ranked,file.path(root,"ranked_gene_list.tsv"))
    mapped_sources<-unique(mapping$original_gene_id[!is.na(mapping$mapped_entrez_id)]);rate<-if(sum(finite)==0)0 else length(mapped_sources)/sum(finite); ranking<-list(all_genes_rows=nrow(source_table),finite_stat_source_genes=sum(finite),nonfinite_or_missing_stat_source_genes=nrow(source_table)-sum(finite),mapped_source_genes=length(mapped_sources),unmapped_source_genes=sum(finite)-length(mapped_sources),mapping_rate=rate,one_to_many_source_ids=length(unique(mapping$original_gene_id[mapping$mapping_status=="ONE_TO_MANY"])),duplicate_target_ids=sum(table(mapping$mapped_entrez_id[!is.na(mapping$mapped_entrez_id)])>1),duplicate_target_rows_collapsed=sum(mapping$rank_status=="COLLAPSED_DUPLICATE_TARGET"),final_ranked_genes=nrow(ranked),positive_stats=sum(ranked$stat>0),negative_stats=sum(ranked$stat<0),zero_stats=sum(ranked$stat==0),stat_ties=sum(duplicated(ranked$stat)))
    reason<-NULL;if(rate<cfg$annotation$minimum_mapping_rate)reason<-paste0("KEGG GSEA blocked: mapped ",sprintf("%.1f%%",100*rate)," of finite tested genes; required minimum is ",sprintf("%.1f%%",100*cfg$annotation$minimum_mapping_rate),".");if(is.null(reason)&&nrow(ranked)<as.integer(cfg$annotation$enrichment$kegg$gsea$minimum_ranked_genes))reason<-paste0("KEGG GSEA blocked: ",nrow(ranked)," unique ranked Entrez genes; minimum is ",as.integer(cfg$annotation$enrichment$kegg$gsea$minimum_ranked_genes),".")
    if(!is.null(reason)){contrast_summaries[[length(contrast_summaries)+1]]<-list(contrast_id=contrast$contrast_id,ranking=ranking,status="BLOCKED",reason=reason);blocked<-c(blocked,reason);next}
    gene_list<-ranked$stat;names(gene_list)<-ranked$entrez_id; result<-tryCatch(suppressMessages(gseKEGG(geneList=gene_list,organism=cfg$kegg$organism_code,keyType="ncbi-geneid",pvalueCutoff=cfg$annotation$enrichment$kegg$gsea$pvalue_cutoff,pAdjustMethod=cfg$annotation$enrichment$kegg$gsea$p_adjust_method,minGSSize=as.integer(cfg$annotation$enrichment$kegg$gsea$min_gs_size),maxGSSize=as.integer(cfg$annotation$enrichment$kegg$gsea$max_gs_size),eps=0,verbose=FALSE,seed=TRUE)),error=function(e)e);if(inherits(result,"error")){reason<-paste0("KEGG service error after successful probe: ",conditionMessage(result));contrast_summaries[[length(contrast_summaries)+1]]<-list(contrast_id=contrast$contrast_id,ranking=ranking,status="NETWORK_UNAVAILABLE",reason=reason);next}
    raw<-if(is.null(result))data.frame()else as.data.frame(result); wanted<-names(empty_gsea());if(nrow(raw)==0)terms<-empty_gsea()else{for(n in setdiff(wanted,names(raw)))raw[[n]]<-NA;terms<-raw[,wanted,drop=FALSE];terms<-terms[order(terms$p.adjust,terms$ID),,drop=FALSE]};significant<-terms[!is.na(terms$p.adjust)&terms$p.adjust<=as.numeric(cfg$annotation$enrichment$kegg$gsea$padj_cutoff)&!is.na(terms$pvalue)&terms$pvalue<=as.numeric(cfg$annotation$enrichment$kegg$gsea$pvalue_cutoff),,drop=FALSE];positive<-significant[significant$NES>0,,drop=FALSE];negative<-significant[significant$NES<0,,drop=FALSE];write_table(terms,file.path(root,"all_terms.tsv"));write_table(significant,file.path(root,"significant.tsv"));write_table(positive,file.path(root,"positive_enrichment.tsv"));write_table(negative,file.path(root,"negative_enrichment.tsv"));core_audit<-write_kegg_member_table(terms,mapping,"core_enrichment",file.path(root,"core_members.tsv"),strict=TRUE,max_members=as.integer(cfg$annotation$enrichment$kegg$gsea$max_gs_size));if(nrow(significant)>0){top<-head(significant[order(-abs(significant$NES),significant$ID),,drop=FALSE],15);top$Description<-factor(top$Description,levels=rev(top$Description));p<-ggplot(top,aes(x=Description,y=NES,size=setSize,color=NES))+geom_point()+coord_flip()+scale_color_gradient2(low="#b2182b",mid="grey90",high="#2166ac")+theme_minimal()+labs(x=NULL,y="Normalized enrichment score");ggsave(file.path(root,"dotplot.png"),p,width=7,height=5,dpi=150);ggsave(file.path(root,"dotplot.tiff"),p,width=7,height=5,dpi=300,compression="lzw")}
    contrast_summaries[[length(contrast_summaries)+1]]<-list(contrast_id=contrast$contrast_id,ranking=ranking,status=if(nrow(significant)==0)"NO_SIGNIFICANT_TERMS"else"SUCCESS",all_terms=nrow(terms),significant_terms=nrow(significant),positive_terms=nrow(positive),negative_terms=nrow(negative))
  }
  overall<-if(any(vapply(contrast_summaries,function(x)x$status=="NETWORK_UNAVAILABLE",logical(1))))"NETWORK_UNAVAILABLE"else if(length(blocked)>0)"BLOCKED"else if(all(vapply(contrast_summaries,function(x)x$status=="NO_SIGNIFICANT_TERMS",logical(1))))"NO_SIGNIFICANT_TERMS"else"SUCCESS";summary<-list(status=overall,reason=if(length(blocked)>0)paste(blocked,collapse=" ")else NULL,resource=resource,clusterProfiler_version=as.character(packageVersion("clusterProfiler")),r_version=R.version.string,contrasts=contrast_summaries);write(toJSON(summary,auto_unbox=TRUE,pretty=TRUE,null="null"),file.path(cfg$output_dir,"gsea_kegg_backend_summary.json"))
}
