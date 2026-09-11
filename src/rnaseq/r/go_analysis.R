args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2 || args[[1]] != "--config") stop("usage: go_analysis.R --config CONFIG.json")
suppressPackageStartupMessages({ library(jsonlite); library(AnnotationDbi); library(clusterProfiler); library(ggplot2) })
cfg <- fromJSON(args[[2]], simplifyVector = FALSE)
dir.create(cfg$output_dir, recursive = TRUE, showWarnings = FALSE)
annotation_dir <- file.path(dirname(dirname(cfg$output_dir)), "annotation"); dir.create(annotation_dir, recursive = TRUE, showWarnings = FALSE)
orgdb <- get(cfg$orgdb_package, envir = asNamespace(cfg$orgdb_package))
input_type <- cfg$annotation$input_id_type
normalize_id <- function(x) if (input_type == "ENSEMBL") sub("\\.[0-9]+$", "", x) else x
all_genes <- read.delim(cfg$contrasts[[1]]$all_genes, check.names = FALSE, stringsAsFactors = FALSE)$gene_id
source_ids <- unique(as.character(all_genes)); normalized <- normalize_id(source_ids)
valid_keys <- AnnotationDbi::keys(orgdb, keytype = input_type)
requested_keys <- intersect(unique(normalized), valid_keys)
selected <- if (length(requested_keys) == 0) data.frame() else suppressMessages(AnnotationDbi::select(orgdb, keys = requested_keys, keytype = input_type, columns = c("ENTREZID", "SYMBOL")))
if (nrow(selected) == 0) selected <- data.frame(key = character(), ENTREZID = character(), SYMBOL = character(), stringsAsFactors = FALSE)
names(selected)[[1]] <- input_type
selected <- selected[!is.na(selected$ENTREZID) & selected$ENTREZID != "", , drop = FALSE]
mapping_rows <- lapply(seq_along(source_ids), function(i) {
  hits <- selected[selected[[input_type]] == normalized[[i]], c("ENTREZID", "SYMBOL"), drop = FALSE]
  targets <- unique(hits$ENTREZID)
  status <- if (length(targets) == 0) "UNMAPPED" else if (length(targets) == 1) "MAPPED_UNIQUE" else "ONE_TO_MANY"
  if (nrow(hits) == 0) data.frame(original_gene_id = source_ids[[i]], normalized_gene_id = normalized[[i]], input_id_type = input_type, mapped_entrez_id = NA_character_, mapped_symbol = NA_character_, mapping_status = status, stringsAsFactors = FALSE) else data.frame(original_gene_id = source_ids[[i]], normalized_gene_id = normalized[[i]], input_id_type = input_type, mapped_entrez_id = hits$ENTREZID, mapped_symbol = hits$SYMBOL, mapping_status = status, stringsAsFactors = FALSE)
})
mapping <- do.call(rbind, mapping_rows)
target_source_count <- table(mapping$mapped_entrez_id[!is.na(mapping$mapped_entrez_id)])
mapping$target_duplicate <- !is.na(mapping$mapped_entrez_id) & target_source_count[mapping$mapped_entrez_id] > 1
write.table(mapping, file.path(annotation_dir, "gene_mapping.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
stats_for <- function(ids) {
  ids <- unique(as.character(ids)); submap <- mapping[mapping$original_gene_id %in% ids, , drop = FALSE]
  mapped_sources <- unique(submap$original_gene_id[!is.na(submap$mapped_entrez_id)])
  targets <- unique(submap$mapped_entrez_id[!is.na(submap$mapped_entrez_id)])
  list(source_genes = length(ids), mapped_source_genes = length(mapped_sources), unmapped_source_genes = length(ids) - length(mapped_sources), one_to_many_source_ids = length(unique(submap$original_gene_id[submap$mapping_status == "ONE_TO_MANY"])), duplicate_target_ids = sum(table(submap$mapped_entrez_id[!is.na(submap$mapped_entrez_id)]) > 1), unique_target_genes = length(targets), mapping_rate = if (length(ids) == 0) 0 else length(mapped_sources) / length(ids), targets = as.list(targets))
}
tested_stats <- stats_for(source_ids)
strip_targets <- function(value) { value$targets <- NULL; value }
mapping_summary <- list(tested = strip_targets(tested_stats))
if (tested_stats$mapping_rate < cfg$annotation$mapping_warning_rate) {
  reason <- paste0("GO enrichment blocked. Mapped ", sprintf("%.1f%%", 100 * tested_stats$mapping_rate), " of tested genes. Required minimum: ", sprintf("%.1f%%", 100 * cfg$annotation$mapping_warning_rate), ". Review annotation organism and input identifier type.")
  summary <- list(status = "BLOCKED", reason = reason, mapping = mapping_summary, annotation_database = cfg$orgdb_package, annotation_database_version = as.character(packageVersion(cfg$orgdb_package)), clusterProfiler_version = as.character(packageVersion("clusterProfiler")), contrasts = list())
  write(toJSON(summary, auto_unbox = TRUE, pretty = TRUE), file.path(cfg$output_dir, "go_backend_summary.json")); quit(status = 0)
}
empty_go <- function(path) write.table(data.frame(ID=character(), Description=character(), GeneRatio=character(), BgRatio=character(), pvalue=numeric(), p.adjust=numeric(), qvalue=numeric(), geneID=character(), Count=integer()), path, sep="\t", quote=FALSE, row.names=FALSE)
member_table <- function(result, path) {
  if (is.null(result) || nrow(as.data.frame(result)) == 0) { write.table(data.frame(GO_ID=character(), GO_term=character(), entrez_id=character(), symbol=character(), original_gene_id=character()), path, sep="\t", quote=FALSE, row.names=FALSE); return() }
  terms <- as.data.frame(result); rows <- list()
  for (i in seq_len(nrow(terms))) for (id in strsplit(terms$geneID[[i]], "/", fixed=TRUE)[[1]]) {
    related <- mapping[mapping$mapped_entrez_id == id, c("mapped_symbol", "original_gene_id"), drop=FALSE]
    if (nrow(related) == 0) related <- data.frame(mapped_symbol=NA_character_, original_gene_id=NA_character_)
    rows[[length(rows)+1]] <- data.frame(GO_ID=terms$ID[[i]], GO_term=terms$Description[[i]], entrez_id=id, symbol=related$mapped_symbol, original_gene_id=related$original_gene_id, stringsAsFactors=FALSE)
  }
  write.table(do.call(rbind, rows), path, sep="\t", quote=FALSE, row.names=FALSE, na="NA")
}
run_ontology <- function(targets, ontology, directory) {
  out <- file.path(directory, paste0("GO_", ontology, ".tsv")); members <- file.path(directory, paste0("GO_", ontology, "_gene_members.tsv"))
  result <- suppressMessages(enrichGO(gene=targets, universe=unlist(tested_stats$targets), OrgDb=orgdb, keyType="ENTREZID", ont=ontology, pAdjustMethod=cfg$annotation$enrichment$go$p_adjust_method, pvalueCutoff=cfg$annotation$enrichment$go$pvalue_cutoff, qvalueCutoff=cfg$annotation$enrichment$go$qvalue_cutoff, readable=FALSE))
  table <- if (is.null(result)) data.frame() else as.data.frame(result)
  if (nrow(table) == 0) { empty_go(out); member_table(NULL, members); return(list(status="NO_SIGNIFICANT_TERMS", terms=0)) }
  table <- table[order(table$p.adjust, table$ID), , drop=FALSE]
  write.table(table, out, sep="\t", quote=FALSE, row.names=FALSE, na="NA"); member_table(result, members)
  top <- head(table, 15); top$Description <- factor(top$Description, levels=rev(top$Description))
  p <- ggplot(top, aes(x=Description, y=-log10(p.adjust), size=Count)) + geom_point(color="#2166ac") + coord_flip() + labs(x=NULL, y="-log10(GO adjusted p-value)") + theme_minimal()
  ggsave(file.path(directory, paste0("dotplot_", ontology, ".png")), p, width=7, height=5, dpi=150)
  ggsave(file.path(directory, paste0("dotplot_", ontology, ".tiff")), p, width=7, height=5, dpi=300, compression="lzw")
  list(status="SUCCESS", terms=nrow(table))
}
contrast_summaries <- list()
for (contrast in cfg$contrasts) {
  result_root <- file.path(cfg$output_dir, contrast$contrast_id); dir.create(result_root, recursive=TRUE, showWarnings=FALSE)
  fgs <- list(significant=contrast$significant, up=contrast$up, down=contrast$down); fg_summary <- list()
  for (name in names(fgs)) {
    ids <- read.delim(fgs[[name]], check.names=FALSE, stringsAsFactors=FALSE)$gene_id
    st <- stats_for(ids); targets <- unique(unlist(st$targets)); directory <- file.path(result_root, name); dir.create(directory, recursive=TRUE, showWarnings=FALSE)
    st$targets <- NULL
    if (length(targets) < cfg$annotation$minimum_mapped_foreground) {
      for (ontology in c("BP","MF","CC")) { empty_go(file.path(directory, paste0("GO_",ontology,".tsv"))); member_table(NULL, file.path(directory, paste0("GO_",ontology,"_gene_members.tsv"))) }
      fg_summary[[name]] <- c(st, list(status="NOT_APPLICABLE", reason=paste0("only ", length(targets), " mapped foreground genes"), term_counts=list(BP=0, MF=0, CC=0)))
    } else {
      statuses <- lapply(c("BP","MF","CC"), function(x) run_ontology(targets, x, directory)); names(statuses) <- c("BP","MF","CC")
      nterms <- lapply(statuses, function(x) x$terms); names(nterms) <- names(statuses)
      fg_status <- if (all(vapply(statuses, function(x) x$status == "NO_SIGNIFICANT_TERMS", logical(1)))) "NO_SIGNIFICANT_TERMS" else "SUCCESS"
      fg_summary[[name]] <- c(st, list(status=fg_status, term_counts=nterms))
    }
  }
  contrast_summaries[[length(contrast_summaries)+1]] <- list(contrast_id=contrast$contrast_id, foregrounds=fg_summary)
}
summary <- list(status="SUCCESS", mapping=mapping_summary, annotation_database=cfg$orgdb_package, annotation_database_version=as.character(packageVersion(cfg$orgdb_package)), clusterProfiler_version=as.character(packageVersion("clusterProfiler")), contrasts=contrast_summaries)
write(toJSON(summary, auto_unbox=TRUE, pretty=TRUE), file.path(cfg$output_dir, "go_backend_summary.json"))
