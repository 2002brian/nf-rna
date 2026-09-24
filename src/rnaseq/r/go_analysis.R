args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 2 || args[[1]] != "--config") stop("usage: go_analysis.R --config CONFIG.json")
script_arg <- commandArgs(trailingOnly=FALSE)
script_file <- sub("^--file=", "", script_arg[grep("^--file=", script_arg)][[1]])
source(file.path(dirname(normalizePath(script_file)), "provenance.R"))
source(file.path(dirname(normalizePath(script_file)), "ora_helpers.R"))
source(file.path(dirname(normalizePath(script_file)), "annotation_mapping_qc.R"))
suppressPackageStartupMessages({ library(jsonlite); library(AnnotationDbi); library(clusterProfiler); library(ggplot2) })
cfg <- fromJSON(args[[2]], simplifyVector=FALSE)
dir.create(cfg$output_dir, recursive=TRUE, showWarnings=FALSE)
orgdb <- get(cfg$orgdb_package, envir=asNamespace(cfg$orgdb_package))
input_type <- cfg$annotation$input_id_type
normalize_id <- function(x) if (input_type == "ENSEMBL") sub("\\.[0-9]+$", "", x) else x
write_table <- function(x, path) write.table(x, path, sep="\t", quote=FALSE, row.names=FALSE, na="NA")
empty_go <- function() data.frame(ID=character(), Description=character(), GeneRatio=character(), BgRatio=character(), pvalue=numeric(), p.adjust=numeric(), qvalue=numeric(), geneID=character(), Count=integer(), stringsAsFactors=FALSE)

map_sources <- function(source) {
  source <- unique(as.character(source[!is.na(source) & source != ""]))
  normalized <- normalize_id(source); valid <- AnnotationDbi::keys(orgdb, keytype=input_type)
  requested <- intersect(unique(normalized), valid)
  selected <- if (!length(requested)) data.frame() else suppressMessages(AnnotationDbi::select(orgdb, keys=requested, keytype=input_type, columns=c("ENTREZID", "SYMBOL")))
  if (!nrow(selected)) selected <- data.frame(key=character(), ENTREZID=character(), SYMBOL=character(), stringsAsFactors=FALSE)
  names(selected)[[1]] <- input_type
  selected <- selected[!is.na(selected$ENTREZID) & selected$ENTREZID != "", , drop=FALSE]
  rows <- lapply(seq_along(source), function(i) {
    hits <- selected[selected[[input_type]] == normalized[[i]], c("ENTREZID", "SYMBOL"), drop=FALSE]
    targets <- unique(hits$ENTREZID); status <- if (!length(targets)) "UNMAPPED" else if (length(targets)==1) "MAPPED_UNIQUE" else "ONE_TO_MANY"
    if (!nrow(hits)) data.frame(original_gene_id=source[[i]], normalized_gene_id=normalized[[i]], input_id_type=input_type, mapped_entrez_id=NA_character_, mapped_symbol=NA_character_, mapping_status=status, stringsAsFactors=FALSE) else data.frame(original_gene_id=source[[i]], normalized_gene_id=normalized[[i]], input_id_type=input_type, mapped_entrez_id=as.character(hits$ENTREZID), mapped_symbol=as.character(hits$SYMBOL), mapping_status=status, stringsAsFactors=FALSE)
  })
  if (!length(rows)) return(data.frame(original_gene_id=character(), normalized_gene_id=character(), input_id_type=character(), mapped_entrez_id=character(), mapped_symbol=character(), mapping_status=character(), stringsAsFactors=FALSE))
  do.call(rbind, rows)
}

mapping_stats <- function(mapping, source) {
  source <- unique(as.character(source[!is.na(source) & source != ""]))
  mapped <- unique(mapping$original_gene_id[!is.na(mapping$mapped_entrez_id)])
  targets <- unique(as.character(mapping$mapped_entrez_id[!is.na(mapping$mapped_entrez_id)]))
  list(source_genes=length(source), mapped_source_genes=length(mapped), unmapped_source_genes=length(source)-length(mapped), unique_target_genes=length(targets), one_to_many_source_ids=length(unique(mapping$original_gene_id[mapping$mapping_status=="ONE_TO_MANY"])), duplicate_target_ids=sum(table(mapping$mapped_entrez_id[!is.na(mapping$mapped_entrez_id)])>1), mapping_rate=if (!length(source)) 0 else length(mapped)/length(source), targets=targets)
}
without_targets <- function(x) { x$targets <- NULL; x }

member_table <- function(terms, mapping, path) {
  if (!nrow(terms)) { write_table(data.frame(GO_ID=character(), GO_term=character(), entrez_id=character(), symbol=character(), original_gene_id=character()), path); return() }
  rows <- list()
  for (i in seq_len(nrow(terms))) for (id in strsplit(terms$geneID[[i]], "/", fixed=TRUE)[[1]]) {
    related <- mapping[mapping$mapped_entrez_id == id, c("mapped_symbol", "original_gene_id"), drop=FALSE]
    if (!nrow(related)) related <- data.frame(mapped_symbol=NA_character_, original_gene_id=NA_character_)
    rows[[length(rows)+1]] <- data.frame(GO_ID=terms$ID[[i]], GO_term=terms$Description[[i]], entrez_id=id, symbol=related$mapped_symbol, original_gene_id=related$original_gene_id, stringsAsFactors=FALSE)
  }
  write_table(do.call(rbind, rows), path)
}

run_ontology <- function(targets, universe, ontology, directory, mapping) {
  all_path <- file.path(directory, paste0("GO_", ontology, "_all_terms.tsv")); out <- file.path(directory, paste0("GO_", ontology, ".tsv")); members <- file.path(directory, paste0("GO_", ontology, "_gene_members.tsv"))
  result <- suppressMessages(enrichGO(gene=targets, universe=universe, OrgDb=orgdb, keyType="ENTREZID", ont=ontology, pAdjustMethod=cfg$annotation$enrichment$go$p_adjust_method, pvalueCutoff=1, qvalueCutoff=1, readable=FALSE))
  terms <- if (is.null(result)) empty_go() else as.data.frame(result)
  if (nrow(terms)) terms <- terms[order(terms$p.adjust, terms$ID), , drop=FALSE]
  significant <- terms[!is.na(terms$pvalue) & terms$pvalue <= cfg$annotation$enrichment$go$pvalue_cutoff & !is.na(terms$qvalue) & terms$qvalue <= cfg$annotation$enrichment$go$qvalue_cutoff, , drop=FALSE]
  write_table(terms, all_path); write_table(significant, out); member_table(significant, mapping, members)
  if (nrow(terms) < nrow(significant)) stop("GO evaluated term count cannot be smaller than significant term count")
  if (!nrow(significant)) return(list(status="NO_SIGNIFICANT_TERMS", evaluated_term_count=nrow(terms), significant_term_count=0))
  top <- head(significant, 15); top$Description <- factor(top$Description, levels=rev(top$Description))
  p <- ggplot(top, aes(x=Description, y=-log10(p.adjust), size=Count)) + geom_point(color="#2166ac") + coord_flip() + labs(x=NULL, y="-log10(GO adjusted p-value)") + theme_minimal()
  ggsave(file.path(directory, paste0("dotplot_", ontology, ".png")), p, width=7, height=5, dpi=150)
  ggsave(file.path(directory, paste0("dotplot_", ontology, ".tiff")), p, width=7, height=5, dpi=300, compression="lzw")
  list(status="SUCCESS", evaluated_term_count=nrow(terms), significant_term_count=nrow(significant))
}

write_empty_ontology_artifacts <- function(directory, mapping) {
  for (ontology in c("BP", "MF", "CC")) {
    write_table(empty_go(), file.path(directory, paste0("GO_", ontology, "_all_terms.tsv")))
    write_table(empty_go(), file.path(directory, paste0("GO_", ontology, ".tsv")))
    member_table(empty_go(), mapping, file.path(directory, paste0("GO_", ontology, "_gene_members.tsv")))
  }
  setNames(lapply(c("BP", "MF", "CC"), function(x) list(status="NOT_APPLICABLE", evaluated_term_count=0, significant_term_count=0)), c("BP", "MF", "CC"))
}

contrast_summaries <- list(); all_mappings <- list()
for (contrast in cfg$contrasts) {
  root <- file.path(cfg$output_dir, contrast$contrast_id); dir.create(root, recursive=TRUE, showWarnings=FALSE)
  all_table <- read.delim(contrast$all_genes, check.names=FALSE, stringsAsFactors=FALSE)
  gene_sets <- nf_rna_ora_gene_sets(all_table)
  retained <- gene_sets$retained; tested <- gene_sets$tested
  mapping <- map_sources(retained); write_table(mapping, file.path(root, "gene_mapping.tsv")); all_mappings[[length(all_mappings)+1]] <- mapping
  tested_mapping <- mapping[mapping$original_gene_id %in% tested, , drop=FALSE]; universe_stats <- mapping_stats(tested_mapping, tested); universe <- universe_stats$targets
  # Shared dual-threshold mapping QC on the mapped tested-gene universe: only a
  # rate below minimum_mapping_rate blocks; a WARNING still runs ORA.
  mapping_qc <- annotation_mapping_qc(universe_stats$mapping_rate, cfg$annotation)
  universe_summary <- c(gene_sets$counts, list(successfully_mapped_tested_genes=length(universe)), without_targets(universe_stats), list(annotation_qc=mapping_qc))
  fg_summary <- list(); fgs <- list(significant=contrast$significant, up=contrast$up, down=contrast$down)
  for (name in names(fgs)) {
    foreground <- unique(as.character(read.delim(fgs[[name]], check.names=FALSE, stringsAsFactors=FALSE)$gene_id)); foreground <- foreground[!is.na(foreground) & foreground != ""]
    if (!all(foreground %in% tested)) stop(paste0("GO ORA foreground ", name, " is not a subset of statistically tested genes for contrast ", contrast$contrast_id))
    fg_mapping <- tested_mapping[tested_mapping$original_gene_id %in% foreground, , drop=FALSE]; fg_stats <- mapping_stats(fg_mapping, foreground); targets <- fg_stats$targets
    directory <- file.path(root, name); dir.create(directory, recursive=TRUE, showWarnings=FALSE)
    if (mapping_qc$status == "BLOCKED") outcome <- list(status="BLOCKED", reason=paste0("GO ORA blocked: ", mapping_qc$reason), evaluated_term_count=0, significant_term_count=0) else if (length(targets) < cfg$annotation$minimum_mapped_foreground) {
      outcome <- list(status="NOT_APPLICABLE", reason=paste0("only ", length(targets), " mapped foreground genes"), ontologies=write_empty_ontology_artifacts(directory, tested_mapping), evaluated_term_count=0, significant_term_count=0)
    } else {
      outcomes <- lapply(c("BP", "MF", "CC"), function(ontology) run_ontology(targets, universe, ontology, directory, tested_mapping)); names(outcomes) <- c("BP", "MF", "CC")
      evaluated_term_count <- sum(vapply(outcomes, function(x) x$evaluated_term_count, numeric(1)))
      significant_term_count <- sum(vapply(outcomes, function(x) x$significant_term_count, numeric(1)))
      if (evaluated_term_count < significant_term_count) stop("GO aggregate evaluated term count cannot be smaller than significant term count")
      outcome <- list(status=if (all(vapply(outcomes, function(x) x$status=="NO_SIGNIFICANT_TERMS", logical(1)))) "NO_SIGNIFICANT_TERMS" else "SUCCESS", ontologies=outcomes, evaluated_term_count=evaluated_term_count, significant_term_count=significant_term_count)
    }
    fg_summary[[name]] <- c(list(foreground_significant_genes=length(foreground), successfully_mapped_foreground_genes=length(targets)), without_targets(fg_stats), outcome)
  }
  contrast_summaries[[length(contrast_summaries)+1]] <- list(contrast_id=contrast$contrast_id, universe=universe_summary, foregrounds=fg_summary)
}
annotation_dir <- file.path(dirname(dirname(cfg$output_dir)), "annotation"); dir.create(annotation_dir, recursive=TRUE, showWarnings=FALSE)
if (length(all_mappings)) write_table(do.call(rbind, all_mappings), file.path(annotation_dir, "gene_mapping.tsv"))
statuses <- unlist(lapply(contrast_summaries, function(x) vapply(x$foregrounds, function(y) y$status, character(1))))
overall <- if (any(statuses == "FAILED")) "FAILED" else if (any(statuses == "BLOCKED")) "BLOCKED" else if (all(statuses == "NOT_APPLICABLE")) "NOT_APPLICABLE" else if (all(statuses != "SUCCESS")) "NO_SIGNIFICANT_TERMS" else "SUCCESS"
annotation_contract <- list(input_id_type=cfg$annotation$input_id_type, target_id_type=cfg$annotation$target_id_type, warning_threshold=as.numeric(cfg$annotation$mapping_warning_rate), blocking_threshold=as.numeric(cfg$annotation$minimum_mapping_rate))
summary <- list(status=overall, annotation_qc_status=annotation_qc_status(contrast_summaries), annotation=annotation_contract, annotation_database=cfg$orgdb_package, annotation_database_version=as.character(packageVersion(cfg$orgdb_package)), clusterProfiler_version=as.character(packageVersion("clusterProfiler")), contrasts=contrast_summaries)
write(toJSON(summary, auto_unbox=TRUE, pretty=TRUE, null="null"), file.path(cfg$output_dir, "go_backend_summary.json"))
nf_rna_write_provenance(cfg, cfg$output_dir, "GO_ORA", overall, c("AnnotationDbi", "clusterProfiler", "ggplot2", "jsonlite", cfg$orgdb_package), summary)
