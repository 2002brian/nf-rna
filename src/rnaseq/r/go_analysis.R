args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 2 || args[[1]] != "--config") stop("usage: go_analysis.R --config CONFIG.json")
script_arg <- commandArgs(trailingOnly=FALSE)
script_file <- sub("^--file=", "", script_arg[grep("^--file=", script_arg)][[1]])
source(file.path(dirname(normalizePath(script_file)), "provenance.R"))
source(file.path(dirname(normalizePath(script_file)), "gene_annotation.R"))
source(file.path(dirname(normalizePath(script_file)), "ora_helpers.R"))
source(file.path(dirname(normalizePath(script_file)), "annotation_mapping_qc.R"))
source(file.path(dirname(normalizePath(script_file)), "figures.R"))
source(file.path(dirname(normalizePath(script_file)), "enrichment_common.R"))
suppressPackageStartupMessages({ library(jsonlite); library(AnnotationDbi); library(clusterProfiler); library(ggplot2) })
cfg <- fromJSON(args[[2]], simplifyVector=FALSE)
dir.create(cfg$output_dir, recursive=TRUE, showWarnings=FALSE)
fig <- nf_rna_figure_context(cfg)
orgdb <- get(cfg$orgdb_package, envir=asNamespace(cfg$orgdb_package))
input_type <- cfg$annotation$input_id_type
normalize_id <- function(x) nf_rna_annotation_key(x, input_type)
write_table <- function(x, path) write.table(x, path, sep="\t", quote=FALSE, row.names=FALSE, na="NA")
empty_go <- function() data.frame(ID=character(), Description=character(), GeneRatio=character(), BgRatio=character(), pvalue=numeric(), p.adjust=numeric(), qvalue=numeric(), geneID=character(), Count=integer(), significant=logical(), stringsAsFactors=FALSE)
# clusterProfiler's enrichGO defaults, now passed explicitly and recorded.
GO_MIN_GS_SIZE <- 10L; GO_MAX_GS_SIZE <- 500L
# Entrez IDs with at least one GO annotation in each ontology: the gene-set
# database the ORA universe is intersected with.
go_database_genes <- setNames(lapply(c("BP", "MF", "CC"), function(ontology) {
  unique(as.character(suppressMessages(AnnotationDbi::select(orgdb, keys=ontology, keytype="ONTOLOGY", columns="ENTREZID"))$ENTREZID))
}), c("BP", "MF", "CC"))

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
    # which() drops the NA comparisons of unmapped genes; logical indexing would
    # emit one all-NA row per unmapped gene for every term-gene pair.
    related <- mapping[which(mapping$mapped_entrez_id == id), c("mapped_symbol", "original_gene_id"), drop=FALSE]
    if (!nrow(related)) related <- data.frame(mapped_symbol=NA_character_, original_gene_id=NA_character_)
    rows[[length(rows)+1]] <- data.frame(GO_ID=terms$ID[[i]], GO_term=terms$Description[[i]], entrez_id=id, symbol=related$mapped_symbol, original_gene_id=related$original_gene_id, stringsAsFactors=FALSE)
  }
  write_table(do.call(rbind, rows), path)
}

run_ontology <- function(targets, universe, ontology, directory, mapping, contrast, foreground, lookup) {
  all_path <- file.path(directory, paste0("GO_", ontology, "_all_terms.tsv")); out <- file.path(directory, paste0("GO_", ontology, ".tsv")); members <- file.path(directory, paste0("GO_", ontology, "_gene_members.tsv"))
  result <- suppressMessages(enrichGO(gene=targets, universe=universe$genes, OrgDb=orgdb, keyType="ENTREZID", ont=ontology, pAdjustMethod=NF_RNA_ENRICHMENT_ADJUST, pvalueCutoff=1, qvalueCutoff=1, minGSSize=GO_MIN_GS_SIZE, maxGSSize=GO_MAX_GS_SIZE, readable=FALSE))
  terms <- if (is.null(result)) empty_go() else nf_rna_flag_significant(as.data.frame(result))
  if (nrow(terms)) terms <- terms[order(terms$p.adjust, terms$ID), , drop=FALSE]
  significant <- terms[terms$significant, , drop=FALSE]
  write_table(terms, all_path); write_table(significant, out); member_table(significant, mapping, members)
  member_record <- nf_rna_write_members(significant, nf_rna_ora_member_pairs(significant, mapping), lookup, file.path(directory, paste0("GO_", ontology, "_members.tsv")))
  universe_record <- c(universe$record, list(foreground_entrez_ids=length(targets), foreground_in_universe=length(intersect(targets, universe$genes)),
                                           clusterprofiler_bgratio_n=nf_rna_bgratio_n(terms)))
  if (nrow(terms) < nrow(significant)) stop("GO evaluated term count cannot be smaller than significant term count")
  base <- list(evaluated_term_count=nrow(terms), significant_term_count=nrow(significant), universe=universe_record, members=member_record,
               min_gs_size=GO_MIN_GS_SIZE, max_gs_size=GO_MAX_GS_SIZE)
  if (!nrow(significant)) return(c(list(status="NO_SIGNIFICANT_TERMS"), base, list(figures=list(dotplot="not drawn: no significant terms (p.adjust <= 0.05)", dotplot_simplified="not drawn: no significant terms (p.adjust <= 0.05)"),
                                                                                    simplified=list(significant_terms=0, simplified_terms=0, reason="no significant terms"))))
  simplified <- nf_rna_simplify_go(result, significant)
  simplified_path <- file.path(directory, paste0("GO_", ontology, "_simplified.tsv")); write_table(simplified$table, simplified_path)
  dotplot_meta <- function(table, path, what) list(ontology=ontology, contrast=list(contrast_id=contrast$contrast_id, numerator=contrast$numerator, denominator=contrast$denominator), foreground=foreground,
    display_subset=sprintf("top %d %s by adjusted P, ties by ID, of %d (p.adjust <= %s)", NF_RNA_DOTPLOT_TOP, what, nrow(table), NF_RNA_ENRICHMENT_PADJ),
    significance_semantics=NF_RNA_SIGNIFICANCE_SEMANTICS, source_tables=list(nf_rna_source_data(path, cfg$output_dir, what)))
  x_label <- sprintf("-log10(GO %s adjusted P, %s)", ontology, NF_RNA_ENRICHMENT_ADJUST)
  nf_rna_ora_dotplot(fig, head(significant, NF_RNA_DOTPLOT_TOP), x_label, file.path(directory, paste0("dotplot_", ontology)), cfg$output_dir,
    dotplot_meta(significant, out, "significant terms"), subtitle=nf_rna_ora_subtitle(fig, contrast, foreground), formats=NF_RNA_ENRICHMENT_FORMATS)
  nf_rna_ora_dotplot(fig, head(simplified$table, NF_RNA_DOTPLOT_TOP), x_label, file.path(directory, paste0("dotplot_", ontology, "_simplified")), cfg$output_dir,
    dotplot_meta(simplified$table, simplified_path, "simplified significant terms"), subtitle=nf_rna_ora_subtitle(fig, contrast, foreground), formats=NF_RNA_ENRICHMENT_FORMATS)
  c(list(status="SUCCESS"), base, list(simplified=nf_rna_simplify_record(significant, simplified), figures=list(dotplot="drawn", dotplot_simplified="drawn")))
}

write_empty_ontology_artifacts <- function(directory, mapping) {
  for (ontology in c("BP", "MF", "CC")) {
    write_table(empty_go(), file.path(directory, paste0("GO_", ontology, "_all_terms.tsv")))
    write_table(empty_go(), file.path(directory, paste0("GO_", ontology, ".tsv")))
    member_table(empty_go(), mapping, file.path(directory, paste0("GO_", ontology, "_gene_members.tsv")))
    nf_rna_write_members(empty_go(), nf_rna_ora_member_pairs(empty_go(), mapping), NULL, file.path(directory, paste0("GO_", ontology, "_members.tsv")))
  }
  setNames(lapply(c("BP", "MF", "CC"), function(x) list(status="NOT_APPLICABLE", evaluated_term_count=0, significant_term_count=0, figures=list(dotplot="not drawn: too few mapped foreground genes"))), c("BP", "MF", "CC"))
}

contrast_summaries <- list(); all_mappings <- list()
for (contrast in cfg$contrasts) {
  root <- file.path(cfg$output_dir, contrast$contrast_id); dir.create(root, recursive=TRUE, showWarnings=FALSE)
  all_table <- read.delim(contrast$all_genes, check.names=FALSE, stringsAsFactors=FALSE)
  lookup <- nf_rna_gene_lookup(all_table)
  gene_sets <- nf_rna_ora_gene_sets(all_table)
  retained <- gene_sets$retained; tested <- gene_sets$tested
  mapping <- map_sources(retained); write_table(mapping, file.path(root, "gene_mapping.tsv")); all_mappings[[length(all_mappings)+1]] <- mapping
  tested_mapping <- mapping[mapping$original_gene_id %in% tested, , drop=FALSE]; universe_stats <- mapping_stats(tested_mapping, tested); universe <- universe_stats$targets
  universes <- lapply(go_database_genes, function(genes) nf_rna_ora_universe(tested, tested_mapping, genes))
  # Shared dual-threshold mapping QC on the mapped tested-gene universe: only a
  # rate below minimum_mapping_rate blocks; a WARNING still runs ORA.
  mapping_qc <- annotation_mapping_qc(universe_stats$mapping_rate, cfg$annotation)
  universe_summary <- c(gene_sets$counts, list(successfully_mapped_tested_genes=length(universe)), without_targets(universe_stats), list(annotation_qc=mapping_qc),
                        list(by_ontology=lapply(universes, function(u) u$record)))
  fg_summary <- list(); fgs <- list(significant=contrast$significant, up=contrast$up, down=contrast$down)
  for (name in names(fgs)) {
    foreground <- unique(as.character(read.delim(fgs[[name]], check.names=FALSE, stringsAsFactors=FALSE)$gene_id)); foreground <- foreground[!is.na(foreground) & foreground != ""]
    if (!all(foreground %in% tested)) stop(paste0("GO ORA foreground ", name, " is not a subset of the padj-non-NA genes for contrast ", contrast$contrast_id))
    fg_mapping <- tested_mapping[tested_mapping$original_gene_id %in% foreground, , drop=FALSE]; fg_stats <- mapping_stats(fg_mapping, foreground); targets <- fg_stats$targets
    directory <- file.path(root, name); dir.create(directory, recursive=TRUE, showWarnings=FALSE)
    if (mapping_qc$status == "BLOCKED") outcome <- list(status="BLOCKED", reason=paste0("GO ORA blocked: ", mapping_qc$reason), evaluated_term_count=0, significant_term_count=0) else if (length(targets) < cfg$annotation$minimum_mapped_foreground) {
      outcome <- list(status="NOT_APPLICABLE", reason=paste0("only ", length(targets), " mapped foreground genes"), ontologies=write_empty_ontology_artifacts(directory, tested_mapping), evaluated_term_count=0, significant_term_count=0)
    } else {
      outcomes <- lapply(c("BP", "MF", "CC"), function(ontology) run_ontology(targets, universes[[ontology]], ontology, directory, tested_mapping[tested_mapping$original_gene_id %in% foreground, , drop=FALSE], contrast, name, lookup)); names(outcomes) <- c("BP", "MF", "CC")
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
summary <- list(status=overall, significance_semantics=NF_RNA_SIGNIFICANCE_SEMANTICS, significance=nf_rna_significance_record(),
                ora=list(foregrounds="significant = significant.tsv, up = upregulated.tsv (log2FoldChange > 0), down = downregulated.tsv (log2FoldChange < 0); the unchanged L2 DEG definition",
                         universe="padj-non-NA genes mapped to Entrez, intersected with the Entrez IDs annotated in the ontology", min_gs_size=GO_MIN_GS_SIZE, max_gs_size=GO_MAX_GS_SIZE),
                gene_set_source=list(kind="GO via OrgDb", orgdb=cfg$orgdb_package, orgdb_version=as.character(packageVersion(cfg$orgdb_package)), go_db_version=as.character(packageVersion("GO.db")),
                                     clusterProfiler_version=as.character(packageVersion("clusterProfiler")), snapshot="none; GO results depend on these package versions"),
                annotation_qc_status=annotation_qc_status(contrast_summaries), annotation=annotation_contract, annotation_database=cfg$orgdb_package, annotation_database_version=as.character(packageVersion(cfg$orgdb_package)), clusterProfiler_version=as.character(packageVersion("clusterProfiler")), contrasts=contrast_summaries)
write(toJSON(summary, auto_unbox=TRUE, pretty=TRUE, null="null"), file.path(cfg$output_dir, "go_backend_summary.json"))
nf_rna_write_figure_manifest(fig, cfg$output_dir)
nf_rna_write_provenance(cfg, cfg$output_dir, "GO_ORA", overall, c("AnnotationDbi", "clusterProfiler", "GO.db", "GOSemSim", "DOSE", "ggplot2", "scales", "systemfonts", "jsonlite", cfg$orgdb_package), summary)
