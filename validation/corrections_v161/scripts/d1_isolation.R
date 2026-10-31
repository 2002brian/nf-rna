# D-1 isolation: the leading-edge correction alone, applied to the exact v1.6.0 GSEA.
#
# For one v1.6.0 run and one contrast, the GSEA of every collection is recomputed
# exactly as v1.6.0 did (set.seed(seed); seed = TRUE) from the run's own ranked
# gene list.  The script then checks that this reproduces the v1.6.0 all_terms.tsv,
# applies nf_rna_correct_leading_edge() and records which columns change, and
# compares the leading edge before and after with the Phase 2 independent fgsea
# reference (leadingEdge depends only on the ranked list and the gene set, so any
# reference seed serves; seed_1 is used).
#
# Usage: Rscript d1_isolation.R <repo> <v160_run_l2_dir> <contrast_id> <orgdb> <kegg_snapshot_dir> <reference_contrast_dir> <out.tsv>
args <- commandArgs(trailingOnly=TRUE)
repo <- args[[1]]; l2 <- args[[2]]; cid <- args[[3]]; orgdb_name <- args[[4]]; snapshot_dir <- args[[5]]; ref_dir <- args[[6]]; out <- args[[7]]
suppressPackageStartupMessages({ library(clusterProfiler); library(orgdb_name, character.only=TRUE) })
source(file.path(repo, "src/rnaseq/r/gsea_leading_edge.R"))
orgdb <- get(orgdb_name)

read_tsv <- function(path) read.delim(path, check.names=FALSE, stringsAsFactors=FALSE, colClasses=c(ID="character"))
ranked_list <- function(path) { r <- read.delim(path, colClasses=c(entrez_id="character")); setNames(r$stat, r$entrez_id) }
as_set <- function(x) sort(strsplit(x, "/", fixed=TRUE)[[1]])
same_set <- function(a, b) mapply(function(x, y) identical(as_set(x), as_set(y)), a, b)

collections <- list(
  GO_BP=list(dir="gsea_go", sub="BP"), GO_MF=list(dir="gsea_go", sub="MF"), GO_CC=list(dir="gsea_go", sub="CC"), KEGG=list(dir="gsea_kegg", sub=NULL))
term2gene <- read.delim(file.path(snapshot_dir, "TERM2GENE.tsv"), colClasses="character")
term2name <- read.delim(file.path(snapshot_dir, "TERM2NAME.tsv"), colClasses="character")
rows <- list()
for (name in names(collections)) {
  spec <- collections[[name]]
  root <- file.path(l2, "enrichment", spec$dir, cid)
  gene_list <- ranked_list(file.path(root, "ranked_gene_list.tsv"))
  v160 <- read_tsv(file.path(if (is.null(spec$sub)) root else file.path(root, spec$sub), "all_terms.tsv"))
  set.seed(1L)
  result <- if (is.null(spec$sub)) {
    suppressMessages(GSEA(geneList=gene_list, TERM2GENE=term2gene, TERM2NAME=term2name, minGSSize=10, maxGSSize=500, eps=0, pvalueCutoff=1, pAdjustMethod="BH", verbose=FALSE, seed=TRUE, nPermSimple=10000L))
  } else {
    suppressMessages(gseGO(geneList=gene_list, OrgDb=orgdb, keyType="ENTREZID", ont=spec$sub, minGSSize=10, maxGSSize=500, pvalueCutoff=1, pAdjustMethod="BH", eps=0, verbose=FALSE, seed=TRUE, nPermSimple=10000L))
  }
  before <- result@result
  at <- match(v160$ID, before$ID)
  numeric_cols <- c("enrichmentScore", "NES", "pvalue", "p.adjust", "qvalue")
  # write.table keeps 15 significant digits, so the comparison with the file is relative 1e-14.
  reproduced <- !anyNA(at) && nrow(before) == nrow(v160) &&
    all(vapply(numeric_cols, function(col) isTRUE(all.equal(before[[col]][at], v160[[col]], tolerance=1e-14)), logical(1))) &&
    identical(before$core_enrichment[at], v160$core_enrichment) && identical(as.integer(before$rank[at]), as.integer(v160$rank)) && identical(before$leading_edge[at], v160$leading_edge)
  corrected <- nf_rna_correct_leading_edge(result)
  after <- corrected$result@result
  stats_cols <- c("ID", "Description", "setSize", numeric_cols)
  statistics_unchanged <- identical(before[, stats_cols], after[, stats_cols])
  changed <- before$core_enrichment != after$core_enrichment
  ref <- read.delim(file.path(ref_dir, "gsea", name, "seed_1.tsv"), colClasses=c(set_id="character", leadingEdge="character"))
  ref_at <- match(after$ID, ref$set_id)
  common <- !is.na(ref_at)
  sig <- !is.na(after$p.adjust) & after$p.adjust <= 0.05
  neg <- after$enrichmentScore < 0
  before_ref <- same_set(before$core_enrichment[common], ref$leadingEdge[ref_at[common]])
  after_ref <- same_set(after$core_enrichment[common], ref$leadingEdge[ref_at[common]])
  rows[[name]] <- data.frame(contrast_id=cid, collection=name, v160_reproduced=reproduced, statistics_unchanged=statistics_unchanged,
    terms=nrow(after), negative_es_terms=sum(neg), core_changed=sum(changed), core_changed_positive_es=sum(changed & !neg), core_changed_negative_es=sum(changed & neg),
    core_changed_significant=sum(changed & sig), rank_changed=sum(before$rank != after$rank), leading_edge_tag_changed=sum(before$leading_edge != after$leading_edge),
    max_abs_es_difference=corrected$audit$max_abs_es_difference, terms_in_reference=sum(common),
    v160_equal_reference=sum(before_ref), corrected_equal_reference=sum(after_ref), corrected_differ_reference=sum(!after_ref), stringsAsFactors=FALSE)
  message(name, ": ", paste(names(rows[[name]]), unlist(rows[[name]]), sep="=", collapse=" "))
}
write.table(do.call(rbind, rows), out, sep="\t", quote=FALSE, row.names=FALSE)
