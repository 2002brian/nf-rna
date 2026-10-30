# Shared OrgDb identifier rules for L2 DEG annotation and the GO/KEGG ORA and
# GSEA mappings.  An ENSEMBL key drops its trailing version suffix for the
# OrgDb query only; the delivered gene_id is never rewritten.
nf_rna_annotation_key <- function(x, input_type) if (identical(input_type, "ENSEMBL")) sub("\\.[0-9]+$", "", x) else x

# Used only when a project has no annotation block: every gene_id must match
# one unambiguous namespace, otherwise the DEG tables stay unannotated.
nf_rna_infer_id_type <- function(ids) {
  if (length(ids) && all(grepl("^ENS[A-Z]*G[0-9]+(\\.[0-9]+)?$", ids))) "ENSEMBL"
  else if (length(ids) && all(grepl("^[0-9]+$", ids))) "ENTREZID"
  else NA_character_
}

# One row per gene_id, in input order.  symbol and description use
# mapIds(multiVals = "first"); a gene_id without a mapping keeps NA.
nf_rna_gene_annotation <- function(gene_ids, settings) {
  table <- data.frame(gene_id = gene_ids, symbol = NA_character_, description = NA_character_, stringsAsFactors = FALSE)
  package <- settings$orgdb_package
  input_type <- settings$input_id_type
  id_source <- "annotation.input_id_type"
  if (is.null(input_type) && !is.null(package)) {
    input_type <- nf_rna_infer_id_type(gene_ids)
    id_source <- "inferred from the gene_id pattern (project has no annotation block)"
  }
  info <- list(status = "NOT_AVAILABLE", orgdb_package = package, input_id_type = input_type, input_id_type_source = id_source,
               version_suffix = "ENSEMBL gene_id version suffixes are removed for the OrgDb query only; gene_id is unchanged",
               one_to_many_rule = "mapIds(multiVals = \"first\"); unmapped values are NA; one row per gene_id")
  if (is.null(package)) return(list(table = table, info = c(info, list(reason = "organism.species has no supported OrgDb"))))
  if (is.null(input_type) || is.na(input_type)) return(list(table = table, info = c(info, list(reason = "gene_id namespace is not declared and could not be inferred"))))
  suppressPackageStartupMessages(library(AnnotationDbi))
  orgdb <- get(package, envir = asNamespace(package))
  keys <- nf_rna_annotation_key(gene_ids, input_type)
  query <- intersect(unique(keys), AnnotationDbi::keys(orgdb, keytype = input_type))
  one_to_many <- 0L
  if (length(query)) {
    symbol <- suppressMessages(AnnotationDbi::mapIds(orgdb, keys = query, column = "SYMBOL", keytype = input_type, multiVals = "first"))
    description <- suppressMessages(AnnotationDbi::mapIds(orgdb, keys = query, column = "GENENAME", keytype = input_type, multiVals = "first"))
    all_symbols <- suppressMessages(AnnotationDbi::mapIds(orgdb, keys = query, column = "SYMBOL", keytype = input_type, multiVals = "list"))
    table$symbol <- unname(symbol[keys]); table$description <- unname(description[keys])
    one_to_many <- sum(lengths(lapply(all_symbols, function(v) unique(v[!is.na(v)])))[keys[keys %in% query]] > 1)
  }
  n <- length(gene_ids)
  info$status <- "SUCCESS"
  info <- c(info, list(orgdb_version = as.character(utils::packageVersion(package)), genes = n,
                       symbol_mapped = sum(!is.na(table$symbol)), description_mapped = sum(!is.na(table$description)),
                       symbol_mapping_rate = if (n) sum(!is.na(table$symbol)) / n else 0,
                       description_mapping_rate = if (n) sum(!is.na(table$description)) / n else 0,
                       one_to_many_symbol_gene_ids = one_to_many))
  list(table = table, info = info)
}
