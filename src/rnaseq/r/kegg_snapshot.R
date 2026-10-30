# KEGG gene-set snapshots: KEGG REST is read once into
# <root>/<snapshot_id>/{TERM2GENE.tsv,TERM2NAME.tsv,manifest.json}; every KEGG
# ORA/GSEA then reads a snapshot and never queries KEGG itself.  Snapshots live
# in the user cache, never in the repository.

nf_rna_sha256 <- function(path) paste(unclass(openssl::sha256(file(path))), collapse="")

nf_rna_snapshot_manifest <- function(directory) {
  path <- file.path(directory, "manifest.json")
  if (!file.exists(path)) return(NULL)
  tryCatch(jsonlite::fromJSON(path, simplifyVector=FALSE), error=function(e) NULL)
}

nf_rna_verify_snapshot <- function(directory) {
  manifest <- nf_rna_snapshot_manifest(directory)
  if (is.null(manifest)) stop("KEGG snapshot has no readable manifest.json: ", directory)
  for (name in c("TERM2GENE.tsv", "TERM2NAME.tsv")) {
    actual <- nf_rna_sha256(file.path(directory, name))
    if (!identical(actual, manifest$files[[name]]$sha256)) stop(sprintf("KEGG snapshot %s: SHA-256 of %s does not match manifest.json", basename(directory), name))
  }
  manifest
}

nf_rna_download_snapshot <- function(root, organism) {
  retrieved <- Sys.time()
  info <- KEGGREST::keggInfo(organism)
  pathway_info <- KEGGREST::keggInfo("pathway")
  links <- KEGGREST::keggLink("pathway", organism)
  names_raw <- KEGGREST::keggList("pathway", organism)
  if (!length(links) || !length(names_raw)) stop("KEGG REST returned no pathway records")
  term2gene <- data.frame(term=sub("^path:", "", unname(links)), gene=sub("^[a-z]+:", "", names(links)), stringsAsFactors=FALSE)
  term2gene <- unique(term2gene[order(term2gene$term, as.numeric(term2gene$gene), term2gene$gene), , drop=FALSE])
  term2name <- data.frame(term=sub("^path:", "", names(names_raw)), name=sub(" - [^-]+\\([^)]*\\)$", "", unname(names_raw)), stringsAsFactors=FALSE)
  term2name <- term2name[order(term2name$term), , drop=FALSE]
  # KEGG REST /info no longer prints a "Release" line; the pathway database
  # "Last update" date is the release information it does provide.
  info_lines <- trimws(strsplit(paste(c(pathway_info, info), collapse="\n"), "\n")[[1]])
  release_line <- head(c(grep("Release", info_lines, value=TRUE), grep("^Last update", info_lines, value=TRUE)), 1)
  release_tag <- if (length(release_line)) gsub("[^0-9.]", "", sub("^.*(Release|Last update) *", "", release_line)) else "unknown"
  id <- sprintf("%s_kegg%s_%s", organism, release_tag, format(retrieved, "%Y%m%dT%H%M%SZ", tz="UTC"))
  staging <- file.path(root, sprintf(".staging-%s-%d", id, Sys.getpid())); dir.create(staging, recursive=TRUE)
  on.exit(unlink(staging, recursive=TRUE), add=TRUE)
  write.table(term2gene, file.path(staging, "TERM2GENE.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
  write.table(term2name, file.path(staging, "TERM2NAME.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
  file_record <- function(name) list(sha256=nf_rna_sha256(file.path(staging, name)), bytes=file.size(file.path(staging, name)))
  manifest <- list(
    schema="nf-rna.kegg-snapshot.v1", snapshot_id=id, organism=organism,
    retrieved_at_utc=format(retrieved, "%Y-%m-%dT%H:%M:%SZ", tz="UTC"),
    source_urls=list(term2gene=sprintf("https://rest.kegg.jp/link/pathway/%s", organism), term2name=sprintf("https://rest.kegg.jp/list/pathway/%s", organism),
                     info=sprintf("https://rest.kegg.jp/info/%s", organism), release=sprintf("https://rest.kegg.jp/info/pathway")),
    retrieval=sprintf("KEGGREST %s (keggLink, keggList, keggInfo)", as.character(utils::packageVersion("KEGGREST"))),
    kegg_release=if (length(release_line)) sprintf("KEGG pathway database: %s", release_line) else NA_character_,
    kegg_info=paste(info, collapse="\n"), kegg_pathway_info=paste(pathway_info, collapse="\n"),
    transformation="pathway IDs without 'path:'; gene IDs without the organism prefix (NCBI Gene / Entrez); pathway names without the trailing ' - <organism>' suffix",
    files=list(TERM2GENE.tsv=file_record("TERM2GENE.tsv"), TERM2NAME.tsv=file_record("TERM2NAME.tsv")),
    pathways=length(unique(term2gene$term)), genes=length(unique(term2gene$gene)), pathway_gene_pairs=nrow(term2gene), named_pathways=nrow(term2name))
  write(jsonlite::toJSON(manifest, auto_unbox=TRUE, pretty=TRUE, null="null"), file.path(staging, "manifest.json"))
  if (!file.rename(staging, file.path(root, id))) stop("could not publish KEGG snapshot ", id)
  id
}

nf_rna_latest_snapshot <- function(root) {
  candidates <- list.dirs(root, recursive=FALSE, full.names=TRUE)
  candidates <- candidates[!startsWith(basename(candidates), ".")]
  manifests <- lapply(candidates, nf_rna_snapshot_manifest)
  keep <- !vapply(manifests, is.null, logical(1))
  if (!any(keep)) return(NULL)
  stamps <- vapply(manifests[keep], function(m) as.character(m$retrieved_at_utc %||% ""), character(1))
  ids <- basename(candidates[keep])
  ids[order(stamps, ids, decreasing=TRUE)][[1]]
}

# Resolve the snapshot of this run: the configured one, else the latest
# existing one, else a new download.  A lock directory serialises the
# KEGG ORA and KEGG GSEA tasks, which may start together.
nf_rna_resolve_kegg_snapshot <- function(root, organism, requested=NULL) {
  dir.create(root, recursive=TRUE, showWarnings=FALSE)
  created <- FALSE
  if (!is.null(requested) && nzchar(requested)) {
    id <- requested; selection <- "configured (annotation.enrichment.kegg.snapshot)"
    if (!dir.exists(file.path(root, id))) stop(sprintf("configured KEGG snapshot %s does not exist under %s", id, root))
  } else {
    lock <- file.path(root, ".lock"); waited <- 0
    while (!dir.create(lock, showWarnings=FALSE)) {
      if (waited >= 1800) stop("timed out waiting for the KEGG snapshot lock ", lock)
      Sys.sleep(2); waited <- waited + 2
    }
    on.exit(unlink(lock, recursive=TRUE), add=TRUE)
    id <- nf_rna_latest_snapshot(root); selection <- "latest existing snapshot"
    if (is.null(id)) { id <- nf_rna_download_snapshot(root, organism); selection <- "created for this run (no snapshot existed)"; created <- TRUE }
  }
  directory <- file.path(root, id)
  manifest <- nf_rna_verify_snapshot(directory)
  list(id=id, directory=directory, selection=selection, created=created, manifest=manifest,
       term2gene=read.delim(file.path(directory, "TERM2GENE.tsv"), colClasses="character"),
       term2name=read.delim(file.path(directory, "TERM2NAME.tsv"), colClasses="character", quote=""))
}

nf_rna_snapshot_record <- function(snapshot) list(
  provider="kegg_rest_snapshot", snapshot_id=snapshot$id, selection=snapshot$selection, created_by_this_run=snapshot$created,
  retrieved_at_utc=snapshot$manifest$retrieved_at_utc, kegg_release=snapshot$manifest$kegg_release, source_urls=snapshot$manifest$source_urls,
  sha256=list(TERM2GENE.tsv=snapshot$manifest$files$TERM2GENE.tsv$sha256, TERM2NAME.tsv=snapshot$manifest$files$TERM2NAME.tsv$sha256),
  sha256_verified=TRUE, pathways=snapshot$manifest$pathways, genes=snapshot$manifest$genes, reproducible=TRUE,
  enrichment_functions="clusterProfiler::enricher() / GSEA() with this TERM2GENE / TERM2NAME; no live KEGG query")
