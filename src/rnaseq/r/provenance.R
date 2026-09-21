# Shared production provenance for every R backend.  Values derived from the
# task filesystem are deliberately omitted from the normalized configuration.

nf_rna_package_versions <- function(packages) {
  values <- lapply(packages, function(package) {
    if (requireNamespace(package, quietly=TRUE)) as.character(utils::packageVersion(package)) else NA_character_
  })
  names(values) <- packages
  c(list(R=R.version.string), values)
}

nf_rna_normalized_config <- function(value, key=NULL) {
  path_keys <- c("output_dir", "counts", "metadata", "tx2gene", "quant_sf", "l1_vst", "all_genes", "significant", "up", "down")
  if (!is.null(key) && key %in% path_keys) return(NULL)
  if (is.list(value)) {
    keys <- names(value)
    if (is.null(keys)) {
      result <- lapply(value, nf_rna_normalized_config)
      return(result[!vapply(result, is.null, logical(1))])
    }
    result <- lapply(keys, function(name) nf_rna_normalized_config(value[[name]], name))
    names(result) <- keys
    return(result[!vapply(result, is.null, logical(1))])
  }
  value
}

nf_rna_contrast_identity <- function(contrasts) {
  lapply(contrasts %||% list(), function(contrast) {
    if (!is.list(contrast)) return(NULL)
    keep <- intersect(c("contrast_id", "factor", "numerator", "denominator"), names(contrast))
    contrast[keep]
  })
}

nf_rna_write_provenance <- function(cfg, output_dir, module, result_state, packages, extra=list(), design_details=list()) {
  dir.create(output_dir, recursive=TRUE, showWarnings=FALSE)
  session_path <- file.path(output_dir, "r_session_info.txt")
  writeLines(capture.output(sessionInfo()), session_path, useBytes=TRUE)
  runtime <- if (is.list(cfg$runtime)) cfg$runtime else list()
  execution_state <- if (identical(result_state, "SUCCESS")) "COMPLETED" else result_state
  document <- list(
    schema_version="nf-rna.scientific-provenance.v1",
    module=module,
    execution_state=execution_state,
    result_state=result_state,
    generated_at_utc=format(Sys.time(), tz="UTC", usetz=TRUE),
    script_identity=list(name=runtime$script %||% NA_character_, source_revision=runtime$source_revision %||% NA_character_),
    container_image=runtime$container_image %||% NA_character_,
    r_version=R.version.string,
    package_versions=nf_rna_package_versions(packages),
    organism=if (is.list(cfg$annotation)) list(organism=cfg$annotation$organism, input_id_type=cfg$annotation$input_id_type, target_id_type=cfg$annotation$target_id_type) else NULL,
    design=c(
      list(formula=cfg$formula %||% NULL, contrasts=nf_rna_contrast_identity(cfg$contrasts)),
      design_details
    ),
    normalized_configuration=nf_rna_normalized_config(cfg),
    output_schema=runtime$output_schema %||% "nf-rna.scientific-provenance.v1",
    result=extra
  )
  write(jsonlite::toJSON(document, auto_unbox=TRUE, pretty=TRUE, null="null", na="null"), file.path(output_dir, "scientific_provenance.json"))
}

`%||%` <- function(value, fallback) if (is.null(value) || length(value) == 0) fallback else value
