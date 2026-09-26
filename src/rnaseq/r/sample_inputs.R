# Shared L1/L2 resolution of per-sample Salmon quant.sf files.  The frozen
# config maps sample_id -> quant.sf path; its key order is not the metadata
# sample order (the control plane writes JSON with sorted keys), so files are
# always selected by sample name and never by position.

nf_rna_quant_sf_files <- function(quant_sf, samples) {
  if (length(samples) == 0 || anyNA(samples) || any(samples == "")) stop("metadata samples must be non-empty sample identifiers")
  if (anyDuplicated(samples)) stop(paste("metadata samples are not unique:", paste(unique(samples[duplicated(samples)]), collapse = ", ")))
  keys <- names(quant_sf)
  if (length(quant_sf) == 0 || is.null(keys)) stop("quant_sf must map each sample_id to its quant.sf path")
  if (anyNA(keys) || any(keys == "")) stop("quant_sf contains an entry without a sample_id")
  if (anyDuplicated(keys)) stop(paste("quant_sf lists a sample_id more than once:", paste(unique(keys[duplicated(keys)]), collapse = ", ")))
  missing <- setdiff(samples, keys)
  if (length(missing) > 0) stop(paste("quant_sf is missing metadata samples:", paste(missing, collapse = ", ")))
  unexpected <- setdiff(keys, samples)
  if (length(unexpected) > 0) stop(paste("quant_sf lists samples absent from metadata:", paste(unexpected, collapse = ", ")))
  for (sample in samples) {
    path <- quant_sf[[sample]]
    if (!is.character(path) || length(path) != 1 || is.na(path) || path == "") stop(paste("quant_sf path is missing or invalid for sample:", sample))
  }
  files <- vapply(samples, function(sample) quant_sf[[sample]], character(1), USE.NAMES = FALSE)
  names(files) <- samples
  files
}

nf_rna_tximport_counts <- function(txi, samples) {
  if (!identical(colnames(txi$counts), samples)) stop("tximport sample columns do not match metadata samples")
  txi$counts
}
