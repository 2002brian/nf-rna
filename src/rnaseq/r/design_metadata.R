# Shared L1/L2 design-metadata loading.  Metadata is read as text so that
# categorical labels keep their identity ("1", "01" and "001" stay distinct);
# the frozen design_variable_types then decide how each variable enters the
# model.  Only a literal "NA" is read as missing, as before.

# Level order matches the previous numeric-then-text parsing, so reference
# levels are unchanged: labels that all parse as numbers are ordered by value
# (ties such as "1"/"01" broken by their text); other labels use sort().
nf_rna_categorical_levels <- function(values) {
  levels <- unique(values[!is.na(values)])
  numbers <- suppressWarnings(as.numeric(levels))
  if (length(levels) > 0 && !anyNA(numbers)) levels[order(numbers, levels)] else sort(levels)
}

nf_rna_design_metadata <- function(cfg, samples) {
  metadata <- read.csv(cfg$metadata, check.names = FALSE, colClasses = "character")
  rownames(metadata) <- metadata$sample_id
  metadata <- metadata[samples, , drop = FALSE]
  types <- cfg$design_variable_types
  for (variable in names(types)) {
    if (!(variable %in% colnames(metadata))) stop(paste("configured design variable is absent from metadata:", variable))
    values <- metadata[[variable]]
    if (types[[variable]] == "categorical") metadata[[variable]] <- factor(values, levels = nf_rna_categorical_levels(values))
    if (types[[variable]] == "continuous") {
      numbers <- suppressWarnings(as.numeric(values))
      invalid <- !is.finite(numbers)
      if (any(invalid)) stop(paste0("continuous design variable is not finite: ", variable, " (sample ", samples[invalid][[1]], ": '", values[invalid][[1]], "')"))
      metadata[[variable]] <- numbers
    }
  }
  # Contracts frozen before typed designs carry no types; keep the previous
  # read.csv inference for their untyped formula variables.
  for (variable in setdiff(intersect(all.vars(as.formula(cfg$formula)), colnames(metadata)), names(types))) {
    metadata[[variable]] <- type.convert(metadata[[variable]], as.is = TRUE)
  }
  if (!is.null(cfg$pair_id)) {
    if (!(cfg$pair_id %in% colnames(metadata))) stop("configured pair_id is absent from metadata")
    if (!is.factor(metadata[[cfg$pair_id]])) metadata[[cfg$pair_id]] <- factor(metadata[[cfg$pair_id]])
  }
  metadata
}
