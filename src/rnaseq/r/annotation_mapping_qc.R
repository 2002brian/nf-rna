# Shared annotation-mapping QC for the two preranked GSEA backends.
#
# Mapping quality is distinct from GSEA rank construction: callers always map
# every finite DESeq2 Wald-statistic source gene before this policy decides
# whether GSEA can proceed.
annotation_mapping_qc <- function(mapping_rate, annotation) {
  # Direct R callers from the pre-dual-threshold contract may not have passed
  # a warning threshold.  Their legacy value remains a hard block; normal
  # nf-rna execution receives the normalized contract from AnnotationConfig.
  legacy_single_threshold <- is.null(annotation$mapping_warning_rate)
  warning_threshold <- as.numeric(if (legacy_single_threshold) annotation$minimum_mapping_rate else annotation$mapping_warning_rate)
  blocking_threshold <- as.numeric(annotation$minimum_mapping_rate)
  if (!is.finite(mapping_rate) || !is.finite(warning_threshold) || !is.finite(blocking_threshold) ||
      warning_threshold < 0 || warning_threshold > 1 || blocking_threshold < 0 ||
      blocking_threshold > 1 || blocking_threshold > warning_threshold) {
    stop("invalid frozen annotation mapping-QC thresholds")
  }
  status <- if (mapping_rate < blocking_threshold) {
    "BLOCKED"
  } else if (mapping_rate < warning_threshold) {
    "WARNING"
  } else {
    "PASS"
  }
  reason <- if (status == "BLOCKED") {
    paste0("mapping rate ", sprintf("%.1f%%", 100 * mapping_rate),
           " is below blocking threshold ", sprintf("%.1f%%", 100 * blocking_threshold), ".")
  } else if (status == "WARNING") {
    paste0("mapping rate ", sprintf("%.1f%%", 100 * mapping_rate),
           " is below warning threshold ", sprintf("%.1f%%", 100 * warning_threshold),
           "; GSEA was executed.")
  } else {
    "mapping rate meets the annotation warning threshold."
  }
  list(
    status=status,
    mapping_rate=mapping_rate,
    warning_threshold=warning_threshold,
    blocking_threshold=blocking_threshold,
    reason=reason
  )
}

annotation_qc_status <- function(contrast_summaries) {
  statuses <- vapply(contrast_summaries, function(item) {
    qc <- item$ranking$annotation_qc
    if (is.null(qc$status)) "BLOCKED" else qc$status
  }, character(1))
  if (any(statuses == "BLOCKED")) "BLOCKED" else if (any(statuses == "WARNING")) "WARNING" else "PASS"
}
