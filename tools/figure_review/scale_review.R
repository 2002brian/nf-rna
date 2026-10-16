#!/usr/bin/env Rscript
# Review of the CONTINUOUS colour scales of the standard preset (review aid only; nothing here edits figures.R).
#
# 1. The expression (diverging) scale of the DEG heatmap versus alternatives, with the approved discrete group colours
#    alongside, because the diverging ends of the standard preset coincide with group colours 1 and 2.
# 2. The proposed enrichment scale (adjusted P in GSEA dotplots) versus the legacy sequential one.
# Output: scale_review.{png,pdf} and scale_metrics.json.  Candidates are NOT approved; no scale, limit or scaling is changed.
#
# usage: Rscript scale_review.R OUTPUT_DIR
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) stop("usage: scale_review.R OUTPUT_DIR")
out <- args[[1]]; dir.create(out, recursive = TRUE, showWarnings = FALSE)
script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[[1]])
repo <- normalizePath(file.path(dirname(normalizePath(script_file)), "..", ".."))
suppressPackageStartupMessages({ library(ggplot2); library(jsonlite); library(colorspace) })
source(file.path(repo, "src", "rnaseq", "r", "figures.R"))
source(file.path(dirname(normalizePath(script_file)), "colour_tools.R"))

std <- NF_RNA_STYLE_PRESETS$standard; leg <- NF_RNA_STYLE_PRESETS$legacy
groups <- unname(std$categorical[1:3])
diverging <- list(
  current_standard = unname(std$diverging),                           # #4E78A7 / #F7F7F7 / #CC5A6A (approved earlier)
  legacy_rdbu_ends = unname(leg$diverging),                            # #2166AC / #F7F7F7 / #B2182B (also the volcano colours)
  candidate_teal_orange = c("#2C7C94", "#F7F7F7", "#C8651B")           # unapproved: avoids the group colours
)
sequential <- list(
  legacy_blues_gsea = rev(leg$sequential[-1]),                         # most significant first
  proposed_enrichment = unname(std$enrichment)
)
ramp <- function(stops, n = 101) colorRampPalette(stops)(n)

metric_rows <- list()
for (name in names(diverging)) {
  ends <- diverging[[name]]; half_lo <- ramp(ends[1:2], 50); half_hi <- ramp(ends[2:3], 50)
  mono <- function(x) { L <- lab(x)[, 1]; all(diff(L) >= -0.01) || all(diff(L) <= 0.01) }
  for (v in names(visions)) {
    L <- lab(visions[[v]](ends)); metric_rows[[length(metric_rows) + 1]] <- data.frame(scale = name, kind = "diverging", vision = v,
      de2000_between_ends = round(de2000(L[1, , drop = FALSE], L[3, , drop = FALSE]), 1),
      de2000_end_to_group1 = round(de2000(L[1, , drop = FALSE], lab(visions[[v]](groups[[1]]))), 1),
      de2000_end_to_group2 = round(de2000(L[3, , drop = FALSE], lab(visions[[v]](groups[[2]]))), 1),
      end_contrast_vs_white = paste(round(contrast_white(ends[c(1, 3)]), 2), collapse = " / "), stringsAsFactors = FALSE)
  }
}
for (name in names(sequential)) {
  stops <- sequential[[name]]; L <- lab(stops)[, 1]
  for (v in names(visions)) {
    Lv <- lab(visions[[v]](ramp(stops, 9)))[, 1]
    metric_rows[[length(metric_rows) + 1]] <- data.frame(scale = name, kind = "sequential", vision = v,
      de2000_between_ends = round(de2000(lab(visions[[v]](stops[1])), lab(visions[[v]](stops[length(stops)]))), 1),
      de2000_end_to_group1 = NA_real_, de2000_end_to_group2 = NA_real_,
      end_contrast_vs_white = paste(round(contrast_white(stops[c(1, length(stops))]), 2), collapse = " / "),
      lightness_monotonic = all(diff(Lv) >= -0.5) || all(diff(Lv) <= 0.5), stringsAsFactors = FALSE)
  }
}
metrics <- do.call(rbind, lapply(metric_rows, function(r) { r$lightness_monotonic <- if (is.null(r$lightness_monotonic)) NA else r$lightness_monotonic; r }))
write(toJSON(list(status = "REVIEW ONLY - no scale approved or changed", note = "The standard diverging ends equal group colours 1 and 2 (de2000 = 0); the annotation bar and the expression extremes therefore share colours.",
                  diverging = diverging, sequential = sequential, group_colours_1_3 = groups, metrics = metrics), auto_unbox = TRUE, pretty = TRUE), file.path(out, "scale_metrics.json"))
write.table(metrics, file.path(out, "scale_metrics.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# ---- figure: each scale as a ramp under normal vision and the three CVD simulations
family <- tryCatch(nf_rna_require_strict_font(std$font$preferred, std$font$bundled)$family, error = function(e) "sans")
strips <- do.call(rbind, c(
  lapply(names(diverging), function(n) do.call(rbind, lapply(names(visions), function(v) data.frame(scale = paste0("expression: ", n), vision = v, x = seq(0, 1, length.out = 101), colour = visions[[v]](ramp(diverging[[n]])))))),
  lapply(names(sequential), function(n) do.call(rbind, lapply(names(visions), function(v) data.frame(scale = paste0("enrichment: ", n), vision = v, x = seq(0, 1, length.out = 101), colour = visions[[v]](ramp(sequential[[n]]))))))))
strips$scale <- factor(strips$scale, unique(strips$scale)); strips$vision <- factor(strips$vision, names(visions))
p <- ggplot(strips, aes(x, 1, fill = colour)) + geom_tile(height = 1, width = 0.01) + scale_fill_identity() +
  facet_grid(scale ~ vision) + labs(title = "Continuous scales of the standard preset: review", subtitle = "Normal vision and simulated deuteranopia / protanopia / tritanopia. Review only: nothing is approved or changed.", x = "expression: low z (left) to high z (right); enrichment: most significant adjusted P (left) to cutoff (right)", y = NULL) +
  theme_minimal(base_family = family, base_size = 7) + theme(axis.text = element_blank(), panel.grid = element_blank(), strip.text.y = element_text(angle = 0, hjust = 0), plot.title.position = "plot")
for (ext in c("png", "pdf")) {
  path <- file.path(out, paste0("scale_review.", ext))
  if (ext == "png") grDevices::png(path, width = 1500, height = 900, res = 150, type = "cairo", bg = "white", family = family) else grDevices::cairo_pdf(path, width = 10, height = 6, family = family, bg = "white", onefile = FALSE)
  print(p); dev.off()
}
cat(sprintf("wrote %s\n", out))
