#!/usr/bin/env Rscript
# Candidate palette review for the standard figure preset (review aid only).
#
# Nothing here changes figures.R: the candidates are NOT approved values.  The
# approved group colours 1-3 are read from the preset table so the review can
# never drift from it.  Output: palette_metrics.{tsv,json}, palette_groups.{png,pdf},
# palette_batch.{png,pdf}, palette_batch_demo.{png,pdf}.
#
# usage: Rscript palette_review.R OUTPUT_DIR

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) stop("usage: palette_review.R OUTPUT_DIR")
out <- args[[1]]; dir.create(out, recursive = TRUE, showWarnings = FALSE)
script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[[1]])
repo <- normalizePath(file.path(dirname(normalizePath(script_file)), "..", ".."))
suppressPackageStartupMessages({ library(ggplot2); library(pheatmap); library(jsonlite); library(colorspace) })
source(file.path(repo, "src", "rnaseq", "r", "figures.R"))

approved <- unname(NF_RNA_STYLE_PRESETS$standard$categorical[1:3])
stopifnot(identical(toupper(approved), c("#4E78A7", "#CC5A6A", "#59A14F")))

# ---- candidates (all unapproved) -------------------------------------------------
OI <- NF_RNA_OKABE_ITO_ORDERED
group_candidates <- list(
  G1_okabe_ito_in_order  = unname(OI[c("orange", "sky_blue", "reddish_purple", "vermillion", "yellow")]),
  G2_okabe_ito_black     = unname(OI[c("orange", "sky_blue", "reddish_purple", "vermillion", "black")]),
  G3_okabe_ito_no_vermillion = unname(OI[c("orange", "sky_blue", "reddish_purple", "black", "yellow")]),
  G4_okabe_ito_black_grey    = c(unname(OI[c("orange", "sky_blue", "reddish_purple", "black")]), "#7F7F7F")
)
batch_candidates <- list(
  B1_muted_neutrals   = c("#8C6D31", "#7F7F7F", "#9467BD", "#17BECF", "#BCBD22", "#E377C2", "#393B79", "#843C39"),
  B2_tol_muted_subset = c("#332288", "#DDCC77", "#88CCEE", "#882255", "#44AA99", "#AA4499", "#999933", "#666666"),
  B3_dark_earth_cool  = c("#6A3D9A", "#B15928", "#00838F", "#666666", "#E7298A", "#A6761D", "#3B4CC0", "#000000")
)

# ---- metrics ---------------------------------------------------------------------
source(file.path(dirname(normalizePath(script_file)), "colour_tools.R"))
rows <- list()
add_group_rows <- function(name, full, role) {
  for (v in names(visions)) {
    mp <- min_pair(full, v); low <- which.min(contrast_white(full))
    rows[[length(rows) + 1]] <<- data.frame(role = role, candidate = name, vision = v, min_de2000_within = round(mp$de, 1), closest_pair = mp$pair, min_de2000_to_groups_1_3 = NA_real_, closest_to_groups = NA_character_,
      min_contrast_vs_white = round(contrast_white(full)[[low]], 2), lowest_contrast_colour = full[[low]], stringsAsFactors = FALSE)
  }
}
add_batch_rows <- function(name, b, role) {
  for (v in names(visions)) {
    mp <- min_pair(b, v); cm <- cross_min(approved, b, v); low <- which.min(contrast_white(b))
    rows[[length(rows) + 1]] <<- data.frame(role = role, candidate = name, vision = v, min_de2000_within = round(mp$de, 1), closest_pair = mp$pair, min_de2000_to_groups_1_3 = round(cm$de, 1), closest_to_groups = cm$pair,
      min_contrast_vs_white = round(contrast_white(b)[[low]], 2), lowest_contrast_colour = b[[low]], stringsAsFactors = FALSE)
  }
}
# Typical designs use few levels, so each candidate is also scored on its first 5 (groups) / 4 (batch) colours.
for (name in names(group_candidates)) {
  full <- c(approved, group_candidates[[name]])
  add_group_rows(name, full, "groups 1-8"); add_group_rows(name, full[1:5], "groups 1-5")
}
for (name in names(batch_candidates)) {
  b <- batch_candidates[[name]]
  add_batch_rows(name, b, "batch 1-8"); add_batch_rows(name, b[1:4], "batch 1-4")
}
metrics <- do.call(rbind, rows)
write.table(metrics, file.path(out, "palette_metrics.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write(toJSON(list(approved_groups_1_3 = approved, group_candidates_4_8 = group_candidates, batch_candidates = batch_candidates,
                  thresholds = list(note = "review guides, not standards", de2000_comfortable = 10, de2000_marginal = 5, wcag_non_text_contrast = 3),
                  metrics = metrics, status = "CANDIDATES ONLY - none approved"), auto_unbox = TRUE, pretty = TRUE, na = "null"), file.path(out, "palette_metrics.json"))

# ---- figures ---------------------------------------------------------------------
font <- nf_rna_resolve_font(); family <- font$family
theme_review <- function() theme_minimal(base_size = 7, base_family = family) + theme(panel.grid = element_blank(), axis.text = element_blank(), axis.title = element_blank(),
  strip.text.y = element_text(angle = 0, hjust = 0, size = 6), strip.text.x = element_text(size = 6.5, face = "bold"), plot.title = element_text(size = 8, face = "bold"), plot.subtitle = element_text(size = 6.5))
# One row per candidate x vision; tiles are wide enough for the hex code.
swatch_plot <- function(sets, title, subtitle, divider = NULL) {
  rows <- do.call(rbind, lapply(names(sets), function(name) do.call(rbind, lapply(names(visions), function(v) {
    hex <- sets[[name]]
    data.frame(candidate = name, vision = v, row = paste0(sub("_.*", "", name), "  ", v), slot = seq_along(hex), hex = hex, drawn = visions[[v]](hex), stringsAsFactors = FALSE)
  }))))
  levels_rows <- unique(rows$row); rows$y <- length(levels_rows) + 1 - match(rows$row, levels_rows)
  labels <- stats::setNames(levels_rows, length(levels_rows) + 1 - seq_along(levels_rows))
  block_edges <- length(levels_rows) + 0.5 - length(visions) * seq_len(length(sets) - 1)
  p <- ggplot(rows, aes(slot, y)) + geom_tile(aes(fill = drawn), width = 0.95, height = 0.9, colour = "grey85", linewidth = 0.2) + scale_fill_identity() +
    geom_text(aes(label = hex), size = 1.9, colour = ifelse(contrast_white(rows$drawn) < 3, "black", "white"), family = family) +
    scale_y_continuous(breaks = as.integer(names(labels)), labels = unname(labels), expand = expansion(add = 0.5)) + scale_x_continuous(expand = expansion(add = 0.3)) +
    geom_hline(yintercept = block_edges, colour = "grey30", linewidth = 0.3) +
    labs(title = title, subtitle = paste(strwrap(subtitle, 120), collapse = "\n")) +
    theme_minimal(base_size = 7, base_family = family) + theme(panel.grid = element_blank(), axis.title = element_blank(), axis.text.x = element_blank(), axis.text.y = element_text(size = 6, hjust = 0, colour = "black"),
                                                              plot.title = element_text(size = 8, face = "bold"), plot.subtitle = element_text(size = 6.5))
  if (!is.null(divider)) p <- p + geom_vline(xintercept = divider + 0.5, colour = "black", linewidth = 0.3)
  p
}
save_review <- function(plot, stem, height_mm) {
  width_mm <- 180
  grDevices::png(file.path(out, paste0(stem, ".png")), width = round(width_mm / 25.4 * 150), height = round(height_mm / 25.4 * 150), res = 150, type = "cairo", bg = "white", family = family); print(plot); dev.off()
  grDevices::cairo_pdf(file.path(out, paste0(stem, ".pdf")), width = width_mm / 25.4, height = height_mm / 25.4, family = family, bg = "white", onefile = FALSE); print(plot); dev.off()
}
group_sets <- lapply(group_candidates, function(x) c(approved, x))
save_review(swatch_plot(group_sets, "Group colours 1-8: approved 1-3 + candidates for 4-8",
                        "Slots 1-3 are the approved colours; slots 4-8 are unapproved candidates. Columns: normal vision and simulated colour-vision deficiencies (Machado et al. 2009, severity 1).", divider = 3),
            "palette_groups", 22 + 6 * length(visions) * length(group_sets))
batch_sets <- lapply(batch_candidates, function(x) c(approved, x))
save_review(swatch_plot(batch_sets, "Batch palette candidates (slots 4-11), shown after the approved group colours (slots 1-3)",
                        "Batch is drawn as colour in the batch PCA while group is the shape; the batch colours should stay apart from group colours 1-3.", divider = 3),
            "palette_batch", 22 + 6 * length(visions) * length(batch_sets))

set.seed(3)
demo <- do.call(rbind, lapply(names(batch_candidates), function(name) {
  n <- 12; data.frame(candidate = name, PC1 = rnorm(n, rep(c(-2, 0, 2), each = 4), 0.7), PC2 = rnorm(n, 0, 1), group = rep(c("A", "B", "C"), each = 4), batch = rep(c("b1", "b2", "b3", "b4"), 3), stringsAsFactors = FALSE)
}))
demo_plot <- do.call(rbind, lapply(names(visions)[c(1, 2)], function(v) do.call(rbind, lapply(names(batch_candidates), function(name) {
  d <- demo[demo$candidate == name, ]; cols <- visions[[v]](batch_candidates[[name]])[1:4]; d$colour <- cols[match(d$batch, c("b1", "b2", "b3", "b4"))]; d$vision <- v; d }))))
demo_plot$vision <- factor(demo_plot$vision, names(visions)[1:2]); demo_plot$candidate <- factor(demo_plot$candidate, names(batch_candidates))
p <- ggplot(demo_plot, aes(PC1, PC2, colour = colour, shape = group)) + geom_point(size = 2) + scale_colour_identity() + scale_shape_manual(values = nf_rna_categorical_shapes(c("A", "B", "C"))) +
  facet_grid(vision ~ candidate) + labs(title = "Synthetic batch-PCA preview (colour = batch b1-b4, shape = group)", subtitle = "Random data for colour review only; not a real analysis.") + theme_review() + theme(legend.position = "bottom", axis.text = element_blank())
save_review(p, "palette_batch_demo", 110)

cat(sprintf("font used for review figures: %s (%s)\n", font$family, font$file))
cat(sprintf("wrote %s\n", out))
