#!/usr/bin/env Rscript
# Audit of the CURRENT palettes of the standard preset, with proposals for the provisional ones (review aid only).
#
# Reads every palette from the preset table in figures.R (so the audit cannot drift), scores it for distinguishability
# (CIEDE2000 under normal vision and simulated deuteranopia / protanopia / tritanopia, contrast against white) and proposes
# replacements for the PROVISIONAL slots by a deterministic maximin selection from a fixed candidate pool.  Nothing here
# edits figures.R: the approved groups 1-3 are fixed inputs and every proposal is labelled UNAPPROVED.
# Output: palette_audit.{json,tsv}, palette_audit.{png,pdf}.
#
# usage: Rscript palette_audit.R OUTPUT_DIR
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) stop("usage: palette_audit.R OUTPUT_DIR")
out <- args[[1]]; dir.create(out, recursive = TRUE, showWarnings = FALSE)
script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[[1]])
repo <- normalizePath(file.path(dirname(normalizePath(script_file)), "..", ".."))
suppressPackageStartupMessages({ library(ggplot2); library(jsonlite); library(colorspace) })
source(file.path(repo, "src", "rnaseq", "r", "figures.R"))
source(file.path(dirname(normalizePath(script_file)), "colour_tools.R"))

std <- NF_RNA_STYLE_PRESETS$standard; leg <- NF_RNA_STYLE_PRESETS$legacy
approved <- unname(std$categorical[1:3]); stopifnot(identical(toupper(approved), c("#4E78A7", "#CC5A6A", "#59A14F")))
current_groups <- unname(std$categorical); current_batch <- unname(std$batch); legacy_groups <- unname(leg$categorical)

# ---- score helpers ------------------------------------------------------------------
score_set <- function(hex) {
  do.call(rbind, lapply(names(visions), function(v) {
    mp <- min_pair(hex, v); low <- which.min(contrast_white(hex))
    data.frame(vision = v, min_de2000 = round(mp$de, 1), closest_pair = mp$pair, min_contrast_vs_white = round(contrast_white(hex)[[low]], 2), lowest_contrast_colour = hex[[low]], stringsAsFactors = FALSE)
  }))
}
worst <- function(hex) min(score_set(hex)$min_de2000)
cross <- function(a, b) do.call(rbind, lapply(names(visions), function(v) { cm <- cross_min(a, b, v); data.frame(vision = v, min_de2000 = round(cm$de, 1), pair = cm$pair, stringsAsFactors = FALSE) }))

# ---- deterministic maximin proposal ---------------------------------------------------
pool <- c("#F28E2B", "#E15759", "#76B7B2", "#EDC948", "#B07AA1", "#FF9DA7", "#9C755F", "#BAB0AC",
          "#E69F00", "#56B4E9", "#CC79A7", "#D55E00", "#F0E442", "#009E73", "#0072B2", "#000000",
          "#7F7F7F", "#8C6D31", "#9467BD", "#17BECF", "#BCBD22", "#393B79", "#843C39", "#6A3D9A",
          "#B15928", "#00838F", "#E7298A", "#A6761D", "#666666", "#E08214", "#7B3294", "#1B9E77",
          "#D95F02", "#7570B3", "#E6AB02", "#2E2E2E")
lch <- function(hex) { L <- lab(hex); cbind(L = L[, 1], C = sqrt(L[, 2]^2 + L[, 3]^2)) }
# Every slot maximises the smallest CIEDE2000 distance (over the four vision types) to the slots already fixed or chosen.
# `usable` restricts the pool: group colours must be chromatic, mid-lightness colours (no near-black or grey that a marker
# would share with axis ink or the "not significant" grey); batch colours may be muted but never near-black.
greedy <- function(fixed, n, usable, avoid = character()) {
  chosen <- character()
  for (k in seq_len(n)) {
    cand <- setdiff(pool[usable(pool)], c(fixed, chosen, avoid))
    anchor <- c(fixed, chosen)
    scores <- vapply(cand, function(h) min(vapply(names(visions), function(v) cross_min(anchor, h, v)$de, numeric(1))), numeric(1))
    chosen <- c(chosen, cand[[which.max(scores)]])   # first maximum in pool order: deterministic
  }
  chosen
}
usable_group <- function(h) { m <- lch(h); contrast_white(h) >= 2.0 & m[, "C"] >= 25 & m[, "L"] >= 30 }
usable_batch <- function(h) { m <- lch(h); contrast_white(h) >= 2.0 & m[, "L"] >= 25 }
proposal_groups <- c(approved, greedy(approved, 5, usable_group))
proposal_batch <- greedy(approved, 8, usable_batch, avoid = proposal_groups[4:8])   # away from groups 1-3 and from the proposed groups 4-8

# ---- group / batch tables -----------------------------------------------------------------------
sets <- list(legacy_okabe_ito_1_8 = legacy_groups, standard_current_groups_1_8 = current_groups, standard_proposed_groups_1_8 = proposal_groups)
rows <- list()
for (name in names(sets)) for (k in c(3, 5, 8)) { s <- score_set(sets[[name]][1:k]); s$set <- name; s$slots <- k; rows[[length(rows) + 1]] <- s }
for (name in c("current", "proposed")) {
  b <- if (name == "current") current_batch else proposal_batch
  for (k in c(4, 8)) { s <- score_set(b[1:k]); s$set <- paste0("batch_", name); s$slots <- k; rows[[length(rows) + 1]] <- s }
}
metrics <- do.call(rbind, rows)[, c("set", "slots", "vision", "min_de2000", "closest_pair", "min_contrast_vs_white", "lowest_contrast_colour")]
batch_vs_groups <- list(current = cross(approved, current_batch), proposed = cross(approved, proposal_batch),
                        current_vs_groups_1_8 = cross(current_groups, current_batch), proposed_vs_groups_1_8 = cross(proposal_groups, proposal_batch))

# ---- continuous scales ---------------------------------------------------------------------------
ramp <- function(stops, n = 9) colorRampPalette(stops)(n)
steps_uniform <- function(hex) { d <- de2000(lab(hex)[-length(hex), , drop = FALSE], lab(hex)[-1, , drop = FALSE]); round(c(min = min(d), max = max(d), ratio = max(d) / min(d)), 2) }
cont <- list(
  deg_heatmap_diverging_standard = list(kind = "diverging, midpoint 0", stops = unname(std$diverging), role = "DEG heatmap z-scores; correlation heatmap when negative correlations occur"),
  deg_heatmap_diverging_legacy = list(kind = "diverging, midpoint 0", stops = unname(leg$diverging), role = "legacy ends (RdBu)"),
  correlation_sequential_standard = list(kind = "sequential (Blues)", stops = unname(std$sequential), role = "correlation heatmap when all r >= 0 (same stops in legacy and standard)"),
  enrichment_standard_proposed = list(kind = "sequential, dark = most significant", stops = unname(std$enrichment), role = "GSEA / KEGG-GSEA adjusted P; ORA mark colour is the first stop"),
  enrichment_legacy = list(kind = "sequential, dark = most significant", stops = rev(leg$sequential[-1]), role = "legacy GSEA scale")
)
cont_rows <- list()
for (name in names(cont)) {
  stops <- cont[[name]]$stops; samples <- ramp(stops, 9)
  for (v in names(visions)) {
    sv <- visions[[v]](samples); Lv <- lab(sv)[, 1]
    to_groups <- min(vapply(1:3, function(i) min(de2000(lab(sv), lab(rep(visions[[v]](approved[[i]]), length(sv))))), numeric(1)))
    cont_rows[[length(cont_rows) + 1]] <- data.frame(scale = name, kind = cont[[name]]$kind, vision = v,
      end_to_end_de2000 = round(de2000(lab(sv[1]), lab(sv[length(sv)])), 1),
      lightness_monotonic = all(diff(Lv) >= -0.5) || all(diff(Lv) <= 0.5),
      min_de2000_to_group_colours_1_3 = round(to_groups, 1),
      darkest_end_contrast_vs_white = round(max(contrast_white(stops)), 2), lightest_end_contrast_vs_white = round(min(contrast_white(stops)), 2), stringsAsFactors = FALSE)
  }
}
cont_metrics <- do.call(rbind, cont_rows)
uniformity <- lapply(cont, function(x) as.list(steps_uniform(ramp(x$stops, 9))))

write.table(rbind(data.frame(table = "discrete", metrics, stringsAsFactors = FALSE)[, c("table", "set", "slots", "vision", "min_de2000", "min_contrast_vs_white")]), file.path(out, "palette_audit.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(cont_metrics, file.path(out, "palette_audit_continuous.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write(toJSON(list(status = "AUDIT ONLY - proposals are UNAPPROVED; no preset value was changed",
                  approved_groups_1_3 = approved, current = list(groups_1_8 = current_groups, batch_1_8 = current_batch, legacy_groups = legacy_groups),
                  proposal = list(groups_1_8 = proposal_groups, batch_1_8 = proposal_batch,
                                  method = "greedy maximin over a fixed candidate pool: each slot maximises the smallest CIEDE2000 distance (normal, deuteranopia, protanopia, tritanopia) to all earlier slots; pool restricted to contrast vs white >= 2.0 and, for groups, chroma >= 25 and L >= 30; batch slots are scored against groups 1-3 and exclude the proposed groups 4-8"),
                  thresholds = list(note = "review guides, not standards", de2000_comfortable = 10, de2000_marginal = 5, wcag_non_text_contrast = 3),
                  discrete_metrics = metrics, batch_vs_groups = batch_vs_groups, continuous = cont, continuous_metrics = cont_metrics, step_uniformity = uniformity),
            auto_unbox = TRUE, pretty = TRUE, na = "null"), file.path(out, "palette_audit.json"))

# ---- figure ------------------------------------------------------------------------------------------
family <- tryCatch(nf_rna_require_strict_font(std$font$preferred, std$font$bundled)$family, error = function(e) "sans")
swatch <- function(sets, title, subtitle, divider = NULL) {
  rows <- do.call(rbind, lapply(names(sets), function(name) do.call(rbind, lapply(names(visions), function(v) {
    hex <- sets[[name]]; data.frame(row = paste0(name, "  ", v), slot = seq_along(hex), hex = hex, drawn = visions[[v]](hex), stringsAsFactors = FALSE) }))))
  lv <- unique(rows$row); rows$y <- length(lv) + 1 - match(rows$row, lv)
  p <- ggplot(rows, aes(slot, y)) + geom_tile(aes(fill = drawn), width = 0.95, height = 0.9, colour = "grey85", linewidth = 0.2) + scale_fill_identity() +
    geom_text(aes(label = hex), size = 1.9, colour = ifelse(contrast_white(rows$drawn) < 3, "black", "white"), family = family) +
    scale_y_continuous(breaks = length(lv) + 1 - seq_along(lv), labels = lv, expand = expansion(add = 0.5)) +
    geom_hline(yintercept = length(lv) + 0.5 - length(visions) * seq_len(length(sets) - 1), colour = "grey30", linewidth = 0.3) +
    labs(title = title, subtitle = subtitle) + theme_minimal(base_size = 7, base_family = family) +
    theme(panel.grid = element_blank(), axis.title = element_blank(), axis.text.x = element_blank(), axis.text.y = element_text(size = 6, hjust = 0, colour = "black"), plot.title = element_text(size = 8, face = "bold"), plot.subtitle = element_text(size = 6.5), plot.title.position = "plot")
  if (!is.null(divider)) p <- p + geom_vline(xintercept = divider + 0.5, colour = "black", linewidth = 0.3)
  p
}
save_both <- function(plot, stem, height_mm) {
  grDevices::png(file.path(out, paste0(stem, ".png")), width = round(180 / 25.4 * 150), height = round(height_mm / 25.4 * 150), res = 150, type = "cairo", bg = "white", family = family); print(plot); dev.off()
  grDevices::cairo_pdf(file.path(out, paste0(stem, ".pdf")), width = 180 / 25.4, height = height_mm / 25.4, family = family, bg = "white", onefile = FALSE); print(plot); dev.off()
}
save_both(swatch(list(legacy_OkabeIto = legacy_groups, standard_current = current_groups, standard_PROPOSED = proposal_groups),
                 "Group colours 1-8: legacy, current standard (4-8 provisional) and unapproved proposal", "Slots 1-3 are the approved colours in both standard rows. Rows per set: normal vision, then simulated deuteranopia / protanopia / tritanopia.", divider = 3), "palette_audit_groups", 105)
save_both(swatch(list(batch_current = current_batch, batch_PROPOSED = proposal_batch), "Batch colours 1-8: current (provisional) and unapproved proposal", "Batch colours appear next to group colours in annotation bars, so they are scored against groups 1-3 as well."), "palette_audit_batch", 75)
strips <- do.call(rbind, lapply(names(cont), function(n) do.call(rbind, lapply(names(visions), function(v) data.frame(scale = n, vision = v, x = seq(0, 1, length.out = 101), colour = visions[[v]](colorRampPalette(cont[[n]]$stops)(101)))))))
strips$scale <- factor(strips$scale, names(cont)); strips$vision <- factor(strips$vision, names(visions))
p <- ggplot(strips, aes(x, 1, fill = colour)) + geom_tile(height = 1, width = 0.01) + scale_fill_identity() + facet_grid(scale ~ vision) +
  labs(title = "Continuous scales: expression (diverging), correlation (sequential) and enrichment (sequential)", subtitle = "Normal vision and simulated deuteranopia / protanopia / tritanopia. Audit only: nothing is approved or changed.", x = NULL, y = NULL) +
  theme_minimal(base_family = family, base_size = 7) + theme(axis.text = element_blank(), panel.grid = element_blank(), strip.text.y = element_text(angle = 0, hjust = 0, size = 5.5), plot.title.position = "plot")
save_both(p, "palette_audit_continuous", 90)
cat(sprintf("wrote %s\n", out))
