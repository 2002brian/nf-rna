#!/usr/bin/env Rscript
# Validation of the standard preset's palettes as implemented (review aid): approved group slots 1-8, the deterministic
# extension from slot 9, and the batch palette candidate.  Scores CIEDE2000 under normal vision and simulated
# deuteranopia / protanopia / tritanopia and contrast against white.  Output: JSON on stdout.
script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[[1]])
repo <- normalizePath(file.path(dirname(normalizePath(script_file)), "..", ".."))
suppressPackageStartupMessages({ library(jsonlite); library(colorspace) })
source(file.path(repo, "src", "rnaseq", "r", "figures.R"))
source(file.path(dirname(normalizePath(script_file)), "colour_tools.R"))
std <- NF_RNA_STYLE_PRESETS$standard
score <- function(hex) sapply(names(visions), function(v) round(min_pair(hex, v)$de, 1))
groups <- unname(std$categorical); batch <- unname(std$batch)
out <- list(groups_1_8 = list(colours = groups, min_de2000_by_vision = as.list(score(groups)), min_contrast_vs_white = round(min(contrast_white(groups)), 2)),
            batch_1_8 = list(colours = batch, min_de2000_by_vision = as.list(score(batch)), min_contrast_vs_white = round(min(contrast_white(batch)), 2),
                             min_de2000_to_groups_1_3_by_vision = as.list(sapply(names(visions), function(v) round(cross_min(groups[1:3], batch, v)$de, 1)))),
            extension = lapply(setNames(9:12, paste0("n", 9:12)), function(n) {
              e <- nf_rna_extend_palette(groups, n)
              list(added = e[9:n], min_de2000_by_vision_all_slots = as.list(score(e)), min_de2000_new_to_slots_1_8 = as.list(sapply(names(visions), function(v) round(cross_min(groups, e[9:n], v)$de, 1))),
                   min_contrast_vs_white_new = round(min(contrast_white(e[9:n])), 2))
            }),
            batch_extension_n9 = {e <- nf_rna_extend_palette(batch, 9); list(added = e[9], min_de2000_to_batch_by_vision = as.list(sapply(names(visions), function(v) round(cross_min(batch, e[9], v)$de, 1))))})
cat(toJSON(out, auto_unbox = TRUE, pretty = TRUE))
