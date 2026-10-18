# Central figure policy for every nf-rna R backend: output profiles,
# typography, colour policies, display-only transformations and export.
# Nothing here computes or alters a scientific result; helpers that transform
# values (z-scores, colour saturation) return new objects for drawing only.

`%||%` <- function(value, fallback) if (is.null(value) || length(value) == 0) fallback else value

# Physical targets per profile.  These are layout targets taken from publisher
# artwork guidance (Nature research figures; PLOS ONE figure requirements;
# Elsevier artwork guidance); they are not a guarantee that any journal will
# accept a figure.  Production uses "general" unless a config names another.
NF_RNA_FIGURE_PROFILES <- list(
  general=list(name="general", single_col_mm=85, double_col_mm=180, max_width_mm=180, max_height_mm=220,
               base_pt=8, small_pt=7, min_pt=6, raster_dpi=300, preview_dpi=150,
               note="nf-rna default; a conservative common layout, not a journal specification"),
  nature=list(name="nature", single_col_mm=89, double_col_mm=183, max_width_mm=183, max_height_mm=170,
              base_pt=7, small_pt=6, min_pt=5, raster_dpi=300, preview_dpi=150,
              note="Nature research figures: 89/183 mm columns, <=170 mm high, 5-7 pt text"),
  plos=list(name="plos", single_col_mm=NA_real_, double_col_mm=190.5, max_width_mm=190.5, max_height_mm=222.3,
            base_pt=9, small_pt=8, min_pt=8, raster_dpi=300, preview_dpi=150,
            note="PLOS ONE: <=190.5 x 222.3 mm, 8-12 pt, 300-600 dpi RGB TIFF (LZW, no alpha)"),
  elsevier=list(name="elsevier", single_col_mm=90, double_col_mm=190, max_width_mm=190, max_height_mm=240,
                base_pt=8, small_pt=7, min_pt=7, raster_dpi=500, preview_dpi=150,
                note="Elsevier artwork: halftone 300, combination 500, line art 1000 dpi; vector preferred")
)
NF_RNA_FONT_PREFERRED <- "Arial"
# DejaVu Sans ships with the locked Conda runtime (fonts-conda-forge); it is the
# explicit, recorded substitute when Arial is not installed.
NF_RNA_FONT_FALLBACK <- "DejaVu Sans"

nf_rna_figure_profile <- function(name=NULL) {
  name <- name %||% "general"
  profile <- NF_RNA_FIGURE_PROFILES[[name]]
  if (is.null(profile)) stop(paste("unknown figure profile:", name))
  profile
}

nf_rna_resolve_font <- function(preferred=NF_RNA_FONT_PREFERRED, fallback=NF_RNA_FONT_FALLBACK) {
  if (!requireNamespace("systemfonts", quietly=TRUE)) {
    return(list(requested=preferred, family=fallback, fallback_used=TRUE, verified=FALSE, file=NA_character_,
                reason="systemfonts unavailable; installed fonts could not be inspected, so the fallback family is used"))
  }
  families <- unique(systemfonts::system_fonts()$family)
  family <- if (preferred %in% families) preferred else if (fallback %in% families) fallback else NA_character_
  if (is.na(family)) stop(paste0("Neither ", preferred, " nor the fallback ", fallback, " is installed; refusing to draw with an unknown font."))
  match <- if (exists("match_fonts", asNamespace("systemfonts"))) systemfonts::match_fonts(family) else systemfonts::match_font(family)
  list(requested=preferred, family=family, fallback_used=!identical(family, preferred), verified=TRUE,
       file=as.character(match$path[[1]]),
       reason=if (identical(family, preferred)) NULL else paste0(preferred, " is not installed in this runtime"))
}

# ------------------------------------------------------------------ style presets

# The preset controls the visual style; the profile (above) controls the
# publication specification.  They are independent: every preset works with every
# profile.  Preset values live only in this table; scripts read them through the
# accessors below and never hard-code a colour, shape or label size.
#
# legacy   : the appearance nf-rna has always produced (kept exactly).
# standard : nf-rna's publication style (its own conventions, defined here and
#            documented in docs/figure-standard.md).  "bioinfo" is a
#            backward-compatible alias that resolves to "standard".  It is applied
#            figure type by figure type: only the types in NF_RNA_STANDARD_FIGURES
#            are drawn with it in this release; every other figure keeps the legacy
#            appearance (and the legacy font) until it is migrated.
NF_RNA_PRESET_ALIASES <- c(bioinfo="standard")
NF_RNA_STANDARD_FIGURES <- c("volcano", "pca", "boxplot", "correlation", "heatmap", "dotplot")
NF_RNA_OKABE_ITO_ORDERED <- c(blue="#0072B2", vermillion="#D55E00", bluish_green="#009E73", reddish_purple="#CC79A7",
                              orange="#E69F00", sky_blue="#56B4E9", yellow="#F0E442", black="#000000")
NF_RNA_STYLE_PRESETS <- list(
  legacy=list(
    name="legacy",
    # Okabe & Ito (2008) colour-universal-design palette, ordered for contrast on white.
    categorical=NF_RNA_OKABE_ITO_ORDERED,
    categorical_name="Okabe-Ito",
    batch=NULL,                                  # NULL: batch uses the categorical palette
    shapes=c(16, 17, 15, 18, 1, 2, 0, 5, 3, 4, 8, 6),
    direction=c(Up="#D55E00", Down="#0072B2", "Not significant"="#BDBDBD"),
    diverging=c(low="#2166AC", mid="#F7F7F7", high="#B2182B"),   # ColorBrewer RdBu
    sequential=c("#F7FBFF", "#C6DBEF", "#6BAED6", "#2171B5", "#08306B"),  # ColorBrewer Blues
    missing="#A0A0A0",
    primary="#0072B2",
    pca_label_pt=NA_real_,                       # NA: the profile's min_pt
    point_alpha=NA_real_,                        # NA: opaque, no alpha channel requested
    volcano_legend_counts=TRUE,
    enrichment=NULL,                             # NULL: reverse of the sequential stops without the lightest (the historic GSEA scale)
    layout=NULL,                                 # NULL: each script's historical size
    font=list(policy="preferred_with_recorded_fallback", preferred="Arial", fallback="DejaVu Sans")
  ),
  standard=list(
    name="standard",
    # Group slots 1-8 are approved (1-3 earlier, 4-8 in Phase 4).  Levels beyond slot 8 keep slots 1-8 and receive
    # deterministic maximin extension colours (nf_rna_extend_palette).
    categorical=c("#4E78A7", "#CC5A6A", "#59A14F", "#843C39", "#56B4E9", "#BCBD22", "#6A3D9A", "#8C6D31"),
    categorical_name="standard groups 1-8 (approved)",
    categorical_extension_name="standard groups 1-8 (approved) + deterministic maximin extension from slot 9",
    # Batch palette: approved as "nf-rna Standard Batch Palette v1" (slot 2 is the darker neutral that replaced #BAB0AC).
    batch=c("#393B79", "#585858", "#17BECF", "#E6AB02", "#CC79A7", "#A6761D", "#7F7F7F", "#E08214"),
    batch_name="nf-rna Standard Batch Palette v1",
    batch_extension_name="nf-rna Standard Batch Palette v1 + deterministic maximin extension from slot 9",
    shapes=c(16, 17, 15, 18, 1, 2, 0, 5, 3, 4, 8, 6),
    direction=c(Up="#B2182B", Down="#2166AC", "Not significant"="#BDBDBD"),
    diverging=c(low="#4E78A7", mid="#F7F7F7", high="#CC5A6A"),
    sequential=c("#F7FBFF", "#C6DBEF", "#6BAED6", "#2171B5", "#08306B"),
    missing="#A0A0A0",
    primary="#4E78A7",
    pca_label_pt=6,
    compact_width_mm=85,                         # single-column width for profiles that have no single column (PLOS)
    point_alpha=0.85,
    volcano_legend_counts=TRUE,
    # Raster dpi is not a preset property: the selected profile's raster_dpi always applies (300, Elsevier 500).
    # Figure widths are columns of the selected profile ("single" uses compact_width_mm where the profile has no
    # single column); heights follow from the aspect ratio or the data.
    # Continuous enrichment scale (adjusted P in the GSEA dotplots; the ORA mark colour is its first stop).  Ordered from the
    # most significant (dark) to the least significant (light) so low P is dark, as in legacy.  It is a sequential scale of its
    # own (Viridis hues), deliberately separate from the discrete group palette and from the diverging expression scale.
    # Approved.
    enrichment=c("#4B0055", "#353E7C", "#007094", "#009B95", "#00BE7D"),
    # Content-aware sizes: widths and heights follow the data and the measured labels, never above the
    # profile's maximum.  Cell and panel sizes are in mm.
    layout=list(
      volcano=list(width="single", aspect=1.0, legend_ncol=1),
      pca=list(width="single", min_height_mm=40),
      boxplot=list(min_width_mm=60, sample_mm=3.4, panel_height_mm=42, pad_mm=4, max_legend_rows=4),
      correlation=list(cell_mm=6, min_cell_mm=1.2, tree_mm=9, matrix_min_mm=32, label_fit=1.1, pad_mm=2),
      heatmap=list(cell_mm=5, min_cell_mm=0.8, tree_mm=9, hidden_row_mm=1.5, matrix_min_mm=32, label_fit=1.25, pad_mm=2),
      dotplot=list(panel_mm=55, min_width_mm=85, row_gap_mm=2.2, line_factor=1.25, min_panel_height_mm=28, pad_mm=2)
    ),
    font=list(policy="strict", preferred="Arimo", fallback=NULL, bundled="arimo")
  )
)
NF_RNA_PRESET_NAMES <- names(NF_RNA_STYLE_PRESETS)

# Canonical preset name: aliases ("bioinfo") resolve to "standard"; anything else is unchanged.
nf_rna_canonical_preset <- function(name) {
  if (is.character(name) && length(name) == 1 && !is.na(name) && name %in% names(NF_RNA_PRESET_ALIASES)) unname(NF_RNA_PRESET_ALIASES[[name]]) else name
}

# Read and validate the resolved figure style from a backend config.  Python
# sends figure_style = {preset, profile, source}.  The older cfg$figures$profile
# hook is still accepted when figure_style is absent; if both are present they
# must agree.  A configuration with neither is the legacy preset, general profile.
nf_rna_read_figure_style <- function(cfg) {
  style <- cfg$figure_style
  hook <- if (is.list(cfg$figures)) cfg$figures$profile else NULL
  if (!is.null(style) && !is.list(style)) stop("figure_style must be an object with preset and profile")
  given <- style$preset %||% "legacy"
  preset <- nf_rna_canonical_preset(given)
  profile <- style$profile %||% hook %||% "general"
  if (!is.null(style$profile) && !is.null(hook) && !identical(style$profile, hook)) {
    stop(sprintf("conflicting figure profiles: figure_style.profile is '%s' but figures.profile is '%s'", style$profile, hook))
  }
  if (!is.character(preset) || length(preset) != 1 || !(preset %in% NF_RNA_PRESET_NAMES)) {
    stop(paste0("unknown figure preset: ", paste(given, collapse=","), " (expected one of ", paste(c(NF_RNA_PRESET_NAMES, names(NF_RNA_PRESET_ALIASES)), collapse=", "), ")"))
  }
  nf_rna_figure_profile(profile)  # stops on an unknown profile; never silently ignored
  list(preset=preset, given=given, profile=profile, source=style$source %||% NA_character_)
}

# ---- strict fonts (presets that forbid substitution: standard -> Arimo) ----
#
# Arimo (SIL Open Font License 1.1) is bundled with nf-rna in rnaseq/fonts/arimo
# with its license and a manifest pinning the SHA-256 of every file.  Nothing is
# looked up elsewhere and nothing is mapped onto another family:
#   1. nf_rna_verify_bundled_font(): the bundled files exist and match the pinned SHA-256.
#   2. nf_rna_activate_bundled_font(): a generated fontconfig file adds ONLY the
#      bundled directory (and hides any other Arimo build) before fontconfig is
#      first used by this process.
#   3. nf_rna_strict_font_report(): systemfonts lists exactly the bundled files.
#   4. nf_rna_verify_drawn_font(): a probe PDF is drawn in all four faces and the
#      embedded font names are read back; fontconfig otherwise maps unknown
#      families to another font silently (Arial -> Liberation Sans, ...).
#   5. nf_rna_require_strict_font(): all of the above must pass or the run stops
#      with an actionable error before any figure is generated.

# Directory of the R scripts (set while this file is sourced).
NF_RNA_R_DIR <- local({
  found <- NA_character_
  for (i in rev(seq_len(sys.nframe()))) {
    ofile <- tryCatch(get("ofile", envir=sys.frame(i), inherits=FALSE), error=function(e) NULL)
    if (is.character(ofile) && length(ofile) == 1 && identical(basename(ofile), "figures.R")) { found <- dirname(normalizePath(ofile)); break }
  }
  if (is.na(found)) {
    script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value=TRUE))
    if (length(script)) found <- dirname(normalizePath(script[[1]]))
  }
  found
})
NF_RNA_FONT_STATE <- new.env(parent=emptyenv())

nf_rna_file_sha256 <- function(path) {
  if (!exists("sha256sum", envir=asNamespace("tools"), inherits=FALSE)) return(NA_character_)  # R < 4.5
  unname(tools::sha256sum(path))
}

nf_rna_bundled_font_dir <- function(name) {
  if (is.na(NF_RNA_R_DIR)) stop("cannot locate the nf-rna R scripts, so the bundled font directory is unknown", call.=FALSE)
  file.path(dirname(NF_RNA_R_DIR), "fonts", name)
}

nf_rna_font_requirement <- function(family, bundled=NULL) {
  if (is.null(bundled)) {
    return(paste0("Install ", family, " so that systemfonts lists the exact family '", family, "'; Arial, Helvetica, Liberation Sans or DejaVu Sans are never substituted."))
  }
  paste0(family, " (SIL Open Font License 1.1) is bundled with nf-rna in rnaseq/fonts/", bundled,
         "; reinstall nf-rna if those files are missing or modified. Arial, Helvetica, Liberation Sans or DejaVu Sans are never substituted.")
}

# Check the bundled files against the manifest; returns the manifest plus the actual hashes.
nf_rna_verify_bundled_font <- function(name) {
  dir <- nf_rna_bundled_font_dir(name)
  manifest_path <- file.path(dir, "FONT_MANIFEST.json")
  if (!file.exists(manifest_path)) stop(paste0("bundled font manifest is missing: ", manifest_path), call.=FALSE)
  manifest <- jsonlite::read_json(manifest_path)
  faces <- unlist(manifest$faces)
  files <- list()
  for (file in c(unname(faces), manifest$license_file)) {
    path <- file.path(dir, file)
    if (!file.exists(path)) stop(paste0("bundled font file is missing: ", path, ". ", nf_rna_font_requirement(manifest$family, name)), call.=FALSE)
    actual <- nf_rna_file_sha256(path)
    if (!is.na(actual) && !identical(actual, manifest$files[[file]]$sha256)) {
      stop(paste0("bundled font file ", file, " does not match its pinned SHA-256 (", actual, " != ", manifest$files[[file]]$sha256, "). ", nf_rna_font_requirement(manifest$family, name)), call.=FALSE)
    }
    files[[file]] <- list(path=normalizePath(path), sha256=actual)
  }
  list(dir=normalizePath(dir), manifest=manifest, files=files, faces=as.list(faces))
}

# Point fontconfig at the bundled directory.  Must run before this R process first
# uses fontconfig (it is read once); a late call is detected by the drawn-font probe.
nf_rna_activate_bundled_font <- function(name) {
  if (!is.null(NF_RNA_FONT_STATE$config)) return(invisible(NF_RNA_FONT_STATE$config))
  bundle <- nf_rna_verify_bundled_font(name)
  prefix <- normalizePath(file.path(R.home(), "..", ".."), mustWork=FALSE)
  candidates <- c(Sys.getenv("FONTCONFIG_FILE"), file.path(Sys.getenv("FONTCONFIG_PATH"), "fonts.conf"), file.path(prefix, "etc", "fonts", "fonts.conf"), "/etc/fonts/fonts.conf")
  base <- Filter(function(p) nzchar(p) && file.exists(p), candidates)
  escape <- function(x) gsub(">", "&gt;", gsub("<", "&lt;", gsub("&", "&amp;", x, fixed=TRUE), fixed=TRUE), fixed=TRUE)
  family <- bundle$manifest$family
  lines <- c('<?xml version="1.0"?>', '<!DOCTYPE fontconfig SYSTEM "urn:fontconfig:fonts.dtd">', "<fontconfig>",
             if (length(base)) sprintf('  <include ignore_missing="yes">%s</include>', escape(normalizePath(base[[1]]))),
             # hide every other build of the family so that exactly the pinned files can be drawn
             sprintf('  <selectfont><rejectfont><pattern><patelt name="family"><string>%s</string></patelt></pattern></rejectfont>', escape(family)),
             paste0("    <acceptfont>", paste(sprintf('<pattern><patelt name="file"><string>%s</string></patelt></pattern>', escape(vapply(bundle$files[unlist(bundle$faces)], function(x) x$path, character(1)))), collapse=""), "</acceptfont></selectfont>"),
             sprintf("  <dir>%s</dir>", escape(bundle$dir)), "</fontconfig>")
  config <- tempfile("nf-rna-fonts-", fileext=".conf")
  writeLines(lines, config)
  NF_RNA_FONT_STATE$previous <- Sys.getenv("FONTCONFIG_FILE", unset=NA)
  Sys.setenv(FONTCONFIG_FILE=config)
  NF_RNA_FONT_STATE$config <- config
  invisible(config)
}

nf_rna_strict_font_report <- function(family, bundle=NULL) {
  requirement <- nf_rna_font_requirement(family, if (is.null(bundle)) NULL else basename(bundle$dir))
  fail <- function(reason) list(requested=family, available=FALSE, family=NA_character_, file=NA_character_, sha256=NA_character_, variable=NA, requirement=requirement, reason=reason)
  if (!requireNamespace("systemfonts", quietly=TRUE)) return(fail("systemfonts unavailable; installed fonts could not be inspected"))
  fonts <- systemfonts::system_fonts()
  rows <- fonts[fonts$family == family, , drop=FALSE]
  if (!nrow(rows)) return(fail(paste0(family, " is not installed")))
  paths <- normalizePath(as.character(rows$path), mustWork=FALSE)
  if (!is.null(bundle)) {
    expected <- vapply(bundle$files[unlist(bundle$faces)], function(x) x$path, character(1))
    stray <- setdiff(paths, expected)
    if (length(stray)) return(fail(paste0("fontconfig lists ", family, " files that are not the bundled build: ", paste(stray, collapse=", "))))
    missing <- setdiff(expected, paths)
    if (length(missing)) return(fail(paste0("fontconfig does not list the bundled ", family, " files: ", paste(basename(missing), collapse=", "))))
  }
  regular <- rows[rows$style %in% c("Regular", "Book", "Roman"), , drop=FALSE]
  path <- as.character((if (nrow(regular)) regular else rows)$path[[1]])
  list(requested=family, available=TRUE, family=family, file=path, sha256=nf_rna_file_sha256(path), variable=grepl("\\[.*\\]", basename(path)), requirement=NULL, reason=NULL)
}

# Font names (subset tag removed) embedded in a PDF, read from raw and from
# Flate-inflated streams (cairo keeps font descriptors in object streams).
nf_rna_pdf_font_names <- function(path) {
  bytes <- readBin(path, "raw", file.size(path))
  chunks <- list(bytes)
  starts <- grepRaw("stream\r?\n", bytes, all=TRUE); ends <- grepRaw("\nendstream", bytes, all=TRUE)
  for (start in starts) {
    end <- ends[ends > start]
    if (!length(end)) next
    body <- bytes[(start + (if (bytes[start + 6L] == as.raw(13)) 8L else 7L)):(end[[1]] - 1L)]
    inflated <- tryCatch(memDecompress(body, type="gzip"), error=function(e) NULL)
    if (!is.null(inflated)) chunks[[length(chunks) + 1L]] <- inflated
  }
  names <- unlist(lapply(chunks, function(chunk) {
    hits <- grepRaw("/(BaseFont|FontName) ?/[A-Za-z0-9+._-]+", chunk, all=TRUE, value=TRUE)
    vapply(hits, function(hit) sub("^/(BaseFont|FontName) ?/", "", rawToChar(hit)), character(1))
  }))
  sort(unique(sub("^[A-Z]{6}\\+", "", names)))
}

nf_rna_font_key <- function(x) tolower(gsub("[^A-Za-z0-9]", "", x))

# Draw plain, bold, italic and bold-italic text with the family and read the
# embedded font names back.  `expected` (optional) lists the exact embedded names
# that must all be present.
nf_rna_verify_drawn_font <- function(family, expected=NULL) {
  probe <- tempfile(fileext=".pdf"); on.exit(unlink(probe), add=TRUE)
  grDevices::cairo_pdf(probe, width=3, height=1.5, family=family, onefile=FALSE)
  tryCatch({
    nf_rna_prime_device_metrics(family)
    grid::grid.newpage()
    for (face in 1:4) grid::grid.text("Hamburgefonstiv 0123", y=1 - face / 5, gp=grid::gpar(fontfamily=family, fontface=face))
  }, finally=grDevices::dev.off())
  embedded <- nf_rna_pdf_font_names(probe)
  want <- nf_rna_font_key(family)
  matched <- length(embedded) > 0 && all(startsWith(nf_rna_font_key(embedded), want))
  if (!is.null(expected)) matched <- matched && all(nf_rna_font_key(expected) %in% nf_rna_font_key(embedded))
  list(requested=family, embedded_fonts=embedded, matched=matched)
}

nf_rna_require_strict_font <- function(family, bundled=NULL) {
  bundle <- NULL
  if (!is.null(bundled)) { nf_rna_activate_bundled_font(bundled); bundle <- nf_rna_verify_bundled_font(bundled) }
  report <- nf_rna_strict_font_report(family, bundle)
  if (!report$available) stop(paste0("Required figure font '", family, "' is not available: ", report$reason, ". ", report$requirement), call.=FALSE)
  drawn <- nf_rna_verify_drawn_font(family, if (is.null(bundle)) NULL else unlist(bundle$manifest$expected_embedded_names))
  if (!drawn$matched) {
    stop(paste0("Figure font '", family, "' is listed as installed but the drawing device embedded [", paste(drawn$embedded_fonts, collapse=", "),
                "] instead (fontconfig substitution, or fontconfig was initialised before the bundled font was activated); refusing to draw with a different font. ", report$requirement), call.=FALSE)
  }
  out <- list(requested=family, family=family, fallback_used=FALSE, verified=TRUE, file=report$file, sha256=report$sha256, variable=report$variable,
              embedded_fonts=as.list(drawn$embedded_fonts), reason=NULL)
  if (!is.null(bundle)) {
    out$bundled <- list(package_dir=paste0("rnaseq/fonts/", bundled), version=bundle$manifest$version, license=bundle$manifest$license,
                        source=bundle$manifest$source$repository, source_commit=bundle$manifest$source$commit,
                        files=lapply(bundle$files, function(x) list(sha256=x$sha256)))
  }
  out
}

# The context describes the run: profile, requested preset and the legacy font used by
# figure types that are not migrated yet.  Migrated types draw through a scoped context
# (nf_rna_figure_scope) that carries the standard style and its strict font.
nf_rna_figure_context <- function(cfg) {
  requested <- nf_rna_read_figure_style(cfg)
  # fontconfig reads its configuration once per process, so the bundled font is activated
  # before anything (including the legacy resolver below) touches it.
  if (identical(requested$preset, "standard")) nf_rna_activate_bundled_font(NF_RNA_STYLE_PRESETS$standard$font$bundled)
  context <- new.env(parent=emptyenv())
  context$profile <- nf_rna_figure_profile(requested$profile)
  context$requested_preset <- requested$preset
  context$requested_name <- requested$given
  context$preset <- "legacy"                       # appearance of figures that are not scoped
  context$style <- NF_RNA_STYLE_PRESETS$legacy
  context$font <- nf_rna_resolve_font()
  context$standard_figures <- if (identical(requested$preset, "standard")) NF_RNA_STANDARD_FIGURES else character(0)
  context$scoped <- character(0)
  context$standard_font <- NULL
  context$records <- list()
  context
}

# The context to draw one figure type with.  Types migrated to the standard preset get a
# derived context (standard style, verified strict font); everything else, and every
# legacy run, returns the context unchanged.  Records from the derived context are
# collected on the run's context.
nf_rna_figure_scope <- function(context, type) {
  if (!(type %in% context$standard_figures)) return(context)
  if (is.null(context$standard_font)) {
    spec <- NF_RNA_STYLE_PRESETS$standard$font
    context$standard_font <- nf_rna_require_strict_font(spec$preferred, spec$bundled)
  }
  scoped <- new.env(parent=emptyenv())
  scoped$profile <- context$profile
  scoped$requested_preset <- context$requested_preset
  scoped$requested_name <- context$requested_name
  scoped$preset <- "standard"
  scoped$style <- NF_RNA_STYLE_PRESETS$standard
  scoped$font <- context$standard_font
  scoped$parent <- context
  scoped$figure_type <- type
  if (!(type %in% context$scoped)) context$scoped <- c(context$scoped, type)
  scoped
}

# Width in mm of a figure type under the context's preset; `default` when the preset sets no layout.
nf_rna_layout_width_mm <- function(context, type, default) {
  layout <- nf_rna_style(context)$layout[[type]]
  if (is.null(layout)) return(default)
  if (identical(layout$width, "single")) nf_rna_single_column_mm(context) else context$profile$double_col_mm
}
# Width of a "single column" figure.  A profile without a single column (PLOS: single_col_mm is NA) gets the preset's compact
# width, which never exceeds the profile's maximum; a preset without one falls back to the double column.
nf_rna_single_column_mm <- function(context) {
  profile <- context$profile
  if (!is.null(profile$single_col_mm) && !is.na(profile$single_col_mm)) return(profile$single_col_mm)
  compact <- nf_rna_style(context)$compact_width_mm
  if (is.null(compact)) profile$double_col_mm else min(compact, profile$max_width_mm)
}
nf_rna_layout_value <- function(context, type, key, default=NULL) nf_rna_style(context)$layout[[type]][[key]] %||% default
nf_rna_layout_aspect <- function(context, type, default=NA_real_) nf_rna_style(context)$layout[[type]]$aspect %||% default

# ggplot2 text geoms take mm by default; theme()/gpar() take pt.
nf_rna_pt_to_mm <- function(pt) pt * 25.4 / 72.27
nf_rna_text_size <- function(pt) {
  if ("size.unit" %in% names(formals(ggplot2::geom_text))) list(size=pt, size.unit="pt") else list(size=nf_rna_pt_to_mm(pt))
}
nf_rna_linewidth <- function(pt) pt / ggplot2::.pt

nf_rna_theme <- function(context) {
  profile <- context$profile
  ggplot2::theme_classic(base_size=profile$base_pt, base_family=context$font$family, base_line_size=nf_rna_linewidth(0.5), base_rect_size=nf_rna_linewidth(0.5)) +
    ggplot2::theme(
      text=ggplot2::element_text(size=profile$base_pt, family=context$font$family, colour="black"),
      axis.text=ggplot2::element_text(size=profile$small_pt, colour="black"),
      axis.title=ggplot2::element_text(size=profile$base_pt),
      legend.text=ggplot2::element_text(size=profile$small_pt),
      legend.title=ggplot2::element_text(size=profile$base_pt),
      plot.title=ggplot2::element_text(size=profile$base_pt, face="bold"),
      plot.subtitle=ggplot2::element_text(size=profile$base_pt, colour="black", lineheight=1.05, margin=ggplot2::margin(b=3)),
      plot.background=ggplot2::element_rect(fill="white", colour=NA),
      panel.background=ggplot2::element_rect(fill="white", colour=NA),
      legend.key.size=grid::unit(3.5, "mm")
    )
}

# ------------------------------------------------------------------ colour

# Legacy-valued aliases kept for existing callers and tests; the preset table is
# the single definition.
NF_RNA_OKABE_ITO <- NF_RNA_STYLE_PRESETS$legacy$categorical
NF_RNA_SHAPES <- NF_RNA_STYLE_PRESETS$legacy$shapes
NF_RNA_DIRECTION_COLOURS <- NF_RNA_STYLE_PRESETS$legacy$direction
NF_RNA_NA_COLOUR <- NF_RNA_STYLE_PRESETS$legacy$missing
NF_RNA_DIVERGING_ENDS <- NF_RNA_STYLE_PRESETS$legacy$diverging
NF_RNA_SEQUENTIAL_STOPS <- NF_RNA_STYLE_PRESETS$legacy$sequential

# `context` is optional everywhere below; NULL means the legacy preset.
nf_rna_style <- function(context=NULL) if (is.null(context)) NF_RNA_STYLE_PRESETS$legacy else context$style

# Colour-vision-deficiency simulation (Machado, Oliveira & Fernandes 2009, severity 1.0, applied in linear RGB).
NF_RNA_CVD_MATRICES <- list(
  deuteranopia=matrix(c(0.367322, 0.860646, -0.227968, 0.280085, 0.672501, 0.047413, -0.011820, 0.042940, 0.968881), 3, byrow=TRUE),
  protanopia=matrix(c(0.152286, 1.052583, -0.204868, 0.114503, 0.786281, 0.099216, -0.003882, -0.048116, 1.051998), 3, byrow=TRUE),
  tritanopia=matrix(c(1.255528, -0.076749, -0.178779, -0.078411, 0.930809, 0.147602, 0.004733, 0.691367, 0.303900), 3, byrow=TRUE))
nf_rna_cvd_lab <- function(hex, matrix=NULL) {
  rgb <- t(grDevices::col2rgb(hex)) / 255
  lin <- ifelse(rgb <= 0.04045, rgb / 12.92, ((rgb + 0.055) / 1.055)^2.4)
  if (!is.null(matrix)) lin <- pmin(pmax(lin %*% t(matrix), 0), 1)
  srgb <- ifelse(lin <= 0.0031308, lin * 12.92, 1.055 * lin^(1 / 2.4) - 0.055)
  grDevices::convertColor(pmin(pmax(srgb, 0), 1), from="sRGB", to="Lab")
}

# Deterministic extension of a categorical palette beyond its own slots: slots 1..length(palette) are kept and every further
# colour is chosen, one at a time, from a fixed HCL pool as the one whose smallest CIE Lab distance to all colours already
# chosen is largest, taking the worst of normal vision and simulated deuteranopia, protanopia and tritanopia (first maximum in
# pool order).  Pool colours must be mid-lightness and chromatic enough to be told apart from the neutral axis ink and have at
# least 2:1 contrast with white.  Nothing is random: the same n always gives the same colours, and a longer palette starts with
# the shorter one.
nf_rna_extend_palette <- function(palette, n) {
  palette <- unname(palette)
  if (n <= length(palette)) return(palette[seq_len(n)])
  pool <- unique(toupper(c(grDevices::hcl.colors(72, "Dark 3"), grDevices::hcl.colors(72, "Set 2"), grDevices::hcl.colors(72, "Dark 2"))))
  lin <- function(v) ifelse(v <= 0.04045, v / 12.92, ((v + 0.055) / 1.055)^2.4)
  rgb <- t(grDevices::col2rgb(pool)) / 255
  contrast <- 1.05 / (0.2126 * lin(rgb[, 1]) + 0.7152 * lin(rgb[, 2]) + 0.0722 * lin(rgb[, 3]) + 0.05)
  lab <- nf_rna_cvd_lab(pool); chroma <- sqrt(lab[, 2]^2 + lab[, 3]^2)
  pool <- pool[contrast >= 2 & chroma >= 25 & lab[, 1] >= 30 & lab[, 1] <= 80]
  views <- c(list(NULL), NF_RNA_CVD_MATRICES)
  chosen <- toupper(palette)
  while (length(chosen) < n) {
    candidates <- setdiff(pool, chosen)
    if (!length(candidates)) stop("the extension pool is exhausted; too many categorical levels")
    score <- rep(Inf, length(candidates))
    for (view in views) {
      anchor <- nf_rna_cvd_lab(chosen, view); cand <- nf_rna_cvd_lab(candidates, view)
      nearest <- vapply(seq_len(nrow(cand)), function(i) min(sqrt(rowSums((anchor - matrix(cand[i, ], nrow(anchor), 3, byrow=TRUE))^2))), numeric(1))
      score <- pmin(score, nearest)
    }
    chosen <- c(chosen, candidates[[which.max(score)]])
  }
  chosen
}

# `extend`: TRUE keeps the palette's own slots and extends deterministically (standard preset); FALSE is the historic behaviour
# (legacy): beyond the palette, every level is recoloured from HCL "Dark 3".
nf_rna_palette_values <- function(levels, palette, extend=FALSE) {
  levels <- unique(as.character(levels))
  if (anyNA(levels)) stop("categorical colour levels must not be NA")
  values <- if (length(levels) <= length(palette)) unname(palette[seq_along(levels)])
            else if (isTRUE(extend)) nf_rna_extend_palette(palette, length(levels))
            else grDevices::hcl.colors(length(levels), "Dark 3")
  stats::setNames(values, levels)
}
nf_rna_extends_palette <- function(context=NULL) !is.null(nf_rna_style(context)$categorical_extension_name)
nf_rna_categorical_colors <- function(levels, context=NULL) nf_rna_palette_values(levels, nf_rna_style(context)$categorical, nf_rna_extends_palette(context))
# Batch colours; a preset without an independent batch palette reuses the group palette.
nf_rna_batch_colors <- function(levels, context=NULL) {
  style <- nf_rna_style(context)
  nf_rna_palette_values(levels, style$batch %||% style$categorical, nf_rna_extends_palette(context))
}
# Name of the palette strategy that produced `n` colours (recorded in the manifests).  `batch`: the batch palette.
nf_rna_categorical_palette_name <- function(n, context=NULL, batch=FALSE) {
  style <- nf_rna_style(context)
  palette <- if (batch) style$batch %||% style$categorical else style$categorical
  if (n <= length(palette)) { if (batch && !is.null(style$batch_name)) style$batch_name else style$categorical_name }
  else if (nf_rna_extends_palette(context)) { if (batch && !is.null(style$batch_extension_name)) style$batch_extension_name else style$categorical_extension_name }
  else "HCL qualitative 'Dark 3'"
}
nf_rna_categorical_shapes <- function(levels, context=NULL) {
  levels <- unique(as.character(levels))
  shapes <- nf_rna_style(context)$shapes
  if (length(levels) > length(shapes)) return(NULL)
  stats::setNames(shapes[seq_along(levels)], levels)
}
nf_rna_direction_colours <- function(context=NULL) nf_rna_style(context)$direction
nf_rna_diverging_ends <- function(context=NULL) nf_rna_style(context)$diverging
nf_rna_sequential_stops <- function(context=NULL) nf_rna_style(context)$sequential
nf_rna_na_colour <- function(context=NULL) nf_rna_style(context)$missing
nf_rna_primary_colour <- function(context=NULL) nf_rna_style(context)$primary
# Continuous enrichment colours (most significant first) and the single mark colour of ORA dotplots.
nf_rna_enrichment_colours <- function(context=NULL) { style <- nf_rna_style(context); style$enrichment %||% rev(style$sequential[-1]) }
nf_rna_enrichment_mark <- function(context=NULL) { style <- nf_rna_style(context); if (is.null(style$enrichment)) style$primary else style$enrichment[[1]] }
nf_rna_diverging_palette <- function(n, context=NULL) grDevices::colorRampPalette(unname(nf_rna_diverging_ends(context)))(n)
nf_rna_sequential_palette <- function(n, context=NULL) grDevices::colorRampPalette(nf_rna_sequential_stops(context))(n)
# Point alpha is NA (ggplot2's "opaque") unless the preset sets one.
nf_rna_point_alpha <- function(context=NULL) nf_rna_style(context)$point_alpha
# PCA label size in pt: the preset's (6 pt), raised to the profile's smallest allowed text where that is larger.
nf_rna_pca_label_pt <- function(context) {
  value <- nf_rna_style(context)$pca_label_pt
  if (is.na(value)) context$profile$min_pt else max(value, context$profile$min_pt)  # never below the profile's smallest allowed text
}
# PCA y axis under equal axis scaling; the scores themselves are never changed.  PC2 is constant when it has
# no variance (two samples, or samples on one line); prcomp centres the scores, so that constant is 0 up to
# rounding.  A constant PC2 gets a fixed extent of +/-20% of the PC1 span and a single break at the constant:
# the former +/-1% sliver left room for none of its default breaks.  A varying PC2 keeps at least +/-1% of the
# PC1 span so the panel never has zero height.  `extent` is the y span the figure height is sized for.
nf_rna_pca_y_axis <- function(pc1, pc2) {
  x_range <- max(diff(range(pc1)), 1e-9)
  centre <- mean(range(pc2))
  tolerance <- sqrt(.Machine$double.eps) * x_range
  if (diff(range(pc2)) > tolerance) {
    return(list(constant = FALSE, limits = centre + c(-1, 1) * 0.01 * x_range, breaks = waiver(), extent = max(diff(range(pc2)), x_range * 0.02)))
  }
  if (abs(centre) <= tolerance) centre <- 0
  list(constant = TRUE, limits = centre + c(-1, 1) * 0.2 * x_range, breaks = centre, extent = 0.4 * x_range)
}

# Stable group/batch colours.  The colour of a level is decided by its position in the
# `universe` (the complete, ordered set of levels of the design variable, for example
# levels(metadata$condition)), never by the levels a particular figure happens to
# contain.  `present` selects the levels a figure needs; a subset, a reordering or a
# missing group therefore never moves a colour.  Group and batch use separate palettes.
nf_rna_stable_colours <- function(universe, present, palette, extend=FALSE) {
  universe <- unique(as.character(universe)); present <- as.character(present)
  if (anyNA(universe) || anyNA(present)) stop("categorical colour levels must not be NA")
  unknown <- setdiff(present, universe)
  if (length(unknown)) stop(paste0("levels not in the colour universe: ", paste(unknown, collapse=", ")))
  nf_rna_palette_values(universe, palette, extend)[present]
}
nf_rna_group_colours <- function(universe, present=universe, context=NULL) nf_rna_stable_colours(universe, present, nf_rna_style(context)$categorical, nf_rna_extends_palette(context))
nf_rna_batch_colours <- function(universe, present=universe, context=NULL) {
  style <- nf_rna_style(context)
  nf_rna_stable_colours(universe, present, style$batch %||% style$categorical, nf_rna_extends_palette(context))
}

# Limits for a sequential scale: fixed ends are kept; a degenerate range (for
# example an all-zero distance matrix) is widened by one unit on the free side.
nf_rna_sequential_limits <- function(values, lower=NULL, upper=NULL) {
  finite <- values[is.finite(values)]
  lo <- lower %||% (if (length(finite)) min(finite) else 0)
  hi <- upper %||% (if (length(finite)) max(finite) else lo + 1)
  degenerate <- !(hi > lo)
  if (degenerate) { if (!is.null(upper) && is.null(lower)) lo <- hi - 1 else hi <- lo + 1 }
  list(limits=c(lo, hi), degenerate=degenerate)
}

# ------------------------------------------------------------------ display-only transformations

# Row z-scores that keep missing values missing.  A row with fewer than two
# finite values or zero variance has no defined z-score and is returned as NA
# (never as 0); its names are listed in attr(, "undefined_rows").  For complete
# rows this equals t(scale(t(x))).
nf_rna_row_zscore <- function(x) {
  x <- as.matrix(x)
  storage.mode(x) <- "double"
  x[!is.finite(x)] <- NA_real_
  n_finite <- rowSums(!is.na(x))
  centre <- rowMeans(x, na.rm=TRUE)
  spread <- apply(x, 1, function(row) if (sum(!is.na(row)) >= 2) stats::sd(row, na.rm=TRUE) else NA_real_)
  z <- (x - centre) / spread
  defined <- n_finite >= 2 & is.finite(spread) & spread > 0
  z[!defined, ] <- NA_real_
  dimnames(z) <- dimnames(x)
  attr(z, "undefined_rows") <- rownames(x)[!defined] %||% character()
  z
}

# Colour saturation for drawing only: returns a clipped copy, NA preserved.
nf_rna_saturate <- function(x, limit) {
  y <- x
  y[!is.na(y) & y > limit] <- limit
  y[!is.na(y) & y < -limit] <- -limit
  attr(y, "undefined_rows") <- NULL
  y
}

# Complete-linkage Euclidean clustering (pheatmap's default), or NULL when
# missing values leave a pairwise distance undefined.
nf_rna_cluster <- function(m) {
  if (nrow(m) < 2) return(NULL)
  d <- stats::dist(m, method="euclidean")
  if (any(!is.finite(d))) return(NULL)
  stats::hclust(d, method="complete")
}

# Wrap long category labels; anything beyond max_lines is truncated with "...".
# The full text always stays in the figure's source-data table.
nf_rna_display_labels <- function(x, width=45, max_lines=2) {
  x <- as.character(x)
  wrapped <- lapply(x, function(text) strwrap(text, width=width))
  truncated <- vapply(wrapped, length, integer(1)) > max_lines
  labels <- vapply(seq_along(wrapped), function(i) {
    lines <- wrapped[[i]]
    if (truncated[[i]]) {
      lines <- lines[seq_len(max_lines)]
      lines[[max_lines]] <- paste0(substr(lines[[max_lines]], 1, max(1, width - 3)), "...")
    }
    paste(lines, collapse="\n")
  }, character(1))
  attr(labels, "truncated") <- truncated
  labels
}

# Compact heatmap title: "<numerator> vs <denominator>" over a detail line.
# pheatmap draws titles bold at 1.3 x fontsize; when the contrast line would
# not fit the figure width (conservative 0.7 em per character) it is broken
# before "vs", so group names are never truncated.
nf_rna_heatmap_title <- function(context, numerator, denominator, detail, width_mm) {
  title_pt <- 1.3 * context$profile$base_pt
  budget <- floor(width_mm / (nf_rna_pt_to_mm(title_pt) * 0.7))
  contrast <- sprintf("%s vs %s", numerator, denominator)
  if (nchar(contrast) > budget) contrast <- sprintf("%s\nvs %s", numerator, denominator)
  paste0(contrast, "\n", detail)
}

# pheatmap places its title over the matrix column only, with clipping off,
# so a title wider than the matrix runs past the page edge.  Span it across
# every column of the figure instead; the figure size is unchanged.
nf_rna_span_pheatmap_title <- function(heat) {
  layout <- heat$gtable$layout
  main <- layout$name == "main"
  if (any(main)) {
    heat$gtable$layout$l[main] <- 1L
    heat$gtable$layout$r[main] <- ncol(heat$gtable)
  }
  heat
}

# ------------------------------------------------------------------ content-aware sizing

# Guard against stale text metrics across graphics devices.  R reuses the metric it last computed for a given font
# specification when the next text height request names the same font, even if a different device (another resolution:
# PDF 72, PNG 150, TIFF 300 or 500 dpi) was opened in between.  The first text measured on the new device then comes
# back in the previous device's units (exactly 1/2, 1/4 or 2 times the right height), which makes row heights of the
# drawn layout wrong: the PDF or TIFF of a figure intermittently differs from the PNG layout, in the ordinary
# png -> tiff -> pdf sequence as well as after the sizing probe.  Measuring one throw-away string in a font that no
# figure uses, immediately after a device is opened and before anything is drawn or measured, replaces that remembered
# entry, so every real measurement is made on the device that is current.  It draws nothing and changes no figure.
# (Evidence and measurements: docs/figure-standard.md, "Intermittent PDF layout".)
nf_rna_prime_device_metrics <- function(family="") {
  invisible(grid::convertHeight(grid::grobHeight(grid::textGrob("M", gp=grid::gpar(fontsize=1.37, fontface=4L, fontfamily=family))), "pt", valueOnly=TRUE))
}

# Run `fun` on a throw-away Cairo PDF device that uses the context's font, so grid
# measures text with the family that will actually be drawn.
nf_rna_with_probe_device <- function(context, fun) {
  probe <- tempfile(fileext=".pdf")
  grDevices::cairo_pdf(probe, width=6, height=6, family=context$font$family, onefile=FALSE)
  on.exit({ grDevices::dev.off(); unlink(probe) }, add=TRUE)
  nf_rna_prime_device_metrics(context$font$family)
  fun()
}

# Absolute width and height (mm) of a gtable, ignoring the null/npc-sized panel or matrix
# cells named by `flexible` (prefix match on layout names).  Call inside nf_rna_with_probe_device.
nf_rna_gtable_extent_mm <- function(gtable, flexible=c("panel", "matrix")) {
  layout <- gtable$layout
  flex <- layout[Reduce(`|`, lapply(flexible, function(name) startsWith(layout$name, name)), FALSE), , drop=FALSE]
  skip_w <- unique(unlist(Map(seq, flex$l, flex$r))); skip_h <- unique(unlist(Map(seq, flex$t, flex$b)))
  size <- function(units, skip, convert) sum(vapply(setdiff(seq_along(units), skip), function(i) {
    u <- units[i]
    if (grid::unitType(u) == "null") 0 else convert(u, "mm", valueOnly=TRUE)
  }, numeric(1)))
  list(width=size(gtable$widths, skip_w, grid::convertWidth), height=size(gtable$heights, skip_h, grid::convertHeight))
}

# Width (mm) of a named grob of a gtable (for example a legend or a title), or 0.
nf_rna_grob_width_mm <- function(gtable, name) {
  index <- which(startsWith(gtable$layout$name, name))
  if (!length(index)) return(0)
  max(vapply(index, function(i) grid::convertWidth(grid::grobWidth(gtable$grobs[[i]]), "mm", valueOnly=TRUE), numeric(1)))
}

# Pheatmap cell size in bigpts for a size in mm.
nf_rna_mm_to_bigpts <- function(mm) mm * 72 / 25.4


# Fit absolute cell sizes (mm) for a pheatmap so that the figure follows its content.
# `build(cell_w_mm, cell_h_mm)` returns a silent pheatmap whose cellwidth/cellheight are the
# given sizes.  The gtable is then measured with the real font (labels, trees, legends, annotation)
# and the cells shrink, never below min_cell_mm, only when the profile maximum would be exceeded.
nf_rna_fit_cells <- function(context, spec, build, n_cols, n_rows, cw, ch, square=FALSE) {
  profile <- context$profile
  measure <- function(cw, ch) {
    heat <- build(cw, ch)
    m <- nf_rna_with_probe_device(context, function() {
      e <- nf_rna_gtable_extent_mm(heat$gtable, character(0))
      c(width=e$width, height=e$height, title=nf_rna_grob_width_mm(heat$gtable, "main"))
    })
    list(heat=heat, w=m[["width"]], h=m[["height"]], title=m[["title"]])
  }
  pad <- 2 * spec$pad_mm
  ch <- max(ch, spec$matrix_min_mm / n_rows)
  r <- measure(cw, ch)
  over_w <- r$w + pad - profile$max_width_mm; over_h <- r$h + pad - profile$max_height_mm
  if (over_w > 0 || over_h > 0) {
    new_cw <- if (over_w > 0) max(spec$min_cell_mm, cw - over_w / n_cols) else cw
    new_ch <- if (over_h > 0) max(spec$min_cell_mm, ch - over_h / n_rows) else ch
    if (square) new_cw <- new_ch <- min(new_cw, new_ch)
    cw <- new_cw; ch <- new_ch
    r <- measure(cw, ch)
  }
  c(r, list(cw=cw, ch=ch, width_mm=min(profile$max_width_mm, max(r$w, r$title + 4) + pad), height_mm=min(profile$max_height_mm, r$h + pad),
    clipped=(r$w + pad > profile$max_width_mm + 1e-6) || (r$h + pad > profile$max_height_mm + 1e-6)))
}

# Choose labels and cell sizes for a standard-preset heatmap.  `build(cw, ch, show_rows, show_cols, row_pt)`
# returns the pheatmap.  Row labels keep the largest font (down to the profile minimum) that stays
# legible at the fitted cell height; labels that cannot be legible are hidden and the cells are
# sized for an unlabelled matrix.  Nothing here changes a value, an order or a clustering.
nf_rna_fit_heatmap <- function(context, type, build, n_cols, n_rows, row_pt, col_pt, square=FALSE) {
  spec <- nf_rna_style(context)$layout[[type]]; profile <- context$profile
  col_need <- spec$label_fit * nf_rna_pt_to_mm(col_pt)
  attempt <- function(show_rows, show_cols, cw, ch, pt) {
    r <- nf_rna_fit_cells(context, spec, function(w, h) build(w, h, show_rows, show_cols, pt), n_cols, n_rows, cw, ch, square)
    r <- c(r, list(show_rows=show_rows, show_cols=show_cols, row_pt=pt))
    # Labels that still do not fit the profile maximum are hidden rather than clipped.
    if (r$clipped && show_cols) return(attempt(show_rows, FALSE, cw, ch, pt))
    r
  }
  if (square) {
    r <- attempt(TRUE, TRUE, spec$cell_mm, spec$cell_mm, row_pt)
    return(if (min(r$cw, r$ch) >= col_need - 1e-6) r else attempt(FALSE, FALSE, spec$cell_mm, spec$cell_mm, row_pt))
  }
  for (pt in seq(row_pt, max(profile$min_pt, 1), by=-0.5)) {
    if (pt > row_pt || pt < profile$min_pt) next
    ch_nominal <- spec$label_fit * nf_rna_pt_to_mm(pt)
    r <- attempt(TRUE, TRUE, spec$cell_mm, ch_nominal, pt)
    if (r$ch >= ch_nominal - 1e-6) return(if (r$cw >= col_need - 1e-6) r else attempt(TRUE, FALSE, spec$cell_mm, ch_nominal, pt))
  }
  r <- attempt(FALSE, TRUE, spec$cell_mm, spec$hidden_row_mm, row_pt)
  if (r$cw >= col_need - 1e-6) r else attempt(FALSE, FALSE, spec$cell_mm, spec$hidden_row_mm, row_pt)
}

# Standard-preset DEG heatmap: the matrix, annotation, clustering objects and title are given; this only
# chooses colours, typography and a content-aware size.  Returns the fit (heat, width_mm, height_mm, cell sizes, label state).
nf_rna_deg_heatmap <- function(context, display, annotation, annotation_colours, row_tree, column_tree, title, z_limit) {
  profile <- context$profile; layout <- nf_rna_style(context)$layout$heatmap
  legend_at <- seq(-z_limit, z_limit, by = 1)
  args <- list(mat = display, color = nf_rna_diverging_palette(100, context), breaks = seq(-z_limit, z_limit, length.out = 101), na_col = nf_rna_na_colour(context),
    cluster_rows = if (is.null(row_tree)) FALSE else row_tree, cluster_cols = if (is.null(column_tree)) FALSE else column_tree,
    annotation_col = annotation, annotation_colors = annotation_colours, border_color = NA,
    legend_breaks = legend_at, legend_labels = c(paste0("\u2264", -z_limit), legend_at[-c(1, length(legend_at))], paste0("\u2265", z_limit)),
    main = title, fontsize = profile$base_pt, fontsize_col = profile$small_pt, silent = TRUE)
  fit <- nf_rna_fit_heatmap(context, "heatmap", function(cw, ch, show_rows, show_cols, row_pt) do.call(pheatmap::pheatmap, c(args, list(
    show_rownames = show_rows, show_colnames = show_cols, fontsize_row = row_pt, cellwidth = nf_rna_mm_to_bigpts(cw), cellheight = nf_rna_mm_to_bigpts(ch),
    treeheight_row = nf_rna_mm_to_bigpts(layout$tree_mm), treeheight_col = nf_rna_mm_to_bigpts(layout$tree_mm)))),
    ncol(display), nrow(display), profile$small_pt, profile$small_pt)
  fit$heat <- nf_rna_span_pheatmap_title(fit$heat)
  fit
}

# Standard-preset sample correlation heatmap (square cells, group/batch annotation in the PCA colours).
nf_rna_correlation_heatmap <- function(context, corr, colours, breaks, annotation, annotation_colours, title) {
  profile <- context$profile; layout <- nf_rna_style(context)$layout$correlation
  fit <- nf_rna_fit_heatmap(context, "correlation", function(cw, ch, show_rows, show_cols, row_pt) {
    pheatmap::pheatmap(corr, color = colours, breaks = breaks, border_color = NA, main = title,
      annotation_col = annotation, annotation_colors = annotation_colours,
      cellwidth = nf_rna_mm_to_bigpts(cw), cellheight = nf_rna_mm_to_bigpts(ch),
      treeheight_row = nf_rna_mm_to_bigpts(layout$tree_mm), treeheight_col = nf_rna_mm_to_bigpts(layout$tree_mm),
      show_rownames = show_rows, show_colnames = show_cols,
      fontsize = profile$base_pt, fontsize_row = profile$small_pt, fontsize_col = profile$small_pt, silent = TRUE)
  }, ncol(corr), nrow(corr), profile$small_pt, profile$small_pt, square = TRUE)
  fit$heat <- nf_rna_span_pheatmap_title(fit$heat)
  fit
}

# Size of a standard-preset sample boxplot: width follows the number of samples, height follows the
# measured labels; the legend is wrapped into more rows until it fits the width.  `make(legend_rows)`
# returns the ggplot.  Returns the plot actually to be drawn and its size in mm.
nf_rna_fit_boxplot <- function(context, make, n_samples) {
  spec <- nf_rna_style(context)$layout$boxplot; profile <- context$profile
  measure <- function(p) nf_rna_with_probe_device(context, function() {
    g <- ggplot2::ggplotGrob(p)
    e <- nf_rna_gtable_extent_mm(g, "panel")
    c(width=e$width, height=e$height, legend=nf_rna_grob_width_mm(g, "guide-box"))
  })
  pad <- 2 * spec$pad_mm
  rows <- 1L; p <- make(rows); m <- measure(p)
  wanted <- max(spec$min_width_mm, m[["width"]] + n_samples * spec$sample_mm + pad)
  wanted <- min(wanted, profile$max_width_mm)
  while (m[["legend"]] + pad > wanted && rows < spec$max_legend_rows) { rows <- rows + 1L; p <- make(rows); m <- measure(p) }
  width <- min(profile$max_width_mm, max(wanted, m[["legend"]] + pad))
  height <- min(profile$max_height_mm, m[["height"]] + spec$panel_height_mm + pad)
  list(plot=p, width_mm=width, height_mm=height, legend_rows=rows, panel_width_mm=width - m[["width"]] - pad)
}

# Greedy word wrap of `text` to `width_mm` using the real font metrics (call inside nf_rna_with_probe_device).
nf_rna_wrap_measured <- function(context, text, width_mm, pt) {
  size <- function(x) grid::convertWidth(grid::grobWidth(grid::textGrob(x, gp=grid::gpar(fontfamily=context$font$family, fontsize=pt))), "mm", valueOnly=TRUE)
  out <- character()
  for (paragraph in strsplit(text, "\n", fixed=TRUE)[[1]]) {
    line <- ""
    for (word in strsplit(paragraph, " ", fixed=TRUE)[[1]]) {
      candidate <- if (nzchar(line)) paste(line, word) else word
      if (nzchar(line) && size(candidate) > width_mm) { out <- c(out, line); line <- word } else line <- candidate
    }
    out <- c(out, line)
  }
  paste(out, collapse="\n")
}

# Size of a standard-preset enrichment dotplot.  The y labels are measured (not estimated): the width is the
# measured label/axis/margin width plus the panel, the legend moves to a vertical stack if it is wider than the
# figure, the subtitle is wrapped to the figure width with the real font, and the height is the measured fixed
# height plus one row per label line and term.  `make(subtitle, legend_box, x_label)` returns the ggplot.
nf_rna_fit_dotplot <- function(context, make, labels, subtitle=NULL, x_label=NULL) {
  spec <- nf_rna_style(context)$layout$dotplot; profile <- context$profile
  pad <- 2 * spec$pad_mm
  measure <- function(p) nf_rna_with_probe_device(context, function() {
    g <- ggplot2::ggplotGrob(p)
    e <- nf_rna_gtable_extent_mm(g, "panel")
    c(width=e$width, height=e$height, legend=nf_rna_grob_width_mm(g, "guide-box"))
  })
  lines <- sum(lengths(strsplit(labels, "\n", fixed=TRUE)))
  box <- "horizontal"; p <- make(NULL, box, x_label); m <- measure(p)
  width <- min(profile$max_width_mm, max(spec$min_width_mm, m[["width"]] + spec$panel_mm + pad))
  if (m[["legend"]] + pad > width) { box <- "vertical"; p <- make(NULL, box, x_label); m <- measure(p); width <- min(profile$max_width_mm, max(width, m[["legend"]] + pad)) }
  wrapped <- if (is.null(subtitle)) NULL else nf_rna_with_probe_device(context, function() nf_rna_wrap_measured(context, gsub("\n", " ", subtitle, fixed=TRUE), width - pad - 2, profile$base_pt))
  # The x-axis title is centred under the panel: wrap it (same words) so it never runs past the figure edge.
  x_wrapped <- if (is.null(x_label)) NULL else nf_rna_with_probe_device(context, function() nf_rna_wrap_measured(context, x_label, max(20, width - m[["width"]] - pad), profile$base_pt))
  p <- make(wrapped, box, x_wrapped); m <- measure(p)
  panel_h <- max(spec$min_panel_height_mm, lines * nf_rna_pt_to_mm(profile$small_pt) * spec$line_factor + length(labels) * spec$row_gap_mm)
  list(plot=p, width_mm=width, height_mm=min(profile$max_height_mm, m[["height"]] + panel_h + pad), subtitle=wrapped,
       layout=list(label_axis_margin_mm=m[["width"]], legend_mm=m[["legend"]], legend_box=box, panel_height_mm=panel_h, label_lines=lines,
                   subtitle_lines=if (is.null(wrapped)) 0L else length(strsplit(wrapped, "\n", fixed=TRUE)[[1]])))
}

# ------------------------------------------------------------------ export

nf_rna_draw <- function(figure) {
  if (inherits(figure, "ggplot")) print(figure)
  else if (inherits(figure, "pheatmap")) { grid::grid.newpage(); grid::grid.draw(figure$gtable) }
  else if (inherits(figure, "gtable") || inherits(figure, "grob")) { grid::grid.newpage(); grid::grid.draw(figure) }
  else if (is.function(figure)) figure()
  else stop("unsupported figure object")
}

# Draw one figure directly at its final physical size on each device: a PNG
# preview for the HTML report, a raster TIFF and a Cairo PDF.  Nothing is
# resampled from another output.  Returns (and records) the export metadata.
nf_rna_save_figure <- function(context, figure, path_stem, width_mm, height_mm, output_root, meta=list(), formats=c("png", "tiff", "pdf")) {
  profile <- context$profile; family <- context$font$family
  if (width_mm > profile$max_width_mm + 1e-9 || height_mm > profile$max_height_mm + 1e-9) {
    warning(sprintf("%s is %.1f x %.1f mm, beyond the %s profile limit of %.1f x %.1f mm", basename(path_stem), width_mm, height_mm, profile$name, profile$max_width_mm, profile$max_height_mm))
  }
  # Cairo PDF pages are whole points (R truncates), so every output uses the
  # same size snapped to whole points and records exactly what was drawn.
  width_pt <- round(width_mm / 25.4 * 72); height_pt <- round(height_mm / 25.4 * 72)
  width_mm <- width_pt * 25.4 / 72; height_mm <- height_pt * 25.4 / 72
  files <- list()
  for (format in formats) {
    path <- paste0(path_stem, ".", format)
    dpi <- if (format == "png") profile$preview_dpi else profile$raster_dpi
    pixels <- c(round(width_pt / 72 * dpi), round(height_pt / 72 * dpi))
    if (format == "png") {
      grDevices::png(path, width=pixels[[1]], height=pixels[[2]], units="px", res=dpi, type="cairo", bg="white", family=family)
    } else if (format == "tiff") {
      grDevices::tiff(path, width=pixels[[1]], height=pixels[[2]], units="px", res=dpi, type="cairo", compression="lzw", bg="white", family=family)
    } else if (format == "pdf") {
      grDevices::cairo_pdf(path, width=(width_pt + 1e-4) / 72, height=(height_pt + 1e-4) / 72, family=family, bg="white", onefile=FALSE)
    } else stop(paste("unsupported figure format:", format))
    tryCatch({ nf_rna_prime_device_metrics(family); nf_rna_draw(figure) }, finally=grDevices::dev.off())
    entry <- list(path=nf_rna_relative_path(path, output_root), format=format, width_mm=width_mm, height_mm=height_mm)
    if (format %in% c("png", "tiff")) {
      entry <- c(entry, list(dpi=dpi, width_px=pixels[[1]], height_px=pixels[[2]], role=if (format == "png") "report preview" else "raster submission",
                             colour_model="RGB, opaque white background", compression=if (format == "tiff") "LZW" else "PNG deflate"))
    } else {
      entry <- c(entry, list(width_pt=width_pt, height_pt=height_pt, role="vector submission", note="Cairo PDF: text and lines are vector; a .pdf suffix alone does not guarantee all content is vector"))
    }
    files[[length(files) + 1]] <- entry
  }
  record <- c(list(id=nf_rna_relative_path(path_stem, output_root), width_mm=width_mm, height_mm=height_mm, files=files), meta)
  if (!is.null(context$parent)) {
    # drawn through a scoped (standard preset) context: say so, with the font actually embedded
    record$figure_style <- list(preset=context$preset, profile=profile$name, font=family, font_sha256=context$font$sha256, font_embedded=context$font$embedded_fonts)
  }
  target <- context$parent %||% context  # records of a scoped context are collected on the run's context
  target$records[[length(target$records) + 1]] <- record
  invisible(record)
}

nf_rna_relative_path <- function(path, root) {
  path <- normalizePath(path, mustWork=FALSE); root <- normalizePath(root, mustWork=FALSE)
  if (startsWith(path, paste0(root, "/"))) substring(path, nchar(root) + 2) else basename(path)
}

nf_rna_source_data <- function(path, output_root, description) {
  list(path=nf_rna_relative_path(path, output_root), md5=if (file.exists(path)) unname(tools::md5sum(path)) else NA_character_, description=description)
}

nf_rna_write_figure_manifest <- function(context, output_dir) {
  packages <- c("ggplot2", "ggrepel", "pheatmap", "scales", "systemfonts")
  versions <- lapply(packages, function(p) if (requireNamespace(p, quietly=TRUE)) as.character(utils::packageVersion(p)) else NA_character_)
  names(versions) <- packages
  software <- grDevices::grSoftVersion()
  # What was actually drawn: when the standard preset was requested and its figures were drawn, the manifest-level font and
  # colour policy describe the standard preset (its strict Arimo font and its palettes), not the legacy resolver's font.
  applied_standard <- !identical(context$requested_preset, "legacy") && length(context$scoped) > 0 && !is.null(context$standard_font)
  style <- if (applied_standard) NF_RNA_STYLE_PRESETS$standard else nf_rna_style(context)
  categorical_text <- paste0("deterministic named mapping in sorted level order; ", style$categorical_name, " for up to ", length(style$categorical), " levels, otherwise ",
                             if (!is.null(style$categorical_extension_name)) paste0(style$categorical_extension_name, " (slots 1-", length(style$categorical), " unchanged)") else "HCL qualitative 'Dark 3'",
                             "; group identity is also encoded by shape where available")
  colour_policy <- list(
    categorical=categorical_text,
    direction=as.list(style$direction),
    diverging=list(colours=as.list(style$diverging), use="only for values with a meaningful zero midpoint (row z-score; correlation when negative values occur)"),
    sequential=list(colours=as.list(style$sequential), use="magnitudes without a signed midpoint"),
    missing=style$missing
  )
  if (applied_standard) {
    colour_policy <- append(colour_policy, list(
      categorical_colours=as.list(style$categorical),
      batch=list(colours=as.list(style$batch), name=style$batch_name, use="batch variable colours; independent of the group palette"),
      enrichment=list(colours=as.list(style$enrichment), name="nf-rna Standard Enrichment Scale v1", use="continuous enrichment scale (adjusted P; ORA mark colour is the first stop)"),
      applied_preset="standard"), after=1)
  }
  document <- list(
    schema_version="nf-rna.figure-manifest.v1",
    profile=context$profile,
    font=if (applied_standard) context$standard_font else context$font,
    typography="theme and pheatmap sizes in pt; ggplot2 text geoms use size.unit='pt' (or pt converted to mm)",
    colour_policy=colour_policy,
    devices=list(png="grDevices::png(type='cairo')", tiff="grDevices::tiff(type='cairo', compression='lzw')", pdf="grDevices::cairo_pdf",
                 cairo=unname(software[["cairo"]] %||% NA_character_), libtiff=unname(software[["libtiff"]] %||% NA_character_)),
    r_version=R.version.string,
    package_versions=versions,
    figures=context$records
  )
  # Present only when a non-legacy preset was requested, so legacy manifests are unchanged.
  if (!identical(context$requested_preset, "legacy")) {
    others <- "legacy appearance and legacy font until migrated"
    document <- append(document, list(figure_style=list(
      requested_preset=context$requested_name, preset=context$requested_preset, profile=context$profile$name,
      standard_figure_types=as.list(context$standard_figures), standard_figures_drawn=as.list(context$scoped),
      other_figure_types=others, font=context$standard_font)), after=3)
  }
  write(jsonlite::toJSON(document, auto_unbox=TRUE, pretty=TRUE, null="null", na="null", digits=NA), file.path(output_dir, "figure_manifest.json"))
}

# ------------------------------------------------------------------ sample distribution boxplots

# Per-sample Tukey box statistics of a genes x samples matrix.  Columns are
# matched to `samples` by name and the result follows `samples` exactly; a
# missing or duplicated column is an error, never a silent reorder or drop.
nf_rna_sample_box_stats <- function(values, samples, matrix_label) {
  if (anyDuplicated(colnames(values)) || !all(samples %in% colnames(values))) {
    stop(paste0("sample distribution matrix '", matrix_label, "' does not contain exactly the expected samples"))
  }
  values <- values[, samples, drop = FALSE]
  rows <- lapply(seq_along(samples), function(i) {
    x <- as.numeric(values[, i]); x <- x[is.finite(x)]
    s <- if (length(x)) grDevices::boxplot.stats(x, coef = 1.5)$stats else rep(NA_real_, 5)
    data.frame(sample_id = samples[[i]], matrix = matrix_label, genes = length(x), whisker_low = s[[1]], q1 = s[[2]], median = s[[3]], q3 = s[[4]], whisker_high = s[[5]], stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}

# Metadata values for `samples`, looked up by sample_id.  Absent samples or
# missing values stop the plot instead of producing an NA/shifted annotation.
nf_rna_sample_annotation <- function(metadata, samples, variable) {
  if (!(variable %in% colnames(metadata))) stop(paste("annotation variable is absent from metadata:", variable))
  if (!all(samples %in% rownames(metadata))) stop("metadata does not contain every sample being plotted")
  values <- metadata[samples, variable]
  if (anyNA(values)) stop(paste("annotation variable has missing values for plotted samples:", variable))
  factor(as.character(values), levels = if (is.factor(metadata[[variable]])) levels(metadata[[variable]]) else sort(unique(as.character(values))))
}

# One box per sample, in the order of `stats$sample_id`.  Fill = group; an
# optional black marker row below the boxes = batch (shape).  Outliers are not
# drawn: with thousands of genes per sample they carry no per-sample message.
nf_rna_sample_boxplot <- function(context, stats, group, group_name, y_label, batch = NULL, batch_name = NULL) {
  profile <- context$profile
  samples <- stats$sample_id
  standard <- identical(context$preset, "standard")
  stats$sample_id <- factor(samples, levels = samples)
  stats$group <- group
  colours <- nf_rna_group_colours(levels(group), context = context)
  shapes <- NULL
  if (!is.null(batch)) shapes <- nf_rna_categorical_shapes(levels(batch), context)
  show_labels <- length(samples) <= 100
  rotation <- if (length(samples) > 8) 90 else 45
  make <- function(legend_rows = NULL) {
    p <- ggplot2::ggplot(stats, ggplot2::aes(x = sample_id)) +
      ggplot2::geom_boxplot(ggplot2::aes(ymin = whisker_low, lower = q1, middle = median, upper = q3, ymax = whisker_high, fill = group),
                            stat = "identity", width = 0.7, linewidth = nf_rna_linewidth(0.5), colour = "black") +
      ggplot2::scale_fill_manual(values = colours, name = group_name, drop = FALSE)
    if (!is.null(shapes)) {
      span <- max(stats$whisker_high, na.rm = TRUE) - min(stats$whisker_low, na.rm = TRUE)
      stats$batch <- batch; stats$batch_y <- min(stats$whisker_low, na.rm = TRUE) - 0.08 * max(span, 1e-9)
      p <- p + ggplot2::geom_point(data = stats, ggplot2::aes(x = sample_id, y = batch_y, shape = batch), size = 1.8, colour = "black", inherit.aes = FALSE) +
        ggplot2::scale_shape_manual(values = shapes, name = batch_name, drop = FALSE)
    }
    p <- p + ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = 0.05)) +
      ggplot2::labs(x = NULL, y = y_label) + nf_rna_theme(context) +
      ggplot2::theme(legend.position = "bottom", axis.text.x = if (show_labels) ggplot2::element_text(angle = rotation, hjust = 1, vjust = if (rotation == 90) 0.5 else 1) else ggplot2::element_blank(),
                     axis.ticks.x = if (show_labels) ggplot2::element_line() else ggplot2::element_blank())
    if (standard) {
      p <- p + ggplot2::theme(legend.box = "vertical", legend.box.just = "left", legend.margin = ggplot2::margin(0, 0, 0, 0), legend.box.spacing = grid::unit(1.5, "mm"))
      if (!is.null(legend_rows)) {
        p <- p + ggplot2::guides(fill = ggplot2::guide_legend(nrow = legend_rows, byrow = TRUE, order = 1))
        if (!is.null(shapes)) p <- p + ggplot2::guides(shape = ggplot2::guide_legend(nrow = legend_rows, byrow = TRUE, order = 2))
      }
    }
    p
  }
  result <- list(colours = colours, shapes = shapes, labels_shown = show_labels, rotation = if (show_labels) rotation else NA_real_)
  if (standard) {
    fit <- nf_rna_fit_boxplot(context, make, length(samples))
    return(c(list(plot = fit$plot), result, list(width_mm = fit$width_mm, height_mm = fit$height_mm, legend_rows = fit$legend_rows)))
  }
  label_mm <- if (show_labels) nf_rna_pt_to_mm(profile$small_pt) * 0.6 * max(nchar(samples)) else 0
  c(list(plot = make()), result, list(height_mm = min(profile$max_height_mm, 65 + label_mm)))
}

# ------------------------------------------------------------------ enrichment dotplots

nf_rna_dotplot_height <- function(context, labels, subtitle_lines=0) {
  lines <- sum(lengths(strsplit(labels, "\n", fixed=TRUE)))
  body <- max(60, 25 + lines * nf_rna_pt_to_mm(context$profile$small_pt) * 1.6 + length(labels) * 2)
  min(context$profile$max_height_mm, body + subtitle_lines * nf_rna_pt_to_mm(context$profile$base_pt) * 1.5)
}

# Break text to at most `width` characters per line; a single token longer than
# the width is broken as well, so no group name can run past the figure edge.
nf_rna_wrap_text <- function(text, width) {
  width <- max(10L, as.integer(width))
  unlist(lapply(strsplit(text, "\n", fixed=TRUE)[[1]], function(line) {
    pieces <- strwrap(line, width=width)
    if (!length(pieces)) return("")
    unlist(lapply(pieces, function(piece) if (nchar(piece) <= width) piece else regmatches(piece, gregexpr(sprintf(".{1,%d}", width), piece))[[1]]))
  }))
}

NF_RNA_ORA_FOREGROUNDS <- c(significant="Significant DE genes", up="Upregulated genes", down="Downregulated genes")

# Caption that makes a standalone ORA dotplot self-identifying.  The comparison
# and its direction come from the contrast's numerator/denominator (never from
# the contrast ID); "up" is higher in the numerator, "down" in the denominator.
nf_rna_ora_subtitle <- function(context, contrast, foreground, width_mm=context$profile$double_col_mm, max_lines=4) {
  numerator <- contrast$numerator; denominator <- contrast$denominator
  if (!is.character(numerator) || !is.character(denominator) || !nzchar(numerator) || !nzchar(denominator)) {
    stop("ORA figure subtitle requires the contrast numerator and denominator")
  }
  label <- NF_RNA_ORA_FOREGROUNDS[foreground]
  if (length(foreground) != 1 || is.na(label)) stop(paste("unknown ORA foreground:", foreground))
  direction <- switch(foreground, up=sprintf(" (higher in %s)", numerator), down=sprintf(" (higher in %s)", denominator), "")
  budget <- floor(0.8 * width_mm / (nf_rna_pt_to_mm(context$profile$base_pt) * 0.6))
  comparison <- sprintf("%s vs %s", numerator, denominator)
  if (nchar(comparison) > budget) comparison <- sprintf("%s\nvs %s", numerator, denominator)
  lines <- nf_rna_wrap_text(paste0(comparison, " \u2014 ", unname(label), direction), budget)
  if (length(lines) > max_lines) {
    lines <- lines[seq_len(max_lines)]
    lines[[max_lines]] <- paste0(substr(lines[[max_lines]], 1, max(1, budget - 3)), "...")
  }
  paste(lines, collapse="\n")
}

nf_rna_write_source_table <- function(table, path) write.table(table, path, sep="\t", quote=FALSE, row.names=FALSE, na="NA")

# Over-representation dotplot of an already selected, ordered term subset.
nf_rna_ora_dotplot <- function(context, top, x_label, path_stem, output_root, meta=list(), subtitle=NULL) {
  context <- nf_rna_figure_scope(context, "dotplot")  # standard preset when requested, else the context unchanged
  standard <- identical(context$preset, "standard")
  labels <- nf_rna_display_labels(top$Description)
  source <- data.frame(display_rank=seq_len(nrow(top)), ID=top$ID, Description=top$Description, display_label=gsub("\n", " ", labels, fixed=TRUE),
                       label_truncated=attr(labels, "truncated"), top[, intersect(c("GeneRatio", "BgRatio", "pvalue", "p.adjust", "qvalue", "Count"), names(top)), drop=FALSE],
                       check.names=FALSE, stringsAsFactors=FALSE)
  source_path <- paste0(path_stem, "_source_data.tsv")
  nf_rna_write_source_table(source, source_path)
  top$axis <- factor(top$ID, levels=rev(top$ID))
  # Display floor as for the volcano: an adjusted P of 0 is drawn at the
  # smallest positive double and marked; the tables keep the reported value.
  top$capped <- !is.na(top$p.adjust) & top$p.adjust < .Machine$double.xmin
  top$plot_padj <- pmax(top$p.adjust, .Machine$double.xmin)
  make <- function(subtitle_text=subtitle, legend_box="horizontal", x_text=x_label) {
    p <- ggplot2::ggplot(top, ggplot2::aes(x=-log10(plot_padj), y=axis, size=Count, shape=capped)) +
      ggplot2::geom_point(colour=nf_rna_enrichment_mark(context)) +
      ggplot2::scale_shape_manual(values=c("FALSE"=16, "TRUE"=17), breaks="TRUE", labels="adjusted P = 0, drawn at display floor", name=NULL, guide=if (any(top$capped)) "legend" else "none") +
      ggplot2::scale_y_discrete(labels=stats::setNames(labels, top$ID)) +
      ggplot2::scale_size_area(max_size=5, name="Genes") + ggplot2::expand_limits(x=0) +
      ggplot2::labs(x=x_text, y=NULL, subtitle=subtitle_text) + nf_rna_theme(context)
    if (!is.null(subtitle_text)) p <- p + ggplot2::theme(plot.title.position="plot")
    if (standard) p <- p + ggplot2::theme(legend.position="bottom", legend.location="plot", legend.justification="left", legend.box=legend_box, legend.box.just="left", legend.margin=ggplot2::margin(0, 0, 0, 0), legend.box.spacing=grid::unit(1.5, "mm"))
    p
  }
  if (standard) {
    fit <- nf_rna_fit_dotplot(context, make, labels, subtitle, x_label)
    p <- fit$plot; subtitle <- fit$subtitle; width_mm <- fit$width_mm; height_mm <- fit$height_mm
  } else {
    p <- make()
    subtitle_lines <- if (is.null(subtitle)) 0 else length(strsplit(subtitle, "\n", fixed=TRUE)[[1]])
    width_mm <- context$profile$double_col_mm; height_mm <- nf_rna_dotplot_height(context, labels, subtitle_lines)
  }
  nf_rna_save_figure(context, p, path_stem, width_mm, height_mm, output_root, c(list(
    kind="ORA dotplot", subtitle=subtitle, x="-log10(adjusted P) as reported by clusterProfiler", size="number of foreground genes in the term (Count)",
    display=list(padj_zero_floor=.Machine$double.xmin, capped_points=sum(top$capped)),
    labels=sprintf("term descriptions wrapped to 45 characters, at most 2 lines; %d truncated; full ID/Description in the source-data table", sum(attr(labels, "truncated"))),
    source_data=list(nf_rna_source_data(source_path, output_root, "plotted terms with full IDs and descriptions"))),
    if (standard) list(layout=fit$layout, term_order=as.list(top$ID), mark_colour=nf_rna_enrichment_mark(context)), meta))
}

# Preranked GSEA dotplot.  The colour scale for adjusted P has fixed limits
# [0, padj_cutoff] so that plots with the same cutoff are directly comparable.
nf_rna_gsea_dotplot <- function(context, top, contrast, padj_cutoff, p_adjust_method, path_stem, output_root, meta=list()) {
  context <- nf_rna_figure_scope(context, "dotplot")  # standard preset when requested, else the context unchanged
  standard <- identical(context$preset, "standard")
  top <- top[order(-top$NES, top$ID), , drop=FALSE]
  labels <- nf_rna_display_labels(top$Description)
  source <- data.frame(display_rank=seq_len(nrow(top)), ID=top$ID, Description=top$Description, display_label=gsub("\n", " ", labels, fixed=TRUE),
                       label_truncated=attr(labels, "truncated"), top[, intersect(c("setSize", "enrichmentScore", "NES", "pvalue", "p.adjust", "qvalue", "rank", "leading_edge"), names(top)), drop=FALSE],
                       check.names=FALSE, stringsAsFactors=FALSE)
  source_path <- paste0(path_stem, "_source_data.tsv")
  nf_rna_write_source_table(source, source_path)
  numerator <- contrast$numerator %||% NULL; denominator <- contrast$denominator %||% NULL
  x_label <- if (!is.null(numerator) && !is.null(denominator)) sprintf("NES (> 0: higher in %s; < 0: higher in %s)", numerator, denominator) else "NES (> 0: enriched among positively ranked genes)"
  top$axis <- factor(top$ID, levels=rev(top$ID))
  make <- function(subtitle_text=NULL, legend_box="horizontal", x_text=x_label) {
    p <- ggplot2::ggplot(top, ggplot2::aes(x=NES, y=axis, size=setSize, colour=p.adjust)) +
      ggplot2::geom_vline(xintercept=0, linewidth=nf_rna_linewidth(0.5), colour="grey40") +
      ggplot2::geom_point() +
      ggplot2::scale_y_discrete(labels=stats::setNames(labels, top$ID)) +
      ggplot2::scale_colour_gradientn(colours=nf_rna_enrichment_colours(context), limits=c(0, as.numeric(padj_cutoff)), oob=scales::squish, name=sprintf("Adjusted P (%s)", p_adjust_method)) +
      ggplot2::scale_size_area(max_size=5, name="Gene set size") +
      ggplot2::labs(x=x_text, y=NULL) + nf_rna_theme(context)
    if (standard) p <- p + ggplot2::theme(legend.position="bottom", legend.location="plot", legend.justification="left", legend.box=legend_box, legend.box.just="left", legend.margin=ggplot2::margin(0, 0, 0, 0), legend.box.spacing=grid::unit(1.5, "mm")) + ggplot2::guides(colour=ggplot2::guide_colourbar(theme=ggplot2::theme(legend.key.width=grid::unit(40, "mm"))))
    p
  }
  if (standard) {
    fit <- nf_rna_fit_dotplot(context, make, labels, NULL, x_label)
    p <- fit$plot; width_mm <- fit$width_mm; height_mm <- fit$height_mm
  } else {
    p <- make(); width_mm <- context$profile$double_col_mm; height_mm <- nf_rna_dotplot_height(context, labels)
  }
  nf_rna_save_figure(context, p, path_stem, width_mm, height_mm, output_root, c(list(
    kind="GSEA dotplot", contrast=contrast[intersect(c("contrast_id", "numerator", "denominator"), names(contrast))],
    x="normalized enrichment score; sign gives the direction", rank_metric="DESeq2 Wald statistic (stat) for numerator vs denominator; ties ordered by Entrez ID",
    colour=list(value=sprintf("adjusted P (%s)", p_adjust_method), limits=c(0, as.numeric(padj_cutoff)), kind="sequential", harmonized="fixed limits shared by every GSEA dotplot with this cutoff"),
    size="gene set size (setSize)", display_order="NES descending",
    labels=sprintf("term descriptions wrapped to 45 characters, at most 2 lines; %d truncated; full ID/Description in the source-data table", sum(attr(labels, "truncated"))),
    source_data=list(nf_rna_source_data(source_path, output_root, "plotted terms with full IDs, NES, adjusted P and leading edge"))),
    if (standard) list(layout=fit$layout, term_order=as.list(top$ID), enrichment_colours=as.list(nf_rna_enrichment_colours(context))), meta))
}
