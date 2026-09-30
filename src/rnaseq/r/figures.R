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

nf_rna_figure_context <- function(cfg) {
  figures <- if (is.list(cfg$figures)) cfg$figures else list()
  context <- new.env(parent=emptyenv())
  context$profile <- nf_rna_figure_profile(figures$profile)
  context$font <- nf_rna_resolve_font()
  context$records <- list()
  context
}

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
      plot.background=ggplot2::element_rect(fill="white", colour=NA),
      panel.background=ggplot2::element_rect(fill="white", colour=NA),
      legend.key.size=grid::unit(3.5, "mm")
    )
}

# ------------------------------------------------------------------ colour

# Okabe & Ito (2008) colour-universal-design palette, ordered for contrast on white.
NF_RNA_OKABE_ITO <- c(blue="#0072B2", vermillion="#D55E00", bluish_green="#009E73", reddish_purple="#CC79A7",
                      orange="#E69F00", sky_blue="#56B4E9", yellow="#F0E442", black="#000000")
NF_RNA_SHAPES <- c(16, 17, 15, 18, 1, 2, 0, 5, 3, 4, 8, 6)
NF_RNA_DIRECTION_COLOURS <- c(Up="#D55E00", Down="#0072B2", "Not significant"="#BDBDBD")
NF_RNA_NA_COLOUR <- "#A0A0A0"
NF_RNA_DIVERGING_ENDS <- c(low="#2166AC", mid="#F7F7F7", high="#B2182B")  # ColorBrewer RdBu
NF_RNA_SEQUENTIAL_STOPS <- c("#F7FBFF", "#C6DBEF", "#6BAED6", "#2171B5", "#08306B")  # ColorBrewer Blues

nf_rna_categorical_colors <- function(levels) {
  levels <- unique(as.character(levels))
  if (anyNA(levels)) stop("categorical colour levels must not be NA")
  values <- if (length(levels) <= length(NF_RNA_OKABE_ITO)) unname(NF_RNA_OKABE_ITO[seq_along(levels)]) else grDevices::hcl.colors(length(levels), "Dark 3")
  stats::setNames(values, levels)
}
nf_rna_categorical_palette_name <- function(n) if (n <= length(NF_RNA_OKABE_ITO)) "Okabe-Ito" else "HCL qualitative 'Dark 3'"
nf_rna_categorical_shapes <- function(levels) {
  levels <- unique(as.character(levels))
  if (length(levels) > length(NF_RNA_SHAPES)) return(NULL)
  stats::setNames(NF_RNA_SHAPES[seq_along(levels)], levels)
}
nf_rna_diverging_palette <- function(n) grDevices::colorRampPalette(unname(NF_RNA_DIVERGING_ENDS))(n)
nf_rna_sequential_palette <- function(n) grDevices::colorRampPalette(NF_RNA_SEQUENTIAL_STOPS)(n)

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
    tryCatch(nf_rna_draw(figure), finally=grDevices::dev.off())
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
  context$records[[length(context$records) + 1]] <- record
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
  packages <- c("ggplot2", "pheatmap", "scales", "systemfonts")
  versions <- lapply(packages, function(p) if (requireNamespace(p, quietly=TRUE)) as.character(utils::packageVersion(p)) else NA_character_)
  names(versions) <- packages
  software <- grDevices::grSoftVersion()
  document <- list(
    schema_version="nf-rna.figure-manifest.v1",
    profile=context$profile,
    font=context$font,
    typography="theme and pheatmap sizes in pt; ggplot2 text geoms use size.unit='pt' (or pt converted to mm)",
    colour_policy=list(
      categorical=paste0("deterministic named mapping in sorted level order; Okabe-Ito for up to ", length(NF_RNA_OKABE_ITO), " levels, otherwise HCL qualitative 'Dark 3'; group identity is also encoded by shape where available"),
      direction=as.list(NF_RNA_DIRECTION_COLOURS),
      diverging=list(colours=as.list(NF_RNA_DIVERGING_ENDS), use="only for values with a meaningful zero midpoint (row z-score; correlation when negative values occur)"),
      sequential=list(colours=as.list(NF_RNA_SEQUENTIAL_STOPS), use="magnitudes without a signed midpoint"),
      missing=NF_RNA_NA_COLOUR
    ),
    devices=list(png="grDevices::png(type='cairo')", tiff="grDevices::tiff(type='cairo', compression='lzw')", pdf="grDevices::cairo_pdf",
                 cairo=unname(software[["cairo"]] %||% NA_character_), libtiff=unname(software[["libtiff"]] %||% NA_character_)),
    r_version=R.version.string,
    package_versions=versions,
    figures=context$records
  )
  write(jsonlite::toJSON(document, auto_unbox=TRUE, pretty=TRUE, null="null", na="null", digits=NA), file.path(output_dir, "figure_manifest.json"))
}

# ------------------------------------------------------------------ enrichment dotplots

nf_rna_dotplot_height <- function(context, labels) {
  lines <- sum(lengths(strsplit(labels, "\n", fixed=TRUE)))
  min(context$profile$max_height_mm, max(60, 25 + lines * nf_rna_pt_to_mm(context$profile$small_pt) * 1.6 + length(labels) * 2))
}

nf_rna_write_source_table <- function(table, path) write.table(table, path, sep="\t", quote=FALSE, row.names=FALSE, na="NA")

# Over-representation dotplot of an already selected, ordered term subset.
nf_rna_ora_dotplot <- function(context, top, x_label, path_stem, output_root, meta=list()) {
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
  p <- ggplot2::ggplot(top, ggplot2::aes(x=-log10(plot_padj), y=axis, size=Count, shape=capped)) +
    ggplot2::geom_point(colour=NF_RNA_OKABE_ITO[["blue"]]) +
    ggplot2::scale_shape_manual(values=c("FALSE"=16, "TRUE"=17), breaks="TRUE", labels="adjusted P = 0, drawn at display floor", name=NULL, guide=if (any(top$capped)) "legend" else "none") +
    ggplot2::scale_y_discrete(labels=stats::setNames(labels, top$ID)) +
    ggplot2::scale_size_area(max_size=5, name="Genes") + ggplot2::expand_limits(x=0) +
    ggplot2::labs(x=x_label, y=NULL) + nf_rna_theme(context)
  nf_rna_save_figure(context, p, path_stem, context$profile$double_col_mm, nf_rna_dotplot_height(context, labels), output_root, c(list(
    kind="ORA dotplot", x="-log10(adjusted P) as reported by clusterProfiler", size="number of foreground genes in the term (Count)",
    display=list(padj_zero_floor=.Machine$double.xmin, capped_points=sum(top$capped)),
    labels=sprintf("term descriptions wrapped to 45 characters, at most 2 lines; %d truncated; full ID/Description in the source-data table", sum(attr(labels, "truncated"))),
    source_data=list(nf_rna_source_data(source_path, output_root, "plotted terms with full IDs and descriptions"))), meta))
}

# Preranked GSEA dotplot.  The colour scale for adjusted P has fixed limits
# [0, padj_cutoff] so that plots with the same cutoff are directly comparable.
nf_rna_gsea_dotplot <- function(context, top, contrast, padj_cutoff, p_adjust_method, path_stem, output_root, meta=list()) {
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
  p <- ggplot2::ggplot(top, ggplot2::aes(x=NES, y=axis, size=setSize, colour=p.adjust)) +
    ggplot2::geom_vline(xintercept=0, linewidth=nf_rna_linewidth(0.5), colour="grey40") +
    ggplot2::geom_point() +
    ggplot2::scale_y_discrete(labels=stats::setNames(labels, top$ID)) +
    ggplot2::scale_colour_gradientn(colours=rev(NF_RNA_SEQUENTIAL_STOPS[-1]), limits=c(0, as.numeric(padj_cutoff)), oob=scales::squish, name=sprintf("Adjusted P (%s)", p_adjust_method)) +
    ggplot2::scale_size_area(max_size=5, name="Gene set size") +
    ggplot2::labs(x=x_label, y=NULL) + nf_rna_theme(context)
  nf_rna_save_figure(context, p, path_stem, context$profile$double_col_mm, nf_rna_dotplot_height(context, labels), output_root, c(list(
    kind="GSEA dotplot", contrast=contrast[intersect(c("contrast_id", "numerator", "denominator"), names(contrast))],
    x="normalized enrichment score; sign gives the direction", rank_metric="DESeq2 Wald statistic (stat) for numerator vs denominator; ties ordered by Entrez ID",
    colour=list(value=sprintf("adjusted P (%s)", p_adjust_method), limits=c(0, as.numeric(padj_cutoff)), kind="sequential", harmonized="fixed limits shared by every GSEA dotplot with this cutoff"),
    size="gene set size (setSize)", display_order="NES descending",
    labels=sprintf("term descriptions wrapped to 45 characters, at most 2 lines; %d truncated; full ID/Description in the source-data table", sum(attr(labels, "truncated"))),
    source_data=list(nf_rna_source_data(source_path, output_root, "plotted terms with full IDs, NES, adjusted P and leading edge"))), meta))
}
