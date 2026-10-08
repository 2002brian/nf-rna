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
      plot.subtitle=ggplot2::element_text(size=profile$base_pt, colour="black", lineheight=1.05, margin=ggplot2::margin(b=3)),
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
  packages <- c("ggplot2", "ggrepel", "pheatmap", "scales", "systemfonts")
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
  stats$sample_id <- factor(samples, levels = samples)
  stats$group <- group
  colours <- nf_rna_categorical_colors(levels(group))
  p <- ggplot2::ggplot(stats, ggplot2::aes(x = sample_id)) +
    ggplot2::geom_boxplot(ggplot2::aes(ymin = whisker_low, lower = q1, middle = median, upper = q3, ymax = whisker_high, fill = group),
                          stat = "identity", width = 0.7, linewidth = nf_rna_linewidth(0.5), colour = "black") +
    ggplot2::scale_fill_manual(values = colours, name = group_name, drop = FALSE)
  shapes <- NULL
  if (!is.null(batch)) {
    shapes <- nf_rna_categorical_shapes(levels(batch))
    if (!is.null(shapes)) {
      span <- max(stats$whisker_high, na.rm = TRUE) - min(stats$whisker_low, na.rm = TRUE)
      stats$batch <- batch; stats$batch_y <- min(stats$whisker_low, na.rm = TRUE) - 0.08 * max(span, 1e-9)
      p <- p + ggplot2::geom_point(data = stats, ggplot2::aes(x = sample_id, y = batch_y, shape = batch), size = 1.8, colour = "black", inherit.aes = FALSE) +
        ggplot2::scale_shape_manual(values = shapes, name = batch_name, drop = FALSE)
    }
  }
  show_labels <- length(samples) <= 100
  rotation <- if (length(samples) > 8) 90 else 45
  p <- p + ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = 0.05)) +
    ggplot2::labs(x = NULL, y = y_label) + nf_rna_theme(context) +
    ggplot2::theme(legend.position = "bottom", axis.text.x = if (show_labels) ggplot2::element_text(angle = rotation, hjust = 1, vjust = if (rotation == 90) 0.5 else 1) else ggplot2::element_blank(),
                   axis.ticks.x = if (show_labels) ggplot2::element_line() else ggplot2::element_blank())
  label_mm <- if (show_labels) nf_rna_pt_to_mm(profile$small_pt) * 0.6 * max(nchar(samples)) else 0
  list(plot = p, colours = colours, shapes = shapes, labels_shown = show_labels, rotation = if (show_labels) rotation else NA_real_,
       height_mm = min(profile$max_height_mm, 65 + label_mm))
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
    ggplot2::labs(x=x_label, y=NULL, subtitle=subtitle) + nf_rna_theme(context)
  if (!is.null(subtitle)) p <- p + ggplot2::theme(plot.title.position="plot")
  subtitle_lines <- if (is.null(subtitle)) 0 else length(strsplit(subtitle, "\n", fixed=TRUE)[[1]])
  nf_rna_save_figure(context, p, path_stem, context$profile$double_col_mm, nf_rna_dotplot_height(context, labels, subtitle_lines), output_root, c(list(
    kind="ORA dotplot", subtitle=subtitle, x="-log10(adjusted P) as reported by clusterProfiler", size="number of foreground genes in the term (Count)",
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
