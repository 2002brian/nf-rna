"""Central R figure policy: display-only transformations, colour policy and physical export."""

from __future__ import annotations

import hashlib
import json
import re
import struct
import subprocess
import zlib
from pathlib import Path

import pytest

from conftest import pdf_text_origins, require_r_packages


FIGURES_R = Path(__file__).resolve().parents[1] / "src" / "rnaseq" / "r" / "figures.R"


def run_r(code: str, cwd: Path) -> dict:
    executable = require_r_packages("ggplot2", "pheatmap", "jsonlite")
    script = cwd / "check.R"
    script.write_text(
        f'suppressPackageStartupMessages({{library(ggplot2); library(pheatmap); library(jsonlite)}})\nsource("{FIGURES_R}")\n'
        + code + '\n',
        encoding="utf-8",
    )
    result = subprocess.run([executable, str(script)], cwd=cwd, capture_output=True, text=True, check=False)
    assert result.returncode == 0, result.stderr
    return json.loads(result.stdout.strip().splitlines()[-1])


def tiff_tags(path: Path) -> dict[int, object]:
    data = path.read_bytes()
    order = "<" if data[:2] == b"II" else ">"
    assert struct.unpack(order + "H", data[2:4])[0] == 42
    offset = struct.unpack(order + "I", data[4:8])[0]
    count = struct.unpack(order + "H", data[offset:offset + 2])[0]
    sizes = {1: 1, 2: 1, 3: 2, 4: 4, 5: 8}
    tags: dict[int, object] = {}
    for index in range(count):
        entry = data[offset + 2 + 12 * index: offset + 14 + 12 * index]
        tag, kind, number = struct.unpack(order + "HHI", entry[:8])
        size = sizes[kind] * number
        raw = entry[8:12] if size <= 4 else data[struct.unpack(order + "I", entry[8:12])[0]:][:size]
        if kind == 3:
            values = list(struct.unpack(order + "H" * number, raw[:2 * number]))
        elif kind == 4:
            values = list(struct.unpack(order + "I" * number, raw[:4 * number]))
        elif kind == 5:
            pairs = struct.unpack(order + "I" * (2 * number), raw)
            values = [pairs[i] / pairs[i + 1] for i in range(0, len(pairs), 2)]
        else:
            values = raw
        tags[tag] = values[0] if isinstance(values, list) and len(values) == 1 else values
    tags["_data"] = data
    return tags


def first_lzw_bytes(strip: bytes, wanted: int) -> bytes:
    """Decode the first bytes of a TIFF LZW strip (MSB-first codes, early change)."""

    table = [bytes([i]) for i in range(256)] + [b"", b""]
    out, previous, position, width = bytearray(), None, 0, 9
    while len(out) < wanted:
        code = int.from_bytes(strip[position // 8: position // 8 + 3], "big") >> (24 - width - position % 8) & ((1 << width) - 1)
        position += width
        if code == 257:
            break
        if code == 256:
            table, previous, width = table[:258], None, 9
            continue
        entry = table[code] if code < len(table) else previous + previous[:1]
        if previous is not None:
            table.append(previous + entry[:1])
        out += entry
        previous = entry
        if len(table) >= (1 << width) - 1 and width < 12:
            width += 1
    return bytes(out[:wanted])


def pdf_pages(path: Path) -> list[tuple[float, ...]]:
    data = path.read_bytes()
    chunks = [data]
    for match in re.finditer(rb"stream\r?\n(.*?)endstream", data, re.S):
        try:
            chunks.append(zlib.decompress(match.group(1)))
        except zlib.error:
            pass
    boxes = []
    for chunk in chunks:
        boxes += [tuple(float(value) for value in box.split()) for box in re.findall(rb"/MediaBox\s*\[\s*([^\]]+)\]", chunk)]
    return boxes


def test_row_zscore_preserves_missing_values_and_marks_constant_rows_undefined(tmp_path):
    result = run_r(r'''
x <- matrix(c(1, 2, 3, 4,   5, 5, 5, 5,   NA, 1, 2, 3,   7, NA, NA, NA,   2, 4, 6, 8), 5, byrow=TRUE, dimnames=list(c("a","b","c","d","e"), c("S1","S2","S3","S4")))
original <- x
z <- nf_rna_row_zscore(x)
complete <- x[c("a","e"), ]
cat(toJSON(list(
  unchanged=identical(x, original),
  equals_scale=isTRUE(all.equal(unname(z[c("a","e"), ]), unname(t(scale(t(complete)))), check.attributes=FALSE)),
  constant_row_na=all(is.na(z["b", ])), single_value_row_na=all(is.na(z["d", ])),
  missing_cell_na=is.na(z["c", "S1"]), partial_row_defined=all(is.finite(z["c", 2:4])),
  no_zero_imputation=!any(z == 0 & is.na(x), na.rm=TRUE),
  undefined=attr(z, "undefined_rows")), auto_unbox=TRUE))
''', tmp_path)
    assert result["unchanged"] and result["equals_scale"]
    assert result["constant_row_na"] and result["single_value_row_na"]
    assert result["missing_cell_na"] and result["partial_row_defined"] and result["no_zero_imputation"]
    assert result["undefined"] == ["b", "d"]


def test_colour_saturation_is_visual_only(tmp_path):
    result = run_r(r'''
z <- matrix(c(-3.5, -1, 0, NA, 2.5, 1.9), 2, dimnames=list(c("g1","g2"), c("A","B","C")))
before <- z
shown <- nf_rna_saturate(z, 2)
cat(toJSON(list(source_unchanged=identical(z, before), shown=as.vector(shown), source=as.vector(z)), auto_unbox=TRUE, na="null"))
''', tmp_path)
    assert result["source_unchanged"]
    assert result["source"] == [-3.5, -1, 0, None, 2.5, 1.9]
    assert result["shown"] == [-2, -1, 0, None, 2, 1.9]


def test_sequential_scale_for_distances_starts_at_zero_and_handles_all_zero_matrix(tmp_path):
    result = run_r(r'''
d <- as.matrix(dist(matrix(c(0, 0, 0, 3, 4, 0), 3)))
zero <- matrix(0, 3, 3)
palette <- nf_rna_sequential_palette(9)
lum <- apply(grDevices::col2rgb(palette), 2, function(v) sum(v * c(0.2126, 0.7152, 0.0722)))
cat(toJSON(list(distance=nf_rna_sequential_limits(d, lower=0), zero=nf_rna_sequential_limits(zero, lower=0),
  correlation=nf_rna_sequential_limits(c(0.97, 1), upper=1), identical_corr=nf_rna_sequential_limits(1, upper=1),
  monotonic=all(diff(lum) < 0), breaks_valid=all(diff(seq(0, nf_rna_sequential_limits(zero, lower=0)$limits[2], length.out=10)) > 0)), auto_unbox=TRUE))
''', tmp_path)
    assert result["distance"]["limits"][0] == 0 and result["distance"]["degenerate"] is False
    assert result["zero"] == {"limits": [0, 1], "degenerate": True}
    assert result["correlation"]["limits"] == [0.97, 1]
    assert result["identical_corr"] == {"limits": [0, 1], "degenerate": True}
    # One-directional lightness ramp: sequential, not diverging.
    assert result["monotonic"] and result["breaks_valid"]


def test_categorical_colours_are_named_deterministic_and_not_interpolated(tmp_path):
    result = run_r(r'''
a <- nf_rna_categorical_colors(c("Control", "Treatment"))
b <- nf_rna_categorical_colors(c("Control", "Treatment"))
many <- nf_rna_categorical_colors(sprintf("L%02d", 1:10))
cat(toJSON(list(a=as.list(a), same=identical(a, b), okabe=unname(a) %in% unname(NF_RNA_OKABE_ITO), many_unique=length(unique(many)) == 10,
  shapes=as.list(nf_rna_categorical_shapes(c("Control", "Treatment"))), too_many_shapes=is.null(nf_rna_categorical_shapes(sprintf("L%02d", 1:13))),
  palette_many=nf_rna_categorical_palette_name(10)), auto_unbox=TRUE))
''', tmp_path)
    assert result["a"] == {"Control": "#0072B2", "Treatment": "#D55E00"}
    assert result["same"] and all(result["okabe"])
    assert result["many_unique"] and result["palette_many"] == "HCL qualitative 'Dark 3'"
    assert result["shapes"] == {"Control": 16, "Treatment": 17} and result["too_many_shapes"]


def test_point_size_helpers(tmp_path):
    result = run_r(r'''
cat(toJSON(list(mm=nf_rna_pt_to_mm(72.27), text=nf_rna_text_size(7), line=nf_rna_linewidth(ggplot2::.pt)), auto_unbox=TRUE, digits=NA))
''', tmp_path)
    assert result["mm"] == pytest.approx(25.4)
    assert result["text"] == {"size": 7, "size.unit": "pt"}
    assert result["line"] == pytest.approx(1)


def test_display_labels_wrap_and_truncate_but_keep_full_text_elsewhere(tmp_path):
    result = run_r(r'''
labels <- nf_rna_display_labels(c("short", paste(rep("very long pathway description", 6), collapse=" ")))
cat(toJSON(list(labels=as.vector(labels), truncated=attr(labels, "truncated")), auto_unbox=TRUE))
''', tmp_path)
    assert result["labels"][0] == "short"
    assert result["labels"][1].count("\n") == 1 and result["labels"][1].endswith("...")
    assert result["truncated"] == [False, True]


def test_export_draws_each_format_at_the_physical_size_and_records_it(tmp_path):
    result = run_r(r'''
ctx <- nf_rna_figure_context(list())
p <- ggplot(data.frame(x=c(-2, 3), y=c(1, 4)), aes(x, y)) + geom_point() + nf_rna_theme(ctx)
heat <- pheatmap(matrix(c(1, 2, 3, 4), 2), silent=TRUE, fontsize=ctx$profile$base_pt)
nf_rna_save_figure(ctx, p, file.path(getwd(), "plot"), 85, 60, getwd(), list(kind="fixture"))
dir.create("sub"); nf_rna_save_figure(ctx, heat, file.path(getwd(), "sub", "heat"), 180, 150, getwd())
nf_rna_write_figure_manifest(ctx, getwd())
cat(toJSON(list(font=ctx$font), auto_unbox=TRUE))
''', tmp_path)
    font = result["font"]
    assert font["requested"] == "Arial"
    assert font["family"] in {"Arial", "DejaVu Sans"}
    assert font["fallback_used"] == (font["family"] != "Arial")
    manifest = json.loads((tmp_path / "figure_manifest.json").read_text())
    assert manifest["schema_version"] == "nf-rna.figure-manifest.v1"
    assert manifest["profile"]["name"] == "general" and manifest["font"]["family"] == font["family"]
    records = {item["id"]: item for item in manifest["figures"]}
    assert set(records) == {"plot", "sub/heat"} and records["plot"]["kind"] == "fixture"
    for stem, (width_pt, height_pt) in {"plot": (241, 170), "sub/heat": (510, 425)}.items():
        files = {item["format"]: item for item in records[stem]["files"]}
        assert set(files) == {"png", "tiff", "pdf"}
        assert records[stem]["width_mm"] == pytest.approx(width_pt * 25.4 / 72)
        # Every format is drawn directly; no output is a resampled copy of another.
        for fmt, dpi in (("tiff", 300), ("png", 150)):
            assert (files[fmt]["width_px"], files[fmt]["height_px"]) == (round(width_pt / 72 * dpi), round(height_pt / 72 * dpi))
        tags = tiff_tags(tmp_path / f"{stem}.tiff")
        assert (tags[256], tags[257]) == (files["tiff"]["width_px"], files["tiff"]["height_px"])
        assert tags[259] == 5  # LZW
        assert tags[262] == 2  # RGB
        assert tags[277] == 3 and 338 not in tags  # no alpha / extra samples
        assert tags[282] == tags[283] == 300 and tags[296] == 2  # 300 pixels per inch
        strip = tags["_data"][tags[273] if isinstance(tags[273], int) else tags[273][0]:]
        assert first_lzw_bytes(strip, 3) == b"\xff\xff\xff"  # opaque white background
        png = (tmp_path / f"{stem}.png").read_bytes()
        assert struct.unpack(">II", png[16:24]) == (files["png"]["width_px"], files["png"]["height_px"])
        assert pdf_pages(tmp_path / f"{stem}.pdf") == [(0.0, 0.0, float(width_pt), float(height_pt))]


def test_pheatmap_title_spans_the_figure_and_stays_on_the_page(tmp_path):
    # A narrow matrix with long row labels and a legend: pheatmap centres its
    # title over the matrix column, so a long title starts left of the page.
    result = run_r(r'''
ctx <- nf_rna_figure_context(list())
m <- matrix(c(1, 2, 3, 4, 5, 6, 2, 1, 3), 3, dimnames=list(paste0("ENSMUSG0000000000", 1:3, "_long_label"), c("S1", "S2", "S3")))
title <- nf_rna_heatmap_title(ctx, "Resistant", "Parental", "Top 50 DEGs \u00b7 row z-score of VST", width_mm=85)
make <- function() pheatmap(m, main=title, fontsize=ctx$profile$base_pt, silent=TRUE)
nf_rna_save_figure(ctx, make(), file.path(getwd(), "unspanned"), 85, 80, getwd(), formats="pdf")
spanned <- nf_rna_span_pheatmap_title(make())
nf_rna_save_figure(ctx, spanned, file.path(getwd(), "spanned"), 85, 80, getwd(), formats="pdf")
main <- spanned$gtable$layout[spanned$gtable$layout$name == "main", ]
grDevices::pdf(NULL, width=85 / 25.4, height=80 / 25.4)
title_width <- grid::convertWidth(grid::grobWidth(spanned$gtable$grobs[[which(spanned$gtable$layout$name == "main")]]), "mm", valueOnly=TRUE)
grDevices::dev.off()
cat(toJSON(list(l=main$l, r=main$r, columns=ncol(spanned$gtable), title_width_mm=title_width), auto_unbox=TRUE))
''', tmp_path)
    assert (result["l"], result["r"]) == (1, result["columns"])
    assert result["title_width_mm"] < 85
    assert min(x for x, _ in pdf_text_origins(tmp_path / "unspanned.pdf")) < 0  # the defect this guards against
    spanned = pdf_text_origins(tmp_path / "spanned.pdf")
    assert spanned and min(x for x, _ in spanned) >= 0
    # Spanning the title never changes the figure size.
    assert pdf_pages(tmp_path / "spanned.pdf") == pdf_pages(tmp_path / "unspanned.pdf")


def test_heatmap_title_breaks_long_group_names_instead_of_truncating(tmp_path):
    result = run_r(r'''
ctx <- nf_rna_figure_context(list())
detail <- "Top 50 DEGs \u00b7 row z-score of VST"
long_num <- "Irradiated_resistant_subline_133Gy_passage_12"; long_den <- "Irradiated_parental_line_22Gy_passage_3"
cat(toJSON(list(short=nf_rna_heatmap_title(ctx, "133Gy", "22Gy", detail, width_mm=180),
  long=nf_rna_heatmap_title(ctx, long_num, long_den, detail, width_mm=180)), auto_unbox=TRUE))
''', tmp_path)
    assert result["short"] == "133Gy vs 22Gy\nTop 50 DEGs \u00b7 row z-score of VST"
    assert result["long"] == (
        "Irradiated_resistant_subline_133Gy_passage_12\nvs Irradiated_parental_line_22Gy_passage_3\nTop 50 DEGs \u00b7 row z-score of VST"
    )


def test_tiff_export_is_byte_deterministic(tmp_path):
    code = r'''
ctx <- nf_rna_figure_context(list())
set.seed(1)
p <- ggplot(data.frame(x=rnorm(200), y=rnorm(200)), aes(x, y)) + geom_point() + nf_rna_theme(ctx)
nf_rna_save_figure(ctx, p, file.path(getwd(), "plot"), 85, 60, getwd())
cat("{}")
'''
    digests = []
    for name in ("one", "two"):
        (tmp_path / name).mkdir()
        run_r(code, tmp_path / name)
        digests.append({fmt: hashlib.sha256((tmp_path / name / f"plot.{fmt}").read_bytes()).hexdigest() for fmt in ("tiff", "png")})
    assert digests[0] == digests[1]


def test_profiles_carry_the_reviewed_targets(tmp_path):
    result = run_r("cat(toJSON(NF_RNA_FIGURE_PROFILES, auto_unbox=TRUE, digits=NA))", tmp_path)
    assert result["general"]["raster_dpi"] == 300 and result["general"]["max_width_mm"] == 180
    assert (result["nature"]["single_col_mm"], result["nature"]["double_col_mm"], result["nature"]["max_height_mm"]) == (89, 183, 170)
    assert 5 <= result["nature"]["min_pt"] <= result["nature"]["base_pt"] <= 7
    assert (result["plos"]["max_width_mm"], result["plos"]["max_height_mm"], result["plos"]["min_pt"]) == (190.5, 222.3, 8)
    assert result["elsevier"]["raster_dpi"] == 500
    for profile in result.values():
        assert profile["min_pt"] <= profile["small_pt"] <= profile["base_pt"]
        assert profile["double_col_mm"] <= profile["max_width_mm"]
