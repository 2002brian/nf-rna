"""PDF/TIFF figure downloads: pipeline files unchanged, no enlarged rasters, correct TIFF dpi and pixel size."""

from __future__ import annotations

import io
import shutil

import pytest
from PIL import Image, ImageChops, ImageDraw

from rnaseq.gui import figure_export

needs_poppler = pytest.mark.skipif(shutil.which("pdftocairo") is None or shutil.which("pdfinfo") is None, reason="poppler not installed")


def _drawing(size: tuple[int, int]) -> Image.Image:
    image = Image.new("RGB", size, "white")
    draw = ImageDraw.Draw(image)
    draw.rectangle([0, 0, size[0] // 2, size[1] - 1], fill=(200, 30, 40))  # left half red, right half white
    return image


def _figure(directory, stem="volcano", tiff_dpi=300, raster_png_dpi=None, pdf=True):
    """A figure as the pipeline writes it: 150 dpi preview, optional TIFF / raster PNG / PDF of a 144 x 72 pt page."""

    page = (144, 72)  # points: 2 x 1 inch
    _drawing((300, 150)).save(directory / f"{stem}.png", dpi=(150, 150))
    if tiff_dpi:
        _drawing((2 * tiff_dpi, tiff_dpi)).save(directory / f"{stem}.tiff", compression="tiff_lzw", dpi=(tiff_dpi, tiff_dpi))
    if raster_png_dpi:
        _drawing((2 * raster_png_dpi, raster_png_dpi)).save(directory / f"{stem}_{raster_png_dpi}dpi.png", dpi=(raster_png_dpi, raster_png_dpi))
    if pdf:
        _drawing(page).save(directory / f"{stem}.pdf", resolution=72)  # page size in points = pixels at 72 dpi
    return directory / stem


def _open(data: bytes) -> Image.Image:
    image = Image.open(io.BytesIO(data))
    image.load()
    return image


def test_pipeline_tiff_and_pdf_are_offered_unchanged(tmp_path):
    stem = _figure(tmp_path)
    option = figure_export.tiff_option(stem, 300)
    assert option.method == figure_export.PIPELINE_FILE and option.file_name == "volcano.tiff"
    assert figure_export.tiff_bytes(option) == stem.with_suffix(".tiff").read_bytes()


def test_lossless_conversion_of_the_raster_png_keeps_every_pixel(tmp_path):
    stem = _figure(tmp_path, stem="positive_1_hsa04110", tiff_dpi=None, raster_png_dpi=300)  # a GSEA running-score plot
    option = figure_export.tiff_option(stem, 300)
    assert option.method.startswith("lossless conversion") and option.file_name == "positive_1_hsa04110_300dpi.tiff"
    tiff = _open(figure_export.tiff_bytes(option))
    assert tiff.format == "TIFF" and tiff.mode == "RGB" and tiff.info["compression"] == "tiff_lzw"
    assert tiff.info["dpi"] == (300, 300) and tiff.size == (600, 300)
    with Image.open(stem.with_name("positive_1_hsa04110_300dpi.png")) as source:
        assert ImageChops.difference(tiff, source.convert("RGB")).getbbox() is None


@needs_poppler
def test_600_dpi_is_rendered_from_the_pdf_at_the_physical_size(tmp_path):
    stem = _figure(tmp_path)
    option = figure_export.tiff_option(stem, 600)
    assert option.method.startswith("rendered from the pipeline's vector PDF")
    tiff = _open(figure_export.tiff_bytes(option))
    assert tiff.size == (1200, 600) and tiff.info["dpi"] == (600, 600)  # 2 x 1 inch at 600 dpi
    assert tiff.mode == "RGB" and tiff.info["compression"] == "tiff_lzw"
    # (the fixture PDF holds a JPEG-compressed image, hence the small tolerance; Cairo PDFs are vector)
    assert max(abs(a - b) for a, b in zip(tiff.getpixel((100, 300)), (200, 30, 40))) <= 2
    assert min(tiff.getpixel((1100, 300))) >= 253
    # The pipeline TIFF is not used for 600 dpi, and is never enlarged.
    assert figure_export.tiff_option(stem, 300).method == figure_export.PIPELINE_FILE


@needs_poppler
def test_a_tiff_at_another_dpi_is_not_relabelled(tmp_path):
    stem = _figure(tmp_path, tiff_dpi=500)  # the elsevier profile draws at 500 dpi
    option = figure_export.tiff_option(stem, 300)
    assert option.method.startswith("rendered from the pipeline's vector PDF")
    assert _open(figure_export.tiff_bytes(option)).size == (600, 300)


def test_without_a_pdf_or_raster_tiff_is_unavailable(tmp_path):
    stem = _figure(tmp_path, tiff_dpi=None, pdf=False)  # only the 150 dpi preview: never enlarged to 300/600
    for dpi, option in figure_export.tiff_options(stem).items():
        assert not option.available and "no PDF" in option.reason
        with pytest.raises(ValueError):
            figure_export.tiff_bytes(option)


def test_pixel_rule_matches_r_devices():
    assert figure_export.pixels_for((306, 577), 300) == (1275, 2404)  # the airway DEG heatmap
    assert figure_export.pixels_for((306, 577), 600) == (2550, 4808)
    assert figure_export.pixels_for((339, 255), 300) == (1412, 1062)  # 1412.5 -> 1412, as R's round()
