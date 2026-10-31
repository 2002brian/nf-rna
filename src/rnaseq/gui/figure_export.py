"""PDF and TIFF downloads of the pipeline's own figures (no Streamlit here).

Nothing is redrawn from data and nothing is written into a run.  For each
figure the GUI offers:

* PDF: the pipeline's Cairo PDF, byte for byte.
* TIFF at 300 dpi (default) or 600 dpi, from the best source that really has
  that resolution, in this order:
  1. the pipeline's own TIFF when its recorded resolution is the requested dpi
     (offered byte for byte);
  2. the pipeline's own ``<stem>_<dpi>dpi.png`` (drawn by R at that dpi and at
     the figure's final physical size), converted losslessly to TIFF: every
     pixel is unchanged, only the container changes;
  3. the pipeline's own vector PDF, rasterized with poppler's ``pdftocairo`` at
     exactly ``dpi`` and cropped to ``round(page_pt / 72 * dpi)`` pixels, the
     pixel rule R uses for its own PNG/TIFF devices.

A raster file is never enlarged to reach a higher dpi.  When none of the
sources exists the option is unavailable and says why.  TIFF output is RGB on
an opaque white background with LZW compression and the dpi recorded in its
resolution tags, the same encoding as the pipeline's TIFF.
"""

from __future__ import annotations

import io
import re
import shutil
import subprocess
import tempfile
from dataclasses import dataclass
from pathlib import Path

from PIL import Image

TIFF_DPIS = (300, 600)
DEFAULT_TIFF_DPI = 300
RENDER_TIMEOUT_S = 120

PIPELINE_FILE = "pipeline TIFF, unchanged"
PNG_CONVERSION = "lossless conversion of the pipeline's {dpi} dpi PNG (pixels unchanged)"
PDF_RENDER = "rendered from the pipeline's vector PDF at {dpi} dpi (text layout as in the PDF)"


@dataclass(frozen=True)
class TiffOption:
    """How a TIFF at ``dpi`` would be obtained for one figure, or why it cannot."""

    dpi: int
    file_name: str
    source: Path | None
    method: str | None
    reason: str | None = None

    @property
    def available(self) -> bool:
        return self.source is not None


def file_dpi(path: Path) -> int | None:
    """The horizontal resolution recorded in a raster file, rounded; ``None`` if absent."""

    try:
        with Image.open(path) as image:
            dpi = image.info.get("dpi")
    except (OSError, ValueError):
        return None
    return int(round(float(dpi[0]))) if dpi else None


def pdf_page_points(pdf: Path) -> tuple[float, float] | None:
    """Page size of the first page in points, from ``pdfinfo``; ``None`` if unreadable."""

    tool = shutil.which("pdfinfo")
    if tool is None:
        return None
    try:
        out = subprocess.run([tool, str(pdf)], capture_output=True, text=True, timeout=30, check=False).stdout
    except (OSError, subprocess.TimeoutExpired):
        return None
    match = re.search(r"^Page size:\s+([\d.]+) x ([\d.]+) pts", out, re.MULTILINE)
    return (float(match.group(1)), float(match.group(2))) if match else None


def pixels_for(points: tuple[float, float], dpi: int) -> tuple[int, int]:
    """Pixel size of a page drawn at ``dpi``, rounded as ``nf_rna_save_figure`` does (R ``round``: half to even)."""

    return tuple(round(value / 72 * dpi) for value in points)  # type: ignore[return-value]


def tiff_option(stem: Path, dpi: int, pdf: Path | None = None) -> TiffOption:
    """The source for a ``dpi`` TIFF of the figure ``stem`` (path without suffix)."""

    pdf = pdf if pdf is not None else stem.with_suffix(".pdf")
    pipeline_tiff = stem.with_suffix(".tiff")
    if pipeline_tiff.is_file() and file_dpi(pipeline_tiff) == dpi:
        return TiffOption(dpi, pipeline_tiff.name, pipeline_tiff, PIPELINE_FILE)
    name = f"{stem.name}_{dpi}dpi.tiff"
    raster_png = stem.with_name(f"{stem.name}_{dpi}dpi.png")
    if raster_png.is_file() and file_dpi(raster_png) == dpi:
        return TiffOption(dpi, name, raster_png, PNG_CONVERSION.format(dpi=dpi))
    if not pdf.is_file():
        return TiffOption(dpi, name, None, None, f"no {dpi} dpi raster and no PDF for this figure")
    if shutil.which("pdftocairo") is None or shutil.which("pdfinfo") is None:
        return TiffOption(dpi, name, None, None, f"no {dpi} dpi raster from the pipeline, and poppler (pdftocairo, pdfinfo) is not installed to render the PDF")
    return TiffOption(dpi, name, pdf, PDF_RENDER.format(dpi=dpi))


def tiff_options(stem: Path, pdf: Path | None = None) -> dict[int, TiffOption]:
    return {dpi: tiff_option(stem, dpi, pdf) for dpi in TIFF_DPIS}


def _encode(image: Image.Image, dpi: int) -> bytes:
    if image.mode != "RGB":
        background = Image.new("RGB", image.size, "white")
        rgba = image.convert("RGBA")
        background.paste(rgba, mask=rgba.getchannel("A"))
        image = background
    buffer = io.BytesIO()
    image.save(buffer, format="TIFF", compression="tiff_lzw", dpi=(dpi, dpi))
    return buffer.getvalue()


def _render_pdf(pdf: Path, dpi: int) -> bytes:
    points = pdf_page_points(pdf)
    if points is None:
        raise RuntimeError(f"cannot read the page size of {pdf}")
    width, height = pixels_for(points, dpi)
    with tempfile.TemporaryDirectory(prefix="nf-rna-tiff-") as tmp:
        target = Path(tmp) / "page"
        # Draw at exactly ``dpi`` and crop to R's pixel count: poppler rounds the page up, so this drops at most
        # one pixel row/column of the (white) page edge and never rescales the drawing.
        subprocess.run(["pdftocairo", "-png", "-singlefile", "-f", "1", "-l", "1", "-r", str(dpi), "-x", "0", "-y", "0", "-W", str(width), "-H", str(height),
                        str(pdf), str(target)], check=True, capture_output=True, timeout=RENDER_TIMEOUT_S)
        with Image.open(target.with_suffix(".png")) as image:
            image.load()
            if image.size != (width, height):
                raise RuntimeError(f"rendered {image.size[0]}x{image.size[1]} px, expected {width}x{height}")
            return _encode(image, dpi)


def tiff_bytes(option: TiffOption) -> bytes:
    """The TIFF file for an available option (read from disk or encoded in memory)."""

    if option.source is None:
        raise ValueError(option.reason or "TIFF unavailable")
    if option.method == PIPELINE_FILE:
        return option.source.read_bytes()
    if option.source.suffix == ".png":
        with Image.open(option.source) as image:
            image.load()
            return _encode(image, option.dpi)
    return _render_pdf(option.source, option.dpi)
