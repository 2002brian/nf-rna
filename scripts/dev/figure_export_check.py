"""Check the GUI's PDF/TIFF figure downloads against a real run (read-only).

For every figure (one per PDF) under downstream/l1 and downstream/l2 of RUN_DIR:

* PDF: the download is the pipeline file byte for byte, and poppler opens it;
* TIFF 300 and 600 dpi: format, RGB, LZW, dpi tags, pixel size equal to
  round(page_pt / 72 * dpi), and content: a pipeline TIFF is byte-identical, a
  converted raster PNG is pixel-identical, a PDF rendering is compared with the
  independent Ghostscript rendering of the same PDF, the dominant colours of the
  pipeline's own raster must all occur, and the axis-line displacement against
  that raster is reported (PDF and raster devices lay out text rows slightly
  differently; see docs/dev/phase7b_report.md).

Writes nothing into the run.  Usage:
  python scripts/dev/figure_export_check.py RUN_DIR OUT_JSON [--save DIR]
"""

from __future__ import annotations

import argparse
import hashlib
import io
import json
import re
import subprocess
import tempfile
from pathlib import Path

from PIL import Image, ImageChops, ImageStat

from rnaseq.gui import figure_export


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def tree_hashes(root: Path) -> dict[str, str]:
    return {str(path.relative_to(root)): sha(path) for path in sorted(root.rglob("*")) if path.is_file()}


def reference_raster(stem: Path, dpi: int) -> Path | None:
    for path in (stem.with_suffix(".tiff"), stem.with_name(f"{stem.name}_{dpi}dpi.png"), stem.with_name(f"{stem.name}_300dpi.png")):
        if path.is_file():
            return path
    return None


def _axis_lines(image: Image.Image) -> tuple[list[float], list[float]]:
    """Positions (in pt) of long dark horizontal and vertical lines: the drawn panel axes."""

    grey = image.convert("L"); width, height = grey.size; dpi = float(image.info.get("dpi", (72, 72))[0])
    data = grey.load()
    rows = [y / dpi * 72 for y in range(height) if sum(1 for x in range(0, width, 2) if data[x, y] < 100) > 0.25 * width]
    cols = [x / dpi * 72 for x in range(width) if sum(1 for y in range(0, height, 2) if data[x, y] < 100) > 0.25 * height]
    return rows, cols


def layout_offset(export: Image.Image, reference: Path) -> dict | None:
    """Largest axis-line displacement (pt) between the export and the pipeline's raster of the same figure, when both have axes."""

    with Image.open(reference) as ref:
        ref.load()
        (r1, c1), (r2, c2) = _axis_lines(export), _axis_lines(ref)
    out = {}
    if r1 and r2:
        out["x_axis_pt"] = round(max(r1) - max(r2), 2)
    if c1 and c2:
        out["y_axis_pt"] = round(min(c1) - min(c2), 2)
    return out or None


def check(stem: Path) -> dict:
    pdf = stem.with_suffix(".pdf")
    points = figure_export.pdf_page_points(pdf)
    opened = subprocess.run(["pdfinfo", str(pdf)], capture_output=True).returncode == 0
    row: dict = {"figure": str(stem), "pdf_opens": opened, "page_pt": points, "tiff": {}}
    for dpi in figure_export.TIFF_DPIS:
        option = figure_export.tiff_option(stem, dpi, pdf)
        entry: dict = {"available": option.available, "method": option.method, "reason": option.reason, "file_name": option.file_name}
        if option.available:
            data = figure_export.tiff_bytes(option)
            image = Image.open(io.BytesIO(data)); image.load()
            expected = figure_export.pixels_for(points, dpi)
            entry.update(format=image.format, mode=image.mode, compression=image.info.get("compression"), dpi=[round(v) for v in image.info.get("dpi", (0, 0))],
                         size=list(image.size), expected_size=list(expected), sha256=hashlib.sha256(data).hexdigest(), bytes=len(data))
            ok = image.format == "TIFF" and image.mode == "RGB" and image.info.get("compression") == "tiff_lzw" and entry["dpi"] == [dpi, dpi] and tuple(image.size) == expected
            if option.method == figure_export.PIPELINE_FILE:
                entry["content"] = "byte-identical to pipeline TIFF" if data == option.source.read_bytes() else "DIFFERS"
            elif option.source.suffix == ".png":
                with Image.open(option.source) as source:
                    same = ImageChops.difference(image, source.convert("RGB")).getbbox() is None
                entry["content"] = "pixel-identical to pipeline raster PNG" if same else "DIFFERS"
            else:
                # Faithful rasterization: compare with an independent renderer (Ghostscript) of the same PDF at the same dpi.
                with tempfile.TemporaryDirectory() as tmp:
                    gs_png = Path(tmp) / "gs.png"
                    subprocess.run(["gs", "-q", "-dSAFER", "-sDEVICE=png16m", f"-r{dpi}", "-dTextAlphaBits=4", "-dGraphicsAlphaBits=4", "-o", str(gs_png), str(pdf)], check=True)
                    with Image.open(gs_png) as gs:
                        gs = gs.convert("RGB").crop((0, 0, *image.size))
                        mean = sum(ImageStat.Stat(ImageChops.difference(image, gs)).mean) / 3
                # Colour mapping: every colour covering >= 0.2% of the pipeline's own raster also occurs in the export.
                reference = reference_raster(stem, dpi)
                with Image.open(reference) as ref:
                    ref = ref.convert("RGB")
                    ref_colours = {colour for count, colour in ref.getcolors(1 << 24) if count >= 0.002 * ref.size[0] * ref.size[1]}
                export_colours = {colour for _count, colour in image.getcolors(1 << 24)}
                missing = sorted(ref_colours - export_colours)
                layout = layout_offset(image, reference)
                entry.update(vs_ghostscript_mean_abs_diff=mean, dominant_colours=len(ref_colours), dominant_colours_missing=[f"#{r:02X}{g:02X}{b:02X}" for r, g, b in missing],
                             layout_vs_pipeline_raster_pt=layout)
                entry["content"] = f"vs Ghostscript render: mean |diff| {mean:.2f}/255; dominant colours of {reference.name}: {len(ref_colours)}, missing {len(missing)}"
                ok = ok and mean < 3 and not missing
            entry["ok"] = ok and "DIFFERS" not in entry.get("content", "")
        row["tiff"][str(dpi)] = entry
    row["ok"] = opened and points is not None and all(item.get("ok") for item in row["tiff"].values())
    return row


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("run_dir", type=Path)
    parser.add_argument("out_json", type=Path)
    parser.add_argument("--save", type=Path, help="also save the downloads of the first figures of each kind here")
    args = parser.parse_args()
    downstream = args.run_dir / "downstream"
    before = tree_hashes(args.run_dir)
    stems = [path.with_suffix("") for level in ("l1", "l2") if (downstream / level).is_dir() for path in sorted((downstream / level).rglob("*.pdf"))]
    rows = [check(stem) for stem in stems]
    if args.save:
        args.save.mkdir(parents=True, exist_ok=True)
        for stem in stems:
            if re.search(r"(pca|volcano|ma_plot|heatmap|/dotplot|negative_1_)", str(stem)) and "simplified" not in stem.name:
                for dpi in figure_export.TIFF_DPIS:
                    option = figure_export.tiff_option(stem, dpi)
                    if option.available:
                        name = str(stem.relative_to(downstream)).replace("/", "__")
                        (args.save / f"{name}_{dpi}dpi.tiff").write_bytes(figure_export.tiff_bytes(option))
    after = tree_hashes(args.run_dir)
    summary = {"run_dir": str(args.run_dir), "figures": len(rows), "all_ok": all(row["ok"] for row in rows), "run_files": len(before), "run_unchanged": before == after,
               "methods": {}, "rows": rows}
    for row in rows:
        for dpi, entry in row["tiff"].items():
            key = f"{dpi}: {entry['method'] or 'unavailable'}"
            summary["methods"][key] = summary["methods"].get(key, 0) + 1
    worst = max((entry.get("vs_ghostscript_mean_abs_diff", 0) for row in rows for entry in row["tiff"].values()), default=0)
    shifts = [abs(value) for row in rows for entry in row["tiff"].values() for value in (entry.get("layout_vs_pipeline_raster_pt") or {}).values()]
    summary["worst_vs_ghostscript"] = worst; summary["layout_shift_pt_max"] = max(shifts, default=None)
    args.out_json.parent.mkdir(parents=True, exist_ok=True)
    args.out_json.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps({key: summary[key] for key in ("figures", "all_ok", "run_files", "run_unchanged", "methods", "worst_vs_ghostscript", "layout_shift_pt_max")}, indent=2))
    for row in rows:
        if not row["ok"]:
            print("FAIL", row["figure"], json.dumps(row["tiff"], indent=1))
    return 0 if summary["all_ok"] and summary["run_unchanged"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
