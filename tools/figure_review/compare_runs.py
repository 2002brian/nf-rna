#!/usr/bin/env python3
"""Side-by-side comparison of two nf-rna figure runs (stdlib only).

Typical use: ``--legacy`` is a legacy-preset run (for example the Phase 0
baseline) and ``--candidate`` a run of the same project with another preset.
Figures are matched through their ``figure_manifest.json`` files, so the tool
never guesses by file name.  For every figure it reports the physical size and
pixel/point sizes recorded in the manifest *and* re-read from the files, the font
named in the manifest *and* the font actually embedded in the PDF, whether the
PNG pixels are identical (with a highlighted difference image when not), and
whether the source-data tables are byte-identical.

The output is deterministic (no timestamps), so the report itself can be diffed.

    compare_runs.py --legacy RUN_A --candidate RUN_B --out REPORT_DIR

Comparing two releases, ``--expected-changes FILE`` (JSON, see
``expectations/``) gates on a declared list of intended changes, each with its
reason: figure kinds whose content may change, new figures matched by ID pattern
and kind, and added output formats.  Everything not declared must be identical,
and a figure missing from the candidate always fails.
"""

from __future__ import annotations

import argparse
import hashlib
import html
import json
import re
import shutil
import struct
import subprocess
import sys
import zlib
from pathlib import Path


# ---------------------------------------------------------------------- file readers


def read_png(path: Path) -> tuple[int, int, int, list[bytes]]:
    """Decode an 8-bit, non-interlaced RGB/RGBA PNG into scanlines of channel bytes."""

    data = path.read_bytes()
    if data[:8] != b"\x89PNG\r\n\x1a\n":
        raise ValueError(f"not a PNG: {path}")
    pos, idat, header = 8, [], None
    while pos < len(data):
        length, kind = struct.unpack(">I4s", data[pos:pos + 8])
        body = data[pos + 8:pos + 8 + length]
        if kind == b"IHDR":
            header = struct.unpack(">IIBBBBB", body)
        elif kind == b"IDAT":
            idat.append(body)
        pos += 12 + length
    assert header is not None
    width, height, depth, colour, _, _, interlace = header
    if depth != 8 or colour not in (2, 6) or interlace:
        raise ValueError(f"unsupported PNG layout in {path} (depth {depth}, colour type {colour}, interlace {interlace})")
    channels = 3 if colour == 2 else 4
    stride = width * channels
    raw = zlib.decompress(b"".join(idat))
    rows: list[bytes] = []
    previous = bytearray(stride)
    for y in range(height):
        base = y * (stride + 1)
        filter_type, line = raw[base], bytearray(raw[base + 1:base + 1 + stride])
        for i in range(stride):
            left = line[i - channels] if i >= channels else 0
            up = previous[i]
            upleft = previous[i - channels] if i >= channels else 0
            if filter_type == 1:
                line[i] = (line[i] + left) & 255
            elif filter_type == 2:
                line[i] = (line[i] + up) & 255
            elif filter_type == 3:
                line[i] = (line[i] + (left + up) // 2) & 255
            elif filter_type == 4:
                p = left + up - upleft
                pa, pb, pc = abs(p - left), abs(p - up), abs(p - upleft)
                line[i] = (line[i] + (left if pa <= pb and pa <= pc else up if pb <= pc else upleft)) & 255
        rows.append(bytes(line))
        previous = line
    return width, height, channels, rows


def write_png(path: Path, width: int, height: int, rows: list[bytes]) -> None:
    def chunk(kind: bytes, body: bytes) -> bytes:
        return struct.pack(">I", len(body)) + kind + body + struct.pack(">I", zlib.crc32(kind + body) & 0xFFFFFFFF)

    raw = b"".join(b"\x00" + row for row in rows)
    path.write_bytes(b"\x89PNG\r\n\x1a\n" + chunk(b"IHDR", struct.pack(">IIBBBBB", width, height, 8, 2, 0, 0, 0)) + chunk(b"IDAT", zlib.compress(raw, 9)) + chunk(b"IEND", b""))


def tiff_info(path: Path) -> dict:
    data = path.read_bytes()
    order = "<" if data[:2] == b"II" else ">"
    offset = struct.unpack(order + "I", data[4:8])[0]
    count = struct.unpack(order + "H", data[offset:offset + 2])[0]
    sizes = {1: 1, 2: 1, 3: 2, 4: 4, 5: 8}
    tags: dict[int, object] = {}
    for index in range(count):
        entry = data[offset + 2 + 12 * index:offset + 14 + 12 * index]
        tag, kind, number = struct.unpack(order + "HHI", entry[:8])
        size = sizes.get(kind, 1) * number
        raw = entry[8:12] if size <= 4 else data[struct.unpack(order + "I", entry[8:12])[0]:][:size]
        if kind == 3:
            tags[tag] = struct.unpack(order + "H", raw[:2])[0]
        elif kind == 4:
            tags[tag] = struct.unpack(order + "I", raw[:4])[0]
        elif kind == 5:
            numerator, denominator = struct.unpack(order + "II", raw[:8])
            tags[tag] = numerator / denominator
    return {"width_px": tags.get(256), "height_px": tags.get(257), "dpi": tags.get(282), "compression": {5: "LZW"}.get(tags.get(259), tags.get(259))}


def _pdf_chunks(data: bytes, *, skip_xref: bool = False) -> list[bytes]:
    chunks = [data]
    for match in re.finditer(rb"stream\r?\n(.*?)\r?\nendstream", data, re.S):
        if skip_xref:
            head = data[max(0, match.start() - 300):match.start()]
            if b"/Type /XRef" in head[head.rfind(b" obj"):]:
                continue  # byte offsets of objects: they shift when the creation date compresses differently
        try:
            chunks.append(zlib.decompress(match.group(1)))
        except zlib.error:
            pass
    return chunks


def pdf_info(path: Path) -> dict:
    chunks = _pdf_chunks(path.read_bytes())
    boxes = [re.search(rb"/MediaBox \[\s*([\d.\s-]+)\]", chunk) for chunk in chunks]
    box = next((m for m in boxes if m), None)
    fonts = sorted({re.sub(r"^[A-Z]{6}\+", "", name.decode("latin1")) for chunk in chunks for name in re.findall(rb"/(?:BaseFont|FontName) ?/([A-Za-z0-9+._-]+)", chunk)})
    size = [float(v) for v in box.group(1).split()] if box else None
    return {"width_pt": size[2] if size else None, "height_pt": size[3] if size else None, "embedded_fonts": fonts}


def pdf_content_hash(path: Path) -> str:
    """Hash of a PDF with creation dates, lengths and offsets (which shift with the date) removed."""

    parts = []
    for chunk in _pdf_chunks(path.read_bytes(), skip_xref=True)[1:]:
        chunk = re.sub(rb"/(?:CreationDate|ModDate) ?\([^)]*\)", b"", chunk)
        chunk = re.sub(rb"\d+ 0 obj\n\s+\d+\nendobj", b"", chunk)
        chunk = re.sub(rb"/Length \d+", b"", chunk)
        parts.append(hashlib.sha256(chunk).digest())
    return hashlib.sha256(b"".join(sorted(parts))).hexdigest()


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


# ---------------------------------------------------------------------- comparison


def pixel_compare(legacy: Path, candidate: Path, diff_path: Path | None) -> dict:
    la, ca = read_png(legacy), read_png(candidate)
    if la[:3] != ca[:3]:
        return {"same_size": False, "legacy_px": [la[0], la[1]], "candidate_px": [ca[0], ca[1]], "identical": False}
    width, height, channels, lrows = la
    crows = ca[3]
    differing, worst = 0, 0
    out = []
    for lrow, crow in zip(lrows, crows):
        line = bytearray()
        for x in range(width):
            a, b = lrow[x * channels:x * channels + 3], crow[x * channels:x * channels + 3]
            delta = max(abs(a[0] - b[0]), abs(a[1] - b[1]), abs(a[2] - b[2]))
            if delta:
                differing += 1
                worst = max(worst, delta)
                line += b"\xe0\x00\x00"
            else:
                grey = (a[0] + a[1] + a[2]) // 3
                line += bytes([128 + grey // 2] * 3)  # faded context
        out.append(bytes(line))
    if differing and diff_path is not None:
        write_png(diff_path, width, height, out)
    return {"same_size": True, "identical": differing == 0, "differing_pixels": differing, "fraction": differing / (width * height), "max_channel_delta": worst}


def load_run(root: Path) -> dict[str, dict]:
    figures: dict[str, dict] = {}
    for manifest_path in sorted(root.rglob("figure_manifest.json")):
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        base = manifest_path.parent.relative_to(root).as_posix()
        for figure in manifest["figures"]:
            key = f"{base}/{figure['id']}"
            figures[key] = {"manifest_path": manifest_path, "manifest": manifest, "figure": figure, "dir": manifest_path.parent}
    return figures


def describe(entry: dict) -> dict:
    figure, folder = entry["figure"], entry["dir"]
    files = {item["format"]: item for item in figure["files"]}
    actual = {}
    for fmt, item in files.items():
        path = folder / item["path"]
        if fmt == "png":
            width, height, _channels, _rows = read_png(path)
            actual["png"] = {"width_px": width, "height_px": height, "dpi": item.get("dpi")}
        elif fmt == "tiff":
            actual["tiff"] = tiff_info(path)
        elif fmt == "pdf":
            actual["pdf"] = pdf_info(path)
    return {
        "width_mm": round(figure["width_mm"], 2), "height_mm": round(figure["height_mm"], 2),
        "formats": sorted(files), "kind": figure.get("kind"),
        "manifest": {fmt: {k: item.get(k) for k in ("width_px", "height_px", "dpi", "width_pt", "height_pt", "compression")} for fmt, item in files.items()},
        "actual": actual,
    }


def edge_ink(path: Path) -> int:
    """Number of non-white pixels on the outermost one-pixel frame of a white-background PNG.

    Text or marks cut off by the figure edge leave ink on the frame, so a non-zero count is a clipping signal."""

    width, height, channels, rows = read_png(path)
    def inked(row: bytes, x: int) -> bool:
        return min(row[x * channels:x * channels + 3]) < 250
    count = sum(inked(rows[0], x) for x in range(width)) + sum(inked(rows[-1], x) for x in range(width))
    count += sum(inked(row, 0) + inked(row, width - 1) for row in rows[1:-1])
    return count


def pdf_words(path: Path) -> set[str] | None:
    """Words drawn in a PDF (poppler's pdftotext), or None when pdftotext is not installed."""

    try:
        result = subprocess.run(["pdftotext", "-layout", str(path), "-"], capture_output=True, text=True, check=False)
    except OSError:
        return None
    return set(re.findall(r"\S+", result.stdout)) if result.returncode == 0 else None


LEGEND_TITLES = {"ORA dotplot": ["Genes"], "GSEA dotplot": ["Gene set size", "Adjusted P"]}


def _source_ids(entry: dict) -> list[str] | None:
    for item in entry["figure"].get("source_data", []):
        if item["path"].endswith("_source_data.tsv") and (entry["dir"] / item["path"]).is_file():
            lines = (entry["dir"] / item["path"]).read_text(encoding="utf-8").splitlines()
            header = lines[0].split("\t")
            return [line.split("\t")[header.index("ID")] for line in lines[1:]] if "ID" in header else None
    return None


def term_check(entry: dict) -> dict | None:
    """For enrichment dotplots: the term IDs in drawing order (source data) and, for standard figures, the order the manifest says was drawn."""

    ids = _source_ids(entry)
    if ids is None:
        return None
    drawn = entry["figure"].get("term_order")
    return {"ids": ids, "drawn_matches_source": None if drawn is None else list(drawn) == ids}


def legend_check(entry: dict) -> dict | None:
    """Legend completeness: every legend title and every label word of the plotted terms is drawn in the PDF."""

    if entry["figure"].get("kind") not in LEGEND_TITLES:
        return None
    pdf = next((item["path"] for item in entry["figure"]["files"] if item["format"] == "pdf"), None)
    words = pdf_words(entry["dir"] / pdf) if pdf else None
    if words is None:
        return {"verified": False}
    labels = set()
    for item in entry["figure"].get("source_data", []):
        if item["path"].endswith("_source_data.tsv") and (entry["dir"] / item["path"]).is_file():
            lines = (entry["dir"] / item["path"]).read_text(encoding="utf-8").splitlines()
            header = lines[0].split("\t")
            if "display_label" in header:
                for line in lines[1:]:
                    labels.update(line.split("\t")[header.index("display_label")].split())
    missing = sorted(w for title in LEGEND_TITLES[entry["figure"]["kind"]] for w in title.split() if w not in words)
    missing += sorted(w for w in labels if w not in words and not any(w in drawn for drawn in words))
    return {"verified": True, "complete": not missing, "missing": missing}


def run_context(entry: dict) -> dict:
    manifest = entry["manifest"]
    style = manifest.get("figure_style") or {}
    font = manifest.get("font") or {}
    standard_font = style.get("font") or {}
    return {
        "profile": manifest["profile"]["name"],
        "requested_preset": style.get("requested_preset", "legacy (no figure_style recorded)"),
        "preset": style.get("preset", "legacy"),
        "standard_figures_drawn": style.get("standard_figures_drawn"),
        "other_figure_types": style.get("other_figure_types"),
        "font_requested": font.get("requested"), "font_family": font.get("family"), "font_fallback_used": font.get("fallback_used"),
        "font_file": font.get("file"), "font_sha256": font.get("sha256"),
        "standard_font_family": standard_font.get("family"), "standard_font_file": standard_font.get("file"),
        "standard_font_sha256": standard_font.get("sha256"), "standard_font_verified": standard_font.get("verified"),
        "standard_font_embedded_in_probe": standard_font.get("embedded_fonts"),
        "standard_font_license": (standard_font.get("bundled") or {}).get("license"),
        "r_version": manifest.get("r_version"),
    }


def _font_key(name: str) -> str:
    return re.sub(r"[^a-z0-9]", "", name.lower())


def font_verdict(entry: dict, described: dict) -> dict:
    """The font a figure claims (its own record, else the run's) against the font names embedded in its PDF."""

    figure_style = entry["figure"].get("figure_style") or {}
    declared = figure_style.get("font") or (entry["manifest"].get("font") or {}).get("family")
    embedded = (described["actual"].get("pdf") or {}).get("embedded_fonts") or []
    verified = bool(declared) and bool(embedded) and all(_font_key(name).startswith(_font_key(declared)) for name in embedded)
    return {"declared": declared, "preset": figure_style.get("preset", "legacy"), "embedded_in_pdf": embedded, "pdf_matches_declared": verified,
            "sha256": figure_style.get("font_sha256") or (entry["manifest"].get("font") or {}).get("sha256")}


def side_by_side(legacy: Path, candidate: Path, out: Path, gap: int = 24) -> None:
    """One PNG with both previews at their true relative scale (both are 150 dpi previews)."""

    parts = []
    for path in (legacy, candidate):
        width, height, channels, rows = read_png(path)
        parts.append((width, height, [b"".join(row[x * channels:x * channels + 3] for x in range(width)) if channels == 4 else row for row in rows]))
    height = max(p[1] for p in parts)
    total = sum(p[0] for p in parts) + gap
    blank = b"\xff\xff\xff"
    rows = []
    for y in range(height):
        line = []
        for index, (width, h, data) in enumerate(parts):
            line.append(data[y] if y < h else blank * width)
            if index == 0:
                line.append(b"\x90\x90\x90" * gap if gap else b"")
        rows.append(b"".join(line))
    write_png(out, total, height, rows)


def source_data_status(legacy: dict, candidate: dict) -> dict:
    status = {}
    for name, entry in (("legacy", legacy), ("candidate", candidate)):
        status[name] = {item["path"]: sha256(entry["dir"] / item["path"]) for item in entry["figure"].get("source_data", []) if (entry["dir"] / item["path"]).is_file()}
    return {"identical": status["legacy"] == status["candidate"], "files": sorted(status["legacy"])}


def compare(legacy_root: Path, candidate_root: Path, out: Path) -> dict:
    legacy, candidate = load_run(legacy_root), load_run(candidate_root)
    assets = out / "assets"
    if assets.exists():
        shutil.rmtree(assets)
    assets.mkdir(parents=True)
    report = {"legacy_root": str(legacy_root), "candidate_root": str(candidate_root),
              "only_in_legacy": sorted(set(legacy) - set(candidate)), "only_in_candidate": sorted(set(candidate) - set(legacy)), "figures": []}
    report["only_in_candidate_kinds"] = {key: str(candidate[key]["figure"].get("kind") or "") for key in report["only_in_candidate"]}
    if legacy:
        report["legacy_context"] = run_context(next(iter(legacy.values())))
    if candidate:
        report["candidate_context"] = run_context(next(iter(candidate.values())))
    for key in sorted(set(legacy) & set(candidate)):
        a, b = legacy[key], candidate[key]
        record = {"figure": key, "legacy": describe(a), "candidate": describe(b)}
        stem = key.replace("/", "__")
        copies = {}
        for label, entry in (("legacy", a), ("candidate", b)):
            png = entry["dir"] / next(item["path"] for item in entry["figure"]["files"] if item["format"] == "png")
            copies[label] = assets / f"{stem}.{label}.png"
            shutil.copyfile(png, copies[label])
            record[f"{label}_png"] = copies[label].name
        diff_path = assets / f"{stem}.diff.png"
        record["pixels"] = pixel_compare(copies["legacy"], copies["candidate"], diff_path)
        side_path = assets / f"{stem}.side_by_side.png"
        side_by_side(copies["legacy"], copies["candidate"], side_path)
        record["side_by_side_png"] = side_path.name
        record["font"] = {"legacy": font_verdict(a, record["legacy"]), "candidate": font_verdict(b, record["candidate"])}
        record["clipping"] = {"legacy_edge_ink_px": edge_ink(copies["legacy"]), "candidate_edge_ink_px": edge_ink(copies["candidate"])}
        terms_a, terms_b = term_check(a), term_check(b)
        record["terms"] = None if terms_a is None or terms_b is None else {"order_identical": terms_a["ids"] == terms_b["ids"], "terms": len(terms_b["ids"]), "drawn_matches_source": terms_b["drawn_matches_source"]}
        record["legend"] = {"legacy": legend_check(a), "candidate": legend_check(b)} if legend_check(b) is not None else None
        record["diff_png"] = diff_path.name if diff_path.exists() else None
        formats = {item["format"]: item for item in a["figure"]["files"]}
        other = {item["format"]: item for item in b["figure"]["files"]}
        # A format absent from both runs is not a difference; absent from only one is.
        record["tiff_identical"] = ("tiff" not in formats and "tiff" not in other) or ("tiff" in formats and "tiff" in other and sha256(a["dir"] / formats["tiff"]["path"]) == sha256(b["dir"] / other["tiff"]["path"]))
        record["pdf_content_identical"] = ("pdf" not in formats and "pdf" not in other) or ("pdf" in formats and "pdf" in other and pdf_content_hash(a["dir"] / formats["pdf"]["path"]) == pdf_content_hash(b["dir"] / other["pdf"]["path"]))
        record["same_formats"] = record["legacy"]["formats"] == record["candidate"]["formats"]
        record["same_size_mm"] = (record["legacy"]["width_mm"], record["legacy"]["height_mm"]) == (record["candidate"]["width_mm"], record["candidate"]["height_mm"])
        record["source_data"] = source_data_status(a, b)
        same_files = record["pixels"]["identical"] and record["tiff_identical"] and record["pdf_content_identical"] and record["same_formats"] and record["source_data"]["identical"]
        # A different physical size or format set is a layout change even when pixels happen to match.
        record["status"] = "layout differs" if not (record["same_size_mm"] and record["same_formats"]) else ("identical" if same_files else "appearance differs")
        report["figures"].append(record)
    counts: dict[str, int] = {}
    for record in report["figures"]:
        counts[record["status"]] = counts.get(record["status"], 0) + 1
    report["summary"] = {"matched_figures": len(report["figures"]), "by_status": counts, "unmatched": len(report["only_in_legacy"]) + len(report["only_in_candidate"])}
    return report


# ---------------------------------------------------------------------- HTML


def _dims(side: dict) -> str:
    actual = side["actual"]
    parts = [f"{side['width_mm']} × {side['height_mm']} mm"]
    if "png" in actual:
        parts.append(f"PNG {actual['png']['width_px']}×{actual['png']['height_px']} px @ {actual['png']['dpi']} dpi")
    if "tiff" in actual:
        parts.append(f"TIFF {actual['tiff']['width_px']}×{actual['tiff']['height_px']} px @ {actual['tiff']['dpi']:.0f} dpi, {actual['tiff']['compression']}")
    if "pdf" in actual:
        parts.append(f"PDF {actual['pdf']['width_pt']:g}×{actual['pdf']['height_pt']:g} pt; fonts: {', '.join(actual['pdf']['embedded_fonts']) or 'none'}")
    return "<br>".join(html.escape(p) for p in parts)


def _fonts(fonts: dict) -> str:
    def side(label: str, item: dict) -> str:
        return (f"{label}: declared {html.escape(str(item['declared']))} ({html.escape(item['preset'])}); embedded in PDF: {html.escape(', '.join(item['embedded_in_pdf']) or 'none')}; "
                f"PDF matches declared font: {item['pdf_matches_declared']}")

    return side("legacy", fonts["legacy"]) + "<br>" + side("candidate", fonts["candidate"])


def _checks(record: dict) -> str:
    parts = [f"edge ink (clipping signal), legacy/candidate: {record['clipping']['legacy_edge_ink_px']}/{record['clipping']['candidate_edge_ink_px']} px"]
    if record.get("terms"):
        terms = record["terms"]
        parts.append(f"{terms['terms']} terms; term order identical to legacy: {terms['order_identical']}; drawn order = source order: {terms['drawn_matches_source']}")
    if record.get("legend"):
        legend = record["legend"]["candidate"]
        parts.append("legend/labels in PDF: " + ("not verified (pdftotext missing)" if not legend["verified"] else f"complete={legend['complete']} missing={legend['missing']}"))
    return html.escape("; ".join(parts))


def render_html(report: dict) -> str:
    def context_table() -> str:
        keys = ["profile", "requested_preset", "preset", "standard_figures_drawn", "other_figure_types", "font_requested", "font_family", "font_fallback_used", "font_file", "font_sha256",
                "standard_font_family", "standard_font_file", "standard_font_sha256", "standard_font_verified", "standard_font_embedded_in_probe", "standard_font_license", "r_version"]
        legacy, candidate = report.get("legacy_context", {}), report.get("candidate_context", {})
        rows = "".join(f"<tr><th>{k}</th><td>{html.escape(str(legacy.get(k)))}</td><td>{html.escape(str(candidate.get(k)))}</td></tr>" for k in keys)
        return f"<table class='ctx'><tr><th></th><th>legacy run</th><th>candidate run</th></tr>{rows}</table>"

    body = []
    for record in report["figures"]:
        diff = f"<img src='assets/{record['diff_png']}'>" if record["diff_png"] else "<p class='ok'>no pixel differences</p>"
        pixels = record["pixels"]
        detail = "identical" if pixels["identical"] else (f"{pixels.get('differing_pixels')} px differ ({pixels.get('fraction', 0):.2%}), max channel delta {pixels.get('max_channel_delta')}" if pixels["same_size"] else f"different pixel size: {pixels['legacy_px']} vs {pixels['candidate_px']}")
        body.append(
            f"<section class='fig {record['status'].replace(' ', '-')}'><h3>{html.escape(record['figure'])} <span class='badge'>{html.escape(record['status'])}</span></h3>"
            f"<p class='meta'>PNG pixels: {html.escape(detail)}; TIFF bytes identical: {record['tiff_identical']}; PDF content identical: {record['pdf_content_identical']}; "
            f"same formats: {record['same_formats']}; source data identical: {record['source_data']['identical']} ({len(record['source_data']['files'])} table(s))</p>"
            f"<p class='meta'>{_fonts(record['font'])}</p><p class='meta'>{_checks(record)}</p>"
            f"<div class='row'><figure><img src='assets/{record['legacy_png']}' style='max-width:{record['legacy']['width_mm']}mm'><figcaption>legacy<br>{_dims(record['legacy'])}</figcaption></figure>"
            f"<figure><img src='assets/{record['candidate_png']}' style='max-width:{record['candidate']['width_mm']}mm'><figcaption>candidate<br>{_dims(record['candidate'])}</figcaption></figure>"
            f"<figure>{diff}<figcaption>difference (red = differing pixel)</figcaption></figure></div>"
            f"<p class='meta'>side-by-side file: <code>assets/{record['side_by_side_png']}</code> (both previews at true relative scale)</p></section>")
    summary = report["summary"]
    unmatched = "".join(f"<li>{html.escape(k)}</li>" for k in report["only_in_legacy"] + report["only_in_candidate"])
    return (
        "<!doctype html><html lang='en'><head><meta charset='utf-8'><title>nf-rna figure comparison</title><style>"
        "body{font-family:sans-serif;margin:1.5rem;max-width:1500px}table.ctx{border-collapse:collapse;font-size:12px}table.ctx td,table.ctx th{border:1px solid #bbb;padding:2px 6px;text-align:left;vertical-align:top}"
        ".row{display:flex;gap:12px;align-items:flex-start}.row figure{margin:0;flex:1;min-width:0}.row img{width:100%;border:1px solid #ccc}figcaption{font-size:11px}"
        ".badge{font-size:12px;padding:1px 6px;border-radius:4px;background:#ddd}.identical .badge{background:#cfe8cf}.appearance-differs .badge,.layout-differs .badge{background:#f3d0a8}.meta{font-size:12px}.ok{color:#2a6e2a;font-size:12px}"
        "</style></head><body><h1>nf-rna figure comparison</h1>"
        f"<p>Legacy root: <code>{html.escape(report['legacy_root'])}</code><br>Candidate root: <code>{html.escape(report['candidate_root'])}</code></p>"
        f"<p><b>{summary['matched_figures']}</b> figures matched; status counts: {html.escape(json.dumps(summary['by_status']))}; unmatched: {summary['unmatched']}</p>"
        f"{('<ul>' + unmatched + '</ul>') if unmatched else ''}{context_table()}{''.join(body)}</body></html>")


def gate_failures(report: dict, kinds: list[str]) -> list[str]:
    """Why the report fails the "only these figure kinds may change" gate (empty list: it passes)."""

    allowed = {kind.lower() for kind in kinds}
    failures = [f"unmatched figure: {key}" for key in report["only_in_legacy"] + report["only_in_candidate"]]
    for record in report["figures"]:
        kind = str(record["candidate"].get("kind") or "").lower()
        if not record["source_data"]["identical"]:
            failures.append(f"{record['figure']}: source data differs")
        if kind in allowed:
            if not record["font"]["candidate"]["pdf_matches_declared"]:
                failures.append(f"{record['figure']}: embedded PDF font does not match the declared font")
            if record["clipping"]["candidate_edge_ink_px"]:
                failures.append(f"{record['figure']}: ink on the figure edge (possible clipping)")
            terms = record.get("terms")
            if terms and (not terms["order_identical"] or terms["drawn_matches_source"] is False):
                failures.append(f"{record['figure']}: term order differs")
            legend = (record.get("legend") or {}).get("candidate")
            if legend and legend["verified"] and not legend["complete"]:
                failures.append(f"{record['figure']}: legend/labels incomplete in the PDF: {legend['missing']}")
        elif record["status"] != "identical":
            failures.append(f"{record['figure']}: unexpected change ({record['status']})")
    return failures


def expectation_failures(report: dict, expected: dict) -> tuple[list[str], list[str]]:
    """Gate a release-to-release comparison on declared, reasoned changes.

    Returns (failures, notes).  Nothing outside the declaration may change:
    * a figure only in the legacy run fails (removed figure);
    * a figure only in the candidate passes only if a ``new_figures`` entry matches its ID (regex) and kind;
    * a matched figure of a ``changed_kinds`` kind may differ in content, but must still embed its declared
      font, keep ink off the edge and keep a complete legend;
    * any other matched figure must be identical, or differ only by ``added_formats`` (identical pixels,
      PDF content, TIFF, size and source data).
    """

    changed = {kind.lower(): reason for kind, reason in (expected.get("changed_kinds") or {}).items()}
    added = expected.get("added_formats") or {}
    new_rules = expected.get("new_figures") or []
    failures = [f"removed figure: {key}" for key in report["only_in_legacy"]]
    notes: list[str] = []
    for key in report["only_in_candidate"]:
        kind = report.get("only_in_candidate_kinds", {}).get(key, "")
        rule = next((r for r in new_rules if re.search(r["pattern"], key) and r.get("kind", kind) == kind), None)
        if rule is None:
            failures.append(f"unexpected new figure: {key} ({kind})")
        else:
            notes.append(f"new figure {key} ({kind}): {rule['reason']}")
    for record in report["figures"]:
        key, kind = record["figure"], str(record["candidate"].get("kind") or "")
        if record["status"] == "identical":
            continue
        if kind.lower() in changed:
            notes.append(f"changed {key} ({kind}): {changed[kind.lower()]}")
            if not record["font"]["candidate"]["pdf_matches_declared"]:
                failures.append(f"{key}: embedded PDF font does not match the declared font")
            if record["clipping"]["candidate_edge_ink_px"]:
                failures.append(f"{key}: ink on the figure edge (possible clipping)")
            legend = (record.get("legend") or {}).get("candidate")
            if legend and legend["verified"] and not legend["complete"]:
                failures.append(f"{key}: legend/labels incomplete in the PDF: {legend['missing']}")
            continue
        legacy_formats, candidate_formats = set(record["legacy"]["formats"]), set(record["candidate"]["formats"])
        extra = candidate_formats - legacy_formats
        format_only = (record["pixels"]["identical"] and record.get("pdf_content_identical") is True and record.get("tiff_identical") is not False
                       and record["source_data"]["identical"] and record.get("same_size_mm") and legacy_formats <= candidate_formats and extra <= set(added))
        if format_only and extra:
            notes.append(f"formats added to {key} ({kind}): " + "; ".join(f"{name}: {added[name]}" for name in sorted(extra)))
        else:
            failures.append(f"{key} ({kind}): unexpected change ({record['status']})")
    return failures, notes


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--legacy", type=Path, required=True, help="root of the legacy-preset run (searched for figure_manifest.json)")
    parser.add_argument("--candidate", type=Path, required=True, help="root of the candidate run")
    parser.add_argument("--out", type=Path, required=True, help="report directory (created; its assets/ folder is replaced)")
    parser.add_argument("--require-identical", action="store_true", help="exit 1 unless every figure is identical and nothing is unmatched")
    parser.add_argument("--require-identical-except", metavar="KIND", nargs="+", help="exit 1 unless every figure whose manifest kind is NOT listed is identical, nothing is unmatched, "
                        "every listed figure carries identical source data and an embedded PDF font that matches its declared font")
    parser.add_argument("--expected-changes", type=Path, metavar="FILE", help="JSON declaration of intended changes between two releases (see expectations/); "
                        "exit 1 on any undeclared difference or removed figure")
    args = parser.parse_args(argv)
    args.out.mkdir(parents=True, exist_ok=True)
    report = compare(args.legacy.resolve(), args.candidate.resolve(), args.out)
    (args.out / "report.json").write_text(json.dumps(report, indent=1, sort_keys=True, default=str) + "\n", encoding="utf-8")
    (args.out / "report.html").write_text(render_html(report), encoding="utf-8")
    if args.require_identical_except:
        for line in gate_failures(report, args.require_identical_except):
            print("GATE FAILURE:", line)
    summary = report["summary"]
    print(f"matched {summary['matched_figures']} figures: {summary['by_status']}; unmatched {summary['unmatched']}; report: {args.out / 'report.html'}")
    if args.expected_changes:
        failures, notes = expectation_failures(report, json.loads(args.expected_changes.read_text(encoding="utf-8")))
        for line in notes:
            print("EXPECTED:", line)
        for line in failures:
            print("GATE FAILURE:", line)
        print(f"expected-changes gate: {'PASS' if not failures else 'FAIL'} ({len(notes)} declared changes, {len(failures)} failures)")
        if failures:
            return 1
    if args.require_identical_except:
        return 1 if gate_failures(report, args.require_identical_except) else 0
    return 1 if args.require_identical and (summary["by_status"].keys() - {"identical"} or summary["unmatched"]) else 0


if __name__ == "__main__":
    sys.exit(main())
