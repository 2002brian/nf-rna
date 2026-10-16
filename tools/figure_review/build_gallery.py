#!/usr/bin/env python3
"""Build ONE self-contained HTML gallery of legacy and standard figures (stdlib only; review aid).

Each ``--section NAME LEGACY_RUN STANDARD_RUN`` pairs two runs of the same project (for example produced by
``generate_runs.py`` with ``--figures`` unset / ``{preset: standard}``).  ``compare_runs.py`` matches the figures through
their manifests; this tool embeds the PNG previews as data URIs and prints, per figure: the standard and the legacy
rendering drawn side by side at a COMMON physical scale, physical dimensions, dpi, declared and embedded fonts (with the
bundled-font SHA-256), the palette recorded in the manifest, source-data identity, term-order identity (enrichment) and
automatic visual warnings.  Sections are labelled synthetic or existing example data.

    build_gallery.py --out gallery.html --work WORKDIR \\
        --section "Synthetic 12-sample project" LEGACY STANDARD --synthetic "counts simulated by make_demo_project.py"
"""

from __future__ import annotations

import argparse
import base64
import html
import json
import re
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from profile_audit import text_sizes  # noqa: E402  (exact PDF font sizes)
HEX = re.compile(r"#[0-9A-Fa-f]{6}\b")
ORDER = [
    ("PCA", re.compile(r"(^|/)pca$")), ("PCA by batch", re.compile(r"(^|/)pca_by_batch$")),
    ("Raw-count boxplot", re.compile(r"boxplot_log2_raw_counts$")), ("VST boxplot", re.compile(r"boxplot_vst$")),
    ("Sample correlation heatmap", re.compile(r"sample_correlation$")), ("Volcano", re.compile(r"volcano")),
    ("DEG heatmap", re.compile(r"(^|/)heatmap")), ("GO dotplot (ORA)", re.compile(r"^enrichment/go/")),
    ("GSEA dotplot (GO)", re.compile(r"^enrichment/gsea_go/")), ("KEGG ORA dotplot", re.compile(r"^enrichment/kegg/")),
    ("KEGG GSEA dotplot", re.compile(r"^enrichment/gsea_kegg/")),
]


def label_of(figure: str) -> tuple[int, str]:
    for index, (label, pattern) in enumerate(ORDER):
        if pattern.search(figure):
            return index, label
    return len(ORDER), "Other"


def data_uri(path: Path) -> str:
    return "data:image/png;base64," + base64.b64encode(path.read_bytes()).decode("ascii")


def manifest_figure(run: Path, figure: str) -> dict | None:
    """The manifest record of ``figure`` (a path relative to the run, as compare_runs prints it)."""

    parts = figure.split("/")
    for cut in range(len(parts) - 1, 0, -1):
        candidate = run.joinpath(*parts[:cut]) / "figure_manifest.json"
        if candidate.exists():
            data = json.loads(candidate.read_text(encoding="utf-8"))
            rest = "/".join(parts[cut:])
            for record in data["figures"]:
                if record["id"] == rest:
                    return {"record": record, "manifest": data, "dir": candidate.parent}
    return None


def palette_of(found: dict | None) -> dict:
    if not found:
        return {"swatches": [], "note": "no manifest record"}
    record = found["record"]
    swatches: list[tuple[str, str]] = []

    def walk(node, path):
        if isinstance(node, dict):
            for key, value in node.items():
                walk(value, f"{path}.{key}" if path else key)
        elif isinstance(node, list):
            for item in node:
                walk(item, path)
        elif isinstance(node, str) and HEX.fullmatch(node.strip()):
            swatches.append((path, node.strip().upper()))

    for key in ("encoding", "colour_scale", "colours", "palette", "group_colours", "enrichment_colours", "mark_colour", "annotation", "legend"):
        if key in record:
            walk(record[key], key)
    names = []
    encoding = record.get("encoding")
    if isinstance(encoding, dict) and encoding.get("palette"):
        names.append(str(encoding["palette"]))
    if record.get("palette"):
        names.append(str(record["palette"]))
    return {"swatches": swatches, "names": names}


def size_note(width, height, profile: dict | None) -> str:
    return f"{width:.1f} x {height:.1f} mm"


def warnings_for(figure: dict, found_std: dict | None, palette: dict) -> list[str]:
    out = []
    cand, leg = figure["candidate"], figure["legacy"]
    clip = figure.get("clipping") or {}
    if clip.get("candidate_edge_ink_px"):
        out.append(f"standard PNG has {clip['candidate_edge_ink_px']} px of ink on its border (possible clipping)")
    if clip.get("legacy_edge_ink_px"):
        out.append(f"legacy PNG has {clip['legacy_edge_ink_px']} px of ink on its border (legacy behaviour, informational)")
    legend = (figure.get("legend") or {}).get("candidate")
    if legend and legend.get("verified") and not legend.get("complete"):
        out.append("standard legend incomplete in the PDF: " + ", ".join(legend.get("missing", [])))
    terms = figure.get("terms")
    if terms and not (terms.get("order_identical") and terms.get("drawn_matches_source")):
        out.append("term order differs from legacy or from the source data")
    font = figure["font"]["candidate"]
    if not font.get("pdf_matches_declared"):
        out.append(f"standard PDF embeds {font.get('embedded_in_pdf')} but declares {font.get('declared')}")
    if not figure["source_data"]["identical"]:
        out.append("source-data tables differ between legacy and standard")
    if found_std:
        profile = found_std["manifest"].get("profile") or {}
        if profile.get("max_width_mm") and cand["width_mm"] > profile["max_width_mm"] + 0.05:
            out.append(f"standard width {cand['width_mm']:.1f} mm exceeds the profile maximum {profile['max_width_mm']} mm")
        if profile.get("max_height_mm") and cand["height_mm"] > profile["max_height_mm"] + 0.05:
            out.append(f"standard height {cand['height_mm']:.1f} mm exceeds the profile maximum {profile['max_height_mm']} mm")
    names = palette.get("names", [])
    record = found_std["record"] if found_std else {}
    encoding = record.get("encoding") if isinstance(record.get("encoding"), dict) else {}
    levels = len(encoding.get("colours") or {})
    if any(name.startswith("HCL") for name in names) and found_std and found_std["manifest"].get("figure_style", {}).get("preset") == "standard":
        out.append(f"{levels} levels exceed the 8 palette slots: HCL 'Dark 3' is used for ALL levels, so approved group colours 1-3 do not appear")
    if profile_obj := (found_std["manifest"].get("profile") if found_std else None):
        if profile_obj.get("single_col_mm") is None and cand["width_mm"] >= profile_obj["double_col_mm"] - 0.5 and figure["candidate"]["kind"] in ("PCA", "volcano"):
            out.append(f"profile {profile_obj['name']} has no single column: this {figure['candidate']['kind']} fills the full {profile_obj['double_col_mm']} mm width")
        pdf_files = [item for item in record.get("files", []) if item["format"] == "pdf"]
        if pdf_files and found_std:
            sizes = text_sizes(found_std["dir"] / pdf_files[0]["path"])
            if sizes.get("runs") and sizes["implied_min_pt"] < profile_obj["min_pt"] - 0.05:
                out.append(f"smallest text is {sizes['implied_min_pt']} pt, below the {profile_obj['name']} profile minimum of {profile_obj['min_pt']} pt")
    if cand["width_mm"] > leg["width_mm"] + 0.05:
        out.append("standard figure is wider than the legacy one")
    return out


def card(figure: dict, label: str, compare_dir: Path, legacy_run: Path, standard_run: Path, scale: float) -> tuple[str, list[str]]:
    found_std = manifest_figure(standard_run, figure["figure"])
    palette = palette_of(found_std)
    warns = warnings_for(figure, found_std, palette)
    cand, leg = figure["candidate"], figure["legacy"]
    assets = compare_dir / "assets"
    legacy_png, standard_png = assets / figure["legacy_png"], assets / figure["candidate_png"]
    font_c, font_l = figure["font"]["candidate"], figure["font"]["legacy"]
    sha = (font_c.get("sha256") or "")[:12]
    swatches = "".join(f'<span class="sw" style="background:{hexv}" title="{html.escape(path)} {hexv}"></span>' for path, hexv in palette["swatches"][:24])
    terms = figure.get("terms")
    term_txt = f"{terms['terms']} terms; order identical to legacy: {terms['order_identical']}; drawn = source order: {terms['drawn_matches_source']}" if terms else "n/a (not an enrichment dotplot)"
    legend = (figure.get("legend") or {}).get("candidate")
    legend_txt = "n/a" if not legend else ("complete" if legend.get("complete") else "INCOMPLETE") + ("" if legend.get("verified") else " (not verifiable)")
    def smallest(found):
        files = [item for item in (found["record"].get("files", []) if found else []) if item["format"] == "pdf"]
        sizes = text_sizes(found["dir"] / files[0]["path"]) if files else {}
        return f"{sizes['implied_min_pt']} pt (largest {sizes['implied_max_pt']} pt)" if sizes.get("runs") else "n/a"

    rows = [
        ("Physical size (legacy)", f"{leg['width_mm']:.1f} x {leg['height_mm']:.1f} mm"),
        ("Physical size (standard)", f"{cand['width_mm']:.1f} x {cand['height_mm']:.1f} mm"),
        ("Raster dpi (legacy / standard)", f"{leg['actual']['png']['dpi']} png, {leg['actual']['tiff']['dpi']:.0f} tiff / {cand['actual']['png']['dpi']} png, {cand['actual']['tiff']['dpi']:.0f} tiff"),
        ("Font (legacy)", f"{font_l.get('declared')} (PDF embeds {', '.join(font_l.get('embedded_in_pdf') or [])})"),
        ("Font (standard)", f"{font_c.get('declared')} (PDF embeds {', '.join(font_c.get('embedded_in_pdf') or [])}); sha256 {sha or 'n/a'}"),
        ("Smallest text (legacy / standard PDF)", f"{smallest(manifest_figure(legacy_run, figure['figure']))} / {smallest(found_std)}"),
        ("Palette (standard)", html.escape("; ".join(palette.get("names", [])) or "see swatches") + " " + swatches),
        ("Source data identical", "yes" if figure["source_data"]["identical"] else "NO"),
        ("Term order", term_txt),
        ("Legend (standard)", legend_txt),
        ("Status", figure["status"]),
    ]
    table = "".join(f"<tr><th>{html.escape(a)}</th><td>{b if a == 'Palette (standard)' else html.escape(b)}</td></tr>" for a, b in rows)
    warn_html = "".join(f"<li>{html.escape(w)}</li>" for w in warns) or "<li class='ok'>none</li>"
    images = (
        f'<figure><img style="width:{leg["width_mm"] * scale:.0f}px" src="{data_uri(legacy_png)}" alt="legacy {html.escape(figure["figure"])}"><figcaption>legacy</figcaption></figure>'
        f'<figure><img style="width:{cand["width_mm"] * scale:.0f}px" src="{data_uri(standard_png)}" alt="standard {html.escape(figure["figure"])}"><figcaption>standard</figcaption></figure>'
    )
    flagged = any(not w.startswith("legacy PNG") for w in warns)
    html_card = (f'<article class="card{" flagged" if flagged else ""}"><h3>{html.escape(label)} <code>{html.escape(figure["figure"])}</code></h3>'
                 f'<div class="pair">{images}</div><table>{table}</table><div class="warn"><b>Visual warnings</b><ul>{warn_html}</ul></div></article>')
    return html_card, warns


CSS = """
:root{--bg:#fff;--fg:#1d1d1f;--mut:#5b5b62;--line:#d7d7dc;--card:#fafafa;--warn:#a4161a;--ok:#1b7f3b;--syn:#fff4d6}
@media (prefers-color-scheme:dark){:root{--bg:#161618;--fg:#ececf0;--mut:#a0a0aa;--line:#3a3a40;--card:#1e1e21;--warn:#ff8a8d;--ok:#6fd48d;--syn:#3b3217}}
body{font:14px/1.45 system-ui,sans-serif;margin:0 auto;max-width:1400px;padding:16px;background:var(--bg);color:var(--fg)}
h1{font-size:1.5rem}h2{margin-top:2rem;border-bottom:1px solid var(--line);padding-bottom:.3rem}h3{font-size:1rem;margin:.2rem 0}code{font-size:.8rem;color:var(--mut)}
.syn{background:var(--syn);border:1px solid var(--line);padding:.5rem .8rem;border-radius:6px;margin:.6rem 0}
.card{border:1px solid var(--line);background:var(--card);border-radius:8px;padding:.8rem;margin:1rem 0}.card.flagged{border-color:var(--warn)}
.pair{display:flex;gap:16px;align-items:flex-start;flex-wrap:wrap;background:#fff;padding:8px;border-radius:6px}
.pair figure{margin:0}.pair img{max-width:100%;height:auto;border:1px solid #ddd}.pair figcaption{font-size:.75rem;color:#444;text-align:center}
table{border-collapse:collapse;margin:.6rem 0;width:100%}th{text-align:left;width:14rem;color:var(--mut);font-weight:500;vertical-align:top}td,th{padding:2px 6px;border-bottom:1px solid var(--line)}
.sw{display:inline-block;width:14px;height:14px;border:1px solid #888;margin-right:2px;vertical-align:middle}
.warn ul{margin:.2rem 0 0 1.1rem;padding:0;color:var(--warn)}.warn li.ok{color:var(--ok);list-style:none;margin-left:-1.1rem}
.summary td,.summary th{width:auto;font-size:.85rem}.note{color:var(--mut);font-size:.85rem}
@media (max-width:700px){th{width:auto;display:block}.pair{display:block}}
"""


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--work", type=Path, required=True, help="directory for the compare_runs reports")
    parser.add_argument("--section", nargs=3, action="append", metavar=("NAME", "LEGACY_RUN", "STANDARD_RUN"), required=True)
    parser.add_argument("--synthetic", action="append", default=[], help="per section (same order): data description; use 'EXAMPLE:...' for existing example data")
    parser.add_argument("--scale", type=float, default=4.0, help="screen px per mm for the common physical scale")
    args = parser.parse_args(argv)
    body, summary_rows, total_warn = [], [], 0
    for index, (name, legacy, standard) in enumerate(args.section):
        legacy_run, standard_run = Path(legacy), Path(standard)
        compare_dir = args.work / f"section_{index}"
        subprocess.run([sys.executable, str(HERE / "compare_runs.py"), "--legacy", str(legacy_run), "--candidate", str(standard_run), "--out", str(compare_dir)], capture_output=True, text=True, check=False)
        report = json.loads((compare_dir / "report.json").read_text(encoding="utf-8"))
        note = args.synthetic[index] if index < len(args.synthetic) else "SYNTHETIC: description not supplied"
        kind = "Existing example data" if note.startswith("EXAMPLE:") else "SYNTHETIC DATA"
        ctx = report["candidate_context"]
        body.append(f'<h2>{html.escape(name)}</h2><div class="syn"><b>{kind}.</b> {html.escape(note.removeprefix("EXAMPLE:").strip())}<br>'
                    f'Standard run: preset <b>{html.escape(str(ctx.get("preset")))}</b>, profile <b>{html.escape(str(ctx.get("profile")))}</b>; bundled font {html.escape(str(ctx.get("standard_font_family")))} '
                    f'sha256 <code>{html.escape(str(ctx.get("standard_font_sha256")))}</code>, verified: {ctx.get("standard_font_verified")}. Pairs are drawn at one common physical scale ({args.scale:g} px per mm).</div>')
        figures = sorted(report["figures"], key=lambda f: (label_of(f["figure"])[0], f["figure"]))
        for unmatched in report.get("only_in_legacy", []) + report.get("only_in_candidate", []):
            body.append(f'<p class="warn"><b>Unmatched figure:</b> {html.escape(str(unmatched))}</p>')
        cards = []
        for figure in figures:
            html_card, warns = card(figure, label_of(figure["figure"])[1], compare_dir, legacy_run, standard_run, args.scale)
            cards.append(html_card)
            real = [w for w in warns if not w.startswith("legacy PNG")]
            total_warn += len(real)
            summary_rows.append((name, label_of(figure["figure"])[1], figure["figure"], figure["legacy"], figure["candidate"], figure["source_data"]["identical"], len(real)))
        body.extend(cards)
    rows = "".join(f"<tr><td>{html.escape(n)}</td><td>{html.escape(l)}</td><td><code>{html.escape(f)}</code></td><td>{lg['width_mm']:.1f} x {lg['height_mm']:.1f}</td><td>{st['width_mm']:.1f} x {st['height_mm']:.1f}</td><td>{'yes' if same else 'NO'}</td><td>{w}</td></tr>" for n, l, f, lg, st, same, w in summary_rows)
    summary = (f'<table class="summary"><tr><th>Section</th><th>Figure</th><th>Id</th><th>Legacy mm</th><th>Standard mm</th><th>Source data identical</th><th>Warnings</th></tr>{rows}</table>')
    page = (f'<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>nf-rna figure gallery</title><style>{CSS}</style></head><body>'
            f'<h1>nf-rna figure gallery: standard versus legacy</h1><p class="note">Generated by tools/figure_review/build_gallery.py from figure manifests and the PNG previews. Review aid; nothing here is a result. {len(summary_rows)} figure pairs, {total_warn} warnings (informational legacy-only notes excluded).</p>'
            f'{summary}{"".join(body)}</body></html>')
    args.out.write_text(page, encoding="utf-8")
    print(f"wrote {args.out} ({args.out.stat().st_size / 1e6:.1f} MB), {len(summary_rows)} figure pairs, {total_warn} warnings")
    return 0


if __name__ == "__main__":
    sys.exit(main())
