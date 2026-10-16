#!/usr/bin/env python3
"""Publication-profile audit of legacy and standard runs (stdlib only; review aid).

For every profile given as ``--profile NAME LEGACY_RUN STANDARD_RUN`` the figures of both runs are matched through
their manifests (``compare_runs.py``) and checked against the profile recorded in the standard manifest:

* physical width and height against ``max_width_mm`` / ``max_height_mm``;
* the raster dpi actually written to the TIFF (and PNG preview) against ``raster_dpi`` / ``preview_dpi``;
* the smallest text actually drawn in the PDF (read from the Tm/Tf operators of the PDF content streams) against the profile's ``min_pt``;
* clipping (ink on the PNG border), legend completeness and PDF fonts as reported by ``compare_runs.py``.

    profile_audit.py --out DIR --profile general LEG STD --profile nature LEG STD ...
"""

from __future__ import annotations

import argparse
import json
import re
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent


TEXT_MATRIX = re.compile(rb"([-\d.e]+) ([-\d.e]+) ([-\d.e]+) ([-\d.e]+) [-\d.e]+ [-\d.e]+ Tm\s*/\S+ ([-\d.e]+) Tf")


def text_sizes(pdf: Path) -> dict:
    """Font sizes (pt) of every text run drawn in the PDF, read from its content streams.

    Cairo writes each run as ``a b c d x y Tm /font size Tf``; the size is hypot(a, b) * size, which is exact for horizontal and
    rotated text alike (no bounding-box heuristics).
    """

    import zlib

    sizes = []
    data = pdf.read_bytes()
    for stream in re.finditer(rb"stream\r?\n(.*?)endstream", data, re.S):
        try:
            content = zlib.decompress(stream.group(1))
        except zlib.error:
            continue
        for m in TEXT_MATRIX.finditer(content):
            a, b, size = float(m.group(1)), float(m.group(2)), float(m.group(5))
            sizes.append(round((a * a + b * b) ** 0.5 * size, 2))
    if not sizes:
        return {"runs": 0}
    return {"runs": len(sizes), "implied_min_pt": min(sizes), "implied_max_pt": max(sizes)}


def pdfs(run: Path) -> dict[str, Path]:
    """Figure key (as compare_runs.py prints it) -> PDF path, read from the manifests."""

    found = {}
    for manifest in sorted(run.rglob("figure_manifest.json")):
        base = manifest.parent.relative_to(run).as_posix()
        for figure in json.loads(manifest.read_text(encoding="utf-8"))["figures"]:
            for item in figure["files"]:
                if item["format"] == "pdf":
                    found[f"{base}/{figure['id']}"] = manifest.parent / item["path"]
    return found


def audit_profile(name: str, legacy: Path, standard: Path, out: Path) -> dict:
    compare = out / name / "compare"
    subprocess.run([sys.executable, str(HERE / "compare_runs.py"), "--legacy", str(legacy), "--candidate", str(standard), "--out", str(compare)], capture_output=True, text=True, check=False)
    report = json.loads((compare / "report.json").read_text(encoding="utf-8"))
    first = next(iter(standard.rglob("figure_manifest.json")))
    profile = json.loads(first.read_text(encoding="utf-8"))["profile"]
    rows = []
    for figure in report["figures"]:
        cand, leg = figure["candidate"], figure["legacy"]
        pdf_std, pdf_leg = pdfs(standard).get(figure["figure"]), pdfs(legacy).get(figure["figure"])
        entry = {
            "figure": figure["figure"], "kind": cand["kind"],
            "legacy_mm": [leg["width_mm"], leg["height_mm"]], "standard_mm": [cand["width_mm"], cand["height_mm"]],
            "standard_within_max": cand["width_mm"] <= profile["max_width_mm"] + 0.05 and cand["height_mm"] <= profile["max_height_mm"] + 0.05,
            "legacy_within_max": leg["width_mm"] <= profile["max_width_mm"] + 0.05 and leg["height_mm"] <= profile["max_height_mm"] + 0.05,
            "standard_tiff_dpi": cand["actual"]["tiff"]["dpi"], "standard_png_dpi": cand["actual"]["png"]["dpi"],
            "legacy_tiff_dpi": leg["actual"]["tiff"]["dpi"], "dpi_ok": round(cand["actual"]["tiff"]["dpi"]) == profile["raster_dpi"] and cand["actual"]["png"]["dpi"] == profile["preview_dpi"],
            "edge_ink_px": figure["clipping"]["candidate_edge_ink_px"], "legend": (figure.get("legend") or {}).get("candidate"),
            "pdf_fonts": cand["actual"]["pdf"]["embedded_fonts"], "source_data_identical": figure["source_data"]["identical"],
            "terms": figure.get("terms"),
        }
        if pdf_std is not None and pdf_std.exists():
            sizes = text_sizes(pdf_std)
            entry["text_standard"] = sizes
            entry["text_ok"] = bool(sizes.get("runs")) and sizes["implied_min_pt"] >= profile["min_pt"] - 0.35
        if pdf_leg is not None and pdf_leg.exists():
            entry["text_legacy"] = text_sizes(pdf_leg)
        rows.append(entry)
    return {"profile": profile, "figures": rows, "unmatched": len(report.get("only_in_legacy", [])) + len(report.get("only_in_candidate", []))}


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--profile", nargs=3, action="append", metavar=("NAME", "LEGACY_RUN", "STANDARD_RUN"), required=True)
    args = parser.parse_args(argv)
    summary = {name: audit_profile(name, Path(leg), Path(std), args.out) for name, leg, std in args.profile}
    (args.out / "profile_audit.json").write_text(json.dumps(summary, indent=1, sort_keys=True) + "\n", encoding="utf-8")
    failed = False
    for name, data in summary.items():
        bad = [f for f in data["figures"] if not (f["standard_within_max"] and f["dpi_ok"] and f["edge_ink_px"] == 0 and f.get("text_ok", True) and (f["legend"] or {}).get("complete", True) and f["source_data_identical"])]
        failed |= bool(bad) or bool(data["unmatched"])
        print(f"{name}: {len(data['figures'])} figures, unmatched {data['unmatched']}, issues {len(bad)}")
        for f in bad:
            print("   ", f["figure"], {k: f.get(k) for k in ("standard_mm", "standard_within_max", "dpi_ok", "edge_ink_px", "text_ok", "legend", "source_data_identical")})
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
