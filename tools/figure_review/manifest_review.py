#!/usr/bin/env python3
"""Manifest review (stdlib only; review aid).

1. ``--compat BASELINE_RUN CURRENT_RUN``: backward compatibility of ``figure_manifest.json`` records.  For every figure present in
   both runs, every key of the baseline record must still exist; values that differ and keys that were added are listed.
2. ``--audit STANDARD_RUN``: which provenance fields a manifest of that run carries (applied preset, profile, actual font, font
   SHA-256, embedded fonts, palette strategy, dimensions, dpi) at the manifest level and per figure, and whether the manifest-level
   ``colour_policy`` describes the preset that was actually applied.

    manifest_review.py --compat BASELINE CURRENT --audit STANDARD_RUN [--out report.json]
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path


def records(root: Path) -> dict[str, tuple[dict, dict]]:
    found = {}
    for manifest in sorted(root.rglob("figure_manifest.json")):
        data = json.loads(manifest.read_text(encoding="utf-8"))
        base = manifest.parent.relative_to(root).as_posix()
        for record in data["figures"]:
            found[f"{base}/{record['id']}"] = (record, data)
    return found


def compat(baseline: Path, current: Path) -> dict:
    a, b = records(baseline), records(current)
    result = {"schema_versions": sorted({d["schema_version"] for _, d in {**a, **b}.values()}), "figures_compared": 0, "removed_keys": {}, "changed_values": {}, "added_keys": {}, "only_in_baseline": sorted(set(a) - set(b)), "only_in_current": sorted(set(b) - set(a))}
    for key in sorted(set(a) & set(b)):
        old, new = a[key][0], b[key][0]
        result["figures_compared"] += 1
        if gone := sorted(set(old) - set(new)):
            result["removed_keys"][key] = gone
        if added := sorted(set(new) - set(old)):
            result["added_keys"][key] = added
        changed = sorted(k for k in set(old) & set(new) if old[k] != new[k])
        if changed:
            result["changed_values"][key] = changed
    top_old = next(iter(a.values()))[1] if a else {}
    top_new = next(iter(b.values()))[1] if b else {}
    result["top_level_removed_keys"] = sorted(set(top_old) - set(top_new))
    result["top_level_added_keys"] = sorted(set(top_new) - set(top_old))
    return result


def audit(run: Path) -> dict:
    recs = records(run)
    sample_key, (sample, top) = next(iter(recs.items()))
    style = top.get("figure_style", {})
    legacy_font_family = top.get("font", {}).get("family")
    standard_font = (style.get("font") or {})
    out = {
        "manifest_level": {
            "preset": style.get("preset"), "requested_preset": style.get("requested_preset"), "profile": (top.get("profile") or {}).get("name"),
            "profile_limits_recorded": sorted((top.get("profile") or {}).keys()),
            "font_block_describes": legacy_font_family, "standard_font_block": {k: standard_font.get(k) for k in ("family", "sha256", "embedded_fonts", "variable")},
            "standard_figures_drawn": style.get("standard_figures_drawn"),
        },
        "per_figure_fields": sorted({k for r, _ in recs.values() for k in r}),
        "per_figure_style_block": sorted({k for r, _ in recs.values() for k in (r.get("figure_style") or {})}),
        "figures_missing_style_block": sorted(k for k, (r, _) in recs.items() if not r.get("figure_style")),
        "figures_with_palette_name": sorted(k for k, (r, _) in recs.items() if "palette" in json.dumps(r.get("encoding", "")) or "palette" in r),
        "figures_with_dimensions_and_dpi": all(r.get("width_mm") and r.get("height_mm") and all("dpi" in f or f["format"] == "pdf" for f in r["files"]) for r, _ in recs.values()),
    }
    policy = top.get("colour_policy", {})
    applied = style.get("preset")
    direction = (policy.get("direction") or {}).get("Up")
    out["colour_policy"] = {"manifest_level_direction_up": direction, "describes_applied_preset": None}
    standard_up = "#B2182B"
    if applied == "standard" and direction:
        out["colour_policy"]["describes_applied_preset"] = direction.upper() == standard_up
    volcano = [r for k, (r, _) in recs.items() if k.endswith("/volcano")]
    if volcano:
        out["colour_policy"]["volcano_record_colours"] = volcano[0].get("colours")
    out["font_block_is_legacy_font_under_standard"] = bool(applied == "standard" and legacy_font_family and legacy_font_family != standard_font.get("family"))
    palette_names = {}
    for key, (r, _) in recs.items():
        enc = r.get("encoding") or {}
        if isinstance(enc, dict) and enc.get("palette"):
            palette_names[key] = enc["palette"]
    out["palette_strategy_names"] = sorted(set(palette_names.values()))
    return out


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--compat", nargs=2, metavar=("BASELINE", "CURRENT"), type=Path)
    parser.add_argument("--audit", type=Path)
    parser.add_argument("--out", type=Path)
    args = parser.parse_args(argv)
    report = {}
    if args.compat:
        report["compat"] = compat(*args.compat)
    if args.audit:
        report["audit"] = audit(args.audit)
    text = json.dumps(report, indent=1, sort_keys=True)
    if args.out:
        args.out.write_text(text + "\n", encoding="utf-8")
    print(text)
    c = report.get("compat")
    return 1 if c and (c["removed_keys"] or c["top_level_removed_keys"] or c["only_in_baseline"]) else 0


if __name__ == "__main__":
    sys.exit(main())
