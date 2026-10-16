#!/usr/bin/env python3
"""Legacy-versus-standard review of the enrichment dotplots (GO, GSEA, KEGG ORA, KEGG GSEA).

Every scenario runs the real enrichment backends once with the legacy selection and once with the
standard preset on the synthetic enrichment fixture (``generate_runs.py``).  KEGG runs against the
OFFLINE synthetic fixture (``kegg_offline_profile.R``): the statistics are synthetic and exist only to
exercise the plotting and backend paths.  ``compare_runs.py`` then reports, per figure, the legacy and
standard dimensions, embedded PDF fonts, source-data identity, term-order identity, legend
completeness and a clipping check, and fails unless only the dotplots changed.

    review_enrichment.py --out DIR [--only SCENARIO ...]
"""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
KINDS = ["ORA dotplot", "GSEA dotplot"]
SCENARIOS = {
    "full_many_terms": ["--kegg-scenario", "many"],
    "kegg_few_terms": ["--kegg-scenario", "few", "--only-kegg"],
    "kegg_long_names": ["--kegg-scenario", "long", "--only-kegg"],
    "kegg_short_names": ["--kegg-scenario", "short", "--only-kegg"],
    "kegg_nonsignificant": ["--kegg-scenario", "nonsig", "--only-kegg"],
    "empty_results": ["--kegg-scenario", "empty", "--empty-foreground"],
}


def run(*args: str) -> subprocess.CompletedProcess:
    return subprocess.run([sys.executable, *args], capture_output=True, text=True, check=False)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--example", type=Path, default=HERE.parents[1] / "examples" / "simple_two_group")
    parser.add_argument("--profile", default="general")
    parser.add_argument("--only", nargs="*")
    args = parser.parse_args(argv)
    summary, failed = {}, False
    for name, extra in SCENARIOS.items():
        if args.only and name not in args.only:
            continue
        base = args.out / name
        for label, figures in (("legacy", None), ("standard", f"{{preset: standard, profile: {args.profile}}}")):
            if not (base / label).exists():
                result = run(str(HERE / "generate_runs.py"), "--example", str(args.example), "--out", str(base / label), "--skip-analysis", *extra, *(["--figures", figures] if figures else []))
                if result.returncode:
                    print(name, label, "FAILED", result.stdout[-1500:], result.stderr[-1500:])
                    return 1
        gate = run(str(HERE / "compare_runs.py"), "--legacy", str(base / "legacy"), "--candidate", str(base / "standard"), "--out", str(base / "compare"), "--require-identical-except", *KINDS)
        data = json.loads((base / "compare" / "report.json").read_text(encoding="utf-8"))
        figures = {}
        for f in data["figures"]:
            figures[f["figure"]] = {
                "legacy_mm": [f["legacy"]["width_mm"], f["legacy"]["height_mm"]], "standard_mm": [f["candidate"]["width_mm"], f["candidate"]["height_mm"]],
                "pdf_fonts": f["font"]["candidate"]["embedded_in_pdf"], "source_data_identical": f["source_data"]["identical"], "terms": f.get("terms"),
                "legend": (f.get("legend") or {}).get("candidate"), "edge_ink_px": f["clipping"]["candidate_edge_ink_px"], "legacy_edge_ink_px": f["clipping"]["legacy_edge_ink_px"]}
        summary[name] = {"gate": "pass" if gate.returncode == 0 else "FAIL", "unmatched": data["summary"]["unmatched"], "figures": figures}
        failed |= gate.returncode != 0 or data["summary"]["unmatched"] > 0
        print(name, summary[name]["gate"], gate.stdout.strip().splitlines()[-1] if gate.stdout.strip() else gate.stderr[-300:])
        for line in gate.stdout.splitlines():
            if line.startswith("GATE FAILURE"):
                print("  ", line)
    (args.out / "enrichment_summary.json").write_text(json.dumps(summary, indent=1, sort_keys=True) + "\n", encoding="utf-8")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
