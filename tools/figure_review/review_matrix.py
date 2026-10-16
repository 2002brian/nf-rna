#!/usr/bin/env python3
"""Legacy-versus-standard figure review over a matrix of representative datasets.

For every scenario (sample count, grouping, label length) a synthetic project is written
(``make_demo_project.py``), L1 + L2 are run once with the legacy selection and once with the
standard preset (``generate_runs.py``), and ``compare_runs.py`` writes the side-by-side report.
The gate allows only the migrated figure kinds to change.

    review_matrix.py --out DIR [--profile general]
"""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
KINDS = ["volcano", "PCA", "sample expression distribution boxplot", "sample correlation heatmap", "DEG heatmap"]
SCENARIOS = {
    "n06": ["--samples-per-group", "3"],
    "n12": ["--samples-per-group", "6"],
    "n26_groups4": ["--samples-per-group", "13", "--dose-levels", "4"],
    "n12_long_labels": ["--samples-per-group", "6", "--sample-prefix", "Organoid_Patient_Derived_Line_"],
}


def run(*args: str) -> subprocess.CompletedProcess:
    return subprocess.run([sys.executable, *args], capture_output=True, text=True, check=False)


def order_identity(standard_run: Path) -> dict:
    """Standard-run display order (from its manifests) against the legacy code path recomputed from the source tables."""

    result = subprocess.run(["Rscript", str(HERE / "display_order.R"), str(standard_run)], capture_output=True, text=True, check=False)
    if result.returncode:
        return {"error": result.stderr[-500:]}
    expected = json.loads(result.stdout)
    verdicts = {}
    for manifest in sorted(standard_run.rglob("figure_manifest.json")):
        for figure in json.loads(manifest.read_text(encoding="utf-8"))["figures"]:
            if "display_order" not in figure:
                continue
            key = (manifest.parent / figure["id"]).relative_to(standard_run).as_posix()
            want = expected.get(key)
            drawn = figure["display_order"]
            verdicts[key] = bool(want) and drawn["rows"] == want["rows"] and drawn["columns"] == want["columns"]
    return verdicts


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--profile", default="general")
    parser.add_argument("--only", nargs="*", help="scenario names to run (default: all)")
    args = parser.parse_args(argv)
    summary = {}
    for name, extra in SCENARIOS.items():
        if args.only and name not in args.only:
            continue
        base = args.out / name
        project = base / "project"
        if not project.exists():
            run(str(HERE / "make_demo_project.py"), str(project), *extra)
        for label, figures in (("legacy", None), ("standard", f"{{preset: standard, profile: {args.profile}}}")):
            target = base / label
            if not target.exists():
                result = run(str(HERE / "generate_runs.py"), "--example", str(project), "--out", str(target), "--skip-enrichment", *(["--figures", figures] if figures else []))
                if result.returncode:
                    print(result.stdout[-2000:], result.stderr[-2000:])
                    return 1
        report = base / "compare"
        gate = run(str(HERE / "compare_runs.py"), "--legacy", str(base / "legacy"), "--candidate", str(base / "standard"), "--out", str(report), "--require-identical-except", *KINDS)
        data = json.loads((report / "report.json").read_text(encoding="utf-8"))
        summary[name] = {
            "gate": "pass" if gate.returncode == 0 else "FAIL",
            "order_identical_to_legacy_path": order_identity(base / "standard"),
            "figures": {f["figure"].split("/")[-1]: {"legacy_mm": [f["legacy"]["width_mm"], f["legacy"]["height_mm"]], "standard_mm": [f["candidate"]["width_mm"], f["candidate"]["height_mm"]],
                                                     "status": f["status"], "pdf_fonts": f["font"]["candidate"]["embedded_in_pdf"], "source_data_identical": f["source_data"]["identical"]} for f in data["figures"]},
        }
        print(name, summary[name]["gate"], gate.stdout.strip().splitlines()[-1] if gate.stdout.strip() else "")
    (args.out / "matrix_summary.json").write_text(json.dumps(summary, indent=1, sort_keys=True) + "\n", encoding="utf-8")
    return 0 if all(item["gate"] == "pass" and all(item["order_identical_to_legacy_path"].values()) for item in summary.values()) else 1


if __name__ == "__main__":
    sys.exit(main())
