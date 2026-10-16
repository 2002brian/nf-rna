#!/usr/bin/env python3
"""Repeat one offline KEGG render many times and report how many distinct PDF layouts come out (review aid).

The standard KEGG GSEA dotplot was observed to be drawn, in about 4% of runs under CPU load, with a different PDF layout
(the axis and legend rows 15 pt taller each, the panel shorter) while its PNG and TIFF were identical.  This tool reproduces
that measurement: it runs ``generate_runs.py --only-kegg`` N times with a private HOME each, optionally under CPU load, and
groups the resulting PDFs by content (``compare_runs.pdf_content_hash`` ignores dates).

    pdf_layout_repeat.py --out DIR [--runs 100] [--parallel 24] [--busy 16] [--preset standard] [--profile general]
"""

from __future__ import annotations

import argparse
import collections
import concurrent.futures
import importlib.util
import os
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location("compare_runs", HERE / "compare_runs.py")
compare_runs = importlib.util.module_from_spec(spec)
spec.loader.exec_module(compare_runs)


def group_by_content(pdfs: list[Path]) -> dict[str, list[Path]]:
    groups: dict[str, list[Path]] = collections.defaultdict(list)
    for pdf in pdfs:
        groups[compare_runs.pdf_content_hash(pdf)].append(pdf)
    return dict(groups)


def one(job: tuple[int, Path, str, str, Path]) -> Path | None:
    index, out, preset, profile, example = job
    home = out / f"home{index}"
    home.mkdir(parents=True)
    env = {**os.environ, "HOME": str(home)}
    subprocess.run([sys.executable, str(HERE / "generate_runs.py"), "--example", str(example), "--out", str(out / f"run{index}"), "--skip-analysis", "--only-kegg",
                    "--kegg-scenario", "many", "--figures", f"{{preset: {preset}, profile: {profile}}}"], env=env, capture_output=True, check=False)
    pdf = out / f"run{index}" / "enrichment" / "gsea_kegg" / "Treatment_vs_Control" / "dotplot.pdf"
    return pdf if pdf.exists() else None


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--runs", type=int, default=100)
    parser.add_argument("--parallel", type=int, default=24)
    parser.add_argument("--busy", type=int, default=0, help="number of busy-loop processes to run during the batch")
    parser.add_argument("--preset", default="standard")
    parser.add_argument("--profile", default="general")
    parser.add_argument("--example", type=Path, default=HERE.parents[1] / "examples" / "simple_two_group")
    args = parser.parse_args(argv)
    args.out.mkdir(parents=True, exist_ok=False)
    busy = [subprocess.Popen([sys.executable, "-c", "while True: pass"]) for _ in range(args.busy)]
    try:
        with concurrent.futures.ThreadPoolExecutor(args.parallel) as pool:
            pdfs = [p for p in pool.map(one, [(i, args.out, args.preset, args.profile, args.example) for i in range(args.runs)]) if p]
    finally:
        for process in busy:
            process.kill()
    groups = group_by_content(pdfs)
    print(f"{len(pdfs)} of {args.runs} runs produced a PDF; distinct PDF layouts: {len(groups)}")
    for digest, files in sorted(groups.items(), key=lambda item: -len(item[1])):
        print(f"  {digest[:12]}  {len(files)} runs  e.g. {files[0]}")
    return 0 if len(groups) <= 1 else 1


if __name__ == "__main__":
    sys.exit(main())
