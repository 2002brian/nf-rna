#!/usr/bin/env python3
"""Packaging audit of a built wheel and sdist (stdlib only; review aid).

    package_audit.py --repo REPO --wheel A.whl --sdist X.tar.gz [--wheel-from-sdist B.whl]

Checks that: the bundled Arimo files and OFL.txt in the wheel match the SHA-256 values pinned in FONT_MANIFEST.json and the
working tree; the OFL text is the SIL Open Font License 1.1; every packaged file under rnaseq/ is byte-identical to the tree;
no packaged file mentions the design-reference toolkit or a user home path; the R scripts and fonts are shipped; the sdist
carries the review tools, tests and skill; and the wheel built from the sdist is byte-identical to the direct wheel.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
import tarfile
import zipfile
from pathlib import Path

FORBIDDEN = [b"bioinfo_" + b"toolkit", b"bioinfo-" + b"toolkit", b"/home/" + b"brian"]  # split so this file does not match itself


def sha(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--repo", type=Path, required=True)
    parser.add_argument("--wheel", type=Path, required=True)
    parser.add_argument("--sdist", type=Path, required=True)
    parser.add_argument("--wheel-from-sdist", type=Path)
    args = parser.parse_args(argv)
    results: dict[str, object] = {}
    wheel = zipfile.ZipFile(args.wheel)
    names = wheel.namelist()
    manifest = json.loads((args.repo / "src/rnaseq/fonts/arimo/FONT_MANIFEST.json").read_text(encoding="utf-8"))
    pinned = {name: info["sha256"] for name, info in manifest["files"].items()}
    font_ok = {}
    for name, digest in pinned.items():
        member = f"rnaseq/fonts/arimo/{name}"
        font_ok[name] = member in names and sha(wheel.read(member)) == digest and sha((args.repo / "src" / member).read_bytes()) == digest
    results["pinned_files_match_manifest_wheel_and_tree"] = font_ok
    licence = wheel.read("rnaseq/fonts/arimo/OFL.txt").decode("utf-8", "replace") if "rnaseq/fonts/arimo/OFL.txt" in names else ""
    results["ofl_text"] = {"present": bool(licence), "sil_open_font_license_1_1": "SIL OPEN FONT LICENSE Version 1.1" in licence, "copyright_notice": "Copyright" in licence}
    package_files = [n for n in names if n.startswith("rnaseq/") and not n.endswith("/")]
    def source_of(name: str) -> Path:
        # pyproject force-includes the top-level workflow/ directory as rnaseq/workflows/
        mapped = args.repo / "workflow" / name.removeprefix("rnaseq/workflows/")
        return mapped if name.startswith("rnaseq/workflows/") and mapped.exists() else args.repo / "src" / name

    differ = [n for n in package_files if source_of(n).exists() and wheel.read(n) != source_of(n).read_bytes()]
    missing = [n for n in package_files if not source_of(n).exists()]
    results["wheel_files"] = {"count": len(package_files), "differing_from_tree": differ, "not_in_tree": missing}
    results["forbidden_strings"] = {n: [f.decode() for f in FORBIDDEN if f in wheel.read(n)] for n in names if any(f in wheel.read(n) for f in FORBIDDEN)}
    results["r_scripts"] = sorted(Path(n).name for n in names if n.endswith(".R"))
    results["font_files_shipped"] = sorted(n for n in names if "/fonts/" in n)
    results["tools_in_wheel"] = [n for n in names if n.startswith(("tools/", ".claude/"))]
    sdist = tarfile.open(args.sdist)
    sd_names = sdist.getnames()
    results["sdist"] = {"has_tools": any("/tools/figure_review/" in n for n in sd_names), "has_skill": any(n.endswith(".claude/skills/nf-rna-figures/SKILL.md") for n in sd_names),
                        "has_tests": any("/tests/test_figure_" in n for n in sd_names), "has_figure_standard_doc": any(n.endswith("docs/figure-standard.md") for n in sd_names),
                        "has_font_licence": any(n.endswith("fonts/arimo/OFL.txt") for n in sd_names)}
    sd_forbidden = {}
    for member in sdist.getmembers():
        if member.isfile():
            data = sdist.extractfile(member).read()
            hits = [f.decode() for f in FORBIDDEN if f in data]
            if hits and not member.name.endswith("tests/test_figure_review_tools.py"):  # that test lists the forbidden strings it checks for
                sd_forbidden[member.name] = hits
    results["sdist_forbidden_strings"] = sd_forbidden
    if args.wheel_from_sdist:
        results["wheel_from_sdist_byte_identical"] = sha(args.wheel.read_bytes()) == sha(args.wheel_from_sdist.read_bytes())
    ok = (not sd_forbidden and all(font_ok.values()) and results["ofl_text"]["sil_open_font_license_1_1"] and not differ and not missing and not results["forbidden_strings"]
          and results["wheel_from_sdist_byte_identical" if args.wheel_from_sdist else "r_scripts"] and not results["tools_in_wheel"])
    print(json.dumps(results, indent=1, sort_keys=True))
    print("PACKAGE AUDIT:", "PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
