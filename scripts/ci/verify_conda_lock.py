"""Verify that a Conda prefix contains exactly the packages of an nf-rna lock file.

``conda env create --file LOCK`` (and micromamba) can exit successfully while a
Bioconductor *data* package's post-link download failed: the conda record is
installed but the R library directory is missing.  This check therefore compares
the prefix's ``conda-meta`` records with the lock (name, version and build, both
directions) and, for R prefixes, requires every locked Bioconductor/CRAN package
to be loadable as an installed R library directory.

Usage: python scripts/ci/verify_conda_lock.py LOCK.yml PREFIX [--r-library]
"""

from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path

SPEC = re.compile(r"^(?:(?P<channel>[^:]+)::)?(?P<name>[^=\s]+)==?(?P<version>[^=\s]+)=(?P<build>\S+)$")


def locked(lock: Path) -> dict[str, tuple[str, str]]:
    packages: dict[str, tuple[str, str]] = {}
    in_dependencies = False
    for raw in lock.read_text(encoding="utf-8").splitlines():
        line = raw.strip()
        if line.startswith("dependencies:"):
            in_dependencies = True
            continue
        if not in_dependencies or not line.startswith("- "):
            continue
        match = SPEC.match(line[2:].strip())
        if not match:
            raise SystemExit(f"unparseable lock entry: {line}")
        packages[match["name"]] = (match["version"], match["build"])
    if not packages:
        raise SystemExit(f"no dependencies found in {lock}")
    return packages


def installed(prefix: Path) -> dict[str, tuple[str, str]]:
    records = {}
    for path in (prefix / "conda-meta").glob("*.json"):
        record = json.loads(path.read_text(encoding="utf-8"))
        records[record["name"]] = (record["version"], record["build"])
    return records


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("lock", type=Path)
    parser.add_argument("prefix", type=Path)
    parser.add_argument("--r-library", action="store_true", help="also require an R library directory for every r-*/bioconductor-* record")
    args = parser.parse_args(argv)
    expected, actual = locked(args.lock), installed(args.prefix)
    problems = [f"missing from prefix: {name} {v}={b}" for name, (v, b) in sorted(expected.items()) if name not in actual]
    problems += [f"version/build differs: {name} locked {expected[name]} installed {actual[name]}" for name in sorted(expected) if name in actual and actual[name] != expected[name]]
    problems += [f"not in lock: {name} {actual[name]}" for name in sorted(set(actual) - set(expected))]
    if args.r_library:
        library = args.prefix / "lib" / "R" / "library"
        present = {path.name.lower() for path in library.iterdir()} if library.is_dir() else set()
        for name in sorted(expected):
            if (name.startswith("bioconductor-") or name.startswith("r-")) and name not in {"r-base", "bioconductor-data-packages"}:
                package = name.split("-", 1)[1]
                if package not in present:
                    problems.append(f"R library directory absent (failed post-link download?): {name}")
    print(f"lock packages: {len(expected)}; installed records: {len(actual)}; problems: {len(problems)}")
    for problem in problems:
        print("  " + problem)
    return 1 if problems else 0


if __name__ == "__main__":
    sys.exit(main())
