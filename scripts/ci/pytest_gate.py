"""Release gate over a pytest JUnit report: count outcomes and refuse unapproved skips.

A scientific CI run must not turn missing coverage into a green result.  Every
skipped test needs a matching entry in the approved-skip policy (a skip that is
legitimately platform specific, with a written justification); any other skip,
any failure or error, too few collected tests, or a scientific test module that
did not execute in full fails the gate.

Usage::

    python scripts/ci/pytest_gate.py REPORT.xml --policy ci/approved_skips.json \
        [--min-collected N] [--require-module tests/science/test_x.py ...] [--summary PATH]

Only the standard library is used, so the gate itself has no dependencies.
"""

from __future__ import annotations

import argparse
import json
import re
import sys
import xml.etree.ElementTree as ET
from collections import Counter
from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class Case:
    nodeid: str
    module: str
    outcome: str  # passed | failed | error | skipped
    message: str


def _module(classname: str, file: str | None) -> str:
    if file:
        return file.replace("\\", "/")
    # pytest writes classname as dotted module path (plus class); convert to a file path.
    parts = classname.split(".")
    while parts and parts[-1][:1].isupper():
        parts.pop()
    return "/".join(parts) + ".py"


def read_report(path: Path) -> list[Case]:
    root = ET.parse(path).getroot()
    cases: list[Case] = []
    for element in root.iter("testcase"):
        classname = element.get("classname", "")
        module = _module(classname, element.get("file"))
        nodeid = f"{module}::{element.get('name', '')}"
        outcome, message = "passed", ""
        for child in element:
            if child.tag in {"failure", "error", "skipped"}:
                outcome = {"failure": "failed", "error": "error", "skipped": "skipped"}[child.tag]
                message = (child.get("message") or child.text or "").strip()
                if child.tag == "skipped" and child.get("type") == "pytest.xfail":
                    outcome = "xfailed"
                break
        cases.append(Case(nodeid, module, outcome, message))
    return cases


def load_policy(path: Path) -> list[dict[str, str]]:
    entries = json.loads(path.read_text(encoding="utf-8"))["approved_skips"]
    for entry in entries:
        missing = {"reason_regex", "category", "justification"} - set(entry)
        if missing:
            raise SystemExit(f"approved skip entry lacks {sorted(missing)}: {entry}")
        re.compile(entry["reason_regex"])
    return entries


def approved(case: Case, policy: list[dict[str, str]]) -> dict[str, str] | None:
    for entry in policy:
        if re.search(entry["reason_regex"], case.message) and re.search(entry.get("module_regex", ".*"), case.module):
            return entry
    return None


def evaluate(cases: list[Case], policy: list[dict[str, str]], *, min_collected: int, required_modules: list[str]) -> tuple[list[str], dict[str, object]]:
    problems: list[str] = []
    counts = Counter(case.outcome for case in cases)
    approved_skips: list[tuple[Case, dict[str, str]]] = []
    for case in cases:
        if case.outcome in {"failed", "error"}:
            problems.append(f"{case.outcome.upper()}: {case.nodeid}: {case.message[:300]}")
        elif case.outcome == "skipped":
            entry = approved(case, policy)
            if entry is None:
                problems.append(f"UNAPPROVED SKIP (missing coverage): {case.nodeid}: {case.message}")
            else:
                approved_skips.append((case, entry))
    if len(cases) < min_collected:
        problems.append(f"only {len(cases)} tests collected; at least {min_collected} are required")
    for module in required_modules:
        mine = [case for case in cases if case.module.endswith(module)]
        if not mine:
            problems.append(f"required scientific module did not run: {module}")
        elif any(case.outcome != "passed" for case in mine):
            problems.append(f"required scientific module did not pass in full: {module} ({Counter(c.outcome for c in mine)})")
    summary = {
        "collected": len(cases), "passed": counts.get("passed", 0), "failed": counts.get("failed", 0),
        "errors": counts.get("error", 0), "skipped": counts.get("skipped", 0), "xfailed": counts.get("xfailed", 0),
        "approved_skips": [{"test": case.nodeid, "reason": case.message, "category": entry["category"]} for case, entry in approved_skips],
        "required_modules": {module: Counter(c.outcome for c in cases if c.module.endswith(module)) for module in required_modules},
        "problems": problems,
    }
    return problems, summary


def markdown(summary: dict[str, object]) -> str:
    lines = ["## Scientific test gate", "", "| collected | passed | failed | errors | skipped (approved) | xfailed |", "|---|---|---|---|---|---|"]
    lines.append(f"| {summary['collected']} | {summary['passed']} | {summary['failed']} | {summary['errors']} | {summary['skipped']} | {summary['xfailed']} |")
    lines += ["", "### Required scientific modules", ""]
    for module, outcomes in summary["required_modules"].items():  # type: ignore[union-attr]
        lines.append(f"- `{module}`: {dict(outcomes) or 'NOT RUN'}")
    skips = summary["approved_skips"]
    lines += ["", f"### Approved platform skips ({len(skips)})", ""]  # type: ignore[arg-type]
    lines += [f"- `{item['test']}` [{item['category']}]: {item['reason']}" for item in skips]  # type: ignore[union-attr]
    problems = summary["problems"]
    lines += ["", "### Result: " + ("FAIL" if problems else "PASS"), ""]
    lines += [f"- {problem}" for problem in problems]  # type: ignore[union-attr]
    return "\n".join(lines) + "\n"


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("report", type=Path)
    parser.add_argument("--policy", type=Path, required=True)
    parser.add_argument("--min-collected", type=int, default=1)
    parser.add_argument("--require-module", action="append", default=[])
    parser.add_argument("--summary", type=Path, help="also append a Markdown summary here (e.g. $GITHUB_STEP_SUMMARY)")
    parser.add_argument("--json", type=Path, help="write the machine-readable summary here")
    args = parser.parse_args(argv)
    problems, summary = evaluate(read_report(args.report), load_policy(args.policy), min_collected=args.min_collected, required_modules=args.require_module)
    text = markdown(summary)
    print(text)
    if args.summary:
        with args.summary.open("a", encoding="utf-8") as handle:
            handle.write(text)
    if args.json:
        args.json.write_text(json.dumps(summary, indent=2, default=dict) + "\n", encoding="utf-8")
    return 1 if problems else 0


if __name__ == "__main__":
    sys.exit(main())
