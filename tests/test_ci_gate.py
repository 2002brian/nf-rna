"""The scientific CI gate itself must refuse missing coverage (scripts/ci)."""

from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]


def _load(name: str):
    spec = importlib.util.spec_from_file_location(f"nf_rna_ci_{name}", ROOT / "scripts" / "ci" / f"{name}.py")
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    sys.modules[f"nf_rna_ci_{name}"] = module  # dataclasses resolve their module through sys.modules
    spec.loader.exec_module(module)
    return module


gate = _load("pytest_gate")
lock = _load("verify_conda_lock")


def _report(tmp_path: Path, cases: list[tuple[str, str, str | None, str]]) -> Path:
    body = []
    for classname, name, outcome, message in cases:
        child = "" if outcome is None else f'<{outcome} message="{message}" type="pytest.skip"/>'
        body.append(f'<testcase classname="{classname}" name="{name}">{child}</testcase>')
    path = tmp_path / "report.xml"
    path.write_text(f'<testsuites><testsuite name="pytest">{"".join(body)}</testsuite></testsuites>', encoding="utf-8")
    return path


POLICY = [{"reason_regex": "^real SJ project not present$", "module_regex": "tests/test_run_status\\.py$", "category": "local-real-data", "justification": "x"}]


def test_gate_passes_only_with_approved_skips_and_full_required_modules(tmp_path):
    cases = gate.read_report(_report(tmp_path, [
        ("tests.science.test_deseq2_designs", "test_a", None, ""),
        ("tests.test_run_status", "test_real", "skipped", "real SJ project not present"),
    ]))
    problems, summary = gate.evaluate(cases, POLICY, min_collected=2, required_modules=["tests/science/test_deseq2_designs.py"])
    assert problems == [] and summary["passed"] == 1 and summary["skipped"] == 1


@pytest.mark.parametrize(("case", "expected"), [
    (("tests.test_l2", "test_r", "skipped", "R packages unavailable: DESeq2"), "UNAPPROVED SKIP"),
    (("tests.test_l2", "test_r", "failure", "boom"), "FAILED"),
    (("tests.test_l2", "test_r", "error", "boom"), "ERROR"),
    # An approved reason in the wrong module is still missing coverage.
    (("tests.test_l2", "test_r", "skipped", "real SJ project not present"), "UNAPPROVED SKIP"),
])
def test_gate_refuses_failures_errors_and_unapproved_skips(tmp_path, case, expected):
    cases = gate.read_report(_report(tmp_path, [case]))
    problems, _ = gate.evaluate(cases, POLICY, min_collected=1, required_modules=[])
    assert len(problems) == 1 and problems[0].startswith(expected)


def test_gate_refuses_a_required_module_that_did_not_run_or_skipped(tmp_path):
    cases = gate.read_report(_report(tmp_path, [("tests.science.test_x", "test_a", "skipped", "real SJ project not present")]))
    problems, _ = gate.evaluate(cases, POLICY, min_collected=1, required_modules=["tests/science/test_x.py", "tests/science/test_y.py"])
    assert any("did not pass in full: tests/science/test_x.py" in p for p in problems)
    assert any("did not run: tests/science/test_y.py" in p for p in problems)


def test_gate_refuses_too_few_collected_tests(tmp_path):
    cases = gate.read_report(_report(tmp_path, [("tests.test_a", "test_a", None, "")]))
    problems, _ = gate.evaluate(cases, POLICY, min_collected=10, required_modules=[])
    assert problems == ["only 1 tests collected; at least 10 are required"]


def test_repository_skip_policy_is_valid_and_justified():
    entries = gate.load_policy(ROOT / "ci" / "approved_skips.json")
    assert entries and all(entry["justification"].strip() for entry in entries)
    # No R, Bioconductor, tool or Nextflow skip may ever be approved: those are missing scientific coverage.
    for entry in entries:
        assert not any(word in entry["reason_regex"].lower() for word in ("rscript", "r package", "runtime", "nextflow", "featurecounts", "samtools", "poppler"))


def test_lock_verifier_parses_both_lock_formats_and_finds_drift(tmp_path):
    downstream = lock.locked(ROOT / "workflow" / "envs" / "locks" / "nf-rna-downstream-linux-64.lock.yml")
    hisat2 = lock.locked(ROOT / "workflow" / "envs" / "hisat2-featurecounts-linux-64.yml")
    assert downstream["bioconductor-deseq2"] == ("1.50.2", "r45ha27e39d_0")
    assert hisat2["subread"] == ("2.0.6", "he4a0461_2")
    prefix = tmp_path / "prefix"
    (prefix / "conda-meta").mkdir(parents=True)
    for name, (version, build) in hisat2.items():
        version = "9.9" if name == "subread" else version
        (prefix / "conda-meta" / f"{name}-{version}-{build}.json").write_text(json.dumps({"name": name, "version": version, "build": build}), encoding="utf-8")
    assert lock.main([str(ROOT / "workflow" / "envs" / "hisat2-featurecounts-linux-64.yml"), str(prefix)]) == 1


def test_every_file_the_scientific_gate_needs_is_visible_to_git():
    """A clean checkout must contain the oracles and fixtures (.gitignore excludes **/reference(s)/ and quant.sf)."""

    import shutil
    import subprocess

    if shutil.which("git") is None or not (ROOT / ".git").exists():
        pytest.fail("the scientific gate runs from a Git checkout; git and .git are required")
    needed = [path.relative_to(ROOT).as_posix() for path in (ROOT / "tests" / "science").rglob("*")
              if path.is_file() and "__pycache__" not in path.parts]
    needed += ["scripts/ci/pytest_gate.py", "scripts/ci/verify_conda_lock.py", "ci/approved_skips.json", "ci/requirements-ci.txt"]
    assert any(name.startswith("tests/science/oracles/") for name in needed)
    ignored = subprocess.run(["git", "check-ignore", "--no-index", *needed], cwd=ROOT, capture_output=True, text=True).stdout.split()
    assert ignored == [], f"ignored by .gitignore, so absent from a clean checkout: {ignored}"
