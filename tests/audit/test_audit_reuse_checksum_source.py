"""Audit regression: L1 -> L2 upstream reuse vs. the FASTQ checksum cache.

``_attach_fastq_checksums`` records *how* each SHA256 was obtained
(``sha256_source``: ``computed`` on the first run, ``cache`` afterwards) inside
the frozen ``input`` section.  ``reuse_upstream_if_compatible`` compares that
whole section, so an earlier run that computed its checksums and a later run
that reused them from the cache are judged incompatible although every FASTQ
identity (path, size, mtime, SHA256) is identical.

Marked ``xfail(strict=True)``: it documents the defect today and will start
failing (forcing removal of the marker) once production code is fixed.
"""

from __future__ import annotations

import pytest
import yaml

from rnaseq.planner import generate_plan
from rnaseq.service import _upstream_identity
from test_fastq_identity import _start_run, linux_host  # noqa: F401  (fixture re-export)
from test_runtime_dispatch import _salmon_report, mocked_downstream_runtime  # noqa: F401


pytestmark = pytest.mark.usefixtures("production_capable_execution_capacity")


def _input(run_dir):
    return yaml.safe_load((run_dir / "frozen" / "input_manifest.yaml").read_text(encoding="utf-8"))["input"]


def test_fastq_identity_is_identical_between_first_and_cached_run(monkeypatch, tmp_path, linux_host):
    report = _salmon_report(tmp_path / "project")
    generate_plan(report)
    first = _start_run(monkeypatch, report.project_dir, "CASE-L1")
    second = _start_run(monkeypatch, report.project_dir, "CASE-L2")
    strip = lambda files: [{k: v for k, v in item.items() if k != "sha256_source"} for item in files]  # noqa: E731
    # The scientific FASTQ identity is the same in both runs ...
    assert strip(_input(first)["files"]) == strip(_input(second)["files"])
    # ... only the bookkeeping field differs.
    assert {item["sha256_source"] for item in _input(first)["files"]} == {"computed"}
    assert {item["sha256_source"] for item in _input(second)["files"]} == {"cache"}


@pytest.mark.xfail(strict=True, reason="AUDIT D-01: sha256_source is part of the upstream reuse identity")
def test_upstream_reuse_identity_ignores_checksum_provenance(monkeypatch, tmp_path, linux_host):
    report = _salmon_report(tmp_path / "project")
    generate_plan(report)
    first = _start_run(monkeypatch, report.project_dir, "CASE-L1")
    second = _start_run(monkeypatch, report.project_dir, "CASE-L2")
    assert _upstream_identity(first / "frozen" / "input_manifest.yaml") == _upstream_identity(
        second / "frozen" / "input_manifest.yaml"
    )
