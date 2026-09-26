"""Raw FASTQs are identified without re-reading them before Nextflow starts.

Regression for 4T1 run grcm39-v130-salmon-r4: before launching Nextflow,
``rnaseq run`` hashed every FASTQ twice (plan-freshness check and frozen input
manifest) and copied it into the run, reading ~145 GB for ~37 GB of input over
the WSL 9p mount.  Planning and validation now identify FASTQs by path, size
and mtime; a run computes each SHA256 at most once per file content and never
copies the raw data.
"""

from __future__ import annotations

import collections
import csv
import json
import os
from pathlib import Path

import pytest
import yaml

from rnaseq.errors import ExecutionPreflightError
from rnaseq.execution import _reference_runtime_check
from rnaseq.planner import generate_plan
from rnaseq.service import execute_retry_service_run, execute_service_run
from rnaseq.validators import validate_project
from test_runtime_dispatch import (  # noqa: F401  (fixture re-export)
    BACKEND_CONDA, LINUX, _fake_nextflow, _guard_backend, _host, _salmon_report, mocked_downstream_runtime,
)


pytestmark = pytest.mark.usefixtures("production_capable_execution_capacity")


class _StopBeforeNextflow(Exception):
    pass


def _fastqs(root: Path) -> list[Path]:
    return sorted((root / "input" / "fastq").glob("*.fastq.gz"))


@pytest.fixture
def fastq_reads(monkeypatch):
    """Count every binary open of a FASTQ, whichever code path does it."""

    reads: collections.Counter[str] = collections.Counter()
    real_open = Path.open

    def counting_open(self, mode="r", *args, **kwargs):
        if self.name.endswith(".fastq.gz") and "r" in mode and "b" in mode:
            reads[self.name] += 1
        return real_open(self, mode, *args, **kwargs)

    monkeypatch.setattr(Path, "open", counting_open)
    return reads


@pytest.fixture
def linux_host(monkeypatch, tmp_path, mocked_downstream_runtime):
    monkeypatch.setenv("RNASEQ_EXECUTION_ROOT", str(tmp_path / "execution"))
    _host(monkeypatch, *LINUX)
    _guard_backend(monkeypatch, BACKEND_CONDA)


def _start_run(monkeypatch, root: Path, case_id: str):
    """Run everything before the first Nextflow launch, then stop."""

    def first_nextflow(*_args, **_kwargs):
        raise _StopBeforeNextflow

    monkeypatch.setattr("rnaseq.service._run_command", first_nextflow)
    with pytest.raises(_StopBeforeNextflow):
        execute_service_run(validate_project(root), case_id=case_id)
    (run_dir,) = (root / "runs" / case_id).iterdir()
    return run_dir


# ---------------------------------------------------------------- A: routine checks read no FASTQ content


@pytest.mark.skipif(os.geteuid() == 0, reason="root ignores file permissions")
def test_validate_plan_and_doctor_never_read_fastq_contents(tmp_path):
    report = _salmon_report(tmp_path / "project")
    root = report.project_dir
    for path in _fastqs(root):
        path.chmod(0)  # stat still works; any read (open, sendfile, copy) would fail
    try:
        assert validate_project(root).is_valid
        generate_plan(validate_project(root))
        assert _reference_runtime_check(root).verdict != "FAIL"
    finally:
        for path in _fastqs(root):
            path.chmod(0o644)
    manifest = yaml.safe_load((root / "planning" / "manifest.preview.yaml").read_text(encoding="utf-8"))
    for entry in manifest["input"]["files"]:
        assert set(entry) == {"relative_path", "size_bytes", "mtime_ns"}


# ---------------------------------------------------------------- B: one pass per file per content, no copy


def test_fresh_run_reads_each_fastq_once_never_copies_and_reuses_checksums(monkeypatch, tmp_path, linux_host, fastq_reads):
    report = _salmon_report(tmp_path / "project")
    root = report.project_dir
    generate_plan(report)
    assert sum(fastq_reads.values()) == 0

    first = _start_run(monkeypatch, root, "CASE-FIRST")
    assert dict(fastq_reads) == {path.name: 1 for path in _fastqs(root)}
    assert not (first / "frozen" / "input" / "fastq").exists()

    fastq_reads.clear()
    second = _start_run(monkeypatch, root, "CASE-SECOND")
    assert sum(fastq_reads.values()) == 0
    first_files = yaml.safe_load((first / "frozen" / "input_manifest.yaml").read_text(encoding="utf-8"))["input"]
    second_files = yaml.safe_load((second / "frozen" / "input_manifest.yaml").read_text(encoding="utf-8"))["input"]
    assert [item["sha256"] for item in first_files["files"]] == [item["sha256"] for item in second_files["files"]]
    assert {item["sha256_source"] for item in first_files["files"]} == {"computed"}
    assert {item["sha256_source"] for item in second_files["files"]} == {"cache"}
    assert second_files["identity_contract"]["checksum_reuse_key"] == "resolved_path + size_bytes + mtime_ns"


# ---------------------------------------------------------------- C: a changed FASTQ is detected


def test_changed_fastq_requires_replanning_and_gets_a_new_checksum(monkeypatch, tmp_path, linux_host, fastq_reads):
    report = _salmon_report(tmp_path / "project")
    root = report.project_dir
    generate_plan(report)
    first = _start_run(monkeypatch, root, "CASE-BEFORE")
    changed = _fastqs(root)[0]
    status = changed.stat()
    changed.write_bytes(b"y")  # same size, different content
    os.utime(changed, ns=(status.st_atime_ns, status.st_mtime_ns + 1_000_000_000))

    with pytest.raises(ExecutionPreflightError, match="changed since planning"):
        execute_service_run(validate_project(root), case_id="CASE-STALE")
    generate_plan(validate_project(root))
    fastq_reads.clear()
    after = _start_run(monkeypatch, root, "CASE-AFTER")
    assert dict(fastq_reads) == {changed.name: 1}
    before_sha = {item["relative_path"]: item["sha256"] for item in yaml.safe_load((first / "frozen" / "input_manifest.yaml").read_text(encoding="utf-8"))["input"]["files"]}
    after_sha = {item["relative_path"]: item["sha256"] for item in yaml.safe_load((after / "frozen" / "input_manifest.yaml").read_text(encoding="utf-8"))["input"]["files"]}
    relative = f"input/fastq/{changed.name}"
    assert before_sha[relative] != after_sha[relative]
    assert {key: value for key, value in before_sha.items() if key != relative} == {key: value for key, value in after_sha.items() if key != relative}


@pytest.mark.parametrize("source_status", ["FAILED", "INTERRUPTED"])
def test_retry_refuses_a_fastq_that_no_longer_matches_its_frozen_identity(monkeypatch, tmp_path, linux_host, source_status):
    report = _salmon_report(tmp_path / "project")
    root = report.project_dir
    generate_plan(report)
    monkeypatch.setattr("rnaseq.service._run_command", lambda *_args, **_kwargs: 1)
    with pytest.raises(Exception):
        execute_service_run(validate_project(root), case_id="CASE-FAILED")
    (source,) = (root / "runs" / "CASE-FAILED").iterdir()
    state = json.loads((source / "run_state.json").read_text(encoding="utf-8"))
    assert state["status"] == "FAILED"
    (source / "run_state.json").write_text(json.dumps({**state, "status": source_status}), encoding="utf-8")
    changed = _fastqs(root)[0]
    changed.write_bytes(b"different content")
    with pytest.raises(ExecutionPreflightError, match="no longer matches its frozen identity"):
        execute_retry_service_run(root, retry_of=f"CASE-FAILED/{source.name}")


# ---------------------------------------------------------------- D / E: pairing and provenance are intact


def test_frozen_samplesheet_pairs_the_validated_fastqs_in_place_with_full_provenance(monkeypatch, tmp_path, linux_host):
    report = _salmon_report(tmp_path / "project")
    root = report.project_dir
    generate_plan(report)
    observed: list[list[str]] = []
    monkeypatch.setattr("rnaseq.service._run_command", _fake_nextflow(observed))
    run = execute_service_run(validate_project(root), case_id="CASE-PAIRS")
    with (run.run_dir / "frozen" / "samplesheet.csv").open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle))
    fastq = (root / "input" / "fastq").resolve()
    assert [(row["sample"], row["fastq_1"], row["fastq_2"]) for row in rows] == [
        (sample, str(fastq / f"{sample}_R1.fastq.gz"), str(fastq / f"{sample}_R2.fastq.gz")) for sample in ("C1", "C2", "T1", "T2")
    ]
    manifest = yaml.safe_load((run.run_dir / "frozen" / "input_manifest.yaml").read_text(encoding="utf-8"))
    assert len(manifest["input"]["files"]) == 8
    assert all(len(item["sha256"]) == 64 and item["resolved_path"].startswith(str(fastq)) for item in manifest["input"]["files"])
    provenance = yaml.safe_load((run.run_dir / "provenance" / "run_provenance.yaml").read_text(encoding="utf-8"))
    assert provenance["source_revision"] == "test-revision"
    assert provenance["input_manifest_sha256"]
    upstream = observed[0]
    assert str(run.run_dir / "frozen" / "samplesheet.csv") in " ".join(upstream)
