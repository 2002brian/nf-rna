"""L1 -> ``rnaseq design`` -> L2 with ``--reuse-upstream`` on real frozen manifests (audit D-01).

The L1 run computes every FASTQ SHA256 (``sha256_source: computed``); the L2
run gets the same digests from the checksum cache (``sha256_source: cache``).
Only Nextflow is faked; validation, planning, freezing, reuse checks and
delivery are the production code paths.
"""

from __future__ import annotations

import json

import pytest
import yaml
from typer.testing import CliRunner

from rnaseq.cli import app
from rnaseq.errors import ExecutionPreflightError
from rnaseq.planner import generate_plan
from rnaseq.service import execute_service_run
from rnaseq.validators import validate_project
from test_runtime_dispatch import LINUX, _fake_nextflow, _guard_conda, _host, _salmon_report, mocked_downstream_runtime  # noqa: F401


pytestmark = pytest.mark.usefixtures("production_capable_execution_capacity")
runner = CliRunner()


def _input_files(run_dir):
    return yaml.safe_load((run_dir / "frozen" / "input_manifest.yaml").read_text(encoding="utf-8"))["input"]["files"]


@pytest.fixture
def l1_project(monkeypatch, tmp_path, mocked_downstream_runtime):
    monkeypatch.setenv("RNASEQ_EXECUTION_ROOT", str(tmp_path / "execution"))
    _host(monkeypatch, *LINUX)
    _guard_conda(monkeypatch)
    root = tmp_path / "project"
    _salmon_report(root)
    config = yaml.safe_load((root / "project.yaml").read_text(encoding="utf-8"))
    config["project"]["preset"] = "L1"
    (root / "project.yaml").write_text(yaml.safe_dump(config, sort_keys=False), encoding="utf-8")
    (root / "contrasts.csv").write_text("contrast_id,factor,numerator,denominator\n", encoding="utf-8")
    report = validate_project(root)
    assert report.is_valid, report.errors
    generate_plan(report)
    observed: list[list[str]] = []
    monkeypatch.setattr("rnaseq.service._run_command", _fake_nextflow(observed))
    l1 = execute_service_run(report, case_id="L1-001")
    assert json.loads(l1.state_path.read_text(encoding="utf-8"))["status"] == "SUCCESS"
    assert len(observed) == 2  # upstream + downstream
    result = runner.invoke(app, ["design", str(root), "--no-enrichment", "--yes"])
    assert result.exit_code == 0, result.output
    return root, l1


def test_l2_reuses_l1_upstream_although_checksum_provenance_differs(monkeypatch, l1_project):
    root, l1 = l1_project
    report = validate_project(root)
    assert report.is_valid, report.errors
    assert report.config.project.preset.value == "L2"
    generate_plan(report)
    observed: list[list[str]] = []
    monkeypatch.setattr("rnaseq.service._run_command", _fake_nextflow(observed))
    l2 = execute_service_run(report, case_id="L2-001", reuse_upstream=f"L1-001/{l1.run_id}")

    state = json.loads(l2.state_path.read_text(encoding="utf-8"))
    assert state["status"] == "SUCCESS"
    # Only the downstream workflow ran; the upstream came from L1.
    assert len(observed) == 1 and "nf-core/rnaseq" not in observed[0]
    assert (l2.run_dir / "upstream" / "nfcore_rnaseq" / "salmon" / "salmon.merged.gene_counts.tsv").is_file()
    # The realistic difference that used to make reuse fail ...
    assert {item["sha256_source"] for item in _input_files(l1.run_dir)} == {"computed"}
    assert {item["sha256_source"] for item in _input_files(l2.run_dir)} == {"cache"}
    # ... while every scientific FASTQ identity is identical.
    strip = lambda files: [{k: v for k, v in item.items() if k != "sha256_source"} for item in files]  # noqa: E731
    assert strip(_input_files(l1.run_dir)) == strip(_input_files(l2.run_dir))


def test_l2_reuse_is_still_refused_when_a_fastq_identity_changed(monkeypatch, l1_project):
    """Ignoring checksum provenance must not ignore the checksum itself."""

    root, l1 = l1_project
    manifest_path = l1.run_dir / "frozen" / "input_manifest.yaml"
    manifest = yaml.safe_load(manifest_path.read_text(encoding="utf-8"))
    manifest["input"]["files"][0]["sha256"] = "0" * 64
    manifest_path.write_text(yaml.safe_dump(manifest, sort_keys=False), encoding="utf-8")
    report = validate_project(root)
    generate_plan(report)
    observed: list[list[str]] = []
    monkeypatch.setattr("rnaseq.service._run_command", _fake_nextflow(observed))
    with pytest.raises(ExecutionPreflightError, match="input_manifest.yaml differs"):
        execute_service_run(report, case_id="L2-001", reuse_upstream=f"L1-001/{l1.run_id}")
    assert observed == []
