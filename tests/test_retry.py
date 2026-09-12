"""Focused immutable retry regression coverage without biological workflow runs."""

from __future__ import annotations

import json
from pathlib import Path

import pytest
import yaml
from typer.testing import CliRunner

from rnaseq.cli import app
from rnaseq.errors import ExecutionPreflightError, UpstreamExecutionError
from rnaseq.execution import RuntimeCheck
from rnaseq.planner import generate_plan
from rnaseq.service import (
    _write_state,
    create_case_run,
    execute_retry_service_run,
    freeze_case_inputs,
    verify_delivery_manifest,
)
from rnaseq.validators import validate_project


runner = CliRunner()
pytestmark = pytest.mark.usefixtures("production_capable_execution_capacity")


def _failed_frozen_run(project_factory):
    root = project_factory()
    report = validate_project(root)
    assert report.is_valid
    generate_plan(report)
    run = create_case_run(report, "CASE-RETRY")
    freeze_case_inputs(report, run, profile="local", command=["rnaseq", "run", str(root), "--case-id", "CASE-RETRY"])
    (run.run_dir / "provenance" / "run_provenance.yaml").write_text("design: {}\n", encoding="utf-8")
    _write_state(run, "FAILED", phase="downstream", error="synthetic downstream failure")
    return root, run


def _mock_runtime(monkeypatch, root: Path):
    monkeypatch.setenv("RNASEQ_EXECUTION_ROOT", str(root.parent / "retry-nextflow-cache"))
    monkeypatch.setattr("rnaseq.service.check_nextflow", lambda: RuntimeCheck("Nextflow", "FOUND", "test"))
    monkeypatch.setattr("rnaseq.service.check_docker", lambda: RuntimeCheck("Docker", "FOUND", "test"))
    monkeypatch.setattr("rnaseq.service.check_container_runtime", lambda *_args: RuntimeCheck("Control-plane container", "FOUND", "test"))


def _successful_downstream(command, *, cwd, stdout_path, stderr_path):
    assert command[0:2] == ["nextflow", "run"]
    outdir = Path(command[command.index("--outdir") + 1])
    (outdir / "report").mkdir(parents=True)
    (outdir / "report" / "report.html").write_text("<html>retry</html>", encoding="utf-8")
    stdout_path.write_text("retry stdout\n", encoding="utf-8")
    stderr_path.write_text("retry stderr\n", encoding="utf-8")
    return 0


def test_retry_creates_new_run_and_preserves_frozen_scientific_contract(monkeypatch, project_factory):
    root, source = _failed_frozen_run(project_factory)
    _mock_runtime(monkeypatch, root)
    monkeypatch.setattr("rnaseq.service._run_command", _successful_downstream)
    original = {
        name: (source.run_dir / "frozen" / name).read_bytes()
        for name in ("project.yaml", "metadata.csv", "contrasts.csv", "input_manifest.yaml", "downstream_contract.json")
    }

    # These mutations are deliberately never read by retry.
    (root / "project.yaml").write_text("not: a retry contract\n", encoding="utf-8")
    (root / "metadata.csv").write_text("sample_id,condition\nchanged,Changed\n", encoding="utf-8")
    (root / "contrasts.csv").write_text("contrast_id,factor,numerator,denominator\nchanged,condition,A,B\n", encoding="utf-8")
    (root / "input" / "counts.csv").write_text("gene_id,changed\nGeneA,1\n", encoding="utf-8")

    retry = execute_retry_service_run(root, retry_of=f"{source.case_id}/{source.run_id}")

    assert retry.run_id != source.run_id
    assert retry.run_dir.parent == source.run_dir.parent
    assert json.loads(source.state_path.read_text(encoding="utf-8"))["status"] == "FAILED"
    assert all((source.run_dir / "frozen" / name).read_bytes() == value for name, value in original.items())
    state = json.loads(retry.state_path.read_text(encoding="utf-8"))
    assert state["status"] == "SUCCESS"
    assert state["attempt_type"] == "RETRY"
    assert state["retry_of"] == {"case_id": source.case_id, "run_id": source.run_id, "status": "FAILED"}
    for name in ("project.yaml", "metadata.csv", "contrasts.csv", "input_manifest.yaml"):
        assert (retry.run_dir / "frozen" / name).read_bytes() == original[name]
    contract = json.loads((retry.run_dir / "frozen" / "downstream_contract.json").read_text(encoding="utf-8"))
    assert contract["case"]["run_id"] == retry.run_id
    assert contract["project_config"] == str((retry.run_dir / "frozen" / "project.yaml").resolve())
    provenance = yaml.safe_load((retry.run_dir / "provenance" / "run_provenance.yaml").read_text(encoding="utf-8"))
    assert provenance["retry"]["retry_of"] == {"case_id": source.case_id, "run_id": source.run_id}
    assert provenance["retry"]["nextflow_resume_requested"] is False
    assert (retry.run_dir / "delivery" / "delivery_manifest.yaml").is_file()
    assert verify_delivery_manifest(retry.run_dir / "delivery") == ()
    assert not (source.run_dir / "delivery" / "delivery_manifest.yaml").exists()


def test_retry_rejects_success_and_malformed_source_before_creating_attempt(monkeypatch, project_factory):
    root, source = _failed_frozen_run(project_factory)
    _mock_runtime(monkeypatch, root)
    before = sorted((root / "runs" / source.case_id).iterdir())
    _write_state(source, "SUCCESS")
    with pytest.raises(ExecutionPreflightError, match="Only FAILED"):
        execute_retry_service_run(root, retry_of=f"{source.case_id}/{source.run_id}")
    assert sorted((root / "runs" / source.case_id).iterdir()) == before

    _write_state(source, "FAILED")
    (source.run_dir / "frozen" / "project.yaml").unlink()
    with pytest.raises((ExecutionPreflightError, UpstreamExecutionError), match="frozen"):
        execute_retry_service_run(root, retry_of=f"{source.case_id}/{source.run_id}")
    assert sorted((root / "runs" / source.case_id).iterdir()) == before


def test_retry_failure_preserves_both_attempts_and_resume_is_opt_in(monkeypatch, project_factory):
    root, source = _failed_frozen_run(project_factory)
    _mock_runtime(monkeypatch, root)
    seen: list[list[str]] = []

    def failing_downstream(command, *, cwd, stdout_path, stderr_path):
        seen.append(command)
        stdout_path.write_text("retry stdout\n", encoding="utf-8")
        stderr_path.write_text("retry stderr\n", encoding="utf-8")
        return 23

    monkeypatch.setattr("rnaseq.service._run_command", failing_downstream)
    with pytest.raises(UpstreamExecutionError, match="return code 23"):
        execute_retry_service_run(root, retry_of=f"{source.case_id}/{source.run_id}", nextflow_resume=True)

    attempts = sorted((root / "runs" / source.case_id).glob("*/run_state.json"))
    assert len(attempts) == 2
    retry_state = json.loads(next(path for path in attempts if path.parent != source.run_dir).read_text(encoding="utf-8"))
    assert retry_state["status"] == "FAILED"
    assert json.loads(source.state_path.read_text(encoding="utf-8"))["status"] == "FAILED"
    assert (source.run_dir / "logs").is_dir()
    assert all("-resume" in command for command in seen)
    assert "-resume" not in (json.loads(source.state_path.read_text(encoding="utf-8")).get("command") or [])


def test_retry_restarts_failed_upstream_from_frozen_fastqs_and_only_then_uses_resume(monkeypatch, project_factory):
    from conftest import base_config

    config = base_config()
    config["input"] = {"type": "fastq", "path": "input/fastq", "layout": "paired_end"}
    config["upstream"] = {
        "engine": "nfcore_rnaseq", "pipeline_version": "3.26.0", "aligner": None,
        "strandedness": "auto", "quantification": {"method": "salmon"},
    }
    config["reference"] = {"source": "igenomes", "genome": "test_reference"}
    root = project_factory(config=config)
    fastq = root / "input" / "fastq"
    fastq.mkdir()
    for sample in ("C1", "C2", "T1", "T2"):
        for read in ("R1", "R2"):
            (fastq / f"{sample}_{read}.fastq.gz").write_bytes(b"frozen fixture")
    (root / "metadata.csv").write_text(
        "sample_id,condition\nC1,Control\nC2,Control\nT1,Treatment\nT2,Treatment\n", encoding="utf-8"
    )
    report = validate_project(root)
    assert report.is_valid, report.errors
    generate_plan(report)
    source = create_case_run(report, "CASE-UPSTREAM-RETRY")
    freeze_case_inputs(report, source, profile="local", command=["rnaseq", "run"])
    (source.run_dir / "provenance" / "run_provenance.yaml").write_text("design: {}\n", encoding="utf-8")
    _write_state(source, "FAILED", phase="upstream", error="synthetic upstream failure")
    _mock_runtime(monkeypatch, root)
    seen: list[list[str]] = []

    def successful_retry(command, *, cwd, stdout_path, stderr_path):
        seen.append(command)
        stdout_path.write_text("retry stdout\n", encoding="utf-8")
        stderr_path.write_text("retry stderr\n", encoding="utf-8")
        outdir = Path(command[command.index("--outdir") + 1])
        if "nf-core/rnaseq" in command:
            salmon = outdir / "salmon"
            salmon.mkdir(parents=True)
            (salmon / "salmon.merged.gene_counts.tsv").write_text("gene_id\tC1\tC2\tT1\tT2\nGeneA\t1\t2\t3\t4\n", encoding="utf-8")
            (salmon / "salmon.merged.tx2gene_augmented.tsv").write_text("transcript_id\tgene_id\nTx1\tGeneA\n", encoding="utf-8")
            for sample in ("C1", "C2", "T1", "T2"):
                sample_dir = salmon / sample
                sample_dir.mkdir()
                (sample_dir / "quant.sf").write_text("Name\tLength\tEffectiveLength\tTPM\tNumReads\nTx1\t100\t80\t1\t1\n", encoding="utf-8")
            (outdir / "multiqc" / "multiqc_data").mkdir(parents=True)
            (outdir / "multiqc" / "multiqc_report.html").write_text("<html></html>", encoding="utf-8")
        else:
            (outdir / "report").mkdir(parents=True)
            (outdir / "report" / "report.html").write_text("<html></html>", encoding="utf-8")
        return 0

    monkeypatch.setattr("rnaseq.service._run_command", successful_retry)
    retry = execute_retry_service_run(root, retry_of=f"{source.case_id}/{source.run_id}", nextflow_resume=True)
    assert json.loads(retry.state_path.read_text(encoding="utf-8"))["status"] == "SUCCESS"
    assert json.loads(source.state_path.read_text(encoding="utf-8"))["status"] == "FAILED"
    assert len(seen) == 2
    assert "nf-core/rnaseq" in seen[0]
    assert all("-resume" in command for command in seen)
    assert (retry.run_dir / "logs" / "upstream.stdout.log").is_file()
    assert (retry.run_dir / "logs" / "downstream.stdout.log").is_file()


def test_status_labels_retry_attempt(monkeypatch, project_factory):
    root, source = _failed_frozen_run(project_factory)
    _mock_runtime(monkeypatch, root)
    monkeypatch.setattr("rnaseq.service._run_command", _successful_downstream)
    execute_retry_service_run(root, retry_of=f"{source.case_id}/{source.run_id}")

    result = runner.invoke(app, ["status", str(root)])
    assert result.exit_code == 0, result.output
    assert "Status: FAILED (original run)" in result.output
    assert "Status: SUCCESS (retry attempt)" in result.output
    assert f"Retry of: {source.case_id}/{source.run_id} (source status: FAILED)" in result.output
