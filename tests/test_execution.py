from __future__ import annotations

import csv
import json
from copy import deepcopy
from pathlib import Path
from types import SimpleNamespace

import yaml
import pytest
from typer.testing import CliRunner

from conftest import BASE_CONTRASTS, base_config
from rnaseq.cli import app
from rnaseq.errors import ExecutionPreflightError, UpstreamExecutionError
from rnaseq.execution import (
    PreparedRun,
    RESOURCE_CONTRACTS,
    RuntimeCheck,
    RuntimeSnapshot,
    ResourceContract,
    LocalResourceCapacity,
    build_nextflow_command,
    classify_execution_failure,
    detect_local_resource_capacity,
    _host_memory_bytes,
    _reference_runtime_check,
    execute_prepared_run,
    load_run_states,
    prepare_run,
    render_local_resource_config,
    render_upstream_conda_config,
    suggested_local_resources,
    validate_local_execution_budget,
    resolve_execution_workspace,
    effective_resource_budget,
    validate_effective_resource_budget,
)
from rnaseq.planner import generate_plan
from rnaseq.validators import validate_project

runner = CliRunner()


def test_experimental_osx_arm64_nfcore_conda_overrides_are_withdrawn(monkeypatch, tmp_path):
    """macOS is unsupported; the experimental per-process Conda overrides must not return."""
    monkeypatch.setattr("rnaseq.execution.native_platform", lambda: "osx-arm64")
    config = render_upstream_conda_config(tmp_path / "cache")
    assert "EAUTILS_GTF2BED" not in config and "TXIMETA_TXIMPORT" not in config
    assert "withName:" not in config


@pytest.mark.parametrize("platform", ["linux-64", "osx-64"])
def test_upstream_conda_config_has_no_process_overrides(monkeypatch, tmp_path, platform):
    monkeypatch.setattr("rnaseq.execution.native_platform", lambda: platform)
    config = render_upstream_conda_config(tmp_path / "cache")
    assert "withName:" not in config
    assert "process {" not in config
    assert "conda.enabled = true" in config


def _ready_fastq_project(tmp_path: Path) -> Path:
    root = tmp_path / "ready_fastq"
    fastq = root / "input" / "fastq"
    fastq.mkdir(parents=True)
    for sample in ("C1", "C2", "T1", "T2"):
        for read in ("R1", "R2"):
            (fastq / f"{sample}_{read}.fastq.gz").write_bytes(b"fixture")
    (root / "metadata.csv").write_text(
        "sample_id,condition\nC1,Control\nC2,Control\nT1,Treatment\nT2,Treatment\n",
        encoding="utf-8",
    )
    (root / "contrasts.csv").write_text(BASE_CONTRASTS, encoding="utf-8")
    config = deepcopy(base_config())
    config["input"] = {"type": "fastq", "path": "input/fastq", "layout": "paired_end"}
    config["upstream"] = {
        "engine": "nfcore_rnaseq",
        "pipeline_version": "3.26.0",
        "aligner": None,
        "strandedness": "auto",
        "quantification": {"method": "salmon"},
    }
    config["reference"] = {"source": "igenomes", "genome": "test_reference"}
    (root / "project.yaml").write_text(yaml.safe_dump(config, sort_keys=False), encoding="utf-8")
    report = validate_project(root)
    assert report.is_valid
    generate_plan(report)
    return root


def _prepared(monkeypatch, root: Path, production_capable_execution_capacity) -> PreparedRun:
    monkeypatch.setenv("RNASEQ_EXECUTION_ROOT", str(root.parent / "local-execution-root"))
    monkeypatch.setattr("rnaseq.execution.check_nextflow", lambda: RuntimeCheck("Nextflow", "FOUND", "25.10.4"))
    return prepare_run(validate_project(root), "local")


def test_raw_count_run_is_rejected(project_factory):
    with pytest.raises(ExecutionPreflightError, match="Raw-count execution"):
        prepare_run(validate_project(project_factory()), "local")


def test_mocked_execution_capacity_fixture_supplies_a_production_capable_runtime_snapshot(
    production_capable_execution_capacity,
):
    snapshot = production_capable_execution_capacity
    resources = effective_resource_budget(snapshot, ResourceContract("PROJECT_LOCAL", 8, 12, 12))
    assert snapshot.logical_cpus == 16
    assert (resources.effective_cpus, resources.effective_memory_gib) == (8, 12)
    validate_effective_resource_budget(resources)


def test_preflight_rejects_server_and_stale_plan(monkeypatch, tmp_path, production_capable_execution_capacity):
    root = _ready_fastq_project(tmp_path)
    report = validate_project(root)
    with pytest.raises(ExecutionPreflightError, match="server"):
        prepare_run(report, "server")
    (root / "metadata.csv").write_text(
        "sample_id,condition\nC1,Control\nC2,Control\nT1,Treatment\nT2,Treatment\n\n", encoding="utf-8"
    )
    # Byte changes alone invalidate the checksum-bearing frozen plan.
    with pytest.raises(ExecutionPreflightError, match="changed since planning"):
        prepare_run(validate_project(root), "local")


def test_preflight_requires_reference_and_nextflow(monkeypatch, tmp_path, production_capable_execution_capacity):
    root = _ready_fastq_project(tmp_path)
    config = yaml.safe_load((root / "project.yaml").read_text())
    config["reference"]["genome"] = None
    (root / "project.yaml").write_text(yaml.safe_dump(config, sort_keys=False), encoding="utf-8")
    with pytest.raises(ExecutionPreflightError, match="reference genome not configured"):
        prepare_run(validate_project(root), "local")

    root = _ready_fastq_project(tmp_path / "nextflow")
    monkeypatch.setattr(
        "rnaseq.execution.check_nextflow",
        lambda: RuntimeCheck("Nextflow", "NOT FOUND", "Nextflow executable was not found on PATH."),
    )
    with pytest.raises(ExecutionPreflightError, match="Nextflow is required"):
        prepare_run(validate_project(root), "local")


def test_fastq_change_invalidates_plan(monkeypatch, tmp_path, production_capable_execution_capacity):
    root = _ready_fastq_project(tmp_path)
    with (root / "input" / "fastq" / "C1_R1.fastq.gz").open("ab") as handle:
        handle.write(b"changed")
    with pytest.raises(ExecutionPreflightError, match="changed since planning"):
        prepare_run(validate_project(root), "local")


def test_safe_pinned_command(monkeypatch, tmp_path, production_capable_execution_capacity):
    root = _ready_fastq_project(tmp_path)
    prepared = _prepared(monkeypatch, root, production_capable_execution_capacity)
    command = build_nextflow_command(
        prepared.report,
        samplesheet=root / "planning" / "samplesheet.csv",
        output_dir=root / "runs" / "out",
        profile="local",
        work_dir=root.parent / "local-execution-root" / "work" / "upstream",
    )
    assert command[:7] == ["nextflow", "run", "nf-core/rnaseq", "-r", "3.26.0", "-profile", "conda"]
    assert "--pseudo_aligner" in command and "salmon" in command
    assert "--skip_alignment" not in command
    assert command[command.index("-work-dir") + 1].endswith("local-execution-root/work/upstream")
    assert "test_reference" in command
    assert prepared.container_runtime == "conda"


def test_successful_mocked_execution_freezes_state_and_handoff(monkeypatch, tmp_path, production_capable_execution_capacity):
    root = _ready_fastq_project(tmp_path)
    prepared = _prepared(monkeypatch, root, production_capable_execution_capacity)

    def successful(arguments, **_kwargs):
        outdir = Path(arguments[arguments.index("--outdir") + 1])
        (outdir / "salmon").mkdir(parents=True)
        (outdir / "salmon" / "salmon.merged.gene_counts.tsv").write_text("gene_id\tC1\nGeneA\t1.0\n")
        (outdir / "salmon" / "salmon.merged.tx2gene_augmented.tsv").write_text("transcript_id\tgene_id\nTx1\tGeneA\n")
        for sample in ("C1", "C2", "T1", "T2"):
            sample_dir = outdir / "salmon" / sample
            sample_dir.mkdir()
            (sample_dir / "quant.sf").write_text("Name\tLength\tEffectiveLength\tTPM\tNumReads\nTx1\t100\t80\t1\t1\n")
        (outdir / "multiqc" / "salmon" / "multiqc_data").mkdir(parents=True)
        (outdir / "multiqc" / "salmon" / "multiqc_report.html").write_text("<html></html>")
        return SimpleNamespace(returncode=0)

    monkeypatch.setattr("rnaseq.execution.subprocess.run", successful)
    result = execute_prepared_run(prepared)
    state = json.loads(result.state_path.read_text())
    handoff = yaml.safe_load(result.handoff_path.read_text())
    assert state["status"] == "SUCCESS"
    assert (result.run_dir / "frozen" / "samplesheet.csv").is_file()
    assert json.loads((result.run_dir / "frozen" / "nextflow.params.json").read_text()) == {
        "skip_alignment": True
    }
    runtime_params = json.loads((result.run_dir / "frozen" / "nextflow.params.json").read_text())
    assert runtime_params["skip_alignment"] is True
    assert "skip_trimming" not in runtime_params
    resource_config = (result.run_dir / "frozen" / "local.nextflow.config").read_text()
    effective = yaml.safe_load((result.run_dir / "provenance" / "run_provenance.yaml").read_text())["runtime_resources"]["effective"]
    assert f"memory: '{effective['memory_gib']}.GB'" in resource_config
    assert "SALMON_QUANT" not in resource_config and "maxForks" not in resource_config
    assert handoff["gene_level_counts"]["format"] == "TSV"
    assert handoff["gene_level_counts"]["identifier_column"] == "gene_id"
    assert sorted(handoff["salmon"]["quant_sf"]) == ["C1", "C2", "T1", "T2"]
    assert handoff["salmon"]["tx2gene"]["mapping_type"] == "nfcore_tx2gene_augmented"
    assert handoff["salmon"]["tx2gene"]["path"].endswith("salmon.merged.tx2gene_augmented.tsv")
    assert len(handoff["salmon"]["tx2gene"]["sha256"]) == 64
    assert handoff["multiqc"]["html"].endswith("multiqc_report.html")
    assert load_run_states(root)[0]["handoff_available"] is True


def test_handoff_accepts_modern_multiqc_report_data(monkeypatch, tmp_path, production_capable_execution_capacity):
    root = _ready_fastq_project(tmp_path)
    prepared = _prepared(monkeypatch, root, production_capable_execution_capacity)

    def successful(arguments, **_kwargs):
        outdir = Path(arguments[arguments.index("--outdir") + 1])
        (outdir / "salmon").mkdir(parents=True)
        (outdir / "salmon" / "salmon.merged.gene_counts.tsv").write_text("gene_id\tC1\nGeneA\t1.0\n")
        (outdir / "salmon" / "salmon.merged.tx2gene_augmented.tsv").write_text("transcript_id\tgene_id\nTx1\tGeneA\n")
        for sample in ("C1", "C2", "T1", "T2"):
            sample_dir = outdir / "salmon" / sample
            sample_dir.mkdir()
            (sample_dir / "quant.sf").write_text("Name\tLength\tEffectiveLength\tTPM\tNumReads\nTx1\t100\t80\t1\t1\n")
        (outdir / "multiqc" / "multiqc_report_data").mkdir(parents=True)
        (outdir / "multiqc" / "multiqc_report.html").write_text("<html></html>")
        return SimpleNamespace(returncode=0)

    monkeypatch.setattr("rnaseq.execution.subprocess.run", successful)
    result = execute_prepared_run(prepared)
    handoff = yaml.safe_load(result.handoff_path.read_text())
    assert handoff["multiqc"]["data_directory"].endswith("multiqc_report_data")


def test_augmented_tx2gene_fixture_retains_nfcore_self_mapping():
    fixture = Path(__file__).parent / "fixtures" / "salmon_augmented_tx2gene"
    with (fixture / "S1" / "quant.sf").open() as handle:
        quant = {row["Name"]: float(row["NumReads"]) for row in csv.DictReader(handle, delimiter="\t")}
    with (fixture / "ordinary.tsv").open() as handle:
        ordinary = {row["transcript_id"] for row in csv.DictReader(handle, delimiter="\t")}
    with (fixture / "augmented.tsv").open() as handle:
        augmented = {row["transcript_id"]: row["gene_id"] for row in csv.DictReader(handle, delimiter="\t")}
    assert sum(quant[name] for name in ordinary) == 10
    assert augmented["TxOrphan"] == "TxOrphan"
    assert sum(quant[name] for name in augmented) == 15


def test_failed_subprocess_preserves_run_and_marks_failed(monkeypatch, tmp_path, production_capable_execution_capacity):
    root = _ready_fastq_project(tmp_path)
    prepared = _prepared(monkeypatch, root, production_capable_execution_capacity)
    monkeypatch.setattr("rnaseq.execution.subprocess.run", lambda *_args, **_kwargs: SimpleNamespace(returncode=23))
    with pytest.raises(UpstreamExecutionError, match="return code 23"):
        execute_prepared_run(prepared)
    state_path = next((root / "runs").glob("*/run_state.json"))
    assert json.loads(state_path.read_text())["status"] == "FAILED"


def test_execution_root_override_is_local_and_portable(monkeypatch, tmp_path):
    root = tmp_path / "local-cache"
    monkeypatch.setenv("RNASEQ_EXECUTION_ROOT", str(root))
    workspace = resolve_execution_workspace("CASE-001", "20260828-120000+0800")
    assert workspace.root == root / "CASE-001" / "20260828-120000+0800"
    assert workspace.launch_dir == workspace.root / "launch"
    assert workspace.work_dir == workspace.root / "work"


def test_project_doctor_surfaces_a_missing_adopted_reference_error(monkeypatch, tmp_path):
    failed_report = SimpleNamespace(
        is_valid=False,
        config=None,
        errors=[SimpleNamespace(message="Local Salmon index is missing or incomplete.")],
    )
    monkeypatch.setattr("rnaseq.validators.validate_project", lambda _project: failed_report)
    check = _reference_runtime_check(tmp_path)
    assert check.verdict == "FAIL"
    assert "Local Salmon index is missing or incomplete." in check.detail


def test_runtime_failure_classification_is_actionable_and_does_not_change_resources(tmp_path):
    stderr = tmp_path / "nextflow.stderr.log"
    stderr.write_text("Rscript: Killed\n", encoding="utf-8")
    diagnostic = classify_execution_failure("SALMON_QUANT", 137, stderr, resource=RESOURCE_CONTRACTS["MEDIUM"])
    assert "LIKELY_OOM" in diagnostic
    assert "exit_code=137" in diagnostic
    assert "4 CPUs, 8 GiB, 8 h" in diagnostic
    assert RESOURCE_CONTRACTS["MEDIUM"].memory_gib == 8


def test_rendered_resource_contract_uses_only_the_aggregate_ceiling(tmp_path):
    rendered = render_local_resource_config()
    assert "executor { cpus = 8; memory = '12.GB' }" in rendered
    assert "resourceLimits = [cpus: 8, memory: '12.GB', time: '12.h']" in rendered
    assert "withLabel" not in rendered and "withName" not in rendered and "maxForks" not in rendered
    assert "skip_alignment" not in rendered and "skip_trimming" not in rendered


def test_local_resource_suggestion_and_validation_are_conservative():
    capacity = LocalResourceCapacity(20, 64, 62)
    assert suggested_local_resources(capacity) == (16, 48)
    validate_local_execution_budget(16, 48, capacity)
    with pytest.raises(ExecutionPreflightError, match="positive"):
        validate_local_execution_budget(0, 48, capacity)
    validate_local_execution_budget(21, 48, capacity)
    assert suggested_local_resources(LocalResourceCapacity(None, None, None)) == (8, 12)


def test_effective_budget_clamps_the_project_budget_to_host_capacity():
    requested = ResourceContract("PROJECT_LOCAL", 24, 48, 12)
    linux = effective_resource_budget(RuntimeSnapshot("Linux", "amd64", 20, 64 * 1024**3), requested)
    assert (linux.effective_cpus, linux.effective_memory_gib) == (20, 48)
    assert linux.clamped is True
    # The retired Docker VM ceiling stays in the record schema, always inert.
    assert linux.as_dict()["container_runtime"] == {"cpus": None, "memory_gib": None, "ceiling_applies": False}


def test_effective_budget_keeps_smaller_project_budget_and_handles_missing_host_detection():
    project = ResourceContract("PROJECT_LOCAL", 8, 12, 12)
    resources = effective_resource_budget(RuntimeSnapshot("Linux", "amd64", None, None), project)
    assert (resources.effective_cpus, resources.effective_memory_gib) == (8, 12)
    assert resources.clamped is False
    assert len(resources.warnings) == 2
    validate_effective_resource_budget(resources)


def test_effective_budget_fails_when_largest_process_cannot_fit():
    resources = effective_resource_budget(
        RuntimeSnapshot("Linux", "amd64", 12, 10 * 1024**3), ResourceContract("PROJECT_LOCAL", 16, 32, 12)
    )
    with pytest.raises(ExecutionPreflightError, match="at least 8 CPUs/12 GiB"):
        validate_effective_resource_budget(resources)


def test_host_resource_detection_uses_linux_procfs_and_darwin_sysctl_without_procps(monkeypatch):
    gib = 1024**3
    monkeypatch.setattr("rnaseq.execution.platform.system", lambda: "Linux")
    monkeypatch.setattr("rnaseq.execution._linux_memory_bytes", lambda: (64 * gib, 48 * gib))
    linux = detect_local_resource_capacity()
    assert linux.total_memory_gib == 64
    assert linux.available_memory_gib == 48
    assert _host_memory_bytes() == 64 * gib

    monkeypatch.setattr("rnaseq.execution.platform.system", lambda: "Darwin")
    monkeypatch.setattr("rnaseq.execution._darwin_memory_bytes", lambda: 32 * gib)
    darwin = detect_local_resource_capacity()
    assert darwin.total_memory_gib == 32
    assert darwin.available_memory_gib == 32
    assert _host_memory_bytes() == 32 * gib


def test_cli_run_requires_explicit_confirmation(monkeypatch, tmp_path):
    root = _ready_fastq_project(tmp_path)
    monkeypatch.setattr("rnaseq.cli.prepare_service_run", lambda report, profile: None)
    monkeypatch.setattr("rnaseq.cli.execute_service_run", lambda *_args, **_kwargs: (_ for _ in ()).throw(AssertionError("must not execute")))
    result = runner.invoke(app, ["run", str(root), "--case-id", "CASE-001"], input="n\n")
    assert result.exit_code == 0, result.output
    assert "Execution cancelled" in result.output
