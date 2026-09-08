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
    LocalResourceCapacity,
    build_nextflow_command,
    classify_execution_failure,
    check_container_runtime,
    downstream_docker_user_mapping,
    downstream_docker_user_mapping_check,
    doctor_checks,
    detect_local_resource_capacity,
    _host_memory_bytes,
    _reference_runtime_check,
    execute_prepared_run,
    load_run_states,
    prepare_run,
    render_local_resource_config,
    suggested_local_resources,
    validate_local_execution_budget,
    resolve_execution_workspace,
    runtime_resource_checks,
)
from rnaseq.planner import generate_plan
from rnaseq.validators import validate_project

runner = CliRunner()


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


def _prepared(monkeypatch, root: Path) -> PreparedRun:
    monkeypatch.setenv("RNASEQ_EXECUTION_ROOT", str(root.parent / "local-execution-root"))
    monkeypatch.setattr("rnaseq.execution.check_nextflow", lambda: RuntimeCheck("Nextflow", "FOUND", "25.10.4"))
    monkeypatch.setattr("rnaseq.execution.check_docker", lambda: RuntimeCheck("Docker", "FOUND", "Docker daemon is available."))
    return prepare_run(validate_project(root), "local")


def test_raw_count_run_is_rejected(project_factory):
    with pytest.raises(ExecutionPreflightError, match="Raw-count execution"):
        prepare_run(validate_project(project_factory()), "local")


def test_preflight_rejects_server_and_stale_plan(monkeypatch, tmp_path):
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


def test_preflight_requires_reference_and_nextflow(monkeypatch, tmp_path):
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


def test_fastq_change_invalidates_plan(monkeypatch, tmp_path):
    root = _ready_fastq_project(tmp_path)
    with (root / "input" / "fastq" / "C1_R1.fastq.gz").open("ab") as handle:
        handle.write(b"changed")
    with pytest.raises(ExecutionPreflightError, match="changed since planning"):
        prepare_run(validate_project(root), "local")


def test_safe_pinned_command(monkeypatch, tmp_path):
    root = _ready_fastq_project(tmp_path)
    prepared = _prepared(monkeypatch, root)
    command = build_nextflow_command(
        prepared.report,
        samplesheet=root / "planning" / "samplesheet.csv",
        output_dir=root / "runs" / "out",
        profile="local",
        work_dir=root.parent / "local-execution-root" / "work" / "upstream",
    )
    assert command[:7] == ["nextflow", "run", "nf-core/rnaseq", "-r", "3.26.0", "-profile", "docker"]
    assert "--pseudo_aligner" in command and "salmon" in command
    assert "--skip_alignment" not in command
    assert command[command.index("-work-dir") + 1].endswith("local-execution-root/work/upstream")
    assert "test_reference" in command


def test_successful_mocked_execution_freezes_state_and_handoff(monkeypatch, tmp_path):
    root = _ready_fastq_project(tmp_path)
    prepared = _prepared(monkeypatch, root)

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
    monkeypatch.setattr(
        "rnaseq.execution.runtime_snapshot",
        lambda *_args: RuntimeSnapshot("Darwin", "arm64", 12, 24 * 1024**3, "arm64", 15 * 1024**3, "test", "arm64"),
    )
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
    assert "memory: '12.GB'" in resource_config
    assert "SALMON_QUANT" in resource_config and "maxForks = 1" in resource_config
    assert handoff["gene_level_counts"]["format"] == "TSV"
    assert handoff["gene_level_counts"]["identifier_column"] == "gene_id"
    assert sorted(handoff["salmon"]["quant_sf"]) == ["C1", "C2", "T1", "T2"]
    assert handoff["salmon"]["tx2gene"]["mapping_type"] == "nfcore_tx2gene_augmented"
    assert handoff["salmon"]["tx2gene"]["path"].endswith("salmon.merged.tx2gene_augmented.tsv")
    assert len(handoff["salmon"]["tx2gene"]["sha256"]) == 64
    assert handoff["multiqc"]["html"].endswith("multiqc_report.html")
    assert load_run_states(root)[0]["handoff_available"] is True


def test_handoff_accepts_modern_multiqc_report_data(monkeypatch, tmp_path):
    root = _ready_fastq_project(tmp_path)
    prepared = _prepared(monkeypatch, root)

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
    monkeypatch.setattr(
        "rnaseq.execution.runtime_snapshot",
        lambda *_args: RuntimeSnapshot("Darwin", "arm64", 12, 24 * 1024**3, "arm64", 15 * 1024**3, "test", "arm64"),
    )
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


def test_failed_subprocess_preserves_run_and_marks_failed(monkeypatch, tmp_path):
    root = _ready_fastq_project(tmp_path)
    prepared = _prepared(monkeypatch, root)
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


def test_container_runtime_probe_requires_python_r_and_r_packages(monkeypatch):
    calls: list[list[str]] = []

    def successful(arguments):
        calls.append(arguments)
        return SimpleNamespace(returncode=0, stdout="", stderr="")

    monkeypatch.setattr("rnaseq.execution.check_docker", lambda: RuntimeCheck("Docker", "FOUND", "available"))
    monkeypatch.setattr("rnaseq.execution._run_capture", successful)
    result = check_container_runtime()
    assert result.state == "FOUND"
    assert calls[0][:3] == ["docker", "image", "inspect"]
    probe = calls[1]
    assert probe[:6] == ["docker", "run", "--rm", "rnaseq-control-plane:latest", "sh", "-c"]
    assert "--entrypoint" not in probe
    assert "-lc" not in probe
    assert "for executable in python Rscript" in probe[-1]
    assert "command -v \"$executable\"" in probe[-1]
    assert "DESeq2" in probe[-1] and "org.Mm.eg.db" in probe[-1]
    assert "python -m rnaseq.workflow_support report --help" in probe[-1]
    assert "--enrichment" in probe[-1]


def test_downstream_docker_user_mapping_is_dynamic_on_linux_wsl_and_absent_on_macos():
    assert downstream_docker_user_mapping(host_os="Linux", uid=24701, gid=24703) == "24701:24703"
    assert downstream_docker_user_mapping(host_os="Darwin", uid=24701, gid=24703) is None
    assert downstream_docker_user_mapping(host_os="Linux", uid=-1, gid=24703) is None
    source = Path("src/rnaseq/execution.py").read_text(encoding="utf-8")
    assert "1000:1000" not in source


def test_doctor_explains_linux_wsl_downstream_user_mapping(monkeypatch):
    monkeypatch.setattr("rnaseq.execution.downstream_docker_user_mapping", lambda: "24701:24703")
    check = downstream_docker_user_mapping_check()
    assert check.verdict == "PASS"
    assert "--user 24701:24703" in check.detail


def test_container_runtime_probe_reports_the_failed_prerequisite(monkeypatch):
    results = iter((
        SimpleNamespace(returncode=0, stdout="", stderr=""),
        SimpleNamespace(returncode=1, stdout="", stderr="missing executable: ps\n"),
    ))
    monkeypatch.setattr("rnaseq.execution.check_docker", lambda: RuntimeCheck("Docker", "FOUND", "available"))
    monkeypatch.setattr("rnaseq.execution._run_capture", lambda _arguments: next(results))
    result = check_container_runtime()
    assert result.state == "NOT FOUND"
    assert result.detail == "container prerequisite probe exited 1: stderr: missing executable: ps"


def test_doctor_reports_a_successful_container_probe(monkeypatch):
    monkeypatch.setattr("rnaseq.execution.check_nextflow", lambda: RuntimeCheck("Nextflow", "FOUND", "available"))
    monkeypatch.setattr("rnaseq.execution.check_docker", lambda: RuntimeCheck("Docker", "FOUND", "available"))
    monkeypatch.setattr("rnaseq.execution.check_container_runtime", lambda *_args: RuntimeCheck("Control-plane container", "FOUND", "available"))
    monkeypatch.setattr("rnaseq.downstream.r_runtime_checks", lambda: ())
    checks = doctor_checks()
    assert RuntimeCheck("Control-plane container", "FOUND", "available") in checks


def test_doctor_reports_requested_and_observed_image_identity(monkeypatch, project_factory):
    root = project_factory()
    config_path = root / "project.yaml"
    config = yaml.safe_load(config_path.read_text(encoding="utf-8"))
    config["runtime"] = {"control_plane_image": "rnaseq-control-plane:0.5.1"}
    config_path.write_text(yaml.safe_dump(config, sort_keys=False), encoding="utf-8")
    monkeypatch.setattr("rnaseq.execution.check_nextflow", lambda: RuntimeCheck("Nextflow", "FOUND", "available"))
    monkeypatch.setattr("rnaseq.execution.check_docker", lambda: RuntimeCheck("Docker", "FOUND", "available"))
    monkeypatch.setattr("rnaseq.execution.check_container_runtime", lambda *_args: RuntimeCheck("Control-plane container", "FOUND", "available"))
    monkeypatch.setattr(
        "rnaseq.execution.inspect_container_image",
        lambda image: {"reference": image, "image_id": "sha256:" + "a" * 64, "repo_digests": ["repo@sha256:" + "b" * 64]},
    )
    monkeypatch.setattr(
        "rnaseq.execution.runtime_snapshot",
        lambda *_args: RuntimeSnapshot("Darwin", "arm64", 12, 24 * 1024**3, "arm64", 15 * 1024**3, "test", "arm64"),
    )
    monkeypatch.setattr("rnaseq.downstream.r_runtime_checks", lambda: ())
    identity = {check.name: check for check in doctor_checks(root)}["Control-plane image identity"]
    assert identity.verdict == "PASS"
    assert "requested=rnaseq-control-plane:0.5.1" in identity.detail
    assert "observed_image_id=sha256:" in identity.detail


def test_doctor_distinguishes_missing_nextflow_and_docker_from_architecture_warnings(monkeypatch):
    monkeypatch.setattr("rnaseq.execution.check_nextflow", lambda: RuntimeCheck("Nextflow", "NOT FOUND", "not installed"))
    monkeypatch.setattr("rnaseq.execution.check_docker", lambda: RuntimeCheck("Docker", "NOT FOUND", "daemon unavailable"))
    monkeypatch.setattr("rnaseq.execution.check_container_runtime", lambda *_args: RuntimeCheck("Control-plane container", "NOT FOUND", "daemon unavailable"))
    monkeypatch.setattr(
        "rnaseq.execution.runtime_snapshot",
        lambda *_args: RuntimeSnapshot("Darwin", "arm64", 12, 24 * 1024**3, None, None, None, "amd64"),
    )
    monkeypatch.setattr("rnaseq.downstream.r_runtime_checks", lambda: ())

    by_name = {check.name: check for check in doctor_checks()}
    assert by_name["Nextflow"].verdict == "FAIL"
    assert by_name["Docker"].verdict == "FAIL"
    assert by_name["Docker runtime"].verdict == "FAIL"
    assert by_name["Control-plane image architecture"].verdict == "WARN"
    assert by_name["Selected local ceiling"].verdict == "WARN"
    assert "cannot be confirmed" in by_name["Selected local ceiling"].detail
    assert "architecture=arm64" in by_name["Host runtime"].detail


def test_runtime_doctor_warns_for_amd64_image_on_arm64_and_low_docker_memory():
    checks = runtime_resource_checks(RuntimeSnapshot(
        host_os="Darwin", host_architecture="arm64", logical_cpus=12, host_memory_bytes=24 * 1024**3,
        docker_architecture="arm64", docker_memory_bytes=8 * 1024**3, docker_version="28.0.1",
        control_plane_image_architecture="amd64", docker_cpus=12,
    ))
    by_name = {item.name: item for item in checks}
    assert by_name["Control-plane image architecture"].verdict == "WARN"
    assert "Rosetta" in by_name["Control-plane image architecture"].detail
    assert by_name["Selected local ceiling"].verdict == "WARN"


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


def test_rendered_resource_contract_bounds_salmon_concurrency_without_touching_nfcore_params(tmp_path):
    rendered = render_local_resource_config()
    assert "SMALL" in rendered and "MEDIUM" in rendered and "LARGE=8/12 GiB" in rendered
    assert "executor { cpus = 8; memory = '12.GB' }" in rendered
    assert "resourceLimits = [cpus: 8, memory: '12.GB', time: '12.h']" in rendered
    assert "SALMON_QUANT" in rendered and "maxForks = 1" in rendered
    assert "skip_alignment" not in rendered and "skip_trimming" not in rendered


def test_local_resource_suggestion_and_validation_are_conservative():
    capacity = LocalResourceCapacity(20, 64, 62)
    assert suggested_local_resources(capacity) == (16, 48)
    validate_local_execution_budget(16, 48, capacity)
    with pytest.raises(ExecutionPreflightError, match="positive"):
        validate_local_execution_budget(0, 48, capacity)
    with pytest.raises(ExecutionPreflightError, match="exceeds detected host"):
        validate_local_execution_budget(21, 48, capacity)
    assert suggested_local_resources(LocalResourceCapacity(None, None, None)) == (8, 12)


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
