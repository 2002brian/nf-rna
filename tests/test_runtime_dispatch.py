"""Runtime policy (v1.3.0+): linux-64/WSL2 -> Nextflow + Conda; every other host is refused."""

from __future__ import annotations

import json
from copy import deepcopy
from pathlib import Path
from types import SimpleNamespace

import pytest
import yaml

from conftest import base_config
from rnaseq.downstream_runtime import DownstreamRuntime
from rnaseq.errors import ExecutionPreflightError
from rnaseq import execution as execution_module
from rnaseq.execution import (
    BACKEND_CONDA,
    RuntimeCheck,
    build_hisat2_featurecounts_command,
    build_nextflow_command,
    execution_backend,
    render_upstream_conda_config,
)
from rnaseq.planner import generate_plan
from rnaseq.service import (
    _frozen_backend,
    build_downstream_nextflow_command,
    create_case_run,
    execute_service_run,
    freeze_case_inputs,
    prepare_service_run,
    write_downstream_runtime_config,
)
from rnaseq.validators import validate_project


pytestmark = pytest.mark.usefixtures("production_capable_execution_capacity")


@pytest.fixture
def mocked_downstream_runtime(monkeypatch, tmp_path):
    runtime = DownstreamRuntime(
        prefix=tmp_path / "native-runtime", platform="linux-64",
        lock_filename="nf-rna-downstream-linux-64.lock.yml", lock_sha256="a" * 64,
        wheel_filename="nf_rna-1.2.0-py3-none-any.whl", wheel_sha256="b" * 64,
        nf_rna_version="1.2.0", source_revision="test-revision",
        r_scripts=({"name": "l2_analysis.R", "sha256": "c" * 64},), r_scripts_sha256="d" * 64,
    )
    monkeypatch.setattr("rnaseq.service.downstream_runtime_preflight", lambda: None)
    monkeypatch.setattr("rnaseq.service.ensure_downstream_runtime", lambda: runtime)
    return runtime


LINUX = ("Linux", "x86_64")
DARWIN = ("Darwin", "arm64")


def _host(monkeypatch, system: str, machine: str, *, wsl: bool = False) -> None:
    monkeypatch.setattr("platform.system", lambda: system)
    monkeypatch.setattr("platform.machine", lambda: machine)
    monkeypatch.setattr("rnaseq.execution._is_wsl", lambda: wsl)


def _fail(message: str):
    def raiser(*_args, **_kwargs):
        raise AssertionError(message)
    return raiser


def _guard_conda(monkeypatch) -> None:
    """Mock the Conda backend's prerequisites and fail on any Docker invocation."""

    real_capture = execution_module._run_capture

    def no_docker(arguments):
        if arguments and Path(arguments[0]).name == "docker":
            raise AssertionError(f"Conda backend invoked Docker: {arguments}")
        return real_capture(arguments)

    monkeypatch.setattr("rnaseq.execution._run_capture", no_docker)
    monkeypatch.setattr("rnaseq.service.check_upstream_conda", lambda: RuntimeCheck("Conda", "FOUND", "conda 26.5.3"))
    monkeypatch.setattr("rnaseq.service.check_nextflow", lambda: RuntimeCheck("Nextflow", "FOUND", "26.04.6"))


# ---------------------------------------------------------------- policy


@pytest.mark.parametrize(("system", "machine", "wsl", "backend"), [
    ("Linux", "x86_64", False, BACKEND_CONDA),
    ("Linux", "x86_64", True, BACKEND_CONDA),
    ("Linux", "amd64", True, BACKEND_CONDA),
])
def test_backend_is_selected_automatically_from_the_host(monkeypatch, system, machine, wsl, backend):
    _host(monkeypatch, system, machine, wsl=wsl)
    assert execution_backend() == backend


@pytest.mark.parametrize("system", ["Windows", "CYGWIN_NT-10.0", "MSYS_NT-10.0"])
def test_native_windows_is_unsupported_and_directs_to_wsl2(monkeypatch, system):
    _host(monkeypatch, system, "AMD64")
    with pytest.raises(ExecutionPreflightError, match="Native Windows is not a supported nf-rna runtime.*WSL2"):
        execution_backend()


@pytest.mark.parametrize(("system", "machine"), [("Linux", "aarch64"), ("Darwin", "x86_64"), ("Darwin", "arm64")])
def test_other_platforms_are_unsupported(monkeypatch, system, machine):
    # macOS Apple Silicon Docker dispatch (last released in v1.2.1, unqualified in v1.3.0) is retired.
    _host(monkeypatch, system, machine)
    with pytest.raises(ExecutionPreflightError, match=r"Unsupported nf-rna runtime platform.*Linux x86-64 including WSL2 \(Nextflow \+ Conda\)\.$"):
        execution_backend()


def test_macos_run_preflight_is_refused(monkeypatch, project_factory):
    report = validate_project(project_factory())
    generate_plan(report)
    _host(monkeypatch, *DARWIN)
    monkeypatch.setattr("rnaseq.service.check_nextflow", lambda: RuntimeCheck("Nextflow", "FOUND", "26.04.6"))
    with pytest.raises(ExecutionPreflightError, match="Unsupported nf-rna runtime platform darwin-arm64"):
        prepare_service_run(report, profile="local")


def test_native_windows_run_preflight_directs_to_wsl2(monkeypatch, project_factory):
    report = validate_project(project_factory())
    generate_plan(report)
    _host(monkeypatch, "Windows", "AMD64")
    with pytest.raises(ExecutionPreflightError, match="WSL2"):
        prepare_service_run(report, profile="local")


# ---------------------------------------------------------------- upstream commands


def _salmon_report(root: Path):
    fastq = root / "input" / "fastq"
    fastq.mkdir(parents=True)
    for sample in ("C1", "C2", "T1", "T2"):
        for read in ("R1", "R2"):
            (fastq / f"{sample}_{read}.fastq.gz").write_bytes(b"x")
    (root / "metadata.csv").write_text("sample_id,condition\nC1,Control\nC2,Control\nT1,Treatment\nT2,Treatment\n", encoding="utf-8")
    (root / "contrasts.csv").write_text("contrast_id,factor,numerator,denominator\nTreatment_vs_Control,condition,Treatment,Control\n", encoding="utf-8")
    config = base_config()
    config["input"] = {"type": "fastq", "path": "input/fastq", "layout": "paired_end"}
    config["upstream"] = {"engine": "nfcore_rnaseq", "pipeline_version": "3.26.0", "aligner": None, "strandedness": "auto", "quantification": {"method": "salmon"}}
    config["reference"] = {"source": "igenomes", "genome": "test_reference"}
    (root / "project.yaml").write_text(yaml.safe_dump(config, sort_keys=False), encoding="utf-8")
    report = validate_project(root)
    assert report.is_valid
    return report


def _hisat2_report(root: Path):
    fastq = root / "input" / "fastq"
    fastq.mkdir(parents=True)
    for name in ("C1_R1.fastq.gz", "C1_R2.fastq.gz", "T1_R1.fastq.gz", "T1_R2.fastq.gz"):
        (fastq / name).write_bytes(b"x")
    (root / "metadata.csv").write_text("sample_id,condition\nC1,C\nT1,T\n", encoding="utf-8")
    (root / "contrasts.csv").write_text("contrast_id,factor,numerator,denominator\nT_vs_C,condition,T,C\n", encoding="utf-8")
    reference = root / "reference"
    (reference / "index").mkdir(parents=True)
    (reference / "genome.fa").write_text(">1\nACGT\n", encoding="utf-8")
    (reference / "genes.gtf").write_text("1\tt\texon\t1\t4\t.\t+\t.\tgene_id \"GeneA\";\n", encoding="utf-8")
    for number in range(1, 9):
        (reference / "index" / f"genome.{number}.ht2").write_bytes(b"x")
    config = deepcopy(base_config())
    config.update({
        "schema_version": "1.2", "input": {"type": "fastq", "path": "input/fastq", "layout": "paired_end", "preprocessing": "pretrimmed"},
        "upstream": {"engine": "nfcore_rnaseq", "pipeline_version": "3.26.0", "strandedness": "reverse", "quantification": {"method": "hisat2_featurecounts"}},
        "reference": {"source": "custom", "fasta": "reference/genome.fa", "gtf": "reference/genes.gtf", "hisat2_index": "reference/index"},
        "analysis": {"enrichment": []},
    })
    (root / "project.yaml").write_text(yaml.safe_dump(config, sort_keys=False), encoding="utf-8")
    report = validate_project(root)
    assert report.is_valid and report.execution_ready
    return report


def test_salmon_nfcore_uses_the_conda_profile(monkeypatch, tmp_path):
    report = _salmon_report(tmp_path / "project")
    _host(monkeypatch, *LINUX)
    conda_config = tmp_path / "conda.config"
    conda_config.write_text(render_upstream_conda_config(tmp_path / "cache"), encoding="utf-8")
    command = build_nextflow_command(report, samplesheet=tmp_path / "s.csv", output_dir=tmp_path / "out", profile="local", conda_config_file=conda_config)
    assert command[command.index("-profile") + 1] == "conda"
    assert str(conda_config.resolve()) in command
    assert command[command.index("-r") + 1] == "3.26.0"
    assert "docker" not in command


def test_hisat2_uses_the_conda_profile(monkeypatch, tmp_path):
    report = _hisat2_report(tmp_path / "project")
    _host(monkeypatch, *LINUX)
    conda_config = tmp_path / "conda.config"
    conda_config.write_text(render_upstream_conda_config(tmp_path / "cache"), encoding="utf-8")
    command = build_hisat2_featurecounts_command(report, samplesheet=tmp_path / "s.csv", output_dir=tmp_path / "out", profile="local", conda_config_file=conda_config)
    assert command[command.index("-profile") + 1] == "conda"
    assert str(conda_config.resolve()) in command
    assert command[command.index("--strandedness") + 1] == "reverse"
    assert "docker" not in command


# ---------------------------------------------------------------- downstream + preflight


def _frozen(monkeypatch, project_factory):
    report = validate_project(project_factory())
    generate_plan(report)
    _host(monkeypatch, *LINUX)
    _guard_conda(monkeypatch)
    run = create_case_run(report, "CASE-DISPATCH")
    freeze_case_inputs(report, run, profile="local", command=["rnaseq", "run"])
    return report, run


def test_linux_downstream_freezes_the_conda_prefix_and_local_profile(monkeypatch, project_factory, mocked_downstream_runtime):
    _report, run = _frozen(monkeypatch, project_factory)
    contract = json.loads((run.run_dir / "frozen" / "downstream_contract.json").read_text(encoding="utf-8"))
    assert contract["execution"] == {"backend": "conda", "downstream_runtime": None}
    assert _frozen_backend(run) == BACKEND_CONDA
    config = write_downstream_runtime_config(run, mocked_downstream_runtime).read_text(encoding="utf-8")
    assert "params.downstream_runtime_prefix" in config and "first_party_image" not in config
    command = build_downstream_nextflow_command(run)
    assert command[command.index("-profile") + 1] == "local"


def test_linux_preflight_never_touches_docker(monkeypatch, project_factory):
    report = validate_project(project_factory())
    generate_plan(report)
    _host(monkeypatch, *LINUX, wsl=True)
    _guard_conda(monkeypatch)
    prepare_service_run(report, profile="local")


# ---------------------------------------------------------------- end-to-end orchestration


def _fake_nextflow(observed):
    def run(command, *, cwd, stdout_path, stderr_path):
        observed.append(command)
        stdout_path.write_text("", encoding="utf-8")
        stderr_path.write_text("", encoding="utf-8")
        outdir = Path(command[command.index("--outdir") + 1])
        if "nf-core/rnaseq" in command:
            (outdir / "salmon").mkdir(parents=True)
            (outdir / "salmon" / "salmon.merged.gene_counts.tsv").write_text("gene_id\tC1\tC2\tT1\tT2\nGeneA\t1.0\t2.0\t3.0\t4.0\n", encoding="utf-8")
            (outdir / "salmon" / "salmon.merged.tx2gene_augmented.tsv").write_text("transcript_id\tgene_id\nTx1\tGeneA\n", encoding="utf-8")
            for sample in ("C1", "C2", "T1", "T2"):
                (outdir / "salmon" / sample).mkdir()
                (outdir / "salmon" / sample / "quant.sf").write_text("Name\tLength\tEffectiveLength\tTPM\tNumReads\nTx1\t100\t80\t1\t1\n", encoding="utf-8")
            (outdir / "multiqc" / "multiqc_data").mkdir(parents=True)
            (outdir / "multiqc" / "multiqc_report.html").write_text("<html></html>", encoding="utf-8")
        else:
            (outdir / "report").mkdir(parents=True)
            (outdir / "report" / "report.html").write_text("<html></html>", encoding="utf-8")
        return 0
    return run


def test_salmon_run_dispatches_every_stage_to_conda(monkeypatch, tmp_path, mocked_downstream_runtime):
    report = _salmon_report(tmp_path / "project")
    generate_plan(report)
    monkeypatch.setenv("RNASEQ_EXECUTION_ROOT", str(tmp_path / "execution"))
    _host(monkeypatch, *LINUX)
    _guard_conda(monkeypatch)
    observed: list[list[str]] = []
    monkeypatch.setattr("rnaseq.service._run_command", _fake_nextflow(observed))
    run = execute_service_run(report, case_id="CASE-SALMON")
    upstream, downstream = observed
    runtime_config = (run.run_dir / "frozen" / "downstream.runtime.config").read_text(encoding="utf-8")
    provenance = yaml.safe_load((run.run_dir / "provenance" / "run_provenance.yaml").read_text(encoding="utf-8"))
    assert provenance["execution_backend"] == BACKEND_CONDA
    assert not any("docker" in token.lower() for token in upstream + downstream)
    assert upstream[upstream.index("-profile") + 1] == "conda"
    assert downstream[downstream.index("-profile") + 1] == "local"
    assert "params.downstream_runtime_prefix" in runtime_config and "first_party_image" not in runtime_config
    assert provenance["upstream_runtime"]["kind"] == "conda"
    assert provenance["downstream_runtime"]["kind"] == "conda"
    # Retired Docker-backend fields keep the record schema and are always null.
    assert provenance["container_runtime"] is None and provenance["upstream_container_runtime"] is None
    assert provenance["execution_image"] is None and provenance["container_image"] is None
    assert all(provenance["runtime_resources"][key] is None for key in (
        "docker_architecture", "docker_memory_bytes", "docker_version", "first_party_image_architecture"))
    # The one source commit that produced the run, as frozen with its Conda runtime.
    assert provenance["source_revision"] == provenance["downstream_runtime"]["source_revision"] == "test-revision"


def test_raw_count_run_dispatches_downstream_to_conda(monkeypatch, tmp_path, project_factory, mocked_downstream_runtime):
    report = validate_project(project_factory())
    generate_plan(report)
    monkeypatch.setenv("RNASEQ_EXECUTION_ROOT", str(tmp_path / "execution"))
    _host(monkeypatch, *LINUX)
    _guard_conda(monkeypatch)
    observed: list[list[str]] = []
    monkeypatch.setattr("rnaseq.service._run_command", _fake_nextflow(observed))
    run = execute_service_run(report, case_id="CASE-RAW")
    (downstream,) = observed
    assert not any("docker" in token.lower() for token in downstream)
    assert downstream[downstream.index("-profile") + 1] == "local"
    provenance = yaml.safe_load((run.run_dir / "provenance" / "run_provenance.yaml").read_text(encoding="utf-8"))
    assert provenance["execution_backend"] == BACKEND_CONDA
    assert provenance["downstream_runtime"] is not None


def test_frozen_backend_recognizes_legacy_docker_contracts(monkeypatch, project_factory):
    _report, run = _frozen(monkeypatch, project_factory)
    assert _frozen_backend(run) == BACKEND_CONDA
    contract_path = run.run_dir / "frozen" / "downstream_contract.json"
    legacy = json.loads(contract_path.read_text(encoding="utf-8"))
    # v1.3.0's unqualified macOS Docker dispatch named its backend and image.
    legacy["execution"] = {"backend": "docker", "downstream_runtime": None, "image": "ghcr.io/2002brian/nf-rna:1.3.0", "source_revision": "0123abc"}
    contract_path.write_text(json.dumps(legacy), encoding="utf-8")
    assert _frozen_backend(run) == "docker"
    # Pre-dual-runtime Docker contracts carried only image/source_revision.
    legacy["execution"] = {"image": "ghcr.io/2002brian/nf-rna:1.1.2", "source_revision": "0123abc"}
    contract_path.write_text(json.dumps(legacy), encoding="utf-8")
    assert _frozen_backend(run) == "docker"
    legacy["execution"] = {"downstream_runtime": {"kind": "conda"}}
    contract_path.write_text(json.dumps(legacy), encoding="utf-8")
    assert _frozen_backend(run) == BACKEND_CONDA


def _incompatible_channels(monkeypatch):
    check = RuntimeCheck("Conda channels for nf-core", "NOT FOUND", "observed channels=['conda-forge']; missing bioconda.")
    monkeypatch.setattr("rnaseq.service.check_conda_channels", lambda: check)


def test_salmon_run_preflight_refuses_incompatible_conda_channels_before_any_run(monkeypatch, tmp_path, mocked_downstream_runtime):
    report = _salmon_report(tmp_path / "project")
    generate_plan(report)
    monkeypatch.setenv("RNASEQ_EXECUTION_ROOT", str(tmp_path / "execution"))
    _host(monkeypatch, *LINUX, wsl=True)
    _guard_conda(monkeypatch)
    _incompatible_channels(monkeypatch)
    with pytest.raises(ExecutionPreflightError, match="nf-core/rnaseq -profile conda needs a compatible Conda channel configuration.*missing bioconda"):
        execute_service_run(report, case_id="CASE-CHANNELS")
    assert not (tmp_path / "project" / "runs").exists()


def test_routes_without_nf_core_do_not_require_its_channels(monkeypatch, tmp_path, project_factory, mocked_downstream_runtime):
    _host(monkeypatch, *LINUX, wsl=True)
    _guard_conda(monkeypatch)
    monkeypatch.setattr("rnaseq.service.check_conda_channels", _fail("non-nf-core route checked nf-core channels"))
    for report in (_hisat2_report(tmp_path / "hisat2"), validate_project(project_factory())):
        generate_plan(report)
        prepare_service_run(report, profile="local")


@pytest.mark.real_conda_channels
@pytest.mark.parametrize(("channels", "allowed"), [(["conda-forge", "bioconda"], True), (["bioconda", "conda-forge"], False)])
def test_salmon_run_preflight_queries_the_effective_conda_channels(monkeypatch, tmp_path, mocked_downstream_runtime, channels, allowed):
    report = _salmon_report(tmp_path / "project")
    generate_plan(report)
    _host(monkeypatch, *LINUX, wsl=True)
    _guard_conda(monkeypatch)
    queried: list[list[str]] = []

    def conda(arguments):
        queried.append(list(arguments))
        return SimpleNamespace(returncode=0, stdout=json.dumps({"channel_priority": "strict", "channels": channels}), stderr="")

    monkeypatch.setattr("rnaseq.execution.shutil.which", lambda name: f"/opt/conda/bin/{name}")
    monkeypatch.setattr("rnaseq.execution._run_capture", conda)
    if allowed:
        prepare_service_run(report, profile="local")
    else:
        with pytest.raises(ExecutionPreflightError, match="conda-forge must come before bioconda"):
            prepare_service_run(report, profile="local")
    assert queried == [["/opt/conda/bin/conda", "config", "--show", "channels", "channel_priority", "--json"]]
