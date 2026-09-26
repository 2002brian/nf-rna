"""Native linux-64 ``rnaseq doctor``: Nextflow + Conda prerequisites, never Docker."""

from __future__ import annotations

import hashlib
import json
import os
import platform
import shutil
from pathlib import Path
from types import SimpleNamespace

import pytest

from rnaseq import downstream_runtime
from rnaseq.execution import (
    HISAT2_LINUX_CONDA_ENV,
    HISAT2_LINUX_CONDA_ENV_SHA256,
    _writable_location_check,
    check_conda_channels,
    check_conda_functional,
    check_java,
    check_nextflow_suitability,
    doctor_checks,
    doctor_readiness,
    downstream_runtime_checks,
    hisat2_conda_env_check,
)


ROOT = Path(__file__).parents[1]


def _completed(stdout: str = "", stderr: str = "", returncode: int = 0) -> SimpleNamespace:
    return SimpleNamespace(returncode=returncode, stdout=stdout, stderr=stderr)


def _native_host(
    monkeypatch, tmp_path: Path, *, java: str = '17.0.20.1', nextflow: str = "26.04.6", conda_platform: str = "linux-64",
    channels: object = ("conda-forge", "bioconda"),
) -> list[list[str]]:
    """A qualified linux-64 host whose Docker is absent and must never be queried."""

    calls: list[list[str]] = []

    def run(arguments):
        calls.append(list(arguments))
        executable = Path(arguments[0]).name
        if executable == "docker":
            raise AssertionError(f"Linux doctor invoked Docker: {arguments}")
        if executable == "java":
            return _completed(stderr=f'openjdk version "{java}" 2026-08-18\nOpenJDK Runtime Environment\n')
        if executable == "nextflow":
            return _completed(stdout=f"      version {nextflow} build 12646\n")
        if executable == "conda" and arguments[1:3] == ["config", "--show"]:
            if isinstance(channels, SimpleNamespace):
                return channels
            return _completed(stdout=json.dumps({"channel_priority": "strict", "channels": list(channels)}))
        if executable == "conda":
            return _completed(stdout=json.dumps({"platform": conda_platform, "conda_version": "26.5.3"}))
        raise AssertionError(f"unexpected doctor subprocess: {arguments}")

    def no_docker(*_args, **_kwargs):
        raise AssertionError("Linux doctor used a Docker-oriented check")

    monkeypatch.setattr(platform, "system", lambda: "Linux")
    monkeypatch.setattr(platform, "machine", lambda: "x86_64")
    monkeypatch.setattr("rnaseq.execution._run_capture", run)
    monkeypatch.setattr(shutil, "which", lambda name: None if name == "docker" else f"/opt/native/bin/{name}")
    for guarded in ("check_docker", "check_container_runtime", "inspect_container_image", "runtime_snapshot", "downstream_docker_user_mapping_check"):
        monkeypatch.setattr(f"rnaseq.execution.{guarded}", no_docker)
    monkeypatch.setattr("rnaseq.downstream.r_runtime_checks", no_docker)
    monkeypatch.setattr("importlib.util.find_spec", lambda name, *args: object() if name == "build" else None)
    for variable in ("NXF_JAVA_HOME", "JAVA_HOME"):
        monkeypatch.delenv(variable, raising=False)
    monkeypatch.setenv("RNASEQ_EXECUTION_ROOT", str(tmp_path / "execution"))
    monkeypatch.setenv("RNASEQ_RUNTIME_ROOT", str(tmp_path / "runtime"))
    monkeypatch.setenv("NXF_HOME", str(tmp_path / "nxf-home"))
    return calls


def test_linux_doctor_never_calls_docker_and_passes_without_it(monkeypatch, tmp_path):
    calls = _native_host(monkeypatch, tmp_path)
    checks = doctor_checks()
    by_name = {check.name: check for check in checks}

    assert not any(call[0].endswith("docker") for call in calls)
    assert not any("docker" in check.name.lower() or "image" in check.name.lower() for check in checks)
    assert [check for check in checks if check.verdict == "FAIL"] == []
    for required in ("Runtime platform", "Execution backend", "Java", "Nextflow", "Conda", "nf-core/rnaseq pin", "Execution root",
                     "Upstream Conda cache", "HISAT2/featureCounts Conda environment", "Downstream Conda lock",
                     "Downstream runtime location", "Downstream nf-rna wheel"):
        assert by_name[required].verdict == "PASS" or (required == "nf-core/rnaseq pin" and by_name[required].verdict == "WARN")
    assert "container-runtime ceiling applies=False" in by_name["Effective local budget"].detail
    # Doctor is non-mutating: no execution, cache, runtime, or Nextflow directories were created.
    assert list(tmp_path.iterdir()) == []


def test_linux_doctor_does_not_require_host_bioinformatics_or_r_tools(monkeypatch, tmp_path):
    calls = _native_host(monkeypatch, tmp_path)
    doctor_checks()
    executables = {Path(call[0]).name for call in calls}
    assert executables <= {"java", "nextflow", "conda"}


def test_linux_doctor_names_the_platform_and_conda_backend(monkeypatch, tmp_path):
    _native_host(monkeypatch, tmp_path)
    by_name = {check.name: check for check in doctor_checks()}
    assert by_name["Runtime platform"].verdict == "PASS"
    assert by_name["Runtime platform"].detail.startswith("linux-64")
    assert by_name["Execution backend"].detail == "Nextflow + Conda"


def test_macos_doctor_selects_docker_checks(monkeypatch):
    monkeypatch.setattr(platform, "system", lambda: "Darwin")
    monkeypatch.setattr(platform, "machine", lambda: "arm64")
    monkeypatch.setattr("rnaseq.execution.container_doctor_checks", lambda _project=None: ("docker-checks",))
    monkeypatch.setattr("rnaseq.execution.native_linux_doctor_checks", lambda *_a: pytest.fail("macOS ran native Conda doctor"))
    platform_check, backend, *rest = doctor_checks()
    assert (platform_check.name, platform_check.detail) == ("Runtime platform", "darwin-arm64")
    assert (backend.name, backend.detail) == ("Execution backend", "Docker")
    assert rest == ["docker-checks"]


def test_native_windows_doctor_directs_to_wsl2(monkeypatch):
    monkeypatch.setattr(platform, "system", lambda: "Windows")
    monkeypatch.setattr(platform, "machine", lambda: "AMD64")
    (check,) = doctor_checks()
    assert check.verdict == "FAIL"
    assert "WSL2" in check.detail and "Native Windows is not a supported" in check.detail


@pytest.mark.parametrize(("conda_platform", "verdict"), [("linux-64", "PASS"), ("osx-arm64", "FAIL")])
def test_conda_must_be_functional_for_the_native_platform(monkeypatch, tmp_path, conda_platform, verdict):
    _native_host(monkeypatch, tmp_path, conda_platform=conda_platform)
    assert check_conda_functional().verdict == verdict


def test_missing_or_broken_conda_fails_with_an_actionable_message(monkeypatch, tmp_path):
    _native_host(monkeypatch, tmp_path)
    monkeypatch.setattr("rnaseq.execution._run_capture", lambda _arguments: _completed(stderr="boom", returncode=1))
    assert "not functional" in check_conda_functional().detail
    monkeypatch.setattr(shutil, "which", lambda _name: None)
    missing = check_conda_functional()
    assert missing.verdict == "FAIL"
    assert "Conda executable was not found on PATH" in missing.detail


@pytest.mark.parametrize(("version", "verdict"), [("17.0.20.1", "PASS"), ("21.0.4", "PASS"), ("11.0.2", "FAIL"), ("1.8.0_392", "FAIL")])
def test_java_must_be_17_or_later(monkeypatch, tmp_path, version, verdict):
    _native_host(monkeypatch, tmp_path, java=version)
    check = check_java()
    assert check.verdict == verdict
    if verdict == "FAIL":
        assert "requires Java 17 or later" in check.detail


def test_java_resolution_follows_the_nextflow_launcher(monkeypatch, tmp_path):
    calls = _native_host(monkeypatch, tmp_path)
    monkeypatch.setenv("JAVA_HOME", "/opt/jdk")
    monkeypatch.setenv("NXF_JAVA_HOME", "/opt/nxf-jdk")
    assert "from NXF_JAVA_HOME" in check_java().detail
    assert calls[-1][0] == "/opt/nxf-jdk/bin/java"
    monkeypatch.setattr(shutil, "which", lambda _name: None)
    monkeypatch.delenv("NXF_JAVA_HOME")
    monkeypatch.delenv("JAVA_HOME")
    assert check_java().verdict == "FAIL"


@pytest.mark.parametrize(("version", "verdict"), [("26.04.6", "PASS"), ("25.04.3", "PASS"), ("25.04.2", "FAIL"), ("24.10.0", "FAIL")])
def test_nextflow_must_satisfy_the_pinned_nfcore_minimum(monkeypatch, tmp_path, version, verdict):
    _native_host(monkeypatch, tmp_path, nextflow=version)
    check = check_nextflow_suitability()
    assert check.verdict == verdict
    if verdict == "FAIL":
        assert ">=25.04.3" in check.detail


def test_missing_nextflow_fails(monkeypatch, tmp_path):
    _native_host(monkeypatch, tmp_path)

    def missing(arguments):
        raise FileNotFoundError(arguments[0])

    monkeypatch.setattr("rnaseq.execution._run_capture", missing)
    check = check_nextflow_suitability()
    assert check.verdict == "FAIL"
    assert "not found on PATH" in check.detail


def test_downstream_linux_lock_identity_is_verified(monkeypatch, tmp_path):
    _native_host(monkeypatch, tmp_path)
    lock = downstream_runtime.lock_path_for_platform("linux-64")
    passed = downstream_runtime_checks()[0]
    assert passed.verdict == "PASS"
    assert "18197d49b7087d39c3a8e3aef8a16d7bd65af7e56e07d87aa47be656b219617d" in passed.detail

    tampered = tmp_path / "locks" / lock.name
    tampered.parent.mkdir()
    tampered.write_bytes(lock.read_bytes() + b"# drift\n")
    monkeypatch.setattr(downstream_runtime, "lock_path_for_platform", lambda *_args: tampered)
    failed = downstream_runtime_checks()
    assert failed[0].verdict == "FAIL"
    assert "checksum mismatch" in failed[0].detail


def test_hisat2_environment_identity_is_the_reviewed_file():
    assert hashlib.sha256(HISAT2_LINUX_CONDA_ENV.read_bytes()).hexdigest() == HISAT2_LINUX_CONDA_ENV_SHA256
    check = hisat2_conda_env_check()
    assert check.verdict == "PASS"
    assert "hisat2=2.2.3 samtools=1.21 subread=2.0.6" in check.detail


def test_hisat2_environment_drift_fails(monkeypatch, tmp_path):
    changed = tmp_path / HISAT2_LINUX_CONDA_ENV.name
    changed.write_text(HISAT2_LINUX_CONDA_ENV.read_text(encoding="utf-8") + "  - conda-forge::extra=1.0=0\n", encoding="utf-8")
    monkeypatch.setattr("rnaseq.execution.HISAT2_LINUX_CONDA_ENV", changed)
    assert "sha256 mismatch" in hisat2_conda_env_check().detail

    def reviewed(text: str) -> None:
        changed.write_text(text, encoding="utf-8")
        monkeypatch.setattr("rnaseq.execution.HISAT2_LINUX_CONDA_ENV_SHA256", hashlib.sha256(changed.read_bytes()).hexdigest())

    original = HISAT2_LINUX_CONDA_ENV.read_text(encoding="utf-8")
    reviewed(original.replace("bioconda::hisat2=2.2.3=h8471819_0", "bioconda::hisat2"))
    assert "Not exactly pinned" in hisat2_conda_env_check().detail
    reviewed(original.replace("bioconda::subread=2.0.6=he4a0461_2", "bioconda::subread=2.0.8=h577a1d6_0"))
    drift = hisat2_conda_env_check()
    assert drift.verdict == "FAIL"
    assert "subread=2.0.8 (expected 2.0.6)" in drift.detail


def test_writable_location_check_is_non_mutating(tmp_path):
    target = tmp_path / "cache" / "nested"
    check = _writable_location_check("Cache", target, "test")
    assert check.verdict == "PASS"
    assert f"will be created under {tmp_path}" in check.detail
    assert not (tmp_path / "cache").exists()


@pytest.mark.skipif(hasattr(os, "geteuid") and os.geteuid() == 0, reason="root bypasses directory permissions")
def test_unwritable_location_fails(tmp_path):
    locked = tmp_path / "locked"
    locked.mkdir()
    locked.chmod(0o555)
    try:
        check = _writable_location_check("Cache", locked / "cache", "test")
    finally:
        locked.chmod(0o755)
    assert check.verdict == "FAIL"
    assert "not a writable directory" in check.detail


def test_downstream_prefix_without_identity_marker_fails(monkeypatch, tmp_path):
    _native_host(monkeypatch, tmp_path)
    sha = downstream_runtime._lock_checksum(downstream_runtime.lock_path_for_platform("linux-64"))
    (tmp_path / "runtime" / "prefixes" / f"nf-rna-downstream-linux-64-{sha[:16]}").mkdir(parents=True)
    location = {check.name: check for check in downstream_runtime_checks()}["Downstream runtime location"]
    assert location.verdict == "FAIL"
    assert "without its wheel identity marker" in location.detail


def test_source_checkout_without_build_tool_fails_before_run(monkeypatch, tmp_path):
    _native_host(monkeypatch, tmp_path)
    monkeypatch.setattr("importlib.util.find_spec", lambda _name, *args: None)
    wheel = {check.name: check for check in downstream_runtime_checks()}["Downstream nf-rna wheel"]
    assert wheel.verdict == "FAIL"
    assert "not importable" in wheel.detail and "install nf-rna normally" in wheel.detail


def test_installed_distribution_needs_no_build_tool(monkeypatch, tmp_path):
    _native_host(monkeypatch, tmp_path)
    monkeypatch.setattr("importlib.util.find_spec", lambda _name, *args: None)
    monkeypatch.setattr(downstream_runtime, "_source_checkout", lambda: None)
    monkeypatch.setattr(downstream_runtime, "_installed_distribution", lambda: object())
    monkeypatch.setattr(downstream_runtime, "_installed_wheel_members", lambda _dist: ("nf_rna-1.2.0.dist-info", {"a": b""}))
    wheel = {check.name: check for check in downstream_runtime_checks()}["Downstream nf-rna wheel"]
    assert wheel.verdict == "PASS"
    assert "without build tools" in wheel.detail


def test_unowned_installation_fails_with_install_guidance(monkeypatch, tmp_path):
    _native_host(monkeypatch, tmp_path)
    monkeypatch.setattr(downstream_runtime, "_source_checkout", lambda: None)
    monkeypatch.setattr(downstream_runtime.metadata, "distributions", lambda **_kwargs: [])
    wheel = {check.name: check for check in downstream_runtime_checks()}["Downstream nf-rna wheel"]
    assert wheel.verdict == "FAIL"
    assert "install nf-rna with pip" in wheel.detail


# ---------------------------------------------------------------- nf-core Conda channels


@pytest.mark.parametrize(
    ("channels", "verdict"),
    [
        (["conda-forge", "bioconda"], "PASS"),
        (["conda-forge"], "FAIL"),
        (["bioconda"], "FAIL"),
        (["bioconda", "conda-forge"], "FAIL"),
        ([], "FAIL"),
        (["defaults", "conda-forge", "my-lab", "bioconda", "r"], "PASS"),
        (["conda-forge", "defaults", "bioconda"], "PASS"),
    ],
    ids=["exact", "missing-bioconda", "missing-conda-forge", "wrong-order", "empty", "extra-channels", "extra-between"],
)
def test_nfcore_conda_channels_require_conda_forge_before_bioconda(monkeypatch, tmp_path, channels, verdict):
    calls = _native_host(monkeypatch, tmp_path, channels=channels)
    check = check_conda_channels()
    assert check.verdict == verdict, check.detail
    # The effective, merged configuration is queried from Conda itself, as JSON.
    assert calls[-1][1:] == ["config", "--show", "channels", "channel_priority", "--json"]
    assert f"observed channels={channels}" in check.detail
    assert "conda-forge before bioconda" in check.detail
    if verdict == "FAIL":
        assert "conda config --add channels bioconda && conda config --add channels conda-forge" in check.detail


@pytest.mark.parametrize(
    "response",
    [
        SimpleNamespace(returncode=0, stdout="channels:\n  - conda-forge\n", stderr=""),
        SimpleNamespace(returncode=1, stdout="", stderr="CondaError: bad .condarc"),
        SimpleNamespace(returncode=0, stdout=json.dumps({"channels": "conda-forge,bioconda"}), stderr=""),
        SimpleNamespace(returncode=0, stdout=json.dumps({"channel_priority": "strict"}), stderr=""),
    ],
    ids=["not-json", "conda-error", "channels-not-a-list", "channels-absent"],
)
def test_unreadable_conda_channel_configuration_is_not_ready(monkeypatch, tmp_path, response):
    _native_host(monkeypatch, tmp_path, channels=response)
    check = check_conda_channels()
    assert check.verdict == "FAIL"
    assert "Could not read the effective channel list" in check.detail


def test_doctor_is_not_ready_when_bioconda_is_missing(monkeypatch, tmp_path):
    _native_host(monkeypatch, tmp_path, channels=["conda-forge"])
    checks = doctor_checks()
    by_name = {check.name: check for check in checks}
    assert by_name["Conda channels for nf-core"].verdict == "FAIL"
    assert doctor_readiness(checks) == (False, ("Conda channels for nf-core",))


def test_doctor_is_ready_with_nfcore_compatible_channels(monkeypatch, tmp_path):
    _native_host(monkeypatch, tmp_path)
    assert doctor_readiness(doctor_checks()) == (True, ())


def test_doctor_cli_reports_overall_readiness_and_exit_code(monkeypatch, tmp_path):
    from typer.testing import CliRunner

    from rnaseq.cli import app

    _native_host(monkeypatch, tmp_path, channels=["conda-forge"])
    result = CliRunner().invoke(app, ["doctor"])
    assert result.exit_code == 1
    assert "Conda channels for nf-core: FAIL" in result.output
    assert "Overall: NOT READY — failed: Conda channels for nf-core" in result.output

    _native_host(monkeypatch, tmp_path, channels=["conda-forge", "bioconda"])
    result = CliRunner().invoke(app, ["doctor"])
    assert result.exit_code == 0, result.output
    assert "Overall: READY" in result.output


def test_channel_check_only_warns_for_a_project_that_does_not_use_nf_core(monkeypatch, tmp_path, project_factory):
    raw_counts_project = project_factory()
    _native_host(monkeypatch, tmp_path, channels=["conda-forge"])
    by_name = {check.name: check for check in doctor_checks(raw_counts_project)}
    assert by_name["Conda channels for nf-core"].verdict == "WARN"
    assert "not used by this project's route" in by_name["Conda channels for nf-core"].detail
