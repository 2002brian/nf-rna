"""Hardware-aware local resource policy: detection, auto sizing, explicit limits, fallback."""

from __future__ import annotations

import json
import os
from pathlib import Path

import pytest
import yaml

from rnaseq.errors import ExecutionPreflightError
from rnaseq.execution import (
    RuntimeSnapshot,
    detect_local_resource_capacity,
    effective_resource_budget,
    project_execution_budget,
    resolve_project_resources,
    validate_effective_resource_budget,
    validate_local_execution_budget,
    LocalResourceCapacity,
)
from rnaseq.models import ExecutionConfig
from rnaseq.resource_policy import (
    AUTO,
    GIB,
    DetectedResources,
    ProcessLimits,
    describe_policy,
    detect_process_limits,
    render_process_tuning_config,
    resolve_policy,
    salmon_quant_tuning,
)

from test_run_status import salmon_case_project, fake_nextflow_factory


def _detected(cpus: int | None, memory_gib: float | None, available_gib: float | None = None) -> DetectedResources:
    return DetectedResources(
        usable_cpus=cpus,
        usable_memory_bytes=int(memory_gib * GIB) if memory_gib is not None else None,
        available_memory_bytes=int(available_gib * GIB) if available_gib is not None else None,
    )


# ------------------------------------------------------------------ auto sizing

def test_auto_uses_detected_capacity_minus_os_reserve():
    policy = resolve_policy(AUTO, AUTO, _detected(28, 94.3, 92))
    assert (policy.requested_cpus, policy.requested_memory_gib) == (24, 79)
    assert (policy.reserve_cpus, policy.reserve_memory_gib) == (4, 15)
    assert policy.cpu_mode == policy.memory_mode == "auto"
    assert policy.fallback is False
    recorded = policy.as_dict()
    assert recorded["selected"] == {"cpus": 24, "memory_gib": 79}
    assert recorded["detected"]["usable_cpus"] == 28


def test_auto_reserve_scales_with_machine_size():
    large = resolve_policy(AUTO, AUTO, _detected(128, 512))
    assert (large.reserve_cpus, large.reserve_memory_gib) == (16, 77)
    assert (large.requested_cpus, large.requested_memory_gib) == (112, 435)


def test_auto_warns_when_currently_available_memory_is_lower_than_selected():
    policy = resolve_policy(AUTO, AUTO, _detected(28, 94, available_gib=30))
    assert policy.requested_memory_gib == 79
    assert any("Only 30.0 GiB is currently available" in note for note in policy.notes)


# ------------------------------------------------------------------ explicit settings

def test_explicit_project_limits_are_never_replaced_by_auto_sizing():
    policy = resolve_policy(22, 20, _detected(28, 94))
    assert (policy.requested_cpus, policy.requested_memory_gib) == (22, 20)
    assert policy.cpu_mode == policy.memory_mode == "explicit"
    assert policy.reserve_cpus is None and policy.reserve_memory_gib is None
    assert describe_policy(policy) == "CPUs 22 explicit in project.yaml; memory 20 GiB explicit in project.yaml"
    assert describe_policy(resolve_policy(AUTO, AUTO, _detected(28, 94))) == "CPUs 24 auto (usable 28, OS reserve 4); memory 79 GiB auto (usable 94 GiB, OS reserve 15 GiB)"


def test_mixed_explicit_cpus_and_auto_memory():
    policy = resolve_policy(10, AUTO, _detected(28, 94))
    assert (policy.requested_cpus, policy.requested_memory_gib) == (10, 79)
    assert (policy.cpu_mode, policy.memory_mode) == ("explicit", "auto")


def test_explicit_limits_survive_resolution_against_a_larger_host(monkeypatch):
    config = ExecutionConfig(max_cpus=22, max_memory_gb=20)
    holder = type("Config", (), {"execution": config})()
    snapshot = RuntimeSnapshot("Linux", "x86_64", 28, 94 * GIB, None, None, None, None)
    resources = resolve_project_resources(holder, snapshot)
    assert (resources.effective_cpus, resources.effective_memory_gib) == (22, 20)
    assert resources.clamped is False
    assert resources.as_dict()["policy"]["cpu_mode"] == "explicit"


def test_explicit_limits_larger_than_the_machine_are_clamped_visibly_not_rewritten():
    config = type("Config", (), {"execution": ExecutionConfig(max_cpus=64, max_memory_gb=256)})()
    snapshot = RuntimeSnapshot("Linux", "x86_64", 16, 64 * GIB, None, None, None, None)
    resources = resolve_project_resources(config, snapshot)
    assert (resources.requested_cpus, resources.requested_memory_gib) == (64, 256)
    assert (resources.effective_cpus, resources.effective_memory_gib) == (16, 64)
    assert resources.clamped is True and any("clamped" in warning for warning in resources.warnings)


def test_execution_config_accepts_auto_and_defaults_to_it():
    assert ExecutionConfig().max_cpus == AUTO and ExecutionConfig().max_memory_gb == AUTO
    assert ExecutionConfig(max_cpus=12).max_cpus == 12
    with pytest.raises(ValueError):
        ExecutionConfig(max_cpus="lots")
    with pytest.raises(ValueError):
        ExecutionConfig(max_memory_gb=0)
    capacity = LocalResourceCapacity(4, 8, 8)
    validate_local_execution_budget(AUTO, AUTO, capacity)  # resolved at run time instead
    with pytest.raises(ExecutionPreflightError, match="at least 8 CPUs"):
        validate_local_execution_budget(4, AUTO, capacity)


# ------------------------------------------------------------------ small machines

def test_small_host_auto_keeps_the_contract_floor_instead_of_a_reserve():
    policy = resolve_policy(AUTO, AUTO, _detected(8, 16))
    assert (policy.requested_cpus, policy.requested_memory_gib) == (8, 12)
    assert any("8-CPU local contract floor" in note for note in policy.notes)


def test_tiny_host_auto_never_exceeds_the_machine_and_existing_preflight_explains():
    policy = resolve_policy(AUTO, AUTO, _detected(4, 8))
    assert (policy.requested_cpus, policy.requested_memory_gib) == (4, 8)
    snapshot = RuntimeSnapshot("Linux", "x86_64", 4, 8 * GIB, None, None, None, None)
    config = type("Config", (), {"execution": ExecutionConfig()})()
    resources = resolve_project_resources(config, snapshot)
    assert (resources.effective_cpus, resources.effective_memory_gib) == (4, 8)
    with pytest.raises(ExecutionPreflightError, match="at least 8 CPUs/12 GiB"):
        validate_effective_resource_budget(resources)


def test_single_cpu_small_memory_host_is_handled():
    policy = resolve_policy(AUTO, AUTO, _detected(1, 1.5))
    assert (policy.requested_cpus, policy.requested_memory_gib) == (1, 1)


# ------------------------------------------------------------------ fallback

def test_undetectable_capacity_falls_back_to_the_historical_default():
    policy = resolve_policy(AUTO, AUTO, _detected(None, None))
    assert (policy.requested_cpus, policy.requested_memory_gib) == (8, 12)
    assert policy.fallback is True
    assert "auto fallback" in describe_policy(policy)
    snapshot = RuntimeSnapshot("Linux", "x86_64", None, None, None, None, None, None)
    config = type("Config", (), {"execution": ExecutionConfig()})()
    resources = resolve_project_resources(config, snapshot)
    assert (resources.effective_cpus, resources.effective_memory_gib) == (8, 12)
    validate_effective_resource_budget(resources)  # detection failure never blocks a run


def test_policy_resolution_errors_fall_back_instead_of_blocking(monkeypatch):
    def broken():
        raise RuntimeError("probe exploded")

    monkeypatch.setattr("rnaseq.execution.detect_local_resource_capacity", broken)
    config = type("Config", (), {"execution": ExecutionConfig()})()
    snapshot = RuntimeSnapshot("Linux", "x86_64", 28, 94 * GIB, None, None, None, None)
    resources = resolve_project_resources(config, snapshot)
    assert (resources.requested_cpus, resources.requested_memory_gib) == (8, 12)
    assert resources.policy.fallback is True
    assert any("probe exploded" in note for note in resources.policy.notes)
    validate_effective_resource_budget(resources)


def test_detection_failures_never_raise(monkeypatch, tmp_path):
    monkeypatch.setattr("rnaseq.execution.platform.system", lambda: "Linux")
    monkeypatch.setattr("rnaseq.execution.os.cpu_count", lambda: None)

    def broken():
        raise OSError("procfs unavailable")

    monkeypatch.setattr("rnaseq.execution._linux_memory_bytes", broken)
    monkeypatch.setattr("rnaseq.execution.detect_process_limits", lambda: ProcessLimits())
    assert detect_local_resource_capacity() == LocalResourceCapacity(None, None, None)
    # An unreadable/empty filesystem root yields "no limits", not an error.
    assert detect_process_limits(tmp_path).cgroup_cpu_limit is None


# ------------------------------------------------------------------ cgroup / affinity

def _cgroup_root(tmp_path: Path, files: dict[str, str]) -> Path:
    for relative, content in files.items():
        path = tmp_path / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(content, encoding="utf-8")
    return tmp_path


def test_cgroup_v2_limits_follow_this_process_cgroup_and_its_ancestors(tmp_path):
    root = _cgroup_root(tmp_path, {
        "proc/self/cgroup": "0::/user.slice/job.scope\n",
        "sys/fs/cgroup/user.slice/job.scope/cpu.max": "250000 100000\n",
        "sys/fs/cgroup/user.slice/job.scope/memory.max": "max\n",
        "sys/fs/cgroup/user.slice/memory.max": str(16 * GIB) + "\n",
    })
    limits = detect_process_limits(root)
    assert limits.cgroup_cpu_limit == 2.5
    assert limits.cgroup_memory_limit_bytes == 16 * GIB
    assert limits.apply_cpus(28) == 2
    assert limits.apply_memory(94 * GIB) == 16 * GIB


def test_cgroup_v1_limits_and_unlimited_values(tmp_path):
    root = _cgroup_root(tmp_path, {
        "sys/fs/cgroup/cpu/cpu.cfs_quota_us": "400000\n",
        "sys/fs/cgroup/cpu/cpu.cfs_period_us": "100000\n",
        "sys/fs/cgroup/memory/memory.limit_in_bytes": "9223372036854771712\n",
    })
    limits = detect_process_limits(root)
    assert limits.cgroup_cpu_limit == 4.0
    assert limits.cgroup_memory_limit_bytes is None
    unlimited = _cgroup_root(tmp_path / "unlimited", {"sys/fs/cgroup/cpu/cpu.cfs_quota_us": "-1\n", "sys/fs/cgroup/cpu/cpu.cfs_period_us": "100000\n"})
    assert detect_process_limits(unlimited).cgroup_cpu_limit is None


def test_capacity_respects_affinity_and_cgroup_limits(monkeypatch):
    monkeypatch.setattr("rnaseq.execution.platform.system", lambda: "Linux")
    monkeypatch.setattr("rnaseq.execution.os.cpu_count", lambda: 28)
    monkeypatch.setattr("rnaseq.execution._linux_memory_bytes", lambda: (94 * GIB, 90 * GIB))
    monkeypatch.setattr("rnaseq.execution.detect_process_limits", lambda: ProcessLimits(affinity_cpus=12, cgroup_cpu_limit=6.5, cgroup_memory_limit_bytes=32 * GIB))
    capacity = detect_local_resource_capacity()
    assert capacity == LocalResourceCapacity(6, 32, 32)


# ------------------------------------------------------------------ Salmon scheduling request

def test_salmon_quant_request_is_sized_from_the_prebuilt_index(tmp_path):
    index = tmp_path / "index"
    index.mkdir()
    with (index / "seq.bin").open("wb") as handle:
        handle.truncate(int(14.7 * GIB))  # sparse: size without disk use
    tuning, notes = salmon_quant_tuning(index)
    assert notes == ()
    assert [(item.selector, item.memory_gib) for item in tuning] == [(".*:SALMON_QUANT", 19)]
    rendered = render_process_tuning_config(tuning)
    assert "withName: '.*:SALMON_QUANT' { memory = '19.GB' }" in rendered
    assert "cpus" not in rendered.split("process {", 1)[1]  # Salmon --threads is never changed


def test_salmon_quant_request_is_capped_at_nfcore_default_and_skipped_without_an_index(tmp_path):
    index = tmp_path / "huge"
    index.mkdir()
    with (index / "seq.bin").open("wb") as handle:
        handle.truncate(60 * GIB)
    tuning, _ = salmon_quant_tuning(index)
    assert tuning[0].memory_gib == 36
    assert salmon_quant_tuning(None)[0] == ()
    assert salmon_quant_tuning(tmp_path / "missing")[0] == ()


# ------------------------------------------------------------------ recorded in a real (mocked) run

def test_run_freezes_and_records_the_resource_policy_and_salmon_tuning(monkeypatch, tmp_path):
    report = salmon_case_project(tmp_path, monkeypatch)
    from rnaseq.resource_policy import ProcessTuning
    monkeypatch.setattr(
        "rnaseq.service.salmon_quant_tuning",
        lambda index: ((ProcessTuning(".*:SALMON_QUANT", 19, "test index"),), ()),
    )
    observed: list[list[str]] = []
    monkeypatch.setattr("rnaseq.service._run_command", fake_nextflow_factory(observed))
    from rnaseq.service import execute_service_run

    run = execute_service_run(report, case_id="CASE-POLICY")
    tuning = run.run_dir / "frozen" / "nfcore.tuning.config"
    assert "memory = '19.GB'" in tuning.read_text(encoding="utf-8")
    upstream, downstream = observed
    assert str(tuning.resolve()) in upstream and str(tuning.resolve()) not in downstream
    local = (run.run_dir / "frozen" / "nfcore.local.config").read_text(encoding="utf-8")
    assert "executor { cpus = 14; memory = '54.GB' }" in local  # conftest host 16 CPUs / 64 GiB, auto
    provenance = yaml.safe_load((run.run_dir / "provenance" / "run_provenance.yaml").read_text(encoding="utf-8"))
    policy = provenance["runtime_resources"]["policy"]
    assert policy["selected"] == {"cpus": 14, "memory_gib": 54}
    assert policy["process_tuning"][0]["memory_gib"] == 19
    assert provenance["frozen_nfcore_tuning_config"]["path"] == "frozen/nfcore.tuning.config"
    manifest = yaml.safe_load((run.run_dir / "frozen" / "execution_manifest.yaml").read_text(encoding="utf-8"))
    assert manifest["resource_policy"]["effective"] == {"cpus": 14, "memory_gib": 54}
    log = (run.run_dir / "logs" / "rnaseq.log").read_text(encoding="utf-8")
    assert "resources: 14 CPUs / 54 GiB effective" in log


def test_run_with_explicit_limits_freezes_them_unchanged(monkeypatch, tmp_path):
    report = salmon_case_project(tmp_path, monkeypatch, execution={"profile": "local", "max_cpus": 10, "max_memory_gb": 20})
    observed: list[list[str]] = []
    monkeypatch.setattr("rnaseq.service._run_command", fake_nextflow_factory(observed))
    from rnaseq.service import execute_service_run

    run = execute_service_run(report, case_id="CASE-EXPLICIT")
    local = (run.run_dir / "frozen" / "nfcore.local.config").read_text(encoding="utf-8")
    assert "executor { cpus = 10; memory = '20.GB' }" in local
    provenance = yaml.safe_load((run.run_dir / "provenance" / "run_provenance.yaml").read_text(encoding="utf-8"))
    assert provenance["runtime_resources"]["policy"]["cpu_mode"] == "explicit"
    assert provenance["runtime_resources"]["effective"] == {"cpus": 10, "memory_gib": 20}
    # No prebuilt index for this iGenomes fixture: nf-core's default request is kept and said so.
    assert not (run.run_dir / "frozen" / "nfcore.tuning.config").exists()
    assert any("keeps nf-core's default request" in note for note in provenance["runtime_resources"]["policy"]["notes"])
