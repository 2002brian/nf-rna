"""Minimal, safe M2 launcher for one local nf-core/rnaseq execution route."""

from __future__ import annotations

import csv
import hashlib
import json
import math
import os
import platform
import re
import shutil
import subprocess
import sys
import tempfile
from dataclasses import dataclass
from datetime import UTC, datetime
from io import StringIO
from pathlib import Path
from typing import Any

import yaml

from rnaseq.errors import ExecutionPreflightError, UpstreamExecutionError
from rnaseq.downstream_runtime import native_platform
from rnaseq.models import FastqPreprocessing, InputType, NFCORE_RNASEQ_VERSION, PIPELINE_VERSION, ReferenceConfig
from rnaseq.hisat2_featurecounts import COUNTING_POLICY, HISAT2_VERSION, SAMTOOLS_VERSION, SUBREAD_VERSION
from rnaseq.planner import render_manifest, render_samplesheet
from rnaseq.references import LocalReferenceError, load_local_reference, sha256_file
from rnaseq.resource_policy import AUTO, GIB, DetectedResources, ResourcePolicy, describe_policy, detect_process_limits, resolve_policy
from rnaseq.validators import ValidationReport
from rnaseq.workflow_assets import workflow_asset_path

LOCAL_PROFILE = "local"
NFCORE_CONDA_PROFILE = "conda"
# Production runtime policy (v1.3.0+): Linux x86-64 including WSL2 with
# Nextflow + Conda is the only backend.  v1.2.1 was the last Docker release;
# the unqualified macOS Docker dispatch still present in v1.3.0 was removed.
BACKEND_CONDA = "conda"
BACKEND_LABELS = {BACKEND_CONDA: "Nextflow + Conda"}
WINDOWS_UNSUPPORTED_MESSAGE = (
    "Native Windows is not a supported nf-rna runtime. Install WSL2 with an Ubuntu distribution, "
    "then install and run nf-rna inside WSL2 (Linux x86-64, Nextflow + Conda)."
)
NFCORE_RNASEQ_REVISION = "e7ca46272c8f9d5ceee3f71759f4ba551d3217a4"
RUN_STATES = {"CREATED", "RUNNING", "SUCCESS", "FAILED", "INTERRUPTED"}
EXECUTION_ROOT_ENV = "RNASEQ_EXECUTION_ROOT"
WORK_ROOT_ENV = "RNASEQ_WORK_ROOT"
HISAT2_WORKFLOW = workflow_asset_path("hisat2_featurecounts.nf")
HISAT2_LINUX_CONDA_ENV = HISAT2_WORKFLOW.parent / "envs" / "hisat2-featurecounts-linux-64.yml"
# Reviewed identity of the qualified linux-64 HISAT2/featureCounts environment.
HISAT2_LINUX_CONDA_ENV_SHA256 = "2c407fb2b37529db8d1b70b7c757140124868e3085b2c83a7f57100249e0aa0c"
# nf-core/rnaseq 3.26.0 manifest: nextflowVersion = '!>=25.04.3'; Nextflow 25+ needs Java 17+.
NEXTFLOW_MINIMUM_VERSION = "25.04.3"
JAVA_MINIMUM_MAJOR = 17


@dataclass(frozen=True)
class RuntimeCheck:
    name: str
    state: str
    detail: str
    level: str | None = None

    @property
    def verdict(self) -> str:
        if self.level is not None:
            return self.level
        return "PASS" if self.state == "FOUND" else "FAIL"


@dataclass(frozen=True)
class ResourceContract:
    """Explicit, portable task resource declaration expressed in GiB."""

    name: str
    cpus: int
    memory_gib: int
    time_hours: int
    max_forks: int | None = None


RESOURCE_CONTRACTS = {
    "SMALL": ResourceContract("SMALL", cpus=1, memory_gib=2, time_hours=2),
    "MEDIUM": ResourceContract("MEDIUM", cpus=4, memory_gib=8, time_hours=8),
    "LARGE": ResourceContract("LARGE", cpus=8, memory_gib=12, time_hours=12),
}
LOCAL_RESOURCE_CEILING = RESOURCE_CONTRACTS["LARGE"]


@dataclass(frozen=True)
class LocalResourceCapacity:
    logical_cpus: int | None
    total_memory_gib: int | None
    available_memory_gib: int | None


def detect_local_resource_capacity() -> LocalResourceCapacity:
    """Read the capacity this process may use; Linux/WSL and macOS are supported.

    On Linux the CPU count respects CPU affinity and cgroup v1/v2 CPU quotas,
    and memory respects a cgroup memory limit.  Any failure leaves a value
    unavailable rather than raising.
    """

    cpus = os.cpu_count()
    total_bytes: int | None = None
    available_bytes: int | None = None
    system = platform.system().lower()
    try:
        if system == "darwin":
            total_bytes = _darwin_memory_bytes()
            available_bytes = total_bytes
        elif system == "linux":
            total_bytes, available_bytes = _linux_memory_bytes()
    except (OSError, ValueError, KeyError):
        pass
    cpus = cpus if cpus and cpus > 0 else None
    if system == "linux":
        limits = detect_process_limits()
        cpus = limits.apply_cpus(cpus)
        total_bytes = limits.apply_memory(total_bytes)
        if available_bytes is not None and total_bytes is not None:
            available_bytes = min(available_bytes, total_bytes)
    total = total_bytes // GIB if total_bytes is not None else None
    available = available_bytes // GIB if available_bytes is not None else total
    return LocalResourceCapacity(cpus, total, available)


def suggested_local_resources(capacity: LocalResourceCapacity) -> tuple[int, int]:
    """Reserve roughly one fifth of CPUs and memory, rounded conservatively."""

    fallback = (LOCAL_RESOURCE_CEILING.cpus, LOCAL_RESOURCE_CEILING.memory_gib)
    if capacity.logical_cpus is None or capacity.total_memory_gib is None:
        return fallback
    cpus = max(1, math.floor(capacity.logical_cpus * 0.8))
    usable_memory = min(capacity.total_memory_gib, capacity.available_memory_gib or capacity.total_memory_gib)
    memory = max(1, (math.floor(usable_memory * 0.8) // 4) * 4)
    return cpus, memory


def validate_local_execution_budget(cpus: int | str, memory_gb: int | str, capacity: LocalResourceCapacity) -> None:
    """Check explicit project limits; ``auto`` is resolved and checked at runtime."""

    for value in (cpus, memory_gb):
        if value == AUTO:
            continue
        if isinstance(value, bool) or not isinstance(value, int) or value < 1:
            raise ExecutionPreflightError("Execution CPU and memory limits must be positive integers or 'auto'.")
    if (cpus != AUTO and cpus < LOCAL_RESOURCE_CEILING.cpus) or (memory_gb != AUTO and memory_gb < LOCAL_RESOURCE_CEILING.memory_gib):
        raise ExecutionPreflightError("Execution budget must be at least 8 CPUs and 12 GiB to satisfy enabled local process contracts.")
    # A project budget is a user-selected upper bound, not a claim about this
    # machine.  Runtime preflight computes and records a visible effective
    # budget instead of rejecting portable project configuration here.


def project_resource_policy(config: Any, snapshot: "RuntimeSnapshot | None" = None) -> ResourcePolicy:
    """Resolve the project's explicit or ``auto`` limits against this machine.

    The runtime snapshot is the capacity authority (affinity/cgroup aware on
    Linux); raw observations are recorded alongside it for provenance.
    """

    try:
        return _project_resource_policy(config, snapshot)
    except Exception as exc:  # resource tuning must never prevent a run
        policy = resolve_policy(config.execution.max_cpus, config.execution.max_memory_gb, DetectedResources())
        return policy.with_process_tuning((), (f"Resource detection failed ({exc}); used the safe fallback.",))


def _project_resource_policy(config: Any, snapshot: "RuntimeSnapshot | None") -> ResourcePolicy:
    snapshot = snapshot or native_runtime_snapshot()
    capacity = detect_local_resource_capacity()
    details: dict[str, object] = {"logical_cpus": os.cpu_count(), "host_memory_gib": (
        round(snapshot.host_memory_bytes / GIB, 1) if snapshot.host_memory_bytes else None
    )}
    if platform.system().lower() == "linux":
        details.update(detect_process_limits().as_dict())
    warnings = tuple(
        message for condition, message in (
            (snapshot.logical_cpus is None, "Usable CPU count is unavailable."),
            (snapshot.host_memory_bytes is None, "Usable memory is unavailable."),
        ) if condition
    )
    detected = DetectedResources(
        usable_cpus=snapshot.logical_cpus,
        usable_memory_bytes=snapshot.host_memory_bytes,
        available_memory_bytes=capacity.available_memory_gib * GIB if capacity.available_memory_gib is not None else None,
        details=details,
        warnings=warnings,
    )
    return resolve_policy(config.execution.max_cpus, config.execution.max_memory_gb, detected)


def project_execution_budget(config: Any, policy: ResourcePolicy | None = None) -> ResourceContract:
    policy = policy or project_resource_policy(config)
    return ResourceContract("PROJECT_LOCAL", policy.requested_cpus, policy.requested_memory_gib, LOCAL_RESOURCE_CEILING.time_hours)


def resolve_project_resources(config: Any, snapshot: "RuntimeSnapshot") -> "EffectiveResourceBudget":
    """Policy plus the effective budget for this machine, resolved once per run."""

    policy = project_resource_policy(config, snapshot)
    return effective_resource_budget(snapshot, project_execution_budget(config, policy), policy=policy)


@dataclass(frozen=True)
class RuntimeSnapshot:
    host_os: str
    host_architecture: str
    logical_cpus: int | None
    host_memory_bytes: int | None


@dataclass(frozen=True)
class EffectiveResourceBudget:
    requested_cpus: int
    requested_memory_gib: int
    host_cpus: int | None
    host_memory_gib: int | None
    effective_cpus: int
    effective_memory_gib: int
    clamped: bool
    warnings: tuple[str, ...] = ()
    policy: ResourcePolicy | None = None

    def as_dict(self) -> dict[str, object]:
        payload: dict[str, object] = {
            "requested": {"cpus": self.requested_cpus, "memory_gib": self.requested_memory_gib},
            "host": {"cpus": self.host_cpus, "memory_gib": self.host_memory_gib},
            # Record-schema compatibility: the retired Docker backend's VM
            # ceiling; the native Conda runtime has none.
            "container_runtime": {"cpus": None, "memory_gib": None, "ceiling_applies": False},
            "effective": {"cpus": self.effective_cpus, "memory_gib": self.effective_memory_gib},
            "clamped": self.clamped,
            "warnings": list(self.warnings),
        }
        if self.policy is not None:
            payload["policy"] = self.policy.as_dict()
        return payload


@dataclass(frozen=True)
class PreparedRun:
    report: ValidationReport
    nextflow_version: str
    container_runtime: str
    resource_budget: EffectiveResourceBudget | None = None


@dataclass(frozen=True)
class RunResult:
    run_id: str
    run_dir: Path
    state_path: Path
    handoff_path: Path


@dataclass(frozen=True)
class ExecutionWorkspace:
    """Local, operational Nextflow paths kept separate from immutable results."""

    root: Path
    launch_dir: Path
    work_dir: Path


def _execution_base() -> Path:
    configured = os.environ.get(EXECUTION_ROOT_ENV)
    if configured:
        base = Path(configured).expanduser()
    elif sys.platform == "darwin":
        base = Path.home() / "Library" / "Caches" / "nf-rna"
    else:
        base = Path(os.environ.get("XDG_CACHE_HOME", str(Path.home() / ".cache"))) / "nf-rna"
    if not base.is_absolute():
        raise ExecutionPreflightError(f"{EXECUTION_ROOT_ENV} must be an absolute path when configured.")
    return base


def _work_base() -> Path | None:
    configured = os.environ.get(WORK_ROOT_ENV)
    if not configured:
        return None
    try:
        base = Path(configured).expanduser()
    except RuntimeError:  # an unknown ~user
        base = Path(configured)
    if not base.is_absolute():
        raise ExecutionPreflightError(f"{WORK_ROOT_ENV} must be an absolute path when configured: {configured!r}")
    return base


def resolve_execution_workspace(case_id: str, run_id: str) -> ExecutionWorkspace:
    """Return the local execution paths for one logical case/run without writing them.

    Nextflow creates ``.nextflow/cache`` below its current working directory.  That
    LevelDB state is operational scratch data, not a durable project artifact, so
    it must never be created on a client or shared project filesystem.  Only the
    task work directory, which holds the large transient files, may be relocated
    with ``RNASEQ_WORK_ROOT``; the launch directory stays under the execution root.
    """

    root = (_execution_base() / case_id / run_id).resolve()
    work_base = _work_base()
    work_dir = (work_base / case_id / run_id / "work").resolve() if work_base is not None else root / "work"
    return ExecutionWorkspace(root=root, launch_dir=root / "launch", work_dir=work_dir)


def execution_backend() -> str:
    """Select the production backend: only linux-64 (including WSL2) with Nextflow + Conda is supported."""

    system = platform.system().lower()
    machine = platform.machine().lower()
    if system == "linux" and machine in {"x86_64", "amd64"}:
        return BACKEND_CONDA
    if system == "windows" or system.startswith(("cygwin", "msys", "mingw")):
        raise ExecutionPreflightError(WINDOWS_UNSUPPORTED_MESSAGE)
    raise ExecutionPreflightError(
        f"Unsupported nf-rna runtime platform {system}-{machine}. Supported: Linux x86-64 including "
        "WSL2 (Nextflow + Conda)."
    )


def runtime_platform_label() -> str:
    system, machine = platform.system().lower(), platform.machine().lower()
    if system == "linux":
        return "linux-64 / WSL2" if _is_wsl() else "linux-64"
    if system == "darwin":
        return f"darwin-{'arm64' if machine in {'arm64', 'aarch64'} else machine}"
    return f"{system}-{machine}"


def execution_cache_dir(name: str) -> Path:
    """A shared cache under the execution root; never relocated by ``RNASEQ_WORK_ROOT``."""

    return (_execution_base() / "cache" / name).resolve()


def upstream_conda_cache() -> Path:
    """Shared nf-core process environments outside all per-run work directories."""

    return execution_cache_dir("upstream-conda")


def check_upstream_conda() -> RuntimeCheck:
    try:
        result = _run_capture(["conda", "--version"])
    except FileNotFoundError:
        return RuntimeCheck("Conda", "NOT FOUND", "Conda executable was not found on PATH.")
    if result.returncode != 0:
        return RuntimeCheck("Conda", "NOT FOUND", (result.stderr or result.stdout).strip())
    return RuntimeCheck("Conda", "FOUND", result.stdout.strip())


def prepare_upstream_conda_cache() -> Path:
    cache = upstream_conda_cache()
    try:
        cache.mkdir(parents=True, exist_ok=True)
        with tempfile.NamedTemporaryFile(prefix=".nf-rna-write-check-", dir=cache):
            pass
    except OSError as exc:
        raise ExecutionPreflightError(f"Upstream Conda cache is not writable: {cache}: {exc}") from exc
    return cache


def render_upstream_conda_config(cache: Path) -> str:
    # JSON string quoting is valid Groovy syntax and handles paths with spaces.
    # The former experimental osx-arm64 per-process Conda overrides were
    # withdrawn; macOS is not a supported runtime.
    return f"conda.enabled = true\ndocker.enabled = false\nconda.cacheDir = {json.dumps(str(cache))}\n"


def native_runtime_snapshot() -> RuntimeSnapshot:
    """Collect host capacity for the native nf-core Conda path."""

    return RuntimeSnapshot(
        host_os=platform.system(), host_architecture=platform.machine().lower(),
        logical_cpus=detect_local_resource_capacity().logical_cpus,
        host_memory_bytes=_native_usable_memory_bytes(),
    )


def _native_usable_memory_bytes() -> int | None:
    total = _host_memory_bytes()
    return detect_process_limits().apply_memory(total) if platform.system().lower() == "linux" else total


def prepare_execution_workspace(workspace: ExecutionWorkspace) -> None:
    """Create local operational directories, never deleting existing scratch data."""

    workspace.launch_dir.mkdir(parents=True, exist_ok=True)
    workspace.work_dir.mkdir(parents=True, exist_ok=True)


def _utc_now() -> str:
    return datetime.now(UTC).replace(microsecond=0).isoformat().replace("+00:00", "Z")


def hisat2_workflow_identity() -> dict[str, str]:
    """Return the immutable identity of the first-party H2 workflow source."""

    digest = sha256_file(HISAT2_WORKFLOW)
    return {
        "name": "nf-rna/hisat2_featurecounts",
        "nf_rna_version": PIPELINE_VERSION,
        "workflow_path": "workflow/hisat2_featurecounts.nf",
        "workflow_sha256": digest,
    }


def resolved_upstream_implementation(report: ValidationReport) -> dict[str, object]:
    """Separate accepted legacy config metadata from the executable route."""

    assert report.config is not None
    method = report.config.upstream.quantification.method if report.config.upstream.quantification else "salmon"
    if method == "salmon":
        return {"implementation": {"name": "nf-core/rnaseq", "version": report.config.upstream.pipeline_version}}
    return {
        "implementation": hisat2_workflow_identity(),
        "legacy_fastq_config": {
            "engine": report.config.upstream.engine,
            "pipeline_version": report.config.upstream.pipeline_version,
            "meaning": "accepted compatibility metadata; not executed by the HISAT2 + featureCounts backend",
        },
    }


def _run_capture(arguments: list[str]) -> subprocess.CompletedProcess[str]:
    return subprocess.run(arguments, capture_output=True, text=True, check=False)


def _linux_memory_bytes() -> tuple[int | None, int | None]:
    """Read Linux/WSL host memory from procfs without an external utility."""

    fields = {
        pieces[0].rstrip(":"): int(pieces[1]) * 1024
        for line in Path("/proc/meminfo").read_text(encoding="utf-8").splitlines()
        if (pieces := line.split()) and len(pieces) >= 2
    }
    total = fields.get("MemTotal")
    return total, fields.get("MemAvailable", total)


def _darwin_memory_bytes() -> int | None:
    """Read macOS physical memory through its native sysctl interface."""

    result = subprocess.run(["sysctl", "-n", "hw.memsize"], capture_output=True, text=True, check=False)
    return int(result.stdout.strip()) if result.returncode == 0 else None


def _host_memory_bytes() -> int | None:
    system = platform.system().lower()
    try:
        if system == "linux":
            total, _available = _linux_memory_bytes()
            return total
        if system == "darwin":
            return _darwin_memory_bytes()
    except (OSError, ValueError, KeyError):
        return None
    try:
        return os.sysconf("SC_PAGE_SIZE") * os.sysconf("SC_PHYS_PAGES")
    except (AttributeError, OSError, ValueError):
        return None


def _is_wsl() -> bool:
    try:
        release_identity = f"{platform.release()} {Path('/proc/version').read_text(encoding='utf-8')}".lower()
    except OSError:
        return False
    return "microsoft" in release_identity or "wsl" in release_identity


def effective_resource_budget(
    snapshot: RuntimeSnapshot, budget: ResourceContract, *, policy: ResourcePolicy | None = None,
) -> EffectiveResourceBudget:
    """Resolve the portable project ceiling against this host's usable capacity."""

    host_memory = snapshot.host_memory_bytes // 1024**3 if snapshot.host_memory_bytes else None
    cpu_limits = [budget.cpus]
    memory_limits = [budget.memory_gib]
    warnings: list[str] = []
    if snapshot.logical_cpus is not None:
        cpu_limits.append(snapshot.logical_cpus)
    else:
        warnings.append("Host CPU capacity is unavailable; it could not constrain the project budget.")
    if host_memory is not None:
        memory_limits.append(host_memory)
    else:
        warnings.append("Host memory capacity is unavailable; it could not constrain the project budget.")
    effective_cpus = min(cpu_limits)
    effective_memory = min(memory_limits)
    clamped = effective_cpus < budget.cpus or effective_memory < budget.memory_gib
    if clamped:
        warnings.append(
            f"Project budget {budget.cpus} CPUs/{budget.memory_gib} GiB was clamped to "
            f"{effective_cpus} CPUs/{effective_memory} GiB for this run."
        )
    return EffectiveResourceBudget(
        requested_cpus=budget.cpus,
        requested_memory_gib=budget.memory_gib,
        host_cpus=snapshot.logical_cpus,
        host_memory_gib=host_memory,
        effective_cpus=effective_cpus,
        effective_memory_gib=effective_memory,
        clamped=clamped,
        warnings=tuple(warnings),
        policy=policy,
    )


def validate_effective_resource_budget(resources: EffectiveResourceBudget) -> None:
    """Require enough effective capacity for the largest enabled local task."""

    if resources.effective_cpus < LOCAL_RESOURCE_CEILING.cpus or resources.effective_memory_gib < LOCAL_RESOURCE_CEILING.memory_gib:
        raise ExecutionPreflightError(
            "Effective local capacity is "
            f"{resources.effective_cpus} CPUs/{resources.effective_memory_gib} GiB, but enabled local "
            "process contracts require at least 8 CPUs/12 GiB. Increase host capacity."
        )


def _gib(value: int | None) -> str:
    return "unavailable" if value is None else f"{value / (1024 ** 3):.1f} GiB"


def _budget_checks(
    snapshot: RuntimeSnapshot, budget: ResourceContract, policy: ResourcePolicy | None = None,
) -> tuple[RuntimeCheck, ...]:
    resources = effective_resource_budget(snapshot, budget, policy=policy)
    level = "WARN" if resources.warnings or resources.effective_cpus < 8 or resources.effective_memory_gib < 12 else None
    return (
        RuntimeCheck(
            "Project resource budget", "FOUND",
            f"requested aggregate ceiling={resources.requested_cpus} CPUs/{resources.requested_memory_gib} GiB"
            + (f" ({describe_policy(policy)})." if policy is not None else "."),
        ),
        RuntimeCheck(
            "Effective local budget", "NOT FOUND" if level else "FOUND",
            f"effective aggregate ceiling={resources.effective_cpus} CPUs/{resources.effective_memory_gib} GiB; "
            + (" ".join(resources.warnings) if resources.warnings else "independent ready tasks may run concurrently within this ceiling."),
            level,
        ),
    )


def render_local_resource_config(budget: ResourceContract = LOCAL_RESOURCE_CEILING) -> str:
    """Render only the aggregate local ceiling, preserving process requests."""

    return (
        f"// Effective aggregate local ceiling={budget.cpus} CPUs/{budget.memory_gib} GiB.\n"
        "// Nextflow schedules independent ready tasks within this shared budget; process-specific requests remain intact.\n"
        f"executor {{ cpus = {budget.cpus}; memory = '{budget.memory_gib}.GB' }}\n"
        "process {\n"
        f"  resourceLimits = [cpus: {budget.cpus}, memory: '{budget.memory_gib}.GB', time: '{budget.time_hours}.h']\n"
        "}\n"
    )


def classify_execution_failure(stage: str, returncode: int, stderr_path: Path, *, resource: ResourceContract) -> str:
    """Return a concise, deterministic infrastructure diagnosis without retrying."""

    try:
        stderr = stderr_path.read_text(encoding="utf-8", errors="replace")[-16000:]
    except OSError:
        stderr = ""
    requested = f"{resource.cpus} CPUs, {resource.memory_gib} GiB, {resource.time_hours} h"
    if returncode in {137, -9} or "killed" in stderr.lower():
        return (
            f"LIKELY_OOM: process={stage}; exit_code={returncode}; requested_resources={requested}; "
            "suggestion=check host memory and the frozen local resource contract; do not change scientific parameters automatically."
        )
    if "no space left on device" in stderr.lower():
        return f"LIKELY_DISK_EXHAUSTION: process={stage}; exit_code={returncode}; suggestion=free space in the local execution workspace."
    if "no such file" in stderr.lower() and ("staged" in stderr.lower() or "downstream_inputs" in stderr.lower()):
        return f"MISSING_STAGED_INPUT: process={stage}; exit_code={returncode}; suggestion=inspect frozen handoff and downstream_inputs staging."
    return f"{stage} failed with return code {returncode}. Logs: {stderr_path.parent}"


def check_nextflow() -> RuntimeCheck:
    """Check for Nextflow without installing or invoking a pipeline."""

    try:
        result = _run_capture(["nextflow", "-version"])
    except FileNotFoundError:
        return RuntimeCheck("Nextflow", "NOT FOUND", "Nextflow executable was not found on PATH.")
    if result.returncode != 0:
        detail = (result.stderr or result.stdout).strip() or "nextflow -version failed."
        return RuntimeCheck("Nextflow", "NOT FOUND", detail)
    output = (result.stdout + "\n" + result.stderr).strip()
    match = re.search(r"version\s+([0-9][^\s]*)", output, flags=re.IGNORECASE)
    return RuntimeCheck("Nextflow", "FOUND", match.group(1) if match else output)


def _reference_runtime_check(project_dir: Path | None) -> RuntimeCheck:
    if project_dir is None:
        return RuntimeCheck("Reference readiness", "NOT FOUND", "Project-specific; run 'rnaseq doctor PROJECT' to inspect selected backend index readiness.", "WARN")
    try:
        from rnaseq.validators import validate_project

        report = validate_project(project_dir)
    except OSError as exc:
        return RuntimeCheck("Reference readiness", "NOT FOUND", f"Cannot inspect project reference: {exc}", "FAIL")
    if not report.is_valid or report.config is None:
        details = "; ".join(issue.message for issue in report.errors[:3]) or "reference readiness is unavailable"
        return RuntimeCheck("Reference readiness", "NOT FOUND", f"Project validation failed: {details}", "FAIL")
    if report.config.input.type is InputType.RAW_COUNTS:
        return RuntimeCheck("Reference readiness", "FOUND", "Raw-count route does not require a FASTQ index.")
    method = report.config.upstream.quantification.method if report.config.upstream.quantification else "salmon"
    reference = report.local_reference
    if reference is not None:
        status = reference.salmon_status if method == "salmon" else reference.hisat2_status
        detail = f"backend={method}; local index_status={status}"
        if status == "built":
            return RuntimeCheck("Reference readiness", "FOUND", detail)
        return RuntimeCheck("Reference readiness", "NOT FOUND", detail + "; execution would require an unavailable index.", "FAIL")
    return RuntimeCheck("Reference readiness", "FOUND", f"source={report.config.reference.source}; no dynamic reference preparation is performed by doctor.", "WARN")


def _version_tuple(text: str) -> tuple[int, ...] | None:
    match = re.match(r"(\d+)(?:\.(\d+))?(?:\.(\d+))?", text.strip())
    return tuple(int(part) for part in match.groups() if part is not None) if match else None


def runtime_policy_checks() -> tuple[RuntimeCheck, ...]:
    """Name the host platform and the automatically selected production backend."""

    try:
        backend = execution_backend()
    except ExecutionPreflightError as exc:
        return (RuntimeCheck("Runtime platform", "NOT FOUND", f"{runtime_platform_label()}; {exc}"),)
    return (
        RuntimeCheck("Runtime platform", "FOUND", runtime_platform_label()),
        RuntimeCheck("Execution backend", "FOUND", BACKEND_LABELS[backend]),
    )


def check_java() -> RuntimeCheck:
    """Resolve Java the way the Nextflow launcher does, then require Java 17+."""

    home = os.environ.get("NXF_JAVA_HOME") or os.environ.get("JAVA_HOME")
    java = str(Path(home) / "bin" / "java") if home else shutil.which("java")
    source = "NXF_JAVA_HOME" if os.environ.get("NXF_JAVA_HOME") else "JAVA_HOME" if home else "PATH"
    if java is None:
        return RuntimeCheck("Java", "NOT FOUND", f"No java executable on PATH; Nextflow requires Java {JAVA_MINIMUM_MAJOR} or later.")
    try:
        result = _run_capture([java, "-version"])
    except FileNotFoundError:
        return RuntimeCheck("Java", "NOT FOUND", f"{java} (from {source}) does not exist; Nextflow requires Java {JAVA_MINIMUM_MAJOR} or later.")
    output = (result.stderr + "\n" + result.stdout).strip()
    match = re.search(r'version "(\d+)(?:\.(\d+))?', output)
    if result.returncode != 0 or match is None:
        return RuntimeCheck("Java", "NOT FOUND", f"{java} -version failed: {output or 'no output'}")
    major = int(match.group(2)) if match.group(1) == "1" and match.group(2) else int(match.group(1))
    detail = f"{output.splitlines()[0]}; executable={java} (from {source})"
    if major < JAVA_MINIMUM_MAJOR:
        return RuntimeCheck("Java", "NOT FOUND", f"{detail}; Nextflow requires Java {JAVA_MINIMUM_MAJOR} or later.")
    return RuntimeCheck("Java", "FOUND", detail)


def check_nextflow_suitability() -> RuntimeCheck:
    found = check_nextflow()
    if found.state != "FOUND":
        return found
    minimum = NEXTFLOW_MINIMUM_VERSION
    version = _version_tuple(found.detail)
    if version is None or version < _version_tuple(minimum):
        return RuntimeCheck("Nextflow", "NOT FOUND", f"version={found.detail}; nf-core/rnaseq {NFCORE_RNASEQ_VERSION} requires Nextflow >={minimum}.")
    return RuntimeCheck("Nextflow", "FOUND", f"version={found.detail}; required >={minimum} by nf-core/rnaseq {NFCORE_RNASEQ_VERSION}.")


def check_conda_functional() -> RuntimeCheck:
    """Ask Conda itself for its platform; this is small and never solves or installs."""

    executable = shutil.which("conda")
    if executable is None:
        return RuntimeCheck("Conda", "NOT FOUND", "Conda executable was not found on PATH; native upstream and downstream environments are Conda-managed.")
    try:
        result = _run_capture([executable, "info", "--json"])
        info = json.loads(result.stdout) if result.returncode == 0 else None
    except (FileNotFoundError, json.JSONDecodeError):
        info = None
    if not isinstance(info, dict):
        return RuntimeCheck("Conda", "NOT FOUND", f"{executable} info --json failed; Conda is not functional.")
    subdir = info.get("platform")
    detail = f"conda {info.get('conda_version', 'unknown')}; executable={executable}; platform={subdir}"
    try:
        expected = native_platform()
    except ExecutionPreflightError:
        expected = None
    if expected is not None and subdir != expected:
        return RuntimeCheck("Conda", "NOT FOUND", f"{detail}; expected Conda platform {expected}.")
    return RuntimeCheck("Conda", "FOUND", detail)


# nf-core/rnaseq -profile conda requires both channels, conda-forge first
# (utils_nextflow_pipeline checkCondaChannels); other channels may also exist.
NFCORE_CONDA_CHANNELS = ("conda-forge", "bioconda")
NFCORE_CONDA_CHANNEL_REMEDY = (
    "conda config --add channels bioconda && conda config --add channels conda-forge "
    "&& conda config --set channel_priority strict"
)


def check_conda_channels() -> RuntimeCheck:
    """Apply nf-core's channel rule to Conda's effective, merged configuration."""

    name = "Conda channels for nf-core"
    required = "conda-forge before bioconda"
    executable = shutil.which("conda")
    if executable is None:
        return RuntimeCheck(name, "NOT FOUND", "Conda executable was not found on PATH; channel configuration cannot be verified.")
    command = [executable, "config", "--show", "channels", "channel_priority", "--json"]
    try:
        result = _run_capture(command)
        payload = json.loads(result.stdout) if result.returncode == 0 else None
    except (OSError, json.JSONDecodeError):
        payload = None
    channels = payload.get("channels") if isinstance(payload, dict) else None
    if not isinstance(channels, list) or not all(isinstance(channel, str) for channel in channels):
        return RuntimeCheck(
            name, "NOT FOUND",
            f"Could not read the effective channel list from '{' '.join(command[1:])}'; required channels: {required}. "
            "Inspect with 'conda config --show-sources'.",
        )
    detail = (
        f"observed channels={channels}; channel_priority={payload.get('channel_priority', 'unknown')}; "
        f"required by nf-core/rnaseq {NFCORE_RNASEQ_VERSION} -profile conda: {required}"
    )
    missing = [channel for channel in NFCORE_CONDA_CHANNELS if channel not in channels]
    if missing:
        return RuntimeCheck(name, "NOT FOUND", f"{detail}; missing {', '.join(missing)}. Fix with: {NFCORE_CONDA_CHANNEL_REMEDY}")
    if [channel for channel in channels if channel in NFCORE_CONDA_CHANNELS] != list(NFCORE_CONDA_CHANNELS):
        return RuntimeCheck(name, "NOT FOUND", f"{detail}; conda-forge must come before bioconda. Fix with: {NFCORE_CONDA_CHANNEL_REMEDY}")
    return RuntimeCheck(name, "FOUND", detail)


def _project_uses_nfcore(project_dir: Path | None) -> bool:
    """Whether the nf-core Conda upstream applies; without a project it may."""

    if project_dir is None:
        return True
    try:
        from rnaseq.project import load_project
        config = load_project(project_dir).config
    except (OSError, ValueError):
        return True
    method = config.upstream.quantification.method if config.upstream.quantification else "salmon"
    return config.input.type is InputType.FASTQ and method == "salmon"


def _nfcore_conda_channel_check(project_dir: Path | None) -> RuntimeCheck:
    check = check_conda_channels()
    if check.state == "FOUND" or _project_uses_nfcore(project_dir):
        return check
    return RuntimeCheck(check.name, check.state, f"{check.detail} (not used by this project's route)", "WARN")


def doctor_readiness(checks: tuple[RuntimeCheck, ...]) -> tuple[bool, tuple[str, ...]]:
    """Ready only when no check failed; WARN never blocks."""

    failed = tuple(check.name for check in checks if check.verdict == "FAIL")
    return not failed, failed


def _writable_location_check(name: str, path: Path, purpose: str) -> RuntimeCheck:
    """Check writability without creating anything: the nearest existing ancestor must be writable."""

    anchor = path
    while not anchor.exists() and anchor != anchor.parent:
        anchor = anchor.parent
    if not anchor.is_dir() or not os.access(anchor, os.W_OK | os.X_OK):
        return RuntimeCheck(name, "NOT FOUND", f"path={path}; {anchor} is not a writable directory; {purpose} cannot be provisioned.")
    state = "exists" if anchor == path else f"will be created under {anchor}"
    return RuntimeCheck(name, "FOUND", f"path={path}; writable ({state}); {purpose}.")


def nextflow_home() -> Path:
    return Path(os.environ.get("NXF_HOME") or Path.home() / ".nextflow").expanduser()


def nfcore_pin_check() -> RuntimeCheck:
    detail = f"nf-core/rnaseq {NFCORE_RNASEQ_VERSION} revision={NFCORE_RNASEQ_REVISION}; -profile {NFCORE_CONDA_PROFILE}"
    assets = nextflow_home() / "assets"
    cached = (assets / ".repos" / "nf-core" / "rnaseq" / "clones" / NFCORE_RNASEQ_REVISION, assets / "nf-core" / "rnaseq")
    for clone in cached:
        if (clone / "main.nf").is_file():
            return RuntimeCheck("nf-core/rnaseq pin", "FOUND", f"{detail}; cached at {clone}.")
    return RuntimeCheck("nf-core/rnaseq pin", "NOT FOUND", f"{detail}; not cached under {assets}; Nextflow pulls it on the first FASTQ run (network required).", "WARN")


def hisat2_conda_env_check() -> RuntimeCheck:
    """Verify the reviewed HISAT2/featureCounts environment identity and exact pins."""

    name = "HISAT2/featureCounts Conda environment"
    if not HISAT2_LINUX_CONDA_ENV.is_file():
        return RuntimeCheck(name, "NOT FOUND", f"Missing {HISAT2_LINUX_CONDA_ENV}.")
    digest = sha256_file(HISAT2_LINUX_CONDA_ENV)
    if digest != HISAT2_LINUX_CONDA_ENV_SHA256:
        return RuntimeCheck(name, "NOT FOUND", f"sha256 mismatch for {HISAT2_LINUX_CONDA_ENV.name}: expected {HISAT2_LINUX_CONDA_ENV_SHA256}, observed {digest}.")
    try:
        dependencies = yaml.safe_load(HISAT2_LINUX_CONDA_ENV.read_text(encoding="utf-8"))["dependencies"]
    except (yaml.YAMLError, KeyError, TypeError) as exc:
        return RuntimeCheck(name, "NOT FOUND", f"{HISAT2_LINUX_CONDA_ENV.name} is not a valid Conda environment: {exc}")
    loose = [item for item in dependencies if not (isinstance(item, str) and re.fullmatch(r"[\w.-]+::[^=\s]+=[^=\s]+=[^=\s]+", item))]
    if loose:
        return RuntimeCheck(name, "NOT FOUND", f"Not exactly pinned (channel::name=version=build): {', '.join(map(str, loose[:3]))}")
    versions = {item.split("::", 1)[1].split("=")[0]: item.split("=")[1] for item in dependencies}
    expected = {"hisat2": HISAT2_VERSION, "samtools": SAMTOOLS_VERSION, "subread": SUBREAD_VERSION}
    drift = [f"{tool}={versions.get(tool)} (expected {version})" for tool, version in expected.items() if versions.get(tool) != version]
    if drift:
        return RuntimeCheck(name, "NOT FOUND", "Tool versions disagree with the scientific contract: " + ", ".join(drift))
    reused = sorted(upstream_conda_cache().glob(f"env-*/conda-meta/hisat2-{HISAT2_VERSION}-*.json"))
    cache_state = f"provisioned env reused from {reused[0].parents[1]}" if reused else "Nextflow provisions it on the first HISAT2 run"
    return RuntimeCheck(
        name, "FOUND",
        f"{HISAT2_LINUX_CONDA_ENV.name} sha256={digest}; {len(dependencies)} exact pins; "
        f"hisat2={HISAT2_VERSION} samtools={SAMTOOLS_VERSION} subread={SUBREAD_VERSION}; {cache_state}.",
    )


def downstream_runtime_checks() -> tuple[RuntimeCheck, ...]:
    """Verify the reviewed downstream lock and that its runtime can be provisioned or reused."""

    from rnaseq import downstream_runtime

    try:
        lock = downstream_runtime.lock_path_for_platform()
        lock_sha = downstream_runtime._lock_checksum(lock)
    except ExecutionPreflightError as exc:
        return (RuntimeCheck("Downstream Conda lock", "NOT FOUND", str(exc)),)
    checks = [RuntimeCheck("Downstream Conda lock", "FOUND", f"{lock.name} sha256={lock_sha} matches SHA256SUMS.")]
    root = downstream_runtime._runtime_root()
    # One prefix per (lock, wheel); the lock-only name is the pre-1.3.0 layout.
    stem = f"nf-rna-downstream-{downstream_runtime.native_platform()}-{lock_sha[:16]}"
    prefixes = sorted(path for path in (root / "prefixes").glob(f"{stem}*") if path.is_dir()) if (root / "prefixes").is_dir() else []
    location = _writable_location_check("Downstream runtime location", root, "the locked downstream prefix and wheel cache")
    if location.state == "FOUND" and prefixes:
        unmarked = [path for path in prefixes if not (path / "runtime" / "nf-rna-wheel.json").is_file()]
        if unmarked:
            location = RuntimeCheck(location.name, "NOT FOUND", f"{unmarked[0]} exists without its wheel identity marker; remove it manually before running.")
        else:
            location = RuntimeCheck(location.name, "FOUND", f"{location.detail} Provisioned prefixes for this lock: {len(prefixes)}.")
    checks.append(location)
    ready, detail = downstream_runtime.wheel_source_status()
    checks.append(RuntimeCheck("Downstream nf-rna wheel", "FOUND" if ready else "NOT FOUND", detail))
    try:
        revision = downstream_runtime.require_identified_source()
        checks.append(RuntimeCheck("nf-rna source revision", "FOUND", f"{revision} (recorded in every run's provenance)."))
    except ExecutionPreflightError as exc:
        checks.append(RuntimeCheck("nf-rna source revision", "NOT FOUND", str(exc)))
    return tuple(checks)


def _work_root_check() -> RuntimeCheck:
    name = "Execution work root"
    try:
        work_root = resolve_execution_workspace("doctor", "resource-check").work_dir.parents[2]
    except ExecutionPreflightError as exc:
        return RuntimeCheck(name, "NOT FOUND", str(exc))
    source = WORK_ROOT_ENV if os.environ.get(WORK_ROOT_ENV) else f"execution root; set {WORK_ROOT_ENV} to relocate"
    return _writable_location_check(name, work_root, f"Nextflow task work directories ({source})")


def _disk_free_check(name: str, path: Path) -> RuntimeCheck:
    disk_probe = path
    while not disk_probe.exists() and disk_probe != disk_probe.parent:
        disk_probe = disk_probe.parent
    try:
        free = shutil.disk_usage(disk_probe).free
        return RuntimeCheck(name, "FOUND", f"path={path}; available_at={disk_probe}; free={_gib(free)}")
    except OSError as exc:
        return RuntimeCheck(name, "NOT FOUND", f"path={path}; unable to inspect free space: {exc}", "WARN")


def _free_space_checks() -> tuple[RuntimeCheck, ...]:
    """Free space where task work lands and, when relocated, where the launch state and caches stay."""

    try:
        workspace = resolve_execution_workspace("doctor", "resource-check").work_dir
    except ExecutionPreflightError as exc:
        return (RuntimeCheck("Execution work-directory free space", "NOT FOUND", str(exc), "WARN"),)
    checks = [_disk_free_check("Execution work-directory free space", workspace)]
    if os.environ.get(WORK_ROOT_ENV):
        checks.append(_disk_free_check("Execution root free space", _execution_base().resolve()))
    return tuple(checks)


def _doctor_budget(project_dir: Path | None) -> tuple[ResourceContract, ResourcePolicy | None]:
    budget = LOCAL_RESOURCE_CEILING
    policy: ResourcePolicy | None = None
    if project_dir is not None:
        try:
            from rnaseq.project import load_project
            config = load_project(project_dir).config
            policy = project_resource_policy(config)
            budget = project_execution_budget(config, policy)
        except (OSError, ValueError):
            pass
    return budget, policy


def native_linux_doctor_checks(project_dir: Path | None = None) -> tuple[RuntimeCheck, ...]:
    """Checks for the qualified linux-64 architecture: Nextflow + Conda, no Docker."""

    budget, policy = _doctor_budget(project_dir)
    snapshot = native_runtime_snapshot()
    probe = Path.cwd()
    return (
        RuntimeCheck("Python", "FOUND", f"{sys.executable} ({platform.python_version()})"),
        check_java(),
        check_nextflow_suitability(),
        check_conda_functional(),
        _nfcore_conda_channel_check(project_dir),
        nfcore_pin_check(),
        _writable_location_check("Execution root", _execution_base().resolve(), "Nextflow launch directories and runtime caches"),
        _work_root_check(),
        _writable_location_check("Upstream Conda cache", upstream_conda_cache(), "nf-core and HISAT2/featureCounts process environments"),
        hisat2_conda_env_check(),
        *downstream_runtime_checks(),
        RuntimeCheck(
            "Host runtime", "FOUND",
            f"OS={snapshot.host_os}; architecture={snapshot.host_architecture}; logical_cpus={snapshot.logical_cpus or 'unavailable'}; memory={_gib(snapshot.host_memory_bytes)}",
        ),
        *_budget_checks(snapshot, budget, policy),
        _reference_runtime_check(project_dir),
        RuntimeCheck("Disk write access", "FOUND" if os.access(probe, os.W_OK) else "NOT FOUND", str(probe)),
        *_free_space_checks(),
    )


def doctor_checks(project_dir: Path | None = None) -> tuple[RuntimeCheck, ...]:
    """Return inexpensive, non-mutating execution prerequisite checks."""

    policy = runtime_policy_checks()
    if len(policy) == 1:
        return policy
    return (*policy, *native_linux_doctor_checks(project_dir))


def _require_fastq_execution_report(report: ValidationReport) -> None:
    if not report.is_valid:
        raise ExecutionPreflightError("Execution is blocked because validation failed.")
    if report.config is None:
        raise ExecutionPreflightError("Execution is blocked because project configuration is unavailable.")
    if report.config.input.type is InputType.RAW_COUNTS:
        raise ExecutionPreflightError(
            "Raw-count execution is not implemented in Milestone 2. "
            "Raw-count projects remain valid for planning and future downstream analysis."
        )
    if not report.execution_ready:
        blockers = "; ".join(report.execution_blockers) or "reference strategy is unresolved"
        raise ExecutionPreflightError(f"Execution is blocked: {blockers}.")
    if report.config.upstream.quantification is None:
        raise ExecutionPreflightError(
            "FASTQ execution configuration is incomplete: upstream.quantification.method is required."
        )
    if report.config.upstream.quantification.method not in {"salmon", "hisat2_featurecounts"}:
        raise ExecutionPreflightError("Unsupported FASTQ quantification method.")


def _load_current_plan_manifest(report: ValidationReport) -> dict[str, Any]:
    path = report.project_dir / "planning" / "manifest.preview.yaml"
    if not path.exists():
        raise ExecutionPreflightError("No planning manifest found. Run 'rnaseq plan PROJECT' before execution.")
    try:
        manifest = yaml.safe_load(path.read_text(encoding="utf-8"))
    except yaml.YAMLError as exc:
        raise ExecutionPreflightError(f"Planning manifest is invalid YAML: {exc}") from exc
    if not isinstance(manifest, dict):
        raise ExecutionPreflightError("Planning manifest must contain a YAML mapping.")
    return manifest


def require_fresh_plan(report: ValidationReport) -> None:
    """Block runs when source config or input identities diverge from the frozen plan."""

    saved = _load_current_plan_manifest(report)
    current = yaml.safe_load(render_manifest(report))
    for key in ("configuration", "input", "metadata", "contrasts", "upstream", "reference", "planned_downstream"):
        if saved.get(key) != current.get(key):
            raise ExecutionPreflightError(
                "Project inputs or configuration changed since planning. "
                "Run 'rnaseq plan PROJECT' again before execution."
            )


def _validate_custom_reference_files(report: ValidationReport) -> None:
    """Recheck reference assets immediately before launching Nextflow.

    The name is retained for the service import used by earlier releases.  Local
    references are checked here too so their checksum guarantee covers the gap
    between planning and execution.
    """

    assert report.config is not None
    reference = report.config.reference
    if reference.source == "local":
        try:
            report.local_reference = load_local_reference(reference, report.config.organism.species.value)
        except LocalReferenceError as exc:
            raise ExecutionPreflightError(str(exc)) from exc
        return
    if reference.source != "custom":
        return
    assert reference.fasta is not None and reference.gtf is not None
    required = (("fasta", reference.fasta), ("gtf", reference.gtf))
    optional = tuple(
        (label, value)
        for label, value in (("transcript_fasta", reference.transcript_fasta), ("salmon_index", reference.salmon_index), ("hisat2_index", reference.hisat2_index), ("hisat2_splice_sites", reference.hisat2_splice_sites))
        if value is not None
    )
    for label, value in required + optional:
        assert value is not None
        candidate = (report.project_dir / value).resolve()
        try:
            candidate.relative_to(report.project_dir.resolve())
        except ValueError as exc:
            raise ExecutionPreflightError(
                f"Custom reference {label} must use a relative path inside the project."
            ) from exc
        if label == "hisat2_index":
            from rnaseq.references import validate_hisat2_index
            try:
                validate_hisat2_index(candidate)
            except LocalReferenceError as exc:
                raise ExecutionPreflightError(str(exc)) from exc
        elif not candidate.is_file():
            raise ExecutionPreflightError(f"Custom reference file not found: {candidate}")


def prepare_run(report: ValidationReport, profile: str) -> PreparedRun:
    """Perform all non-mutating preflight checks required before confirmation."""

    if profile != LOCAL_PROFILE:
        raise ExecutionPreflightError(
            "Profile 'server' is not supported for execution in Milestone 2. "
            "Only '--profile local' is available."
        )
    _require_fastq_execution_report(report)
    assert report.config is not None
    validate_local_execution_budget(
        report.config.execution.max_cpus, report.config.execution.max_memory_gb, detect_local_resource_capacity()
    )
    resources = resolve_project_resources(report.config, native_runtime_snapshot())
    validate_effective_resource_budget(resources)
    _validate_custom_reference_files(report)
    require_fresh_plan(report)
    # Also prove the generated samplesheet remains renderable before any run directory exists.
    render_samplesheet(report)
    nextflow = check_nextflow()
    if nextflow.state != "FOUND":
        raise ExecutionPreflightError(
            "Nextflow is required for FASTQ execution but was not found. " + nextflow.detail
        )
    conda = check_upstream_conda()
    if conda.state != "FOUND":
        raise ExecutionPreflightError(
            "Conda is required for upstream FASTQ execution: " + conda.detail
        )
    prepare_upstream_conda_cache()
    resolve_execution_workspace("preflight", "paths")  # rejects an invalid RNASEQ_WORK_ROOT before a run exists
    return PreparedRun(report, nextflow.detail, "conda", resources)


def _render_execution_samplesheet(report: ValidationReport) -> str:
    """Freeze an absolute-path sheet so execution does not depend on mutable planning files."""

    assert report.fastq is not None and report.config is not None
    output = StringIO(newline="")
    writer = csv.writer(output, lineterminator="\n")
    writer.writerow(["sample", "fastq_1", "fastq_2", "strandedness"])
    for record in report.fastq.records:
        writer.writerow(
            [
                record.sample_id,
                str(record.fastq_1.resolve()),
                str(record.fastq_2.resolve()) if record.fastq_2 else "",
                report.config.upstream.strandedness,
            ]
        )
    return output.getvalue()


def _run_fingerprint(report: ValidationReport) -> str:
    return hashlib.sha256(render_manifest(report).encode("utf-8")).hexdigest()[:12]


def create_run_directory(report: ValidationReport) -> tuple[str, Path]:
    """Create a new immutable run location; never reuse a completed run directory."""

    stamp = datetime.now(UTC).strftime("%Y%m%dT%H%M%SZ")
    stem = f"run-{stamp}-{_run_fingerprint(report)}"
    runs_dir = report.project_dir / "runs"
    runs_dir.mkdir(exist_ok=True)
    for index in range(1, 1000):
        run_id = stem if index == 1 else f"{stem}-{index:02d}"
        run_dir = runs_dir / run_id
        try:
            run_dir.mkdir()
        except FileExistsError:
            continue
        for child in ("frozen", "upstream", "logs", "provenance", "handoff"):
            (run_dir / child).mkdir()
        return run_id, run_dir
    raise UpstreamExecutionError("Unable to allocate a unique run directory.")


def _write_text(path: Path, content: str) -> None:
    path.write_text(content, encoding="utf-8", newline="\n")


def _write_state(path: Path, **state: Any) -> None:
    _write_text(path, json.dumps(state, indent=2, sort_keys=True) + "\n")


def nfcore_runtime_params(report: ValidationReport) -> dict[str, bool]:
    """Return JSON-native nf-core booleans for the frozen params file.

    Nextflow's direct CLI parameter syntax can stringify a valueless option.
    Keep booleans in the JSON params file so nf-core schema validation receives
    actual JSON booleans rather than the string ``"true"``.
    """

    assert report.config is not None and report.config.input.type is InputType.FASTQ
    if report.config.upstream.quantification and report.config.upstream.quantification.method == "hisat2_featurecounts":
        return {}
    params = {"skip_alignment": True}
    if report.config.input.preprocessing is FastqPreprocessing.PRETRIMMED:
        params["skip_trimming"] = True
    return params


def _freeze_inputs(prepared: PreparedRun, run_dir: Path) -> tuple[Path, Path, Path, Path]:
    report = prepared.report
    assert report.loaded is not None and report.config is not None
    frozen = run_dir / "frozen"
    for source, destination in (
        (report.loaded.config_path, frozen / "project.yaml"),
        (report.loaded.metadata_path, frozen / "metadata.csv"),
        (report.loaded.contrasts_path, frozen / "contrasts.csv"),
    ):
        shutil.copy2(source, destination)
    samplesheet = frozen / "samplesheet.csv"
    _write_text(samplesheet, _render_execution_samplesheet(report))
    upstream = {
        "engine": "nextflow",
        **resolved_upstream_implementation(report),
        "execution_profile": LOCAL_PROFILE,
        "container_runtime": prepared.container_runtime,
        "reference": (
            report.local_reference.provenance()
            if report.local_reference is not None
            else report.config.reference.model_dump()
        ),
        "preprocessing": {
            "declared": report.config.input.preprocessing.value,
            "skip_trimming": report.config.input.preprocessing is FastqPreprocessing.PRETRIMMED,
            "nfcore_arguments": (
                ["--skip_trimming"]
                if report.config.input.preprocessing is FastqPreprocessing.PRETRIMMED else []
            ),
        },
        "quantification": report.config.upstream.quantification.model_dump(),
    }
    upstream_path = frozen / "upstream_run.yaml"
    _write_text(upstream_path, yaml.safe_dump(upstream, sort_keys=False, allow_unicode=True))
    params_path = frozen / "nextflow.params.json"
    _write_text(params_path, json.dumps(nfcore_runtime_params(report), sort_keys=True) + "\n")
    # Resource declarations are frozen separately from scientific parameters.
    runtime_config = frozen / "local.nextflow.config"
    effective = prepared.resource_budget or resolve_project_resources(report.config, native_runtime_snapshot())
    _write_text(runtime_config, render_local_resource_config(ResourceContract(
        "EFFECTIVE_LOCAL", effective.effective_cpus, effective.effective_memory_gib, LOCAL_RESOURCE_CEILING.time_hours
    )))
    _write_text(frozen / "nfcore.conda.config", render_upstream_conda_config(upstream_conda_cache()))
    return samplesheet, upstream_path, params_path, runtime_config


def build_nextflow_command(
    report: ValidationReport, *, samplesheet: Path, output_dir: Path, profile: str,
    params_file: Path | None = None, config_file: Path | None = None,
    conda_config_file: Path | None = None,
    reference_paths: dict[str, Path] | None = None, work_dir: Path | None = None,
    tuning_config_file: Path | None = None,
) -> list[str]:
    """Build the exact argument vector; it is never run through a shell."""

    _require_fastq_execution_report(report)
    if profile != LOCAL_PROFILE:
        raise ExecutionPreflightError("Only the local profile is executable in Milestone 2.")
    assert report.config is not None
    reference: ReferenceConfig = report.config.reference
    execution_backend()  # refuses unsupported platforms before any command exists
    command = ["nextflow", "run"]
    if config_file is not None:
        command.extend(["-c", str(config_file.resolve())])
    if conda_config_file is not None:
        command.extend(["-c", str(conda_config_file.resolve())])
    if tuning_config_file is not None:
        command.extend(["-c", str(tuning_config_file.resolve())])
    command.extend([
        "nf-core/rnaseq", "-r", NFCORE_RNASEQ_VERSION,
        "-profile", NFCORE_CONDA_PROFILE,
    ])
    if work_dir is not None:
        command.extend(["-work-dir", str(work_dir.resolve())])
    if params_file is not None:
        command.extend(["-params-file", str(params_file.resolve())])
    command.extend([
        "--input", str(samplesheet.resolve()), "--outdir", str(output_dir.resolve()),
        "--pseudo_aligner", "salmon",
    ])
    if reference.source == "igenomes":
        assert reference.genome is not None
        command.extend(["--genome", reference.genome])
    elif reference.source == "local":
        local_reference = report.local_reference
        if local_reference is None:
            try:
                local_reference = load_local_reference(reference, report.config.organism.species.value)
            except LocalReferenceError as exc:
                raise ExecutionPreflightError(str(exc)) from exc
        for option, path in local_reference.nfcore_arguments():
            command.extend([option, str((reference_paths or {}).get(option.removeprefix("--"), path))])
    else:
        assert reference.fasta is not None and reference.gtf is not None
        for option, value in (("--fasta", reference.fasta), ("--gtf", reference.gtf)):
            key = option.removeprefix("--")
            candidate = (reference_paths or {}).get(key, (report.project_dir / value).resolve())
            try:
                candidate.relative_to(report.project_dir.resolve())
            except ValueError as exc:
                raise ExecutionPreflightError(f"Custom reference {option} must remain inside the project.") from exc
            if not candidate.is_file():
                raise ExecutionPreflightError(f"Custom reference file not found: {candidate}")
            command.extend([option, str(candidate)])
        for option, value in (("--transcript_fasta", reference.transcript_fasta), ("--salmon_index", reference.salmon_index)):
            if value is None:
                continue
            key = option.removeprefix("--")
            candidate = (reference_paths or {}).get(key, (report.project_dir / value).resolve())
            try:
                candidate.relative_to(report.project_dir.resolve())
            except ValueError as exc:
                raise ExecutionPreflightError(f"Custom reference {option} must remain inside the project.") from exc
            if not candidate.is_file():
                raise ExecutionPreflightError(f"Custom reference file not found: {candidate}")
            command.extend([option, str(candidate)])
    return command


def build_hisat2_featurecounts_command(
    report: ValidationReport, *, samplesheet: Path, output_dir: Path, profile: str,
    reference_paths: dict[str, Path] | None = None, work_dir: Path | None = None,
    config_file: Path | None = None, conda_config_file: Path | None = None,
    observer_config_file: Path | None = None,
) -> list[str]:
    """Build the first-party alignment/counting command without a shell."""

    _require_fastq_execution_report(report)
    if profile != LOCAL_PROFILE:
        raise ExecutionPreflightError("Only '--profile local' is implemented.")
    assert report.config is not None and report.fastq is not None
    if report.config.upstream.quantification.method != "hisat2_featurecounts":
        raise ExecutionPreflightError("HISAT2 command requested for a non-HISAT2 project.")
    reference = report.config.reference
    paths = reference_paths or {}
    if reference.source == "local":
        local = report.local_reference
        if local is None:
            local = load_local_reference(reference, report.config.organism.species.value)
        reference_arguments = {
            "fasta": paths.get("fasta", local.genome_fasta.path),
            "gtf": paths.get("gtf", local.annotation_gtf.path),
            "hisat2_index": paths.get("hisat2_index", local.hisat2_index),
            "hisat2_splice_sites": paths.get("hisat2_splice_sites", local.hisat2_splice_sites.path if local.hisat2_splice_sites else None),
            "hisat2_index_basename": local.hisat2_index_prefix.name if local.hisat2_index_prefix else "genome",
        }
        use_runtime_splices = local.hisat2_strategy in {
            "genome_only_runtime_splicesites",
            "genome_only_runtime_splices",  # legacy manifest spelling
        }
    else:
        if not reference.fasta or not reference.gtf or not reference.hisat2_index:
            raise ExecutionPreflightError("HISAT2 requires reference.fasta, reference.gtf, and reference.hisat2_index.")
        reference_arguments = {
            key: paths.get(key, (report.project_dir / value).resolve())
            for key, value in (("fasta", reference.fasta), ("gtf", reference.gtf), ("hisat2_index", reference.hisat2_index))
        }
        reference_arguments["hisat2_splice_sites"] = paths.get("hisat2_splice_sites")
        reference_arguments["hisat2_index_basename"] = "genome"
        use_runtime_splices = False
    if any(reference_arguments[key] is None for key in ("fasta", "gtf", "hisat2_index")):
        raise ExecutionPreflightError("HISAT2 reference is not prepared.")
    execution_backend()  # refuses unsupported platforms before any command exists
    command = ["nextflow", "run"]
    if config_file is not None:
        command.extend(["-c", str(config_file.resolve())])
    if conda_config_file is None:
        raise ExecutionPreflightError("Linux HISAT2 execution requires a frozen Conda runtime config.")
    if not HISAT2_LINUX_CONDA_ENV.is_file():
        raise ExecutionPreflightError(f"Linux HISAT2 Conda environment is missing: {HISAT2_LINUX_CONDA_ENV}")
    command.extend(["-c", str(conda_config_file.resolve())])
    if observer_config_file is not None:
        command.extend(["-c", str(observer_config_file.resolve())])
    command.extend([
        str(HISAT2_WORKFLOW), "-profile", NFCORE_CONDA_PROFILE,
        "--input", str(samplesheet.resolve()), "--outdir", str(output_dir.resolve()),
        "--fasta", str(reference_arguments["fasta"]), "--gtf", str(reference_arguments["gtf"]),
        "--hisat2_index", str(reference_arguments["hisat2_index"]),
        "--hisat2_index_basename", str(reference_arguments["hisat2_index_basename"]),
        "--assembly_script", str((Path(__file__).resolve().parent / "hisat2_featurecounts.py")),
        "--layout", report.fastq.layout.value, "--strandedness", report.config.upstream.strandedness,
        "--pretrimmed", str(report.config.input.preprocessing is FastqPreprocessing.PRETRIMMED).lower(),
    ])
    if use_runtime_splices:
        splice_sites = reference_arguments["hisat2_splice_sites"]
        assert splice_sites is not None
        command.extend(["--hisat2_splice_sites", str(splice_sites), "--hisat2_use_runtime_splices", "true"])
    if work_dir is not None:
        command.extend(["-work-dir", str(work_dir.resolve())])
    return command


def _relative_to_run(run_dir: Path, path: Path) -> str:
    return path.resolve().relative_to(run_dir.resolve()).as_posix()


def generate_handoff_manifest(prepared: PreparedRun, run_dir: Path) -> Path:
    """Record documented nf-core output locations only after a successful run."""

    report = prepared.report
    assert report.config is not None and report.fastq is not None
    output = run_dir / "upstream" / "nfcore_rnaseq"
    counts = output / "salmon" / "salmon.merged.gene_counts.tsv"
    multiqc_reports = sorted((output / "multiqc").glob("**/multiqc_report.html"))
    if not counts.is_file():
        raise UpstreamExecutionError(
            "nf-core completed but the documented Salmon gene-count matrix was not found: " + str(counts)
        )
    if not multiqc_reports:
        raise UpstreamExecutionError("nf-core completed but no documented MultiQC report was found.")
    multiqc = multiqc_reports[0]
    multiqc_data = next(
        (candidate for candidate in (multiqc.parent / "multiqc_data", multiqc.parent / "multiqc_report_data") if candidate.is_dir()),
        None,
    )
    if multiqc_data is None:
        raise UpstreamExecutionError("nf-core completed but MultiQC data directory was not found.")
    with counts.open("r", encoding="utf-8", newline="") as handle:
        header = next(csv.reader(handle, delimiter="\t"), [])
    if not header:
        raise UpstreamExecutionError("The documented Salmon gene-count matrix is empty.")
    quant_files = {
        sample: output / "salmon" / sample / "quant.sf"
        for sample in report.fastq.sample_ids
    }
    tx2gene = output / "salmon" / "salmon.merged.tx2gene_augmented.tsv"
    missing_salmon = [str(path) for path in (*quant_files.values(), tx2gene) if not path.is_file()]
    if missing_salmon:
        raise UpstreamExecutionError(
            "nf-core completed but required Salmon transcript-level handoff artifacts were not found: "
            + ", ".join(missing_salmon)
        )
    manifest = {
        "pipeline": {"name": "nf-core/rnaseq", "version": report.config.upstream.pipeline_version},
        "quantification": {
            "method": "salmon",
            "route": "--pseudo_aligner salmon --skip_alignment",
        },
        "preprocessing": {
            "declared": report.config.input.preprocessing.value,
            "skip_trimming": report.config.input.preprocessing is FastqPreprocessing.PRETRIMMED,
            "nfcore_arguments": (
                ["--skip_trimming"]
                if report.config.input.preprocessing is FastqPreprocessing.PRETRIMMED else []
            ),
        },
        "reference": (
            report.local_reference.provenance()
            if report.local_reference is not None
            else report.config.reference.model_dump()
        ),
        "samples": list(report.fastq.sample_ids),
        "gene_level_counts": {
            "path": _relative_to_run(run_dir, counts),
            "format": "TSV",
            "identifier_column": header[0],
            "description": "nf-core/rnaseq Salmon merged gene-level estimated-count matrix.",
            "downstream_note": "Milestone 3 must choose its DESeq2 import strategy; no downstream analysis was run here.",
        },
        "salmon": {
            "quant_sf": {sample: _relative_to_run(run_dir, path) for sample, path in quant_files.items()},
            "tx2gene": {
                "path": _relative_to_run(run_dir, tx2gene),
                "sha256": sha256_file(tx2gene),
                "mapping_type": "nfcore_tx2gene_augmented",
                "role": (
                    "Mapping used by nf-core/rnaseq 3.26.0 tximport; includes self-mappings "
                    "for quantified transcripts absent from the GTF-derived mapping."
                ),
            },
            "transcript_counts": _relative_to_run(
                run_dir, output / "salmon" / "salmon.merged.transcript_counts.tsv"
            ),
            "contract_note": (
                "L1 imports per-sample quant.sf with salmon.merged.tx2gene_augmented.tsv via tximport; "
                "estimated counts are not silently rounded by Python."
            ),
        },
        "multiqc": {
            "html": _relative_to_run(run_dir, multiqc),
            "data_directory": _relative_to_run(run_dir, multiqc_data),
        },
    }
    path = run_dir / "handoff" / "upstream_manifest.yaml"
    _write_text(path, yaml.safe_dump(manifest, sort_keys=False, allow_unicode=True))
    return path


def generate_hisat2_featurecounts_handoff(report: ValidationReport, run_dir: Path) -> Path:
    """Validate the first-party H2/featureCounts boundary before downstream use."""

    assert report.config is not None and report.fastq is not None
    output = run_dir / "upstream" / "hisat2_featurecounts"
    matrix = output / "counts" / "canonical_counts.csv"
    sample_map = output / "counts" / "sample_map.csv"
    multiqc = output / "multiqc" / "multiqc_report.html"
    required = (matrix, sample_map, multiqc)
    missing = [str(path) for path in required if not path.is_file()]
    if missing:
        raise UpstreamExecutionError("HISAT2 + featureCounts completed without required artifacts: " + ", ".join(missing))
    with matrix.open(encoding="utf-8", newline="") as handle:
        header = next(csv.reader(handle), [])
    samples = list(report.fastq.sample_ids)
    if header != ["gene_id", *samples]:
        raise UpstreamExecutionError("Canonical featureCounts matrix columns do not match declared sample IDs.")
    manifest = {
        "pipeline": {**hisat2_workflow_identity(), "versions": {"hisat2": HISAT2_VERSION, "samtools": SAMTOOLS_VERSION, "subread": SUBREAD_VERSION}},
        "quantification": {"method": "hisat2_featurecounts", "route": "first_party_alignment_counting"},
        "counting": COUNTING_POLICY,
        "strandedness": report.config.upstream.strandedness,
        "samples": samples,
        "featurecounts": {
            "canonical_matrix": _relative_to_run(run_dir, matrix),
            "sample_map": _relative_to_run(run_dir, sample_map),
            "per_sample_directory": _relative_to_run(run_dir, output / "counts" / "per_sample"),
            "assignment_summary_directory": _relative_to_run(run_dir, output / "counts" / "per_sample"),
        },
        "alignment": {
            "bam_directory": _relative_to_run(run_dir, output / "bam"),
            "count_only_bam_directory": _relative_to_run(run_dir, output / "bam" / "count_only"),
            "count_only_transformation": "samtools view -bh -F 0x900; retain the original coordinate-sorted BAM and its tags",
        },
        "multiqc": {"html": _relative_to_run(run_dir, multiqc)},
        "reference": report.local_reference.provenance() if report.local_reference is not None else report.config.reference.model_dump(),
    }
    path = run_dir / "handoff" / "upstream_manifest.yaml"
    _write_text(path, yaml.safe_dump(manifest, sort_keys=False, allow_unicode=True))
    return path


def execute_prepared_run(prepared: PreparedRun) -> RunResult:
    """Freeze inputs, run Nextflow, preserve logs, and write a standard handoff."""

    run_id, run_dir = create_run_directory(prepared.report)
    state_path = run_dir / "run_state.json"
    started = _utc_now()
    _write_state(state_path, run_id=run_id, status="CREATED", started_at=None, completed_at=None, return_code=None)
    samplesheet, _, params_file, runtime_config = _freeze_inputs(prepared, run_dir)
    workspace = resolve_execution_workspace("legacy", run_id)
    command = build_nextflow_command(
        prepared.report,
        samplesheet=samplesheet,
        output_dir=run_dir / "upstream" / "nfcore_rnaseq",
        profile=LOCAL_PROFILE,
        params_file=params_file,
        config_file=runtime_config,
        conda_config_file=run_dir / "frozen" / "nfcore.conda.config",
        work_dir=workspace.work_dir / "upstream",
    )
    runtime = native_runtime_snapshot()
    resources = prepared.resource_budget or resolve_project_resources(prepared.report.config, runtime)
    source_root = Path(__file__).resolve().parents[2]
    source_checkout = source_root if (source_root / ".git").exists() else None
    git_commit: str | None = None
    try:
        if source_checkout is not None:
            git_result = _run_capture(["git", "-C", str(source_checkout), "rev-parse", "HEAD"])
            if git_result.returncode == 0 and isinstance(getattr(git_result, "stdout", None), str):
                git_commit = git_result.stdout.strip() or None
    except (OSError, ValueError, AttributeError):
        pass
    provenance = {
        "nextflow_version": prepared.nextflow_version,
        "nfcore_rnaseq_version": NFCORE_RNASEQ_VERSION,
        "execution_profile": LOCAL_PROFILE,
        "container_runtime": prepared.container_runtime,
        "upstream_runtime": {"kind": "conda", "profile": NFCORE_CONDA_PROFILE,
                             "revision": NFCORE_RNASEQ_REVISION, "cache_dir": str(upstream_conda_cache()),
                             "platform": native_platform()},
        "production_intended": prepared.report.config.reference.acceptance == "production",
        "git_commit": git_commit,
        "source_checkout": str(source_checkout) if source_checkout is not None else None,
        "workflow_sha256": {
            "workflow/main.nf": sha256_file(workflow_asset_path("main.nf")),
            "workflow/hisat2_featurecounts.nf": sha256_file(HISAT2_WORKFLOW),
        },
        "execution_root": str(workspace.root),
        "execution_launch_dir": str(workspace.launch_dir),
        "execution_work_dir": str(workspace.work_dir),
        "command": command,
        "input_identity": yaml.safe_load(render_manifest(prepared.report))["input"],
        "fastq_preprocessing": prepared.report.config.input.preprocessing.value,
        "skip_trimming": prepared.report.config.input.preprocessing is FastqPreprocessing.PRETRIMMED,
        "runtime_resources": {
            **resources.as_dict(),
            "host_architecture": runtime.host_architecture,
            # Record-schema compatibility with the retired Docker backend.
            "docker_architecture": None,
            "docker_memory_bytes": None,
            "first_party_image_architecture": None,
            "resource_profile": "M5_LOCAL_SMALL_MEDIUM_LARGE",
        },
        "frozen_local_nextflow_config": {
            "path": "frozen/local.nextflow.config",
            "sha256": sha256_file(run_dir / "frozen" / "local.nextflow.config"),
        },
    }
    _write_text(run_dir / "provenance" / "run_provenance.yaml", yaml.safe_dump(provenance, sort_keys=False))
    _write_state(state_path, run_id=run_id, status="RUNNING", started_at=started, completed_at=None, return_code=None)
    stdout_path, stderr_path = run_dir / "logs" / "nextflow.stdout.log", run_dir / "logs" / "nextflow.stderr.log"
    try:
        prepare_execution_workspace(workspace)
        with stdout_path.open("w", encoding="utf-8", newline="\n") as stdout, stderr_path.open("w", encoding="utf-8", newline="\n") as stderr:
            result = subprocess.run(command, cwd=workspace.launch_dir, stdout=stdout, stderr=stderr, check=False)
    except OSError as exc:
        _write_state(state_path, run_id=run_id, status="FAILED", started_at=started, completed_at=_utc_now(), return_code=None, error=str(exc))
        raise UpstreamExecutionError(f"Unable to launch Nextflow: {exc}. Logs: {run_dir / 'logs'}") from exc
    if result.returncode != 0:
        _write_state(state_path, run_id=run_id, status="FAILED", started_at=started, completed_at=_utc_now(), return_code=result.returncode)
        raise UpstreamExecutionError(
            classify_execution_failure(
                "nf-core/rnaseq", result.returncode, stderr_path, resource=RESOURCE_CONTRACTS["MEDIUM"]
            )
        )
    try:
        handoff = generate_handoff_manifest(prepared, run_dir)
    except UpstreamExecutionError as exc:
        _write_state(state_path, run_id=run_id, status="FAILED", started_at=started, completed_at=_utc_now(), return_code=0, error=str(exc))
        raise
    handoff_payload = yaml.safe_load(handoff.read_text(encoding="utf-8"))
    mapping = handoff_payload["salmon"]["tx2gene"]
    provenance["salmon_tx2gene"] = {
        "path": str((run_dir / mapping["path"]).resolve()),
        "sha256": mapping["sha256"],
        "mapping_type": mapping["mapping_type"],
        "role": mapping["role"],
    }
    _write_text(run_dir / "provenance" / "run_provenance.yaml", yaml.safe_dump(provenance, sort_keys=False))
    _write_state(state_path, run_id=run_id, status="SUCCESS", started_at=started, completed_at=_utc_now(), return_code=0)
    return RunResult(run_id, run_dir, state_path, handoff)


def load_run_states(project_dir: Path) -> list[dict[str, Any]]:
    states: list[dict[str, Any]] = []
    runs_dir = project_dir / "runs"
    if not runs_dir.is_dir():
        return states
    # Support both legacy flat runs and the production case/run hierarchy.
    paths = set(runs_dir.glob("*/run_state.json")) | set(runs_dir.glob("*/*/run_state.json"))
    for path in sorted(paths, reverse=True):
        try:
            state = json.loads(path.read_text(encoding="utf-8"))
        except (json.JSONDecodeError, OSError):
            continue
        if isinstance(state, dict) and state.get("status") in RUN_STATES:
            state["run_dir"] = str(path.parent)
            state["handoff_available"] = (path.parent / "handoff" / "upstream_manifest.yaml").is_file()
            delivery = path.parent / "delivery"
            dated_readme = any(
                candidate.is_file() and re.fullmatch(r"README_[0-9]{8}\.md", candidate.name)
                for candidate in delivery.iterdir()
            ) if delivery.is_dir() else False
            state["delivery_available"] = (delivery / "README.md").is_file() or dated_readme
            states.append(state)
    return states
