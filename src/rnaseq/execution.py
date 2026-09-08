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
from dataclasses import dataclass
from datetime import UTC, datetime
from io import StringIO
from pathlib import Path
from typing import Any

import yaml

from rnaseq.errors import ExecutionPreflightError, UpstreamExecutionError
from rnaseq.models import FastqPreprocessing, InputType, NFCORE_RNASEQ_VERSION, PIPELINE_VERSION, ReferenceConfig
from rnaseq.hisat2_featurecounts import COUNTING_POLICY, HISAT2_VERSION, SAMTOOLS_VERSION, SUBREAD_VERSION
from rnaseq.planner import render_manifest, render_samplesheet
from rnaseq.references import LocalReferenceError, load_local_reference, sha256_file
from rnaseq.validators import ValidationReport

LOCAL_PROFILE = "local"
CONTAINER_PROFILE = "docker"
RUN_STATES = {"CREATED", "RUNNING", "SUCCESS", "FAILED"}
EXECUTION_ROOT_ENV = "RNASEQ_EXECUTION_ROOT"
CONTROL_PLANE_IMAGE = "rnaseq-control-plane:latest"
HISAT2_WORKFLOW = Path(__file__).resolve().parents[2] / "workflow" / "hisat2_featurecounts.nf"
CONTAINER_R_PACKAGES = (
    "DESeq2", "tximport", "ggplot2", "pheatmap", "yaml", "jsonlite",
    "clusterProfiler", "AnnotationDbi", "org.Hs.eg.db", "org.Mm.eg.db",
)


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
    """Read host capacity without Docker; Linux/WSL and macOS are supported."""

    cpus = os.cpu_count()
    total: int | None = None
    available: int | None = None
    try:
        if sys.platform == "darwin":
            result = subprocess.run(["sysctl", "-n", "hw.memsize"], capture_output=True, text=True, check=False)
            total = int(result.stdout.strip()) // 1024**3 if result.returncode == 0 else None
            available = total
        else:
            fields = {
                pieces[0].rstrip(":"): pieces[1]
                for line in Path("/proc/meminfo").read_text(encoding="utf-8").splitlines()
                if (pieces := line.split()) and len(pieces) >= 2
            }
            total = int(fields["MemTotal"]) // 1024**2
            available = int(fields.get("MemAvailable", fields["MemTotal"])) // 1024**2
    except (OSError, ValueError, KeyError):
        pass
    return LocalResourceCapacity(cpus if cpus and cpus > 0 else None, total, available)


def suggested_local_resources(capacity: LocalResourceCapacity) -> tuple[int, int]:
    """Reserve roughly one fifth of CPUs and memory, rounded conservatively."""

    fallback = (LOCAL_RESOURCE_CEILING.cpus, LOCAL_RESOURCE_CEILING.memory_gib)
    if capacity.logical_cpus is None or capacity.total_memory_gib is None:
        return fallback
    cpus = max(1, math.floor(capacity.logical_cpus * 0.8))
    usable_memory = min(capacity.total_memory_gib, capacity.available_memory_gib or capacity.total_memory_gib)
    memory = max(1, (math.floor(usable_memory * 0.8) // 4) * 4)
    return cpus, memory


def validate_local_execution_budget(cpus: int, memory_gb: int, capacity: LocalResourceCapacity) -> None:
    if isinstance(cpus, bool) or isinstance(memory_gb, bool) or cpus < 1 or memory_gb < 1:
        raise ExecutionPreflightError("Execution CPU and memory limits must be positive integers.")
    if cpus < LOCAL_RESOURCE_CEILING.cpus or memory_gb < LOCAL_RESOURCE_CEILING.memory_gib:
        raise ExecutionPreflightError("Execution budget must be at least 8 CPUs and 12 GiB to satisfy enabled local process contracts.")
    if capacity.logical_cpus is not None and cpus > capacity.logical_cpus:
        raise ExecutionPreflightError(f"Requested {cpus} CPUs exceeds detected host capacity of {capacity.logical_cpus}.")
    memory_capacity = min(capacity.total_memory_gib or memory_gb, capacity.available_memory_gib or memory_gb)
    if capacity.total_memory_gib is not None and memory_gb > memory_capacity:
        raise ExecutionPreflightError(f"Requested {memory_gb} GiB exceeds detected available memory of {memory_capacity} GiB.")


def project_execution_budget(config: Any) -> ResourceContract:
    return ResourceContract("PROJECT_LOCAL", config.execution.max_cpus, config.execution.max_memory_gb, LOCAL_RESOURCE_CEILING.time_hours)


@dataclass(frozen=True)
class RuntimeSnapshot:
    host_os: str
    host_architecture: str
    logical_cpus: int | None
    host_memory_bytes: int | None
    docker_architecture: str | None
    docker_memory_bytes: int | None
    docker_version: str | None
    control_plane_image_architecture: str | None
    docker_cpus: int | None = None


@dataclass(frozen=True)
class PreparedRun:
    report: ValidationReport
    nextflow_version: str
    container_runtime: str


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


def resolve_execution_workspace(case_id: str, run_id: str) -> ExecutionWorkspace:
    """Return the local execution paths for one logical case/run without writing them.

    Nextflow creates ``.nextflow/cache`` below its current working directory.  That
    LevelDB state is operational scratch data, not a durable project artifact, so
    it must never be created on a client or shared project filesystem.
    """

    configured = os.environ.get(EXECUTION_ROOT_ENV)
    if configured:
        base = Path(configured).expanduser()
    elif sys.platform == "darwin":
        base = Path.home() / "Library" / "Caches" / "rnaseq-control-plane"
    else:
        base = Path(os.environ.get("XDG_CACHE_HOME", str(Path.home() / ".cache"))) / "rnaseq-control-plane"
    if not base.is_absolute():
        raise ExecutionPreflightError(f"{EXECUTION_ROOT_ENV} must be an absolute path when configured.")
    root = (base / case_id / run_id).resolve()
    return ExecutionWorkspace(root=root, launch_dir=root / "launch", work_dir=root / "work")


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


def inspect_container_image(image: str) -> dict[str, object]:
    """Return Docker's observed immutable identity for an already-resolved image.

    A build tag is mutable and a locally built image has no registry digest, so
    provenance must preserve Docker's content ID as well as any observed
    RepoDigests.  This helper never pulls an image and is deliberately best
    effort: a failed inspection is recorded as unavailable rather than guessed.
    """

    observed: dict[str, object] = {
        "reference": image,
        "image_id": None,
        "repo_digests": [],
        "architecture": None,
    }
    try:
        result = _run_capture(["docker", "image", "inspect", image, "--format", "{{json .}}"])
        if result.returncode != 0:
            return observed
        payload = json.loads(result.stdout)
    except (FileNotFoundError, json.JSONDecodeError, TypeError, ValueError, AttributeError):
        return observed
    if not isinstance(payload, dict):
        return observed
    identifier = payload.get("Id")
    if isinstance(identifier, str) and identifier:
        observed["image_id"] = identifier
    digests = payload.get("RepoDigests")
    if isinstance(digests, list):
        observed["repo_digests"] = sorted(item for item in digests if isinstance(item, str) and item)
    architecture = _normalise_architecture(str(payload.get("Architecture") or ""))
    if architecture:
        observed["architecture"] = architecture
    return observed


def _normalise_architecture(value: str | None) -> str | None:
    if not value:
        return None
    normalized = value.strip().lower()
    aliases = {"aarch64": "arm64", "arm64/v8": "arm64", "x86_64": "amd64"}
    return aliases.get(normalized, normalized)


def _host_memory_bytes() -> int | None:
    try:
        return os.sysconf("SC_PAGE_SIZE") * os.sysconf("SC_PHYS_PAGES")
    except (AttributeError, OSError, ValueError):
        pass
    if sys.platform == "darwin":
        try:
            result = _run_capture(["sysctl", "-n", "hw.memsize"])
            return int(result.stdout.strip()) if result.returncode == 0 else None
        except (FileNotFoundError, ValueError):
            return None
    return None


def runtime_snapshot(image: str = CONTROL_PLANE_IMAGE) -> RuntimeSnapshot:
    """Collect cheap host/Docker facts without launching workflow containers."""

    host_architecture = _normalise_architecture(platform.machine()) or "unknown"
    docker_architecture: str | None = None
    docker_cpus: int | None = None
    docker_memory_bytes: int | None = None
    docker_version: str | None = None
    image_architecture: str | None = None
    try:
        info = _run_capture(["docker", "info", "--format", "{{json .}}"])
        if info.returncode == 0:
            payload = json.loads(info.stdout)
            if isinstance(payload, dict):
                docker_architecture = _normalise_architecture(str(payload.get("Architecture") or ""))
                memory = payload.get("MemTotal")
                docker_memory_bytes = memory if isinstance(memory, int) and memory > 0 else None
                cpus = payload.get("NCPU")
                docker_cpus = cpus if isinstance(cpus, int) and cpus > 0 else None
                version = payload.get("ServerVersion")
                docker_version = str(version) if version else None
            image_result = _run_capture(["docker", "image", "inspect", image, "--format", "{{json .}}"])
            if image_result.returncode == 0:
                image_payload = json.loads(image_result.stdout)
                if isinstance(image_payload, dict):
                    image_architecture = _normalise_architecture(str(image_payload.get("Architecture") or ""))
    except (FileNotFoundError, json.JSONDecodeError, TypeError):
        pass
    return RuntimeSnapshot(
        host_os=platform.system() or "unknown",
        host_architecture=host_architecture,
        logical_cpus=os.cpu_count(),
        host_memory_bytes=_host_memory_bytes(),
        docker_architecture=docker_architecture,
        docker_memory_bytes=docker_memory_bytes,
        docker_version=docker_version,
        control_plane_image_architecture=image_architecture,
        docker_cpus=docker_cpus,
    )


def downstream_docker_user_mapping(
    *, host_os: str | None = None, uid: int | None = None, gid: int | None = None,
) -> str | None:
    """Return the invoking Linux user's Docker identity for downstream tasks.

    Nextflow creates each task work directory on the host before Docker starts.
    On Linux (including WSL), Docker bind mounts preserve numeric ownership, so
    the image's fixed mamba user cannot necessarily create ``.command.*`` task
    files.  Docker Desktop on macOS already virtualizes shared-file ownership;
    preserving its default container user avoids changing the validated macOS
    path.  Values are resolved locally and constrained to non-negative integer
    IDs before they are rendered into a Nextflow config.
    """

    if (host_os or platform.system()).strip().lower() != "linux":
        return None
    try:
        resolved_uid = os.getuid() if uid is None else uid
        resolved_gid = os.getgid() if gid is None else gid
    except AttributeError:
        return None
    if type(resolved_uid) is not int or type(resolved_gid) is not int:
        return None
    if resolved_uid < 0 or resolved_gid < 0:
        return None
    return f"{resolved_uid}:{resolved_gid}"


def downstream_docker_user_mapping_check() -> RuntimeCheck:
    """Explain the downstream task-user policy in ``rnaseq doctor`` output."""

    mapping = downstream_docker_user_mapping()
    if mapping is not None:
        return RuntimeCheck(
            "Downstream Docker user mapping", "FOUND",
            f"Linux/WSL downstream tasks will use --user {mapping} for host-mounted Nextflow work directories.",
        )
    return RuntimeCheck(
        "Downstream Docker user mapping", "FOUND",
        "Not applied outside Linux/WSL; Docker Desktop macOS shared-file behavior keeps the image default user.",
    )


def _gib(value: int | None) -> str:
    return "unavailable" if value is None else f"{value / (1024 ** 3):.1f} GiB"


def runtime_resource_checks(snapshot: RuntimeSnapshot, budget: ResourceContract = LOCAL_RESOURCE_CEILING) -> tuple[RuntimeCheck, ...]:
    """Turn portable runtime facts into doctor PASS/WARN/FAIL checks."""

    host = RuntimeCheck(
        "Host runtime", "FOUND",
        f"OS={snapshot.host_os}; architecture={snapshot.host_architecture}; logical_cpus={snapshot.logical_cpus or 'unavailable'}; memory={_gib(snapshot.host_memory_bytes)}",
    )
    if snapshot.docker_architecture is None:
        docker = RuntimeCheck("Docker runtime", "NOT FOUND", "Docker runtime details are unavailable.", "FAIL")
    else:
        docker = RuntimeCheck(
            "Docker runtime", "FOUND",
            f"architecture={snapshot.docker_architecture}; logical_cpus={snapshot.docker_cpus or 'unavailable'}; memory={_gib(snapshot.docker_memory_bytes)}; version={snapshot.docker_version or 'unavailable'}",
        )
    checks: list[RuntimeCheck] = [host, docker]
    if snapshot.control_plane_image_architecture is None:
        checks.append(RuntimeCheck("Control-plane image architecture", "NOT FOUND", "Image architecture is unavailable; build or inspect rnaseq-control-plane:latest.", "WARN"))
    elif snapshot.host_architecture == "arm64" and snapshot.control_plane_image_architecture == "amd64":
        checks.append(RuntimeCheck("Control-plane image architecture", "FOUND", "amd64 image on arm64 host; Docker/Rosetta emulation may reduce throughput.", "WARN"))
    else:
        checks.append(RuntimeCheck("Control-plane image architecture", "FOUND", f"image={snapshot.control_plane_image_architecture}; host={snapshot.host_architecture}"))
    required_memory = budget.memory_gib * 1024 ** 3
    requested = f"{budget.cpus} CPUs / {_gib(required_memory)}"
    unavailable = snapshot.docker_memory_bytes is None or snapshot.docker_cpus is None
    exceeds = (
        (snapshot.docker_memory_bytes is not None and snapshot.docker_memory_bytes < required_memory)
        or (snapshot.docker_cpus is not None and snapshot.docker_cpus < budget.cpus)
    )
    if unavailable:
        checks.append(RuntimeCheck("Selected local ceiling", "NOT FOUND", f"Requested local ceiling is {requested}; Docker CPU or memory allocation is unavailable, so capacity cannot be confirmed.", "WARN"))
    elif exceeds:
        checks.append(RuntimeCheck("Selected local ceiling", "NOT FOUND", f"Requested local ceiling is {requested}; Docker exposes {snapshot.docker_cpus} CPUs / {_gib(snapshot.docker_memory_bytes)}. Reduce the requested profile or increase Docker allocation.", "WARN"))
    else:
        checks.append(RuntimeCheck("Selected local ceiling", "FOUND", f"Requested aggregate ceiling={budget.cpus} CPUs/{budget.memory_gib} GiB; per-process declarations and maxForks remain in effect; one project at a time."))
    checks.append(RuntimeCheck("nf-core upstream image architecture", "FOUND", "nf-core/rnaseq resolves process images dynamically; inspect the frozen Nextflow trace for per-process image architecture.", "WARN"))
    return tuple(checks)


def render_local_resource_config(budget: ResourceContract = LOCAL_RESOURCE_CEILING) -> str:
    """Render the frozen, auditable local nf-core resource policy."""

    small, medium, large = (RESOURCE_CONTRACTS[name] for name in ("SMALL", "MEDIUM", "LARGE"))
    return (
        f"// Local research resource contract: SMALL=1/2 GiB, MEDIUM=4/8 GiB, LARGE=8/12 GiB; aggregate ceiling={budget.cpus} CPUs/{budget.memory_gib} GiB; one project at a time.\n"
        "// The local executor treats these as the aggregate per-run budget.\n"
        f"executor {{ cpus = {budget.cpus}; memory = '{budget.memory_gib}.GB' }}\n"
        "process {\n"
        f"  resourceLimits = [cpus: {budget.cpus}, memory: '{budget.memory_gib}.GB', time: '{budget.time_hours}.h']\n"
        f"  withLabel:process_low {{ cpus = {small.cpus}; memory = '{small.memory_gib}.GB'; time = '{small.time_hours}.h' }}\n"
        f"  withLabel:process_medium {{ cpus = {medium.cpus}; memory = '{medium.memory_gib}.GB'; time = '{medium.time_hours}.h' }}\n"
        f"  withLabel:process_high {{ cpus = {large.cpus}; memory = '{large.memory_gib}.GB'; time = '{large.time_hours}.h' }}\n"
        "  // One 8 GiB Salmon task at a time avoids local Docker oversubscription.\n"
        f"  withName: '.*:SALMON_QUANT' {{ cpus = {medium.cpus}; memory = '{medium.memory_gib}.GB'; time = '{medium.time_hours}.h'; maxForks = 1 }}\n"
        "  // First-party workflow process names are deliberately explicit and match its task directives.\n"
        "  withName: '.*HISAT2_ALIGN' { maxForks = 1 }\n"
        "  withName: '.*(SORT_LANE_BAM|MERGE_AND_INDEX|PREPARE_COUNT_BAM)' { maxForks = 2 }\n"
        "  withName: '.*FEATURECOUNTS' { maxForks = 2 }\n"
        "  withName: '.*(ASSEMBLE_COUNTS|MULTIQC|L1_ANALYSIS|L2_ANALYSIS|ENRICHMENT_ANALYSIS|TECHNICAL_REPORT.*)' { maxForks = 1 }\n"
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
            "suggestion=check Docker memory allocation and the frozen local resource contract; do not change scientific parameters automatically."
        )
    if "no space left on device" in stderr.lower():
        return f"LIKELY_DISK_EXHAUSTION: process={stage}; exit_code={returncode}; suggestion=free space in the local execution workspace or Docker runtime."
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


def check_docker() -> RuntimeCheck:
    """Confirm that Docker's daemon is usable, not just that its CLI exists."""

    try:
        result = _run_capture(["docker", "info"])
    except FileNotFoundError:
        return RuntimeCheck("Docker", "NOT FOUND", "Docker executable was not found on PATH.")
    if result.returncode != 0:
        detail = (result.stderr or result.stdout).strip() or "docker info failed."
        return RuntimeCheck("Docker", "NOT FOUND", detail)
    return RuntimeCheck("Docker", "FOUND", "Docker daemon is available.")


def check_container_runtime(image: str = CONTROL_PLANE_IMAGE) -> RuntimeCheck:
    """Verify the built control-plane image has Nextflow task prerequisites.

    This intentionally runs only a short shell/R package probe; it does not run a
    workflow, access project inputs, or pull an image implicitly.
    """

    docker = check_docker()
    if docker.state != "FOUND":
        return RuntimeCheck("Control-plane container", "NOT FOUND", "Docker daemon is unavailable.")
    try:
        present = _run_capture(["docker", "image", "inspect", image])
    except FileNotFoundError:
        return RuntimeCheck("Control-plane container", "NOT FOUND", "Docker executable was not found on PATH.")
    if present.returncode != 0:
        return RuntimeCheck(
            "Control-plane container", "NOT FOUND",
            f"Required image {image} is not available locally; build or resolve it before execution.",
        )
    packages = ", ".join(repr(package) for package in CONTAINER_R_PACKAGES)
    probe = (
        "for executable in ps python Rscript; do "
        "command -v \"$executable\" >/dev/null || { echo \"missing executable: $executable\" >&2; exit 1; }; "
        "done; "
        f"Rscript -e \"packages <- c({packages}); missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly=TRUE)]; if (length(missing)) {{ cat('missing R package(s): ', paste(missing, collapse=', '), '\\n', file=stderr()); quit(status=1) }}\"; "
        "python -m rnaseq.workflow_support report --help | grep -F -- '--enrichment' >/dev/null "
        "|| { echo 'missing report CLI option: --enrichment' >&2; exit 1; }"
    )
    result = _run_capture([
        # Keep the image entrypoint and use a non-login shell so the probe sees
        # the same activated micromamba PATH as a Nextflow task container.
        "docker", "run", "--rm", image, "sh", "-c", probe,
    ])
    if result.returncode != 0:
        diagnostic_parts = [
            f"stdout: {result.stdout.strip()}" for result.stdout in (result.stdout,) if result.stdout.strip()
        ] + [
            f"stderr: {result.stderr.strip()}" for result.stderr in (result.stderr,) if result.stderr.strip()
        ]
        output = "; ".join(diagnostic_parts)
        detail = (
            f"container prerequisite probe exited {result.returncode}: {output}"
            if output
            else f"container prerequisite probe exited {result.returncode} without diagnostic output."
        )
        return RuntimeCheck("Control-plane container", "NOT FOUND", detail)
    return RuntimeCheck(
        "Control-plane container", "FOUND",
        f"requested={image}; ps, python, Rscript, required R packages, and the final-report CLI contract are available.",
    )


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


def doctor_checks(project_dir: Path | None = None) -> tuple[RuntimeCheck, ...]:
    """Return inexpensive, non-mutating execution prerequisite checks."""

    probe = Path.cwd()
    writable = probe.exists() and probe.is_dir() and probe.stat().st_mode != 0
    from rnaseq.downstream import r_runtime_checks
    requested_image = CONTROL_PLANE_IMAGE
    budget = LOCAL_RESOURCE_CEILING
    if project_dir is not None:
        try:
            from rnaseq.project import load_project
            config = load_project(project_dir).config
            requested_image = config.runtime.control_plane_image
            budget = project_execution_budget(config)
        except (OSError, ValueError):
            pass
    observed_image = inspect_container_image(requested_image)
    snapshot = runtime_snapshot(requested_image)

    workspace = resolve_execution_workspace("doctor", "resource-check").work_dir
    disk_probe = workspace
    while not disk_probe.exists() and disk_probe != disk_probe.parent:
        disk_probe = disk_probe.parent
    try:
        free = shutil.disk_usage(disk_probe).free
        disk_check = RuntimeCheck("Execution work-directory free space", "FOUND", f"path={workspace}; available_at={disk_probe}; free={_gib(free)}")
    except OSError as exc:
        disk_check = RuntimeCheck("Execution work-directory free space", "NOT FOUND", f"path={workspace}; unable to inspect free space: {exc}", "WARN")

    return (
        RuntimeCheck("Python", "FOUND", "Python runtime is active."),
        check_nextflow(),
        check_docker(),
        check_container_runtime(requested_image),
        RuntimeCheck(
            "Control-plane image identity",
            "FOUND" if observed_image.get("image_id") else "NOT FOUND",
            f"requested={requested_image}; observed_image_id={observed_image.get('image_id') or 'unavailable'}; "
            f"observed_repo_digests={observed_image.get('repo_digests') or []}",
            None if observed_image.get("image_id") else "WARN",
        ),
        downstream_docker_user_mapping_check(),
        *runtime_resource_checks(snapshot, budget),
        _reference_runtime_check(project_dir),
        RuntimeCheck("Disk write access", "FOUND" if writable else "NOT FOUND", str(probe)),
        disk_check,
        *r_runtime_checks(),
    )


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
    _validate_custom_reference_files(report)
    require_fresh_plan(report)
    # Also prove the generated samplesheet remains renderable before any run directory exists.
    render_samplesheet(report)
    nextflow = check_nextflow()
    if nextflow.state != "FOUND":
        raise ExecutionPreflightError(
            "Nextflow is required for FASTQ execution but was not found. " + nextflow.detail
        )
    docker = check_docker()
    if docker.state != "FOUND":
        raise ExecutionPreflightError(
            "Docker is required for the verified local execution profile but is unavailable. "
            + docker.detail
        )
    return PreparedRun(report, nextflow.detail, "docker")


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
    _write_text(runtime_config, render_local_resource_config(project_execution_budget(report.config)))
    return samplesheet, upstream_path, params_path, runtime_config


def build_nextflow_command(
    report: ValidationReport, *, samplesheet: Path, output_dir: Path, profile: str,
    params_file: Path | None = None, config_file: Path | None = None,
    reference_paths: dict[str, Path] | None = None, work_dir: Path | None = None,
) -> list[str]:
    """Build the exact argument vector; it is never run through a shell."""

    _require_fastq_execution_report(report)
    if profile != LOCAL_PROFILE:
        raise ExecutionPreflightError("Only the local profile is executable in Milestone 2.")
    assert report.config is not None
    reference: ReferenceConfig = report.config.reference
    command = ["nextflow", "run"]
    if config_file is not None:
        command.extend(["-c", str(config_file.resolve())])
    command.extend([
        "nf-core/rnaseq", "-r", NFCORE_RNASEQ_VERSION,
        "-profile", CONTAINER_PROFILE,
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
    config_file: Path | None = None,
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
        use_runtime_splices = local.hisat2_strategy == "genome_only_runtime_splices"
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
    command = ["nextflow", "run"]
    if config_file is not None:
        command.extend(["-c", str(config_file.resolve())])
    command.extend([
        str(HISAT2_WORKFLOW), "-profile", CONTAINER_PROFILE,
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
        work_dir=workspace.work_dir / "upstream",
    )
    requested_image = prepared.report.config.runtime.control_plane_image
    runtime = runtime_snapshot(requested_image)
    source_root = Path(__file__).resolve().parents[2]
    git_commit: str | None = None
    try:
        git_result = _run_capture(["git", "-C", str(source_root), "rev-parse", "HEAD"])
        if git_result.returncode == 0 and isinstance(getattr(git_result, "stdout", None), str):
            git_commit = git_result.stdout.strip() or None
    except (OSError, ValueError, AttributeError):
        pass
    provenance = {
        "nextflow_version": prepared.nextflow_version,
        "nfcore_rnaseq_version": NFCORE_RNASEQ_VERSION,
        "execution_profile": LOCAL_PROFILE,
        "container_runtime": prepared.container_runtime,
        "container_image": inspect_container_image(requested_image),
        "production_intended": prepared.report.config.reference.acceptance == "production",
        "git_commit": git_commit,
        "source_checkout": str(source_root),
        "workflow_sha256": {
            "workflow/main.nf": sha256_file(source_root / "workflow" / "main.nf"),
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
            "selected_local_ceiling": {
                "cpus": LOCAL_RESOURCE_CEILING.cpus,
                "memory_gib": LOCAL_RESOURCE_CEILING.memory_gib,
                "one_project_at_a_time": True,
            },
            "host_architecture": runtime.host_architecture,
            "docker_architecture": runtime.docker_architecture,
            "docker_memory_bytes": runtime.docker_memory_bytes,
            "control_plane_image_architecture": runtime.control_plane_image_architecture,
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
