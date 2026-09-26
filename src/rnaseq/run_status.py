"""Read-only progress view of nf-rna case runs for ``rnaseq status``.

Everything here reads durable artifacts only: ``run_state.json``, the run's
own ``logs/`` and Nextflow trace files, and frozen provenance.  Nothing is
written, and Nextflow does not need to be running.  A run is reported as
SUCCESS only when its durable state says so; a vanished process is never
taken as success.
"""

from __future__ import annotations

import csv
import json
import os
import re
import socket
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path
from typing import Any

import yaml

RUN_ID_PATTERN = re.compile(r"^[0-9]{8}-[0-9]{6}\+[0-9]{4}(?:-[0-9]{2,3})?$")
PROCESS_RECORD = Path("logs") / "rnaseq.process.json"
EXECUTION_LOG = Path("logs") / "rnaseq.log"
PHASES = ("validation", "upstream", "downstream", "delivery", "completed")
ACTIVE_STATES = {"CREATED", "RUNNING"}
_SUBMITTED = re.compile(r"^\[([0-9a-f]{2}/[0-9a-f]{6})\] (Submitted|Cached) process > (.+?)\s*$")
_TRACE_COLUMNS = ("hash", "name", "status")


# --------------------------------------------------------------------- process identity

def _boot_id() -> str | None:
    try:
        return Path("/proc/sys/kernel/random/boot_id").read_text(encoding="utf-8").strip() or None
    except OSError:
        return None


def _process_start_ticks(pid: int) -> int | None:
    """Kernel start time of ``pid`` (clock ticks since boot); guards against PID reuse."""

    try:
        stat = Path(f"/proc/{pid}/stat").read_text(encoding="utf-8")
    except OSError:
        return None
    try:
        # Field 22; the command name (field 2) may contain spaces and ')' so split after the last ')'.
        return int(stat[stat.rindex(")") + 2:].split()[19])
    except (ValueError, IndexError):
        return None


def current_process_identity() -> dict[str, Any]:
    pid = os.getpid()
    return {"pid": pid, "hostname": socket.gethostname(), "boot_id": _boot_id(), "start_ticks": _process_start_ticks(pid)}


def process_liveness(record: dict[str, Any] | None) -> tuple[bool | None, str]:
    """Return (alive, explanation); ``None`` means it cannot be determined here."""

    if not record or not isinstance(record.get("pid"), int):
        return None, "no process identity was recorded (run predates per-run process records)"
    pid = record["pid"]
    host = record.get("hostname")
    if host and host != socket.gethostname():
        return None, f"launched on host {host}; liveness cannot be checked from {socket.gethostname()}"
    boot = record.get("boot_id")
    if boot and _boot_id() and boot != _boot_id():
        return False, f"rnaseq process {pid} belonged to an earlier boot of this machine"
    try:
        os.kill(pid, 0)
    except ProcessLookupError:
        return False, f"rnaseq process {pid} no longer exists"
    except PermissionError:
        pass  # exists, owned by someone else
    except OSError:
        return None, f"cannot check rnaseq process {pid}"
    expected = record.get("start_ticks")
    if isinstance(expected, int):
        observed = _process_start_ticks(pid)
        if observed is not None and observed != expected:
            return False, f"rnaseq process {pid} no longer exists (PID reused by another process)"
    return True, f"rnaseq process {pid} is alive"


# --------------------------------------------------------------------- cached reads

@dataclass
class _Cache:
    """(path -> (size, mtime_ns, value)); lets --watch skip unchanged files."""

    entries: dict[Path, tuple[int, int, Any]] = field(default_factory=dict)

    def load(self, path: Path, parser: Any) -> Any:
        try:
            status = path.stat()
        except OSError:
            self.entries.pop(path, None)
            return None
        key = (status.st_size, status.st_mtime_ns)
        cached = self.entries.get(path)
        if cached is not None and cached[:2] == key:
            return cached[2]
        try:
            value = parser(path)
        except (OSError, UnicodeError, ValueError, yaml.YAMLError, csv.Error):
            value = None
        self.entries[path] = (*key, value)
        return value


def _json(path: Path) -> Any:
    return json.loads(path.read_text(encoding="utf-8"))


def _yaml(path: Path) -> Any:
    return yaml.safe_load(path.read_text(encoding="utf-8"))


def _trace_rows(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", errors="replace", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None or not set(_TRACE_COLUMNS) <= set(reader.fieldnames):
            return []
        return [{key: row.get(key) or "" for key in _TRACE_COLUMNS} for row in reader]


def _submissions(path: Path) -> dict[str, str]:
    """hash -> task name from Nextflow's plain log ('Submitted process > NAME (tag)')."""

    found: dict[str, str] = {}
    with path.open(encoding="utf-8", errors="replace") as handle:
        for line in handle:
            match = _SUBMITTED.match(line)
            if match:
                found[match.group(1)] = match.group(3)
    return found


# --------------------------------------------------------------------- model

@dataclass(frozen=True)
class TaskProgress:
    completed: int = 0
    cached: int = 0
    failed: int = 0
    running: tuple[str, ...] = ()
    failed_names: tuple[str, ...] = ()
    trace_available: bool = False
    submitted: int = 0

    @property
    def done(self) -> int:
        return self.completed + self.cached


@dataclass(frozen=True)
class StageStatus:
    state: str
    tasks: TaskProgress | None = None
    note: str | None = None


@dataclass(frozen=True)
class RunStatus:
    case_id: str
    run_id: str
    run_dir: Path
    recorded_state: str
    state: str
    phase: str
    started_at: str | None
    completed_at: str | None
    upstream: StageStatus
    downstream: StageStatus
    delivery: StageStatus
    liveness: str | None = None
    error: str | None = None
    resources: dict[str, Any] | None = None
    last_activity: float | None = None
    logs: tuple[Path, ...] = ()
    attempt: str | None = None

    @property
    def is_final(self) -> bool:
        return self.state in {"SUCCESS", "FAILED", "INTERRUPTED"}


# --------------------------------------------------------------------- discovery

def discover_runs(project_dir: Path) -> list[tuple[Path, dict[str, Any]]]:
    """All recorded runs (case/run layout and legacy flat runs), newest first."""

    runs_dir = project_dir / "runs"
    if not runs_dir.is_dir():
        return []
    found: list[tuple[Path, dict[str, Any]]] = []
    for path in set(runs_dir.glob("*/run_state.json")) | set(runs_dir.glob("*/*/run_state.json")):
        try:
            state = _json(path)
        except (OSError, UnicodeError, ValueError):
            continue
        if isinstance(state, dict):
            found.append((path.parent, state))

    def key(item: tuple[Path, dict[str, Any]]) -> tuple[float, float]:
        run_dir, state = item
        started = _parse_time(state.get("started_at"))
        try:
            mtime = (run_dir / "run_state.json").stat().st_mtime
        except OSError:
            mtime = 0.0
        stamp = started.timestamp() if started and started.tzinfo else (started.replace().timestamp() if started else mtime)
        return (stamp, mtime)

    return sorted(found, key=key, reverse=True)


def select_run(project_dir: Path, case_id: str | None = None, run_id: str | None = None) -> tuple[Path, dict[str, Any]] | None:
    for run_dir, state in discover_runs(project_dir):
        if case_id is not None and state.get("case_id", run_dir.parent.name) != case_id:
            continue
        if run_id is not None and state.get("run_id", run_dir.name) != run_id:
            continue
        return run_dir, state
    return None


# --------------------------------------------------------------------- derivation

def _parse_time(value: object) -> datetime | None:
    if not isinstance(value, str) or not value:
        return None
    try:
        return datetime.fromisoformat(value.replace("Z", "+00:00"))
    except ValueError:
        return None


def _phase_index(phase: str) -> int:
    return PHASES.index(phase) if phase in PHASES else 0


def _recorded_phase(state: dict[str, Any]) -> str:
    phase = state.get("phase")
    if phase in {"upstream", "downstream", "delivery"}:
        return str(phase)
    return "validation"  # CREATED, freeze, retry_freeze, or unknown


def _tasks(trace: list[dict[str, str]] | None, submitted: dict[str, str] | None, *, running_possible: bool) -> TaskProgress | None:
    if trace is None and submitted is None:
        return None
    rows = trace or []
    completed = sum(1 for row in rows if row["status"] == "COMPLETED")
    cached = sum(1 for row in rows if row["status"] == "CACHED")
    failed_rows = [row for row in rows if row["status"] in {"FAILED", "ABORTED"}]
    running: tuple[str, ...] = ()
    if running_possible and submitted:
        finished = {row["hash"] for row in rows}
        running = tuple(name for task_hash, name in submitted.items() if task_hash not in finished)
    return TaskProgress(
        completed, cached, len(failed_rows), running, tuple(row["name"] for row in failed_rows), trace is not None,
        len(submitted or {}),
    )


def _latest(paths: list[Path]) -> Path | None:
    candidates = [path for path in paths if path.is_file()]
    return max(candidates, key=lambda path: (path.name, path.stat().st_mtime)) if candidates else None


def _upstream_trace(run_dir: Path) -> Path | None:
    nfcore = sorted((run_dir / "upstream" / "nfcore_rnaseq" / "pipeline_info").glob("execution_trace_*.txt"))
    return _latest([*nfcore, run_dir / "provenance" / "upstream.trace.txt"])


def _stage_state(stage: str, phase: str, overall: str) -> str:
    """PENDING / RUNNING / SUCCESS / FAILED / INTERRUPTED / NOT STARTED for one stage.

    A stage is SUCCESS only when the durable state moved past it (or the whole
    run recorded SUCCESS); the stage the run stopped in inherits its outcome.
    """

    if overall == "SUCCESS":
        return "SUCCESS"
    order, current = _phase_index(stage), _phase_index(phase)
    if current > order:
        return "SUCCESS"
    if current < order:
        return "PENDING" if overall in {"RUNNING", "PLANNED"} else "NOT STARTED"
    if overall in {"FAILED", "INTERRUPTED", "RUNNING"}:
        return overall
    return "PENDING"


def build_status(run_dir: Path, state: dict[str, Any], cache: _Cache | None = None) -> RunStatus:
    cache = cache or _Cache()
    recorded = str(state.get("status") or "UNKNOWN")
    phase = _recorded_phase(state)
    process = cache.load(run_dir / PROCESS_RECORD, _json)
    liveness_note: str | None = None
    overall = recorded
    if recorded in ACTIVE_STATES:
        alive, liveness_note = process_liveness(process if isinstance(process, dict) else None)
        if alive is False:
            overall = "INTERRUPTED"
            liveness_note = f"stale: recorded {recorded}, but {liveness_note}; the run did not finish"
        elif alive is None:
            overall = "RUNNING" if recorded == "RUNNING" else "PLANNED"
            liveness_note = f"unverified: {liveness_note}"
        else:
            overall = "RUNNING"
    elif recorded not in {"SUCCESS", "FAILED", "INTERRUPTED"}:
        overall = "UNKNOWN"

    display_phase = "completed" if overall == "SUCCESS" else phase
    active = overall == "RUNNING"

    # Upstream
    frozen_manifest = cache.load(run_dir / "frozen" / "input_manifest.yaml", _yaml)
    input_type = (frozen_manifest or {}).get("input", {}).get("type") if isinstance(frozen_manifest, dict) else None
    contract = cache.load(run_dir / "frozen" / "downstream_contract.json", _json)
    source = contract.get("source") if isinstance(contract, dict) and isinstance(contract.get("source"), dict) else {}
    if input_type == "raw_counts" or source.get("type") == "raw_counts":
        upstream = StageStatus("NOT APPLICABLE", note="raw-count input")
    else:
        trace_path = _upstream_trace(run_dir)
        trace = cache.load(trace_path, _trace_rows) if trace_path else None
        submitted = cache.load(run_dir / "logs" / "upstream.stdout.log", _submissions)
        upstream_state = _stage_state("upstream", phase, overall)
        if source.get("reused_from") and upstream_state == "SUCCESS":
            upstream = StageStatus("SUCCESS", note=f"reused from {source['reused_from']}")
        else:
            upstream = StageStatus(upstream_state, _tasks(trace, submitted, running_possible=active and phase == "upstream"))

    # Downstream
    if state.get("downstream_skipped"):
        downstream = StageStatus("SKIPPED", note="technical QC only")
    else:
        trace_path = run_dir / "provenance" / "downstream.trace.txt"
        trace = cache.load(trace_path, _trace_rows) if trace_path.is_file() else None
        submitted = cache.load(run_dir / "logs" / "downstream.stdout.log", _submissions)
        downstream = StageStatus(
            _stage_state("downstream", phase, overall),
            _tasks(trace, submitted, running_possible=active and phase == "downstream"),
        )
    delivery = StageStatus(_stage_state("delivery", phase, overall))

    provenance = cache.load(run_dir / "provenance" / "run_provenance.yaml", _yaml)
    resources = provenance.get("runtime_resources") if isinstance(provenance, dict) else None
    launch_dir = provenance.get("execution_launch_dir") if isinstance(provenance, dict) else None
    if isinstance(provenance, dict) and isinstance(resources, dict):
        resources = {**resources, "_tuning_config": provenance.get("frozen_nfcore_tuning_config")}
    logs = tuple(path for path in (
        run_dir / EXECUTION_LOG, run_dir / "logs" / "upstream.stdout.log", run_dir / "logs" / "downstream.stdout.log",
    ) if path.is_file())
    activity = [run_dir / "run_state.json", *logs, run_dir / "provenance" / "downstream.trace.txt"]
    if isinstance(launch_dir, str) and launch_dir:
        activity.append(Path(launch_dir) / ".nextflow.log")  # Nextflow's own heartbeat; stat only
    trace_path = _upstream_trace(run_dir)
    if trace_path:
        activity.append(trace_path)
    mtimes = []
    for path in activity:
        try:
            mtimes.append(path.stat().st_mtime)
        except OSError:
            pass
    attempt = None
    retry = state.get("retry_of")
    if isinstance(retry, dict):
        attempt = f"retry of {retry.get('case_id', '?')}/{retry.get('run_id', '?')}"
    return RunStatus(
        case_id=str(state.get("case_id") or run_dir.parent.name),
        run_id=str(state.get("run_id") or run_dir.name),
        run_dir=run_dir,
        recorded_state=recorded,
        state=overall,
        phase=display_phase,
        started_at=state.get("started_at"),
        completed_at=state.get("completed_at"),
        upstream=upstream,
        downstream=downstream,
        delivery=delivery,
        liveness=liveness_note,
        error=state.get("error") if overall in {"FAILED", "INTERRUPTED"} else None,
        resources=resources if isinstance(resources, dict) else None,
        last_activity=max(mtimes) if mtimes else None,
        logs=logs,
        attempt=attempt,
    )


# --------------------------------------------------------------------- rendering

def _format_time(value: str | None) -> str:
    moment = _parse_time(value)
    return moment.strftime("%Y-%m-%d %H:%M:%S") if moment else "-"


def format_duration(seconds: float) -> str:
    seconds = max(0, int(seconds))
    days, remainder = divmod(seconds, 86400)
    hours, remainder = divmod(remainder, 3600)
    minutes, secs = divmod(remainder, 60)
    if days:
        return f"{days}d {hours:02d}h {minutes:02d}m"
    if hours:
        return f"{hours}h {minutes:02d}m"
    return f"{minutes}m {secs:02d}s"


def _elapsed(status: RunStatus, now: datetime) -> str:
    started = _parse_time(status.started_at)
    if started is None:
        return "-"
    finished = _parse_time(status.completed_at)
    if finished is None:
        if status.state != "RUNNING":
            return "-"
        finished = now.astimezone(started.tzinfo) if started.tzinfo else now
    return format_duration((finished - started).total_seconds())


def _task_line(stage: StageStatus) -> str | None:
    tasks = stage.tasks
    if tasks is None:
        return None
    if not tasks.trace_available:
        running = f", {len(tasks.running)} running" if tasks.running else ""
        return f"{tasks.submitted} submitted{running}; no completed-task trace yet"
    if stage.state == "SUCCESS":
        return f"{tasks.done}/{tasks.done} completed, {tasks.cached} cached, {tasks.failed} failed"
    running = f", {len(tasks.running)} running" if tasks.running else ""
    return f"{tasks.done} completed ({tasks.cached} cached){running}, {tasks.failed} failed; total not yet known"


def _running_lines(stage: StageStatus, limit: int = 8) -> list[str]:
    if stage.tasks is None or not stage.tasks.running:
        return []
    names = list(stage.tasks.running)
    shown = names[:limit]
    lines = [f"  - {name}" for name in shown]
    if len(names) > limit:
        lines.append(f"  - ... {len(names) - limit} more")
    return lines


def _resource_lines(resources: dict[str, Any]) -> list[str]:
    effective = resources.get("effective") if isinstance(resources.get("effective"), dict) else None
    if not effective:
        return []
    policy = resources.get("policy") if isinstance(resources.get("policy"), dict) else None
    lines = [f"Resources:  {effective.get('cpus')} CPUs / {effective.get('memory_gib')} GiB (Nextflow local executor ceiling)"]
    if policy:
        modes = {policy.get("cpu_mode"), policy.get("memory_mode")}
        mode = "auto" if modes == {"auto"} else "explicit" if modes == {"explicit"} else "mixed auto/explicit"
        detected = policy.get("detected") if isinstance(policy.get("detected"), dict) else {}
        reserve = policy.get("os_reserve") if isinstance(policy.get("os_reserve"), dict) else {}
        detail = f"Policy:     {mode}; usable {detected.get('usable_cpus', '?')} CPUs / {detected.get('usable_memory_gib', '?')} GiB"
        if reserve.get("cpus") is not None or reserve.get("memory_gib") is not None:
            detail += f"; OS reserve {reserve.get('cpus') or 0} CPUs / {reserve.get('memory_gib') or 0} GiB"
        lines.append(detail)
        for item in policy.get("process_tuning") or []:
            if isinstance(item, dict):
                lines.append(f"Tuning:     {item.get('selector')} memory {item.get('memory_gib')} GiB/task; CPUs unchanged")
    else:
        requested = resources.get("requested") if isinstance(resources.get("requested"), dict) else {}
        lines.append(f"Requested:  {requested.get('cpus', '?')} CPUs / {requested.get('memory_gib', '?')} GiB (no policy recorded; pre-policy run)")
    if resources.get("clamped"):
        lines.append("Note:       project limits were clamped to this machine's capacity")
    return lines


def render_status(status: RunStatus, *, now: datetime | None = None) -> str:
    now = now or datetime.now().astimezone()
    lines = [
        f"Case:       {status.case_id}",
        f"Run:        {status.run_id}" + (f" ({status.attempt})" if status.attempt else ""),
        f"State:      {status.state}",
        f"Phase:      {status.phase}",
    ]
    if status.liveness and (status.state != status.recorded_state or status.liveness.startswith(("unverified", "stale"))):
        lines.append(f"Note:       {status.liveness}")
    lines += ["", f"Started:    {_format_time(status.started_at)}"]
    if status.completed_at:
        lines.append(f"Finished:   {_format_time(status.completed_at)}")
    lines.append(f"Elapsed:    {_elapsed(status, now)}")
    if status.state == "RUNNING" and status.last_activity:
        lines.append(f"Activity:   last write {format_duration(now.timestamp() - status.last_activity)} ago")
    lines.append("")

    def stage_lines(label: str, stage: StageStatus) -> list[str]:
        out = [f"{label + ':':<12}{stage.state}" + (f" ({stage.note})" if stage.note else "")]
        task_line = _task_line(stage)
        if task_line:
            out.append(f"{'Tasks:':<12}{task_line}")
        if stage.state == "RUNNING":
            running = _running_lines(stage)
            if running:
                out.append("Running:")
                out.extend(running)
        if stage.tasks and stage.tasks.failed_names:
            out.append("Failed tasks:")
            out.extend(f"  - {name}" for name in stage.tasks.failed_names[:8])
        return out

    lines += stage_lines("Upstream", status.upstream)
    lines += stage_lines("Downstream", status.downstream)
    lines += stage_lines("Delivery", status.delivery)
    if status.error:
        first = str(status.error).strip().splitlines()[0] if str(status.error).strip() else ""
        lines += ["", f"Error:      {first[:300]}"]
    if status.resources:
        resource_lines = _resource_lines(status.resources)
        if resource_lines:
            lines += ["", *resource_lines]
    lines += ["", f"Run dir:    {status.run_dir}"]
    execution_log = status.run_dir / EXECUTION_LOG
    if execution_log in status.logs:
        lines.append(f"Log:        {execution_log}")
    elif status.logs:
        lines.append(f"Logs:       {execution_log.parent} (no per-run rnaseq.log; run predates it)")
    return "\n".join(lines)
