"""Read-only execution dashboard of nf-rna case runs for ``rnaseq status``.

Everything here reads durable artifacts only: ``run_state.json``, the run's own
``logs/``, the Nextflow trace files, frozen provenance, and the start/exit
markers Nextflow writes into each task's work directory.  Nothing is written,
and Nextflow does not need to be running.  A run is reported as SUCCESS only
when its durable state says so; a vanished process is never taken as success.

Status is observability only: it never reads or summarizes scientific results
(DESeq2 tables, enrichment, PCA); those belong to the report and delivery.
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
from typing import Any, Callable

import yaml

RUN_ID_PATTERN = re.compile(r"^[0-9]{8}-[0-9]{6}\+[0-9]{4}(?:-[0-9]{2,3})?$")
PROCESS_RECORD = Path("logs") / "rnaseq.process.json"
EXECUTION_LOG = Path("logs") / "rnaseq.log"
PHASES = ("validation", "upstream", "downstream", "delivery", "completed")
ACTIVE_STATES = {"CREATED", "RUNNING"}
FINAL_STATES = {"SUCCESS", "FAILED", "INTERRUPTED"}
_SUBMITTED = re.compile(r"^\[([0-9a-f]{2}/[0-9a-f]{6})\] (Submitted|Cached) process > (.+?)\s*$")
_TASK_NAME = re.compile(r"^(?P<path>[^ (]+)(?: \((?P<tag>.*)\))?$")
_TRACE_COLUMNS = ("hash", "name", "status")

# Stage states and their dashboard symbols.
DONE, RUNNING, WAITING, FAILED, INTERRUPTED, NOT_STARTED, SKIPPED = (
    "DONE", "RUNNING", "WAITING", "FAILED", "INTERRUPTED", "NOT_STARTED", "SKIPPED",
)
SYMBOLS = {DONE: "✓", RUNNING: "▶", WAITING: "○", FAILED: "!", INTERRUPTED: "⏸", NOT_STARTED: "–", SKIPPED: "–"}
ASCII_SYMBOLS = {DONE: "+", RUNNING: ">", WAITING: "o", FAILED: "!", INTERRUPTED: "=", NOT_STARTED: "-", SKIPPED: "-"}
_RULE = "━" * 44


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


def _samples(path: Path) -> dict[str, int]:
    """sample -> number of rows (FASTQ lanes) in the frozen samplesheet or metadata."""

    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        column = next((name for name in ("sample", "sample_id") if name in (reader.fieldnames or ())), None)
        lanes: dict[str, int] = {}
        for row in reader if column else ():
            if row.get(column):
                lanes[row[column]] = lanes.get(row[column], 0) + 1
        return lanes


# --------------------------------------------------------------------- model

@dataclass(frozen=True)
class Task:
    """One Nextflow task, from the trace (finished) or the plain log (submitted)."""

    hash: str
    name: str
    status: str  # COMPLETED CACHED FAILED ABORTED | RUNNING QUEUED SUBMITTED | STOPPED
    started: float | None = None

    @property
    def path(self) -> str:
        match = _TASK_NAME.match(self.name)
        return match.group("path") if match else self.name

    @property
    def process(self) -> str:
        return self.path.rsplit(":", 1)[-1]

    @property
    def tag(self) -> str | None:
        match = _TASK_NAME.match(self.name)
        return match.group("tag") if match else None

    @property
    def sample(self) -> str | None:
        tag = self.tag
        return tag.split()[-1] if tag else None

    @property
    def done(self) -> bool:
        return self.status in {"COMPLETED", "CACHED"}

    @property
    def active(self) -> bool:
        return self.status in {"RUNNING", "QUEUED", "SUBMITTED"}


@dataclass(frozen=True)
class StageSpec:
    """A user-facing pipeline stage and the Nextflow processes that make it up."""

    label: str
    match: Callable[[Task], bool]
    per_sample: tuple[str, ...] = ()  # process(es) whose completed tasks count finished samples
    final: str | None = None  # process whose completion marks the stage done
    expected: int | None = None  # known number of tasks, when the stage is not per-sample


@dataclass(frozen=True)
class Stage:
    label: str
    state: str
    detail: str = ""


@dataclass(frozen=True)
class TaskCounts:
    completed: int = 0
    cached: int = 0
    running: int = 0
    queued: int = 0
    failed: int = 0
    retried: int = 0
    evidence: bool = False


@dataclass(frozen=True)
class RunStatus:
    project: str
    case_id: str
    run_id: str
    run_dir: Path
    recorded_state: str
    state: str
    phase: str
    route: str
    started_at: str | None
    completed_at: str | None
    stages: tuple[Stage, ...]
    running: tuple[Task, ...] = ()
    failed: tuple[Task, ...] = ()
    tasks: TaskCounts = TaskCounts()
    liveness: str | None = None
    error: str | None = None
    resources: dict[str, Any] | None = None
    last_activity: float | None = None
    logs: tuple[Path, ...] = ()
    attempt: str | None = None
    delivery: Path | None = None

    @property
    def is_final(self) -> bool:
        return self.state in FINAL_STATES


# --------------------------------------------------------------------- discovery and selection

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


def _case_of(run_dir: Path, state: dict[str, Any]) -> str:
    if state.get("case_id"):
        return str(state["case_id"])
    return "legacy" if run_dir.parent.name == "runs" else run_dir.parent.name  # flat runs/RUN layout


def _project_of(run_dir: Path) -> str:
    runs = next((parent for parent in run_dir.parents if parent.name == "runs"), None)
    return runs.parent.name if runs is not None else run_dir.parent.name


def _attempt(state: dict[str, Any]) -> str | None:
    retry = state.get("retry_of")
    if not isinstance(retry, dict):
        return None
    source = f" (source {retry['status']})" if retry.get("status") else ""
    return f"retry of {retry.get('case_id', '?')}/{retry.get('run_id', '?')}{source}"


def _run_of(run_dir: Path, state: dict[str, Any]) -> str:
    return str(state.get("run_id") or run_dir.name)


def matching_runs(project_dir: Path, case_id: str | None = None, run_id: str | None = None) -> list[tuple[Path, dict[str, Any]]]:
    """Recorded runs filtered by exact case and/or run ID, newest first."""

    return [
        (run_dir, state) for run_dir, state in discover_runs(project_dir)
        if (case_id is None or _case_of(run_dir, state) == case_id) and (run_id is None or _run_of(run_dir, state) == run_id)
    ]


def select_run(project_dir: Path, case_id: str | None = None, run_id: str | None = None) -> tuple[Path, dict[str, Any]] | None:
    """The newest run matching the selection (the latest run when nothing is given)."""

    found = matching_runs(project_dir, case_id, run_id)
    return found[0] if found else None


def load_run(run_dir: Path, cache: _Cache) -> tuple[Path, dict[str, Any]] | None:
    """Re-read one already selected run; --watch stays pinned to it."""

    state = cache.load(run_dir / "run_state.json", _json)
    return (run_dir, state) if isinstance(state, dict) else None


# --------------------------------------------------------------------- derivation helpers

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


def _group_state(group: str, phase: str, overall: str) -> str:
    """Coarse state of one execution phase from the durable run state alone.

    A phase is DONE only when the durable state moved past it (or the whole run
    recorded SUCCESS); the phase the run stopped in inherits its outcome.
    """

    if overall == "SUCCESS":
        return DONE
    order, current = _phase_index(group), _phase_index(phase)
    if current > order:
        return DONE
    if current < order:
        return WAITING if overall in {"RUNNING", "PLANNED"} else NOT_STARTED
    if overall == "RUNNING":
        return RUNNING
    if overall == "FAILED":
        return FAILED
    if overall == "INTERRUPTED":
        return INTERRUPTED
    return WAITING


def effective_state(run_dir: Path, state: dict[str, Any], cache: _Cache) -> tuple[str, str | None]:
    """The state to display, reconciling a recorded RUNNING with process liveness."""

    recorded = str(state.get("status") or "UNKNOWN")
    if recorded in ACTIVE_STATES:
        process = cache.load(run_dir / PROCESS_RECORD, _json)
        alive, note = process_liveness(process if isinstance(process, dict) else None)
        if alive is False:
            return "INTERRUPTED", f"stale: recorded {recorded}, but {note}; the run did not finish"
        if alive is None:
            return ("RUNNING" if recorded == "RUNNING" else "PLANNED"), f"unverified: {note}"
        return "RUNNING", note
    if recorded in FINAL_STATES:
        return recorded, None
    return "UNKNOWN", None


def _work_dir(state: dict[str, Any], key: str) -> Path | None:
    command = state.get(key)
    if isinstance(command, list) and "-work-dir" in command:
        index = command.index("-work-dir")
        if index + 1 < len(command) and isinstance(command[index + 1], str):
            return Path(command[index + 1])
    return None


def _task_start(work_dir: Path | None, task_hash: str) -> tuple[str, float | None]:
    """RUNNING/QUEUED/SUBMITTED for a submitted task, from Nextflow's own task markers.

    ``.command.begin`` is touched when the task starts executing and
    ``.exitcode`` when it ends; without a resolvable work directory the task is
    only known to have been submitted.
    """

    if work_dir is None:
        return "SUBMITTED", None
    prefix, rest = task_hash.split("/", 1)
    try:
        candidates = [path for path in (work_dir / prefix).glob(f"{rest}*") if path.is_dir()]
    except OSError:
        candidates = []
    if len(candidates) != 1:
        return "SUBMITTED", None
    task_dir = candidates[0]
    try:
        begin = (task_dir / ".command.begin").stat().st_mtime
    except OSError:
        return "QUEUED", None
    if (task_dir / ".exitcode").exists():
        return "SUBMITTED", None  # finished; the trace row follows shortly
    return "RUNNING", begin


def _collect_tasks(
    trace: list[dict[str, str]] | None, submitted: dict[str, str] | None, *, live: bool, work_dir: Path | None,
) -> list[Task]:
    rows = trace or []
    # ABORTED tasks were killed by Nextflow because of another failure or an interrupt: stopped, not failed.
    tasks = [Task(row["hash"], row["name"], "STOPPED" if row["status"] == "ABORTED" else row["status"]) for row in rows]
    finished = {row["hash"] for row in rows}
    for task_hash, name in (submitted or {}).items():
        if task_hash in finished:
            continue
        if live:
            status, started = _task_start(work_dir, task_hash)
            tasks.append(Task(task_hash, name, status, started))
        else:
            tasks.append(Task(task_hash, name, "STOPPED"))
    return tasks


# HISAT2-route processes that run once per FASTQ lane and share the name 'PROCESS (sample)'.
PER_LANE_PROCESSES = frozenset({"FASTQC_RAW", "FASTP_PREPARE", "FASTQC_PROCESSED", "HISAT2_ALIGN", "SORT_LANE_BAM"})


def _expected_instances(task: Task, lanes: dict[str, int]) -> int:
    return lanes.get(task.sample or "", 1) if task.process in PER_LANE_PROCESSES else 1


def _unresolved_failures(tasks: list[Task], lanes: dict[str, int]) -> tuple[list[Task], int]:
    """(failed attempts not made good, failed attempts retried or being retried).

    Tasks are grouped by name. A failure is resolved when enough attempts of that
    name succeeded or are active again (a Nextflow retry in progress) to cover the
    expected instances: one, or one per lane for per-lane processes.
    """

    by_name: dict[str, list[Task]] = {}
    for task in tasks:
        by_name.setdefault(task.name, []).append(task)
    unresolved: list[Task] = []
    retried = 0
    for same in by_name.values():
        failed = [task for task in same if task.status == "FAILED"]
        if not failed:
            continue
        covered = sum(1 for task in same if task.done or task.active)
        shortfall = max(0, _expected_instances(same[0], lanes) - covered)
        unresolved += failed[:shortfall]
        retried += len(failed) - min(len(failed), shortfall)
    return unresolved, retried


# --------------------------------------------------------------------- route definitions

def _process_in(*names: str) -> Callable[[Task], bool]:
    wanted = set(names)
    return lambda task: task.process in wanted


def _path_contains(fragment: str) -> Callable[[Task], bool]:
    return lambda task: fragment in task.path


def _either(*predicates: Callable[[Task], bool]) -> Callable[[Task], bool]:
    return lambda task: any(predicate(task) for predicate in predicates)


_NFCORE_QC = "FASTQ_QC_TRIM_FILTER_SETSTRANDEDNESS"
_NFCORE_QC_PROCESSES = ("FASTQC", "TRIMGALORE", "FASTP", "FQ_LINT", "FQ_LINT_AFTER_TRIMMING", "CAT_FASTQ", "SORTMERNA", "BBMAP_BBSPLIT", "FQ_SUBSAMPLE", "UMITOOLS_EXTRACT")

HISAT2_STAGES = (
    StageSpec("FASTQ QC / preprocessing", _process_in("FASTQC_RAW", "FASTP_PREPARE", "FASTQC_PROCESSED"), per_sample=("FASTQC_PROCESSED",)),
    StageSpec("HISAT2 alignment", _process_in("HISAT2_ALIGN", "SORT_LANE_BAM", "MERGE_AND_INDEX"), per_sample=("MERGE_AND_INDEX",)),
    StageSpec("featureCounts", _process_in("PREPARE_COUNT_BAM", "FEATURECOUNTS", "ASSEMBLE_COUNTS"), per_sample=("FEATURECOUNTS",), final="ASSEMBLE_COUNTS"),
    StageSpec("MultiQC", _process_in("MULTIQC"), final="MULTIQC"),
)

SALMON_STAGES = (
    StageSpec("Reference preparation", _path_contains(":PREPARE_GENOME:")),
    StageSpec("FASTQ QC / trimming", _either(_path_contains(_NFCORE_QC), _process_in(*_NFCORE_QC_PROCESSES)), per_sample=("TRIMGALORE", "FASTP")),
    StageSpec("Salmon quantification", _process_in("SALMON_QUANT"), per_sample=("SALMON_QUANT",)),
    StageSpec("tximport / gene summary", _either(_path_contains("QUANT_TXIMPORT_SUMMARIZEDEXPERIMENT"), _process_in("CUSTOM_TX2GENE", "TXIMETA_TXIMPORT", "TXIMPORT"))),
    StageSpec("MultiQC", _process_in("MULTIQC", "DESEQ2_QC_PSEUDO", "DESEQ2_QC_STAR_SALMON"), final="MULTIQC"),
)

GENERIC_UPSTREAM = (StageSpec("Upstream", lambda task: True),)


def _downstream_stages(preset: str | None, enrichment: tuple[str, ...], backends: int) -> tuple[StageSpec, ...]:
    stages = [StageSpec("L1 expression QC", _process_in("L1_ANALYSIS"), final="L1_ANALYSIS")]
    if preset != "L1":
        stages.append(StageSpec("L2 / DESeq2", _process_in("L2_ANALYSIS"), final="L2_ANALYSIS"))
    if enrichment:
        stages.append(StageSpec(f"Enrichment ({', '.join(enrichment)})", _process_in("ENRICHMENT_ANALYSIS"), expected=backends or None))
    stages.append(StageSpec("Technical report", lambda task: task.process.startswith("TECHNICAL_REPORT")))
    return tuple(stages)


def _route(manifest: dict[str, Any]) -> tuple[str, str]:
    """(route key, human label) from the frozen input manifest."""

    input_type = (manifest.get("input") or {}).get("type") if isinstance(manifest.get("input"), dict) else None
    if input_type == "raw_counts":
        return "raw_counts", "raw counts → DESeq2"
    upstream = manifest.get("upstream") if isinstance(manifest.get("upstream"), dict) else {}
    quantification = upstream.get("quantification") if isinstance(upstream.get("quantification"), dict) else {}
    method = quantification.get("method") or upstream.get("quantification_method")
    if method == "hisat2_featurecounts":
        return "hisat2", "HISAT2 → featureCounts"
    if method == "salmon" or upstream.get("engine") == "nfcore_rnaseq":
        return "salmon", "nf-core/rnaseq Salmon → tximport"
    return "unknown", "FASTQ upstream"


def _enrichment(manifest: dict[str, Any]) -> tuple[tuple[str, ...], int]:
    """(selected public methods, number of backend tasks) from the frozen plan."""

    from rnaseq.models import normalize_enrichment_selection, production_enrichment_backends

    planned = manifest.get("planned_downstream") if isinstance(manifest.get("planned_downstream"), dict) else {}
    try:
        selected = normalize_enrichment_selection(planned.get("enrichment") or ())
        return tuple(selected), len(production_enrichment_backends(selected))
    except (ValueError, TypeError):
        return (), 0


# --------------------------------------------------------------------- stage evaluation

def _samples_done(spec: StageSpec, tasks: list[Task], lanes: dict[str, int]) -> tuple[int, str | None]:
    """Samples whose counting process finished for every lane (per-lane processes) or once."""

    for process in spec.per_sample:
        mine = [task for task in tasks if task.process == process]
        if mine:
            finished: dict[str, int] = {}
            for task in mine:
                if task.done and task.sample:
                    finished[task.sample] = finished.get(task.sample, 0) + 1
            per_lane = process in PER_LANE_PROCESSES
            return sum(1 for sample, count in finished.items() if count >= (lanes.get(sample, 1) if per_lane else 1)), process
    return 0, None


def _stage_complete(index: int, specs: tuple[StageSpec, ...], grouped: list[list[Task]], lanes: dict[str, int]) -> bool:
    spec, tasks = specs[index], grouped[index]
    if not tasks or any(task.active or task.status == "STOPPED" for task in tasks):
        return False
    if spec.final is not None:
        return any(task.process == spec.final and task.done for task in tasks)
    if spec.per_sample and lanes:
        done, process = _samples_done(spec, tasks, lanes)
        if process is not None:
            return done >= len(lanes)
    if spec.expected is not None:
        return len({task.name for task in tasks if task.done}) >= spec.expected
    # No completion criterion: done once a later stage of the same phase has started.
    return any(task.done for task in tasks) and any(grouped[later] for later in range(index + 1, len(specs)))


def _progress(spec: StageSpec, tasks: list[Task], lanes: dict[str, int]) -> str:
    if spec.per_sample and lanes:
        done, process = _samples_done(spec, tasks, lanes)
        if process is not None or not any(task.done for task in tasks):
            return f"{done} / {len(lanes)} samples"
    running = sum(1 for task in tasks if task.status == "RUNNING")
    queued = sum(1 for task in tasks if task.status in {"QUEUED", "SUBMITTED"})
    parts = [f"{sum(1 for task in tasks if task.done)} done"] if any(task.done for task in tasks) else []
    if running:
        parts.append(f"{running} running")
    if queued:
        parts.append(f"{queued} queued")
    return ", ".join(parts)


def _complete_samples(spec: StageSpec, tasks: list[Task], lanes: dict[str, int]) -> str:
    """'n / n samples' for a finished per-sample stage; blank when the trace cannot confirm every sample."""

    if not (spec.per_sample and lanes and tasks):
        return ""
    done, process = _samples_done(spec, tasks, lanes)
    return f"{done} / {len(lanes)} samples" if process is not None and done >= len(lanes) else ""


def _evaluate_group(specs: tuple[StageSpec, ...], tasks: list[Task], group: str, lanes: dict[str, int]) -> list[Stage]:
    grouped: list[list[Task]] = [[] for _ in specs]
    for task in tasks:
        for index, spec in enumerate(specs):
            if spec.match(task):
                grouped[index].append(task)
                break
    stages: list[Stage] = []
    for index, spec in enumerate(specs):
        mine = grouped[index]
        failed, _retried = _unresolved_failures(mine, lanes)
        complete = _stage_complete(index, specs, grouped, lanes)
        progress = _progress(spec, mine, lanes)
        if group == DONE:
            stages.append(Stage(spec.label, DONE, _complete_samples(spec, mine, lanes)))
        elif group == WAITING:
            stages.append(Stage(spec.label, WAITING, "waiting"))
        elif group == NOT_STARTED:
            stages.append(Stage(spec.label, NOT_STARTED, "not started"))
        elif failed:
            stages.append(Stage(spec.label, FAILED, f"{len(failed)} failed task{'s' if len(failed) != 1 else ''}"))
        elif complete:
            stages.append(Stage(spec.label, DONE, _complete_samples(spec, mine, lanes)))
        elif group == RUNNING:
            # "Running" needs a started or finished task; submitted-but-not-started tasks only queue the stage.
            started = any(task.done or task.status == "RUNNING" for task in mine)
            queued = sum(1 for task in mine if task.status in {"QUEUED", "SUBMITTED"})
            stages.append(Stage(spec.label, RUNNING, progress) if started else Stage(spec.label, WAITING, f"{queued} queued" if queued else "waiting"))
        elif any(task.status == "STOPPED" for task in mine):
            # Tasks still in flight when the run ended: interrupted with it, or stopped because another task failed.
            if group == INTERRUPTED:
                stages.append(Stage(spec.label, INTERRUPTED, progress))
            else:
                stages.append(Stage(spec.label, NOT_STARTED, f"stopped at {progress}" if progress else "stopped"))
        elif mine:
            stages.append(Stage(spec.label, NOT_STARTED, f"stopped at {progress}" if progress else "stopped"))
        else:
            stages.append(Stage(spec.label, NOT_STARTED, "not started"))
    if group in {FAILED, INTERRUPTED} and not any(stage.state in {FAILED, INTERRUPTED} for stage in stages):
        # The phase ended badly but no task carries the outcome (e.g. Nextflow
        # failed before or between tasks): mark the first unfinished stage.
        for index, stage in enumerate(stages):
            if stage.state != DONE:
                stages[index] = Stage(stage.label, group, stage.detail if stage.detail not in {"not started", "stopped"} else "")
                break
    return stages


def build_status(run_dir: Path, state: dict[str, Any], cache: _Cache | None = None, *, project: str | None = None) -> RunStatus:
    cache = cache or _Cache()
    recorded = str(state.get("status") or "UNKNOWN")
    phase = _recorded_phase(state)
    overall, liveness_note = effective_state(run_dir, state, cache)
    live = overall == "RUNNING"

    manifest = cache.load(run_dir / "frozen" / "input_manifest.yaml", _yaml)
    manifest = manifest if isinstance(manifest, dict) else {}
    contract = cache.load(run_dir / "frozen" / "downstream_contract.json", _json)
    source = contract.get("source") if isinstance(contract, dict) and isinstance(contract.get("source"), dict) else {}
    route, route_label = _route(manifest)
    if source.get("type") == "raw_counts":
        route, route_label = "raw_counts", "raw counts → DESeq2"
    lanes = cache.load(run_dir / "frozen" / "samplesheet.csv", _samples) or cache.load(run_dir / "frozen" / "metadata.csv", _samples) or {}
    planned = manifest.get("planned_downstream") if isinstance(manifest.get("planned_downstream"), dict) else {}
    preset = planned.get("preset") if isinstance(planned.get("preset"), str) else None
    enrichment, backends = _enrichment(manifest)

    stages: list[Stage] = [Stage("Input / preflight", _group_state("validation", phase, overall))]
    all_tasks: list[Task] = []

    # Upstream
    upstream_trace = _upstream_trace(run_dir)
    upstream_group = _group_state("upstream", phase, overall)
    if route == "raw_counts":
        stages.append(Stage("Upstream", SKIPPED, "not applicable (raw-count input)"))
    elif source.get("reused_from") and upstream_group == DONE:
        stages.append(Stage("Upstream", DONE, f"reused from {source['reused_from']}"))
    else:
        trace = cache.load(upstream_trace, _trace_rows) if upstream_trace else None
        submitted = cache.load(run_dir / "logs" / "upstream.stdout.log", _submissions)
        tasks = _collect_tasks(trace, submitted, live=live and upstream_group == RUNNING, work_dir=_work_dir(state, "upstream_command"))
        all_tasks += tasks
        specs = HISAT2_STAGES if route == "hisat2" else SALMON_STAGES if route == "salmon" else GENERIC_UPSTREAM
        if not tasks and upstream_group not in {DONE, WAITING, NOT_STARTED}:
            stages.append(Stage("Upstream", upstream_group, "starting; no task submitted yet" if upstream_group == RUNNING else "no task trace recorded"))
        else:
            stages += _evaluate_group(specs, tasks, upstream_group, lanes)

    # Downstream
    downstream_group = _group_state("downstream", phase, overall)
    if state.get("downstream_skipped") or preset == "QC":
        stages.append(Stage("L1 / L2", SKIPPED, "skipped (technical QC only)"))
    else:
        trace_path = run_dir / "provenance" / "downstream.trace.txt"
        trace = cache.load(trace_path, _trace_rows) if trace_path.is_file() else None
        submitted = cache.load(run_dir / "logs" / "downstream.stdout.log", _submissions)
        tasks = _collect_tasks(trace, submitted, live=live and downstream_group == RUNNING, work_dir=_work_dir(state, "downstream_command"))
        all_tasks += tasks
        # Without a readable frozen plan the downstream modules are unknown; do not guess them.
        specs = _downstream_stages(preset, enrichment, backends) if manifest else (StageSpec("Downstream", lambda task: True),)
        if not tasks and downstream_group not in {DONE, WAITING, NOT_STARTED}:
            stages.append(Stage("Downstream", downstream_group, "starting; no task submitted yet" if downstream_group == RUNNING else "no task trace recorded"))
        else:
            stages += _evaluate_group(specs, tasks, downstream_group, lanes)

    stages.append(Stage("Delivery", _group_state("delivery", phase, overall)))

    failed, retried = _unresolved_failures(all_tasks, lanes)
    counts = TaskCounts(
        completed=sum(1 for task in all_tasks if task.done),
        cached=sum(1 for task in all_tasks if task.status == "CACHED"),
        running=sum(1 for task in all_tasks if task.status == "RUNNING"),
        queued=sum(1 for task in all_tasks if task.status in {"QUEUED", "SUBMITTED"}),
        failed=len(failed),
        retried=retried,
        evidence=bool(all_tasks),
    )
    running = tuple(sorted((task for task in all_tasks if task.status == "RUNNING"), key=lambda task: task.started or 0.0))

    provenance = cache.load(run_dir / "provenance" / "run_provenance.yaml", _yaml)
    resources = provenance.get("runtime_resources") if isinstance(provenance, dict) else None
    launch_dir = provenance.get("execution_launch_dir") if isinstance(provenance, dict) else None
    logs = tuple(path for path in (
        run_dir / EXECUTION_LOG, run_dir / "logs" / "upstream.stdout.log", run_dir / "logs" / "downstream.stdout.log",
    ) if path.is_file())
    activity = [run_dir / "run_state.json", *logs, run_dir / "provenance" / "downstream.trace.txt"]
    if isinstance(launch_dir, str) and launch_dir:
        activity.append(Path(launch_dir) / ".nextflow.log")  # Nextflow's own heartbeat; stat only
    if upstream_trace:
        activity.append(upstream_trace)
    mtimes = []
    for path in activity:
        try:
            mtimes.append(path.stat().st_mtime)
        except OSError:
            pass
    attempt = _attempt(state)
    delivery = state.get("delivery")
    delivery_path = Path(delivery) if overall == "SUCCESS" and isinstance(delivery, str) and delivery else None
    project_block = manifest.get("project") if isinstance(manifest.get("project"), dict) else {}
    project_name = project or project_block.get("id") or _project_of(run_dir)
    return RunStatus(
        project=str(project_name),
        case_id=_case_of(run_dir, state),
        run_id=_run_of(run_dir, state),
        run_dir=run_dir,
        recorded_state=recorded,
        state=overall,
        phase="completed" if overall == "SUCCESS" else phase,
        route=route_label,
        started_at=state.get("started_at"),
        completed_at=state.get("completed_at"),
        stages=tuple(stages),
        running=running,
        failed=tuple(failed),
        tasks=counts,
        liveness=liveness_note,
        error=state.get("error") if overall in {"FAILED", "INTERRUPTED"} else None,
        resources=resources if isinstance(resources, dict) else None,
        last_activity=max(mtimes) if mtimes else None,
        logs=logs,
        attempt=attempt,
        delivery=delivery_path,
    )


def _latest(paths: list[Path]) -> Path | None:
    candidates = [path for path in paths if path.is_file()]
    return max(candidates, key=lambda path: (path.name, path.stat().st_mtime)) if candidates else None


def _upstream_trace(run_dir: Path) -> Path | None:
    nfcore = sorted((run_dir / "upstream" / "nfcore_rnaseq" / "pipeline_info").glob("execution_trace_*.txt"))
    return _latest([*nfcore, run_dir / "provenance" / "upstream.trace.txt"])


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
        return f"{hours}h {minutes:02d}m {secs:02d}s"
    return f"{minutes}m {secs:02d}s"


def _clock(seconds: float) -> str:
    seconds = max(0, int(seconds))
    hours, remainder = divmod(seconds, 3600)
    minutes, secs = divmod(remainder, 60)
    return f"{hours}:{minutes:02d}:{secs:02d}" if hours else f"{minutes:02d}:{secs:02d}"


def _elapsed_seconds(started_at: str | None, completed_at: str | None, state: str, now: datetime) -> float | None:
    started = _parse_time(started_at)
    if started is None:
        return None
    finished = _parse_time(completed_at)
    if finished is None:
        if state != "RUNNING":
            return None
        finished = now.astimezone(started.tzinfo) if started.tzinfo else now.replace(tzinfo=None)
    return (finished - started).total_seconds()


def _task_summary(status: RunStatus) -> str | None:
    counts = status.tasks
    if not counts.evidence:
        return None
    completed = f"{counts.completed} completed" + (f" ({counts.cached} cached)" if counts.cached else "")
    if status.state == "RUNNING":
        text = f"{completed} / {counts.running} running / {counts.queued} queued / {counts.failed} failed"
        text += "\n(further tasks are created as Nextflow progresses; total not yet known)"
    else:
        text = f"{completed} / {counts.failed} failed"
    if counts.retried:
        text += f"\n{counts.retried} failed attempt{'s' if counts.retried != 1 else ''} later succeeded on retry"
    return text


def _resource_lines(resources: dict[str, Any]) -> list[str]:
    effective = resources.get("effective") if isinstance(resources.get("effective"), dict) else None
    if not effective:
        return []
    lines = [
        "Resources (configured Nextflow ceiling, not live usage)",
        f"{'CPU':<12}{effective.get('cpus')} max",
        f"{'Memory':<12}{effective.get('memory_gib')} GiB max",
    ]
    policy = resources.get("policy") if isinstance(resources.get("policy"), dict) else None
    if policy:
        modes = {policy.get("cpu_mode"), policy.get("memory_mode")}
        mode = "auto" if modes == {"auto"} else "explicit" if modes == {"explicit"} else "mixed auto/explicit"
        detected = policy.get("detected") if isinstance(policy.get("detected"), dict) else {}
        reserve = policy.get("os_reserve") if isinstance(policy.get("os_reserve"), dict) else {}
        detail = f"{'Policy':<12}{mode}; usable {detected.get('usable_cpus', '?')} CPUs / {detected.get('usable_memory_gib', '?')} GiB"
        if reserve.get("cpus") is not None or reserve.get("memory_gib") is not None:
            detail += f"; OS reserve {reserve.get('cpus') or 0} CPUs / {reserve.get('memory_gib') or 0} GiB"
        lines.append(detail)
        for item in policy.get("process_tuning") or []:
            if isinstance(item, dict):
                lines.append(f"{'Tuning':<12}{item.get('selector')} memory {item.get('memory_gib')} GiB/task; CPUs unchanged")
    else:
        requested = resources.get("requested") if isinstance(resources.get("requested"), dict) else {}
        lines.append(f"{'Requested':<12}{requested.get('cpus', '?')} CPUs / {requested.get('memory_gib', '?')} GiB (pre-policy run)")
    if resources.get("clamped"):
        lines.append("Note        project limits were clamped to this machine's capacity")
    return lines


def render_status(status: RunStatus, *, now: datetime | None = None, ascii_only: bool = False) -> str:
    now = now or datetime.now().astimezone()
    symbols = ASCII_SYMBOLS if ascii_only else SYMBOLS
    elapsed = _elapsed_seconds(status.started_at, status.completed_at, status.state, now)
    headline = f"{status.state:<9} " + (format_duration(elapsed) if elapsed is not None else "")
    if status.state == "RUNNING" and elapsed is not None:
        headline += " elapsed"
    lines = [f"{status.project} / {status.case_id}", "-" * 44 if ascii_only else _RULE, headline.rstrip()]
    if status.liveness and status.liveness.startswith(("unverified", "stale")):
        lines.append(f"Note: {status.liveness}")

    lines += ["", f"Pipeline ({status.route})"]
    width = max(len(stage.label) for stage in status.stages) + 2
    for stage in status.stages:
        detail = stage.detail
        if stage.state == INTERRUPTED and not detail:
            detail = "interrupted"
        elif stage.state == FAILED and not detail:
            detail = "failed"
        lines.append(f"{symbols[stage.state]} {stage.label:<{width}}{detail}".rstrip())

    if status.state == "RUNNING" and status.running:
        lines += ["", "Running now"]
        shown = status.running[:10]
        sample_width = max(len(task.tag or "-") for task in shown) + 2
        process_width = max(len(task.process) for task in shown) + 2
        for task in shown:
            since = _clock(now.timestamp() - task.started) if task.started else "-"
            lines.append(f"{task.tag or '-':<{sample_width}}{task.process:<{process_width}}{since}")
        if len(status.running) > len(shown):
            lines.append(f"... {len(status.running) - len(shown)} more")

    if status.failed:
        lines += ["", "Failed tasks"]
        lines += [f"{task.name}" for task in status.failed[:8]]
        if len(status.failed) > 8:
            lines.append(f"... {len(status.failed) - 8} more")

    summary = _task_summary(status)
    if summary:
        lines += ["", "Tasks", summary]

    if status.error:
        first = str(status.error).strip().splitlines()[0] if str(status.error).strip() else ""
        lines += ["", f"Error: {first[:300]}"]

    if status.resources:
        resource_lines = _resource_lines(status.resources)
        if resource_lines:
            lines += ["", *resource_lines]

    lines += ["", f"{'Started':<12}{_format_time(status.started_at)}"]
    if status.completed_at:
        lines.append(f"{'Completed':<12}{_format_time(status.completed_at)}")
    lines.append(f"{'Run':<12}{status.run_id}" + (f" ({status.attempt})" if status.attempt else ""))
    lines.append(f"{'Run dir':<12}{status.run_dir}")
    if status.delivery is not None:
        lines.append(f"{'Delivery':<12}{status.delivery}")
    execution_log = status.run_dir / EXECUTION_LOG
    if execution_log in status.logs:
        lines.append(f"{'Log':<12}{execution_log}")
    elif status.logs:
        lines.append(f"{'Logs':<12}{execution_log.parent} (no per-run rnaseq.log; run predates it)")
    if status.state == "RUNNING":
        if status.last_activity:
            lines.append(f"{'Activity':<12}last write {format_duration(now.timestamp() - status.last_activity)} ago")
        lines.append(f"{'Last update':<12}{now.strftime('%H:%M:%S')}")
    text = "\n".join(lines)
    if ascii_only:
        text = text.replace("→", "->").encode("ascii", "replace").decode("ascii")
    return text


def render_history(project_dir: Path, runs: list[tuple[Path, dict[str, Any]]], cache: _Cache, *, latest: Path | None = None, now: datetime | None = None) -> str:
    """Compact one-line-per-run history; ``latest`` (the run the default view shows) is starred."""

    now = now or datetime.now().astimezone()
    rows = [("", "CASE", "RUN", "STATUS", "PHASE", "STARTED", "ELAPSED", "DELIVERY", "ATTEMPT")]
    for run_dir, state in runs:
        shown, note = effective_state(run_dir, state, cache)
        label = shown + (" (stale)" if note and note.startswith("stale") else " (unverified)" if note and note.startswith("unverified") else "")
        elapsed = _elapsed_seconds(state.get("started_at"), state.get("completed_at"), shown, now)
        delivery = state.get("delivery")
        available = shown == "SUCCESS" and isinstance(delivery, str) and Path(delivery).is_dir()
        rows.append((
            "*" if latest == run_dir else "",
            _case_of(run_dir, state), _run_of(run_dir, state), label,
            "completed" if shown == "SUCCESS" else _recorded_phase(state),
            _format_time(state.get("started_at")), format_duration(elapsed) if elapsed is not None else "-",
            "available" if available else "-", _attempt(state) or "original",
        ))
    widths = [max(len(row[index]) for row in rows) for index in range(len(rows[0]))]
    lines = ["  ".join(value.ljust(width) for value, width in zip(row, widths)).rstrip() for row in rows]
    lines += ["", f"Run directories: {project_dir / 'runs'}/CASE/RUN",
              "* shown by default ('rnaseq status PROJECT'); select another with --case/--run."]
    return "\n".join(lines)

