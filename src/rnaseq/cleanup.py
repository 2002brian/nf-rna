"""Explicit, fail-closed removal of one finished run's disposable Nextflow work directory.

The target is never built from user-supplied strings.  The case and run IDs only
select a durable run record under ``PROJECT/runs/CASE/RUN``; the path to delete
is the ``execution_work_dir`` that run recorded in its provenance (or, for a run
without that key, the ``-work-dir`` values of its recorded Nextflow commands).
That recorded path must be canonical, shaped ``<base>/<case>/<run>/work``, hold
only Nextflow's ``upstream``/``downstream`` work trees, lie outside the project,
and contain none of the durable or shared locations nf-rna depends on.
"""

from __future__ import annotations

import json
import os
import shutil
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path, PurePosixPath
from typing import Any

import yaml

from rnaseq.errors import ExecutionPreflightError
from rnaseq.execution import EXECUTION_ROOT_ENV, WORK_ROOT_ENV, execution_cache_dir
from rnaseq.run_status import PROCESS_RECORD, RUN_ID_PATTERN, _work_dir, process_liveness
from rnaseq.service import validate_case_id


WORK_STAGES = ("upstream", "downstream")
CLEANUP_RECORD = Path("logs") / "work_cleanup.json"
FORCE_STATES = frozenset({"FAILED", "INTERRUPTED"})


class CleanupRefused(ExecutionPreflightError):
    """The request is well formed but cleaning this target is not allowed."""


class AlreadyAbsent(Exception):
    """The work directory vanished between planning and deletion (e.g. a concurrent clean)."""


@dataclass(frozen=True)
class WorkUsage:
    bytes: int
    files: int
    unreadable: int


@dataclass(frozen=True)
class CleanPlan:
    project: Path
    run_dir: Path
    case_id: str
    run_id: str
    recorded_state: str
    effective_state: str
    state_note: str | None
    work_dir: Path
    source: str
    execution_root: Path | None

    @property
    def exists(self) -> bool:
        return os.path.lexists(self.work_dir)

    def refusal(self, *, force: bool) -> str | None:
        """Why cleaning is not allowed now, or ``None`` when it is."""

        if self.effective_state == "SUCCESS":
            return None
        if self.effective_state in FORCE_STATES:
            if force:
                return None
            return (
                f"the run is {self.effective_state}; its work directory is retained by default for "
                "diagnosis and retry. Re-run with --force to remove it anyway."
            )
        if self.effective_state == "RUNNING":
            return (f"the run is recorded as {self.recorded_state} and cannot be shown to have stopped ({self.state_note or 'recorded as active'}); "
                "a running run's work is never removed, even with --force. If it is certainly dead, check with 'rnaseq status' on the machine that launched it.")
        return f"the run state {self.recorded_state!r} is not a finished state; nothing is removed."


def _fail(message: str) -> CleanupRefused:
    return CleanupRefused(message)


def _real_directory(path: Path, label: str) -> None:
    if path.is_symlink() or not path.is_dir():
        raise _fail(f"{label} is not a real directory: {path}")


def _load_state(run_dir: Path) -> dict[str, Any]:
    path = run_dir / "run_state.json"
    try:
        state = json.loads(path.read_text(encoding="utf-8"))
    except FileNotFoundError as exc:
        raise _fail(f"run record has no run_state.json: {run_dir}") from exc
    except (OSError, UnicodeError, ValueError) as exc:
        raise _fail(f"run_state.json is unreadable: {exc}") from exc
    if not isinstance(state, dict):
        raise _fail("run_state.json is not a mapping.")
    return state


def _load_provenance(run_dir: Path) -> dict[str, Any]:
    path = run_dir / "provenance" / "run_provenance.yaml"
    if not path.exists():
        return {}
    try:
        payload = yaml.safe_load(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, yaml.YAMLError) as exc:
        raise _fail(f"provenance/run_provenance.yaml is unreadable: {exc}") from exc
    if not isinstance(payload, dict):
        raise _fail("provenance/run_provenance.yaml is not a mapping.")
    return payload


def _canonical_recorded_path(value: object, label: str) -> Path:
    """Accept only an absolute, already-normalized path string; never repair one."""

    if not isinstance(value, str) or not value or "\x00" in value:
        raise _fail(f"recorded {label} is missing or not a path string.")
    pure = PurePosixPath(value)
    if not pure.is_absolute():
        raise _fail(f"recorded {label} is not absolute: {value!r}")
    if any(part in {".", ".."} for part in pure.parts) or os.path.normpath(value) != value:
        raise _fail(f"recorded {label} is not a normalized path: {value!r}")
    return Path(value)


def _recorded_work_dir(state: dict[str, Any], provenance: dict[str, Any]) -> tuple[Path, str]:
    commands = {stage: _work_dir(state, f"{stage}_command") for stage in WORK_STAGES}
    recorded = provenance.get("execution_work_dir")
    if recorded is not None:
        work = _canonical_recorded_path(recorded, "execution_work_dir")
        source = "provenance/run_provenance.yaml (execution_work_dir)"
    else:
        parents = {stage: path for stage, path in commands.items() if path is not None}
        if not parents:
            raise _fail(
                "this run recorded no execution work directory (no execution_work_dir in provenance and no "
                "-work-dir in its Nextflow commands); nf-rna does not guess a path to delete."
            )
        candidates = {str(path.parent) for path in parents.values()}
        if len(candidates) != 1:
            raise _fail(f"recorded Nextflow commands disagree on the work directory: {sorted(candidates)}")
        work = _canonical_recorded_path(candidates.pop(), "Nextflow -work-dir parent")
        source = "run_state.json (Nextflow -work-dir)"
    for stage, path in commands.items():
        if path is not None and path != work / stage:
            raise _fail(f"recorded {stage} -work-dir {path} does not match the recorded work directory {work}.")
    return work, source


def _protected_locations(project: Path, run_dir: Path, provenance: dict[str, Any]) -> list[Path]:
    from rnaseq.downstream_runtime import _runtime_root  # local: keeps CLI start-up light
    from rnaseq.execution import _execution_base

    protected = [Path("/"), Path.home(), Path.home() / ".nextflow", project, project / "runs", run_dir]
    for key in ("execution_root", "execution_launch_dir"):
        value = provenance.get(key)
        if isinstance(value, str) and PurePosixPath(value).is_absolute():
            protected.append(Path(value))
    try:
        protected.extend([_execution_base(), execution_cache_dir("upstream-conda"), execution_cache_dir("fastq-sha256")])
    except ExecutionPreflightError:
        pass  # an invalid RNASEQ_EXECUTION_ROOT cannot also be the recorded target
    configured = os.environ.get(WORK_ROOT_ENV)
    if configured and PurePosixPath(configured).is_absolute():
        protected.append(Path(configured))
    try:
        protected.append(_runtime_root())
    except (OSError, RuntimeError):
        pass
    expanded: list[Path] = []
    for path in protected:
        expanded.append(Path(os.path.normpath(path)))
        try:
            expanded.append(path.resolve())
        except (OSError, RuntimeError):
            pass
    return expanded


def _shared_locations(protected: list[Path]) -> list[Path]:
    """Protected locations whose *contents* are shared state, never a run's work tree."""

    from rnaseq.downstream_runtime import _runtime_root

    shared = [Path.home() / ".nextflow"]
    for factory in (lambda: execution_cache_dir("upstream-conda"), lambda: execution_cache_dir("fastq-sha256"), _runtime_root):
        try:
            shared.append(factory())
        except (ExecutionPreflightError, OSError, RuntimeError):
            pass
    return [path for path in shared if path in protected or path.resolve() in protected]


def _live_processes_using(paths: list[Path]) -> list[str]:
    """Processes on this machine whose working directory or arguments reference ``paths``.

    Liveness of the recorded rnaseq process is not enough: Nextflow runs in its
    own session and its task processes can outlive a killed rnaseq.
    """

    prefixes = [str(path) for path in paths]
    found: list[str] = []
    proc = Path("/proc")
    if not proc.is_dir():
        return found
    for entry in proc.iterdir():
        if not entry.name.isdigit() or int(entry.name) == os.getpid():
            continue
        references: list[str] = []
        try:
            references.append(os.readlink(entry / "cwd"))
        except OSError:
            pass
        try:
            references += [item.decode("utf-8", "replace") for item in (entry / "cmdline").read_bytes().split(b"\0") if item]
        except OSError:
            pass
        for reference in references:
            if any(reference == prefix or reference.startswith(prefix + "/") or f"{prefix}/" in reference for prefix in prefixes):
                try:
                    name = (entry / "comm").read_text(encoding="utf-8").strip()
                except OSError:
                    name = "?"
                found.append(f"pid {entry.name} ({name})")
                break
    return found


def _open_without_symlinks(path: Path) -> int:
    """Open a directory walking from '/' one component at a time, refusing any symbolic link."""

    flags = os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW | getattr(os, "O_CLOEXEC", 0)
    fd = os.open("/", flags)
    try:
        for part in path.parts[1:]:
            child = os.open(part, flags, dir_fd=fd)
            os.close(fd)
            fd = child
    except OSError:
        os.close(fd)
        raise
    return fd


def _validate_target(work: Path, case_id: str, run_id: str, project: Path, run_dir: Path, provenance: dict[str, Any]) -> None:
    if work.name != "work" or work.parent.name != run_id or work.parent.parent.name != case_id:
        raise _fail(f"recorded work directory {work} is not shaped <base>/{case_id}/{run_id}/work for this run.")
    if len(work.parts) < 5:  # '/', at least one base component, case, run, 'work'
        raise _fail(f"recorded work directory {work} has no base directory above the case.")
    try:
        resolved = work.resolve()
    except (OSError, RuntimeError) as exc:  # e.g. a symbolic-link loop
        raise _fail(f"recorded work directory {work} cannot be resolved safely: {exc}") from exc
    if work.is_relative_to(project) or resolved.is_relative_to(project):
        raise _fail(f"recorded work directory {work} lies inside the project; durable project data is never cleaned.")
    protected_paths = _protected_locations(project, run_dir, provenance)
    for protected in protected_paths:
        if protected == work or protected.is_relative_to(work):
            raise _fail(f"recorded work directory {work} is, or contains, protected location {protected}.")
    for shared in _shared_locations(protected_paths):
        if work.is_relative_to(shared) or resolved.is_relative_to(shared):
            raise _fail(f"recorded work directory {work} lies inside shared location {shared}.")
    if not os.path.lexists(work):
        return
    if work.is_symlink() or resolved != work:
        raise _fail(f"recorded work directory {work} is, or is reached through, a symbolic link; refusing to follow it.")
    if not work.is_dir():
        raise _fail(f"recorded work directory {work} is not a directory.")
    unexpected = []
    for child in sorted(work.iterdir()):
        if child.name not in WORK_STAGES or child.is_symlink() or not child.is_dir():
            unexpected.append(child.name)
    if unexpected:
        raise _fail(
            f"{work} contains entries that are not Nextflow upstream/downstream work trees: "
            f"{', '.join(unexpected)}. Inspect it manually; nf-rna removes only a recognisable run work directory."
        )


def _durable_links_into(run_dir: Path, work: Path) -> list[Path]:
    """Symbolic links inside the durable run record that point into the work tree."""

    found: list[Path] = []
    for root, directories, files in os.walk(run_dir, followlinks=False):
        for name in [*directories, *files]:
            path = Path(root) / name
            if path.is_symlink():
                try:
                    target = path.resolve()
                except (OSError, RuntimeError):
                    continue
                if target == work or target.is_relative_to(work):
                    found.append(path)
    return found


def _effective_state(run_dir: Path, recorded: str) -> tuple[str, str | None]:
    if recorded in {"SUCCESS", *FORCE_STATES}:
        return recorded, None
    if recorded in {"CREATED", "RUNNING"}:
        try:
            process = json.loads((run_dir / PROCESS_RECORD).read_text(encoding="utf-8"))
        except (OSError, UnicodeError, ValueError):
            process = None
        alive, note = process_liveness(process if isinstance(process, dict) else None)
        if alive is False:
            # The recording process is verifiably gone: the run can never finish.
            return "INTERRUPTED", f"stale {recorded}: {note}"
        return "RUNNING", note
    return recorded or "UNKNOWN", None


def plan_work_cleanup(project_dir: Path, case_id: str, run_id: str) -> CleanPlan:
    """Resolve and validate one run's recorded work directory without changing anything."""

    validate_case_id(case_id)
    if not isinstance(run_id, str) or not RUN_ID_PATTERN.fullmatch(run_id):
        raise ExecutionPreflightError(f"run ID must look like YYYYMMDD-HHMMSS+ZZZZ: {run_id!r}")
    if not project_dir.is_dir():
        raise ExecutionPreflightError(f"project directory not found: {project_dir}")
    project = project_dir.resolve(strict=True)
    runs = project / "runs"
    run_dir = runs / case_id / run_id
    if not os.path.lexists(run_dir):
        raise ExecutionPreflightError(f"no recorded run {case_id}/{run_id} in {project}; list runs with: rnaseq status PROJECT --all")
    for path, label in ((runs, "PROJECT/runs"), (runs / case_id, "case directory"), (run_dir, "run directory")):
        _real_directory(path, label)
    state = _load_state(run_dir)
    for key, expected in (("case_id", case_id), ("run_id", run_id)):
        if state.get(key) not in (None, expected):
            raise _fail(f"run_state.json records {key}={state.get(key)!r}, not {expected!r}.")
    provenance = _load_provenance(run_dir)
    work, source = _recorded_work_dir(state, provenance)
    _validate_target(work, case_id, run_id, project, run_dir, provenance)
    recorded = str(state.get("status") or "UNKNOWN")
    effective, note = _effective_state(run_dir, recorded)
    execution_root = provenance.get("execution_root")
    return CleanPlan(
        project=project, run_dir=run_dir, case_id=case_id, run_id=run_id,
        recorded_state=recorded, effective_state=effective, state_note=note,
        work_dir=work, source=source,
        execution_root=Path(execution_root) if isinstance(execution_root, str) and execution_root else None,
    )


def check_no_live_processes(plan: CleanPlan) -> str | None:
    """A refusal message when any live process still uses this run's work or launch directory."""

    paths = [plan.work_dir]
    if plan.execution_root is not None:
        paths.append(plan.execution_root / "launch")
    busy = _live_processes_using(paths)
    if busy:
        return (f"{len(busy)} live process(es) still use this run's work or launch directory ({', '.join(busy[:5])}); "
                "nothing is removed, even with --force. Stop them first.")
    return None


def measure_work(path: Path) -> WorkUsage:
    """Allocated size of a tree without following symbolic links (Nextflow stages inputs as links)."""

    total = files = unreadable = 0

    def _error(_exc: OSError) -> None:
        nonlocal unreadable
        unreadable += 1

    for root, directories, names in os.walk(path, followlinks=False, onerror=_error):
        for name, is_file in [*((item, False) for item in directories), *((item, True) for item in names)]:
            try:
                status = os.lstat(os.path.join(root, name))
            except OSError:
                unreadable += 1
                continue
            total += getattr(status, "st_blocks", 0) * 512 or status.st_size
            files += is_file
    return WorkUsage(total, files, unreadable)


def format_bytes(value: int) -> str:
    for unit in ("B", "KiB", "MiB", "GiB"):
        if value < 1024 or unit == "GiB":
            return f"{value:.0f} {unit}" if unit == "B" else f"{value:.1f} {unit}"
        value /= 1024
    raise AssertionError("unreachable")


def _record_cleanup(plan: CleanPlan, usage: WorkUsage, *, forced: bool) -> str | None:
    from rnaseq.models import PIPELINE_VERSION

    path = plan.run_dir / CLEANUP_RECORD
    try:
        events = json.loads(path.read_text(encoding="utf-8")) if path.exists() else []
        if not isinstance(events, list):
            events = []
        events.append({
            "cleaned_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
            "work_dir": str(plan.work_dir), "source": plan.source,
            "recorded_state": plan.recorded_state, "effective_state": plan.effective_state,
            "forced": forced, "approximate_bytes": usage.bytes, "files": usage.files,
            "nf_rna_version": PIPELINE_VERSION,
        })
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(events, indent=2) + "\n", encoding="utf-8")
    except (OSError, UnicodeError, ValueError) as exc:
        return f"the cleanup was not recorded in {path}: {exc}"
    return None


def execute_work_cleanup(plan: CleanPlan, *, force: bool, usage: WorkUsage) -> str | None:
    """Delete exactly the planned work directory after re-validating it; returns a warning, if any."""

    if not shutil.rmtree.avoids_symlink_attacks:
        raise _fail("this platform's recursive delete is not symlink-attack resistant; refusing to clean.")
    current = plan_work_cleanup(plan.project, plan.case_id, plan.run_id)
    if current.work_dir != plan.work_dir or current.effective_state != plan.effective_state:
        raise _fail("the run record changed while cleanup was being confirmed; nothing was removed.")
    refusal = current.refusal(force=force)
    if refusal:
        raise _fail(refusal)
    if not current.exists:
        raise AlreadyAbsent(str(current.work_dir))
    busy = check_no_live_processes(current)
    if busy:
        raise _fail(busy)
    links = _durable_links_into(current.run_dir, current.work_dir)
    if links:
        raise _fail(f"durable run files link into the work directory (for example {links[0]}); nothing was removed.")
    # Delete relative to a parent directory opened without following any
    # symbolic link, so a component swapped after validation cannot redirect it.
    parent_fd = _open_without_symlinks(current.work_dir.parent)
    try:
        if os.stat(current.work_dir.parent, follow_symlinks=False).st_ino != os.fstat(parent_fd).st_ino:
            raise _fail("the work directory's parent changed during cleanup; nothing was removed.")
        shutil.rmtree(current.work_dir.name, dir_fd=parent_fd)
    finally:
        os.close(parent_fd)
    for directory in (current.work_dir.parent, current.work_dir.parent.parent):
        protected = _protected_locations(current.project, current.run_dir, _load_provenance(current.run_dir))
        if directory in protected or (current.execution_root is not None and directory == current.execution_root):
            break
        try:
            directory.rmdir()  # only an empty per-run / per-case shell under the work root
        except OSError:
            break
    return _record_cleanup(current, usage, forced=force and current.effective_state != "SUCCESS")
