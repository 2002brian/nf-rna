"""Background ``rnaseq run`` submissions from the GUI.

Analysis state is never kept here: a run's state is its own ``run_state.json``
and ``logs/rnaseq.process.json``, read through :mod:`rnaseq.run_status`.  The
GUI keeps only what it needs to find the process it started:

* ``active_run.lock`` -- the single active submission, created atomically
  (``O_CREAT | O_EXCL``) so two clicks or two browser tabs cannot both submit;
* ``launches/<id>.json`` -- what was launched, how, and which process owns it;
* ``logs/<id>.log`` -- the launcher's stdout/stderr, which exists before the
  run's own log directory does (and is the only log of a run refused by
  ``rnaseq run``'s preflight).

A run is started either as a transient ``systemd --user`` unit, which survives
a stop or restart of the GUI service, or as a detached process in a new
session, which survives the GUI process when no user service manager is used.
"""

from __future__ import annotations

import getpass
import json
import os
import shutil
import socket
import subprocess
import sys
import threading
import uuid
from dataclasses import asdict, dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from rnaseq.gui.config import GuiConfig, gui_state_dir
from rnaseq.run_status import _Cache, _boot_id, _process_start_ticks, effective_state, process_liveness

LOCK_NAME = "active_run.lock"
PROCESS_RECORD = Path("logs") / "rnaseq.process.json"
# A lock without a process identity belongs to a submission that is still
# starting; after this many seconds it is treated as abandoned.
STARTING_GRACE_SECONDS = 120
# Environment the analysis needs; a transient unit starts from the user
# manager's environment, so these are passed on explicitly.
_PASSED_ENV = ("PATH", "HOME", "LANG", "LC_ALL", "TZ", "TMPDIR", "CONDARC", "JAVA_HOME", "JAVA_CMD")
_PASSED_ENV_PREFIXES = ("RNASEQ_", "NXF_", "NF_RNA_", "CONDA_", "XDG_")


class ActiveRunError(RuntimeError):
    """Another GUI submission is still active."""


class LaunchError(RuntimeError):
    """The background process could not be started."""


@dataclass
class Launch:
    job_id: str
    project: str
    case_id: str
    command: list[str]
    launched_at: str
    launcher: str
    log: str
    unit: str | None = None
    process: dict[str, Any] = field(default_factory=dict)
    run_dir: str | None = None

    @property
    def log_path(self) -> Path:
        return Path(self.log)


def _now() -> str:
    return datetime.now(timezone.utc).astimezone().replace(microsecond=0).isoformat()


def _write_json(path: Path, data: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    temporary.write_text(json.dumps(data, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    os.replace(temporary, path)


def _read_json(path: Path) -> dict[str, Any] | None:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return None
    return value if isinstance(value, dict) else None


def lock_path() -> Path:
    return gui_state_dir() / LOCK_NAME


def _launch_path(job_id: str) -> Path:
    return gui_state_dir() / "launches" / f"{job_id}.json"


def save_launch(launch: Launch) -> None:
    _write_json(_launch_path(launch.job_id), asdict(launch))


def load_launch(job_id: str) -> Launch | None:
    data = _read_json(_launch_path(job_id))
    if data is None:
        return None
    try:
        return Launch(**data)
    except TypeError:
        return None


def list_launches(project: Path | None = None) -> list[Launch]:
    """Launch records, newest first (optionally only those of one project)."""

    directory = gui_state_dir() / "launches"
    launches = [load_launch(path.stem) for path in directory.glob("*.json")] if directory.is_dir() else []
    wanted = str(project.resolve()) if project is not None else None
    found = [item for item in launches if item is not None and (wanted is None or item.project == wanted)]
    return sorted(found, key=lambda item: item.launched_at, reverse=True)


# --------------------------------------------------------------------- the single active run

def read_lock() -> dict[str, Any] | None:
    return _read_json(lock_path())


def _lock_is_live(lock: dict[str, Any]) -> bool:
    process = lock.get("process")
    if isinstance(process, dict) and process.get("pid"):
        alive, _note = process_liveness(process)
        return alive is not False
    try:
        started = datetime.fromisoformat(str(lock.get("launched_at")))
    except ValueError:
        return False
    return (datetime.now(timezone.utc) - started).total_seconds() < STARTING_GRACE_SECONDS


def clear_stale_lock() -> bool:
    """Remove the lock when the process it names no longer exists; the run's own files are untouched."""

    lock = read_lock()
    path = lock_path()
    if lock is None:
        if path.exists():  # unreadable or half-written lock of a crashed GUI
            try:
                age = datetime.now().timestamp() - path.stat().st_mtime
            except OSError:
                return False
            if age > STARTING_GRACE_SECONDS:
                path.unlink(missing_ok=True)
                return True
        return False
    if _lock_is_live(lock):
        return False
    path.unlink(missing_ok=True)
    return True


def _acquire_lock(record: dict[str, Any]) -> None:
    path = lock_path()
    path.parent.mkdir(parents=True, exist_ok=True)
    for attempt in range(2):
        try:
            descriptor = os.open(path, os.O_CREAT | os.O_EXCL | os.O_WRONLY, 0o600)
        except FileExistsError:
            if attempt == 0 and clear_stale_lock():
                continue
            lock = read_lock() or {}
            raise ActiveRunError(
                f"A run is already active: case {lock.get('case_id', '?')} of {lock.get('project', '?')} "
                f"(submitted {lock.get('launched_at', '?')}). Only one GUI run can be active at a time."
            ) from None
        with os.fdopen(descriptor, "w", encoding="utf-8") as handle:
            json.dump(record, handle, indent=2, sort_keys=True)
        return


def _release_lock(job_id: str) -> None:
    lock = read_lock()
    if lock is not None and lock.get("job_id") == job_id:
        lock_path().unlink(missing_ok=True)


def active_launch() -> Launch | None:
    """The launch holding the lock, after clearing a lock whose process is gone."""

    clear_stale_lock()
    lock = read_lock()
    if lock is None:
        return None
    return load_launch(str(lock.get("job_id"))) or Launch(
        job_id=str(lock.get("job_id")), project=str(lock.get("project")), case_id=str(lock.get("case_id")),
        command=[], launched_at=str(lock.get("launched_at")), launcher=str(lock.get("launcher", "?")), log="",
    )


# --------------------------------------------------------------------- launching

def rnaseq_command(project: Path, case_id: str, reuse_upstream: str | None = None) -> list[str]:
    """The ``rnaseq run`` argument list (never a shell string) for this installation."""

    executable = Path(sys.executable).with_name("rnaseq")
    base = [str(executable)] if executable.is_file() else [sys.executable, "-m", "rnaseq.cli"]
    command = [*base, "run", str(project), "--case-id", case_id, "--yes"]
    if reuse_upstream:
        command.extend(["--reuse-upstream", reuse_upstream])
    return command


def child_environment(config: GuiConfig) -> dict[str, str]:
    """This process's environment, plus the configured CONDARC unless one is already set."""

    env = dict(os.environ)
    if config.condarc is not None and "CONDARC" not in env and config.condarc.is_file():
        env["CONDARC"] = str(config.condarc)
    return env


def _systemd_user_available() -> bool:
    if shutil.which("systemd-run") is None or shutil.which("systemctl") is None:
        return False
    try:
        result = subprocess.run(["systemctl", "--user", "is-system-running"], capture_output=True, text=True, timeout=10, check=False)
    except (OSError, subprocess.TimeoutExpired):
        return False
    return result.stdout.strip() in {"running", "degraded"}


def linger_enabled() -> bool:
    return Path("/var/lib/systemd/linger", getpass.getuser()).exists()


def launcher_mode(config: GuiConfig) -> str:
    """``systemd`` (transient user unit) or ``process`` (detached session).

    Under a systemd service, children share the service's cgroup and would be
    stopped with it, so a transient unit is required.  Outside a service a
    transient unit is used only with lingering enabled: without it the user
    manager, and every unit it runs, stops at the last logout.
    """

    if config.launcher in {"systemd", "process"}:
        return config.launcher
    if _systemd_user_available() and (os.environ.get("INVOCATION_ID") or linger_enabled()):
        return "systemd"
    return "process"


def _identity(pid: int) -> dict[str, Any]:
    return {"pid": pid, "hostname": socket.gethostname(), "boot_id": _boot_id(), "start_ticks": _process_start_ticks(pid)}


def _start_systemd(command: list[str], env: dict[str, str], log: Path, cwd: Path, unit: str) -> dict[str, Any]:
    passed = {
        key: value for key, value in env.items()
        if key in _PASSED_ENV or key.startswith(_PASSED_ENV_PREFIXES)
    }
    arguments = [
        "systemd-run", "--user", "--unit", unit, "--collect", "--quiet",
        f"--working-directory={cwd}",
        f"--property=StandardOutput=append:{log}", f"--property=StandardError=append:{log}",
        *(f"--setenv={key}={value}" for key, value in sorted(passed.items())),
        "--", *command,
    ]
    result = subprocess.run(arguments, capture_output=True, text=True, timeout=60, check=False)
    if result.returncode != 0:
        raise LaunchError(f"systemd-run failed: {(result.stderr or result.stdout).strip()}")
    shown = subprocess.run(
        ["systemctl", "--user", "show", "--property=MainPID", "--value", unit],
        capture_output=True, text=True, timeout=30, check=False,
    )
    try:
        pid = int(shown.stdout.strip() or 0)
    except ValueError:
        pid = 0
    # MainPID 0: the command already exited (for example a refused preflight).
    return {**_identity(pid), "unit": unit} if pid else {"pid": None, "unit": unit, "exited_immediately": True}


def _start_process(command: list[str], env: dict[str, str], log: Path, cwd: Path) -> dict[str, Any]:
    with log.open("ab") as handle:
        process = subprocess.Popen(
            command, stdin=subprocess.DEVNULL, stdout=handle, stderr=subprocess.STDOUT,
            cwd=cwd, env=env, start_new_session=True, close_fds=True,
        )
    # Reap the child when it exits so its PID is not kept alive as a zombie.
    threading.Thread(target=process.wait, daemon=True).start()
    return _identity(process.pid)


def launch_run(project: Path, case_id: str, config: GuiConfig, *, reuse_upstream: str | None = None) -> Launch:
    """Start ``rnaseq run --yes`` in the background; raises :class:`ActiveRunError` if one is active."""

    project = project.resolve()
    job_id = datetime.now().strftime("%Y%m%d-%H%M%S") + "-" + uuid.uuid4().hex[:6]
    log = gui_state_dir() / "logs" / f"{job_id}.log"
    mode = launcher_mode(config)
    command = rnaseq_command(project, case_id, reuse_upstream)
    launch = Launch(job_id, str(project), case_id, command, _now(), mode, str(log), unit=f"nf-rna-run-{job_id}" if mode == "systemd" else None)
    _acquire_lock({"job_id": job_id, "project": launch.project, "case_id": case_id, "launched_at": launch.launched_at, "launcher": mode, "gui_pid": os.getpid()})
    try:
        log.parent.mkdir(parents=True, exist_ok=True)
        with log.open("a", encoding="utf-8") as handle:
            handle.write(f"{launch.launched_at} nf-rna GUI launch {job_id} ({mode}): {' '.join(command)}\n")
        env = child_environment(config)
        if mode == "systemd":
            assert launch.unit is not None
            launch.process = _start_systemd(command, env, log, project, launch.unit)
        else:
            launch.process = _start_process(command, env, log, project)
        save_launch(launch)
        lock = read_lock() or {}
        _write_json(lock_path(), {**lock, "process": launch.process, "log": launch.log, "unit": launch.unit})
    except BaseException:
        _release_lock(job_id)
        raise
    return launch


# --------------------------------------------------------------------- state of a launch

def launch_alive(launch: Launch) -> tuple[bool | None, str]:
    if not launch.process.get("pid"):
        return False, "the launched rnaseq process exited immediately"
    return process_liveness(launch.process)


def launch_run_dir(launch: Launch) -> Path | None:
    """The run directory created by this launch, matched by the process that owns it."""

    if launch.run_dir:
        return Path(launch.run_dir)
    pid = launch.process.get("pid")
    case_dir = Path(launch.project) / "runs" / launch.case_id
    if not pid or not case_dir.is_dir():
        return None
    for run_dir in sorted(case_dir.iterdir(), reverse=True):
        record = _read_json(run_dir / PROCESS_RECORD)
        if record and record.get("pid") == pid and record.get("start_ticks") == launch.process.get("start_ticks"):
            launch.run_dir = str(run_dir)
            save_launch(launch)
            return run_dir
    return None


def run_display_state(run_dir: Path, state: dict[str, Any], cache: _Cache | None = None) -> tuple[str, str | None]:
    """The state shown in the GUI: the existing effective state, with a dead RUNNING shown as INCOMPLETE.

    ``rnaseq status`` reports a recorded RUNNING/CREATED run whose process is
    gone as INTERRUPTED (stale).  The GUI labels it INCOMPLETE so it is not
    mistaken for a run that recorded its own interruption; nothing is written.
    """

    shown, note = effective_state(run_dir, state, cache or _Cache())
    if shown == "INTERRUPTED" and note and note.startswith("stale"):
        return "INCOMPLETE", f"process no longer exists ({note})"
    return shown, note


def tail(path: Path, lines: int) -> str:
    """The last ``lines`` lines of a text file, read from its end."""

    try:
        with path.open("rb") as handle:
            handle.seek(0, os.SEEK_END)
            size = handle.tell()
            block, data = 65536, b""
            while size > 0 and data.count(b"\n") <= lines:
                step = min(block, size)
                size -= step
                handle.seek(size)
                data = handle.read(step) + data
    except OSError as exc:
        return f"(cannot read {path}: {exc})"
    return "\n".join(data.decode("utf-8", errors="replace").splitlines()[-lines:])
