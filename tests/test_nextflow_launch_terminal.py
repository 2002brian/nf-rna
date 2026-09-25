"""`rnaseq run` must be safe as a background job of an interactive terminal.

Regression for 4T1 run grcm39-v130-salmon-r3: launched as
``nohup rnaseq run ... &``, Nextflow's ANSI progress log initialised JLine,
which ran ``stty -icanon min 1 -icrnl -inlcr < /dev/tty`` from a background
process group.  The kernel answered with SIGTTOU and the whole workflow stayed
stopped before any task ran.
"""

from __future__ import annotations

import fcntl
import os
import shutil
import signal
import subprocess
import termios
import time
from pathlib import Path

import pytest

from rnaseq import service


# Exactly what JLine 2.9 (bundled in Nextflow 26.04.6) runs to initialise a Unix terminal.
JLINE_STTY = "stty -icanon min 1 -icrnl -inlcr < /dev/tty"
AGENT_ENVIRONMENT = ("CLAUDECODE", "AI_AGENT", "NXF_AGENT_MODE", "NXF_ANSI_LOG")


def _in_background_job_of_a_terminal(action) -> tuple[str, int]:
    """Run ``action`` like ``nohup CMD &`` from an interactive shell.

    A session leader owns a pty as its controlling terminal; the job runs in
    another (background) process group of that session, with stdin
    /dev/null, SIGHUP ignored and default SIGTTOU handling.  Returns
    ("stopped", signal) if job control stopped it, else ("exited", status).
    """

    master, slave = os.openpty()
    read_end, write_end = os.pipe()
    leader = os.fork()
    if leader == 0:
        try:
            os.close(read_end)
            os.setsid()
            fcntl.ioctl(slave, termios.TIOCSCTTY, 0)
            job = os.fork()
            if job == 0:
                try:
                    os.setpgid(0, 0)
                    signal.signal(signal.SIGTTOU, signal.SIG_DFL)
                    signal.signal(signal.SIGTTIN, signal.SIG_DFL)
                    signal.signal(signal.SIGHUP, signal.SIG_IGN)
                    null = os.open(os.devnull, os.O_RDONLY)
                    os.dup2(null, 0)
                    os._exit(action())
                except BaseException:
                    os._exit(99)
            deadline = time.monotonic() + 120
            outcome = "timeout:0"
            while time.monotonic() < deadline:
                pid, status = os.waitpid(job, os.WUNTRACED | os.WNOHANG)
                if pid and os.WIFSTOPPED(status):
                    outcome = f"stopped:{os.WSTOPSIG(status)}"
                    os.killpg(job, signal.SIGKILL)
                    os.waitpid(job, 0)
                    break
                if pid and os.WIFEXITED(status):
                    outcome = f"exited:{os.WEXITSTATUS(status)}"
                    break
                time.sleep(0.05)
            else:
                os.killpg(job, signal.SIGKILL)
            os.write(write_end, outcome.encode())
        finally:
            os._exit(0)
    os.close(write_end)
    os.close(slave)
    os.close(master)
    os.waitpid(leader, 0)
    kind, value = os.read(read_end, 64).decode().split(":")
    os.close(read_end)
    return kind, int(value)


def _fake_nextflow(tmp_path: Path) -> Path:
    script = tmp_path / "nextflow"
    script.write_text(
        "#!/bin/sh\n"
        f"{JLINE_STTY} 2>/dev/null\n"
        'echo "stty_exit=$?"\n'
        'echo "ansi=${NXF_ANSI_LOG:-unset}"\n'
        'echo "stdin=$(readlink /proc/$$/fd/0)"\n'
        'echo "controlling_tty=$(ps -o tty= -p $$ | tr -d " ")"\n',
        encoding="utf-8",
    )
    script.chmod(0o755)
    return script


def _launch(tmp_path: Path, command: list[str]) -> int:
    return service._run_command(
        command, cwd=tmp_path, stdout_path=tmp_path / "stdout.log", stderr_path=tmp_path / "stderr.log",
    )


def test_harness_reproduces_the_r3_stop_with_an_attached_launcher(tmp_path):
    """Control: launching as before (inheriting the terminal) is stopped by SIGTTOU."""

    script = _fake_nextflow(tmp_path)

    def legacy_launch() -> int:
        with (tmp_path / "legacy.log").open("w") as out:
            return subprocess.run([str(script)], cwd=tmp_path, stdout=out, stderr=out, check=False).returncode

    assert _in_background_job_of_a_terminal(legacy_launch) == ("stopped", signal.SIGTTOU)


def test_background_launch_is_not_stopped_by_terminal_job_control(tmp_path):
    script = _fake_nextflow(tmp_path)
    assert _in_background_job_of_a_terminal(lambda: _launch(tmp_path, [str(script)])) == ("exited", 0)
    report = dict(line.split("=", 1) for line in (tmp_path / "stdout.log").read_text(encoding="utf-8").split())
    assert report["ansi"] == "false"
    assert report["stdin"] == "/dev/null"
    # No controlling terminal: JLine's stty fails harmlessly instead of stopping the job.
    assert report["controlling_tty"] == "?"
    assert report["stty_exit"] != "0"


def test_interactive_foreground_launch_still_works(tmp_path):
    script = _fake_nextflow(tmp_path)
    assert _launch(tmp_path, [str(script)]) == 0
    assert "ansi=false" in (tmp_path / "stdout.log").read_text(encoding="utf-8")


@pytest.mark.parametrize("signum", [signal.SIGINT, signal.SIGTERM])
def test_interrupting_rnaseq_stops_the_detached_nextflow_cleanly(tmp_path, signum):
    script = tmp_path / "nextflow"
    script.write_text(
        "#!/bin/sh\n"
        "trap 'echo terminated > terminated; exit 143' TERM\n"
        "echo $$ > started\n"
        "while :; do sleep 0.1; done\n",
        encoding="utf-8",
    )
    script.chmod(0o755)
    child = os.fork()
    if child == 0:
        try:
            signal.signal(signal.SIGINT, signal.default_int_handler)
            signal.signal(signal.SIGTERM, signal.SIG_DFL)
            _launch(tmp_path, [str(script)])
            os._exit(0)
        except KeyboardInterrupt:
            os._exit(130)
        except service._Termination:
            os._exit(143)
        except BaseException:
            os._exit(99)
    deadline = time.monotonic() + 30
    while not (tmp_path / "started").exists() and time.monotonic() < deadline:
        time.sleep(0.05)
    nextflow_pid = int((tmp_path / "started").read_text(encoding="utf-8"))
    os.kill(child, signum)
    _pid, status = os.waitpid(child, 0)
    assert os.WEXITSTATUS(status) == (130 if signum == signal.SIGINT else 143)
    assert (tmp_path / "terminated").read_text(encoding="utf-8").strip() == "terminated"
    with pytest.raises(ProcessLookupError):
        os.kill(nextflow_pid, 0)


@pytest.mark.skipif(shutil.which("nextflow") is None or shutil.which("java") is None, reason="real Nextflow is not installed")
def test_real_nextflow_ansi_log_does_not_stop_a_background_run(tmp_path, monkeypatch):
    """End to end with Nextflow's own JLine: the r3 failure mode, under the launcher."""

    for variable in AGENT_ENVIRONMENT:
        # Nextflow 26 switches to agent logging when it detects an AI agent;
        # the ANSI path that failed in r3 only runs without these variables.
        monkeypatch.delenv(variable, raising=False)
    workflow = tmp_path / "main.nf"
    workflow.write_text('process HELLO {\n  output: stdout\n  script: "echo hello"\n}\nworkflow { HELLO() | view }\n', encoding="utf-8")
    command = ["nextflow", "run", str(workflow), "-work-dir", str(tmp_path / "work")]
    assert _in_background_job_of_a_terminal(lambda: _launch(tmp_path, command)) == ("exited", 0)
    assert "hello" in (tmp_path / "stdout.log").read_text(encoding="utf-8")
