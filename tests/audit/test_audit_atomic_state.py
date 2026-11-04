"""Control-plane files are replaced atomically and damaged state is preserved (audit R-01)."""

from __future__ import annotations

import json
import os
import stat
import threading

import pytest

from rnaseq.atomic_io import atomic_write_text
from rnaseq.service import CaseRun, _write_state


def _run(tmp_path) -> CaseRun:
    run_dir = tmp_path / "runs" / "CASE" / "20261104-101010+0800"
    (run_dir / "logs").mkdir(parents=True)
    return CaseRun(case_id="CASE", run_id="20261104-101010+0800", run_dir=run_dir, started_at="2026-11-04T10:10:10+08:00")


def test_failed_write_keeps_the_previous_complete_file(tmp_path, monkeypatch):
    target = tmp_path / "run_state.json"
    atomic_write_text(target, '{"status": "RUNNING"}\n')

    def crash(*_args, **_kwargs):
        raise OSError("simulated crash before rename")

    monkeypatch.setattr("rnaseq.atomic_io.os.replace", crash)
    with pytest.raises(OSError, match="simulated crash"):
        atomic_write_text(target, '{"status": "SUCCESS"}\n')
    assert json.loads(target.read_text(encoding="utf-8")) == {"status": "RUNNING"}
    assert sorted(path.name for path in tmp_path.iterdir()) == ["run_state.json"]  # no temporary left behind


def test_atomic_write_keeps_ordinary_permissions(tmp_path):
    target = tmp_path / "new.json"
    atomic_write_text(target, "{}\n")
    mask = os.umask(0)
    os.umask(mask)
    assert stat.S_IMODE(target.stat().st_mode) == 0o666 & ~mask
    os.chmod(target, 0o640)
    atomic_write_text(target, "{}\n")
    assert stat.S_IMODE(target.stat().st_mode) == 0o640


def test_concurrent_readers_never_see_a_partial_state(tmp_path):
    target = tmp_path / "run_state.json"
    payload = json.dumps({"status": "RUNNING", "padding": "x" * 200_000}) + "\n"
    atomic_write_text(target, payload)
    stop, problems = threading.Event(), []

    def reader():
        while not stop.is_set():
            try:
                json.loads(target.read_text(encoding="utf-8"))
            except ValueError as exc:  # pragma: no cover - the failure being guarded against
                problems.append(exc)

    thread = threading.Thread(target=reader)
    thread.start()
    try:
        for _ in range(200):
            atomic_write_text(target, payload)
    finally:
        stop.set()
        thread.join()
    assert problems == []


def test_run_state_writes_go_through_the_atomic_writer(tmp_path, monkeypatch):
    calls = []
    real = atomic_write_text
    monkeypatch.setattr("rnaseq.service.atomic_write_text", lambda path, value: (calls.append(path), real(path, value)))
    run = _run(tmp_path)
    _write_state(run, "RUNNING", phase="upstream")
    assert run.state_path in calls


def test_damaged_state_is_preserved_not_silently_discarded(tmp_path):
    run = _run(tmp_path)
    _write_state(run, "RUNNING", phase="upstream", command=["nextflow", "run"])
    run.state_path.write_text('{"status": "RUNNING", "comm', encoding="utf-8")  # external truncation
    _write_state(run, "FAILED", error="boom")
    state = json.loads(run.state_path.read_text(encoding="utf-8"))
    assert state["status"] == "FAILED" and state["started_at"] == run.started_at
    (preserved,) = run.run_dir.glob("run_state.json.corrupt-*")
    assert preserved.read_text(encoding="utf-8") == '{"status": "RUNNING", "comm'
    log = (run.run_dir / "logs").rglob("*")
    assert any("unreadable" in path.read_text(encoding="utf-8") for path in log if path.is_file())
