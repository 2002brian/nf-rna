"""``rnaseq status`` and per-run execution logging."""

from __future__ import annotations

import hashlib
import json
import os
import subprocess
import sys
from datetime import datetime
from pathlib import Path

import pytest
import yaml
from typer.testing import CliRunner

from conftest import base_config
from rnaseq.cli import app
from rnaseq.downstream_runtime import DownstreamRuntime
from rnaseq.execution import LocalResourceCapacity, RuntimeCheck
from rnaseq.planner import generate_plan
from rnaseq.run_status import _Cache, build_status, current_process_identity, render_status, select_run
from rnaseq.validators import validate_project

runner = CliRunner()
RUN_ID = "20260925-203032+0800"
TRACE_HEADER = "task_id\thash\tnative_id\tname\tstatus\texit\tsubmit\tduration\trealtime\t%cpu\tpeak_rss\tpeak_vmem\trchar\twchar\n"


# ------------------------------------------------------------------ shared helpers (also used by test_resource_policy)

def salmon_case_project(tmp_path: Path, monkeypatch, *, execution: dict | None = None):
    """A validated, planned Salmon FASTQ project whose runtime checks are mocked."""

    root = tmp_path / "project"
    (root / "input" / "fastq").mkdir(parents=True)
    for sample in ("C1", "C2", "T1", "T2"):
        for read in ("R1", "R2"):
            (root / "input" / "fastq" / f"{sample}_{read}.fastq.gz").write_bytes(b"fixture")
    (root / "metadata.csv").write_text("sample_id,condition\nC1,Control\nC2,Control\nT1,Treatment\nT2,Treatment\n", encoding="utf-8")
    (root / "contrasts.csv").write_text("contrast_id,factor,numerator,denominator\nTreatment_vs_Control,condition,Treatment,Control\n", encoding="utf-8")
    config = base_config()
    config["input"] = {"type": "fastq", "path": "input/fastq", "layout": "paired_end"}
    config["upstream"] = {"engine": "nfcore_rnaseq", "pipeline_version": "3.26.0", "aligner": None, "strandedness": "auto", "quantification": {"method": "salmon"}}
    config["reference"] = {"source": "igenomes", "genome": "test_reference"}
    if execution is not None:
        config["execution"] = execution
    (root / "project.yaml").write_text(yaml.safe_dump(config, sort_keys=False), encoding="utf-8")
    report = validate_project(root)
    assert report.is_valid
    generate_plan(report)
    monkeypatch.setenv("RNASEQ_EXECUTION_ROOT", str(tmp_path / "execution-root"))
    capacity = LocalResourceCapacity(16, 64, 60)
    monkeypatch.setattr("rnaseq.execution.detect_local_resource_capacity", lambda: capacity)
    monkeypatch.setattr("rnaseq.service.detect_local_resource_capacity", lambda: capacity)
    monkeypatch.setattr("rnaseq.service.check_nextflow", lambda: RuntimeCheck("Nextflow", "FOUND", "25.10.4"))
    monkeypatch.setattr("rnaseq.service.check_upstream_conda", lambda: RuntimeCheck("Conda", "FOUND", "conda 25.3.1"))
    runtime = DownstreamRuntime(
        prefix=tmp_path / "native-runtime", platform="linux-64",
        lock_filename="nf-rna-downstream-linux-64.lock.yml", lock_sha256="a" * 64,
        wheel_filename="nf_rna-1.3.0-py3-none-any.whl", wheel_sha256="b" * 64,
        nf_rna_version="1.3.0", source_revision="test-revision",
        r_scripts=({"name": "l2_analysis.R", "sha256": "c" * 64},), r_scripts_sha256="d" * 64,
    )
    monkeypatch.setattr("rnaseq.service.downstream_runtime_preflight", lambda: None)
    monkeypatch.setattr("rnaseq.service.ensure_downstream_runtime", lambda: runtime)
    return validate_project(root)


def fake_nextflow_factory(observed: list[list[str]], *, on_upstream=None):
    """A successful fake Nextflow that writes the documented outputs and plain-log lines."""

    def fake(command, *, cwd, stdout_path, stderr_path):
        observed.append(command)
        outdir = Path(command[command.index("--outdir") + 1])
        if "nf-core/rnaseq" in command:
            if on_upstream is not None:
                on_upstream(command)
            stdout_path.write_text("[ab/cdef01] Submitted process > NFCORE_RNASEQ:RNASEQ:QUANTIFY_PSEUDO_ALIGNMENT:SALMON_QUANT (C1)\n", encoding="utf-8")
            stderr_path.write_text("", encoding="utf-8")
            (outdir / "salmon").mkdir(parents=True)
            (outdir / "salmon" / "salmon.merged.gene_counts.tsv").write_text("gene_id\tC1\tC2\tT1\tT2\nGeneA\t1.0\t2.0\t3.0\t4.0\n", encoding="utf-8")
            (outdir / "salmon" / "salmon.merged.tx2gene_augmented.tsv").write_text("transcript_id\tgene_id\nTx1\tGeneA\n", encoding="utf-8")
            for sample in ("C1", "C2", "T1", "T2"):
                (outdir / "salmon" / sample).mkdir()
                (outdir / "salmon" / sample / "quant.sf").write_text("Name\tLength\tEffectiveLength\tTPM\tNumReads\nTx1\t100\t80\t1\t1\n", encoding="utf-8")
            (outdir / "multiqc" / "multiqc_data").mkdir(parents=True)
            (outdir / "multiqc" / "multiqc_report.html").write_text("<html></html>", encoding="utf-8")
            (outdir / "pipeline_info").mkdir()
            (outdir / "pipeline_info" / "execution_trace_2026-09-25_20-34-41.txt").write_text(
                TRACE_HEADER + "1\tab/cdef01\t1\tNFCORE_RNASEQ:RNASEQ:QUANTIFY_PSEUDO_ALIGNMENT:SALMON_QUANT (C1)\tCOMPLETED\t0\n", encoding="utf-8",
            )
        else:
            stdout_path.write_text("downstream\n", encoding="utf-8")
            stderr_path.write_text("", encoding="utf-8")
            (outdir / "report").mkdir(parents=True)
            (outdir / "report" / "report.html").write_text("<html></html>", encoding="utf-8")
        return 0

    return fake


def _run_dir(project: Path, case: str = "grcm39-v130-salmon-r6", run_id: str = RUN_ID) -> Path:
    run_dir = project / "runs" / case / run_id
    for child in ("frozen", "logs", "provenance", "upstream", "downstream", "delivery"):
        (run_dir / child).mkdir(parents=True, exist_ok=True)
    return run_dir


def _state(run_dir: Path, **values) -> None:
    state = {"case_id": run_dir.parent.name, "run_id": run_dir.name, "timezone": "Asia/Taipei", "started_at": "2026-09-25T20:30:32+08:00", "completed_at": None}
    state.update(values)
    (run_dir / "run_state.json").write_text(json.dumps(state, indent=2), encoding="utf-8")


def _trace(path: Path, rows: list[tuple[str, str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(TRACE_HEADER + "".join(f"{index}\t{task_hash}\t{index}\t{name}\t{status}\t0\n" for index, (task_hash, name, status) in enumerate(rows, 1)), encoding="utf-8")


def _upstream_trace(run_dir: Path) -> Path:
    return run_dir / "upstream" / "nfcore_rnaseq" / "pipeline_info" / "execution_trace_2026-09-25_20-34-41.txt"


def _tree_fingerprint(root: Path) -> dict[str, tuple[int, str]]:
    return {
        str(path.relative_to(root)): (path.stat().st_mtime_ns, hashlib.sha256(path.read_bytes()).hexdigest() if path.is_file() else "dir")
        for path in sorted(root.rglob("*"))
    }


def _dead_pid() -> int:
    process = subprocess.Popen([sys.executable, "-c", "pass"])
    process.wait()
    return process.pid


def _status(project: Path, *arguments: str) -> str:
    result = runner.invoke(app, ["status", str(project), *arguments])
    assert result.exit_code == 0, result.output
    return result.output


# ------------------------------------------------------------------ completed / failed / interrupted

def test_completed_success_status_matches_the_documented_summary(tmp_path):
    run_dir = _run_dir(tmp_path)
    _state(run_dir, status="SUCCESS", phase="delivery", completed_at="2026-09-26T00:39:36+08:00")
    _trace(_upstream_trace(run_dir), [(f"{i:02x}/000000", f"NFCORE_RNASEQ:TASK_{i}", "COMPLETED") for i in range(65)])
    _trace(run_dir / "provenance" / "downstream.trace.txt", [("aa/000000", "L2_ANALYSIS", "COMPLETED")])
    (run_dir / "provenance" / "run_provenance.yaml").write_text(yaml.safe_dump({"runtime_resources": {
        "requested": {"cpus": 22, "memory_gib": 20}, "effective": {"cpus": 22, "memory_gib": 20}, "clamped": False,
        "policy": {"cpu_mode": "explicit", "memory_mode": "explicit", "detected": {"usable_cpus": 28, "usable_memory_gib": 94.3},
                   "os_reserve": {"cpus": None, "memory_gib": None},
                   "process_tuning": [{"selector": ".*:SALMON_QUANT", "memory_gib": 19}]},
    }}), encoding="utf-8")

    output = _status(tmp_path)

    assert "Case:       grcm39-v130-salmon-r6\nRun:        20260925-203032+0800\nState:      SUCCESS\nPhase:      completed" in output
    assert "Started:    2026-09-25 20:30:32\nFinished:   2026-09-26 00:39:36\nElapsed:    4h 09m" in output
    assert "Upstream:   SUCCESS\nTasks:      65/65 completed, 0 cached, 0 failed" in output
    assert "Downstream: SUCCESS\nTasks:      1/1 completed, 0 cached, 0 failed" in output
    assert "Delivery:   SUCCESS" in output
    assert "Resources:  22 CPUs / 20 GiB" in output
    assert "Policy:     explicit; usable 28 CPUs / 94.3 GiB" in output
    assert "Tuning:     .*:SALMON_QUANT memory 19 GiB/task; CPUs unchanged" in output


def test_failed_run_reports_the_failed_stage_and_tasks(tmp_path):
    run_dir = _run_dir(tmp_path)
    _state(run_dir, status="FAILED", phase="downstream", completed_at="2026-09-25T23:00:00+08:00",
           error="downstream Nextflow failed with return code 1. Logs: x")
    _trace(_upstream_trace(run_dir), [("aa/000001", "SALMON_QUANT (C1)", "COMPLETED"), ("aa/000002", "MULTIQC", "CACHED")])
    _trace(run_dir / "provenance" / "downstream.trace.txt", [("bb/000001", "L1_ANALYSIS", "COMPLETED"), ("bb/000002", "L2_ANALYSIS", "FAILED")])

    output = _status(tmp_path)

    assert "State:      FAILED\nPhase:      downstream" in output
    assert "Upstream:   SUCCESS\nTasks:      2/2 completed, 1 cached, 0 failed" in output
    assert "Downstream: FAILED" in output and "Failed tasks:\n  - L2_ANALYSIS" in output
    assert "Delivery:   NOT STARTED" in output
    assert "Error:      downstream Nextflow failed with return code 1." in output
    assert "Elapsed:    2h 29m" in output


def test_recorded_interrupted_run(tmp_path):
    run_dir = _run_dir(tmp_path)
    _state(run_dir, status="INTERRUPTED", phase="upstream", completed_at="2026-09-25T20:40:00+08:00",
           interrupted_by="SIGTERM", error="Interrupted by SIGTERM; the run did not finish.")

    output = _status(tmp_path)

    assert "State:      INTERRUPTED\nPhase:      upstream" in output
    assert "Upstream:   INTERRUPTED" in output and "Downstream: NOT STARTED" in output
    assert "Error:      Interrupted by SIGTERM" in output


# ------------------------------------------------------------------ live and stale RUNNING

def test_active_running_status_lists_running_tasks(tmp_path):
    run_dir = _run_dir(tmp_path)
    _state(run_dir, status="RUNNING", phase="upstream")
    (run_dir / "logs" / "rnaseq.process.json").write_text(json.dumps(current_process_identity()), encoding="utf-8")
    (run_dir / "logs" / "rnaseq.log").write_text("2026-09-25T20:30:32+08:00 state RUNNING phase=upstream\n", encoding="utf-8")
    (run_dir / "logs" / "upstream.stdout.log").write_text(
        "N E X T F L O W\n"
        "[aa/000001] Submitted process > NFCORE_RNASEQ:RNASEQ:FASTQC (S1)\n"
        "[aa/000002] Submitted process > NFCORE_RNASEQ:RNASEQ:QUANTIFY_PSEUDO_ALIGNMENT:SALMON_QUANT (S1)\n"
        "[aa/000003] Submitted process > NFCORE_RNASEQ:RNASEQ:QUANTIFY_PSEUDO_ALIGNMENT:SALMON_QUANT (S2)\n",
        encoding="utf-8",
    )
    _trace(_upstream_trace(run_dir), [("aa/000001", "NFCORE_RNASEQ:RNASEQ:FASTQC (S1)", "COMPLETED")])

    output = _status(tmp_path)

    assert "State:      RUNNING\nPhase:      upstream" in output
    assert "Upstream:   RUNNING\nTasks:      1 completed (0 cached), 2 running, 0 failed; total not yet known" in output
    assert "Running:\n  - NFCORE_RNASEQ:RNASEQ:QUANTIFY_PSEUDO_ALIGNMENT:SALMON_QUANT (S1)\n  - NFCORE_RNASEQ:RNASEQ:QUANTIFY_PSEUDO_ALIGNMENT:SALMON_QUANT (S2)" in output
    assert "Downstream: PENDING" in output and "Delivery:   PENDING" in output
    assert "Activity:   last write" in output
    assert f"Log:        {run_dir / 'logs' / 'rnaseq.log'}" in output


def test_stale_running_state_whose_process_is_gone_is_interrupted_not_success(tmp_path):
    run_dir = _run_dir(tmp_path)
    _state(run_dir, status="RUNNING", phase="downstream")
    identity = {**current_process_identity(), "pid": _dead_pid()}
    (run_dir / "logs" / "rnaseq.process.json").write_text(json.dumps(identity), encoding="utf-8")
    # Even complete-looking outputs never make a dead RUNNING run a SUCCESS.
    (run_dir / "delivery" / "README.md").write_text("x", encoding="utf-8")

    output = _status(tmp_path)

    assert "State:      INTERRUPTED\nPhase:      downstream" in output
    assert "Note:       stale: recorded RUNNING, but rnaseq process" in output and "no longer exists" in output
    assert "SUCCESS" not in output.split("Downstream:")[1].split("\n")[0]


def test_reused_pid_is_detected_by_process_start_time(tmp_path):
    run_dir = _run_dir(tmp_path)
    _state(run_dir, status="RUNNING", phase="upstream")
    identity = current_process_identity()
    if identity["start_ticks"] is None:
        pytest.skip("procfs start times unavailable")
    identity["start_ticks"] -= 1
    (run_dir / "logs" / "rnaseq.process.json").write_text(json.dumps(identity), encoding="utf-8")
    assert "State:      INTERRUPTED" in _status(tmp_path)


def test_legacy_running_state_without_process_record_is_unverified(tmp_path):
    run_dir = _run_dir(tmp_path)
    _state(run_dir, status="RUNNING", phase="upstream")
    output = _status(tmp_path)
    assert "State:      RUNNING" in output
    assert "Note:       unverified: no process identity was recorded" in output


def test_running_state_from_another_host_is_unverified(tmp_path):
    run_dir = _run_dir(tmp_path)
    _state(run_dir, status="RUNNING", phase="upstream")
    (run_dir / "logs" / "rnaseq.process.json").write_text(json.dumps({"pid": 1, "hostname": "some-other-host"}), encoding="utf-8")
    assert "unverified: launched on host some-other-host" in _status(tmp_path)


# ------------------------------------------------------------------ missing / partial artifacts, read-only

def test_status_with_missing_or_corrupt_optional_artifacts(tmp_path):
    run_dir = _run_dir(tmp_path)
    _state(run_dir, status="RUNNING", phase="upstream")
    (run_dir / "provenance" / "run_provenance.yaml").write_text(": not [valid yaml", encoding="utf-8")
    _upstream_trace(run_dir).parent.mkdir(parents=True)
    _upstream_trace(run_dir).write_text("garbage without header\n", encoding="utf-8")
    (run_dir / "logs" / "upstream.stdout.log").write_bytes(b"\xff\xfe binary \x00\n")
    # A malformed state of another run is skipped, not fatal.
    broken = _run_dir(tmp_path, case="OTHER")
    (broken / "run_state.json").write_text("{", encoding="utf-8")

    output = _status(tmp_path)
    assert "Case:       grcm39-v130-salmon-r6" in output
    assert "Upstream:   RUNNING" in output


def test_status_with_only_a_run_state_and_no_runs(tmp_path):
    run_dir = _run_dir(tmp_path)
    _state(run_dir, status="CREATED")
    output = _status(tmp_path)
    assert "Phase:      validation" in output and "Upstream:   PENDING" in output

    empty = tmp_path / "empty"
    (empty / "planning").mkdir(parents=True)
    (empty / "planning" / "manifest.preview.yaml").write_text("{}", encoding="utf-8")
    assert "State:      PLANNED" in _status(empty)
    nothing = tmp_path / "nothing-yet"
    nothing.mkdir()
    assert "No recorded case runs." in _status(nothing)


def test_status_never_modifies_the_project(tmp_path):
    run_dir = _run_dir(tmp_path)
    _state(run_dir, status="RUNNING", phase="upstream")
    (run_dir / "logs" / "rnaseq.process.json").write_text(json.dumps({**current_process_identity(), "pid": _dead_pid()}), encoding="utf-8")
    _trace(_upstream_trace(run_dir), [("aa/000001", "FASTQC", "COMPLETED")])
    before = _tree_fingerprint(tmp_path)
    _status(tmp_path)
    _status(tmp_path, "--all")
    runner.invoke(app, ["status", str(tmp_path), "--watch", "--interval", "2"])  # final (stale) state: prints once
    assert _tree_fingerprint(tmp_path) == before


# ------------------------------------------------------------------ --watch

def test_watch_stops_cleanly_on_ctrl_c(monkeypatch, tmp_path):
    run_dir = _run_dir(tmp_path)
    _state(run_dir, status="RUNNING", phase="upstream")
    (run_dir / "logs" / "rnaseq.process.json").write_text(json.dumps(current_process_identity()), encoding="utf-8")

    def interrupt(_seconds):
        raise KeyboardInterrupt

    monkeypatch.setattr("rnaseq.cli.time.sleep", interrupt)
    result = runner.invoke(app, ["status", str(tmp_path), "--watch"])
    assert result.exit_code == 0, result.output
    assert "State:      RUNNING" in result.output
    assert "Stopped watching; the run itself was not affected." in result.output
    assert "Traceback" not in result.output


def test_watch_refreshes_until_the_run_reaches_a_final_state(monkeypatch, tmp_path):
    run_dir = _run_dir(tmp_path)
    _state(run_dir, status="RUNNING", phase="downstream")
    (run_dir / "logs" / "rnaseq.process.json").write_text(json.dumps(current_process_identity()), encoding="utf-8")
    sleeps: list[float] = []

    def finish(seconds):
        sleeps.append(seconds)
        _state(run_dir, status="SUCCESS", phase="delivery", completed_at="2026-09-25T22:00:00+08:00")

    monkeypatch.setattr("rnaseq.cli.time.sleep", finish)
    result = runner.invoke(app, ["status", str(tmp_path), "--watch", "--interval", "5"])
    assert result.exit_code == 0, result.output
    assert sleeps == [5.0]
    assert result.output.count("Case:") == 2
    assert result.output.rstrip().split("State:")[-1].startswith("      SUCCESS")


def test_watch_interval_has_a_floor(tmp_path):
    result = runner.invoke(app, ["status", str(tmp_path), "--watch", "--interval", "0.1"])
    assert result.exit_code != 0


# ------------------------------------------------------------------ per-run logs; same case ID twice

def test_two_runs_with_the_same_case_id_stay_distinguishable_and_logs_do_not_mix(monkeypatch, tmp_path):
    """The r6 scenario: a duplicate launch is SIGTERM'd while the official run succeeds."""

    from rnaseq.service import _Termination, execute_service_run

    report = salmon_case_project(tmp_path, monkeypatch)
    moments = iter([datetime(2026, 9, 25, 20, 30, 32), datetime(2026, 9, 25, 20, 31, 5)])
    from rnaseq import service
    original_create = service.create_case_run
    monkeypatch.setattr("rnaseq.service.create_case_run", lambda report, case_id: original_create(report, case_id, moment=next(moments)))

    def terminated(command, *, cwd, stdout_path, stderr_path):
        stdout_path.write_text("duplicate nextflow output\n", encoding="utf-8")
        raise _Termination("SIGTERM")

    monkeypatch.setattr("rnaseq.service._run_command", terminated)
    with pytest.raises(_Termination):
        execute_service_run(report, case_id="SAME-CASE")
    observed: list[list[str]] = []
    monkeypatch.setattr("rnaseq.service._run_command", fake_nextflow_factory(observed))
    official = execute_service_run(report, case_id="SAME-CASE")

    runs = sorted((tmp_path / "project" / "runs" / "SAME-CASE").iterdir())
    assert [run.name for run in runs] == ["20260925-203032+0800", "20260925-203105+0800"]
    duplicate, official_dir = runs
    assert official_dir == official.run_dir
    assert json.loads((duplicate / "run_state.json").read_text())["status"] == "INTERRUPTED"
    assert json.loads((duplicate / "run_state.json").read_text())["interrupted_by"] == "SIGTERM"
    assert json.loads((official_dir / "run_state.json").read_text())["status"] == "SUCCESS"

    duplicate_log = (duplicate / "logs" / "rnaseq.log").read_text(encoding="utf-8")
    official_log = (official_dir / "logs" / "rnaseq.log").read_text(encoding="utf-8")
    assert "SAME-CASE/20260925-203032+0800 started" in duplicate_log and "state INTERRUPTED" in duplicate_log
    assert "20260925-203105+0800" not in duplicate_log and "state SUCCESS" not in duplicate_log
    assert "SAME-CASE/20260925-203105+0800 started" in official_log and "run finished: SUCCESS" in official_log
    assert "20260925-203032+0800" not in official_log and "INTERRUPTED" not in official_log
    assert (duplicate / "logs" / "upstream.stdout.log").read_text() == "duplicate nextflow output\n"
    assert "duplicate" not in (official_dir / "logs" / "upstream.stdout.log").read_text()
    # Nothing nf-rna-owned is written at the project root.
    assert not list((tmp_path / "project").glob("*.log"))

    latest = _status(tmp_path / "project")
    assert "Run:        20260925-203105+0800" in latest and "State:      SUCCESS" in latest
    assert "Upstream:   SUCCESS\nTasks:      1/1 completed, 0 cached, 0 failed" in latest
    older = _status(tmp_path / "project", "--case", "SAME-CASE", "--run", "20260925-203032+0800")
    assert "Run:        20260925-203032+0800" in older and "State:      INTERRUPTED" in older
    # Delivery still excludes internal logs.
    assert not any(path.name.startswith("rnaseq") for path in (official_dir / "delivery").rglob("*"))


def test_ctrl_c_during_a_run_is_recorded_as_interrupted(monkeypatch, tmp_path):
    report = salmon_case_project(tmp_path, monkeypatch)

    def interrupted(command, *, cwd, stdout_path, stderr_path):
        raise KeyboardInterrupt

    monkeypatch.setattr("rnaseq.service._run_command", interrupted)
    result = runner.invoke(app, ["run", str(tmp_path / "project"), "--case-id", "CASE-CTRL-C", "--yes"])
    assert result.exit_code == 130, result.output
    assert "INTERRUPTED: rnaseq received SIGINT" in result.output
    state_path = next((tmp_path / "project" / "runs" / "CASE-CTRL-C").glob("*/run_state.json"))
    state = json.loads(state_path.read_text(encoding="utf-8"))
    assert state["status"] == "INTERRUPTED" and state["interrupted_by"] == "SIGINT" and state["completed_at"]


def test_unexpected_error_is_recorded_as_failed_with_a_traceback_in_the_run_log(monkeypatch, tmp_path):
    from rnaseq.service import execute_service_run

    report = salmon_case_project(tmp_path, monkeypatch)

    def broken(command, *, cwd, stdout_path, stderr_path):
        raise KeyError("unexpected")

    monkeypatch.setattr("rnaseq.service._run_command", broken)
    with pytest.raises(KeyError):
        execute_service_run(report, case_id="CASE-BUG")
    run_dir = next((tmp_path / "project" / "runs" / "CASE-BUG").iterdir())
    state = json.loads((run_dir / "run_state.json").read_text(encoding="utf-8"))
    assert state["status"] == "FAILED" and "Unexpected KeyError" in state["error"]
    log = (run_dir / "logs" / "rnaseq.log").read_text(encoding="utf-8")
    assert "Traceback" in log and "KeyError: 'unexpected'" in log


def test_successful_run_log_records_lifecycle_and_process_identity(monkeypatch, tmp_path):
    from rnaseq.service import execute_service_run

    report = salmon_case_project(tmp_path, monkeypatch)
    monkeypatch.setattr("rnaseq.service._run_command", fake_nextflow_factory([]))
    run = execute_service_run(report, case_id="CASE-LOG")
    log = (run.run_dir / "logs" / "rnaseq.log").read_text(encoding="utf-8")
    for expected in (
        f"run CASE-LOG/{run.run_id} started by pid {os.getpid()}",
        "state RUNNING phase=upstream", "upstream Nextflow launched", "upstream Nextflow exited with code 0",
        "state RUNNING phase=downstream", "downstream Nextflow exited with code 0",
        "state RUNNING phase=delivery", "state SUCCESS phase=delivery", "run finished: SUCCESS",
    ):
        assert expected in log
    record = json.loads((run.run_dir / "logs" / "rnaseq.process.json").read_text(encoding="utf-8"))
    assert record["pid"] == os.getpid() and record["command"][:2] == ["rnaseq", "run"]
    # After the run, status reports SUCCESS from durable state even though this process is alive.
    status = build_status(*select_run(tmp_path / "project"), cache=_Cache())
    assert status.state == "SUCCESS" and status.phase == "completed"
    assert "Resources:  14 CPUs / 54 GiB" in render_status(status)


def test_hisat2_route_records_an_upstream_trace_that_status_reads(tmp_path):
    from rnaseq.service import CaseRun, write_upstream_observer_config

    run_dir = _run_dir(tmp_path, case="H2-CASE")
    config = write_upstream_observer_config(CaseRun("H2-CASE", RUN_ID, run_dir, "2026-09-25T20:30:32+08:00"))
    text = config.read_text(encoding="utf-8")
    assert "trace {" in text and str((run_dir / "provenance" / "upstream.trace.txt").resolve()) in text
    assert "report" not in text and "timeline" not in text
    _state(run_dir, status="SUCCESS", phase="delivery", completed_at="2026-09-25T21:00:00+08:00")
    _trace(run_dir / "provenance" / "upstream.trace.txt", [("aa/000001", "HISAT2_ALIGN (S1)", "COMPLETED"), ("aa/000002", "FEATURECOUNTS (S1)", "COMPLETED")])
    assert "Upstream:   SUCCESS\nTasks:      2/2 completed, 0 cached, 0 failed" in _status(tmp_path)


def test_raw_count_runs_report_upstream_not_applicable(tmp_path):
    run_dir = _run_dir(tmp_path, case="RAW")
    _state(run_dir, status="SUCCESS", phase="delivery", completed_at="2026-09-25T20:40:00+08:00")
    (run_dir / "frozen" / "input_manifest.yaml").write_text(yaml.safe_dump({"input": {"type": "raw_counts"}}), encoding="utf-8")
    assert "Upstream:   NOT APPLICABLE (raw-count input)" in _status(tmp_path)


def test_real_sigterm_during_upstream_records_interrupted_and_stops_nextflow(monkeypatch, tmp_path):
    import signal
    import time

    from rnaseq.service import _Termination, execute_service_run

    report = salmon_case_project(tmp_path, monkeypatch)
    script = tmp_path / "fake-nextflow"
    script.write_text(
        "#!/bin/sh\n"
        f"trap 'echo terminated > {tmp_path}/terminated; exit 143' TERM\n"
        f"echo $$ > {tmp_path}/started\n"
        "while :; do sleep 0.1; done\n",
        encoding="utf-8",
    )
    script.chmod(0o755)
    monkeypatch.setattr("rnaseq.service.build_nextflow_command", lambda *_args, **_kwargs: [str(script)])
    child = os.fork()
    if child == 0:
        try:
            signal.signal(signal.SIGTERM, signal.SIG_DFL)
            execute_service_run(report, case_id="CASE-SIGTERM")
            os._exit(0)
        except _Termination:
            os._exit(143)
        except BaseException:
            os._exit(99)
    deadline = time.monotonic() + 30
    while not (tmp_path / "started").exists() and time.monotonic() < deadline:
        time.sleep(0.05)
    nextflow_pid = int((tmp_path / "started").read_text(encoding="utf-8"))
    os.kill(child, signal.SIGTERM)
    _pid, status = os.waitpid(child, 0)
    assert os.WEXITSTATUS(status) == 143
    assert (tmp_path / "terminated").read_text(encoding="utf-8").strip() == "terminated"
    with pytest.raises(ProcessLookupError):
        os.kill(nextflow_pid, 0)
    run_dir = next((tmp_path / "project" / "runs" / "CASE-SIGTERM").iterdir())
    state = json.loads((run_dir / "run_state.json").read_text(encoding="utf-8"))
    assert state["status"] == "INTERRUPTED" and state["interrupted_by"] == "SIGTERM" and state["phase"] == "upstream"
    assert "upstream Nextflow launched" in (run_dir / "logs" / "rnaseq.log").read_text(encoding="utf-8")
    output = _status(tmp_path / "project")
    assert "State:      INTERRUPTED\nPhase:      upstream" in output
