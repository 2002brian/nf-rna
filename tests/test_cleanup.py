"""`rnaseq clean`: explicit, provenance-resolved removal of one run's Nextflow work."""

from __future__ import annotations

import json
import os
from pathlib import Path

import pytest
import yaml
from typer.testing import CliRunner

from rnaseq.cleanup import plan_work_cleanup
from rnaseq.cli import app
from rnaseq.errors import ExecutionPreflightError

runner = CliRunner()
CASE = "CASE-001"
RUN = "20260929-123856+0800"
DURABLE = ("delivery/report.html", "logs/rnaseq.log", "provenance/upstream.trace.txt", "frozen/project.yaml", "downstream/l2/all_genes.tsv")


@pytest.fixture
def roots(tmp_path, monkeypatch):
    execution = tmp_path / "exec"
    monkeypatch.setenv("RNASEQ_EXECUTION_ROOT", str(execution))
    monkeypatch.setenv("RNASEQ_RUNTIME_ROOT", str(tmp_path / "runtime"))
    return tmp_path, execution


def _write_tree(root: Path) -> None:
    for relative in ("upstream/ab/cdef12/.command.sh", "upstream/ab/cdef12/out.sam", "downstream/0f/99aa/vst.tsv"):
        path = root / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text("x" * 4096, encoding="utf-8")


def make_run(
    tmp_path: Path, execution: Path, *, status: str = "SUCCESS", work_root: Path | None = None,
    provenance: dict | None = None, commands: bool = True, case: str = CASE, run: str = RUN, create_work: bool = True,
) -> tuple[Path, Path, Path]:
    project = tmp_path / "project"
    run_dir = project / "runs" / case / run
    for relative in DURABLE:
        (run_dir / relative).parent.mkdir(parents=True, exist_ok=True)
        (run_dir / relative).write_text("durable", encoding="utf-8")
    exec_root = (execution / case / run).resolve()
    work = (work_root / case / run / "work").resolve() if work_root else exec_root / "work"
    (exec_root / "launch").mkdir(parents=True, exist_ok=True)
    (exec_root / "launch" / ".nextflow.log").write_text("log", encoding="utf-8")
    if create_work:
        _write_tree(work)
    state = {"case_id": case, "run_id": run, "status": status}
    if commands:
        state["upstream_command"] = ["nextflow", "run", "x", "-work-dir", str(work / "upstream")]
        state["downstream_command"] = ["nextflow", "run", "y", "-work-dir", str(work / "downstream")]
    (run_dir / "run_state.json").write_text(json.dumps(state), encoding="utf-8")
    record = provenance if provenance is not None else {
        "execution_root": str(exec_root), "execution_launch_dir": str(exec_root / "launch"), "execution_work_dir": str(work),
    }
    if record:
        (run_dir / "provenance" / "run_provenance.yaml").write_text(yaml.safe_dump(record), encoding="utf-8")
    return project, run_dir, work


def _snapshot(path: Path) -> dict[str, str]:
    return {str(item.relative_to(path)): item.read_text(encoding="utf-8") for item in sorted(path.rglob("*")) if item.is_file()}


def clean(project: Path, *extra: str, case: str = CASE, run: str = RUN):
    return runner.invoke(app, ["clean", str(project), "--case-id", case, "--run", run, "--yes", *extra])


def test_success_run_work_in_relocated_root_is_removed_and_durable_record_untouched(roots):
    tmp_path, execution = roots
    work_root = tmp_path / "nf-rna-work"
    project, run_dir, work = make_run(tmp_path, execution, work_root=work_root)
    before = _snapshot(run_dir)
    result = clean(project)
    assert result.exit_code == 0, result.output
    assert not work.exists()
    assert "Work size: about" in result.output and "Cleaned: removed" in result.output
    assert f"Durable run record preserved: {run_dir.resolve()}" in result.output
    after = _snapshot(run_dir)
    record = json.loads(after.pop("logs/work_cleanup.json"))
    assert after == before
    assert record[0]["work_dir"] == str(work) and record[0]["recorded_state"] == "SUCCESS" and record[0]["forced"] is False
    # Only empty per-run/per-case shells are tidied; the work root itself and the execution root stay.
    assert work_root.is_dir() and not (work_root / CASE).exists()
    assert (execution / CASE / RUN / "launch" / ".nextflow.log").is_file()


def test_default_legacy_work_location_keeps_execution_root_and_launch_dir(roots):
    tmp_path, execution = roots
    project, _run_dir, work = make_run(tmp_path, execution)
    result = clean(project)
    assert result.exit_code == 0, result.output
    assert not work.exists()
    assert (execution / CASE / RUN).is_dir()
    assert (execution / CASE / RUN / "launch" / ".nextflow.log").is_file()


def test_missing_work_directory_is_an_idempotent_no_op(roots):
    tmp_path, execution = roots
    project, _run_dir, _work = make_run(tmp_path, execution, create_work=False)
    for _ in range(2):
        result = clean(project)
        assert result.exit_code == 0, result.output
        assert "Nothing to clean" in result.output
    project2, _run_dir2, _work2 = make_run(tmp_path / "second", execution)
    assert clean(project2).exit_code == 0
    again = clean(project2)
    assert again.exit_code == 0 and "Nothing to clean" in again.output


@pytest.mark.parametrize("status", ["FAILED", "INTERRUPTED"])
def test_failed_and_interrupted_work_is_retained_unless_forced(roots, status):
    tmp_path, execution = roots
    project, run_dir, work = make_run(tmp_path, execution, status=status)
    refused = clean(project)
    assert refused.exit_code == 1
    assert "retained by default" in refused.output and "--force" in refused.output
    assert work.is_dir()
    forced = clean(project, "--force")
    assert forced.exit_code == 0, forced.output
    assert not work.exists()
    assert json.loads((run_dir / "logs" / "work_cleanup.json").read_text())[0]["forced"] is True


@pytest.mark.parametrize("status", ["RUNNING", "CREATED"])
def test_running_run_is_refused_even_with_force(roots, status):
    tmp_path, execution = roots
    project, run_dir, work = make_run(tmp_path, execution, status=status)
    from rnaseq.run_status import current_process_identity

    (run_dir / "logs" / "rnaseq.process.json").write_text(json.dumps(current_process_identity()), encoding="utf-8")
    result = clean(project, "--force")
    assert result.exit_code == 1
    assert "RUNNING" in result.output and "never removed" in result.output
    assert work.is_dir()


def test_running_without_verifiable_process_is_refused(roots):
    tmp_path, execution = roots
    project, _run_dir, work = make_run(tmp_path, execution, status="RUNNING")
    result = clean(project, "--force")
    assert result.exit_code == 1 and work.is_dir()


def test_stale_running_run_counts_as_interrupted_and_needs_force(roots):
    tmp_path, execution = roots
    project, run_dir, work = make_run(tmp_path, execution, status="RUNNING")
    (run_dir / "logs" / "rnaseq.process.json").write_text(
        json.dumps({"pid": 2**22 + 12345, "hostname": os.uname().nodename}), encoding="utf-8",
    )
    refused = clean(project)
    assert refused.exit_code == 1 and "INTERRUPTED (recorded RUNNING)" in refused.output
    assert clean(project, "--force").exit_code == 0
    assert not work.exists()


def test_dry_run_and_declined_confirmation_remove_nothing(roots):
    tmp_path, execution = roots
    project, _run_dir, work = make_run(tmp_path, execution)
    dry = clean(project, "--dry-run")
    assert dry.exit_code == 0 and "Dry run" in dry.output and work.is_dir()
    declined = runner.invoke(app, ["clean", str(project), "--case", CASE, "--run", RUN], input="n\n")
    assert declined.exit_code == 0 and "cancelled" in declined.output and work.is_dir()


@pytest.mark.parametrize(("case", "run"), [
    ("..", RUN), ("../CASE-001", RUN), ("CASE-001/..", RUN), (CASE, "../" + RUN), (CASE, "not-a-run"), (CASE, RUN + "/.."),
])
def test_crafted_case_or_run_arguments_are_rejected(roots, case, run):
    tmp_path, execution = roots
    project, _run_dir, work = make_run(tmp_path, execution)
    result = clean(project, case=case, run=run)
    assert result.exit_code == 1, result.output
    assert work.is_dir()


@pytest.mark.parametrize("recorded", [
    "relative/CASE-001/20260929-123856+0800/work",
    "/tmp/../etc/CASE-001/20260929-123856+0800/work",
    "/tmp/x/CASE-001/20260929-123856+0800/./work",
    "/tmp/x//CASE-001/20260929-123856+0800/work",
    "/tmp/x/OTHER-CASE/20260929-123856+0800/work",
    "/tmp/x/CASE-001/20260101-000000+0800/work",
    "/tmp/x/CASE-001/20260929-123856+0800/workdir",
    "/CASE-001/20260929-123856+0800/work",
    "",
])
def test_malformed_recorded_work_paths_are_refused(roots, recorded):
    tmp_path, execution = roots
    project, _run_dir, work = make_run(tmp_path, execution, provenance={"execution_work_dir": recorded}, commands=False)
    with pytest.raises(ExecutionPreflightError):
        plan_work_cleanup(project, CASE, RUN)
    assert clean(project).exit_code == 1
    assert work.is_dir()


def test_command_work_dirs_must_agree_with_provenance(roots):
    tmp_path, execution = roots
    project, run_dir, work = make_run(tmp_path, execution)
    state = json.loads((run_dir / "run_state.json").read_text())
    state["downstream_command"][-1] = str(tmp_path / "elsewhere" / "downstream")
    (run_dir / "run_state.json").write_text(json.dumps(state), encoding="utf-8")
    result = clean(project)
    assert result.exit_code == 1 and "does not match" in result.output
    assert work.is_dir()


def test_symlinked_work_directory_is_not_followed(roots):
    tmp_path, execution = roots
    victim = tmp_path / "victim"
    _write_tree(victim)
    work = tmp_path / "wr" / CASE / RUN / "work"
    work.parent.mkdir(parents=True)
    work.symlink_to(victim, target_is_directory=True)
    project, _run_dir, _ = make_run(
        tmp_path, execution, create_work=False, commands=False, provenance={"execution_work_dir": str(work)},
    )
    result = clean(project)
    assert result.exit_code == 1 and "symbolic link" in result.output
    assert (victim / "upstream" / "ab" / "cdef12" / "out.sam").is_file()


def test_symlinked_ancestor_of_work_directory_is_not_followed(roots):
    tmp_path, execution = roots
    real = tmp_path / "real-root"
    _write_tree(real / CASE / RUN / "work")
    (tmp_path / "linked-root").symlink_to(real, target_is_directory=True)
    work = tmp_path / "linked-root" / CASE / RUN / "work"
    project, _run_dir, _ = make_run(
        tmp_path, execution, create_work=False, commands=False, provenance={"execution_work_dir": str(work)},
    )
    assert clean(project).exit_code == 1
    assert (real / CASE / RUN / "work" / "upstream").is_dir()


def test_links_inside_work_are_removed_without_touching_their_targets(roots):
    tmp_path, execution = roots
    fastq = tmp_path / "input" / "S1_R1.fastq.gz"
    fastq.parent.mkdir()
    fastq.write_text("reads", encoding="utf-8")
    project, _run_dir, work = make_run(tmp_path, execution, work_root=tmp_path / "wr")
    (work / "upstream" / "ab" / "cdef12" / "S1_R1.fastq.gz").symlink_to(fastq)
    (work / "upstream" / "ab" / "cdef12" / "inputs").symlink_to(fastq.parent, target_is_directory=True)
    assert clean(project).exit_code == 0
    assert not work.exists() and fastq.read_text() == "reads"


def test_unrecognised_work_contents_are_refused(roots):
    tmp_path, execution = roots
    project, _run_dir, work = make_run(tmp_path, execution, work_root=tmp_path / "wr")
    (work / "important.txt").write_text("keep", encoding="utf-8")
    result = clean(project)
    assert result.exit_code == 1 and "important.txt" in result.output
    assert (work / "important.txt").is_file()


def test_durable_links_into_work_block_cleanup(roots):
    tmp_path, execution = roots
    project, run_dir, work = make_run(tmp_path, execution, work_root=tmp_path / "wr")
    (run_dir / "upstream").mkdir()
    (run_dir / "upstream" / "result.bam").symlink_to(work / "upstream" / "ab" / "cdef12" / "out.sam")
    result = clean(project)
    assert result.exit_code == 1 and "link into the work directory" in result.output
    assert work.is_dir()


def test_work_root_itself_cannot_be_the_target(roots, monkeypatch):
    tmp_path, execution = roots
    # A crafted record whose base is shaped like a work root that *is* <case>/<run>/work.
    work_root = tmp_path / "wr" / CASE / RUN / "work"
    monkeypatch.setenv("RNASEQ_WORK_ROOT", str(work_root))
    _write_tree(work_root)
    project, _run_dir, _ = make_run(tmp_path, execution, create_work=False, commands=False, provenance={"execution_work_dir": str(work_root)})
    result = clean(project)
    assert result.exit_code == 1 and "protected location" in result.output
    assert work_root.is_dir()


@pytest.mark.parametrize("target", ["execution_root", "project", "run_dir"])
def test_execution_root_project_and_run_dir_cannot_be_targets(roots, target):
    tmp_path, execution = roots
    project = tmp_path / "project"
    run_dir = project / "runs" / CASE / RUN
    locations = {
        # An execution root shaped exactly like <base>/<case>/<run>/work.
        "execution_root": execution / "x" / CASE / RUN / "work",
        "project": run_dir / CASE / RUN / "work",
        "run_dir": project / "runs" / CASE / RUN / "work",
    }
    recorded = locations[target].resolve()
    _write_tree(recorded)
    provenance = {"execution_work_dir": str(recorded)}
    if target == "execution_root":
        provenance["execution_root"] = str(recorded)
    make_run(tmp_path, execution, create_work=False, commands=False, provenance=provenance)
    result = clean(project)
    assert result.exit_code == 1, result.output
    assert recorded.is_dir()


def test_run_without_any_recorded_work_directory_is_refused(roots):
    tmp_path, execution = roots
    project, _run_dir, work = make_run(tmp_path, execution, provenance={}, commands=False)
    result = clean(project)
    assert result.exit_code == 1 and "does not guess" in result.output
    assert work.is_dir()


def test_old_run_without_provenance_key_uses_recorded_nextflow_commands(roots):
    tmp_path, execution = roots
    project, _run_dir, work = make_run(tmp_path, execution, provenance={"execution_root": "/unused"})
    plan = plan_work_cleanup(project, CASE, RUN)
    assert plan.work_dir == work and "run_state.json" in plan.source
    assert clean(project).exit_code == 0 and not work.exists()


def test_run_state_identity_must_match_the_requested_run(roots):
    tmp_path, execution = roots
    project, run_dir, work = make_run(tmp_path, execution)
    state = json.loads((run_dir / "run_state.json").read_text())
    state["case_id"] = "OTHER"
    (run_dir / "run_state.json").write_text(json.dumps(state), encoding="utf-8")
    assert clean(project).exit_code == 1 and work.is_dir()


def test_symlinked_run_directory_is_refused(roots):
    tmp_path, execution = roots
    project, run_dir, work = make_run(tmp_path, execution, case="REAL")
    (project / "runs" / CASE).symlink_to(project / "runs" / "REAL", target_is_directory=True)
    result = clean(project)
    assert result.exit_code == 1 and "not a real directory" in result.output
    assert work.is_dir()


def test_unknown_run_is_reported_without_traceback(roots):
    tmp_path, execution = roots
    project, _run_dir, _work = make_run(tmp_path, execution)
    result = clean(project, run="20250101-000000+0800")
    assert result.exit_code == 1 and "no recorded run" in result.output
    assert result.exception is None or isinstance(result.exception, SystemExit)


@pytest.mark.parametrize("status", ["SUCCESS", "INTERRUPTED"])
def test_live_process_inside_work_blocks_cleanup_even_with_force(roots, status):
    import subprocess

    tmp_path, execution = roots
    project, _run_dir, work = make_run(tmp_path, execution, status=status, work_root=tmp_path / "wr")
    task = work / "upstream" / "ab" / "cdef12"
    orphan = subprocess.Popen(["sleep", "60"], cwd=task)  # e.g. a Nextflow task that outlived a killed rnaseq
    try:
        result = clean(project, "--force")
        assert result.exit_code == 1 and "live process" in result.output and str(orphan.pid) in result.output
        assert work.is_dir()
    finally:
        orphan.kill()
        orphan.wait()
    assert clean(project, "--force").exit_code == 0 and not work.exists()


def test_symlink_loop_is_refused_without_traceback(roots):
    tmp_path, execution = roots
    work = tmp_path / "wr" / CASE / RUN / "work"
    work.parent.mkdir(parents=True)
    work.symlink_to(work)
    project, _run_dir, _ = make_run(tmp_path, execution, create_work=False, commands=False, provenance={"execution_work_dir": str(work)})
    result = clean(project)
    assert result.exit_code == 1 and "cannot be resolved safely" in result.output
    assert work.is_symlink()


def test_work_inside_shared_cache_is_refused(roots):
    tmp_path, execution = roots
    work = (execution / "cache" / "upstream-conda" / CASE / RUN / "work").resolve()
    _write_tree(work)
    project, _run_dir, _ = make_run(tmp_path, execution, create_work=False, commands=False, provenance={"execution_work_dir": str(work)})
    result = clean(project)
    assert result.exit_code == 1 and "shared location" in result.output
    assert work.is_dir()


def test_work_removed_concurrently_reports_nothing_to_clean(roots, monkeypatch):
    import shutil as _shutil

    import rnaseq.cleanup as cleanup

    tmp_path, execution = roots
    project, _run_dir, work = make_run(tmp_path, execution, work_root=tmp_path / "wr")
    original = cleanup.measure_work

    def measure_then_vanish(path):
        usage = original(path)
        _shutil.rmtree(work)
        return usage

    monkeypatch.setattr(cleanup, "measure_work", measure_then_vanish)
    result = clean(project)
    assert result.exit_code == 0 and "disappeared before deletion" in result.output and "Cleaned:" not in result.output
