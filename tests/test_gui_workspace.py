"""Local Research Workspace (v1.7.0 Phase 7a): datasets, project overview, execution status and scientific QC.

Service tests need no GUI extra; the page tests use Streamlit's AppTest and
are skipped without it. Run fixtures are small synthetic run directories with
the records the pipeline writes; nothing here runs an analysis.
"""

from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path

import pytest
import yaml

from rnaseq.gui import datasets, workspace
from rnaseq.gui.config import PathNotAllowed, load_config
from rnaseq.models import DesignType, InputType, Preset, SequencingLayout, Species
from rnaseq.project import create_project, import_count_matrix
from rnaseq.run_status import current_process_identity
from rnaseq.validators import validate_project

from conftest import BASE_COUNTS


@pytest.fixture
def gui_home(tmp_path, monkeypatch):
    data = tmp_path / "data"
    (data / "projects").mkdir(parents=True)
    settings = tmp_path / "gui.yaml"
    settings.write_text(yaml.safe_dump({"allowed_roots": [str(data)], "project_parent": str(data / "projects"), "condarc": None, "launcher": "process"}), encoding="utf-8")
    monkeypatch.setenv("NF_RNA_GUI_CONFIG", str(settings))
    monkeypatch.setenv("XDG_CACHE_HOME", str(tmp_path / "cache"))
    return data


def _write(path: Path, text: str | bytes) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    (path.write_bytes if isinstance(text, bytes) else path.write_text)(text)
    return path


def _raw_project(parent: Path, name: str, *, metadata: str | None = None, contrasts: str | None = None) -> Path:
    project = create_project(
        project_name=name, destination=parent, species=Species.HUMAN, preset=Preset.L2, design_type=DesignType.TWO_GROUP,
        input_type=InputType.RAW_COUNTS, layout=None, scaffold=True, gene_id_type="ENSEMBL",
    )
    counts = parent / f"{name}-counts.csv"
    counts.write_text(BASE_COUNTS, encoding="utf-8")
    import_count_matrix(project, counts)
    if metadata is not None:
        (project / "metadata.csv").write_text(metadata, encoding="utf-8")
    if contrasts is not None:
        (project / "contrasts.csv").write_text(contrasts, encoding="utf-8")
    return project


def _run(project: Path, case: str, *, status: str = "SUCCESS", input_type: str = "raw_counts", strand: dict | None = None,
         flags: list[str] | None = None, r_warnings: dict | None = None, started: str = "2026-10-30T12:00:00+08:00", level: str = "L2") -> Path:
    run = project / "runs" / case / "20261030-120000+0800"
    _write(run / "run_state.json", json.dumps({"case_id": case, "run_id": run.name, "status": status, "started_at": started, **({"error": "R exited 1"} if status == "FAILED" else {})}))
    contract: dict = {"analysis_level": level}
    if strand is not None:
        contract["upstream_qc"] = {"strandedness": strand}
    _write(run / "frozen" / "downstream_contract.json", json.dumps(contract))
    _write(run / "frozen" / "project.yaml", yaml.safe_dump({"input": {"type": input_type}}))
    if flags is not None:
        _write(run / "downstream" / "l1" / "backend_summary.json", json.dumps({"flags": flags}))
    if r_warnings is not None:
        metadata = {"shrinkage": {"method": "apeglm", "status": "success"}, "r_warnings": r_warnings}
        _write(run / "downstream" / "l2" / "contrasts" / "B_vs_A" / "contrast_metadata.json", json.dumps(metadata))
    return run


def _digest(directory: Path) -> dict[str, str]:
    return {str(path.relative_to(directory)): hashlib.sha256(path.read_bytes()).hexdigest() for path in sorted(directory.rglob("*")) if path.is_file()}


F1_WARNING = {"status": "WARNING", "declared": "forward", "samples": {
    "D5_R_A": {"status": "CONTRADICTS_DECLARATION", "forward_over_reverse": 0.096},
    "D5_R_B": {"status": "CONTRADICTS_DECLARATION", "forward_over_reverse": 0.075},
}}
F3_WARNINGS = {"status": "WARNINGS", "fit": [], "contrast": [
    {"step": "shrinkage: lfcShrink(type = apeglm)", "message": "the line search routine failed, unable to sufficiently decrease the function value", "count": 6},
]}


# --------------------------------------------------------------------- datasets

def test_dataset_scan_recognises_kinds_sizes_and_permissions(gui_home):
    folder = gui_home / "fastq" / "batch1"
    _write(folder / "S1_R1.fastq.gz", b"x" * 10)
    _write(folder / "S1_R2.fastq.gz", b"x" * 10)
    _write(folder / "S2_R1.fastq.gz", b"x" * 12)
    _write(folder / "S2_R2.fastq.gz", b"x" * 12)
    _write(folder / "counts.csv", "gene_id,S1,S2\nG1,1,2\n")
    _write(folder / "metadata.csv", "sample_id,condition\nS1,A\nS2,B\n")
    _write(folder / "samplesheet.csv", "sample,fastq_1,fastq_2,strandedness\n")
    _write(folder / "notes.txt", "free text\n")
    locked = _write(folder / "locked.csv", "sample_id\n")
    locked.chmod(0)
    try:
        scan = datasets.scan_dataset(load_config(), folder)
    finally:
        locked.chmod(0o644)
    kinds = {item.relative: item.kind for item in scan.files}
    assert kinds["S1_R1.fastq.gz"] == "fastq" and kinds["counts.csv"] == "count_matrix"
    assert kinds["metadata.csv"] == "metadata" and kinds["samplesheet.csv"] == "samplesheet" and kinds["notes.txt"] == "table"
    assert len(scan.of_kind("fastq")) == 4 and scan.fastq_folders == [folder.resolve()]
    if os.geteuid() != 0:
        assert [item.relative for item in scan.unreadable] == ["locked.csv"]
    summary, errors = datasets.recognise_fastq(folder, datasets.guess_layout(scan.of_kind("fastq")))
    assert not errors and summary.sample_ids == ("S1", "S2") and summary.layout is SequencingLayout.PAIRED_END
    counts, issues = datasets.recognise_counts(folder / "counts.csv")
    assert counts is not None and counts.sample_ids == ("S1", "S2")


def test_dataset_scan_never_follows_a_symlink_out_of_the_allowed_roots(gui_home, tmp_path):
    outside = tmp_path / "private"
    _write(outside / "secret_R1.fastq.gz", b"x")
    folder = gui_home / "linked"
    folder.mkdir()
    (folder / "escape").symlink_to(outside)
    (folder / "secret_R1.fastq.gz").symlink_to(outside / "secret_R1.fastq.gz")
    config = load_config()
    scan = datasets.scan_dataset(config, folder)
    assert scan.files == [] and sorted(scan.outside_links) == ["escape", "secret_R1.fastq.gz"]
    with pytest.raises(PathNotAllowed):
        datasets.scan_dataset(config, outside)
    with pytest.raises(PathNotAllowed):
        datasets.scan_dataset(config, folder / "escape")
    sheet = _write(gui_home / "sheet.csv", f"sample,fastq_1,fastq_2,strandedness\nsecret,{outside / 'secret_R1.fastq.gz'},,auto\n")
    view = datasets.recognise_samplesheet(config, sheet)
    assert view.error is None and view.outside == (str((outside / "secret_R1.fastq.gz").resolve()),)


def test_fastq_copy_plan_reports_size_and_free_space(gui_home):
    folder = gui_home / "fq"
    _write(folder / "A_R1.fastq.gz", b"x" * 100)
    _write(folder / "A_R2.fastq.gz", b"x" * 50)
    summary, _ = datasets.recognise_fastq(folder, SequencingLayout.PAIRED_END)
    plan = datasets.fastq_copy_plan(summary.records, gui_home / "projects" / "p" / "input" / "fastq")
    assert (plan.files, plan.bytes) == (2, 150) and plan.free_bytes and plan.fits


# --------------------------------------------------------------------- execution status and scientific QC

def test_execution_status_has_four_states():
    assert [workspace.execution_status(state) for state in ("RUNNING", "PLANNED", "SUCCESS", "FAILED", "INCOMPLETE", "INTERRUPTED", "UNKNOWN")] == [
        "RUNNING", "RUNNING", "SUCCESS", "FAILED", "INCOMPLETE", "INCOMPLETE", "INCOMPLETE"]


def test_scientific_qc_reads_recorded_f1_and_f3_verdicts_and_writes_nothing(tmp_path):
    project = tmp_path / "proj"
    f1 = _run(project, "F1", input_type="fastq", strand=F1_WARNING, flags=[], level="L1")
    f3 = _run(project, "F3", flags=[], r_warnings=F3_WARNINGS)
    clean = _run(project, "CLEAN", input_type="fastq", strand={"status": "CONSISTENT", "declared": "reverse", "samples": {}}, flags=[], r_warnings={"status": "NONE", "fit": [], "contrast": []})
    legacy = _run(project, "LEGACY", flags=[], r_warnings=None)
    _write(legacy / "downstream" / "l2" / "contrasts" / "B_vs_A" / "contrast_metadata.json", json.dumps({"shrinkage": {"status": "success"}}))
    salmon = _run(project, "SALMON", input_type="fastq", strand={"status": "NOT_EVALUATED", "reason": "Salmon route"}, flags=[])
    inconclusive = _run(project, "LOWREADS", input_type="fastq", strand={"status": "INCONCLUSIVE", "samples": {"S": {"status": "INCONCLUSIVE"}}}, flags=[])
    failed_early = _run(project, "FAILED", status="FAILED", input_type="fastq")
    before = _digest(project)

    qc = workspace.run_qc(f1)
    assert qc.level == "WARNING" and qc.items[0].recorded == "WARNING" and "D5_R_A, D5_R_B" in qc.items[0].detail
    qc = workspace.run_qc(f3)
    assert qc.level == "WARNING" and qc.items[-1].recorded == "WARNINGS: 6" and "apeglm" in qc.items[-1].detail
    assert workspace.run_qc(clean).level == "PASS"
    assert [item.level for item in workspace.run_qc(legacy).items] == ["PASS", "NOT_EVALUATED"] and workspace.run_qc(legacy).level == "PASS"
    assert workspace.run_qc(salmon).items[0].level == "NOT_EVALUATED"
    assert workspace.run_qc(inconclusive).level == "INCONCLUSIVE"
    assert workspace.run_qc(failed_early).level == "NOT_EVALUATED"
    assert _digest(project) == before


def test_run_history_separates_execution_status_from_qc(tmp_path):
    project = tmp_path / "proj"
    _run(project, "OK", flags=[], r_warnings=F3_WARNINGS, started="2026-10-30T10:00:00+08:00")
    _run(project, "BROKEN", status="FAILED", started="2026-10-30T11:00:00+08:00")
    stale = _run(project, "KILLED", status="RUNNING", flags=[], started="2026-10-30T12:00:00+08:00")
    _write(stale / "logs" / "rnaseq.process.json", json.dumps({**current_process_identity(), "start_ticks": -1}))
    before = (stale / "run_state.json").read_bytes()
    rows = workspace.project_runs(project)
    assert [(row.case_id, row.execution, row.qc.level) for row in rows] == [("KILLED", "INCOMPLETE", "PASS"), ("BROKEN", "FAILED", "NOT_EVALUATED"), ("OK", "SUCCESS", "WARNING")]
    assert (stale / "run_state.json").read_bytes() == before


# --------------------------------------------------------------------- project summary and workflow

def test_project_summary_comes_from_the_project_files(gui_home):
    project = _raw_project(gui_home / "projects", "alpha", metadata="sample_id,condition\nS1,A\nS2,A\nS3,B\nS4,B\n",
                           contrasts="contrast_id,factor,numerator,denominator\nB_vs_A,condition,B,A\n")
    _run(project, "CASE-1", flags=[], r_warnings=F3_WARNINGS)
    summary = workspace.project_summary(project)
    assert (summary.project_id, summary.species, summary.input_type, summary.scope, summary.samples, summary.contrasts) == (
        "alpha", "Homo sapiens", "raw_counts", "L2", 4, ("B_vs_A",))
    assert summary.quantification.startswith("none") and summary.reference.startswith("not used")
    assert summary.latest.execution == "SUCCESS" and summary.latest.qc.headline() == "WARNING (1)"
    broken = gui_home / "projects" / "broken"
    _write(broken / "project.yaml", "schema_version: '9'\n")
    assert workspace.project_summary(broken).error


def test_workflow_steps_point_to_the_next_unfinished_step(gui_home):
    project = _raw_project(gui_home / "projects", "beta")
    steps = workspace.workflow_steps(validate_project(project), [])
    assert [step.name for step in steps] == list(workspace.STEP_NAMES)
    assert steps[0].status == "done" and workspace.next_step(steps).name == "Metadata"
    (project / "metadata.csv").write_text("sample_id,condition\nC1,control\nC2,control\nC3,control\nT1,treated\nT2,treated\nT3,treated\n", encoding="utf-8")
    steps = workspace.workflow_steps(validate_project(project), [])
    assert workspace.next_step(steps).name == "Contrasts"


# --------------------------------------------------------------------- pages (AppTest)

def _page(name: str):
    pytest.importorskip("streamlit")
    from streamlit.testing.v1 import AppTest

    return AppTest.from_string(f"from rnaseq.gui import pages\npages.{name}()\n", default_timeout=60)


def test_datasets_page_lists_and_recognises_a_fastq_folder(gui_home):
    folder = gui_home / "fastq" / "run1"
    for sample in ("C1", "C2"):
        for read in ("R1", "R2"):
            _write(folder / f"{sample}_{read}.fastq.gz", b"x" * 20)
    page = _page("datasets_page").run()
    page.text_input(key="ds_text").set_value(str(folder)).run()
    assert not page.exception, page.exception
    assert [metric.value for metric in page.metric][:3] == ["4", "80 B", "4"]
    sample_tables = [frame.value for frame in page.dataframe if "sample_id" in frame.value.columns]
    assert sample_tables and sample_tables[0]["sample_id"].tolist() == ["C1", "C2"]
    assert any("choose a project" in item.value for item in page.info)


def test_project_page_and_switching_keep_each_projects_own_state(gui_home):
    alpha = _raw_project(gui_home / "projects", "alpha", metadata="sample_id,condition\nS1,A\nS2,A\nS3,B\nS4,B\n")
    beta = _raw_project(gui_home / "projects", "beta")
    _run(alpha, "ALPHA-RUN", flags=[], r_warnings=F3_WARNINGS)
    page = _page("project_page")
    page.session_state["project"] = str(alpha)
    page.run()
    assert not page.exception, page.exception
    assert page.title[0].value == "alpha"
    history = [frame.value for frame in page.dataframe if "execution" in frame.value.columns][0]
    assert history["run"].tolist() == ["ALPHA-RUN/20261030-120000+0800"] and "WARNING" in history["scientific QC"].iloc[0]
    assert any("Scientific QC: **WARNING**" in item.value for item in page.warning)
    first, second = "ALPHA-RUN/20261030-120000+0800", "ALPHA-OLD/20261030-120000+0800"
    _run(alpha, "ALPHA-OLD", status="FAILED", started="2026-10-29T12:00:00+08:00")
    page.run()
    page.selectbox(key=f"project_run::{alpha}").set_value(second).run()
    assert any("NOT EVALUATED" in item.value for item in page.caption)

    page.session_state["project"] = str(beta)
    page.run()
    assert not page.exception and page.title[0].value == "beta"
    assert not [frame for frame in page.dataframe if "execution" in frame.value.columns]
    assert any("Next: **Metadata**" in item.value for item in page.info)

    page.session_state["project"] = str(alpha)
    page.run()
    assert page.title[0].value == "alpha" and page.selectbox(key=f"project_run::{alpha}").value == second
    assert first != second


def test_run_page_shows_validation_by_step_and_qc_of_runs(gui_home):
    project = _raw_project(gui_home / "projects", "gamma")
    other = _raw_project(gui_home / "projects", "delta")
    _run(project, "OLD", status="FAILED")
    page = _page("runner_page")
    page.session_state["project"] = str(project)
    page.run()
    assert not page.exception, page.exception
    assert any("cannot run yet" in item.value for item in page.error)
    assert any(item.value == "**Metadata**" for item in page.markdown)
    assert next(button for button in page.button if button.key == "run_start").disabled
    assert any("FAILED: R exited 1" in item.value for item in page.error)
    page.text_input(key=f"case::{project}").set_value("GAMMA-1").run()
    page.session_state["project"] = str(other)
    page.run()
    assert page.text_input(key=f"case::{other}").value == ""
    page.session_state["project"] = str(project)
    page.run()
    assert page.text_input(key=f"case::{project}").value == "GAMMA-1"


def test_a_project_opened_by_url_survives_a_page_switch(gui_home):
    project = _raw_project(gui_home / "projects", "epsilon")
    page = _page("project_page")
    page.query_params["project"] = str(project)
    page.run()
    assert not page.exception and page.title[0].value == "epsilon"
    assert page.session_state["project"] == str(project)  # kept when switch_page drops the query string


def test_datasets_page_assigns_a_fastq_folder_to_the_current_project_by_copy(gui_home):
    folder = gui_home / "fastq" / "batch2"
    for sample in ("Ctrl1", "Trt1"):
        for read in ("R1", "R2"):
            _write(folder / f"{sample}_{read}.fastq.gz", b"@r\nACGT\n+\nIIII\n")
    _write(folder / "metadata.csv", "sample_id,condition\nCtrl1,control\nTrt1,treated\n")
    project = create_project(project_name="fq", destination=gui_home / "projects", species=Species.HUMAN, preset=Preset.L1,
                             design_type=DesignType.TWO_GROUP, scaffold=True, gene_id_type="ENSEMBL")
    before = _digest(folder)
    page = _page("datasets_page")
    page.session_state["project"] = str(project)
    page.run()
    page.text_input(key="ds_text").set_value(str(folder)).run()
    assert not page.exception, page.exception
    assert any("Importing copies 4 FASTQ file(s)" in item.value for item in page.info)
    page.button(key="ds_fq_import").click().run()
    assert not page.exception and not page.error, [item.value for item in page.error]
    assert sorted(path.name for path in (project / "input" / "fastq").iterdir()) == ["Ctrl1_R1.fastq.gz", "Ctrl1_R2.fastq.gz", "Trt1_R1.fastq.gz", "Trt1_R2.fastq.gz"]
    assert workspace.input_origin(project) == [str(folder.resolve())]
    assert _digest(folder) == before  # the originals are not moved or changed


def test_reuse_hint_refuses_a_source_run_of_another_nf_rna_version(tmp_path):
    from types import SimpleNamespace

    from rnaseq.models import PIPELINE_VERSION

    project = tmp_path / "fq"
    for case, version in (("OLD", "1.6.1"), ("NEW", PIPELINE_VERSION)):
        run = project / "runs" / case / "20261031-120000+0800"
        (run / "handoff").mkdir(parents=True)
        (run / "frozen").mkdir()
        (run / "handoff" / "upstream_manifest.yaml").write_text("{}\n", encoding="utf-8")
        (run / "frozen" / "input_manifest.yaml").write_text(yaml.safe_dump({"pipeline": {"name": "bulk_rnaseq", "version": version}}), encoding="utf-8")
        (run / "run_state.json").write_text(json.dumps({"case_id": case, "run_id": run.name, "status": "SUCCESS"}), encoding="utf-8")
    report = SimpleNamespace(config=SimpleNamespace(input=SimpleNamespace(type=InputType.FASTQ)), project_dir=project)
    options = {option.reference.split("/")[0]: option for option in workspace.reuse_options(report)}
    assert options["NEW"].accepted
    assert not options["OLD"].accepted and "made by nf-rna 1.6.1" in options["OLD"].reason
