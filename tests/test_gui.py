"""GUI service layer and pages (v1.6.0 Phase 6a).

The service tests need no GUI extra; the page tests use Streamlit's AppTest
and are skipped when Streamlit is not installed.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest
import yaml
from typer.testing import CliRunner

from rnaseq.cli import app
from rnaseq.gui import GUI_ADDRESS, streamlit_command
from rnaseq.gui import jobs, services
from rnaseq.gui.config import PathNotAllowed, ensure_data_layout, load_config, root_status
from rnaseq.models import DesignType, InputType, Preset, SequencingLayout, Species
from rnaseq.project import create_project, fastq_samplesheet_text, import_count_matrix, import_fastq_samplesheet
from rnaseq.run_status import current_process_identity
from rnaseq.validators import validate_project

from conftest import BASE_COUNTS

runner = CliRunner()


@pytest.fixture
def gui_home(tmp_path, monkeypatch):
    """GUI settings and state in temporary directories, with one allowed data root."""

    data = tmp_path / "data"
    (data / "projects").mkdir(parents=True)
    settings = tmp_path / "gui.yaml"
    settings.write_text(yaml.safe_dump({"allowed_roots": [str(data)], "project_parent": str(data / "projects"), "condarc": None, "launcher": "process"}), encoding="utf-8")
    monkeypatch.setenv("NF_RNA_GUI_CONFIG", str(settings))
    monkeypatch.setenv("XDG_CACHE_HOME", str(tmp_path / "cache"))
    return data


def _raw_project(parent: Path, name: str = "demo", preset: Preset = Preset.L2) -> Path:
    project = create_project(
        project_name=name, destination=parent, species=Species.HUMAN, preset=preset, design_type=DesignType.TWO_GROUP,
        input_type=InputType.RAW_COUNTS, layout=None, scaffold=True, gene_id_type="ENSEMBL",
    )
    counts = parent / f"{name}-counts.csv"
    counts.write_text(BASE_COUNTS, encoding="utf-8")
    import_count_matrix(project, counts)
    return project


def test_paths_outside_the_allowed_roots_are_refused(gui_home, tmp_path):
    config = load_config()
    inside = gui_home / "projects"
    assert config.resolve_allowed(str(inside)) == inside.resolve()
    (gui_home / "escape").symlink_to(tmp_path)
    with pytest.raises(PathNotAllowed):
        config.resolve_allowed(gui_home / "escape" / "gui.yaml")
    with pytest.raises(PathNotAllowed):
        config.resolve_allowed(tmp_path)


def test_an_unmounted_media_root_is_unavailable_and_nothing_is_created(tmp_path):
    missing = Path("/media/nf-rna-test-not-mounted/brian/nf-rna")
    status = root_status(missing)
    assert not status.available and "not mounted" in status.message
    config = load_config(tmp_path / "gui.yaml")
    config = type(config)((missing,), missing / "projects", None)
    assert "not created" in ensure_data_layout(config)[0]
    assert not missing.exists()


def test_metadata_can_be_saved_incomplete_but_never_structurally_invalid(gui_home):
    project = _raw_project(gui_home / "projects")
    columns = ["sample_id", "condition"]
    incomplete = [{"sample_id": sample, "condition": ""} for sample in ("C1", "C2", "C3", "T1", "T2", "T3")]
    result = services.save_metadata(project, columns, incomplete)
    assert result.saved and {issue.code for issue in result.incomplete} == {"missing_design_value"}
    saved = (project / "metadata.csv").read_bytes()

    duplicate = [{"sample_id": "C1", "condition": "Control"}, {"sample_id": "C1", "condition": "Control"}]
    result = services.save_metadata(project, columns, duplicate)
    assert not result.saved and result.errors[0].code == "duplicate_sample_id"
    assert services.issue_column(result.errors[0], columns) == "sample_id"
    assert (project / "metadata.csv").read_bytes() == saved
    assert not list(project.glob(".metadata.csv.gui-*"))


def test_pasted_excel_rows_are_tab_separated(gui_home):
    header, rows = services.parse_table_text("sample_id\tcondition\r\nC1\tControl\r\nT1\tTreatment\r\n")
    assert header == ["sample_id", "condition"]
    assert rows == [{"sample_id": "C1", "condition": "Control"}, {"sample_id": "T1", "condition": "Treatment"}]


def test_gui_contrasts_write_what_rnaseq_design_writes(gui_home, tmp_path):
    import shutil

    from rnaseq.design_setup import apply_design, resolve_design_enrichment

    project = _raw_project(gui_home / "projects", preset=Preset.L1)
    (project / "metadata.csv").write_text("sample_id,condition\nC1,Control\nC2,Control\nC3,Control\nT1,Treatment\nT2,Treatment\nT3,Treatment\n", encoding="utf-8")
    cli_copy = tmp_path / "cli-copy"
    shutil.copytree(project, cli_copy)
    assert runner.invoke(app, ["design", str(cli_copy), "--yes"]).exit_code == 0

    report = validate_project(project)
    view = services.design_view(report)
    assert view.levels == ("Control", "Treatment")
    definition = services.contrast_definition(view.factor, "Treatment", "Control")
    selection, annotation = resolve_design_enrichment(report.config, None, None)
    after = apply_design(report, services.design_proposal(report, view.levels, [definition]), enrichment=selection, annotation=annotation)
    assert after.is_valid
    for name in ("project.yaml", "contrasts.csv", "metadata.csv"):
        assert (project / name).read_bytes() == (cli_copy / name).read_bytes(), name


def test_fastq_files_are_recognized_and_imported_by_samplesheet(gui_home):
    folder = gui_home / "fastq" / "batch1"
    folder.mkdir(parents=True)
    for sample in ("Ctrl1_S1", "Trt1_S2"):
        for read in ("R1", "R2"):
            (folder / f"{sample}_L001_{read}_001.fastq.gz").write_bytes(b"")
    summary, errors = services.discover_fastq(folder, SequencingLayout.PAIRED_END)
    assert not errors and summary.sample_ids == ("Ctrl1_S1", "Trt1_S2")
    sheet = gui_home / "samplesheet.csv"
    sheet.write_text(fastq_samplesheet_text(summary.records, "auto"), encoding="utf-8")
    assert sheet.read_text().splitlines()[1] == f"Ctrl1_S1,{folder}/Ctrl1_S1_L001_R1_001.fastq.gz,{folder}/Ctrl1_S1_L001_R2_001.fastq.gz,auto"

    project = create_project(project_name="fq", destination=gui_home / "projects", species=Species.HUMAN, preset=Preset.L1, design_type=DesignType.TWO_GROUP, scaffold=True, gene_id_type="ENSEMBL")
    import_fastq_samplesheet(project, sheet)
    assert sorted(path.name for path in (project / "input" / "fastq").iterdir()) == sorted(path.name for path in folder.iterdir())
    assert validate_project(project).fastq.sample_ids == ("Ctrl1_S1", "Trt1_S2")
    with pytest.raises(Exception, match="already contains"):
        import_fastq_samplesheet(project, sheet)


def _fake_start(identity):
    def start(command, env, log, cwd):
        log.write_text("started\n", encoding="utf-8")
        return identity
    return start


def test_only_one_gui_run_can_be_active(gui_home, monkeypatch):
    project = _raw_project(gui_home / "projects")
    config = load_config()
    alive = current_process_identity()  # this test process: certainly alive
    monkeypatch.setattr(jobs, "_start_process", _fake_start(alive))
    first = jobs.launch_run(project, "CASE-1", config)
    assert jobs.read_lock()["job_id"] == first.job_id and jobs.active_launch().job_id == first.job_id
    with pytest.raises(jobs.ActiveRunError):
        jobs.launch_run(project, "CASE-2", config)
    assert len(jobs.list_launches(project)) == 1


def test_a_lock_whose_process_is_gone_is_cleared(gui_home, monkeypatch):
    project = _raw_project(gui_home / "projects")
    config = load_config()
    gone = {**current_process_identity(), "start_ticks": -1}  # same PID, different process start: not ours
    monkeypatch.setattr(jobs, "_start_process", _fake_start(gone))
    jobs.launch_run(project, "CASE-1", config)
    assert jobs.active_launch() is None and jobs.read_lock() is None
    jobs.launch_run(project, "CASE-2", config)  # a new run can start


def test_a_running_run_without_its_process_is_shown_incomplete_and_not_rewritten(tmp_path):
    run_dir = tmp_path / "proj" / "runs" / "CASE" / "20261030-120000+0800"
    (run_dir / "logs").mkdir(parents=True)
    state = {"case_id": "CASE", "run_id": run_dir.name, "status": "RUNNING", "phase": "downstream", "started_at": "2026-10-30T12:00:00+08:00"}
    (run_dir / "run_state.json").write_text(json.dumps(state), encoding="utf-8")
    (run_dir / "logs" / "rnaseq.process.json").write_text(json.dumps({**current_process_identity(), "start_ticks": -1}), encoding="utf-8")
    before = (run_dir / "run_state.json").read_bytes()
    shown, note = jobs.run_display_state(run_dir, state)
    assert shown == "INCOMPLETE" and "no longer exists" in note
    assert (run_dir / "run_state.json").read_bytes() == before


def test_rnaseq_command_is_an_argument_list_and_gui_listens_on_loopback(tmp_path):
    command = jobs.rnaseq_command(tmp_path / "my project", "CASE-1", "L1-001/20261030-120000+0800")
    assert command[-7:] == ["run", str(tmp_path / "my project"), "--case-id", "CASE-1", "--yes", "--reuse-upstream", "L1-001/20261030-120000+0800"]
    server = streamlit_command(8600)
    assert server[server.index("--server.address") + 1] == GUI_ADDRESS == "127.0.0.1"
    assert server[server.index("--server.port") + 1] == "8600"


def test_gui_command_without_the_extra_names_the_install_command(monkeypatch):
    monkeypatch.setattr("rnaseq.gui.streamlit_available", lambda: False)
    result = runner.invoke(app, ["gui"])
    assert result.exit_code == 1 and "pip install" in result.output and "[gui]" in result.output


# --------------------------------------------------------------------- pages (AppTest)

def _page(name: str):
    pytest.importorskip("streamlit")
    from streamlit.testing.v1 import AppTest

    script = f"from rnaseq.gui import pages\npages.{name}()\n"
    return AppTest.from_string(script, default_timeout=60)


def test_dashboard_lists_projects_and_runs(gui_home, monkeypatch):
    project = _raw_project(gui_home / "projects")
    monkeypatch.setattr("rnaseq.execution.doctor_checks", lambda project_dir=None: ())
    page = _page("dashboard_page").run()
    assert not page.exception
    tables = [frame.value for frame in page.dataframe]
    assert any("demo" in table.to_string() for table in tables)
    assert any(str(project) in table.to_string() for table in tables)


def test_new_project_page_creates_the_wizard_scaffold(gui_home):
    page = _page("new_project_page").run()
    page.text_input(key="np_name").set_value("gui_made").run()
    page.radio(key="np_input").set_value("raw_counts").run()
    page.radio(key="np_preset_raw_counts").set_value("L2").run()
    page.selectbox(key="np_design").set_value("paired_two_group").run()
    page.text_input(key="np_pair").set_value("donor").run()
    page.button(key="np_create").click().run()
    assert not page.exception and not page.error, [item.value for item in page.error]
    created = gui_home / "projects" / "gui_made"
    config = yaml.safe_load((created / "project.yaml").read_text())
    assert config["design"] == {"type": "paired_two_group", "formula": "~ donor + condition", "variables": {"donor": "categorical", "condition": "categorical"}, "pair_id": "donor"}
    assert config["analysis"]["enrichment"] == ["go", "kegg", "gsea"] and config["annotation"]["input_id_type"] == "ENSEMBL"
    assert (created / "metadata.csv").read_text() == "sample_id,donor,condition\n"


def test_run_page_disables_run_while_a_gui_run_is_active(gui_home, monkeypatch):
    project = _raw_project(gui_home / "projects")
    monkeypatch.setattr(jobs, "_start_process", _fake_start(current_process_identity()))
    jobs.launch_run(project, "CASE-1", load_config())
    page = _page("runner_page")
    page.session_state["project"] = str(project)
    page.run()
    assert not page.exception
    run = next(button for button in page.button if button.key == "run_start")
    assert run.disabled
    assert any("A GUI run is active" in item.value for item in page.warning)
