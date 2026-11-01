"""v1.7.0 pre-release GUI hardening: path picker, project URL guard, metadata integrity,
numeric table display, run status refresh and disabled-control reasons.

Service tests need no GUI extra; page tests use Streamlit's AppTest and are
skipped without it.
"""

from __future__ import annotations

import json
import os
from pathlib import Path

import pytest
import yaml

from rnaseq.gui import results, services
from rnaseq.models import DesignType, InputType, Preset, Species
from rnaseq.project import create_project, import_count_matrix
from rnaseq.run_status import current_process_identity

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


def _write(path: Path, text: str) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")
    return path


def _paired_project(parent: Path, name: str = "paired", metadata: str | None = None) -> Path:
    """A raw-count L2 project with a pairing column and a technical batch (~ donor + batch + condition)."""

    project = create_project(
        project_name=name, destination=parent, species=Species.HUMAN, preset=Preset.L2, design_type=DesignType.PAIRED_TWO_GROUP,
        input_type=InputType.RAW_COUNTS, layout=None, scaffold=True, pair_id="donor", batch=True, gene_id_type="ENSEMBL",
    )
    counts = _write(parent / f"{name}-counts.csv", BASE_COUNTS)
    import_count_matrix(project, counts)
    if metadata is not None:
        _write(project / "metadata.csv", metadata)
    return project


def _page(name: str):
    pytest.importorskip("streamlit")
    from streamlit.testing.v1 import AppTest

    return AppTest.from_string(f"from rnaseq.gui import pages\npages.{name}()\n", default_timeout=60)


# A page whose st.data_editor also applies the cell edits stored for its key in
# session_state["_test_edits"], the way Streamlit applies a user's typed edits
# kept in the editor's own widget state (lost when the editor key changes).
EDITOR_SCRIPT = """
import streamlit
from rnaseq.gui import pages

class _Streamlit:
    def __getattr__(self, name):
        return getattr(streamlit, name)

    @staticmethod
    def data_editor(data, **kwargs):
        result = streamlit.data_editor(data, **kwargs)
        edits = streamlit.session_state.get("_test_edits", {}).get(kwargs.get("key"))
        if edits:
            result = result.copy()
            for sample, column, value in edits.get("cells", ()):
                result.loc[result["sample_id"] == sample, column] = value
            if edits.get("delete"):
                result = result[~result["sample_id"].isin(edits["delete"])]
        return result

pages.st = _Streamlit()
try:
    pages.{page}()
finally:
    pages.st = streamlit
"""


def _editor_page(page: str = "inputs_page"):
    pytest.importorskip("streamlit")
    from streamlit.testing.v1 import AppTest

    return AppTest.from_string(EDITOR_SCRIPT.replace("{page}", page), default_timeout=60)


def _type_cells(page, project: Path, cells=(), delete=()) -> None:
    key = f"metadata::{project}"
    editor = f"{key}::editor::{page.session_state[f'{key}::version']}"
    page.session_state["_test_edits"] = {editor: {"cells": list(cells), "delete": list(delete)}}


def _working(page, project: Path):
    return page.session_state[f"metadata::{project}::working"]


def _value(table, sample: str, column: str) -> str:
    return services.clean_cell(table.loc[table["sample_id"] == sample, column].iloc[0])


# --------------------------------------------------------------------- P0-1 path picker

def test_use_this_folder_selects_the_browsed_folder_without_an_exception(gui_home):
    folder = gui_home / "fastq" / "run1"
    for sample in ("C1", "C2"):
        for read in ("R1", "R2"):
            _write(folder / f"{sample}_{read}.fastq.gz", "x")
    page = _page("datasets_page").run()
    page.button(key="ds_goroot").click().run()
    assert page.session_state["ds_cwd"] == str(gui_home.resolve())
    page.selectbox(key="ds_sub").set_value("fastq").run()
    page.button(key="ds_open").click().run()
    page.selectbox(key="ds_sub").set_value("run1").run()
    page.button(key="ds_open").click().run()
    page.button(key="ds_usedir").click().run()
    assert not page.exception, page.exception
    assert page.text_input(key="ds_text").value == str(folder.resolve())
    assert [metric.value for metric in page.metric][:3] == ["4", "4 B", "4"]


def test_use_this_file_selects_the_browsed_file_without_an_exception(gui_home):
    project = _paired_project(gui_home / "projects")
    page = _page("inputs_page")
    page.session_state["project"] = str(project)
    page.run()
    page.selectbox(key="counts_sub").set_value("projects").run()
    page.button(key="counts_open").click().run()
    page.selectbox(key="counts_file").set_value("paired-counts.csv").run()
    page.button(key="counts_usefile").click().run()
    assert not page.exception, page.exception
    assert page.text_input(key="counts_text").value == str((gui_home / "projects" / "paired-counts.csv").resolve())
    assert any(button.key == "counts_import" for button in page.button)


def test_the_path_picker_never_browses_or_selects_outside_the_allowed_roots(gui_home, tmp_path):
    outside = tmp_path / "outside"
    _write(outside / "private" / "secret.csv", "a,b\n")
    (gui_home / "escape").symlink_to(outside, target_is_directory=True)
    (gui_home / "escape.csv").symlink_to(outside / "private" / "secret.csv")
    page = _page("datasets_page").run()
    assert "escape" not in page.selectbox(key="ds_sub").options
    assert any("not shown" in item.value for item in page.caption)
    # A forged or stale browse location falls back to the first allowed root.
    page.session_state["ds_cwd"] = str(gui_home / "escape" / "private")
    page.run()
    assert not page.exception and page.session_state["ds_cwd"] == str(gui_home.resolve())
    page.text_input(key="ds_text").set_value(str(gui_home / "escape")).run()
    assert any("outside the allowed roots" in item.value for item in page.error)

    project = _paired_project(gui_home / "projects")
    page = _page("inputs_page")
    page.session_state["project"] = str(project)
    page.session_state["counts_cwd"] = str(gui_home)
    page.run()
    assert "escape.csv" not in page.selectbox(key="counts_file").options


# --------------------------------------------------------------------- project URL guard

def test_a_url_project_outside_the_allowed_roots_is_not_opened(gui_home, tmp_path):
    elsewhere = tmp_path / "elsewhere"
    elsewhere.mkdir()
    outside = create_project(project_name="outsider", destination=elsewhere, species=Species.HUMAN, preset=Preset.L1,
                             design_type=DesignType.TWO_GROUP, input_type=InputType.RAW_COUNTS, layout=None, scaffold=True, gene_id_type="ENSEMBL")
    (gui_home / "projects" / "linked").symlink_to(outside, target_is_directory=True)
    for value in (str(outside), str(gui_home / "projects" / "linked"), str(gui_home / "projects" / ".." / ".." / "elsewhere" / "outsider")):
        for page_name in ("project_page", "inputs_page", "runner_page"):
            page = _page(page_name)
            page.query_params["project"] = value
            page.run()
            assert not page.exception, page.exception
            assert not page.title or page.title[0].value != "outsider"
            assert any("was not opened" in item.value for item in page.error), (value, page_name)
            assert "project" not in page.session_state and "project" not in page.query_params


def test_a_url_project_inside_the_allowed_roots_still_opens(gui_home):
    project = _paired_project(gui_home / "projects", "inside")
    page = _page("project_page")
    page.query_params["project"] = str(project)
    page.run()
    assert not page.exception and page.title[0].value == "inside"
    assert page.session_state["project"] == str(project.resolve())


# --------------------------------------------------------------------- P0-2 metadata integrity (services)

COLUMNS = ["sample_id", "donor", "batch", "condition"]


def _rows(*values: tuple[str, str, str, str]) -> list[dict[str, str]]:
    return [dict(zip(COLUMNS, row)) for row in values]


def test_merging_a_pasted_table_keeps_every_current_value_and_lists_replacements():
    current = _rows(("C1", "P1", "b1", "Control"), ("C2", "P2", "", "Control"), ("T1", "P1", "b1", ""), ("T2", "", "", "Treatment"))
    header, pasted = services.parse_table_text("sample_id\tbatch\tcondition\tRIN\nC1\tb2\t\t9.1\nC2\tb1\tControl\t8.7\nT1\tb1\tTreatment\t\nX9\tb3\tTreatment\t7\n")
    change = services.merge_metadata_table(COLUMNS, current, header, pasted)
    assert not change.errors
    assert change.columns == [*COLUMNS, "RIN"] and change.new_columns == ("RIN",)
    assert change.overwritten == (("C1", "batch", "b1", "b2"),)  # the only non-empty value replaced
    assert change.added == ("X9",)
    by_id = {row["sample_id"]: row for row in change.rows}
    assert by_id["C1"] == {"sample_id": "C1", "donor": "P1", "batch": "b2", "condition": "Control", "RIN": "9.1"}  # empty pasted condition kept
    assert by_id["T1"]["condition"] == "Treatment" and by_id["T1"]["donor"] == "P1"
    assert by_id["T2"] == {"sample_id": "T2", "donor": "", "batch": "", "condition": "Treatment", "RIN": ""}  # absent from the paste: unchanged
    assert [row["sample_id"] for row in change.rows] == ["C1", "C2", "T1", "T2", "X9"]  # order kept, new samples last


def test_merging_refuses_ambiguous_sample_identity():
    current = _rows(("C1", "P1", "", ""), ("C1", "P2", "", ""))
    header, pasted = services.parse_table_text("sample_id,condition\nC1,Control\n")
    assert "repeated in the table being edited" in services.merge_metadata_table(COLUMNS, current, header, pasted).errors[0]
    header, pasted = services.parse_table_text("sample_id,condition\nC2,Control\nC2,Treatment\n,Control\n")
    errors = services.merge_metadata_table(COLUMNS, _rows(("C2", "", "", "")), header, pasted).errors
    assert any("no sample_id" in error for error in errors) and any("repeated in the new table" in error for error in errors)
    assert services.merge_metadata_table(COLUMNS, [], ["condition"], [{"condition": "x"}]).errors


def test_bulk_assignment_changes_only_the_chosen_cells():
    current = _rows(("C1", "P1", "b1", "Control"), ("T1", "P1", "b1", "Treated"), ("T2", "P2", "", ""))
    change = services.assign_metadata_value(COLUMNS, current, ["T1", "T2"], "condition", "Treatment")
    assert change.overwritten == (("T1", "condition", "Treated", "Treatment"),) and change.filled == 1
    assert [row["condition"] for row in change.rows] == ["Control", "Treatment", "Treatment"]
    assert [row["donor"] for row in change.rows] == ["P1", "P1", "P2"]


def test_metadata_difference_reports_unsaved_cells_rows_and_columns():
    saved = _rows(("C1", "P1", "b1", "Control"), ("C2", "P2", "b1", "Control"))
    same = services.metadata_difference(COLUMNS, saved, COLUMNS, [dict(row) for row in saved])
    assert not same.unsaved
    edited = [dict(saved[0], condition="Treatment"), {"sample_id": "T9", "donor": "", "batch": "", "condition": "", "RIN": "7"}]
    difference = services.metadata_difference(COLUMNS, saved, [*COLUMNS, "RIN"], edited)
    assert difference.unsaved and difference.changed_cells == 1
    assert difference.added_rows == ("T9",) and difference.removed_rows == ("C2",) and difference.columns_added == ("RIN",)
    assert "1 changed value(s)" in difference.describe()


# --------------------------------------------------------------------- P0-2 metadata integrity (page)

def test_manual_metadata_edits_survive_bulk_assignment_and_are_saved(gui_home):
    project = _paired_project(gui_home / "projects")
    page = _editor_page()
    page.session_state["project"] = str(project)
    page.run()
    assert not page.exception, page.exception
    assert list(_working(page, project).columns) == ["sample_id", "donor", "batch", "condition"]
    assert any("No unsaved changes" in item.value for item in page.caption) is False  # prefilled sample IDs are not saved yet

    # The user types condition, pair and batch values into the table ...
    _type_cells(page, project, cells=[("C1", "condition", "Control"), ("C2", "condition", "Control"), ("C3", "condition", "Control"),
                                      ("C1", "donor", "P1"), ("T1", "donor", "P1"), ("C1", "batch", "b1"), ("T1", "batch", "b2")])
    page.run()
    # ... then assigns the treated samples in bulk.
    version = page.session_state[f"metadata::{project}::version"]
    bulk = f"metadata::{project}::bulk::{version}"
    page.multiselect(key=f"{bulk}::samples").set_value(["T1", "T2", "T3"]).run()
    page.selectbox(key=f"{bulk}::column").set_value("condition").run()
    page.text_input(key=f"{bulk}::value").set_value("Treatment").run()
    page.button(key=f"{bulk}::apply").click().run()
    assert not page.exception, page.exception
    table = _working(page, project)
    for sample in ("C1", "C2", "C3"):
        assert _value(table, sample, "condition") == "Control"  # typed before the bulk assignment: kept
    for sample in ("T1", "T2", "T3"):
        assert _value(table, sample, "condition") == "Treatment"
    assert (_value(table, "C1", "donor"), _value(table, "T1", "donor"), _value(table, "T1", "batch")) == ("P1", "P1", "b2")
    assert any("condition = 'Treatment' applied" in item.value for item in page.success)
    assert any("Unsaved changes" in item.value for item in page.warning)

    page.button(key=f"metadata::{project}::save").click().run()
    assert not page.exception, page.exception
    saved = (project / "metadata.csv").read_text().splitlines()
    assert saved[0] == "sample_id,donor,batch,condition"
    assert "C1,P1,b1,Control" in saved and "T1,P1,b2,Treatment" in saved and "T3,,,Treatment" in saved
    assert any("No unsaved changes" in item.value for item in page.caption)


def test_pasting_merges_into_unsaved_edits_and_replacing_needs_a_tick(gui_home):
    project = _paired_project(gui_home / "projects")
    page = _editor_page()
    page.session_state["project"] = str(project)
    page.run()
    _type_cells(page, project, cells=[("C1", "condition", "Control"), ("C1", "batch", "b1"), ("T1", "donor", "P1")], delete=["T3"])
    page.run()
    key = f"metadata::{project}"
    page.text_area(key=f"{key}::paste").set_value("sample_id\tbatch\tcondition\nC1\tb2\t\nC2\tb1\tControl\nX9\tb3\tTreatment\n").run()
    assert not page.exception, page.exception
    apply = next(button for button in page.button if button.key == f"{key}::paste_apply")
    assert apply.disabled  # C1 batch b1 -> b2 replaces a value the user typed
    assert any("existing value(s) would be replaced" in item.value for item in page.warning)
    assert any("X9" in item.value and "not samples of the input" in item.value for item in page.warning)
    page.checkbox(key=f"{key}::paste_confirm").check().run()
    page.button(key=f"{key}::paste_apply").click().run()
    assert not page.exception, page.exception
    table = _working(page, project)
    assert _value(table, "C1", "condition") == "Control" and _value(table, "C1", "batch") == "b2"
    assert _value(table, "T1", "donor") == "P1" and _value(table, "C2", "condition") == "Control"
    assert "T3" not in table["sample_id"].tolist()  # the row the user deleted stays deleted
    assert table["sample_id"].tolist() == ["C1", "C2", "C3", "T1", "T2", "X9"]
    assert any("Pasted table applied" in item.value and "replaced 1 value(s)" in item.value for item in page.success)


def test_reload_discards_unsaved_edits_only_after_a_tick(gui_home):
    project = _paired_project(gui_home / "projects")
    page = _editor_page()
    page.session_state["project"] = str(project)
    page.run()
    _type_cells(page, project, cells=[("C1", "condition", "Control")])
    page.run()
    key = f"metadata::{project}"
    assert next(button for button in page.button if button.key == f"{key}::reload").disabled
    page.checkbox(key=f"{key}::reload_confirm").check().run()
    page.button(key=f"{key}::reload").click().run()
    assert not page.exception
    assert _value(_working(page, project), "C1", "condition") == ""


def test_datasets_metadata_load_merges_with_unsaved_edits(gui_home):
    project = _paired_project(gui_home / "projects")
    folder = gui_home / "sheets"
    _write(folder / "metadata.csv", "sample_id,condition\nC1,Treatment\nC2,Control\nT1,Treatment\n")
    page = _editor_page()
    page.session_state["project"] = str(project)
    page.run()
    _type_cells(page, project, cells=[("C1", "condition", "Control"), ("T2", "donor", "P2")])
    page.run()
    working = _working(page, project)

    datasets = _page("datasets_page")
    datasets.session_state["project"] = str(project)
    datasets.session_state[f"metadata::{project}"] = working
    datasets.session_state[f"metadata::{project}::working"] = working
    datasets.session_state[f"metadata::{project}::version"] = page.session_state[f"metadata::{project}::version"]
    datasets.run()
    datasets.text_input(key="ds_text").set_value(str(folder)).run()
    assert not datasets.exception, datasets.exception
    load = next(button for button in datasets.button if button.key == "ds_meta_load")
    assert load.disabled and any("C1" in frame.value.to_string() and "Treatment" in frame.value.to_string() for frame in datasets.dataframe)
    confirm = next(box for box in datasets.checkbox if box.key.startswith("ds_meta_confirm::"))
    confirm.check().run()
    datasets.button(key="ds_meta_load").click().run()
    assert not datasets.exception, datasets.exception
    merged = datasets.session_state[f"metadata::{project}::working"]
    assert _value(merged, "C1", "condition") == "Treatment" and _value(merged, "T2", "donor") == "P2"
    assert any("not saved yet" in item.value for item in datasets.success)


# --------------------------------------------------------------------- numeric display

def test_numeric_columns_sort_numerically_with_na_last():
    table = results.Table(("gene_id", "symbol", "baseMean", "log2FoldChange", "pvalue", "padj", "flag"), tuple(
        dict(zip(("gene_id", "symbol", "baseMean", "log2FoldChange", "pvalue", "padj", "flag"), values)) for values in (
            ("ENSG1", "A", "9.5", "-10.2", "1e-18", "0.002", "TRUE"),
            ("ENSG2", "B", "100", "2.1", "3.03342369782533e-129", "NA", "FALSE"),
            ("ENSG3", "C", "20", "-Inf", "0.84", "1.8e-125", "TRUE"),
        )))
    assert results.numeric_columns(table) == {"baseMean": "float", "log2FoldChange": "float", "pvalue": "float", "padj": "float"}
    pytest.importorskip("streamlit")
    from rnaseq.gui import pages

    shown = pages.frame(table)
    assert shown.sort_values("pvalue")["symbol"].tolist() == ["B", "A", "C"]  # lexical order would be C, A, B
    assert shown.sort_values("log2FoldChange")["symbol"].tolist() == ["C", "A", "B"]
    assert shown.sort_values("padj")["symbol"].tolist() == ["C", "A", "B"]  # NA last
    assert shown["padj"].isna().tolist() == [False, True, False] and shown.loc[1, "pvalue"] == float("3.03342369782533e-129")
    assert shown["gene_id"].tolist() == ["ENSG1", "ENSG2", "ENSG3"] and shown["flag"].tolist() == ["TRUE", "FALSE", "TRUE"]


def test_enrichment_columns_are_numeric_and_identifiers_and_ratios_stay_text():
    columns = ("ID", "Description", "GeneRatio", "BgRatio", "pvalue", "p.adjust", "qvalue", "Count", "setSize", "NES", "enrichmentScore", "rank", "geneID", "entrez_id")
    rows = (("GO:0001", "a", "5/100", "50/18000", "1e-6", "0.001", "0.0009", "5", "120", "2.0", "0.54", "2634", "1/2/3", "7105"),
            ("GO:0002", "b", "12/100", "500/18000", "0.04", "NA", "NA", "12", "30", "-1.5", "-0.3", "120", "4/5", "8404"))
    table = results.Table(columns, tuple(dict(zip(columns, values)) for values in rows))
    found = results.numeric_columns(table)
    assert found == {"pvalue": "float", "p.adjust": "float", "qvalue": "float", "Count": "int", "setSize": "int", "NES": "float", "enrichmentScore": "float", "rank": "int"}
    pytest.importorskip("streamlit")
    from rnaseq.gui import pages

    shown = pages.frame(table)
    assert str(shown["Count"].dtype) == "Int64" and shown.sort_values("rank")["ID"].tolist() == ["GO:0002", "GO:0001"]
    assert shown["GeneRatio"].tolist() == ["5/100", "12/100"] and shown["entrez_id"].tolist() == ["7105", "8404"]


def test_numeric_display_leaves_the_tsv_files_unchanged(tmp_path):
    path = _write(tmp_path / "all_genes.tsv", "gene_id\tpadj\nENSG1\t1.82328987063955e-125\nENSG2\tNA\n")
    before = path.read_bytes()
    pytest.importorskip("streamlit")
    from rnaseq.gui import pages

    pages.frame(results.read_table(path))
    assert path.read_bytes() == before and results.read_table(path).rows[0]["padj"] == "1.82328987063955e-125"


# --------------------------------------------------------------------- run status refresh

def _running_run(project: Path, case: str = "LIVE") -> Path:
    run = project / "runs" / case / "20261101-120000+0800"
    _write(run / "run_state.json", json.dumps({"case_id": case, "run_id": run.name, "status": "RUNNING", "started_at": "2026-11-01T12:00:00+08:00"}))
    _write(run / "logs" / "rnaseq.process.json", json.dumps(current_process_identity()))
    _write(run / "logs" / "rnaseq.log", "2026-11-01 12:00:00 run started\n")
    return run


def test_an_active_run_refreshes_its_status_section_only(gui_home, monkeypatch):
    from rnaseq.gui import pages

    project = _paired_project(gui_home / "projects")
    assert not pages.live_run_in_project(project, None)
    run = _running_run(project)
    assert pages.live_run_in_project(project, None)
    page = _page("runner_page")
    page.session_state["project"] = str(project)
    page.run()
    assert not page.exception, page.exception
    assert any(f"every {pages.STATUS_REFRESH_SECONDS} s" in item.value for item in page.caption)
    assert any("How to stop this run" in item.label for item in page.expander)
    assert any(f"kill -TERM {os.getpid()}" in item.value for item in page.code)

    # The refreshed part re-reads run files only: validating the project there is an error.
    script = (
        "from pathlib import Path\nimport rnaseq.validators as validators\nfrom rnaseq.gui import pages\n"
        "def refuse(*args, **kwargs):\n    raise AssertionError('validated during a status refresh')\n"
        "validators.validate_project = refuse\n"
        f"pages._run_status_body(pages.gui_config(), Path({str(project)!r}), True)\n"
    )
    from streamlit.testing.v1 import AppTest

    original = __import__("rnaseq.validators", fromlist=["validate_project"]).validate_project
    try:
        body = AppTest.from_string(script, default_timeout=60).run()
    finally:
        import rnaseq.validators

        rnaseq.validators.validate_project = original
    assert not body.exception, body.exception

    # Terminal state: no automatic refresh any more.
    _write(run / "run_state.json", json.dumps({"case_id": "LIVE", "run_id": run.name, "status": "SUCCESS", "started_at": "2026-11-01T12:00:00+08:00"}))
    assert not pages.live_run_in_project(project, None)
    page.run()
    assert not page.exception
    assert any("press Refresh to update" in item.value for item in page.caption)
    assert not any("How to stop this run" in item.label for item in page.expander)


# --------------------------------------------------------------------- disabled controls and stale messages

def test_disabled_run_and_contrast_controls_say_why(gui_home):
    project = _paired_project(gui_home / "projects", metadata=(
        "sample_id,donor,batch,condition\nC1,P1,b1,Control\nC2,P2,b1,Control\nC3,P3,b2,Control\nT1,P1,b2,Treatment\nT2,P2,b2,Treatment\nT3,P3,b1,Treatment\n"))
    page = _page("contrasts_page")
    page.session_state["project"] = str(project)
    page.run()
    assert not page.exception, page.exception
    key = f"contrasts::{project}"
    assert any("*Save contrasts* is disabled: Add at least one contrast" in item.value for item in page.caption)
    page.selectbox(key=f"{key}::num").set_value("Control").run()
    page.selectbox(key=f"{key}::den").set_value("Control").run()
    assert next(button for button in page.button if button.key == f"{key}::add").disabled
    assert any("must be different levels" in item.value for item in page.caption)
    page.selectbox(key=f"{key}::num").set_value("Treatment").run()
    page.button(key=f"{key}::add").click().run()
    assert any("already in the list" in item.value for item in page.caption)
    page.button(key=f"{key}::save").click().run()
    assert not page.exception, page.exception
    assert any("Wrote 1 contrast(s)" in item.value for item in page.success)
    assert any("Treatment_vs_Control" in frame.value.to_string() for frame in page.dataframe)  # the saved table is redrawn
    assert any("Tick *Replace the saved contrasts*" in item.value for item in page.caption)

    runner = _page("runner_page")
    runner.session_state["project"] = str(project)
    runner.run()
    assert any("*Run L2* is disabled: enter a case ID" in item.value for item in runner.caption)
    runner.text_input(key=f"case::{project}").set_value("bad/../case").run()
    assert any("case ID is not usable" in item.value for item in runner.caption)
    assert next(button for button in runner.button if button.key == "run_start").disabled
    runner.text_input(key=f"case::{project}").set_value("GOOD-1").run()
    assert not next(button for button in runner.button if button.key == "run_start").disabled


def test_a_full_validation_report_of_an_older_project_state_is_not_shown(gui_home):
    project = _paired_project(gui_home / "projects")
    page = _page("runner_page")
    page.session_state["project"] = str(project)
    page.run()
    page.button(key="run_validate").click().run()
    assert any(item.label.startswith("Full validation report") for item in page.expander)
    metadata = project / "metadata.csv"
    metadata.write_text("sample_id,donor,batch,condition\nC1,P1,b1,Control\n", encoding="utf-8")
    os.utime(metadata, ns=(metadata.stat().st_atime_ns, metadata.stat().st_mtime_ns + 5_000_000_000))
    page.run()
    assert not any(item.label.startswith("Full validation report") for item in page.expander)
    assert any("full validation report was removed" in item.value for item in page.caption)


def test_other_reference_fields_are_visible_and_not_an_error_before_they_are_filled(gui_home):
    page = _page("new_project_page").run()
    page.text_input(key="np_name").set_value("ref_other").run()
    page.selectbox(key="np_reference").set_value("Other (advanced)").run()
    assert not page.exception and not page.error, [item.value for item in page.error]
    assert any(item.key == "np_ref_root" for item in page.text_input)
    assert any("Enter the managed reference folder" in item.value for item in page.caption)
    assert not any("--reference" in item.value for item in page.error)
    page.radio(key="np_ref_source").set_value("custom").run()
    assert any("Enter the genome FASTA and the annotation GTF" in item.value for item in page.caption)
    page.button(key="np_create").click().run()
    assert any("fill in the required fields" in item.value for item in page.error)
    assert not (gui_home / "projects" / "ref_other").exists()
