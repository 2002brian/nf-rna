"""New Project reference selection per species and backend (v1.7.0 Human HISAT2 reference fix).

References come only from the workstation registry; an absent or incomplete
reference is reported as not prepared, and a managed reference is re-checked
against the project's species and backend before the project is created.
Page tests use Streamlit's AppTest and are skipped without it.
"""

from __future__ import annotations

from pathlib import Path

import pytest
import yaml

from rnaseq.gui import services
from rnaseq.references import load_local_reference_root, register_local_reference
from rnaseq.validators import validate_project

from test_local_reference import _local_fastq_project, _mark_genome_only_hisat2, _promote_reference_for_production

HUMAN, MOUSE = "Homo sapiens", "Mus musculus"


def _reference(base: Path, name: str, species: str, *, hisat2: str | None = "validated", salmon: bool = True) -> Path:
    """A synthetic managed reference: production manifest, optional Salmon index, optional HISAT2 index."""

    project, reference = _local_fastq_project(base / name, species=species)
    if salmon:
        _promote_reference_for_production(project, reference)
    else:
        manifest_path = reference / "reference_manifest.yaml"
        manifest = yaml.safe_load(manifest_path.read_text(encoding="utf-8"))
        manifest.update({"schema_version": "1.1", "purpose": "production"})
        manifest_path.write_text(yaml.safe_dump(manifest, sort_keys=False), encoding="utf-8")
    if hisat2 is not None:
        _mark_genome_only_hisat2(reference, runtime_compatibility=hisat2)
    return reference


@pytest.fixture
def home(tmp_path, monkeypatch):
    data = tmp_path / "data"
    (data / "projects").mkdir(parents=True)
    settings = tmp_path / "gui.yaml"
    settings.write_text(yaml.safe_dump({"allowed_roots": [str(data)], "project_parent": str(data / "projects"), "condarc": None, "launcher": "process"}), encoding="utf-8")
    monkeypatch.setenv("NF_RNA_GUI_CONFIG", str(settings))
    monkeypatch.setenv("XDG_CACHE_HOME", str(tmp_path / "cache"))
    monkeypatch.setenv("XDG_CONFIG_HOME", str(tmp_path / "config"))  # the registry used by the GUI
    return tmp_path


def _register(*roots: Path) -> None:
    for root in roots:
        register_local_reference(root)


def _page():
    pytest.importorskip("streamlit")
    from streamlit.testing.v1 import AppTest

    page = AppTest.from_string("from rnaseq.gui import pages\npages.new_project_page()\n", default_timeout=60).run()
    page.radio(key="np_input").set_value("fastq").run()
    return page


def _set(page, species: str | None = None, method: str | None = None):
    if species is not None:
        page.selectbox(key="np_species").set_value(species).run()
    if method is not None:
        page.selectbox(key="np_method").set_value(method).run()
    assert not page.exception, page.exception
    return _ref(page).options


def _ref(page):
    """The reference selectbox of the current species and backend (its key includes both)."""

    return page.selectbox(key=f"np_reference::{page.selectbox(key='np_species').value}::{page.selectbox(key='np_method').value}")


def _create(page, name: str, *, strandedness: str = "reverse"):
    page.text_input(key="np_name").set_value(name).run()
    if page.selectbox(key="np_method").value == "hisat2_featurecounts":
        page.selectbox(key="np_strand_hisat2").set_value(strandedness).run()
    page.button(key="np_create").click().run()
    assert not page.exception, page.exception
    return page


# --------------------------------------------------------------------- services

def test_availability_follows_species_and_backend_and_explains_unusable_references(home):
    human = _reference(home / "refs", "human", HUMAN)
    mouse = _reference(home / "refs", "mouse", MOUSE)
    human_salmon_only = _reference(home / "refs", "human_salmon_only", HUMAN, hisat2=None)
    human_unvalidated = _reference(home / "refs", "human_unvalidated", HUMAN, hisat2="requires_smoke_validation")
    _register(human, mouse, human_salmon_only, human_unvalidated)

    hisat2 = services.reference_availability(HUMAN, "hisat2_featurecounts")
    assert [item.root for item in hisat2.compatible] == [human.resolve()]
    reasons = {item.root: reason for item, reason in hisat2.unusable}
    assert "no HISAT2 index" in reasons[human_salmon_only.resolve()] and "rnaseq reference prepare-hisat2" in reasons[human_salmon_only.resolve()]
    assert "no validated runtime compatibility" in reasons[human_unvalidated.resolve()]
    assert mouse.resolve() not in reasons  # another species is never listed for Human

    salmon = services.reference_availability(HUMAN, "salmon")
    assert {item.root for item in salmon.compatible} == {human.resolve(), human_salmon_only.resolve(), human_unvalidated.resolve()}
    assert [item.root for item in services.reference_availability(MOUSE, "hisat2_featurecounts").compatible] == [mouse.resolve()]


def test_a_salmon_index_is_never_taken_for_a_hisat2_index(home):
    salmon_only = _reference(home / "refs", "human_salmon_only", HUMAN, hisat2=None)
    loaded = load_local_reference_root(salmon_only)
    assert loaded.salmon_index is not None and loaded.hisat2_index is None
    problem = services.check_project_reference({"source": "local", "root": str(salmon_only), "manifest": "reference_manifest.yaml"}, HUMAN, "hisat2_featurecounts")
    assert problem and "no HISAT2 index" in problem
    assert services.check_project_reference({"source": "local", "root": str(salmon_only), "manifest": "reference_manifest.yaml"}, HUMAN, "salmon") is None


def test_the_creation_check_refuses_another_species_and_a_broken_root(home):
    mouse = _reference(home / "refs", "mouse", MOUSE)
    wrong = services.check_project_reference({"source": "local", "root": str(mouse), "manifest": "reference_manifest.yaml"}, HUMAN, "hisat2_featurecounts")
    assert wrong and "Mus musculus reference, not Homo sapiens" in wrong
    missing = services.check_project_reference({"source": "local", "root": str(home / "nowhere"), "manifest": "reference_manifest.yaml"}, HUMAN, "salmon")
    assert missing and "cannot be used" in missing
    assert services.check_project_reference({"source": "igenomes", "genome": "GRCh38"}, HUMAN, "salmon") is None


# --------------------------------------------------------------------- New Project page

def test_each_species_and_backend_offers_only_its_own_registered_reference(home):
    human = _reference(home / "refs", "human", HUMAN)
    mouse = _reference(home / "refs", "mouse", MOUSE)
    _register(mouse, human)
    page = _page()
    cases = [("human", "salmon", human, mouse), ("human", "hisat2_featurecounts", human, mouse),
             ("mouse", "salmon", mouse, human), ("mouse", "hisat2_featurecounts", mouse, human)]
    for species, method, own, other in cases:
        options = _set(page, species, method)
        assert any(str(own.resolve()) in option for option in options), (species, method, options)
        assert not any(str(other.resolve()) in option for option in options), (species, method, options)
        assert str(own.resolve()) in _ref(page).value  # preselected, never the other species
        assert ("iGenomes (Salmon only)" in options) == (method == "salmon")
        assert not any("not prepared" in item.value for item in page.warning)


def test_switching_species_and_backend_never_keeps_an_incompatible_selection(home):
    human = _reference(home / "refs", "human", HUMAN)
    mouse = _reference(home / "refs", "mouse", MOUSE)
    _register(human, mouse)
    page = _page()
    _set(page, "human", "hisat2_featurecounts")
    assert str(human.resolve()) in _ref(page).value
    _set(page, "mouse")  # Human -> Mouse
    assert str(mouse.resolve()) in _ref(page).value
    _set(page, "human")  # Mouse -> Human
    assert str(human.resolve()) in _ref(page).value
    _set(page, method="salmon")
    _ref(page).set_value("iGenomes (Salmon only)").run()
    _set(page, method="hisat2_featurecounts")  # Salmon -> HISAT2: iGenomes is not offered and not kept
    assert _ref(page).value != "iGenomes (Salmon only)"
    assert str(human.resolve()) in _ref(page).value


def test_other_chosen_for_lack_of_a_human_reference_is_not_carried_over_to_mouse(home):
    mouse = _reference(home / "refs", "mouse", MOUSE)
    _register(mouse)  # Human has no registered reference: "Other" is its only option
    page = _page()
    _set(page, "human", "hisat2_featurecounts")
    assert _ref(page).value == "Other (advanced)"
    _set(page, "mouse")  # Human -> Mouse: the registered Mouse reference, not the Human page's "Other"
    assert str(mouse.resolve()) in _ref(page).value
    assert not any("not prepared" in item.value for item in page.warning)
    _set(page, "human")  # Mouse -> Human: "not prepared", never the Mouse reference
    assert _ref(page).value == "Other (advanced)" and any("Homo sapiens HISAT2 reference not prepared" in item.value for item in page.warning)


def test_a_missing_human_hisat2_reference_is_reported_not_filled_in(home):
    mouse = _reference(home / "refs", "mouse", MOUSE)
    _register(mouse)  # the server state found on 2026-11-01: only Mouse registered
    page = _page()
    options = _set(page, "human", "hisat2_featurecounts")
    assert options == ["Other (advanced)"]
    assert any("Homo sapiens HISAT2 reference not prepared" in item.value for item in page.warning)
    assert any("rnaseq reference register" in item.value for item in page.caption)
    assert page.text_input(key="np_ref_root").value == ""  # no path is invented
    _create(page, "human_missing")
    assert any("Choose a reference" in item.value for item in page.error)
    assert not (home / "data" / "projects" / "human_missing").exists()
    # The Mouse root typed under Other for a Human project is refused, not substituted.
    page.text_input(key="np_ref_root").set_value(str(mouse)).run()
    assert any("Mus musculus reference, not Homo sapiens" in item.value for item in page.error)
    _create(page, "human_with_mouse_ref")
    assert not (home / "data" / "projects" / "human_with_mouse_ref").exists()


def test_an_incomplete_human_reference_names_what_is_missing(home):
    salmon_only = _reference(home / "refs", "human_salmon_only", HUMAN, hisat2=None)
    unvalidated = _reference(home / "refs", "human_unvalidated", HUMAN, hisat2="requires_smoke_validation")
    _register(salmon_only, unvalidated)
    page = _page()
    _set(page, "human", "hisat2_featurecounts")
    warning = next(item.value for item in page.warning if "not prepared" in item.value)
    assert str(salmon_only.resolve()) in warning and "rnaseq reference prepare-hisat2" in warning
    assert str(unvalidated.resolve()) in warning and "no validated runtime compatibility" in warning
    assert _ref(page).options == ["Other (advanced)"]
    page.text_input(key="np_ref_root").set_value(str(unvalidated)).run()
    assert any("no validated runtime compatibility" in item.value for item in page.error)


def test_registered_and_valid_custom_references_create_valid_human_hisat2_projects(home):
    human = _reference(home / "refs", "human", HUMAN)
    _register(human)
    page = _page()
    _set(page, "human", "hisat2_featurecounts")
    _create(page, "human_registered")
    assert not page.error, [item.value for item in page.error]
    config = yaml.safe_load((home / "data" / "projects" / "human_registered" / "project.yaml").read_text())
    assert config["reference"] == {"source": "local", "root": str(human.resolve()), "manifest": "reference_manifest.yaml"}
    assert config["organism"]["species"] == HUMAN and config["upstream"]["quantification"]["method"] == "hisat2_featurecounts"
    report = validate_project(home / "data" / "projects" / "human_registered")
    assert report.local_reference is not None and report.local_reference.hisat2_runtime_ready
    assert not [issue for issue in report.errors if "reference" in issue.code or "hisat2" in issue.code], report.errors

    # The same reference entered by hand under Other (advanced).
    page = _page()
    _set(page, "human", "hisat2_featurecounts")
    _ref(page).set_value("Other (advanced)").run()
    page.text_input(key="np_ref_root").set_value(str(human)).run()
    assert not page.error, [item.value for item in page.error]
    _create(page, "human_custom")
    assert not page.error, [item.value for item in page.error]
    assert yaml.safe_load((home / "data" / "projects" / "human_custom" / "project.yaml").read_text())["reference"]["root"] == str(human.resolve())


def test_the_validator_still_refuses_a_mismatched_reference_at_run_time(home):
    human = _reference(home / "refs", "human", HUMAN)
    mouse = _reference(home / "refs", "mouse", MOUSE)
    _register(human, mouse)
    page = _page()
    _set(page, "human", "hisat2_featurecounts")
    _create(page, "edited_later")
    project = home / "data" / "projects" / "edited_later"
    config = yaml.safe_load((project / "project.yaml").read_text())
    config["reference"]["root"] = str(mouse.resolve())  # project.yaml edited by hand after creation
    (project / "project.yaml").write_text(yaml.safe_dump(config, sort_keys=False))
    issues = [issue for issue in validate_project(project).errors if issue.code == "invalid_local_reference"]
    assert issues and "does not match project organism 'Homo sapiens'" in issues[0].message


def test_quantification_backends_use_the_same_labels_in_new_project_and_overview(home):
    """Display labels only: project.yaml keeps the backend identifiers, and Salmon is the default."""

    from streamlit.testing.v1 import AppTest

    from rnaseq.gui import workspace

    labels = {"salmon": "Salmon (pseudoalignment)", "hisat2_featurecounts": "HISAT2 + featureCounts (genome alignment)"}
    human = _reference(home / "refs", "human", HUMAN)
    _register(human)
    page = _page()
    method = page.selectbox(key="np_method")
    assert method.value == "salmon" and method.options == list(labels.values())
    for value, label in labels.items():
        page = _page()
        _set(page, "human", value)
        _create(page, f"label_{value}")
        assert not page.error, [item.value for item in page.error]
        project = home / "data" / "projects" / f"label_{value}"
        assert yaml.safe_load((project / "project.yaml").read_text())["upstream"]["quantification"]["method"] == value
        assert workspace.project_summary(project).quantification == label
        overview = AppTest.from_string("from rnaseq.gui import pages\npages.project_page()\n", default_timeout=60)
        overview.session_state["project"] = str(project)
        overview.run()
        assert not overview.exception, overview.exception
        shown = [frame.value.to_csv() for frame in overview.dataframe]
        assert any(f"quantification,{label}" in table for table in shown), shown
