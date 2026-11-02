"""GUI default settings are portable: they derive from the current user's home, never from a workstation."""

from __future__ import annotations

from pathlib import Path

import pytest
import yaml

from rnaseq.gui.config import PathNotAllowed, config_path, ensure_data_layout, load_config

ROOT = Path(__file__).resolve().parents[1]


@pytest.fixture
def new_user(tmp_path, monkeypatch):
    """A user with an empty home, no gui.yaml and no workspace."""

    home = tmp_path / "home" / "alice"
    home.mkdir(parents=True)
    monkeypatch.setenv("HOME", str(home))
    monkeypatch.setenv("NF_RNA_GUI_CONFIG", str(home / ".config" / "nf-rna" / "gui.yaml"))
    monkeypatch.setenv("XDG_CACHE_HOME", str(home / ".cache"))
    return home


def test_defaults_follow_the_current_home_and_are_written_once(new_user):
    config = load_config()
    workspace = new_user / "nf-rna-work"
    assert config.allowed_roots == (workspace,)
    assert config.project_parent == workspace / "projects"
    assert config.condarc is None  # no condarc file in the new home
    written = yaml.safe_load(config.source.read_text(encoding="utf-8"))
    assert written["allowed_roots"] == [str(workspace)] and written["project_parent"] == str(workspace / "projects")
    assert "/media/" not in config.source.read_text(encoding="utf-8") and "~" not in config.source.read_text(encoding="utf-8")
    assert config.notes and "Wrote default GUI settings" in config.notes[0]


def test_two_users_get_their_own_workspaces(tmp_path, monkeypatch):
    for name in ("alice", "bob"):
        home = tmp_path / name
        home.mkdir()
        monkeypatch.setenv("HOME", str(home))
        config = load_config(home / "gui.yaml")
        assert config.allowed_roots == (home / "nf-rna-work",) and config.project_parent == home / "nf-rna-work" / "projects"


def test_an_existing_default_condarc_is_used(new_user):
    condarc = new_user / "nf-rna-work" / "conda" / "condarc.yml"
    condarc.parent.mkdir(parents=True)
    condarc.write_text("channels: [conda-forge]\n", encoding="utf-8")
    assert load_config().condarc == condarc


def test_a_missing_default_workspace_is_created_on_first_use(new_user):
    config = load_config()
    workspace = new_user / "nf-rna-work"
    assert not workspace.exists() and config.available_roots() == []
    messages = ensure_data_layout(config)
    assert (workspace / "projects").is_dir() and (workspace / "fastq").is_dir()
    assert any("default workspace" in message for message in messages)
    assert config.available_roots() == [workspace]
    assert config.resolve_allowed(workspace / "projects") == (workspace / "projects").resolve()
    assert ensure_data_layout(config) == []  # idempotent


def test_a_dangling_default_workspace_link_is_not_replaced(new_user, tmp_path):
    workspace = new_user / "nf-rna-work"
    workspace.symlink_to(tmp_path / "gone")
    messages = ensure_data_layout(load_config())
    assert "not created" in messages[0]
    assert workspace.is_symlink() and not (tmp_path / "gone").exists()


def test_a_custom_gui_yaml_takes_precedence_and_is_never_rewritten(new_user, tmp_path):
    data = tmp_path / "server-disk" / "nf-rna"
    testdata = new_user / "nf-rna-work" / "testdata"
    (data / "projects").mkdir(parents=True)
    testdata.mkdir(parents=True)
    condarc = tmp_path / "condarc.yml"
    settings = config_path()
    settings.parent.mkdir(parents=True)
    settings.write_text(yaml.safe_dump({"allowed_roots": [str(data), str(testdata)], "project_parent": str(data / "projects"), "condarc": str(condarc), "launcher": "process"}), encoding="utf-8")
    before = settings.read_bytes()
    config = load_config()
    assert config.allowed_roots == (data, testdata) and config.project_parent == data / "projects" and config.condarc == condarc
    assert config.launcher == "process" and config.notes == ()
    assert ensure_data_layout(config) == [f"Created {data / 'fastq'}."]
    assert settings.read_bytes() == before
    assert not (new_user / "nf-rna-work" / "projects").exists()  # the default workspace is not created for a custom data root


def test_a_tilde_in_a_custom_gui_yaml_expands_to_the_current_home(new_user):
    settings = new_user / "custom.yaml"
    settings.write_text(yaml.safe_dump({"allowed_roots": ["~/data"], "project_parent": "~/data/projects", "condarc": "~/condarc.yml"}), encoding="utf-8")
    config = load_config(settings)
    assert config.allowed_roots == (new_user / "data",) and config.project_parent == new_user / "data" / "projects" and config.condarc == new_user / "condarc.yml"


def test_a_partial_gui_yaml_falls_back_to_the_portable_defaults(new_user):
    settings = new_user / "partial.yaml"
    settings.write_text(yaml.safe_dump({"log_tail_lines": 50}), encoding="utf-8")
    config = load_config(settings)
    assert config.allowed_roots == (new_user / "nf-rna-work",) and config.project_parent == new_user / "nf-rna-work" / "projects"
    assert config.log_tail_lines == 50


def test_a_project_parent_outside_the_allowed_roots_creates_nothing(new_user, tmp_path):
    root = tmp_path / "root"
    root.mkdir()
    outside = tmp_path / "elsewhere"
    outside.mkdir()
    settings = new_user / "outside.yaml"
    settings.write_text(yaml.safe_dump({"allowed_roots": [str(root)], "project_parent": str(outside / "projects")}), encoding="utf-8")
    config = load_config(settings)
    messages = ensure_data_layout(config)
    assert "outside the allowed roots" in messages[0] and "not created" in messages[0]
    assert not (outside / "projects").exists() and not (outside / "fastq").exists()
    with pytest.raises(PathNotAllowed):
        config.resolve_allowed(config.project_parent)


def test_the_default_workspace_keeps_the_symlink_and_traversal_guard(new_user, tmp_path):
    config = load_config()
    ensure_data_layout(config)
    workspace = new_user / "nf-rna-work"
    secret = tmp_path / "secret"
    secret.mkdir()
    (workspace / "escape").symlink_to(secret)
    for value in (workspace / "escape", workspace / "escape" / "file", f"{workspace}/projects/../../.config", new_user, "/etc", "~/.ssh"):
        with pytest.raises(PathNotAllowed):
            config.resolve_allowed(value)
    assert config.resolve_allowed("~/nf-rna-work/projects") == (workspace / "projects").resolve()


def test_shipped_gui_code_and_docs_name_no_workstation_path():
    needles = ("/media/" + "duobook02", "/home/" + "brian")
    files = [*(ROOT / "src" / "rnaseq").rglob("*.py"), ROOT / "docs" / "gui.md"]
    offenders = [f"{path.relative_to(ROOT)}: {needle}" for path in files for needle in needles if needle in path.read_text(encoding="utf-8")]
    assert offenders == []
