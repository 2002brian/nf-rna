"""Release automation and metadata guards for the Conda runtime (v1.3.0+)."""

import re
from pathlib import Path

import yaml

from rnaseq import __version__


ROOT = Path(__file__).parents[1]
WORKFLOWS = sorted((ROOT / ".github" / "workflows").glob("*.y*ml"))


def _triggers(path: Path) -> dict:
    document = yaml.safe_load(path.read_text(encoding="utf-8"))
    # PyYAML reads the bare ``on:`` key as boolean True.
    triggers = document.get("on", document.get(True))
    if isinstance(triggers, str):
        return {triggers: None}
    if isinstance(triggers, list):
        return dict.fromkeys(triggers)
    return triggers or {}


def test_no_workflow_publishes_a_docker_image_or_reacts_to_a_release():
    # GHCR publication ended with the historical v1.2.1 Docker release.
    assert WORKFLOWS
    assert not (ROOT / ".github" / "workflows" / "publish-ghcr.yml").exists()
    for path in WORKFLOWS:
        triggers = _triggers(path)
        assert "release" not in triggers, path.name
        assert "tags" not in (triggers.get("push") or {}), path.name
        text = path.read_text(encoding="utf-8")
        for marker in ("docker/build-push-action", "docker/login-action", "ghcr.io", "packages: write"):
            assert marker not in text, f"{path.name}: {marker}"


def test_release_metadata_matches_the_package_version():
    citation = yaml.safe_load((ROOT / "CITATION.cff").read_text(encoding="utf-8"))
    assert citation["version"] == __version__
    changelog = (ROOT / "CHANGELOG.md").read_text(encoding="utf-8")
    first = re.search(r"^## (.+)$", changelog, flags=re.MULTILINE).group(1)
    assert first == f"{__version__} — {citation['date-released']}"
