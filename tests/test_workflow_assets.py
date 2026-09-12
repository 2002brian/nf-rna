"""Regression coverage for first-party Nextflow workflow package assets."""

from pathlib import Path

from rnaseq.workflow_assets import REQUIRED_WORKFLOW_ASSETS, required_workflow_assets


ROOT = Path(__file__).parents[1]


def test_required_workflow_assets_are_discoverable_and_match_source_workflows():
    discovered = required_workflow_assets()

    assert tuple(discovered) == REQUIRED_WORKFLOW_ASSETS
    for name, asset in discovered.items():
        assert asset.is_file()
        assert asset.read_bytes() == (ROOT / "workflow" / name).read_bytes()
