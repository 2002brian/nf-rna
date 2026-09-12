"""Locate bundled first-party Nextflow assets without assuming a source tree."""

from __future__ import annotations

from importlib import resources
from pathlib import Path


WORKFLOW_RESOURCE_PACKAGE = "rnaseq.workflows"
REQUIRED_WORKFLOW_ASSETS = (
    "main.nf",
    "hisat2_featurecounts.nf",
    "nextflow.config",
)


def workflow_asset_path(name: str) -> Path:
    """Return one installed workflow asset as a real filesystem path.

    Wheels are installed as normal filesystem trees, which Nextflow requires
    because it receives paths rather than Python streams.  ``importlib.resources``
    keeps lookup tied to the installed package instead of a repository-relative
    assumption.
    """

    if name not in REQUIRED_WORKFLOW_ASSETS:
        raise ValueError(f"Unsupported workflow asset: {name!r}.")
    asset = resources.files(WORKFLOW_RESOURCE_PACKAGE).joinpath(name)
    try:
        path = Path(asset)
    except TypeError as exc:  # pragma: no cover - wheels install as filesystems.
        raise RuntimeError("Installed workflow assets are not filesystem-backed.") from exc
    if path.is_file():
        return path

    # Hatch maps the repository-level workflow directory into the wheel.  An
    # editable/source checkout intentionally keeps that directory in place, so
    # retain this narrowly-scoped development fallback without using it for an
    # installed distribution.
    development_asset = Path(__file__).resolve().parents[2] / "workflow" / name
    if development_asset.is_file():
        return development_asset
    raise RuntimeError(f"Installed workflow asset is missing: {name}.")


def required_workflow_assets() -> dict[str, Path]:
    """Return every first-party executable workflow/config asset."""

    return {name: workflow_asset_path(name) for name in REQUIRED_WORKFLOW_ASSETS}
