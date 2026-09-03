"""Project creation and strict configuration loading."""

from __future__ import annotations

import shutil
import tempfile
from dataclasses import dataclass
from importlib import resources
from pathlib import Path
from typing import Any

import yaml
from pydantic import ValidationError

from rnaseq.errors import ProjectConfigError, ProjectCreationError
from rnaseq.models import (
    DesignType,
    InputType,
    NFCORE_RNASEQ_VERSION,
    Preset,
    ProjectConfig,
    ProjectInfo,
    SequencingLayout,
    Species,
    SUPPORTED_SCHEMA_VERSION,
)


@dataclass(frozen=True)
class LoadedProject:
    """A validated configuration and its resolved project paths."""

    root: Path
    config_path: Path
    config: ProjectConfig
    input_path: Path
    metadata_path: Path
    contrasts_path: Path

    @property
    def counts_path(self) -> Path:
        """Backward-compatible alias used only by the raw-count validator."""

        return self.input_path


def _format_pydantic_errors(error: ValidationError) -> str:
    messages: list[str] = []
    for item in error.errors(include_url=False):
        location = ".".join(str(part) for part in item["loc"])
        message = str(item["msg"])
        if message.startswith("Value error, "):
            message = message.removeprefix("Value error, ")
        messages.append(f"{location}: {message}")
    return "\n".join(messages)


def _resolve_project_path(root: Path, configured: str, label: str) -> Path:
    if not configured:
        raise ProjectConfigError(f"{label} must not be blank.")
    relative = Path(configured)
    if relative.is_absolute():
        raise ProjectConfigError(f"{label} must be relative to the project directory.")
    root_resolved = root.resolve()
    resolved = (root / relative).resolve()
    try:
        resolved.relative_to(root_resolved)
    except ValueError as exc:
        raise ProjectConfigError(
            f"{label} must remain inside the project directory."
        ) from exc
    return resolved


def load_project(project_dir: Path | str) -> LoadedProject:
    """Load a project configuration without reading or modifying its data files."""

    root = Path(project_dir)
    if not root.exists():
        raise ProjectConfigError(f"Project directory does not exist: {root}")
    if not root.is_dir():
        raise ProjectConfigError(f"Project path is not a directory: {root}")

    config_path = root / "project.yaml"
    if not config_path.exists():
        raise ProjectConfigError(f"Missing project configuration: {config_path}")
    if not config_path.is_file():
        raise ProjectConfigError(f"Project configuration is not a file: {config_path}")

    try:
        raw: Any = yaml.safe_load(config_path.read_text(encoding="utf-8"))
    except yaml.YAMLError as exc:
        raise ProjectConfigError(f"Invalid YAML in project.yaml: {exc}") from exc

    if not isinstance(raw, dict):
        raise ProjectConfigError("project.yaml must contain a YAML mapping at its root.")

    try:
        config = ProjectConfig.model_validate(raw)
    except ValidationError as exc:
        raise ProjectConfigError(_format_pydantic_errors(exc)) from exc

    return LoadedProject(
        root=root.resolve(),
        config_path=config_path.resolve(),
        config=config,
        input_path=_resolve_project_path(root, config.input.path, "input.path"),
        metadata_path=_resolve_project_path(root, config.metadata_file, "metadata_file"),
        contrasts_path=_resolve_project_path(root, config.contrasts_file, "contrasts_file"),
    )


def _template_text(name: str) -> str:
    return resources.files("rnaseq.templates").joinpath(name).read_text(encoding="utf-8")


def _project_yaml(
    project_id: str,
    species: Species,
    preset: Preset,
    design_type: DesignType,
    input_type: InputType,
    layout: SequencingLayout | None,
) -> str:
    formula = "~ subject_id + condition" if design_type is DesignType.PAIRED else "~ condition"
    content = {
        "schema_version": SUPPORTED_SCHEMA_VERSION,
        "project": {
            "id": project_id,
            "pipeline": "bulk_rnaseq",
            "preset": preset.value,
        },
        "organism": {"species": species.value},
        "input": {
            "type": input_type.value,
            "path": "input/fastq" if input_type is InputType.FASTQ else "input/counts.csv",
        },
        "design": {"type": design_type.value, "formula": formula},
        "metadata_file": "metadata.csv",
        "contrasts_file": "contrasts.csv",
        "upstream": (
            {
                "engine": "nfcore_rnaseq",
                "pipeline_version": NFCORE_RNASEQ_VERSION,
                "aligner": None,
                "strandedness": "auto",
                "quantification": {"method": "salmon"},
            }
            if input_type is InputType.FASTQ
            else {
                "engine": "external",
                "provider": "external_provider",
                "quantification_method": "unknown",
            }
        ),
        "reference": {"source": "igenomes", "genome": None},
        "thresholds": {"padj": 0.05, "abs_log2fc": 1.0},
        "analysis": {"enrichment": []},
    }
    if input_type is InputType.FASTQ:
        content["input"]["layout"] = (layout or SequencingLayout.PAIRED_END).value
        content["input"]["preprocessing"] = "raw"
    return yaml.safe_dump(content, sort_keys=False, allow_unicode=True)


def create_project(
    *,
    project_name: str,
    destination: Path | str,
    species: Species,
    preset: Preset,
    design_type: DesignType,
    input_type: InputType = InputType.FASTQ,
    layout: SequencingLayout | None = SequencingLayout.PAIRED_END,
) -> Path:
    """Create a new project atomically and return its final path."""

    try:
        ProjectInfo(id=project_name, pipeline="bulk_rnaseq", preset=preset)
    except ValidationError as exc:
        raise ProjectCreationError(_format_pydantic_errors(exc)) from exc

    destination_path = Path(destination).expanduser().resolve()
    destination_path.mkdir(parents=True, exist_ok=True)
    # `rnaseq new` invoked from an empty directory already named after the
    # project must not create the surprising project/project nesting.
    use_destination_directly = destination_path.name == project_name
    target = destination_path if use_destination_directly else destination_path / project_name
    if target.exists() and (not use_destination_directly or any(target.iterdir())):
        raise ProjectCreationError(f"Target project already exists: {target}")

    staging = Path(tempfile.mkdtemp(prefix=f".{project_name}.tmp-", dir=destination_path.parent))
    try:
        (staging / "input").mkdir()
        if input_type is InputType.FASTQ:
            (staging / "input" / "fastq").mkdir()
        (staging / "planning").mkdir()
        metadata_template = (
            "metadata_paired.csv" if design_type is DesignType.PAIRED else "metadata.csv"
        )
        (staging / "metadata.csv").write_text(
            _template_text(metadata_template), encoding="utf-8", newline="\n"
        )
        (staging / "contrasts.csv").write_text(
            _template_text("contrasts.csv"), encoding="utf-8", newline="\n"
        )
        (staging / "project.yaml").write_text(
            _project_yaml(project_name, species, preset, design_type, input_type, layout),
            encoding="utf-8",
            newline="\n",
        )
        if target.exists():
            target.rmdir()
        staging.replace(target)
    except Exception:
        shutil.rmtree(staging, ignore_errors=True)
        raise
    return target
