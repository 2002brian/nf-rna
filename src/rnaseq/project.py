"""Project creation and strict configuration loading."""

from __future__ import annotations

import shutil
import tempfile
import csv
import re
from dataclasses import dataclass
from importlib import resources
from pathlib import Path
from typing import Any

import yaml
from pydantic import ValidationError

from rnaseq.errors import ProjectConfigError, ProjectCreationError
from rnaseq.models import (
    DesignType,
    FastqPreprocessing,
    InputType,
    NFCORE_RNASEQ_VERSION,
    Preset,
    ProjectConfig,
    ProjectInfo,
    SequencingLayout,
    Species,
    SUPPORTED_SCHEMA_VERSION,
)

FASTQ_SAMPLESHEET_HEADER = ("sample", "fastq_1", "fastq_2", "strandedness")
FASTQ_NAME_PATTERN = re.compile(
    r"^(?P<sample>.+?)(?:_L(?P<lane>\d{3}))?_R(?P<read>[12])(?:_\d{3})?\.(?:fastq|fq)\.gz$"
)


@dataclass(frozen=True)
class ImportedFastq:
    """Validated row from the portable, nf-core-compatible FASTQ samplesheet."""

    sample: str
    fastq_1: Path
    fastq_2: Path | None
    strandedness: str


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
    preprocessing: FastqPreprocessing = FastqPreprocessing.RAW,
    strandedness: str = "auto",
    quantification_method: str = "salmon",
    reference: dict[str, object] | None = None,
    execution: dict[str, object] | None = None,
    formula: str | None = None,
    pair_id: str | None = None,
) -> str:
    formula = formula or (
        f"~ {pair_id or 'patient'} + condition"
        if design_type is DesignType.PAIRED_TWO_GROUP else "~ condition"
    )
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
        "design": {
            "type": design_type.value,
            "formula": formula,
            **({"pair_id": pair_id or "patient"} if design_type is DesignType.PAIRED_TWO_GROUP else {}),
        },
        "metadata_file": "metadata.csv",
        "contrasts_file": "contrasts.csv",
        "upstream": (
            {
                "engine": "nfcore_rnaseq",
                "pipeline_version": NFCORE_RNASEQ_VERSION,
                "aligner": None,
                "strandedness": strandedness,
                "quantification": {"method": quantification_method},
            }
            if input_type is InputType.FASTQ
            else {
                "engine": "external",
                "provider": "external_provider",
                "quantification_method": "unknown",
            }
        ),
        "reference": reference or {"source": "igenomes", "genome": None},
        "runtime": {"control_plane_image": "rnaseq-control-plane:latest"},
        "execution": execution or {"profile": "local", "max_cpus": 8, "max_memory_gb": 12},
        "thresholds": {"padj": 0.05, "abs_log2fc": 1.0},
        "analysis": {"enrichment": []},
    }
    if input_type is InputType.FASTQ:
        content["input"]["layout"] = (layout or SequencingLayout.PAIRED_END).value
        content["input"]["preprocessing"] = preprocessing.value
    return yaml.safe_dump(content, sort_keys=False, allow_unicode=True)


def _copy_file(source: Path, destination: Path) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(source, destination)


def _read_fastq_samplesheet(path: Path) -> tuple[tuple[ImportedFastq, ...], SequencingLayout, str]:
    """Read a strict, portable FASTQ samplesheet before creating anything.

    The project model has one run-wide strandedness setting.  A samplesheet
    with mixed values is therefore rejected rather than silently flattening
    biological metadata into a default.
    """

    if not path.is_file():
        raise ProjectCreationError(f"FASTQ samplesheet is not a file: {path}")
    try:
        with path.open(encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle)
            if tuple(reader.fieldnames or ()) != FASTQ_SAMPLESHEET_HEADER:
                raise ProjectCreationError(
                    "FASTQ samplesheet columns must be exactly: sample,fastq_1,fastq_2,strandedness."
                )
            rows = list(reader)
    except (OSError, UnicodeError, csv.Error) as exc:
        raise ProjectCreationError(f"Cannot read FASTQ samplesheet {path}: {exc}") from exc
    if not rows:
        raise ProjectCreationError("FASTQ samplesheet has no sample rows.")

    imported: list[ImportedFastq] = []
    names: set[str] = set()
    sample_strands: dict[str, str] = {}
    layout: SequencingLayout | None = None
    for number, row in enumerate(rows, start=2):
        sample = (row.get("sample") or "").strip()
        first_text, second_text = (row.get("fastq_1") or "").strip(), (row.get("fastq_2") or "").strip()
        strand = (row.get("strandedness") or "").strip()
        if not sample or not first_text or not strand:
            raise ProjectCreationError(f"FASTQ samplesheet row {number} requires sample, fastq_1, and strandedness.")
        if strand not in {"auto", "unstranded", "forward", "reverse"}:
            raise ProjectCreationError(f"FASTQ samplesheet row {number} has unsupported strandedness: {strand!r}.")
        first = Path(first_text).expanduser()
        second = Path(second_text).expanduser() if second_text else None
        if not first.is_absolute():
            first = path.parent / first
        if second is not None and not second.is_absolute():
            second = path.parent / second
        first, second = first.resolve(), second.resolve() if second else None
        if not first.is_file() or (second is not None and not second.is_file()):
            raise ProjectCreationError(f"FASTQ samplesheet row {number} refers to a missing FASTQ file.")
        first_match = FASTQ_NAME_PATTERN.fullmatch(first.name)
        second_match = FASTQ_NAME_PATTERN.fullmatch(second.name) if second else None
        if first_match is None or first_match.group("read") != "1":
            raise ProjectCreationError(f"FASTQ samplesheet row {number} fastq_1 must have an unambiguous _R1 .fastq.gz/.fq.gz name.")
        row_layout = SequencingLayout.PAIRED_END if second else SequencingLayout.SINGLE_END
        if second and (second_match is None or second_match.group("read") != "2" or second_match.group("sample") != first_match.group("sample") or second_match.group("lane") != first_match.group("lane")):
            raise ProjectCreationError(f"FASTQ samplesheet row {number} fastq_2 must be the matching _R2 file for fastq_1.")
        if first_match.group("sample") != sample:
            raise ProjectCreationError(f"FASTQ samplesheet row {number} sample {sample!r} does not match FASTQ filename sample {first_match.group('sample')!r}.")
        if layout is not None and layout is not row_layout:
            raise ProjectCreationError("FASTQ samplesheet mixes single-end and paired-end rows; create separate projects.")
        layout = row_layout
        previous = sample_strands.setdefault(sample, strand)
        if previous != strand:
            raise ProjectCreationError(f"FASTQ samplesheet gives sample {sample!r} conflicting strandedness values: {previous!r} and {strand!r}.")
        for item in (first, second):
            if item is not None:
                if item.name in names:
                    raise ProjectCreationError(f"Imported FASTQ basename is duplicated: {item.name}. Rename source files before import.")
                names.add(item.name)
        imported.append(ImportedFastq(sample, first, second, strand))
    strands = set(sample_strands.values())
    if len(strands) != 1:
        raise ProjectCreationError("FASTQ samplesheet has mixed strandedness; this schema supports one explicit run-wide value.")
    assert layout is not None
    return tuple(imported), layout, strands.pop()


def create_project(
    *,
    project_name: str,
    destination: Path | str,
    species: Species,
    preset: Preset,
    design_type: DesignType,
    input_type: InputType = InputType.FASTQ,
    layout: SequencingLayout | None = SequencingLayout.PAIRED_END,
    preprocessing: FastqPreprocessing = FastqPreprocessing.RAW,
    strandedness: str = "auto",
    quantification_method: str = "salmon",
    reference: dict[str, object] | None = None,
    fastq_samplesheet: Path | str | None = None,
    counts_file: Path | str | None = None,
    metadata_file: Path | str | None = None,
    contrasts_file: Path | str | None = None,
    scaffold: bool = False,
    formula: str | None = None,
    pair_id: str | None = None,
    execution: dict[str, object] | None = None,
) -> Path:
    """Create a new project atomically and return its final path."""

    if design_type is DesignType.PAIRED_TWO_GROUP and (pair_id is None or not pair_id.strip()):
        raise ProjectCreationError(
            "paired_two_group project creation requires an explicit biological pair_id column."
        )

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

    imported_fastq: tuple[ImportedFastq, ...] = ()
    source_fastq_samplesheet: Path | None = None
    if input_type is InputType.FASTQ and fastq_samplesheet is not None:
        source_fastq_samplesheet = Path(fastq_samplesheet).expanduser().resolve()
        imported_fastq, imported_layout, imported_strandedness = _read_fastq_samplesheet(source_fastq_samplesheet)
        if layout is not None and layout is not imported_layout:
            raise ProjectCreationError(f"Requested layout {layout.value} conflicts with imported {imported_layout.value} FASTQs.")
        layout, strandedness = imported_layout, imported_strandedness
    if input_type is InputType.FASTQ and quantification_method == "hisat2_featurecounts" and strandedness == "auto":
        raise ProjectCreationError("HISAT2 + featureCounts requires explicit strandedness: unstranded, forward, or reverse.")
    source_counts = Path(counts_file).expanduser().resolve() if counts_file else None
    source_metadata = Path(metadata_file).expanduser().resolve() if metadata_file else None
    source_contrasts = Path(contrasts_file).expanduser().resolve() if contrasts_file else None
    for label, source in (("count matrix", source_counts), ("metadata", source_metadata), ("contrasts", source_contrasts)):
        if source is not None and not source.is_file():
            raise ProjectCreationError(f"Imported {label} is not a file: {source}")
    if input_type is InputType.RAW_COUNTS and not scaffold and source_counts is None:
        raise ProjectCreationError("Raw-count import requires --counts, or use --scaffold to create an incomplete project.")
    resolved_reference = dict(reference) if reference is not None else {"source": "igenomes", "genome": None}
    custom_reference_assets: list[tuple[Path, str, str]] = []
    if resolved_reference.get("source") == "custom":
        for key in ("fasta", "gtf", "transcript_fasta", "salmon_index", "hisat2_index", "hisat2_splice_sites"):
            configured = resolved_reference.get(key)
            if configured is None:
                continue
            source = Path(str(configured)).expanduser().resolve()
            if not source.exists() or not (source.is_file() or source.is_dir()):
                raise ProjectCreationError(f"Custom reference {key} does not exist: {source}")
            relative = f"reference/{key}/{source.name}"
            custom_reference_assets.append((source, key, relative))
            resolved_reference[key] = relative

    staging = Path(tempfile.mkdtemp(prefix=f".{project_name}.tmp-", dir=destination_path.parent))
    try:
        (staging / "input").mkdir()
        if input_type is InputType.FASTQ:
            (staging / "input" / "fastq").mkdir()
        (staging / "planning").mkdir()
        metadata_template = (
            "metadata_paired.csv" if design_type is DesignType.PAIRED_TWO_GROUP else "metadata.csv"
        )
        if source_metadata is not None:
            _copy_file(source_metadata, staging / "metadata.csv")
        elif design_type is DesignType.PAIRED_TWO_GROUP:
            (staging / "metadata.csv").write_text(
                f"sample_id,{pair_id or 'patient'},condition\n", encoding="utf-8", newline="\n"
            )
        else:
            (staging / "metadata.csv").write_text(_template_text(metadata_template), encoding="utf-8", newline="\n")
        if source_contrasts is not None:
            _copy_file(source_contrasts, staging / "contrasts.csv")
        else:
            (staging / "contrasts.csv").write_text(_template_text("contrasts.csv"), encoding="utf-8", newline="\n")
        if source_counts is not None:
            _copy_file(source_counts, staging / "input" / "counts.csv")
        for record in imported_fastq:
            _copy_file(record.fastq_1, staging / "input" / "fastq" / record.fastq_1.name)
            if record.fastq_2 is not None:
                _copy_file(record.fastq_2, staging / "input" / "fastq" / record.fastq_2.name)
        if source_fastq_samplesheet is not None:
            _copy_file(source_fastq_samplesheet, staging / "planning" / "imported_fastq_samplesheet.csv")
        for source, _key, relative in custom_reference_assets:
            target_reference = staging / relative
            if source.is_dir():
                shutil.copytree(source, target_reference)
            else:
                _copy_file(source, target_reference)
        (staging / "project.yaml").write_text(
            _project_yaml(
                project_name, species, preset, design_type, input_type, layout,
                preprocessing, strandedness, quantification_method,
                resolved_reference, execution, formula, pair_id,
            ),
            encoding="utf-8",
            newline="\n",
        )
        if not scaffold:
            # Validate the complete imported project while it is still private
            # staging state.  An invalid import therefore never leaves a
            # project-shaped partial directory at the requested destination.
            from rnaseq.validators import validate_project

            report = validate_project(staging)
            if not report.is_valid:
                details = "; ".join(issue.message for issue in report.errors)
                raise ProjectCreationError("Imported project failed strict validation: " + details)
        if target.exists():
            target.rmdir()
        staging.replace(target)
    except Exception:
        shutil.rmtree(staging, ignore_errors=True)
        raise
    return target
