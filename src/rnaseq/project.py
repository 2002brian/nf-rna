"""Project creation and strict configuration loading."""

from __future__ import annotations

import shutil
import tempfile
import csv
import io
import re
from collections.abc import Iterable
from dataclasses import dataclass
from importlib import resources
from pathlib import Path
from typing import Any

import yaml
from pydantic import ValidationError

from rnaseq.errors import ProjectConfigError, ProjectCreationError
from rnaseq.models import (
    DEFAULT_L2_ENRICHMENT,
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
    default_annotation,
)

FASTQ_SAMPLESHEET_HEADER = ("sample", "fastq_1", "fastq_2", "strandedness")
# Standardised column names used by the template-first lifecycle.
CONDITION_COLUMN = "condition"
BATCH_COLUMN = "batch"
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


def derive_design(
    design_type: DesignType, pair_id: str | None = None, batch: bool = False
) -> tuple[str, dict[str, str], tuple[str, ...]]:
    """Mechanically derive the provisional design from the planning choices.

    Returns ``(formula, variable types, metadata columns after sample_id)``.
    The formula lists the pairing block first, then batch, with the contrast
    factor ``condition`` last; the metadata header lists ``condition`` before
    ``batch`` so it reads as the user's primary annotation first.  Nothing here
    inspects data: levels are discovered later from the completed metadata.
    """

    paired = design_type is DesignType.PAIRED_TWO_GROUP
    pair = (pair_id or "patient") if paired else None
    formula_terms = [*([pair] if pair else []), *([BATCH_COLUMN] if batch else []), CONDITION_COLUMN]
    columns = (*([pair] if pair else []), CONDITION_COLUMN, *([BATCH_COLUMN] if batch else []))
    return "~ " + " + ".join(formula_terms), {term: "categorical" for term in formula_terms}, columns


def metadata_template_header(design_type: DesignType, pair_id: str | None = None, batch: bool = False) -> str:
    """Return the header-only metadata template for a design (with newline)."""

    return ",".join(("sample_id", *derive_design(design_type, pair_id, batch)[2])) + "\n"


def resolve_enrichment_defaults(
    preset: Preset, species: Species, enrichment: tuple[str, ...] | None, gene_id_type: str | None,
) -> tuple[tuple[str, ...], dict[str, str] | None]:
    """Return ``(analysis.enrichment, annotation)`` for a new project.

    L2 selects every supported method unless the caller opts out with ``()``.
    Enrichment is never dropped to make creation succeed: a missing gene ID
    type or an unsupported organism is an error naming the remedy.
    """

    if preset is Preset.QC and (enrichment or gene_id_type):
        raise ProjectCreationError("Enrichment and --gene-id-type do not apply to the qc preset.")
    if preset is Preset.L1 and enrichment is not None:
        raise ProjectCreationError(
            "Enrichment is selected when the project is configured for L2 (rnaseq new --preset L2, "
            "or rnaseq design after L1); --enrichment/--no-enrichment does not apply to L1."
        )
    selected = (DEFAULT_L2_ENRICHMENT if enrichment is None else tuple(enrichment)) if preset is Preset.L2 else ()
    if selected and gene_id_type is None:
        raise ProjectCreationError(
            "L2 enables GO ORA, KEGG ORA, and preranked GSEA (GO and KEGG) by default, which need the "
            "gene identifier type used by the count matrix. Pass --gene-id-type ENSEMBL, ENTREZID, or SYMBOL "
            "(Ensembl GTF/managed references use ENSEMBL), or opt out with --no-enrichment."
        )
    try:
        annotation = default_annotation(species, gene_id_type) if gene_id_type is not None else None
    except ValueError as exc:
        raise ProjectCreationError(str(exc)) from exc
    return selected, annotation


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
    design_variables: dict[str, str] | None = None,
    batch: bool = False,
    enrichment: tuple[str, ...] = (),
    annotation: dict[str, str] | None = None,
) -> str:
    if formula is None:
        formula, derived_variables, _columns = derive_design(design_type, pair_id, batch)
        design_variables = design_variables or derived_variables
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
            "variables": design_variables or {
                variable: "categorical"
                for variable in re.findall(r"[A-Za-z_][A-Za-z0-9_.]*", formula.removeprefix("~"))
            },
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
        "execution": execution or {"profile": "local", "max_cpus": "auto", "max_memory_gb": "auto"},
        "thresholds": {"padj": 0.05, "abs_log2fc": 1.0},
        **({"annotation": dict(annotation)} if annotation else {}),
        "analysis": {"enrichment": list(enrichment)},
        # New projects select the current figure style explicitly; projects without
        # this section keep rendering the legacy style.
        "figures": {"preset": "standard", "profile": "general"},
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


def _copy_imported_fastq(records: tuple[ImportedFastq, ...], samplesheet: Path, root: Path, fastq_dir: Path) -> None:
    for record in records:
        _copy_file(record.fastq_1, fastq_dir / record.fastq_1.name)
        if record.fastq_2 is not None:
            _copy_file(record.fastq_2, fastq_dir / record.fastq_2.name)
    _copy_file(samplesheet, root / "planning" / "imported_fastq_samplesheet.csv")


def fastq_samplesheet_text(records: Iterable[Any], strandedness: str) -> str:
    """A portable samplesheet (``sample,fastq_1,fastq_2,strandedness``) for discovered FASTQs.

    ``records`` are :class:`rnaseq.validators.FastqRecord` values, so sample IDs
    and R1/R2 pairing follow the validator's FASTQ naming rule.
    """

    buffer = io.StringIO()
    writer = csv.writer(buffer, lineterminator="\n")
    writer.writerow(FASTQ_SAMPLESHEET_HEADER)
    for record in records:
        writer.writerow([record.sample_id, str(record.fastq_1), str(record.fastq_2) if record.fastq_2 else "", strandedness])
    return buffer.getvalue()


def import_fastq_samplesheet(project_dir: Path | str, samplesheet: Path | str) -> tuple[ImportedFastq, ...]:
    """Copy the FASTQs of a samplesheet into an existing project's FASTQ input directory.

    The same rules as ``rnaseq new --fastq-samplesheet`` apply; in addition the
    layout and strandedness must match project.yaml, and existing FASTQ files
    are never overwritten.
    """

    loaded = load_project(project_dir)
    config = loaded.config
    if config.input.type is not InputType.FASTQ:
        raise ProjectCreationError("This project takes a raw-count matrix, not FASTQ files.")
    source = Path(samplesheet).expanduser().resolve()
    records, layout, strandedness = _read_fastq_samplesheet(source)
    if config.input.layout is not None and layout is not config.input.layout:
        raise ProjectCreationError(f"The samplesheet is {layout.value}, but project.yaml declares {config.input.layout.value}.")
    if strandedness != config.upstream.strandedness:
        raise ProjectCreationError(
            f"The samplesheet strandedness {strandedness!r} differs from project.yaml upstream.strandedness {config.upstream.strandedness!r}."
        )
    target = loaded.input_path
    existing = sorted(item.name for item in target.iterdir() if item.is_file() and not item.name.startswith("._")) if target.is_dir() else []
    if existing:
        raise ProjectCreationError(f"{target} already contains {len(existing)} file(s); remove them before importing another samplesheet.")
    _copy_imported_fastq(records, source, loaded.root, target)
    return records


def import_count_matrix(project_dir: Path | str, counts_file: Path | str, *, replace: bool = False) -> Path:
    """Copy a raw-count matrix to an existing project's configured input path (as ``rnaseq new --counts``)."""

    loaded = load_project(project_dir)
    if loaded.config.input.type is not InputType.RAW_COUNTS:
        raise ProjectCreationError("This project takes FASTQ input, not a count matrix.")
    source = Path(counts_file).expanduser().resolve()
    if not source.is_file():
        raise ProjectCreationError(f"Imported count matrix is not a file: {source}")
    if loaded.input_path.exists() and not replace:
        raise ProjectCreationError(f"{loaded.input_path} already exists; replace it explicitly.")
    _copy_file(source, loaded.input_path)
    return loaded.input_path


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
    design_variables: dict[str, str] | None = None,
    execution: dict[str, object] | None = None,
    batch: bool = False,
    enrichment: tuple[str, ...] | None = None,
    gene_id_type: str | None = None,
) -> Path:
    """Create a new project atomically and return its final path.

    ``batch`` declares a known technical batch for template-first projects: the
    metadata template gains a ``batch`` column and the provisional formula
    adjusts for it.  It cannot be combined with an explicit ``formula`` (the
    advanced path states its own model) or an imported metadata file.

    ``enrichment`` ``None`` selects the preset default: every supported method
    for L2 and none otherwise; ``()`` is an explicit opt-out.  ``gene_id_type``
    records the annotation contract that enrichment requires.
    """

    selected_enrichment, annotation = resolve_enrichment_defaults(preset, species, enrichment, gene_id_type)

    if batch and (formula is not None or metadata_file is not None):
        raise ProjectCreationError(
            "A declared batch generates the template and formula; it cannot be combined with an explicit "
            "formula or imported metadata (name 'batch' as a covariate for imported metadata instead)."
        )

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
        if source_metadata is not None:
            _copy_file(source_metadata, staging / "metadata.csv")
        else:
            (staging / "metadata.csv").write_text(
                metadata_template_header(design_type, pair_id, batch), encoding="utf-8", newline="\n"
            )
        if source_contrasts is not None:
            _copy_file(source_contrasts, staging / "contrasts.csv")
        else:
            (staging / "contrasts.csv").write_text(_template_text("contrasts.csv"), encoding="utf-8", newline="\n")
        if source_counts is not None:
            _copy_file(source_counts, staging / "input" / "counts.csv")
        if source_fastq_samplesheet is not None:
            _copy_imported_fastq(imported_fastq, source_fastq_samplesheet, staging, staging / "input" / "fastq")
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
                resolved_reference, execution, formula, pair_id, design_variables, batch,
                selected_enrichment, annotation,
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
