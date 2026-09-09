"""Command-line interface for the FASTQ and raw-count Milestone 1 control plane."""

from __future__ import annotations

import csv
import sys
from collections.abc import Iterator
from contextlib import contextmanager
from pathlib import Path
from typing import Any

import typer

from rnaseq.errors import ExecutionPreflightError, ProjectCreationError, UpstreamExecutionError
from rnaseq.execution import (
    detect_local_resource_capacity, doctor_checks, load_run_states,
    suggested_local_resources, validate_local_execution_budget,
)
from rnaseq.service import execute_service_run, prepare_service_run, sanitize_completed_delivery, validate_case_id
from rnaseq.models import DesignType, FastqPreprocessing, InputType, PIPELINE_VERSION, Preset, SequencingLayout, Species
from rnaseq.planner import generate_plan
from rnaseq.project import create_project
from rnaseq.references import (
    LocalReferenceError,
    LocalReference,
    ReferenceAdoptionError,
    ReferencePreparationError,
    adopt_local_salmon_index,
    compatible_registered_references,
    prepare_local_reference,
    prepare_local_hisat2_reference,
    reference_registry_path,
    register_local_reference,
)
from rnaseq.validators import ValidationReport, validate_project

app = typer.Typer(
    name="rnaseq",
    help="Create, validate, and plan versioned FASTQ or raw-count bulk RNA-seq projects.",
    no_args_is_help=True,
)
reference_app = typer.Typer(help="Prepare and inspect managed local references.", no_args_is_help=True)
app.add_typer(reference_app, name="reference")


def render_validation_report(report: ValidationReport) -> str:
    lines = ["Bulk RNA-seq Project Validation", "=" * 32]
    config = report.config
    if config is not None:
        lines.extend(
            [
                "Project",
                "-------",
                f"ID: {config.project.id}",
                f"Schema version: {config.schema_version}",
                f"Pipeline version: {PIPELINE_VERSION}",
                f"Preset: {config.project.preset.value}",
                f"Species: {config.organism.species.value}",
            ]
        )
    if report.counts is not None:
        lines.extend(
            [
                "Input",
                "-----",
                "Type: raw_counts",
                f"Genes: {report.counts.gene_count}",
                f"Samples: {len(report.counts.sample_ids)}",
                f"All-zero genes: {report.counts.all_zero_genes}",
            ]
        )
    if report.fastq is not None:
        lines.extend(
            [
                "Input",
                "-----",
                "Type: fastq",
                f"Layout: {report.fastq.layout.value}",
                f"Biological samples: {len(report.fastq.sample_ids)}",
                f"FASTQ assignments: {len(report.fastq.records)}",
            ]
        )
    if report.local_reference is not None:
        reference = report.local_reference
        lines.extend(
            [
                "Reference",
                "---------",
                "Source: local",
                f"Identity: {reference.species}; {reference.provider} release {reference.release}; {reference.assembly_identity}",
                f"Manifest SHA256: {reference.manifest_sha256}",
                f"Salmon index: {reference.salmon_status}",
            ]
        )
    if report.groups:
        lines.extend(["Groups", "------"])
        for factor, levels in report.groups.items():
            if len(report.groups) > 1:
                lines.append(f"Factor: {factor}")
            for level, count in levels.items():
                lines.append(f"{level:<20} n={count}")
    if config is not None:
        lines.extend(["Design", "------", config.design.formula])
    if report.contrasts is not None and report.contrasts.contrasts:
        lines.extend(["Contrasts", "---------"])
        for contrast in report.contrasts.contrasts:
            lines.extend(
                [
                    contrast.contrast_id,
                    f"{contrast.numerator} vs {contrast.denominator}",
                    "Direction:",
                    f"log2FC > 0 = higher in {contrast.numerator}",
                    f"log2FC < 0 = higher in {contrast.denominator}",
                ]
            )
    lines.extend(["Warnings", "--------"])
    lines.extend(issue.message for issue in report.warnings)
    if not report.warnings:
        lines.append("None")
    if report.errors:
        lines.extend(["Errors", "------"])
        lines.extend(issue.message for issue in report.errors)
    lines.extend(["Validation result:", report.state])
    lines.extend([
        "Execution readiness:",
        "READY" if report.execution_ready else "NOT READY",
    ])
    if report.execution_blockers:
        lines.append("Blocking requirements:")
        lines.extend(f"- {item}" for item in report.execution_blockers)
    return "\n".join(lines)


@reference_app.command("prepare")
def prepare_reference_command(
    reference_root: Path,
    threads: int = typer.Option(4, "--threads", min=1, help="Host-native builder threads."),
) -> None:
    """Optionally build a checksum-bound, host-native Salmon index."""

    try:
        reference = prepare_local_reference(reference_root, threads=threads)
    except (LocalReferenceError, ReferencePreparationError) as exc:
        typer.echo(f"ERROR: {exc}", err=True)
        raise typer.Exit(code=1) from exc
    except (OSError, UnicodeError) as exc:
        typer.echo(f"SYSTEM ERROR: {exc}", err=True)
        raise typer.Exit(code=2) from exc
    assert reference.salmon_index is not None
    typer.echo(f"Prepared local reference: {reference.root}")
    typer.echo(f"Salmon index: {reference.salmon_index}")
    typer.echo("Salmon status: built")


@reference_app.command("register")
def register_reference_command(reference_root: Path) -> None:
    """Validate and register one managed reference for this workstation."""

    try:
        reference = register_local_reference(reference_root)
    except LocalReferenceError as exc:
        typer.echo(f"ERROR: {exc}", err=True)
        raise typer.Exit(code=1) from exc
    except (OSError, UnicodeError) as exc:
        typer.echo(f"SYSTEM ERROR: {exc}", err=True)
        raise typer.Exit(code=2) from exc
    typer.echo(f"Registered managed reference: {reference.root}")
    typer.echo(
        f"Identity: {reference.species}; {reference.provider} release {reference.release}; "
        f"{reference.assembly_identity}"
    )
    typer.echo(f"Registry: {reference_registry_path()}")


@reference_app.command("prepare-hisat2")
def prepare_hisat2_reference_command(
    reference_root: Path,
    threads: int = typer.Option(4, "--threads", min=1, help="Host-native builder threads."),
) -> None:
    """Optionally build a checksum-bound, host-native HISAT2 index."""

    try:
        reference = prepare_local_hisat2_reference(reference_root, threads=threads)
    except (LocalReferenceError, ReferencePreparationError) as exc:
        typer.echo(f"ERROR: {exc}", err=True)
        raise typer.Exit(code=1) from exc
    except (OSError, UnicodeError) as exc:
        typer.echo(f"SYSTEM ERROR: {exc}", err=True)
        raise typer.Exit(code=2) from exc
    assert reference.hisat2_index is not None
    typer.echo(f"Prepared local reference: {reference.root}")
    typer.echo(f"HISAT2 index: {reference.hisat2_index}")
    typer.echo("HISAT2 status: built")


@reference_app.command("adopt-salmon-index")
def adopt_salmon_index_command(
    reference_root: Path,
    index: str = typer.Option(..., "--index", help="Reference-root-relative prebuilt Salmon index directory."),
    transcriptome: str | None = typer.Option(None, "--transcriptome", help="Optional root-relative transcriptome; it must match files.transcript_fasta."),
    strategy: str = typer.Option(..., "--strategy", help="Declared Salmon strategy: transcriptome_only or decoy_aware."),
    validation_artifact: str | None = typer.Option(None, "--validation-artifact", help="Optional root-relative external validation record."),
) -> None:
    """Atomically register an explicitly selected prebuilt Salmon index."""

    try:
        reference = adopt_local_salmon_index(
            reference_root,
            index=index,
            transcriptome=transcriptome,
            strategy=strategy,
            validation_artifact=validation_artifact,
        )
    except (LocalReferenceError, ReferenceAdoptionError) as exc:
        typer.echo(f"ERROR: {exc}", err=True)
        raise typer.Exit(code=1) from exc
    except (OSError, UnicodeError) as exc:
        typer.echo(f"SYSTEM ERROR: {exc}", err=True)
        raise typer.Exit(code=2) from exc
    assert reference.salmon_index is not None
    typer.echo(f"Adopted local reference: {reference.root}")
    typer.echo(f"Salmon strategy: {reference.salmon_strategy_type}")
    typer.echo(f"Salmon index: {reference.salmon_index}")


def _wizard_completion_candidates(prefix: str, choices: list[str]) -> list[str]:
    """Return canonical wizard values that begin with the typed prefix.

    This deliberately implements shell-style prefix matching only: it never
    normalizes, fuzzily matches, or otherwise changes a user's input.
    """

    return [value for value in choices if value.startswith(prefix)]


def _readline_module() -> Any | None:
    """Load optional readline support without making it a CLI dependency."""

    try:
        import readline
    except ImportError:
        return None
    return readline


@contextmanager
def _wizard_tab_completion(choices: list[str]) -> Iterator[None]:
    """Temporarily expose canonical choice values to a terminal completer.

    Python's standard GNU readline and libedit-compatible readline bindings
    normally bind Tab to completion.  We intentionally do not alter that
    binding because readline provides no portable way to recover a user's
    previous binding afterwards.  The completer and its delimiters are restored
    immediately after the prompt instead.
    """

    readline = _readline_module()
    if readline is None or not _is_interactive_terminal():
        yield
        return

    get_completer = getattr(readline, "get_completer", None)
    set_completer = getattr(readline, "set_completer", None)
    get_delimiters = getattr(readline, "get_completer_delims", None)
    set_delimiters = getattr(readline, "set_completer_delims", None)
    if not callable(get_completer) or not callable(set_completer):
        yield
        return

    previous_completer = get_completer()
    previous_delimiters = get_delimiters() if callable(get_delimiters) else None

    def complete(text: str, state: int) -> str | None:
        candidates = _wizard_completion_candidates(text, choices)
        return candidates[state] if state < len(candidates) else None

    set_completer(complete)
    # Keep snake_case canonical values as one word for readline/libedit.
    if previous_delimiters is not None and callable(set_delimiters):
        set_delimiters(previous_delimiters.replace("_", ""))
    try:
        yield
    finally:
        set_completer(previous_completer)
        if previous_delimiters is not None and callable(set_delimiters):
            set_delimiters(previous_delimiters)


def _wizard_choice(label: str, choices: list[tuple[str, str]], *, default: str | None = None) -> str:
    """Prompt one canonical wizard choice, retrying ordinary typing mistakes."""

    typer.echo(label)
    for value, description in choices:
        typer.echo(f"  {value}: {description}")
    allowed = [value for value, _ in choices]
    while True:
        # Keep this a string prompt instead of click.Choice: click raises a
        # BadParameter exception before this interactive wizard can retry just
        # the current question.
        with _wizard_tab_completion(allowed):
            selected = typer.prompt("Choose", default=default)
        if selected in allowed:
            return selected
        typer.echo(f"Invalid choice {selected!r}.")
        typer.echo("Please choose one of: " + ", ".join(allowed) + ".")


def _wizard_positive_integer(label: str, *, default: int) -> int:
    """Prompt one positive integer without exposing a conversion traceback."""

    while True:
        raw = typer.prompt(label, default=str(default))
        try:
            value = int(raw)
        except (TypeError, ValueError):
            value = 0
        if value > 0:
            return value
        typer.echo(f"Invalid value {raw!r}. Please enter a positive integer.")


def _is_interactive_terminal() -> bool:
    """Keep the TTY gate explicit and testable without changing CLI semantics."""

    return sys.stdin.isatty()


def _managed_reference_identity(reference: LocalReference) -> str:
    return f"{reference.provider} {reference.release} / {reference.assembly_identity}"


def _registered_reference_choice(species: Species, backend: str) -> LocalReference | None:
    """Offer only manifest-validated, production registered references."""

    try:
        candidates = compatible_registered_references(species.value, backend)
    except LocalReferenceError:
        # A malformed or stale local registry must not block the established
        # manual-managed, custom, and iGenomes routes below.
        return None
    if not candidates:
        return None
    asset_label = "Salmon index" if backend == "salmon" else "HISAT2 index"
    if len(candidates) == 1:
        reference = candidates[0]
        typer.echo("Reference\n---------")
        typer.echo("Detected managed reference:")
        typer.echo(f"  Species: {reference.species}")
        typer.echo(f"  Provider: {reference.provider}")
        typer.echo(f"  Release: {reference.release}")
        typer.echo(f"  Assembly: {reference.assembly_identity}")
        typer.echo(f"  Root: {reference.root}")
        typer.echo(f"  {asset_label}: available")
        return reference if typer.confirm("Use this reference?", default=True) else None

    typer.echo("Reference\n---------")
    typer.echo("Multiple compatible managed references are registered:")
    for number, reference in enumerate(candidates, start=1):
        typer.echo(f"  {number}. {_managed_reference_identity(reference)}")
    while True:
        selected = typer.prompt("Choose reference number")
        try:
            position = int(selected)
        except (TypeError, ValueError):
            position = 0
        if 1 <= position <= len(candidates):
            return candidates[position - 1]
        typer.echo(f"Please choose a number from 1 to {len(candidates)}.")


def _reference_options(
    *, source: str, species: Species, method: str, local_root: Path | None,
    local_manifest: str | None, fasta: Path | None, gtf: Path | None,
    transcript_fasta: Path | None, salmon_index: Path | None, hisat2_index: Path | None,
) -> dict[str, object]:
    if source == "igenomes":
        if method != "salmon":
            raise ProjectCreationError("iGenomes is supported only by the Salmon backend; choose local or custom for HISAT2 + featureCounts.")
        return {"source": "igenomes", "genome": "GRCh38" if species is Species.HUMAN else "GRCm39"}
    if source == "local":
        if local_root is None:
            raise ProjectCreationError("Managed local reference requires --reference-root.")
        return {"source": "local", "root": str(local_root.expanduser().resolve()), "manifest": local_manifest or "reference_manifest.yaml"}
    if fasta is None or gtf is None:
        raise ProjectCreationError("Custom reference requires --reference-fasta and --reference-gtf.")
    reference: dict[str, object] = {"source": "custom", "fasta": str(fasta), "gtf": str(gtf)}
    if transcript_fasta is not None:
        reference["transcript_fasta"] = str(transcript_fasta)
    if salmon_index is not None:
        reference["salmon_index"] = str(salmon_index)
    if hisat2_index is not None:
        reference["hisat2_index"] = str(hisat2_index)
    return reference


def _new_summary(values: dict[str, Any]) -> None:
    typer.echo("\nProject review")
    typer.echo("--------------")
    for label, value in values.items():
        typer.echo(f"{label}: {value}")


def _metadata_fields(path: Path) -> tuple[dict[str, list[str]], list[str]]:
    """Return actual metadata columns/levels for wizard design choices."""

    try:
        with path.open(encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle)
            columns = list(reader.fieldnames or ())
            if "sample_id" not in columns:
                raise ProjectCreationError("Imported metadata must contain sample_id before design fields can be selected.")
            values = {column: [] for column in columns if column != "sample_id"}
            for row in reader:
                for column in values:
                    value = (row.get(column) or "").strip()
                    if value and value not in values[column]:
                        values[column].append(value)
    except (OSError, UnicodeError, csv.Error) as exc:
        raise ProjectCreationError(f"Cannot read metadata for design choices: {exc}") from exc
    fields = [column for column in columns if column != "sample_id"]
    if not fields:
        raise ProjectCreationError("Imported metadata has no fields beyond sample_id.")
    return values, fields


@app.command("new")
def new_project(
    name: str | None = typer.Option(None, "--name", help="Project identifier (enables non-interactive creation)."),
    destination: Path | None = typer.Option(None, "--destination", help="Parent directory, or an empty directory named after the project."),
    species: str | None = typer.Option(None, "--species", help="human or mouse."),
    input_type: str | None = typer.Option(None, "--input-type", help="fastq or raw_counts."),
    fastq_samplesheet: Path | None = typer.Option(None, "--fastq-samplesheet", help="CSV: sample,fastq_1,fastq_2,strandedness."),
    counts: Path | None = typer.Option(None, "--counts", help="Raw integer count matrix CSV."),
    metadata: Path | None = typer.Option(None, "--metadata", help="Existing metadata CSV to preserve."),
    contrasts: Path | None = typer.Option(None, "--contrasts", help="Existing contrasts CSV to preserve."),
    layout: str | None = typer.Option(None, "--layout", help="paired_end or single_end; inferred from imported samplesheets."),
    preprocessing: str | None = typer.Option(None, "--preprocessing", help="raw or pretrimmed FASTQ."),
    method: str | None = typer.Option(None, "--method", "--backend", help="salmon or hisat2_featurecounts."),
    strandedness: str | None = typer.Option(None, "--strandedness", help="auto, unstranded, forward, or reverse."),
    reference_source: str | None = typer.Option(None, "--reference-source", help="igenomes, local, or custom."),
    reference_root: Path | None = typer.Option(None, "--reference-root", help="Managed local reference root."),
    reference_manifest: str | None = typer.Option(None, "--reference-manifest", help="Managed local manifest relative to its root."),
    reference_fasta: Path | None = typer.Option(None, "--reference-fasta", help="Custom genome FASTA."),
    reference_gtf: Path | None = typer.Option(None, "--reference-gtf", help="Custom annotation GTF."),
    reference_transcript_fasta: Path | None = typer.Option(None, "--reference-transcript-fasta", help="Optional custom transcript FASTA."),
    reference_salmon_index: Path | None = typer.Option(None, "--reference-salmon-index", help="Optional custom Salmon index."),
    reference_hisat2_index: Path | None = typer.Option(None, "--reference-hisat2-index", help="Optional prepared custom HISAT2 index."),
    preset: str | None = typer.Option(None, "--preset", help="qc, L1, or L2."),
    design_type: str | None = typer.Option(None, "--design-type", help="two_group, multi_group, or paired."),
    condition_column: str | None = typer.Option(None, "--condition-column", help="Imported metadata factor used for contrasts and the design formula."),
    covariate: list[str] | None = typer.Option(None, "--covariate", help="Additional imported metadata field; repeat as needed."),
    pairing_column: str | None = typer.Option(None, "--pairing-column", help="Imported metadata field for paired designs."),
    execution_profile: str | None = typer.Option(None, "--execution-profile", help="Only local is supported."),
    cpus: int | None = typer.Option(None, "--cpus", help="Total local Nextflow CPU ceiling, not per-task CPUs."),
    memory_gb: int | None = typer.Option(None, "--memory-gb", help="Total local Nextflow memory ceiling in GiB, not per-task memory."),
    scaffold: bool = typer.Option(False, "--scaffold", help="Create templates only; the project remains incomplete until inputs are supplied."),
    yes: bool = typer.Option(False, "--yes", help="Create without an interactive confirmation."),
) -> None:
    """Create a reviewed RNA-seq project interactively or from explicit flags."""

    noninteractive = any(value is not None for value in (name, destination, species, input_type, fastq_samplesheet, counts, metadata, contrasts, layout, preprocessing, method, strandedness, reference_source, reference_root, reference_manifest, reference_fasta, reference_gtf, reference_transcript_fasta, reference_salmon_index, reference_hisat2_index, preset, design_type, condition_column, covariate, pairing_column, execution_profile, cpus, memory_gb)) or scaffold or yes
    selected_managed_reference: LocalReference | None = None
    try:
        if not noninteractive and not _is_interactive_terminal():
            raise ProjectCreationError(
                "rnaseq new needs an interactive terminal, or explicit non-interactive flags (use --scaffold for templates)."
            )
        if not noninteractive:
            name = typer.prompt("Project name")
            destination = Path(typer.prompt("Destination directory", default="."))
            species = _wizard_choice("Species", [("human", "Homo sapiens"), ("mouse", "Mus musculus")], default="human")
            input_type = _wizard_choice("Input type", [("fastq", "Reads requiring upstream processing"), ("raw_counts", "Integer gene-level raw-count matrix")], default="fastq")
            # Ordinary interactive creation is deliberately scaffold-first.
            # Importing is an explicit, reproducible flags-only operation.
            scaffold = True
            if input_type == "fastq":
                layout = _wizard_choice("Sequencing layout", [("paired_end", "R1 and R2 per lane"), ("single_end", "R1 only")], default="paired_end")
                preprocessing = _wizard_choice("FASTQ preprocessing", [("raw", "Run adapter/quality trimming"), ("pretrimmed", "Keep supplied reads; skip trimming")], default="raw")
                method = _wizard_choice("Quantification backend", [("salmon", "nf-core/rnaseq pseudoalignment (default)"), ("hisat2_featurecounts", "HISAT2 alignment plus gene-level featureCounts")], default="salmon")
                strand_choices = [("auto", "infer with Salmon"), ("unstranded", "no stranded protocol"), ("forward", "forward stranded"), ("reverse", "reverse stranded")]
                if method == "hisat2_featurecounts":
                    strand_choices = strand_choices[1:]
                strandedness = _wizard_choice("Strandedness", strand_choices, default="auto" if method == "salmon" else "unstranded")
                selected_managed_reference = _registered_reference_choice(
                    Species.HUMAN if species == "human" else Species.MOUSE, method
                )
                if selected_managed_reference is not None:
                    reference_source = "local"
                    reference_root = selected_managed_reference.root
                    reference_manifest = str(
                        selected_managed_reference.manifest_path.relative_to(selected_managed_reference.root)
                    )
                else:
                    reference_source = _wizard_choice("Reference", [("igenomes", "managed iGenomes convenience route (Salmon only)"), ("local", "checksum-bound managed local reference"), ("custom", "copy supported custom assets into the project")], default="igenomes" if method == "salmon" else "local")
                if reference_source == "local" and selected_managed_reference is None:
                    reference_root = Path(typer.prompt("Managed reference root"))
                    reference_manifest = typer.prompt("Reference manifest", default="reference_manifest.yaml")
                elif reference_source == "custom":
                    reference_fasta = Path(typer.prompt("Custom genome FASTA"))
                    reference_gtf = Path(typer.prompt("Custom annotation GTF"))
                    if method == "hisat2_featurecounts":
                        reference_hisat2_index = Path(typer.prompt("Prepared HISAT2 index directory (blank if not built)", default="")) if typer.confirm("Is a HISAT2 index already prepared?", default=False) else None
            if input_type == "fastq":
                preset = _wizard_choice("Analysis scope", [("qc", "Quantification + technical QC only"), ("L1", "Expression-level QC and exploratory analysis"), ("L2", "L1 plus differential expression and enrichment")], default="L1")
            else:
                preset = _wizard_choice("Analysis scope", [("L1", "Expression-level QC and exploratory analysis"), ("L2", "L1 plus differential expression and enrichment")], default="L2")
            design_type = "two_group" if preset == "qc" else _wizard_choice("Experimental design", [("two_group", "One two-level condition"), ("multi_group", "At least three condition levels"), ("paired", "Paired subject and condition design")], default="two_group")
            if preset != "qc" and metadata is not None:
                fields_by_level, fields = _metadata_fields(metadata)
                typer.echo("Imported metadata fields: " + "; ".join(f"{field}=[{', '.join(fields_by_level[field])}]" for field in fields))
                condition_column = _wizard_choice("Condition/contrast field", [(field, f"levels: {', '.join(fields_by_level[field]) or 'none'}") for field in fields], default=fields[0])
                remaining = [field for field in fields if field != condition_column]
                covariate = []
                if remaining and typer.confirm("Add a covariate to the design?", default=False):
                    covariate.append(_wizard_choice("Covariate", [(field, f"levels: {', '.join(fields_by_level[field]) or 'none'}") for field in remaining], default=remaining[0]))
                if design_type == "paired":
                    candidates = [field for field in fields if field != condition_column]
                    if not candidates:
                        raise ProjectCreationError("Paired design needs a pairing field besides the condition field.")
                    pairing_column = _wizard_choice("Pairing field", [(field, f"levels: {', '.join(fields_by_level[field]) or 'none'}") for field in candidates], default=candidates[0])
            capacity = detect_local_resource_capacity()
            suggested_cpus, suggested_memory = suggested_local_resources(capacity)
            typer.echo("\nLocal execution resources\n-------------------------")
            typer.echo(f"Detected: {capacity.logical_cpus or 'unavailable'} logical CPUs / {capacity.available_memory_gib or capacity.total_memory_gib or 'unavailable'} GiB memory")
            typer.echo(f"Suggested: {suggested_cpus} CPUs / {suggested_memory} GiB memory")
            cpus = _wizard_positive_integer("CPU limit", default=suggested_cpus)
            memory_gb = _wizard_positive_integer("Memory limit in GiB", default=suggested_memory)
            execution_profile = "local"

        if preset is not None and preset.lower() == "qc" and design_type is None:
            # QC-only does not use a statistical design; keep a schema-valid
            # inert value without asking the user an unrelated question.
            design_type = "two_group"
        if name is None or destination is None or species is None or input_type is None or preset is None or design_type is None:
            raise ProjectCreationError("Non-interactive creation requires --name, --destination, --species, --input-type, --preset, and --design-type (or --scaffold with explicit choices).")
        normalized_species = {"human": Species.HUMAN, "mouse": Species.MOUSE}.get(species.lower())
        if normalized_species is None:
            raise ProjectCreationError("--species must be human or mouse.")
        normalized_input = InputType(input_type)
        normalized_preset = Preset.QC if preset.lower() == "qc" else Preset(preset)
        normalized_design = DesignType(design_type)
        formula: str | None = None
        if normalized_preset is not Preset.QC and metadata is not None and condition_column is not None:
            fields_by_level, available_fields = _metadata_fields(metadata)
            selected_covariates = covariate or []
            required_fields = [*selected_covariates, condition_column]
            if normalized_design is DesignType.PAIRED:
                if pairing_column is None:
                    raise ProjectCreationError("Paired imported design requires --pairing-column.")
                required_fields.insert(0, pairing_column)
            missing_fields = [field for field in required_fields if field not in available_fields]
            if missing_fields:
                raise ProjectCreationError("Selected design field(s) are not in imported metadata: " + ", ".join(missing_fields))
            # Formula order is deterministic: pairing first, then covariates,
            # with the contrast factor last. Pairing identity itself is stored
            # explicitly and never inferred from this ordering.
            formula = "~ " + " + ".join(dict.fromkeys(required_fields))
        normalized_layout = SequencingLayout(layout) if layout else (None if fastq_samplesheet else SequencingLayout.PAIRED_END)
        normalized_preprocessing = FastqPreprocessing(preprocessing or "raw")
        normalized_method = method or "salmon"
        if normalized_method not in {"salmon", "hisat2_featurecounts"}:
            raise ProjectCreationError("--method must be salmon or hisat2_featurecounts.")
        normalized_strand = strandedness or ("auto" if normalized_method == "salmon" else "unstranded")
        if normalized_strand not in {"auto", "unstranded", "forward", "reverse"}:
            raise ProjectCreationError("--strandedness must be auto, unstranded, forward, or reverse.")
        if normalized_preset is Preset.QC and normalized_input is not InputType.FASTQ:
            raise ProjectCreationError("The qc preset is available only with --input-type fastq.")
        if noninteractive and not scaffold:
            if normalized_input is InputType.FASTQ and fastq_samplesheet is None:
                raise ProjectCreationError("FASTQ import requires --fastq-samplesheet, or use --scaffold.")
            if normalized_input is InputType.RAW_COUNTS and (counts is None or metadata is None or contrasts is None):
                raise ProjectCreationError("Raw-count import requires --counts, --metadata, and --contrasts, or use --scaffold.")
            if normalized_input is InputType.FASTQ and normalized_preset is not Preset.QC and (metadata is None or contrasts is None):
                raise ProjectCreationError("Analytical FASTQ import requires --metadata and --contrasts, or use --scaffold.")
        if normalized_input is InputType.RAW_COUNTS:
            normalized_method, normalized_strand, normalized_layout = "salmon", "auto", None
            reference = {"source": "igenomes", "genome": None}
        else:
            source = reference_source or ("igenomes" if normalized_method == "salmon" else "local")
            if source not in {"igenomes", "local", "custom"}:
                raise ProjectCreationError("--reference-source must be igenomes, local, or custom.")
            reference = _reference_options(source=source, species=normalized_species, method=normalized_method, local_root=reference_root, local_manifest=reference_manifest, fasta=reference_fasta, gtf=reference_gtf, transcript_fasta=reference_transcript_fasta, salmon_index=reference_salmon_index, hisat2_index=reference_hisat2_index)
        if execution_profile not in {None, "local"}:
            raise ProjectCreationError("Only --execution-profile local is supported.")
        selected_cpus, selected_memory = cpus or 8, memory_gb or 12
        validate_local_execution_budget(selected_cpus, selected_memory, detect_local_resource_capacity())
        execution = {"profile": "local", "max_cpus": selected_cpus, "max_memory_gb": selected_memory}
        review_reference = (
            _managed_reference_identity(selected_managed_reference)
            if selected_managed_reference is not None
            else reference.get("source")
        )
        review = {"Project": name, "Destination": destination, "Species": normalized_species.value, "Input": normalized_input.value, "Input handling": "scaffold — add data after project creation" if scaffold else "import supplied inputs", "Scope": normalized_preset.value, "Design": "not applicable for QC" if normalized_preset is Preset.QC else normalized_design.value, "Backend": normalized_method if normalized_input is InputType.FASTQ else "external raw counts", "Reference": review_reference, "Execution": f"local, {selected_cpus} CPUs / {selected_memory} GiB"}
        _new_summary(review)
        if not yes and not typer.confirm("Create this project?", default=True):
            if not noninteractive and typer.confirm("Revise choices?", default=True):
                # No filesystem action has occurred yet; restart the wizard
                # rather than trying to mutate a partially written project.
                return new_project(
                    name=None, destination=None, species=None, input_type=None,
                    fastq_samplesheet=None, counts=None, metadata=None, contrasts=None,
                    layout=None, preprocessing=None, method=None, strandedness=None,
                    reference_source=None, reference_root=None, reference_manifest=None,
                    reference_fasta=None, reference_gtf=None, reference_transcript_fasta=None,
                    reference_salmon_index=None, reference_hisat2_index=None, preset=None,
                    design_type=None, condition_column=None, covariate=None, pairing_column=None,
                    execution_profile=None, cpus=None, memory_gb=None,
                    scaffold=False, yes=False,
                )
            typer.echo("Project creation cancelled; no project was written.")
            return
        target = create_project(project_name=name, destination=destination, species=normalized_species, preset=normalized_preset, design_type=normalized_design, input_type=normalized_input, layout=normalized_layout, preprocessing=normalized_preprocessing, strandedness=normalized_strand, quantification_method=normalized_method, reference=reference, fastq_samplesheet=fastq_samplesheet, counts_file=counts, metadata_file=metadata, contrasts_file=contrasts, scaffold=scaffold, formula=formula, pairing_column=pairing_column, execution=execution)
    except (ValueError, ProjectCreationError) as exc:
        typer.echo(f"ERROR: {exc}", err=True)
        raise typer.Exit(code=1) from exc
    except (OSError, UnicodeError) as exc:
        typer.echo(f"SYSTEM ERROR: {exc}", err=True)
        raise typer.Exit(code=2) from exc

    typer.echo(f"Created project: {target}")
    if scaffold:
        typer.echo("This is an incomplete scaffold. Add valid inputs, metadata, and contrasts, then run: rnaseq validate PROJECT")
    elif normalized_input is InputType.FASTQ:
        typer.echo("Imported FASTQs were copied into input/fastq; run rnaseq validate PROJECT to check readiness.")
    else:
        typer.echo("Imported raw counts and supplied metadata/contrasts were preserved; run rnaseq validate PROJECT.")


@app.command("validate")
def validate_command(project_dir: Path) -> None:
    """Validate a FASTQ or raw-count project without running analysis."""

    try:
        report = validate_project(project_dir)
    except (OSError, UnicodeError) as exc:
        typer.echo(f"SYSTEM ERROR: {exc}", err=True)
        raise typer.Exit(code=2) from exc
    typer.echo(render_validation_report(report))
    if not report.is_valid:
        raise typer.Exit(code=1)


@app.command("plan")
def plan_command(project_dir: Path) -> None:
    """Validate a project and generate deterministic planning artifacts."""

    try:
        report = validate_project(project_dir)
        typer.echo(render_validation_report(report))
        if not report.is_valid:
            raise typer.Exit(code=1)
        artifacts = generate_plan(report)
    except typer.Exit:
        raise
    except (OSError, UnicodeError) as exc:
        typer.echo(f"SYSTEM ERROR: {exc}", err=True)
        raise typer.Exit(code=2) from exc
    for artifact in artifacts:
        typer.echo(f"Planning artifact: {artifact}")


@app.command("run")
def run_command(
    project_dir: Path,
    case_id: str = typer.Option(..., "--case-id", help="Required filesystem-safe client case identifier."),
    profile: str = typer.Option("local", "--profile", help="Execution profile; only local is implemented."),
    reuse_upstream: str | None = typer.Option(None, "--reuse-upstream", help="Optional compatible CASE-ID/RUN-ID upstream reuse source."),
    yes: bool = typer.Option(False, "--yes", help="Authorize execution without an interactive prompt."),
) -> None:
    """Execute one immutable client case through upstream and downstream workflows."""

    try:
        report = validate_project(project_dir)
        typer.echo(render_validation_report(report))
        validate_case_id(case_id)
        prepare_service_run(report, profile=profile)
        assert report.config is not None
        if not yes:
            typer.echo(f"Project: {report.config.project.id}")
            typer.echo(f"Case ID: {case_id}")
            typer.echo(f"Input type: {report.config.input.type.value}")
            if report.fastq is not None:
                method = report.config.upstream.quantification.method if report.config.upstream.quantification else "salmon"
                typer.echo(
                    f"Pipeline: nf-core/rnaseq {report.config.upstream.pipeline_version}"
                    if method == "salmon" else "Pipeline: HISAT2 + featureCounts"
                )
                typer.echo(f"Samples: {len(report.fastq.sample_ids)}")
            typer.echo(f"Profile: {profile}")
            typer.echo(
                "An immutable case run and technical-QC delivery package will be created."
                if report.config.project.preset is Preset.QC
                else "An immutable case run, downstream workflow, report, and delivery package will be created."
            )
            if not typer.confirm("Proceed?", default=False):
                typer.echo("Execution cancelled; no run directory was created.")
                raise typer.Exit(code=0)
        result = execute_service_run(report, case_id=case_id, profile=profile, reuse_upstream=reuse_upstream)
    except typer.Exit:
        raise
    except ExecutionPreflightError as exc:
        typer.echo(f"ERROR: {exc}", err=True)
        raise typer.Exit(code=1) from exc
    except UpstreamExecutionError as exc:
        typer.echo(f"ERROR: {exc}", err=True)
        raise typer.Exit(code=1) from exc
    except (OSError, UnicodeError) as exc:
        typer.echo(f"SYSTEM ERROR: {exc}", err=True)
        raise typer.Exit(code=2) from exc
    typer.echo(f"Case run: SUCCESS ({result.case_id}/{result.run_id})")
    typer.echo(f"Run directory: {result.run_dir}")
    typer.echo(f"Delivery package: {result.run_dir / 'delivery'}")


@app.command("status")
def status_command(project_dir: Path) -> None:
    """Show persisted case/run states without inspecting live processes."""

    states = load_run_states(project_dir.resolve())
    if not states:
        typer.echo("No recorded case runs.")
        return
    for state in states:
        typer.echo(f"Run ID: {state.get('run_id', 'unknown')}")
        typer.echo(f"Status: {state.get('status', 'unknown')}")
        typer.echo(f"Case: {state.get('case_id', 'legacy')}")
        typer.echo(f"Profile: {state.get('profile', 'local')}")
        typer.echo(f"Started: {state.get('started_at')}")
        typer.echo(f"Completed: {state.get('completed_at')}")
        typer.echo(f"Upstream outputs: {'available' if state.get('handoff_available') else 'unavailable'}")
        typer.echo(f"Delivery package: {'available' if state.get('delivery_available') else 'unavailable'}")
        typer.echo(f"Run directory: {state.get('run_dir')}")


@app.command("sanitize-delivery")
def sanitize_delivery_command(run_dir: Path) -> None:
    """Remove AppleDouble sidecars from exactly one completed run's delivery tree."""

    try:
        delivery = sanitize_completed_delivery(run_dir)
    except ExecutionPreflightError as exc:
        typer.echo(f"ERROR: {exc}", err=True)
        raise typer.Exit(code=1) from exc
    except (OSError, UnicodeError, UpstreamExecutionError) as exc:
        typer.echo(f"SYSTEM ERROR: {exc}", err=True)
        raise typer.Exit(code=2) from exc
    typer.echo(f"Delivery package sanitized: {delivery}")


@app.command("analyze")
def analyze_command(
    project_dir: Path,
    level: str = typer.Option(..., "--level", help="Legacy option retained for a controlled migration."),
    run_id: str | None = typer.Option(None, "--run-id", help="Legacy option retained for a controlled migration."),
    enrichment: str | None = typer.Option(None, "--enrichment", help="Legacy option retained for a controlled migration."),
) -> None:
    """Compatibility command; direct Python-to-R execution is retired."""

    _ = (project_dir, level, run_id, enrichment)
    typer.echo(
        "DEPRECATED: 'rnaseq analyze' no longer executes statistical backends directly. "
        "Use 'rnaseq run PROJECT --case-id CASE-ID' so downstream R tasks are orchestrated by Nextflow.",
        err=True,
    )
    raise typer.Exit(code=2)


@app.command("doctor")
def doctor_command(project_dir: Path | None = typer.Argument(None, help="Optional project for adopted-reference readiness.")) -> None:
    """Report non-mutating execution prerequisites; nothing is installed automatically."""

    for item in doctor_checks(project_dir):
        typer.echo(f"{item.name}: {item.verdict} — {item.detail}")


if __name__ == "__main__":
    app()
