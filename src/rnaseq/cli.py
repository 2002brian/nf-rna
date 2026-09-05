"""Command-line interface for the FASTQ and raw-count Milestone 1 control plane."""

from __future__ import annotations

import csv
import sys
from pathlib import Path
from typing import Any

import click
import typer

from rnaseq.errors import ExecutionPreflightError, ProjectCreationError, UpstreamExecutionError
from rnaseq.execution import doctor_checks, load_run_states
from rnaseq.service import execute_service_run, prepare_service_run, sanitize_completed_delivery, validate_case_id
from rnaseq.models import DesignType, FastqPreprocessing, InputType, PIPELINE_VERSION, Preset, SequencingLayout, Species
from rnaseq.planner import generate_plan
from rnaseq.project import create_project
from rnaseq.references import (
    LocalReferenceError,
    ReferenceAdoptionError,
    ReferencePreparationError,
    adopt_local_salmon_index,
    prepare_local_reference,
    prepare_local_hisat2_reference,
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
                f"Identity: {reference.species}; {reference.provider} release {reference.release}; {reference.assembly} {reference.assembly_patch}",
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
def prepare_reference_command(reference_root: Path) -> None:
    """Build the checksum-bound decoy-aware Salmon index for one reference root."""

    try:
        reference = prepare_local_reference(reference_root)
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


@reference_app.command("prepare-hisat2")
def prepare_hisat2_reference_command(reference_root: Path) -> None:
    """Build the checksum-bound HISAT2 index for one local reference root."""

    try:
        reference = prepare_local_hisat2_reference(reference_root)
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
    transcriptome: str = typer.Option(..., "--transcriptome", help="Reference-root-relative GTF-derived transcriptome FASTA."),
    strategy: str = typer.Option(..., "--strategy", help="Declared Salmon strategy; currently transcriptome_only."),
    validation_artifact: str = typer.Option(..., "--validation-artifact", help="Reference-root-relative transcript validation JSON."),
) -> None:
    """Atomically adopt an explicitly selected, validated existing Salmon index."""

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


def _wizard_choice(label: str, choices: list[tuple[str, str]], *, default: str | None = None) -> str:
    """Prompt with visible descriptions while retaining scriptable enum values."""

    typer.echo(label)
    for value, description in choices:
        typer.echo(f"  {value}: {description}")
    return typer.prompt("Choose", type=click.Choice([value for value, _ in choices], case_sensitive=True), default=default)


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
    method: str | None = typer.Option(None, "--method", help="salmon or hisat2_featurecounts."),
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
    scaffold: bool = typer.Option(False, "--scaffold", help="Create templates only; the project remains incomplete until inputs are supplied."),
    yes: bool = typer.Option(False, "--yes", help="Create without an interactive confirmation."),
) -> None:
    """Create a reviewed RNA-seq project interactively or from explicit flags."""

    noninteractive = any(value is not None for value in (name, destination, species, input_type, fastq_samplesheet, counts, metadata, contrasts, layout, preprocessing, method, strandedness, reference_source, reference_root, reference_manifest, reference_fasta, reference_gtf, reference_transcript_fasta, reference_salmon_index, reference_hisat2_index, preset, design_type, condition_column, covariate, pairing_column)) or scaffold or yes
    try:
        if not noninteractive and not sys.stdin.isatty():
            raise ProjectCreationError(
                "rnaseq new needs an interactive terminal, or explicit non-interactive flags (use --scaffold for templates)."
            )
        if not noninteractive:
            name = typer.prompt("Project name")
            destination = Path(typer.prompt("Destination directory", default="."))
            species = _wizard_choice("Species", [("human", "Homo sapiens"), ("mouse", "Mus musculus")], default="human")
            input_type = _wizard_choice("Input type", [("fastq", "Reads requiring upstream processing"), ("raw_counts", "Integer gene-level raw-count matrix")], default="fastq")
            importing = typer.confirm("Import existing inputs now?", default=False)
            scaffold = not importing
            if input_type == "fastq":
                if importing:
                    fastq_samplesheet = Path(typer.prompt("FASTQ samplesheet path"))
                else:
                    layout = _wizard_choice("Sequencing layout", [("paired_end", "R1 and R2 per lane"), ("single_end", "R1 only")], default="paired_end")
                preprocessing = _wizard_choice("FASTQ preprocessing", [("raw", "Run adapter/quality trimming"), ("pretrimmed", "Keep supplied reads; skip trimming")], default="raw")
                method = _wizard_choice("Quantification backend", [("salmon", "nf-core/rnaseq pseudoalignment (default)"), ("hisat2_featurecounts", "HISAT2 alignment plus gene-level featureCounts")], default="salmon")
                if not importing:
                    strand_choices = [("auto", "infer with Salmon"), ("unstranded", "no stranded protocol"), ("forward", "forward stranded"), ("reverse", "reverse stranded")]
                    if method == "hisat2_featurecounts":
                        strand_choices = strand_choices[1:]
                    strandedness = _wizard_choice("Strandedness", strand_choices, default="auto" if method == "salmon" else "unstranded")
                reference_source = _wizard_choice("Reference", [("igenomes", "managed iGenomes convenience route (Salmon only)"), ("local", "checksum-bound managed local reference"), ("custom", "copy supported custom assets into the project")], default="igenomes" if method == "salmon" else "local")
                if reference_source == "local":
                    reference_root = Path(typer.prompt("Managed reference root"))
                    reference_manifest = typer.prompt("Reference manifest", default="reference_manifest.yaml")
                elif reference_source == "custom":
                    reference_fasta = Path(typer.prompt("Custom genome FASTA"))
                    reference_gtf = Path(typer.prompt("Custom annotation GTF"))
                    if method == "hisat2_featurecounts":
                        reference_hisat2_index = Path(typer.prompt("Prepared HISAT2 index directory (blank if not built)", default="")) if typer.confirm("Is a HISAT2 index already prepared?", default=False) else None
            else:
                if importing:
                    counts = Path(typer.prompt("Raw-count matrix path"))
            if importing:
                metadata_text = typer.prompt("Metadata CSV path (blank keeps a template)", default="")
                contrasts_text = typer.prompt("Contrasts CSV path (blank keeps a template)", default="")
                metadata, contrasts = (Path(metadata_text) if metadata_text else None), (Path(contrasts_text) if contrasts_text else None)
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
            required_fields = [condition_column, *selected_covariates]
            if normalized_design is DesignType.PAIRED:
                if pairing_column is None:
                    raise ProjectCreationError("Paired imported design requires --pairing-column.")
                required_fields.insert(0, pairing_column)
            missing_fields = [field for field in required_fields if field not in available_fields]
            if missing_fields:
                raise ProjectCreationError("Selected design field(s) are not in imported metadata: " + ", ".join(missing_fields))
            # Formula order is deterministic: pairing first, then covariates,
            # with the contrast factor last.
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
        review = {"Project": name, "Destination": destination, "Species": normalized_species.value, "Input": normalized_input.value, "Scope": normalized_preset.value, "Design": normalized_design.value, "Formula": formula or "template default", "Backend": normalized_method if normalized_input is InputType.FASTQ else "external raw counts", "Reference": reference.get("source"), "Mode": "scaffold (incomplete)" if scaffold else "import"}
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
                    scaffold=False, yes=False,
                )
            typer.echo("Project creation cancelled; no project was written.")
            return
        target = create_project(project_name=name, destination=destination, species=normalized_species, preset=normalized_preset, design_type=normalized_design, input_type=normalized_input, layout=normalized_layout, preprocessing=normalized_preprocessing, strandedness=normalized_strand, quantification_method=normalized_method, reference=reference, fastq_samplesheet=fastq_samplesheet, counts_file=counts, metadata_file=metadata, contrasts_file=contrasts, scaffold=scaffold, formula=formula)
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
