"""Command-line interface for the FASTQ and raw-count Milestone 1 control plane."""

from __future__ import annotations

from pathlib import Path

import click
import typer

from rnaseq.errors import ExecutionPreflightError, ProjectCreationError, UpstreamExecutionError
from rnaseq.execution import doctor_checks, load_run_states
from rnaseq.service import execute_service_run, prepare_service_run, sanitize_completed_delivery, validate_case_id
from rnaseq.models import DesignType, InputType, PIPELINE_VERSION, Preset, SequencingLayout, Species
from rnaseq.planner import generate_plan
from rnaseq.project import create_project
from rnaseq.references import (
    LocalReferenceError,
    ReferenceAdoptionError,
    ReferencePreparationError,
    adopt_local_salmon_index,
    prepare_local_reference,
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


@app.command("new")
def new_project() -> None:
    """Interactively create a FASTQ or raw-count RNA-seq project."""

    try:
        project_name = typer.prompt("Project name")
        destination = typer.prompt("Destination directory", default=".")
        species = typer.prompt(
            "Species",
            type=click.Choice([item.value for item in Species], case_sensitive=True),
        )
        input_type = typer.prompt(
            "Input type",
            type=click.Choice([item.value for item in InputType], case_sensitive=True),
            default=InputType.FASTQ.value,
        )
        layout = None
        if input_type == InputType.FASTQ.value:
            layout = typer.prompt(
                "Sequencing layout",
                type=click.Choice([item.value for item in SequencingLayout], case_sensitive=True),
                default=SequencingLayout.PAIRED_END.value,
            )
        preset = typer.prompt(
            "Analysis preset",
            type=click.Choice([item.value for item in Preset], case_sensitive=True),
        )
        design_type = typer.prompt(
            "Experimental design",
            type=click.Choice([item.value for item in DesignType], case_sensitive=True),
        )
        target = create_project(
            project_name=project_name,
            destination=destination,
            species=Species(species),
            preset=Preset(preset),
            design_type=DesignType(design_type),
            input_type=InputType(input_type),
            layout=SequencingLayout(layout) if layout else None,
        )
    except ProjectCreationError as exc:
        typer.echo(f"ERROR: {exc}", err=True)
        raise typer.Exit(code=1) from exc
    except (OSError, UnicodeError) as exc:
        typer.echo(f"SYSTEM ERROR: {exc}", err=True)
        raise typer.Exit(code=2) from exc

    typer.echo(f"Created project: {target}")
    if input_type == InputType.FASTQ.value:
        typer.echo(f"Place gzipped FASTQ files at: {target / 'input' / 'fastq'}")
    else:
        typer.echo(f"Place the raw count matrix at: {target / 'input' / 'counts.csv'}")
    typer.echo("Then complete metadata.csv and contrasts.csv before validation.")


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
                typer.echo(f"Pipeline: nf-core/rnaseq {report.config.upstream.pipeline_version}")
                typer.echo(f"Samples: {len(report.fastq.sample_ids)}")
            typer.echo(f"Profile: {profile}")
            typer.echo("An immutable case run, downstream workflow, report, and delivery package will be created.")
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
