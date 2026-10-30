"""Command-line interface for the FASTQ and raw-count Milestone 1 control plane."""

from __future__ import annotations

import csv
import json
import shlex
import sys
import time
from collections.abc import Iterator
from contextlib import contextmanager
from pathlib import Path
from typing import Any

import typer

from rnaseq.errors import ExecutionPreflightError, ProjectCreationError, UpstreamExecutionError
from rnaseq.execution import (
    detect_local_resource_capacity, doctor_checks, doctor_readiness,
    validate_local_execution_budget,
)
from rnaseq.service import RunTermination, execute_retry_service_run, execute_service_run, preflight_reuse_upstream, prepare_service_run, sanitize_completed_delivery, validate_case_id
from rnaseq.models import ANNOTATION_DATABASES, DEFAULT_L2_ENRICHMENT, GENE_ID_TYPES, DesignType, FastqPreprocessing, InputType, PIPELINE_VERSION, Preset, SequencingLayout, Species, default_annotation, describe_enrichment
from rnaseq.planner import generate_plan
from rnaseq.resource_policy import AUTO
from rnaseq.design_setup import (
    DesignRejected, DesignSetupError, apply_design, build_proposal, design_blockers, design_enrichment_selection,
    design_needs_gene_id_type, residual_df_blocker, resolve_design_enrichment,
)
from rnaseq.project import create_project, derive_design, resolve_enrichment_defaults
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
from rnaseq.validators import ValidationReport, contingency_table, design_dimensions, validate_project
from rnaseq.wizard import (
    DEFAULT_PAIR_ID, DESIGN_CHOICES, GENE_ID_TYPE_CHOICES, INPUT_TYPE_CHOICES, LAYOUT_CHOICES, METHOD_CHOICES,
    PREPROCESSING_CHOICES, REFERENCE_SOURCE_CHOICES, SPECIES_CHOICES, default_reference_source, default_strandedness,
    managed_reference_settings, preset_choices, reference_options, resource_suggestion, species_value,
    strandedness_choices,
)

app = typer.Typer(
    name="rnaseq",
    help="Create, validate, and plan versioned FASTQ or raw-count bulk RNA-seq projects.",
    no_args_is_help=True,
)
reference_app = typer.Typer(help="Prepare and inspect managed local references.", no_args_is_help=True)
app.add_typer(reference_app, name="reference")


def _show_version(value: bool) -> None:
    if value:
        typer.echo(PIPELINE_VERSION)
        raise typer.Exit()


@app.callback()
def main(
    version: bool = typer.Option(
        False,
        "--version",
        help="Show the nf-rna version and exit.",
        is_eager=True,
        callback=_show_version,
    ),
) -> None:
    """nf-rna command-line control plane."""


def _render_batch_contingency(report: ValidationReport) -> list[str]:
    """Show batch x condition sample counts when the design adjusts for a batch."""

    if report.metadata is None or "batch" not in report.formula_variables or not report.metadata.rows:
        return []
    factor = next((name for name in report.groups if name != "batch"), None)
    if factor is None or factor not in report.metadata.columns or "batch" not in report.metadata.columns:
        return []
    table = contingency_table(report.metadata, "batch", factor)
    columns = list(report.groups[factor])
    lines = ["Batch x condition (samples)", "---------------------------", f"{'batch':<20}" + "".join(f"{column:<14}" for column in columns)]
    for batch_level, counts in table.items():
        lines.append(f"{batch_level:<20}" + "".join(f"{counts.get(column, 0):<14}" for column in columns))
    return lines


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
        lines.extend(["Design", "------", f"Type: {config.design.type.value}"])
        if config.design.pair_id is not None:
            lines.append(f"Pairing variable: {config.design.pair_id}")
        lines.append(f"DESeq2 formula: {config.design.formula}")
        if config.design.condition_order:
            lines.append(f"Condition order: {', '.join(config.design.condition_order)}")
        if config.project.preset is not Preset.QC and report.contrasts is not None and not report.contrasts.contrasts:
            lines.append(
                "Contrasts: not yet defined (required for L2: run `rnaseq design PROJECT` to configure them)"
                if config.project.preset is Preset.L2
                else "Contrasts: not yet defined (run `rnaseq design PROJECT` after reviewing L1)"
            )
    lines.extend(_render_batch_contingency(report))
    if report.pairing is not None:
        lines.extend(
            [
                f"Total samples: {report.pairing.total_samples}",
                f"Unique pair IDs: {report.pairing.unique_pair_ids}",
            ]
        )
        for summary in report.pairing.contrasts:
            lines.extend(
                [
                    f"Comparison: {summary.numerator} vs {summary.denominator} ({summary.contrast_id})",
                    f"Complete pairs: {summary.complete_pairs}",
                    f"Samples in comparison: {summary.samples_in_complete_pairs}",
                    f"Incomplete pairs: {len(summary.incomplete_pairs)}",
                    "Duplicated/invalid pair structures: "
                    f"{len(set(summary.duplicate_numerator_pairs) | set(summary.duplicate_denominator_pairs))}",
                ]
            )
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
    rebuild_salmon: bool = typer.Option(
        False, "--rebuild-salmon",
        help="Replace an existing Salmon declaration with a GTF-derived index; the old manifest is archived and the old index kept.",
    ),
) -> None:
    """Build a checksum-bound, host-native decoy-aware Salmon index from a GTF-derived transcriptome."""

    try:
        reference = prepare_local_reference(reference_root, threads=threads, rebuild_salmon=rebuild_salmon)
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


def _wizard_prompt_completer(choices: list[str]) -> Any | None:
    """Build the terminal completer used by the live wizard prompt."""

    try:
        from prompt_toolkit.completion import WordCompleter
    except ImportError:
        return None
    return WordCompleter(choices, sentence=True, match_middle=False)


def _prompt_toolkit_choice_prompt(choices: list[str], default: str | None) -> str | None:
    """Read one value with completion active while the user edits the line."""

    # Click's test/non-TTY streams must keep the established plain prompt path.
    if not sys.stdin.isatty():
        return None
    completer = _wizard_prompt_completer(choices)
    if completer is None:
        return None
    try:
        from prompt_toolkit import prompt
    except ImportError:
        return None
    selected = prompt(
        f"Choose [{default}]: " if default is not None else "Choose: ",
        completer=completer,
        complete_while_typing=False,
        default="",
    )
    return default if selected == "" and default is not None else selected


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
        selected = _prompt_toolkit_choice_prompt(allowed, default) if _is_interactive_terminal() else None
        if selected is None:
            with _wizard_tab_completion(allowed):
                selected = typer.prompt("Choose", default=default)
        if selected in allowed:
            return selected
        typer.echo(f"Invalid choice {selected!r}.")
        typer.echo("Please choose one of: " + ", ".join(allowed) + ".")


def _wizard_positive_integer(label: str, *, default: int | str, allow_auto: bool = False) -> int | str:
    """Prompt one positive integer (or ``auto``) without exposing a conversion traceback."""

    while True:
        raw = typer.prompt(label, default=str(default))
        if allow_auto and raw.strip().lower() == AUTO:
            return AUTO
        try:
            value = int(raw)
        except (TypeError, ValueError):
            value = 0
        if value > 0:
            return value
        typer.echo(f"Invalid value {raw!r}. Please enter a positive integer" + (" or 'auto'." if allow_auto else "."))


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


def _reference_source(explicit: str | None, method: str, *, local_flags: list[str], custom_flags: list[str]) -> str:
    """Resolve the FASTQ reference source without ever ignoring a reference flag.

    Without --reference-source, local or custom flags select that source; with
    none, the documented default applies (iGenomes for Salmon, local for HISAT2).
    The machine-local registry is never consulted here: only the wizard offers it.
    """

    if explicit is not None and explicit not in {"igenomes", "local", "custom"}:
        raise ProjectCreationError("--reference-source must be igenomes, local, or custom.")
    if explicit is None and local_flags and custom_flags:
        raise ProjectCreationError(f"Conflicting reference flags: {', '.join(local_flags)} select a local reference, {', '.join(custom_flags)} a custom one; set --reference-source.")
    source = explicit or ("local" if local_flags else "custom" if custom_flags else "igenomes" if method == "salmon" else "local")
    ignored = custom_flags if source == "local" else local_flags if source == "custom" else local_flags + custom_flags
    if ignored:
        raise ProjectCreationError(f"{', '.join(ignored)} cannot be used with --reference-source {source}.")
    return source


def _reference_review(reference: dict[str, object], registered: LocalReference | None, *, defaulted: bool) -> str:
    """State the selected reference and how it was chosen in the creation review."""

    if registered is not None:
        return f"{_managed_reference_identity(registered)} (registered local reference, root {registered.root})"
    if reference["source"] == "igenomes":
        how = "Salmon default, no --reference-source given; registered local references are used only with --reference-source local --reference-root" if defaulted else "selected"
        return f"igenomes {reference['genome']} ({how})"
    if reference["source"] == "local":
        return f"local {reference['root']} (manifest {reference['manifest']})"
    return f"custom {reference['fasta']} + {reference['gtf']}"


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


def _enrichment_selection(enabled: bool | None) -> tuple[str, ...] | None:
    """Map --enrichment/--no-enrichment onto the project-level selection."""

    return None if enabled is None else (DEFAULT_L2_ENRICHMENT if enabled else ())


def _annotation_review(species: Species, gene_id_type: str | None) -> str:
    if gene_id_type is None:
        return "not configured"
    database = ANNOTATION_DATABASES.get(species.value, "unsupported organism")
    return f"{species.value}; {gene_id_type} -> ENTREZID via {database}"


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
    design_type: str | None = typer.Option(None, "--design-type", help="two_group, paired_two_group, or multi_group."),
    condition_column: str | None = typer.Option(None, "--condition-column", help="Advanced: imported metadata factor used for contrasts and the design formula (requires --metadata; default 'condition')."),
    covariate: list[str] | None = typer.Option(None, "--covariate", help="Advanced: additional imported categorical metadata field (requires --metadata); repeat as needed."),
    continuous_covariate: list[str] | None = typer.Option(None, "--continuous-covariate", help="Advanced: imported numeric adjustment variable (requires --metadata); repeat as needed."),
    pair_id: str | None = typer.Option(None, "--pair-id", "--pairing-column", help="Biological pair column (paired_two_group only); with --metadata it must exist there."),
    batch: bool | None = typer.Option(None, "--batch/--no-batch", help="Declare a known technical batch: the template gets a 'batch' column and the model adjusts for it. With --metadata, equivalent to --covariate batch."),
    gene_id_type: str | None = typer.Option(None, "--gene-id-type", help="Gene identifier type of the count matrix for enrichment annotation: ENSEMBL, ENTREZID, or SYMBOL. Required for L2 unless --no-enrichment."),
    enrichment: bool | None = typer.Option(None, "--enrichment/--no-enrichment", help="L2 only. Default: GO ORA, KEGG ORA, and preranked GO + KEGG GSEA (analysis.enrichment lists go, kegg, gsea). --no-enrichment records an explicit opt-out."),
    execution_profile: str | None = typer.Option(None, "--execution-profile", help="Only local is supported."),
    cpus: int | None = typer.Option(None, "--cpus", help="Total local Nextflow CPU ceiling, not per-task CPUs."),
    memory_gb: int | None = typer.Option(None, "--memory-gb", help="Total local Nextflow memory ceiling in GiB, not per-task memory."),
    scaffold: bool = typer.Option(False, "--scaffold", help="Create templates only; the project remains incomplete until inputs are supplied."),
    yes: bool = typer.Option(False, "--yes", help="Create without an interactive confirmation."),
) -> None:
    """Create a reviewed RNA-seq project interactively or from explicit flags."""

    noninteractive = any(value is not None for value in (name, destination, species, input_type, fastq_samplesheet, counts, metadata, contrasts, layout, preprocessing, method, strandedness, reference_source, reference_root, reference_manifest, reference_fasta, reference_gtf, reference_transcript_fasta, reference_salmon_index, reference_hisat2_index, preset, design_type, condition_column, covariate, continuous_covariate, pair_id, batch, gene_id_type, enrichment, execution_profile, cpus, memory_gb)) or scaffold or yes
    selected_managed_reference: LocalReference | None = None
    try:
        if not noninteractive and not _is_interactive_terminal():
            raise ProjectCreationError(
                "rnaseq new needs an interactive terminal, or explicit non-interactive flags (use --scaffold for templates)."
            )
        if not noninteractive:
            name = typer.prompt("Project name")
            destination = Path(typer.prompt("Destination directory", default="."))
            species = _wizard_choice("Species", SPECIES_CHOICES, default="human")
            input_type = _wizard_choice("Input type", INPUT_TYPE_CHOICES, default="fastq")
            # Ordinary interactive creation is deliberately scaffold-first.
            # Importing is an explicit, reproducible flags-only operation.
            scaffold = True
            if input_type == "fastq":
                layout = _wizard_choice("Sequencing layout", LAYOUT_CHOICES, default="paired_end")
                preprocessing = _wizard_choice("FASTQ preprocessing", PREPROCESSING_CHOICES, default="raw")
                method = _wizard_choice("Quantification backend", METHOD_CHOICES, default="salmon")
                strandedness = _wizard_choice("Strandedness", strandedness_choices(method), default=default_strandedness(method))
                selected_managed_reference = _registered_reference_choice(species_value(species), method)
                if selected_managed_reference is not None:
                    reference_source = "local"
                    reference_root, reference_manifest = managed_reference_settings(selected_managed_reference)
                else:
                    reference_source = _wizard_choice("Reference", REFERENCE_SOURCE_CHOICES, default=default_reference_source(method))
                if reference_source == "local" and selected_managed_reference is None:
                    reference_root = Path(typer.prompt("Managed reference root"))
                    reference_manifest = typer.prompt("Reference manifest", default="reference_manifest.yaml")
                elif reference_source == "custom":
                    reference_fasta = Path(typer.prompt("Custom genome FASTA"))
                    reference_gtf = Path(typer.prompt("Custom annotation GTF"))
                    if method == "hisat2_featurecounts":
                        reference_hisat2_index = Path(typer.prompt("Prepared HISAT2 index directory (blank if not built)", default="")) if typer.confirm("Is a HISAT2 index already prepared?", default=False) else None
            preset = _wizard_choice("Analysis scope", preset_choices(input_type), default="L1")
            design_type = "two_group" if preset == "qc" else _wizard_choice("Experimental design", DESIGN_CHOICES, default="two_group")
            if design_type == "paired_two_group":
                pair_id = typer.prompt("Pairing metadata column", default=DEFAULT_PAIR_ID)
            if preset != "qc":
                # Only what cannot be read from the completed metadata is asked
                # here: condition levels, batch levels and contrasts are
                # discovered later by `rnaseq design`.
                batch = typer.confirm("Is there a known technical batch (adds a 'batch' column and adjusts the model for it)?", default=False)
                if preset == "L2":
                    typer.echo(f"Functional enrichment (default for L2): {describe_enrichment(DEFAULT_L2_ENRICHMENT)}")
                    enrichment = typer.confirm("Run these enrichment analyses?", default=True)
                if preset == "L1" or enrichment:
                    # Asked, never assumed: the annotation contract maps these
                    # identifiers to Entrez IDs for GO/KEGG.
                    organism = species_value(species).value
                    typer.echo(f"Annotation database for {organism}: {ANNOTATION_DATABASES[organism]}")
                    gene_id_type = _wizard_choice("Gene identifier type in the count matrix", GENE_ID_TYPE_CHOICES, default="ENSEMBL")
            suggestion = resource_suggestion(detect_local_resource_capacity())
            typer.echo("\nLocal execution resources\n-------------------------")
            typer.echo(f"Detected: {suggestion.logical_cpus or 'unavailable'} logical CPUs / {suggestion.memory_gib or 'unavailable'} GiB memory")
            typer.echo(f"auto on this machine: {suggestion.auto_cpus} CPUs / {suggestion.auto_memory_gib} GiB memory (resolved again on each run's machine)")
            typer.echo("Enter 'auto' to size limits from the running machine, or a number to fix a deliberate limit.")
            cpus = _wizard_positive_integer("CPU limit", default=AUTO, allow_auto=True)
            memory_gb = _wizard_positive_integer("Memory limit in GiB", default=AUTO, allow_auto=True)
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
        # Keep the interim CLI spelling readable while all newly written
        # projects use the M5A contract name.
        normalized_design = DesignType("paired_two_group" if design_type == "paired" else design_type)
        if normalized_design is DesignType.PAIRED_TWO_GROUP and (pair_id is None or not pair_id.strip()):
            raise ProjectCreationError("paired_two_group requires an explicit --pair-id metadata column.")
        formula: str | None = None
        design_variables: dict[str, str] | None = None
        template_batch = False
        # No design option may be silently ignored: each one either takes
        # effect or is rejected with the reason.
        design_flags = [
            flag for flag, value in (
                ("--condition-column", condition_column), ("--covariate", covariate),
                ("--continuous-covariate", continuous_covariate),
            ) if value
        ]
        if normalized_preset is Preset.QC:
            inapplicable = design_flags + (["--pair-id"] if pair_id else []) + (["--batch"] if batch else [])
            if inapplicable:
                raise ProjectCreationError(
                    f"{', '.join(inapplicable)} {'does' if len(inapplicable) == 1 else 'do'} not apply to the qc preset, which has no statistical design."
                )
        if pair_id is not None and normalized_design is not DesignType.PAIRED_TWO_GROUP:
            raise ProjectCreationError("--pair-id applies only to --design-type paired_two_group.")
        if metadata is None and design_flags and normalized_preset is not Preset.QC:
            hint = " For a technical batch use --batch (adds a 'batch' column to the template)." if "--covariate" in design_flags else ""
            raise ProjectCreationError(
                f"{', '.join(design_flags)} {'requires' if len(design_flags) == 1 else 'require'} --metadata "
                f"(they describe an imported metadata file); the generated template always uses the standard 'condition' column.{hint}"
            )
        if normalized_preset is not Preset.QC and metadata is None:
            template_batch = batch is True
        if normalized_preset is not Preset.QC and metadata is not None:
            fields_by_level, available_fields = _metadata_fields(metadata)
            selected_condition = condition_column
            if selected_condition is None:
                if "condition" not in available_fields:
                    raise ProjectCreationError(
                        "Imported metadata has no 'condition' column; name the contrast factor with --condition-column. "
                        "Available fields: " + ", ".join(available_fields)
                    )
                selected_condition = "condition"
            categorical_covariates = list(covariate or [])
            continuous_covariates = list(continuous_covariate or [])
            if batch is True and "batch" not in categorical_covariates:
                categorical_covariates.append("batch")
            if batch is False and "batch" in categorical_covariates:
                raise ProjectCreationError("--no-batch conflicts with --covariate batch.")
            selected_covariates = [*categorical_covariates, *continuous_covariates]
            if selected_condition in selected_covariates:
                raise ProjectCreationError(f"{selected_condition!r} is the condition column and cannot also be a covariate.")
            if selected_condition in continuous_covariates:
                raise ProjectCreationError("--condition-column is a categorical contrast factor and cannot be continuous.")
            if pair_id is not None and pair_id in selected_covariates:
                raise ProjectCreationError(f"--pair-id {pair_id!r} is already part of the design and cannot also be a covariate.")
            if pair_id is not None and pair_id == selected_condition:
                raise ProjectCreationError("--pair-id cannot be the condition column.")
            if normalized_design is DesignType.PAIRED_TWO_GROUP and pair_id is None:
                raise ProjectCreationError("paired_two_group imported design requires --pair-id.")
            required_fields = [*selected_covariates, selected_condition]
            if pair_id is not None:
                required_fields.insert(0, pair_id)
            missing_fields = [field for field in required_fields if field not in available_fields]
            if missing_fields:
                raise ProjectCreationError(
                    "Selected design field(s) are not in imported metadata: " + ", ".join(missing_fields)
                    + ". Available fields: " + ", ".join(available_fields)
                )
            # Formula order is deterministic: pairing first, then covariates,
            # with the contrast factor last. Pairing identity itself is stored
            # explicitly and never inferred from this ordering.
            formula = "~ " + " + ".join(dict.fromkeys(required_fields))
            design_variables = {field: "categorical" for field in required_fields}
            design_variables.update({field: "continuous" for field in continuous_covariates})
            typer.echo(f"Design from imported metadata: {formula}")
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
            source = _reference_source(
                reference_source, normalized_method,
                local_flags=[flag for flag, value in (("--reference-root", reference_root), ("--reference-manifest", reference_manifest)) if value is not None],
                custom_flags=[flag for flag, value in (("--reference-fasta", reference_fasta), ("--reference-gtf", reference_gtf), ("--reference-transcript-fasta", reference_transcript_fasta), ("--reference-salmon-index", reference_salmon_index), ("--reference-hisat2-index", reference_hisat2_index)) if value is not None],
            )
            reference = reference_options(source=source, species=normalized_species, method=normalized_method, local_root=reference_root, local_manifest=reference_manifest, fasta=reference_fasta, gtf=reference_gtf, transcript_fasta=reference_transcript_fasta, salmon_index=reference_salmon_index, hisat2_index=reference_hisat2_index)
        if execution_profile not in {None, "local"}:
            raise ProjectCreationError("Only --execution-profile local is supported.")
        selected_cpus, selected_memory = cpus or AUTO, memory_gb or AUTO
        validate_local_execution_budget(selected_cpus, selected_memory, detect_local_resource_capacity())
        execution = {"profile": "local", "max_cpus": selected_cpus, "max_memory_gb": selected_memory}
        normalized_gene_id_type = gene_id_type.upper() if gene_id_type else None
        selected_enrichment, _annotation = resolve_enrichment_defaults(
            normalized_preset, normalized_species, _enrichment_selection(enrichment), normalized_gene_id_type
        )
        review_reference = (
            "not used (external raw counts)"
            if normalized_input is InputType.RAW_COUNTS
            else _reference_review(reference, selected_managed_reference, defaulted=reference_source is None)
        )
        review = {"Project": name, "Destination": destination, "Species": normalized_species.value, "Input": normalized_input.value, "Input handling": "scaffold — add data after project creation" if scaffold else "import supplied inputs", "Scope": normalized_preset.value, "Design": "not applicable for QC" if normalized_preset is Preset.QC else normalized_design.value, **({} if normalized_preset is Preset.QC else {"Technical batch": "yes (metadata template gets a 'batch' column)" if template_batch else ("yes (in model)" if formula is not None and "batch" in formula.replace("~", " ").replace("+", " ").split() else "no"), "Provisional model": formula or derive_design(normalized_design, pair_id, template_batch)[0], "Contrasts": "imported" if contrasts is not None else ("set later by `rnaseq design` (required before L2 can run)" if normalized_preset is Preset.L2 else "set later by `rnaseq design`"), "Enrichment": describe_enrichment(selected_enrichment) if normalized_preset is Preset.L2 else "selected when L2 is configured by `rnaseq design`", "Annotation": _annotation_review(normalized_species, normalized_gene_id_type)}), "Backend": normalized_method if normalized_input is InputType.FASTQ else "external raw counts", "Reference": review_reference, "Execution": f"local, {selected_cpus} CPUs / {selected_memory} GiB" if AUTO not in (selected_cpus, selected_memory) else f"local, CPUs={selected_cpus} / memory GiB={selected_memory} (auto = sized from the running machine)"}
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
                    design_type=None, condition_column=None, covariate=None, continuous_covariate=None, pair_id=None, batch=None,
                    gene_id_type=None, enrichment=None, execution_profile=None, cpus=None, memory_gb=None,
                    scaffold=False, yes=False,
                )
            typer.echo("Project creation cancelled; no project was written.")
            return
        target = create_project(project_name=name, destination=destination, species=normalized_species, preset=normalized_preset, design_type=normalized_design, input_type=normalized_input, layout=normalized_layout, preprocessing=normalized_preprocessing, strandedness=normalized_strand, quantification_method=normalized_method, reference=reference, fastq_samplesheet=fastq_samplesheet, counts_file=counts, metadata_file=metadata, contrasts_file=contrasts, scaffold=scaffold, formula=formula, pair_id=pair_id, design_variables=design_variables, execution=execution, batch=template_batch, enrichment=_enrichment_selection(enrichment), gene_id_type=normalized_gene_id_type)
    except (ValueError, ProjectCreationError) as exc:
        typer.echo(f"ERROR: {exc}", err=True)
        raise typer.Exit(code=1) from exc
    except (OSError, UnicodeError) as exc:
        typer.echo(f"SYSTEM ERROR: {exc}", err=True)
        raise typer.Exit(code=2) from exc

    typer.echo(f"Created project: {target}")
    project_arg = shlex.quote(str(target))
    if scaffold:
        typer.echo(f"This is an incomplete scaffold. Add valid inputs and complete metadata.csv, then run: rnaseq validate {project_arg}")
        if normalized_preset is not Preset.QC:
            header = ",".join(("sample_id", *derive_design(normalized_design, pair_id, template_batch)[2])) if metadata is None else "(imported)"
            for line in _scaffold_next_steps(target, normalized_preset, normalized_input, header):
                typer.echo(line)
    elif normalized_input is InputType.FASTQ:
        typer.echo(f"Imported FASTQs were copied into input/fastq; run rnaseq validate {project_arg} to check readiness.")
    else:
        typer.echo(f"Imported raw counts and supplied metadata/contrasts were preserved; run rnaseq validate {project_arg}.")


def _scaffold_next_steps(target: Path, preset: Preset, input_type: InputType, header: str) -> list[str]:
    """The ordered commands that take a new L1 or L2 scaffold to its first run."""

    project = shlex.quote(str(target))
    inputs = (
        f"Copy the FASTQ files into {target / 'input' / 'fastq'} and fill {target / 'metadata.csv'} (columns: {header}).\n"
        "     Each sample_id must equal the FASTQ sample ID: the file name before _Lnnn/_R1/_R2, "
        "e.g. Ctrl1_S1_L001_R1_001.fastq.gz -> Ctrl1_S1."
        if input_type is InputType.FASTQ
        else f"Put the count matrix in {target / 'input' / 'counts.csv'} and fill {target / 'metadata.csv'} (columns: {header})."
    )
    if preset is Preset.L2:
        return [
            "Next steps (L2: differential expression and enrichment; the L2 run also performs the L1 QC):",
            f"  1. {inputs}",
            f"  2. Configure the contrasts: rnaseq design {project}",
            "     (writes contrasts.csv and refreshes the planning artifacts; L2 cannot be validated or run without it)",
            f"  3. Confirm readiness: rnaseq validate {project}  and  rnaseq doctor {project}",
            f"  4. Run L2: rnaseq run {project} --case-id CASE-ID",
            "     One L2 run produces the L1 QC, differential expression and enrichment; no separate L1 run is needed.",
        ]
    reuse = " --reuse-upstream CASE-ID/RUN-ID (of the successful L1 run)" if input_type is InputType.FASTQ else ""
    return [
        "Next steps (L1: QC and exploratory analysis):",
        f"  1. {inputs}",
        f"  2. Validate: rnaseq validate {project}",
        f"  3. Plan and check the machine: rnaseq plan {project}  and  rnaseq doctor {project}",
        f"  4. Run L1: rnaseq run {project} --case-id CASE-ID",
        f"  5. Review the L1 QC in {target / 'runs'}/CASE-ID/RUN-ID/delivery/",
        f"  6. Configure L2 (differential expression): rnaseq design {project}",
        f"     then run L2 with a new case ID: rnaseq run {project} --case-id NEW-CASE-ID{reuse}",
    ]


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


def _latest_successful_run(project: Path) -> str | None:
    """Return ``CASE-ID/RUN-ID`` of the newest successful case run, if any."""

    candidates: list[tuple[str, str]] = []
    for state in (project / "runs").glob("*/*/run_state.json"):
        try:
            if json.loads(state.read_text(encoding="utf-8")).get("status") == "SUCCESS":
                candidates.append((state.parent.name, state.parent.parent.name))
        except (OSError, ValueError):
            continue
    if not candidates:
        return None
    run_id, case_id = max(candidates)
    return f"{case_id}/{run_id}"


def _show_design_summary(report: ValidationReport, proposal: Any) -> None:
    config = report.config
    assert config is not None and report.metadata is not None
    typer.echo(f"Design setup for {config.project.id}")
    typer.echo("-" * 24)
    typer.echo(f"Design type: {config.design.type.value}")
    origin = "derived from the design choices" if proposal.standard_formula else "custom formula kept as written"
    typer.echo(f"Model: {proposal.formula}  ({origin})")
    ordering = "frozen order" if config.design.condition_order else "order of first appearance in metadata.csv"
    typer.echo(f"{proposal.factor.capitalize()} levels ({ordering}): " + ", ".join(f"{level} (n={proposal.counts[level]})" for level in proposal.levels))
    for line in _render_batch_contingency(report):
        typer.echo(line)
    dimensions = design_dimensions(report)
    if dimensions is not None:
        samples, coefficients = dimensions
        typer.echo(f"Residual degrees of freedom: {samples} samples - {coefficients} model coefficients = {samples - coefficients}")
    for issue in report.warnings:
        typer.echo(f"WARNING: {issue.message}")


def _show_contrasts(proposal: Any) -> None:
    for contrast in proposal.contrasts:
        typer.echo(f"  {contrast.contrast_id}: {contrast.numerator} vs {contrast.denominator}")
        typer.echo(f"      log2FC > 0 = higher in {contrast.numerator}; log2FC < 0 = higher in {contrast.denominator}")


def _design_enrichment(
    config: Any, enabled: bool | None, gene_id_type: str | None, *, interactive: bool,
) -> tuple[tuple[str, ...] | None, dict[str, str] | None]:
    """Decide the enrichment selection and annotation written by `rnaseq design`.

    Interactively, a missing gene identifier type is asked for instead of
    failing; the decision itself is :func:`resolve_design_enrichment`.
    """

    selection = design_enrichment_selection(config, enabled)
    if interactive and gene_id_type is None and design_needs_gene_id_type(config, selection):
        typer.echo(f"Functional enrichment (default for L2): {describe_enrichment(DEFAULT_L2_ENRICHMENT)}")
        typer.echo(f"Annotation database for {config.organism.species.value}: {ANNOTATION_DATABASES[config.organism.species.value]}")
        gene_id_type = _wizard_choice("Gene identifier type in the count matrix", GENE_ID_TYPE_CHOICES, default="ENSEMBL")
    return resolve_design_enrichment(config, enabled, gene_id_type)


@app.command("design")
def design_command(
    project_dir: Path,
    condition_order: str | None = typer.Option(None, "--condition-order", help="Comma-separated condition levels, reference (earliest) first. Default: order of first appearance in metadata.csv."),
    swap: bool = typer.Option(False, "--swap", help="Two-level designs: reverse the proposed contrast direction."),
    force: bool = typer.Option(False, "--force", help="Replace an existing non-empty contrasts.csv."),
    yes: bool = typer.Option(False, "--yes", help="Accept the proposal without prompting."),
    gene_id_type: str | None = typer.Option(None, "--gene-id-type", help="Gene identifier type of the count matrix (ENSEMBL, ENTREZID, or SYMBOL); needed when the project has no annotation contract yet."),
    enrichment: bool | None = typer.Option(None, "--enrichment/--no-enrichment", help="Default when moving to L2: GO ORA, KEGG ORA, and preranked GO + KEGG GSEA. --no-enrichment records an explicit opt-out; an existing L2 selection is otherwise kept."),
) -> None:
    """Propose and confirm L2 contrasts from the completed metadata, then prepare the project for L2."""

    try:
        report = validate_project(project_dir)
        config = report.config
        if config is None:
            typer.echo(render_validation_report(report), err=True)
            raise typer.Exit(code=1)
        if config.project.preset is Preset.QC:
            typer.echo("ERROR: This project's scope is QC-only; there is no statistical design to configure.", err=True)
            raise typer.Exit(code=1)
        blocking = design_blockers(report)
        if blocking:
            typer.echo("ERROR: The metadata/design is not valid yet; fix the following and re-run `rnaseq design`:", err=True)
            for issue in blocking:
                typer.echo(f"  - {issue.message}", err=True)
            raise typer.Exit(code=1)
        existing = report.contrasts.contrasts if report.contrasts is not None else ()
        if existing and not force:
            typer.echo(f"ERROR: contrasts.csv already defines {len(existing)} contrast(s); it was not changed:", err=True)
            for contrast in existing:
                typer.echo(f"  {contrast.contrast_id}: {contrast.numerator} vs {contrast.denominator}", err=True)
            typer.echo("Edit it directly (advanced), or re-run with --force to replace it with the proposal.", err=True)
            raise typer.Exit(code=1)
        if not yes and not _is_interactive_terminal():
            raise ProjectCreationError("rnaseq design needs an interactive terminal, or --yes to accept the proposal.")
        order = tuple(level.strip() for level in condition_order.split(",")) if condition_order else None
        proposal = build_proposal(report, condition_order=order, swap=swap)
        selected_enrichment, new_annotation = _design_enrichment(config, enrichment, gene_id_type, interactive=not yes)
        _show_design_summary(report, proposal)
        no_df = residual_df_blocker(report)
        if no_df is not None:
            raise DesignSetupError(no_df + " Add replicates or simplify the design before configuring L2.")
        if existing:
            typer.echo("These existing contrasts will be replaced:")
            for contrast in existing:
                typer.echo(f"  {contrast.contrast_id}: {contrast.numerator} vs {contrast.denominator}")
        two_level = len(proposal.levels) == 2
        while True:
            typer.echo("Proposed contrast:" if two_level else f"Proposed contrasts ({len(proposal.contrasts)}, every pair; later level vs earlier level):")
            _show_contrasts(proposal)
            if yes:
                break
            if two_level:
                answer = typer.prompt("Accept? (y = write it, swap = reverse the direction, n = cancel)", default="y").strip().lower()
                if answer in {"y", "yes"}:
                    break
                if answer in {"s", "swap"}:
                    proposal = build_proposal(report, condition_order=proposal.levels, swap=not proposal.swapped)
                    continue
                if answer in {"n", "no"}:
                    typer.echo("No changes were written.")
                    return
                typer.echo("Please answer y, n, or swap.")
            else:
                if typer.confirm(f"Write these {len(proposal.contrasts)} contrasts and prepare L2?", default=True):
                    break
                typer.echo("No changes were written. Use --condition-order to change the order, or edit contrasts.csv.")
                return
        after = apply_design(report, proposal, enrichment=selected_enrichment, annotation=new_annotation)
    except typer.Exit:
        raise
    except DesignRejected as exc:
        typer.echo(f"ERROR: {exc.args[0].removesuffix('.')}:", err=True)
        for issue in exc.issues:
            typer.echo(f"  - {issue.message}", err=True)
        raise typer.Exit(code=1) from exc
    except (DesignSetupError, ProjectCreationError) as exc:
        typer.echo(f"ERROR: {exc}", err=True)
        raise typer.Exit(code=1) from exc
    except (OSError, UnicodeError) as exc:
        typer.echo(f"SYSTEM ERROR: {exc}", err=True)
        raise typer.Exit(code=2) from exc
    typer.echo(f"Wrote {len(proposal.contrasts)} contrast(s) to contrasts.csv, froze design.condition_order ({', '.join(proposal.levels)}), and set the scope to L2.")
    final_enrichment = after.config.analysis.enrichment if after.config and after.config.analysis else ()
    typer.echo(f"Functional enrichment: {describe_enrichment(tuple(final_enrichment))}")
    typer.echo("L2 validation: " + after.state)
    typer.echo("Planning artifacts refreshed for the new configuration (equivalent to `rnaseq plan`).")
    # Only FASTQ projects have an upstream stage whose results a later run can reuse.
    reuse = _latest_successful_run(Path(project_dir)) if config.input.type is InputType.FASTQ else None
    hint = f" --reuse-upstream {reuse}" if reuse else ""
    typer.echo(f"Next: rnaseq run {project_dir} --case-id CASE-ID{hint}")
    if reuse:
        typer.echo("(--reuse-upstream is accepted only when the frozen upstream inputs match.)")


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
    reuse_upstream: str | None = typer.Option(None, "--reuse-upstream", help="FASTQ projects: reuse the upstream results of a successful run, given as CASE-ID/RUN-ID (e.g. L1-001/20261029-101010+0800; see rnaseq status PROJECT --all)."),
    yes: bool = typer.Option(False, "--yes", help="Authorize execution without an interactive prompt."),
) -> None:
    """Execute one immutable client case through upstream and downstream workflows."""

    try:
        report = validate_project(project_dir)
        typer.echo(render_validation_report(report))
        validate_case_id(case_id)
        if reuse_upstream and report.config is not None and report.is_valid:
            # Fail on an unusable reuse source before the confirmation prompt
            # and before any run directory exists.
            preflight_reuse_upstream(report, reuse_upstream)
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
    except (KeyboardInterrupt, RunTermination) as exc:
        _report_interruption(exc)
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


@app.command("retry")
def retry_command(
    project_dir: Path,
    retry_of: str = typer.Option(..., "--retry-of", help="Required FAILED or INTERRUPTED CASE-ID/RUN-ID source run."),
    nextflow_resume: bool = typer.Option(False, "--nextflow-resume", help="Opt in to Nextflow -resume for this new retry attempt."),
    yes: bool = typer.Option(False, "--yes", help="Authorize retry without an interactive prompt."),
) -> None:
    """Create a new immutable attempt from one FAILED or INTERRUPTED run's frozen contract."""

    try:
        if not yes:
            typer.echo(f"Retry source: {retry_of}")
            typer.echo("A new immutable run will use the source run's frozen scientific contract; the original run will not change.")
            if nextflow_resume:
                typer.echo("Nextflow cache reuse is requested for this new retry attempt; it does not define retry identity.")
            if not typer.confirm("Proceed?", default=False):
                typer.echo("Retry cancelled; no run directory was created.")
                raise typer.Exit(code=0)
        result = execute_retry_service_run(project_dir, retry_of=retry_of, nextflow_resume=nextflow_resume)
    except typer.Exit:
        raise
    except (KeyboardInterrupt, RunTermination) as exc:
        _report_interruption(exc)
    except (ExecutionPreflightError, UpstreamExecutionError) as exc:
        typer.echo(f"ERROR: {exc}", err=True)
        raise typer.Exit(code=1) from exc
    except (OSError, UnicodeError) as exc:
        typer.echo(f"SYSTEM ERROR: {exc}", err=True)
        raise typer.Exit(code=2) from exc
    typer.echo(f"Case retry: SUCCESS ({result.case_id}/{result.run_id})")
    typer.echo(f"Run directory: {result.run_dir}")
    typer.echo(f"Delivery package: {result.run_dir / 'delivery'}")


def _report_interruption(exc: BaseException) -> None:
    name = str(exc) if isinstance(exc, RunTermination) and str(exc) else "SIGINT"
    typer.echo(f"INTERRUPTED: rnaseq received {name}. A run that had started was stopped (Nextflow included) and recorded as INTERRUPTED; see 'rnaseq status'.", err=True)
    raise typer.Exit(code=130 if name == "SIGINT" else 143)


@app.command("status")
def status_command(
    project_dir: Path,
    watch: bool = typer.Option(False, "--watch", "-w", help="Refresh the dashboard until the run reaches SUCCESS, FAILED or INTERRUPTED. Ctrl-C stops watching only; the run is never signalled."),
    interval: float = typer.Option(10.0, "--interval", min=2.0, help="Seconds between --watch refreshes (default 10, minimum 2)."),
    case_id: str | None = typer.Option(None, "--case", "--case-id", help="Show the latest run of this case ID (or, with --run, that exact run)."),
    run_id: str | None = typer.Option(None, "--run", help="Show this exact run ID; add --case if the run ID exists in more than one case."),
    all_runs: bool = typer.Option(False, "--all", help="One-line summary of every recorded run (optionally only --case), newest first."),
) -> None:
    """Execution dashboard of the latest run; read-only, Nextflow need not be running.

    Shows pipeline stages for the run's actual route, running tasks, task counts
    and the configured resource ceiling. Scientific results are in the report
    and delivery package, not here.
    """

    from rnaseq.run_status import _Cache, build_status, load_run, matching_runs, render_history, render_status

    project = project_dir.resolve()
    if not project.is_dir():
        typer.echo(f"ERROR: project directory not found: {project}", err=True)
        raise typer.Exit(code=1)
    cache = _Cache()
    if all_runs:
        if watch or run_id:
            typer.echo("ERROR: --all is a history summary; it cannot be combined with --watch or --run.", err=True)
            raise typer.Exit(code=2)
        runs = matching_runs(project, case_id)
        if not runs:
            typer.echo(f"No recorded runs for case {case_id}." if case_id else "No recorded case runs.")
            raise typer.Exit(code=1 if case_id else 0)
        latest = matching_runs(project)[:1] if case_id else runs[:1]
        typer.echo(render_history(project, runs, cache, latest=latest[0][0] if latest else None))
        return

    matches = matching_runs(project, case_id, run_id)
    if run_id and not case_id and len({path.parent for path, _state in matches}) > 1:
        cases = ", ".join(sorted({str(state.get("case_id") or path.parent.name) for path, state in matches}))
        typer.echo(f"ERROR: run ID {run_id} exists in more than one case ({cases}); add --case.", err=True)
        raise typer.Exit(code=2)
    if not matches:
        if case_id or run_id:
            wanted = " / ".join(item for item in (case_id and f"case {case_id}", run_id and f"run {run_id}") if item)
            typer.echo(f"ERROR: no recorded run matches {wanted}; list runs with: rnaseq status PROJECT --all", err=True)
            raise typer.Exit(code=1)
        planned = (project / "planning" / "manifest.preview.yaml").is_file()
        typer.echo("PLANNED: no case runs recorded yet; start one with: rnaseq run PROJECT --case-id CASE-ID"
                   if planned else "No recorded case runs.")
        return
    pinned = matches[0][0]  # --watch follows this run even if a newer one starts
    ascii_only = not (getattr(sys.stdout, "encoding", None) or "").lower().replace("-", "").startswith("utf")

    def snapshot() -> tuple[str, bool]:
        selected = load_run(pinned, cache)
        if selected is None:
            return f"ERROR: run state of {pinned} is no longer readable.", True
        status = build_status(*selected, cache=cache)
        return render_status(status, ascii_only=ascii_only), status.is_final

    text, final = snapshot()
    if not watch:
        typer.echo(text)
        return
    refresh = sys.stdout.isatty()
    try:
        while True:
            if refresh:
                typer.echo("\033[H\033[2J", nl=False)  # redraw in place instead of scrolling copies
            typer.echo(text)
            if final:
                return
            typer.echo(f"\nRefreshing every {interval:g}s; Ctrl-C stops watching (the run is not affected).")
            time.sleep(interval)
            text, final = snapshot()
            if not refresh:
                typer.echo("")
    except KeyboardInterrupt:
        typer.echo("\nStopped watching; the run itself was not affected.")
        raise typer.Exit(code=0)


@app.command("clean")
def clean_command(
    project_dir: Path,
    case_id: str = typer.Option(..., "--case-id", "--case", help="Case ID of the run whose work directory is removed."),
    run_id: str = typer.Option(..., "--run", help="Exact run ID (YYYYMMDD-HHMMSS+ZZZZ) of that case."),
    force: bool = typer.Option(False, "--force", help="Also remove the work of a FAILED or INTERRUPTED run (never of a running run)."),
    dry_run: bool = typer.Option(False, "--dry-run", help="Show what would be removed; delete nothing."),
    yes: bool = typer.Option(False, "--yes", help="Delete without an interactive confirmation."),
) -> None:
    """Remove one finished run's disposable Nextflow work directory; the durable run record is kept.

    The path comes from the run's own provenance, never from the arguments.
    SUCCESS runs are eligible; FAILED and INTERRUPTED runs need --force because
    their work is what diagnosis and retry use; a running run is always refused.
    """

    from rnaseq.cleanup import AlreadyAbsent, check_no_live_processes, execute_work_cleanup, format_bytes, measure_work, plan_work_cleanup

    try:
        plan = plan_work_cleanup(project_dir, case_id, run_id)
    except ExecutionPreflightError as exc:
        typer.echo(f"ERROR: {exc}", err=True)
        raise typer.Exit(code=1) from exc
    except (OSError, UnicodeError) as exc:
        typer.echo(f"SYSTEM ERROR: {exc}", err=True)
        raise typer.Exit(code=2) from exc
    state = plan.effective_state if plan.effective_state == plan.recorded_state else f"{plan.effective_state} (recorded {plan.recorded_state})"
    typer.echo(f"Run: {plan.case_id}/{plan.run_id}")
    typer.echo(f"State: {state}" + (f" - {plan.state_note}" if plan.state_note else ""))
    typer.echo(f"Work directory: {plan.work_dir}")
    typer.echo(f"  recorded in: {plan.source}")
    typer.echo(f"Durable run record (never removed): {plan.run_dir}")
    if not plan.exists:
        typer.echo("Nothing to clean: the work directory does not exist (already cleaned or never created).")
        return
    refusal = plan.refusal(force=force) or check_no_live_processes(plan)
    if refusal:
        typer.echo(f"REFUSED: {refusal}", err=True)
        raise typer.Exit(code=1)
    typer.echo("Measuring work directory size ...")
    usage = measure_work(plan.work_dir)
    unreadable = f"; {usage.unreadable} entries could not be read" if usage.unreadable else ""
    typer.echo(f"Work size: about {format_bytes(usage.bytes)} in {usage.files:,} files{unreadable}")
    if dry_run:
        typer.echo("Dry run: nothing was removed.")
        return
    if not yes and not typer.confirm(f"Permanently remove {plan.work_dir}?", default=False):
        typer.echo("Cleanup cancelled; nothing was removed.")
        return
    try:
        warning = execute_work_cleanup(plan, force=force, usage=usage)
    except AlreadyAbsent:
        typer.echo("Nothing to clean: the work directory disappeared before deletion (another cleanup may have removed it).")
        return
    except ExecutionPreflightError as exc:
        typer.echo(f"REFUSED: {exc}", err=True)
        raise typer.Exit(code=1) from exc
    except OSError as exc:
        typer.echo(f"SYSTEM ERROR: cleanup stopped part-way: {exc}. The durable run record is untouched; re-run to finish.", err=True)
        raise typer.Exit(code=2) from exc
    typer.echo(f"Cleaned: removed {plan.work_dir} (about {format_bytes(usage.bytes)} freed).")
    typer.echo(f"Durable run record preserved: {plan.run_dir}")
    if warning:
        typer.echo(f"WARNING: {warning}", err=True)


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


@app.command("gui")
def gui_command(
    port: int = typer.Option(8501, "--port", min=1, max=65535, help="Local port; the GUI listens on 127.0.0.1 only."),
) -> None:
    """Start the browser GUI on 127.0.0.1 (needs the optional 'gui' extra)."""

    from rnaseq.gui import GUI_ADDRESS, INSTALL_HINT, launch, streamlit_available

    if not streamlit_available():
        typer.echo(f"ERROR: the GUI needs the optional 'gui' extra (Streamlit). Install it with: {INSTALL_HINT}", err=True)
        raise typer.Exit(code=1)
    typer.echo(f"nf-rna GUI on http://{GUI_ADDRESS}:{port}  (from another machine: ssh -L {port}:127.0.0.1:{port} <server>)")
    launch(port)


@app.command("doctor")
def doctor_command(project_dir: Path | None = typer.Argument(None, help="Optional project for adopted-reference readiness.")) -> None:
    """Report non-mutating execution prerequisites; nothing is installed automatically."""

    checks = doctor_checks(project_dir)
    for item in checks:
        typer.echo(f"{item.name}: {item.verdict} — {item.detail}")
    ready, failed = doctor_readiness(checks)
    if ready:
        typer.echo("Overall: READY")
    else:
        typer.echo(f"Overall: NOT READY — failed: {', '.join(failed)}")
        raise typer.Exit(code=1)


if __name__ == "__main__":
    app()
