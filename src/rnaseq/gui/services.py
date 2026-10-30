"""GUI-facing helpers composed from the existing nf-rna services (no Streamlit here).

Every scientific decision -- schema, validation, design, contrasts, run
preflight and run state -- is delegated to the modules the CLI uses.
"""

from __future__ import annotations

import csv
import io
import os
import tempfile
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from rnaseq.design_setup import CONTRAST_ID_PATTERN, DesignProposal, build_proposal, discover_levels
from rnaseq.enrichment_config import KEGG_CODES, latest_kegg_snapshot
from rnaseq.errors import ProjectConfigError
from rnaseq.gui.config import GuiConfig
from rnaseq.gui.jobs import run_display_state
from rnaseq.models import InputType, Preset, SequencingLayout, production_enrichment_backends
from rnaseq.project import load_project
from rnaseq.references import LocalReference, LocalReferenceError, registered_local_references
from rnaseq.run_status import _Cache, discover_runs
from rnaseq.validators import (
    ContrastDefinition, FastqSummary, ValidationIssue, ValidationReport, condition_factor, parse_formula,
    validate_fastq, validate_metadata, validate_project,
)

_SKIPPED_DIRS = {"runs", "input", "planning", "downstream", "delivery", "work", "reference", "fastq", "raw-data", "archives", "references"}
# Metadata findings that make the table unusable; anything else (an unfilled
# design value) may be saved and completed later.
_STRUCTURAL_METADATA_CODES = {
    "empty_metadata", "blank_metadata_column", "duplicate_metadata_column", "missing_sample_id_column",
    "missing_design_variables", "blank_sample_id", "duplicate_sample_id", "malformed_metadata_row", "missing_metadata_rows",
}
SAMPLE_AGREEMENT_CODES = {"count_only_samples", "metadata_only_samples", "fastq_only_samples", "metadata_only_fastq_samples"}


# --------------------------------------------------------------------- projects and runs

def discover_projects(config: GuiConfig) -> list[Path]:
    """Projects (directories with project.yaml) up to two levels below each available root."""

    found: set[Path] = set()
    for root in config.available_roots():
        level = [root]
        for _depth in range(2):
            children: list[Path] = []
            for directory in level:
                try:
                    entries = [entry for entry in directory.iterdir() if entry.is_dir() and not entry.name.startswith(".")]
                except OSError:
                    continue
                for entry in entries:
                    if (entry / "project.yaml").is_file():
                        found.add(entry.resolve())
                    elif entry.name not in _SKIPPED_DIRS:
                        children.append(entry)
            level = children
    return sorted(found)


@dataclass(frozen=True)
class ProjectRow:
    path: Path
    project_id: str
    preset: str
    input_type: str
    latest_run: str
    latest_state: str
    error: str | None = None


def project_row(path: Path, cache: _Cache | None = None) -> ProjectRow:
    try:
        config = load_project(path).config
    except ProjectConfigError as exc:
        return ProjectRow(path, path.name, "?", "?", "", "", error=str(exc))
    runs = discover_runs(path)
    latest, shown = "", "no runs"
    if runs:
        run_dir, state = runs[0]
        latest = f"{run_dir.parent.name}/{run_dir.name}"
        shown = run_display_state(run_dir, state, cache)[0]
    return ProjectRow(path, config.project.id, config.project.preset.value, config.input.type.value, latest, shown)


def recent_runs(projects: list[Path], limit: int = 15, cache: _Cache | None = None) -> list[dict[str, str]]:
    cache = cache or _Cache()
    rows = []
    for project in projects:
        for run_dir, state in discover_runs(project)[:limit]:
            shown, note = run_display_state(run_dir, state, cache)
            rows.append({
                "started": str(state.get("started_at") or ""), "project": project.name,
                "case": run_dir.parent.name, "run": run_dir.name, "state": shown, "note": note or "",
                "path": str(run_dir),
            })
    return sorted(rows, key=lambda row: row["started"], reverse=True)[:limit]


def registered_references() -> tuple[list[LocalReference], str | None]:
    try:
        return registered_local_references(), None
    except LocalReferenceError as exc:
        return [], str(exc)


# --------------------------------------------------------------------- metadata

def _csv_header(path: Path) -> list[str]:
    try:
        with path.open(encoding="utf-8", newline="") as handle:
            return next(csv.reader(handle), [])
    except (OSError, UnicodeError, csv.Error):
        return []


def design_variables(report: ValidationReport) -> tuple[str, ...]:
    assert report.config is not None
    return parse_formula(report.config.design.formula, ValidationReport(project_dir=report.project_dir))


def metadata_columns(report: ValidationReport) -> list[str]:
    """sample_id, the design variables, then any further columns already in metadata.csv."""

    assert report.loaded is not None
    existing = [column for column in _csv_header(report.loaded.metadata_path) if column]
    columns = ["sample_id", *design_variables(report)]
    return list(dict.fromkeys([*columns, *existing]))


def read_metadata(report: ValidationReport) -> list[dict[str, str]]:
    assert report.loaded is not None
    try:
        with report.loaded.metadata_path.open(encoding="utf-8", newline="") as handle:
            return [dict(row) for row in csv.DictReader(handle) if any((value or "").strip() for value in row.values())]
    except (OSError, UnicodeError, csv.Error):
        return []


def input_sample_ids(report: ValidationReport) -> list[str]:
    """Sample IDs from the count-matrix header or the FASTQ file names (never their conditions)."""

    if report.counts is not None:
        return list(report.counts.sample_ids)
    if report.fastq is not None:
        return list(report.fastq.sample_ids)
    return []


def merge_sample_ids(rows: list[dict[str, str]], columns: list[str], sample_ids: list[str]) -> list[dict[str, str]]:
    present = {row.get("sample_id", "") for row in rows}
    merged = [{column: row.get(column, "") or "" for column in columns} for row in rows]
    merged.extend({column: (sample if column == "sample_id" else "") for column in columns} for sample in sample_ids if sample not in present)
    return merged


def parse_table_text(text: str) -> tuple[list[str], list[dict[str, str]]]:
    """A pasted table: tab-separated (as copied from Excel) or comma-separated, header first."""

    lines = [line for line in text.replace("\r\n", "\n").replace("\r", "\n").split("\n") if line.strip()]
    if not lines:
        return [], []
    delimiter = "\t" if "\t" in lines[0] else ","
    parsed = list(csv.reader(lines, delimiter=delimiter))
    header = [name.strip() for name in parsed[0]]
    rows = [{name: (values[index].strip() if index < len(values) else "") for index, name in enumerate(header)} for values in parsed[1:]]
    return header, rows


@dataclass
class MetadataSave:
    saved: bool
    errors: list[ValidationIssue] = field(default_factory=list)
    incomplete: list[ValidationIssue] = field(default_factory=list)
    agreement: list[ValidationIssue] = field(default_factory=list)


def save_metadata(project: Path, columns: list[str], rows: list[dict[str, Any]]) -> MetadataSave:
    """Write metadata.csv when its structure is valid; unfilled design values may be completed later.

    The table is checked with the project's metadata validator before it
    replaces metadata.csv; a structurally invalid table is not written.
    """

    loaded = load_project(project)
    check = ValidationReport(project_dir=loaded.root)
    variables = parse_formula(loaded.config.design.formula, check)
    buffer = io.StringIO()
    writer = csv.writer(buffer, lineterminator="\n")
    writer.writerow(columns)
    for row in rows:
        values = ["" if row.get(column) is None else str(row.get(column)).strip() for column in columns]
        if any(values):
            writer.writerow(values)
    target = loaded.metadata_path
    descriptor, temporary = tempfile.mkstemp(prefix=f".{target.name}.gui-", dir=target.parent)
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8", newline="\n") as handle:
            handle.write(buffer.getvalue())
        validate_metadata(Path(temporary), variables, loaded.config.design.variables, check)
        errors = [issue for issue in check.errors if issue.code in _STRUCTURAL_METADATA_CODES]
        if errors:
            return MetadataSave(False, errors=errors)
        os.replace(temporary, target)
    finally:
        Path(temporary).unlink(missing_ok=True)
    after = validate_project(loaded.root)
    return MetadataSave(
        True,
        incomplete=[issue for issue in check.issues if issue.code not in _STRUCTURAL_METADATA_CODES],
        agreement=[issue for issue in after.issues if issue.code in SAMPLE_AGREEMENT_CODES],
    )


def issue_column(issue: ValidationIssue, columns: list[str]) -> str:
    """The metadata column a validator finding refers to, for display beside that column."""

    if issue.code in {"blank_sample_id", "duplicate_sample_id", "missing_sample_id_column"} | SAMPLE_AGREEMENT_CODES:
        return "sample_id"
    for column in columns:
        if f"'{column}'" in issue.message or f"variable {column!r}" in issue.message:
            return column
    return "table"


# --------------------------------------------------------------------- FASTQ input

def discover_fastq(directory: Path, layout: SequencingLayout) -> tuple[FastqSummary | None, list[ValidationIssue]]:
    """Sample IDs and R1/R2 pairs by the validator's FASTQ naming rule (nothing is read or copied)."""

    report = ValidationReport(project_dir=directory)
    summary = validate_fastq(directory, layout, report)
    return summary, report.errors


# --------------------------------------------------------------------- design and contrasts

@dataclass(frozen=True)
class DesignView:
    factor: str | None
    levels: tuple[str, ...]
    counts: dict[str, int]
    formula: str
    paired: bool
    pair_id: str | None


def design_view(report: ValidationReport) -> DesignView:
    assert report.config is not None
    design = report.config.design
    factor = condition_factor(report.formula_variables, report.design_variable_types)
    levels: tuple[str, ...] = ()
    counts: dict[str, int] = {}
    if factor and report.metadata is not None and factor in report.metadata.columns:
        observed = discover_levels(report.metadata, factor)
        frozen = tuple(design.condition_order or ())
        levels = frozen if frozen and set(frozen) == set(observed) else observed
        counts = {level: sum(1 for row in report.metadata.rows if row[factor] == level) for level in levels}
    return DesignView(factor, levels, counts, design.formula, design.pair_id is not None, design.pair_id)


def contrast_definition(factor: str, numerator: str, denominator: str) -> ContrastDefinition:
    """The contrast ID rule used by ``rnaseq design``: ``<numerator>_vs_<denominator>``."""

    return ContrastDefinition(f"{numerator}_vs_{denominator}", factor, numerator, denominator)


def unsafe_levels(levels: tuple[str, ...]) -> list[str]:
    return [level for level in levels if not CONTRAST_ID_PATTERN.fullmatch(level)]


def design_proposal(report: ValidationReport, levels: tuple[str, ...], contrasts: list[ContrastDefinition]) -> DesignProposal:
    """The proposal ``rnaseq design`` would write, with the user's chosen contrasts."""

    base = build_proposal(report, condition_order=levels)
    return DesignProposal(base.factor, base.levels, base.counts, tuple(contrasts), base.formula, base.standard_formula, base.has_batch)


@dataclass(frozen=True)
class KeggPin:
    applies: bool
    pinned: str | None
    latest: str | None
    warning: str | None


def kegg_pin_status(report: ValidationReport) -> KeggPin:
    config = report.config
    if config is None or config.annotation is None:
        return KeggPin(False, None, None, None)
    backends = production_enrichment_backends(config.analysis.enrichment if config.analysis else ())
    code = KEGG_CODES.get(config.annotation.organism)
    if code is None or not ({"kegg", "gsea-kegg"} & set(backends)):
        return KeggPin(False, None, None, None)
    latest = latest_kegg_snapshot(code)
    warning = next((issue.message for issue in report.issues if issue.code in {"kegg_snapshot_not_pinned", "kegg_snapshot_missing"}), None)
    return KeggPin(True, config.annotation.enrichment.kegg.snapshot, latest[0] if latest else None, warning)


# --------------------------------------------------------------------- running

def reuse_candidates(report: ValidationReport) -> list[str]:
    """``CASE-ID/RUN-ID`` of SUCCESS runs accepted by the existing ``--reuse-upstream`` preflight.

    The byte-level comparison of frozen upstream inputs happens in the run
    itself, exactly as on the command line.
    """

    from rnaseq.errors import ExecutionPreflightError
    from rnaseq.service import preflight_reuse_upstream

    if report.config is None or report.config.input.type is not InputType.FASTQ:
        return []
    candidates = []
    for run_dir, state in discover_runs(report.project_dir):
        if state.get("status") != "SUCCESS":
            continue
        reference = f"{run_dir.parent.name}/{run_dir.name}"
        try:
            preflight_reuse_upstream(report, reference)
        except ExecutionPreflightError:
            continue
        candidates.append(reference)
    return candidates


def run_preflight(project: Path, case_id: str, reuse_upstream: str | None) -> ValidationReport:
    """The checks ``rnaseq run`` performs before it creates a run directory; raises on refusal."""

    from rnaseq.service import prepare_service_run, preflight_reuse_upstream, validate_case_id

    report = validate_project(project)
    validate_case_id(case_id)
    if reuse_upstream and report.config is not None and report.is_valid:
        preflight_reuse_upstream(report, reuse_upstream)
    prepare_service_run(report, profile="local")
    return report


def scope_label(report: ValidationReport) -> str:
    assert report.config is not None
    return {Preset.QC: "QC (technical FASTQ QC only)", Preset.L1: "L1 (QC and exploratory analysis)", Preset.L2: "L2 (L1 + differential expression and enrichment)"}[report.config.project.preset]
