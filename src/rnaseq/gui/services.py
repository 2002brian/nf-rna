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


BACKEND_LABELS = {"salmon": "Salmon", "hisat2_featurecounts": "HISAT2"}


def reference_unusable_reason(reference: LocalReference, species: str, method: str, *, require_production: bool = True) -> str | None:
    """Why a managed reference cannot serve a project of ``species`` with ``method`` (None when it can).

    The same rule as ``compatible_registered_references`` and the project
    validator: same species, a production manifest, and the backend's own
    index (a Salmon index never stands in for a HISAT2 genome index).
    """

    if reference.species != species:
        return f"it is a {reference.species} reference, not {species}"
    if require_production and reference.purpose != "production":
        return f"its manifest purpose is {reference.purpose!r}, not production"
    if method == "salmon" and reference.salmon_index is None:
        return f"it has no Salmon index (build one with `rnaseq reference prepare {reference.root}`)"
    if method == "hisat2_featurecounts":
        if reference.hisat2_index is None:
            return f"it has no HISAT2 index (build one with `rnaseq reference prepare-hisat2 {reference.root}`)"
        if not reference.hisat2_runtime_ready:
            return "its HISAT2 index has no validated runtime compatibility (manifest hisat2.runtime_compatibility)"
    return None


@dataclass(frozen=True)
class ReferenceAvailability:
    compatible: list[LocalReference]
    unusable: list[tuple[LocalReference, str]]  # registered for this species, with the reason it cannot be used
    error: str | None = None


def reference_availability(species: str, method: str) -> ReferenceAvailability:
    """Registered references for a species and backend, from the workstation registry only (nothing is scanned)."""

    from rnaseq.references import compatible_registered_references

    try:
        registered = registered_local_references()
        compatible = compatible_registered_references(species, method)
    except LocalReferenceError as exc:
        return ReferenceAvailability([], [], str(exc))
    usable = {(item.root, item.manifest_path) for item in compatible}
    unusable = []
    for item in registered:
        if item.species != species or (item.root, item.manifest_path) in usable:
            continue
        unusable.append((item, reference_unusable_reason(item, species, method) or "not selected by the registry rule"))
    return ReferenceAvailability(compatible, unusable)


def check_project_reference(reference: dict[str, Any] | None, species: str, method: str) -> str | None:
    """Re-check a managed (``source: local``) reference against the project before it is created.

    Loads the manifest with the same loader as the project validator and
    applies :func:`reference_unusable_reason`; returns the problem, or None.
    Custom FASTA/GTF and iGenomes references are checked by the validator.
    """

    from rnaseq.references import _load_local_reference_root

    if not reference or reference.get("source") != "local":
        return None
    root = Path(str(reference.get("root") or ""))
    manifest = str(reference.get("manifest") or "reference_manifest.yaml")
    try:
        loaded, _manifest = _load_local_reference_root(root.expanduser().resolve(), manifest, None)
    except LocalReferenceError as exc:
        return f"The reference at {root} cannot be used: {exc}"
    reason = reference_unusable_reason(loaded, species, method, require_production=False)
    if reason:
        return f"The reference at {root} cannot be used for a {species} {BACKEND_LABELS.get(method, method)} project: {reason}."
    return None


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


def metadata_file_columns(report: ValidationReport) -> list[str]:
    """The header of metadata.csv as written (empty when the file is missing)."""

    assert report.loaded is not None
    return _csv_header(report.loaded.metadata_path)


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


@dataclass(frozen=True)
class MetadataChange:
    """The table that results from applying pasted, imported or assigned values to the table being edited.

    Nothing is dropped: samples absent from the new values keep every value,
    columns absent from them are kept, and an empty new cell leaves the current
    value unchanged.  Non-empty values that would be replaced are listed in
    ``overwritten`` so the GUI can ask before applying them.
    """

    columns: list[str]
    rows: list[dict[str, str]]
    filled: int = 0
    overwritten: tuple[tuple[str, str, str, str], ...] = ()  # sample_id, column, current value, new value
    added: tuple[str, ...] = ()
    new_columns: tuple[str, ...] = ()
    errors: tuple[str, ...] = ()

    @property
    def changes(self) -> int:
        return self.filled + len(self.overwritten) + len(self.added) + len(self.new_columns)


def clean_cell(value: Any) -> str:
    """A table cell as stripped text; None and NaN (an empty data_editor cell) become ""."""

    return _clean(value)


def _clean(value: Any) -> str:
    if value is None:
        return ""
    try:
        if value != value:  # NaN from an empty data_editor cell
            return ""
    except TypeError:
        pass
    return str(value).strip()


def _table_rows(columns: list[str], rows: list[dict[str, Any]]) -> list[dict[str, str]]:
    return [{column: _clean(row.get(column)) for column in columns} for row in rows]


def merge_metadata_table(columns: list[str], rows: list[dict[str, Any]], new_header: list[str], new_rows: list[dict[str, Any]]) -> MetadataChange:
    """Merge a pasted or imported table into the current one by ``sample_id``."""

    current = _table_rows(columns, rows)
    if "sample_id" not in new_header:
        return MetadataChange(list(columns), current, errors=("The new table needs a header row with a sample_id column.",))
    incoming = _table_rows(new_header, new_rows)
    errors: list[str] = []
    blank = sum(1 for row in incoming if not row["sample_id"] and any(row.values()))
    if blank:
        errors.append(f"{blank} row(s) of the new table have no sample_id.")
    seen: dict[str, int] = {}
    for row in incoming:
        if row["sample_id"]:
            seen[row["sample_id"]] = seen.get(row["sample_id"], 0) + 1
    duplicated = sorted(sample for sample, count in seen.items() if count > 1)
    if duplicated:
        errors.append("sample_id repeated in the new table: " + ", ".join(duplicated[:10]))
    index: dict[str, list[int]] = {}
    for position, row in enumerate(current):
        if row.get("sample_id"):
            index.setdefault(row["sample_id"], []).append(position)
    ambiguous = sorted(sample for sample in seen if len(index.get(sample, [])) > 1)
    if ambiguous:
        errors.append("sample_id repeated in the table being edited (fix these rows first): " + ", ".join(ambiguous[:10]))
    if errors:
        return MetadataChange(list(columns), current, errors=tuple(errors))
    new_columns = [column for column in new_header if column and column not in columns]
    merged_columns = [*columns, *new_columns]
    merged = [{column: row.get(column, "") for column in merged_columns} for row in current]
    filled, overwritten, added = 0, [], []
    for row in incoming:
        sample = row["sample_id"]
        if not sample:
            continue
        if sample not in index:
            merged.append({column: row.get(column, "") for column in merged_columns})
            added.append(sample)
            continue
        target = merged[index[sample][0]]
        for column in new_header:
            if column == "sample_id" or not column:
                continue
            value, old = row[column], target.get(column, "")
            if not value or value == old:
                continue
            if old:
                overwritten.append((sample, column, old, value))
            else:
                filled += 1
            target[column] = value
    return MetadataChange(merged_columns, merged, filled, tuple(overwritten), tuple(added), tuple(new_columns))


def assign_metadata_value(columns: list[str], rows: list[dict[str, Any]], samples: list[str], column: str, value: str) -> MetadataChange:
    """Set ``column`` to ``value`` for the rows of ``samples``; every other cell is unchanged."""

    current = _table_rows(columns, rows)
    if column not in columns:
        return MetadataChange(list(columns), current, errors=(f"Unknown column {column!r}.",))
    value = _clean(value)
    wanted = set(samples)
    filled, overwritten = 0, []
    for row in current:
        if row.get("sample_id") not in wanted or row[column] == value:
            continue
        if row[column]:
            overwritten.append((row["sample_id"], column, row[column], value))
        else:
            filled += 1
        row[column] = value
    return MetadataChange(list(columns), current, filled, tuple(overwritten))


@dataclass(frozen=True)
class MetadataDifference:
    changed_cells: int
    added_rows: tuple[str, ...]
    removed_rows: tuple[str, ...]
    columns_added: tuple[str, ...]
    columns_removed: tuple[str, ...]

    @property
    def unsaved(self) -> bool:
        return bool(self.changed_cells or self.added_rows or self.removed_rows or self.columns_added or self.columns_removed)

    def describe(self) -> str:
        parts = []
        if self.changed_cells:
            parts.append(f"{self.changed_cells} changed value(s)")
        if self.added_rows:
            parts.append(f"{len(self.added_rows)} row(s) not in metadata.csv ({', '.join(self.added_rows[:8])}{' …' if len(self.added_rows) > 8 else ''})")
        if self.removed_rows:
            parts.append(f"{len(self.removed_rows)} row(s) of metadata.csv removed ({', '.join(self.removed_rows[:8])}{' …' if len(self.removed_rows) > 8 else ''})")
        if self.columns_added:
            parts.append("new column(s) " + ", ".join(self.columns_added))
        if self.columns_removed:
            parts.append("column(s) of metadata.csv missing: " + ", ".join(self.columns_removed))
        return "; ".join(parts)


def metadata_difference(saved_columns: list[str], saved_rows: list[dict[str, Any]], columns: list[str], rows: list[dict[str, Any]]) -> MetadataDifference:
    """How the table being edited differs from metadata.csv (rows matched by sample_id, in order for blank IDs).

    Rows and columns that are entirely empty are ignored, as :func:`save_metadata` does not write them.
    """

    def keyed(names: list[str], table: list[dict[str, Any]]) -> list[tuple[str, dict[str, str]]]:
        cleaned = [row for row in _table_rows(names, table) if any(row.values())]
        return [(row.get("sample_id", "") or f"(row {number})", row) for number, row in enumerate(cleaned, 1)]

    saved_names = [column for column in saved_columns if column]
    names = [column for column in columns if column]
    before, after = dict(keyed(saved_names, saved_rows)), dict(keyed(names, rows))
    used = [column for column in names if any(row.get(column) for row in after.values()) or column in saved_names]
    columns_added = tuple(column for column in used if column not in saved_names)
    columns_removed = tuple(column for column in saved_names if column not in names)
    shared = [column for column in names if column in saved_names]
    changed = sum(1 for sample, row in after.items() if sample in before for column in shared if row.get(column, "") != before[sample].get(column, ""))
    changed += sum(1 for sample, row in after.items() if sample in before for column in columns_added if row.get(column))
    return MetadataDifference(
        changed, tuple(sample for sample in after if sample not in before), tuple(sample for sample in before if sample not in after),
        columns_added, columns_removed,
    )


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


# --------------------------------------------------------------------- paired-design preview

@dataclass(frozen=True)
class PairRow:
    pair: str
    samples: dict[str, tuple[str, ...]]  # condition level -> sample IDs
    status: str  # "complete" or a short description of the problem


@dataclass(frozen=True)
class PairingPreview:
    pair_column: str
    condition_column: str
    levels: tuple[str, ...]
    pairs: tuple[PairRow, ...]
    blank_pair: tuple[str, ...]  # samples without a pair ID
    blank_condition: tuple[str, ...]  # samples without a condition

    @property
    def complete(self) -> int:
        return sum(row.status == "complete" for row in self.pairs)

    @property
    def problems(self) -> int:
        return len(self.pairs) - self.complete + len(self.blank_pair)


def pairing_preview(columns: list[str], rows: list[dict[str, Any]], pair_column: str, condition_column: str) -> PairingPreview | None:
    """How the samples of a (possibly unsaved) metadata table form pairs.

    Display only, for the table being edited: it applies the validator's rule
    for paired_two_group designs (each pair ID has exactly one sample of each
    of the two condition levels; no blank pair ID) so mistakes show before
    saving. ``rnaseq validate`` remains the authority. None when either column
    is missing from the table.
    """

    if pair_column not in columns or condition_column not in columns:
        return None
    levels: list[str] = []
    by_pair: dict[str, dict[str, list[str]]] = {}
    blank_pair: list[str] = []
    blank_condition: list[str] = []
    for row in rows:
        sample = clean_cell(row.get("sample_id"))
        if not sample:
            continue
        pair, level = clean_cell(row.get(pair_column)), clean_cell(row.get(condition_column))
        if not level:
            blank_condition.append(sample)
        elif level not in levels:
            levels.append(level)
        if not pair:
            blank_pair.append(sample)
            continue
        by_pair.setdefault(pair, {}).setdefault(level, []).append(sample)
    out = []
    for pair in sorted(by_pair):
        members = by_pair[pair]
        problems = []
        for level in levels:
            count = len(members.get(level, ()))
            if count == 0:
                problems.append(f"no {level} sample")
            elif count > 1:
                problems.append(f"{count} {level} samples (exactly one needed)")
        if members.get(""):
            problems.append(f"{len(members[''])} sample(s) without a condition")
        if len(levels) != 2:
            problems.append(f"a paired design needs exactly 2 condition levels (found {len(levels)})")
        out.append(PairRow(pair, {level: tuple(members.get(level, ())) for level in levels}, "complete" if not problems else "; ".join(problems)))
    return PairingPreview(pair_column, condition_column, tuple(levels), tuple(out), tuple(blank_pair), tuple(blank_condition))
