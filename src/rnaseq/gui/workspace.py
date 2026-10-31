"""Project workspace views for the GUI (no Streamlit here).

Everything is read from the files nf-rna already writes: ``project.yaml``,
``metadata.csv``, ``contrasts.csv``, each run's ``run_state.json``, frozen
contract, provenance and downstream records. Nothing is written and no
statistic is recomputed.

Two separate notions are shown for a run:

* **Execution status** (RUNNING, SUCCESS, FAILED, INCOMPLETE) -- whether the
  run finished, from ``run_state.json`` and process liveness
  (:func:`rnaseq.gui.jobs.run_display_state`).
* **Scientific QC** (PASS, WARNING, INCONCLUSIVE, NOT_EVALUATED) -- a display
  roll-up of verdicts the pipeline itself recorded: the strandedness check
  (v1.6.1 F-1), the L1 sample flags, the R warnings of the DESeq2 fit and
  shrinkage (v1.6.1 F-3) and the shrinkage status. Each recorded verdict is
  shown verbatim next to the roll-up; the roll-up never changes a verdict and
  never changes the execution status.
"""

from __future__ import annotations

import csv
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from rnaseq.errors import ProjectConfigError
from rnaseq.gui.jobs import run_display_state
from rnaseq.gui.results import read_json, read_yaml
from rnaseq.models import PIPELINE_VERSION, InputType, Preset, ProjectConfig, describe_enrichment
from rnaseq.project import load_project
from rnaseq.run_status import _Cache, discover_runs
from rnaseq.validators import ValidationIssue, ValidationReport

EXECUTION_STATES = ("RUNNING", "SUCCESS", "FAILED", "INCOMPLETE")
QC_LEVELS = ("PASS", "WARNING", "INCONCLUSIVE", "NOT_EVALUATED")
QC_ICONS = {"PASS": "🟢", "WARNING": "🟠", "INCONCLUSIVE": "⚪", "NOT_EVALUATED": "·"}
EXECUTION_ICONS = {"RUNNING": "🔵", "SUCCESS": "✅", "FAILED": "❌", "INCOMPLETE": "⏸"}

# Run-level strandedness status (hisat2_featurecounts.evaluate_strandedness) -> QC level.
# ADVISORY (declared unstranded, reads look stranded; unstranded counting is valid)
# is surfaced as WARNING so it is looked at; its own text says counting is valid.
STRAND_LEVELS = {"CONSISTENT": "PASS", "WARNING": "WARNING", "ADVISORY": "WARNING", "INCONCLUSIVE": "INCONCLUSIVE", "NOT_EVALUATED": "NOT_EVALUATED"}

METHOD_LABELS = {"salmon": "Salmon (+ tximport)", "hisat2_featurecounts": "HISAT2 + featureCounts"}


# --------------------------------------------------------------------- execution status

def execution_status(shown: str) -> str:
    """The four execution states of the GUI, from the displayed run state.

    CREATED/PLANNED/RUNNING with a live (or unverifiable) process is RUNNING;
    INTERRUPTED (recorded or stale) and UNKNOWN are INCOMPLETE.
    """

    if shown in {"RUNNING", "PLANNED", "CREATED"}:
        return "RUNNING"
    if shown in {"SUCCESS", "FAILED"}:
        return shown
    return "INCOMPLETE"


# --------------------------------------------------------------------- scientific QC

@dataclass(frozen=True)
class QcItem:
    check: str
    level: str
    recorded: str
    detail: str


@dataclass(frozen=True)
class RunQc:
    level: str
    items: tuple[QcItem, ...]

    @property
    def warnings(self) -> list[QcItem]:
        return [item for item in self.items if item.level == "WARNING"]

    def headline(self) -> str:
        if self.level == "WARNING":
            return f"WARNING ({len(self.warnings)})"
        return self.level


def roll_up(levels: list[str]) -> str:
    for level in ("WARNING", "INCONCLUSIVE", "PASS"):
        if level in levels:
            return level
    return "NOT_EVALUATED"


def _strandedness_item(run_dir: Path, contract: dict[str, Any], input_type: str | None) -> QcItem | None:
    if input_type != InputType.FASTQ.value:
        return None  # raw counts: no reads, nothing to check
    upstream_qc = contract.get("upstream_qc") if isinstance(contract.get("upstream_qc"), dict) else {}
    strand = upstream_qc.get("strandedness") if isinstance(upstream_qc.get("strandedness"), dict) else None
    if strand is None:
        provenance = read_yaml(run_dir / "provenance" / "run_provenance.yaml") or {}
        strand = provenance.get("strandedness_qc") if isinstance(provenance.get("strandedness_qc"), dict) else None
    if strand is None:
        return QcItem("Strandedness check", "NOT_EVALUATED", "not recorded", "No strandedness record in this run (runs before v1.6.1, or the upstream stage did not finish).")
    status = str(strand.get("status") or "NOT_EVALUATED")
    samples = strand.get("samples") if isinstance(strand.get("samples"), dict) else {}
    flagged = [name for name, item in samples.items() if isinstance(item, dict) and item.get("status") not in (None, "CONSISTENT")]
    detail = {
        "CONSISTENT": f"Reads are consistent with the declared `{strand.get('declared')}` for every sample.",
        "WARNING": f"Reads contradict the declared `{strand.get('declared')}` for: {', '.join(flagged) or 'at least one sample'}. Counts used the declaration unchanged; check the library protocol.",
        "ADVISORY": "Declared unstranded, but the reads look stranded. Unstranded counting is valid; check the library protocol.",
        "INCONCLUSIVE": f"No clear strand evidence for: {', '.join(flagged) or 'at least one sample'}; no judgement was made.",
    }.get(status, str(strand.get("reason") or "not evaluated on this route"))
    return QcItem("Strandedness check", STRAND_LEVELS.get(status, "NOT_EVALUATED"), status, detail)


def _l1_item(run_dir: Path) -> QcItem:
    summary = read_json(run_dir / "downstream" / "l1" / "backend_summary.json")
    if summary is None:
        return QcItem("L1 sample QC flags", "NOT_EVALUATED", "not produced", "downstream/l1/backend_summary.json was not written.")
    flags = [str(flag) for flag in summary.get("flags") or []]
    if not flags:
        return QcItem("L1 sample QC flags", "PASS", "none", "No sample was flagged by the L1 QC rules.")
    return QcItem("L1 sample QC flags", "WARNING", f"{len(flags)} flag(s)", "; ".join(flags))


def _contrast_items(run_dir: Path) -> list[QcItem]:
    root = run_dir / "downstream" / "l2" / "contrasts"
    if not root.is_dir():
        return []
    items = []
    for directory in sorted(path for path in root.iterdir() if path.is_dir()):
        metadata = read_json(directory / "contrast_metadata.json")
        name = directory.name
        if metadata is None:
            items.append(QcItem(f"DE warnings: {name}", "NOT_EVALUATED", "not produced", "contrast_metadata.json was not written."))
            continue
        record = metadata.get("r_warnings")
        if not isinstance(record, dict):
            items.append(QcItem(f"DE warnings: {name}", "NOT_EVALUATED", "not recorded", "R warnings are recorded from nf-rna v1.6.1."))
        else:
            warnings = [item for part in ("fit", "contrast") for item in (record.get(part) or []) if isinstance(item, dict)]
            if warnings:
                total = sum(int(item.get("count") or 1) for item in warnings)
                steps = sorted({str(item.get("step")) for item in warnings})
                items.append(QcItem(f"DE warnings: {name}", "WARNING", f"{record.get('status') or 'WARNINGS'}: {total}",
                                    f"{total} R warning(s) during {', '.join(steps)}. The analysis completed; see Results → Differential expression."))
            else:
                items.append(QcItem(f"DE warnings: {name}", "PASS", str(record.get("status") or "none"), "No R warning during the fit, results() or shrinkage."))
        shrinkage = metadata.get("shrinkage") if isinstance(metadata.get("shrinkage"), dict) else None
        if shrinkage is not None and shrinkage.get("status") != "success":
            items.append(QcItem(f"LFC shrinkage: {name}", "WARNING", str(shrinkage.get("status")),
                                f"{shrinkage.get('method')}: {shrinkage.get('error') or 'no error recorded'}. log2FoldChange_shrunk is not a valid shrunken estimate."))
    return items


def run_qc(run_dir: Path) -> RunQc:
    """The recorded scientific QC verdicts of one run, with their roll-up."""

    contract = read_json(run_dir / "frozen" / "downstream_contract.json") or {}
    project = read_yaml(run_dir / "frozen" / "project.yaml") or {}
    input_type = (project.get("input") or {}).get("type")
    items: list[QcItem] = []
    strand = _strandedness_item(run_dir, contract, input_type)
    if strand is not None:
        items.append(strand)
    if (run_dir / "downstream").is_dir():
        items.append(_l1_item(run_dir))
        items.extend(_contrast_items(run_dir))
    return RunQc(roll_up([item.level for item in items]), tuple(items))


# --------------------------------------------------------------------- runs of a project

@dataclass(frozen=True)
class RunRow:
    run_dir: Path
    case_id: str
    run_id: str
    shown: str
    execution: str
    note: str | None
    level: str | None
    started: str
    completed: str
    qc: RunQc
    has_results: bool
    has_report: bool

    @property
    def label(self) -> str:
        return f"{self.case_id}/{self.run_id}"


def project_runs(project: Path, cache: _Cache | None = None) -> list[RunRow]:
    """Every recorded run of the project, newest first, with execution status and QC roll-up."""

    cache = cache or _Cache()
    rows = []
    for run_dir, state in discover_runs(project):
        shown, note = run_display_state(run_dir, state, cache)
        contract = read_json(run_dir / "frozen" / "downstream_contract.json") or {}
        has_report = any((run_dir / "delivery").glob("report_*.html")) or (run_dir / "downstream" / "report" / "report.html").is_file()
        rows.append(RunRow(
            run_dir, run_dir.parent.name, run_dir.name, shown, execution_status(shown), note, contract.get("analysis_level"),
            str(state.get("started_at") or ""), str(state.get("completed_at") or ""), run_qc(run_dir),
            (run_dir / "downstream").is_dir(), has_report,
        ))
    return rows


# --------------------------------------------------------------------- project summary

def _rows(path: Path) -> list[dict[str, str]]:
    try:
        with path.open(encoding="utf-8", newline="") as handle:
            return [row for row in csv.DictReader(handle) if any((value or "").strip() for value in row.values() if isinstance(value, str))]
    except (OSError, UnicodeError, csv.Error):
        return []


def quantification_label(config: ProjectConfig) -> str:
    if config.input.type is InputType.RAW_COUNTS:
        return "none (imported raw counts)"
    method = (config.upstream.quantification.method if config.upstream.quantification else None) or config.upstream.quantification_method or "salmon"
    return METHOD_LABELS.get(method, method)


def reference_label(config: ProjectConfig) -> str:
    if config.input.type is InputType.RAW_COUNTS:
        return "not used (counts imported)"
    reference = config.reference
    if reference.source == "local":
        return f"managed: {reference.root}"
    if reference.source == "igenomes":
        return f"iGenomes {reference.genome or '(genome not set)'}"
    parts = [value for value in (reference.fasta, reference.gtf, reference.salmon_index, reference.hisat2_index) if value]
    return "custom: " + (", ".join(parts) or "(no files set)")


@dataclass(frozen=True)
class ProjectSummary:
    path: Path
    project_id: str
    species: str = "?"
    input_type: str = "?"
    layout: str | None = None
    strandedness: str | None = None
    quantification: str = "?"
    reference: str = "?"
    scope: str = "?"
    design: str = ""
    samples: int = 0
    contrasts: tuple[str, ...] = ()
    enrichment: str = ""
    thresholds: str = ""
    runs: tuple[RunRow, ...] = ()
    error: str | None = None

    @property
    def latest(self) -> RunRow | None:
        return self.runs[0] if self.runs else None

    def analysis_settings(self) -> str:
        if self.scope == Preset.QC.value:
            return "QC only (no statistical design)"
        text = f"{self.scope}: {self.design}"
        if self.scope == Preset.L2.value:
            text += f"; DEG {self.thresholds}; enrichment {self.enrichment}"
        return text


def project_summary(path: Path, cache: _Cache | None = None) -> ProjectSummary:
    try:
        loaded = load_project(path)
    except ProjectConfigError as exc:
        return ProjectSummary(path, path.name, error=str(exc))
    config = loaded.config
    metadata = _rows(loaded.metadata_path)
    contrasts = tuple(row.get("contrast_id", "") for row in _rows(loaded.contrasts_path) if row.get("contrast_id"))
    enrichment = tuple(config.analysis.enrichment) if config.analysis else ()
    return ProjectSummary(
        path=path, project_id=config.project.id, species=config.organism.species.value, input_type=config.input.type.value,
        layout=config.input.layout.value if config.input.layout else None, strandedness=config.upstream.strandedness,
        quantification=quantification_label(config), reference=reference_label(config), scope=config.project.preset.value,
        design=config.design.formula, samples=sum(1 for row in metadata if (row.get("sample_id") or "").strip()),
        contrasts=contrasts, enrichment=describe_enrichment(enrichment) if enrichment else "none",
        thresholds=f"padj < {config.thresholds.padj}, |log2FC| ≥ {config.thresholds.abs_log2fc}",
        runs=tuple(project_runs(path, cache)),
    )


def input_origin(path: Path) -> list[str]:
    """Where imported FASTQ files came from (the importer's own copy of the samplesheet)."""

    sheet = path / "planning" / "imported_fastq_samplesheet.csv"
    folders: list[str] = []
    for row in _rows(sheet):
        for column in ("fastq_1", "fastq_2"):
            value = (row.get(column) or "").strip()
            if value:
                folder = str(Path(value).parent)
                if folder not in folders:
                    folders.append(folder)
    return folders


# --------------------------------------------------------------------- workflow steps

STEP_NAMES = ("Dataset", "Project", "Metadata", "Contrasts", "Validate", "Run", "Results")
_INPUT = re.compile(r"count|fastq|gene|orphan|unexpected_r2|lane")
_METADATA = re.compile(r"metadata|sample|design|pair|batch|formula|variable|level|continuous|confounded|rank_deficient|residual|replication")
_CONTRAST = re.compile(r"contrast|condition_order")


def issue_area(issue: ValidationIssue) -> str:
    code = issue.code
    if _CONTRAST.search(code):
        return "Contrasts"
    if code.endswith("_only_samples") or code == "missing_count_samples":
        return "Metadata"
    if _INPUT.search(code):
        return "Dataset"
    if _METADATA.search(code):
        return "Metadata"
    return "Project"


@dataclass(frozen=True)
class Step:
    name: str
    status: str  # done, todo, problem, not needed
    detail: str
    issues: tuple[ValidationIssue, ...] = field(default_factory=tuple)


def workflow_steps(report: ValidationReport, runs: list[RunRow]) -> list[Step]:
    """Dataset → Project → Metadata → Contrasts → Validate → Run → Results for one project."""

    assert report.config is not None and report.loaded is not None
    config = report.config
    by_area: dict[str, list[ValidationIssue]] = {}
    for issue in report.issues:
        by_area.setdefault(issue_area(issue), []).append(issue)

    def errors(area: str) -> tuple[ValidationIssue, ...]:
        return tuple(issue for issue in by_area.get(area, []) if issue.severity.value == "ERROR")

    steps: list[Step] = []
    if report.counts is not None:
        steps.append(Step("Dataset", "done" if not errors("Dataset") else "problem", f"{len(report.counts.sample_ids)} samples × {report.counts.gene_count} genes", errors("Dataset")))
    elif report.fastq is not None and report.fastq.records:
        steps.append(Step("Dataset", "done" if not errors("Dataset") else "problem", f"{len(report.fastq.sample_ids)} FASTQ samples", errors("Dataset")))
    else:
        steps.append(Step("Dataset", "todo", "no input data yet", errors("Dataset")))
    steps.append(Step("Project", "problem" if errors("Project") else "done", f"{config.project.preset.value}, {config.organism.species.value}", errors("Project")))
    if config.project.preset is Preset.QC:
        steps.append(Step("Metadata", "not needed", "QC-only project"))
        steps.append(Step("Contrasts", "not needed", "QC-only project"))
    else:
        rows = report.metadata.rows if report.metadata else ()
        meta_errors = errors("Metadata")
        steps.append(Step("Metadata", "problem" if meta_errors else ("done" if rows else "todo"), f"{len(rows)} samples described" if rows else "metadata.csv has no rows", meta_errors))
        contrasts = report.contrasts.contrasts if report.contrasts else ()
        if config.project.preset is Preset.L1:
            steps.append(Step("Contrasts", "not needed", "L1 has no contrasts (saving contrasts moves the project to L2)"))
        else:
            con_errors = errors("Contrasts")
            steps.append(Step("Contrasts", "problem" if con_errors else ("done" if contrasts else "todo"), f"{len(contrasts)} contrast(s)" if contrasts else "no contrast yet", con_errors))
    steps.append(Step("Validate", "done" if report.is_valid else "problem", "valid" + (f", {len(report.warnings)} warning(s)" if report.warnings else "") if report.is_valid else f"{len(report.errors)} error(s)", tuple(report.errors)))
    if runs:
        latest = runs[0]
        steps.append(Step("Run", "done" if latest.execution == "SUCCESS" else ("todo" if latest.execution == "RUNNING" else "problem"), f"latest run {latest.execution}"))
    else:
        steps.append(Step("Run", "todo", "no run yet"))
    with_results = [run for run in runs if run.has_results]
    steps.append(Step("Results", "done" if with_results else "todo", f"{len(with_results)} run(s) with results" if with_results else "no results yet"))
    return steps


def next_step(steps: list[Step]) -> Step | None:
    for step in steps:
        if step.status in {"problem", "todo"}:
            return step
    return None


# --------------------------------------------------------------------- upstream reuse

@dataclass(frozen=True)
class ReuseOption:
    reference: str
    accepted: bool
    reason: str


def reuse_options(report: ValidationReport) -> list[ReuseOption]:
    """SUCCESS runs of a FASTQ project and whether ``--reuse-upstream``'s own preflight accepts them."""

    from rnaseq.errors import ExecutionPreflightError
    from rnaseq.service import preflight_reuse_upstream

    if report.config is None or report.config.input.type is not InputType.FASTQ:
        return []
    options = []
    for run_dir, state in discover_runs(report.project_dir):
        if state.get("status") != "SUCCESS":
            continue
        reference = f"{run_dir.parent.name}/{run_dir.name}"
        try:
            preflight_reuse_upstream(report, reference)
        except ExecutionPreflightError as exc:
            options.append(ReuseOption(reference, False, str(exc)))
            continue
        made_by = _frozen_pipeline_version(run_dir)
        if made_by != PIPELINE_VERSION:
            # The run-time check compares the upstream identity, which includes the nf-rna version
            # (service._upstream_identity), so this run would be refused after it had started.
            options.append(ReuseOption(reference, False, f"made by nf-rna {made_by or 'unknown'}; --reuse-upstream needs a source run made by this "
                                                          f"version ({PIPELINE_VERSION}), because the frozen upstream identity includes the nf-rna version"))
            continue
        options.append(ReuseOption(reference, True, "accepted by the reuse preflight; the frozen upstream inputs are compared byte for byte when the run starts"))
    return options


def _frozen_pipeline_version(run_dir: Path) -> str | None:
    manifest = read_yaml(run_dir / "frozen" / "input_manifest.yaml")
    pipeline = manifest.get("pipeline") if isinstance(manifest, dict) else None
    version = pipeline.get("version") if isinstance(pipeline, dict) else None
    return str(version) if version is not None else None


# --------------------------------------------------------------------- logs

_LOG_ALERT = re.compile(r"\b(ERROR|WARN(ING)?|FAILED|Exception|Error)\b")


def alert_lines(text: str) -> str:
    """Only the lines of a log that report a warning or an error."""

    return "\n".join(line for line in text.splitlines() if _LOG_ALERT.search(line))
