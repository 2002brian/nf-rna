"""Read-only access to a run's finished outputs for the Results and Reports pages (no Streamlit here).

Everything shown comes from files the pipeline wrote: figures are the R
pipeline's PNG/PDF files, tables are read as text so values are displayed
exactly as written, and significance is the backend's own significant table or
``significant`` column. Searching, filtering and sorting only select rows.
"""

from __future__ import annotations

import csv
import getpass
import json
import os
import shlex
import socket
import zipfile
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import yaml

from rnaseq.gui.config import gui_state_dir
from rnaseq.gui.jobs import run_display_state
from rnaseq.models import normalize_enrichment_selection
from rnaseq.run_status import _Cache, discover_runs
from rnaseq.service import DELIVERY_MANIFEST_FILENAME, verify_delivery_manifest
from rnaseq.workflow_support import _figure_reason, _l2_versions

ONTOLOGIES = ("BP", "MF", "CC")
ORA_FOREGROUNDS = (("significant", "all significant DEGs"), ("up", "upregulated DEGs"), ("down", "downregulated DEGs"))
LEGACY_SIGNIFICANCE = {
    "ora": "pvalue <= 0.05 and qvalue <= 0.2 (p.adjust not used)",
    "gsea": "p.adjust <= 0.05 and pvalue <= 0.05",
}


def read_json(path: Path) -> dict[str, Any] | None:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return None
    return value if isinstance(value, dict) else None


def read_yaml(path: Path) -> dict[str, Any] | None:
    try:
        value = yaml.safe_load(path.read_text(encoding="utf-8"))
    except (OSError, yaml.YAMLError):
        return None
    return value if isinstance(value, dict) else None


@dataclass(frozen=True)
class Table:
    """A TSV as text: every cell is the string written by the pipeline."""

    columns: tuple[str, ...]
    rows: tuple[dict[str, str], ...]

    def column(self, name: str) -> list[str]:
        return [row.get(name, "") for row in self.rows]

    def where(self, column: str, value: str) -> "Table":
        return Table(self.columns, tuple(row for row in self.rows if row.get(column) == value))


def read_table(path: Path) -> Table:
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows = tuple(reader)
        return Table(tuple(reader.fieldnames or ()), rows)


def _number(value: str) -> float | None:
    try:
        return float(value)
    except ValueError:
        return None  # NA and empty cells


def tsv_rows(path: Path) -> int | None:
    if not path.is_file():
        return None
    with path.open(encoding="utf-8", newline="") as handle:
        return sum(1 for _ in csv.DictReader(handle, delimiter="\t"))


# --------------------------------------------------------------------- runs

@dataclass(frozen=True)
class ResultRun:
    run_dir: Path
    case_id: str
    run_id: str
    state: str
    level: str | None

    @property
    def label(self) -> str:
        return f"{self.case_id}/{self.run_id}"

    @property
    def l1(self) -> Path:
        return self.run_dir / "downstream" / "l1"

    @property
    def l2(self) -> Path:
        return self.run_dir / "downstream" / "l2"

    @property
    def delivery(self) -> Path:
        return self.run_dir / "delivery"

    @property
    def contract(self) -> dict[str, Any]:
        return read_json(self.run_dir / "frozen" / "downstream_contract.json") or {}

    @property
    def project(self) -> dict[str, Any]:
        return read_yaml(self.run_dir / "frozen" / "project.yaml") or {}

    @property
    def is_l2(self) -> bool:
        return self.level == "L2"


def result_runs(project: Path, cache: _Cache | None = None) -> list[ResultRun]:
    """Runs with downstream outputs, newest first (any state: a failed run shows what it produced)."""

    found = []
    for run_dir, state in discover_runs(project):
        if not (run_dir / "downstream").is_dir():
            continue
        shown, _note = run_display_state(run_dir, state, cache)
        contract = read_json(run_dir / "frozen" / "downstream_contract.json") or {}
        found.append(ResultRun(run_dir, run_dir.parent.name, run_dir.name, shown, contract.get("analysis_level")))
    return found


# --------------------------------------------------------------------- figures

@dataclass(frozen=True)
class Figure:
    label: str
    png: Path | None
    pdf: Path | None
    reason: str | None = None


def figure(directory: Path, stem: str, label: str, reason: str) -> Figure:
    png, pdf = directory / f"{stem}.png", directory / f"{stem}.pdf"
    return Figure(label, png if png.is_file() else None, pdf if pdf.is_file() else None, None if png.is_file() else reason)


def running_score_figures(directory: Path, record: Any, label: str) -> list[Figure]:
    """Every running-score plot the backend drew, or one unavailable entry per NES sign with its reason."""

    running = (record.get("figures") or {}).get("running_score") if isinstance(record, dict) else None
    running = running if isinstance(running, dict) else {}
    figures = []
    for sign in ("positive", "negative"):
        found = sorted((path for path in directory.glob(f"running_score/{sign}_*.png") if not path.stem.endswith("_300dpi")), key=lambda path: int(path.name.split("_")[1]))
        for path in found:
            _sign, rank, *rest = path.stem.split("_")
            figures.append(Figure(f"{label} running score, {sign} NES, rank {rank}: {'_'.join(rest).replace('GO_', 'GO:')}", path, path.with_suffix(".pdf") if path.with_suffix(".pdf").is_file() else None))
        if not found:
            if running.get("reason"):
                reason = str(running["reason"])
            elif running:
                reason = f"no significant gene set with {sign} NES (p.adjust <= 0.05)"
            else:
                reason = "no running-score plots in these results (produced from nf-rna v1.5)"
            figures.append(Figure(f"{label} running score, {sign} NES", None, None, reason))
    return figures


# --------------------------------------------------------------------- L1

@dataclass(frozen=True)
class QcView:
    summary: dict[str, Any] | None
    tables: dict[str, Path]
    figures: list[Figure]
    multiqc: Path | None
    upstream_qc: dict[str, Any] | None = None

    @property
    def strandedness(self) -> dict[str, Any]:
        value = (self.upstream_qc or {}).get("strandedness")
        return value if isinstance(value, dict) else {}

    @property
    def execution_threads(self) -> dict[str, Any]:
        value = (self.upstream_qc or {}).get("execution_threads")
        return value if isinstance(value, dict) else {}


def qc_view(run: ResultRun) -> QcView:
    l1 = run.l1
    variables = (run.project.get("design") or {}).get("variables") or {}
    figures = [
        figure(l1, "pca", "PCA (blind VST)", "pca.png was not produced"),
        *([figure(l1, "pca_by_batch", "PCA coloured by batch", "pca_by_batch.png was not produced")] if "batch" in variables or (l1 / "pca_by_batch.png").is_file() else []),
        figure(l1, "sample_correlation", "Sample correlation (blind VST)", "sample_correlation.png was not produced"),
        figure(l1, "boxplot_log2_raw_counts", "Per-sample log2(raw count + 1), before normalization", "boxplot_log2_raw_counts.png was not produced"),
        figure(l1, "boxplot_vst", "Per-sample blind VST expression", "boxplot_vst.png was not produced"),
    ]
    tables = {name: l1 / name for name in ("library_size_qc.tsv", "sample_distribution_qc.tsv", "pca_variance.tsv", "sample_correlation.tsv") if (l1 / name).is_file()}
    multiqc = next(iter(sorted(run.delivery.glob("multiqc/*.html"))), None)
    upstream_qc = run.contract.get("upstream_qc")
    return QcView(read_json(l1 / "backend_summary.json"), tables, figures, multiqc, upstream_qc if isinstance(upstream_qc, dict) else None)


# --------------------------------------------------------------------- L2 / DEG

def contrast_ids(run: ResultRun) -> list[str]:
    """Contrasts in the frozen contrasts.csv order that have L2 outputs."""

    frozen = run.run_dir / "frozen" / "contrasts.csv"
    ordered: list[str] = []
    if frozen.is_file():
        with frozen.open(encoding="utf-8", newline="") as handle:
            ordered = [row["contrast_id"] for row in csv.DictReader(handle) if row.get("contrast_id")]
    present = {path.name for path in (run.l2 / "contrasts").iterdir() if path.is_dir()} if (run.l2 / "contrasts").is_dir() else set()
    return [item for item in ordered if item in present] + sorted(present - set(ordered))


@dataclass(frozen=True)
class DegView:
    contrast_id: str
    root: Path
    metadata: dict[str, Any] | None
    summary: dict[str, Any] | None
    counts: dict[str, int | None]
    padj_below: object
    abs_log2fc_at_least: object
    rule: str
    figures: list[Figure]

    @property
    def shrinkage(self) -> dict[str, Any]:
        value = (self.metadata or {}).get("shrinkage")
        return value if isinstance(value, dict) else {}

    @property
    def r_warnings(self) -> list[dict[str, Any]] | None:
        """Recorded R warnings of the shared fit and this contrast; None for runs before v1.6.1."""

        record = (self.metadata or {}).get("r_warnings")
        if not isinstance(record, dict):
            return None
        return [item for part in ("fit", "contrast") for item in (record.get(part) or []) if isinstance(item, dict)]


def deg_view(run: ResultRun, contrast_id: str) -> DegView:
    root = run.l2 / "contrasts" / contrast_id
    metadata = read_json(root / "contrast_metadata.json")
    definition = (metadata or {}).get("deg_definition") if isinstance((metadata or {}).get("deg_definition"), dict) else None
    thresholds = run.project.get("thresholds") or {}
    if definition is not None:
        padj, lfc, rule = definition.get("padj_below"), definition.get("abs_log2fc_at_least"), str(definition.get("rule") or "")
    else:
        padj, lfc = thresholds.get("padj"), thresholds.get("abs_log2fc")
        rule = "not recorded (L2 results from before nf-rna v1.5); thresholds from the run's frozen project.yaml"
    counts = {name: tsv_rows(root / f"{name}.tsv") for name in ("significant", "upregulated", "downregulated")}
    figures = [
        figure(root, "volcano", f"Volcano plot: {contrast_id}", "volcano.png was not produced"),
        figure(root, "ma_plot", f"MA plot: {contrast_id}", "ma_plot.png was not produced (L2 results from before nf-rna v1.5)"),
        figure(root, "heatmap", f"DEG heatmap: {contrast_id}", str((read_json(root / "backend_summary.json") or {}).get("heatmap_status") or "heatmap.png was not produced")),
    ]
    return DegView(contrast_id, root, metadata, read_json(root / "backend_summary.json"), counts, padj, lfc, rule, figures)


DEG_SORT_COLUMNS = ("padj", "pvalue", "log2FoldChange", "log2FoldChange_shrunk", "baseMean", "stat")


def significant_ids(view: DegView) -> set[str]:
    path = view.root / "significant.tsv"
    return set(read_table(path).column("gene_id")) if path.is_file() else set()


def filter_deg(
    table: Table, significant: set[str], *, query: str = "", only_significant: bool = False,
    padj_max: float | None = None, abs_lfc_min: float | None = None, sort_by: str = "padj", descending: bool = False,
) -> Table:
    """Select and order rows of all_genes.tsv; the cells stay the strings of the file."""

    text = query.strip().lower()

    def keep(row: dict[str, str]) -> bool:
        if text and text not in row.get("gene_id", "").lower() and text not in row.get("symbol", "").lower():
            return False
        if only_significant and row.get("gene_id") not in significant:
            return False
        if padj_max is not None:
            padj = _number(row.get("padj", ""))
            if padj is None or padj > padj_max:
                return False
        if abs_lfc_min is not None:
            lfc = _number(row.get("log2FoldChange", ""))
            if lfc is None or abs(lfc) < abs_lfc_min:
                return False
        return True

    rows = [row for row in table.rows if keep(row)]
    if sort_by in table.columns:
        numbered = [(_number(row.get(sort_by, "")), row) for row in rows]
        present = sorted((item for item in numbered if item[0] is not None), key=lambda item: item[0], reverse=descending)
        rows = [row for _value, row in present] + [row for value, row in numbered if value is None]  # NA last, stable
    return Table(table.columns, tuple(rows))


# --------------------------------------------------------------------- enrichment

@dataclass(frozen=True)
class EnrichmentResult:
    """One result set: a foreground x ontology of ORA, an ontology of GO GSEA, or KEGG GSEA."""

    label: str
    directory: Path
    all_terms: Path | None
    significant: Path | None
    members: Path | None
    members_legacy: bool
    status: str
    reason: str | None
    facts: dict[str, object]
    figures: list[Figure]


@dataclass(frozen=True)
class EnrichmentModule:
    key: str
    label: str
    root: Path
    summary: dict[str, Any] | None
    reason: str | None = None
    facts: dict[str, object] = field(default_factory=dict)

    @property
    def legacy(self) -> bool:
        return self.summary is not None and self.summary.get("significance_semantics") != "v1.5"


# (key, label, directory, summary file, selection name)
MODULES = (
    ("go", "GO ORA", "go", "go_backend_summary.json", "go"),
    ("kegg", "KEGG ORA", "kegg", "kegg_backend_summary.json", "kegg"),
    ("gsea_go", "GO GSEA", "gsea_go", "gsea_backend_summary.json", "gsea"),
    ("gsea_kegg", "KEGG GSEA", "gsea_kegg", "gsea_kegg_backend_summary.json", "gsea"),
)


def _existing(*paths: Path) -> tuple[Path | None, bool]:
    """The first existing path, and whether it was a fallback (legacy) name."""

    for index, path in enumerate(paths):
        if path.is_file():
            return path, index > 0
    return None, False


def enrichment_modules(run: ResultRun) -> list[EnrichmentModule]:
    """Selected enrichment analyses of the run, with their significance and resource facts."""

    selected = normalize_enrichment_selection((run.contract.get("analysis") or {}).get("enrichment") or [])
    pinned = ((((run.project.get("annotation") or {}).get("enrichment") or {}).get("kegg") or {}).get("snapshot"))
    modules = []
    for key, label, directory, summary_name, selection in MODULES:
        if selection not in selected:
            continue
        root = run.l2 / "enrichment" / directory
        summary = read_json(root / summary_name)
        if summary is None:
            modules.append(EnrichmentModule(key, label, root, None, f"{summary_name} was not produced (the backend did not finish)"))
            continue
        kind = "gsea" if key.startswith("gsea") else "ora"
        facts: dict[str, object] = {"status": summary.get("status"), "reason": summary.get("reason")}
        if summary.get("significance_semantics") == "v1.5":
            significance = summary.get("significance") or {}
            facts["significance semantics"] = "v1.5"
            facts["significance rule"] = significance.get("rule")
            facts["p-value adjustment"] = (significance.get("clusterprofiler_call") or {}).get("pAdjustMethod")
        else:
            facts["significance semantics"] = "legacy (before v1.5; not recorded in the summary)"
            facts["significance rule"] = f"as written by that nf-rna version: {LEGACY_SIGNIFICANCE[kind]}"
        if kind == "ora":
            facts["ORA universe"] = (summary.get("ora") or {}).get("universe")
            facts["gene-set size"] = f"{(summary.get('ora') or {}).get('min_gs_size')}..{(summary.get('ora') or {}).get('max_gs_size')}"
        else:
            gsea = summary.get("gsea") or {}
            facts["GSEA ranking metric"] = gsea.get("ranking_metric")
            facts["GSEA seed"] = gsea.get("seed")
            facts["nPermSimple / eps"] = f"{gsea.get('nPermSimple')} / {gsea.get('eps')}" if gsea else None
        source = summary.get("gene_set_source") or summary.get("resource") or {}
        if key in {"kegg", "gsea_kegg"}:
            facts["KEGG snapshot"] = source.get("snapshot_id")
            facts["KEGG release"] = source.get("kegg_release")
            facts["KEGG snapshot pinned"] = f"yes ({pinned})" if pinned else "no (latest local snapshot at run time)"
        else:
            facts["GO annotation"] = f"{source.get('orgdb')} {source.get('orgdb_version')}, GO.db {source.get('go_db_version')}" if source.get("orgdb") else None
        modules.append(EnrichmentModule(key, label, root, summary, facts={name: value for name, value in facts.items() if value not in (None, "")}))
    return modules


def _contrast_record(summary: dict[str, Any] | None, contrast_id: str) -> dict[str, Any]:
    for item in (summary or {}).get("contrasts", []) or []:
        if isinstance(item, dict) and item.get("contrast_id") == contrast_id:
            return item
    return {}


def _status(record: dict[str, Any], fallback: dict[str, Any]) -> tuple[str, str | None]:
    status = record.get("status") or fallback.get("status") or "NOT_AVAILABLE"
    return str(status), (record.get("reason") or fallback.get("reason") or None)


def enrichment_results(module: EnrichmentModule, contrast_id: str) -> list[EnrichmentResult]:
    item = _contrast_record(module.summary, contrast_id)
    base = module.root / contrast_id
    results = []
    if module.key == "go":
        for foreground, words in ORA_FOREGROUNDS:
            outcome = (item.get("foregrounds") or {}).get(foreground) or {}
            for ontology in ONTOLOGIES:
                record = (outcome.get("ontologies") or {}).get(ontology) or {}
                directory = base / foreground
                status, reason = _status(record, outcome)
                members, legacy = _existing(directory / f"GO_{ontology}_members.tsv", directory / f"GO_{ontology}_gene_members.tsv")
                universe = (record.get("universe") or {})
                results.append(EnrichmentResult(
                    f"{words} · GO {ontology}", directory, _existing(directory / f"GO_{ontology}_all_terms.tsv")[0], _existing(directory / f"GO_{ontology}.tsv")[0],
                    members, legacy, status, reason,
                    {"ORA universe genes": universe.get("universe_genes"), "foreground genes in universe": universe.get("foreground_in_universe"),
                     "evaluated terms": record.get("evaluated_term_count"), "significant terms": record.get("significant_term_count")},
                    [figure(directory, f"dotplot_{ontology}", f"GO {ontology} dotplot ({words})", _figure_reason(record, "dotplot", str(reason or "not drawn"))),
                     figure(directory, f"dotplot_{ontology}_simplified", f"GO {ontology} dotplot, simplified ({words})", _figure_reason(record, "dotplot_simplified", str((record.get("simplified") or {}).get("reason") or reason or "not produced by this GO ORA backend")))],
                ))
    elif module.key == "kegg":
        for foreground, words in ORA_FOREGROUNDS:
            record = (item.get("foregrounds") or {}).get(foreground) or {}
            directory = base / foreground
            status, reason = _status(record, module.summary or {})
            members, legacy = _existing(directory / "kegg_members.tsv", directory / "kegg_gene_members.tsv")
            universe = record.get("universe") or {}
            results.append(EnrichmentResult(
                f"{words} · KEGG", directory, _existing(directory / "kegg_all_terms.tsv")[0], _existing(directory / "kegg.tsv")[0], members, legacy, status, reason,
                {"ORA universe genes": universe.get("universe_genes"), "foreground genes in universe": universe.get("foreground_in_universe"),
                 "evaluated pathways": record.get("tested_terms"), "significant pathways": record.get("significant_terms")},
                [figure(directory, "dotplot", f"KEGG dotplot ({words})", _figure_reason(record, "dotplot", str(reason or "not drawn")))],
            ))
    elif module.key == "gsea_go":
        ranking = item.get("ranking") or {}
        for ontology in ONTOLOGIES:
            record = (item.get("ontologies") or {}).get(ontology) or {}
            directory = base / ontology
            status, reason = _status(record, item)
            members, legacy = _existing(directory / "members.tsv", directory / "gsea_core_members.tsv")
            results.append(EnrichmentResult(
                f"GO {ontology}", directory, _existing(directory / "all_terms.tsv")[0], _existing(directory / "significant.tsv")[0], members, legacy, status, reason,
                {"ranked genes": ranking.get("final_ranked_genes"), "evaluated terms": record.get("evaluated_terms", record.get("all_terms")),
                 "significant terms": record.get("significant_terms"), "positive / negative NES": f"{record.get('positive_terms')} / {record.get('negative_terms')}" if record else None},
                [figure(directory, "dotplot", f"GO {ontology} GSEA dotplot", _figure_reason(record, "dotplot", str(reason or "not drawn"))),
                 figure(directory, "dotplot_simplified", f"GO {ontology} GSEA dotplot, simplified", _figure_reason(record, "dotplot_simplified", str((record.get("simplified") or {}).get("reason") or reason or "not produced by this GO GSEA backend"))),
                 *running_score_figures(directory, record, f"GO {ontology}")],
            ))
    else:
        status, reason = _status(item, module.summary or {})
        members, legacy = _existing(base / "members.tsv", base / "core_members.tsv")
        ranking = (_contrast_record(read_json(module.root.parent / "gsea_go" / "gsea_backend_summary.json"), contrast_id).get("ranking") or {})
        results.append(EnrichmentResult(
            "KEGG", base, _existing(base / "all_terms.tsv")[0], _existing(base / "significant.tsv")[0], members, legacy, status, reason,
            {"ranked genes": ranking.get("final_ranked_genes"), "evaluated pathways": item.get("evaluated_terms", item.get("all_terms")),
             "significant pathways": item.get("significant_terms"), "positive / negative NES": f"{item.get('positive_terms')} / {item.get('negative_terms')}" if item else None},
            [figure(base, "dotplot", "KEGG GSEA dotplot", _figure_reason(item, "dotplot", str(reason or "not drawn"))), *running_score_figures(base, item, "KEGG")],
        ))
    return results


MEMBER_ID_COLUMNS = ("pathway_id", "GO_ID", "ID")


def term_id_column(table: Table) -> str | None:
    return next((column for column in MEMBER_ID_COLUMNS if column in table.columns), None)


def pathway_members(path: Path, pathway_id: str) -> Table:
    table = read_table(path)
    column = term_id_column(table)
    return table.where(column, pathway_id) if column else Table(table.columns, ())


# --------------------------------------------------------------------- provenance

def provenance(run: ResultRun) -> list[tuple[str, str]]:
    """nf-rna build, lock, reference and package versions, as recorded by the run."""

    record = read_yaml(run.run_dir / "provenance" / "run_provenance.yaml") or {}
    runtime = record.get("downstream_runtime") or {}
    project = run.project
    reference = project.get("reference") or {}
    rows: list[tuple[str, object]] = [
        ("nf-rna version", record.get("pipeline_version")),
        ("git commit", record.get("git_commit")),
        ("Nextflow", record.get("nextflow_version")),
        ("Python", record.get("python_version")),
        ("downstream lock", f"{(runtime.get('lock') or {}).get('filename')} sha256 {(runtime.get('lock') or {}).get('sha256')}" if runtime.get("lock") else None),
        ("runtime prefix", runtime.get("prefix")),
        ("input type", (project.get("input") or {}).get("type")),
        ("reference", ", ".join(f"{name}: {value}" for name, value in reference.items() if value) or None),
        ("upstream", record.get("upstream_implementation") or record.get("fastq_backend")),
        ("strandedness check", (record.get("strandedness_qc") or {}).get("status")),
    ]
    threads = record.get("execution_threads") if isinstance(record.get("execution_threads"), dict) else {}
    for name, item in (threads.get("processes") or {}).items():
        if isinstance(item, dict):
            rows.append((f"threads: {name}", f"{', '.join(map(str, item.get('cpus') or [])) or 'not recorded'} CPU(s) per task, {item.get('tasks')} task(s)"))
    for sample, value in (threads.get("salmon_quant_threads_cmd_info") or {}).items():
        rows.append((f"Salmon --threads: {sample}", value))
    if threads and threads.get("status") != "RECORDED":
        rows.append(("upstream threads", f"not recorded: {threads.get('reason')}"))
    l1_summary = read_json(run.l1 / "backend_summary.json") or {}
    if run.is_l2:
        metadata = {item: read_json(run.l2 / "contrasts" / item / "contrast_metadata.json") for item in contrast_ids(run)}
        indexed = {path.name: path for path in (run.l2 / "enrichment").iterdir() if path.is_dir()} if (run.l2 / "enrichment").is_dir() else {}
        versions, snapshot = _l2_versions(l1_summary.get("package_versions") or {}, metadata, indexed)
        rows.extend(versions.items())
        if snapshot:
            rows.append(("KEGG snapshot", snapshot))
    else:
        rows.extend((name, value) for name, value in (l1_summary.get("package_versions") or {}).items())
    return [(name, str(value)) for name, value in rows if value]


# --------------------------------------------------------------------- downloads

@dataclass(frozen=True)
class Download:
    group: str
    label: str
    path: Path


def report_path(run: ResultRun) -> Path | None:
    delivered = sorted(run.delivery.glob("report_*.html"))
    if delivered:
        return delivered[-1]
    built = run.run_dir / "downstream" / "report" / "report.html"
    return built if built.is_file() else None


def downloads(run: ResultRun) -> list[Download]:
    """Single files offered for download, grouped; only files that exist, and no L2 file for an L1 run."""

    items: list[Download] = []
    report = report_path(run)
    if report is not None:
        items.append(Download("Report", report.name, report))
    normalized = run.l2 / "normalized_counts.tsv" if run.is_l2 and (run.l2 / "normalized_counts.tsv").is_file() else run.l1 / "normalized_counts.tsv"
    for path in (normalized, run.l1 / "vst.tsv"):
        if path.is_file():
            items.append(Download("Counts", f"{path.parent.name}/{path.name}", path))
    for path in sorted(run.l1.glob("*.tsv")):
        if path.name not in {"normalized_counts.tsv", "vst.tsv"}:
            items.append(Download("QC tables", f"l1/{path.name}", path))
    if run.is_l2:
        for contrast in contrast_ids(run):
            for name in ("all_genes", "significant", "upregulated", "downregulated"):
                path = run.l2 / "contrasts" / contrast / f"{name}.tsv"
                if path.is_file():
                    items.append(Download("DEG tables", f"{contrast}/{name}.tsv", path))
        enrichment = run.l2 / "enrichment"
        for path in sorted(enrichment.rglob("*.tsv")) if enrichment.is_dir() else ():
            items.append(Download("Enrichment tables and members", str(path.relative_to(enrichment)), path))
    roots = [run.l1, *([run.l2] if run.is_l2 else [])]
    for root in roots:
        for path in sorted(root.rglob("*")):
            if path.suffix in {".png", ".pdf"} and not path.stem.endswith("_300dpi"):
                items.append(Download("Figures (PNG / PDF)", str(path.relative_to(run.run_dir / "downstream")), path))
    return items


def too_large(path: Path, limit_mb: int) -> bool:
    return path.stat().st_size > limit_mb * 1024 * 1024


def scp_command(path: Path) -> str:
    """The command to copy a server file or folder to the current directory of another machine."""

    flag = "-r " if path.is_dir() else ""
    return f"scp {flag}{getpass.getuser()}@{socket.gethostname()}:{shlex.quote(str(path))} ."


def delivery_size(run: ResultRun) -> int | None:
    manifest = read_yaml(run.delivery / DELIVERY_MANIFEST_FILENAME)
    if manifest is None:
        return None
    return sum(int(entry.get("size_bytes") or 0) for entry in manifest.get("files") or [] if isinstance(entry, dict))


def delivery_archive_path(run: ResultRun) -> Path:
    return gui_state_dir() / "packages" / f"{run.case_id}_{run.run_id.replace('+', '_')}_delivery.zip"


def build_delivery_archive(run: ResultRun) -> Path:
    """A zip of the run's verified delivery folder (written once to the GUI cache, file by file).

    The delivery folder is checked against its own ``delivery_manifest.yaml``
    first; the run directory is never written.
    """

    failures = verify_delivery_manifest(run.delivery)
    if failures:
        raise OSError("The delivery package does not match its manifest: " + "; ".join(failures[:5]))
    target = delivery_archive_path(run)
    manifest = run.delivery / DELIVERY_MANIFEST_FILENAME
    if target.is_file() and target.stat().st_mtime >= manifest.stat().st_mtime:
        return target
    entries = [entry["relative_path"] for entry in (read_yaml(manifest) or {}).get("files") or []]
    target.parent.mkdir(parents=True, exist_ok=True)
    partial = target.with_suffix(f".{os.getpid()}.partial")
    prefix = f"{run.case_id}_{run.run_id.replace('+', '_')}_delivery"
    with zipfile.ZipFile(partial, "w", compression=zipfile.ZIP_DEFLATED, allowZip64=True) as archive:
        for relative in [DELIVERY_MANIFEST_FILENAME, *entries]:
            archive.write(run.delivery / relative, f"{prefix}/{relative}")
    partial.replace(target)
    return target
