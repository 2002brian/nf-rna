"""Milestone 4A orchestration for simple-design DESeq2 inference only."""

from __future__ import annotations

import csv
import json
import subprocess
from dataclasses import dataclass
from datetime import UTC, datetime
from importlib import resources
from pathlib import Path
from typing import Any

import yaml

from rnaseq.downstream import L1_FILTER, PreparedL1, _read_yaml, _require_r, _write, execute_l1, prepare_l1
from rnaseq.errors import DownstreamExecutionError
from rnaseq.models import PIPELINE_VERSION, DesignType, InputType, ProjectConfig
from rnaseq.validators import CONTRAST_HEADER, ContrastDefinition, ValidationReport

RESULT_HEADER = ("gene_id", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")
L2_STATES = {"RUNNING", "SUCCESS", "PARTIAL", "FAILED"}
HEATMAP_TOP_N = 50


@dataclass(frozen=True)
class PreparedL2:
    l1: PreparedL1
    output_dir: Path
    config: ProjectConfig
    contrasts: tuple[ContrastDefinition, ...]
    metadata_rows: tuple[dict[str, str], ...]


@dataclass(frozen=True)
class L2Result:
    output_dir: Path
    state_path: Path
    status: str
    summaries: tuple[dict[str, Any], ...]


def _state(path: Path, **value: Any) -> None:
    _write(path, json.dumps(value, indent=2, sort_keys=True) + "\n")


def _read_metadata(path: Path) -> tuple[tuple[str, ...], tuple[dict[str, str], ...]]:
    try:
        with path.open(encoding="utf-8", newline="") as handle:
            reader = csv.reader(handle)
            header = next(reader)
            rows = tuple(dict(zip(header, row, strict=True)) for row in reader if row and any(row))
    except (OSError, StopIteration, ValueError) as exc:
        raise DownstreamExecutionError(f"Cannot read frozen metadata for L2: {exc}") from exc
    return tuple(header), rows


def _read_contrasts(path: Path, columns: tuple[str, ...], rows: tuple[dict[str, str], ...], formula: str) -> tuple[ContrastDefinition, ...]:
    try:
        with path.open(encoding="utf-8", newline="") as handle:
            reader = csv.reader(handle)
            header = next(reader)
            values = [row for row in reader if row and any(row)]
    except (OSError, StopIteration) as exc:
        raise DownstreamExecutionError(f"Cannot read contrasts for L2: {exc}") from exc
    if header != CONTRAST_HEADER:
        raise DownstreamExecutionError("L2 requires contrasts columns exactly: " + ",".join(CONTRAST_HEADER))
    formula_variables = {part.strip() for part in formula.removeprefix("~").split("+")}
    seen: set[str] = set()
    contrasts: list[ContrastDefinition] = []
    for number, row in enumerate(values, start=2):
        if len(row) != 4 or any(not value.strip() for value in row):
            raise DownstreamExecutionError(f"L2 contrast row {number} must contain four nonblank values.")
        contrast_id, factor, numerator, denominator = row
        if contrast_id in seen:
            raise DownstreamExecutionError(f"L2 contrast_id is duplicated: {contrast_id}")
        seen.add(contrast_id)
        if factor not in columns or factor not in formula_variables:
            raise DownstreamExecutionError(f"L2 contrast factor {factor!r} is not available in the validated design.")
        levels = {row[factor] for row in rows}
        if numerator not in levels or denominator not in levels:
            raise DownstreamExecutionError(f"L2 contrast {contrast_id} has a numerator or denominator level absent from metadata.")
        if numerator == denominator:
            raise DownstreamExecutionError(f"L2 contrast {contrast_id} numerator and denominator must differ.")
        contrasts.append(ContrastDefinition(contrast_id, factor, numerator, denominator))
    if not contrasts:
        raise DownstreamExecutionError("L2 requires at least one validated contrast.")
    return tuple(contrasts)


def _guard_replicates(config: ProjectConfig, metadata_rows: tuple[dict[str, str], ...], contrasts: tuple[ContrastDefinition, ...]) -> None:
    """Block n=1 inference and recheck complete pairs before invoking DESeq2."""
    for contrast in contrasts:
        numerator_rows = [row for row in metadata_rows if row[contrast.factor] == contrast.numerator]
        denominator_rows = [row for row in metadata_rows if row[contrast.factor] == contrast.denominator]
        if len(numerator_rows) < 2 or len(denominator_rows) < 2:
            raise DownstreamExecutionError(
                "Differential-expression inference requires biological replication. "
                f"Current contrast {contrast.contrast_id}: {contrast.denominator} n={len(denominator_rows)}, "
                f"{contrast.numerator} n={len(numerator_rows)}. L1 QC remains available; L2 inference was not executed."
            )
        if config.design.type is DesignType.PAIRED_TWO_GROUP:
            pair_id = config.design.pair_id
            if pair_id is None or pair_id not in metadata_rows[0]:
                raise DownstreamExecutionError("Paired L2 inference requires the declared design.pair_id metadata column.")
            pairs: dict[str, list[str]] = {}
            for row in metadata_rows:
                if row[contrast.factor] in (contrast.numerator, contrast.denominator):
                    pairs.setdefault(row[pair_id], []).append(row[contrast.factor])
            incomplete = sorted(
                pair for pair, levels in pairs.items()
                if sorted(levels) != sorted([contrast.numerator, contrast.denominator])
            )
            if incomplete or len(pairs) < 2:
                detail = ", ".join(incomplete) if incomplete else "fewer than two complete pairs"
                raise DownstreamExecutionError(f"Paired L2 inference requires at least two complete pairs; invalid pair IDs: {detail}")


def prepare_l2(report: ValidationReport, *, run_id: str | None) -> PreparedL2:
    """Select immutable inputs and repeat critical contrast checks immediately before inference."""
    l1 = prepare_l1(report, run_id=run_id)
    assert report.config is not None and report.loaded is not None
    if l1.source_type == "raw_counts":
        config = report.config
        metadata_path = report.loaded.metadata_path
        contrasts_path = report.loaded.contrasts_path
        output = report.project_dir / "downstream" / "l2"
    else:
        run_dir = l1.output_dir.parents[1]
        raw = _read_yaml(run_dir / "frozen" / "project.yaml", "frozen project configuration")
        try:
            config = ProjectConfig.model_validate(raw)
        except Exception as exc:  # Pydantic's rendered error remains useful at the CLI boundary.
            raise DownstreamExecutionError(f"Frozen project configuration is invalid for L2: {exc}") from exc
        metadata_path = run_dir / "frozen" / "metadata.csv"
        contrasts_path = run_dir / "frozen" / "contrasts.csv"
        output = run_dir / "downstream" / "l2"
    columns, metadata_rows = _read_metadata(metadata_path)
    contrasts = _read_contrasts(contrasts_path, columns, metadata_rows, config.design.formula)
    _guard_replicates(config, metadata_rows, contrasts)
    return PreparedL2(l1, output, config, contrasts, metadata_rows)


def _parse_table(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise DownstreamExecutionError(f"Empty L2 result table: {path}")
        return reader.fieldnames, list(reader)


def _as_float(value: str) -> float | None:
    return None if value in {"", "NA"} else float(value)


def _validate_contrast_outputs(directory: Path, contrast: ContrastDefinition, thresholds: dict[str, float]) -> dict[str, Any]:
    required = ("all_genes.tsv", "significant.tsv", "upregulated.tsv", "downregulated.tsv", "volcano.png", "backend_summary.json")
    missing = [name for name in required if not (directory / name).is_file()]
    if missing:
        raise DownstreamExecutionError(f"Contrast {contrast.contrast_id} is missing required output(s): " + ", ".join(missing))
    header, all_rows = _parse_table(directory / "all_genes.tsv")
    if tuple(header) != RESULT_HEADER or not all_rows:
        raise DownstreamExecutionError(f"Contrast {contrast.contrast_id} master table has an invalid schema.")
    _, significant = _parse_table(directory / "significant.tsv")
    _, up = _parse_table(directory / "upregulated.tsv")
    _, down = _parse_table(directory / "downregulated.tsv")
    all_ids = {row["gene_id"] for row in all_rows}
    sig_ids, up_ids, down_ids = ({row["gene_id"] for row in rows} for rows in (significant, up, down))
    if not sig_ids <= all_ids or up_ids | down_ids != sig_ids or up_ids & down_ids:
        raise DownstreamExecutionError(f"Contrast {contrast.contrast_id} has inconsistent significant/up/down result tables.")
    for row in significant:
        padj, lfc = _as_float(row["padj"]), _as_float(row["log2FoldChange"])
        if padj is None or lfc is None or not (padj < thresholds["padj"] and abs(lfc) >= thresholds["abs_log2fc"]):
            raise DownstreamExecutionError(f"Contrast {contrast.contrast_id} significant table violates configured thresholds.")
    for row in up:
        if (_as_float(row["log2FoldChange"]) or 0) <= 0:
            raise DownstreamExecutionError(f"Contrast {contrast.contrast_id} upregulated table has nonpositive LFC.")
    for row in down:
        if (_as_float(row["log2FoldChange"]) or 0) >= 0:
            raise DownstreamExecutionError(f"Contrast {contrast.contrast_id} downregulated table has nonnegative LFC.")
    summary = json.loads((directory / "backend_summary.json").read_text(encoding="utf-8"))
    if not isinstance(summary, dict):
        raise DownstreamExecutionError(f"Contrast {contrast.contrast_id} backend summary is invalid.")
    summary.update({"contrast_id": contrast.contrast_id, "factor": contrast.factor, "numerator": contrast.numerator, "denominator": contrast.denominator, "table_genes": len(all_rows), "significant": len(significant), "upregulated": len(up), "downregulated": len(down)})
    return summary


def _contrast_report(summary: dict[str, Any], prepared: PreparedL2) -> str:
    lines = [
        f"# Contrast: {summary['contrast_id']}", "", "## Statistical definition", "",
        f"- Design: `{prepared.config.design.formula}`", f"- Factor: `{summary['factor']}`",
        f"- Positive log2FoldChange: higher in `{summary['numerator']}`", f"- Negative log2FoldChange: higher in `{summary['denominator']}`",
        f"- Samples: {summary['denominator']} n={summary['denominator_samples']}; {summary['numerator']} n={summary['numerator_samples']}",
        f"- Significance: `padj < {prepared.config.thresholds.padj}` and `abs(log2FoldChange) >= {prepared.config.thresholds.abs_log2fc}`", "",
        "## Results", "", f"- Input genes: {summary.get('input_genes')}", f"- Explicitly filtered genes: {summary.get('filtered_genes')}",
        f"- Tested genes: {summary.get('tested_genes')}", f"- p-value NA: {summary.get('pvalue_na')}", f"- padj NA: {summary.get('padj_na')}",
        f"- Significant: {summary['significant']}", f"- Upregulated: {summary['upregulated']}", f"- Downregulated: {summary['downregulated']}",
        f"- DESeq2 independent filtering: {summary.get('independent_filtering')}", f"- Multiple-testing method: {summary.get('multiple_testing_method')}",
        f"- Volcano: available", f"- Heatmap: {summary.get('heatmap_status')}", "",
        "This report contains statistical facts only; no biological interpretation or enrichment analysis was performed.", "",
    ]
    return "\n".join(lines)


def _overall_report(summaries: list[dict[str, Any]]) -> str:
    lines = ["# L2 differential-expression summary", "", "| contrast | tested | significant | up | down | status |", "|---|---:|---:|---:|---:|---|"]
    lines.extend(f"| {item['contrast_id']} | {item.get('tested_genes')} | {item['significant']} | {item['upregulated']} | {item['downregulated']} | {item.get('status', 'SUCCESS')} |" for item in summaries)
    lines.extend(["", "Positive log2FoldChange is numerator-enriched; negative log2FoldChange is denominator-enriched.", "No biological interpretation or functional enrichment was performed.", ""])
    return "\n".join(lines)


def _provenance(prepared: PreparedL2, l1_summary: dict[str, Any]) -> dict[str, Any]:
    thresholds = prepared.config.thresholds
    pairing: dict[str, Any] | None = None
    if prepared.config.design.type is DesignType.PAIRED_TWO_GROUP:
        pair_id = prepared.config.design.pair_id
        assert pair_id is not None
        contrast_stats = []
        for contrast in prepared.contrasts:
            grouped: dict[str, list[str]] = {}
            for row in prepared.metadata_rows:
                grouped.setdefault(row[pair_id], []).append(row[contrast.factor])
            complete = sum(
                sorted(levels) == sorted([contrast.numerator, contrast.denominator])
                for levels in grouped.values()
            )
            contrast_stats.append(
                {
                    "contrast_id": contrast.contrast_id,
                    "complete_pairs": complete,
                    "analyzed_samples": complete * 2,
                }
            )
        pairing = {
            "pair_id": pair_id,
            "unique_pair_ids": len({row[pair_id] for row in prepared.metadata_rows}),
            "contrasts": contrast_stats,
        }
    return {
        "pipeline": {"version": PIPELINE_VERSION, "analysis_level": "L2"},
        "input": prepared.l1.config,
        "design": {
            "type": prepared.config.design.type.value,
            "formula": prepared.config.design.formula,
            **({"pair_id": prepared.config.design.pair_id} if prepared.config.design.pair_id else {}),
            **({"pairing": pairing} if pairing is not None else {}),
        },
        "contrasts": [{"contrast_id": item.contrast_id, "factor": item.factor, "numerator": item.numerator, "denominator": item.denominator} for item in prepared.contrasts],
        "filtering": L1_FILTER,
        "deseq2": {"model_fit_count": 1, "independent_filtering": True, "multiple_testing_method": "Benjamini-Hochberg (DESeq2 default)"},
        "thresholds": {"padj": thresholds.padj, "abs_log2fc": thresholds.abs_log2fc, "significance": "padj < threshold AND abs(log2FoldChange) >= threshold"},
        "plots": {"vst_source": str(prepared.l1.output_dir / "vst.tsv"), "heatmap_gene_selection": f"top {HEATMAP_TOP_N} significant genes by padj, then gene_id"},
        "software": l1_summary.get("package_versions", {}),
        "automatic_sample_removal": False,
    }


def execute_l2(prepared: PreparedL2) -> L2Result:
    """Run deterministic L1 preparation then one DESeq2 fit for all simple contrasts."""
    _require_r()
    l1_result = execute_l1(prepared.l1)  # Explicit deterministic re-preparation; no L1 cache engine is introduced.
    output = prepared.output_dir
    (output / "contrasts").mkdir(parents=True, exist_ok=True)
    (output / "logs").mkdir(exist_ok=True)
    state_path = output / "l2_state.json"
    started = datetime.now(UTC).replace(microsecond=0).isoformat().replace("+00:00", "Z")
    _state(state_path, status="RUNNING", started_at=started, completed_at=None)
    config = {
        **prepared.l1.config, "metadata": str(prepared.l1.metadata_path.resolve()), "formula": prepared.config.design.formula,
        "pair_id": prepared.config.design.pair_id,
        "samples": list(prepared.l1.sample_ids), "output_dir": str(output.resolve()), "filter": L1_FILTER,
        "thresholds": {"padj": prepared.config.thresholds.padj, "abs_log2fc": prepared.config.thresholds.abs_log2fc},
        "contrasts": [{"contrast_id": item.contrast_id, "factor": item.factor, "numerator": item.numerator, "denominator": item.denominator} for item in prepared.contrasts],
        "l1_vst": str(prepared.l1.output_dir / "vst.tsv"), "heatmap_top_n": HEATMAP_TOP_N,
    }
    _write(output / "backend_config.json", json.dumps(config, indent=2, sort_keys=True) + "\n")
    provenance = _provenance(prepared, l1_result.summary)
    _write(output / "statistical_provenance.yaml", yaml.safe_dump(provenance, sort_keys=False))
    script = resources.files("rnaseq.r").joinpath("l2_analysis.R")
    command = ["Rscript", str(script), "--config", str(output / "backend_config.json")]
    try:
        with (output / "logs" / "r.stdout.log").open("w", encoding="utf-8", newline="\n") as stdout, (output / "logs" / "r.stderr.log").open("w", encoding="utf-8", newline="\n") as stderr:
            result = subprocess.run(command, stdout=stdout, stderr=stderr, check=False)
        if result.returncode != 0:
            raise DownstreamExecutionError(f"L2 R backend failed with return code {result.returncode}. Logs: {output / 'logs'}")
        summaries: list[dict[str, Any]] = []
        for contrast in prepared.contrasts:
            directory = output / "contrasts" / contrast.contrast_id
            summary = _validate_contrast_outputs(directory, contrast, config["thresholds"])
            summary["numerator_samples"] = sum(row[contrast.factor] == contrast.numerator for row in prepared.metadata_rows)
            summary["denominator_samples"] = sum(row[contrast.factor] == contrast.denominator for row in prepared.metadata_rows)
            summary["status"] = "SUCCESS"
            _write(directory / "summary.yaml", yaml.safe_dump(summary, sort_keys=False))
            _write(directory / "contrast_report.md", _contrast_report(summary, prepared))
            summaries.append(summary)
    except (OSError, ValueError, json.JSONDecodeError, DownstreamExecutionError) as exc:
        _state(state_path, status="FAILED", started_at=started, completed_at=datetime.now(UTC).replace(microsecond=0).isoformat().replace("+00:00", "Z"), error=str(exc))
        if isinstance(exc, DownstreamExecutionError):
            raise
        raise DownstreamExecutionError(f"Unable to execute L2 backend: {exc}. Logs: {output / 'logs'}") from exc
    provenance["deseq2"]["fit_method"] = summaries[0].get("fit_method")
    _write(output / "statistical_provenance.yaml", yaml.safe_dump(provenance, sort_keys=False))
    _write(output / "l2_summary.md", _overall_report(summaries))
    _state(state_path, status="SUCCESS", started_at=started, completed_at=datetime.now(UTC).replace(microsecond=0).isoformat().replace("+00:00", "Z"))
    return L2Result(output, state_path, "SUCCESS", tuple(summaries))
