"""Offline GO over-representation analysis after a successful M4A L2 run."""

from __future__ import annotations

import json
import subprocess
from dataclasses import dataclass
from datetime import UTC, datetime
from importlib import resources
from pathlib import Path
from typing import Any

import yaml

from rnaseq.downstream import _write
from rnaseq.errors import DownstreamExecutionError
from rnaseq.l2 import PreparedL2, prepare_l2
from rnaseq.models import AnnotationConfig
from rnaseq.validators import ValidationReport

GO_STATES = {"SUCCESS", "BLOCKED", "NOT_APPLICABLE", "NO_SIGNIFICANT_TERMS", "FAILED"}
ORGDB = {"Homo sapiens": "org.Hs.eg.db", "Mus musculus": "org.Mm.eg.db"}


@dataclass(frozen=True)
class PreparedGo:
    l2: PreparedL2
    annotation: AnnotationConfig
    output_dir: Path


@dataclass(frozen=True)
class GoResult:
    output_dir: Path
    state_path: Path
    status: str
    summary: dict[str, Any]


def _state(path: Path, **value: Any) -> None:
    _write(path, json.dumps(value, indent=2, sort_keys=True) + "\n")


def _check_go_runtime(organism: str) -> None:
    package = ORGDB[organism]
    script = (
        "required <- c('clusterProfiler','AnnotationDbi','" + package + "'); "
        "quit(status=if (all(vapply(required, requireNamespace, logical(1), quietly=TRUE))) 0 else 1)"
    )
    try:
        result = subprocess.run(["Rscript", "-e", script], capture_output=True, text=True, check=False)
    except FileNotFoundError as exc:
        raise DownstreamExecutionError("GO enrichment requires Rscript, which was not found.") from exc
    if result.returncode != 0:
        raise DownstreamExecutionError(f"GO enrichment requires local packages clusterProfiler, AnnotationDbi, and {package}.")


def prepare_go(report: ValidationReport, *, run_id: str | None) -> PreparedGo:
    l2 = prepare_l2(report, run_id=run_id)
    annotation = l2.config.annotation
    if annotation is None:
        organism = l2.config.organism.species.value
        if organism not in ORGDB:
            raise DownstreamExecutionError(f"GO enrichment is not currently supported for organism {organism}.")
        raise DownstreamExecutionError(
            "GO enrichment requires an explicit annotation contract in project.yaml (annotation.organism and annotation.input_id_type)."
        )
    if annotation.organism not in ORGDB:
        raise DownstreamExecutionError(f"GO enrichment is not currently supported for organism {annotation.organism}.")
    state_path = l2.output_dir / "l2_state.json"
    if not state_path.is_file() or json.loads(state_path.read_text(encoding="utf-8")).get("status") != "SUCCESS":
        raise DownstreamExecutionError("GO enrichment requires an existing successful L2 DE analysis; run 'rnaseq analyze PROJECT --level L2' first.")
    for contrast in l2.contrasts:
        required = l2.output_dir / "contrasts" / contrast.contrast_id / "all_genes.tsv"
        if not required.is_file():
            raise DownstreamExecutionError(f"Successful L2 output is missing {required}.")
    return PreparedGo(l2, annotation, l2.output_dir / "enrichment" / "go")


def _mapping_report(summary: dict[str, Any], annotation: AnnotationConfig) -> str:
    tested = summary["mapping"]["tested"]
    lines = [
        "# Gene identifier mapping audit", "", f"- Organism: `{annotation.organism}`",
        f"- Input identifier type: `{annotation.input_id_type}`", "- Target identifier type: `ENTREZID`",
        f"- Tested genes: {tested['source_genes']}", f"- Mapped source genes: {tested['mapped_source_genes']}",
        f"- Unique mapped target genes: {tested['unique_target_genes']}", f"- Unmapped source genes: {tested['unmapped_source_genes']}",
        f"- One-to-many source IDs: {tested['one_to_many_source_ids']}", f"- Duplicate target IDs: {tested['duplicate_target_ids']}",
        f"- Mapping rate: {tested['mapping_rate']:.1%}", f"- Required minimum mapping rate: {annotation.mapping_warning_rate:.1%}", "",
        "", "## Contrast foregrounds", "",
    ]
    for contrast in summary.get("contrasts", []):
        lines.append(f"### {contrast['contrast_id']}")
        for name, values in contrast["foregrounds"].items():
            lines.append(f"- {name}: source={values['source_genes']}, mapped={values['mapped_source_genes']}, unique targets={values['unique_target_genes']}, unmapped={values['unmapped_source_genes']}")
    lines.extend(["", "One-to-many mappings retain every valid Entrez target for ORA; repeated target IDs are deduplicated before ORA. Source DE tables are not modified.", ""])
    return "\n".join(lines)


def _go_report(summary: dict[str, Any], annotation: AnnotationConfig) -> str:
    lines = ["# GO over-representation analysis", "", f"Status: {summary['status']}", ""]
    if summary["status"] == "BLOCKED":
        lines.extend([summary["reason"], ""])
        return "\n".join(lines)
    lines.extend([
        f"- Universe: {summary['mapping']['tested']['unique_target_genes']} mapped, tested Entrez genes",
        f"- Minimum mapped foreground size: {annotation.minimum_mapped_foreground}",
        f"- GO thresholds: pvalue <= {annotation.enrichment.go.pvalue_cutoff}; qvalue <= {annotation.enrichment.go.qvalue_cutoff}; adjustment {annotation.enrichment.go.p_adjust_method}",
        "- BP, MF, and CC are evaluated separately for significant, upregulated, and downregulated foregrounds.",
        "- This report contains no biological interpretation.", "",
    ])
    for item in summary.get("contrasts", []):
        lines.append(f"## {item['contrast_id']}")
        for foreground, values in item["foregrounds"].items():
            lines.append(f"- {foreground}: {values['status']}; mapped foreground={values['unique_target_genes']}; BP/MF/CC terms={values['term_counts']}")
        lines.append("")
    return "\n".join(lines)


def execute_go(prepared: PreparedGo) -> GoResult:
    """Map once and run local clusterProfiler GO ORA without altering L2 outputs."""
    _check_go_runtime(prepared.annotation.organism)
    output = prepared.output_dir
    (output / "logs").mkdir(parents=True, exist_ok=True)
    state_path = output / "go_state.json"
    started = datetime.now(UTC).replace(microsecond=0).isoformat().replace("+00:00", "Z")
    _state(state_path, status="RUNNING", started_at=started, completed_at=None)
    contrasts = []
    for contrast in prepared.l2.contrasts:
        root = prepared.l2.output_dir / "contrasts" / contrast.contrast_id
        contrasts.append({
            "contrast_id": contrast.contrast_id,
            "all_genes": str(root / "all_genes.tsv"),
            "significant": str(root / "significant.tsv"),
            "up": str(root / "upregulated.tsv"),
            "down": str(root / "downregulated.tsv"),
        })
    cfg = {
        "output_dir": str(output.resolve()), "annotation": prepared.annotation.model_dump(),
        "orgdb_package": ORGDB[prepared.annotation.organism], "contrasts": contrasts,
    }
    _write(output / "backend_config.json", json.dumps(cfg, indent=2, sort_keys=True) + "\n")
    script = resources.files("rnaseq.r").joinpath("go_analysis.R")
    try:
        with (output / "logs" / "r.stdout.log").open("w", encoding="utf-8", newline="\n") as stdout, (output / "logs" / "r.stderr.log").open("w", encoding="utf-8", newline="\n") as stderr:
            result = subprocess.run(["Rscript", str(script), "--config", str(output / "backend_config.json")], stdout=stdout, stderr=stderr, check=False)
        if result.returncode != 0:
            raise DownstreamExecutionError(f"GO R backend failed with return code {result.returncode}. Logs: {output / 'logs'}")
        summary = json.loads((output / "go_backend_summary.json").read_text(encoding="utf-8"))
        if not isinstance(summary, dict) or summary.get("status") not in GO_STATES:
            raise DownstreamExecutionError("GO backend summary is invalid.")
        if not (prepared.l2.output_dir / "annotation" / "gene_mapping.tsv").is_file():
            raise DownstreamExecutionError("GO backend did not produce annotation/gene_mapping.tsv.")
    except (OSError, json.JSONDecodeError, DownstreamExecutionError) as exc:
        _state(state_path, status="FAILED", started_at=started, completed_at=datetime.now(UTC).replace(microsecond=0).isoformat().replace("+00:00", "Z"), error=str(exc))
        if isinstance(exc, DownstreamExecutionError):
            raise
        raise DownstreamExecutionError(f"Unable to run GO enrichment: {exc}") from exc
    annotation_dir = prepared.l2.output_dir / "annotation"
    annotation_dir.mkdir(exist_ok=True)
    _write(annotation_dir / "mapping_summary.yaml", yaml.safe_dump({"tested": summary["mapping"]["tested"], "contrasts": summary.get("contrasts", [])}, sort_keys=False))
    _write(annotation_dir / "mapping_report.md", _mapping_report(summary, prepared.annotation))
    _write(output / "go_summary.yaml", yaml.safe_dump(summary, sort_keys=False))
    _write(output / "go_report.md", _go_report(summary, prepared.annotation))
    provenance = {
        "project_id": prepared.l2.config.project.id,
        "input": prepared.l2.l1.config,
        "annotation": prepared.annotation.model_dump(),
        "annotation_database": summary.get("annotation_database"),
        "annotation_database_version": summary.get("annotation_database_version"),
        "clusterProfiler_version": summary.get("clusterProfiler_version"),
        "mapping_policy": "retain all valid Entrez targets for one-to-many source IDs; deduplicate targets before ORA",
        "mapping_rate_guardrail": prepared.annotation.mapping_warning_rate,
        "minimum_mapped_foreground": prepared.annotation.minimum_mapped_foreground,
        "universe": summary["mapping"]["tested"],
        "thresholds": prepared.annotation.enrichment.go.model_dump(),
        "ontologies": ["BP", "MF", "CC"],
        "automatic_sample_removal": False,
    }
    _write(output / "go_provenance.yaml", yaml.safe_dump(provenance, sort_keys=False))
    _state(state_path, status=summary["status"], started_at=started, completed_at=datetime.now(UTC).replace(microsecond=0).isoformat().replace("+00:00", "Z"))
    return GoResult(output, state_path, summary["status"], summary)
