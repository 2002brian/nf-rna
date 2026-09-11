"""Offline GO preranked GSEA after an existing successful M4A L2 run."""

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
from rnaseq.go import ORGDB, _check_go_runtime
from rnaseq.l2 import PreparedL2, prepare_l2
from rnaseq.models import AnnotationConfig
from rnaseq.validators import ValidationReport

GSEA_STATES = {"SUCCESS", "BLOCKED", "NO_SIGNIFICANT_TERMS", "FAILED"}


@dataclass(frozen=True)
class PreparedGsea:
    l2: PreparedL2
    annotation: AnnotationConfig
    output_dir: Path


@dataclass(frozen=True)
class GseaResult:
    output_dir: Path
    state_path: Path
    status: str
    summary: dict[str, Any]


def _state(path: Path, **value: Any) -> None:
    _write(path, json.dumps(value, indent=2, sort_keys=True) + "\n")


def prepare_gsea(report: ValidationReport, *, run_id: str | None) -> PreparedGsea:
    """Require L2 outputs but never re-run DESeq2 or invoke ORA."""
    l2 = prepare_l2(report, run_id=run_id)
    annotation = l2.config.annotation
    if annotation is None:
        organism = l2.config.organism.species.value
        if organism not in ORGDB:
            raise DownstreamExecutionError(f"GO preranked GSEA is not currently supported for organism {organism}.")
        raise DownstreamExecutionError(
            "GO preranked GSEA requires an explicit annotation contract in project.yaml "
            "(annotation.organism and annotation.input_id_type)."
        )
    if annotation.organism not in ORGDB:
        raise DownstreamExecutionError(f"GO preranked GSEA is not currently supported for organism {annotation.organism}.")
    state_path = l2.output_dir / "l2_state.json"
    if not state_path.is_file() or json.loads(state_path.read_text(encoding="utf-8")).get("status") != "SUCCESS":
        raise DownstreamExecutionError(
            "GO preranked GSEA requires an existing successful L2 DE analysis; "
            "run 'rnaseq analyze PROJECT --level L2' first."
        )
    for contrast in l2.contrasts:
        table = l2.output_dir / "contrasts" / contrast.contrast_id / "all_genes.tsv"
        if not table.is_file():
            raise DownstreamExecutionError(f"Successful L2 output is missing {table}.")
    return PreparedGsea(l2, annotation, l2.output_dir / "enrichment" / "gsea_go")


def _gsea_report(summary: dict[str, Any], annotation: AnnotationConfig) -> str:
    lines = ["# GO preranked GSEA", "", f"Status: {summary['status']}", ""]
    lines.extend([
        "- Ranking source: successful M4A `all_genes.tsv` tested genes with finite DESeq2 `stat`.",
        "- Ranking metric: DESeq2 Wald statistic; positive values favor the contrast numerator and negative values favor the denominator.",
        "- No p-value, adjusted-p-value, significance, or log2-fold-change prefilter was applied.",
        "- One-to-many source IDs propagate the same statistic; duplicate Entrez targets retain the largest absolute statistic (lexical source-ID tie-break).",
        "- Tied statistics are not perturbed; the ranked vector uses Entrez ID as a deterministic secondary sort key before fgsea.",
        f"- Minimum ranked genes: {annotation.enrichment.gsea.minimum_ranked_genes}",
        f"- Gene-set size: {annotation.enrichment.gsea.min_gs_size}–{annotation.enrichment.gsea.max_gs_size}",
        f"- Thresholds: pvalue <= {annotation.enrichment.gsea.pvalue_cutoff}; adjusted p <= {annotation.enrichment.gsea.padj_cutoff}; adjustment {annotation.enrichment.gsea.p_adjust_method}",
        "- `all_terms.tsv` contains every structurally eligible term successfully evaluated and returned by gseGO with calculation pvalueCutoff=1; `significant.tsv` applies nf-rna's configured p-value and adjusted-p filters to those rows.",
        "- BP, MF, and CC are run separately with completed ontology objects released before the next ontology. No biological interpretation is generated.", "",
    ])
    if summary["status"] == "BLOCKED":
        lines.extend([summary.get("reason", "GSEA ranking guardrail blocked execution."), ""])
    for item in summary.get("contrasts", []):
        rank = item.get("ranking", {})
        annotation_qc = rank.get("annotation_qc", {})
        lines.extend([
            f"## {item['contrast_id']}",
            f"- Ranked Entrez genes: {rank.get('final_ranked_genes', 0)} (positive={rank.get('positive_stats', 0)}, negative={rank.get('negative_stats', 0)}, zero={rank.get('zero_stats', 0)})",
            f"- Mapping rate: {rank.get('mapping_rate', 0):.1%}",
            f"- Annotation mapping QC: {annotation_qc.get('status', 'not available')} (warning threshold={annotation_qc.get('warning_threshold', 'not available')}; blocking threshold={annotation_qc.get('blocking_threshold', 'not available')})",
        ])
        if annotation_qc.get("status") == "WARNING":
            lines.append(f"- Annotation mapping warning: {annotation_qc.get('reason', 'not available')}")
        for ontology, values in item.get("ontologies", {}).items():
            lines.append(f"- {ontology}: {values['status']}; all terms={values['all_terms']}; significant={values['significant_terms']}; NA pathways={values.get('na_pathways', 0)}")
        if item.get("failed_ontology"):
            lines.append(f"- Failed ontology: {item['failed_ontology']}; completed: {', '.join(item.get('completed_ontologies', []))}")
        lines.append("")
    return "\n".join(lines)


def execute_gsea(prepared: PreparedGsea) -> GseaResult:
    """Run local clusterProfiler gseGO with a deterministic, audited ranked vector."""
    _check_go_runtime(prepared.annotation.organism)
    output = prepared.output_dir
    (output / "logs").mkdir(parents=True, exist_ok=True)
    state_path = output / "gsea_state.json"
    started = datetime.now(UTC).replace(microsecond=0).isoformat().replace("+00:00", "Z")
    _state(state_path, status="RUNNING", started_at=started, completed_at=None)
    contrasts = [
        {
            "contrast_id": contrast.contrast_id,
            "all_genes": str(prepared.l2.output_dir / "contrasts" / contrast.contrast_id / "all_genes.tsv"),
        }
        for contrast in prepared.l2.contrasts
    ]
    cfg = {
        "output_dir": str(output.resolve()),
        "annotation": prepared.annotation.model_dump(),
        "orgdb_package": ORGDB[prepared.annotation.organism],
        "contrasts": contrasts,
    }
    _write(output / "backend_config.json", json.dumps(cfg, indent=2, sort_keys=True) + "\n")
    script = resources.files("rnaseq.r").joinpath("gsea_analysis.R")
    try:
        with (output / "logs" / "r.stdout.log").open("w", encoding="utf-8", newline="\n") as stdout, (output / "logs" / "r.stderr.log").open("w", encoding="utf-8", newline="\n") as stderr:
            result = subprocess.run(["Rscript", str(script), "--config", str(output / "backend_config.json")], stdout=stdout, stderr=stderr, check=False)
        if result.returncode != 0:
            raise DownstreamExecutionError(f"GO preranked GSEA R backend failed with return code {result.returncode}. Logs: {output / 'logs'}")
        summary = json.loads((output / "gsea_backend_summary.json").read_text(encoding="utf-8"))
        if not isinstance(summary, dict) or summary.get("status") not in GSEA_STATES:
            raise DownstreamExecutionError("GO preranked GSEA backend summary is invalid.")
    except (OSError, json.JSONDecodeError, DownstreamExecutionError) as exc:
        _state(state_path, status="FAILED", started_at=started, completed_at=datetime.now(UTC).replace(microsecond=0).isoformat().replace("+00:00", "Z"), error=str(exc))
        if isinstance(exc, DownstreamExecutionError):
            raise
        raise DownstreamExecutionError(f"Unable to run GO preranked GSEA: {exc}") from exc
    _write(output / "gsea_summary.yaml", yaml.safe_dump(summary, sort_keys=False))
    for contrast in summary.get("contrasts", []):
        ranking_path = output / contrast["contrast_id"] / "gsea_ranking_summary.yaml"
        ranking_path.parent.mkdir(parents=True, exist_ok=True)
        ranking_summary = {
            "contrast_id": contrast["contrast_id"],
            "ranking_source": "M4A all_genes.tsv; finite DESeq2 stat rows only",
            "ranking_metric": "DESeq2 Wald stat",
            "direction": "positive = numerator-enriched; negative = denominator-enriched",
            "target_collapse": "largest absolute stat, then lexical original gene ID",
            "target_sort": "decreasing stat, then ascending Entrez ID",
            **contrast.get("ranking", {}),
            "status": contrast.get("status", summary["status"]),
            "reason": contrast.get("reason"),
        }
        _write(ranking_path, yaml.safe_dump(ranking_summary, sort_keys=False))
    _write(output / "gsea_report.md", _gsea_report(summary, prepared.annotation))
    provenance = {
        "project_id": prepared.l2.config.project.id,
        "ranking_source": "M4A all_genes.tsv; finite DESeq2 stat rows only",
        "ranking_metric": "DESeq2 Wald stat",
        "direction": "positive = numerator-enriched; negative = denominator-enriched",
        "annotation": prepared.annotation.model_dump(),
        "annotation_database": summary.get("annotation_database"),
        "annotation_database_version": summary.get("annotation_database_version"),
        "annotation_qc_status": summary.get("annotation_qc_status"),
        "clusterProfiler_version": summary.get("clusterProfiler_version"),
        "mapping_policy": "one-to-many source IDs retain all valid Entrez targets; duplicate targets retain largest absolute stat, then lexical original gene ID",
        "tie_handling": "DESeq2 stat is unmodified; ties are ordered by ascending Entrez ID before fgsea.",
        "memory_lifecycle": "BP, MF, and CC run sequentially; result objects and intermediates are released and garbage-collected after each ontology",
        "gsea_parameters": prepared.annotation.enrichment.gsea.model_dump(),
        "clusterprofiler_calculation_pvalue_cutoff": 1,
        "nf_rna_significance_policy": "finite p.adjust <= configured padj_cutoff and finite pvalue <= configured pvalue_cutoff",
        "all_terms_semantics": "terms successfully evaluated and returned by gseGO under configured structural gene-set constraints before nf-rna significance filtering",
        "ontologies": ["BP", "MF", "CC"],
        "automatic_sample_removal": False,
    }
    _write(output / "gsea_provenance.yaml", yaml.safe_dump(provenance, sort_keys=False))
    ontology_progress = [
        {
            "contrast_id": contrast["contrast_id"],
            "completed_ontologies": contrast.get("completed_ontologies", []),
            "failed_ontology": contrast.get("failed_ontology"),
        }
        for contrast in summary.get("contrasts", [])
    ]
    _state(
        state_path,
        status=summary["status"],
        started_at=started,
        completed_at=datetime.now(UTC).replace(microsecond=0).isoformat().replace("+00:00", "Z"),
        ontology_progress=ontology_progress,
    )
    return GseaResult(output, state_path, summary["status"], summary)
