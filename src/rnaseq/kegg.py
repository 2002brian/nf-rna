"""Explicit online KEGG ORA and preranked GSEA after a successful M4A L2 run."""

from __future__ import annotations

import json
import subprocess
from dataclasses import dataclass
from datetime import UTC, datetime
from importlib import resources
from pathlib import Path
from typing import Any, Literal

import yaml

from rnaseq.downstream import _write
from rnaseq.errors import DownstreamExecutionError
from rnaseq.go import ORGDB
from rnaseq.l2 import PreparedL2, prepare_l2
from rnaseq.models import AnnotationConfig
from rnaseq.validators import ValidationReport

KEGG_CODES = {"Homo sapiens": "hsa", "Mus musculus": "mmu"}
KEGG_STATES = {"SUCCESS", "BLOCKED", "NETWORK_UNAVAILABLE", "NOT_APPLICABLE", "NO_SIGNIFICANT_TERMS", "FAILED"}
KeggMode = Literal["ora", "gsea"]


@dataclass(frozen=True)
class KeggResourceAdapter:
    """The sole M4B-3 provider: live KEGG REST, consumed by clusterProfiler."""

    organism_code: str
    provider: str = "online_kegg_rest_via_clusterprofiler"

    @property
    def probe_endpoint(self) -> str:
        return f"https://rest.kegg.jp/list/pathway/{self.organism_code}"


@dataclass(frozen=True)
class PreparedKegg:
    l2: PreparedL2
    annotation: AnnotationConfig
    adapter: KeggResourceAdapter
    mode: KeggMode
    output_dir: Path


@dataclass(frozen=True)
class KeggResult:
    output_dir: Path
    state_path: Path
    status: str
    summary: dict[str, Any]


def _state(path: Path, **value: Any) -> None:
    _write(path, json.dumps(value, indent=2, sort_keys=True) + "\n")


def _check_kegg_runtime() -> None:
    script = "quit(status=if (requireNamespace('clusterProfiler', quietly=TRUE) && requireNamespace('AnnotationDbi', quietly=TRUE)) 0 else 1)"
    try:
        result = subprocess.run(["Rscript", "-e", script], capture_output=True, text=True, check=False)
    except FileNotFoundError as exc:
        raise DownstreamExecutionError("KEGG enrichment requires Rscript, which was not found.") from exc
    if result.returncode != 0:
        raise DownstreamExecutionError("KEGG enrichment requires local clusterProfiler and AnnotationDbi packages.")


def prepare_kegg(report: ValidationReport, *, run_id: str | None, mode: KeggMode) -> PreparedKegg:
    """Prepare one isolated KEGG action without running DE, GO, or a network probe."""
    l2 = prepare_l2(report, run_id=run_id)
    annotation = l2.config.annotation
    if annotation is None:
        organism = l2.config.organism.species.value
        if organism not in KEGG_CODES:
            raise DownstreamExecutionError(f"KEGG enrichment is not currently supported for organism {organism}.")
        raise DownstreamExecutionError(
            "KEGG enrichment requires an explicit annotation contract in project.yaml "
            "(annotation.organism and annotation.input_id_type)."
        )
    if annotation.organism not in KEGG_CODES:
        raise DownstreamExecutionError(f"KEGG enrichment is not currently supported for organism {annotation.organism}.")
    state_path = l2.output_dir / "l2_state.json"
    if not state_path.is_file() or json.loads(state_path.read_text(encoding="utf-8")).get("status") != "SUCCESS":
        raise DownstreamExecutionError(
            "KEGG enrichment requires an existing successful L2 DE analysis; "
            "run 'rnaseq analyze PROJECT --level L2' first."
        )
    for contrast in l2.contrasts:
        table = l2.output_dir / "contrasts" / contrast.contrast_id / "all_genes.tsv"
        if not table.is_file():
            raise DownstreamExecutionError(f"Successful L2 output is missing {table}.")
    adapter = KeggResourceAdapter(KEGG_CODES[annotation.organism])
    name = "kegg" if mode == "ora" else "gsea_kegg"
    return PreparedKegg(l2, annotation, adapter, mode, l2.output_dir / "enrichment" / name)


def _report(mode: KeggMode, summary: dict[str, Any], adapter: KeggResourceAdapter) -> str:
    title = "# KEGG over-representation analysis" if mode == "ora" else "# KEGG preranked GSEA"
    lines = [title, "", f"Status: {summary['status']}", "", f"- KEGG organism code: `{adapter.organism_code}`", f"- Resource provider: `{adapter.provider}`", "- KEGG resources are externally retrieved and may change over time.", "- No biological interpretation is generated.", ""]
    if summary["status"] in {"BLOCKED", "NETWORK_UNAVAILABLE"}:
        lines.extend([summary.get("reason", "KEGG resource or mapping guardrail prevented execution."), ""])
    for contrast in summary.get("contrasts", []):
        lines.append(f"## {contrast['contrast_id']}")
        if mode == "ora":
            lines.append(f"- Tested mapped KEGG universe: {contrast.get('universe_size', 0)}")
            for name, values in contrast.get("foregrounds", {}).items():
                lines.append(f"- {name}: {values['status']}; foreground={values['unique_target_genes']}; pathways={values['pathway_count']}")
        else:
            rank = contrast.get("ranking", {})
            lines.append(f"- Ranked genes: {rank.get('final_ranked_genes', 0)}; positive={rank.get('positive_stats', 0)}; negative={rank.get('negative_stats', 0)}")
            annotation_qc = rank.get("annotation_qc", {})
            lines.append(f"- Annotation mapping QC: {annotation_qc.get('status', 'not available')} (warning threshold={annotation_qc.get('warning_threshold', 'not available')}; blocking threshold={annotation_qc.get('blocking_threshold', 'not available')})")
            if annotation_qc.get("status") == "WARNING":
                lines.append(f"- Annotation mapping warning: {annotation_qc.get('reason', 'not available')}")
            lines.append(f"- Pathways: evaluated={contrast.get('evaluated_terms', contrast.get('all_terms', 0))}; significant={contrast.get('significant_terms', 0)}")
        lines.append("")
    return "\n".join(lines)


def execute_kegg(prepared: PreparedKegg) -> KeggResult:
    """Run one KEGG action through its explicit live-resource adapter."""
    _check_kegg_runtime()
    output = prepared.output_dir
    (output / "logs").mkdir(parents=True, exist_ok=True)
    state_name = "kegg_state.json" if prepared.mode == "ora" else "gsea_kegg_state.json"
    state_path = output / state_name
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
            "go_gsea_ranked": str(prepared.l2.output_dir / "enrichment" / "gsea_go" / contrast.contrast_id / "ranked_gene_list.tsv"),
        })
    cfg = {
        "mode": prepared.mode,
        "output_dir": str(output.resolve()),
        "annotation": prepared.annotation.model_dump(),
        "orgdb_package": ORGDB[prepared.annotation.organism],
        "kegg": {"organism_code": prepared.adapter.organism_code, "provider": prepared.adapter.provider, "probe_endpoint": prepared.adapter.probe_endpoint},
        "contrasts": contrasts,
    }
    _write(output / "backend_config.json", json.dumps(cfg, indent=2, sort_keys=True) + "\n")
    script = resources.files("rnaseq.r").joinpath("kegg_analysis.R")
    summary_name = "kegg_backend_summary.json" if prepared.mode == "ora" else "gsea_kegg_backend_summary.json"
    try:
        with (output / "logs" / "r.stdout.log").open("w", encoding="utf-8", newline="\n") as stdout, (output / "logs" / "r.stderr.log").open("w", encoding="utf-8", newline="\n") as stderr:
            result = subprocess.run(["Rscript", str(script), "--config", str(output / "backend_config.json")], stdout=stdout, stderr=stderr, check=False)
        if result.returncode != 0:
            raise DownstreamExecutionError(f"KEGG {prepared.mode} R backend failed with return code {result.returncode}. Logs: {output / 'logs'}")
        summary = json.loads((output / summary_name).read_text(encoding="utf-8"))
        if not isinstance(summary, dict) or summary.get("status") not in KEGG_STATES:
            raise DownstreamExecutionError("KEGG backend summary is invalid.")
    except (OSError, json.JSONDecodeError, DownstreamExecutionError) as exc:
        _state(state_path, status="FAILED", started_at=started, completed_at=datetime.now(UTC).replace(microsecond=0).isoformat().replace("+00:00", "Z"), error=str(exc))
        if isinstance(exc, DownstreamExecutionError):
            raise
        raise DownstreamExecutionError(f"Unable to run KEGG {prepared.mode}: {exc}") from exc
    prefix = "kegg" if prepared.mode == "ora" else "gsea_kegg"
    _write(output / f"{prefix}_summary.yaml", yaml.safe_dump(summary, sort_keys=False))
    _write(output / f"{prefix}_report.md", _report(prepared.mode, summary, prepared.adapter))
    provenance = {
        "project_id": prepared.l2.config.project.id,
        "mode": prepared.mode,
        "organism": prepared.annotation.organism,
        "kegg_organism_code": prepared.adapter.organism_code,
        "input_id_namespace": prepared.annotation.input_id_type,
        "target_id_namespace": "ENTREZID / ncbi-geneid",
        "resource_provider": prepared.adapter.provider,
        "resource_probe_endpoint": prepared.adapter.probe_endpoint,
        "retrieval": summary.get("resource"),
        "clusterProfiler_version": summary.get("clusterProfiler_version"),
        "r_version": summary.get("r_version"),
        "annotation_database": summary.get("annotation_database"),
        "annotation_database_version": summary.get("annotation_database_version"),
        "annotation_qc_status": summary.get("annotation_qc_status"),
        "annotation": prepared.annotation.model_dump(),
        "clusterprofiler_calculation_pvalue_cutoff": 1 if prepared.mode == "gsea" else None,
        "nf_rna_significance_policy": (
            "finite p.adjust <= configured padj_cutoff and finite pvalue <= configured pvalue_cutoff"
            if prepared.mode == "gsea" else None
        ),
        "all_terms_semantics": (
            "terms successfully evaluated and returned by gseKEGG under configured structural gene-set constraints before nf-rna significance filtering"
            if prepared.mode == "gsea" else None
        ),
        "mapping_policy": "local AnnotationDbi Entrez mapping; one-to-many retained; duplicate KEGG targets deduplicated deterministically",
        "warning": "KEGG results may vary when the external KEGG resource changes.",
        "contrasts": summary.get("contrasts", []),
        "automatic_sample_removal": False,
    }
    _write(output / f"{prefix}_provenance.yaml", yaml.safe_dump(provenance, sort_keys=False))
    _state(state_path, status=summary["status"], started_at=started, completed_at=datetime.now(UTC).replace(microsecond=0).isoformat().replace("+00:00", "Z"))
    return KeggResult(output, state_path, summary["status"], summary)
