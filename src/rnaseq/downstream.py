"""Bounded, provenance-aware Milestone 3 L1 expression-QC orchestration.

Python owns validation, input selection, state, and output checks.  The
mathematical work is delegated to a small R/Bioconductor backend; this module
deliberately contains no differential-expression testing implementation.
"""

from __future__ import annotations

import csv
import json
import shutil
import subprocess
from dataclasses import dataclass
from datetime import UTC, datetime
from importlib import resources
from pathlib import Path
from typing import Any

import yaml

from rnaseq.errors import DownstreamExecutionError
from rnaseq.execution import RuntimeCheck
from rnaseq.models import InputType
from rnaseq.validators import ValidationReport

L1_FILTER = {"rule": "remove genes with total count < 10 after all-zero removal", "minimum_total_count": 10}
L1_OUTPUTS = (
    "normalized_counts.tsv", "vst.tsv", "library_size_qc.tsv", "pca_scores.tsv",
    "pca_variance.tsv", "sample_correlation.tsv", "pca.png", "sample_correlation.png",
    "backend_summary.json",
)


@dataclass(frozen=True)
class PreparedL1:
    project_dir: Path
    output_dir: Path
    source_type: str
    sample_ids: tuple[str, ...]
    metadata_path: Path
    formula: str
    config: dict[str, Any]


@dataclass(frozen=True)
class L1Result:
    output_dir: Path
    state_path: Path
    summary: dict[str, Any]


def _capture(args: list[str]) -> subprocess.CompletedProcess[str]:
    return subprocess.run(args, capture_output=True, text=True, check=False)


def r_runtime_checks() -> tuple[RuntimeCheck, ...]:
    """Report R and the explicit M3 Bioconductor dependencies without installing them."""
    try:
        result = _capture(["Rscript", "--version"])
    except FileNotFoundError:
        return (
            RuntimeCheck("Rscript", "NOT FOUND", "Rscript executable was not found on PATH."),
            RuntimeCheck("DESeq2", "NOT FOUND", "Rscript is unavailable."),
            RuntimeCheck("tximport", "NOT FOUND", "Rscript is unavailable."),
        )
    if result.returncode != 0:
        detail = (result.stderr or result.stdout).strip() or "Rscript --version failed."
        return (RuntimeCheck("Rscript", "NOT FOUND", detail),)
    checks = [RuntimeCheck("Rscript", "FOUND", (result.stdout or result.stderr).strip())]
    for package in ("DESeq2", "tximport", "ggplot2", "pheatmap", "yaml", "jsonlite", "clusterProfiler", "AnnotationDbi", "org.Hs.eg.db", "org.Mm.eg.db"):
        probe = _capture(["Rscript", "-e", f"quit(status=if (requireNamespace('{package}', quietly=TRUE)) 0 else 1)"])
        checks.append(RuntimeCheck(package, "FOUND" if probe.returncode == 0 else "NOT FOUND", "available" if probe.returncode == 0 else "required for L1"))
    return tuple(checks)


def _require_r() -> None:
    missing = [item.name for item in r_runtime_checks() if item.state != "FOUND"]
    if missing:
        raise DownstreamExecutionError("L1 requires unavailable runtime component(s): " + ", ".join(missing))


def _read_yaml(path: Path, label: str) -> dict[str, Any]:
    try:
        value = yaml.safe_load(path.read_text(encoding="utf-8"))
    except (OSError, yaml.YAMLError) as exc:
        raise DownstreamExecutionError(f"Cannot read {label}: {exc}") from exc
    if not isinstance(value, dict):
        raise DownstreamExecutionError(f"{label} must contain a YAML mapping.")
    return value


def _relative_existing(run_dir: Path, value: Any, label: str) -> Path:
    if not isinstance(value, str) or not value:
        raise DownstreamExecutionError(f"Handoff is missing {label}.")
    candidate = (run_dir / value).resolve()
    try:
        candidate.relative_to(run_dir.resolve())
    except ValueError as exc:
        raise DownstreamExecutionError(f"Handoff {label} must remain inside its run directory.") from exc
    if not candidate.is_file():
        raise DownstreamExecutionError(f"Handoff {label} does not exist: {candidate}")
    return candidate


def _upgrade_salmon_handoff(run_dir: Path, manifest: dict[str, Any]) -> dict[str, Any]:
    """Add only real M2 Salmon artifact references to a legacy successful handoff."""
    if isinstance(manifest.get("salmon"), dict):
        return manifest
    samples = manifest.get("samples")
    if not isinstance(samples, list) or not all(isinstance(item, str) and item for item in samples):
        raise DownstreamExecutionError("Legacy handoff has no valid sample list for Salmon import.")
    output = run_dir / "upstream" / "nfcore_rnaseq" / "salmon"
    quant = {sample: output / sample / "quant.sf" for sample in samples}
    tx2gene = output / "salmon.merged.tx2gene.tsv"
    if not tx2gene.is_file() or any(not path.is_file() for path in quant.values()):
        raise DownstreamExecutionError("Legacy handoff cannot be upgraded: required real Salmon quant.sf/tx2gene artifacts are absent.")
    manifest["salmon"] = {
        "quant_sf": {sample: path.relative_to(run_dir).as_posix() for sample, path in quant.items()},
        "tx2gene": tx2gene.relative_to(run_dir).as_posix(),
        "transcript_counts": (output / "salmon.merged.transcript_counts.tsv").relative_to(run_dir).as_posix(),
        "contract_note": "L1 imports per-sample quant.sf with salmon.merged.tx2gene.tsv via tximport; estimated counts are not silently rounded by Python.",
    }
    manifest_path = run_dir / "handoff" / "upstream_manifest.yaml"
    manifest_path.write_text(yaml.safe_dump(manifest, sort_keys=False, allow_unicode=True), encoding="utf-8", newline="\n")
    return manifest


def prepare_l1(report: ValidationReport, *, run_id: str | None) -> PreparedL1:
    if not report.is_valid or report.config is None or report.loaded is None or report.metadata is None:
        raise DownstreamExecutionError("L1 is blocked because project validation failed.")
    config = report.config
    if config.input.type is InputType.RAW_COUNTS:
        if run_id is not None:
            raise DownstreamExecutionError("--run-id is only valid for FASTQ projects.")
        if report.counts is None:
            raise DownstreamExecutionError("Validated raw-count matrix is unavailable.")
        return PreparedL1(
            report.project_dir, report.project_dir / "downstream" / "l1", "raw_counts",
            report.counts.sample_ids, report.loaded.metadata_path, config.design.formula,
            {"source_type": "raw_counts", "counts": str(report.counts.path.resolve())},
        )
    if not run_id:
        raise DownstreamExecutionError("FASTQ L1 analysis requires explicit --run-id; a latest run is never selected implicitly.")
    run_dir = report.project_dir / "runs" / run_id
    state = _read_yaml(run_dir / "run_state.json", "upstream run state") if (run_dir / "run_state.json").is_file() else None
    if state is None or state.get("status") != "SUCCESS":
        raise DownstreamExecutionError(f"Selected upstream run is not a successful persisted run: {run_id}")
    manifest = _read_yaml(run_dir / "handoff" / "upstream_manifest.yaml", "upstream handoff manifest")
    manifest = _upgrade_salmon_handoff(run_dir, manifest)
    salmon = manifest.get("salmon")
    if not isinstance(salmon, dict):
        raise DownstreamExecutionError("Handoff does not provide a Salmon input contract.")
    quant = salmon.get("quant_sf")
    samples = manifest.get("samples")
    if not isinstance(quant, dict) or not isinstance(samples, list):
        raise DownstreamExecutionError("Handoff Salmon quant.sf mapping is invalid.")
    if set(quant) != set(samples):
        raise DownstreamExecutionError("Handoff Salmon quant.sf sample set disagrees with its recorded sample list.")
    quant_paths = {sample: str(_relative_existing(run_dir, value, f"salmon.quant_sf.{sample}")) for sample, value in sorted(quant.items())}
    tx2gene = _relative_existing(run_dir, salmon.get("tx2gene"), "salmon.tx2gene")
    frozen_metadata = run_dir / "frozen" / "metadata.csv"
    if not frozen_metadata.is_file():
        raise DownstreamExecutionError("Selected run has no frozen metadata.csv.")
    # A run is immutable: downstream uses its frozen project formula and metadata.
    frozen_project = _read_yaml(run_dir / "frozen" / "project.yaml", "frozen project configuration")
    formula = frozen_project.get("design", {}).get("formula") if isinstance(frozen_project.get("design"), dict) else None
    if not isinstance(formula, str):
        raise DownstreamExecutionError("Selected run frozen project configuration has no design.formula.")
    return PreparedL1(
        report.project_dir, run_dir / "downstream" / "l1", "salmon_tximport", tuple(sorted(samples)),
        frozen_metadata, formula,
        {"source_type": "salmon_tximport", "run_id": run_id, "quant_sf": quant_paths, "tx2gene": str(tx2gene)},
    )


def _write(path: Path, content: str) -> None:
    path.write_text(content, encoding="utf-8", newline="\n")


def _state(path: Path, **value: Any) -> None:
    _write(path, json.dumps(value, indent=2, sort_keys=True) + "\n")


def _validate_outputs(prepared: PreparedL1) -> dict[str, Any]:
    missing = [name for name in L1_OUTPUTS if not (prepared.output_dir / name).is_file()]
    if missing:
        raise DownstreamExecutionError("L1 backend did not produce required output(s): " + ", ".join(missing))
    with (prepared.output_dir / "normalized_counts.tsv").open(encoding="utf-8", newline="") as handle:
        rows = list(csv.reader(handle, delimiter="\t"))
    if len(rows) < 2 or rows[0][0] != "gene_id" or tuple(rows[0][1:]) != prepared.sample_ids:
        raise DownstreamExecutionError("normalized_counts.tsv has an invalid feature/sample structure.")
    if any(len(row) != len(rows[0]) for row in rows[1:]):
        raise DownstreamExecutionError("normalized_counts.tsv has malformed rows.")
    summary = json.loads((prepared.output_dir / "backend_summary.json").read_text(encoding="utf-8"))
    if not isinstance(summary, dict):
        raise DownstreamExecutionError("L1 backend summary is invalid.")
    return summary


def _render_report(prepared: PreparedL1, summary: dict[str, Any]) -> str:
    flags = summary.get("flags", [])
    lines = [
        "# L1 expression-level QC", "", "Status: completed", "", "## Input", "",
        f"- Source: `{prepared.source_type}`", f"- Samples retained: {len(prepared.sample_ids)} (no automatic sample removal)",
        f"- Formula: `{prepared.formula}`", f"- Filter: {L1_FILTER['rule']}", "", "## Outputs", "",
        "- DESeq2 size-factor normalization and normalized counts", "- blind VST, PCA, and sample Pearson correlation", "- No differential-expression hypothesis testing, DEG calling, or enrichment was run.", "",
        "## QC flags", "",
    ]
    lines.extend(f"- {flag}" for flag in flags) if flags else lines.append("- None")
    lines.append("")
    return "\n".join(lines)


def execute_l1(prepared: PreparedL1) -> L1Result:
    """Execute the bounded R backend and preserve logs/state on any failure."""
    _require_r()
    output = prepared.output_dir
    output.mkdir(parents=True, exist_ok=True)
    logs = output / "logs"
    logs.mkdir(exist_ok=True)
    state_path = output / "l1_state.json"
    started = datetime.now(UTC).replace(microsecond=0).isoformat().replace("+00:00", "Z")
    _state(state_path, status="RUNNING", started_at=started, completed_at=None, source_type=prepared.source_type)
    backend_config = {
        **prepared.config, "metadata": str(prepared.metadata_path.resolve()), "formula": prepared.formula,
        "samples": list(prepared.sample_ids), "output_dir": str(output.resolve()), "filter": L1_FILTER,
    }
    config_path = output / "backend_config.json"
    _write(config_path, json.dumps(backend_config, indent=2, sort_keys=True) + "\n")
    _write(output / "input_provenance.yaml", yaml.safe_dump({"input": prepared.config, "metadata": str(prepared.metadata_path), "formula": prepared.formula, "filter": L1_FILTER}, sort_keys=False))
    script = resources.files("rnaseq.r").joinpath("l1_analysis.R")
    command = ["Rscript", str(script), "--config", str(config_path)]
    try:
        with (logs / "r.stdout.log").open("w", encoding="utf-8", newline="\n") as stdout, (logs / "r.stderr.log").open("w", encoding="utf-8", newline="\n") as stderr:
            result = subprocess.run(command, stdout=stdout, stderr=stderr, check=False)
        if result.returncode != 0:
            raise DownstreamExecutionError(f"L1 R backend failed with return code {result.returncode}. Logs: {logs}")
        summary = _validate_outputs(prepared)
    except (OSError, json.JSONDecodeError, DownstreamExecutionError) as exc:
        _state(state_path, status="FAILED", started_at=started, completed_at=datetime.now(UTC).replace(microsecond=0).isoformat().replace("+00:00", "Z"), error=str(exc))
        if isinstance(exc, DownstreamExecutionError):
            raise
        raise DownstreamExecutionError(f"Unable to execute L1 R backend: {exc}. Logs: {logs}") from exc
    _write(
        output / "filtering_summary.yaml",
        yaml.safe_dump(
            {
                "rule": L1_FILTER,
                "genes_input": summary.get("genes_input"),
                "genes_removed_all_zero": summary.get("genes_removed_all_zero"),
                "genes_removed_low_total": summary.get("genes_removed_low_total"),
                "genes_retained": summary.get("genes_retained"),
                "automatic_sample_removal": False,
            },
            sort_keys=False,
        ),
    )
    _write(output / "l1_report.md", _render_report(prepared, summary))
    _state(state_path, status="SUCCESS", started_at=started, completed_at=datetime.now(UTC).replace(microsecond=0).isoformat().replace("+00:00", "Z"), source_type=prepared.source_type)
    return L1Result(output, state_path, summary)
