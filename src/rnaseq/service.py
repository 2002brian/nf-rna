"""Case/run lifecycle and client-facing delivery for the production workflow.

This module deliberately contains no statistics.  It freezes validated inputs,
starts independent Nextflow workflows, and curates their artifacts for a client.
"""

from __future__ import annotations

import csv
import hashlib
import json
import os
import platform
import re
import shutil
import subprocess
import sys
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Any
from zoneinfo import ZoneInfo

import yaml

from rnaseq.errors import ExecutionPreflightError, UpstreamExecutionError
from rnaseq.execution import (
    CONTROL_PLANE_IMAGE,
    HISAT2_WORKFLOW,
    CONTAINER_PROFILE,
    LOCAL_PROFILE,
    PreparedRun,
    EffectiveResourceBudget,
    ExecutionWorkspace,
    _validate_custom_reference_files,
    build_nextflow_command,
    build_hisat2_featurecounts_command,
    classify_execution_failure,
    check_container_runtime,
    check_docker,
    check_nextflow,
    downstream_docker_user_mapping,
    generate_handoff_manifest,
    generate_hisat2_featurecounts_handoff,
    inspect_container_image,
    effective_resource_budget,
    resolved_upstream_implementation,
    nfcore_runtime_params,
    prepare_execution_workspace,
    project_execution_budget,
    RESOURCE_CONTRACTS,
    render_local_resource_config,
    require_fresh_plan,
    resolve_execution_workspace,
    runtime_snapshot,
    validate_effective_resource_budget,
    validate_local_execution_budget,
    detect_local_resource_capacity,
    ResourceContract,
    LOCAL_RESOURCE_CEILING,
)
from rnaseq.models import FastqPreprocessing, InputType, PIPELINE_VERSION, Preset, production_enrichment_backends
from rnaseq.hisat2_featurecounts import FASTP_IMAGE, FASTP_VERSION, FASTQC_IMAGE, FASTQC_VERSION, HISAT2_IMAGE, HISAT2_VERSION, MULTIQC_IMAGE, MULTIQC_VERSION, SAMTOOLS_IMAGE, SAMTOOLS_VERSION, SUBREAD_IMAGE, SUBREAD_VERSION
from rnaseq.planner import pairing_contract, render_manifest
from rnaseq.validators import ValidationReport

CASE_ID_PATTERN = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]{0,127}$")
RUN_ID_PATTERN = re.compile(r"^[0-9]{8}-[0-9]{6}\+[0-9]{4}(?:-[0-9]{2,3})?$")
TAIPEI = ZoneInfo("Asia/Taipei")
FINAL_STATES = {"SUCCESS", "FAILED"}


@dataclass(frozen=True)
class CaseRun:
    case_id: str
    run_id: str
    run_dir: Path
    started_at: str

    @property
    def state_path(self) -> Path:
        return self.run_dir / "run_state.json"


@dataclass(frozen=True)
class FrozenInputs:
    samplesheet: Path | None
    contract: Path
    manifest: Path
    upstream_params: Path | None
    upstream_config: Path | None
    reference_paths: dict[str, Path]


@dataclass(frozen=True)
class ResolvedDownstreamInputs:
    """Container-stageable analysis inputs, distinct from frozen provenance paths."""

    root: Path
    manifest: Path
    source_type: str
    samples: tuple[str, ...]


def validate_case_id(value: str) -> str:
    """Validate a human supplied logical case identifier before using it in paths."""

    if not isinstance(value, str) or not CASE_ID_PATTERN.fullmatch(value):
        raise ExecutionPreflightError(
            "case ID must be 1-128 filesystem-safe characters: letters, digits, '.', '_' or '-'."
        )
    if value in {".", ".."} or value.endswith(".") or ".." in value or any(ord(char) < 32 for char in value):
        raise ExecutionPreflightError("case ID must not contain path traversal or control characters.")
    return value


def taipei_run_timestamp(moment: datetime | None = None) -> str:
    """Return the documented, sortable execution identifier in Asia/Taipei time."""

    value = moment or datetime.now(TAIPEI)
    if value.tzinfo is None:
        value = value.replace(tzinfo=TAIPEI)
    else:
        value = value.astimezone(TAIPEI)
    value = value.replace(microsecond=0)
    return value.strftime("%Y%m%d-%H%M%S%z")


def delivery_filename(name: str, run_id: str) -> str:
    """Return a deterministic client filename containing only the run date."""

    if not RUN_ID_PATTERN.fullmatch(run_id):
        raise ExecutionPreflightError(f"Cannot derive a delivery date from unsafe run ID: {run_id!r}.")
    source = Path(name)
    return f"{source.stem}_{run_id[:8]}{source.suffix}"


def _write_text(path: Path, value: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(value, encoding="utf-8", newline="\n")


def _write_yaml(path: Path, value: dict[str, Any]) -> None:
    _write_text(path, yaml.safe_dump(value, sort_keys=False, allow_unicode=True))


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _salmon_mapping_contract(run_dir: Path, value: object) -> tuple[Path, dict[str, str]]:
    """Resolve both current structured and historical string tx2gene contracts."""

    if isinstance(value, str):
        path = _safe_existing_under(
            run_dir, value, "salmon.tx2gene", allowed_root=run_dir / "upstream" / "nfcore_rnaseq"
        )
        return path, {
            "mapping_type": "historical_ordinary",
            "role": "Historical ordinary GTF-derived tx2gene mapping; no augmented self-mapping claim.",
            "sha256": _sha256(path),
        }
    if not isinstance(value, dict):
        raise UpstreamExecutionError("Frozen Salmon handoff has no valid tx2gene mapping contract.")
    path = _safe_existing_under(
        run_dir, value.get("path"), "salmon.tx2gene.path", allowed_root=run_dir / "upstream" / "nfcore_rnaseq"
    )
    mapping_type = value.get("mapping_type")
    role = value.get("role")
    expected = value.get("sha256")
    if mapping_type not in {"nfcore_tx2gene_augmented", "historical_ordinary"}:
        raise UpstreamExecutionError(f"Unsupported Salmon tx2gene mapping type: {mapping_type!r}.")
    if not isinstance(role, str) or not role.strip():
        raise UpstreamExecutionError("Salmon tx2gene mapping contract must record its role.")
    observed = _sha256(path)
    if expected != observed:
        raise UpstreamExecutionError(
            f"Salmon tx2gene checksum mismatch: expected {expected!r}, observed {observed}."
        )
    return path, {"mapping_type": mapping_type, "role": role, "sha256": observed}


def _is_appledouble(path: Path) -> bool:
    return path.name.startswith("._")


def _read_json_mapping(path: Path, label: str) -> dict[str, Any]:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError) as exc:
        raise UpstreamExecutionError(f"Cannot read {label}: {exc}") from exc
    if not isinstance(value, dict):
        raise UpstreamExecutionError(f"{label} must be a JSON object.")
    return value


def _read_yaml_mapping(path: Path, label: str) -> dict[str, Any]:
    try:
        value = yaml.safe_load(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, yaml.YAMLError) as exc:
        raise UpstreamExecutionError(f"Cannot read {label}: {exc}") from exc
    if not isinstance(value, dict):
        raise UpstreamExecutionError(f"{label} must be a YAML mapping.")
    return value


def _safe_existing_under(root: Path, configured: object, label: str, *, allowed_root: Path | None = None) -> Path:
    if not isinstance(configured, str) or not configured:
        raise UpstreamExecutionError(f"Downstream handoff is missing {label}.")
    candidate = (root / configured).resolve()
    try:
        candidate.relative_to(root.resolve())
    except ValueError as exc:
        raise UpstreamExecutionError(f"Downstream handoff {label} escapes immutable run directory.") from exc
    if allowed_root is not None:
        try:
            candidate.relative_to(allowed_root.resolve())
        except ValueError as exc:
            raise UpstreamExecutionError(
                f"Downstream handoff {label} is outside the immutable upstream output: {candidate}"
            ) from exc
    if _is_appledouble(candidate):
        raise UpstreamExecutionError(f"Downstream handoff {label} selected an AppleDouble artifact: {candidate}")
    if not candidate.is_file():
        raise UpstreamExecutionError(f"Downstream handoff {label} is unavailable: {candidate}")
    try:
        with candidate.open("rb") as handle:
            handle.read(1)
    except OSError as exc:
        raise UpstreamExecutionError(f"Downstream handoff {label} is unreadable: {candidate}: {exc}") from exc
    return candidate


def _frozen_sample_ids(metadata: Path) -> tuple[str, ...]:
    try:
        with metadata.open(encoding="utf-8", newline="") as handle:
            rows = list(csv.DictReader(handle))
    except (OSError, UnicodeError, csv.Error) as exc:
        raise UpstreamExecutionError(f"Frozen metadata is unreadable: {metadata}: {exc}") from exc
    samples = [row.get("sample_id") for row in rows]
    if not samples or any(not isinstance(sample, str) or not sample for sample in samples) or len(set(samples)) != len(samples):
        raise UpstreamExecutionError("Frozen metadata must contain unique, nonblank sample_id values for downstream staging.")
    return tuple(samples)


def resolve_downstream_inputs(run: CaseRun) -> ResolvedDownstreamInputs:
    """Resolve immutable provenance into narrow, stageable Nextflow file inputs.

    Paths in the frozen contract are provenance paths. Processes consume this
    copied bundle through a Nextflow ``path`` input rather than dereferencing
    arbitrary host paths from the contract.
    """

    frozen = run.run_dir / "frozen"
    contract_path = frozen / "downstream_contract.json"
    contract = _read_json_mapping(contract_path, "frozen downstream contract")
    source = contract.get("source")
    if not isinstance(source, dict):
        raise UpstreamExecutionError("Frozen downstream contract has no source mapping.")
    source_type = source.get("type")
    if source_type not in {"raw_counts", "salmon_tximport", "featurecounts_raw_counts"}:
        raise UpstreamExecutionError(f"Unsupported downstream source type: {source_type!r}.")
    project = frozen / "project.yaml"
    metadata = frozen / "metadata.csv"
    contrasts = frozen / "contrasts.csv"
    for label, path in (("frozen project configuration", project), ("frozen metadata", metadata), ("frozen contrasts", contrasts)):
        if not path.is_file():
            raise UpstreamExecutionError(f"Selected run has no {label}: {path}")
    samples = _frozen_sample_ids(metadata)
    destination = run.run_dir / "downstream_inputs"
    if destination.exists():
        raise UpstreamExecutionError(f"Refusing to overwrite immutable downstream execution inputs: {destination}")
    temporary = run.run_dir / f".downstream-inputs-{os.urandom(8).hex()}"
    try:
        temporary.mkdir()
        for source_path, target_name in ((project, "project.yaml"), (metadata, "metadata.csv"), (contrasts, "contrasts.csv")):
            _copy_snapshot(source_path, temporary / target_name)
        execution_manifest = frozen / "execution_manifest.yaml"
        if execution_manifest.is_file():
            _copy_snapshot(execution_manifest, temporary / "execution_manifest.yaml")
        execution_source: dict[str, Any]
        if source_type == "raw_counts":
            expected = frozen / "input" / "counts.csv"
            observed = Path(source.get("counts", "")).resolve() if isinstance(source.get("counts"), str) else None
            if observed != expected.resolve():
                raise UpstreamExecutionError("Frozen raw-count provenance does not point to this run's frozen count matrix.")
            if not expected.is_file():
                raise UpstreamExecutionError(f"Frozen count matrix is unavailable: {expected}")
            _copy_snapshot(expected, temporary / "source" / "counts.csv")
            execution_source = {"type": "raw_counts", "counts": "source/counts.csv"}
        elif source_type == "featurecounts_raw_counts":
            expected_handoff = frozen / "upstream_handoff_manifest.yaml"
            observed_handoff = Path(source.get("upstream_handoff", "")).resolve() if isinstance(source.get("upstream_handoff"), str) else None
            if observed_handoff != expected_handoff.resolve():
                raise UpstreamExecutionError("Frozen featureCounts provenance does not point to this run's frozen upstream handoff manifest.")
            handoff = _read_yaml_mapping(expected_handoff, "frozen upstream handoff manifest")
            featurecounts = handoff.get("featurecounts")
            if not isinstance(featurecounts, dict):
                raise UpstreamExecutionError("Frozen upstream handoff has no featureCounts artifact mapping.")
            matrix = _safe_existing_under(run.run_dir, featurecounts.get("canonical_matrix"), "featurecounts.canonical_matrix", allowed_root=run.run_dir / "upstream" / "hisat2_featurecounts")
            with matrix.open(encoding="utf-8", newline="") as handle:
                header = next(csv.reader(handle), [])
            if header != ["gene_id", *samples]:
                raise UpstreamExecutionError("Frozen featureCounts matrix sample IDs disagree with frozen metadata.")
            _copy_snapshot(matrix, temporary / "source" / "canonical_counts.csv")
            execution_source = {"type": "featurecounts_raw_counts", "counts": "source/canonical_counts.csv"}
        else:
            expected_handoff = frozen / "upstream_handoff_manifest.yaml"
            observed_handoff = Path(source.get("upstream_handoff", "")).resolve() if isinstance(source.get("upstream_handoff"), str) else None
            if observed_handoff != expected_handoff.resolve():
                raise UpstreamExecutionError("Frozen Salmon provenance does not point to this run's frozen upstream handoff manifest.")
            handoff = _read_yaml_mapping(expected_handoff, "frozen upstream handoff manifest")
            handoff_samples = handoff.get("samples")
            salmon = handoff.get("salmon")
            if not isinstance(handoff_samples, list) or not all(isinstance(item, str) and item for item in handoff_samples):
                raise UpstreamExecutionError("Frozen upstream handoff has no valid Salmon sample list.")
            if tuple(sorted(handoff_samples)) != samples:
                raise UpstreamExecutionError(
                    "Frozen metadata sample IDs disagree with upstream Salmon handoff samples: "
                    f"metadata={list(samples)!r}; handoff={sorted(handoff_samples)!r}."
                )
            if not isinstance(salmon, dict) or not isinstance(salmon.get("quant_sf"), dict):
                raise UpstreamExecutionError("Frozen upstream handoff has no Salmon quant.sf mapping.")
            quant = salmon["quant_sf"]
            if set(quant) != set(samples):
                raise UpstreamExecutionError(
                    "Frozen upstream handoff Salmon quant.sf sample set disagrees with metadata samples."
                )
            upstream_root = run.run_dir / "upstream" / "nfcore_rnaseq"
            tx2gene, mapping = _salmon_mapping_contract(run.run_dir, salmon.get("tx2gene"))
            staged_name = (
                "salmon.merged.tx2gene_augmented.tsv"
                if mapping["mapping_type"] == "nfcore_tx2gene_augmented"
                else "salmon.merged.tx2gene.tsv"
            )
            _copy_snapshot(tx2gene, temporary / "source" / staged_name)
            staged_quant: dict[str, str] = {}
            for index, sample in enumerate(samples, start=1):
                original = _safe_existing_under(
                    run.run_dir, quant[sample], f"salmon.quant_sf.{sample}", allowed_root=upstream_root
                )
                relative = f"source/quants/{index:04d}/quant.sf"
                _copy_snapshot(original, temporary / relative)
                staged_quant[sample] = relative
            execution_source = {
                "type": "salmon_tximport",
                "tx2gene": f"source/{staged_name}",
                "tx2gene_mapping": mapping,
                "quant_sf": staged_quant,
            }
        execution_manifest = {
            "schema_version": "1.0",
            "project_config": "project.yaml",
            "metadata": "metadata.csv",
            "contrasts": "contrasts.csv",
            "samples": list(samples),
            "source": execution_source,
        }
        _write_text(temporary / "execution_inputs.json", json.dumps(execution_manifest, indent=2, sort_keys=True) + "\n")
        temporary.replace(destination)
    except Exception:
        shutil.rmtree(temporary, ignore_errors=True)
        raise
    return ResolvedDownstreamInputs(destination, destination / "execution_inputs.json", source_type, samples)


def _resolved_delivery_root(delivery: Path, *, run_dir: Path | None = None) -> Path:
    """Resolve exactly one delivery root without allowing symlink escape."""

    if delivery.is_symlink():
        raise OSError(f"delivery root must not be a symlink: {delivery}")
    root = delivery.resolve(strict=True)
    if not root.is_dir():
        raise OSError(f"delivery root is not a directory: {delivery}")
    if run_dir is not None:
        try:
            root.relative_to(run_dir.resolve(strict=True))
        except ValueError as exc:
            raise OSError(f"delivery root escapes its immutable run directory: {delivery}") from exc
    return root


def _delivery_appledouble_entries(root: Path) -> list[Path]:
    """List AppleDouble entries under one already-resolved delivery root."""

    entries: list[Path] = []
    for current, directories, files in os.walk(root, followlinks=False):
        parent = Path(current)
        parent.relative_to(root)
        for name in [*directories, *files]:
            if name.startswith("._"):
                candidate = parent / name
                candidate.relative_to(root)
                entries.append(candidate.relative_to(root))
    return sorted(entries)


def _assert_delivery_appledouble_free(root: Path) -> None:
    remaining = _delivery_appledouble_entries(root)
    if remaining:
        raise OSError(
            "Delivery finalization failed: AppleDouble entries remain after sanitization: "
            + ", ".join(path.as_posix() for path in remaining)
        )


def _sanitize_delivery_appledouble(delivery: Path) -> None:
    """Remove only AppleDouble entries from one resolved delivery tree.

    Copying allowed artifacts to a macOS or external filesystem can create
    sidecars *after* discovery filtered their source counterparts. This final
    pass walks the resolved delivery root without following symlinks; a
    sidecar symlink is unlinked rather than traversed.
    """

    root = _resolved_delivery_root(delivery)

    def fail_walk(error: OSError) -> None:
        raise error

    for current, directories, files in os.walk(root, topdown=False, followlinks=False, onerror=fail_walk):
        parent = Path(current)
        parent.relative_to(root)
        for name in sorted([*files, *directories]):
            if not name.startswith("._"):
                continue
            candidate = parent / name
            candidate.relative_to(root)
            if candidate.is_symlink() or candidate.is_file():
                candidate.unlink()
            elif candidate.is_dir():
                # shutil.rmtree does not follow symlinks; the explicit check
                # above makes that boundary clear for a directory sidecar.
                shutil.rmtree(candidate)

    _assert_delivery_appledouble_free(root)


def create_case_run(report: ValidationReport, case_id: str, *, moment: datetime | None = None) -> CaseRun:
    """Allocate a new immutable case/run directory without ever reusing one."""

    if not report.is_valid:
        raise ExecutionPreflightError("Execution is blocked because validation failed.")
    case_id = validate_case_id(case_id)
    stamp = taipei_run_timestamp(moment)
    case_dir = report.project_dir / "runs" / case_id
    case_dir.mkdir(parents=True, exist_ok=True)
    for index in range(0, 1000):
        run_id = stamp if index == 0 else f"{stamp}-{index:02d}"
        run_dir = case_dir / run_id
        try:
            run_dir.mkdir()
        except FileExistsError:
            continue
        for child in ("frozen", "upstream", "downstream", "logs", "provenance", "handoff", "delivery"):
            (run_dir / child).mkdir()
        started_at = datetime.now(TAIPEI).replace(microsecond=0).isoformat()
        run = CaseRun(case_id, run_id, run_dir, started_at)
        _write_state(run, "CREATED", command=None)
        return run
    raise ExecutionPreflightError("Unable to allocate a unique run directory without overwriting an existing run.")


def _write_state(run: CaseRun, status: str, **values: Any) -> None:
    current: dict[str, Any] = {}
    if run.state_path.is_file():
        try:
            current = json.loads(run.state_path.read_text(encoding="utf-8"))
        except json.JSONDecodeError:
            current = {}
    current.update(
        {
            "case_id": run.case_id,
            "run_id": run.run_id,
            "timezone": "Asia/Taipei",
            "started_at": current.get("started_at") or run.started_at,
            "status": status,
            **values,
        }
    )
    if status in FINAL_STATES:
        current["completed_at"] = datetime.now(TAIPEI).replace(microsecond=0).isoformat()
    else:
        current.setdefault("completed_at", None)
    _write_text(run.state_path, json.dumps(current, indent=2, sort_keys=True) + "\n")


def _copy_snapshot(source: Path, destination: Path) -> None:
    if _is_appledouble(source):
        return
    destination.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(source, destination)


def _copy_delivery_artifact(source: Path, destination: Path) -> None:
    """Copy delivery payload bytes without propagating filesystem metadata.

    AppleDouble sidecars on external macOS volumes commonly represent source
    extended attributes. Delivery does not need those attributes, so use a
    data-only copy and keep the final sanitizer as the acceptance backstop.
    """

    if _is_appledouble(source):
        return
    destination.parent.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(source, destination)


def _staged_fastq_samplesheet(report: ValidationReport, fastq_root: Path) -> str:
    assert report.fastq is not None and report.config is not None
    output: list[list[str]] = [["sample", "fastq_1", "fastq_2", "strandedness"]]
    for record in report.fastq.records:
        first = fastq_root / record.fastq_1.name
        second = fastq_root / record.fastq_2.name if record.fastq_2 else None
        output.append([record.sample_id, str(first.resolve()), str(second.resolve()) if second else "", report.config.upstream.strandedness or "auto"])
    from io import StringIO

    text = StringIO(newline="")
    writer = csv.writer(text, lineterminator="\n")
    writer.writerows(output)
    return text.getvalue()


def freeze_case_inputs(
    report: ValidationReport, run: CaseRun, *, profile: str, command: list[str],
    resources: EffectiveResourceBudget | None = None,
) -> FrozenInputs:
    """Copy every mutable analysis input into the run before any compute starts."""

    assert report.loaded is not None and report.config is not None
    frozen = run.run_dir / "frozen"
    for source, target in (
        (report.loaded.config_path, frozen / "project.yaml"),
        (report.loaded.metadata_path, frozen / "metadata.csv"),
        (report.loaded.contrasts_path, frozen / "contrasts.csv"),
    ):
        _copy_snapshot(source, target)

    manifest = yaml.safe_load(render_manifest(report))
    assert isinstance(manifest, dict)
    manifest_path = frozen / "input_manifest.yaml"
    _write_yaml(manifest_path, manifest)
    staged_input = frozen / "input"
    reference_paths: dict[str, Path] = {}
    if report.config.reference.source == "local":
        assert report.local_reference is not None
        # The managed reference can be large.  Freeze its manifest and resolved
        # checksum identity rather than copying multi-GB immutable assets into a
        # client run; execution re-verifies those checksum-bound source files.
        manifest_snapshot = frozen / "reference" / "reference_manifest.yaml"
        _copy_snapshot(report.local_reference.manifest_path, manifest_snapshot)
        method = report.config.upstream.quantification.method if report.config.input.type is InputType.FASTQ and report.config.upstream.quantification else "salmon"
        reference_paths = {option.removeprefix("--"): path for option, path in (
            report.local_reference.nfcore_arguments() if method == "salmon" else report.local_reference.hisat2_arguments()
        )}
    elif report.config.reference.source == "custom":
        for key in ("fasta", "gtf", "transcript_fasta", "salmon_index", "hisat2_index", "hisat2_splice_sites"):
            configured = getattr(report.config.reference, key)
            if configured is None:
                continue
            source_reference = (report.project_dir / configured).resolve()
            if source_reference.is_dir():
                # Indexes are immutable, potentially large assets. They stay at
                # their checksum-validated project path; the frozen manifest
                # records their exact configured identity.
                reference_paths[key] = source_reference
            else:
                target_reference = frozen / "reference" / source_reference.name
                _copy_snapshot(source_reference, target_reference)
                reference_paths[key] = target_reference.resolve()
    samplesheet: Path | None = None
    source: dict[str, Any]
    if report.config.input.type is InputType.RAW_COUNTS:
        assert report.counts is not None
        counts = staged_input / "counts.csv"
        _copy_snapshot(report.counts.path, counts)
        source = {"type": "raw_counts", "construction_method": "DESeqDataSetFromMatrix", "counts": str(counts.resolve())}
    else:
        assert report.fastq is not None
        fastq_dir = staged_input / "fastq"
        staged: set[str] = set()
        for record in report.fastq.records:
            for item in (record.fastq_1, record.fastq_2):
                if item is None or _is_appledouble(item):
                    continue
                if item.name in staged:
                    continue
                staged.add(item.name)
                _copy_snapshot(item, fastq_dir / item.name)
        samplesheet = frozen / "samplesheet.csv"
        _write_text(samplesheet, _staged_fastq_samplesheet(report, fastq_dir))
        method = report.config.upstream.quantification.method if report.config.upstream.quantification else "salmon"
        source = {
            "type": "salmon_tximport" if method == "salmon" else "featurecounts_raw_counts",
            "construction_method": "DESeqDataSetFromTximport" if method == "salmon" else "DESeqDataSetFromMatrix",
            "upstream_handoff": None,
            "fastq_preprocessing": report.config.input.preprocessing.value,
            "skip_trimming": report.config.input.preprocessing is FastqPreprocessing.PRETRIMMED,
        }

    params: Path | None = None
    runtime: Path | None = None
    if report.config.input.type is InputType.FASTQ:
        params = frozen / "nfcore.params.json"
        _write_text(params, json.dumps(nfcore_runtime_params(report), sort_keys=True) + "\n")
        # This single path-free config is passed to Salmon, HISAT2/featureCounts,
        # and the first-party downstream workflow for every local FASTQ run.
        runtime = frozen / "nfcore.local.config"
        resolved = resources or effective_resource_budget(
            runtime_snapshot(report.config.runtime.control_plane_image), project_execution_budget(report.config)
        )
        _write_text(runtime, render_local_resource_config(ResourceContract(
            "EFFECTIVE_LOCAL", resolved.effective_cpus, resolved.effective_memory_gib, LOCAL_RESOURCE_CEILING.time_hours
        )))

    execution = {
        "case_id": run.case_id,
        "run_id": run.run_id,
        "timezone": "Asia/Taipei",
        "profile": profile,
        "execution_budget": report.config.execution.model_dump(),
        "command": command,
        "pipeline": {"name": report.config.project.pipeline, "version": PIPELINE_VERSION},
        "upstream_implementation": (
            resolved_upstream_implementation(report)
            if report.config.input.type is InputType.FASTQ else None
        ),
        "reference": (
            report.local_reference.provenance()
            if report.local_reference is not None
            else report.config.reference.model_dump()
        ),
        "fastq_preprocessing": (
            report.config.input.preprocessing.value
            if report.config.input.type is InputType.FASTQ else None
        ),
        "skip_trimming": (
            report.config.input.preprocessing is FastqPreprocessing.PRETRIMMED
            if report.config.input.type is InputType.FASTQ else None
        ),
        "nfcore_preprocessing_arguments": (
            ["--skip_trimming"]
            if report.config.input.type is InputType.FASTQ
            and report.config.input.preprocessing is FastqPreprocessing.PRETRIMMED else []
        ),
        "nfcore_runtime_params": (
            nfcore_runtime_params(report)
            if report.config.input.type is InputType.FASTQ else None
        ),
        "local_nextflow_config": (
            {"path": runtime.relative_to(run.run_dir).as_posix(), "sha256": _sha256(runtime)}
            if runtime is not None else None
        ),
    }
    _write_yaml(frozen / "execution_manifest.yaml", execution)
    contract = {
        "schema_version": "1.0",
        "case": {"id": run.case_id, "run_id": run.run_id, "timezone": "Asia/Taipei"},
        "source": source,
        "project_config": str((frozen / "project.yaml").resolve()),
        "metadata": str((frozen / "metadata.csv").resolve()),
        "contrasts": str((frozen / "contrasts.csv").resolve()),
        "input_manifest": str(manifest_path.resolve()),
        # The immutable project preset is the sole execution-level selector.
        # Downstream tasks must never infer L2 from contrasts or metadata.
        "analysis_level": report.config.project.preset.value,
        "analysis": report.config.analysis.model_dump(mode="json") if report.config.analysis else {"enrichment": []},
        "design": {
            "type": report.config.design.type.value,
            "formula": report.config.design.formula,
            **({"pair_id": report.config.design.pair_id, "pairing": pairing_contract(report)} if report.config.design.pair_id else {}),
        },
        # This is the immutable, normalized enrichment input for downstream
        # tasks.  Do not make those tasks re-read mutable project.yaml.
        "annotation": report.config.annotation.model_dump(mode="json") if report.config.annotation else None,
        "output_dir": str((run.run_dir / "downstream").resolve()),
    }
    contract_path = frozen / "downstream_contract.json"
    _write_text(contract_path, json.dumps(contract, indent=2, sort_keys=True) + "\n")
    return FrozenInputs(samplesheet, contract_path, manifest_path, params, runtime, reference_paths)


def _update_contract(contract_path: Path, mutate: dict[str, Any]) -> None:
    contract = json.loads(contract_path.read_text(encoding="utf-8"))
    contract.update(mutate)
    _write_text(contract_path, json.dumps(contract, indent=2, sort_keys=True) + "\n")


def _provenance(
    report: ValidationReport, run: CaseRun, *, profile: str, command: list[str], workspace: ExecutionWorkspace,
    resources: EffectiveResourceBudget | None = None,
) -> dict[str, Any]:
    git_commit: str | None = None
    try:
        source_root = Path(__file__).resolve().parents[2]
        result = subprocess.run(["git", "rev-parse", "HEAD"], cwd=source_root, capture_output=True, text=True, check=False)
        if result.returncode == 0:
            git_commit = result.stdout.strip() or None
    except OSError:
        pass
    nextflow = check_nextflow()
    requested_image = report.config.runtime.control_plane_image if report.config else CONTROL_PLANE_IMAGE
    runtime = runtime_snapshot(requested_image)
    resolved_resources = resources or effective_resource_budget(runtime, project_execution_budget(report.config))
    method = report.config.upstream.quantification.method if report.config and report.config.upstream.quantification else None
    def tool_identity(version: str, image: str) -> dict[str, object]:
        return {"version": version, **inspect_container_image(image)}

    control_plane = inspect_container_image(requested_image)
    source_root = Path(__file__).resolve().parents[2]
    workflow_hashes = {
        "workflow/main.nf": _sha256(source_root / "workflow" / "main.nf"),
        "workflow/hisat2_featurecounts.nf": _sha256(HISAT2_WORKFLOW),
    }
    local_config = run.run_dir / "frozen" / "nfcore.local.config"
    return {
        "case_id": run.case_id,
        "run_id": run.run_id,
        "timezone": "Asia/Taipei",
        "started_at": run.started_at,
        "pipeline_version": PIPELINE_VERSION,
        "git_commit": git_commit,
        "source_checkout": str(source_root),
        "workflow_sha256": workflow_hashes,
        "python_version": sys.version.split()[0],
        "nextflow_version": nextflow.detail if nextflow.state == "FOUND" else None,
        "upstream_implementation": (resolved_upstream_implementation(report) if report.config and report.config.input.type is InputType.FASTQ else None),
        "nfcore_rnaseq_version": (report.config.upstream.pipeline_version if method == "salmon" else None),
        "fastq_backend": method,
        "hisat2_featurecounts_tools": (
            {
                "hisat2": tool_identity(HISAT2_VERSION, HISAT2_IMAGE),
                "samtools": tool_identity(SAMTOOLS_VERSION, SAMTOOLS_IMAGE),
                "subread": tool_identity(SUBREAD_VERSION, SUBREAD_IMAGE),
                "fastqc": tool_identity(FASTQC_VERSION, FASTQC_IMAGE),
                "fastp": tool_identity(FASTP_VERSION, FASTP_IMAGE),
                "multiqc": tool_identity(MULTIQC_VERSION, MULTIQC_IMAGE),
            } if method == "hisat2_featurecounts" else None
        ),
        "container_runtime": "docker",
        "container_image": control_plane,
        "production_intended": bool(report.config and report.config.reference.acceptance == "production"),
        "runtime_resources": {
            **resolved_resources.as_dict(),
            "host_os": runtime.host_os,
            "host_architecture": runtime.host_architecture,
            "logical_cpus": runtime.logical_cpus,
            "host_memory_bytes": runtime.host_memory_bytes,
            "docker_architecture": runtime.docker_architecture,
            "docker_memory_bytes": runtime.docker_memory_bytes,
            "docker_version": runtime.docker_version,
            "control_plane_image_architecture": runtime.control_plane_image_architecture,
            "resource_profile": "M5_LOCAL_SMALL_MEDIUM_LARGE",
        },
        "frozen_local_nextflow_config": (
            {"path": local_config.relative_to(run.run_dir).as_posix(), "sha256": _sha256(local_config)}
            if local_config.is_file() else None
        ),
        "profile": profile,
        "command": command,
        "execution_root": str(workspace.root),
        "execution_launch_dir": str(workspace.launch_dir),
        "execution_work_dir": str(workspace.work_dir),
        "input_manifest_sha256": _sha256(run.run_dir / "frozen" / "input_manifest.yaml"),
        "design": {
            "type": report.config.design.type.value,
            "formula": report.config.design.formula,
            **({"pair_id": report.config.design.pair_id, "pairing": pairing_contract(report)} if report.config.design.pair_id else {}),
        } if report.config is not None else None,
        "frozen_project_sha256": _sha256(run.run_dir / "frozen" / "project.yaml"),
        "reference": (
            report.local_reference.provenance()
            if report.local_reference is not None
            else report.config.reference.model_dump() if report.config is not None else None
        ),
        "fastq_preprocessing": (
            report.config.input.preprocessing.value
            if report.config and report.config.input.type is InputType.FASTQ else None
        ),
        "skip_trimming": (
            report.config.input.preprocessing is FastqPreprocessing.PRETRIMMED
            if report.config and report.config.input.type is InputType.FASTQ else None
        ),
        "nfcore_preprocessing_arguments": (
            ["--skip_trimming"]
            if report.config and report.config.input.type is InputType.FASTQ
            and report.config.input.preprocessing is FastqPreprocessing.PRETRIMMED else []
        ),
        "nfcore_runtime_params": (
            nfcore_runtime_params(report)
            if report.config and report.config.input.type is InputType.FASTQ else None
        ),
    }


def write_downstream_observer_config(run: CaseRun) -> Path:
    """Freeze per-run Nextflow observer paths with safe same-run overwrites.

    Nextflow can invoke observers again while finalising a failed workflow.  Each
    immutable run has unique provenance paths, so allowing an observer to replace
    only its own partial file prevents secondary shutdown warnings without
    permitting any cross-run overwrite.
    """

    targets = {
        "trace": run.run_dir / "provenance" / "downstream.trace.txt",
        "report": run.run_dir / "provenance" / "downstream.report.html",
        "timeline": run.run_dir / "provenance" / "downstream.timeline.html",
        "dag": run.run_dir / "provenance" / "downstream.dag.html",
    }
    lines: list[str] = []
    for observer, target in targets.items():
        lines.extend([
            f"{observer} {{",
            "  enabled = true",
            f"  file = {json.dumps(str(target.resolve()))}",
            "  overwrite = true",
            "}",
            "",
        ])
    path = run.run_dir / "frozen" / "downstream.observers.config"
    _write_text(path, "\n".join(lines))
    return path


def write_downstream_docker_user_config(
    run: CaseRun, *, host_os: str | None = None, uid: int | None = None, gid: int | None = None,
) -> Path | None:
    """Freeze a Linux/WSL Docker user override for one downstream run.

    The override is deliberately a separate Nextflow config rather than an
    image change or a broad filesystem permission change.  It affects only the
    downstream task containers and only on Linux/WSL.  Numeric values are
    produced by :func:`downstream_docker_user_mapping`, so the rendered Groovy
    string cannot interpolate a user-controlled shell value.
    """

    mapping = downstream_docker_user_mapping(
        host_os=host_os or platform.system(), uid=uid, gid=gid,
    )
    if mapping is None:
        return None
    path = run.run_dir / "frozen" / "downstream.docker-user.config"
    _write_text(path, "docker {\n  runOptions = '--user " + mapping + "'\n}\n")
    return path


def write_downstream_runtime_config(run: CaseRun, image: str) -> Path:
    """Freeze the requested per-run downstream image instead of inheriting latest."""

    path = run.run_dir / "frozen" / "downstream.runtime.config"
    _write_text(path, f"process.container = {json.dumps(image)}\n")
    return path


def finalize_fastq_handoff(report: ValidationReport, run: CaseRun, contract: Path, *, reused_from: str | None = None) -> Path:
    """Validate and freeze the stable nf-core handoff boundary for downstream."""

    prepared = PreparedRun(report, "unavailable", "docker")
    method = report.config.upstream.quantification.method if report.config and report.config.upstream.quantification else "salmon"
    handoff = generate_handoff_manifest(prepared, run.run_dir) if method == "salmon" else generate_hisat2_featurecounts_handoff(report, run.run_dir)
    frozen_handoff = run.run_dir / "frozen" / "upstream_handoff_manifest.yaml"
    _copy_snapshot(handoff, frozen_handoff)
    _update_contract(
        contract,
        {"source": {"type": "salmon_tximport" if method == "salmon" else "featurecounts_raw_counts", "construction_method": "DESeqDataSetFromTximport" if method == "salmon" else "DESeqDataSetFromMatrix", "upstream_handoff": str(frozen_handoff.resolve()), "reused_from": reused_from}},
    )
    if method == "salmon":
        payload = _read_yaml_mapping(frozen_handoff, "frozen upstream handoff manifest")
        salmon = payload.get("salmon")
        if isinstance(salmon, dict):
            mapping_path, mapping = _salmon_mapping_contract(run.run_dir, salmon.get("tx2gene"))
            provenance_path = run.run_dir / "provenance" / "run_provenance.yaml"
            if provenance_path.is_file():
                provenance = _read_yaml_mapping(provenance_path, "run provenance")
                provenance["salmon_tx2gene"] = {
                    "path": str(mapping_path),
                    **mapping,
                }
                _write_yaml(provenance_path, provenance)
    return frozen_handoff


def build_downstream_nextflow_command(
    run: CaseRun, *, profile: str = LOCAL_PROFILE, work_dir: Path | None = None,
    observer_config: Path | None = None, docker_user_config: Path | None = None,
    execution_inputs: ResolvedDownstreamInputs | None = None, runtime_config: Path | None = None,
    local_resource_config: Path | None = None,
) -> list[str]:
    root = Path(__file__).resolve().parents[2]
    contract = json.loads((run.run_dir / "frozen" / "downstream_contract.json").read_text(encoding="utf-8"))
    analysis_level = contract.get("analysis_level")
    if analysis_level not in {"L1", "L2"}:
        raise ValueError("Frozen downstream contract analysis_level must be L1 or L2.")
    modules = production_enrichment_backends(contract.get("analysis", {}).get("enrichment", []))
    if analysis_level == "L1" and modules:
        raise ValueError("Frozen downstream contract cannot enable enrichment for analysis_level L1.")
    resolved_work_dir = work_dir or (resolve_execution_workspace(run.case_id, run.run_id).work_dir / "downstream")
    command = [
        "nextflow", "run", str(root / "workflow" / "main.nf"), "-c", str(root / "workflow" / "nextflow.config"), "-profile", CONTAINER_PROFILE if profile == LOCAL_PROFILE else profile,
    ]
    if observer_config is not None:
        command.extend(["-c", str(observer_config.resolve())])
    if runtime_config is not None:
        command.extend(["-c", str(runtime_config.resolve())])
    if local_resource_config is not None:
        command.extend(["-c", str(local_resource_config.resolve())])
    if docker_user_config is not None:
        command.extend(["-c", str(docker_user_config.resolve())])
    command.extend([
        "-work-dir", str(resolved_work_dir.resolve()),
        "--contract", str((run.run_dir / "frozen" / "downstream_contract.json").resolve()),
        "--inputs", str((execution_inputs.root if execution_inputs is not None else run.run_dir / "downstream_inputs").resolve()),
        "--outdir", str((run.run_dir / "downstream").resolve()),
        "--analysis_level", analysis_level,
    ])
    # Semantic absence must remain absent at the Nextflow boundary.  Passing an
    # empty CLI value becomes a literal invalid module under Nextflow parsing.
    if modules:
        command.extend(["--enrichment", ",".join(modules)])
    return command


def _run_command(command: list[str], *, cwd: Path, stdout_path: Path, stderr_path: Path) -> int:
    with stdout_path.open("w", encoding="utf-8", newline="\n") as stdout, stderr_path.open("w", encoding="utf-8", newline="\n") as stderr:
        return subprocess.run(command, cwd=cwd, stdout=stdout, stderr=stderr, check=False).returncode


def _validate_delivery_count_matrix(
    path: Path, samples: tuple[str, ...], *, integer_required: bool, nonnegative: bool = True, delimiter: str = ",",
) -> tuple[int, int]:
    """Validate a source matrix before exposing it as a delivery count artifact."""

    try:
        with path.open(encoding="utf-8", newline="") as handle:
            rows = list(csv.reader(handle, delimiter=delimiter))
    except (OSError, UnicodeError, csv.Error) as exc:
        raise UpstreamExecutionError(f"Delivery count matrix is unreadable: {path}: {exc}") from exc
    if len(rows) < 2 or rows[0] != ["gene_id", *samples]:
        raise UpstreamExecutionError(f"Delivery count matrix has unexpected gene/sample columns: {path}")
    genes: set[str] = set()
    for row in rows[1:]:
        if len(row) != len(rows[0]) or not row[0] or row[0] in genes:
            raise UpstreamExecutionError(f"Delivery count matrix has invalid gene identifiers: {path}")
        genes.add(row[0])
        for value in row[1:]:
            try:
                numeric = float(value)
            except ValueError as exc:
                raise UpstreamExecutionError(f"Delivery count matrix contains a non-numeric value: {path}") from exc
            if numeric != numeric or numeric in {float("inf"), float("-inf")} or (nonnegative and numeric < 0):
                raise UpstreamExecutionError(f"Delivery matrix contains a non-finite or invalid value: {path}")
            if integer_required and not numeric.is_integer():
                raise UpstreamExecutionError(f"Delivery raw-count matrix contains a non-integer value: {path}")
    return len(rows) - 1, len(samples)


def _delivery_count_source(run: CaseRun) -> tuple[Path, str, bool, str, str, str] | None:
    """Resolve the exact untransformed matrix used by this newly-created run."""

    contract = _read_json_mapping(run.run_dir / "frozen" / "downstream_contract.json", "frozen downstream contract")
    source = contract.get("source")
    if not isinstance(source, dict):
        return None
    source_type = source.get("type")
    construction = source.get("construction_method")
    if not isinstance(construction, str):
        raise UpstreamExecutionError("Frozen downstream contract has no source construction method.")
    if source_type == "raw_counts":
        return run.run_dir / "frozen" / "input" / "counts.csv", "raw_counts.csv", True, "integer_raw_counts", construction, ","
    if source_type == "featurecounts_raw_counts":
        handoff = _read_yaml_mapping(run.run_dir / "frozen" / "upstream_handoff_manifest.yaml", "frozen upstream handoff")
        featurecounts = handoff.get("featurecounts")
        if not isinstance(featurecounts, dict):
            raise UpstreamExecutionError("Frozen upstream handoff has no featureCounts count matrix.")
        matrix = _safe_existing_under(run.run_dir, featurecounts.get("canonical_matrix"), "featurecounts.canonical_matrix", allowed_root=run.run_dir / "upstream" / "hisat2_featurecounts")
        return matrix, "raw_counts.csv", True, "integer_raw_counts", construction, ","
    if source_type == "salmon_tximport":
        l1 = run.run_dir / "downstream" / "l1" / "source_counts.csv"
        if l1.is_file():
            return l1, "estimated_counts.csv", False, "salmon_estimated_counts", construction, ","
        # QC-only has no L1/tximport task. Deliver the canonical upstream
        # estimate as explicitly upstream-only, never as raw counts.
        handoff = _read_yaml_mapping(run.run_dir / "frozen" / "upstream_handoff_manifest.yaml", "frozen upstream handoff")
        salmon = handoff.get("salmon")
        if not isinstance(salmon, dict):
            return None
        gene_counts = handoff.get("gene_level_counts")
        if not isinstance(gene_counts, dict):
            return None
        matrix = _safe_existing_under(run.run_dir, gene_counts.get("path"), "salmon.gene_level_counts", allowed_root=run.run_dir / "upstream" / "nfcore_rnaseq")
        return matrix, "estimated_counts.csv", False, "salmon_estimated_counts", construction, "\t"
    raise UpstreamExecutionError(f"Unsupported delivery count source: {source_type!r}")


def _write_delivery_count_manifest(delivery: Path, artifacts: list[dict[str, object]]) -> None:
    _write_text(delivery / "counts" / "artifact_manifest.json", json.dumps({"schema_version": "1.0", "artifacts": artifacts}, indent=2, sort_keys=True) + "\n")


def _copy_matrix_as_csv(source: Path, destination: Path, *, delimiter: str) -> None:
    if delimiter == ",":
        _copy_delivery_artifact(source, destination)
        return
    with source.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.reader(handle, delimiter=delimiter))
    destination.parent.mkdir(parents=True, exist_ok=True)
    with destination.open("w", encoding="utf-8", newline="") as handle:
        csv.writer(handle, lineterminator="\n").writerows(rows)


def _delivery_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def assemble_delivery(run: CaseRun) -> Path:
    """Build a client-safe package from an explicit allowlist of final artifacts."""

    delivery = _resolved_delivery_root(run.run_dir / "delivery", run_dir=run.run_dir)
    figures = {".png": delivery / "figures" / "png", ".tif": delivery / "figures" / "tiff_300dpi", ".tiff": delivery / "figures" / "tiff_300dpi"}
    tables = delivery / "tables"
    counts_dir = delivery / "counts"
    for directory in [*figures.values(), tables, counts_dir, delivery / "methods_and_versions"]:
        directory.mkdir(parents=True, exist_ok=True)
    downstream = run.run_dir / "downstream"
    report = downstream / "report" / "report.html"
    if report.is_file():
        _copy_delivery_artifact(report, delivery / delivery_filename("report.html", run.run_id))
    artifacts: list[dict[str, object]] = []
    source = _delivery_count_source(run)
    samples = _frozen_sample_ids(run.run_dir / "frozen" / "metadata.csv")
    if source is not None:
        matrix, filename, integer_required, semantics, construction, delimiter = source
        rows, columns = _validate_delivery_count_matrix(matrix, samples, integer_required=integer_required, delimiter=delimiter)
        target = counts_dir / filename
        _copy_matrix_as_csv(matrix, target, delimiter=delimiter)
        artifacts.append({
            "role": "analysis_input" if matrix.is_relative_to(run.run_dir / "downstream") else "frozen_or_upstream_input",
            "filename": f"counts/{filename}", "source_type": _read_json_mapping(run.run_dir / "frozen" / "downstream_contract.json", "frozen downstream contract")["source"]["type"],
            "value_semantics": semantics, "normalized": False, "integer_required": integer_required,
            "deseq2_construction_method": construction, "sha256": _delivery_sha256(target), "rows": rows, "columns": columns,
            "ordered_sample_ids": list(samples), "gene_identifier_namespace": "gene_id",
        })
    vst = run.run_dir / "downstream" / "l1" / "vst.csv"
    if vst.is_file():
        rows, columns = _validate_delivery_count_matrix(vst, samples, integer_required=False, nonnegative=False)
        target = counts_dir / "vst.csv"
        _copy_delivery_artifact(vst, target)
        artifacts.append({
            "role": "visualization", "filename": "counts/vst.csv", "source_type": "downstream_l1",
            "value_semantics": "variance_stabilized_expression", "normalized": True, "integer_required": False,
            "deseq2_construction_method": "varianceStabilizingTransformation_or_vst", "sha256": _delivery_sha256(target),
            "rows": rows, "columns": columns, "ordered_sample_ids": list(samples), "gene_identifier_namespace": "gene_id",
        })
    if artifacts:
        _write_delivery_count_manifest(delivery, artifacts)
    for path in downstream.rglob("*") if downstream.is_dir() else ():
        if not path.is_file() or _is_appledouble(path):
            continue
        # These have a single, source-labelled canonical location under
        # delivery/counts.  Keeping an unlabelled duplicate under tables
        # invites using VST as a count matrix or treating Salmon estimates as
        # raw counts.
        if path in {downstream / "l1" / "source_counts.csv", downstream / "l1" / "vst.csv"}:
            continue
        suffix = path.suffix.lower()
        relative = path.relative_to(downstream)
        if suffix in figures:
            _copy_delivery_artifact(path, figures[suffix] / relative)
        elif suffix in {".tsv", ".csv"}:
            _copy_delivery_artifact(path, tables / relative)
    handoff = run.run_dir / "frozen" / "upstream_handoff_manifest.yaml"
    if handoff.is_file():
        payload = yaml.safe_load(handoff.read_text(encoding="utf-8"))
        multiqc = payload.get("multiqc", {}) if isinstance(payload, dict) else {}
        html = multiqc.get("html")
        if isinstance(html, str):
            source = run.run_dir / html
            if source.is_file():
                _copy_delivery_artifact(
                    source,
                    delivery / "multiqc" / delivery_filename("multiqc_report.html", run.run_id),
                )
    for name in ("execution_manifest.yaml", "input_manifest.yaml", "upstream_handoff_manifest.yaml"):
        source = run.run_dir / "frozen" / name
        if source.is_file():
            _copy_delivery_artifact(
                source,
                delivery / "methods_and_versions" / delivery_filename(name, run.run_id),
            )
    for source in (run.state_path, run.run_dir / "provenance" / "run_provenance.yaml"):
        if source.is_file():
            _copy_delivery_artifact(
                source,
                delivery / "methods_and_versions" / delivery_filename(source.name, run.run_id),
            )
    _write_text(
        delivery / delivery_filename("README.md", run.run_id),
        f"# RNA-seq delivery package\n\nCase: `{run.case_id}`  \nRun: `{run.run_id}`  \nTimezone: `Asia/Taipei`\n\nThis directory contains only curated client deliverables. Internal logs, work directories, caches and temporary artifacts remain outside this package.\n",
    )
    # This is deliberately the last mutation of delivery.  SUCCESS is recorded
    # only after this hard check confirms there are no AppleDouble entries.
    _sanitize_delivery_appledouble(delivery)
    _assert_delivery_appledouble_free(delivery)
    return delivery


def sanitize_completed_delivery(run_dir: Path) -> Path:
    """Repair only a completed run's delivery tree; never rerun analysis."""

    root = run_dir.resolve(strict=True)
    state_path = root / "run_state.json"
    state = _read_json_mapping(state_path, "run state")
    if state.get("status") != "SUCCESS":
        raise ExecutionPreflightError("Delivery sanitization is allowed only for a successful immutable run.")
    delivery = _resolved_delivery_root(root / "delivery", run_dir=root)
    _sanitize_delivery_appledouble(delivery)
    _assert_delivery_appledouble_free(delivery)
    return delivery


def prepare_service_run(report: ValidationReport, *, profile: str) -> EffectiveResourceBudget:
    if profile != LOCAL_PROFILE:
        raise ExecutionPreflightError("Only '--profile local' is implemented; server profiles are configuration placeholders.")
    if not report.is_valid or report.config is None:
        raise ExecutionPreflightError("Execution is blocked because validation failed.")
    validate_local_execution_budget(
        report.config.execution.max_cpus, report.config.execution.max_memory_gb, detect_local_resource_capacity()
    )
    resources = effective_resource_budget(
        runtime_snapshot(report.config.runtime.control_plane_image), project_execution_budget(report.config)
    )
    validate_effective_resource_budget(resources)
    require_fresh_plan(report)
    if report.config.input.type is InputType.FASTQ:
        if not report.execution_ready:
            raise ExecutionPreflightError("Execution is blocked: " + "; ".join(report.execution_blockers))
        _validate_custom_reference_files(report)
    nextflow, docker = check_nextflow(), check_docker()
    if nextflow.state != "FOUND":
        raise ExecutionPreflightError("Nextflow is required: " + nextflow.detail)
    if docker.state != "FOUND":
        raise ExecutionPreflightError("Docker is required: " + docker.detail)
    requested_image = report.config.runtime.control_plane_image
    container = check_container_runtime(requested_image)
    if container.state != "FOUND":
        raise ExecutionPreflightError("Control-plane container is required: " + container.detail)
    if report.config.reference.acceptance == "production":
        observed = inspect_container_image(requested_image)
        if not observed.get("image_id"):
            raise ExecutionPreflightError(
                "Production-intended execution requires an observed immutable control-plane image ID/digest; "
                f"Docker could not establish one for {requested_image}."
            )
    return resources


def reuse_upstream_if_compatible(run: CaseRun, frozen: FrozenInputs, reference: str) -> str:
    """Copy a prior nf-core output only when its frozen upstream contract matches."""

    pieces = reference.split("/")
    if len(pieces) != 2:
        raise ExecutionPreflightError("--reuse-upstream must be formatted as CASE-ID/RUN-ID.")
    source_case, source_run = pieces
    validate_case_id(source_case)
    if not RUN_ID_PATTERN.fullmatch(source_run):
        raise ExecutionPreflightError("--reuse-upstream has an unsafe run ID.")
    runs_dir = run.run_dir.parent.parent
    source = runs_dir / source_case / source_run
    try:
        source.resolve().relative_to(runs_dir.resolve())
    except ValueError as exc:
        raise ExecutionPreflightError("--reuse-upstream must refer to a run in this project.") from exc
    state_path = source / "run_state.json"
    if not state_path.is_file() or json.loads(state_path.read_text(encoding="utf-8")).get("status") != "SUCCESS":
        raise ExecutionPreflightError("--reuse-upstream must refer to a successful case run.")
    for name in ("input_manifest.yaml", "nfcore.params.json", "nfcore.local.config"):
        old, current = source / "frozen" / name, run.run_dir / "frozen" / name
        if not old.is_file() or not current.is_file() or old.read_bytes() != current.read_bytes():
            raise ExecutionPreflightError(f"--reuse-upstream is incompatible: frozen {name} differs.")
    if not (source / "handoff" / "upstream_manifest.yaml").is_file():
        raise ExecutionPreflightError("--reuse-upstream has no validated upstream handoff manifest.")
    shutil.copytree(source / "upstream", run.run_dir / "upstream", dirs_exist_ok=True, ignore=shutil.ignore_patterns("._*"))
    return f"{source_case}/{source_run}"


def execute_service_run(
    report: ValidationReport, *, case_id: str, profile: str = LOCAL_PROFILE, reuse_upstream: str | None = None,
) -> CaseRun:
    """Run the complete service lifecycle; final success requires downstream completion."""

    resources = prepare_service_run(report, profile=profile)
    run = create_case_run(report, case_id)
    command = ["rnaseq", "run", str(report.project_dir), "--case-id", case_id, "--profile", profile]
    if reuse_upstream:
        command.extend(["--reuse-upstream", reuse_upstream])
    try:
        workspace = resolve_execution_workspace(run.case_id, run.run_id)
        frozen = freeze_case_inputs(report, run, profile=profile, command=command, resources=resources)
        _write_yaml(
            run.run_dir / "provenance" / "run_provenance.yaml",
            _provenance(report, run, profile=profile, command=command, workspace=workspace, resources=resources),
        )
        prepare_execution_workspace(workspace)
        _write_state(run, "RUNNING", phase="freeze", command=command)
        assert report.config is not None
        if report.config.input.type is InputType.FASTQ:
            assert frozen.samplesheet and frozen.upstream_params and frozen.upstream_config
            reused_from = reuse_upstream_if_compatible(run, frozen, reuse_upstream) if reuse_upstream else None
            if reused_from is None:
                method = report.config.upstream.quantification.method if report.config.upstream.quantification else "salmon"
                upstream = (build_nextflow_command(
                    report, samplesheet=frozen.samplesheet, output_dir=run.run_dir / "upstream" / "nfcore_rnaseq",
                    profile=profile, params_file=frozen.upstream_params, config_file=frozen.upstream_config,
                    reference_paths=frozen.reference_paths, work_dir=workspace.work_dir / "upstream",
                ) if method == "salmon" else build_hisat2_featurecounts_command(
                    report, samplesheet=frozen.samplesheet, output_dir=run.run_dir / "upstream" / "hisat2_featurecounts",
                    profile=profile, reference_paths=frozen.reference_paths, work_dir=workspace.work_dir / "upstream",
                    config_file=frozen.upstream_config,
                ))
                _write_state(run, "RUNNING", phase="upstream", upstream_command=upstream)
                result = _run_command(upstream, cwd=workspace.launch_dir, stdout_path=run.run_dir / "logs" / "upstream.stdout.log", stderr_path=run.run_dir / "logs" / "upstream.stderr.log")
                if result != 0:
                    raise UpstreamExecutionError(
                        classify_execution_failure(
                            "nf-core/rnaseq" if method == "salmon" else "HISAT2 + featureCounts", result, run.run_dir / "logs" / "upstream.stderr.log",
                            resource=RESOURCE_CONTRACTS["MEDIUM"],
                        )
                    )
            finalize_fastq_handoff(report, run, frozen.contract, reused_from=reused_from)
        if report.config.project.preset is Preset.QC:
            # QC projects intentionally stop at the immutable FASTQ backend.
            # MultiQC and frozen provenance are still assembled into the normal
            # client delivery package, but no metadata-dependent downstream
            # statistical workflow is launched.
            delivery = assemble_delivery(run)
            _write_state(run, "SUCCESS", phase="delivery", delivery=str(delivery), downstream_skipped="technical_qc_only")
            return run
        execution_inputs = resolve_downstream_inputs(run)
        observer_config = write_downstream_observer_config(run)
        runtime_config = write_downstream_runtime_config(run, report.config.runtime.control_plane_image)
        docker_user_config = write_downstream_docker_user_config(run)
        downstream = build_downstream_nextflow_command(
            run, profile=profile, work_dir=workspace.work_dir / "downstream", observer_config=observer_config,
            docker_user_config=docker_user_config, execution_inputs=execution_inputs, runtime_config=runtime_config,
            local_resource_config=frozen.upstream_config,
        )
        _write_state(run, "RUNNING", phase="downstream", downstream_command=downstream)
        result = _run_command(downstream, cwd=workspace.launch_dir, stdout_path=run.run_dir / "logs" / "downstream.stdout.log", stderr_path=run.run_dir / "logs" / "downstream.stderr.log")
        if result != 0:
            raise UpstreamExecutionError(
                classify_execution_failure(
                    "downstream Nextflow", result, run.run_dir / "logs" / "downstream.stderr.log",
                    resource=RESOURCE_CONTRACTS["LARGE"],
                )
            )
        delivery = assemble_delivery(run)
        _write_state(run, "SUCCESS", phase="delivery", delivery=str(delivery))
    except (OSError, UpstreamExecutionError) as exc:
        _write_state(run, "FAILED", error=str(exc))
        raise
    return run
