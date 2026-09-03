"""Case/run lifecycle and client-facing delivery for the production workflow.

This module deliberately contains no statistics.  It freezes validated inputs,
starts independent Nextflow workflows, and curates their artifacts for a client.
"""

from __future__ import annotations

import csv
import hashlib
import json
import os
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
    CONTAINER_PROFILE,
    LOCAL_PROFILE,
    PreparedRun,
    ExecutionWorkspace,
    _validate_custom_reference_files,
    build_nextflow_command,
    classify_execution_failure,
    check_container_runtime,
    check_docker,
    check_nextflow,
    generate_handoff_manifest,
    nfcore_runtime_params,
    prepare_execution_workspace,
    RESOURCE_CONTRACTS,
    render_local_resource_config,
    require_fresh_plan,
    resolve_execution_workspace,
    runtime_snapshot,
)
from rnaseq.models import FastqPreprocessing, InputType, production_enrichment_backends
from rnaseq.planner import render_manifest
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
    return tuple(sorted(samples))


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
    if source_type not in {"raw_counts", "salmon_tximport"}:
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
            tx2gene = _safe_existing_under(
                run.run_dir, salmon.get("tx2gene"), "salmon.tx2gene", allowed_root=upstream_root
            )
            _copy_snapshot(tx2gene, temporary / "source" / "salmon.merged.tx2gene.tsv")
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
                "tx2gene": "source/salmon.merged.tx2gene.tsv",
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


def freeze_case_inputs(report: ValidationReport, run: CaseRun, *, profile: str, command: list[str]) -> FrozenInputs:
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
        reference_paths = {
            option.removeprefix("--"): path
            for option, path in report.local_reference.nfcore_arguments()
        }
    elif report.config.reference.source == "custom":
        for key in ("fasta", "gtf", "transcript_fasta", "salmon_index"):
            configured = getattr(report.config.reference, key)
            if configured is None:
                continue
            source_reference = (report.project_dir / configured).resolve()
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
        source = {
            "type": "salmon_tximport",
            "construction_method": "DESeqDataSetFromTximport",
            "upstream_handoff": None,
            "fastq_preprocessing": report.config.input.preprocessing.value,
            "skip_trimming": report.config.input.preprocessing is FastqPreprocessing.PRETRIMMED,
        }

    execution = {
        "case_id": run.case_id,
        "run_id": run.run_id,
        "timezone": "Asia/Taipei",
        "profile": profile,
        "command": command,
        "pipeline": {"name": report.config.project.pipeline, "version": "0.4.3"},
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
        # This is the immutable, normalized enrichment input for downstream
        # tasks.  Do not make those tasks re-read mutable project.yaml.
        "annotation": report.config.annotation.model_dump(mode="json") if report.config.annotation else None,
        "output_dir": str((run.run_dir / "downstream").resolve()),
    }
    contract_path = frozen / "downstream_contract.json"
    _write_text(contract_path, json.dumps(contract, indent=2, sort_keys=True) + "\n")
    params: Path | None = None
    runtime: Path | None = None
    if report.config.input.type is InputType.FASTQ:
        params = frozen / "nfcore.params.json"
        _write_text(params, json.dumps(nfcore_runtime_params(report), sort_keys=True) + "\n")
        runtime = frozen / "nfcore.local.config"
        _write_text(runtime, render_local_resource_config())
    return FrozenInputs(samplesheet, contract_path, manifest_path, params, runtime, reference_paths)


def _update_contract(contract_path: Path, mutate: dict[str, Any]) -> None:
    contract = json.loads(contract_path.read_text(encoding="utf-8"))
    contract.update(mutate)
    _write_text(contract_path, json.dumps(contract, indent=2, sort_keys=True) + "\n")


def _provenance(
    report: ValidationReport, run: CaseRun, *, profile: str, command: list[str], workspace: ExecutionWorkspace,
) -> dict[str, Any]:
    git_commit: str | None = None
    try:
        result = subprocess.run(["git", "rev-parse", "HEAD"], cwd=report.project_dir, capture_output=True, text=True, check=False)
        if result.returncode == 0:
            git_commit = result.stdout.strip() or None
    except OSError:
        pass
    nextflow = check_nextflow()
    runtime = runtime_snapshot()
    return {
        "case_id": run.case_id,
        "run_id": run.run_id,
        "timezone": "Asia/Taipei",
        "started_at": run.started_at,
        "pipeline_version": "0.4.3",
        "git_commit": git_commit,
        "python_version": sys.version.split()[0],
        "nextflow_version": nextflow.detail if nextflow.state == "FOUND" else None,
        "nfcore_rnaseq_version": report.config.upstream.pipeline_version if report.config and report.config.input.type is InputType.FASTQ else None,
        "container_runtime": "docker",
        "container_image": "rnaseq-control-plane:latest",
        "container_digest": None,
        "runtime_resources": {
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
        "profile": profile,
        "command": command,
        "execution_root": str(workspace.root),
        "execution_launch_dir": str(workspace.launch_dir),
        "execution_work_dir": str(workspace.work_dir),
        "input_manifest_sha256": _sha256(run.run_dir / "frozen" / "input_manifest.yaml"),
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


def finalize_fastq_handoff(report: ValidationReport, run: CaseRun, contract: Path, *, reused_from: str | None = None) -> Path:
    """Validate and freeze the stable nf-core handoff boundary for downstream."""

    prepared = PreparedRun(report, "unavailable", "docker")
    handoff = generate_handoff_manifest(prepared, run.run_dir)
    frozen_handoff = run.run_dir / "frozen" / "upstream_handoff_manifest.yaml"
    _copy_snapshot(handoff, frozen_handoff)
    _update_contract(
        contract,
        {"source": {"type": "salmon_tximport", "construction_method": "DESeqDataSetFromTximport", "upstream_handoff": str(frozen_handoff.resolve()), "reused_from": reused_from}},
    )
    return frozen_handoff


def build_downstream_nextflow_command(
    run: CaseRun, *, profile: str = LOCAL_PROFILE, work_dir: Path | None = None,
    observer_config: Path | None = None, execution_inputs: ResolvedDownstreamInputs | None = None,
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


def assemble_delivery(run: CaseRun) -> Path:
    """Build a client-safe package from an explicit allowlist of final artifacts."""

    delivery = _resolved_delivery_root(run.run_dir / "delivery", run_dir=run.run_dir)
    figures = {".png": delivery / "figures" / "png", ".tif": delivery / "figures" / "tiff_300dpi", ".tiff": delivery / "figures" / "tiff_300dpi"}
    tables = delivery / "tables"
    for directory in [*figures.values(), tables, delivery / "methods_and_versions"]:
        directory.mkdir(parents=True, exist_ok=True)
    downstream = run.run_dir / "downstream"
    report = downstream / "report" / "report.html"
    if report.is_file():
        _copy_delivery_artifact(report, delivery / delivery_filename("report.html", run.run_id))
    for path in downstream.rglob("*") if downstream.is_dir() else ():
        if not path.is_file() or _is_appledouble(path):
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


def prepare_service_run(report: ValidationReport, *, profile: str) -> None:
    if profile != LOCAL_PROFILE:
        raise ExecutionPreflightError("Only '--profile local' is implemented; server profiles are configuration placeholders.")
    if not report.is_valid or report.config is None:
        raise ExecutionPreflightError("Execution is blocked because validation failed.")
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
    container = check_container_runtime()
    if container.state != "FOUND":
        raise ExecutionPreflightError("Control-plane container is required: " + container.detail)


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

    prepare_service_run(report, profile=profile)
    run = create_case_run(report, case_id)
    command = ["rnaseq", "run", str(report.project_dir), "--case-id", case_id, "--profile", profile]
    if reuse_upstream:
        command.extend(["--reuse-upstream", reuse_upstream])
    try:
        workspace = resolve_execution_workspace(run.case_id, run.run_id)
        frozen = freeze_case_inputs(report, run, profile=profile, command=command)
        _write_yaml(
            run.run_dir / "provenance" / "run_provenance.yaml",
            _provenance(report, run, profile=profile, command=command, workspace=workspace),
        )
        prepare_execution_workspace(workspace)
        _write_state(run, "RUNNING", phase="freeze", command=command)
        assert report.config is not None
        if report.config.input.type is InputType.FASTQ:
            assert frozen.samplesheet and frozen.upstream_params and frozen.upstream_config
            reused_from = reuse_upstream_if_compatible(run, frozen, reuse_upstream) if reuse_upstream else None
            if reused_from is None:
                upstream = build_nextflow_command(
                    report, samplesheet=frozen.samplesheet, output_dir=run.run_dir / "upstream" / "nfcore_rnaseq",
                    profile=profile, params_file=frozen.upstream_params, config_file=frozen.upstream_config,
                    reference_paths=frozen.reference_paths, work_dir=workspace.work_dir / "upstream",
                )
                _write_state(run, "RUNNING", phase="upstream", upstream_command=upstream)
                result = _run_command(upstream, cwd=workspace.launch_dir, stdout_path=run.run_dir / "logs" / "upstream.stdout.log", stderr_path=run.run_dir / "logs" / "upstream.stderr.log")
                if result != 0:
                    raise UpstreamExecutionError(
                        classify_execution_failure(
                            "nf-core/rnaseq", result, run.run_dir / "logs" / "upstream.stderr.log",
                            resource=RESOURCE_CONTRACTS["MEDIUM"],
                        )
                    )
            finalize_fastq_handoff(report, run, frozen.contract, reused_from=reused_from)
        execution_inputs = resolve_downstream_inputs(run)
        observer_config = write_downstream_observer_config(run)
        downstream = build_downstream_nextflow_command(
            run, profile=profile, work_dir=workspace.work_dir / "downstream", observer_config=observer_config,
            execution_inputs=execution_inputs,
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
