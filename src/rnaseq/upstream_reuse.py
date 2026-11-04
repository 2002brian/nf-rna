"""Copy a reused upstream output into a new run without its large alignment artifacts.

``--reuse-upstream`` (and a retry of a run that reused its upstream) gives the
new run its own copy of the source run's ``upstream/`` tree, so the new run
stays self-contained for everything downstream and delivery actually read: the
count matrix, the sample map, featureCounts summaries, strandedness evidence,
MultiQC and QC reports.

The HISAT2 + featureCounts route also publishes coordinate-sorted BAMs (and
runs before nf-rna 1.7.0's fix also published uncompressed lane SAMs and the
prepared FASTQs, for pretrimmed input byte copies of the source FASTQ).  No
downstream, delivery or status code reads them; they are alignment provenance
of the run that produced them.  Copying them multiplied storage and I/O (about
240 GB per reuse for a real 8-sample project) without any analytical use, so
they are *referenced* in place instead: :data:`REFERENCE_RECORD` lists each one
with its location in the producing run and its size.

Deliberately not used:

* hard links -- they fail across devices, and a shared inode means a later
  write through either path silently changes the other run's evidence;
* symbolic links -- they dangle when the source run is archived or deleted and
  would let a later path resolution escape the new run directory.

The copy is staged in a sibling directory and moved into place only after every
file has been copied and its size verified, so an interrupted or failed copy
never leaves a partial ``upstream/`` that looks complete.
"""

from __future__ import annotations

import json
import os
import shutil
import uuid
from pathlib import Path, PurePosixPath
from typing import Any

import yaml

#: Upstream route whose alignment artifacts are referenced, not copied.  Salmon
#: (``nfcore_rnaseq``) trees are copied exactly as before.
REFERENCED_ROUTE = "hisat2_featurecounts"
#: Alignment artifact suffixes that downstream never reads.
ALIGNMENT_SUFFIXES = (".sam", ".bam", ".bai", ".cram", ".crai", ".csi")
#: Prepared reads that runs before the FASTP_PREPARE publish fix published
#: (for pretrimmed input, byte copies of the source FASTQ); never read downstream.
PREPARED_READS_DIRECTORY = ("qc", "fastp", "prepared")
PREPARED_READS_SUFFIXES = (".fastq.gz", ".fq.gz", ".fastq", ".fq")
#: Record written at ``upstream/hisat2_featurecounts/<REFERENCE_RECORD>``.
REFERENCE_RECORD = "referenced_alignment_artifacts.json"
REFERENCE_RECORD_SCHEMA = "nf-rna.referenced-alignment-artifacts.v1"
_REFERENCE_POLICY = (
    "Alignment artifacts (SAM/BAM/CRAM and their indexes) and published prepared reads (qc/fastp/prepared) of a "
    "reused HISAT2 + featureCounts upstream are not copied: downstream analysis and delivery never read them. "
    "Each entry names the run that holds the file."
)


class UpstreamReuseCopyError(OSError):
    """The reused upstream could not be copied completely; nothing was put in place."""


def is_referenced_artifact(relative: PurePosixPath) -> bool:
    """Whether ``relative`` (a path under ``upstream/``) is a large artifact left in its source run."""

    return artifact_kind(relative) is not None


def artifact_kind(relative: PurePosixPath) -> str | None:
    """``alignment`` or ``prepared_reads`` for a referenced HISAT2-route artifact, otherwise ``None``."""

    if len(relative.parts) < 2 or relative.parts[0] != REFERENCED_ROUTE:
        return None
    if relative.name.endswith(ALIGNMENT_SUFFIXES):
        return "alignment"
    if relative.parts[1:-1] == PREPARED_READS_DIRECTORY and relative.name.endswith(PREPARED_READS_SUFFIXES):
        return "prepared_reads"
    return None


def _ignored(name: str) -> bool:
    # Same filter as the previous shutil.copytree(..., ignore=ignore_patterns("._*")).
    return name.startswith("._")


def _load_record(path: Path) -> dict[str, Any] | None:
    if not path.is_file():
        return None
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError) as exc:
        raise UpstreamReuseCopyError(f"Reused upstream has an unreadable {REFERENCE_RECORD}: {exc}") from exc
    if not isinstance(payload, dict) or payload.get("schema_version") != REFERENCE_RECORD_SCHEMA or not isinstance(payload.get("artifacts"), list):
        raise UpstreamReuseCopyError(f"Reused upstream has an invalid {REFERENCE_RECORD}.")
    return payload


def _source_pipeline(source_run_dir: Path) -> object:
    try:
        manifest = yaml.safe_load((source_run_dir / "handoff" / "upstream_manifest.yaml").read_text(encoding="utf-8"))
    except (OSError, yaml.YAMLError):
        return None
    return manifest.get("pipeline") if isinstance(manifest, dict) else None


def copy_reusable_upstream(source_run_dir: Path, destination: Path, *, source_label: str) -> dict[str, Any]:
    """Copy ``source_run_dir/upstream`` to ``destination`` (an absent or empty directory).

    Returns a summary with the copied and referenced file counts and bytes.
    Raises :class:`UpstreamReuseCopyError` (an ``OSError``) on any failure, after
    removing the staging directory; ``destination`` is then left as it was.
    """

    source = source_run_dir / "upstream"
    if not source.is_dir():
        raise UpstreamReuseCopyError(f"Reused run {source_label} has no upstream output directory.")
    if destination.exists() and (not destination.is_dir() or any(destination.iterdir())):
        raise UpstreamReuseCopyError(f"Refusing to reuse an upstream into a non-empty directory: {destination}")
    inherited = _load_record(source / REFERENCED_ROUTE / REFERENCE_RECORD)
    staging = destination.parent / f".{destination.name}.partial-{uuid.uuid4().hex}"
    copied_files = copied_bytes = newly_referenced = 0
    referenced: list[dict[str, Any]] = list(inherited["artifacts"]) if inherited else []
    try:
        staging.parent.mkdir(parents=True, exist_ok=True)
        staging.mkdir()
        for directory, subdirectories, files in os.walk(source):
            subdirectories[:] = sorted(name for name in subdirectories if not _ignored(name))
            current = Path(directory)
            relative_dir = current.relative_to(source)
            (staging / relative_dir).mkdir(parents=True, exist_ok=True)
            for name in sorted(files):
                if _ignored(name):
                    continue
                origin = current / name
                relative = PurePosixPath(*(relative_dir / name).parts)
                if inherited and relative == PurePosixPath(REFERENCED_ROUTE, REFERENCE_RECORD):
                    continue  # rewritten below with the inherited entries kept
                size = origin.stat().st_size
                if is_referenced_artifact(relative):
                    referenced.append({
                        "path": relative.as_posix(),
                        "kind": artifact_kind(relative),
                        "size_bytes": size,
                        "source_run": source_label,
                        "source_path": str(origin.resolve()),
                    })
                    newly_referenced += 1
                    continue
                target = staging / relative_dir / name
                shutil.copy2(origin, target)
                if target.stat().st_size != size:
                    raise UpstreamReuseCopyError(f"Incomplete copy of reused upstream file {relative} ({target.stat().st_size} of {size} bytes).")
                copied_files += 1
                copied_bytes += size
        if referenced:
            record = {
                "schema_version": REFERENCE_RECORD_SCHEMA,
                "policy": _REFERENCE_POLICY,
                "reused_from": inherited.get("reused_from") if inherited else source_label,
                "source_pipeline": inherited.get("source_pipeline") if inherited else _source_pipeline(source_run_dir),
                "artifacts": sorted(referenced, key=lambda item: str(item.get("path"))),
                "referenced_bytes": sum(int(item.get("size_bytes") or 0) for item in referenced),
            }
            (staging / REFERENCED_ROUTE).mkdir(parents=True, exist_ok=True)
            (staging / REFERENCED_ROUTE / REFERENCE_RECORD).write_text(json.dumps(record, indent=2) + "\n", encoding="utf-8")
        if destination.exists():
            destination.rmdir()  # verified empty above
        os.replace(staging, destination)
    except BaseException as exc:
        shutil.rmtree(staging, ignore_errors=True)
        if isinstance(exc, OSError) and not isinstance(exc, UpstreamReuseCopyError):
            raise UpstreamReuseCopyError(f"Could not copy the reused upstream of {source_label}: {exc}") from exc
        raise
    return {
        "copied_files": copied_files,
        "copied_bytes": copied_bytes,
        "referenced_files": len(referenced),
        "referenced_bytes": sum(int(item.get("size_bytes") or 0) for item in referenced),
        "newly_referenced_files": newly_referenced,
    }


def referenced_alignment_record(upstream_output: Path) -> dict[str, Any] | None:
    """The reference record of a reused HISAT2 upstream, or ``None`` when its artifacts are local.

    An unreadable or invalid record raises :class:`UpstreamReuseCopyError`.
    """

    return _load_record(upstream_output / REFERENCE_RECORD)
