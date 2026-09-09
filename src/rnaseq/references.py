"""Versioned local-reference validation and one-time Salmon index preparation."""

from __future__ import annotations

import hashlib
import json
import os
import platform
import re
import shutil
import subprocess
import uuid
from copy import deepcopy
from dataclasses import dataclass
from datetime import UTC, datetime
from pathlib import Path
from typing import Any, Callable

import yaml

from rnaseq.models import ReferenceConfig


SHA256_PATTERN = re.compile(r"^[0-9a-fA-F]{64}$")
MD5_PATTERN = re.compile(r"^[0-9a-fA-F]{32}$")
SEMVER_PATTERN = re.compile(r"^[0-9]+\.[0-9]+(?:\.[0-9]+)?$")
LEGACY_LOCAL_REFERENCE_MANIFEST_VERSION = "1.0"
LOCAL_REFERENCE_MANIFEST_VERSION = "1.2"
LOCAL_REFERENCE_MANIFEST_VERSIONS = frozenset(
    (LEGACY_LOCAL_REFERENCE_MANIFEST_VERSION, "1.1", LOCAL_REFERENCE_MANIFEST_VERSION)
)
REFERENCE_REGISTRY_SCHEMA_VERSION = "1.0"
REFERENCE_PURPOSES = frozenset(("synthetic_test", "production"))
SALMON_NOT_BUILT = "not_built"
SALMON_BUILT = "built"
SALMON_VERSION = "1.10.3"
SALMON_KMER_SIZE = 31
# Kept only to read the provenance of historic nf-rna-built indexes.  New
# manifests bind an index directly to ``files.transcript_fasta`` and do not
# require RSEM to have created that file.
TRANSCRIPTOME_STRATEGY = "nfcore_rnaseq_3.26.0_rsem_from_genome_fasta_and_gtf"
DECOY_STRATEGY = "nfcore_rnaseq_3.26.0_gentrome"
SALMON_STRATEGY_TRANSCRIPTOME_ONLY = "transcriptome_only"
SALMON_STRATEGY_DECOY_AWARE = "decoy_aware"
SALMON_STRATEGIES = frozenset((SALMON_STRATEGY_TRANSCRIPTOME_ONLY, SALMON_STRATEGY_DECOY_AWARE))
ADOPTED_EXISTING_INDEX = "adopted_existing_index"
HISAT2_VERSION = "2.2.3"
HISAT2_RUNTIME_VERSION = HISAT2_VERSION
HISAT2_NOT_BUILT = "not_built"
HISAT2_BUILT = "built"
HISAT2_STRATEGY_GRAPH_EMBEDDED = "graph_embedded_splice_sites"
HISAT2_STRATEGY_LEGACY_ANNOTATION_AWARE = "annotation_aware_genome_index_with_splice_sites"
HISAT2_STRATEGY_GENOME_ONLY_RUNTIME_SPLICES = "genome_only_runtime_splicesites"
# 1.2 manifests written during the earlier compatibility work used this
# spelling.  It remains readable, but new manifests use the public spelling
# above.
HISAT2_STRATEGY_GENOME_ONLY_RUNTIME_SPLICES_LEGACY = "genome_only_runtime_splices"
HISAT2_STRATEGIES = frozenset((
    HISAT2_STRATEGY_GRAPH_EMBEDDED,
    HISAT2_STRATEGY_LEGACY_ANNOTATION_AWARE,
    HISAT2_STRATEGY_GENOME_ONLY_RUNTIME_SPLICES,
    HISAT2_STRATEGY_GENOME_ONLY_RUNTIME_SPLICES_LEGACY,
))
HISAT2_RUNTIME_COMPATIBILITY_VALIDATED = "validated"
HISAT2_RUNTIME_COMPATIBILITY_SMOKE_REQUIRED = "requires_smoke_validation"
REQUIRED_HISAT2_INDEX_SUFFIXES = tuple(f".{index}.ht2" for index in range(1, 9))
REQUIRED_SALMON_INDEX_FILES = (
    "complete_ref_lens.bin", "ctable.bin", "ctg_offsets.bin", "duplicate_clusters.tsv",
    "info.json", "mphf.bin", "pos.bin", "rank.bin", "refAccumLengths.bin", "reflengths.bin",
    "refseq.bin", "seq.bin", "versionInfo.json",
)


class LocalReferenceError(ValueError):
    """Raised when a configured local reference cannot be used safely."""


class ReferencePreparationError(RuntimeError):
    """Raised when one-time local-reference preparation cannot complete."""


class ReferenceAdoptionError(RuntimeError):
    """Raised when an existing Salmon index cannot be adopted safely."""


@dataclass(frozen=True)
class LocalReferenceAsset:
    name: str
    path: Path
    relative_path: str
    sha256: str


@dataclass(frozen=True)
class LocalReference:
    """Resolved immutable identity of one local FASTQ reference release."""

    root: Path
    manifest_path: Path
    manifest_sha256: str
    manifest_schema_version: str
    purpose: str | None
    species: str
    provider: str
    release: int
    assembly: str
    assembly_patch: str | None
    genome_fasta: LocalReferenceAsset
    annotation_gtf: LocalReferenceAsset
    transcript_fasta: LocalReferenceAsset | None
    salmon_index: Path | None
    salmon_status: str
    salmon_strategy_type: str | None
    salmon_provenance: dict[str, object] | None
    salmon_transcriptome: LocalReferenceAsset | None = None
    salmon_index_metadata: dict[str, object] | None = None
    salmon_validation: dict[str, object] | None = None
    hisat2_index: Path | None = None
    hisat2_index_prefix: Path | None = None
    hisat2_splice_sites: LocalReferenceAsset | None = None
    hisat2_status: str = HISAT2_NOT_BUILT
    hisat2_strategy: str | None = None
    hisat2_provenance: dict[str, object] | None = None
    source_provenance: dict[str, object] | None = None

    @property
    def salmon_strategy(self) -> str:
        if self.salmon_index is None:
            return "no validated pre-built Salmon index is registered in the manifest"
        assert self.salmon_strategy_type is not None
        return f"use the manifest-declared pre-built {self.salmon_strategy_type.replace('_', '-')} Salmon index"

    @property
    def transcriptome_strategy(self) -> str:
        if self.transcript_fasta is None:
            return "not_applicable_without_salmon"
        if self.transcript_fasta is not None and self.salmon_transcriptome == self.transcript_fasta:
            return "manifest_registered_transcript_fasta"
        if self.salmon_transcriptome is not None:
            return "gtf_derived_exact_transcript_id_contract"
        return TRANSCRIPTOME_STRATEGY

    @property
    def external_transcript_fasta_used(self) -> bool:
        """Whether the selected Salmon index is checksum-bound to this FASTA."""

        return self.transcript_fasta is not None and self.salmon_transcriptome == self.transcript_fasta

    @property
    def hisat2_runtime_ready(self) -> bool:
        """Whether the declared index/runtime pair has acceptance evidence."""

        if self.hisat2_index is None:
            return False
        if not _uses_runtime_splice_sites(self.hisat2_strategy):
            return True
        return bool(self.hisat2_provenance) and (
            self.hisat2_provenance.get("runtime_compatibility") == HISAT2_RUNTIME_COMPATIBILITY_VALIDATED
        )

    @property
    def assembly_identity(self) -> str:
        """Render a biological assembly without inventing an absent patch."""

        return f"{self.assembly}.{self.assembly_patch}" if self.assembly_patch else self.assembly

    def assets(self) -> tuple[LocalReferenceAsset, ...]:
        return tuple(asset for asset in (self.genome_fasta, self.annotation_gtf, self.transcript_fasta) if asset is not None)

    def hisat2_arguments(self) -> list[tuple[str, Path]]:
        if self.hisat2_index is None:
            raise LocalReferenceError(
                "Local reference manifest has no validated HISAT2 index. "
                "Register a prebuilt index and matching splice_sites asset, or run the optional host-native builder."
            )
        arguments = [("--fasta", self.genome_fasta.path), ("--gtf", self.annotation_gtf.path), ("--hisat2_index", self.hisat2_index)]
        if _uses_runtime_splice_sites(self.hisat2_strategy):
            assert self.hisat2_splice_sites is not None
            arguments.append(("--hisat2_splice_sites", self.hisat2_splice_sites.path))
        return arguments

    def nfcore_arguments(self) -> list[tuple[str, Path]]:
        if self.salmon_index is None:
            raise LocalReferenceError(
                "Local reference manifest has no validated Salmon index. "
                "Register a prebuilt index, or run the optional host-native builder."
            )
        # Never pass downloaded Ensembl cDNA as --transcript_fasta: nf-core
        # derives tx2gene-compatible transcripts from the frozen GTF itself.
        return [
            ("--fasta", self.genome_fasta.path),
            ("--gtf", self.annotation_gtf.path),
            ("--salmon_index", self.salmon_index),
        ]

    def provenance(self) -> dict[str, object]:
        return {
            "source": "local",
            "root": str(self.root),
            "purpose": self.purpose,
            "species": self.species,
            "provider": self.provider,
            "release": self.release,
            "assembly": self.assembly,
            "assembly_patch": self.assembly_patch,
            "manifest": {
                "path": str(self.manifest_path),
                "sha256": self.manifest_sha256,
                "schema_version": self.manifest_schema_version,
            },
            "source_provenance": self.source_provenance,
            "assets": {
                asset.name: {"path": str(asset.path), "sha256": asset.sha256}
                for asset in self.assets()
            },
            "transcriptome": {
                "strategy": self.transcriptome_strategy,
                "genome_fasta": {"path": str(self.genome_fasta.path), "sha256": self.genome_fasta.sha256},
                "annotation_gtf": {"path": str(self.annotation_gtf.path), "sha256": self.annotation_gtf.sha256},
                "external_transcript_fasta": {
                    "path": str(self.transcript_fasta.path), "sha256": self.transcript_fasta.sha256,
                } if self.transcript_fasta is not None else None,
                "external_transcript_fasta_used": self.external_transcript_fasta_used,
                "adopted_transcriptome": (
                    {"path": str(self.salmon_transcriptome.path), "sha256": self.salmon_transcriptome.sha256}
                    if self.salmon_transcriptome is not None else None
                ),
            },
            "salmon": {
                "status": self.salmon_status,
                "index": str(self.salmon_index) if self.salmon_index is not None else None,
                "strategy": self.salmon_strategy_type,
                "provenance": self.salmon_provenance,
                "index_metadata": self.salmon_index_metadata,
                "validation": self.salmon_validation,
            },
            "hisat2": {
                "status": self.hisat2_status,
                "index": str(self.hisat2_index) if self.hisat2_index is not None else None,
                "index_prefix": str(self.hisat2_index_prefix) if self.hisat2_index_prefix is not None else None,
                "strategy": self.hisat2_strategy,
                "splice_sites": ({"path": str(self.hisat2_splice_sites.path), "sha256": self.hisat2_splice_sites.sha256} if self.hisat2_splice_sites else None),
                "index_builder_version": (
                    self.hisat2_provenance.get("index_builder_version", self.hisat2_provenance.get("hisat2_version"))
                    if self.hisat2_provenance else None
                ),
                "runtime_aligner_version": HISAT2_RUNTIME_VERSION,
                "runtime_arguments": (
                    ["hisat2", "-x", str(self.hisat2_index_prefix or (self.hisat2_index / "genome")), "--known-splicesite-infile", str(self.hisat2_splice_sites.path)]
                    if _uses_runtime_splice_sites(self.hisat2_strategy) and self.hisat2_index and self.hisat2_splice_sites
                    else ["hisat2", "-x", str(self.hisat2_index_prefix or (self.hisat2_index / "genome"))] if self.hisat2_index else None
                ),
                "provenance": self.hisat2_provenance,
            },
        }


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _digest_file(path: Path, algorithm: str) -> str:
    digest = hashlib.new(algorithm)
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _require_mapping(value: object, label: str) -> dict[str, Any]:
    if not isinstance(value, dict):
        raise LocalReferenceError(f"Local reference manifest {label} must be a mapping.")
    return value


def _require_string(value: object, label: str) -> str:
    if not isinstance(value, str) or not value.strip():
        raise LocalReferenceError(f"Local reference manifest {label} must be a non-blank string.")
    return value


def _resolve_under(root: Path, configured: str, label: str, *, directory: bool = False) -> Path:
    candidate = Path(configured)
    if candidate.is_absolute():
        raise LocalReferenceError(f"Local reference manifest {label} must be relative to reference.root.")
    resolved = (root / candidate).resolve()
    try:
        resolved.relative_to(root)
    except ValueError as exc:
        raise LocalReferenceError(f"Local reference manifest {label} escapes reference.root.") from exc
    if directory:
        if not resolved.is_dir():
            raise LocalReferenceError(f"Local reference {label} directory not found: {resolved}")
    elif not resolved.is_file():
        raise LocalReferenceError(f"Local reference {label} file not found: {resolved}")
    return resolved


def _asset(root: Path, files: dict[str, Any], key: str) -> LocalReferenceAsset:
    return _asset_from_payload(root, _require_mapping(files.get(key), f"files.{key}"), f"files.{key}", key)


def _asset_from_payload(root: Path, payload: dict[str, Any], label: str, name: str) -> LocalReferenceAsset:
    relative_path = _require_string(payload.get("path"), f"{label}.path")
    expected_sha256 = _require_string(payload.get("sha256"), f"{label}.sha256")
    if not SHA256_PATTERN.fullmatch(expected_sha256):
        raise LocalReferenceError(f"Local reference manifest {label}.sha256 must be a SHA256 hex digest.")
    path = _resolve_under(root, relative_path, label)
    observed_sha256 = sha256_file(path)
    if observed_sha256.lower() != expected_sha256.lower():
        raise LocalReferenceError(
            f"Local reference checksum mismatch for {name}: expected {expected_sha256}, observed {observed_sha256}."
        )
    return LocalReferenceAsset(name, path, relative_path, expected_sha256.lower())


def _validated_source_provenance(
    value: object,
    assets: dict[str, LocalReferenceAsset],
) -> dict[str, object] | None:
    """Validate only supplied upstream identities and prove supplied checksums."""

    if value is None:
        return None
    sources = _require_mapping(value, "sources")
    unknown = sorted(set(sources) - set(assets))
    if unknown:
        raise LocalReferenceError("Local reference manifest sources has unknown asset(s): " + ", ".join(unknown) + ".")
    result: dict[str, object] = {}
    for name, raw in sources.items():
        item = _require_mapping(raw, f"sources.{name}")
        recorded: dict[str, object] = {}
        for field in ("url", "accession"):
            if field in item:
                recorded[field] = _require_string(item[field], f"sources.{name}.{field}")
        if "upstream_checksum" in item:
            checksum = _require_mapping(item["upstream_checksum"], f"sources.{name}.upstream_checksum")
            algorithm = _require_string(
                checksum.get("algorithm"), f"sources.{name}.upstream_checksum.algorithm"
            ).lower()
            expected = _require_string(
                checksum.get("value"), f"sources.{name}.upstream_checksum.value"
            ).lower()
            pattern = SHA256_PATTERN if algorithm == "sha256" else MD5_PATTERN if algorithm == "md5" else None
            if pattern is None or not pattern.fullmatch(expected):
                raise LocalReferenceError(
                    f"Local reference manifest sources.{name}.upstream_checksum must be a valid md5 or sha256 digest."
                )
            observed = _digest_file(assets[name].path, algorithm)
            if observed != expected:
                raise LocalReferenceError(
                    f"Local reference upstream checksum mismatch for {name}: expected {expected}, observed {observed}."
                )
            recorded["upstream_checksum"] = {
                "algorithm": algorithm,
                "value": expected,
                "verification": "matched_local_asset",
            }
        unknown_fields = sorted(set(item) - {"url", "accession", "upstream_checksum"})
        if unknown_fields:
            raise LocalReferenceError(
                f"Local reference manifest sources.{name} has unsupported field(s): " + ", ".join(unknown_fields) + "."
            )
        if not recorded:
            raise LocalReferenceError(
                f"Local reference manifest sources.{name} must supply a URL, accession, or upstream checksum."
            )
        result[name] = recorded
    return result


def _validated_salmon_provenance(
    value: object, genome_fasta: LocalReferenceAsset, annotation_gtf: LocalReferenceAsset
) -> dict[str, object]:
    provenance = _require_mapping(value, "salmon.provenance")
    expected_strings = {
        "salmon_version": SALMON_VERSION,
        "genome_fasta_sha256": genome_fasta.sha256,
        "annotation_gtf_sha256": annotation_gtf.sha256,
        "transcriptome_strategy": TRANSCRIPTOME_STRATEGY,
        "decoy_strategy": DECOY_STRATEGY,
    }
    for key, expected in expected_strings.items():
        observed = _require_string(provenance.get(key), f"salmon.provenance.{key}")
        if observed != expected:
            raise LocalReferenceError(
                f"Local reference salmon.provenance.{key} is incompatible: expected {expected!r}, observed {observed!r}."
            )
    if provenance.get("kmer_size") != SALMON_KMER_SIZE or isinstance(provenance.get("kmer_size"), bool):
        raise LocalReferenceError(f"Local reference salmon.provenance.kmer_size must be {SALMON_KMER_SIZE}.")
    _require_string(provenance.get("built_at"), "salmon.provenance.built_at")
    return dict(provenance)


def validate_salmon_index(index: Path) -> None:
    """Reject incomplete or clearly corrupt Salmon index directories before execution."""

    if not index.is_dir():
        raise LocalReferenceError(f"Local reference salmon.index directory not found: {index}")
    missing = [name for name in REQUIRED_SALMON_INDEX_FILES if not (index / name).is_file()]
    empty = [
        name for name in REQUIRED_SALMON_INDEX_FILES
        if (index / name).is_file() and (index / name).stat().st_size == 0
    ]
    if missing or empty:
        details = []
        if missing:
            details.append("missing " + ", ".join(missing))
        if empty:
            details.append("empty " + ", ".join(empty))
        raise LocalReferenceError("Local reference Salmon index is incomplete: " + "; ".join(details) + ".")


def validate_hisat2_index(index: Path, *, prefix: str = "genome") -> None:
    """Reject incomplete HISAT2 indexes before they can become executable."""

    if not index.is_dir():
        raise LocalReferenceError(f"Local reference hisat2.index directory not found: {index}")
    # HISAT2 supports both small (.ht2) and large (.ht2l) indexes.  A valid
    # index has one complete numbered family, never a partial mixture.
    if not prefix or Path(prefix).name != prefix:
        raise LocalReferenceError("Local reference HISAT2 index prefix must be a simple non-blank basename.")
    small = [index / f"{prefix}{suffix}" for suffix in REQUIRED_HISAT2_INDEX_SUFFIXES]
    large = [index / f"{prefix}.{number}.ht2l" for number in range(1, 9)]
    family = small if all(path.is_file() for path in small) else large
    if not all(path.is_file() and path.stat().st_size > 0 for path in family):
        raise LocalReferenceError(f"Local reference HISAT2 index is incomplete; expected {prefix}.1..8.ht2 or .ht2l files.")


def _load_json_mapping(path: Path, label: str) -> dict[str, Any]:
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError) as exc:
        raise LocalReferenceError(f"Local reference {label} is not valid JSON: {exc}") from exc
    return _require_mapping(payload, label)


def _require_nonnegative_integer(value: object, label: str) -> int:
    if not isinstance(value, int) or isinstance(value, bool) or value < 0:
        raise LocalReferenceError(f"Local reference manifest {label} must be a non-negative integer.")
    return value


def _require_equal(value: object, expected: object, label: str) -> None:
    if value != expected:
        raise LocalReferenceError(
            f"Local reference {label} is incompatible: expected {expected!r}, observed {value!r}."
        )


def _uses_runtime_splice_sites(strategy: str | None) -> bool:
    return strategy in {
        HISAT2_STRATEGY_GENOME_ONLY_RUNTIME_SPLICES,
        HISAT2_STRATEGY_GENOME_ONLY_RUNTIME_SPLICES_LEGACY,
    }


def _salmon_index_metadata(index: Path) -> dict[str, object]:
    """Read the stable metadata emitted by the pinned Salmon 1.10.3 index format."""

    validate_salmon_index(index)
    info = _load_json_mapping(index / "info.json", "salmon.index/info.json")
    version = _load_json_mapping(index / "versionInfo.json", "salmon.index/versionInfo.json")
    required_info = {
        "index_version": "info_index_version",
        "k": "kmer_size",
        "num_decoys": "num_decoys",
        "SeqHash": "seq_hash",
        "NameHash": "name_hash",
        "seq_length": "seq_length",
        "num_kmers": "num_kmers",
        "num_contigs": "num_contigs",
    }
    result: dict[str, object] = {}
    for source, destination in required_info.items():
        if source not in info:
            raise LocalReferenceError(f"Local reference salmon.index/info.json is missing {source}.")
        result[destination] = info[source]
    if "salmonVersion" not in version or "indexVersion" not in version:
        raise LocalReferenceError("Local reference salmon.index/versionInfo.json is missing salmonVersion or indexVersion.")
    result["salmon_version"] = version["salmonVersion"]
    result["salmon_index_version"] = version["indexVersion"]
    return result


def _validated_prebuilt_salmon(
    root: Path,
    salmon: dict[str, Any],
    genome_fasta: LocalReferenceAsset,
    transcript_fasta: LocalReferenceAsset,
) -> tuple[Path, dict[str, object], dict[str, object]]:
    """Validate a first-class externally built Salmon index declaration.

    The index is validated from its own metadata and is tied to the source
    transcriptome by checksum.  No claim is made about which tool created it.
    """

    strategy = _require_string(salmon.get("strategy"), "salmon.strategy")
    if strategy not in SALMON_STRATEGIES:
        raise LocalReferenceError(
            "Local reference salmon.strategy must be one of: "
            f"{SALMON_STRATEGY_TRANSCRIPTOME_ONLY}, {SALMON_STRATEGY_DECOY_AWARE}."
        )
    index = _resolve_under(root, _require_string(salmon.get("index"), "salmon.index"), "salmon.index", directory=True)
    metadata = _salmon_index_metadata(index)
    version = _require_string(salmon.get("version", salmon.get("salmon_version")), "salmon.version")
    _require_equal(version, metadata["salmon_version"], "salmon.version versus index metadata")
    _require_equal(version, SALMON_VERSION, "salmon.version")
    _require_equal(metadata["kmer_size"], SALMON_KMER_SIZE, "salmon index k")

    source_transcriptome = _require_string(
        salmon.get("source_transcriptome_sha256"), "salmon.source_transcriptome_sha256"
    ).lower()
    if not SHA256_PATTERN.fullmatch(source_transcriptome):
        raise LocalReferenceError("Local reference salmon.source_transcriptome_sha256 must be a SHA256 hex digest.")
    _require_equal(
        source_transcriptome, transcript_fasta.sha256, "salmon.source_transcriptome_sha256 versus files.transcript_fasta.sha256"
    )
    expected_decoys = 0 if strategy == SALMON_STRATEGY_TRANSCRIPTOME_ONLY else None
    if expected_decoys is not None:
        _require_equal(metadata["num_decoys"], expected_decoys, "salmon index num_decoys")
    elif not isinstance(metadata["num_decoys"], int) or metadata["num_decoys"] <= 0:
        raise LocalReferenceError("Local reference decoy_aware Salmon index must declare at least one decoy in info.json.")
    if strategy == SALMON_STRATEGY_DECOY_AWARE:
        genome_sha256 = _require_string(salmon.get("source_genome_sha256"), "salmon.source_genome_sha256").lower()
        if not SHA256_PATTERN.fullmatch(genome_sha256):
            raise LocalReferenceError("Local reference salmon.source_genome_sha256 must be a SHA256 hex digest.")
        _require_equal(genome_sha256, genome_fasta.sha256, "salmon.source_genome_sha256 versus files.genome_fasta.sha256")
    provenance = dict(_require_mapping(salmon.get("provenance", {}), "salmon.provenance"))
    provenance.setdefault("mode", "prebuilt")
    provenance.setdefault("source_transcriptome_sha256", source_transcriptome)
    return index, metadata, provenance


def _is_prebuilt_salmon_declaration(salmon: dict[str, Any]) -> bool:
    return "source_transcriptome_sha256" in salmon


def _validation_artifact(
    root: Path,
    salmon: dict[str, Any],
    genome_fasta: LocalReferenceAsset,
    annotation_gtf: LocalReferenceAsset,
    transcriptome: LocalReferenceAsset,
    identity: dict[str, Any],
) -> tuple[Path, dict[str, object]]:
    validation = _require_mapping(salmon.get("validation"), "salmon.validation")
    artifact_path = _resolve_under(root, _require_string(validation.get("artifact"), "salmon.validation.artifact"), "salmon.validation.artifact")
    artifact = _load_json_mapping(artifact_path, "salmon.validation.artifact")
    _require_equal(artifact.get("status"), "PASS", "salmon validation artifact status")
    artifact_reference = _require_mapping(artifact.get("reference"), "salmon.validation.artifact.reference")
    for key, expected in {
        "species": identity["species"],
        "provider": identity["provider"],
        "release": identity["release"],
        "assembly": identity["assembly"],
    }.items():
        _require_equal(artifact_reference.get(key), expected, f"salmon validation artifact reference.{key}")
    inputs = _require_mapping(artifact.get("inputs"), "salmon.validation.artifact.inputs")
    for key, asset in (("genome", genome_fasta), ("gtf", annotation_gtf)):
        item = _require_mapping(inputs.get(key), f"salmon.validation.artifact.inputs.{key}")
        _require_equal(item.get("path"), asset.relative_path, f"salmon validation artifact inputs.{key}.path")
        _require_equal(item.get("sha256"), asset.sha256, f"salmon validation artifact inputs.{key}.sha256")
    generated = _require_mapping(artifact.get("generated_transcriptome"), "salmon.validation.artifact.generated_transcriptome")
    _require_equal(generated.get("path"), transcriptome.relative_path, "salmon validation artifact generated_transcriptome.path")
    _require_equal(generated.get("sha256"), transcriptome.sha256, "salmon validation artifact generated_transcriptome.sha256")
    transcript_count = _require_nonnegative_integer(generated.get("transcript_count"), "salmon.validation.artifact.generated_transcriptome.transcript_count")
    if transcript_count == 0:
        raise LocalReferenceError("Local reference salmon validation artifact has no generated transcripts.")
    compatibility = _require_mapping(artifact.get("identifier_compatibility"), "salmon.validation.artifact.identifier_compatibility")
    _require_equal(compatibility.get("fasta_only"), 0, "salmon validation artifact identifier_compatibility.fasta_only")
    _require_equal(compatibility.get("fasta_mapping_rate"), 1.0, "salmon validation artifact identifier_compatibility.fasta_mapping_rate")
    _require_equal(compatibility.get("gtf_mapping_rate"), 1.0, "salmon validation artifact identifier_compatibility.gtf_mapping_rate")
    _require_equal(compatibility.get("fasta_unique_ids"), transcript_count, "salmon validation artifact identifier_compatibility.fasta_unique_ids")
    tx2gene = _require_mapping(artifact.get("tx2gene"), "salmon.validation.artifact.tx2gene")
    _require_equal(tx2gene.get("zero_gene_mappings"), 0, "salmon validation artifact tx2gene.zero_gene_mappings")
    _require_equal(tx2gene.get("multi_gene_mappings"), 0, "salmon validation artifact tx2gene.multi_gene_mappings")
    _require_equal(tx2gene.get("one_gene_mappings"), transcript_count, "salmon validation artifact tx2gene.one_gene_mappings")
    fasta_validation = _require_mapping(artifact.get("fasta_validation"), "salmon.validation.artifact.fasta_validation")
    _require_equal(fasta_validation.get("duplicate_transcript_ids"), 0, "salmon validation artifact fasta_validation.duplicate_transcript_ids")
    return artifact_path, validation


def _validated_transcriptome_only_salmon(
    root: Path,
    salmon: dict[str, Any],
    genome_fasta: LocalReferenceAsset,
    annotation_gtf: LocalReferenceAsset,
    *,
    identity: dict[str, Any],
) -> tuple[Path, LocalReferenceAsset, dict[str, object], dict[str, object], dict[str, object]]:
    _require_equal(salmon.get("strategy"), SALMON_STRATEGY_TRANSCRIPTOME_ONLY, "salmon.strategy")
    _require_equal(salmon.get("status"), SALMON_BUILT, "salmon.status")
    _require_equal(salmon.get("version"), SALMON_VERSION, "salmon.version")
    _require_equal(salmon.get("kmer_size"), SALMON_KMER_SIZE, "salmon.kmer_size")
    _require_equal(salmon.get("num_decoys"), 0, "salmon.num_decoys")
    _require_equal(salmon.get("decoy_aware"), False, "salmon.decoy_aware")
    transcriptome = _asset_from_payload(
        root, _require_mapping(salmon.get("transcriptome"), "salmon.transcriptome"), "salmon.transcriptome", "salmon_transcriptome"
    )
    index = _resolve_under(root, _require_string(salmon.get("index"), "salmon.index"), "salmon.index", directory=True)
    artifact_path, validation = _validation_artifact(root, salmon, genome_fasta, annotation_gtf, transcriptome, identity)
    index_metadata = _salmon_index_metadata(index)
    expected_metadata = _require_mapping(salmon.get("index_metadata"), "salmon.index_metadata")
    for key in ("seq_hash", "name_hash", "info_index_version", "salmon_index_version", "seq_length", "num_kmers", "num_contigs"):
        _require_equal(expected_metadata.get(key), index_metadata[key], f"salmon.index_metadata.{key}")
    _require_equal(index_metadata["salmon_version"], SALMON_VERSION, "salmon index Salmon version")
    _require_equal(index_metadata["kmer_size"], SALMON_KMER_SIZE, "salmon index k")
    _require_equal(index_metadata["num_decoys"], 0, "salmon index num_decoys")
    provenance = _require_mapping(salmon.get("provenance"), "salmon.provenance")
    _require_equal(provenance.get("mode"), ADOPTED_EXISTING_INDEX, "salmon.provenance.mode")
    _require_string(provenance.get("adopted_at"), "salmon.provenance.adopted_at")
    return index, transcriptome, dict(provenance), dict(index_metadata), validation


def _load_local_reference_root(
    root: Path, manifest_name: str, expected_species: str | None
) -> tuple[LocalReference, dict[str, Any]]:
    if not root.is_dir():
        raise LocalReferenceError(f"Local reference root directory not found: {root}")
    manifest_path = _resolve_under(root, manifest_name, "manifest")
    try:
        loaded = yaml.safe_load(manifest_path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, yaml.YAMLError) as exc:
        raise LocalReferenceError(f"Unable to parse local reference manifest {manifest_path}: {exc}") from exc
    manifest = _require_mapping(loaded, "root")
    manifest_schema_version = manifest.get("schema_version")
    if manifest_schema_version not in LOCAL_REFERENCE_MANIFEST_VERSIONS:
        raise LocalReferenceError(
            "Unsupported local reference manifest schema version: "
            f"{manifest_schema_version!r}; supported versions: "
            f"{LEGACY_LOCAL_REFERENCE_MANIFEST_VERSION}, {LOCAL_REFERENCE_MANIFEST_VERSION}."
        )
    purpose = manifest.get("purpose")
    if manifest_schema_version != LEGACY_LOCAL_REFERENCE_MANIFEST_VERSION:
        purpose = _require_string(purpose, "purpose")
        if purpose not in REFERENCE_PURPOSES:
            raise LocalReferenceError(
                "Local reference manifest purpose must be synthetic_test or production."
            )
    elif purpose is not None:
        raise LocalReferenceError(
                "Legacy local reference manifest schema 1.0 must be migrated explicitly to 1.1 or later before adding purpose."
        )
    identity = _require_mapping(manifest.get("reference"), "reference")
    species = _require_string(identity.get("species"), "reference.species")
    if expected_species is not None and species != expected_species:
        raise LocalReferenceError(
            f"Local reference species {species!r} does not match project organism {expected_species!r}."
        )
    provider = _require_string(identity.get("provider"), "reference.provider")
    release = identity.get("release")
    if not isinstance(release, int) or isinstance(release, bool):
        raise LocalReferenceError("Local reference manifest reference.release must be an integer.")
    assembly = _require_string(identity.get("assembly"), "reference.assembly")
    assembly_patch_value = identity.get("assembly_patch")
    if assembly_patch_value is None:
        assembly_patch = None
    elif not isinstance(assembly_patch_value, str):
        raise LocalReferenceError("Local reference manifest reference.assembly_patch must be a string or null.")
    else:
        # Historic manifests occasionally used an empty placeholder for an
        # assembly without a named patch.  Interpret it as absent without
        # changing the operator-owned manifest on disk.
        assembly_patch = assembly_patch_value.strip() or None
    files = _require_mapping(manifest.get("files"), "files")
    genome_fasta = _asset(root, files, "genome_fasta")
    annotation_gtf = _asset(root, files, "annotation_gtf")
    transcript_fasta = _asset(root, files, "transcript_fasta") if files.get("transcript_fasta") is not None else None
    source_assets: dict[str, LocalReferenceAsset] = {
        "genome_fasta": genome_fasta,
        "annotation_gtf": annotation_gtf,
    }
    if transcript_fasta is not None:
        source_assets["transcript_fasta"] = transcript_fasta
    source_provenance = _validated_source_provenance(
        manifest.get("sources"),
        source_assets,
    )
    salmon = _require_mapping(manifest.get("salmon"), "salmon")
    salmon_status = _require_string(salmon.get("status"), "salmon.status")
    index_value = salmon.get("index")
    salmon_strategy_type: str | None = None
    salmon_transcriptome: LocalReferenceAsset | None = None
    salmon_index_metadata: dict[str, object] | None = None
    salmon_validation: dict[str, object] | None = None
    if index_value is None:
        if salmon_status != SALMON_NOT_BUILT:
            raise LocalReferenceError("Local reference salmon.status must be 'not_built' when salmon.index is null.")
        salmon_index = None
        salmon_provenance = None
    else:
        index_path = _require_string(index_value, "salmon.index")
        salmon_index = _resolve_under(root, index_path, "salmon.index", directory=True)
        if salmon_status != SALMON_BUILT:
            raise LocalReferenceError("Local reference salmon.status must be 'built' when salmon.index is configured.")
        strategy = salmon.get("strategy")
        if strategy is not None and strategy not in SALMON_STRATEGIES:
            raise LocalReferenceError(
                "Local reference salmon.strategy must be one of: "
                f"{SALMON_STRATEGY_TRANSCRIPTOME_ONLY}, {SALMON_STRATEGY_DECOY_AWARE}."
            )
        if _is_prebuilt_salmon_declaration(salmon):
            if transcript_fasta is None:
                raise LocalReferenceError(
                    "A built Salmon declaration requires files.transcript_fasta and its matching source checksum."
                )
            salmon_index, salmon_index_metadata, salmon_provenance = _validated_prebuilt_salmon(
                root, salmon, genome_fasta, transcript_fasta
            )
            salmon_transcriptome = transcript_fasta
            salmon_strategy_type = _require_string(salmon.get("strategy"), "salmon.strategy")
        elif strategy == SALMON_STRATEGY_TRANSCRIPTOME_ONLY:
            if transcript_fasta is None:
                raise LocalReferenceError("A built Salmon declaration requires files.transcript_fasta.")
            (
                salmon_index,
                salmon_transcriptome,
                salmon_provenance,
                salmon_index_metadata,
                salmon_validation,
            ) = _validated_transcriptome_only_salmon(
                root,
                salmon,
                genome_fasta,
                annotation_gtf,
                identity={"species": species, "provider": provider, "release": release, "assembly": assembly},
            )
            salmon_strategy_type = SALMON_STRATEGY_TRANSCRIPTOME_ONLY
        else:
            # Existing manifests predate explicit strategy declaration. They are
            # deterministic legacy decoy-aware manifests when their immutable
            # prepare provenance validates; new prepare output is explicit.
            validate_salmon_index(salmon_index)
            salmon_provenance = _validated_salmon_provenance(
                salmon.get("provenance"), genome_fasta, annotation_gtf
            )
            salmon_strategy_type = SALMON_STRATEGY_DECOY_AWARE
    # HISAT2 was added after the original manifest contract.  Its absence is a
    # valid historical state and means only that this reference is not ready
    # for the alignment/counting backend.
    hisat2 = manifest.get("hisat2", {"status": HISAT2_NOT_BUILT})
    hisat2 = _require_mapping(hisat2, "hisat2")
    hisat2_status = _require_string(hisat2.get("status"), "hisat2.status")
    hisat2_index: Path | None = None
    hisat2_index_prefix: Path | None = None
    hisat2_splice_sites: LocalReferenceAsset | None = None
    hisat2_provenance: dict[str, object] | None = None
    hisat2_strategy: str | None = None
    if hisat2_status == HISAT2_BUILT:
        # ``index_prefix`` is the portable HISAT2 identity.  Retain the
        # directory field for old manifests, but do not require nf-rna to have
        # created an ``hisat2/index`` layout before an existing index can be
        # used.
        index_prefix_value = hisat2.get("index_prefix")
        if index_prefix_value is not None:
            configured_prefix = Path(_require_string(index_prefix_value, "hisat2.index_prefix"))
            if configured_prefix.is_absolute():
                raise LocalReferenceError("Local reference manifest hisat2.index_prefix must be relative to reference.root.")
            hisat2_index_prefix = (root / configured_prefix).resolve()
            try:
                hisat2_index_prefix.relative_to(root)
            except ValueError as exc:
                raise LocalReferenceError("Local reference hisat2.index_prefix escapes reference.root.") from exc
            hisat2_index = hisat2_index_prefix.parent
        else:
            hisat2_index = _resolve_under(root, _require_string(hisat2.get("index"), "hisat2.index"), "hisat2.index", directory=True)
            hisat2_index_prefix = hisat2_index / "genome"
        if hisat2.get("index") is not None:
            declared_index = _resolve_under(root, _require_string(hisat2.get("index"), "hisat2.index"), "hisat2.index", directory=True)
            if declared_index != hisat2_index:
                raise LocalReferenceError("Local reference hisat2.index must be the parent directory of hisat2.index_prefix.")
        splice_payload = _require_mapping(hisat2.get("splice_sites"), "hisat2.splice_sites")
        splice_path = _resolve_under(root, _require_string(splice_payload.get("path"), "hisat2.splice_sites.path"), "hisat2.splice_sites")
        hisat2_splice_sites = _asset_from_payload(root, splice_payload, "hisat2.splice_sites", "hisat2_splice_sites")
        hisat2_provenance = dict(_require_mapping(hisat2.get("provenance", {}), "hisat2.provenance"))
        declared_strategy = hisat2.get("strategy", hisat2_provenance.get("index_strategy"))
        if declared_strategy is None:
            declared_strategy = HISAT2_STRATEGY_LEGACY_ANNOTATION_AWARE
        hisat2_strategy = _require_string(declared_strategy, "hisat2.strategy")
        if hisat2_strategy not in HISAT2_STRATEGIES:
            raise LocalReferenceError(f"Local reference hisat2.strategy is unsupported: {hisat2_strategy!r}.")
        prefix = hisat2_index_prefix.name
        validate_hisat2_index(hisat2_index, prefix=prefix)
        # New prebuilt declarations keep scientific asset identity at the
        # backend level; historic builder manifests keep it in provenance.
        genome_checksum = hisat2.get("genome_fasta_sha256", hisat2_provenance.get("genome_fasta_sha256"))
        gtf_checksum = hisat2.get("source_gtf_sha256", hisat2_provenance.get("annotation_gtf_sha256"))
        _require_equal(genome_checksum, genome_fasta.sha256, "hisat2.genome_fasta_sha256")
        _require_equal(gtf_checksum, annotation_gtf.sha256, "hisat2.source_gtf_sha256")
        hisat2_provenance.setdefault("genome_fasta_sha256", genome_fasta.sha256)
        hisat2_provenance.setdefault("annotation_gtf_sha256", annotation_gtf.sha256)
        builder_version = _require_string(
            hisat2.get("version", hisat2_provenance.get("index_builder_version", hisat2_provenance.get("hisat2_version"))),
            "hisat2.provenance.index_builder_version",
        )
        if not SEMVER_PATTERN.fullmatch(builder_version):
            raise LocalReferenceError("Local reference HISAT2 version must be a numeric release, for example 2.2.3.")
        if builder_version != HISAT2_VERSION:
            raise LocalReferenceError(
                f"Local reference HISAT2 index builder requires exactly {HISAT2_VERSION}; "
                f"observed {builder_version!r}. Create or update the documented reference-builder environment."
            )
        hisat2_provenance.setdefault("index_builder_version", builder_version)
        if _uses_runtime_splice_sites(hisat2_strategy):
            _require_equal(
                hisat2.get("splice_sites_gtf_sha256", hisat2_provenance.get("splice_sites_derived_from_gtf_sha256", gtf_checksum)), annotation_gtf.sha256,
                "hisat2.provenance.splice_sites_derived_from_gtf_sha256",
            )
            # A version match proves only the pinned binary contract.  It does
            # not replace a recorded FASTQ acceptance/smoke validation for an
            # existing prebuilt index.
            default_compatibility = HISAT2_RUNTIME_COMPATIBILITY_SMOKE_REQUIRED
            compatibility = _require_string(hisat2.get("runtime_compatibility", hisat2_provenance.get("runtime_compatibility", default_compatibility)), "hisat2.runtime_compatibility")
            if compatibility not in {HISAT2_RUNTIME_COMPATIBILITY_VALIDATED, HISAT2_RUNTIME_COMPATIBILITY_SMOKE_REQUIRED}:
                raise LocalReferenceError("Local reference hisat2.provenance.runtime_compatibility is unsupported.")
            hisat2_provenance.setdefault("runtime_compatibility", compatibility)
    elif hisat2_status != HISAT2_NOT_BUILT:
        raise LocalReferenceError("Local reference hisat2.status must be 'not_built' or 'built'.")
    return LocalReference(
        root=root,
        manifest_path=manifest_path,
        manifest_sha256=sha256_file(manifest_path),
        manifest_schema_version=manifest_schema_version,
        purpose=purpose,
        species=species,
        provider=provider,
        release=release,
        assembly=assembly,
        assembly_patch=assembly_patch,
        genome_fasta=genome_fasta,
        annotation_gtf=annotation_gtf,
        transcript_fasta=transcript_fasta,
        salmon_index=salmon_index,
        salmon_status=salmon_status,
        salmon_strategy_type=salmon_strategy_type,
        salmon_provenance=salmon_provenance,
        salmon_transcriptome=salmon_transcriptome,
        salmon_index_metadata=salmon_index_metadata,
        salmon_validation=salmon_validation,
        hisat2_index=hisat2_index,
        hisat2_index_prefix=hisat2_index_prefix,
        hisat2_splice_sites=hisat2_splice_sites,
        hisat2_status=hisat2_status,
        hisat2_strategy=hisat2_strategy,
        hisat2_provenance=hisat2_provenance,
        source_provenance=source_provenance,
    ), manifest


def load_local_reference(reference: ReferenceConfig, project_species: str) -> LocalReference:
    """Resolve a project-configured local reference without modifying it."""

    if reference.source != "local":
        raise LocalReferenceError("Local reference loading requires reference.source: local.")
    assert reference.root is not None and reference.manifest is not None
    root_configured = Path(reference.root).expanduser()
    if not root_configured.is_absolute():
        raise LocalReferenceError("reference.root must be an absolute path for source=local.")
    local_reference, _manifest = _load_local_reference_root(
        root_configured.resolve(), reference.manifest, project_species
    )
    return local_reference


def load_local_reference_root(reference_root: Path | str) -> LocalReference:
    """Load a root for the standalone ``rnaseq reference prepare`` command."""

    root = Path(reference_root).expanduser()
    if not root.is_absolute():
        raise LocalReferenceError("Reference root must be an absolute path.")
    local_reference, _manifest = _load_local_reference_root(root.resolve(), "reference_manifest.yaml", None)
    return local_reference


def reference_registry_path() -> Path:
    """Return the machine-local registry path without creating it.

    The registry deliberately records only locations and manifest-derived
    identity metadata.  The manifest remains the source of truth for all
    scientific assets and checksums when a reference is later selected.
    """

    xdg_config_home = os.environ.get("XDG_CONFIG_HOME")
    if xdg_config_home:
        return Path(xdg_config_home).expanduser() / "nf-rna" / "references.yaml"
    if platform.system() == "Darwin":
        return Path.home() / "Library" / "Application Support" / "nf-rna" / "references.yaml"
    return Path.home() / ".config" / "nf-rna" / "references.yaml"


def _resolved_registry_path(registry_path: Path | str | None) -> Path:
    path = Path(registry_path).expanduser() if registry_path is not None else reference_registry_path()
    return path.resolve()


def _read_reference_registry(path: Path) -> list[dict[str, Any]]:
    if not path.exists():
        return []
    if not path.is_file():
        raise LocalReferenceError(f"Managed-reference registry is not a file: {path}")
    try:
        loaded = yaml.safe_load(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, yaml.YAMLError) as exc:
        raise LocalReferenceError(f"Unable to parse managed-reference registry {path}: {exc}") from exc
    payload = _require_mapping(loaded, "registry")
    if payload.get("schema_version") != REFERENCE_REGISTRY_SCHEMA_VERSION:
        raise LocalReferenceError(
            "Unsupported managed-reference registry schema version: "
            f"{payload.get('schema_version')!r}."
        )
    entries = payload.get("references")
    if not isinstance(entries, list):
        raise LocalReferenceError("Managed-reference registry references must be a list.")
    return [_require_mapping(entry, "registry.references entry") for entry in entries]


def _registry_entry(reference: LocalReference) -> dict[str, object]:
    """Persist a pointer and display identity, never a second asset manifest."""

    return {
        "root": str(reference.root),
        "manifest": str(reference.manifest_path.relative_to(reference.root)),
        "identity": {
            "species": reference.species,
            "provider": reference.provider,
            "release": reference.release,
            "assembly": reference.assembly,
            "assembly_patch": reference.assembly_patch,
        },
        "purpose": reference.purpose,
        "manifest_sha256": reference.manifest_sha256,
    }


def _write_reference_registry(path: Path, entries: list[dict[str, Any]]) -> None:
    """Atomically publish a private, user-local registry file."""

    path.parent.mkdir(mode=0o700, parents=True, exist_ok=True)
    payload = {
        "schema_version": REFERENCE_REGISTRY_SCHEMA_VERSION,
        "references": entries,
    }
    temporary = path.with_name(f".{path.name}.{uuid.uuid4().hex}.tmp")
    try:
        temporary.write_text(yaml.safe_dump(payload, sort_keys=False, allow_unicode=True), encoding="utf-8")
        os.chmod(temporary, 0o600)
        os.replace(temporary, path)
    finally:
        temporary.unlink(missing_ok=True)


def register_local_reference(
    reference_root: Path | str, *, registry_path: Path | str | None = None
) -> LocalReference:
    """Validate a managed reference and record it once for this workstation."""

    reference = load_local_reference_root(reference_root)
    path = _resolved_registry_path(registry_path)
    entries = _read_reference_registry(path)
    entry = _registry_entry(reference)
    existing_index = next(
        (
            index for index, existing in enumerate(entries)
            if existing.get("root") == entry["root"] and existing.get("manifest") == entry["manifest"]
        ),
        None,
    )
    if existing_index is None:
        entries.append(entry)
    else:
        # Re-registration is also the intentional way to refresh the cached
        # display identity after an index is prepared or a manifest is revised.
        entries[existing_index] = entry
    _write_reference_registry(path, entries)
    return reference


def registered_local_references(
    *, registry_path: Path | str | None = None
) -> list[LocalReference]:
    """Load every still-valid registered reference through the manifest loader.

    A moved, deleted, or now-invalid registration is ignored here so that an
    ordinary interactive project creation can still use its manual and custom
    reference escape hatches.  Registration itself remains strict.
    """

    entries = _read_reference_registry(_resolved_registry_path(registry_path))
    references: list[LocalReference] = []
    seen: set[tuple[Path, Path]] = set()
    for entry in entries:
        try:
            root_value = _require_string(entry.get("root"), "registry.references.root")
            manifest = _require_string(entry.get("manifest"), "registry.references.manifest")
            root = Path(root_value).expanduser()
            if not root.is_absolute():
                raise LocalReferenceError("Managed-reference registry roots must be absolute paths.")
            reference, _manifest = _load_local_reference_root(root.resolve(), manifest, None)
            key = (reference.root, reference.manifest_path)
            if key not in seen:
                references.append(reference)
                seen.add(key)
        except LocalReferenceError:
            continue
    return references


def compatible_registered_references(
    species: str, backend: str, *, registry_path: Path | str | None = None
) -> list[LocalReference]:
    """Return deterministic, production-ready registry matches for one backend."""

    if backend not in {"salmon", "hisat2_featurecounts"}:
        raise LocalReferenceError(f"Unsupported managed-reference backend: {backend!r}.")
    compatible = [
        reference
        for reference in registered_local_references(registry_path=registry_path)
        if reference.purpose == "production"
        and reference.species == species
        and (
            reference.salmon_index is not None
            if backend == "salmon"
            else reference.hisat2_index is not None and reference.hisat2_runtime_ready
        )
    ]
    return sorted(
        compatible,
        key=lambda reference: (
            reference.provider.casefold(),
            -reference.release,
            reference.assembly.casefold(),
            (reference.assembly_patch or "").casefold(),
            str(reference.root),
        ),
    )


def _run_reference_command(command: list[str], runner: Callable[..., subprocess.CompletedProcess[str]], *, cwd: Path) -> subprocess.CompletedProcess[str]:
    """Run one host-native reference-builder argv without a shell."""

    try:
        completed = runner(command, check=False, text=True, capture_output=True, cwd=cwd)
    except OSError as exc:
        raise ReferencePreparationError(f"Unable to start host-native reference preparation command {command[0]!r}: {exc}") from exc
    if completed.returncode != 0:
        detail = (completed.stderr or completed.stdout or "no diagnostic output").strip()
        raise ReferencePreparationError(
            f"Reference preparation command failed (exit {completed.returncode}): {detail}"
        )
    return completed


def _host_tool(
    name: str, *, route: str, version_pattern: re.Pattern[str] | None,
    expected: str, runner: Callable[..., subprocess.CompletedProcess[str]],
    resolver: Callable[[str], str | None], cwd: Path,
) -> dict[str, object]:
    """Resolve, version-check and describe one required host executable."""

    configured = resolver(name)
    if not configured:
        raise ReferencePreparationError(
            f"Missing required host executable {name!r} for {route} preparation; expected {expected}. "
            "Create the documented environment with: mamba env create -f environment.reference-builder.yml"
        )
    executable = str(Path(configured).resolve())
    try:
        result = runner([executable, "--version"], check=False, text=True, capture_output=True, cwd=cwd)
    except OSError as exc:
        raise ReferencePreparationError(f"Cannot obtain version for host executable {name!r}: {exc}") from exc
    raw = ((result.stdout or "") + ("\n" + result.stderr if result.stderr else "")).strip()
    if result.returncode != 0:
        raise ReferencePreparationError(
            f"Host executable {name!r} for {route} does not support a successful --version check; "
            f"expected {expected}. Output: {raw or 'none'}"
        )
    match = version_pattern.search(raw) if version_pattern else None
    if version_pattern is not None and match is None:
        raise ReferencePreparationError(
            f"Unable to parse version for host executable {name!r} used by {route}; expected {expected}. Output: {raw or 'none'}"
        )
    version = match.group(1) if match is not None else raw
    required_versions = {
        "salmon": SALMON_VERSION,
        "hisat2-build": HISAT2_VERSION,
    }
    if name in required_versions and version != required_versions[name]:
        expected_version = required_versions[name]
        raise ReferencePreparationError(
            f"Unsupported {name} version {version!r} for {route}; expected {expected_version}. "
            "Create the documented environment with: mamba env create -f environment.reference-builder.yml"
        )
    return {"name": name, "executable": executable, "version": version, "version_output": raw}


def _hisat2_splice_site_helper(
    *, route: str, hisat2_build: dict[str, object], runner: Callable[..., subprocess.CompletedProcess[str]],
    resolver: Callable[[str], str | None], cwd: Path,
) -> dict[str, object]:
    """Validate the packaged splice-site helper without inventing a version API.

    ``hisat2_extract_splice_sites.py`` exposes ``-h`` and ``-v`` (verbose),
    but no semantic-version option.  Tie it to the exact checked
    ``hisat2-build`` installation by location and verify its documented help
    contract instead of treating a verbose flag as a version report.
    """

    name = "hisat2_extract_splice_sites.py"
    configured = resolver(name)
    if not configured:
        raise ReferencePreparationError(
            f"Missing required host executable {name!r} for {route}. "
            "Create the documented environment with: mamba env create -f environment.reference-builder.yml"
        )
    executable = str(Path(configured).resolve())
    build_executable = Path(str(hisat2_build["executable"])).resolve()
    if Path(executable).parent != build_executable.parent:
        raise ReferencePreparationError(
            f"Host executable {name!r} must resolve beside the validated hisat2-build executable "
            f"for {route}; observed {executable}, hisat2-build {build_executable}."
        )
    try:
        result = runner([executable, "-h"], check=False, text=True, capture_output=True, cwd=cwd)
    except OSError as exc:
        raise ReferencePreparationError(f"Unable to inspect host executable {name!r}: {exc}") from exc
    raw = ((result.stdout or "") + ("\n" + result.stderr if result.stderr else "")).strip()
    if result.returncode != 0 or "Extract splice junctions from a GTF file" not in raw or "gtf_file" not in raw:
        raise ReferencePreparationError(
            f"Host executable {name!r} does not provide the expected HISAT2 splice-site helper help contract "
            f"for {route}. Output: {raw or 'none'}"
        )
    return {
        "name": name,
        "executable": executable,
        "validation": "help_contract",
        "help_output": raw,
        "associated_hisat2_build_version": hisat2_build["version"],
    }


def _reference_build_preflight(reference: LocalReference, *, route: str, threads: int) -> dict[str, object]:
    """Validate cheap host/filesystem facts before creating builder staging."""

    if threads < 1:
        raise ReferencePreparationError("Reference preparation threads must be at least 1.")
    for asset in reference.assets():
        if not os.access(asset.path, os.R_OK):
            raise ReferencePreparationError(f"Reference source asset is not readable: {asset.path}")
    if not os.access(reference.root, os.W_OK):
        raise ReferencePreparationError(f"Reference destination is not writable: {reference.root}")
    target_parent = reference.root / ("salmon" if route.startswith("Salmon") else "hisat2")
    # Staging is deliberately a sibling below reference.root. An existing
    # backend directory may itself be a mount point, however, in which case a
    # rename into it would cross filesystems and must be refused before build.
    if target_parent.exists() and reference.root.stat().st_dev != target_parent.stat().st_dev:
        raise ReferencePreparationError(
            f"Reference staging root and final {target_parent.name} destination are on different filesystems; "
            "cannot guarantee atomic publication."
        )
    usage = shutil.disk_usage(target_parent if target_parent.exists() else reference.root)
    return {
        "route": route,
        "destination": str(reference.root),
        "threads": threads,
        "logical_cpus": os.cpu_count(),
        "available_memory_bytes": _available_memory_bytes(),
        "available_disk_bytes": usage.free,
        "operating_system": platform.system() or "unknown",
        "architecture": platform.machine() or "unknown",
    }


def _assert_host_owned(paths: tuple[Path, ...]) -> None:
    """Reject Linux/WSL outputs not owned by the process that built them.

    Native executables should create files as the invoking account.  Checking
    this explicitly catches an accidental privilege wrapper without weakening
    permissions or attempting to repair ownership after publication.
    """

    if platform.system().lower() != "linux" or not hasattr(os, "getuid"):
        return
    expected_uid = os.getuid()
    unexpected: list[Path] = []
    for path in paths:
        candidates = (path, *path.rglob("*")) if path.is_dir() else (path,)
        for candidate in candidates:
            if candidate.is_file() and candidate.stat().st_uid != expected_uid:
                unexpected.append(candidate)
    if unexpected:
        raise ReferencePreparationError(
            "Host-native reference preparation produced files not owned by the invoking user: "
            + ", ".join(str(path) for path in unexpected[:3])
        )


def _available_memory_bytes() -> int | None:
    try:
        return os.sysconf("SC_PAGE_SIZE") * os.sysconf("SC_AVPHYS_PAGES")
    except (AttributeError, OSError, ValueError):
        return None


def _write_decoys(genome: Path, destination: Path) -> None:
    """Write exact FASTA record IDs for the decoy-aware gentrome without a shell."""

    with genome.open(encoding="utf-8") as source, destination.open("w", encoding="utf-8", newline="\n") as target:
        for line in source:
            if line.startswith(">"):
                identifier = line[1:].strip().split(maxsplit=1)[0]
                if not identifier:
                    raise ReferencePreparationError(f"Genome FASTA contains a blank record identifier: {genome}")
                target.write(identifier + "\n")


def _reference_source_assets(reference: LocalReference) -> dict[str, dict[str, str]]:
    """Return biological input identity separately from host-builder facts."""

    return {
        "genome_fasta": {
            "path": reference.genome_fasta.path.relative_to(reference.root).as_posix(),
            "sha256": reference.genome_fasta.sha256,
        },
        "annotation_gtf": {
            "path": reference.annotation_gtf.path.relative_to(reference.root).as_posix(),
            "sha256": reference.annotation_gtf.sha256,
        },
    }


def _host_builder_provenance(
    tools: dict[str, dict[str, object]], preflight: dict[str, object], arguments: dict[str, list[str]],
) -> dict[str, object]:
    """Keep host-specific reproducibility facts out of biological identity."""

    return {
        "mode": "host_native",
        "tools": tools,
        "arguments": arguments,
        "threads": preflight["threads"],
        "operating_system": preflight["operating_system"],
        "architecture": preflight["architecture"],
        "preflight": preflight,
    }


def _concat_files(destination: Path, sources: tuple[Path, ...]) -> None:
    with destination.open("wb") as target:
        for source in sources:
            with source.open("rb") as handle:
                shutil.copyfileobj(handle, target)


def _write_manifest_atomically(path: Path, manifest: dict[str, Any]) -> None:
    temporary = path.with_name(f".{path.name}.{uuid.uuid4().hex}.tmp")
    try:
        temporary.write_text(
            yaml.safe_dump(manifest, sort_keys=False, allow_unicode=True), encoding="utf-8", newline="\n"
        )
        temporary.replace(path)
    finally:
        if temporary.exists():
            temporary.unlink()


def adopt_local_salmon_index(
    reference_root: Path | str,
    *,
    index: str,
    transcriptome: str | None = None,
    strategy: str,
    validation_artifact: str | None = None,
    now: Callable[[], datetime] = lambda: datetime.now(UTC),
) -> LocalReference:
    """Validate and atomically register an already-built Salmon index.

    This is deliberately explicit: the caller selects every asset and strategy;
    no index is discovered, rebuilt, or inferred from its contents.
    """

    root = Path(reference_root).expanduser()
    if not root.is_absolute():
        raise ReferenceAdoptionError("Reference root must be an absolute path.")
    root = root.resolve()
    try:
        reference, manifest = _load_local_reference_root(root, "reference_manifest.yaml", None)
    except LocalReferenceError as exc:
        raise ReferenceAdoptionError(str(exc)) from exc
    if strategy not in SALMON_STRATEGIES:
        raise ReferenceAdoptionError(
            "Unsupported Salmon strategy: "
            f"{strategy!r}. Supported strategies: {SALMON_STRATEGY_TRANSCRIPTOME_ONLY}, {SALMON_STRATEGY_DECOY_AWARE}."
        )
    if reference.salmon_status != SALMON_NOT_BUILT:
        raise ReferenceAdoptionError(
            "Refusing to replace an existing local reference Salmon declaration; "
            f"current status is {reference.salmon_status!r}."
        )
    if reference.transcript_fasta is None:
        raise ReferenceAdoptionError(
            "Adopting a Salmon index requires files.transcript_fasta so its source SHA256 can be verified."
        )
    try:
        index_path = _resolve_under(root, index, "adoption index", directory=True)
        transcriptome_path = (
            _resolve_under(root, transcriptome, "adoption transcriptome")
            if transcriptome is not None else reference.transcript_fasta.path
        )
    except LocalReferenceError as exc:
        raise ReferenceAdoptionError(str(exc)) from exc
    try:
        index_metadata = _salmon_index_metadata(index_path)
    except LocalReferenceError as exc:
        raise ReferenceAdoptionError(str(exc)) from exc
    transcriptome_sha256 = sha256_file(transcriptome_path)
    # Keep the former, stricter GTF-derived adoption record readable for
    # already scripted users.  New adoption (no validation artifact) is the
    # lightweight prebuilt-index contract below.
    if validation_artifact is not None and transcriptome is not None and transcriptome_sha256 != reference.transcript_fasta.sha256:
        try:
            artifact_path = _resolve_under(root, validation_artifact, "adoption validation artifact")
        except LocalReferenceError as exc:
            raise ReferenceAdoptionError(str(exc)) from exc
        candidate = deepcopy(manifest)
        candidate["salmon"] = {
            "status": SALMON_BUILT, "strategy": SALMON_STRATEGY_TRANSCRIPTOME_ONLY,
            "index": index_path.relative_to(root).as_posix(),
            "transcriptome": {"path": transcriptome_path.relative_to(root).as_posix(), "sha256": transcriptome_sha256, "source": "gtf_derived"},
            "version": index_metadata["salmon_version"], "kmer_size": index_metadata["kmer_size"],
            "num_decoys": index_metadata["num_decoys"], "decoy_aware": False,
            "index_metadata": {key: index_metadata[key] for key in ("seq_hash", "name_hash", "info_index_version", "salmon_index_version", "seq_length", "num_kmers", "num_contigs")},
            "validation": {"artifact": artifact_path.relative_to(root).as_posix(), "transcript_id_contract": "exact", "fasta_only_transcripts": 0, "zero_gene_mappings": 0, "multi_gene_mappings": 0},
            "provenance": {"mode": ADOPTED_EXISTING_INDEX, "adopted_at": now().astimezone(UTC).replace(microsecond=0).isoformat()},
        }
        try:
            _validated_transcriptome_only_salmon(root, _require_mapping(candidate["salmon"], "salmon"), reference.genome_fasta, reference.annotation_gtf, identity={"species": reference.species, "provider": reference.provider, "release": reference.release, "assembly": reference.assembly})
        except LocalReferenceError as exc:
            raise ReferenceAdoptionError(str(exc)) from exc
        _write_manifest_atomically(reference.manifest_path, candidate)
        return load_local_reference_root(root)
    if transcriptome_sha256 != reference.transcript_fasta.sha256:
        raise ReferenceAdoptionError(
            "The adopted Salmon transcriptome must match files.transcript_fasta by SHA256; "
            "update and validate the manifest source asset first."
        )
    expected_decoys = 0 if strategy == SALMON_STRATEGY_TRANSCRIPTOME_ONLY else None
    if expected_decoys is not None and index_metadata["num_decoys"] != expected_decoys:
        raise ReferenceAdoptionError("The adopted transcriptome_only Salmon index contains decoys.")
    if strategy == SALMON_STRATEGY_DECOY_AWARE and (
        not isinstance(index_metadata["num_decoys"], int) or index_metadata["num_decoys"] <= 0
    ):
        raise ReferenceAdoptionError("The adopted decoy_aware Salmon index declares no decoys.")
    candidate = deepcopy(manifest)
    candidate["salmon"] = {
        "status": SALMON_BUILT,
        "strategy": strategy,
        "index": index_path.relative_to(root).as_posix(),
        "version": index_metadata["salmon_version"],
        "source_transcriptome_sha256": transcriptome_sha256,
        "provenance": {
            "mode": ADOPTED_EXISTING_INDEX,
            "adopted_at": now().astimezone(UTC).replace(microsecond=0).isoformat(),
        },
    }
    if strategy == SALMON_STRATEGY_DECOY_AWARE:
        candidate["salmon"]["source_genome_sha256"] = reference.genome_fasta.sha256
    if validation_artifact is not None:
        try:
            artifact_path = _resolve_under(root, validation_artifact, "adoption validation artifact")
        except LocalReferenceError as exc:
            raise ReferenceAdoptionError(str(exc)) from exc
        candidate["salmon"]["provenance"]["validation_artifact"] = artifact_path.relative_to(root).as_posix()
    try:
        _validated_prebuilt_salmon(root, _require_mapping(candidate["salmon"], "salmon"), reference.genome_fasta, reference.transcript_fasta)
    except LocalReferenceError as exc:
        raise ReferenceAdoptionError(str(exc)) from exc
    _write_manifest_atomically(reference.manifest_path, candidate)
    return load_local_reference_root(root)


def prepare_local_reference(
    reference_root: Path | str,
    *,
    threads: int = 4,
    runner: Callable[..., subprocess.CompletedProcess[str]] = subprocess.run,
    tool_resolver: Callable[[str], str | None] = shutil.which,
    now: Callable[[], datetime] = lambda: datetime.now(UTC),
) -> LocalReference:
    """Optionally build a host-native decoy-aware index from registered assets."""

    reference, manifest = _load_local_reference_root(
        Path(reference_root).expanduser().resolve(), "reference_manifest.yaml", None
    )
    if reference.salmon_status == SALMON_BUILT:
        raise ReferencePreparationError(
            f"Local reference Salmon index is already built: {reference.salmon_index}."
        )
    if reference.salmon_status != SALMON_NOT_BUILT:
        raise ReferencePreparationError(f"Unsupported local reference Salmon status: {reference.salmon_status}.")
    if reference.transcript_fasta is None:
        raise ReferencePreparationError(
            "Host-native Salmon preparation requires files.transcript_fasta; "
            "register the source transcriptome asset and SHA256 first."
        )

    salmon_root = reference.root / "salmon"
    final_index = salmon_root / "index"
    if final_index.exists():
        raise ReferencePreparationError(
            f"Refusing to overwrite existing local reference index directory: {final_index}."
        )
    preflight = _reference_build_preflight(reference, route="Salmon decoy-aware", threads=threads)
    tools = {
        "salmon": _host_tool(
            "salmon", route="Salmon decoy-aware", version_pattern=re.compile(r"(?:salmon\s+|version\s+)([0-9]+\.[0-9]+\.[0-9]+)", re.I),
            expected=f"Salmon {SALMON_VERSION}", runner=runner, resolver=tool_resolver, cwd=reference.root,
        ),
    }
    temporary_root = reference.root / f".rnaseq-reference-prepare-{uuid.uuid4().hex}"
    salmon_dir = temporary_root / "salmon"
    temporary_index = salmon_dir / "index"
    temporary_root.mkdir()
    salmon_dir.mkdir()
    decoys = salmon_dir / "decoys.txt"
    gentrome = salmon_dir / "gentrome.fa"
    salmon_command = [
        str(tools["salmon"]["executable"]), "index", "--threads", str(threads), "-t", str(gentrome),
        "-d", str(decoys), "-i", str(temporary_index), "-k", str(SALMON_KMER_SIZE),
    ]
    try:
        _write_decoys(reference.genome_fasta.path, decoys)
        _concat_files(gentrome, (reference.transcript_fasta.path, reference.genome_fasta.path))
        _run_reference_command(salmon_command, runner, cwd=temporary_root)
        validate_salmon_index(temporary_index)
        _assert_host_owned((temporary_index,))
        salmon_root.mkdir(exist_ok=True)
        temporary_index.replace(final_index)
        manifest["salmon"] = {
            "index": "salmon/index",
            "status": SALMON_BUILT,
            "strategy": SALMON_STRATEGY_DECOY_AWARE,
            "version": SALMON_VERSION,
            "source_transcriptome_sha256": reference.transcript_fasta.sha256,
            "source_genome_sha256": reference.genome_fasta.sha256,
            "provenance": {
                "builder": _host_builder_provenance(
                    tools, preflight, {"salmon_index": salmon_command}
                ),
                "salmon_version": SALMON_VERSION,
                "index_path": "salmon/index",
                "source_assets": _reference_source_assets(reference),
                "source_transcriptome_sha256": reference.transcript_fasta.sha256,
                "genome_fasta_sha256": reference.genome_fasta.sha256,
                "annotation_gtf_sha256": reference.annotation_gtf.sha256,
                "transcriptome_strategy": "manifest_registered_transcript_fasta",
                "decoy_strategy": DECOY_STRATEGY,
                "kmer_size": SALMON_KMER_SIZE,
                "commands": {"salmon_index": salmon_command},
                "index_validation": "required Salmon index artifacts present and non-empty",
                "built_at": now().astimezone(UTC).replace(microsecond=0).isoformat(),
            },
        }
        try:
            _write_manifest_atomically(reference.manifest_path, manifest)
        except Exception:
            # A manifest write failure must not strand an unregistered index in
            # the final location. Put the newly built artifact back into the
            # retained staging tree, leaving the old manifest untouched.
            if final_index.exists() and not temporary_index.exists():
                final_index.replace(temporary_index)
            raise
    except Exception:
        # Keep failed staging for an operator to inspect; publication has not
        # happened because final_index is renamed only after validation.
        raise
    else:
        shutil.rmtree(temporary_root, ignore_errors=True)
    return load_local_reference_root(reference.root)


def prepare_local_hisat2_reference(
    reference_root: Path | str,
    *,
    threads: int = 4,
    runner: Callable[..., subprocess.CompletedProcess[str]] = subprocess.run,
    tool_resolver: Callable[[str], str | None] = shutil.which,
    now: Callable[[], datetime] = lambda: datetime.now(UTC),
) -> LocalReference:
    """Build and atomically register the HISAT2 genome index for one reference.

    A complete index is moved into place only after every numbered file and
    generated splice-site file has been checked.  The final-directory refusal
    is also the concurrency guard: a competing completed preparation is never
    overwritten.
    """

    reference, manifest = _load_local_reference_root(Path(reference_root).expanduser().resolve(), "reference_manifest.yaml", None)
    if reference.hisat2_status == HISAT2_BUILT:
        raise ReferencePreparationError(f"Local reference HISAT2 index is already built: {reference.hisat2_index}.")
    final_root = reference.root / "hisat2"
    final_index = final_root / "index"
    if final_index.exists():
        raise ReferencePreparationError(f"Refusing to overwrite existing local HISAT2 index directory: {final_index}.")
    preflight = _reference_build_preflight(reference, route="HISAT2 genome-only runtime splice-sites", threads=threads)
    tools = {
        "hisat2-build": _host_tool(
            "hisat2-build", route="HISAT2 genome-only runtime splice-sites", version_pattern=re.compile(r"(?:version\s+)?([0-9]+\.[0-9]+\.[0-9]+)", re.I),
            expected=f"HISAT2 {HISAT2_VERSION}", runner=runner, resolver=tool_resolver, cwd=reference.root,
        ),
    }
    tools["hisat2_extract_splice_sites.py"] = _hisat2_splice_site_helper(
        route="HISAT2 genome-only runtime splice-sites",
        hisat2_build=tools["hisat2-build"],
        runner=runner,
        resolver=tool_resolver,
        cwd=reference.root,
    )
    temporary = reference.root / f".rnaseq-hisat2-prepare-{uuid.uuid4().hex}"
    temporary_index = temporary / "index"
    temporary.mkdir()
    temporary_index.mkdir()
    try:
        splice_command = [str(tools["hisat2_extract_splice_sites.py"]["executable"]), str(reference.annotation_gtf.path)]
        try:
            splice_result = runner(splice_command, check=False, text=True, capture_output=True, cwd=temporary)
        except OSError as exc:
            raise ReferencePreparationError(f"Unable to start host-native HISAT2 splice-site helper: {exc}") from exc
        if splice_result.returncode != 0:
            raise ReferencePreparationError(f"HISAT2 splice-site command failed: {(splice_result.stderr or splice_result.stdout).strip()}")
        splice = temporary / "splice_sites.txt"
        splice.write_text(splice_result.stdout or "", encoding="utf-8", newline="\n")
        index_command = [
            str(tools["hisat2-build"]["executable"]), "--threads", str(threads),
            str(reference.genome_fasta.path), str(temporary_index / "genome"),
        ]
        _run_reference_command(index_command, runner, cwd=temporary)
        validate_hisat2_index(temporary_index)
        # A valid single-exon annotation has no junctions.  HISAT2 accepts an
        # empty --ss file, so preserve it as the explicit, reproducible result
        # rather than rejecting an otherwise valid reference.
        if not splice.is_file():
            raise ReferencePreparationError("HISAT2 splice-site extraction produced no artifact.")
        _assert_host_owned((temporary_index, splice))
        final_root.mkdir(exist_ok=True)
        temporary_index.replace(final_index)
        final_splice = final_root / "splice_sites.txt"
        splice.replace(final_splice)
        manifest["hisat2"] = {
            "status": HISAT2_BUILT,
            "index": "hisat2/index",
            "index_prefix": "hisat2/index/genome",
            "strategy": HISAT2_STRATEGY_GENOME_ONLY_RUNTIME_SPLICES,
            "version": HISAT2_VERSION,
            "genome_fasta_sha256": reference.genome_fasta.sha256,
            "source_gtf_sha256": reference.annotation_gtf.sha256,
            "splice_sites_gtf_sha256": reference.annotation_gtf.sha256,
            "runtime_compatibility": HISAT2_RUNTIME_COMPATIBILITY_VALIDATED,
            "splice_sites": {"path": "hisat2/splice_sites.txt", "sha256": sha256_file(final_splice)},
            "provenance": {
                "builder": _host_builder_provenance(
                    tools, preflight, {"splice_sites": splice_command, "hisat2_build": index_command}
                ),
                "hisat2_version": HISAT2_VERSION,
                "index_builder_version": HISAT2_VERSION,
                "runtime_aligner_version": HISAT2_RUNTIME_VERSION,
                "source_assets": _reference_source_assets(reference),
                "genome_fasta_sha256": reference.genome_fasta.sha256,
                "annotation_gtf_sha256": reference.annotation_gtf.sha256,
                "threads": threads,
                "index_strategy": HISAT2_STRATEGY_GENOME_ONLY_RUNTIME_SPLICES,
                "commands": {"splice_sites": splice_command, "hisat2_build": index_command},
                "index_validation": "complete numbered .ht2 or .ht2l family",
                "built_at": now().astimezone(UTC).replace(microsecond=0).isoformat(),
            },
        }
        try:
            _write_manifest_atomically(reference.manifest_path, manifest)
        except Exception:
            if final_splice.exists() and not splice.exists():
                final_splice.replace(splice)
            if final_index.exists() and not temporary_index.exists():
                final_index.replace(temporary_index)
            raise
    except Exception:
        raise
    else:
        shutil.rmtree(temporary, ignore_errors=True)
    return load_local_reference_root(reference.root)
