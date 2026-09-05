"""Versioned local-reference validation and one-time Salmon index preparation."""

from __future__ import annotations

import hashlib
import json
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
LOCAL_REFERENCE_MANIFEST_VERSION = "1.0"
SALMON_NOT_BUILT = "not_built"
SALMON_BUILT = "built"
SALMON_VERSION = "1.10.3"
SALMON_KMER_SIZE = 31
TRANSCRIPTOME_STRATEGY = "nfcore_rnaseq_3.26.0_rsem_from_genome_fasta_and_gtf"
DECOY_STRATEGY = "nfcore_rnaseq_3.26.0_gentrome"
SALMON_STRATEGY_TRANSCRIPTOME_ONLY = "transcriptome_only"
SALMON_STRATEGY_DECOY_AWARE = "decoy_aware"
SALMON_STRATEGIES = frozenset((SALMON_STRATEGY_TRANSCRIPTOME_ONLY, SALMON_STRATEGY_DECOY_AWARE))
ADOPTED_EXISTING_INDEX = "adopted_existing_index"
RSEM_IMAGE = "community.wave.seqera.io/library/rsem_star:5acb4e8c03239c32"
SALMON_IMAGE = "quay.io/biocontainers/salmon:1.10.3--h6dccd9a_2"
HISAT2_VERSION = "2.2.1"
HISAT2_IMAGE = "quay.io/biocontainers/hisat2:2.2.1--h87f3376_4"
HISAT2_NOT_BUILT = "not_built"
HISAT2_BUILT = "built"
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
    species: str
    provider: str
    release: int
    assembly: str
    assembly_patch: str
    genome_fasta: LocalReferenceAsset
    annotation_gtf: LocalReferenceAsset
    transcript_fasta: LocalReferenceAsset
    salmon_index: Path | None
    salmon_status: str
    salmon_strategy_type: str | None
    salmon_provenance: dict[str, object] | None
    salmon_transcriptome: LocalReferenceAsset | None = None
    salmon_index_metadata: dict[str, object] | None = None
    salmon_validation: dict[str, object] | None = None
    hisat2_index: Path | None = None
    hisat2_splice_sites: LocalReferenceAsset | None = None
    hisat2_status: str = HISAT2_NOT_BUILT
    hisat2_provenance: dict[str, object] | None = None

    @property
    def salmon_strategy(self) -> str:
        if self.salmon_index is None:
            return "no pre-built Salmon index; prepare it before production execution"
        assert self.salmon_strategy_type is not None
        return f"use the manifest-declared pre-built {self.salmon_strategy_type.replace('_', '-')} Salmon index"

    @property
    def transcriptome_strategy(self) -> str:
        if self.salmon_transcriptome is not None:
            return "gtf_derived_exact_transcript_id_contract"
        return TRANSCRIPTOME_STRATEGY

    @property
    def external_transcript_fasta_used(self) -> bool:
        """The downloaded transcript FASTA is provenance only, never a runtime input."""

        return False

    def assets(self) -> tuple[LocalReferenceAsset, ...]:
        return (self.genome_fasta, self.annotation_gtf, self.transcript_fasta)

    def hisat2_arguments(self) -> list[tuple[str, Path]]:
        if self.hisat2_index is None:
            raise LocalReferenceError(
                f"Local reference HISAT2 index is not built. Run: rnaseq reference prepare-hisat2 {self.root}"
            )
        return [("--fasta", self.genome_fasta.path), ("--gtf", self.annotation_gtf.path), ("--hisat2_index", self.hisat2_index)]

    def nfcore_arguments(self) -> list[tuple[str, Path]]:
        if self.salmon_index is None:
            raise LocalReferenceError(
                f"Local reference Salmon index is not built. Run: rnaseq reference prepare {self.root}"
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
            "species": self.species,
            "provider": self.provider,
            "release": self.release,
            "assembly": self.assembly,
            "assembly_patch": self.assembly_patch,
            "manifest": {
                "path": str(self.manifest_path),
                "sha256": self.manifest_sha256,
                "schema_version": LOCAL_REFERENCE_MANIFEST_VERSION,
            },
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
                },
                "external_transcript_fasta_used": False,
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
                "splice_sites": ({"path": str(self.hisat2_splice_sites.path), "sha256": self.hisat2_splice_sites.sha256} if self.hisat2_splice_sites else None),
                "provenance": self.hisat2_provenance,
            },
        }


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
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


def validate_hisat2_index(index: Path) -> None:
    """Reject incomplete HISAT2 indexes before they can become executable."""

    if not index.is_dir():
        raise LocalReferenceError(f"Local reference hisat2.index directory not found: {index}")
    # HISAT2 supports both small (.ht2) and large (.ht2l) indexes.  A valid
    # index has one complete numbered family, never a partial mixture.
    small = [index / f"genome{suffix}" for suffix in REQUIRED_HISAT2_INDEX_SUFFIXES]
    large = [index / f"genome.{number}.ht2l" for number in range(1, 9)]
    family = small if all(path.is_file() for path in small) else large
    if not all(path.is_file() and path.stat().st_size > 0 for path in family):
        raise LocalReferenceError("Local reference HISAT2 index is incomplete; expected genome.1..8.ht2 or .ht2l files.")


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
    if manifest.get("schema_version") != LOCAL_REFERENCE_MANIFEST_VERSION:
        raise LocalReferenceError(
            "Unsupported local reference manifest schema version: "
            f"{manifest.get('schema_version')!r}; supported version: {LOCAL_REFERENCE_MANIFEST_VERSION}."
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
    assembly_patch = _require_string(identity.get("assembly_patch"), "reference.assembly_patch")
    files = _require_mapping(manifest.get("files"), "files")
    genome_fasta = _asset(root, files, "genome_fasta")
    annotation_gtf = _asset(root, files, "annotation_gtf")
    transcript_fasta = _asset(root, files, "transcript_fasta")
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
        if strategy == SALMON_STRATEGY_TRANSCRIPTOME_ONLY:
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
    hisat2_splice_sites: LocalReferenceAsset | None = None
    hisat2_provenance: dict[str, object] | None = None
    if hisat2_status == HISAT2_BUILT:
        hisat2_index = _resolve_under(root, _require_string(hisat2.get("index"), "hisat2.index"), "hisat2.index", directory=True)
        validate_hisat2_index(hisat2_index)
        hisat2_splice_sites = _asset_from_payload(root, _require_mapping(hisat2.get("splice_sites"), "hisat2.splice_sites"), "hisat2.splice_sites", "hisat2_splice_sites")
        hisat2_provenance = _require_mapping(hisat2.get("provenance"), "hisat2.provenance")
        for key, expected in (("hisat2_version", HISAT2_VERSION), ("genome_fasta_sha256", genome_fasta.sha256), ("annotation_gtf_sha256", annotation_gtf.sha256)):
            _require_equal(hisat2_provenance.get(key), expected, f"hisat2.provenance.{key}")
    elif hisat2_status != HISAT2_NOT_BUILT:
        raise LocalReferenceError("Local reference hisat2.status must be 'not_built' or 'built'.")
    return LocalReference(
        root=root,
        manifest_path=manifest_path,
        manifest_sha256=sha256_file(manifest_path),
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
        hisat2_splice_sites=hisat2_splice_sites,
        hisat2_status=hisat2_status,
        hisat2_provenance=hisat2_provenance,
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


def _container_path(reference: LocalReference, path: Path) -> str:
    return "/reference/" + path.resolve().relative_to(reference.root).as_posix()


def _run_reference_container(command: list[str], runner: Callable[..., subprocess.CompletedProcess[str]]) -> None:
    try:
        completed = runner(command, check=False, text=True, capture_output=True)
    except OSError as exc:
        raise ReferencePreparationError(f"Unable to start Docker for local reference preparation: {exc}") from exc
    if completed.returncode != 0:
        detail = (completed.stderr or completed.stdout or "no diagnostic output").strip()
        raise ReferencePreparationError(
            f"Reference preparation command failed (exit {completed.returncode}): {detail}"
        )


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
    transcriptome: str,
    strategy: str,
    validation_artifact: str,
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
    if strategy != SALMON_STRATEGY_TRANSCRIPTOME_ONLY:
        raise ReferenceAdoptionError(
            "Adopting existing indexes currently supports only strategy 'transcriptome_only'. "
            "Use 'rnaseq reference prepare' for decoy-aware construction."
        )
    if reference.salmon_status != SALMON_NOT_BUILT:
        raise ReferenceAdoptionError(
            "Refusing to replace an existing local reference Salmon declaration; "
            f"current status is {reference.salmon_status!r}."
        )
    try:
        index_path = _resolve_under(root, index, "adoption index", directory=True)
        transcriptome_path = _resolve_under(root, transcriptome, "adoption transcriptome")
        artifact_path = _resolve_under(root, validation_artifact, "adoption validation artifact")
    except LocalReferenceError as exc:
        raise ReferenceAdoptionError(str(exc)) from exc
    try:
        index_metadata = _salmon_index_metadata(index_path)
    except LocalReferenceError as exc:
        raise ReferenceAdoptionError(str(exc)) from exc
    transcriptome_sha256 = sha256_file(transcriptome_path)
    candidate = deepcopy(manifest)
    candidate["salmon"] = {
        "status": SALMON_BUILT,
        "strategy": SALMON_STRATEGY_TRANSCRIPTOME_ONLY,
        "index": index_path.relative_to(root).as_posix(),
        "transcriptome": {
            "path": transcriptome_path.relative_to(root).as_posix(),
            "sha256": transcriptome_sha256,
            "source": "gtf_derived",
        },
        "version": index_metadata["salmon_version"],
        "kmer_size": index_metadata["kmer_size"],
        "num_decoys": index_metadata["num_decoys"],
        "decoy_aware": False,
        "index_metadata": {
            key: index_metadata[key]
            for key in ("seq_hash", "name_hash", "info_index_version", "salmon_index_version", "seq_length", "num_kmers", "num_contigs")
        },
        "validation": {
            "artifact": artifact_path.relative_to(root).as_posix(),
            "transcript_id_contract": "exact",
            "fasta_only_transcripts": 0,
            "zero_gene_mappings": 0,
            "multi_gene_mappings": 0,
        },
        "provenance": {
            "mode": ADOPTED_EXISTING_INDEX,
            "adopted_at": now().astimezone(UTC).replace(microsecond=0).isoformat(),
        },
    }
    try:
        _validated_transcriptome_only_salmon(
            root,
            _require_mapping(candidate["salmon"], "salmon"),
            reference.genome_fasta,
            reference.annotation_gtf,
            identity={
                "species": reference.species,
                "provider": reference.provider,
                "release": reference.release,
                "assembly": reference.assembly,
            },
        )
    except LocalReferenceError as exc:
        raise ReferenceAdoptionError(str(exc)) from exc
    _write_manifest_atomically(reference.manifest_path, candidate)
    return load_local_reference_root(root)


def prepare_local_reference(
    reference_root: Path | str,
    *,
    threads: int = 4,
    runner: Callable[..., subprocess.CompletedProcess[str]] = subprocess.run,
    now: Callable[[], datetime] = lambda: datetime.now(UTC),
) -> LocalReference:
    """Build and atomically register a decoy-aware Salmon index once per reference."""

    if threads < 1:
        raise ReferencePreparationError("Reference preparation threads must be at least 1.")
    reference, manifest = _load_local_reference_root(
        Path(reference_root).expanduser().resolve(), "reference_manifest.yaml", None
    )
    if reference.salmon_status == SALMON_BUILT:
        raise ReferencePreparationError(
            f"Local reference Salmon index is already built: {reference.salmon_index}."
        )
    if reference.salmon_status != SALMON_NOT_BUILT:
        raise ReferencePreparationError(f"Unsupported local reference Salmon status: {reference.salmon_status}.")

    salmon_root = reference.root / "salmon"
    final_index = salmon_root / "index"
    if final_index.exists():
        raise ReferencePreparationError(
            f"Refusing to overwrite existing local reference index directory: {final_index}."
        )
    temporary_root = reference.root / f".rnaseq-reference-prepare-{uuid.uuid4().hex}"
    transcript_dir = temporary_root / "transcriptome"
    salmon_dir = temporary_root / "salmon"
    temporary_index = salmon_dir / "index"
    temporary_root.mkdir()
    transcript_dir.mkdir()
    salmon_dir.mkdir()
    mount = f"type=bind,src={reference.root},dst=/reference"
    transcript_container_dir = _container_path(reference, transcript_dir)
    genome = _container_path(reference, reference.genome_fasta.path)
    gtf = _container_path(reference, reference.annotation_gtf.path)
    transcript_prefix = f"{transcript_container_dir}/genome"
    transcript_fasta = f"{transcript_prefix}.transcripts.fa"
    salmon_container_dir = _container_path(reference, salmon_dir)
    try:
        _run_reference_container(
            [
                "docker", "run", "--rm", "--mount", mount, "-w", transcript_container_dir,
                RSEM_IMAGE, "rsem-prepare-reference", "--gtf", gtf, "--num-threads", str(threads),
                genome, transcript_prefix,
            ],
            runner,
        )
        script = (
            f"grep '^>' {genome} | cut -d ' ' -f 1 | cut -d $'\\t' -f 1 | sed 's/>//g' > decoys.txt && "
            f"cat {transcript_fasta} {genome} > gentrome.fa && "
            f"salmon index --threads {threads} -t gentrome.fa -d decoys.txt -i index -k {SALMON_KMER_SIZE}"
        )
        _run_reference_container(
            [
                "docker", "run", "--rm", "--mount", mount, "-w", salmon_container_dir,
                SALMON_IMAGE, "sh", "-c", script,
            ],
            runner,
        )
        validate_salmon_index(temporary_index)
        salmon_root.mkdir(exist_ok=True)
        temporary_index.replace(final_index)
        manifest["salmon"] = {
            "index": "salmon/index",
            "status": SALMON_BUILT,
            "strategy": SALMON_STRATEGY_DECOY_AWARE,
            "provenance": {
                "salmon_version": SALMON_VERSION,
                "index_path": "salmon/index",
                "genome_fasta_sha256": reference.genome_fasta.sha256,
                "annotation_gtf_sha256": reference.annotation_gtf.sha256,
                "transcriptome_strategy": TRANSCRIPTOME_STRATEGY,
                "decoy_strategy": DECOY_STRATEGY,
                "kmer_size": SALMON_KMER_SIZE,
                "built_at": now().astimezone(UTC).replace(microsecond=0).isoformat(),
            },
        }
        _write_manifest_atomically(reference.manifest_path, manifest)
    except Exception:
        if final_index.exists() and reference.salmon_status == SALMON_NOT_BUILT:
            shutil.rmtree(final_index)
        raise
    finally:
        shutil.rmtree(temporary_root, ignore_errors=True)
    return load_local_reference_root(reference.root)


def prepare_local_hisat2_reference(
    reference_root: Path | str,
    *,
    threads: int = 4,
    runner: Callable[..., subprocess.CompletedProcess[str]] = subprocess.run,
    now: Callable[[], datetime] = lambda: datetime.now(UTC),
) -> LocalReference:
    """Build and atomically register the HISAT2 genome index for one reference.

    A complete index is moved into place only after every numbered file and
    generated splice-site file has been checked.  The final-directory refusal
    is also the concurrency guard: a competing completed preparation is never
    overwritten.
    """

    if threads < 1:
        raise ReferencePreparationError("Reference preparation threads must be at least 1.")
    reference, manifest = _load_local_reference_root(Path(reference_root).expanduser().resolve(), "reference_manifest.yaml", None)
    if reference.hisat2_status == HISAT2_BUILT:
        raise ReferencePreparationError(f"Local reference HISAT2 index is already built: {reference.hisat2_index}.")
    final_root = reference.root / "hisat2"
    final_index = final_root / "index"
    if final_index.exists():
        raise ReferencePreparationError(f"Refusing to overwrite existing local HISAT2 index directory: {final_index}.")
    temporary = reference.root / f".rnaseq-hisat2-prepare-{uuid.uuid4().hex}"
    temporary_index = temporary / "index"
    temporary.mkdir()
    # hisat2-build creates the numbered files but not their parent directory.
    # Create this private staging directory before entering the container so a
    # successful build can still be atomically moved into the final location.
    temporary_index.mkdir()
    try:
        genome = _container_path(reference, reference.genome_fasta.path)
        gtf = _container_path(reference, reference.annotation_gtf.path)
        mount = f"type=bind,src={reference.root},dst=/reference"
        work = "/reference/" + temporary.relative_to(reference.root).as_posix()
        script = (
            f"hisat2_extract_splice_sites.py {gtf} > splice_sites.txt && "
            f"hisat2-build --threads {threads} --ss splice_sites.txt {genome} index/genome"
        )
        _run_reference_container(["docker", "run", "--rm", "--mount", mount, "-w", work, HISAT2_IMAGE, "sh", "-ec", script], runner)
        validate_hisat2_index(temporary_index)
        splice = temporary / "splice_sites.txt"
        # A valid single-exon annotation has no junctions.  HISAT2 accepts an
        # empty --ss file, so preserve it as the explicit, reproducible result
        # rather than rejecting an otherwise valid reference.
        if not splice.is_file():
            raise ReferencePreparationError("HISAT2 splice-site extraction produced no artifact.")
        final_root.mkdir(exist_ok=True)
        temporary_index.replace(final_index)
        final_splice = final_root / "splice_sites.txt"
        splice.replace(final_splice)
        manifest["hisat2"] = {
            "status": HISAT2_BUILT,
            "index": "hisat2/index",
            "splice_sites": {"path": "hisat2/splice_sites.txt", "sha256": sha256_file(final_splice)},
            "provenance": {
                "hisat2_version": HISAT2_VERSION,
                "container": HISAT2_IMAGE,
                "genome_fasta_sha256": reference.genome_fasta.sha256,
                "annotation_gtf_sha256": reference.annotation_gtf.sha256,
                "threads": threads,
                "splice_site_command": "hisat2_extract_splice_sites.py GTF",
                "index_command": "hisat2-build --ss splice_sites.txt FASTA index/genome",
                "built_at": now().astimezone(UTC).replace(microsecond=0).isoformat(),
            },
        }
        _write_manifest_atomically(reference.manifest_path, manifest)
    except Exception:
        if final_index.exists() and reference.hisat2_status == HISAT2_NOT_BUILT:
            shutil.rmtree(final_index)
        raise
    finally:
        shutil.rmtree(temporary, ignore_errors=True)
    return load_local_reference_root(reference.root)
