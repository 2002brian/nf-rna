"""Choices and defaults shared by the interactive ``rnaseq new`` wizard and the GUI.

Both front ends offer exactly these canonical values with these descriptions
and defaults, and both turn a reference choice into the same ``reference``
mapping.  Validation of the final combination stays in
:func:`rnaseq.project.create_project`.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from rnaseq.errors import ProjectCreationError
from rnaseq.execution import LocalResourceCapacity, detect_local_resource_capacity
from rnaseq.models import Species
from rnaseq.references import LocalReference
from rnaseq.resource_policy import AUTO, GIB, DetectedResources, resolve_policy

Choice = tuple[str, str]

SPECIES_CHOICES: list[Choice] = [("human", "Homo sapiens"), ("mouse", "Mus musculus")]
INPUT_TYPE_CHOICES: list[Choice] = [
    ("fastq", "Reads requiring upstream processing"),
    ("raw_counts", "Integer gene-level raw-count matrix"),
]
LAYOUT_CHOICES: list[Choice] = [("paired_end", "R1 and R2 per lane"), ("single_end", "R1 only")]
PREPROCESSING_CHOICES: list[Choice] = [
    ("raw", "Run adapter/quality trimming"),
    ("pretrimmed", "Keep supplied reads; skip trimming"),
]
METHOD_CHOICES: list[Choice] = [
    ("salmon", "nf-core/rnaseq pseudoalignment (default)"),
    ("hisat2_featurecounts", "HISAT2 alignment plus gene-level featureCounts"),
]
REFERENCE_SOURCE_CHOICES: list[Choice] = [
    ("igenomes", "managed iGenomes convenience route (Salmon only)"),
    ("local", "checksum-bound managed local reference"),
    ("custom", "copy supported custom assets into the project"),
]
DESIGN_CHOICES: list[Choice] = [
    ("two_group", "Unpaired two-level condition"),
    ("paired_two_group", "Biologically paired two-level condition"),
    ("multi_group", "At least three condition levels"),
]
GENE_ID_TYPE_CHOICES: list[Choice] = [
    ("ENSEMBL", "Ensembl gene IDs, e.g. ENSMUSG00000000001 (Ensembl GTF / managed references)"),
    ("ENTREZID", "NCBI Entrez gene IDs (numeric)"),
    ("SYMBOL", "Official gene symbols"),
]
DEFAULT_PAIR_ID = "patient"

_QC_CHOICE: Choice = ("qc", "Quantification + technical QC only")
_L1_CHOICE: Choice = ("L1", "Expression-level QC and exploratory analysis (differential expression is set up later with `rnaseq design`)")
_L2_CHOICE: Choice = ("L2", "L1 plus differential expression and functional enrichment (GO ORA, KEGG ORA, GO + KEGG GSEA)")


def strandedness_choices(method: str) -> list[Choice]:
    """Salmon may infer strandedness; HISAT2 + featureCounts needs an explicit library type."""

    choices = [("auto", "infer with Salmon"), ("unstranded", "no stranded protocol"), ("forward", "forward stranded"), ("reverse", "reverse stranded")]
    return choices[1:] if method == "hisat2_featurecounts" else choices


def default_strandedness(method: str) -> str:
    return "auto" if method == "salmon" else "unstranded"


def default_reference_source(method: str) -> str:
    return "igenomes" if method == "salmon" else "local"


def preset_choices(input_type: str) -> list[Choice]:
    """The qc scope exists only for FASTQ input."""

    return [_QC_CHOICE, _L1_CHOICE, _L2_CHOICE] if input_type == "fastq" else [_L1_CHOICE, _L2_CHOICE]


def species_value(choice: str) -> Species:
    return Species.HUMAN if choice == "human" else Species.MOUSE


@dataclass(frozen=True)
class ResourceSuggestion:
    logical_cpus: int | None
    memory_gib: int | None
    auto_cpus: int
    auto_memory_gib: int


def resource_suggestion(capacity: LocalResourceCapacity | None = None) -> ResourceSuggestion:
    """What ``auto`` resolves to on this machine (it is resolved again on each run's machine)."""

    capacity = capacity or detect_local_resource_capacity()
    preview = resolve_policy(AUTO, AUTO, DetectedResources(
        capacity.logical_cpus,
        capacity.total_memory_gib * GIB if capacity.total_memory_gib is not None else None,
    ))
    return ResourceSuggestion(
        capacity.logical_cpus, capacity.available_memory_gib or capacity.total_memory_gib,
        preview.requested_cpus, preview.requested_memory_gib,
    )


def managed_reference_settings(reference: LocalReference) -> tuple[Path, str]:
    """``(reference root, manifest relative to the root)`` of a registered reference."""

    return reference.root, str(reference.manifest_path.relative_to(reference.root))


def reference_options(
    *, source: str, species: Species, method: str, local_root: Path | None,
    local_manifest: str | None, fasta: Path | None, gtf: Path | None,
    transcript_fasta: Path | None, salmon_index: Path | None, hisat2_index: Path | None,
) -> dict[str, object]:
    """The ``reference`` mapping of project.yaml for one reference choice."""

    if source == "igenomes":
        if method != "salmon":
            raise ProjectCreationError("iGenomes is supported only by the Salmon backend; choose local or custom for HISAT2 + featureCounts.")
        return {"source": "igenomes", "genome": "GRCh38" if species is Species.HUMAN else "GRCm39"}
    if source == "local":
        if local_root is None:
            raise ProjectCreationError("Managed local reference requires --reference-root.")
        return {"source": "local", "root": str(local_root.expanduser().resolve()), "manifest": local_manifest or "reference_manifest.yaml"}
    if fasta is None or gtf is None:
        raise ProjectCreationError("Custom reference requires --reference-fasta and --reference-gtf.")
    reference: dict[str, object] = {"source": "custom", "fasta": str(fasta), "gtf": str(gtf)}
    if transcript_fasta is not None:
        reference["transcript_fasta"] = str(transcript_fasta)
    if salmon_index is not None:
        reference["salmon_index"] = str(salmon_index)
    if hisat2_index is not None:
        reference["hisat2_index"] = str(hisat2_index)
    return reference
