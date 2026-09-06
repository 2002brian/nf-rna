#!/usr/bin/env python3
"""Generate the deterministic paired-end HISAT2 acceptance fixture.

This creates intentionally tiny synthetic Mus musculus assets.  It is a
software-contract fixture, not a biological reference.  Run this script, then
prepare its managed HISAT2 index with the public CLI:

  python -c 'from rnaseq.cli import app; app()' reference prepare-hisat2 \\
    tests/fixtures/hisat2_paired_raw/reference
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import random
from pathlib import Path


ROOT = Path(__file__).resolve().parent
ADAPTER_R1 = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"
ADAPTER_R2 = "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"
READ_LENGTH = 75


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _reverse_complement(sequence: str) -> str:
    return sequence.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def _contig(seed: int, length: int = 520) -> str:
    rng = random.Random(seed)
    return "".join(rng.choice("ACGT") for _ in range(length))


def _pair_for_gene(
    contig: str, gene_strand: str, fragment_start: int, name: str, *, adapter: bool, fragment_length: int
) -> tuple[tuple[str, str], tuple[str, str]]:
    """Return an inward-facing FR pair from one defined genomic fragment.

    Coordinates here are zero-based.  A forward-stranded RNA library has R1
    on the transcript strand: R1 maps plus for GeneA and minus for GeneB.
    """
    if fragment_length < READ_LENGTH:
        fragment = contig[fragment_start : fragment_start + fragment_length]
        left, right = fragment, fragment
    else:
        left = contig[fragment_start : fragment_start + READ_LENGTH]
        right = contig[fragment_start + fragment_length - READ_LENGTH : fragment_start + fragment_length]
    if gene_strand == "+":
        r1, r2 = left, _reverse_complement(right)
    else:
        r1, r2 = _reverse_complement(right), left
    if adapter:
        r1 += ADAPTER_R1
        r2 += ADAPTER_R2
    return (f"@{name}/1", r1), (f"@{name}/2", r2)


def _write_fastq(path: Path, reads: list[tuple[str, str]]) -> None:
    with gzip.open(path, "wt", encoding="ascii", newline="\n") as handle:
        for name, sequence in reads:
            handle.write(f"{name}\n{sequence}\n+\n{'I' * len(sequence)}\n")


def _write_reference(root: Path, genome_a: str, genome_b: str, *, assembly_patch: str) -> None:
    reference = root / "reference"
    reference.mkdir(parents=True)
    genome = reference / "genome.fa"
    genome.write_text(f">chrFixtureA\n{genome_a}\n>chrFixtureB\n{genome_b}\n", encoding="ascii")
    gtf = reference / "genes.gtf"
    gtf.write_text(
        "chrFixtureA\tfixture\texon\t51\t470\t.\t+\t.\tgene_id \"GeneA\"; transcript_id \"TxA\";\n"
        "chrFixtureB\tfixture\texon\t51\t470\t.\t-\t.\tgene_id \"GeneB\"; transcript_id \"TxB\";\n",
        encoding="ascii",
    )
    transcript = reference / "transcripts.fa"
    transcript.write_text(f">TxA\n{genome_a[50:470]}\n>TxB\n{_reverse_complement(genome_b[50:470])}\n", encoding="ascii")
    manifest = reference / "reference_manifest.yaml"
    manifest.write_text(
        "schema_version: '1.1'\n"
        "purpose: synthetic_test\n"
        "reference:\n"
        "  species: Mus musculus\n"
        "  provider: synthetic-paired-uat\n"
        "  release: 1\n"
        "  assembly: paired-fixture\n"
        f"  assembly_patch: {assembly_patch}\n"
        "files:\n"
        "  genome_fasta:\n"
        "    path: genome.fa\n"
        f"    sha256: {_sha256(genome)}\n"
        "  annotation_gtf:\n"
        "    path: genes.gtf\n"
        f"    sha256: {_sha256(gtf)}\n"
        "  transcript_fasta:\n"
        "    path: transcripts.fa\n"
        f"    sha256: {_sha256(transcript)}\n"
        "salmon:\n"
        "  status: not_built\n"
        "hisat2:\n"
        "  status: not_built\n",
        encoding="utf-8",
    )


def generate(root: Path) -> None:
    existing = [path for path in root.iterdir() if path.resolve() != Path(__file__).resolve()] if root.exists() else []
    if existing:
        raise SystemExit(f"Refusing to overwrite non-empty fixture root: {root}")
    root.mkdir(parents=True, exist_ok=True)
    genome_a, genome_b = _contig(7331), _contig(7332)
    _write_reference(root, genome_a, genome_b, assembly_patch="v2" if root.name.endswith("_v2") else "v1")
    fastq = root / "fastq"
    fastq.mkdir()
    lanes = {
        "PairAlpha_L001": [
            ("GeneA", "+", 110, "PairAlphaA1", False, 180),
            ("GeneA", "+", 130, "PairAlphaA2_adapter", True, 60),
        ],
        "PairAlpha_L002": [
            ("GeneA", "+", 150, "PairAlphaA3", False, 180),
            ("GeneB", "-", 110, "PairAlphaB1_adapter", True, 60),
        ],
        "PairBeta_L001": [
            ("GeneB", "-", 130, "PairBetaB1", False, 180),
            ("GeneB", "-", 150, "PairBetaB2_adapter", True, 60),
        ],
    }
    sample_sheet = ["sample,fastq_1,fastq_2,strandedness"]
    for lane, records in lanes.items():
        r1_reads, r2_reads = [], []
        for gene, strand, start, name, adapter, fragment_length in records:
            r1, r2 = _pair_for_gene(
                genome_a if gene == "GeneA" else genome_b, strand, start, name,
                adapter=adapter, fragment_length=fragment_length,
            )
            r1_reads.append(r1)
            r2_reads.append(r2)
        r1_path, r2_path = fastq / f"{lane}_R1.fastq.gz", fastq / f"{lane}_R2.fastq.gz"
        _write_fastq(r1_path, r1_reads)
        _write_fastq(r2_path, r2_reads)
        sample = "PairAlpha" if lane.startswith("PairAlpha") else "PairBeta"
        sample_sheet.append(f"{sample},{r1_path.relative_to(root)},{r2_path.relative_to(root)},forward")
    (root / "samplesheet.csv").write_text("\n".join(sample_sheet) + "\n", encoding="utf-8")
    (root / "expected_counts.csv").write_text(
        "gene_id,PairAlpha,PairBeta\nGeneA,3,0\nGeneB,1,2\n", encoding="utf-8"
    )
    (root / "design.md").write_text(
        "# Paired-end HISAT2 UAT fixture\n\n"
        "The library is forward stranded: R1 is on the transcript strand, so HISAT2 uses `--rna-strandness FR` and featureCounts uses `-s 1`.\n\n"
        "`PairAlpha` has two technical lanes from one library: lane 001 has two GeneA fragments and lane 002 has one GeneA plus one GeneB fragment. `PairBeta` is a separate sample with two GeneB fragments. All six fragments are uniquely mappable. Three ordinary fragments have 180 bp inserts; three adapter-bearing fragments have 60 bp inserts, allowing paired-end overlap detection to trim standard Illumina adapters while retaining 60 high-quality biological bases.\n\n"
        "Expected counts are independently specified by this construction, before execution: PairAlpha = GeneA 3 / GeneB 1; PairBeta = GeneA 0 / GeneB 2. A fragment contributes once under `-p --countReadPairs -B -C`; BAM alignment-record counts are therefore twice the eligible fragment count.\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=ROOT, help="empty fixture root to populate")
    arguments = parser.parse_args()
    generate(arguments.root.resolve())
