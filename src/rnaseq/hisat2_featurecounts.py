"""Stable contracts shared by the HISAT2/featureCounts workflow and control plane.

This module intentionally does not attempt to infer strandedness or sample
identity.  The workflow supplies an explicit library declaration and this
module turns one featureCounts file per *declared sample* into the canonical
matrix consumed by the existing count-based DESeq2 path.
"""

from __future__ import annotations

import argparse
import csv
from collections import OrderedDict
from pathlib import Path


HISAT2_VERSION = "2.2.3"
SAMTOOLS_VERSION = "1.21"
SUBREAD_VERSION = "2.0.6"
FASTQC_VERSION = "0.12.1"
FASTP_VERSION = "0.24.0"
MULTIQC_VERSION = "1.33"

# Image *references* are pinned to immutable Bioconda build tags.  The service
# records a digest only after Docker has actually resolved it; a source tree
# must not claim a digest it has not observed.
HISAT2_IMAGE = "quay.io/biocontainers/hisat2:2.2.3--h8471819_0"
SAMTOOLS_IMAGE = "quay.io/biocontainers/samtools:1.21--h50ea8bc_0"
SUBREAD_IMAGE = "quay.io/biocontainers/subread:2.0.6--he4a0461_2"
FASTQC_IMAGE = "quay.io/biocontainers/fastqc:0.12.1--hdfd78af_0"
FASTP_IMAGE = "quay.io/biocontainers/fastp:0.24.0--h125f33a_0"
MULTIQC_IMAGE = "community.wave.seqera.io/library/multiqc:1.33--ee7739d47738383b"

COUNTING_POLICY = {
    "feature_type": "exon",
    "grouping_attribute": "gene_id",
    "multimappers": "excluded (featureCounts default; no -M)",
    "multi_gene_overlap": "excluded (featureCounts default; no -O)",
    "fractional_counting": "disabled (no --fraction)",
    "mapq": 0,
    "secondary_supplementary": "excluded in a count-only BAM using samtools -F 0x900; original diagnostic BAM retained",
    "pcr_duplicates": "not removed",
    "paired_end": "fragment; require both mapped mates and exclude chimeras",
    "single_end": "read",
}


def hisat2_strand_option(strandedness: str, layout: str) -> str | None:
    """Translate the project contract into HISAT2's layout-specific syntax."""

    if strandedness == "unstranded":
        return None
    if strandedness not in {"forward", "reverse"}:
        raise ValueError("HISAT2 requires strandedness: unstranded, forward, or reverse.")
    if layout == "paired_end":
        return "FR" if strandedness == "forward" else "RF"
    if layout == "single_end":
        return "F" if strandedness == "forward" else "R"
    raise ValueError(f"Unsupported sequencing layout: {layout!r}.")


def featurecounts_arguments(*, layout: str, strandedness: str) -> list[str]:
    """Return the complete scientific counting policy, excluding input paths."""

    strand = {"unstranded": "0", "forward": "1", "reverse": "2"}.get(strandedness)
    if strand is None:
        raise ValueError("featureCounts requires strandedness: unstranded, forward, or reverse.")
    arguments = ["-t", "exon", "-g", "gene_id", "-s", strand, "-Q", "0", "--primary"]
    if layout == "paired_end":
        return [*arguments, "-p", "--countReadPairs", "-B", "-C"]
    if layout == "single_end":
        return arguments
    raise ValueError(f"Unsupported sequencing layout: {layout!r}.")


def _parse_featurecounts(path: Path) -> OrderedDict[str, int]:
    """Parse one pinned featureCounts output without guessing its BAM column."""

    rows: OrderedDict[str, int] = OrderedDict()
    header: list[str] | None = None
    with path.open(encoding="utf-8", newline="") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if header is None:
                header = fields
                if len(header) < 7 or header[0] != "Geneid":
                    raise ValueError(f"featureCounts output has no supported Geneid header: {path}")
                continue
            if len(fields) != len(header):
                raise ValueError(f"featureCounts output has an inconsistent row width: {path}")
            gene_id = fields[0].strip()
            if not gene_id or gene_id in rows:
                raise ValueError(f"featureCounts output has a blank or duplicate gene ID: {path}")
            # A per-sample invocation must produce exactly one count column;
            # keeping this strict prevents basename-order accidents.
            if len(fields[6:]) != 1:
                raise ValueError(f"featureCounts output must contain one declared BAM count column: {path}")
            try:
                count = int(fields[6])
            except ValueError as exc:
                raise ValueError(f"featureCounts count is not an integer for {gene_id!r}: {path}") from exc
            if count < 0:
                raise ValueError(f"featureCounts count is negative for {gene_id!r}: {path}")
            rows[gene_id] = count
    if header is None or not rows:
        raise ValueError(f"featureCounts output is empty: {path}")
    return rows


def assemble_count_matrix(sample_files: dict[str, Path], output: Path) -> None:
    """Assemble deterministic gene-by-sample CSV using explicit sample IDs."""

    if not sample_files or any(not sample for sample in sample_files):
        raise ValueError("featureCounts assembly requires nonblank explicit sample IDs.")
    parsed = {sample: _parse_featurecounts(path) for sample, path in sorted(sample_files.items())}
    first_sample = next(iter(parsed))
    genes = list(parsed[first_sample])
    expected = set(genes)
    for sample, values in parsed.items():
        if set(values) != expected:
            raise ValueError(f"featureCounts gene set for {sample!r} differs from {first_sample!r}.")
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        samples = sorted(parsed)
        writer.writerow(["gene_id", *samples])
        for gene in genes:
            writer.writerow([gene, *(parsed[sample][gene] for sample in samples)])


def _main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--sample", action="append", required=True, metavar="ID=PATH")
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    mapping: dict[str, Path] = {}
    for item in args.sample:
        sample, separator, configured = item.partition("=")
        if not separator or not sample or not configured or sample in mapping:
            parser.error("--sample must be a unique SAMPLE_ID=FEATURECOUNTS_PATH value")
        mapping[sample] = Path(configured)
    assemble_count_matrix(mapping, args.out)


if __name__ == "__main__":
    _main()
