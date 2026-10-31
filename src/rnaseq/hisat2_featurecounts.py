"""Stable contracts shared by the HISAT2/featureCounts workflow and control plane.

This module intentionally does not infer strandedness or sample identity.  The
workflow supplies an explicit library declaration and this module turns one
featureCounts file per *declared sample* into the canonical matrix consumed by
the existing count-based DESeq2 path.  It also judges whether the declaration
is consistent with the reads (strandedness QC), which is reported but never
changes the declaration or the counts.
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


# Strandedness consistency QC (Phase 3 finding F-1).  The evidence is the same
# count-only BAM counted again with featureCounts -s 1 and -s 2 (all other
# options unchanged), so the ratio compares two strand hypotheses on the same
# reads and does not depend on species, depth or mapping rate.  Annotated
# antisense overlaps move the ratio towards 1 (a perfectly reverse simulated
# yeast library gave forward/reverse 0.075-0.096; a real one 0.17); they cannot
# make the opposite direction dominate, so a contradiction needs >= 3x.  The
# check only reports: the declared strandedness is never changed and nothing is
# recounted.
STRAND_QC_MIN_ASSIGNED = 1000
STRAND_QC_DOMINANT_RATIO = 3.0
STRAND_QC_BALANCED_RATIO = 1.5
STRAND_QC_POLICY = {
    "evidence": "assigned fragments (reads for single-end) of the count-only BAM under featureCounts -s 1 and -s 2, other counting options unchanged",
    "ratio": "assigned(-s 1) / assigned(-s 2)",
    "minimum_strand_assigned": STRAND_QC_MIN_ASSIGNED,
    "stranded_if_ratio_at_least": STRAND_QC_DOMINANT_RATIO,
    "unstranded_if_ratio_within": [round(1 / STRAND_QC_BALANCED_RATIO, 4), STRAND_QC_BALANCED_RATIO],
    "action": "report only: the declared strandedness is used unchanged and nothing is recounted",
}
STRAND_QC_SEVERITY = {
    "CONSISTENT": "pass",
    "CONTRADICTS_DECLARATION": "warning",
    "APPEARS_UNSTRANDED": "warning",
    "APPEARS_STRANDED": "advisory",
    "INCONCLUSIVE": "inconclusive",
}


def featurecounts_assigned(summary: Path) -> int:
    """The ``Assigned`` total of one featureCounts ``.summary`` file."""

    with summary.open(encoding="utf-8") as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if fields[0] == "Assigned" and len(fields) == 2:
                return int(fields[1])
    raise ValueError(f"featureCounts summary has no Assigned row: {summary}")


def strandedness_verdict(declared: str, assigned_forward: int, assigned_reverse: int) -> dict[str, object]:
    """Compare the declared strandedness with the -s 1 / -s 2 assignment evidence of one sample."""

    if declared not in {"unstranded", "forward", "reverse"}:
        raise ValueError("strandedness QC requires declared strandedness: unstranded, forward, or reverse.")
    total = assigned_forward + assigned_reverse
    ratio = assigned_forward / assigned_reverse if assigned_reverse else (float("inf") if assigned_forward else None)
    if total < STRAND_QC_MIN_ASSIGNED or ratio is None:
        observed = "insufficient_evidence"
    elif ratio >= STRAND_QC_DOMINANT_RATIO:
        observed = "forward"
    elif ratio <= 1 / STRAND_QC_DOMINANT_RATIO:
        observed = "reverse"
    elif 1 / STRAND_QC_BALANCED_RATIO <= ratio <= STRAND_QC_BALANCED_RATIO:
        observed = "unstranded"
    else:
        observed = "ambiguous"
    if observed in {"insufficient_evidence", "ambiguous"}:
        status = "INCONCLUSIVE"
        message = (f"fewer than {STRAND_QC_MIN_ASSIGNED} strand-assignable fragments" if observed == "insufficient_evidence"
                   else "neither strand dominates clearly and the two are not balanced")
    elif observed == declared:
        status, message = "CONSISTENT", f"observed {observed}, as declared"
    elif declared == "unstranded":
        status = "APPEARS_STRANDED"
        message = (f"declared unstranded, but the data look {observed}-stranded; counting as unstranded is valid, "
                   "but check the library preparation protocol")
    elif observed == "unstranded":
        status = "APPEARS_UNSTRANDED"
        message = (f"declared {declared}, but the data look unstranded; strand-specific counting may discard or "
                   "misassign about half of the fragments")
    else:
        status = "CONTRADICTS_DECLARATION"
        message = (f"declared {declared}, but the data look {observed}-stranded; counts made with the declared "
                   "strandedness are likely wrong")
    return {
        "declared": declared,
        "assigned_forward_s1": assigned_forward,
        "assigned_reverse_s2": assigned_reverse,
        "forward_over_reverse": None if ratio is None else (ratio if ratio != float("inf") else "inf"),
        "observed": observed,
        "status": status,
        "severity": STRAND_QC_SEVERITY[status],
        "message": message,
    }


def evaluate_strandedness(declared: str, evidence: dict[str, tuple[Path, Path]]) -> dict[str, object]:
    """Run-level strandedness QC from per-sample (-s 1, -s 2) featureCounts summaries."""

    samples = {sample: strandedness_verdict(declared, featurecounts_assigned(s1), featurecounts_assigned(s2))
               for sample, (s1, s2) in sorted(evidence.items())}
    severities = {item["severity"] for item in samples.values()}
    status = next((name for name, severity in (("WARNING", "warning"), ("ADVISORY", "advisory"), ("INCONCLUSIVE", "inconclusive"))
                   if severity in severities), "CONSISTENT")
    return {"status": status, "declared": declared, "policy": STRAND_QC_POLICY, "samples": samples}


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
