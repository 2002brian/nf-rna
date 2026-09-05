"""Optional real-tool checks for featureCounts semantics.

These run only when HISAT2/SAMtools/featureCounts are on PATH.  They are kept
separate from the pinned-container smoke test: a host environment can provide
useful behavioural evidence but cannot substitute for recorded container
versions and digests.
"""

from __future__ import annotations

import csv
import shutil
import subprocess
from pathlib import Path

import pytest


TOOLS = {name: shutil.which(name) for name in ("samtools", "featureCounts")}
pytestmark = pytest.mark.skipif(not all(TOOLS.values()), reason="real SAMtools and featureCounts are required")


def _run(arguments: list[str], cwd: Path) -> None:
    result = subprocess.run(arguments, cwd=cwd, capture_output=True, text=True, check=False)
    assert result.returncode == 0, result.stderr


def _bam(tmp_path: Path, name: str, records: list[str]) -> Path:
    sam = tmp_path / f"{name}.sam"
    sam.write_text("@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:1\tLN:1000\n@SQ\tSN:2\tLN:1000\n" + "\n".join(records) + "\n", encoding="utf-8")
    unsorted, bam = tmp_path / f"{name}.unsorted.bam", tmp_path / f"{name}.bam"
    _run([str(TOOLS["samtools"]), "view", "-bS", "-o", str(unsorted), str(sam)], tmp_path)
    _run([str(TOOLS["samtools"]), "sort", "-o", str(bam), str(unsorted)], tmp_path)
    _run([str(TOOLS["samtools"]), "index", str(bam)], tmp_path)
    return bam


def _counts(tmp_path: Path, bam: Path, name: str, arguments: list[str]) -> dict[str, int]:
    output = tmp_path / f"{name}.counts.txt"
    _run([str(TOOLS["featureCounts"]), "-a", str(tmp_path / "genes.gtf"), "-o", str(output), *arguments, str(bam)], tmp_path)
    with output.open(encoding="utf-8", newline="") as handle:
        rows = [row for row in csv.reader((line for line in handle if not line.startswith("#")), delimiter="\t")]
    return {row[0]: int(row[6]) for row in rows[1:]}


def _countable_bam(tmp_path: Path, bam: Path) -> Path:
    """Mirror the workflow's non-destructive secondary/supplementary filter."""

    output = tmp_path / f"{bam.stem}.countable.bam"
    _run([str(TOOLS["samtools"]), "view", "-bh", "-F", "0x900", "-o", str(output), str(bam)], tmp_path)
    _run([str(TOOLS["samtools"]), "index", str(output)], tmp_path)
    return output


def test_real_featurecounts_counting_policy(tmp_path: Path):
    """Expected counts derive from the hand-written alignment records below."""

    (tmp_path / "genes.gtf").write_text(
        "1\tfixture\texon\t100\t200\t.\t+\t.\tgene_id \"GeneA\";\n"
        "1\tfixture\texon\t150\t250\t.\t+\t.\tgene_id \"GeneOverlap\";\n"
        "1\tfixture\texon\t300\t400\t.\t+\t.\tgene_id \"GeneB\";\n",
        encoding="utf-8",
    )
    seq, quality = "A" * 20, "I" * 20
    # unique A; ambiguous A/overlap; NH=2 multimapper; secondary and
    # supplementary A.  Only the first record is eligible by default.
    single = _bam(tmp_path, "single", [
        f"unique\t0\t1\t110\t60\t20M\t*\t0\t0\t{seq}\t{quality}\tNH:i:1",
        f"overlap\t0\t1\t160\t60\t20M\t*\t0\t0\t{seq}\t{quality}\tNH:i:1",
        f"multi\t0\t1\t310\t60\t20M\t*\t0\t0\t{seq}\t{quality}\tNH:i:2",
        f"secondary\t256\t1\t110\t60\t20M\t*\t0\t0\t{seq}\t{quality}\tNH:i:1",
        f"supplementary\t2048\t1\t110\t60\t20M\t*\t0\t0\t{seq}\t{quality}\tNH:i:1",
    ])
    observed = _counts(tmp_path, _countable_bam(tmp_path, single), "single", ["-t", "exon", "-g", "gene_id", "-s", "0", "-Q", "0", "--primary"])
    assert observed == {"GeneA": 1, "GeneOverlap": 0, "GeneB": 0}

    strands = _bam(tmp_path, "strands", [
        f"forward\t0\t1\t110\t60\t20M\t*\t0\t0\t{seq}\t{quality}\tNH:i:1",
        f"reverse\t16\t1\t310\t60\t20M\t*\t0\t0\t{seq}\t{quality}\tNH:i:1",
    ])
    countable_strands = _countable_bam(tmp_path, strands)
    assert _counts(tmp_path, countable_strands, "forward", ["-t", "exon", "-g", "gene_id", "-s", "1", "-Q", "0", "--primary"]) == {"GeneA": 1, "GeneOverlap": 0, "GeneB": 0}
    assert _counts(tmp_path, countable_strands, "reverse", ["-t", "exon", "-g", "gene_id", "-s", "2", "-Q", "0", "--primary"]) == {"GeneA": 0, "GeneOverlap": 0, "GeneB": 1}

    # One proper pair, one R1-only fragment, and one inter-chromosomal pair.
    paired = _bam(tmp_path, "paired", [
        f"pair\t99\t1\t110\t60\t20M\t=\t150\t60\t{seq}\t{quality}\tNH:i:1",
        f"pair\t147\t1\t150\t60\t20M\t=\t110\t-60\t{seq}\t{quality}\tNH:i:1",
        f"orphan\t73\t1\t110\t60\t20M\t*\t0\t0\t{seq}\t{quality}\tNH:i:1",
        f"chimera\t97\t1\t110\t60\t20M\t2\t110\t0\t{seq}\t{quality}\tNH:i:1",
        f"chimera\t145\t2\t110\t60\t20M\t1\t110\t0\t{seq}\t{quality}\tNH:i:1",
    ])
    paired_counts = _counts(tmp_path, _countable_bam(tmp_path, paired), "paired", ["-t", "exon", "-g", "gene_id", "-s", "0", "-Q", "0", "--primary", "-p", "--countReadPairs", "-B", "-C"])
    assert paired_counts == {"GeneA": 1, "GeneOverlap": 0, "GeneB": 0}

    lane_a = _bam(tmp_path, "lane_a", [f"lanea\t0\t1\t110\t60\t20M\t*\t0\t0\t{seq}\t{quality}\tNH:i:1"])
    lane_b = _bam(tmp_path, "lane_b", [f"laneb\t0\t1\t110\t60\t20M\t*\t0\t0\t{seq}\t{quality}\tNH:i:1"])
    merged = tmp_path / "sample.bam"
    _run([str(TOOLS["samtools"]), "merge", "-o", str(merged), str(lane_a), str(lane_b)], tmp_path)
    _run([str(TOOLS["samtools"]), "index", str(merged)], tmp_path)
    assert _counts(tmp_path, _countable_bam(tmp_path, merged), "merged", ["-t", "exon", "-g", "gene_id", "-s", "0", "-Q", "0", "--primary"])["GeneA"] == 2
