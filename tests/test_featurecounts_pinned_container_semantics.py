"""Production-container acceptance checks for the featureCounts contract.

These fixtures deliberately run the same count-only BAM transformation and
featureCounts arguments as ``workflow/hisat2_featurecounts.nf``.  They are
skipped only when the Docker daemon or the already-pinned process images are
unavailable, so an explicit acceptance run proves the actual production tool
versions rather than substituting host executables.
"""

from __future__ import annotations

import csv
import shutil
import subprocess
from pathlib import Path

import pytest

from rnaseq.hisat2_featurecounts import SAMTOOLS_IMAGE, SUBREAD_IMAGE


DOCKER = shutil.which("docker")
pytestmark = pytest.mark.skipif(DOCKER is None, reason="Docker is required for pinned-container semantics")


def _run(image: str, arguments: list[str], cwd: Path) -> None:
    result = subprocess.run(
        [str(DOCKER), "run", "--rm", "-v", f"{cwd}:/work", "-w", "/work", image, *arguments],
        capture_output=True,
        text=True,
        check=False,
    )
    assert result.returncode == 0, result.stderr


def _samtools(arguments: list[str], cwd: Path) -> None:
    _run(SAMTOOLS_IMAGE, ["samtools", *arguments], cwd)


def _featurecounts(arguments: list[str], cwd: Path) -> None:
    _run(SUBREAD_IMAGE, ["featureCounts", *arguments], cwd)


def _bam(tmp_path: Path, name: str, records: list[str]) -> Path:
    sam = tmp_path / f"{name}.sam"
    sam.write_text("@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:1\tLN:1000\n@SQ\tSN:2\tLN:1000\n" + "\n".join(records) + "\n", encoding="utf-8")
    _samtools(["view", "-bS", "-o", f"{name}.unsorted.bam", f"{name}.sam"], tmp_path)
    _samtools(["sort", "-o", f"{name}.bam", f"{name}.unsorted.bam"], tmp_path)
    _samtools(["index", f"{name}.bam"], tmp_path)
    return tmp_path / f"{name}.bam"


def _countable_bam(tmp_path: Path, bam: Path) -> Path:
    """Make the separately retained production count input without losing NH."""

    output = tmp_path / f"{bam.stem}.countable.bam"
    _samtools(["view", "-bh", "-F", "0x900", "-o", output.name, bam.name], tmp_path)
    _samtools(["index", output.name], tmp_path)
    return output


def _counts(tmp_path: Path, bam: Path, name: str, arguments: list[str]) -> dict[str, int]:
    output = tmp_path / f"{name}.counts.txt"
    _featurecounts(["-a", "genes.gtf", "-o", output.name, *arguments, bam.name], tmp_path)
    with output.open(encoding="utf-8", newline="") as handle:
        rows = [row for row in csv.reader((line for line in handle if not line.startswith("#")), delimiter="\t")]
    return {row[0]: int(row[6]) for row in rows[1:]}


def test_pinned_container_featurecounts_counting_policy(tmp_path: Path):
    (tmp_path / "genes.gtf").write_text(
        "1\tfixture\texon\t100\t200\t.\t+\t.\tgene_id \"GeneA\";\n"
        "1\tfixture\texon\t150\t250\t.\t+\t.\tgene_id \"GeneOverlap\";\n"
        "1\tfixture\texon\t300\t400\t.\t+\t.\tgene_id \"GeneB\";\n",
        encoding="utf-8",
    )
    sequence, quality = "A" * 20, "I" * 20
    single = _bam(tmp_path, "single", [
        f"unique\t0\t1\t110\t60\t20M\t*\t0\t0\t{sequence}\t{quality}\tNH:i:1",
        f"overlap\t0\t1\t160\t60\t20M\t*\t0\t0\t{sequence}\t{quality}\tNH:i:1",
        # Its secondary alignment is removed by -F 0x900, but NH:i:2 remains
        # on this primary alignment and must keep the fragment excluded.
        f"multi\t0\t1\t310\t60\t20M\t*\t0\t0\t{sequence}\t{quality}\tNH:i:2",
        f"multi\t256\t2\t310\t60\t20M\t*\t0\t0\t{sequence}\t{quality}\tNH:i:2",
        f"secondary\t256\t1\t110\t60\t20M\t*\t0\t0\t{sequence}\t{quality}\tNH:i:1",
        f"supplementary\t2048\t1\t110\t60\t20M\t*\t0\t0\t{sequence}\t{quality}\tNH:i:1",
    ])
    countable = _countable_bam(tmp_path, single)
    # Prove the original and transformed artifacts retain NH as their
    # reproducible lineage; only flags 0x100/0x800 are removed.
    original = subprocess.run([str(DOCKER), "run", "--rm", "-v", f"{tmp_path}:/work", "-w", "/work", SAMTOOLS_IMAGE, "samtools", "view", single.name], capture_output=True, text=True, check=False)
    filtered = subprocess.run([str(DOCKER), "run", "--rm", "-v", f"{tmp_path}:/work", "-w", "/work", SAMTOOLS_IMAGE, "samtools", "view", countable.name], capture_output=True, text=True, check=False)
    assert original.returncode == filtered.returncode == 0
    assert "multi\t0" in filtered.stdout and "NH:i:2" in filtered.stdout
    assert "multi\t256" in original.stdout and "multi\t256" not in filtered.stdout
    assert _counts(tmp_path, countable, "single", ["-t", "exon", "-g", "gene_id", "-s", "0", "-Q", "0", "--primary"]) == {"GeneA": 1, "GeneOverlap": 0, "GeneB": 0}

    strands = _bam(tmp_path, "strands", [
        f"forward\t0\t1\t110\t60\t20M\t*\t0\t0\t{sequence}\t{quality}\tNH:i:1",
        f"reverse\t16\t1\t310\t60\t20M\t*\t0\t0\t{sequence}\t{quality}\tNH:i:1",
    ])
    strand_bam = _countable_bam(tmp_path, strands)
    assert _counts(tmp_path, strand_bam, "forward", ["-t", "exon", "-g", "gene_id", "-s", "1", "-Q", "0", "--primary"]) == {"GeneA": 1, "GeneOverlap": 0, "GeneB": 0}
    assert _counts(tmp_path, strand_bam, "reverse", ["-t", "exon", "-g", "gene_id", "-s", "2", "-Q", "0", "--primary"]) == {"GeneA": 0, "GeneOverlap": 0, "GeneB": 1}

    paired = _bam(tmp_path, "paired", [
        f"pair\t99\t1\t110\t60\t20M\t=\t150\t60\t{sequence}\t{quality}\tNH:i:1",
        f"pair\t147\t1\t150\t60\t20M\t=\t110\t-60\t{sequence}\t{quality}\tNH:i:1",
        f"orphan\t73\t1\t110\t60\t20M\t*\t0\t0\t{sequence}\t{quality}\tNH:i:1",
        f"chimera\t97\t1\t110\t60\t20M\t2\t110\t0\t{sequence}\t{quality}\tNH:i:1",
        f"chimera\t145\t2\t110\t60\t20M\t1\t110\t0\t{sequence}\t{quality}\tNH:i:1",
    ])
    assert _counts(tmp_path, _countable_bam(tmp_path, paired), "paired", ["-t", "exon", "-g", "gene_id", "-s", "0", "-Q", "0", "--primary", "-p", "--countReadPairs", "-B", "-C"]) == {"GeneA": 1, "GeneOverlap": 0, "GeneB": 0}

    lane_a = _bam(tmp_path, "lane_a", [f"lanea\t0\t1\t110\t60\t20M\t*\t0\t0\t{sequence}\t{quality}\tNH:i:1"])
    lane_b = _bam(tmp_path, "lane_b", [f"laneb\t0\t1\t110\t60\t20M\t*\t0\t0\t{sequence}\t{quality}\tNH:i:1"])
    _samtools(["merge", "-o", "sample.bam", lane_a.name, lane_b.name], tmp_path)
    _samtools(["index", "sample.bam"], tmp_path)
    assert _counts(tmp_path, _countable_bam(tmp_path, tmp_path / "sample.bam"), "merged", ["-t", "exon", "-g", "gene_id", "-s", "0", "-Q", "0", "--primary"])["GeneA"] == 2
