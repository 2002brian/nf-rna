"""Pinned-tool acceptance checks for the featureCounts contract.

These fixtures deliberately run the same count-only BAM transformation and
featureCounts arguments as ``workflow/hisat2_featurecounts.nf`` with the exact
SAMtools and Subread builds of the qualified linux-64 HISAT2/featureCounts
Conda environment, so they prove the production tool versions rather than
substituting whatever host executables happen to be on PATH.

The environment is the one Nextflow provisions in the shared upstream Conda
cache on the first HISAT2 run.  Set ``NF_RNA_HISAT2_CONDA_PREFIX`` to point at
another prefix built from ``workflow/envs/hisat2-featurecounts-linux-64.yml``.
"""

from __future__ import annotations

import csv
import os
import subprocess
from pathlib import Path

import pytest
import yaml

from rnaseq.execution import HISAT2_LINUX_CONDA_ENV, upstream_conda_cache
from rnaseq.hisat2_featurecounts import featurecounts_arguments


PREFIX_ENV = "NF_RNA_HISAT2_CONDA_PREFIX"


def _pinned_records() -> set[str]:
    """Conda-meta record names of the exactly pinned SAMtools and Subread builds."""

    dependencies = yaml.safe_load(HISAT2_LINUX_CONDA_ENV.read_text(encoding="utf-8"))["dependencies"]
    records = set()
    for spec in dependencies:
        name, version, build = spec.split("::", 1)[1].split("=")
        if name in {"samtools", "subread"}:
            records.add(f"{name}-{version}-{build}.json")
    assert len(records) == 2, dependencies
    return records


def _is_pinned_prefix(prefix: Path, records: set[str]) -> bool:
    return records <= {path.name for path in (prefix / "conda-meta").glob("*.json")} and all(
        (prefix / "bin" / tool).is_file() for tool in ("samtools", "featureCounts")
    )


def _pinned_prefix() -> Path | None:
    records = _pinned_records()
    configured = os.environ.get(PREFIX_ENV)
    if configured:
        prefix = Path(configured).expanduser()
        if not _is_pinned_prefix(prefix, records):
            pytest.fail(f"{PREFIX_ENV}={prefix} does not contain the pinned {sorted(records)} builds.")
        return prefix
    cache = upstream_conda_cache()
    if not cache.is_dir():
        return None
    return next((prefix for prefix in sorted(cache.iterdir()) if _is_pinned_prefix(prefix, records)), None)


PREFIX = _pinned_prefix()
pytestmark = pytest.mark.skipif(
    PREFIX is None,
    reason=f"The pinned HISAT2/featureCounts Conda environment is not provisioned; run a HISAT2 project once or set {PREFIX_ENV}",
)


def _run(tool: str, arguments: list[str], cwd: Path) -> subprocess.CompletedProcess[str]:
    assert PREFIX is not None
    result = subprocess.run([str(PREFIX / "bin" / tool), *arguments], cwd=cwd, capture_output=True, text=True, check=False)
    assert result.returncode == 0, result.stderr
    return result


def _samtools(arguments: list[str], cwd: Path) -> subprocess.CompletedProcess[str]:
    return _run("samtools", arguments, cwd)


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
    _run("featureCounts", ["-a", "genes.gtf", "-o", output.name, *arguments, bam.name], tmp_path)
    with output.open(encoding="utf-8", newline="") as handle:
        rows = [row for row in csv.reader((line for line in handle if not line.startswith("#")), delimiter="\t")]
    return {row[0]: int(row[6]) for row in rows[1:]}


def test_pinned_tool_versions_are_the_production_builds(tmp_path: Path):
    assert _samtools(["--version"], tmp_path).stdout.splitlines()[0] == "samtools 1.21"
    assert PREFIX is not None
    version = subprocess.run([str(PREFIX / "bin" / "featureCounts"), "-v"], capture_output=True, text=True, check=False)
    assert "featureCounts v2.0.6" in version.stdout + version.stderr


def test_pinned_conda_featurecounts_counting_policy(tmp_path: Path):
    single_end = ["-t", "exon", "-g", "gene_id", "-s", "0", "-Q", "0", "--primary"]
    forward = ["-t", "exon", "-g", "gene_id", "-s", "1", "-Q", "0", "--primary"]
    reverse = ["-t", "exon", "-g", "gene_id", "-s", "2", "-Q", "0", "--primary"]
    paired_end = [*single_end, "-p", "--countReadPairs", "-B", "-C"]
    # The literal policy under test is exactly what the production workflow receives.
    assert featurecounts_arguments(layout="single_end", strandedness="unstranded") == single_end
    assert featurecounts_arguments(layout="single_end", strandedness="forward") == forward
    assert featurecounts_arguments(layout="single_end", strandedness="reverse") == reverse
    assert featurecounts_arguments(layout="paired_end", strandedness="unstranded") == paired_end

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
    original = _samtools(["view", single.name], tmp_path)
    filtered = _samtools(["view", countable.name], tmp_path)
    assert "multi\t0" in filtered.stdout and "NH:i:2" in filtered.stdout
    assert "multi\t256" in original.stdout and "multi\t256" not in filtered.stdout
    assert "secondary" not in filtered.stdout and "supplementary" not in filtered.stdout
    assert _counts(tmp_path, countable, "single", single_end) == {"GeneA": 1, "GeneOverlap": 0, "GeneB": 0}

    strands = _bam(tmp_path, "strands", [
        f"forward\t0\t1\t110\t60\t20M\t*\t0\t0\t{sequence}\t{quality}\tNH:i:1",
        f"reverse\t16\t1\t310\t60\t20M\t*\t0\t0\t{sequence}\t{quality}\tNH:i:1",
    ])
    strand_bam = _countable_bam(tmp_path, strands)
    assert _counts(tmp_path, strand_bam, "forward", forward) == {"GeneA": 1, "GeneOverlap": 0, "GeneB": 0}
    assert _counts(tmp_path, strand_bam, "reverse", reverse) == {"GeneA": 0, "GeneOverlap": 0, "GeneB": 1}

    paired = _bam(tmp_path, "paired", [
        f"pair\t99\t1\t110\t60\t20M\t=\t150\t60\t{sequence}\t{quality}\tNH:i:1",
        f"pair\t147\t1\t150\t60\t20M\t=\t110\t-60\t{sequence}\t{quality}\tNH:i:1",
        f"orphan\t73\t1\t110\t60\t20M\t*\t0\t0\t{sequence}\t{quality}\tNH:i:1",
        f"chimera\t97\t1\t110\t60\t20M\t2\t110\t0\t{sequence}\t{quality}\tNH:i:1",
        f"chimera\t145\t2\t110\t60\t20M\t1\t110\t0\t{sequence}\t{quality}\tNH:i:1",
    ])
    assert _counts(tmp_path, _countable_bam(tmp_path, paired), "paired", paired_end) == {"GeneA": 1, "GeneOverlap": 0, "GeneB": 0}

    lane_a = _bam(tmp_path, "lane_a", [f"lanea\t0\t1\t110\t60\t20M\t*\t0\t0\t{sequence}\t{quality}\tNH:i:1"])
    lane_b = _bam(tmp_path, "lane_b", [f"laneb\t0\t1\t110\t60\t20M\t*\t0\t0\t{sequence}\t{quality}\tNH:i:1"])
    _samtools(["merge", "-o", "sample.bam", lane_a.name, lane_b.name], tmp_path)
    _samtools(["index", "sample.bam"], tmp_path)
    assert _counts(tmp_path, _countable_bam(tmp_path, tmp_path / "sample.bam"), "merged", single_end)["GeneA"] == 2
