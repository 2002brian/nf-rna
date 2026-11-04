"""HISAT2 artifact lifecycle and L1 -> L2 upstream reuse (v1.7.0 storage/I/O defect).

Defect: HISAT2_ALIGN published every uncompressed lane SAM (about 242 GB for a
real 8-sample project), and ``--reuse-upstream`` / retry copied the whole
upstream tree, so an L2 freeze re-read and re-wrote all SAM and BAM files that
no downstream component reads.

Covered here with small synthetic fixtures:

* workflow: only the lane summary is published; alignment/counting commands unchanged;
* a real Nextflow run of the workflow with stand-in tools (when Nextflow is available);
* the selective copy: what is copied, what is referenced, atomicity, failure cleanup;
* production L1 -> design -> L2 reuse and retry paths with only Nextflow faked.
"""

from __future__ import annotations

import errno
import hashlib
import json
import os
import re
import shutil
import subprocess
import sys
import time
from copy import deepcopy
from pathlib import Path, PurePosixPath

import pytest
import yaml
from typer.testing import CliRunner

from conftest import base_config
from rnaseq import upstream_reuse
from rnaseq.cli import app
from rnaseq.errors import ExecutionPreflightError, UpstreamExecutionError
from rnaseq.planner import generate_plan
from rnaseq.service import execute_retry_service_run, execute_service_run
from rnaseq.upstream_reuse import REFERENCE_RECORD, UpstreamReuseCopyError, copy_reusable_upstream, is_referenced_artifact
from rnaseq.validators import validate_project
from test_runtime_dispatch import LINUX, _guard_conda, _host, mocked_downstream_runtime  # noqa: F401

ROOT = Path(__file__).resolve().parents[1]
WORKFLOW = ROOT / "workflow" / "hisat2_featurecounts.nf"
SAMPLES = ("C1", "C2", "T1", "T2")
runner = CliRunner()


# --------------------------------------------------------------------------- workflow contract


def _process(name: str) -> str:
    text = WORKFLOW.read_text(encoding="utf-8")
    match = re.search(rf"^process {name} \{{\n(.*?)^\}}", text, flags=re.S | re.M)
    assert match, name
    return match.group(1)


def test_hisat2_align_publishes_only_the_lane_summary_not_the_lane_sam():
    block = _process("HISAT2_ALIGN")
    (publish,) = [line.strip() for line in block.splitlines() if line.strip().startswith("publishDir")]
    assert publish == "publishDir \"${params.outdir}/alignment/lane_summaries\", mode: 'copy', overwrite: false, pattern: '*.hisat2.summary'"
    # The SAM is still an output: it is the channel hand-off to SORT_LANE_BAM.
    assert 'path("*.lane.sam"), path("*.hisat2.summary"), emit: aligned' in block
    assert "publishDir" not in _process("SORT_LANE_BAM")


def test_fastp_prepare_publishes_only_its_reports_not_the_prepared_reads():
    block = _process("FASTP_PREPARE")
    (publish,) = [line.strip() for line in block.splitlines() if line.strip().startswith("publishDir")]
    assert publish == "publishDir \"${params.outdir}/qc/fastp\", mode: 'copy', overwrite: false, pattern: 'fastp_*'"
    # The prepared reads are still the channel input of HISAT2_ALIGN and FASTQC_PROCESSED.
    assert "tuple val(sample), val(strandedness), path('prepared/*.fastq.gz'), emit: prepared" in block


def test_no_process_publishes_an_unfiltered_sam():
    text = WORKFLOW.read_text(encoding="utf-8")
    for name in re.findall(r"^process (\w+) \{", text, flags=re.M):
        block = _process(name)
        outputs = block.split("output:", 1)[-1].split("script:", 1)[0]
        if re.search(r'\.sam["\)]', outputs) and "publishDir" in block:
            assert "pattern:" in block, f"{name} publishes its SAM output"


def test_documented_alignment_and_count_artifacts_are_still_published():
    assert "publishDir \"${params.outdir}/bam\", mode: 'copy', overwrite: false" in _process("MERGE_AND_INDEX")
    assert "publishDir \"${params.outdir}/bam/count_only\", mode: 'copy', overwrite: false" in _process("PREPARE_COUNT_BAM")
    assert "publishDir \"${params.outdir}/counts/per_sample\", mode: 'copy', overwrite: false" in _process("FEATURECOUNTS")
    assert "publishDir \"${params.outdir}/qc/strandedness\", mode: 'copy', overwrite: false" in _process("STRAND_EVIDENCE")
    assert "publishDir \"${params.outdir}/counts\", mode: 'copy', overwrite: false" in _process("ASSEMBLE_COUNTS")
    assert "publishDir \"${params.outdir}/multiqc\", mode: 'copy', overwrite: false" in _process("MULTIQC")


def test_alignment_and_counting_commands_are_scientifically_unchanged():
    """Byte-level guard on every scientific command line of the route (scientific invariance)."""

    assert "hisat2 -p ${task.cpus} -x ${index}/${params.hisat2_index_basename} ${inputs} ${strand} ${spliceArgument} --summary-file ${lane}.hisat2.summary -S ${lane}.lane.sam" in _process("HISAT2_ALIGN")
    assert "samtools sort -@ ${task.cpus} -o ${sam.baseName}.bam ${sam}" in _process("SORT_LANE_BAM")
    assert "samtools view -@ ${task.cpus} -bh -F 0x900 -o ${sample}.countable.bam ${bam}" in _process("PREPARE_COUNT_BAM")
    assert "featureCounts -T ${task.cpus} -a ${gtf} -o ${sample}.counts.txt -t exon -g gene_id -s ${strand} -Q 0 --primary ${pairArgs} ${bam}" in _process("FEATURECOUNTS")
    assert "def pairArgs = paired ? '-p --countReadPairs -B -C' : ''" in _process("FEATURECOUNTS")
    # The diff of this fix against the workflow is confined to the HISAT2_ALIGN and FASTP_PREPARE publishDirs and their comments.
    original = subprocess.run(["git", "-C", str(ROOT), "show", "5cba156:workflow/hisat2_featurecounts.nf"], capture_output=True, text=True)
    if original.returncode != 0:
        pytest.skip("base revision 5cba156 is not available in this checkout")
    strip = lambda text: [line for line in text.splitlines() if not line.strip().startswith(("publishDir", "//"))]  # noqa: E731
    assert strip(original.stdout) == strip(WORKFLOW.read_text(encoding="utf-8"))


# --------------------------------------------------------------------------- real Nextflow run with stand-in tools

_FAKE_TOOLS = {
    "fastqc": r"""#!/bin/bash
for f in "$@"; do case "$f" in *.fastq.gz) b=$(basename "$f" .fastq.gz); echo html > "${b}_fastqc.html"; echo zip > "${b}_fastqc.zip";; esac; done
""",
    "hisat2": r"""#!/bin/bash
while [ $# -gt 0 ]; do case "$1" in -S) sam=$2; shift;; --summary-file) summary=$2; shift;; esac; shift; done
printf '@HD\tVN:1.0\nr1\t99\t1\t1\t60\t4M\t=\t1\t4\tACGT\tIIII\n' > "$sam"
printf '1 reads; of these:\n100.00%% overall alignment rate\n' > "$summary"
""",
    "samtools": r"""#!/bin/bash
sub=$1; shift
out=""; args=()
while [ $# -gt 0 ]; do case "$1" in -o) out=$2; shift;; -@) shift;; -*) ;; *) args+=("$1");; esac; shift; done
case "$sub" in
  sort|view) cp "${args[-1]}" "$out";;
  merge) cat "${args[@]}" > "$out";;
  index) echo bai > "${args[-1]}.bai";;
  flagstat) echo "1 + 0 in total (QC-passed reads + QC-failed reads)";;
  quickcheck) exit 0;;
esac
""",
    "featureCounts": r"""#!/bin/bash
while [ $# -gt 0 ]; do case "$1" in -o) out=$2; shift;; -s) strand=$2; shift;; -a|-T|-t|-g|-Q) shift;; -*) ;; *) bam=$1;; esac; shift; done
n=$([ "$strand" = 2 ] && echo 5000 || echo 100)
printf '# Program:featureCounts v2.0.6\nGeneid\tChr\tStart\tEnd\tStrand\tLength\t%s\nGeneA\t1\t1\t4\t+\t4\t%s\nGeneB\t1\t5\t8\t+\t4\t7\n' "$bam" "$n" > "$out"
printf 'Status\t%s\nAssigned\t%s\nUnassigned_NoFeatures\t3\n' "$bam" "$n" > "$out.summary"
""",
    "fastp": r"""#!/bin/bash
while [ $# -gt 0 ]; do case "$1" in -i) i1=$2; shift;; -I) i2=$2; shift;; -o) o1=$2; shift;; -O) o2=$2; shift;; -h) h=$2; shift;; -j) j=$2; shift;; esac; shift; done
cp "$i1" "$o1"; [ -n "$o2" ] && cp "$i2" "$o2"; echo '<html></html>' > "$h"; echo '{}' > "$j"
""",
    "multiqc": r"""#!/bin/bash
while [ $# -gt 0 ]; do case "$1" in --filename) name=$2; shift;; esac; shift; done
echo '<html></html>' > "$name"
""",
}


@pytest.mark.parametrize("pretrimmed", (True, False), ids=("pretrimmed", "fastp"))
def test_real_nextflow_run_publishes_no_lane_sam_and_keeps_every_essential_output(tmp_path, pretrimmed):
    if shutil.which("nextflow") is None:
        pytest.skip("Nextflow is not installed")
    tools = tmp_path / "bin"
    tools.mkdir()
    for name, body in {**_FAKE_TOOLS, "python": f"#!/bin/bash\nexec {sys.executable} \"$@\"\n"}.items():
        (tools / name).write_text(body, encoding="utf-8")
        (tools / name).chmod(0o755)
    data = tmp_path / "data"
    (data / "index").mkdir(parents=True)
    for number in range(1, 9):
        (data / "index" / f"genome.{number}.ht2").write_bytes(b"x")
    (data / "genome.fa").write_text(">1\nACGTACGT\n", encoding="utf-8")
    (data / "genes.gtf").write_text("1\tt\texon\t1\t4\t.\t+\t.\tgene_id \"GeneA\";\n", encoding="utf-8")
    rows = ["sample,fastq_1,fastq_2,strandedness"]
    for sample in ("S1", "S2"):
        for read in ("R1", "R2"):
            (data / f"{sample}_{read}.fastq.gz").write_bytes(b"@r\nACGT\n+\nIIII\n")
        rows.append(f"{sample},{data / f'{sample}_R1.fastq.gz'},{data / f'{sample}_R2.fastq.gz'},reverse")
    (data / "samplesheet.csv").write_text("\n".join(rows) + "\n", encoding="utf-8")
    outdir = tmp_path / "out"
    env = {**os.environ, "PATH": f"{tools}:{os.environ['PATH']}", "NXF_ANSI_LOG": "false"}
    result = subprocess.run(
        ["nextflow", "run", str(WORKFLOW), "-work-dir", str(tmp_path / "work"),
         "--input", str(data / "samplesheet.csv"), "--outdir", str(outdir), "--fasta", str(data / "genome.fa"),
         "--gtf", str(data / "genes.gtf"), "--hisat2_index", str(data / "index"),
         "--assembly_script", str(ROOT / "src" / "rnaseq" / "hisat2_featurecounts.py"),
         "--layout", "paired_end", "--strandedness", "reverse", "--pretrimmed", "true" if pretrimmed else "false"],
        cwd=tmp_path, env=env, capture_output=True, text=True, timeout=600,
    )
    assert result.returncode == 0, result.stdout[-4000:] + result.stderr[-4000:]
    published = sorted(path.relative_to(outdir).as_posix() for path in outdir.rglob("*") if path.is_file())
    assert not [path for path in published if path.endswith(".sam")], published
    assert not [path for path in published if path.endswith((".fastq.gz", ".fq.gz"))], published
    fastp_reports = {path for path in published if path.startswith("qc/fastp/")}
    assert fastp_reports == (set() if pretrimmed else {"qc/fastp/fastp_S1_R1.fastq.html", "qc/fastp/fastp_S1_R1.fastq.json", "qc/fastp/fastp_S2_R1.fastq.html", "qc/fastp/fastp_S2_R1.fastq.json"})
    for sample in ("S1", "S2"):
        for expected in (
            f"alignment/lane_summaries/{sample}_R1.fastq.fastq.hisat2.summary",
            f"bam/{sample}.bam", f"bam/{sample}.bam.bai", f"bam/{sample}.flagstat.txt",
            f"bam/count_only/{sample}.countable.bam", f"bam/count_only/{sample}.countable.bam.bai",
            f"counts/per_sample/{sample}.counts.txt", f"counts/per_sample/{sample}.counts.txt.summary",
            f"qc/strandedness/{sample}.strand_s1.txt.summary", f"qc/strandedness/{sample}.strand_s2.txt.summary",
        ):
            assert expected in published, (expected, published)
    assert {"counts/canonical_counts.csv", "counts/sample_map.csv", "multiqc/multiqc_report.html"} <= set(published)
    # The SAM still reached SORT_LANE_BAM: the published BAM is the SAM's content.
    assert (outdir / "bam" / "S1.bam").read_text(encoding="utf-8").startswith("@HD")
    assert (outdir / "counts" / "canonical_counts.csv").read_text(encoding="utf-8") == "gene_id,S1,S2\nGeneA,5000,5000\nGeneB,7,7\n"
    # The lane SAM exists only in the HISAT2_ALIGN task work directories (SORT_LANE_BAM stages it as a symlink).
    sams = list((tmp_path / "work").rglob("*.lane.sam"))
    assert len([path for path in sams if not path.is_symlink()]) == 2
    assert all(path.resolve().is_file() for path in sams)


# --------------------------------------------------------------------------- selective copy (unit level)


def _write(path: Path, content: bytes | str) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(content.encode() if isinstance(content, str) else content)
    return path


def _hisat2_upstream(run_dir: Path, samples=SAMPLES, *, historical_sam: bool = True, alignment_bytes: int = 4096) -> Path:
    """A HISAT2 + featureCounts upstream tree as published; ``historical_sam`` mimics runs before this fix."""

    out = run_dir / "upstream" / "hisat2_featurecounts"
    header = "gene_id," + ",".join(samples)
    _write(out / "counts" / "canonical_counts.csv", header + "\nGeneA," + ",".join(str(10 + i) for i, _ in enumerate(samples)) + "\nGeneB," + ",".join("3" for _ in samples) + "\n")
    _write(out / "counts" / "sample_map.csv", "sample_id,featurecounts_file\n" + "".join(f"{s},{s}.counts.txt\n" for s in samples))
    _write(out / "multiqc" / "multiqc_report.html", "<html></html>")
    _write(out / "qc" / "fastqc" / "raw" / "raw_x_fastqc.html", "<html></html>")
    for index, sample in enumerate(samples):
        _write(out / "counts" / "per_sample" / f"{sample}.counts.txt", f"# featureCounts\nGeneid\tChr\tStart\tEnd\tStrand\tLength\t{sample}.countable.bam\nGeneA\t1\t1\t4\t+\t4\t{10 + index}\n")
        _write(out / "counts" / "per_sample" / f"{sample}.counts.txt.summary", f"Status\t{sample}.countable.bam\nAssigned\t{10 + index}\n")
        _write(out / "qc" / "strandedness" / f"{sample}.strand_s1.txt.summary", f"Status\tx\nAssigned\t100\n")
        _write(out / "qc" / "strandedness" / f"{sample}.strand_s2.txt.summary", f"Status\tx\nAssigned\t5000\n")
        _write(out / "alignment" / "lane_summaries" / f"{sample}_R1.fastq.hisat2.summary", "90.00% overall alignment rate\n")
        if historical_sam:
            _write(out / "alignment" / "lane_summaries" / f"{sample}_R1.fastq.lane.sam", os.urandom(alignment_bytes))
            for read in ("R1", "R2"):
                _write(out / "qc" / "fastp" / "prepared" / f"{sample}_{read}.fastq.fastq.gz", os.urandom(alignment_bytes))
        _write(out / "bam" / f"{sample}.bam", os.urandom(alignment_bytes))
        _write(out / "bam" / f"{sample}.bam.bai", os.urandom(64))
        _write(out / "bam" / f"{sample}.flagstat.txt", "1 + 0 in total\n")
        _write(out / "bam" / "count_only" / f"{sample}.countable.bam", os.urandom(alignment_bytes))
        _write(out / "bam" / "count_only" / f"{sample}.countable.bam.bai", os.urandom(64))
        _write(out / "bam" / "count_only" / f"{sample}.flagstat.txt", "1 + 0 in total\n")
    _write(run_dir / "handoff" / "upstream_manifest.yaml", yaml.safe_dump({"pipeline": {"name": "nf-rna/hisat2_featurecounts", "workflow_sha256": "a" * 64}}))
    _write(out / "._appledouble", "ignored")
    return out


def _tree(root: Path) -> dict[str, str]:
    return {path.relative_to(root).as_posix(): hashlib.sha256(path.read_bytes()).hexdigest() for path in root.rglob("*") if path.is_file()}


@pytest.mark.parametrize("historical_sam", (True, False), ids=("historical-L1-with-lane-SAM", "fixed-L1"))
def test_selective_copy_keeps_every_small_file_and_references_alignment_artifacts(tmp_path, historical_sam):
    source_run = tmp_path / "runs" / "L1" / "R1"
    _hisat2_upstream(source_run, historical_sam=historical_sam)
    destination = tmp_path / "runs" / "L2" / "R2" / "upstream"
    destination.mkdir(parents=True)  # create_case_run makes it empty
    summary = copy_reusable_upstream(source_run, destination, source_label="L1/R1")

    source_tree, copied_tree = _tree(source_run / "upstream"), _tree(destination)
    referenced = {path for path in source_tree if is_referenced_artifact(PurePosixPath(path))}
    expected_referenced = 4 * (7 if historical_sam else 4)
    assert len(referenced) == expected_referenced
    # Every non-alignment file is byte-identical; no alignment artifact was copied.
    small = {path: digest for path, digest in source_tree.items() if path not in referenced and not path.endswith("._appledouble")}
    record_path = f"hisat2_featurecounts/{REFERENCE_RECORD}"
    assert {path: digest for path, digest in copied_tree.items() if path != record_path} == small
    assert not [path for path in copied_tree if path.endswith((".sam", ".bam", ".bai", ".fastq.gz"))]
    # Copies are independent files: no hard links, no symbolic links.
    for path in destination.rglob("*"):
        assert not path.is_symlink()
        if path.is_file():
            assert path.stat().st_nlink == 1
    record = json.loads((destination / record_path).read_text(encoding="utf-8"))
    assert record["schema_version"] == "nf-rna.referenced-alignment-artifacts.v1"
    assert record["reused_from"] == "L1/R1"
    assert record["source_pipeline"] == {"name": "nf-rna/hisat2_featurecounts", "workflow_sha256": "a" * 64}
    assert {item["path"] for item in record["artifacts"]} == referenced
    for item in record["artifacts"]:
        origin = Path(item["source_path"])
        assert origin == (source_run / "upstream" / item["path"]).resolve()
        assert item["size_bytes"] == origin.stat().st_size and item["source_run"] == "L1/R1"
        assert item["kind"] == ("prepared_reads" if item["path"].endswith(".fastq.gz") else "alignment")
    assert summary["referenced_files"] == expected_referenced
    assert summary["copied_files"] == len(small)
    assert summary["copied_bytes"] == sum((source_run / "upstream" / path).stat().st_size for path in small)
    # The source run is untouched.
    assert _tree(source_run / "upstream") == source_tree
    assert not list(destination.parent.glob(".upstream.partial-*"))


def test_salmon_upstream_is_still_copied_completely(tmp_path):
    source_run = tmp_path / "L1"
    out = source_run / "upstream" / "nfcore_rnaseq"
    _write(out / "salmon" / "salmon.merged.gene_counts.tsv", "gene_id\tS1\nGeneA\t1.0\n")
    _write(out / "salmon" / "S1" / "quant.sf", "Name\tLength\nTx1\t100\n")
    _write(out / "star_salmon" / "S1.markdup.sorted.bam", os.urandom(256))  # copied as before: not the HISAT2 route
    _write(out / "multiqc" / "multiqc_report.html", "<html></html>")
    destination = tmp_path / "L2" / "upstream"
    copy_reusable_upstream(source_run, destination, source_label="L1/R")
    assert _tree(destination) == _tree(source_run / "upstream")
    assert not (destination / "hisat2_featurecounts").exists()


def test_chained_reuse_keeps_the_original_reference_and_its_owner(tmp_path):
    l1 = tmp_path / "L1"
    _hisat2_upstream(l1)
    l2 = tmp_path / "L2"
    copy_reusable_upstream(l1, l2 / "upstream", source_label="L1/R1")
    _write(l2 / "handoff" / "upstream_manifest.yaml", yaml.safe_dump({"pipeline": {"name": "nf-rna/hisat2_featurecounts", "workflow_sha256": "b" * 64}}))
    retry = tmp_path / "L2-retry"
    summary = copy_reusable_upstream(l2, retry / "upstream", source_label="L2/R2")
    first = json.loads((l2 / "upstream" / "hisat2_featurecounts" / REFERENCE_RECORD).read_text(encoding="utf-8"))
    second = json.loads((retry / "upstream" / "hisat2_featurecounts" / REFERENCE_RECORD).read_text(encoding="utf-8"))
    assert second["artifacts"] == first["artifacts"]  # still point at L1, which holds the files
    assert second["reused_from"] == "L1/R1"
    assert second["source_pipeline"]["workflow_sha256"] == "a" * 64  # the producing workflow, not the intermediate run's
    assert summary["newly_referenced_files"] == 0


def test_a_pre_fix_full_copy_used_as_retry_source_is_referenced_to_that_run(tmp_path):
    """An L2 made before this fix holds a full copy; a retry references the copy it actually holds."""

    old_l2 = tmp_path / "old-L2"
    _hisat2_upstream(old_l2)
    summary = copy_reusable_upstream(old_l2, tmp_path / "retry" / "upstream", source_label="L2/OLD")
    record = json.loads((tmp_path / "retry" / "upstream" / "hisat2_featurecounts" / REFERENCE_RECORD).read_text(encoding="utf-8"))
    assert {item["source_run"] for item in record["artifacts"]} == {"L2/OLD"}
    assert summary["newly_referenced_files"] == summary["referenced_files"] == 28


@pytest.mark.parametrize("failure", ("cross-device", "disk-full", "short-copy"))
def test_failed_copy_leaves_no_partial_upstream(tmp_path, monkeypatch, failure):
    source_run = tmp_path / "L1"
    _hisat2_upstream(source_run)
    destination = tmp_path / "L2" / "upstream"
    destination.mkdir(parents=True)
    real_copy = shutil.copy2
    calls = {"n": 0}

    def flaky(origin, target, *args, **kwargs):
        calls["n"] += 1
        if calls["n"] == 3:
            if failure == "short-copy":
                Path(target).write_bytes(b"")
                return target
            raise OSError(errno.EXDEV if failure == "cross-device" else errno.ENOSPC, os.strerror(errno.EXDEV if failure == "cross-device" else errno.ENOSPC))
        return real_copy(origin, target, *args, **kwargs)

    monkeypatch.setattr(upstream_reuse.shutil, "copy2", flaky)
    with pytest.raises(UpstreamReuseCopyError):
        copy_reusable_upstream(source_run, destination, source_label="L1/R1")
    assert destination.is_dir() and not any(destination.iterdir())
    assert not list(destination.parent.glob(".upstream.partial-*"))


def test_refuses_a_non_empty_destination_and_a_missing_source(tmp_path):
    source_run = tmp_path / "L1"
    _hisat2_upstream(source_run)
    destination = tmp_path / "L2" / "upstream"
    _write(destination / "stale.txt", "x")
    with pytest.raises(UpstreamReuseCopyError, match="non-empty"):
        copy_reusable_upstream(source_run, destination, source_label="L1/R1")
    assert (destination / "stale.txt").read_text(encoding="utf-8") == "x"
    with pytest.raises(UpstreamReuseCopyError, match="no upstream output"):
        copy_reusable_upstream(tmp_path / "missing", tmp_path / "other" / "upstream", source_label="L1/X")


def test_invalid_inherited_reference_record_is_rejected(tmp_path):
    source_run = tmp_path / "L2"
    out = _hisat2_upstream(source_run)
    _write(out / REFERENCE_RECORD, "{not json")
    with pytest.raises(UpstreamReuseCopyError, match="unreadable"):
        copy_reusable_upstream(source_run, tmp_path / "retry" / "upstream", source_label="L2/R")


# --------------------------------------------------------------------------- production L1 -> L2 reuse and retry


def _hisat2_project(root: Path) -> None:
    fastq = root / "input" / "fastq"
    fastq.mkdir(parents=True)
    for sample in SAMPLES:
        for read in ("R1", "R2"):
            (fastq / f"{sample}_{read}.fastq.gz").write_bytes(f"{sample}{read}".encode())
    (root / "metadata.csv").write_text("sample_id,condition\nC1,C\nC2,C\nT1,T\nT2,T\n", encoding="utf-8")
    (root / "contrasts.csv").write_text("contrast_id,factor,numerator,denominator\n", encoding="utf-8")
    reference = root / "reference"
    (reference / "index").mkdir(parents=True)
    (reference / "genome.fa").write_text(">1\nACGT\n", encoding="utf-8")
    (reference / "genes.gtf").write_text("1\tt\texon\t1\t4\t.\t+\t.\tgene_id \"GeneA\";\n", encoding="utf-8")
    for number in range(1, 9):
        (reference / "index" / f"genome.{number}.ht2").write_bytes(b"x")
    config = deepcopy(base_config())
    config.update({
        "schema_version": "1.2", "input": {"type": "fastq", "path": "input/fastq", "layout": "paired_end", "preprocessing": "pretrimmed"},
        "upstream": {"engine": "nfcore_rnaseq", "pipeline_version": "3.26.0", "strandedness": "reverse", "quantification": {"method": "hisat2_featurecounts"}},
        "reference": {"source": "custom", "fasta": "reference/genome.fa", "gtf": "reference/genes.gtf", "hisat2_index": "reference/index"},
        "analysis": {"enrichment": []},
    })
    config["project"]["preset"] = "L1"
    (root / "project.yaml").write_text(yaml.safe_dump(config, sort_keys=False), encoding="utf-8")


def _fake_nextflow(observed: list[list[str]], *, historical_sam: bool, downstream_rc: int = 0):
    def run(command, *, cwd, stdout_path, stderr_path):
        observed.append(command)
        stdout_path.write_text("", encoding="utf-8")
        stderr_path.write_text("", encoding="utf-8")
        outdir = Path(command[command.index("--outdir") + 1])
        if any(token.endswith("hisat2_featurecounts.nf") for token in command):
            _hisat2_upstream(outdir.parent.parent, historical_sam=historical_sam)
            return 0
        (outdir / "report").mkdir(parents=True, exist_ok=True)
        (outdir / "report" / "report.html").write_text("<html></html>", encoding="utf-8")
        return downstream_rc
    return run


@pytest.fixture
def hisat2_l1(monkeypatch, tmp_path, mocked_downstream_runtime, production_capable_execution_capacity, request):
    historical_sam = getattr(request, "param", True)
    monkeypatch.setenv("RNASEQ_EXECUTION_ROOT", str(tmp_path / "execution"))
    _host(monkeypatch, *LINUX)
    _guard_conda(monkeypatch)
    root = tmp_path / "project"
    _hisat2_project(root)
    report = validate_project(root)
    assert report.is_valid and report.execution_ready, report.errors
    generate_plan(report)
    observed: list[list[str]] = []
    monkeypatch.setattr("rnaseq.service._run_command", _fake_nextflow(observed, historical_sam=historical_sam))
    l1 = execute_service_run(report, case_id="L1-001")
    assert json.loads(l1.state_path.read_text(encoding="utf-8"))["status"] == "SUCCESS"
    result = runner.invoke(app, ["design", str(root), "--no-enrichment", "--yes"])
    assert result.exit_code == 0, result.output
    return root, l1


def _l2(monkeypatch, root, l1, *, downstream_rc: int = 0, case_id: str = "L2-001"):
    report = validate_project(root)
    assert report.is_valid, report.errors
    generate_plan(report)
    observed: list[list[str]] = []
    monkeypatch.setattr("rnaseq.service._run_command", _fake_nextflow(observed, historical_sam=True, downstream_rc=downstream_rc))
    return report, observed


@pytest.mark.parametrize("hisat2_l1", (True, False), ids=("historical-L1", "fixed-L1"), indirect=True)
def test_l2_reuse_succeeds_without_copying_alignment_artifacts(monkeypatch, hisat2_l1):
    root, l1 = hisat2_l1
    report, observed = _l2(monkeypatch, root, l1)
    l2 = execute_service_run(report, case_id="L2-001", reuse_upstream=f"L1-001/{l1.run_id}")
    assert json.loads(l2.state_path.read_text(encoding="utf-8"))["status"] == "SUCCESS"
    assert len(observed) == 1 and not any(token.endswith("hisat2_featurecounts.nf") for token in observed[0])

    l1_out, l2_out = l1.run_dir / "upstream" / "hisat2_featurecounts", l2.run_dir / "upstream" / "hisat2_featurecounts"
    assert not [path for path in l2_out.rglob("*") if path.name.endswith((".sam", ".bam", ".bai", ".fastq.gz"))]
    # Small required files are available to L2, byte-identical.
    for relative in ("counts/canonical_counts.csv", "counts/sample_map.csv", "multiqc/multiqc_report.html",
                     "counts/per_sample/C1.counts.txt.summary", "qc/strandedness/T2.strand_s2.txt.summary",
                     "alignment/lane_summaries/C1_R1.fastq.hisat2.summary", "bam/C1.flagstat.txt"):
        assert (l2_out / relative).read_bytes() == (l1_out / relative).read_bytes(), relative
    # Scientific invariance: the L2 staged and delivered count matrix equals L1's.
    assert (l2.run_dir / "delivery" / "counts" / "raw_counts.csv").read_bytes() == (l1_out / "counts" / "canonical_counts.csv").read_bytes()
    assert (l2.run_dir / "downstream_inputs" / "source" / "canonical_counts.csv").read_bytes() == (l1.run_dir / "downstream_inputs" / "source" / "canonical_counts.csv").read_bytes()
    # The handoff records where the alignment artifacts live and which workflow produced them.
    handoff = yaml.safe_load((l2.run_dir / "frozen" / "upstream_handoff_manifest.yaml").read_text(encoding="utf-8"))
    assert handoff["alignment"]["storage"] == "referenced_in_source_run"
    assert handoff["alignment"]["reused_from"] == f"L1-001/{l1.run_id}"
    assert handoff["alignment"]["referenced_files"] == 4 * (7 if any(l1_out.rglob("*.lane.sam")) else 4)
    l1_handoff = yaml.safe_load((l1.run_dir / "handoff" / "upstream_manifest.yaml").read_text(encoding="utf-8"))
    assert handoff["reuse"]["source_pipeline"] == l1_handoff["pipeline"]
    assert handoff["strandedness_qc"]["status"] == l1_handoff["strandedness_qc"]["status"]
    assert {k: v for k, v in handoff["counting"].items()} == l1_handoff["counting"]
    log = (l2.run_dir / "logs" / "rnaseq.log").read_text(encoding="utf-8")
    assert "alignment artifacts" in log and "referenced in their source run, not copied" in log
    # The source run is unchanged and still holds every alignment artifact.
    assert all(Path(item["source_path"]).is_file() for item in json.loads((l2_out / REFERENCE_RECORD).read_text(encoding="utf-8"))["artifacts"])


def test_identity_checks_still_reject_incompatible_reuse_before_copying(monkeypatch, hisat2_l1):
    root, l1 = hisat2_l1
    manifest_path = l1.run_dir / "frozen" / "input_manifest.yaml"
    manifest = yaml.safe_load(manifest_path.read_text(encoding="utf-8"))
    manifest["input"]["files"][0]["sha256"] = "0" * 64
    manifest_path.write_text(yaml.safe_dump(manifest, sort_keys=False), encoding="utf-8")
    report, observed = _l2(monkeypatch, root, l1)
    with pytest.raises(ExecutionPreflightError, match="input_manifest.yaml differs"):
        execute_service_run(report, case_id="L2-001", reuse_upstream=f"L1-001/{l1.run_id}")
    assert observed == []
    (l2_dir,) = (root / "runs" / "L2-001").iterdir()
    assert json.loads((l2_dir / "run_state.json").read_text(encoding="utf-8"))["status"] == "FAILED"
    assert not any((l2_dir / "upstream").iterdir())


def test_incomplete_source_artifacts_fail_the_run_and_never_report_success(monkeypatch, hisat2_l1):
    root, l1 = hisat2_l1
    (l1.run_dir / "upstream" / "hisat2_featurecounts" / "counts" / "canonical_counts.csv").unlink()
    report, observed = _l2(monkeypatch, root, l1)
    with pytest.raises(UpstreamExecutionError, match="required artifacts"):
        execute_service_run(report, case_id="L2-001", reuse_upstream=f"L1-001/{l1.run_id}")
    assert observed == []
    (l2_dir,) = (root / "runs" / "L2-001").iterdir()
    assert json.loads((l2_dir / "run_state.json").read_text(encoding="utf-8"))["status"] == "FAILED"


def test_interrupted_copy_fails_the_run_without_a_partial_upstream(monkeypatch, hisat2_l1):
    root, l1 = hisat2_l1
    report, observed = _l2(monkeypatch, root, l1)
    real_copy = shutil.copy2
    calls = {"n": 0}

    def failing(origin, target, *args, **kwargs):
        calls["n"] += 1
        if calls["n"] == 5:
            raise OSError(errno.ENOSPC, "No space left on device")
        return real_copy(origin, target, *args, **kwargs)

    monkeypatch.setattr(upstream_reuse.shutil, "copy2", failing)
    with pytest.raises(OSError, match="No space left"):
        execute_service_run(report, case_id="L2-001", reuse_upstream=f"L1-001/{l1.run_id}")
    (l2_dir,) = (root / "runs" / "L2-001").iterdir()
    state = json.loads((l2_dir / "run_state.json").read_text(encoding="utf-8"))
    assert state["status"] == "FAILED" and "No space left" in state["error"]
    assert not any((l2_dir / "upstream").iterdir())
    assert not list(l2_dir.glob(".upstream.partial-*"))
    assert observed == []


def test_retry_of_a_failed_reusing_l2_keeps_the_reference_chain(monkeypatch, hisat2_l1):
    root, l1 = hisat2_l1
    report, _observed = _l2(monkeypatch, root, l1, downstream_rc=1)
    with pytest.raises(Exception):
        execute_service_run(report, case_id="L2-001", reuse_upstream=f"L1-001/{l1.run_id}")
    (failed,) = (root / "runs" / "L2-001").iterdir()
    assert json.loads((failed / "run_state.json").read_text(encoding="utf-8"))["status"] == "FAILED"
    # Pre-existing retry rule, unrelated to this fix: a custom (unmanaged) reference must have each
    # asset copied as a file into frozen/reference; a custom HISAT2 index *directory* never is.  The
    # reuse path never runs HISAT2, so a placeholder satisfies the check for this fixture only.
    (failed / "frozen" / "reference").mkdir(exist_ok=True)
    (failed / "frozen" / "reference" / "index").write_bytes(b"placeholder")
    observed: list[list[str]] = []
    monkeypatch.setattr("rnaseq.service._run_command", _fake_nextflow(observed, historical_sam=True))
    retry = execute_retry_service_run(root, retry_of=f"L2-001/{failed.name}")
    assert json.loads(retry.state_path.read_text(encoding="utf-8"))["status"] == "SUCCESS"
    out = retry.run_dir / "upstream" / "hisat2_featurecounts"
    assert not [path for path in out.rglob("*") if path.name.endswith((".sam", ".bam", ".bai"))]
    record = json.loads((out / REFERENCE_RECORD).read_text(encoding="utf-8"))
    assert record["reused_from"] == f"L1-001/{l1.run_id}"
    assert {item["source_run"] for item in record["artifacts"]} == {f"L1-001/{l1.run_id}"}
    assert (retry.run_dir / "delivery" / "counts" / "raw_counts.csv").read_bytes() == (l1.run_dir / "upstream" / "hisat2_featurecounts" / "counts" / "canonical_counts.csv").read_bytes()


# --------------------------------------------------------------------------- performance (measured on a small fixture)


def _legacy_copy(source_run: Path, destination: Path) -> None:
    shutil.copytree(source_run / "upstream", destination, dirs_exist_ok=True, ignore=shutil.ignore_patterns("._*"))


def test_reuse_copy_volume_before_and_after(tmp_path, capsys):
    """Measured on an 8-sample fixture with 1 MiB alignment files; prints the comparison table."""

    samples = tuple(f"S{i}" for i in range(1, 9))
    source_run = tmp_path / "L1"
    _hisat2_upstream(source_run, samples, alignment_bytes=1 << 20)
    measurements = {}
    for label, action in (("legacy copytree", lambda dest: _legacy_copy(source_run, dest)),
                          ("selective copy", lambda dest: copy_reusable_upstream(source_run, dest, source_label="L1/R"))):
        dest = tmp_path / label.replace(" ", "_") / "upstream"
        started = time.perf_counter()
        action(dest)
        elapsed = time.perf_counter() - started
        files = [path for path in dest.rglob("*") if path.is_file()]
        measurements[label] = (len(files), sum(path.stat().st_size for path in files), elapsed)
    source_bytes = sum(path.stat().st_size for path in (source_run / "upstream").rglob("*") if path.is_file())
    legacy, selective = measurements["legacy copytree"], measurements["selective copy"]
    with capsys.disabled():
        print(f"\nsource upstream bytes={source_bytes}")
        for label, (files, size, elapsed) in measurements.items():
            print(f"{label}: files={files} bytes={size} amplification={size / source_bytes:.4f} seconds={elapsed:.4f}")
    assert legacy[1] >= 8 * 3 * (1 << 20)
    assert selective[1] < 64 * 1024  # only summaries, counts, QC and the reference record
    assert selective[1] / legacy[1] < 0.01
