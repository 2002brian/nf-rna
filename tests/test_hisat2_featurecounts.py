from __future__ import annotations

from copy import deepcopy
from pathlib import Path

import pytest
import yaml

from conftest import base_config
from rnaseq.execution import build_hisat2_featurecounts_command, resolved_upstream_implementation
from rnaseq.hisat2_featurecounts import HISAT2_IMAGE, HISAT2_VERSION, assemble_count_matrix, featurecounts_arguments, hisat2_strand_option
from rnaseq.service import CaseRun, freeze_case_inputs
from rnaseq.validators import validate_project


def _featurecounts(path: Path, values: list[tuple[str, int]]) -> None:
    path.write_text(
        "# Program:featureCounts\nGeneid\tChr\tStart\tEnd\tStrand\tLength\tdeclared.bam\n"
        + "".join(f"{gene}\t1\t1\t10\t+\t10\t{count}\n" for gene, count in values),
        encoding="utf-8",
    )


def test_featurecounts_policy_translates_layout_and_strand():
    assert hisat2_strand_option("forward", "paired_end") == "FR"
    assert hisat2_strand_option("reverse", "paired_end") == "RF"
    assert hisat2_strand_option("forward", "single_end") == "F"
    assert hisat2_strand_option("reverse", "single_end") == "R"
    assert hisat2_strand_option("unstranded", "single_end") is None
    paired = featurecounts_arguments(layout="paired_end", strandedness="reverse")
    assert paired == ["-t", "exon", "-g", "gene_id", "-s", "2", "-Q", "0", "--primary", "-p", "--countReadPairs", "-B", "-C"]
    assert featurecounts_arguments(layout="single_end", strandedness="unstranded") == ["-t", "exon", "-g", "gene_id", "-s", "0", "-Q", "0", "--primary"]


def test_hisat2_is_pinned_to_2_2_3_for_builder_and_runtime():
    root = Path(__file__).parents[1]
    assert HISAT2_VERSION == "2.2.3"
    assert HISAT2_IMAGE == "quay.io/biocontainers/hisat2:2.2.3--h8471819_0"
    assert "hisat2=2.2.3" in (root / "environment.reference-builder.yml").read_text(encoding="utf-8")
    assert "hisat2:2.2.3--h8471819_0" in (root / "workflow" / "hisat2_featurecounts.nf").read_text(encoding="utf-8")


def test_hisat2_workflow_publishes_fastp_reports_and_feeds_upstream_reports_to_multiqc():
    text = (Path(__file__).parents[1] / "workflow" / "hisat2_featurecounts.nf").read_text(encoding="utf-8")
    assert "path 'fastp_*.html', optional: true, emit: fastp_html" in text
    assert "path 'fastp_*.json', optional: true, emit: fastp_json" in text
    assert "-h fastp_${r1}.html -j fastp_${r1}.json" in text
    assert "-I ${reads[1]} --detect_adapter_for_pe -o prepared/${r1}.fastq.gz" in text
    assert ".mix(FASTP_PREPARE.out.fastp_json)" in text
    assert ".mix(HISAT2_ALIGN.out.aligned.map { sample, strandedness, sam, summary -> summary })" in text
    assert "path 'raw_*_fastqc.zip', emit: zip" in text
    assert "path 'processed_*_fastqc.zip', emit: zip" in text
    assert 'mv "$f" "raw_$f"' in text and 'mv "$f" "processed_$f"' in text
    for channel in (
        "FASTQC_RAW.out.html", "FASTQC_RAW.out.zip",
        "FASTQC_PROCESSED.out.html", "FASTQC_PROCESSED.out.zip",
    ):
        assert f".mix({channel})" in text


def test_runtime_splice_argument_is_strategy_gated_in_the_hisat2_workflow():
    text = (Path(__file__).parents[1] / "workflow" / "hisat2_featurecounts.nf").read_text(encoding="utf-8")
    assert "params.hisat2_use_runtime_splices = false" in text
    assert "--known-splicesite-infile ${known_splices}" in text
    assert "-x ${index}/${params.hisat2_index_basename}" in text
    assert "HISAT2_ALIGN(FASTP_PREPARE.out.prepared, Channel.value(file(params.hisat2_index)), Channel.value(knownSplices))" in text


def test_hisat2_workflow_binds_supported_tool_threads_and_local_resources():
    text = (Path(__file__).parents[1] / "workflow" / "hisat2_featurecounts.nf").read_text(encoding="utf-8")
    for process, cpus, memory in (
        ("FASTP_PREPARE", 4, "6 GB"), ("HISAT2_ALIGN", 4, "6 GB"),
        ("SORT_LANE_BAM", 4, "6 GB"), ("MERGE_AND_INDEX", 4, "6 GB"),
        ("PREPARE_COUNT_BAM", 2, "3 GB"), ("FEATURECOUNTS", 4, "6 GB"),
        ("ASSEMBLE_COUNTS", 1, "2 GB"), ("MULTIQC", 1, "2 GB"),
    ):
        section = text.split(f"process {process} {{", 1)[1].split("\nprocess ", 1)[0]
        assert f"cpus {cpus}" in section and f"memory '{memory}'" in section
        assert "maxForks" not in section
    assert "fastp --thread ${task.cpus}" in text
    assert "hisat2 -p ${task.cpus}" in text
    assert "samtools sort -@ ${task.cpus}" in text
    assert "samtools merge -@ ${task.cpus}" in text
    assert "samtools view -@ ${task.cpus}" in text
    assert "samtools index -@ ${task.cpus}" in text
    assert "featureCounts -T ${task.cpus}" in text


def test_featurecounts_assembly_uses_declared_sample_ids_and_rejects_mismatched_genes(tmp_path: Path):
    c1, t1 = tmp_path / "lane_a.txt", tmp_path / "lane_b.txt"
    _featurecounts(c1, [("GeneB", 1), ("GeneA", 2)])
    _featurecounts(t1, [("GeneB", 10), ("GeneA", 20)])
    output = tmp_path / "canonical.csv"
    assemble_count_matrix({"T1": t1, "C1": c1}, output)
    assert output.read_text(encoding="utf-8") == "gene_id,C1,T1\nGeneB,1,10\nGeneA,2,20\n"
    _featurecounts(t1, [("GeneC", 10)])
    with pytest.raises(ValueError, match="gene set"):
        assemble_count_matrix({"C1": c1, "T1": t1}, output)


def test_hisat2_backend_requires_explicit_strandedness(project_factory):
    config = deepcopy(base_config())
    config.update({
        "schema_version": "1.2",
        "input": {"type": "fastq", "path": "input/fastq", "layout": "paired_end", "preprocessing": "raw"},
        "upstream": {"engine": "nfcore_rnaseq", "pipeline_version": "3.26.0", "strandedness": "auto", "quantification": {"method": "hisat2_featurecounts"}},
    })
    report = validate_project(project_factory(config=config))
    assert not report.is_valid
    assert any("explicit upstream.strandedness" in issue.message for issue in report.errors)


def test_hisat2_command_uses_first_party_workflow_and_custom_index(tmp_path: Path):
    root = tmp_path / "project"
    fastq = root / "input" / "fastq"
    fastq.mkdir(parents=True)
    for name in ("C1_R1.fastq.gz", "C1_R2.fastq.gz", "T1_R1.fastq.gz", "T1_R2.fastq.gz"):
        (fastq / name).write_bytes(b"x")
    (root / "metadata.csv").write_text("sample_id,condition\nC1,C\nT1,T\n", encoding="utf-8")
    (root / "contrasts.csv").write_text("contrast_id,factor,numerator,denominator\nT_vs_C,condition,T,C\n", encoding="utf-8")
    reference = root / "reference"
    reference.mkdir()
    (reference / "genome.fa").write_text(">1\nACGT\n", encoding="utf-8")
    (reference / "genes.gtf").write_text("1\tt\texon\t1\t4\t.\t+\t.\tgene_id \"GeneA\";\n", encoding="utf-8")
    index = reference / "index"; index.mkdir()
    for number in range(1, 9):
        (index / f"genome.{number}.ht2").write_bytes(b"x")
    config = deepcopy(base_config())
    config.update({
        "schema_version": "1.2", "input": {"type": "fastq", "path": "input/fastq", "layout": "paired_end", "preprocessing": "pretrimmed"},
        "upstream": {"engine": "nfcore_rnaseq", "pipeline_version": "3.26.0", "strandedness": "reverse", "quantification": {"method": "hisat2_featurecounts"}},
        "reference": {"source": "custom", "fasta": "reference/genome.fa", "gtf": "reference/genes.gtf", "hisat2_index": "reference/index"},
        "analysis": {"enrichment": []},
    })
    (root / "project.yaml").write_text(yaml.safe_dump(config, sort_keys=False), encoding="utf-8")
    report = validate_project(root)
    assert report.is_valid and report.execution_ready
    local_config = root / "local.config"
    local_config.write_text("process { resourceLimits = [cpus: 8, memory: '12.GB'] }\n", encoding="utf-8")
    command = build_hisat2_featurecounts_command(report, samplesheet=root / "samples.csv", output_dir=root / "out", profile="local", config_file=local_config)
    assert any(item.endswith("workflow/hisat2_featurecounts.nf") for item in command)
    assert command[command.index("-c") + 1] == str(local_config.resolve())
    assert command[command.index("--strandedness") + 1] == "reverse"
    assert command[command.index("--hisat2_index") + 1] == str(index)
    resolved = resolved_upstream_implementation(report)
    assert resolved["implementation"]["name"] == "nf-rna/hisat2_featurecounts"
    assert resolved["implementation"]["workflow_path"] == "workflow/hisat2_featurecounts.nf"
    assert len(resolved["implementation"]["workflow_sha256"]) == 64
    assert resolved["legacy_fastq_config"]["engine"] == "nfcore_rnaseq"
    run_dir = root / "runs" / "CASE" / "20260905-120000+0800"
    (run_dir / "frozen").mkdir(parents=True)
    frozen = freeze_case_inputs(report, CaseRun("CASE", "20260905-120000+0800", run_dir, "2026-09-05T12:00:00+08:00"), profile="local", command=["rnaseq", "run"])
    execution = yaml.safe_load((run_dir / "frozen" / "execution_manifest.yaml").read_text(encoding="utf-8"))
    assert execution["upstream_implementation"]["implementation"]["name"] == "nf-rna/hisat2_featurecounts"
    assert execution["upstream_implementation"]["legacy_fastq_config"]["pipeline_version"] == "3.26.0"
    assert yaml.safe_load(frozen.manifest.read_text(encoding="utf-8"))["upstream"]["quantification"]["method"] == "hisat2_featurecounts"
