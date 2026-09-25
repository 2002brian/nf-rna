"""The Salmon transcriptome must share the GTF transcript_id namespace exactly.

nf-core/rnaseq 3.26.0 CUSTOM_TX2GENE matches Salmon transcript names against GTF
attributes by exact equality.  An index built from Ensembl cDNA (versioned
``ENSMUST00000200568.2``) against an Ensembl GTF (``transcript_id
"ENSMUST00000200568"; transcript_version "2"``) matches nothing, so such a
reference must be rejected before any nf-core run.
"""

from __future__ import annotations

import hashlib
import json
from datetime import UTC, datetime
from pathlib import Path

import pytest
import yaml

from rnaseq.errors import ExecutionPreflightError
from rnaseq.planner import generate_plan
from rnaseq.references import (
    GTF_DERIVED_TRANSCRIPTOME_STRATEGY,
    LocalReferenceError,
    ReferenceAdoptionError,
    ReferencePreparationError,
    _load_local_reference_root,
    adopt_local_salmon_index,
    filter_gtf_like_nfcore,
    load_local_reference_root,
    prepare_local_reference,
    transcript_id_contract,
    transcript_id_contract_summary,
)
from rnaseq.service import prepare_service_run
from rnaseq.validators import validate_project
from test_local_reference import _fake_builder, _local_fastq_project, _write_salmon_index


def _sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


GENOME = ">1 dna:chromosome chromosome:GRCm39:1:1:40:1 REF\n" + "ACGT" * 10 + "\n"
ENSEMBL_GTF = "".join([
    "#!genome-build GRCm39\n",
    '1\tensembl\tgene\t1\t40\t.\t+\t.\tgene_id "ENSMUSG00000000001"; gene_version "5";\n',
    '1\tensembl\ttranscript\t1\t40\t.\t+\t.\tgene_id "ENSMUSG00000000001"; gene_version "5"; transcript_id "ENSMUST00000000001"; transcript_version "3";\n',
    '1\tensembl\texon\t1\t40\t.\t+\t.\tgene_id "ENSMUSG00000000001"; gene_version "5"; transcript_id "ENSMUST00000000001"; transcript_version "3";\n',
    '1\tensembl\ttranscript\t1\t20\t.\t+\t.\tgene_id "ENSMUSG00000000001"; gene_version "5"; transcript_id "ENSMUST00000200568"; transcript_version "2";\n',
    '1\tensembl\texon\t1\t20\t.\t+\t.\tgene_id "ENSMUSG00000000001"; gene_version "5"; transcript_id "ENSMUST00000200568"; transcript_version "2";\n',
    'JH584295.1\tensembl\texon\t1\t10\t.\t+\t.\tgene_id "ENSMUSG00000000009"; transcript_id "ENSMUST00000000009"; transcript_version "1";\n',
])
VERSIONED_CDNA = ">ENSMUST00000000001.3 cdna chromosome:GRCm39:1:1:40:1 gene:ENSMUSG00000000001.5\n" + "ACGT" * 10 + "\n" \
    + ">ENSMUST00000200568.2 cdna chromosome:GRCm39:1:1:20:1 gene:ENSMUSG00000000001.5\n" + "ACGT" * 5 + "\n"
GTF_DERIVED = ">ENSMUST00000000001\n" + "ACGT" * 10 + "\n>ENSMUST00000200568\n" + "ACGT" * 5 + "\n"


def _ensembl_reference(root: Path) -> Path:
    """Rewrite a fixture reference with Ensembl conventions and a registered versioned cDNA."""

    for relative, content in (
        ("fasta/Homo_sapiens.GRCh38.dna.primary_assembly.fa", GENOME),
        ("annotation/Homo_sapiens.GRCh38.116.gtf", ENSEMBL_GTF),
        ("salmon/Homo_sapiens.GRCh38.cdna.all.fa", VERSIONED_CDNA),
    ):
        (root / relative).write_text(content, encoding="utf-8")
    manifest_path = root / "reference_manifest.yaml"
    manifest = yaml.safe_load(manifest_path.read_text(encoding="utf-8"))
    for key in ("genome_fasta", "annotation_gtf", "transcript_fasta"):
        manifest["files"][key]["sha256"] = _sha(root / manifest["files"][key]["path"])
    manifest_path.write_text(yaml.safe_dump(manifest, sort_keys=False), encoding="utf-8")
    return manifest_path


def _declare_cdna_prebuilt_index(root: Path) -> Path:
    """The forensic 4T1 shape: a decoy-aware index declared from Ensembl cDNA, without a contract."""

    manifest_path = _ensembl_reference(root)
    index = root / "salmon" / "index"
    _write_salmon_index(index, num_decoys=1)
    manifest = yaml.safe_load(manifest_path.read_text(encoding="utf-8"))
    manifest["salmon"] = {
        "index": "salmon/index", "status": "built", "strategy": "decoy_aware", "version": "1.10.3",
        "source_transcriptome_sha256": manifest["files"]["transcript_fasta"]["sha256"],
        "source_genome_sha256": manifest["files"]["genome_fasta"]["sha256"],
        "provenance": {"transcriptome_strategy": "manifest_registered_transcript_fasta"},
    }
    manifest_path.write_text(yaml.safe_dump(manifest, sort_keys=False), encoding="utf-8")
    return index


def _contract(root: Path, transcriptome: Path) -> dict:
    reference, _ = _load_local_reference_root(root, "reference_manifest.yaml", None, validate_salmon=False)
    return transcript_id_contract(
        identity={"species": reference.species, "provider": reference.provider, "release": reference.release,
                  "assembly": reference.assembly, "assembly_patch": reference.assembly_patch},
        genome_fasta=reference.genome_fasta, annotation_gtf=reference.annotation_gtf,
        transcriptome=transcriptome, transcriptome_relative_path=transcriptome.relative_to(root).as_posix(),
    )


# ---------------------------------------------------------------- A / B: the contract itself


def test_versioned_ensembl_cdna_ids_are_rejected_against_an_unversioned_gtf(tmp_path):
    _root, reference = _local_fastq_project(tmp_path)
    _ensembl_reference(reference)
    payload = _contract(reference, reference / "salmon" / "Homo_sapiens.GRCh38.cdna.all.fa")
    assert payload["status"] == "BLOCKED"
    compatibility = payload["identifier_compatibility"]
    assert compatibility["intersection"] == 0
    assert compatibility["fasta_only"] == 2
    # Diagnosed, but never normalized to make the reference pass.
    assert compatibility["fasta_only_matching_gtf_after_version_strip"] == 2
    assert payload["contract"]["normalization"] == "none"
    summary = transcript_id_contract_summary(payload)
    assert "0 of 2 transcriptome IDs exactly equal a GTF transcript_id" in summary
    assert "'.N' version suffix" in summary


def test_gtf_derived_transcriptome_with_exact_transcript_ids_is_accepted(tmp_path):
    _root, reference = _local_fastq_project(tmp_path)
    _ensembl_reference(reference)
    derived = reference / "derived.fa"
    derived.write_text(GTF_DERIVED, encoding="utf-8")
    payload = _contract(reference, derived)
    assert payload["status"] == "PASS"
    assert payload["identifier_compatibility"]["intersection"] == 2
    assert payload["identifier_compatibility"]["fasta_only"] == 0
    # The GTF transcript on a scaffold absent from the genome is not in the transcriptome.
    assert payload["identifier_compatibility"]["gtf_only"] == 1
    assert payload["identifier_compatibility"]["gtf_on_genome_only"] == 0
    assert payload["tx2gene"] == {**payload["tx2gene"], "one_gene_mappings": 2, "zero_gene_mappings": 0, "multi_gene_mappings": 0}


@pytest.mark.parametrize(
    ("transcriptome", "field"),
    [
        (GTF_DERIVED + ">ENSMUST00000000001\nACGT\n", ("fasta_validation", "duplicate_transcript_ids")),
        (GTF_DERIVED + ">ENSMUST99999999999\nACGT\n", ("identifier_compatibility", "fasta_only")),
    ],
    ids=["duplicate-id", "transcript-absent-from-gtf"],
)
def test_contract_blocks_duplicate_or_unannotated_transcripts(tmp_path, transcriptome, field):
    _root, reference = _local_fastq_project(tmp_path)
    _ensembl_reference(reference)
    derived = reference / "derived.fa"
    derived.write_text(transcriptome, encoding="utf-8")
    payload = _contract(reference, derived)
    assert payload["status"] == "BLOCKED"
    assert payload[field[0]][field[1]] == 1


def test_partial_transcriptome_with_matching_ids_is_rejected_for_scope(tmp_path):
    """Matching IDs are not enough: a cDNA-style subset (no lncRNA) must not pass."""

    _root, reference = _local_fastq_project(tmp_path)
    _ensembl_reference(reference)
    derived = reference / "subset.fa"
    derived.write_text(">ENSMUST00000000001\n" + "ACGT" * 10 + "\n", encoding="utf-8")
    payload = _contract(reference, derived)
    assert payload["status"] == "BLOCKED"
    assert payload["identifier_compatibility"]["fasta_only"] == 0
    assert payload["identifier_compatibility"]["gtf_on_genome_only"] == 1
    assert payload["identifier_compatibility"]["gtf_on_genome_only_examples"] == ["ENSMUST00000200568"]


def test_gtf_filter_reproduces_nfcore_custom_gtffilter(tmp_path):
    genome = tmp_path / "genome.fa"
    gtf = tmp_path / "genes.gtf"
    genome.write_text(GENOME, encoding="utf-8")
    gtf.write_text(ENSEMBL_GTF, encoding="utf-8")
    stats = filter_gtf_like_nfcore(genome, gtf, tmp_path / "filtered.gtf")
    kept = (tmp_path / "filtered.gtf").read_text(encoding="utf-8").splitlines(keepends=True)
    # Comments, gene records (no transcript_id) and records off the genome FASTA are dropped;
    # every kept line is byte-identical to the input.
    assert kept == ENSEMBL_GTF.splitlines(keepends=True)[2:6]
    assert (stats["lines_kept"], stats["lines_removed"]) == (4, 3)


# ---------------------------------------------------------------- C: the production builder


def test_builder_derives_the_transcriptome_from_genome_and_gtf_not_registered_cdna(tmp_path):
    _root, reference = _local_fastq_project(tmp_path)
    _ensembl_reference(reference)
    commands: list[list[str]] = []
    resolver, run = _fake_builder(tmp_path, transcripts=GTF_DERIVED, commands=commands)
    prepared = prepare_local_reference(reference, runner=run, tool_resolver=resolver)
    assert prepared.transcriptome_strategy == GTF_DERIVED_TRANSCRIPTOME_STRATEGY
    assert not prepared.external_transcript_fasta_used
    (rsem,) = [command for command in commands if Path(command[0]).name == "rsem-prepare-reference"]
    assert rsem[1] == "--gtf" and rsem[-2] == str(prepared.genome_fasta.path)
    cdna = str(reference / "salmon" / "Homo_sapiens.GRCh38.cdna.all.fa")
    assert not any(cdna in argument for command in commands for argument in command)
    assert prepared.salmon_transcriptome.path.read_text(encoding="utf-8") == GTF_DERIVED


def test_builder_refuses_to_index_a_transcriptome_that_breaks_the_contract(tmp_path):
    _root, reference = _local_fastq_project(tmp_path)
    _ensembl_reference(reference)
    before = (reference / "reference_manifest.yaml").read_bytes()
    commands: list[list[str]] = []
    resolver, run = _fake_builder(tmp_path, transcripts=VERSIONED_CDNA, commands=commands)
    with pytest.raises(ReferencePreparationError, match="transcript-ID contract BLOCKED.*version suffix"):
        prepare_local_reference(reference, runner=run, tool_resolver=resolver)
    assert (reference / "reference_manifest.yaml").read_bytes() == before
    assert not (reference / "salmon" / "gtf_derived").exists()
    assert not any(Path(command[0]).name == "salmon" and "index" in command for command in commands)
    (blocked,) = reference.glob(".rnaseq-reference-prepare-*/work/transcript_id_contract.BLOCKED.json")
    assert json.loads(blocked.read_text(encoding="utf-8"))["status"] == "BLOCKED"


def test_builder_requires_the_nfcore_rsem_package(tmp_path):
    _root, reference = _local_fastq_project(tmp_path)
    resolver, run = _fake_builder(tmp_path, rsem_version="1.3.1")
    with pytest.raises(ReferencePreparationError, match="bioconda rsem=1.3.3"):
        prepare_local_reference(reference, runner=run, tool_resolver=resolver)
    assert not list(reference.glob(".rnaseq-reference-prepare-*"))


# ---------------------------------------------------------------- D: rejected before nf-core execution


def test_cdna_prebuilt_reference_is_rejected_by_validate_and_run_preflight(tmp_path):
    root, reference = _local_fastq_project(tmp_path)
    _declare_cdna_prebuilt_index(reference)
    with pytest.raises(LocalReferenceError, match="without a transcript-ID contract.*--rebuild-salmon"):
        load_local_reference_root(reference)
    report = validate_project(root)
    assert not report.is_valid and not report.execution_ready
    assert any("transcript-ID contract" in issue.message for issue in report.errors)
    with pytest.raises(ExecutionPreflightError, match="validation failed"):
        prepare_service_run(report, profile="local")


def test_a_blocked_contract_artifact_does_not_make_the_reference_ready(tmp_path):
    root, reference = _local_fastq_project(tmp_path)
    _declare_cdna_prebuilt_index(reference)
    artifact = reference / "salmon" / "index.transcript_id_contract.json"
    artifact.write_text(json.dumps(_contract(reference, reference / "salmon" / "Homo_sapiens.GRCh38.cdna.all.fa")), encoding="utf-8")
    manifest_path = reference / "reference_manifest.yaml"
    manifest = yaml.safe_load(manifest_path.read_text(encoding="utf-8"))
    manifest["salmon"]["validation"] = {"artifact": "salmon/index.transcript_id_contract.json", "transcript_id_contract": "exact"}
    manifest_path.write_text(yaml.safe_dump(manifest, sort_keys=False), encoding="utf-8")
    report = validate_project(root)
    assert not report.is_valid
    assert any("salmon validation artifact status" in issue.message for issue in report.errors)


def test_adopting_a_versioned_cdna_index_is_refused_without_writing_anything(tmp_path):
    _root, reference = _local_fastq_project(tmp_path)
    _ensembl_reference(reference)
    _write_salmon_index(reference / "external" / "salmon", num_decoys=1)
    before = (reference / "reference_manifest.yaml").read_bytes()
    with pytest.raises(ReferenceAdoptionError, match="transcript_id namespace.*version suffix"):
        adopt_local_salmon_index(reference, index="external/salmon", strategy="decoy_aware")
    assert (reference / "reference_manifest.yaml").read_bytes() == before
    assert not (reference / "external" / "salmon.transcript_id_contract.json").exists()


# ---------------------------------------------------------------- rebuild preserves the evidence


def test_rebuild_replaces_a_rejected_cdna_index_and_preserves_the_old_identity(tmp_path):
    root, reference = _local_fastq_project(tmp_path)
    old_index = _declare_cdna_prebuilt_index(reference)
    old_manifest = (reference / "reference_manifest.yaml").read_bytes()
    old_index_files = {path.name: path.read_bytes() for path in old_index.iterdir()}
    resolver, run = _fake_builder(tmp_path, transcripts=GTF_DERIVED)

    # Without the explicit flag the rejected declaration stops preparation, with instructions.
    with pytest.raises(LocalReferenceError, match="--rebuild-salmon"):
        prepare_local_reference(reference, runner=run, tool_resolver=resolver)
    rebuilt = prepare_local_reference(
        reference, rebuild_salmon=True, runner=run, tool_resolver=resolver, now=lambda: datetime(2026, 9, 25, tzinfo=UTC),
    )

    archive = reference / f"reference_manifest.{hashlib.sha256(old_manifest).hexdigest()[:12]}.yaml"
    assert archive.read_bytes() == old_manifest
    assert {path.name: path.read_bytes() for path in old_index.iterdir()} == old_index_files
    assert rebuilt.salmon_index == reference / "salmon" / "gtf_derived" / "index"
    replaces = yaml.safe_load((reference / "reference_manifest.yaml").read_text(encoding="utf-8"))["salmon"]["provenance"]["replaces"]
    assert replaces["manifest_archive"] == archive.name
    assert replaces["index"] == "salmon/index"
    assert replaces["reason"].startswith("rejected Salmon declaration:")
    report = validate_project(root)
    assert report.is_valid and report.execution_ready, report.errors
    generate_plan(report)
    plan = (root / "planning" / "analysis_plan.md").read_text(encoding="utf-8")
    assert GTF_DERIVED_TRANSCRIPTOME_STRATEGY in plan
    assert "Salmon index built from it: `false`" in plan
