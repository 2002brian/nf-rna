from __future__ import annotations

import hashlib
import json
import subprocess
from copy import deepcopy
from datetime import UTC, datetime
from pathlib import Path

import pytest
import yaml
from typer.testing import CliRunner

from conftest import BASE_CONTRASTS, base_config
from rnaseq.errors import ExecutionPreflightError
from rnaseq.execution import RuntimeCheck, _validate_custom_reference_files
from rnaseq.execution import build_hisat2_featurecounts_command, build_nextflow_command
from rnaseq.cli import app
from rnaseq.planner import generate_plan
from rnaseq.references import (
    ADOPTED_EXISTING_INDEX,
    DECOY_STRATEGY,
    REQUIRED_SALMON_INDEX_FILES,
    SALMON_BUILT,
    SALMON_KMER_SIZE,
    SALMON_STRATEGY_DECOY_AWARE,
    SALMON_STRATEGY_TRANSCRIPTOME_ONLY,
    SALMON_VERSION,
    HISAT2_VERSION,
    LocalReferenceError,
    TRANSCRIPTOME_STRATEGY,
    ReferenceAdoptionError,
    ReferencePreparationError,
    adopt_local_salmon_index,
    compatible_registered_references,
    load_local_reference_root,
    prepare_local_reference,
    prepare_local_hisat2_reference,
    register_local_reference,
)
from rnaseq.service import CaseRun, create_case_run, freeze_case_inputs, prepare_service_run
from rnaseq.validators import validate_project


runner = CliRunner()


def _sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _write_reference(root: Path, *, species: str = "Homo sapiens") -> None:
    genome = root / "fasta" / "Homo_sapiens.GRCh38.dna.primary_assembly.fa"
    gtf = root / "annotation" / "Homo_sapiens.GRCh38.116.gtf"
    transcript = root / "salmon" / "Homo_sapiens.GRCh38.cdna.all.fa"
    for path, content in (
        (genome, ">chr1\nACGT\n"),
        (gtf, "chr1\tfixture\texon\t1\t4\t.\t+\t.\tgene_id \"GENE1\"; transcript_id \"TX1\";\n"),
        (transcript, ">TX1 gene:GENE1\nACGT\n"),
    ):
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(content, encoding="utf-8")
    manifest = {
        "schema_version": "1.0",
        "reference": {
            "species": species,
            "provider": "Ensembl",
            "release": 116,
            "assembly": "GRCh38",
            "assembly_patch": "p14",
        },
        "files": {
            "genome_fasta": {"path": str(genome.relative_to(root)), "sha256": _sha(genome)},
            "annotation_gtf": {"path": str(gtf.relative_to(root)), "sha256": _sha(gtf)},
            "transcript_fasta": {"path": str(transcript.relative_to(root)), "sha256": _sha(transcript)},
        },
        "salmon": {"index": None, "status": "not_built"},
    }
    (root / "reference_manifest.yaml").write_text(yaml.safe_dump(manifest, sort_keys=False), encoding="utf-8")


def _write_salmon_index(index: Path, *, num_decoys: int = 0) -> None:
    index.mkdir(parents=True)
    for name in REQUIRED_SALMON_INDEX_FILES:
        if name == "info.json":
            (index / name).write_text(json.dumps({
                "index_version": 4,
                "k": 31,
                "num_decoys": num_decoys,
                "SeqHash": "a" * 64,
                "NameHash": "b" * 64,
                "seq_length": 123,
                "num_kmers": 456,
                "num_contigs": 789,
            }), encoding="utf-8")
        elif name == "versionInfo.json":
            (index / name).write_text(json.dumps({"salmonVersion": "1.10.3", "indexVersion": 5}), encoding="utf-8")
        else:
            (index / name).write_text("fixture\n", encoding="utf-8")


def _write_hisat2_index(index: Path) -> None:
    index.mkdir(parents=True, exist_ok=True)
    for number in range(1, 9):
        (index / f"genome.{number}.ht2").write_text("fixture\n", encoding="utf-8")


def _write_gtf_derived_candidate(reference: Path) -> tuple[Path, Path, Path]:
    manifest_path = reference / "reference_manifest.yaml"
    manifest = yaml.safe_load(manifest_path.read_text(encoding="utf-8"))
    transcriptome = reference / "salmon" / "gtf_derived" / "transcripts.fa"
    transcriptome.parent.mkdir(parents=True, exist_ok=True)
    transcriptome.write_text(">TX1\nACGT\n", encoding="utf-8")
    index = transcriptome.parent / "transcriptome_only_index"
    _write_salmon_index(index)
    artifact = transcriptome.parent / "transcriptome_validation.json"
    artifact.write_text(json.dumps({
        "reference": manifest["reference"],
        "inputs": {
            "genome": manifest["files"]["genome_fasta"],
            "gtf": manifest["files"]["annotation_gtf"],
        },
        "generated_transcriptome": {
            "path": str(transcriptome.relative_to(reference)),
            "sha256": _sha(transcriptome),
            "transcript_count": 1,
            "total_sequence_length": 4,
        },
        "identifier_compatibility": {
            "fasta_unique_ids": 1,
            "gtf_unique_transcript_ids": 1,
            "intersection": 1,
            "fasta_only": 0,
            "gtf_only": 0,
            "fasta_mapping_rate": 1.0,
            "gtf_mapping_rate": 1.0,
        },
        "tx2gene": {"one_gene_mappings": 1, "zero_gene_mappings": 0, "multi_gene_mappings": 0},
        "fasta_validation": {"duplicate_transcript_ids": 0},
        "status": "PASS",
    }), encoding="utf-8")
    return index, transcriptome, artifact


def _adopt_gtf_derived_candidate(reference: Path) -> Path:
    index, transcriptome, artifact = _write_gtf_derived_candidate(reference)
    adopted = adopt_local_salmon_index(
        reference,
        index=str(index.relative_to(reference)),
        transcriptome=str(transcriptome.relative_to(reference)),
        strategy=SALMON_STRATEGY_TRANSCRIPTOME_ONLY,
        validation_artifact=str(artifact.relative_to(reference)),
        now=lambda: datetime(2026, 8, 31, tzinfo=UTC),
    )
    assert adopted.salmon_index == index
    return index


def _mark_salmon_built(reference: Path) -> Path:
    index = reference / "salmon" / "index"
    _write_salmon_index(index)
    manifest_path = reference / "reference_manifest.yaml"
    manifest = yaml.safe_load(manifest_path.read_text(encoding="utf-8"))
    manifest["salmon"] = {
        "index": "salmon/index",
        "status": SALMON_BUILT,
        "provenance": {
            "salmon_version": SALMON_VERSION,
            "index_path": "salmon/index",
            "genome_fasta_sha256": manifest["files"]["genome_fasta"]["sha256"],
            "annotation_gtf_sha256": manifest["files"]["annotation_gtf"]["sha256"],
            "transcriptome_strategy": TRANSCRIPTOME_STRATEGY,
            "decoy_strategy": DECOY_STRATEGY,
            "kmer_size": SALMON_KMER_SIZE,
            "built_at": "2026-08-31T00:00:00+00:00",
        },
    }
    manifest_path.write_text(yaml.safe_dump(manifest, sort_keys=False), encoding="utf-8")
    return index


def _mark_genome_only_hisat2(reference: Path, *, runtime_compatibility: str = "validated") -> Path:
    index = reference / "hisat2" / "index"
    _write_hisat2_index(index)
    splice = reference / "hisat2" / "splice_sites.txt"
    splice.parent.mkdir(parents=True, exist_ok=True)
    splice.write_text("chr1\t1\t4\t+\n", encoding="utf-8")
    manifest_path = reference / "reference_manifest.yaml"
    manifest = yaml.safe_load(manifest_path.read_text(encoding="utf-8"))
    manifest["schema_version"] = "1.2"
    manifest["purpose"] = "production"
    manifest["hisat2"] = {
        "status": "built", "index": "hisat2/index", "index_prefix": "hisat2/index/genome",
        "strategy": "genome_only_runtime_splices",
        "splice_sites": {"path": "hisat2/splice_sites.txt", "sha256": _sha(splice)},
        "provenance": {
            "index_builder_version": "2.2.3", "runtime_aligner_version": "2.2.3",
            "genome_fasta_sha256": manifest["files"]["genome_fasta"]["sha256"],
            "annotation_gtf_sha256": manifest["files"]["annotation_gtf"]["sha256"],
            "splice_sites_derived_from_gtf_sha256": manifest["files"]["annotation_gtf"]["sha256"],
            "runtime_compatibility": runtime_compatibility,
            "index_validation": "complete numbered .ht2 family", "built_at": "2026-09-08T00:00:00+00:00",
        },
    }
    manifest_path.write_text(yaml.safe_dump(manifest, sort_keys=False), encoding="utf-8")
    return index


def _local_fastq_project(tmp_path: Path, *, species: str = "Homo sapiens") -> tuple[Path, Path]:
    root = tmp_path / "project"
    fastq = root / "input" / "fastq"
    fastq.mkdir(parents=True)
    rows = ["sample_id,condition"]
    for condition, prefix in (("Control", "C"), ("Treatment", "T")):
        for number in range(1, 4):
            sample = f"{prefix}{number}"
            rows.append(f"{sample},{condition}")
            for read in ("R1", "R2"):
                (fastq / f"{sample}_{read}.fastq.gz").write_bytes(b"fixture")
    (root / "metadata.csv").write_text("\n".join(rows) + "\n", encoding="utf-8")
    (root / "contrasts.csv").write_text(BASE_CONTRASTS, encoding="utf-8")
    reference = tmp_path / "reference"
    _write_reference(reference, species=species)
    config = deepcopy(base_config())
    config["organism"] = {"species": species}
    config["input"] = {"type": "fastq", "path": "input/fastq", "layout": "paired_end"}
    config["upstream"] = {
        "engine": "nfcore_rnaseq",
        "pipeline_version": "3.26.0",
        "aligner": None,
        "strandedness": "auto",
        "quantification": {"method": "salmon"},
    }
    config["reference"] = {"source": "local", "root": str(reference), "manifest": "reference_manifest.yaml"}
    (root / "project.yaml").write_text(yaml.safe_dump(config, sort_keys=False), encoding="utf-8")
    return root, reference


def _promote_reference_for_production(root: Path, reference: Path, *, purpose: str = "production") -> None:
    _mark_salmon_built(reference)
    manifest_path = reference / "reference_manifest.yaml"
    manifest = yaml.safe_load(manifest_path.read_text(encoding="utf-8"))
    manifest["schema_version"] = "1.1"
    manifest["purpose"] = purpose
    manifest["sources"] = {
        name: {
            "accession": f"fixture:{name}",
            "upstream_checksum": {"algorithm": "sha256", "value": item["sha256"]},
        }
        for name, item in manifest["files"].items()
    }
    manifest_path.write_text(yaml.safe_dump(manifest, sort_keys=False), encoding="utf-8")
    config_path = root / "project.yaml"
    config = yaml.safe_load(config_path.read_text(encoding="utf-8"))
    config["reference"]["acceptance"] = "production"
    config["runtime"] = {"control_plane_image": "rnaseq-control-plane:0.5.1"}
    config_path.write_text(yaml.safe_dump(config, sort_keys=False), encoding="utf-8")


def test_production_acceptance_requires_reviewed_managed_reference_and_verified_sources(tmp_path):
    root, reference = _local_fastq_project(tmp_path)
    _promote_reference_for_production(root, reference)
    report = validate_project(root)
    assert report.is_valid, report.errors
    assert report.local_reference.purpose == "production"
    assert report.local_reference.source_provenance["genome_fasta"]["upstream_checksum"]["verification"] == "matched_local_asset"


def test_machine_local_registry_deduplicates_and_resolves_by_backend(tmp_path):
    root, reference = _local_fastq_project(tmp_path)
    _promote_reference_for_production(root, reference)
    registry = tmp_path / "config" / "nf-rna" / "references.yaml"

    register_local_reference(reference, registry_path=registry)
    register_local_reference(reference, registry_path=registry)

    payload = yaml.safe_load(registry.read_text(encoding="utf-8"))
    assert payload["schema_version"] == "1.0"
    assert len(payload["references"]) == 1
    assert "files" not in payload["references"][0]
    assert [item.root for item in compatible_registered_references(
        "Homo sapiens", "salmon", registry_path=registry
    )] == [reference]
    assert compatible_registered_references("Homo sapiens", "hisat2_featurecounts", registry_path=registry) == []

    _mark_genome_only_hisat2(reference)
    assert [item.root for item in compatible_registered_references(
        "Homo sapiens", "hisat2_featurecounts", registry_path=registry
    )] == [reference]

    other_root, other_reference = _local_fastq_project(tmp_path / "older_release")
    _promote_reference_for_production(other_root, other_reference)
    other_manifest_path = other_reference / "reference_manifest.yaml"
    other_manifest = yaml.safe_load(other_manifest_path.read_text(encoding="utf-8"))
    other_manifest["reference"]["release"] = 115
    other_manifest_path.write_text(yaml.safe_dump(other_manifest, sort_keys=False), encoding="utf-8")
    register_local_reference(other_reference, registry_path=registry)
    assert [item.release for item in compatible_registered_references(
        "Homo sapiens", "salmon", registry_path=registry
    )] == [116, 115]


def test_reference_register_command_uses_xdg_registry(monkeypatch, tmp_path):
    root, reference = _local_fastq_project(tmp_path)
    _promote_reference_for_production(root, reference)
    config_home = tmp_path / "config-home"
    monkeypatch.setenv("XDG_CONFIG_HOME", str(config_home))

    result = runner.invoke(app, ["reference", "register", str(reference)])

    assert result.exit_code == 0, result.output
    assert "Registered managed reference" in result.output
    assert (config_home / "nf-rna" / "references.yaml").is_file()


def test_assembly_identity_accepts_patched_and_unpatched_manifests_and_registry_matches(tmp_path):
    human_root, human_reference = _local_fastq_project(tmp_path / "human")
    _promote_reference_for_production(human_root, human_reference)
    human = load_local_reference_root(human_reference)
    assert human.assembly_patch == "p14"
    assert human.assembly_identity == "GRCh38.p14"

    mouse_root, mouse_reference = _local_fastq_project(tmp_path / "mouse", species="Mus musculus")
    _promote_reference_for_production(mouse_root, mouse_reference)
    mouse_manifest_path = mouse_reference / "reference_manifest.yaml"
    mouse_manifest = yaml.safe_load(mouse_manifest_path.read_text(encoding="utf-8"))
    mouse_manifest["reference"]["assembly"] = "GRCm39"
    mouse_manifest["reference"].pop("assembly_patch")
    mouse_manifest_path.write_text(yaml.safe_dump(mouse_manifest, sort_keys=False), encoding="utf-8")
    mouse = load_local_reference_root(mouse_reference)
    assert mouse.assembly_patch is None
    assert mouse.assembly_identity == "GRCm39"

    mouse_manifest = yaml.safe_load(mouse_manifest_path.read_text(encoding="utf-8"))
    mouse_manifest["reference"]["assembly_patch"] = None
    mouse_manifest_path.write_text(yaml.safe_dump(mouse_manifest, sort_keys=False), encoding="utf-8")
    assert load_local_reference_root(mouse_reference).assembly_identity == "GRCm39"
    mouse_manifest["reference"]["assembly_patch"] = ""
    mouse_manifest_path.write_text(yaml.safe_dump(mouse_manifest, sort_keys=False), encoding="utf-8")
    assert load_local_reference_root(mouse_reference).assembly_identity == "GRCm39"

    registry = tmp_path / "config" / "nf-rna" / "references.yaml"
    register_local_reference(human_reference, registry_path=registry)
    register_local_reference(mouse_reference, registry_path=registry)
    assert [item.assembly_identity for item in compatible_registered_references(
        "Homo sapiens", "salmon", registry_path=registry
    )] == ["GRCh38.p14"]
    assert [item.assembly_identity for item in compatible_registered_references(
        "Mus musculus", "salmon", registry_path=registry
    )] == ["GRCm39"]


def test_synthetic_manifest_works_normally_but_fails_production_acceptance(tmp_path):
    root, reference = _local_fastq_project(tmp_path)
    _promote_reference_for_production(root, reference, purpose="synthetic_test")
    production = validate_project(root)
    assert "synthetic_reference_not_production" in {issue.code for issue in production.errors}
    config_path = root / "project.yaml"
    config = yaml.safe_load(config_path.read_text(encoding="utf-8"))
    config["reference"]["acceptance"] = "standard"
    config_path.write_text(yaml.safe_dump(config, sort_keys=False), encoding="utf-8")
    assert validate_project(root).is_valid


def test_legacy_manifest_is_not_silently_promoted_to_production(tmp_path):
    root, reference = _local_fastq_project(tmp_path)
    _mark_salmon_built(reference)
    config_path = root / "project.yaml"
    config = yaml.safe_load(config_path.read_text(encoding="utf-8"))
    config["reference"]["acceptance"] = "production"
    config["runtime"] = {"control_plane_image": "rnaseq-control-plane:0.5.1"}
    config_path.write_text(yaml.safe_dump(config, sort_keys=False), encoding="utf-8")
    report = validate_project(root)
    assert "legacy_reference_not_production" in {issue.code for issue in report.errors}


def test_production_runtime_requires_observed_immutable_image_identity(monkeypatch, tmp_path, production_capable_execution_capacity):
    root, reference = _local_fastq_project(tmp_path)
    _promote_reference_for_production(root, reference)
    report = validate_project(root)
    generate_plan(report)
    monkeypatch.setattr("rnaseq.service.check_nextflow", lambda: RuntimeCheck("Nextflow", "FOUND", "test"))
    monkeypatch.setattr("rnaseq.service.check_docker", lambda: RuntimeCheck("Docker", "FOUND", "test"))
    monkeypatch.setattr("rnaseq.service.check_container_runtime", lambda *_args: RuntimeCheck("Control-plane container", "FOUND", "test"))
    monkeypatch.setattr("rnaseq.service.inspect_container_image", lambda image: {"reference": image, "image_id": None, "repo_digests": []})
    with pytest.raises(ExecutionPreflightError, match="observed immutable"):
        prepare_service_run(report, profile="local")


def test_production_acceptance_rejects_custom_reference_and_latest_runtime(project_factory):
    config = deepcopy(base_config())
    config["input"] = {"type": "fastq", "path": "input/fastq", "layout": "paired_end"}
    config["upstream"] = {
        "engine": "nfcore_rnaseq", "pipeline_version": "3.26.0", "strandedness": "auto",
        "quantification": {"method": "salmon"},
    }
    config["reference"] = {"source": "custom", "fasta": "reference.fa", "gtf": "genes.gtf", "acceptance": "production"}
    report = validate_project(project_factory(config=config))
    assert "invalid_project_config" in {issue.code for issue in report.errors}


def test_not_built_local_reference_is_valid_but_not_execution_ready(tmp_path):
    root, reference = _local_fastq_project(tmp_path)
    report = validate_project(root)
    assert report.is_valid
    assert not report.execution_ready
    assert report.local_reference is not None
    assert report.local_reference.provider == "Ensembl"
    first = [path.read_bytes() for path in generate_plan(report)]
    second = [path.read_bytes() for path in generate_plan(validate_project(root))]
    assert first == second
    preview = yaml.safe_load((root / "planning" / "upstream_run.preview.yaml").read_text(encoding="utf-8"))
    assert preview["reference"]["source"] == "local"
    assert preview["reference"]["release"] == 116
    assert preview["nfcore_reference_arguments"] == []
    assert preview["transcriptome_strategy"] == TRANSCRIPTOME_STRATEGY
    assert "no validated pre-built Salmon index" in preview["salmon_index_strategy"]
    assert preview["blocking_requirements"] == [
        "Local reference manifest has no validated Salmon index; register a compatible prebuilt index or use the optional host-native builder."
    ]
    assert preview["reference"]["transcriptome"]["external_transcript_fasta_used"] is False
    assert preview["reference"]["assets"]["transcript_fasta"]["sha256"]
    with pytest.raises(ExecutionPreflightError, match="no validated Salmon index"):
        build_nextflow_command(
            report,
            samplesheet=root / "planning" / "samplesheet.csv",
            output_dir=root / "runs" / "out",
            profile="local",
        )


def test_empty_fastq_input_does_not_misreport_a_valid_local_reference(tmp_path):
    root, reference = _local_fastq_project(tmp_path)
    _mark_salmon_built(reference)
    for path in (root / "input" / "fastq").iterdir():
        path.unlink()

    report = validate_project(root)

    assert not report.is_valid
    assert report.local_reference is not None
    assert report.execution_blockers == (
        "FASTQ input validation must pass before execution readiness can be confirmed.",
    )


def test_built_local_reference_passes_index_without_external_transcript_fasta(tmp_path):
    root, reference = _local_fastq_project(tmp_path)
    index = _mark_salmon_built(reference)
    report = validate_project(root)
    generate_plan(report)
    preview = yaml.safe_load((root / "planning" / "upstream_run.preview.yaml").read_text(encoding="utf-8"))
    assert preview["nfcore_reference_arguments"] == [
        "--fasta", str(reference / "fasta" / "Homo_sapiens.GRCh38.dna.primary_assembly.fa"),
        "--gtf", str(reference / "annotation" / "Homo_sapiens.GRCh38.116.gtf"),
        "--salmon_index", str(index),
    ]
    command = build_nextflow_command(
        report,
        samplesheet=root / "planning" / "samplesheet.csv",
        output_dir=root / "runs" / "out",
        profile="local",
    )
    assert "--genome" not in command
    assert command[command.index("--salmon_index") + 1] == str(index)
    assert command[command.index("--fasta") + 1] == str(reference / "fasta" / "Homo_sapiens.GRCh38.dna.primary_assembly.fa")
    assert command[command.index("--gtf") + 1] == str(reference / "annotation" / "Homo_sapiens.GRCh38.116.gtf")
    assert "--transcript_fasta" not in command
    assert command.count("--salmon_index") == 1
    assert "--pseudo_aligner" in command and command[command.index("--pseudo_aligner") + 1] == "salmon"
    assert "--skip_trimming" not in command


def test_transcriptome_only_adoption_is_explicit_atomic_and_ready_for_nfcore(tmp_path):
    root, reference = _local_fastq_project(tmp_path)
    index = _adopt_gtf_derived_candidate(reference)
    manifest = yaml.safe_load((reference / "reference_manifest.yaml").read_text(encoding="utf-8"))
    salmon = manifest["salmon"]
    assert salmon["strategy"] == SALMON_STRATEGY_TRANSCRIPTOME_ONLY
    assert salmon["index"] == "salmon/gtf_derived/transcriptome_only_index"
    assert salmon["transcriptome"]["path"] == "salmon/gtf_derived/transcripts.fa"
    assert salmon["transcriptome"]["sha256"] == _sha(reference / "salmon" / "gtf_derived" / "transcripts.fa")
    assert salmon["version"] == SALMON_VERSION
    assert salmon["kmer_size"] == 31
    assert salmon["num_decoys"] == 0
    assert salmon["index_metadata"]["seq_hash"] == "a" * 64
    assert salmon["index_metadata"]["name_hash"] == "b" * 64
    assert salmon["provenance"]["mode"] == ADOPTED_EXISTING_INDEX

    report = validate_project(root)
    assert report.is_valid and report.execution_ready
    assert report.local_reference is not None
    assert report.local_reference.salmon_strategy_type == SALMON_STRATEGY_TRANSCRIPTOME_ONLY
    generate_plan(report)
    preview = yaml.safe_load((root / "planning" / "upstream_run.preview.yaml").read_text(encoding="utf-8"))
    assert preview["nfcore_reference_arguments"] == [
        "--fasta", str(reference / "fasta" / "Homo_sapiens.GRCh38.dna.primary_assembly.fa"),
        "--gtf", str(reference / "annotation" / "Homo_sapiens.GRCh38.116.gtf"),
        "--salmon_index", str(index),
    ]
    command = build_nextflow_command(
        report,
        samplesheet=root / "planning" / "samplesheet.csv",
        output_dir=root / "runs" / "out",
        profile="local",
    )
    assert command[command.index("--salmon_index") + 1] == str(index)
    assert command[command.index("--fasta") + 1] == str(reference / "fasta" / "Homo_sapiens.GRCh38.dna.primary_assembly.fa")
    assert command[command.index("--gtf") + 1] == str(reference / "annotation" / "Homo_sapiens.GRCh38.116.gtf")
    assert "--transcript_fasta" not in command
    assert "--skip_alignment" not in command
    assert preview["nfcore_runtime_params"]["skip_alignment"] is True
    assert preview["nfcore_dynamic_salmon_index_build"] is False


@pytest.mark.parametrize(
    ("mutation", "message"),
    [
        (lambda data: data.__setitem__("status", "BLOCKED"), "artifact status"),
        (lambda data: data["identifier_compatibility"].__setitem__("fasta_only", 1), "fasta_only"),
        (lambda data: data["tx2gene"].__setitem__("zero_gene_mappings", 1), "zero_gene_mappings"),
        (lambda data: data["tx2gene"].__setitem__("multi_gene_mappings", 1), "multi_gene_mappings"),
        (lambda data: data["fasta_validation"].__setitem__("duplicate_transcript_ids", 1), "duplicate_transcript_ids"),
        (lambda data: data["inputs"]["genome"].__setitem__("sha256", "0" * 64), "inputs.genome.sha256"),
        (lambda data: data["inputs"]["gtf"].__setitem__("sha256", "0" * 64), "inputs.gtf.sha256"),
    ],
)
def test_transcriptome_only_adoption_rejects_invalid_artifacts_atomically(tmp_path, mutation, message):
    _root, reference = _local_fastq_project(tmp_path)
    index, transcriptome, artifact = _write_gtf_derived_candidate(reference)
    payload = json.loads(artifact.read_text(encoding="utf-8"))
    mutation(payload)
    artifact.write_text(json.dumps(payload), encoding="utf-8")
    before = (reference / "reference_manifest.yaml").read_bytes()
    with pytest.raises(ReferenceAdoptionError, match=message):
        adopt_local_salmon_index(
            reference,
            index=str(index.relative_to(reference)),
            transcriptome=str(transcriptome.relative_to(reference)),
            strategy=SALMON_STRATEGY_TRANSCRIPTOME_ONLY,
            validation_artifact=str(artifact.relative_to(reference)),
        )
    assert (reference / "reference_manifest.yaml").read_bytes() == before


@pytest.mark.parametrize(
    ("mutation", "message"),
    [
        (lambda index: (index / "seq.bin").unlink(), "index is incomplete"),
        (lambda index: (index / "info.json").write_text("not-json", encoding="utf-8"), "info.json is not valid JSON"),
        (lambda index: (index / "versionInfo.json").write_text("not-json", encoding="utf-8"), "versionInfo.json is not valid JSON"),
        (lambda index: (index / "info.json").write_text(json.dumps({"index_version": 4}), encoding="utf-8"), "missing k"),
    ],
)
def test_transcriptome_only_adoption_rejects_incomplete_or_malformed_index(tmp_path, mutation, message):
    _root, reference = _local_fastq_project(tmp_path)
    index, transcriptome, artifact = _write_gtf_derived_candidate(reference)
    mutation(index)
    before = (reference / "reference_manifest.yaml").read_bytes()
    with pytest.raises((ReferenceAdoptionError, LocalReferenceError), match=message):
        adopt_local_salmon_index(
            reference,
            index=str(index.relative_to(reference)),
            transcriptome=str(transcriptome.relative_to(reference)),
            strategy=SALMON_STRATEGY_TRANSCRIPTOME_ONLY,
            validation_artifact=str(artifact.relative_to(reference)),
        )
    assert (reference / "reference_manifest.yaml").read_bytes() == before


def test_transcriptome_only_normal_validation_detects_changed_frozen_assets(tmp_path):
    root, reference = _local_fastq_project(tmp_path)
    _adopt_gtf_derived_candidate(reference)
    transcriptome = reference / "salmon" / "gtf_derived" / "transcripts.fa"
    transcriptome.write_text(">TX1\nACGTA\n", encoding="utf-8")
    report = validate_project(root)
    assert not report.is_valid
    assert "checksum mismatch for salmon_transcriptome" in report.errors[0].message


def test_reference_strategy_validation_and_cli_adoption(tmp_path):
    _root, reference = _local_fastq_project(tmp_path)
    index, transcriptome, artifact = _write_gtf_derived_candidate(reference)
    with pytest.raises(ReferenceAdoptionError, match="Unsupported Salmon strategy"):
        adopt_local_salmon_index(
            reference,
            index=str(index.relative_to(reference)),
            transcriptome=str(transcriptome.relative_to(reference)),
            strategy="unknown",
            validation_artifact=str(artifact.relative_to(reference)),
        )
    result = runner.invoke(app, [
        "reference", "adopt-salmon-index", str(reference),
        "--index", str(index.relative_to(reference)),
        "--transcriptome", str(transcriptome.relative_to(reference)),
        "--strategy", SALMON_STRATEGY_TRANSCRIPTOME_ONLY,
        "--validation-artifact", str(artifact.relative_to(reference)),
    ])
    assert result.exit_code == 0
    assert "Salmon strategy: transcriptome_only" in result.output


def test_reference_prepare_atomically_builds_and_registers_expected_nfcore_strategy(tmp_path):
    _root, reference = _local_fastq_project(tmp_path)
    commands: list[list[str]] = []

    def fake_runner(command, **_kwargs):
        commands.append(command)
        if command[-1] == "--version":
            output = "salmon 1.10.3" if Path(command[0]).name == "salmon" else "RSEM v1.3.3"
            return subprocess.CompletedProcess(command, 0, output, "")
        if Path(command[0]).name == "salmon":
            _write_salmon_index(Path(command[command.index("-i") + 1]), num_decoys=1)
        return subprocess.CompletedProcess(command, 0, "", "")

    prepared = prepare_local_reference(
        reference,
        threads=5,
        runner=fake_runner,
        tool_resolver=lambda name: f"/tools/{name}",
        now=lambda: datetime(2026, 8, 31, tzinfo=UTC),
    )
    assert prepared.salmon_status == SALMON_BUILT
    assert prepared.salmon_index == reference / "salmon" / "index"
    manifest = yaml.safe_load((reference / "reference_manifest.yaml").read_text(encoding="utf-8"))
    assert manifest["salmon"]["status"] == SALMON_BUILT
    assert manifest["salmon"]["index"] == "salmon/index"
    assert manifest["salmon"]["provenance"]["kmer_size"] == 31
    assert manifest["salmon"]["provenance"]["decoy_strategy"] == DECOY_STRATEGY
    assert manifest["salmon"]["strategy"] == SALMON_STRATEGY_DECOY_AWARE
    assert manifest["salmon"]["provenance"]["builder"]["mode"] == "host_native"
    assert manifest["salmon"]["provenance"]["builder"]["tools"]["salmon"]["executable"] == "/tools/salmon"
    assert manifest["salmon"]["provenance"]["builder"]["preflight"]["threads"] == 5
    assert manifest["salmon"]["provenance"]["source_assets"]["genome_fasta"]["path"] == "fasta/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
    assert manifest["salmon"]["provenance"]["commands"]["salmon_index"][-1] == "31"
    assert "--threads" in manifest["salmon"]["provenance"]["commands"]["salmon_index"]
    assert all("docker" not in item for command in commands for item in command)


def test_reference_prepare_failure_leaves_manifest_not_built_and_no_index(tmp_path):
    _root, reference = _local_fastq_project(tmp_path)

    def failing_runner(command, **_kwargs):
        if command[-1] == "--version":
            output = "salmon 1.10.3" if Path(command[0]).name == "salmon" else "RSEM v1.3.3"
            return subprocess.CompletedProcess(command, 0, output, "")
        return subprocess.CompletedProcess(command, 42, "", "fixture failure")

    with pytest.raises(ReferencePreparationError, match="exit 42"):
        prepare_local_reference(reference, runner=failing_runner, tool_resolver=lambda name: f"/tools/{name}")
    manifest = yaml.safe_load((reference / "reference_manifest.yaml").read_text(encoding="utf-8"))
    assert manifest["salmon"] == {"index": None, "status": "not_built"}
    assert not (reference / "salmon" / "index").exists()
    assert list(reference.glob(".rnaseq-reference-prepare-*"))


def test_manifest_publish_failure_returns_new_index_to_retained_staging(monkeypatch, tmp_path):
    _root, reference = _local_fastq_project(tmp_path)
    before = (reference / "reference_manifest.yaml").read_bytes()

    def fake_runner(command, **_kwargs):
        executable = Path(command[0]).name
        if command[-1] == "--version":
            return subprocess.CompletedProcess(command, 0, "salmon 1.10.3" if executable == "salmon" else "RSEM v1.3.3", "")
        if executable == "salmon":
            _write_salmon_index(Path(command[command.index("-i") + 1]), num_decoys=1)
        return subprocess.CompletedProcess(command, 0, "", "")

    monkeypatch.setattr("rnaseq.references._write_manifest_atomically", lambda *_args, **_kwargs: (_ for _ in ()).throw(OSError("fixture manifest failure")))
    with pytest.raises(OSError, match="fixture manifest failure"):
        prepare_local_reference(reference, runner=fake_runner, tool_resolver=lambda name: f"/tools/{name}")
    assert (reference / "reference_manifest.yaml").read_bytes() == before
    assert not (reference / "salmon" / "index").exists()
    assert any(path.joinpath("salmon", "index").is_dir() for path in reference.glob(".rnaseq-reference-prepare-*"))


def test_host_native_reference_preflight_fails_before_staging_for_missing_or_unsupported_tools(tmp_path):
    _root, reference = _local_fastq_project(tmp_path)
    with pytest.raises(ReferencePreparationError, match="Missing required host executable 'salmon'"):
        prepare_local_reference(reference, tool_resolver=lambda _name: None)
    assert not list(reference.glob(".rnaseq-reference-prepare-*"))

    with pytest.raises(ReferencePreparationError, match="threads must be at least 1"):
        prepare_local_reference(reference, threads=0, tool_resolver=lambda name: f"/tools/{name}")
    assert not list(reference.glob(".rnaseq-reference-prepare-*"))

    def old_salmon(command, **_kwargs):
        output = "salmon 1.9.0" if Path(command[0]).name == "salmon" else "RSEM v1.3.3"
        return subprocess.CompletedProcess(command, 0, output, "")

    with pytest.raises(ReferencePreparationError, match="Unsupported salmon version"):
        prepare_local_reference(reference, runner=old_salmon, tool_resolver=lambda name: f"/tools/{name}")
    assert not list(reference.glob(".rnaseq-reference-prepare-*"))


def test_host_native_reference_command_is_tokenized_for_paths_with_spaces(tmp_path):
    reference = tmp_path / "reference with spaces"
    _write_reference(reference)
    commands: list[list[str]] = []

    def fake_runner(command, **_kwargs):
        commands.append(command)
        executable = Path(command[0]).name
        if command[-1] == "--version":
            return subprocess.CompletedProcess(command, 0, "salmon 1.10.3" if executable == "salmon" else "RSEM v1.3.3", "")
        if executable == "salmon":
            _write_salmon_index(Path(command[command.index("-i") + 1]), num_decoys=1)
        return subprocess.CompletedProcess(command, 0, "", "")

    prepare_local_reference(reference, runner=fake_runner, tool_resolver=lambda name: f"/tools/{name}")
    assert any(str(reference) in argument for command in commands for argument in command)
    assert all(isinstance(argument, str) for command in commands for argument in command)
    assert all("docker" not in argument.lower() for command in commands for argument in command)


def test_existing_valid_index_is_not_rebuilt_or_relabelled(tmp_path):
    _root, reference = _local_fastq_project(tmp_path)
    _mark_salmon_built(reference)
    before = (reference / "reference_manifest.yaml").read_bytes()
    with pytest.raises(ReferencePreparationError, match="already built"):
        prepare_local_reference(reference, tool_resolver=lambda name: f"/tools/{name}")
    assert (reference / "reference_manifest.yaml").read_bytes() == before


def test_source_checksum_mismatch_blocks_preparation_before_host_tool_lookup(tmp_path):
    _root, reference = _local_fastq_project(tmp_path)
    (reference / "fasta" / "Homo_sapiens.GRCh38.dna.primary_assembly.fa").write_text(">chr1\nACGTA\n", encoding="utf-8")
    with pytest.raises(LocalReferenceError, match="checksum mismatch for genome_fasta"):
        prepare_local_reference(reference, tool_resolver=lambda _name: (_ for _ in ()).throw(AssertionError("must not resolve tools")))


def test_host_native_hisat2_preparation_records_annotation_aware_provenance(tmp_path):
    _root, reference = _local_fastq_project(tmp_path)
    commands: list[list[str]] = []

    def fake_runner(command, **_kwargs):
        commands.append(command)
        executable = Path(command[0]).name
        if command[-1] == "--version":
            if executable == "hisat2_extract_splice_sites.py":
                return subprocess.CompletedProcess(command, 2, "", "unrecognized arguments: --version")
            return subprocess.CompletedProcess(command, 0, "hisat2-build version 2.2.3", "")
        if executable == "hisat2_extract_splice_sites.py" and command[-1] == "-h":
            return subprocess.CompletedProcess(
                command, 0, "usage: hisat2_extract_splice_sites.py [-h] [-v] [gtf_file]\nExtract splice junctions from a GTF file\n", ""
            )
        if executable == "hisat2_extract_splice_sites.py":
            return subprocess.CompletedProcess(command, 0, "chr1\t1\t4\t+\n", "")
        _write_hisat2_index(Path(command[-1]).parent)
        return subprocess.CompletedProcess(command, 0, "", "")

    prepared = prepare_local_hisat2_reference(
        reference, threads=7, runner=fake_runner, tool_resolver=lambda name: f"/tools/{name}",
        now=lambda: datetime(2026, 8, 31, tzinfo=UTC),
    )
    assert prepared.hisat2_index == reference / "hisat2" / "index"
    manifest = yaml.safe_load((reference / "reference_manifest.yaml").read_text(encoding="utf-8"))
    provenance = manifest["hisat2"]["provenance"]
    assert provenance["builder"]["mode"] == "host_native"
    assert provenance["threads"] == 7
    assert provenance["index_strategy"] == "genome_only_runtime_splicesites"
    assert manifest["hisat2"]["version"] == "2.2.3"
    assert provenance["index_builder_version"] == "2.2.3"
    assert provenance["runtime_aligner_version"] == "2.2.3"
    helper = provenance["builder"]["tools"]["hisat2_extract_splice_sites.py"]
    assert helper["validation"] == "help_contract"
    assert helper["associated_hisat2_build_version"] == "2.2.3"
    assert "version" not in helper
    assert ["/tools/hisat2_extract_splice_sites.py", "--version"] not in commands
    assert ["/tools/hisat2_extract_splice_sites.py", "-h"] in commands
    assert "--ss" not in provenance["commands"]["hisat2_build"]
    assert all("docker" not in item for command in commands for item in command)
    assert dict(prepared.hisat2_arguments())["--hisat2_splice_sites"] == reference / "hisat2" / "splice_sites.txt"


@pytest.mark.parametrize("observed", ("2.2.1", "2.2.2", "2.3.0"))
def test_host_native_hisat2_preparation_rejects_every_noncanonical_version(tmp_path, observed):
    _root, reference = _local_fastq_project(tmp_path)

    def fake_runner(command, **_kwargs):
        executable = Path(command[0]).name
        output = f"hisat2-build version {observed}" if executable == "hisat2-build" else f"hisat2 {observed}"
        return subprocess.CompletedProcess(command, 0, output, "")

    with pytest.raises(ReferencePreparationError, match=rf"expected {HISAT2_VERSION}"):
        prepare_local_hisat2_reference(
            reference, runner=fake_runner, tool_resolver=lambda name: f"/tools/{name}"
        )
    assert not list(reference.glob(".rnaseq-hisat2-prepare-*"))
    assert not (reference / "hisat2" / "index").exists()


def test_hisat2_preparation_rejects_missing_or_nonfunctional_splice_site_helper(tmp_path):
    _root, reference = _local_fastq_project(tmp_path)

    def valid_build(command, **_kwargs):
        return subprocess.CompletedProcess(command, 0, "hisat2-build version 2.2.3", "")

    with pytest.raises(ReferencePreparationError, match="Missing required host executable 'hisat2_extract_splice_sites.py'"):
        prepare_local_hisat2_reference(
            reference,
            runner=valid_build,
            tool_resolver=lambda name: "/tools/hisat2-build" if name == "hisat2-build" else None,
        )
    assert not list(reference.glob(".rnaseq-hisat2-prepare-*"))

    def broken_helper(command, **_kwargs):
        executable = Path(command[0]).name
        if executable == "hisat2-build":
            return subprocess.CompletedProcess(command, 0, "hisat2-build version 2.2.3", "")
        return subprocess.CompletedProcess(command, 2, "", "helper failure")

    with pytest.raises(ReferencePreparationError, match="expected HISAT2 splice-site helper help contract"):
        prepare_local_hisat2_reference(
            reference,
            runner=broken_helper,
            tool_resolver=lambda name: f"/tools/{name}",
        )
    assert not list(reference.glob(".rnaseq-hisat2-prepare-*"))


def test_existing_2_2_3_index_is_not_rebuilt_or_promoted_without_smoke_validation(tmp_path):
    _root, reference = _local_fastq_project(tmp_path)
    _mark_genome_only_hisat2(reference, runtime_compatibility="requires_smoke_validation")
    manifest_path = reference / "reference_manifest.yaml"
    before = manifest_path.read_bytes()
    index_before = {
        path.name: path.read_bytes() for path in (reference / "hisat2" / "index").iterdir()
    }

    loaded = load_local_reference_root(reference)
    assert loaded.hisat2_index is not None
    assert not loaded.hisat2_runtime_ready
    with pytest.raises(ReferencePreparationError, match="already built"):
        prepare_local_hisat2_reference(reference, tool_resolver=lambda name: f"/tools/{name}")
    assert manifest_path.read_bytes() == before
    assert {
        path.name: path.read_bytes() for path in (reference / "hisat2" / "index").iterdir()
    } == index_before
    assert not loaded.hisat2_runtime_ready


@pytest.mark.parametrize("observed", ("2.2.1", "2.2.2", "2.3.0"))
def test_manifest_rejects_prebuilt_hisat2_indexes_outside_the_exact_contract(tmp_path, observed):
    _root, reference = _local_fastq_project(tmp_path)
    _mark_genome_only_hisat2(reference)
    manifest_path = reference / "reference_manifest.yaml"
    manifest = yaml.safe_load(manifest_path.read_text(encoding="utf-8"))
    manifest["hisat2"]["provenance"]["index_builder_version"] = observed
    manifest_path.write_text(yaml.safe_dump(manifest, sort_keys=False), encoding="utf-8")

    with pytest.raises(LocalReferenceError, match=r"requires exactly 2\.2\.3"):
        load_local_reference_root(reference)


def test_genome_only_runtime_splices_require_registered_matching_nonempty_provenance(tmp_path):
    _root, reference = _local_fastq_project(tmp_path)
    _mark_genome_only_hisat2(reference)
    loaded = load_local_reference_root(reference)
    assert loaded.hisat2_strategy == "genome_only_runtime_splices"
    assert loaded.hisat2_index_prefix == reference / "hisat2" / "index" / "genome"
    assert loaded.hisat2_arguments()[-1] == ("--hisat2_splice_sites", reference / "hisat2" / "splice_sites.txt")
    assert loaded.provenance()["hisat2"]["index_builder_version"] == "2.2.3"
    assert loaded.provenance()["hisat2"]["runtime_aligner_version"] == "2.2.3"

    splice = reference / "hisat2" / "splice_sites.txt"
    splice.unlink()
    with pytest.raises(LocalReferenceError, match="hisat2.splice_sites file not found"):
        load_local_reference_root(reference)

    _mark_genome_only_hisat2(reference)
    splice.write_text("", encoding="utf-8")
    # A reference with no introns has a valid, empty runtime splice-site file.
    manifest_path = reference / "reference_manifest.yaml"
    manifest = yaml.safe_load(manifest_path.read_text(encoding="utf-8"))
    manifest["hisat2"]["splice_sites"]["sha256"] = _sha(splice)
    manifest_path.write_text(yaml.safe_dump(manifest, sort_keys=False), encoding="utf-8")
    assert load_local_reference_root(reference).hisat2_splice_sites.path == splice

    _mark_genome_only_hisat2(reference)
    manifest_path = reference / "reference_manifest.yaml"
    manifest = yaml.safe_load(manifest_path.read_text(encoding="utf-8"))
    manifest["hisat2"]["provenance"]["splice_sites_derived_from_gtf_sha256"] = "0" * 64
    manifest_path.write_text(yaml.safe_dump(manifest, sort_keys=False), encoding="utf-8")
    with pytest.raises(LocalReferenceError, match="splice_sites_derived_from_gtf_sha256"):
        load_local_reference_root(reference)


def test_manifest_accepts_prebuilt_salmon_without_nf_rna_builder_provenance(tmp_path):
    _root, reference = _local_fastq_project(tmp_path)
    index = reference / "vendor indexes" / "salmon"
    _write_salmon_index(index)
    manifest_path = reference / "reference_manifest.yaml"
    manifest = yaml.safe_load(manifest_path.read_text(encoding="utf-8"))
    manifest["schema_version"] = "1.2"
    manifest["purpose"] = "production"
    manifest["salmon"] = {
        "status": "built",
        "index": "vendor indexes/salmon",
        "version": "1.10.3",
        "strategy": "transcriptome_only",
        "source_transcriptome_sha256": manifest["files"]["transcript_fasta"]["sha256"],
    }
    manifest_path.write_text(yaml.safe_dump(manifest, sort_keys=False), encoding="utf-8")

    loaded = load_local_reference_root(reference)
    assert loaded.salmon_index == index
    assert loaded.salmon_strategy_type == "transcriptome_only"
    assert loaded.external_transcript_fasta_used
    assert loaded.nfcore_arguments()[-1] == ("--salmon_index", index)


def test_prebuilt_salmon_rejects_source_checksum_or_strategy_mismatch(tmp_path):
    _root, reference = _local_fastq_project(tmp_path)
    index = reference / "vendor" / "salmon"
    _write_salmon_index(index, num_decoys=1)
    manifest_path = reference / "reference_manifest.yaml"
    manifest = yaml.safe_load(manifest_path.read_text(encoding="utf-8"))
    manifest["schema_version"] = "1.2"
    manifest["purpose"] = "production"
    manifest["salmon"] = {
        "status": "built", "index": "vendor/salmon", "version": "1.10.3",
        "strategy": "decoy_aware", "source_transcriptome_sha256": "0" * 64,
        "source_genome_sha256": manifest["files"]["genome_fasta"]["sha256"],
    }
    manifest_path.write_text(yaml.safe_dump(manifest, sort_keys=False), encoding="utf-8")
    with pytest.raises(LocalReferenceError, match="source_transcriptome_sha256"):
        load_local_reference_root(reference)


def test_adopt_prebuilt_decoy_aware_salmon_without_rsem_or_validation_artifact(tmp_path):
    _root, reference = _local_fastq_project(tmp_path)
    index = reference / "external" / "salmon"
    _write_salmon_index(index, num_decoys=1)

    adopted = adopt_local_salmon_index(
        reference,
        index="external/salmon",
        strategy="decoy_aware",
        now=lambda: datetime(2026, 9, 8, tzinfo=UTC),
    )
    assert adopted.salmon_index == index
    manifest = yaml.safe_load((reference / "reference_manifest.yaml").read_text(encoding="utf-8"))
    assert manifest["salmon"]["strategy"] == "decoy_aware"
    assert manifest["salmon"]["source_transcriptome_sha256"] == manifest["files"]["transcript_fasta"]["sha256"]
    assert manifest["salmon"]["provenance"]["mode"] == ADOPTED_EXISTING_INDEX


def test_manifest_accepts_prebuilt_genome_only_hisat2_prefix_without_builder_provenance(tmp_path):
    _root, reference = _local_fastq_project(tmp_path)
    index = reference / "third party" / "hisat2"
    _write_hisat2_index(index)
    splice = reference / "third party" / "splice_sites.txt"
    splice.write_text("chr1\t1\t4\t+\n", encoding="utf-8")
    manifest_path = reference / "reference_manifest.yaml"
    manifest = yaml.safe_load(manifest_path.read_text(encoding="utf-8"))
    manifest["schema_version"] = "1.2"
    manifest["purpose"] = "production"
    # A HISAT2-only reference does not need a Salmon transcript FASTA.
    manifest["files"].pop("transcript_fasta")
    manifest["hisat2"] = {
        "status": "built", "index_prefix": "third party/hisat2/genome", "version": "2.2.3",
        "strategy": "genome_only_runtime_splicesites",
        "genome_fasta_sha256": manifest["files"]["genome_fasta"]["sha256"],
        "source_gtf_sha256": manifest["files"]["annotation_gtf"]["sha256"],
        "splice_sites_gtf_sha256": manifest["files"]["annotation_gtf"]["sha256"],
        "splice_sites": {"path": "third party/splice_sites.txt", "sha256": _sha(splice)},
        "runtime_compatibility": "validated",
    }
    manifest_path.write_text(yaml.safe_dump(manifest, sort_keys=False), encoding="utf-8")

    loaded = load_local_reference_root(reference)
    assert loaded.hisat2_index == index
    assert loaded.hisat2_index_prefix == index / "genome"
    assert loaded.hisat2_runtime_ready
    assert loaded.hisat2_arguments()[-1] == ("--hisat2_splice_sites", splice)
    config_path = _root / "project.yaml"
    config = yaml.safe_load(config_path.read_text(encoding="utf-8"))
    config["upstream"]["strandedness"] = "forward"
    config["upstream"]["quantification"] = {"method": "hisat2_featurecounts"}
    config_path.write_text(yaml.safe_dump(config, sort_keys=False), encoding="utf-8")
    report = validate_project(_root)
    assert report.execution_ready
    assert generate_plan(report)


def test_genome_only_runtime_splices_require_compatibility_smoke_before_ready(tmp_path):
    root, reference = _local_fastq_project(tmp_path)
    _mark_genome_only_hisat2(reference, runtime_compatibility="requires_smoke_validation")
    config_path = root / "project.yaml"
    config = yaml.safe_load(config_path.read_text(encoding="utf-8"))
    config["upstream"]["strandedness"] = "forward"
    config["upstream"]["quantification"] = {"method": "hisat2_featurecounts"}
    config_path.write_text(yaml.safe_dump(config, sort_keys=False), encoding="utf-8")
    report = validate_project(root)
    assert report.is_valid and not report.execution_ready
    assert "requires the documented smoke validation" in report.execution_blockers[0]


def test_genome_only_runtime_splices_are_frozen_and_passed_once(tmp_path):
    root, reference = _local_fastq_project(tmp_path)
    _mark_genome_only_hisat2(reference)
    config_path = root / "project.yaml"
    config = yaml.safe_load(config_path.read_text(encoding="utf-8"))
    config["upstream"]["strandedness"] = "forward"
    config["upstream"]["quantification"] = {"method": "hisat2_featurecounts"}
    config_path.write_text(yaml.safe_dump(config, sort_keys=False), encoding="utf-8")
    report = validate_project(root)
    assert report.execution_ready
    command = build_hisat2_featurecounts_command(
        report, samplesheet=root / "samples.csv", output_dir=root / "out", profile="local"
    )
    assert command.count("--hisat2_splice_sites") == 1
    assert command[command.index("--hisat2_splice_sites") + 1] == str(reference / "hisat2" / "splice_sites.txt")
    assert command[command.index("--hisat2_index_basename") + 1] == "genome"
    assert command.count("--hisat2_use_runtime_splices") == 1
    run_dir = root / "runs" / "CASE" / "20260908-120000+0800"
    run_dir.joinpath("frozen").mkdir(parents=True)
    run = CaseRun("CASE", "20260908-120000+0800", run_dir, "2026-09-08T12:00:00+08:00")
    freeze_case_inputs(report, run, profile="local", command=["rnaseq", "run"])
    execution = yaml.safe_load((run_dir / "frozen" / "execution_manifest.yaml").read_text(encoding="utf-8"))
    frozen_hisat2 = execution["reference"]["hisat2"]
    assert frozen_hisat2["strategy"] == "genome_only_runtime_splices"
    assert frozen_hisat2["index_builder_version"] == "2.2.3"
    assert frozen_hisat2["runtime_aligner_version"] == "2.2.3"
    assert frozen_hisat2["splice_sites"]["sha256"] == _sha(reference / "hisat2" / "splice_sites.txt")


def test_reference_prepare_cli_exposes_the_one_time_command(monkeypatch, tmp_path):
    _root, reference = _local_fastq_project(tmp_path)
    index = _mark_salmon_built(reference)

    class Prepared:
        root = reference
        salmon_index = index

    monkeypatch.setattr("rnaseq.cli.prepare_local_reference", lambda root, **_kwargs: Prepared())
    result = runner.invoke(app, ["reference", "prepare", str(reference)])
    assert result.exit_code == 0
    assert f"Salmon index: {index}" in result.output
    assert "Salmon status: built" in result.output


def test_incomplete_or_mismatched_built_index_blocks_execution(tmp_path):
    root, reference = _local_fastq_project(tmp_path)
    index = _mark_salmon_built(reference)
    (index / "seq.bin").unlink()
    report = validate_project(root)
    assert not report.is_valid
    assert "Salmon index is incomplete" in report.errors[0].message

    root, reference = _local_fastq_project(tmp_path / "mismatch")
    _mark_salmon_built(reference)
    manifest_path = reference / "reference_manifest.yaml"
    manifest = yaml.safe_load(manifest_path.read_text(encoding="utf-8"))
    manifest["salmon"]["provenance"]["genome_fasta_sha256"] = "0" * 64
    manifest_path.write_text(yaml.safe_dump(manifest, sort_keys=False), encoding="utf-8")
    report = validate_project(root)
    assert not report.is_valid
    assert "salmon.provenance.genome_fasta_sha256 is incompatible" in report.errors[0].message


def test_local_reference_missing_manifest_or_required_assets_fails(tmp_path):
    root, reference = _local_fastq_project(tmp_path)
    (reference / "reference_manifest.yaml").unlink()
    report = validate_project(root)
    assert not report.is_valid
    assert "manifest file not found" in report.errors[0].message

    root, reference = _local_fastq_project(tmp_path / "missing_genome")
    (reference / "fasta" / "Homo_sapiens.GRCh38.dna.primary_assembly.fa").unlink()
    report = validate_project(root)
    assert not report.is_valid
    assert "genome_fasta file not found" in report.errors[0].message

    root, reference = _local_fastq_project(tmp_path / "missing_gtf")
    (reference / "annotation" / "Homo_sapiens.GRCh38.116.gtf").unlink()
    report = validate_project(root)
    assert not report.is_valid
    assert "annotation_gtf file not found" in report.errors[0].message

    root, reference = _local_fastq_project(tmp_path / "missing_transcript")
    (reference / "salmon" / "Homo_sapiens.GRCh38.cdna.all.fa").unlink()
    report = validate_project(root)
    assert not report.is_valid
    assert "transcript_fasta file not found" in report.errors[0].message


def test_local_reference_checksum_and_species_mismatch_fail(tmp_path):
    root, reference = _local_fastq_project(tmp_path)
    genome = reference / "fasta" / "Homo_sapiens.GRCh38.dna.primary_assembly.fa"
    genome.write_text(">chr1\nACGTA\n", encoding="utf-8")
    report = validate_project(root)
    assert not report.is_valid
    assert "checksum mismatch for genome_fasta" in report.errors[0].message

    root, reference = _local_fastq_project(tmp_path / "mismatch")
    manifest = yaml.safe_load((reference / "reference_manifest.yaml").read_text(encoding="utf-8"))
    manifest["reference"]["species"] = "Mus musculus"
    (reference / "reference_manifest.yaml").write_text(yaml.safe_dump(manifest, sort_keys=False), encoding="utf-8")
    report = validate_project(root)
    assert not report.is_valid
    assert "does not match project organism" in report.errors[0].message


def test_malformed_local_reference_manifest_fails(tmp_path):
    root, reference = _local_fastq_project(tmp_path)
    (reference / "reference_manifest.yaml").write_text("files: []\n", encoding="utf-8")
    report = validate_project(root)
    assert not report.is_valid
    assert "Unsupported local reference manifest schema version" in report.errors[0].message


def test_local_reference_freezes_manifest_identity_not_assets(tmp_path):
    root, reference = _local_fastq_project(tmp_path)
    _mark_salmon_built(reference)
    report = validate_project(root)
    generate_plan(report)
    run = create_case_run(report, "CASE-REFERENCE")
    frozen = freeze_case_inputs(report, run, profile="local", command=["rnaseq", "run"])
    assert (run.run_dir / "frozen" / "reference" / "reference_manifest.yaml").is_file()
    payload = yaml.safe_load(frozen.manifest.read_text(encoding="utf-8"))
    assert payload["reference"]["source"] == "local"
    assert payload["reference"]["assets"]["genome_fasta"]["sha256"]
    assert payload["reference"]["transcriptome"]["external_transcript_fasta_used"] is False
    assert frozen.reference_paths["fasta"] == report.local_reference.genome_fasta.path
    assert frozen.reference_paths["gtf"] == report.local_reference.annotation_gtf.path
    assert frozen.reference_paths["salmon_index"] == report.local_reference.salmon_index
    assert "transcript_fasta" not in frozen.reference_paths


def test_pretrimmed_fastq_state_is_frozen_in_execution_contract(tmp_path):
    root, reference = _local_fastq_project(tmp_path)
    _mark_salmon_built(reference)
    config = yaml.safe_load((root / "project.yaml").read_text(encoding="utf-8"))
    config["input"]["preprocessing"] = "pretrimmed"
    (root / "project.yaml").write_text(yaml.safe_dump(config, sort_keys=False), encoding="utf-8")
    report = validate_project(root)
    generate_plan(report)
    run = create_case_run(report, "CASE-PRETRIMMED")
    frozen = freeze_case_inputs(report, run, profile="local", command=["rnaseq", "run"])
    execution_manifest = yaml.safe_load((run.run_dir / "frozen" / "execution_manifest.yaml").read_text(encoding="utf-8"))
    contract = json.loads((run.run_dir / "frozen" / "downstream_contract.json").read_text(encoding="utf-8"))
    assert execution_manifest["fastq_preprocessing"] == "pretrimmed"
    assert execution_manifest["skip_trimming"] is True
    assert execution_manifest["nfcore_preprocessing_arguments"] == ["--skip_trimming"]
    assert execution_manifest["nfcore_runtime_params"] == {
        "skip_alignment": True,
        "skip_trimming": True,
    }
    assert contract["source"]["fastq_preprocessing"] == "pretrimmed"
    assert contract["source"]["skip_trimming"] is True
    runtime_params_path = run.run_dir / "frozen" / "nfcore.params.json"
    runtime_params = json.loads(runtime_params_path.read_text(encoding="utf-8"))
    assert runtime_params == {"skip_alignment": True, "skip_trimming": True}
    assert runtime_params["skip_alignment"] is True
    assert runtime_params["skip_trimming"] is True
    assert '"skip_trimming": true' in runtime_params_path.read_text(encoding="utf-8")
    assert '"skip_trimming": "true"' not in runtime_params_path.read_text(encoding="utf-8")
    runtime_command = build_nextflow_command(
        report,
        samplesheet=frozen.samplesheet,
        output_dir=run.run_dir / "upstream" / "nfcore_rnaseq",
        profile="local",
        params_file=frozen.upstream_params,
        config_file=frozen.upstream_config,
        reference_paths=frozen.reference_paths,
    )
    assert "--fasta" in runtime_command and "--gtf" in runtime_command
    assert "--transcript_fasta" not in runtime_command


def test_local_reference_checksum_is_rechecked_before_execution(tmp_path):
    root, reference = _local_fastq_project(tmp_path)
    _mark_salmon_built(reference)
    report = validate_project(root)
    assert report.execution_ready
    (reference / "salmon" / "Homo_sapiens.GRCh38.cdna.all.fa").write_text(
        ">TX1 gene:GENE1\nACGTA\n", encoding="utf-8"
    )
    try:
        _validate_custom_reference_files(report)
    except ExecutionPreflightError as exc:
        assert "checksum mismatch for transcript_fasta" in str(exc)
    else:
        raise AssertionError("checksum change must block execution preflight")


def test_igenomes_and_raw_count_contracts_remain_compatible(project_factory, tmp_path):
    raw_report = validate_project(project_factory())
    assert raw_report.is_valid and raw_report.execution_ready
    root, _reference = _local_fastq_project(tmp_path)
    config = yaml.safe_load((root / "project.yaml").read_text(encoding="utf-8"))
    config["reference"] = {"source": "igenomes", "genome": "GRCh38"}
    (root / "project.yaml").write_text(yaml.safe_dump(config, sort_keys=False), encoding="utf-8")
    report = validate_project(root)
    assert report.is_valid and report.execution_ready
