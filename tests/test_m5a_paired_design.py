from __future__ import annotations

import json
from copy import deepcopy
from datetime import datetime
from pathlib import Path
from types import SimpleNamespace

import pytest
import yaml
from typer.testing import CliRunner

from conftest import base_config
from rnaseq.cli import app, render_validation_report
from rnaseq.execution import RuntimeCheck
from rnaseq.l2 import prepare_l2
from rnaseq.planner import generate_plan
from rnaseq import service
from rnaseq.service import create_case_run, freeze_case_inputs
from rnaseq.validators import validate_project

runner = CliRunner()


PAIRED_COUNTS = """gene_id,P01_Primary,P01_Metastasis,P02_Primary,P02_Metastasis,P03_Primary,P03_Metastasis
GeneA,10,20,11,24,9,19
GeneB,30,15,28,14,32,16
"""
PAIRED_METADATA = """sample_id,patient,condition
P01_Primary,P01,Primary
P01_Metastasis,P01,Metastasis
P02_Primary,P02,Primary
P02_Metastasis,P02,Metastasis
P03_Primary,P03,Primary
P03_Metastasis,P03,Metastasis
"""
PAIRED_CONTRAST = """contrast_id,factor,numerator,denominator
Metastasis_vs_Primary,condition,Metastasis,Primary
"""


def paired_config() -> dict:
    config = deepcopy(base_config())
    config["schema_version"] = "1.2"
    config["analysis"] = {"enrichment": []}
    config["design"] = {
        "type": "paired_two_group",
        "formula": "~ patient + condition",
        "pair_id": "patient",
    }
    return config


def paired_project(project_factory, *, metadata: str = PAIRED_METADATA, config: dict | None = None):
    samples = [line.split(",", 1)[0] for line in metadata.strip().splitlines()[1:]]
    counts = "gene_id," + ",".join(samples) + "\nGeneA," + ",".join("1" for _ in samples) + "\n"
    return project_factory(
        config=config or paired_config(), counts=counts, metadata=metadata, contrasts=PAIRED_CONTRAST
    )


def issue_codes(report) -> set[str]:
    return {issue.code for issue in report.errors}


def test_valid_paired_raw_count_contract_and_pair_statistics(project_factory):
    report = validate_project(
        project_factory(
            config=paired_config(), counts=PAIRED_COUNTS,
            metadata=PAIRED_METADATA, contrasts=PAIRED_CONTRAST,
        )
    )
    assert report.is_valid
    assert report.config.design.type.value == "paired_two_group"
    assert report.config.design.pair_id == "patient"
    assert report.pairing is not None
    summary = report.pairing.contrasts[0]
    assert (report.pairing.total_samples, report.pairing.unique_pair_ids) == (6, 3)
    assert (summary.complete_pairs, summary.samples_in_complete_pairs) == (3, 6)
    assert summary.incomplete_pairs == ()
    assert summary.membership[0].pair_id == "P01"
    rendered = render_validation_report(report)
    assert "Type: paired_two_group" in rendered
    assert "Pairing variable: patient" in rendered
    assert "Complete pairs: 3" in rendered


@pytest.mark.parametrize(
    ("metadata", "expected"),
    [
        (PAIRED_METADATA.replace("P03,Primary", ",Primary"), "blank_pair_id"),
        (PAIRED_METADATA.replace("P03_Metastasis,P03,Metastasis\n", ""), "pair_missing_numerator"),
        (PAIRED_METADATA + "P01_Metastasis_2,P01,Metastasis\n", "duplicate_pair_numerator"),
        (PAIRED_METADATA + "P01_Primary_2,P01,Primary\n", "duplicate_pair_denominator"),
    ],
)
def test_invalid_pair_structures_are_rejected_actionably(project_factory, metadata, expected):
    report = validate_project(paired_project(project_factory, metadata=metadata))
    assert not report.is_valid
    assert expected in issue_codes(report)


def test_missing_pair_column_and_pair_factor_identity_are_rejected(project_factory):
    missing = PAIRED_METADATA.replace(",patient,condition", ",subject,condition")
    assert "missing_design_variables" in issue_codes(
        validate_project(paired_project(project_factory, metadata=missing))
    )

    config = paired_config()
    config["design"] = {"type": "paired_two_group", "formula": "~ condition", "pair_id": "condition"}
    report = validate_project(paired_project(project_factory, config=config))
    assert "pair_id_is_contrast_factor" in issue_codes(report)


def test_wrong_pair_formula_and_rank_deficiency_fail_before_execution(project_factory):
    config = paired_config()
    config["design"]["formula"] = "~ condition"
    report = validate_project(paired_project(project_factory, config=config))
    assert "invalid_project_config" in issue_codes(report)
    assert "must be present in design.formula" in report.errors[0].message

    config = paired_config()
    config["design"]["formula"] = "~ patient + batch + condition"
    metadata = PAIRED_METADATA.replace(
        "sample_id,patient,condition", "sample_id,patient,batch,condition"
    ).replace(",P01,", ",P01,B1,").replace(",P02,", ",P02,B2,").replace(",P03,", ",P03,B3,")
    assert "rank_deficient_design" in issue_codes(
        validate_project(paired_project(project_factory, metadata=metadata, config=config))
    )


def test_one_pair_is_insufficient_for_paired_l2(project_factory):
    metadata = "sample_id,patient,condition\nP1P,P1,Primary\nP1M,P1,Metastasis\n"
    report = validate_project(paired_project(project_factory, metadata=metadata))
    assert "insufficient_complete_pairs" in issue_codes(report)


def test_planning_and_frozen_contract_preserve_pairing(project_factory):
    root = project_factory(
        config=paired_config(), counts=PAIRED_COUNTS,
        metadata=PAIRED_METADATA, contrasts=PAIRED_CONTRAST,
    )
    report = validate_project(root)
    generate_plan(report)
    first = (root / "planning" / "manifest.preview.yaml").read_bytes()
    generate_plan(report)
    assert (root / "planning" / "manifest.preview.yaml").read_bytes() == first
    manifest = yaml.safe_load(first)
    assert manifest["design"]["pair_id"] == "patient"
    assert manifest["design"]["pairing"]["contrasts"][0]["complete_pairs"] == 3
    assert manifest["design"]["pairing"]["contrasts"][0]["membership"][0] == {
        "pair_id": "P01", "numerator_sample_id": "P01_Metastasis",
        "denominator_sample_id": "P01_Primary",
    }

    run = create_case_run(report, "M5A", moment=datetime(2026, 9, 10, 10, 0, 0))
    frozen = freeze_case_inputs(report, run, profile="local", command=["rnaseq", "run"])
    contract = json.loads(frozen.contract.read_text())
    assert contract["design"]["type"] == "paired_two_group"
    assert contract["design"]["pairing"]["contrasts"][0]["complete_pairs"] == 3


def test_run_provenance_preserves_paired_design(monkeypatch, project_factory, tmp_path):
    report = validate_project(project_factory(
        config=paired_config(), counts=PAIRED_COUNTS,
        metadata=PAIRED_METADATA, contrasts=PAIRED_CONTRAST,
    ))
    run = create_case_run(report, "M5A-PROVENANCE", moment=datetime(2026, 9, 10, 11, 0, 0))
    freeze_case_inputs(report, run, profile="local", command=["rnaseq", "run"])
    monkeypatch.setattr(service, "check_nextflow", lambda: RuntimeCheck("Nextflow", "FOUND", "test"))
    monkeypatch.setattr(service, "inspect_container_image", lambda _image: {"reference": "test"})
    monkeypatch.setattr(
        service,
        "runtime_snapshot",
        lambda _image: SimpleNamespace(
            host_os="test", host_architecture="test", logical_cpus=1,
            host_memory_bytes=1, docker_architecture="test", docker_memory_bytes=1,
            docker_version="test", control_plane_image_architecture="test",
        ),
    )
    workspace = SimpleNamespace(root=tmp_path / "runtime", launch_dir=tmp_path / "launch", work_dir=tmp_path / "work")
    provenance = service._provenance(
        report, run, profile="local", command=["rnaseq", "run"], workspace=workspace
    )
    assert provenance["design"]["type"] == "paired_two_group"
    assert provenance["design"]["formula"] == "~ patient + condition"
    assert provenance["design"]["pair_id"] == "patient"
    paired = provenance["design"]["pairing"]["contrasts"][0]
    assert (paired["complete_pairs"], paired["samples_in_complete_pairs"]) == (3, 6)


def test_l2_receives_declared_formula_pair_id_and_contrast_direction(project_factory):
    report = validate_project(project_factory(
        config=paired_config(), counts=PAIRED_COUNTS,
        metadata=PAIRED_METADATA, contrasts=PAIRED_CONTRAST,
    ))
    prepared = prepare_l2(report, run_id=None)
    assert prepared.config.design.formula == "~ patient + condition"
    assert prepared.config.design.pair_id == "patient"
    contrast = prepared.contrasts[0]
    assert (contrast.numerator, contrast.denominator) == ("Metastasis", "Primary")
    assert prepared.l1.config["pair_id"] == "patient"


def test_paired_fastq_contract_is_independent_of_paired_end_layout(tmp_path: Path):
    root = tmp_path / "fastq"
    (root / "input" / "fastq").mkdir(parents=True)
    for sample in ("P01_Primary", "P01_Metastasis", "P02_Primary", "P02_Metastasis"):
        (root / "input" / "fastq" / f"{sample}_R1.fastq.gz").write_bytes(b"fixture")
    config = paired_config()
    config["input"] = {"type": "fastq", "path": "input/fastq", "layout": "single_end"}
    config["upstream"] = {
        "engine": "nfcore_rnaseq", "pipeline_version": "3.26.0", "aligner": None,
        "strandedness": "auto", "quantification": {"method": "salmon"},
    }
    config["reference"] = {"source": "igenomes", "genome": "fixture"}
    metadata = "\n".join(PAIRED_METADATA.splitlines()[:5]) + "\n"
    (root / "project.yaml").write_text(yaml.safe_dump(config, sort_keys=False))
    (root / "metadata.csv").write_text(metadata)
    (root / "contrasts.csv").write_text(PAIRED_CONTRAST)
    report = validate_project(root)
    assert report.is_valid
    assert report.fastq.layout.value == "single_end"
    assert report.config.design.type.value == "paired_two_group"


def test_interim_paired_spelling_remains_read_compatible(project_factory):
    config = paired_config()
    config["design"] = {
        "type": "paired", "formula": "~ patient + condition", "pairing_column": "patient"
    }
    report = validate_project(paired_project(project_factory, config=config))
    assert report.is_valid
    assert report.config.design.type.value == "paired_two_group"
    assert report.config.design.pair_id == "patient"


def test_new_scaffolds_explicit_paired_contract_without_biological_data(tmp_path):
    result = runner.invoke(app, [
        "new", "--name", "paired-scaffold", "--destination", str(tmp_path),
        "--species", "human", "--input-type", "raw_counts", "--preset", "L2",
        "--design-type", "paired_two_group", "--pair-id", "donor",
        "--scaffold", "--yes",
    ])
    assert result.exit_code == 0, result.output
    config = yaml.safe_load((tmp_path / "paired-scaffold" / "project.yaml").read_text())
    assert config["design"] == {
        "type": "paired_two_group", "formula": "~ donor + condition", "pair_id": "donor"
    }
    assert (tmp_path / "paired-scaffold" / "metadata.csv").read_text() == "sample_id,donor,condition\n"

    missing = runner.invoke(app, [
        "new", "--name", "implicit-pairing", "--destination", str(tmp_path),
        "--species", "human", "--input-type", "raw_counts", "--preset", "L2",
        "--design-type", "paired_two_group", "--scaffold", "--yes",
    ])
    assert missing.exit_code == 1
    assert "explicit --pair-id" in missing.output
    assert not (tmp_path / "implicit-pairing").exists()
