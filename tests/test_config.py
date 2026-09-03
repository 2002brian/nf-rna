from __future__ import annotations

from copy import deepcopy

from conftest import base_config
from rnaseq.project import load_project
from rnaseq.validators import validate_project


def issue_messages(report):
    return "\n".join(issue.message for issue in report.errors)


def test_schema_version_string_passes(project_factory):
    root = project_factory()
    assert load_project(root).config.schema_version == "1.0"
    assert validate_project(root).is_valid


def test_missing_schema_version_fails(project_factory):
    config = base_config()
    del config["schema_version"]
    report = validate_project(project_factory(config=config))
    assert not report.is_valid
    assert "schema_version" in issue_messages(report)
    assert "Field required" in issue_messages(report)


def test_unsupported_schema_version_fails_clearly(project_factory):
    config = base_config()
    config["schema_version"] = "2.0"
    report = validate_project(project_factory(config=config))
    assert not report.is_valid
    assert "Unsupported project schema version: 2.0" in issue_messages(report)
    assert "Supported schema versions: 1.0, 1.1" in issue_messages(report)


def test_numeric_schema_version_fails(project_factory):
    config = base_config()
    config["schema_version"] = 1.0
    report = validate_project(project_factory(config=config))
    assert not report.is_valid
    assert "schema_version" in issue_messages(report)
    assert "valid string" in issue_messages(report)


def test_unknown_nested_key_fails(project_factory):
    config = deepcopy(base_config())
    config["input"]["mystery"] = True
    report = validate_project(project_factory(config=config))
    assert not report.is_valid
    assert "input.mystery" in issue_messages(report)


def test_invalid_preset_fails(project_factory):
    config = deepcopy(base_config())
    config["project"]["preset"] = "L3"
    assert not validate_project(project_factory(config=config)).is_valid


def test_fastq_requires_fastq_contract(project_factory):
    config = deepcopy(base_config())
    config["input"]["type"] = "fastq"
    assert not validate_project(project_factory(config=config)).is_valid


def test_unsupported_design_type_fails(project_factory):
    config = deepcopy(base_config())
    config["design"]["type"] = "time_course"
    assert not validate_project(project_factory(config=config)).is_valid


def test_path_escape_fails(project_factory):
    config = deepcopy(base_config())
    config["metadata_file"] = "../metadata.csv"
    report = validate_project(project_factory(config=config))
    assert not report.is_valid
    assert "inside the project directory" in issue_messages(report)


def test_public_enrichment_accepts_only_gsea_and_normalizes_the_complete_legacy_scope(project_factory):
    config = base_config()
    config["schema_version"] = "1.1"
    config["annotation"] = {"organism": "Mus musculus", "input_id_type": "ENSEMBL"}
    config["analysis"] = {"enrichment": "gsea"}
    assert load_project(project_factory(config=config)).config.analysis.enrichment == ("gsea",)

    legacy = deepcopy(config)
    legacy["analysis"] = {"enrichment": ["go", "gsea-go", "kegg", "gsea-kegg"]}
    assert load_project(project_factory(config=legacy)).config.analysis.enrichment == ("gsea",)

    partial = deepcopy(config)
    partial["analysis"] = {"enrichment": ["gsea-go"]}
    report = validate_project(project_factory(config=partial))
    assert not report.is_valid
    assert "Legacy analysis.enrichment is accepted only" in issue_messages(report)
