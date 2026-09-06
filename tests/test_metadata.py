from __future__ import annotations

from copy import deepcopy

from conftest import base_config
from rnaseq.validators import validate_project


def codes(report):
    return {issue.code for issue in report.issues}


def test_unused_metadata_columns_are_permitted(project_factory):
    metadata = """sample_id,condition,sex,batch,age
C1,Control,M,B1,8
C2,Control,F,B1,9
C3,Control,M,B1,8
T1,Treatment,M,B2,9
T2,Treatment,F,B2,8
T3,Treatment,M,B2,10
"""
    root = project_factory(metadata=metadata)
    before = (root / "metadata.csv").read_bytes()
    report = validate_project(root)
    after = (root / "metadata.csv").read_bytes()
    assert report.is_valid
    assert report.state == "PASS"
    assert not any("unused" in issue.message.lower() for issue in report.issues)
    assert report.metadata.columns == ("sample_id", "condition", "sex", "batch", "age")
    assert before == after


def test_missing_required_formula_variable(project_factory):
    metadata = "sample_id,batch\nC1,B1\nC2,B1\nC3,B1\nT1,B2\nT2,B2\nT3,B2\n"
    report = validate_project(project_factory(metadata=metadata))
    assert "missing_design_variables" in codes(report)


def test_missing_required_formula_value(project_factory):
    metadata = "sample_id,condition\nC1,Control\nC2,\nC3,Control\nT1,Treatment\nT2,Treatment\nT3,Treatment\n"
    report = validate_project(project_factory(metadata=metadata))
    assert "missing_design_value" in codes(report)


def test_duplicate_additional_column(project_factory):
    metadata = "sample_id,condition,sex,sex\nC1,Control,M,M\n"
    report = validate_project(project_factory(metadata=metadata))
    assert "duplicate_metadata_column" in codes(report)


def test_blank_additional_column(project_factory):
    metadata = "sample_id,condition,,age\nC1,Control,M,8\n"
    report = validate_project(project_factory(metadata=metadata))
    assert "blank_metadata_column" in codes(report)


def test_duplicate_and_blank_sample_ids(project_factory):
    metadata = "sample_id,condition\nC1,Control\nC1,Control\n,Treatment\nT1,Treatment\nT2,Treatment\nT3,Treatment\n"
    report = validate_project(project_factory(metadata=metadata))
    assert "duplicate_sample_id" in codes(report)
    assert "blank_sample_id" in codes(report)


def test_exact_sample_agreement_ignores_order(project_factory):
    metadata = "sample_id,condition\nT3,Treatment\nC2,Control\nT1,Treatment\nC1,Control\nT2,Treatment\nC3,Control\n"
    assert validate_project(project_factory(metadata=metadata)).is_valid


def test_sample_set_mismatches_are_exact(project_factory):
    metadata = "sample_id,condition\nC1,Control\nC2,Control\nC3,Control\nT1,Treatment\nT2,Treatment\nS06,Treatment\n"
    report = validate_project(project_factory(metadata=metadata))
    messages = "\n".join(issue.message for issue in report.errors)
    assert "T3" in messages
    assert "S06" in messages
    assert "count_only_samples" in codes(report)
    assert "metadata_only_samples" in codes(report)


def test_n_two_is_warning_not_error(project_factory):
    counts = "gene_id,C1,C2,T1,T2\nGeneA,1,2,3,4\n"
    metadata = "sample_id,condition\nC1,Control\nC2,Control\nT1,Treatment\nT2,Treatment\n"
    report = validate_project(project_factory(counts=counts, metadata=metadata))
    assert report.is_valid
    assert report.state == "PASS WITH WARNINGS"
    assert sum(issue.code == "limited_replication" for issue in report.warnings) == 2


def test_strict_paired_design_passes(project_factory):
    config = deepcopy(base_config())
    config["design"] = {"type": "paired", "formula": "~ subject_id + condition", "pairing_column": "subject_id"}
    counts = "gene_id,S1_pre,S1_post,S2_pre,S2_post\nGeneA,1,2,3,4\n"
    metadata = "sample_id,subject_id,condition\nS1_pre,S1,Pre\nS1_post,S1,Post\nS2_pre,S2,Pre\nS2_post,S2,Post\n"
    contrasts = "contrast_id,factor,numerator,denominator\nPost_vs_Pre,condition,Post,Pre\n"
    report = validate_project(
        project_factory(config=config, counts=counts, metadata=metadata, contrasts=contrasts)
    )
    assert report.is_valid


def test_incomplete_pair_fails(project_factory):
    config = deepcopy(base_config())
    config["design"] = {"type": "paired", "formula": "~ subject_id + condition", "pairing_column": "subject_id"}
    counts = "gene_id,S1_pre,S1_post,S2_pre\nGeneA,1,2,3\n"
    metadata = "sample_id,subject_id,condition\nS1_pre,S1,Pre\nS1_post,S1,Post\nS2_pre,S2,Pre\n"
    contrasts = "contrast_id,factor,numerator,denominator\nPost_vs_Pre,condition,Post,Pre\n"
    report = validate_project(
        project_factory(config=config, counts=counts, metadata=metadata, contrasts=contrasts)
    )
    assert not report.is_valid
    assert "incomplete_pair" in codes(report)
