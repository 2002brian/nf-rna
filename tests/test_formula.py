from __future__ import annotations

from copy import deepcopy

import pytest

from conftest import base_config
from rnaseq.validators import validate_project


def codes(report):
    return {issue.code for issue in report.issues}


def test_supported_additive_formula(project_factory):
    config = deepcopy(base_config())
    config["design"]["formula"] = "~ batch + condition"
    metadata = """sample_id,batch,condition
C1,B1,Control
C2,B1,Control
C3,B2,Control
T1,B1,Treatment
T2,B2,Treatment
T3,B2,Treatment
"""
    report = validate_project(project_factory(config=config, metadata=metadata))
    assert report.is_valid
    assert report.formula_variables == ("batch", "condition")


@pytest.mark.parametrize(
    "formula",
    [
        "~ condition * batch",
        "~ condition:batch",
        "~ condition / batch",
        "~ condition^2",
        "~ log(age) + condition",
        "~ 0 + condition",
    ],
)
def test_unsupported_formula(project_factory, formula):
    config = deepcopy(base_config())
    config["design"]["formula"] = formula
    report = validate_project(project_factory(config=config))
    assert not report.is_valid
    assert "unsupported_formula" in codes(report)


def test_multi_group_requires_three_levels(project_factory):
    config = deepcopy(base_config())
    config["design"]["type"] = "multi_group"
    report = validate_project(project_factory(config=config))
    assert not report.is_valid
    assert "invalid_multi_group_design" in codes(report)


def test_multi_group_with_three_levels_passes(project_factory):
    config = deepcopy(base_config())
    config["design"]["type"] = "multi_group"
    counts = "gene_id,A1,A2,A3,B1,B2,B3,C1,C2,C3\nGeneA,1,2,3,4,5,6,7,8,9\n"
    metadata = """sample_id,condition
A1,A
A2,A
A3,A
B1,B
B2,B
B3,B
C1,C
C2,C
C3,C
"""
    contrasts = "contrast_id,factor,numerator,denominator\nB_vs_A,condition,B,A\n"
    report = validate_project(
        project_factory(config=config, counts=counts, metadata=metadata, contrasts=contrasts)
    )
    assert report.is_valid


def test_paired_design_uses_explicit_pair_id_first(project_factory):
    config = deepcopy(base_config())
    config["design"] = {"type": "paired_two_group", "formula": "~ subject + condition", "pair_id": "subject"}
    counts = "gene_id,S1C,S1T,S2C,S2T\nGeneA,1,2,3,4\n"
    metadata = "sample_id,subject,condition\nS1C,S1,C\nS1T,S1,T\nS2C,S2,C\nS2T,S2,T\n"
    contrasts = "contrast_id,factor,numerator,denominator\nT_vs_C,condition,T,C\n"
    report = validate_project(project_factory(config=config, counts=counts, metadata=metadata, contrasts=contrasts))
    assert report.is_valid
    assert report.config.design.pair_id == "subject"


def test_paired_design_pairing_identity_is_independent_of_covariate_order(project_factory):
    config = deepcopy(base_config())
    config["design"] = {
        "type": "paired_two_group", "formula": "~ batch + subject + condition", "pair_id": "subject",
    }
    counts = "gene_id,S1C,S1T,S2C,S2T\nGeneA,1,2,3,4\n"
    metadata = (
        "sample_id,batch,subject,condition\n"
        "S1C,B1,S1,C\nS1T,B2,S1,T\nS2C,B2,S2,C\nS2T,B1,S2,T\n"
    )
    contrasts = "contrast_id,factor,numerator,denominator\nT_vs_C,condition,T,C\n"
    report = validate_project(project_factory(config=config, counts=counts, metadata=metadata, contrasts=contrasts))
    assert report.is_valid
    assert report.config.design.pair_id == "subject"


@pytest.mark.parametrize(
    "metadata",
    (
        "sample_id,subject,condition\nS1C,S1,C\nS1T,S1,T\nS2C,S2,C\n",
        "sample_id,subject,condition\nS1C,S1,C\nS1C2,S1,C\nS2C,S2,C\nS2T,S2,T\n",
    ),
)
def test_paired_design_rejects_missing_or_duplicate_pair_condition(project_factory, metadata):
    config = deepcopy(base_config())
    config["design"] = {"type": "paired_two_group", "formula": "~ subject + condition", "pair_id": "subject"}
    samples = [line.split(",", 1)[0] for line in metadata.splitlines()[1:]]
    counts = "gene_id," + ",".join(samples) + "\nGeneA," + ",".join("1" for _ in samples) + "\n"
    contrasts = "contrast_id,factor,numerator,denominator\nT_vs_C,condition,T,C\n"
    report = validate_project(project_factory(config=config, counts=counts, metadata=metadata, contrasts=contrasts))
    assert not report.is_valid
    assert {
        "pair_missing_numerator", "pair_missing_denominator",
        "duplicate_pair_numerator", "duplicate_pair_denominator",
    } & codes(report)


def test_rank_deficient_additive_design_is_rejected_actionably(project_factory):
    config = deepcopy(base_config())
    config["design"] = {
        "type": "paired_two_group", "formula": "~ subject + batch + condition", "pair_id": "subject",
    }
    counts = "gene_id,S1C,S1T,S2C,S2T\nGeneA,1,2,3,4\n"
    metadata = (
        "sample_id,subject,batch,condition\n"
        "S1C,S1,B1,C\nS1T,S1,B1,T\nS2C,S2,B2,C\nS2T,S2,B2,T\n"
    )
    contrasts = "contrast_id,factor,numerator,denominator\nT_vs_C,condition,T,C\n"
    report = validate_project(project_factory(config=config, counts=counts, metadata=metadata, contrasts=contrasts))
    assert "rank_deficient_design" in codes(report)
    assert any("confounded or redundant covariate" in issue.message for issue in report.errors)


def test_valid_unpaired_design_remains_valid(project_factory):
    report = validate_project(project_factory())
    assert report.is_valid
    assert report.config.design.pair_id is None
