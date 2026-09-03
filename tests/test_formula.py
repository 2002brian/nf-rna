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
