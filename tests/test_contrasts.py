from __future__ import annotations

import pytest

from rnaseq.validators import validate_project


def codes(report):
    return {issue.code for issue in report.issues}


def test_valid_contrast(project_factory):
    assert validate_project(project_factory()).is_valid


def test_additional_contrast_column_fails(project_factory):
    contrasts = "contrast_id,factor,numerator,denominator,note\nTx,condition,Treatment,Control,x\n"
    report = validate_project(project_factory(contrasts=contrasts))
    assert "invalid_contrast_schema" in codes(report)


def test_wrong_contrast_column_order_fails(project_factory):
    contrasts = "factor,contrast_id,numerator,denominator\ncondition,Tx,Treatment,Control\n"
    report = validate_project(project_factory(contrasts=contrasts))
    assert "invalid_contrast_schema" in codes(report)


@pytest.mark.parametrize(
    ("contrasts", "expected"),
    [
        (
            "contrast_id,factor,numerator,denominator\nTx,missing,Treatment,Control\n",
            "missing_contrast_factor",
        ),
        (
            "contrast_id,factor,numerator,denominator\nTx,condition,Missing,Control\n",
            "missing_contrast_numerator",
        ),
        (
            "contrast_id,factor,numerator,denominator\nTx,condition,Treatment,Missing\n",
            "missing_contrast_denominator",
        ),
        (
            "contrast_id,factor,numerator,denominator\nTx,condition,Control,Control\n",
            "identical_contrast_levels",
        ),
    ],
)
def test_invalid_contrast_values(project_factory, contrasts, expected):
    report = validate_project(project_factory(contrasts=contrasts))
    assert expected in codes(report)


def test_duplicate_contrast_id(project_factory):
    contrasts = """contrast_id,factor,numerator,denominator
Tx,condition,Treatment,Control
Tx,condition,Control,Treatment
"""
    report = validate_project(project_factory(contrasts=contrasts))
    assert "duplicate_contrast_id" in codes(report)
