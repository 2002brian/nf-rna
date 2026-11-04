"""Python validation must not accept values that R reads differently (audit D-03, D-04, D-05).

R's ``read.csv`` reads a literal ``NA`` as missing and parses numbers with a
plain ASCII grammar; Python's ``Decimal``/``float`` also accept PEP 515
underscores and non-ASCII digits.
"""

from __future__ import annotations

import pytest
import yaml

from conftest import BASE_COUNTS, BASE_METADATA, base_config
from rnaseq.validators import validate_project


def _codes(report) -> set[str]:
    return {issue.code for issue in report.issues}


def test_literal_na_gene_id_is_rejected(project_factory):
    lines = BASE_COUNTS.splitlines()
    lines[1] = "NA," + lines[1].split(",", 1)[1]
    report = validate_project(project_factory(counts="\n".join(lines) + "\n"))
    assert "r_missing_gene_id" in _codes(report)


def test_literal_na_in_an_unused_metadata_column_is_harmless(project_factory):
    metadata = BASE_METADATA.replace("sample_id,condition", "sample_id,condition,note", 1)
    metadata = "\n".join(line + (",NA" if index else "") for index, line in enumerate(metadata.strip().splitlines())) + "\n"
    report = validate_project(project_factory(metadata=metadata))
    assert report.is_valid, report.errors


def test_lowercase_or_padded_na_is_an_ordinary_level_for_r_too(project_factory):
    """R's default na.strings is exactly "NA"; other spellings stay text in both languages."""

    report = validate_project(project_factory(metadata=BASE_METADATA.replace("Treatment", "na"),
                                              contrasts="contrast_id,factor,numerator,denominator\nna_vs_Control,condition,na,Control\n"))
    assert "r_missing_design_value" not in _codes(report)


def test_literal_na_design_value_reports_a_specific_code(project_factory):
    metadata = BASE_METADATA.replace("Treatment", "NA")
    contrasts = "contrast_id,factor,numerator,denominator\nNA_vs_Control,condition,NA,Control\n"
    report = validate_project(project_factory(metadata=metadata, contrasts=contrasts))
    assert "r_missing_design_value" in _codes(report)
    assert not report.is_valid


@pytest.mark.parametrize("value", ["1_000", "4_0", "١٢", "４０", "40 ", " 40", "0x28"])
def test_count_spellings_r_cannot_parse_are_rejected(project_factory, value):
    report = validate_project(project_factory(counts=BASE_COUNTS.replace(",40,", f",{value},")))
    assert "non_numeric_count" in _codes(report)


@pytest.mark.parametrize("value", ["40", "+40", "4e1", "4E+1", "40.0", "40.", "400e-1", " 40", "40 "])
def test_count_spellings_both_languages_parse_identically_stay_valid(project_factory, value):
    report = validate_project(project_factory(counts=BASE_COUNTS.replace(",40,", f",{value},")))
    assert report.is_valid, report.errors


def _continuous_project(project_factory, value):
    config = base_config()
    config["design"] = {"type": "two_group", "formula": "~ age + condition", "variables": {"age": "continuous", "condition": "categorical"}}
    metadata = f"sample_id,age,condition\nC1,30,Control\nC2,{value},Control\nC3,32,Control\nT1,33,Treatment\nT2,34,Treatment\nT3,35,Treatment\n"
    return validate_project(project_factory(config=config, metadata=metadata))


@pytest.mark.parametrize("value", ["1_000", "١٢", "3 "])
def test_continuous_values_r_cannot_parse_are_rejected(project_factory, value):
    assert "invalid_continuous_value" in _codes(_continuous_project(project_factory, value))


@pytest.mark.parametrize("value", ["31", "31.5", "3.1e1", "-2"])
def test_plain_continuous_values_stay_valid(project_factory, value):
    report = _continuous_project(project_factory, value)
    assert "invalid_continuous_value" not in _codes(report), report.errors
