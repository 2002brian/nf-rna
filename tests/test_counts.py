from __future__ import annotations

import pytest

from rnaseq.validators import validate_project


def codes(report):
    return {issue.code for issue in report.issues}


def test_valid_integer_and_integral_decimal_counts(project_factory):
    counts = "gene_id,C1,C2,C3,T1,T2,T3\nGeneA,1,2.0,3,4,5,6\n"
    assert validate_project(project_factory(counts=counts)).is_valid


@pytest.mark.parametrize(
    ("value", "expected_code"),
    [
        ("-1", "negative_count"),
        ("1.5", "fractional_count"),
        ("abc", "non_numeric_count"),
        ("NaN", "non_finite_count"),
        ("Inf", "non_finite_count"),
    ],
)
def test_invalid_count_values(project_factory, value, expected_code):
    counts = f"gene_id,C1,C2,C3,T1,T2,T3\nGeneA,{value},2,3,4,5,6\n"
    report = validate_project(project_factory(counts=counts))
    assert not report.is_valid
    assert expected_code in codes(report)


def test_duplicate_gene_ids(project_factory):
    counts = "gene_id,C1,C2,C3,T1,T2,T3\nGeneA,1,2,3,4,5,6\nGeneA,2,3,4,5,6,7\n"
    report = validate_project(project_factory(counts=counts))
    assert "duplicate_gene_id" in codes(report)


def test_duplicate_sample_headers(project_factory):
    counts = "gene_id,C1,C1,C3,T1,T2,T3\nGeneA,1,2,3,4,5,6\n"
    report = validate_project(project_factory(counts=counts))
    assert "duplicate_count_matrix_column" in codes(report)


def test_blank_sample_header(project_factory):
    counts = "gene_id,C1,,C3,T1,T2,T3\nGeneA,1,2,3,4,5,6\n"
    report = validate_project(project_factory(counts=counts))
    assert "blank_count_matrix_column" in codes(report)


def test_all_zero_genes_warn_without_failure(project_factory):
    counts = "gene_id,C1,C2,C3,T1,T2,T3\nGeneA,0,0,0,0,0,0\nGeneB,1,2,3,4,5,6\n"
    report = validate_project(project_factory(counts=counts))
    assert report.is_valid
    assert report.state == "PASS WITH WARNINGS"
    assert report.counts.all_zero_genes == 1
    assert "all_zero_genes" in codes(report)


def test_malformed_count_row(project_factory):
    counts = "gene_id,C1,C2,C3,T1,T2,T3\nGeneA,1,2\n"
    report = validate_project(project_factory(counts=counts))
    assert "malformed_count_row" in codes(report)
