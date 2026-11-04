"""Audit regressions: inputs that pass Python validation but break later.

Each ``xfail(strict=True)`` test states the *safe* expected behaviour and
currently fails, documenting a confirmed defect without changing production
code.  When the defect is fixed the test XPASSes, which strict mode reports as
a failure so the marker is removed deliberately.  Synthetic fixtures only.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from conftest import BASE_CONTRASTS, BASE_COUNTS, BASE_METADATA
from rnaseq.models import SequencingLayout
from rnaseq.validators import ValidationReport, validate_fastq, validate_project


def _codes(report) -> set[str]:
    return {issue.code for issue in report.issues}


# --------------------------------------------------------------- D-02 contrast_id becomes a path component

@pytest.mark.xfail(strict=True, reason="AUDIT D-02: contrast_id is not restricted; R writes output_dir/contrasts/<contrast_id>")
@pytest.mark.parametrize("contrast_id", ["../../escape", "a/b", "with space"])
def test_unsafe_contrast_id_is_rejected(project_factory, contrast_id):
    contrasts = f"contrast_id,factor,numerator,denominator\n{contrast_id},condition,Treatment,Control\n"
    report = validate_project(project_factory(contrasts=contrasts))
    assert not report.is_valid


def test_generated_contrast_ids_are_already_safe(project_factory):
    """Control: the conventional ID used by `rnaseq design` validates."""

    assert validate_project(project_factory(contrasts=BASE_CONTRASTS)).is_valid


# --------------------------------------------------------------- D-03 / D-04 literal "NA" is missing in R

@pytest.mark.xfail(strict=True, reason='AUDIT D-03: design_metadata.R reads literal "NA" as missing; Python accepts it as a level')
def test_literal_na_design_level_is_flagged(project_factory):
    metadata = BASE_METADATA.replace("Treatment", "NA")
    contrasts = "contrast_id,factor,numerator,denominator\nNA_vs_Control,condition,NA,Control\n"
    report = validate_project(project_factory(metadata=metadata, contrasts=contrasts))
    assert report.errors or report.warnings


@pytest.mark.xfail(strict=True, reason='AUDIT D-04: sample_id "NA" becomes an NA row name in R')
def test_literal_na_sample_id_is_flagged(project_factory):
    report = validate_project(project_factory(
        counts=BASE_COUNTS.replace("T3", "NA"), metadata=BASE_METADATA.replace("T3", "NA"),
    ))
    assert report.errors or report.warnings


# --------------------------------------------------------------- D-05 Python numeric grammar is wider than R's

@pytest.mark.xfail(strict=True, reason="AUDIT D-05: Decimal('1_000') is accepted; R read.csv makes the column character -> NA counts")
@pytest.mark.parametrize("value", ["1_000", "١٢"])  # PEP 515 underscore, Arabic-Indic digits
def test_count_values_r_cannot_parse_are_rejected(project_factory, value):
    report = validate_project(project_factory(counts=BASE_COUNTS.replace(",40,", f",{value},")))
    assert not report.is_valid


def test_ordinary_integral_spellings_remain_valid(project_factory):
    """Control: forms both Python and R parse identically must stay valid."""

    report = validate_project(project_factory(counts=BASE_COUNTS.replace(",40,", ",4e1,")))
    assert report.is_valid


# --------------------------------------------------------------- D-06 sample IDs reach unquoted shell scripts

def _fastq_dir(tmp_path: Path, names: list[str]) -> Path:
    root = tmp_path / "fastq"
    root.mkdir()
    for name in names:
        (root / name).write_bytes(b"")
    return root


@pytest.mark.xfail(strict=True, reason="AUDIT D-06: FASTQ-derived sample IDs are interpolated unquoted into hisat2_featurecounts.nf scripts")
@pytest.mark.parametrize("sample", ["a b", "a$(id)", "a;b", "a'b"])
def test_shell_unsafe_fastq_sample_ids_are_rejected(tmp_path, sample):
    directory = _fastq_dir(tmp_path, [f"{sample}_R1.fastq.gz", f"{sample}_R2.fastq.gz"])
    report = ValidationReport(project_dir=tmp_path)
    validate_fastq(directory, SequencingLayout.PAIRED_END, report)
    assert report.errors


# --------------------------------------------------------------- R-03 merged + per-lane FASTQs for one sample

@pytest.mark.xfail(strict=True, reason="AUDIT R-03: a merged FASTQ next to its lane files is silently counted twice")
def test_merged_and_lane_fastqs_for_one_sample_are_flagged(tmp_path):
    directory = _fastq_dir(tmp_path, [
        "S1_L001_R1_001.fastq.gz", "S1_L001_R2_001.fastq.gz",  # lane file
        "S1_R1.fastq.gz", "S1_R2.fastq.gz",                    # merged file of the same sample
    ])
    report = ValidationReport(project_dir=tmp_path)
    summary = validate_fastq(directory, SequencingLayout.PAIRED_END, report)
    assert summary is not None
    assert report.errors or report.warnings


def test_true_multi_lane_sample_stays_valid(tmp_path):
    """Control: genuine technical lanes of one sample are accepted and grouped."""

    directory = _fastq_dir(tmp_path, [
        "S1_L001_R1_001.fastq.gz", "S1_L001_R2_001.fastq.gz",
        "S1_L002_R1_001.fastq.gz", "S1_L002_R2_001.fastq.gz",
    ])
    report = ValidationReport(project_dir=tmp_path)
    summary = validate_fastq(directory, SequencingLayout.PAIRED_END, report)
    assert not report.issues
    assert [(record.sample_id, record.lane) for record in summary.records] == [("S1", "001"), ("S1", "002")]
