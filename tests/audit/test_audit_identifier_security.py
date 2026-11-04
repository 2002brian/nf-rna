"""Security regressions: sample and contrast identifiers (audit D-02, D-06).

FASTQ-derived sample IDs are interpolated into the HISAT2 + featureCounts task
shell scripts and become output file names; contrast IDs become R output
directories (``output_dir/contrasts/<contrast_id>``).  Both must be restricted
to a path- and shell-inert alphabet before anything is planned or run.
"""

from __future__ import annotations

import re
from pathlib import Path

import pytest

from conftest import BASE_CONTRASTS
from rnaseq.design_setup import CONTRAST_ID_PATTERN
from rnaseq.models import SequencingLayout
from rnaseq.validators import SAFE_IDENTIFIER_PATTERN, ValidationReport, validate_fastq, validate_project

WORKFLOW = Path(__file__).resolve().parents[2] / "workflow" / "hisat2_featurecounts.nf"

MALICIOUS_SAMPLE_IDS = [
    "a b", "a\tb", "a$(id)", "a`id`", "a;b", "a&b", "a|b", "a>b", "a<b", "a'b", 'a"b', "a\\b",
    "a*b", "a?b", "a[b]", "a{b}", "a~b", "a!b", "a#b", "a%b", "a=b", "a,b", "a:b", "aé", "樣本1",
    "-rf", ".hidden", "..", "_S1", "$HOME", "a b",
]
SAFE_SAMPLE_IDS = ["S1", "Ctrl-1", "WT.rep_1", "Undetermined_S0", "4T1_KO_3", "A1.b-2_c"]


def _fastq_dir(tmp_path: Path, names: list[str]) -> Path:
    root = tmp_path / "fastq"
    root.mkdir()
    for name in names:
        (root / name).write_bytes(b"")
    return root


def _codes(report) -> set[str]:
    return {issue.code for issue in report.issues}


@pytest.mark.parametrize("sample", MALICIOUS_SAMPLE_IDS)
@pytest.mark.parametrize("layout", [SequencingLayout.PAIRED_END, SequencingLayout.SINGLE_END])
def test_shell_or_path_unsafe_fastq_sample_ids_are_rejected(tmp_path, sample, layout):
    reads = ["R1", "R2"] if layout is SequencingLayout.PAIRED_END else ["R1"]
    directory = _fastq_dir(tmp_path, [f"{sample}_L001_{read}_001.fastq.gz" for read in reads])
    report = ValidationReport(project_dir=tmp_path)
    validate_fastq(directory, layout, report)
    assert "unsafe_fastq_sample_id" in _codes(report)
    assert not report.is_valid


@pytest.mark.parametrize("sample", SAFE_SAMPLE_IDS)
def test_conventional_fastq_sample_ids_stay_valid(tmp_path, sample):
    directory = _fastq_dir(tmp_path, [f"{sample}_R1.fastq.gz", f"{sample}_R2.fastq.gz"])
    report = ValidationReport(project_dir=tmp_path)
    summary = validate_fastq(directory, SequencingLayout.PAIRED_END, report)
    assert not report.issues
    assert summary.sample_ids == (sample,)


def test_unsafe_fastq_sample_blocks_the_whole_project(tmp_path):
    """End to end: a malicious FASTQ name makes the project invalid, so no plan or run is possible."""

    from test_runtime_dispatch import _salmon_report

    root = tmp_path / "project"
    _salmon_report(root)
    for read in ("R1", "R2"):
        (root / "input" / "fastq" / f"T2_{read}.fastq.gz").rename(root / "input" / "fastq" / f"T2$(touch pwned)_{read}.fastq.gz")
    metadata = (root / "metadata.csv").read_text(encoding="utf-8").replace("T2,", "T2$(touch pwned),")
    (root / "metadata.csv").write_text(metadata, encoding="utf-8")
    report = validate_project(root)
    assert "unsafe_fastq_sample_id" in _codes(report)
    assert not report.is_valid
    from rnaseq.planner import generate_plan

    with pytest.raises(Exception):
        generate_plan(report)
    assert not (root / "pwned").exists()


def test_workflow_refuses_unsafe_samplesheet_ids_with_the_same_rule():
    """Defence in depth: the HISAT2 workflow re-checks every samplesheet sample ID."""

    text = WORKFLOW.read_text(encoding="utf-8")
    match = re.search(r"if \(!\(row\.sample ==~ /(?P<regex>[^/]+)/\)\) \{\s*error ", text)
    assert match is not None, "hisat2_featurecounts.nf must validate row.sample before building the reads channel"
    groovy = re.compile(match.group("regex"))  # Groovy ==~ is a whole-string match
    assert f"^{match.group('regex')}$" == SAFE_IDENTIFIER_PATTERN.pattern
    for sample in MALICIOUS_SAMPLE_IDS:
        assert groovy.fullmatch(sample) is None, sample
    for sample in SAFE_SAMPLE_IDS:
        assert groovy.fullmatch(sample) is not None, sample
    # The guard runs where samplesheet rows become the reads channel, before any process.
    assert text.index("row.sample ==~") < text.index("FASTQC_RAW(reads)")


@pytest.mark.parametrize("contrast_id", [
    "../../escape", "..", ".", "a/b", "/abs", "a\\b", "with space", "a$(id)", "a;b", "-x", ".hidden", "Té_vs_C",
])
def test_path_unsafe_contrast_ids_are_rejected(project_factory, contrast_id):
    contrasts = f"contrast_id,factor,numerator,denominator\n{contrast_id},condition,Treatment,Control\n"
    report = validate_project(project_factory(contrasts=contrasts))
    assert "unsafe_contrast_id" in _codes(report)
    assert not report.is_valid


def test_generated_and_validated_contrast_ids_use_one_rule(project_factory):
    assert CONTRAST_ID_PATTERN is SAFE_IDENTIFIER_PATTERN
    assert validate_project(project_factory(contrasts=BASE_CONTRASTS)).is_valid
