from __future__ import annotations

import gzip
import hashlib
from copy import deepcopy
from pathlib import Path

import yaml
from typer.testing import CliRunner

from conftest import BASE_CONTRASTS, base_config
from rnaseq.cli import app
from rnaseq.execution import build_nextflow_command
from rnaseq.planner import generate_plan
from rnaseq.validators import validate_project

runner = CliRunner()
READ = b"@read1\nACGT\n+\n!!!!\n"


def _fastq_config(*, layout: str = "paired_end") -> dict:
    config = deepcopy(base_config())
    config["input"] = {"type": "fastq", "path": "input/fastq", "layout": layout}
    config["upstream"] = {
        "engine": "nfcore_rnaseq",
        "pipeline_version": "3.26.0",
        "aligner": None,
        "strandedness": "auto",
    }
    config["reference"] = {"source": "igenomes", "genome": None}
    return config


def _make_fastq_project(
    tmp_path: Path, *, layout: str = "paired_end", files: tuple[str, ...] | None = None
) -> Path:
    root = tmp_path / "fastq_project"
    fastq_dir = root / "input" / "fastq"
    fastq_dir.mkdir(parents=True)
    names = files or ("C1_R1.fastq.gz", "C1_R2.fastq.gz", "C2_R1.fastq.gz", "C2_R2.fastq.gz", "T1_R1.fastq.gz", "T1_R2.fastq.gz", "T2_R1.fastq.gz", "T2_R2.fastq.gz")
    for name in names:
        with gzip.open(fastq_dir / name, "wb") as handle:
            handle.write(READ)
    metadata = "sample_id,condition,sex,batch,age\nC1,Control,M,B1,8\nC2,Control,F,B1,9\nT1,Treatment,M,B2,8\nT2,Treatment,F,B2,9\n"
    (root / "metadata.csv").write_text(metadata, encoding="utf-8")
    (root / "contrasts.csv").write_text(BASE_CONTRASTS, encoding="utf-8")
    (root / "project.yaml").write_text(yaml.safe_dump(_fastq_config(layout=layout), sort_keys=False), encoding="utf-8")
    return root


def test_valid_paired_fastq_with_extensible_metadata(tmp_path):
    report = validate_project(_make_fastq_project(tmp_path))
    assert report.is_valid
    assert report.fastq is not None
    assert report.fastq.sample_ids == ("C1", "C2", "T1", "T2")
    assert report.config is not None and report.config.input.preprocessing.value == "raw"
    assert not any("unused" in issue.message.lower() for issue in report.issues)
    assert not report.execution_ready
    assert report.config.input.layout.value == "paired_end"
    assert report.config.design.type.value == "two_group"
    assert report.config.design.pair_id is None


def test_valid_single_end_fastq(tmp_path):
    files = ("C1_R1.fastq.gz", "C2_R1.fastq.gz", "T1_R1.fastq.gz", "T2_R1.fastq.gz")
    report = validate_project(_make_fastq_project(tmp_path, layout="single_end", files=files))
    assert report.is_valid
    assert all(record.fastq_2 is None for record in report.fastq.records)


def test_fastq_pairing_and_metadata_mismatch_fail(tmp_path):
    root = _make_fastq_project(tmp_path, files=("C1_R1.fastq.gz", "C1_R2.fastq.gz", "C2_R1.fastq.gz", "T1_R2.fastq.gz"))
    report = validate_project(root)
    codes = {issue.code for issue in report.errors}
    assert {"orphan_r1", "orphan_r2", "metadata_only_fastq_samples"} <= codes


def test_duplicate_assignment_and_fastq_only_sample_fail(tmp_path):
    root = _make_fastq_project(tmp_path)
    extra = root / "input" / "fastq"
    with gzip.open(extra / "C1_R1.fq.gz", "wb") as handle:
        handle.write(READ)
    with gzip.open(extra / "Extra_R1.fastq.gz", "wb") as handle:
        handle.write(READ)
    with gzip.open(extra / "Extra_R2.fastq.gz", "wb") as handle:
        handle.write(READ)
    codes = {issue.code for issue in validate_project(root).errors}
    assert {"duplicate_fastq_assignment", "fastq_only_samples"} <= codes


def test_unsupported_fastq_layout_fails_config_validation(tmp_path):
    root = _make_fastq_project(tmp_path)
    config = _fastq_config()
    config["input"]["layout"] = "mixed_end"
    (root / "project.yaml").write_text(yaml.safe_dump(config, sort_keys=False), encoding="utf-8")
    assert "invalid_project_config" in {issue.code for issue in validate_project(root).errors}


def test_multiple_lanes_are_sorted_and_generate_nf_core_samplesheet(tmp_path):
    files = (
        "T1_R2_001.fastq.gz", "T1_R1_001.fastq.gz", "C1_L002_R2_001.fastq.gz", "C1_L001_R1_001.fastq.gz",
        "C1_L002_R1_001.fastq.gz", "C1_L001_R2_001.fastq.gz", "C2_R1.fastq.gz", "C2_R2.fastq.gz",
        "T2_R1.fastq.gz", "T2_R2.fastq.gz",
    )
    root = _make_fastq_project(tmp_path, files=files)
    result = runner.invoke(app, ["plan", str(root)])
    assert result.exit_code == 0, result.output
    samplesheet = (root / "planning" / "samplesheet.csv").read_text(encoding="utf-8")
    assert samplesheet.splitlines()[0] == "sample,fastq_1,fastq_2,strandedness"
    assert samplesheet.splitlines()[1:3] == [
        "C1,input/fastq/C1_L001_R1_001.fastq.gz,input/fastq/C1_L001_R2_001.fastq.gz,auto",
        "C1,input/fastq/C1_L002_R1_001.fastq.gz,input/fastq/C1_L002_R2_001.fastq.gz,auto",
    ]
    paths = [root / "planning" / name for name in ("samplesheet.csv", "analysis_plan.md", "manifest.preview.yaml", "upstream_run.preview.yaml")]
    first = [hashlib.sha256(path.read_bytes()).hexdigest() for path in paths]
    assert runner.invoke(app, ["plan", str(root)]).exit_code == 0
    assert first == [hashlib.sha256(path.read_bytes()).hexdigest() for path in paths]
    preview = yaml.safe_load(paths[-1].read_text(encoding="utf-8"))
    assert preview["preprocessing"] == "raw"
    assert preview["skip_trimming"] is False
    assert preview["nfcore_preprocessing_arguments"] == []
    assert preview["nfcore_runtime_params"] == {"skip_alignment": True}
    assert preview["status"]["execution_ready"] is False
    assert preview["blocking_requirements"] == ["reference genome not configured"]


def test_fastq_preprocessing_controls_skip_trimming_and_plan(tmp_path):
    root = _make_fastq_project(tmp_path)
    config = yaml.safe_load((root / "project.yaml").read_text(encoding="utf-8"))
    config["upstream"]["quantification"] = {"method": "salmon"}
    config["reference"]["genome"] = "test_reference"
    config["input"]["preprocessing"] = "pretrimmed"
    (root / "project.yaml").write_text(yaml.safe_dump(config, sort_keys=False), encoding="utf-8")
    report = validate_project(root)
    assert report.is_valid and report.execution_ready
    generate_plan(report)
    command = build_nextflow_command(
        report,
        samplesheet=root / "planning" / "samplesheet.csv",
        output_dir=root / "runs" / "out",
        profile="local",
    )
    # Runtime booleans live in nfcore.params.json; a direct valueless CLI flag
    # would be serialized by Nextflow as the string "true".
    assert "--skip_trimming" not in command
    preview = yaml.safe_load((root / "planning" / "upstream_run.preview.yaml").read_text(encoding="utf-8"))
    assert preview["preprocessing"] == "pretrimmed"
    assert preview["skip_trimming"] is True
    assert preview["nfcore_preprocessing_arguments"] == ["--skip_trimming"]
    assert preview["nfcore_runtime_params"] == {"skip_alignment": True, "skip_trimming": True}


def test_invalid_fastq_preprocessing_and_raw_count_preprocessing_fail(project_factory, tmp_path):
    root = _make_fastq_project(tmp_path)
    config = yaml.safe_load((root / "project.yaml").read_text(encoding="utf-8"))
    config["input"]["preprocessing"] = "trimmed_by_filename"
    (root / "project.yaml").write_text(yaml.safe_dump(config, sort_keys=False), encoding="utf-8")
    assert "invalid_project_config" in {issue.code for issue in validate_project(root).errors}

    raw_root = project_factory()
    raw_config = yaml.safe_load((raw_root / "project.yaml").read_text(encoding="utf-8"))
    raw_config["input"]["preprocessing"] = "raw"
    (raw_root / "project.yaml").write_text(yaml.safe_dump(raw_config, sort_keys=False), encoding="utf-8")
    assert "invalid_project_config" in {issue.code for issue in validate_project(raw_root).errors}


def test_invalid_fastq_extension_fails(tmp_path):
    root = _make_fastq_project(tmp_path)
    (root / "input" / "fastq" / "C3_R1.txt").write_text("not fastq", encoding="utf-8")
    assert "invalid_fastq_filename" in {issue.code for issue in validate_project(root).errors}


def test_appledouble_fastq_sidecars_are_ignored(tmp_path):
    root = _make_fastq_project(tmp_path)
    (root / "input" / "fastq" / "._C1_R1.fastq.gz").write_bytes(b"filesystem metadata")
    report = validate_project(root)
    assert report.is_valid
    assert report.fastq is not None
    assert report.fastq.sample_ids == ("C1", "C2", "T1", "T2")
