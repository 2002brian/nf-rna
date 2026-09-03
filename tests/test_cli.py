from __future__ import annotations

import hashlib
import json
from pathlib import Path

import yaml
from typer.testing import CliRunner

from rnaseq.cli import app

runner = CliRunner()


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def test_new_creates_versioned_project_without_placeholder_counts(tmp_path):
    result = runner.invoke(
        app,
        ["new"],
        input=f"created_project\n{tmp_path}\nMus musculus\nraw_counts\nL2\ntwo_group\n",
    )
    assert result.exit_code == 0, result.output
    root = tmp_path / "created_project"
    config = yaml.safe_load((root / "project.yaml").read_text(encoding="utf-8"))
    assert config["schema_version"] == "1.1"
    assert isinstance(config["schema_version"], str)
    assert (root / "input").is_dir()
    assert (root / "planning").is_dir()
    assert not (root / "input" / "counts.csv").exists()
    assert (root / "metadata.csv").read_text(encoding="utf-8") == "sample_id,condition\n"


def test_new_refuses_to_overwrite(tmp_path):
    (tmp_path / "existing").mkdir()
    result = runner.invoke(
        app,
        ["new"],
        input=f"existing\n{tmp_path}\nMus musculus\nraw_counts\nL1\ntwo_group\n",
    )
    assert result.exit_code == 1
    assert "already exists" in result.output


def test_new_uses_an_empty_same_named_destination_without_nesting(tmp_path):
    target = tmp_path / "client_case"
    target.mkdir()
    result = runner.invoke(
        app,
        ["new"],
        input=f"client_case\n{target}\nMus musculus\nraw_counts\nL1\ntwo_group\n",
    )
    assert result.exit_code == 0, result.output
    assert (target / "project.yaml").is_file()
    assert not (target / "client_case").exists()


def test_new_defaults_to_fastq_project(tmp_path):
    result = runner.invoke(
        app,
        ["new"],
        input=f"fastq_project\n{tmp_path}\nMus musculus\n\n\nL1\ntwo_group\n",
    )
    assert result.exit_code == 0, result.output
    config = yaml.safe_load((tmp_path / "fastq_project" / "project.yaml").read_text())
    assert config["input"] == {
        "type": "fastq",
        "path": "input/fastq",
        "layout": "paired_end",
        "preprocessing": "raw",
    }
    assert config["upstream"]["pipeline_version"] == "3.26.0"


def test_validate_example_passes():
    example = Path(__file__).parents[1] / "examples" / "simple_two_group"
    result = runner.invoke(app, ["validate", str(example)])
    assert result.exit_code == 0, result.output
    assert "Schema version: 1.0" in result.output
    assert "Genes: 100" in result.output
    assert "Validation result:\nPASS" in result.output


def test_sanitize_delivery_command_only_cleans_one_successful_run(tmp_path):
    run = tmp_path / "runs" / "CASE" / "20260902-133054+0800"
    delivery = run / "delivery"
    sidecar = delivery / "figures" / "png" / "l1" / "._pca.png"
    sidecar.parent.mkdir(parents=True)
    sidecar.write_bytes(b"sidecar")
    (delivery / "report.html").write_text("report", encoding="utf-8")
    (run / "run_state.json").write_text(json.dumps({"status": "SUCCESS"}), encoding="utf-8")

    result = runner.invoke(app, ["sanitize-delivery", str(run)])

    assert result.exit_code == 0, result.output
    assert "Delivery package sanitized" in result.output
    assert not sidecar.exists()
    assert (delivery / "report.html").read_text(encoding="utf-8") == "report"


def test_plan_is_deterministic_and_records_schema_version():
    example = Path(__file__).parents[1] / "examples" / "simple_two_group"
    plan_path = example / "planning" / "analysis_plan.md"
    manifest_path = example / "planning" / "manifest.preview.yaml"
    plan_path.unlink(missing_ok=True)
    manifest_path.unlink(missing_ok=True)

    first = runner.invoke(app, ["plan", str(example)])
    assert first.exit_code == 0, first.output
    first_digests = (digest(plan_path), digest(manifest_path))
    manifest = yaml.safe_load(manifest_path.read_text(encoding="utf-8"))
    assert manifest["schema"]["project_schema_version"] == "1.0"
    assert isinstance(manifest["schema"]["project_schema_version"], str)
    assert manifest["pipeline"]["version"] == "0.4.3"

    second = runner.invoke(app, ["plan", str(example)])
    assert second.exit_code == 0, second.output
    assert (digest(plan_path), digest(manifest_path)) == first_digests


def test_manifest_checksums_match_sources():
    example = Path(__file__).parents[1] / "examples" / "simple_two_group"
    result = runner.invoke(app, ["plan", str(example)])
    assert result.exit_code == 0
    manifest = yaml.safe_load(
        (example / "planning" / "manifest.preview.yaml").read_text(encoding="utf-8")
    )
    assert manifest["configuration"]["sha256"] == digest(example / "project.yaml")
    assert manifest["input"]["sha256"] == digest(example / "counts.csv")
    assert manifest["metadata"]["sha256"] == digest(example / "metadata.csv")
    assert manifest["contrasts"]["sha256"] == digest(example / "contrasts.csv")


def test_failed_plan_writes_no_artifacts(project_factory):
    root = project_factory(counts="gene_id,C1\nGeneA,-1\n")
    result = runner.invoke(app, ["plan", str(root)])
    assert result.exit_code == 1
    assert not (root / "planning" / "analysis_plan.md").exists()
    assert not (root / "planning" / "manifest.preview.yaml").exists()


def test_validate_system_error_returns_exit_code_two(monkeypatch, tmp_path):
    def fail(_project_dir):
        raise OSError("simulated filesystem failure")

    monkeypatch.setattr("rnaseq.cli.validate_project", fail)
    result = runner.invoke(app, ["validate", str(tmp_path)])
    assert result.exit_code == 2
    assert "SYSTEM ERROR" in result.output
