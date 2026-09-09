from __future__ import annotations

import hashlib
import json
from pathlib import Path

import yaml
from typer.testing import CliRunner

from rnaseq.cli import _wizard_completion_candidates, _wizard_tab_completion, app

runner = CliRunner()


def test_wizard_completion_candidates_match_only_canonical_prefixes():
    assert _wizard_completion_candidates("p", ["paired_end", "single_end"]) == ["paired_end"]
    assert _wizard_completion_candidates("s", ["paired_end", "single_end"]) == ["single_end"]
    assert _wizard_completion_candidates("l", ["igenomes", "local", "custom"]) == ["local"]
    assert _wizard_completion_candidates("h", ["salmon", "hisat2_featurecounts"]) == ["hisat2_featurecounts"]
    assert _wizard_completion_candidates("paried_end", ["paired_end", "single_end"]) == []


def test_wizard_completion_candidates_leave_ambiguous_prefix_unselected():
    assert _wizard_completion_candidates("", ["auto", "unstranded", "forward", "reverse"]) == [
        "auto",
        "unstranded",
        "forward",
        "reverse",
    ]
    assert _wizard_completion_candidates("a", ["auto", "analysis"]) == ["auto", "analysis"]


def test_wizard_tab_completion_is_temporary_and_keeps_ambiguous_candidates(monkeypatch):
    class FakeReadline:
        def __init__(self):
            self.completer = original_completer
            self.delimiters = " \t_"

        def get_completer(self):
            return self.completer

        def set_completer(self, completer):
            self.completer = completer

        def get_completer_delims(self):
            return self.delimiters

        def set_completer_delims(self, delimiters):
            self.delimiters = delimiters

    def original_completer(text, state):
        return None

    fake_readline = FakeReadline()
    monkeypatch.setattr("rnaseq.cli._is_interactive_terminal", lambda: True)
    monkeypatch.setattr("rnaseq.cli._readline_module", lambda: fake_readline)

    with _wizard_tab_completion(["auto", "analysis"]):
        assert fake_readline.delimiters == " \t"
        assert fake_readline.completer("a", 0) == "auto"
        assert fake_readline.completer("a", 1) == "analysis"
        assert fake_readline.completer("a", 2) is None

    assert fake_readline.completer is original_completer
    assert fake_readline.delimiters == " \t_"


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def test_new_creates_versioned_project_without_placeholder_counts(tmp_path):
    result = runner.invoke(
        app,
        ["new", "--name", "created_project", "--destination", str(tmp_path), "--species", "mouse", "--input-type", "raw_counts", "--preset", "L2", "--design-type", "two_group", "--scaffold", "--yes"],
    )
    assert result.exit_code == 0, result.output
    root = tmp_path / "created_project"
    config = yaml.safe_load((root / "project.yaml").read_text(encoding="utf-8"))
    assert config["schema_version"] == "1.2"
    assert isinstance(config["schema_version"], str)
    assert (root / "input").is_dir()
    assert (root / "planning").is_dir()
    assert not (root / "input" / "counts.csv").exists()
    assert (root / "metadata.csv").read_text(encoding="utf-8") == "sample_id,condition\n"


def test_new_records_explicit_total_execution_budget(monkeypatch, tmp_path):
    from rnaseq.execution import LocalResourceCapacity
    monkeypatch.setattr("rnaseq.cli.detect_local_resource_capacity", lambda: LocalResourceCapacity(20, 64, 62))
    result = runner.invoke(
        app,
        ["new", "--name", "resourced", "--destination", str(tmp_path), "--species", "mouse", "--input-type", "raw_counts", "--preset", "L2", "--design-type", "two_group", "--execution-profile", "local", "--cpus", "16", "--memory-gb", "48", "--scaffold", "--yes"],
    )
    assert result.exit_code == 0, result.output
    config = yaml.safe_load((tmp_path / "resourced" / "project.yaml").read_text(encoding="utf-8"))
    assert config["execution"] == {"profile": "local", "max_cpus": 16, "max_memory_gb": 48}


def test_new_refuses_to_overwrite(tmp_path):
    (tmp_path / "existing").mkdir()
    result = runner.invoke(
        app,
        ["new", "--name", "existing", "--destination", str(tmp_path), "--species", "mouse", "--input-type", "raw_counts", "--preset", "L1", "--design-type", "two_group", "--scaffold", "--yes"],
    )
    assert result.exit_code == 1
    assert "already exists" in result.output


def test_new_uses_an_empty_same_named_destination_without_nesting(tmp_path):
    target = tmp_path / "client_case"
    target.mkdir()
    result = runner.invoke(
        app,
        ["new", "--name", "client_case", "--destination", str(target), "--species", "mouse", "--input-type", "raw_counts", "--preset", "L1", "--design-type", "two_group", "--scaffold", "--yes"],
    )
    assert result.exit_code == 0, result.output
    assert (target / "project.yaml").is_file()
    assert not (target / "client_case").exists()


def test_new_defaults_to_fastq_project(tmp_path):
    result = runner.invoke(
        app,
        ["new", "--name", "fastq_project", "--destination", str(tmp_path), "--species", "mouse", "--input-type", "fastq", "--preset", "L1", "--design-type", "two_group", "--scaffold", "--yes"],
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


def test_interactive_fastq_new_retries_only_invalid_choice_and_creates_scaffold(monkeypatch, tmp_path):
    from rnaseq.execution import LocalResourceCapacity

    monkeypatch.setattr("rnaseq.cli._is_interactive_terminal", lambda: True)
    monkeypatch.setattr("rnaseq.cli.detect_local_resource_capacity", lambda: LocalResourceCapacity(20, 64, 62))
    result = runner.invoke(app, ["new"], input=(
        "fastq_scaffold\n"
        f"{tmp_path}\n"
        "mouse\n"
        "fastq\n"
        "paried_end\n"
        "\n"  # Retry the layout question and accept paired_end.
        "\n"  # raw preprocessing
        "\n"  # Salmon
        "\n"  # auto strandedness
        "\n"  # iGenomes
        "qc\n"
        "\n"  # suggested CPUs
        "\n"  # suggested memory
        "y\n"
    ))

    assert result.exit_code == 0, result.output
    assert "Invalid choice 'paried_end'." in result.output
    assert "Please choose one of: paired_end, single_end." in result.output
    assert "Traceback" not in result.output
    assert "Import existing inputs now?" not in result.output
    assert result.output.count("Project name") == 1
    assert "Input handling: scaffold — add data after project creation" in result.output
    root = tmp_path / "fastq_scaffold"
    config = yaml.safe_load((root / "project.yaml").read_text(encoding="utf-8"))
    assert config["organism"]["species"] == "Mus musculus"
    assert config["input"]["layout"] == "paired_end"
    assert config["execution"] == {"profile": "local", "max_cpus": 16, "max_memory_gb": 48}
    assert (root / "input" / "fastq").is_dir()
    assert not any((root / "input" / "fastq").iterdir())


def test_interactive_raw_count_new_is_scaffold_first_and_retries_bad_integer(monkeypatch, tmp_path):
    from rnaseq.execution import LocalResourceCapacity

    monkeypatch.setattr("rnaseq.cli._is_interactive_terminal", lambda: True)
    monkeypatch.setattr("rnaseq.cli.detect_local_resource_capacity", lambda: LocalResourceCapacity(20, 64, 62))
    result = runner.invoke(app, ["new"], input=(
        "count_scaffold\n"
        f"{tmp_path}\n"
        "human\n"
        "raw_counts\n"
        "L2\n"
        "two_group\n"
        "sixteen\n"
        "\n"
        "\n"
        "y\n"
    ))

    assert result.exit_code == 0, result.output
    assert "Invalid value 'sixteen'. Please enter a positive integer." in result.output
    assert "Import existing inputs now?" not in result.output
    root = tmp_path / "count_scaffold"
    assert (root / "input").is_dir()
    assert not (root / "input" / "counts.csv").exists()
    assert (root / "metadata.csv").is_file()
    assert (root / "contrasts.csv").is_file()
    validation = runner.invoke(app, ["validate", str(root)])
    assert validation.exit_code == 1
    assert "Count matrix not found" in validation.output


def test_new_imports_raw_counts_with_metadata_and_contrasts(tmp_path):
    counts = tmp_path / "source_counts.csv"
    metadata = tmp_path / "source_metadata.csv"
    contrasts = tmp_path / "source_contrasts.csv"
    counts.write_text("gene_id,C1,T1\nGeneA,2,5\n", encoding="utf-8")
    metadata.write_text("sample_id,condition\nC1,Control\nT1,Treatment\n", encoding="utf-8")
    contrasts.write_text("contrast_id,factor,numerator,denominator\nT_vs_C,condition,Treatment,Control\n", encoding="utf-8")

    result = runner.invoke(app, [
        "new", "--name", "imported", "--destination", str(tmp_path), "--species", "human",
        "--input-type", "raw_counts", "--counts", str(counts), "--metadata", str(metadata),
        "--contrasts", str(contrasts), "--preset", "L1", "--design-type", "two_group", "--yes",
    ])

    assert result.exit_code == 0, result.output
    root = tmp_path / "imported"
    assert (root / "input" / "counts.csv").read_bytes() == counts.read_bytes()
    assert (root / "metadata.csv").read_bytes() == metadata.read_bytes()
    assert runner.invoke(app, ["validate", str(root)]).exit_code == 0


def test_new_imports_lane_fastqs_and_preserves_explicit_hisat2_backend(tmp_path):
    source = tmp_path / "source"
    source.mkdir()
    for name in ("S1_L001_R1.fastq.gz", "S1_L001_R2.fastq.gz", "S1_L002_R1.fastq.gz", "S1_L002_R2.fastq.gz", "S2_L001_R1.fastq.gz", "S2_L001_R2.fastq.gz"):
        (source / name).write_bytes(b"not-inspected-by-creation")
    sheet = source / "samplesheet.csv"
    sheet.write_text(
        "sample,fastq_1,fastq_2,strandedness\n"
        "S1,S1_L001_R1.fastq.gz,S1_L001_R2.fastq.gz,forward\n"
        "S1,S1_L002_R1.fastq.gz,S1_L002_R2.fastq.gz,forward\n"
        "S2,S2_L001_R1.fastq.gz,S2_L001_R2.fastq.gz,forward\n",
        encoding="utf-8",
    )
    metadata = source / "metadata.csv"
    metadata.write_text("sample_id,condition\nS1,Control\nS2,Treatment\n", encoding="utf-8")
    contrasts = source / "contrasts.csv"
    contrasts.write_text("contrast_id,factor,numerator,denominator\nT_vs_C,condition,Treatment,Control\n", encoding="utf-8")
    fasta, gtf = source / "genome.fa", source / "genes.gtf"
    fasta.write_text(">chr1\nA\n", encoding="utf-8")
    gtf.write_text("# annotation\n", encoding="utf-8")

    result = runner.invoke(app, [
        "new", "--name", "lanes", "--destination", str(tmp_path), "--species", "human",
        "--input-type", "fastq", "--fastq-samplesheet", str(sheet), "--metadata", str(metadata),
        "--contrasts", str(contrasts), "--method", "hisat2_featurecounts", "--strandedness", "forward",
        "--reference-source", "custom", "--reference-fasta", str(fasta), "--reference-gtf", str(gtf),
        "--preset", "L1", "--design-type", "two_group", "--yes",
    ])

    assert result.exit_code == 0, result.output
    root = tmp_path / "lanes"
    config = yaml.safe_load((root / "project.yaml").read_text(encoding="utf-8"))
    assert config["upstream"]["quantification"]["method"] == "hisat2_featurecounts"
    assert config["upstream"]["strandedness"] == "forward"
    assert (root / "input" / "fastq" / "S1_L002_R2.fastq.gz").is_file()
    assert (root / "planning" / "imported_fastq_samplesheet.csv").read_bytes() == sheet.read_bytes()
    validation = runner.invoke(app, ["validate", str(root)])
    assert validation.exit_code == 0, validation.output
    assert "NOT READY" in validation.output
    assert "HISAT2 requires" in validation.output


def test_new_rejects_mixed_fastq_strandedness_without_creating_project(tmp_path):
    source = tmp_path / "source"
    source.mkdir()
    for name in ("S1_R1.fastq.gz", "S1_R2.fastq.gz", "S2_R1.fastq.gz", "S2_R2.fastq.gz"):
        (source / name).write_bytes(b"x")
    sheet = source / "samplesheet.csv"
    sheet.write_text(
        "sample,fastq_1,fastq_2,strandedness\nS1,S1_R1.fastq.gz,S1_R2.fastq.gz,forward\nS2,S2_R1.fastq.gz,S2_R2.fastq.gz,reverse\n",
        encoding="utf-8",
    )
    result = runner.invoke(app, [
        "new", "--name", "bad", "--destination", str(tmp_path), "--species", "mouse", "--input-type", "fastq",
        "--fastq-samplesheet", str(sheet), "--method", "salmon", "--preset", "qc", "--design-type", "two_group", "--yes",
    ])
    assert result.exit_code == 1
    assert "mixed strandedness" in result.output
    assert not (tmp_path / "bad").exists()


def test_fastq_qc_scope_validates_and_plans_without_statistical_metadata(tmp_path):
    source = tmp_path / "source"
    source.mkdir()
    for name in ("S1_R1.fastq.gz", "S1_R2.fastq.gz"):
        (source / name).write_bytes(b"x")
    sheet = source / "samplesheet.csv"
    sheet.write_text("sample,fastq_1,fastq_2,strandedness\nS1,S1_R1.fastq.gz,S1_R2.fastq.gz,auto\n", encoding="utf-8")
    result = runner.invoke(app, [
        "new", "--name", "technical_qc", "--destination", str(tmp_path), "--species", "mouse", "--input-type", "fastq",
        "--fastq-samplesheet", str(sheet), "--method", "salmon", "--reference-source", "igenomes",
        "--preset", "qc", "--yes",
    ])
    assert result.exit_code == 0, result.output
    root = tmp_path / "technical_qc"
    assert runner.invoke(app, ["validate", str(root)]).exit_code == 0
    planned = runner.invoke(app, ["plan", str(root)])
    assert planned.exit_code == 0, planned.output
    assert "technical FASTQ QC only" in (root / "planning" / "analysis_plan.md").read_text(encoding="utf-8")


def test_new_cancellation_writes_no_partial_project(tmp_path):
    result = runner.invoke(app, [
        "new", "--name", "cancelled", "--destination", str(tmp_path), "--species", "mouse",
        "--input-type", "raw_counts", "--preset", "L1", "--design-type", "two_group", "--scaffold",
    ], input="n\n")
    assert result.exit_code == 0, result.output
    assert "cancelled" in result.output
    assert not (tmp_path / "cancelled").exists()


def test_new_requires_flags_when_stdin_is_not_a_terminal():
    result = runner.invoke(app, ["new"])
    assert result.exit_code == 1
    assert "interactive terminal" in result.output


def test_invalid_import_is_rejected_before_project_becomes_visible(tmp_path):
    counts, metadata, contrasts = (tmp_path / "counts.csv", tmp_path / "metadata.csv", tmp_path / "contrasts.csv")
    counts.write_text("gene_id,C1,T1\nGeneA,2,5\n", encoding="utf-8")
    metadata.write_text("sample_id,condition\nC1,Control\nWRONG,Treatment\n", encoding="utf-8")
    contrasts.write_text("contrast_id,factor,numerator,denominator\nT_vs_C,condition,Treatment,Control\n", encoding="utf-8")
    result = runner.invoke(app, [
        "new", "--name", "invalid_import", "--destination", str(tmp_path), "--species", "human",
        "--input-type", "raw_counts", "--counts", str(counts), "--metadata", str(metadata), "--contrasts", str(contrasts),
        "--preset", "L1", "--design-type", "two_group", "--yes",
    ])
    assert result.exit_code == 1
    assert "strict validation" in result.output
    assert not (tmp_path / "invalid_import").exists()


def test_new_uses_imported_metadata_fields_for_noninteractive_formula(tmp_path):
    counts, metadata, contrasts = (tmp_path / "counts.csv", tmp_path / "metadata.csv", tmp_path / "contrasts.csv")
    counts.write_text("gene_id,S1,S2,S3,S4\nGeneA,1,2,3,4\n", encoding="utf-8")
    metadata.write_text("sample_id,subject,condition,batch\nS1,A,Control,B1\nS2,A,Treatment,B2\nS3,B,Control,B2\nS4,B,Treatment,B1\n", encoding="utf-8")
    contrasts.write_text("contrast_id,factor,numerator,denominator\nT_vs_C,condition,Treatment,Control\n", encoding="utf-8")
    result = runner.invoke(app, [
        "new", "--name", "paired", "--destination", str(tmp_path), "--species", "mouse", "--input-type", "raw_counts",
        "--counts", str(counts), "--metadata", str(metadata), "--contrasts", str(contrasts), "--preset", "L2",
        "--design-type", "paired", "--condition-column", "condition", "--pairing-column", "subject", "--covariate", "batch", "--yes",
    ])
    assert result.exit_code == 0, result.output
    config = yaml.safe_load((tmp_path / "paired" / "project.yaml").read_text(encoding="utf-8"))
    assert config["design"]["formula"] == "~ subject + batch + condition"
    assert config["design"]["pairing_column"] == "subject"


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
    assert manifest["pipeline"]["version"] == "0.5.1"

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
