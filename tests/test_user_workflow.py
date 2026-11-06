"""Beginner-facing workflow contract: wizard next steps, validation blockers,
early ``--reuse-upstream`` preflight, and the commands shown in the user docs."""

from __future__ import annotations

import gzip
import json
import re
import shlex
from pathlib import Path

import pytest
import typer
import yaml
from typer.testing import CliRunner

from rnaseq.cli import app
from rnaseq.errors import ExecutionPreflightError
from rnaseq.execution import LocalResourceCapacity
from rnaseq.service import execute_service_run, preflight_reuse_upstream
from rnaseq.validators import FASTQ_NAME_PATTERN, VALIDATION_BLOCKER, validate_project

runner = CliRunner()
SAMPLES = ("Ctrl1_S1", "Ctrl2_S2", "Ctrl3_S3", "Trt1_S4", "Trt2_S5", "Trt3_S6")
METADATA = "sample_id,condition\n" + "".join(
    f"{sample},{'Control' if sample.startswith('Ctrl') else 'Treatment'}\n" for sample in SAMPLES
)
RUN_ID = "20261029-101010+0800"


def _new(tmp_path: Path, name: str, *, preset: str, input_type: str = "fastq", destination: Path | None = None):
    return runner.invoke(app, [
        "new", "--name", name, "--destination", str(destination or tmp_path), "--species", "mouse",
        "--input-type", input_type, "--scaffold", "--preset", preset, "--design-type", "two_group",
        "--gene-id-type", "ENSEMBL", "--yes",
    ])


def _fastq_project(tmp_path: Path, name: str = "proj", *, preset: str = "L1", metadata: str = METADATA) -> Path:
    """An iGenomes GRCm39 Salmon project with Illumina-named placeholder reads (never executed)."""

    result = _new(tmp_path, name, preset=preset)
    assert result.exit_code == 0, result.output
    root = tmp_path / name
    for sample in SAMPLES:
        for read in ("R1", "R2"):
            path = root / "input" / "fastq" / f"{sample}_L001_{read}_001.fastq.gz"
            path.write_bytes(gzip.compress(b"@r1\nACGTACGTAC\n+\nIIIIIIIIII\n"))
    (root / "metadata.csv").write_text(metadata, encoding="utf-8")
    return root


def _record_run(root: Path, case: str = "L1-001", run_id: str = RUN_ID, *, status: str = "SUCCESS", handoff: bool = True) -> Path:
    run_dir = root / "runs" / case / run_id
    (run_dir / "handoff").mkdir(parents=True)
    (run_dir / "run_state.json").write_text(json.dumps({"status": status}), encoding="utf-8")
    if handoff:
        (run_dir / "handoff" / "upstream_manifest.yaml").write_text("pipeline: {}\n", encoding="utf-8")
    return run_dir


# ---------------------------------------------------------------- D1: wizard next steps


def _wizard(monkeypatch, tmp_path: Path, name: str, preset: str):
    monkeypatch.setattr("rnaseq.cli._is_interactive_terminal", lambda: True)
    monkeypatch.setattr("rnaseq.cli.detect_local_resource_capacity", lambda: LocalResourceCapacity(20, 64, 62))
    monkeypatch.setattr("rnaseq.cli.compatible_registered_references", lambda species, backend: [])
    answers = [name, str(tmp_path), "mouse", "fastq", "", "", "", "", "", preset, "", ""]
    answers += [""] if preset == "L2" else []  # accept the default L2 enrichment
    answers += ["", "", "", "y"]  # ENSEMBL, CPU auto, memory auto, create
    return runner.invoke(app, ["new"], input="\n".join(answers) + "\n")


def test_l2_wizard_sends_the_user_to_design_and_never_to_an_independent_l1_run(monkeypatch, tmp_path):
    result = _wizard(monkeypatch, tmp_path, "direct_l2", "L2")
    assert result.exit_code == 0, result.output
    root = tmp_path / "direct_l2"
    assert yaml.safe_load((root / "project.yaml").read_text(encoding="utf-8"))["project"]["preset"] == "L2"
    assert "Scope: L2" in result.output
    assert "Contrasts: set later by `rnaseq design` (required before L2 can run)" in result.output
    output = result.output[result.output.index("Next steps"):]
    assert "run L1" not in output and "Run L1" not in output and "PROJECT" not in output
    assert "the L2 run also performs the L1 QC" in output
    assert "no separate L1 run is needed" in output
    design, validate, run = (output.index(f"rnaseq {command} {root}") for command in ("design", "validate", "run"))
    assert design < validate < run
    assert "refreshes the planning artifacts" in output


def test_l1_wizard_lists_validate_plan_run_review_then_design(monkeypatch, tmp_path):
    result = _wizard(monkeypatch, tmp_path, "staged_l1", "L1")
    assert result.exit_code == 0, result.output
    root = tmp_path / "staged_l1"
    output = result.output[result.output.index("Next steps"):]
    assert "PROJECT" not in output
    steps = [f"rnaseq validate {root}", f"rnaseq plan {root}", f"Run L1: rnaseq run {root} --case-id CASE-ID",
             "Review the L1 QC", f"rnaseq design {root}"]
    positions = [output.index(step) for step in steps]
    assert positions == sorted(positions)
    assert "--reuse-upstream CASE-ID/RUN-ID" in output  # FASTQ: L2 can reuse the L1 quantification
    assert "Ctrl1_S1_L001_R1_001.fastq.gz -> Ctrl1_S1" in output


def test_raw_count_scaffold_steps_never_mention_upstream_reuse_and_quote_the_path(tmp_path):
    destination = tmp_path / "my projects"
    destination.mkdir()
    result = _new(tmp_path, "counts", preset="L1", input_type="raw_counts", destination=destination)
    assert result.exit_code == 0, result.output
    quoted = shlex.quote(str(destination / "counts"))
    assert quoted.startswith("'")
    assert f"rnaseq validate {quoted}" in result.output and f"rnaseq design {quoted}" in result.output
    assert "--reuse-upstream" not in result.output
    assert "input/counts.csv" in result.output


# ---------------------------------------------------------------- D2: validation blockers


def test_input_errors_are_not_reported_as_a_missing_igenomes_reference(tmp_path):
    trimmed = METADATA.replace("_S1", "").replace("_S2", "").replace("_S3", "").replace("_S4", "").replace("_S5", "").replace("_S6", "")
    root = _fastq_project(tmp_path, metadata=trimmed)
    report = validate_project(root)
    assert not report.is_valid
    assert report.execution_blockers == (VALIDATION_BLOCKER,)
    result = runner.invoke(app, ["validate", str(root)])
    assert result.exit_code == 1
    assert "Samples present in FASTQ input but absent from metadata: Ctrl1_S1" in result.output
    assert "reference genome not configured" not in result.output
    assert VALIDATION_BLOCKER in result.output


def test_l2_missing_contrasts_names_design_without_blaming_the_reference(tmp_path):
    root = _fastq_project(tmp_path, preset="L2")
    result = runner.invoke(app, ["validate", str(root)])
    assert result.exit_code == 1
    assert "Contrasts: not yet defined (required for L2: run `rnaseq design PROJECT` to configure them)" in result.output
    assert "after reviewing L1" not in result.output
    assert "reference genome not configured" not in result.output


def test_a_genuinely_missing_reference_is_still_reported(tmp_path):
    root = _fastq_project(tmp_path)
    config = yaml.safe_load((root / "project.yaml").read_text(encoding="utf-8"))
    config["reference"]["genome"] = None
    (root / "project.yaml").write_text(yaml.safe_dump(config, sort_keys=False), encoding="utf-8")
    assert validate_project(root).execution_blockers == ("reference genome not configured",)
    (root / "metadata.csv").write_text("sample_id,condition\n", encoding="utf-8")
    assert validate_project(root).execution_blockers == (VALIDATION_BLOCKER, "reference genome not configured")


# ---------------------------------------------------------------- D3/D4: reuse preflight


def _designed_l2(tmp_path: Path) -> Path:
    root = _fastq_project(tmp_path)
    result = runner.invoke(app, ["design", str(root), "--yes"])
    assert result.exit_code == 0, result.output
    return root


@pytest.mark.parametrize(("value", "message"), (
    (RUN_ID, "must be formatted as CASE-ID/RUN-ID, for example L1-001/20261029-101010+0800"),
    ("L1-001", "must be formatted as CASE-ID/RUN-ID"),
    ("L1-001/latest", "unsafe run ID"),
    ("L1-001/20261029-101010+0800", "no recorded run L1-001/20261029-101010+0800 exists"),
    ("FAILED-1/20261029-101010+0800", "FAILED-1/20261029-101010+0800 is FAILED"),
    ("NOHANDOFF/20261029-101010+0800", "no validated upstream handoff manifest"),
))
def test_unusable_reuse_source_is_refused_before_any_run_record_exists(monkeypatch, tmp_path, value, message):
    root = _designed_l2(tmp_path)
    _record_run(root, "FAILED-1", status="FAILED")
    _record_run(root, "NOHANDOFF", handoff=False)
    before = sorted(path.relative_to(root) for path in (root / "runs").rglob("*"))
    monkeypatch.setattr("rnaseq.cli.prepare_service_run", lambda *_a, **_k: pytest.fail("preflight ran after the reuse check"))
    monkeypatch.setattr("rnaseq.cli.execute_service_run", lambda *_a, **_k: pytest.fail("a run was started"))
    result = runner.invoke(app, ["run", str(root), "--case-id", "L2-001", "--reuse-upstream", value, "--yes"])
    assert result.exit_code == 1
    assert message in result.output
    assert sorted(path.relative_to(root) for path in (root / "runs").rglob("*")) == before
    assert not (root / "runs" / "L2-001").exists()


def test_service_checks_the_reuse_source_before_preflight_or_run_creation(monkeypatch, tmp_path):
    root = _designed_l2(tmp_path)
    monkeypatch.setattr("rnaseq.service.prepare_service_run", lambda *_a, **_k: pytest.fail("prepared before the reuse check"))
    monkeypatch.setattr("rnaseq.service.create_case_run", lambda *_a, **_k: pytest.fail("run created before the reuse check"))
    with pytest.raises(ExecutionPreflightError, match="formatted as CASE-ID/RUN-ID"):
        execute_service_run(validate_project(root), case_id="L2-001", reuse_upstream=RUN_ID)
    assert not (root / "runs").exists()


def test_valid_fastq_reuse_source_passes_the_early_preflight(tmp_path):
    root = _designed_l2(tmp_path)
    _record_run(root)
    preflight_reuse_upstream(validate_project(root), f"L1-001/{RUN_ID}")


def test_raw_count_project_rejects_reuse_upstream_before_running(monkeypatch, tmp_path):
    assert _new(tmp_path, "counts", preset="L1", input_type="raw_counts").exit_code == 0
    root = tmp_path / "counts"
    (root / "metadata.csv").write_text(METADATA, encoding="utf-8")
    rows = "".join(f"G{gene}," + ",".join(str(10 + gene + index) for index in range(len(SAMPLES))) + "\n" for gene in range(5))
    (root / "input" / "counts.csv").write_text("gene_id," + ",".join(SAMPLES) + "\n" + rows, encoding="utf-8")
    _record_run(root)
    report = validate_project(root)
    assert report.is_valid, [issue.message for issue in report.errors]
    with pytest.raises(ExecutionPreflightError, match="applies only to FASTQ projects"):
        preflight_reuse_upstream(report, f"L1-001/{RUN_ID}")
    monkeypatch.setattr("rnaseq.cli.prepare_service_run", lambda *_a, **_k: pytest.fail("raw-count reuse reached preflight"))
    result = runner.invoke(app, ["run", str(root), "--case-id", "L1-002", "--reuse-upstream", f"L1-001/{RUN_ID}", "--yes"])
    assert result.exit_code == 1
    assert "applies only to FASTQ projects" in result.output and "Run without --reuse-upstream" in result.output


# ---------------------------------------------------------------- user documentation


DOCS = ("README.md", "README_zh-TW.md", "docs/quickstart.md", "docs/runtime.md")
ROOT = Path(__file__).parents[1]


def _documented_commands():
    for name in DOCS:
        text = (ROOT / name).read_text(encoding="utf-8")
        for block in re.findall(r"```bash\n(.*?)```", text, flags=re.S):
            for line in block.splitlines():
                line = line.split(" #", 1)[0].strip()
                for segment in re.split(r"\s*(?:&&|;)\s*", line):
                    words = segment.split()
                    while words and words[0] in {"nohup", "setsid"}:
                        words = words[1:]
                    if words and words[0] == "rnaseq":
                        yield name, words


def test_every_documented_rnaseq_command_and_option_exists():
    root_command = typer.main.get_command(app)
    seen = 0
    for name, words in _documented_commands():
        command, position = root_command, 1
        # Typer may vendor its own click, so groups are recognised by their subcommands.
        while hasattr(command, "commands") and position < len(words) and words[position] in command.commands:
            command, position = command.commands[words[position]], position + 1
        assert not hasattr(command, "commands") or "--version" in words, f"{name}: {' '.join(words)}"
        options = {option for parameter in command.params for option in getattr(parameter, "opts", [])}
        for word in words[position:]:
            if word.startswith("--"):
                assert word.split("=", 1)[0] in options, f"{name}: unknown option {word} in: {' '.join(words)}"
        seen += 1
    assert seen > 20


def test_documented_fastq_sample_id_example_matches_the_parser():
    match = FASTQ_NAME_PATTERN.fullmatch("Ctrl1_S1_L001_R1_001.fastq.gz")
    assert match is not None and match.group("sample") == "Ctrl1_S1"
    for name in ("README.md", "README_zh-TW.md", "docs/quickstart.md"):
        text = (ROOT / name).read_text(encoding="utf-8")
        assert "Ctrl1_S1_L001_R1_001.fastq.gz" in text and "`Ctrl1_S1`" in text, name


def test_docs_no_longer_claim_the_l1_to_l2_journey_needs_the_development_version():
    for name in ("README.md", "README_zh-TW.md"):
        text = (ROOT / name).read_text(encoding="utf-8")
        assert "requires the development version" not in text and "需要 `main` 上的開發版本" not in text, name
    assert "v1.3.0" not in (ROOT / "docs" / "runtime.md").read_text(encoding="utf-8").split("## Historical")[0]
