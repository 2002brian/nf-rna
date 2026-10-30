"""L2 enables every supported enrichment method by default, with an explicit opt-out."""

from __future__ import annotations

import csv
import hashlib
import io
import json
from datetime import datetime
from pathlib import Path

import pytest
import yaml
from typer.testing import CliRunner

from conftest import base_config
from rnaseq.cli import app
from rnaseq.errors import ExecutionPreflightError
from rnaseq.models import ProjectConfig, production_enrichment_backends
from rnaseq.planner import generate_plan
from rnaseq.service import (
    _delivery_readme_analysis_summary, _enrichment_preflight, build_downstream_nextflow_command,
    create_case_run, freeze_case_inputs, kegg_rest_preflight, resolve_downstream_inputs,
)
from rnaseq.validators import validate_project
from rnaseq.workflow_support import enrichment_config, report
from test_workflow_support import _frozen_contract, _write_report_artifacts

runner = CliRunner()
ALL_METHODS = ["go", "kegg", "gsea"]
PARTIALLY_CONFOUNDED = (
    "sample_id,condition,batch\n"
    "c1,control,b1\nc2,control,b2\nc3,control,b2\n"
    "h1,high,b1\nh2,high,b2\nh3,high,b3\nh4,high,b3\n"
    "l1,low,b3\nl2,low,b3\nl3,low,b3\n"
)
THREE_CONTRASTS = (
    "contrast_id,factor,numerator,denominator\n"
    "high_vs_control,condition,high,control\nlow_vs_control,condition,low,control\nlow_vs_high,condition,low,high\n"
)


def _new(tmp_path: Path, name: str, *extra: str, preset: str = "L2", design: str = "multi_group", species: str = "mouse"):
    return runner.invoke(app, [
        "new", "--name", name, "--destination", str(tmp_path), "--species", species, "--input-type", "raw_counts",
        "--scaffold", "--preset", preset, "--design-type", design, "--yes", *extra,
    ])


def _config(root: Path) -> dict:
    return yaml.safe_load((root / "project.yaml").read_text(encoding="utf-8"))


def _tree_digest(root: Path) -> dict[str, str]:
    return {str(path.relative_to(root)): hashlib.sha256(path.read_bytes()).hexdigest() for path in sorted(root.rglob("*")) if path.is_file()}


# ------------------------------------------------------------- project creation


@pytest.mark.parametrize("design_flags", ([], ["--batch"]), ids=("multi_group", "batch_adjusted_multi_group"))
def test_noninteractive_l2_selects_every_method_including_batch_adjusted_multi_group(tmp_path, design_flags):
    result = _new(tmp_path, "proj", "--gene-id-type", "ENSEMBL", *design_flags)
    assert result.exit_code == 0, result.output
    config = _config(tmp_path / "proj")
    assert config["analysis"] == {"enrichment": ALL_METHODS}
    assert config["annotation"] == {"organism": "Mus musculus", "input_id_type": "ENSEMBL"}
    assert config["design"]["formula"] == ("~ batch + condition" if design_flags else "~ condition")
    assert "Enrichment: GO ORA, KEGG ORA, GO GSEA, KEGG GSEA" in result.output
    assert "Annotation: Mus musculus; ENSEMBL -> ENTREZID via org.Mm.eg.db" in result.output
    parsed = ProjectConfig.model_validate(config)
    assert production_enrichment_backends(parsed.analysis.enrichment) == ("go", "kegg", "gsea-go", "gsea-kegg")


def test_noninteractive_l2_without_gene_id_type_is_an_actionable_error_and_writes_nothing(tmp_path):
    result = _new(tmp_path, "proj")
    assert result.exit_code == 1
    assert "--gene-id-type ENSEMBL, ENTREZID, or SYMBOL" in result.output
    assert "--no-enrichment" in result.output
    assert not (tmp_path / "proj").exists()


def test_l2_explicit_opt_out_records_an_empty_selection(tmp_path):
    result = _new(tmp_path, "proj", "--no-enrichment")
    assert result.exit_code == 0, result.output
    config = _config(tmp_path / "proj")
    assert config["analysis"] == {"enrichment": []}
    assert "annotation" not in config
    assert "Enrichment: none (analysis.enrichment: [])" in result.output


def test_gene_id_type_is_validated(tmp_path):
    result = _new(tmp_path, "proj", "--gene-id-type", "refseq")
    assert result.exit_code == 1
    assert "ENSEMBL, ENTREZID, SYMBOL" in result.output
    assert not (tmp_path / "proj").exists()


def test_l1_records_annotation_but_selects_no_enrichment_and_rejects_enrichment_flags(tmp_path):
    result = _new(tmp_path, "l1", "--gene-id-type", "SYMBOL", preset="L1")
    assert result.exit_code == 0, result.output
    config = _config(tmp_path / "l1")
    assert config["analysis"] == {"enrichment": []}
    assert config["annotation"] == {"organism": "Mus musculus", "input_id_type": "SYMBOL"}
    assert "selected when L2 is configured by `rnaseq design`" in result.output

    bare = _new(tmp_path, "bare", preset="L1")
    assert bare.exit_code == 0, bare.output
    assert "annotation" not in _config(tmp_path / "bare")

    for flag in ("--no-enrichment", "--enrichment"):
        rejected = _new(tmp_path, f"flag{flag}", flag, preset="L1")
        assert rejected.exit_code == 1
        assert "does not apply to L1" in rejected.output


def test_interactive_l2_shows_methods_and_asks_for_the_gene_id_type(monkeypatch, tmp_path):
    from rnaseq.execution import LocalResourceCapacity

    monkeypatch.setattr("rnaseq.cli._is_interactive_terminal", lambda: True)
    monkeypatch.setattr("rnaseq.cli.detect_local_resource_capacity", lambda: LocalResourceCapacity(20, 64, 62))
    result = runner.invoke(app, ["new"], input=(
        "wiz\n" f"{tmp_path}\n" "mouse\n" "raw_counts\n" "L2\n" "multi_group\n"
        "y\n"  # known technical batch
        "\n"  # run the default enrichment
        "SYMBOL\n"
        "\n" "\n" "y\n"
    ))
    assert result.exit_code == 0, result.output
    assert "Functional enrichment (default for L2): GO ORA, KEGG ORA, GO GSEA, KEGG GSEA" in result.output
    assert "Annotation database for Mus musculus: org.Mm.eg.db" in result.output
    config = _config(tmp_path / "wiz")
    assert config["design"]["formula"] == "~ batch + condition"
    assert config["analysis"] == {"enrichment": ALL_METHODS}
    assert config["annotation"] == {"organism": "Mus musculus", "input_id_type": "SYMBOL"}


def test_interactive_l2_can_decline_enrichment_without_an_annotation_question(monkeypatch, tmp_path):
    from rnaseq.execution import LocalResourceCapacity

    monkeypatch.setattr("rnaseq.cli._is_interactive_terminal", lambda: True)
    monkeypatch.setattr("rnaseq.cli.detect_local_resource_capacity", lambda: LocalResourceCapacity(20, 64, 62))
    result = runner.invoke(app, ["new"], input=(
        "wiz\n" f"{tmp_path}\n" "human\n" "raw_counts\n" "L2\n" "two_group\n" "\n" "n\n" "\n" "\n" "y\n"
    ))
    assert result.exit_code == 0, result.output
    assert "Gene identifier type" not in result.output
    config = _config(tmp_path / "wiz")
    assert config["analysis"] == {"enrichment": []}
    assert "annotation" not in config


def test_qc_preset_rejects_enrichment_options(tmp_path):
    result = runner.invoke(app, [
        "new", "--name", "qc", "--destination", str(tmp_path), "--species", "mouse", "--input-type", "fastq",
        "--scaffold", "--preset", "qc", "--gene-id-type", "ENSEMBL", "--yes",
    ])
    assert result.exit_code == 1
    assert "do not apply to the qc preset" in result.output


# ------------------------------------------------------------- rnaseq design (L1 -> L2)


def _filled_l1(tmp_path: Path, *extra: str) -> Path:
    assert _new(tmp_path, "proj", "--batch", *extra, preset="L1").exit_code == 0
    root = tmp_path / "proj"
    (root / "metadata.csv").write_text(PARTIALLY_CONFOUNDED, encoding="utf-8")
    samples = [line.split(",")[0] for line in PARTIALLY_CONFOUNDED.strip().splitlines()[1:]]
    rows = "".join(f"G{g}," + ",".join(str(10 + 3 * g + 7 * index % 11) for index, _ in enumerate(samples)) + "\n" for g in range(3))
    (root / "input" / "counts.csv").write_text("gene_id," + ",".join(samples) + "\n" + rows, encoding="utf-8")
    return root


def test_design_promotion_selects_every_method_for_a_batch_adjusted_multi_group_project(tmp_path):
    root = _filled_l1(tmp_path, "--gene-id-type", "ENSEMBL")
    result = runner.invoke(app, ["design", str(root), "--yes"])
    assert result.exit_code == 0, result.output
    config = _config(root)
    assert config["project"]["preset"] == "L2"
    assert config["design"]["formula"] == "~ batch + condition"
    assert config["analysis"] == {"enrichment": ALL_METHODS}
    assert "Functional enrichment: GO ORA, KEGG ORA, GO GSEA, KEGG GSEA" in result.output
    assert "occurs in only one batch" in result.output
    assert validate_project(root).is_valid


def test_design_without_annotation_requires_gene_id_type_and_writes_nothing(tmp_path):
    root = _filled_l1(tmp_path)
    before = _tree_digest(root)
    result = runner.invoke(app, ["design", str(root), "--yes"])
    assert result.exit_code == 1
    assert "--gene-id-type" in result.output and "--no-enrichment" in result.output
    assert _tree_digest(root) == before

    added = runner.invoke(app, ["design", str(root), "--yes", "--gene-id-type", "entrezid"])
    assert added.exit_code == 0, added.output
    assert _config(root)["annotation"] == {"organism": "Mus musculus", "input_id_type": "ENTREZID"}
    assert _config(root)["analysis"] == {"enrichment": ALL_METHODS}


def test_design_explicit_opt_out_is_recorded_and_kept_on_later_design_runs(tmp_path):
    root = _filled_l1(tmp_path)
    result = runner.invoke(app, ["design", str(root), "--yes", "--no-enrichment"])
    assert result.exit_code == 0, result.output
    assert _config(root)["analysis"] == {"enrichment": []}
    assert "Functional enrichment: none (analysis.enrichment: [])" in result.output
    # Already L2: re-running design never overrides the recorded opt-out.
    again = runner.invoke(app, ["design", str(root), "--yes", "--force"])
    assert again.exit_code == 0, again.output
    assert _config(root)["analysis"] == {"enrichment": []}


def test_design_rejects_a_conflicting_gene_id_type(tmp_path):
    root = _filled_l1(tmp_path, "--gene-id-type", "ENSEMBL")
    result = runner.invoke(app, ["design", str(root), "--yes", "--gene-id-type", "SYMBOL"])
    assert result.exit_code == 1
    assert "already declares annotation.input_id_type ENSEMBL" in result.output


def test_design_for_an_organism_without_annotation_database_requires_opt_out(tmp_path):
    root = _filled_l1(tmp_path)
    config = _config(root)
    config["organism"]["species"] = "Rattus norvegicus"
    (root / "project.yaml").write_text(yaml.safe_dump(config, sort_keys=False), encoding="utf-8")
    result = runner.invoke(app, ["design", str(root), "--yes", "--gene-id-type", "ENSEMBL"])
    assert result.exit_code == 1
    assert "supports only Homo sapiens, Mus musculus" in result.output and "--no-enrichment" in result.output
    assert runner.invoke(app, ["design", str(root), "--yes", "--no-enrichment"]).exit_code == 0


def test_design_enrichment_promotion_never_touches_run_records(tmp_path):
    root = _filled_l1(tmp_path, "--gene-id-type", "ENSEMBL")
    run = root / "runs" / "L1-001" / "20260101-000000+0000"
    (run / "frozen").mkdir(parents=True)
    (run / "run_state.json").write_text(json.dumps({"status": "SUCCESS"}), encoding="utf-8")
    (run / "frozen" / "project.yaml").write_text((root / "project.yaml").read_text(encoding="utf-8"), encoding="utf-8")
    before = _tree_digest(root / "runs")
    assert runner.invoke(app, ["design", str(root), "--yes"]).exit_code == 0
    assert _tree_digest(root / "runs") == before


# ------------------------------------------------------------- validation and preflight


def test_missing_annotation_is_an_actionable_validation_error(project_factory):
    config = base_config()
    config["schema_version"] = "1.3"
    config["design"]["variables"] = {"condition": "categorical"}
    config["analysis"] = {"enrichment": ALL_METHODS}
    report_ = validate_project(project_factory(config=config))
    assert not report_.is_valid
    message = " ".join(issue.message for issue in report_.errors)
    assert "annotation.input_id_type (ENSEMBL, ENTREZID, SYMBOL)" in message
    assert "analysis.enrichment: []" in message


@pytest.mark.parametrize("config_changes", (
    {"schema_version": "1.0"},  # legacy: no analysis section
    {"schema_version": "1.3", "design": {"type": "two_group", "formula": "~ condition", "variables": {"condition": "categorical"}}, "analysis": {"enrichment": []}},
    {"schema_version": "1.1", "annotation": {"organism": "Mus musculus", "input_id_type": "ENSEMBL"}, "analysis": {"enrichment": ["go", "gsea-go", "kegg", "gsea-kegg"]}},
))
def test_existing_projects_load_with_their_recorded_selection(project_factory, config_changes):
    config = {**base_config(), **config_changes}
    report_ = validate_project(project_factory(config=config))
    assert report_.is_valid, report_.errors
    expected = config.get("analysis", {}).get("enrichment") or []
    selected = report_.config.analysis.enrichment if report_.config.analysis else ()
    assert list(selected) == (ALL_METHODS if expected else [])


def _annotated(selection: list[str]) -> ProjectConfig:
    config = base_config()
    config["annotation"] = {"organism": "Mus musculus", "input_id_type": "ENSEMBL"}
    config["schema_version"] = "1.3"
    config["design"]["variables"] = {"condition": "categorical"}
    config["analysis"] = {"enrichment": selection}
    return ProjectConfig.model_validate(config)


def test_enrichment_preflight_probes_kegg_only_when_a_kegg_backend_is_selected(monkeypatch, tmp_path):
    calls: list[str] = []
    monkeypatch.setenv("XDG_CACHE_HOME", str(tmp_path))  # no local KEGG snapshot
    monkeypatch.setattr("rnaseq.service.kegg_rest_preflight", lambda code, **_: calls.append(code))
    _enrichment_preflight(_annotated([]))
    _enrichment_preflight(_annotated(["go"]))
    assert calls == []
    _enrichment_preflight(_annotated(["gsea"]))
    _enrichment_preflight(_annotated(ALL_METHODS))
    assert calls == ["mmu", "mmu"]
    # With a local snapshot, KEGG no longer needs the network.
    snapshot = tmp_path / "nf-rna" / "genesets" / "kegg" / "mmu" / "mmu_fixture"
    snapshot.mkdir(parents=True)
    (snapshot / "manifest.json").write_text("{}", encoding="utf-8")
    _enrichment_preflight(_annotated(ALL_METHODS))
    assert calls == ["mmu", "mmu"]


@pytest.mark.real_kegg_preflight
def test_unreachable_kegg_fails_preflight_with_the_remedy(monkeypatch, tmp_path):
    from urllib.error import URLError

    monkeypatch.setenv("XDG_CACHE_HOME", str(tmp_path))  # no local KEGG snapshot

    def offline(*_args, **_kwargs):
        raise URLError("no route to host")

    monkeypatch.setattr("urllib.request.urlopen", offline)
    with pytest.raises(ExecutionPreflightError) as error:
        _enrichment_preflight(_annotated(ALL_METHODS))
    message = str(error.value)
    assert "https://rest.kegg.jp/list/pathway/mmu" in message and "no route to host" in message
    assert "remove 'kegg' and 'gsea'" in message


@pytest.mark.real_kegg_preflight
def test_reachable_kegg_passes_preflight(monkeypatch):
    monkeypatch.setattr("urllib.request.urlopen", lambda *_a, **_k: io.BytesIO(b"path:mmu00010\tGlycolysis\n"))
    kegg_rest_preflight("mmu")


# ------------------------------------------------------------- execution


def test_default_selection_launches_four_backends_each_covering_every_contrast(project_factory):
    config = base_config()
    config.update({"schema_version": "1.3", "annotation": {"organism": "Mus musculus", "input_id_type": "ENSEMBL"}, "analysis": {"enrichment": ALL_METHODS}})
    config["design"] = {"type": "multi_group", "formula": "~ batch + condition", "variables": {"batch": "categorical", "condition": "categorical"}}
    counts = "gene_id," + ",".join(row.split(",")[0] for row in PARTIALLY_CONFOUNDED.strip().splitlines()[1:]) + "\nG1," + ",".join(["10"] * 10) + "\nG2," + ",".join(["30"] * 10) + "\n"
    root = project_factory(config=config, counts=counts, metadata=PARTIALLY_CONFOUNDED, contrasts=THREE_CONTRASTS)
    validation = validate_project(root)
    assert validation.is_valid, validation.errors
    generate_plan(validation)
    run = create_case_run(validation, "CASE-FOUR", moment=datetime(2026, 10, 14, 12, 0, 0))
    frozen = freeze_case_inputs(validation, run, profile="local", command=["rnaseq", "run"])
    contract = json.loads(frozen.contract.read_text(encoding="utf-8"))
    assert contract["design"]["formula"] == "~ batch + condition"
    assert [item["code"] for item in contract["design_warnings"]].count("partially_confounded_batch") == 1

    command = build_downstream_nextflow_command(run)
    assert command[command.index("--enrichment") + 1] == "go,kegg,gsea-go,gsea-kegg"
    inputs = resolve_downstream_inputs(run).root
    expected = ["high_vs_control", "low_vs_control", "low_vs_high"]
    for kind in ("go", "kegg", "gsea-go", "gsea-kegg"):
        cfg = enrichment_config(frozen.contract, inputs, Path("l2"), kind, Path("out"))
        assert [item["contrast_id"] for item in cfg["contrasts"]] == expected
        assert cfg["orgdb_package"] == "org.Mm.eg.db"
    main_nf = (Path(__file__).resolve().parents[1] / "workflow" / "main.nf").read_text(encoding="utf-8")
    assert "Channel.fromList(modules).combine(l2)" in main_nf


# ------------------------------------------------------------- reporting


def _status_rows(output: Path) -> list[dict[str, str]]:
    with (output / "enrichment_status.tsv").open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def test_report_without_enrichment_lists_every_method_as_not_requested(project_factory):
    contract, contract_path, l2, inputs = _frozen_contract(project_factory, schema_version="1.1")
    contract["analysis"] = {"enrichment": []}
    contract["design_warnings"] = [{"code": "partially_confounded_batch", "message": "Condition level 'low' occurs in only one batch ('b3')."}]
    contract_path.write_text(json.dumps(contract), encoding="utf-8")
    l1, _go, _kegg = _write_report_artifacts(inputs / "contrasts.csv", l2)
    output = l2.parent / "plain-report"
    report(contract_path, inputs, l1, l2, output, [])
    text = (output / "report.html").read_text(encoding="utf-8")
    rows = _status_rows(output)
    assert [row["method"] for row in rows] == ["GO ORA", "KEGG ORA", "GO GSEA", "KEGG GSEA"]
    assert {row["status"] for row in rows} == {"NOT_REQUESTED"}
    assert "L2 — functional enrichment status" in text
    assert "No functional enrichment method was selected" in text
    assert "GO GSEA evaluates BP, MF, and CC" not in text
    assert "Functional enrichment selected: none (analysis.enrichment: [])" in text
    assert "<strong>WARNING:</strong> Condition level &#x27;low&#x27; occurs in only one batch" in text


def test_report_status_matches_backend_summaries_for_every_method(project_factory):
    contract, contract_path, l2, inputs = _frozen_contract(project_factory, schema_version="1.1")
    contract["analysis"] = {"enrichment": ALL_METHODS}
    contract_path.write_text(json.dumps(contract), encoding="utf-8")
    l1, gsea_go, gsea_kegg = _write_report_artifacts(inputs / "contrasts.csv", l2)
    ora_statuses = {"go": ("SUCCESS", "NO_SIGNIFICANT_TERMS", "NOT_APPLICABLE"), "kegg": ("NOT_APPLICABLE", "NOT_APPLICABLE", "NOT_APPLICABLE")}
    dirs = [gsea_go, gsea_kegg]
    for module, summary_name in (("go", "go_backend_summary.json"), ("kegg", "kegg_backend_summary.json")):
        root = l2 / "enrichment" / module
        root.mkdir(parents=True)
        foregrounds = {name: {"status": status, **({"reason": "only 2 mapped foreground genes"} if status == "NOT_APPLICABLE" else {})} for name, status in zip(("significant", "up", "down"), ora_statuses[module])}
        (root / summary_name).write_text(json.dumps({"status": "SUCCESS", "contrasts": [{"contrast_id": "Treatment_vs_Control", "universe": {}, "foregrounds": foregrounds}]}), encoding="utf-8")
        dirs.append(root)
    output = l2.parent / "full-report"
    report(contract_path, inputs, l1, l2, output, dirs)
    rows = {row["method"]: row for row in _status_rows(output)}
    assert rows["GO ORA"]["status"] == "SUCCESS"
    assert rows["KEGG ORA"]["status"] == "NOT_APPLICABLE" and "only 2 mapped" in rows["KEGG ORA"]["reason"]
    go_gsea = json.loads((gsea_go / "gsea_backend_summary.json").read_text(encoding="utf-8"))
    assert rows["GO GSEA"]["status"] in {"SUCCESS", "NO_SIGNIFICANT_TERMS"}
    assert rows["GO GSEA"]["status"] == (go_gsea["contrasts"][0].get("status") or "SUCCESS")
    assert rows["KEGG GSEA"]["status"] in {"SUCCESS", "NO_SIGNIFICANT_TERMS"}
    assert "NOT_REQUESTED" not in {row["status"] for row in rows.values()}
    text = (output / "report.html").read_text(encoding="utf-8")
    assert "No functional enrichment method was selected" not in text
    assert "GO GSEA evaluates BP, MF, and CC" in text


def test_ora_only_report_marks_gsea_not_requested_without_claiming_nothing_ran(project_factory):
    contract, contract_path, l2, inputs = _frozen_contract(project_factory, schema_version="1.1")
    contract["analysis"] = {"enrichment": ["go"]}
    contract_path.write_text(json.dumps(contract), encoding="utf-8")
    l1, _go, _kegg = _write_report_artifacts(inputs / "contrasts.csv", l2)
    root = l2 / "enrichment" / "go"
    root.mkdir(parents=True)
    (root / "go_backend_summary.json").write_text(json.dumps({"status": "NO_SIGNIFICANT_TERMS", "contrasts": [{"contrast_id": "Treatment_vs_Control", "universe": {}, "foregrounds": {"significant": {"status": "NO_SIGNIFICANT_TERMS"}}}]}), encoding="utf-8")
    output = l2.parent / "go-report"
    report(contract_path, inputs, l1, l2, output, [root])
    rows = {row["method"]: row["status"] for row in _status_rows(output)}
    assert rows == {"GO ORA": "NO_SIGNIFICANT_TERMS", "KEGG ORA": "NOT_REQUESTED", "GO GSEA": "NOT_REQUESTED", "KEGG GSEA": "NOT_REQUESTED"}
    text = (output / "report.html").read_text(encoding="utf-8")
    assert "No functional enrichment method was selected" not in text
    assert "Not requested (gsea is not in analysis.enrichment)" in text
    assert "GO GSEA evaluates BP, MF, and CC" not in text


def test_report_requires_the_summary_of_every_selected_method(project_factory):
    contract, contract_path, l2, inputs = _frozen_contract(project_factory, schema_version="1.1")
    contract["analysis"] = {"enrichment": ["kegg"]}
    contract_path.write_text(json.dumps(contract), encoding="utf-8")
    l1, _go, _kegg = _write_report_artifacts(inputs / "contrasts.csv", l2)
    with pytest.raises(ValueError, match="KEGG ORA"):
        report(contract_path, inputs, l1, l2, l2.parent / "missing", [])


def test_legacy_contract_without_design_warnings_says_so(project_factory):
    contract, contract_path, l2, inputs = _frozen_contract(project_factory, schema_version="1.1")
    contract.pop("design_warnings", None)
    contract["analysis"] = {"enrichment": []}
    contract_path.write_text(json.dumps(contract), encoding="utf-8")
    l1, _go, _kegg = _write_report_artifacts(inputs / "contrasts.csv", l2)
    output = l2.parent / "legacy"
    report(contract_path, inputs, l1, l2, output, [])
    assert "Design warnings: not recorded in this run's frozen contract." in (output / "report.html").read_text(encoding="utf-8")


def test_delivery_readme_summarizes_design_warnings_and_status(tmp_path):
    class Run:
        run_dir = tmp_path

    (tmp_path / "frozen").mkdir()
    (tmp_path / "frozen" / "downstream_contract.json").write_text(json.dumps({
        "analysis_level": "L2", "design_warnings": [{"code": "partially_confounded_batch", "message": "Level 'low' occurs in only one batch."}],
    }), encoding="utf-8")
    assert "Not available" in _delivery_readme_analysis_summary(Run())
    (tmp_path / "downstream" / "report").mkdir(parents=True)
    (tmp_path / "downstream" / "report" / "enrichment_status.tsv").write_text(
        "method\tselection\tcontrast_id\tstatus\treason\nGO ORA\tgo\tlow_vs_control\tSUCCESS\t\nKEGG ORA\tkegg\tlow_vs_control\tNOT_REQUESTED\tx\n",
        encoding="utf-8",
    )
    text = _delivery_readme_analysis_summary(Run())
    assert "- Level 'low' occurs in only one batch." in text
    assert "| GO ORA | low_vs_control | SUCCESS |" in text
    assert "| KEGG ORA | low_vs_control | NOT_REQUESTED |" in text

    (tmp_path / "frozen" / "downstream_contract.json").write_text(json.dumps({"analysis_level": "L1"}), encoding="utf-8")
    assert _delivery_readme_analysis_summary(Run()) == ""


# ------------------------------------------------------------- v1.5 enrichment warnings


def _enrichment_project(project_factory, *, enrichment: dict | None = None, selection: list[str] = ALL_METHODS):
    config = base_config()
    config["schema_version"] = "1.3"
    config["design"]["variables"] = {"condition": "categorical"}
    config["annotation"] = {"organism": "Mus musculus", "input_id_type": "ENTREZID", **({"enrichment": enrichment} if enrichment else {})}
    config["analysis"] = {"enrichment": selection}
    return validate_project(project_factory(config=config))


def test_legacy_cutoff_fields_warn_but_do_not_fail_validation(project_factory, monkeypatch, tmp_path):
    monkeypatch.setenv("XDG_CACHE_HOME", str(tmp_path))
    report = _enrichment_project(project_factory, enrichment={"go": {"qvalue_cutoff": 0.1}, "kegg": {"gsea": {"padj_cutoff": 0.01}}}, selection=["go"])
    warnings = [issue for issue in report.warnings if issue.code == "legacy_enrichment_cutoff"]
    assert report.is_valid and len(warnings) == 1
    assert "annotation.enrichment.go.qvalue_cutoff, annotation.enrichment.kegg.gsea.padj_cutoff" in warnings[0].message
    assert "p.adjust (Benjamini-Hochberg) <= 0.05" in warnings[0].message
    assert not [issue for issue in _enrichment_project(project_factory, selection=["go"]).warnings if issue.code == "legacy_enrichment_cutoff"]


def test_unpinned_kegg_snapshot_is_named_and_a_missing_pinned_one_is_an_error(project_factory, monkeypatch, tmp_path):
    monkeypatch.setenv("XDG_CACHE_HOME", str(tmp_path))
    report = _enrichment_project(project_factory)
    [warning] = [issue for issue in report.warnings if issue.code == "kegg_snapshot_not_pinned"]
    assert "no local snapshot exists" in warning.message
    root = tmp_path / "nf-rna" / "genesets" / "kegg" / "mmu"
    for snapshot_id, stamp in (("mmu_a", "2026-01-01T00:00:00Z"), ("mmu_b", "2026-05-01T00:00:00Z")):
        (root / snapshot_id).mkdir(parents=True)
        (root / snapshot_id / "manifest.json").write_text(json.dumps({"retrieved_at_utc": stamp, "kegg_release": "KEGG pathway database: Last update 2026/04/30"}), encoding="utf-8")
    [warning] = [issue for issue in _enrichment_project(project_factory).warnings if issue.code == "kegg_snapshot_not_pinned"]
    assert warning.message.startswith("KEGG snapshot: mmu_b (not pinned") and "annotation.enrichment.kegg.snapshot: mmu_b" in warning.message
    pinned = _enrichment_project(project_factory, enrichment={"kegg": {"snapshot": "mmu_a"}})
    assert not [issue for issue in pinned.warnings if issue.code == "kegg_snapshot_not_pinned"] and pinned.is_valid
    missing = _enrichment_project(project_factory, enrichment={"kegg": {"snapshot": "mmu_absent"}})
    assert [issue.code for issue in missing.errors] == ["kegg_snapshot_missing"]
    assert not [issue for issue in _enrichment_project(project_factory, selection=["go"]).warnings if issue.code == "kegg_snapshot_not_pinned"]


def test_report_figures_are_shown_or_explicitly_unavailable_with_a_reason(tmp_path):
    from rnaseq.workflow_support import _image_or_unavailable, _running_score_html

    line = _image_or_unavailable(tmp_path / "ma_plot.png", "MA plot: c", "ma_plot.png was not produced")
    assert "MA plot: c: unavailable" in line and "ma_plot.png was not produced" in line and "<img" not in line
    record = {"figures": {"running_score": {"drawn": 0, "reason": "no significant gene sets (p.adjust <= 0.05)"}}}
    parts = _running_score_html(tmp_path, record, "KEGG GSEA")
    assert len(parts) == 2 and all("unavailable" in part and "no significant gene sets" in part for part in parts)
    (tmp_path / "running_score").mkdir()
    for name in ("positive_1_A.png", "positive_1_A_300dpi.png", "positive_2_B.png", "positive_3_C.png"):
        (tmp_path / "running_score" / name).write_bytes(b"\x89PNG")
    parts = _running_score_html(tmp_path, {"figures": {"running_score": {"drawn": 3, "positive": 3, "negative": 0}}}, "KEGG GSEA")
    assert sum("<img" in part for part in parts) == 2 and "positive_3_C" not in "".join(parts)
    assert "no significant gene set with negative NES" in parts[-1]
