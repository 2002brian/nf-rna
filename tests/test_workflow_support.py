from __future__ import annotations

import json
import os
import re
import shutil
import subprocess
import sys
from copy import deepcopy
from datetime import datetime
from pathlib import Path

import pytest

from conftest import base_config
from rnaseq.planner import generate_plan
from rnaseq.service import create_case_run, freeze_case_inputs, resolve_downstream_inputs
from rnaseq.validators import validate_project
from rnaseq.workflow_support import _samples, _source_import_label, enrichment_config, l1_config, l2_config, report


@pytest.mark.parametrize(
    ("source", "label"),
    (
        ("salmon_tximport", "Salmon/tximport import"),
        ("featurecounts_raw_counts", "featureCounts raw-count matrix via DESeqDataSetFromMatrix"),
        ("raw_counts", "Imported raw-count matrix via DESeqDataSetFromMatrix"),
    ),
)
def test_report_source_labels_match_actual_import(source, label):
    assert _source_import_label(source) == label


ROOT = Path(__file__).parents[1]


def _annotation() -> dict[str, object]:
    return {
        "organism": "Mus musculus",
        "input_id_type": "ENSEMBL",
        "target_id_type": "ENTREZID",
        "gene_symbol_output": True,
        "mapping_warning_rate": 0.83,
        "minimum_mapping_rate": 0.50,
        "minimum_mapped_foreground": 7,
        "enrichment": {
            "go": {"pvalue_cutoff": 0.031, "qvalue_cutoff": 0.17, "p_adjust_method": "BH"},
            "gsea": {"minimum_ranked_genes": 61, "min_gs_size": 11, "max_gs_size": 321, "pvalue_cutoff": 0.041, "padj_cutoff": 0.049, "p_adjust_method": "BH", "seed": 9},
            "kegg": {
                "resource_provider": "online_kegg_rest_via_clusterprofiler",
                "ora": {"pvalue_cutoff": 0.021, "qvalue_cutoff": 0.19, "p_adjust_method": "BH", "min_gs_size": 12, "max_gs_size": 322},
                "gsea": {"minimum_ranked_genes": 62, "pvalue_cutoff": 0.91, "padj_cutoff": 0.039, "p_adjust_method": "BH", "min_gs_size": 13, "max_gs_size": 323, "seed": 10},
            },
        },
    }


def _frozen_contract(project_factory, *, schema_version: str) -> tuple[dict[str, object], Path, Path, Path]:
    config = base_config()
    config["schema_version"] = schema_version
    config["annotation"] = _annotation()
    if schema_version == "1.1":
        config["analysis"] = {"enrichment": "gsea"}
    root = project_factory(config=config)
    report = validate_project(root)
    assert report.is_valid, report.errors
    generate_plan(report)
    run = create_case_run(report, "CASE-CONTRACT", moment=datetime(2026, 8, 28, 12, 0, 0))
    frozen = freeze_case_inputs(report, run, profile="local", command=["rnaseq", "run"])
    inputs = resolve_downstream_inputs(run)
    return (
        json.loads(frozen.contract.read_text(encoding="utf-8")),
        frozen.contract,
        run.run_dir / "downstream" / "l2",
        inputs.root,
    )


def _active_l2_contract(project_factory, *, counts: str, metadata: str) -> tuple[Path, Path, object]:
    """Freeze an ordinary L2 project for the production config-builder route."""

    root = project_factory(counts=counts, metadata=metadata)
    validation = validate_project(root)
    assert validation.is_valid, validation.errors
    generate_plan(validation)
    run = create_case_run(validation, "CASE-L2-CONTRACT", moment=datetime(2026, 9, 12, 12, 0, 0))
    frozen = freeze_case_inputs(validation, run, profile="local", command=["rnaseq", "run"])
    return frozen.contract, resolve_downstream_inputs(run).root, validation


@pytest.mark.parametrize(
    ("counts", "metadata", "expected_counts"),
    (
        (
            "gene_id,C1,T1\nGeneA,10,40\nGeneB,100,95\n",
            "sample_id,condition\nC1,Control\nT1,Treatment\n",
            "Control n=1, Treatment n=1",
        ),
        (
            "gene_id,C1,T1,T2\nGeneA,10,40,45\nGeneB,100,95,102\n",
            "sample_id,condition\nC1,Control\nT1,Treatment\nT2,Treatment\n",
            "Control n=1, Treatment n=2",
        ),
    ),
)
def test_active_l2_config_blocks_unreplicated_contrasts_but_l1_remains_available(
    project_factory, tmp_path, counts, metadata, expected_counts,
):
    contract, inputs, validation = _active_l2_contract(project_factory, counts=counts, metadata=metadata)
    assert validation.is_valid
    assert l1_config(contract, inputs, tmp_path / "l1")["output_dir"] == str(tmp_path / "l1")
    with pytest.raises(ValueError, match=expected_counts):
        l2_config(contract, inputs, tmp_path / "l1", tmp_path / "l2")


def test_active_l2_config_accepts_two_by_two_with_existing_low_replication_warnings(project_factory, tmp_path):
    counts = "gene_id,C1,C2,T1,T2\nGeneA,10,12,40,45\nGeneB,100,110,95,102\n"
    metadata = "sample_id,condition\nC1,Control\nC2,Control\nT1,Treatment\nT2,Treatment\n"
    contract, inputs, validation = _active_l2_contract(project_factory, counts=counts, metadata=metadata)
    assert sum(issue.code == "limited_replication" for issue in validation.warnings) == 2
    config = l2_config(contract, inputs, tmp_path / "l1", tmp_path / "l2")
    assert config["contrasts"] == [
        {"contrast_id": "Treatment_vs_Control", "factor": "condition", "numerator": "Treatment", "denominator": "Control"}
    ]


def test_active_nextflow_l2_process_uses_the_guarded_config_builder():
    workflow = (ROOT / "workflow" / "main.nf").read_text(encoding="utf-8")
    assert "python -m rnaseq.workflow_support l2-config" in workflow


def _write_report_artifacts(contrasts_path: Path, l2: Path) -> tuple[Path, Path, Path]:
    l1 = l2.parent / "l1"
    l1.mkdir(parents=True)
    (l1 / "backend_summary.json").write_text(json.dumps({
        "genes_input": 120, "genes_removed_all_zero": 10, "genes_removed_low_total": 5, "genes_retained": 105,
        "filter": {"rule": "fixture filter"}, "normalization": {"method": "fixture normalization"},
        "vst": {"method": "fixture VST"}, "package_versions": {"R": "fixture R", "DESeq2": "fixture DESeq2"},
    }), encoding="utf-8")
    (l1 / "pca_variance.tsv").write_text("component\tproportion_variance\nPC1\t0.4\nPC2\t0.2\n", encoding="utf-8")
    (l1 / "library_size_qc.tsv").write_text("sample_id\tinput_total\tretained_total\tsize_factor\tnormalized_total\nC1\t100\t90\t1\t90\n", encoding="utf-8")
    (l1 / "pca.png").write_bytes(b"png")
    (l1 / "sample_correlation.png").write_bytes(b"png")
    with contrasts_path.open(encoding="utf-8", newline="") as handle:
        contrasts = list(__import__("csv").DictReader(handle))
    for contrast in contrasts:
        root = l2 / "contrasts" / contrast["contrast_id"]
        root.mkdir(parents=True)
        (root / "backend_summary.json").write_text(json.dumps({"tested_genes": 100}), encoding="utf-8")
        (root / "all_genes.tsv").write_text("gene_id\tstat\nGene1\t2\n", encoding="utf-8")
        (root / "significant.tsv").write_text("gene_id\tstat\nGene1\t2\n", encoding="utf-8")
        (root / "upregulated.tsv").write_text("gene_id\tstat\nGene1\t2\n", encoding="utf-8")
        (root / "downregulated.tsv").write_text("gene_id\tstat\n", encoding="utf-8")
        (root / "volcano.png").write_bytes(b"png")
        (root / "heatmap.png").write_bytes(b"png")

    go_root, kegg_root = l2 / "enrichment" / "gsea_go", l2 / "enrichment" / "gsea_kegg"
    go_contrasts, kegg_contrasts = [], []
    for contrast in contrasts:
        contrast_id = contrast["contrast_id"]
        ranking = {
            "finite_stat_source_genes": 100, "mapped_source_genes": 80, "unmapped_source_genes": 20,
            "mapping_rate": 0.8, "final_ranked_genes": 70, "positive_stats": 35, "negative_stats": 35,
            "tie_handling": "fixture tie handling",
        }
        ontologies: dict[str, object] = {}
        for ontology in ("BP", "MF", "CC"):
            root = go_root / contrast_id / ontology
            root.mkdir(parents=True)
            (root / "all_terms.tsv").write_text("ID\tDescription\tNES\tpvalue\tp.adjust\nGO:1\tFixture term\t1.2\t0.01\t0.02\n", encoding="utf-8")
            (root / "dotplot.png").write_bytes(b"png")
            ontologies[ontology] = {"status": "SUCCESS", "all_terms": 4, "significant_terms": 2, "positive_terms": 1, "negative_terms": 1}
        go_contrasts.append({"contrast_id": contrast_id, "status": "SUCCESS", "ranking": ranking, "ontologies": ontologies})
        root = kegg_root / contrast_id
        root.mkdir(parents=True)
        (root / "all_terms.tsv").write_text(
            "ID\tDescription\tsetSize\tenrichmentScore\tNES\tpvalue\tp.adjust\tqvalue\trank\tleading_edge\tcore_enrichment\n"
            "mmu05171\tCoronavirus disease\t143\t0.63\t2.454\t0.0001\t0.002\t0.002\t87\ttags=40%\t101;102\n",
            encoding="utf-8",
        )
        (root / "dotplot.png").write_bytes(b"png")
        kegg_contrasts.append({"contrast_id": contrast_id, "status": "SUCCESS", "ranking": ranking, "all_terms": 3, "significant_terms": 1, "positive_terms": 1, "negative_terms": 0})
    (go_root / "gsea_backend_summary.json").write_text(json.dumps({"status": "SUCCESS", "contrasts": go_contrasts}), encoding="utf-8")
    (kegg_root / "gsea_kegg_backend_summary.json").write_text(json.dumps({"status": "SUCCESS", "resource": {"provider": "fixture KEGG"}, "contrasts": kegg_contrasts}), encoding="utf-8")
    return l1, go_root, kegg_root


def test_schema_11_freezes_complete_annotation_and_enrichment_bridge_is_deterministic(project_factory):
    contract, contract_path, l2, inputs = _frozen_contract(project_factory, schema_version="1.1")
    annotation = contract["annotation"]
    assert annotation == _annotation()

    l1, go_root, kegg_root = _write_report_artifacts(inputs / "contrasts.csv", l2)
    report_dir = l2.parent / "report"
    report(contract_path, inputs, l1, l2, report_dir, [go_root, kegg_root])
    report_text = (report_dir / "report.html").read_text(encoding="utf-8")
    assert "Enrichment method: GSEA" in report_text
    assert "L1 — quality control" in report_text
    assert "Treatment_vs_Control" in report_text
    assert "GO GSEA" in report_text and "BP" in report_text and "MF" in report_text and "CC" in report_text
    assert "KEGG GSEA" in report_text and "fixture KEGG" in report_text
    assert "Genes input: 120" in report_text and "Final ranked genes: 70" in report_text
    assert "tables/l2/enrichment/gsea_go" in report_text and "tables/l2/enrichment/gsea_kegg" in report_text
    assert "GO:1" in report_text and "Fixture term" in report_text
    assert "Significant genes: 1" in report_text
    assert "mmu05171" in report_text and "Coronavirus disease" in report_text
    assert "2.454" in report_text and "0.0001" in report_text and "0.002" in report_text
    assert "max-width: 100%;" in report_text
    assert "height: auto;" in report_text
    assert "object-fit: contain;" in report_text
    assert "overflow: hidden" not in report_text.lower()
    assert re.search(r"height\s*:\s*[0-9]", report_text) is None
    image_tags = re.findall(r"<img\b[^>]*>", report_text)
    assert image_tags
    assert all("class='report-figure'" in tag and "style=" not in tag for tag in image_tags)
    assert report_text.count("class='report-figure-container'") == len(image_tags)
    assert "@media print" in report_text and "page-break-inside: avoid" in report_text

    # The bridge must not consult frozen project.yaml after run creation.
    frozen_project = Path(contract["project_config"])
    frozen_project.write_text("annotation:\n  minimum_mapping_rate: 0.01\n", encoding="utf-8")
    configs = {kind: enrichment_config(contract_path, inputs, l2, kind, l2.parent) for kind in ("gsea-go", "gsea-kegg")}
    assert enrichment_config(contract_path, inputs, l2, "gsea-go", l2.parent) == configs["gsea-go"]
    assert json.dumps(configs["gsea-go"], sort_keys=True) == json.dumps(enrichment_config(contract_path, inputs, l2, "gsea-go", l2.parent), sort_keys=True)
    assert all(cfg["annotation"] == annotation for cfg in configs.values())

    assert configs["gsea-go"]["annotation"]["enrichment"]["gsea"] == _annotation()["enrichment"]["gsea"]
    assert configs["gsea-kegg"]["annotation"]["enrichment"]["kegg"]["gsea"] == _annotation()["enrichment"]["kegg"]["gsea"]
    assert configs["gsea-kegg"]["kegg"]["provider"] == "online_kegg_rest_via_clusterprofiler"
    assert configs["gsea-kegg"]["kegg"]["organism_code"] == "mmu"


def test_final_report_requires_both_enabled_gsea_backends_and_renders_multiple_contrasts(project_factory):
    _contract, contract_path, l2, inputs = _frozen_contract(project_factory, schema_version="1.1")
    contrasts_path = inputs / "contrasts.csv"
    contrasts_path.write_text(
        "contrast_id,factor,numerator,denominator\nTreatment_vs_Control,condition,Treatment,Control\nControl_vs_Treatment,condition,Control,Treatment\n",
        encoding="utf-8",
    )
    l1, go_root, kegg_root = _write_report_artifacts(contrasts_path, l2)
    with pytest.raises(ValueError, match="gsea_kegg"):
        report(contract_path, inputs, l1, l2, l2.parent / "missing-report", [go_root])
    report(contract_path, inputs, l1, l2, l2.parent / "report", [go_root, kegg_root])
    text = (l2.parent / "report" / "report.html").read_text(encoding="utf-8")
    assert text.count("GO GSEA:") == 2
    assert text.count("KEGG GSEA:") == 2
    assert "Control_vs_Treatment" in text


def test_final_report_renders_annotation_warning_without_failing_enabled_gsea(project_factory):
    _contract, contract_path, l2, inputs = _frozen_contract(project_factory, schema_version="1.1")
    l1, go_root, kegg_root = _write_report_artifacts(inputs / "contrasts.csv", l2)
    for root, filename in ((go_root, "gsea_backend_summary.json"), (kegg_root, "gsea_kegg_backend_summary.json")):
        path = root / filename
        summary = json.loads(path.read_text(encoding="utf-8"))
        summary["annotation"] = {
            "input_id_type": "ENSEMBL",
            "target_id_type": "ENTREZID",
            "warning_threshold": 0.70,
            "blocking_threshold": 0.50,
        }
        summary["annotation_database"] = "org.Mm.eg.db"
        summary["annotation_database_version"] = "fixture"
        summary["annotation_qc_status"] = "WARNING"
        for contrast in summary["contrasts"]:
            contrast["ranking"]["mapping_rate"] = 0.697
            contrast["ranking"]["annotation_qc"] = {
                "status": "WARNING",
                "mapping_rate": 0.697,
                "warning_threshold": 0.70,
                "blocking_threshold": 0.50,
                "reason": "mapping rate 69.7% is below warning threshold 70.0%; GSEA was executed.",
            }
        path.write_text(json.dumps(summary), encoding="utf-8")

    report_dir = l2.parent / "warning-report"
    report(contract_path, inputs, l1, l2, report_dir, [go_root, kegg_root])
    report_text = (report_dir / "report.html").read_text(encoding="utf-8")
    assert "Annotation mapping QC: <strong>WARNING</strong>" in report_text
    assert "GSEA was executed." in report_text


def test_final_report_distinguishes_evaluated_from_significant_gsea_terms(project_factory):
    _contract, contract_path, l2, inputs = _frozen_contract(project_factory, schema_version="1.1")
    l1, go_root, kegg_root = _write_report_artifacts(inputs / "contrasts.csv", l2)
    for root, filename in ((go_root, "gsea_backend_summary.json"), (kegg_root, "gsea_kegg_backend_summary.json")):
        path = root / filename
        summary = json.loads(path.read_text(encoding="utf-8"))
        for contrast in summary["contrasts"]:
            contrast["evaluated_terms"] = 100
            contrast["all_terms"] = 100
            contrast["significant_terms"] = 20
            if "ontologies" in contrast:
                for outcome in contrast["ontologies"].values():
                    outcome["evaluated_terms"] = 100
                    outcome["all_terms"] = 100
                    outcome["significant_terms"] = 20
        path.write_text(json.dumps(summary), encoding="utf-8")

    report_dir = l2.parent / "evaluated-terms-report"
    report(contract_path, inputs, l1, l2, report_dir, [go_root, kegg_root])
    report_text = (report_dir / "report.html").read_text(encoding="utf-8")
    assert "Evaluated terms: 100" in report_text
    assert "Evaluated pathways: 100" in report_text
    assert "Significant terms: 20" in report_text
    assert "Significant pathways: 20" in report_text


def test_final_report_remains_fail_closed_for_blocked_enabled_gsea_with_reason(project_factory):
    _contract, contract_path, l2, inputs = _frozen_contract(project_factory, schema_version="1.1")
    l1, go_root, kegg_root = _write_report_artifacts(inputs / "contrasts.csv", l2)
    summary_path = go_root / "gsea_backend_summary.json"
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    summary.update({
        "status": "BLOCKED",
        "reason": "GO preranked GSEA blocked: mapping rate 49.0% is below blocking threshold 50.0%.",
    })
    summary_path.write_text(json.dumps(summary), encoding="utf-8")

    with pytest.raises(ValueError, match="below blocking threshold 50.0%"):
        report(contract_path, inputs, l1, l2, l2.parent / "blocked-report", [go_root, kegg_root])


def test_report_cli_accepts_all_enrichment_paths_and_rejects_incomplete_or_unknown_inputs(project_factory):
    _contract, contract_path, l2, inputs = _frozen_contract(project_factory, schema_version="1.1")
    l1, go_root, kegg_root = _write_report_artifacts(inputs / "contrasts.csv", l2)
    environment = {**os.environ, "PYTHONPATH": str(ROOT / "src")}
    base = [
        sys.executable, "-m", "rnaseq.workflow_support", "report",
        "--contract", str(contract_path), "--inputs", str(inputs), "--l1", str(l1), "--l2", str(l2),
    ]
    complete = subprocess.run(
        [*base, "--enrichment", str(go_root), str(kegg_root), "--out", str(l2.parent / "cli-report")],
        env=environment, capture_output=True, text=True, check=False,
    )
    assert complete.returncode == 0, complete.stderr
    assert (l2.parent / "cli-report" / "report.html").is_file()

    incomplete = subprocess.run(
        [*base, "--enrichment", str(go_root), "--out", str(l2.parent / "incomplete-report")],
        env=environment, capture_output=True, text=True, check=False,
    )
    assert incomplete.returncode != 0
    assert "gsea_kegg" in incomplete.stderr

    unknown = subprocess.run(
        [*base, "--unexpected", "value", "--out", str(l2.parent / "unknown-report")],
        env=environment, capture_output=True, text=True, check=False,
    )
    assert unknown.returncode == 2
    assert "unrecognized arguments" in unknown.stderr


def test_report_cli_without_enrichment_remains_supported(project_factory):
    _contract, contract_path, l2, inputs = _frozen_contract(project_factory, schema_version="1.0")
    l1, _go_root, _kegg_root = _write_report_artifacts(inputs / "contrasts.csv", l2)
    result = subprocess.run(
        [
            sys.executable, "-m", "rnaseq.workflow_support", "report",
            "--contract", str(contract_path), "--inputs", str(inputs), "--l1", str(l1), "--l2", str(l2),
            "--out", str(l2.parent / "no-enrichment-report"),
        ],
        env={**os.environ, "PYTHONPATH": str(ROOT / "src")}, capture_output=True, text=True, check=False,
    )
    assert result.returncode == 0, result.stderr
    assert (l2.parent / "no-enrichment-report" / "report.html").is_file()


def test_l1_report_requires_no_l2_artifacts_and_cli_does_not_require_l2(project_factory):
    contract, contract_path, l2, inputs = _frozen_contract(project_factory, schema_version="1.0")
    contract["analysis_level"] = "L1"
    contract_path.write_text(json.dumps(contract, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    l1, _go_root, _kegg_root = _write_report_artifacts(inputs / "contrasts.csv", l2)

    output = l2.parent / "l1-only-report"
    report(contract_path, inputs, l1, None, output)
    report_text = (output / "report.html").read_text(encoding="utf-8")
    assert "Requested analysis level: L1" in report_text
    assert "Not requested for this L1 project." in report_text
    assert "L2 summary" not in report_text

    result = subprocess.run(
        [
            sys.executable, "-m", "rnaseq.workflow_support", "report",
            "--contract", str(contract_path), "--inputs", str(inputs), "--l1", str(l1),
            "--out", str(l2.parent / "l1-cli-report"),
        ],
        env={**os.environ, "PYTHONPATH": str(ROOT / "src")}, capture_output=True, text=True, check=False,
    )
    assert result.returncode == 0, result.stderr
    assert (l2.parent / "l1-cli-report" / "report.html").is_file()


def test_report_rejects_invalid_frozen_analysis_level_before_reading_r_artifacts(project_factory):
    contract, contract_path, l2, inputs = _frozen_contract(project_factory, schema_version="1.0")
    contract["analysis_level"] = "L3"
    contract_path.write_text(json.dumps(contract, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    with pytest.raises(ValueError, match="frozen analysis_level L1 or L2"):
        report(contract_path, inputs, l2.parent / "missing-l1", None, l2.parent / "invalid-level-report")


def test_report_rejects_tsv_rows_without_any_requested_display_value(project_factory):
    _contract, contract_path, l2, inputs = _frozen_contract(project_factory, schema_version="1.1")
    l1, go_root, kegg_root = _write_report_artifacts(inputs / "contrasts.csv", l2)
    contrast_id = "Treatment_vs_Control"
    (kegg_root / contrast_id / "all_terms.tsv").write_text(
        "ID\tDescription\tNES\tpvalue\tp.adjust\n\t\t\t\t\n",
        encoding="utf-8",
    )
    with pytest.raises(ValueError, match="no values for requested report columns"):
        report(contract_path, inputs, l1, l2, l2.parent / "blank-kegg-report", [go_root, kegg_root])


def test_report_pipeline_provenance_fallback_is_not_duplicated(project_factory):
    _contract, contract_path, l2, inputs = _frozen_contract(project_factory, schema_version="1.1")
    l1, go_root, kegg_root = _write_report_artifacts(inputs / "contrasts.csv", l2)
    execution_manifest = inputs / "execution_manifest.yaml"
    execution_manifest.unlink()
    report(contract_path, inputs, l1, l2, l2.parent / "provenance-fallback", [go_root, kegg_root])
    report_text = (l2.parent / "provenance-fallback" / "report.html").read_text(encoding="utf-8")
    assert "Pipeline: bulk_rnaseq; version: not available." in report_text
    assert "not available not available" not in report_text


def test_enrichment_bridge_fails_before_r_when_a_required_frozen_field_is_missing(project_factory):
    contract, contract_path, l2, inputs = _frozen_contract(project_factory, schema_version="1.1")
    broken = deepcopy(contract)
    del broken["annotation"]["minimum_mapping_rate"]
    contract_path.write_text(json.dumps(broken, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    try:
        enrichment_config(contract_path, inputs, l2, "gsea-go", l2.parent)
    except ValueError as error:
        assert str(error) == "frozen enrichment configuration is missing annotation.minimum_mapping_rate"
    else:
        raise AssertionError("expected missing frozen enrichment configuration to fail before R launch")


def test_schema_10_annotation_contract_remains_compatible(project_factory):
    contract, contract_path, l2, inputs = _frozen_contract(project_factory, schema_version="1.0")
    assert contract["annotation"] == _annotation()
    cfg = enrichment_config(contract_path, inputs, l2, "gsea-go", l2.parent)
    assert cfg["annotation"]["enrichment"]["gsea"]["pvalue_cutoff"] == 0.041


def test_l1_bridge_uses_only_task_staged_execution_paths(project_factory, tmp_path):
    contract, contract_path, _l2, inputs = _frozen_contract(project_factory, schema_version="1.0")
    task_inputs = tmp_path / "nextflow-task" / "inputs"
    shutil.copytree(inputs, task_inputs)
    task_contract = tmp_path / "nextflow-task" / "downstream_contract.json"
    hidden_provenance = deepcopy(contract)
    hidden_provenance["source"]["counts"] = "/Volumes/KOXIA/not-mounted/counts.csv"
    task_contract.write_text(json.dumps(hidden_provenance, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    config = l1_config(task_contract, task_inputs, Path("l1"))

    assert Path(config["counts"]) == task_inputs / "source" / "counts.csv"
    assert Path(config["metadata"]) == task_inputs / "metadata.csv"
    assert "/Volumes/KOXIA" not in json.dumps(config, sort_keys=True)
    assert json.loads(contract_path.read_text(encoding="utf-8"))["source"]["counts"] != hidden_provenance["source"]["counts"]

    escaped_manifest = json.loads((task_inputs / "execution_inputs.json").read_text(encoding="utf-8"))
    escaped_manifest["source"]["counts"] = "../outside-counts.csv"
    (task_inputs.parent / "outside-counts.csv").write_text("gene_id,C1\nGeneA,1\n", encoding="utf-8")
    (task_inputs / "execution_inputs.json").write_text(
        json.dumps(escaped_manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    with pytest.raises(ValueError, match="raw-count matrix escapes its task directory"):
        l1_config(task_contract, task_inputs, Path("l1"))


def test_samples_preserves_non_alphabetical_staged_metadata_order(tmp_path):
    inputs = tmp_path / "downstream_inputs"
    inputs.mkdir()

    samples = ["C1", "T1", "C2", "T2", "C3", "T3"]

    (inputs / "metadata.csv").write_text(
        "sample_id,condition\n"
        "C1,Control\n"
        "T1,Treatment\n"
        "C2,Control\n"
        "T2,Treatment\n"
        "C3,Control\n"
        "T3,Treatment\n",
        encoding="utf-8",
    )

    (inputs / "execution_inputs.json").write_text(
        json.dumps(
            {
                "schema_version": "1.0",
                "metadata": "metadata.csv",
                "samples": samples,
            }
        ),
        encoding="utf-8",
    )

    assert _samples(inputs) == samples


def test_samples_rejects_metadata_order_different_from_staged_sample_order(tmp_path):
    inputs = tmp_path / "downstream_inputs"
    inputs.mkdir()

    samples = ["C1", "T1", "C2", "T2"]

    (inputs / "metadata.csv").write_text(
        "sample_id,condition\n"
        "C1,Control\n"
        "C2,Control\n"
        "T1,Treatment\n"
        "T2,Treatment\n",
        encoding="utf-8",
    )

    (inputs / "execution_inputs.json").write_text(
        json.dumps(
            {
                "schema_version": "1.0",
                "metadata": "metadata.csv",
                "samples": samples,
            }
        ),
        encoding="utf-8",
    )

    with pytest.raises(
        ValueError,
        match="staged metadata sample IDs disagree with staged execution inputs sample list",
    ):
        _samples(inputs)
