from __future__ import annotations

import json
from copy import deepcopy
from pathlib import Path
from types import SimpleNamespace

import pytest
import yaml

from conftest import base_config, require_r_packages
from rnaseq.errors import DownstreamExecutionError
from rnaseq.go import execute_go, prepare_go
from rnaseq.l2 import prepare_l2
from rnaseq.validators import validate_project


def _annotation() -> dict[str, object]:
    return {
        "organism": "Mus musculus", "input_id_type": "ENTREZID", "target_id_type": "ENTREZID",
        "gene_symbol_output": True, "minimum_mapping_rate": 0.70, "minimum_mapped_foreground": 5,
        "enrichment": {"go": {"pvalue_cutoff": 0.05, "qvalue_cutoff": 0.2, "p_adjust_method": "BH"}},
    }


def test_annotation_contract_is_optional_but_must_match_project_organism(project_factory):
    config = base_config()
    config["annotation"] = _annotation()
    assert validate_project(project_factory(config=config)).is_valid
    config = deepcopy(config)
    config["annotation"]["organism"] = "Homo sapiens"
    report = validate_project(project_factory(config=config))
    assert any("annotation.organism must match" in issue.message for issue in report.errors)


def test_go_requires_explicit_contract_and_rejects_unsupported_organism(project_factory):
    with pytest.raises(DownstreamExecutionError, match="explicit annotation contract"):
        prepare_go(validate_project(project_factory()), run_id=None)
    config = base_config()
    config["organism"] = {"species": "Rattus norvegicus"}
    with pytest.raises(DownstreamExecutionError, match="not currently supported"):
        prepare_go(validate_project(project_factory(config=config)), run_id=None)


def test_go_requires_existing_l2_success(project_factory):
    config = base_config()
    config["annotation"] = _annotation()
    with pytest.raises(DownstreamExecutionError, match="existing successful L2"):
        prepare_go(validate_project(project_factory(config=config)), run_id=None)


def test_go_blocked_state_is_separate_from_l2(monkeypatch, project_factory):
    config = base_config()
    config["annotation"] = _annotation()
    root = project_factory(config=config)
    l2 = prepare_l2(validate_project(root), run_id=None)
    l2.output_dir.mkdir(parents=True)
    (l2.output_dir / "l2_state.json").write_text(json.dumps({"status": "SUCCESS"}), encoding="utf-8")
    contrast = l2.output_dir / "contrasts" / "Treatment_vs_Control"
    contrast.mkdir(parents=True)
    for name in ("all_genes.tsv", "significant.tsv", "upregulated.tsv", "downregulated.tsv"):
        (contrast / name).write_text("gene_id\nGeneA\n", encoding="utf-8")
    prepared = prepare_go(validate_project(root), run_id=None)
    monkeypatch.setattr("rnaseq.go._check_go_runtime", lambda _organism: None)

    def fake_backend(arguments, **_kwargs):
        output = Path(arguments[-1]).parent
        annotation = output.parent.parent / "annotation"
        annotation.mkdir(exist_ok=True)
        (annotation / "gene_mapping.tsv").write_text("original_gene_id\nGeneA\n", encoding="utf-8")
        (output / "go_backend_summary.json").write_text(json.dumps({"status": "BLOCKED", "reason": "low mapping", "mapping": {"tested": {"source_genes": 1, "mapped_source_genes": 0, "unmapped_source_genes": 1, "one_to_many_source_ids": 0, "duplicate_target_ids": 0, "unique_target_genes": 0, "mapping_rate": 0}}, "contrasts": [], "annotation_database": "org.Mm.eg.db", "annotation_database_version": "test", "clusterProfiler_version": "test"}), encoding="utf-8")
        return SimpleNamespace(returncode=0)

    monkeypatch.setattr("rnaseq.go.subprocess.run", fake_backend)
    result = execute_go(prepared)
    assert result.status == "BLOCKED"
    assert json.loads((l2.output_dir / "l2_state.json").read_text())["status"] == "SUCCESS"
    assert json.loads(result.state_path.read_text())["status"] == "BLOCKED"


@pytest.mark.parametrize("status", ["NO_SIGNIFICANT_TERMS", "NOT_APPLICABLE"])
def test_go_terminal_nonerror_states_remain_distinguishable(monkeypatch, project_factory, status):
    config = base_config()
    config["annotation"] = _annotation()
    root = project_factory(config=config)
    l2 = prepare_l2(validate_project(root), run_id=None)
    l2.output_dir.mkdir(parents=True)
    (l2.output_dir / "l2_state.json").write_text(json.dumps({"status": "SUCCESS"}), encoding="utf-8")
    contrast = l2.output_dir / "contrasts" / "Treatment_vs_Control"
    contrast.mkdir(parents=True)
    for name in ("all_genes.tsv", "significant.tsv", "upregulated.tsv", "downregulated.tsv"):
        (contrast / name).write_text("gene_id\nGeneA\n", encoding="utf-8")
    prepared = prepare_go(validate_project(root), run_id=None)
    monkeypatch.setattr("rnaseq.go._check_go_runtime", lambda _organism: None)

    def fake_backend(arguments, **_kwargs):
        output = Path(arguments[-1]).parent
        annotation = output.parent.parent / "annotation"
        annotation.mkdir(exist_ok=True)
        (annotation / "gene_mapping.tsv").write_text("original_gene_id\nGeneA\n", encoding="utf-8")
        summary = {"status": status, "mapping": {"tested": {"source_genes": 1, "mapped_source_genes": 1, "unmapped_source_genes": 0, "one_to_many_source_ids": 0, "duplicate_target_ids": 0, "unique_target_genes": 1, "mapping_rate": 1}}, "contrasts": [], "annotation_database": "org.Mm.eg.db", "annotation_database_version": "test", "clusterProfiler_version": "test"}
        (output / "go_backend_summary.json").write_text(json.dumps(summary), encoding="utf-8")
        return SimpleNamespace(returncode=0)

    monkeypatch.setattr("rnaseq.go.subprocess.run", fake_backend)
    result = execute_go(prepared)
    assert result.status == status
    assert json.loads(result.state_path.read_text(encoding="utf-8"))["status"] == status


def test_go_backend_failure_records_failed_state(monkeypatch, project_factory):
    config = base_config()
    config["annotation"] = _annotation()
    root = project_factory(config=config)
    l2 = prepare_l2(validate_project(root), run_id=None)
    l2.output_dir.mkdir(parents=True)
    (l2.output_dir / "l2_state.json").write_text(json.dumps({"status": "SUCCESS"}), encoding="utf-8")
    contrast = l2.output_dir / "contrasts" / "Treatment_vs_Control"
    contrast.mkdir(parents=True)
    for name in ("all_genes.tsv", "significant.tsv", "upregulated.tsv", "downregulated.tsv"):
        (contrast / name).write_text("gene_id\nGeneA\n", encoding="utf-8")
    prepared = prepare_go(validate_project(root), run_id=None)
    monkeypatch.setattr("rnaseq.go._check_go_runtime", lambda _organism: None)
    monkeypatch.setattr("rnaseq.go.subprocess.run", lambda *_args, **_kwargs: SimpleNamespace(returncode=1))
    with pytest.raises(DownstreamExecutionError, match="GO R backend failed"):
        execute_go(prepared)
    assert json.loads((prepared.output_dir / "go_state.json").read_text(encoding="utf-8"))["status"] == "FAILED"


def test_annotation_package_maps_real_entrez_symbol_and_ensembl_suffix():
    import subprocess

    require_r_packages("AnnotationDbi", "org.Mm.eg.db")
    code = "suppressPackageStartupMessages({library(AnnotationDbi);library(org.Mm.eg.db)}); a<-select(org.Mm.eg.db,keys='11303',keytype='ENTREZID',columns=c('SYMBOL','ENSEMBL')); stopifnot(a$SYMBOL[[1]]=='Abca1', a$ENSEMBL[[1]]=='ENSMUSG00000015243'); stopifnot(strsplit(paste0(a$ENSEMBL[[1]],'.7'),'.',fixed=TRUE)[[1]][1]==a$ENSEMBL[[1]])"
    result = subprocess.run(["Rscript", "-e", code], capture_output=True, text=True, check=False)
    assert result.returncode == 0, result.stderr


def test_go_ora_uses_finite_pvalue_mapped_universe_and_emits_not_applicable_schema(tmp_path):
    """ORA excludes every nonfinite p-value but retains finite p-value/NA padj rows."""

    import subprocess

    require_r_packages("AnnotationDbi", "clusterProfiler", "org.Mm.eg.db")
    all_genes = tmp_path / "all_genes.tsv"
    significant = tmp_path / "significant.tsv"
    empty = tmp_path / "empty.tsv"
    all_genes.write_text(
        "gene_id\tpvalue\tpadj\n11303\t0.01\tNA\n11304\t0.02\t0.04\n11305\tNA\tNA\n11306\tNaN\tNA\n11307\tInf\tNA\n11308\t-Inf\tNA\n",
        encoding="utf-8",
    )
    significant.write_text("gene_id\n11303\n", encoding="utf-8")
    empty.write_text("gene_id\n", encoding="utf-8")
    output = tmp_path / "go"
    config = {
        "output_dir": str(output), "orgdb_package": "org.Mm.eg.db",
        "annotation": {"organism": "Mus musculus", "input_id_type": "ENTREZID", "target_id_type": "ENTREZID", "mapping_warning_rate": 0, "minimum_mapped_foreground": 3, "enrichment": {"go": {"pvalue_cutoff": 0.05, "qvalue_cutoff": 0.2, "p_adjust_method": "BH"}}},
        "contrasts": [{"contrast_id": "fixture", "all_genes": str(all_genes), "significant": str(significant), "up": str(empty), "down": str(empty)}],
        "runtime": {"script": "go_analysis.R", "source_revision": "fixture", "container_image": "fixture", "output_schema": "nf-rna.scientific-provenance.v1"},
    }
    config_path = tmp_path / "config.json"
    config_path.write_text(json.dumps(config), encoding="utf-8")
    script = Path(__file__).parents[1] / "src" / "rnaseq" / "r" / "go_analysis.R"
    result = subprocess.run(["Rscript", str(script), "--config", str(config_path)], capture_output=True, text=True, check=False)
    assert result.returncode == 0, result.stderr
    summary = json.loads((output / "go_backend_summary.json").read_text(encoding="utf-8"))
    universe = summary["contrasts"][0]["universe"]
    assert universe["statistically_tested_genes"] == 2
    assert universe["pvalue_na_excluded_genes"] == 1
    assert universe["pvalue_nan_excluded_genes"] == 1
    assert universe["pvalue_infinite_excluded_genes"] == 2
    assert universe["pvalue_nonfinite_excluded_genes"] == 4
    assert universe["successfully_mapped_tested_genes"] == 2
    for foreground in ("significant", "up", "down"):
        outcome = summary["contrasts"][0]["foregrounds"][foreground]
        assert outcome["status"] == "NOT_APPLICABLE"
        assert outcome["evaluated_term_count"] == outcome["significant_term_count"] == 0
        for ontology in ("BP", "MF", "CC"):
            for name in (f"GO_{ontology}_all_terms.tsv", f"GO_{ontology}.tsv", f"GO_{ontology}_gene_members.tsv"):
                lines = (output / "fixture" / foreground / name).read_text(encoding="utf-8").splitlines()
                assert len(lines) == 1
    assert (output / "scientific_provenance.json").is_file()
    assert (output / "r_session_info.txt").is_file()
    provenance = (output / "scientific_provenance.json").read_text(encoding="utf-8")
    assert str(tmp_path) not in provenance
    assert '"source_revision": "fixture"' in provenance
