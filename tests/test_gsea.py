from __future__ import annotations

import json
from pathlib import Path
from types import SimpleNamespace

import pytest

from conftest import base_config
from rnaseq.errors import DownstreamExecutionError
from rnaseq.gsea import execute_gsea, prepare_gsea
from rnaseq.l2 import prepare_l2
from rnaseq.validators import validate_project


def _annotation() -> dict[str, object]:
    return {
        "organism": "Mus musculus", "input_id_type": "ENTREZID", "target_id_type": "ENTREZID",
        "gene_symbol_output": True, "minimum_mapping_rate": 0.70, "minimum_mapped_foreground": 5,
        "enrichment": {"gsea": {"minimum_ranked_genes": 50, "min_gs_size": 10, "max_gs_size": 500, "pvalue_cutoff": 0.05, "padj_cutoff": 0.05, "p_adjust_method": "BH", "seed": 1}},
    }


def test_gsea_requires_annotation_and_successful_l2(project_factory):
    with pytest.raises(DownstreamExecutionError, match="explicit annotation contract"):
        prepare_gsea(validate_project(project_factory()), run_id=None)
    config = base_config()
    config["annotation"] = _annotation()
    with pytest.raises(DownstreamExecutionError, match="existing successful L2"):
        prepare_gsea(validate_project(project_factory(config=config)), run_id=None)


def test_gsea_config_rejects_invalid_gene_set_bounds(project_factory):
    config = base_config()
    config["annotation"] = _annotation()
    config["annotation"]["enrichment"]["gsea"]["min_gs_size"] = 501
    config["annotation"]["enrichment"]["gsea"]["max_gs_size"] = 500
    report = validate_project(project_factory(config=config))
    assert any("max_gs_size must be >=" in issue.message for issue in report.errors)


def test_gsea_blocked_state_is_separate_from_l2(monkeypatch, project_factory):
    config = base_config()
    config["annotation"] = _annotation()
    root = project_factory(config=config)
    l2 = prepare_l2(validate_project(root), run_id=None)
    l2.output_dir.mkdir(parents=True)
    (l2.output_dir / "l2_state.json").write_text(json.dumps({"status": "SUCCESS"}), encoding="utf-8")
    contrast = l2.output_dir / "contrasts" / "Treatment_vs_Control"
    contrast.mkdir(parents=True)
    (contrast / "all_genes.tsv").write_text("gene_id\tstat\n11303\t2.0\n", encoding="utf-8")
    prepared = prepare_gsea(validate_project(root), run_id=None)
    monkeypatch.setattr("rnaseq.gsea._check_go_runtime", lambda _organism: None)

    def fake_backend(arguments, **_kwargs):
        output = Path(arguments[-1]).parent
        (output / "gsea_backend_summary.json").write_text(json.dumps({"status": "BLOCKED", "reason": "minimum ranked genes", "contrasts": [{"contrast_id": "Treatment_vs_Control", "ranking": {"final_ranked_genes": 1, "mapping_rate": 1.0, "positive_stats": 1, "negative_stats": 0, "zero_stats": 0}, "ontologies": {}}], "annotation_database": "org.Mm.eg.db", "annotation_database_version": "test", "clusterProfiler_version": "test"}), encoding="utf-8")
        return SimpleNamespace(returncode=0)

    monkeypatch.setattr("rnaseq.gsea.subprocess.run", fake_backend)
    result = execute_gsea(prepared)
    assert result.status == "BLOCKED"
    assert json.loads((l2.output_dir / "l2_state.json").read_text())["status"] == "SUCCESS"
    assert json.loads(result.state_path.read_text())["status"] == "BLOCKED"


def test_gsea_rejects_unsupported_organism(project_factory):
    config = base_config()
    config["organism"] = {"species": "Rattus norvegicus"}
    with pytest.raises(DownstreamExecutionError, match="not currently supported"):
        prepare_gsea(validate_project(project_factory(config=config)), run_id=None)


def test_gsea_rank_backend_uses_finite_stats_and_deterministic_ties(tmp_path):
    import shutil
    import subprocess

    if shutil.which("Rscript") is None:
        pytest.skip("Rscript unavailable")
    all_genes = tmp_path / "all_genes.tsv"
    all_genes.write_text(
        "gene_id\tstat\n11287\t3\n11298\t3\n11303\t-2\nunknown\t7\n11304\tNA\n",
        encoding="utf-8",
    )
    output = tmp_path / "gsea"
    config = {
        "output_dir": str(output),
        "orgdb_package": "org.Mm.eg.db",
        "annotation": {
            "input_id_type": "ENTREZID", "minimum_mapping_rate": 0.70,
            "enrichment": {"gsea": {"minimum_ranked_genes": 4, "min_gs_size": 10, "max_gs_size": 500, "pvalue_cutoff": 0.05, "padj_cutoff": 0.05, "p_adjust_method": "BH", "seed": 1}},
        },
        "contrasts": [{"contrast_id": "test", "all_genes": str(all_genes)}],
    }
    config_path = tmp_path / "config.json"
    config_path.write_text(json.dumps(config), encoding="utf-8")
    script = Path(__file__).parents[1] / "src" / "rnaseq" / "r" / "gsea_analysis.R"
    result = subprocess.run(["Rscript", str(script), "--config", str(config_path)], capture_output=True, text=True, check=False)
    assert result.returncode == 0, result.stderr
    summary = json.loads((output / "gsea_backend_summary.json").read_text())
    assert summary["status"] == "BLOCKED"
    ranking = summary["contrasts"][0]["ranking"]
    assert ranking["finite_stat_source_genes"] == 4
    assert ranking["nonfinite_or_missing_stat_source_genes"] == 1
    assert ranking["mapped_source_genes"] == 3
    ranked = (output / "test" / "ranked_gene_list.tsv").read_text().splitlines()
    assert ranked[1].split("\t")[1] == "11287"
    assert ranked[2].split("\t")[1] == "11298"


def test_gsea_backend_processes_ontologies_sequentially_with_bounded_lifecycle():
    script = (Path(__file__).parents[1] / "src" / "rnaseq" / "r" / "gsea_analysis.R").read_text()
    helper = (Path(__file__).parents[1] / "src" / "rnaseq" / "r" / "gsea_core_members.R").read_text()
    assert 'for (ontology in c("BP", "MF", "CC"))' in script
    assert 'lapply(c("BP", "MF", "CC"), function(x) run_ontology' not in script
    assert "rm(result, raw, terms, significant, positive, negative, top, p, core_audit)" in script
    assert "gc(verbose=FALSE)" in script
    assert 'write_gsea_core_members' in script
    assert 'append=TRUE' in helper  # core members are streamed, not accumulated in a global rows list
    assert "completed_ontologies" in script and "failed_ontology" in script


def test_gsea_backend_retains_na_pathways_but_excludes_them_from_nes_subsets():
    text = (Path(__file__).parents[1] / "src" / "rnaseq" / "r" / "gsea_analysis.R").read_text()
    assert "na_pathways=sum(is.na(terms$pvalue) | is.na(terms$p.adjust) | is.na(terms$NES))" in text
    assert "!is.na(significant$NES) & significant$NES > 0" in text
    assert "!is.na(significant$NES) & significant$NES < 0" in text


def test_gsea_core_members_have_exact_cardinality_and_no_na_expansion(tmp_path):
    import subprocess

    helper = Path(__file__).parents[1] / "src" / "rnaseq" / "r" / "gsea_core_members.R"
    output = tmp_path / "members.tsv"
    code = f'''source("{helper}")
terms <- data.frame(ID=c("GO:1", "GO:2"), Description=c("first", "second"), setSize=c(3L, 3L), core_enrichment=c("2/1", "2/3"), stringsAsFactors=FALSE)
mapping <- data.frame(mapped_entrez_id=c("1", "2", "3"), original_gene_id=c("gene1", "gene2", "gene3"), mapped_symbol=c("A", "B", "C"), stat=c(1.1, 2.2, -3.3), stringsAsFactors=FALSE)
audit <- write_gsea_core_members(terms, mapping, "{output}", 500L)
tab <- read.delim("{output}", check.names=FALSE, stringsAsFactors=FALSE)
stopifnot(audit$rows == 4L, nrow(tab) == 4L, !anyDuplicated(paste(tab$GO_ID, tab$entrez_id)), !anyNA(tab$original_gene_id), all(is.finite(tab$stat)))
'''
    result = subprocess.run(["Rscript", "-e", code], capture_output=True, text=True, check=False)
    assert result.returncode == 0, result.stderr


def test_gsea_core_member_guard_rejects_nonunique_rank_mapping(tmp_path):
    import subprocess

    helper = Path(__file__).parents[1] / "src" / "rnaseq" / "r" / "gsea_core_members.R"
    output = tmp_path / "members.tsv"
    code = f'''source("{helper}")
terms <- data.frame(ID="GO:1", Description="first", setSize=2L, core_enrichment="1", stringsAsFactors=FALSE)
mapping <- data.frame(mapped_entrez_id=c("1", "1"), original_gene_id=c("gene1", "gene1b"), mapped_symbol=c("A", "A2"), stat=c(1.1, 1.2), stringsAsFactors=FALSE)
err <- try(write_gsea_core_members(terms, mapping, "{output}", 500L), silent=TRUE)
stopifnot(inherits(err, "try-error"))
'''
    result = subprocess.run(["Rscript", "-e", code], capture_output=True, text=True, check=False)
    assert result.returncode == 0, result.stderr


def test_gsea_core_members_ignore_empty_tokens_and_use_stable_pathway_order(tmp_path):
    import subprocess

    helper = Path(__file__).parents[1] / "src" / "rnaseq" / "r" / "gsea_core_members.R"
    output = tmp_path / "members.tsv"
    code = f'''source("{helper}")
terms <- data.frame(ID=c("GO:Z", "GO:A", "GO:EMPTY"), Description=c("z", "a", "empty"), setSize=c(2L, 2L, 2L), core_enrichment=c("2//", "1/", NA), stringsAsFactors=FALSE)
mapping <- data.frame(mapped_entrez_id=c("1", "2"), original_gene_id=c("gene1", "gene2"), mapped_symbol=c("A", "B"), stat=c(1.1, 2.2), stringsAsFactors=FALSE)
audit <- write_gsea_core_members(terms, mapping, "{output}", 500L)
tab <- read.delim("{output}", check.names=FALSE, stringsAsFactors=FALSE)
stopifnot(audit$rows == 2L, identical(as.character(tab$GO_ID), c("GO:A", "GO:Z")), identical(as.character(tab$entrez_id), c("1", "2")))
'''
    result = subprocess.run(["Rscript", "-e", code], capture_output=True, text=True, check=False)
    assert result.returncode == 0, result.stderr
