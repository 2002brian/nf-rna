from __future__ import annotations

import json
from pathlib import Path
from types import SimpleNamespace

import pytest
import yaml

from conftest import base_config, require_rscript
from rnaseq.errors import DownstreamExecutionError
from rnaseq.kegg import KEGG_CODES, KeggResourceAdapter, execute_kegg, prepare_kegg
from rnaseq.l2 import prepare_l2
from rnaseq.models import KeggGseaConfig
from rnaseq.validators import validate_project


def _annotation() -> dict[str, object]:
    return {
        "organism": "Mus musculus", "input_id_type": "ENTREZID", "target_id_type": "ENTREZID",
        "gene_symbol_output": True, "minimum_mapping_rate": 0.70, "minimum_mapped_foreground": 5,
        "enrichment": {"kegg": {"resource_provider": "online_kegg_rest_via_clusterprofiler", "ora": {"pvalue_cutoff": 0.05, "qvalue_cutoff": 0.2, "p_adjust_method": "BH", "min_gs_size": 10, "max_gs_size": 500}, "gsea": {"minimum_ranked_genes": 50, "pvalue_cutoff": 0.05, "padj_cutoff": 0.05, "p_adjust_method": "BH", "min_gs_size": 10, "max_gs_size": 500, "seed": 1}}},
    }


def _ready(root: Path) -> None:
    l2 = prepare_l2(validate_project(root), run_id=None)
    l2.output_dir.mkdir(parents=True)
    (l2.output_dir / "l2_state.json").write_text(json.dumps({"status": "SUCCESS"}), encoding="utf-8")
    contrast = l2.output_dir / "contrasts" / "Treatment_vs_Control"
    contrast.mkdir(parents=True)
    (contrast / "all_genes.tsv").write_text("gene_id\tstat\n11303\t2.0\n", encoding="utf-8")
    for name in ("significant.tsv", "upregulated.tsv", "downregulated.tsv"):
        (contrast / name).write_text("gene_id\n11303\n", encoding="utf-8")


def test_kegg_organism_codes_are_explicit():
    assert KEGG_CODES == {"Homo sapiens": "hsa", "Mus musculus": "mmu"}
    assert KeggResourceAdapter("mmu").probe_endpoint == "https://rest.kegg.jp/list/pathway/mmu"


def test_kegg_gsea_summary_distinguishes_significance_from_calculation_cutoffs():
    assert KeggGseaConfig().pvalue_cutoff == 0.05
    script = (Path(__file__).parents[1] / "src" / "rnaseq" / "r" / "kegg_analysis.R").read_text(encoding="utf-8")
    assert "gseKEGG(geneList=gene_list" in script and "pvalueCutoff=1" in script
    assert "configured_pvalue_cutoff=as.numeric(cfg$annotation$enrichment$kegg$gsea$pvalue_cutoff)" in script
    assert "configured_padj_cutoff=as.numeric(cfg$annotation$enrichment$kegg$gsea$padj_cutoff)" in script
    assert "calculation_pvalue_cutoff=1" in script


def test_kegg_requires_annotation_and_rejects_unsupported_organism(project_factory):
    with pytest.raises(DownstreamExecutionError, match="explicit annotation contract"):
        prepare_kegg(validate_project(project_factory()), run_id=None, mode="ora")
    config = base_config()
    config["organism"] = {"species": "Rattus norvegicus"}
    with pytest.raises(DownstreamExecutionError, match="not currently supported"):
        prepare_kegg(validate_project(project_factory(config=config)), run_id=None, mode="gsea")


@pytest.mark.parametrize("mode,summary_name,state_name", [("ora", "kegg_backend_summary.json", "kegg_state.json"), ("gsea", "gsea_kegg_backend_summary.json", "gsea_kegg_state.json")])
def test_kegg_network_unavailable_isolated_from_l2_and_go(monkeypatch, project_factory, mode, summary_name, state_name):
    config = base_config()
    config["annotation"] = _annotation()
    root = project_factory(config=config)
    _ready(root)
    l2 = prepare_l2(validate_project(root), run_id=None)
    go_state = l2.output_dir / "enrichment" / "go" / "go_state.json"
    go_state.parent.mkdir(parents=True)
    go_state.write_text(json.dumps({"status": "SUCCESS"}), encoding="utf-8")
    prepared = prepare_kegg(validate_project(root), run_id=None, mode=mode)
    monkeypatch.setattr("rnaseq.kegg._check_kegg_runtime", lambda: None)

    def fake_backend(arguments, **_kwargs):
        output = Path(arguments[-1]).parent
        (output / summary_name).write_text(json.dumps({"status": "NETWORK_UNAVAILABLE", "reason": "DNS failed", "resource": {"status": "NETWORK_UNAVAILABLE", "provider": "online_kegg_rest_via_clusterprofiler", "error": "DNS failed"}, "clusterProfiler_version": "test", "r_version": "test", "contrasts": []}), encoding="utf-8")
        return SimpleNamespace(returncode=0)

    monkeypatch.setattr("rnaseq.kegg.subprocess.run", fake_backend)
    result = execute_kegg(prepared)
    assert result.status == "NETWORK_UNAVAILABLE"
    assert json.loads((l2.output_dir / "l2_state.json").read_text())["status"] == "SUCCESS"
    assert json.loads(go_state.read_text())["status"] == "SUCCESS"
    assert json.loads((prepared.output_dir / state_name).read_text())["status"] == "NETWORK_UNAVAILABLE"
    provenance = yaml.safe_load((prepared.output_dir / ("kegg_provenance.yaml" if mode == "ora" else "gsea_kegg_provenance.yaml")).read_text())
    assert provenance["resource_provider"] == "online_kegg_rest_via_clusterprofiler"
    assert "may vary" in provenance["warning"]


def test_kegg_config_rejects_invalid_gene_set_bounds(project_factory):
    config = base_config()
    config["annotation"] = _annotation()
    config["annotation"]["enrichment"]["kegg"]["gsea"]["min_gs_size"] = 501
    report = validate_project(project_factory(config=config))
    assert any("kegg.gsea.max_gs_size" in issue.message for issue in report.errors)


def test_kegg_gsea_core_members_stream_deterministically_without_cartesian_expansion(tmp_path):
    import subprocess

    require_rscript()
    helper = Path(__file__).parents[1] / "src" / "rnaseq" / "r" / "kegg_core_members.R"
    output = tmp_path / "members.tsv"
    code = f'''source("{helper}")
terms <- data.frame(
  ID=c("path:B", "path:A", "path:EMPTY"),
  Description=c("B", "A", "empty"),
  setSize=c(3L, 3L, 3L),
  core_enrichment=c("2/1/2", "3//4", NA),
  stringsAsFactors=FALSE
)
mapping <- data.frame(
  mapped_entrez_id=c("1", "2", "2", "3", "4"),
  mapped_symbol=c("S1", "S2a", "S2b", "S3", "S4"),
  original_gene_id=c("gene1", "gene2a", "gene2b", "gene_shared", "gene_shared"),
  stat=c(1.0, 2.0, 2.0, -3.0, 4.0),
  stringsAsFactors=FALSE
)
audit <- write_kegg_member_table(terms, mapping, "core_enrichment", "{output}", strict=TRUE, max_members=3L)
tab <- read.delim("{output}", check.names=FALSE, stringsAsFactors=FALSE)
stopifnot(
  identical(names(tab), c("pathway_id", "pathway_description", "entrez_id", "symbol", "original_gene_id", "stat")),
  audit$pathway_entrez_pairs == 4L,
  audit$rows == 5L,
  nrow(tab) == 5L,
  identical(as.character(tab$pathway_id), c("path:A", "path:A", "path:B", "path:B", "path:B")),
  identical(as.character(tab$entrez_id), c("3", "4", "1", "2", "2")),
  length(unique(paste(tab$pathway_id, tab$entrez_id, sep="\\r"))) == 4L,
  identical(as.character(tab$original_gene_id[tab$entrez_id %in% c("3", "4")]), c("gene_shared", "gene_shared")),
  !anyNA(tab$original_gene_id),
  all(is.finite(tab$stat))
)
'''
    result = subprocess.run(["Rscript", "-e", code], capture_output=True, text=True, check=False)
    assert result.returncode == 0, result.stderr


def test_kegg_gsea_core_member_writer_keeps_empty_outputs_valid_and_fails_missing_ranked_ids(tmp_path):
    import subprocess

    require_rscript()
    helper = Path(__file__).parents[1] / "src" / "rnaseq" / "r" / "kegg_core_members.R"
    empty_output = tmp_path / "empty.tsv"
    missing_output = tmp_path / "missing.tsv"
    code = f'''source("{helper}")
empty_terms <- data.frame(ID=character(), Description=character(), setSize=integer(), core_enrichment=character(), stringsAsFactors=FALSE)
mapping <- data.frame(mapped_entrez_id="1", mapped_symbol="S1", original_gene_id="gene1", stat=1.0, stringsAsFactors=FALSE)
audit <- write_kegg_member_table(empty_terms, mapping, "core_enrichment", "{empty_output}", strict=TRUE, max_members=10L)
empty_tab <- read.delim("{empty_output}", check.names=FALSE, stringsAsFactors=FALSE)
terms <- data.frame(ID="path:MISSING", Description="missing", setSize=1L, core_enrichment="999", stringsAsFactors=FALSE)
err <- try(write_kegg_member_table(terms, mapping, "core_enrichment", "{missing_output}", strict=TRUE, max_members=10L), silent=TRUE)
stopifnot(audit$rows == 0L, nrow(empty_tab) == 0L, identical(names(empty_tab), c("pathway_id", "pathway_description", "entrez_id", "symbol", "original_gene_id", "stat")), inherits(err, "try-error"), grepl("mapped universe is missing Entrez ID 999", as.character(err), fixed=TRUE))
'''
    result = subprocess.run(["Rscript", "-e", code], capture_output=True, text=True, check=False)
    assert result.returncode == 0, result.stderr


def test_kegg_gsea_core_member_writer_is_indexed_and_streamed():
    root = Path(__file__).parents[1] / "src" / "rnaseq" / "r"
    helper = (root / "kegg_core_members.R").read_text()
    script = (root / "kegg_analysis.R").read_text()
    assert "split(which(valid), mapped_ids[valid]" in helper
    assert "append=TRUE" in helper
    assert "rows <- list()" not in helper
    assert "do.call(rbind, rows)" not in helper
    assert 'write_kegg_member_table(terms,mapping,"core_enrichment"' in script
