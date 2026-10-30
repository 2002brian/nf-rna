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
    assert KeggResourceAdapter("mmu").provider == "kegg_rest_snapshot"
    assert KeggResourceAdapter("mmu").snapshot_root.parts[-4:] == ("nf-rna", "genesets", "kegg", "mmu")


def test_kegg_uses_one_snapshot_with_enricher_and_gsea_and_no_live_query():
    assert KeggGseaConfig().pvalue_cutoff == 0.05  # schema-compatible; not a v1.5 significance input
    script = (Path(__file__).parents[1] / "src" / "rnaseq" / "r" / "kegg_analysis.R").read_text(encoding="utf-8")
    assert "GSEA(geneList=gene_list, TERM2GENE=term2gene, TERM2NAME=term2name" in script
    assert "enricher(gene=targets, universe=universe$genes, TERM2GENE=term2gene, TERM2NAME=term2name" in script
    assert "nPermSimple=NF_RNA_GSEA_NPERM_SIMPLE" in script and "set.seed(as.integer(gsea_cfg$seed))" in script
    assert "enrichKEGG(" not in script and "gseKEGG(" not in script and "use_internal_data" not in script
    assert "calculation_pvalue_cutoff=1" in script


def test_kegg_ora_has_the_same_tested_universe_and_foreground_subset_contract():
    """KEGG is network-dependent at runtime; its deterministic preflight is shared with GO."""

    script = (Path(__file__).parents[1] / "src" / "rnaseq" / "r" / "kegg_analysis.R").read_text(encoding="utf-8")
    go_script = (Path(__file__).parents[1] / "src" / "rnaseq" / "r" / "go_analysis.R").read_text(encoding="utf-8")
    assert '"ora_helpers.R"' in script and '"ora_helpers.R"' in go_script
    assert "gene_sets <- nf_rna_ora_gene_sets(tab)" in script
    assert "gene_sets <- nf_rna_ora_gene_sets(all_table)" in go_script
    assert "is not a subset of the padj-non-NA genes" in script
    assert "pvalueCutoff=1" in script and "qvalueCutoff=1" in script


def test_shared_ora_helper_universe_is_padj_non_na():
    import subprocess

    require_rscript()
    helper = Path(__file__).parents[1] / "src" / "rnaseq" / "r" / "ora_helpers.R"
    code = f'''source("{helper}")
tab <- data.frame(gene_id=c("tested", "filtered", "outlier", "tested"),
                  pvalue=c("0.01", "0.30", NA, "0.02"), padj=c("0.04", NA, NA, "0.05"), stringsAsFactors=FALSE)
sets <- nf_rna_ora_gene_sets(tab)
stopifnot(
  identical(sets$retained, c("tested", "filtered", "outlier")),
  identical(sets$tested, c("tested")),
  sets$counts$finite_pvalue_genes == 2L,
  sets$counts$padj_non_na_genes == 1L,
  sets$counts$padj_na_excluded_genes == 2L
)
'''
    result = subprocess.run(["Rscript", "-e", code], capture_output=True, text=True, check=False)
    assert result.returncode == 0, result.stderr


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
    assert provenance["resource_provider"] == "kegg_rest_snapshot"
    assert "snapshot" in provenance["warning"]


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


def test_kegg_ora_dotplot_is_labelled_from_the_contrast_metadata_and_foreground():
    """KEGG needs the network at runtime; the shared labelling helper is exercised in test_figures and test_go."""

    script = (Path(__file__).parents[1] / "src" / "rnaseq" / "r" / "kegg_analysis.R").read_text(encoding="utf-8")
    assert "run_fg <- function(ids, name, root, mapping, tested, universe, contrast, lookup)" in script
    assert "subtitle=nf_rna_ora_subtitle(fig, contrast, name)" in script
    assert "foreground=name" in script and "numerator=contrast$numerator, denominator=contrast$denominator" in script
    assert "tested_mapping, tested, universe, contrast, lookup)" in script  # the contrast is passed explicitly, not read from the loop scope


def test_kegg_snapshot_resolution_prefers_the_latest_and_rejects_tampered_files(tmp_path):
    import hashlib
    import subprocess

    require_rscript()
    root = tmp_path / "snapshots"
    for snapshot_id, retrieved in (("hsa_old", "2026-01-01T00:00:00Z"), ("hsa_new", "2026-06-01T00:00:00Z")):
        target = root / snapshot_id
        target.mkdir(parents=True)
        (target / "TERM2GENE.tsv").write_text("term\tgene\nhsa00010\t124\n", encoding="utf-8")
        (target / "TERM2NAME.tsv").write_text("term\tname\nhsa00010\tGlycolysis / Gluconeogenesis\n", encoding="utf-8")
        files = {name: {"sha256": hashlib.sha256((target / name).read_bytes()).hexdigest()} for name in ("TERM2GENE.tsv", "TERM2NAME.tsv")}
        (target / "manifest.json").write_text(json.dumps({"snapshot_id": snapshot_id, "retrieved_at_utc": retrieved, "files": files}), encoding="utf-8")
    r_root = Path(__file__).parents[1] / "src" / "rnaseq" / "r"
    code = f'''`%||%` <- function(a, b) if (is.null(a)) b else a
source("{r_root / "kegg_snapshot.R"}")
latest <- nf_rna_resolve_kegg_snapshot("{root}", "hsa")
pinned <- nf_rna_resolve_kegg_snapshot("{root}", "hsa", "hsa_old")
stopifnot(latest$id == "hsa_new", latest$selection == "latest existing snapshot", pinned$id == "hsa_old", nrow(latest$term2gene) == 1L)
cat("x\\n", file="{root / "hsa_old" / "TERM2GENE.tsv"}", append=TRUE)
tampered <- tryCatch(nf_rna_resolve_kegg_snapshot("{root}", "hsa", "hsa_old"), error=function(e) conditionMessage(e))
stopifnot(is.character(tampered), grepl("does not match manifest.json", tampered))
missing <- tryCatch(nf_rna_resolve_kegg_snapshot("{root}", "hsa", "hsa_absent"), error=function(e) conditionMessage(e))
stopifnot(grepl("does not exist", missing), !dir.exists("{root / ".lock"}"))
'''
    result = subprocess.run(["Rscript", "-e", code], capture_output=True, text=True, check=False)
    assert result.returncode == 0, result.stderr
