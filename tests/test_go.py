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
        "annotation": {"organism": "Mus musculus", "input_id_type": "ENTREZID", "target_id_type": "ENTREZID", "mapping_warning_rate": 0, "minimum_mapping_rate": 0, "minimum_mapped_foreground": 3, "enrichment": {"go": {"pvalue_cutoff": 0.05, "qvalue_cutoff": 0.2, "p_adjust_method": "BH"}}},
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


def test_go_ora_dotplots_identify_comparison_and_foreground_from_contrast_metadata(tmp_path):
    """End to end with the real GO annotation: the label comes from numerator/denominator, never the contrast ID."""

    import subprocess

    require_r_packages("AnnotationDbi", "clusterProfiler", "org.Mm.eg.db", "ggplot2", "ggrepel", "pheatmap", "systemfonts", "jsonlite")
    all_genes, significant, empty = tmp_path / "all_genes.tsv", tmp_path / "significant.tsv", tmp_path / "empty.tsv"
    prepare = subprocess.run(["Rscript", "-e", (
        'suppressMessages(library(org.Mm.eg.db)); set.seed(1);'
        'immune <- unique(AnnotationDbi::select(org.Mm.eg.db, keys="GO:0006955", keytype="GOALL", columns="ENTREZID")$ENTREZID);'
        'immune <- head(sort(immune), 150);'
        'background <- head(setdiff(sort(AnnotationDbi::keys(org.Mm.eg.db, "ENTREZID")), immune), 2500);'
        'genes <- c(immune, background);'
        f'write.table(data.frame(gene_id=genes, pvalue=0.5, padj=0.9), "{all_genes}", sep="\t", quote=FALSE, row.names=FALSE);'
        f'write.table(data.frame(gene_id=immune), "{significant}", sep="\t", quote=FALSE, row.names=FALSE)')],
        capture_output=True, text=True, check=False)
    assert prepare.returncode == 0, prepare.stderr
    empty.write_text("gene_id\n", encoding="utf-8")
    output = tmp_path / "go"
    config = {
        "output_dir": str(output), "orgdb_package": "org.Mm.eg.db",
        "annotation": {"organism": "Mus musculus", "input_id_type": "ENTREZID", "target_id_type": "ENTREZID", "mapping_warning_rate": 0, "minimum_mapping_rate": 0, "minimum_mapped_foreground": 5, "enrichment": {"go": {"pvalue_cutoff": 0.05, "qvalue_cutoff": 0.2, "p_adjust_method": "BH"}}},
        "contrasts": [{"contrast_id": "Misleading_vs_Id", "numerator": "133Gy", "denominator": "22Gy", "all_genes": str(all_genes), "significant": str(significant), "up": str(empty), "down": str(empty)}],
        "runtime": {"script": "go_analysis.R", "source_revision": "fixture", "container_image": "fixture", "output_schema": "nf-rna.scientific-provenance.v1"},
    }
    config_path = tmp_path / "config.json"
    config_path.write_text(json.dumps(config), encoding="utf-8")
    script = Path(__file__).parents[1] / "src" / "rnaseq" / "r" / "go_analysis.R"
    result = subprocess.run(["Rscript", str(script), "--config", str(config_path)], capture_output=True, text=True, check=False)
    assert result.returncode == 0, result.stderr
    records = {item["id"]: item for item in json.loads((output / "figure_manifest.json").read_text(encoding="utf-8"))["figures"]}
    assert records, "the immune foreground must produce at least one significant GO ORA dotplot"
    for record in records.values():
        assert record["subtitle"] == "133Gy vs 22Gy \u2014 Significant DE genes"
        assert record["foreground"] == "significant"
        assert record["contrast"] == {"contrast_id": "Misleading_vs_Id", "numerator": "133Gy", "denominator": "22Gy"}
        assert "Misleading" not in record["subtitle"]
        assert {entry["format"] for entry in record["files"]} == {"png", "tiff", "pdf"}
        for entry in record["files"]:
            assert (output / entry["path"]).stat().st_size > 0


def test_go_ora_gene_members_have_no_na_padding_for_unmapped_ensembl_ids(tmp_path):
    """Unmapped tested genes must not add an all-NA member row per term gene.

    Regression: logical indexing on an NA ``mapped_entrez_id`` returned one NA
    row per unmapped gene for every term-gene pair, inflating the member tables
    by orders of magnitude while the term tables stayed correct.
    """

    import csv
    import subprocess

    require_r_packages("AnnotationDbi", "clusterProfiler", "org.Mm.eg.db", "ggplot2", "ggrepel", "pheatmap", "systemfonts", "jsonlite")
    all_genes, significant, empty = tmp_path / "all_genes.tsv", tmp_path / "significant.tsv", tmp_path / "empty.tsv"
    prepare = subprocess.run(["Rscript", "-e", (
        'suppressMessages(library(org.Mm.eg.db));'
        'immune <- unique(AnnotationDbi::select(org.Mm.eg.db, keys="GO:0006955", keytype="GOALL", columns="ENTREZID")$ENTREZID);'
        'immune <- head(sort(immune), 150);'
        'background <- head(setdiff(sort(AnnotationDbi::keys(org.Mm.eg.db, "ENTREZID")), immune), 2500);'
        'ens <- suppressMessages(AnnotationDbi::select(org.Mm.eg.db, keys=c(immune, background), keytype="ENTREZID", columns="ENSEMBL"));'
        'ens <- ens[!is.na(ens$ENSEMBL) & !duplicated(ens$ENTREZID) & !duplicated(ens$ENSEMBL), ];'
        'mapped_fg <- ens$ENSEMBL[ens$ENTREZID %in% immune]; mapped_bg <- ens$ENSEMBL[!(ens$ENTREZID %in% immune)];'
        'unmapped <- sprintf("ENSMUSG9%010d", seq_len(300));'
        'genes <- c(paste0(mapped_fg, ".3"), mapped_bg, unmapped);'
        f'write.table(data.frame(gene_id=genes, pvalue=0.5, padj=0.9), "{all_genes}", sep="\t", quote=FALSE, row.names=FALSE);'
        f'write.table(data.frame(gene_id=c(paste0(mapped_fg, ".3"), head(unmapped, 20))), "{significant}", sep="\t", quote=FALSE, row.names=FALSE)')],
        capture_output=True, text=True, check=False)
    assert prepare.returncode == 0, prepare.stderr
    empty.write_text("gene_id\n", encoding="utf-8")
    output = tmp_path / "go"
    config = {
        "output_dir": str(output), "orgdb_package": "org.Mm.eg.db",
        "annotation": {"organism": "Mus musculus", "input_id_type": "ENSEMBL", "target_id_type": "ENTREZID", "mapping_warning_rate": 0, "minimum_mapping_rate": 0, "minimum_mapped_foreground": 5, "enrichment": {"go": {"pvalue_cutoff": 0.05, "qvalue_cutoff": 0.2, "p_adjust_method": "BH"}}},
        "contrasts": [{"contrast_id": "fixture", "numerator": "T", "denominator": "C", "all_genes": str(all_genes), "significant": str(significant), "up": str(empty), "down": str(empty)}],
        "runtime": {"script": "go_analysis.R", "source_revision": "fixture", "container_image": "fixture", "output_schema": "nf-rna.scientific-provenance.v1"},
    }
    config_path = tmp_path / "config.json"
    config_path.write_text(json.dumps(config), encoding="utf-8")
    script = Path(__file__).parents[1] / "src" / "rnaseq" / "r" / "go_analysis.R"
    result = subprocess.run(["Rscript", str(script), "--config", str(config_path)], capture_output=True, text=True, check=False)
    assert result.returncode == 0, result.stderr

    def rows(path: Path) -> list[dict[str, str]]:
        with path.open(encoding="utf-8", newline="") as handle:
            return list(csv.DictReader(handle, delimiter="\t"))

    mapping = rows(output / "fixture" / "gene_mapping.tsv")
    unmapped = {row["original_gene_id"] for row in mapping if row["mapping_status"] == "UNMAPPED"}
    assert len(unmapped) == 300, "the fixture must exercise unmapped tested genes"
    sources_by_entrez: dict[str, set[str]] = {}
    for row in mapping:
        if row["mapped_entrez_id"] != "NA":
            sources_by_entrez.setdefault(row["mapped_entrez_id"], set()).add(row["original_gene_id"])

    checked = 0
    for ontology in ("BP", "MF", "CC"):
        terms = rows(output / "fixture" / "significant" / f"GO_{ontology}.tsv")
        members = rows(output / "fixture" / "significant" / f"GO_{ontology}_gene_members.tsv")
        assert all(row["original_gene_id"] != "NA" and row["symbol"] != "NA" for row in members)
        assert not {row["original_gene_id"] for row in members} & unmapped
        assert len(members) == len({tuple(row.values()) for row in members}), "member rows must be unique"
        expected = {
            (term["ID"], entrez, source)
            for term in terms for entrez in term["geneID"].split("/") for source in sources_by_entrez[entrez]
        }
        assert {(row["GO_ID"], row["entrez_id"], row["original_gene_id"]) for row in members} == expected
        assert len(members) == len(expected)
        checked += len(members)
    assert checked, "the immune foreground must produce significant GO terms with members"
