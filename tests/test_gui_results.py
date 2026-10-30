"""GUI results, reports and downloads (v1.6.0 Phase 6b) on a small synthetic run.

The fixture is a run directory with the files the R pipeline writes, cut down
to a few genes and one enrichment module; no R runtime is needed. The service
tests need no GUI extra; the page tests use AppTest and skip without Streamlit.
"""

from __future__ import annotations

import hashlib
import json
import zipfile
from pathlib import Path

import pytest
import yaml

from rnaseq.gui import results
from rnaseq.service import write_delivery_manifest

# The smallest valid PNG (1x1); figures are only displayed, never decoded by nf-rna.
PNG = bytes.fromhex("89504e470d0a1a0a0000000d4948445200000001000000010806000000" "1f15c4890000000d49444154789c6360000002000154a24f5d0000000049454e44ae426082")
ALL_GENES = (
    "gene_id\tsymbol\tdescription\tbaseMean\tlog2FoldChange\tlfcSE\tstat\tpvalue\tpadj\tlog2FoldChange_shrunk\n"
    "ENSG01\tDUSP1\tdual specificity phosphatase 1\t3409.03835490143\t2.94784960588663\t0.12\t24.18\t3.03342369782533e-129\t1.82328987063955e-125\t2.93659125491481\n"
    "ENSG02\tDUSP10\tdual specificity phosphatase 10\t512.1\t-1.18489661324826\t0.1\t-11.2\t1e-18\t3.45210226969281e-17\t-1.15865770249413\n"
    "ENSG03\tACTB\tactin beta\t9000.5\t0.01\t0.05\t0.2\t0.84\tNA\t0.005\n"
)


def _write(path: Path, text: str | bytes) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    (path.write_bytes if isinstance(text, bytes) else path.write_text)(text)
    return path


def _run(project: Path, level: str = "L2", case: str = "CASE-1", run_id: str = "20261030-120000+0800") -> Path:
    run = project / "runs" / case / run_id
    _write(run / "run_state.json", json.dumps({"case_id": case, "run_id": run_id, "status": "SUCCESS", "started_at": "2026-10-30T12:00:00+08:00"}))
    _write(run / "frozen" / "downstream_contract.json", json.dumps({"analysis_level": level, "analysis": {"enrichment": ["kegg", "go"] if level == "L2" else []}}))
    _write(run / "frozen" / "project.yaml", yaml.safe_dump({"design": {"variables": {"condition": "categorical"}}, "thresholds": {"padj": 0.05, "abs_log2fc": 1.0},
                                                             "annotation": {"enrichment": {"kegg": {"snapshot": "hsa_kegg20261028_X"}}}, "input": {"type": "raw_counts"}}))
    _write(run / "frozen" / "contrasts.csv", "contrast_id,factor,numerator,denominator\nB_vs_A,condition,B,A\n")
    l1 = run / "downstream" / "l1"
    _write(l1 / "backend_summary.json", json.dumps({"samples": ["S1", "S2"], "genes_input": 3, "genes_retained": 3, "package_versions": {"R": "R 4.5", "DESeq2": "1.50.2"}}))
    _write(l1 / "pca.png", PNG)
    _write(l1 / "pca.pdf", b"%PDF-1.4 fixture")
    _write(l1 / "library_size_qc.tsv", "sample_id\tinput_total\nS1\t100\nS2\t120\n")
    _write(l1 / "normalized_counts.tsv", "gene_id\tS1\tS2\nENSG01\t1.5\t2.5\n")
    if level == "L2":
        contrast = run / "downstream" / "l2" / "contrasts" / "B_vs_A"
        _write(contrast / "all_genes.tsv", ALL_GENES)
        header = ALL_GENES.splitlines()[0]
        rows = ALL_GENES.splitlines()[1:3]
        _write(contrast / "significant.tsv", "\n".join([header, *rows]) + "\n")
        _write(contrast / "upregulated.tsv", "\n".join([header, rows[0]]) + "\n")
        _write(contrast / "downregulated.tsv", "\n".join([header, rows[1]]) + "\n")
        _write(contrast / "volcano.png", PNG)
        _write(contrast / "contrast_metadata.json", json.dumps({
            "deg_definition": {"padj_below": 0.05, "abs_log2fc_at_least": 1, "rule": "padj < threshold AND abs(raw log2FoldChange) >= threshold"},
            "shrinkage": {"method": "ashr", "status": "failed", "error": "apeglm did not converge"},
        }))
        kegg = run / "downstream" / "l2" / "enrichment" / "kegg"
        _write(kegg / "kegg_backend_summary.json", json.dumps({
            "status": "SUCCESS", "significance_semantics": "v1.5",
            "significance": {"rule": "significant = p.adjust (BH) <= 0.05", "clusterprofiler_call": {"pAdjustMethod": "BH"}},
            "ora": {"universe": "padj-non-NA genes mapped to Entrez, intersected with the KEGG snapshot", "min_gs_size": 10, "max_gs_size": 500},
            "gene_set_source": {"snapshot_id": "hsa_kegg20261028_X"},
            "contrasts": [{"contrast_id": "B_vs_A", "foregrounds": {
                "significant": {"status": "SUCCESS", "tested_terms": 2, "significant_terms": 1, "universe": {"universe_genes": 6057}, "figures": {"dotplot": "drawn"}},
                "up": {"status": "NO_SIGNIFICANT_TERMS", "reason": "no pathway with p.adjust <= 0.05", "figures": {"dotplot": "no significant pathway: dotplot not drawn"}},
            }}],
        }))
        significant = kegg / "B_vs_A" / "significant"
        terms = "ID\tDescription\tpvalue\tp.adjust\tsignificant\nhsa04010\tMAPK signaling pathway\t1e-6\t0.001\tTRUE\nhsa00010\tGlycolysis\t0.4\t0.9\tFALSE\n"
        _write(significant / "kegg_all_terms.tsv", terms)
        _write(significant / "kegg.tsv", "\n".join(terms.splitlines()[:2]) + "\n")
        _write(significant / "kegg_members.tsv", "pathway_id\tpathway_description\tgene_id\tsymbol\tlog2FoldChange\tpadj\n"
               "hsa04010\tMAPK signaling pathway\tENSG01\tDUSP1\t2.94784960588663\t1.82328987063955e-125\n"
               "hsa04010\tMAPK signaling pathway\tENSG02\tDUSP10\t-1.18489661324826\t3.45210226969281e-17\n")
        _write(significant / "dotplot.png", PNG)
        # A pre-v1.5 GO ORA summary: no significance semantics recorded.
        _write(run / "downstream" / "l2" / "enrichment" / "go" / "go_backend_summary.json", json.dumps({"status": "SUCCESS", "contrasts": []}))
    _write(run / "downstream" / "report" / "report.html", "<!doctype html><html><body><h1>report</h1></body></html>")
    delivery = run / "delivery"
    _write(delivery / "report_20261030.html", (run / "downstream" / "report" / "report.html").read_text())
    _write(delivery / "tables" / "l1" / "library_size_qc.tsv", (l1 / "library_size_qc.tsv").read_text())
    write_delivery_manifest(delivery)
    return run


@pytest.fixture
def results_home(tmp_path, monkeypatch):
    data = tmp_path / "data"
    project = data / "projects" / "demo"
    _write(project / "project.yaml", "project: {id: demo}\n")
    settings = tmp_path / "gui.yaml"
    settings.write_text(yaml.safe_dump({"allowed_roots": [str(data)], "project_parent": str(data / "projects"), "condarc": None, "launcher": "process"}), encoding="utf-8")
    monkeypatch.setenv("NF_RNA_GUI_CONFIG", str(settings))
    monkeypatch.setenv("XDG_CACHE_HOME", str(tmp_path / "cache"))
    _run(project)
    return project


def test_deg_table_search_shows_the_strings_of_all_genes_tsv(results_home):
    run = results.result_runs(results_home)[0]
    view = results.deg_view(run, "B_vs_A")
    genes = results.read_table(view.root / "all_genes.tsv")
    found = results.filter_deg(genes, results.significant_ids(view), query="dusp1")
    assert found.column("symbol") == ["DUSP1", "DUSP10"]
    row = found.rows[0]
    assert (row["log2FoldChange"], row["padj"], row["log2FoldChange_shrunk"]) == ("2.94784960588663", "1.82328987063955e-125", "2.93659125491481")
    # "Only significant" is membership of significant.tsv, not a recomputed threshold.
    assert results.filter_deg(genes, results.significant_ids(view), only_significant=True).column("gene_id") == ["ENSG01", "ENSG02"]
    ordered = results.filter_deg(genes, set(), sort_by="log2FoldChange", descending=True)
    assert ordered.column("gene_id") == ["ENSG01", "ENSG03", "ENSG02"]
    assert results.filter_deg(genes, set(), padj_max=0.01).column("gene_id") == ["ENSG01", "ENSG02"]  # NA padj is not <= 0.01
    assert results.filter_deg(genes, set(), sort_by="padj").column("gene_id") == ["ENSG01", "ENSG02", "ENSG03"]  # NA sorts last
    assert view.counts == {"significant": 2, "upregulated": 1, "downregulated": 1}
    assert (view.padj_below, view.abs_log2fc_at_least) == (0.05, 1)


def test_missing_outputs_are_unavailable_with_their_recorded_reason(results_home):
    run = results.result_runs(results_home)[0]
    figures = {item.label: item for item in results.deg_view(run, "B_vs_A").figures}
    assert figures["Volcano plot: B_vs_A"].png is not None
    assert figures["MA plot: B_vs_A"].png is None and "before nf-rna v1.5" in figures["MA plot: B_vs_A"].reason
    modules = {module.key: module for module in results.enrichment_modules(run)}
    kegg = modules["kegg"]
    assert not kegg.legacy and kegg.facts["KEGG snapshot pinned"] == "yes (hsa_kegg20261028_X)"
    up = next(item for item in results.enrichment_results(kegg, "B_vs_A") if item.label.startswith("upregulated"))
    assert up.status == "NO_SIGNIFICANT_TERMS" and up.figures[0].png is None
    assert up.figures[0].reason == "no significant pathway: dotplot not drawn"
    assert modules["go"].legacy and "legacy" in modules["go"].facts["significance semantics"]


def test_pathway_members_match_the_members_tsv(results_home):
    run = results.result_runs(results_home)[0]
    kegg = next(module for module in results.enrichment_modules(run) if module.key == "kegg")
    significant = results.enrichment_results(kegg, "B_vs_A")[0]
    members = results.pathway_members(significant.members, "hsa04010")
    source = results.read_table(significant.members)
    assert members.rows == tuple(row for row in source.rows if row["pathway_id"] == "hsa04010") and len(members.rows) == 2


def test_an_l1_run_offers_no_l2_results_or_downloads(results_home):
    _run(results_home, level="L1", case="L1-CASE")
    run = next(item for item in results.result_runs(results_home) if item.case_id == "L1-CASE")
    assert not run.is_l2
    groups = {item.group for item in results.downloads(run)}
    assert groups == {"Report", "Counts", "QC tables", "Figures (PNG / PDF)"}
    assert all("l2" not in str(item.path.relative_to(run.run_dir)) for item in results.downloads(run))


def test_delivery_zip_holds_the_verified_delivery_files(results_home):
    run = results.result_runs(results_home)[0]
    archive = results.build_delivery_archive(run)
    with zipfile.ZipFile(archive) as handle:
        for name in handle.namelist():
            relative = name.split("/", 1)[1]
            assert hashlib.sha256(handle.read(name)).hexdigest() == hashlib.sha256((run.delivery / relative).read_bytes()).hexdigest()
    assert archive.parent == Path(results.gui_state_dir()) / "packages"  # never written into the run
    (run.delivery / "tables" / "l1" / "library_size_qc.tsv").write_text("tampered\n")
    archive.unlink()
    with pytest.raises(OSError, match="does not match its manifest"):
        results.build_delivery_archive(run)
    assert results.scp_command(run.delivery).startswith("scp -r ")


# --------------------------------------------------------------------- pages (AppTest)

def _page(name: str, project: Path):
    pytest.importorskip("streamlit")
    from streamlit.testing.v1 import AppTest

    page = AppTest.from_string(f"from rnaseq.gui import pages\npages.{name}()\n", default_timeout=60)
    page.session_state["project"] = str(project)
    return page


def test_results_page_shows_deg_summary_shrinkage_failure_and_search(results_home):
    page = _page("results_page", results_home).run()
    assert not page.exception
    page.radio(key="results_section::L2").set_value("Differential expression").run()
    assert not page.exception, page.exception
    assert any("Shrinkage **failed**" in item.value for item in page.error)
    assert any("MA plot: B_vs_A: unavailable" in item.value for item in page.warning)
    key = "deg::CASE-1/20261030-120000+0800::B_vs_A"
    page.text_input(key=f"{key}::query").set_value("DUSP1").run()
    table = next(frame.value for frame in page.dataframe if "log2FoldChange_shrunk" in frame.value)
    assert table.iloc[0]["padj"] == "1.82328987063955e-125"
    page.radio(key="results_section::L2").set_value("Enrichment").run()
    assert not page.exception
    assert any("Legacy results" in item.value for item in page.warning)  # GO ORA of the fixture is pre-v1.5
    page.radio(key="enr::CASE-1/20261030-120000+0800::module").set_value("KEGG ORA").run()
    assert not page.exception and not any("Legacy results" in item.value for item in page.warning)
    assert any("hsa_kegg20261028_X" in frame.value.to_string() for frame in page.dataframe)


def test_results_page_of_an_l1_run_has_no_l2_sections(results_home):
    _run(results_home, level="L1", case="L1-CASE", run_id="20261031-120000+0800")
    page = _page("results_page", results_home).run()
    assert not page.exception
    assert page.selectbox(key="results_run").value.startswith("L1-CASE")
    assert page.radio(key="results_section::L1").options == ["QC & exploratory", "Provenance"]


def test_reports_page_offers_report_files_and_delivery(results_home):
    page = _page("reports_page", results_home).run()
    assert not page.exception
    labels = [button.label for button in page.get("download_button")]
    assert "Download report_20261030.html" in labels
    page.button(key="rep::zip").click().run()
    assert not page.exception and not page.error
    assert results.delivery_archive_path(results.result_runs(results_home)[0]).is_file()
