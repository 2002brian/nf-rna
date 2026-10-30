"""Drive the nf-rna GUI Results and Reports pages in headless Chrome (v1.6.0 Phase 6b acceptance).

Usage (the GUI must already be running):

    python scripts/dev/gui_results_acceptance.py STEP [STEP ...] --shots DIR [--project DIR] [--base URL]

Steps:
  browse      QC, DEG search (DUSP1 compared with all_genes.tsv), KEGG pathway members
              (compared with kegg_members.tsv), figures loaded, provenance
  downloads   download one file of every kind and compare its SHA-256 with the server file
  report      open the delivered report.html offline and check that every image renders
  unavailable screenshots of a project whose run lacks outputs, and of an L1-only run

Reuses the helpers of gui_browser_acceptance.py; needs the ``gui-dev`` extra.
This is a development tool, not part of the test suite.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
import zipfile
from pathlib import Path

from playwright.sync_api import Page, sync_playwright

sys.path.insert(0, str(Path(__file__).parent))
import gui_browser_acceptance as gui  # noqa: E402

from rnaseq.gui import results  # noqa: E402

GENE = "DUSP1"
EVIDENCE: dict[str, object] = {}


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def grid_rows(page: Page, index: int) -> list[dict[str, str]]:
    """The rendered rows of a data grid, from its accessibility table (only the rows on screen)."""

    grid = page.locator("[data-testid=stDataFrame]").nth(index)
    headers = grid.locator("thead th").all_inner_texts()
    rows = []
    for row in grid.locator("tbody tr").all():
        rows.append(dict(zip(headers, row.locator("td").all_inner_texts())))
    return rows


def images_loaded(page: Page) -> list[tuple[str, int]]:
    """(caption, naturalWidth) of each figure on the page; naturalWidth 0 means not rendered."""

    return page.evaluate(
        "() => Array.from(document.querySelectorAll('[data-testid=stMain] [data-testid=stImage]')).map(f => {"
        " const i = f.querySelector('img'); const c = f.querySelector('[data-testid=stImageCaption]');"
        " return [c ? c.textContent : '', i && i.complete ? i.naturalWidth : 0]; })"
    )


def section(page: Page, name: str) -> None:
    page.get_by_text(name, exact=True).first.click()
    gui.settle(page, 1500)


def step_browse(page: Page, shots: Path, project: Path, prefix: str) -> None:
    run = results.result_runs(project)[0]
    gui.goto(page, "results", str(project))
    gui.shot(page, shots, f"{prefix}1_results_qc")
    assert all(width > 0 for _alt, width in images_loaded(page)), images_loaded(page)

    section(page, "Differential expression")
    gui.fill(page, "Search symbol or gene_id", GENE)
    gui.shot(page, shots, f"{prefix}2_deg_search_{GENE}")
    shown = next(row for row in grid_rows(page, 0) if row.get("symbol") == GENE)
    source = results.read_table(run.l2 / "contrasts" / "trt_vs_untrt" / "all_genes.tsv")
    expected = source.where("symbol", GENE).rows[0]
    compared = {column: (shown[column], expected[column]) for column in ("gene_id", "log2FoldChange", "padj", "log2FoldChange_shrunk")}
    assert all(a == b for a, b in compared.values()), compared
    EVIDENCE[f"{prefix}_deg_{GENE}"] = {column: value[0] for column, value in compared.items()}
    loaded = images_loaded(page)
    assert len(loaded) >= 3 and all(width > 0 for _alt, width in loaded), loaded
    EVIDENCE[f"{prefix}_deg_images"] = loaded
    body = gui.text_of(page)
    assert "Shrinkage method: apeglm" in body and "status success" in body, body[:2000]

    section(page, "Enrichment")
    page.get_by_text("KEGG ORA", exact=True).click()
    gui.settle(page, 1500)
    grids = page.locator("[data-testid=stDataFrame]")
    terms_index = grids.count() - 1  # after the module and result-set facts tables
    grids.nth(terms_index).click(position={"x": 12, "y": 35 + 17})  # row-selection column of the first term
    gui.settle(page, 1500)
    gui.shot(page, shots, f"{prefix}3_kegg_ora_members")
    pathway = grid_rows(page, terms_index)[0]["ID"]
    members = grid_rows(page, terms_index + 1)
    table = results.read_table(run.l2 / "enrichment" / "kegg" / "trt_vs_untrt" / "significant" / "kegg_members.tsv")
    expected_members = table.where("pathway_id", pathway).rows
    assert members and all(row["pathway_id"] == pathway for row in members), members[:2]
    shown_ids = [row["gene_id"] for row in members]
    assert shown_ids == [row["gene_id"] for row in expected_members], (shown_ids, [row["gene_id"] for row in expected_members])
    for row, source_row in zip(members, expected_members):
        assert all(row[column] == source_row[column] for column in row), (row, source_row)
    EVIDENCE[f"{prefix}_kegg_members"] = {"pathway": pathway, "rows_in_tsv": len(expected_members), "rows_compared_on_screen": len(members)}

    page.get_by_text("GO ORA", exact=True).click()
    gui.settle(page, 1500)
    loaded = images_loaded(page)
    assert any("simplified" in caption for caption, width in loaded if width > 0), loaded
    gui.shot(page, shots, f"{prefix}4_go_ora_simplified_dotplot")
    page.get_by_text("GO GSEA", exact=True).click()
    gui.settle(page, 2500)
    loaded = images_loaded(page)
    running = [caption for caption, width in loaded if "running score" in caption and width > 0]
    assert len(running) >= 2, loaded
    EVIDENCE[f"{prefix}_gsea_running_score_rendered"] = running
    gui.shot(page, shots, f"{prefix}5_go_gsea_running_score")
    page.get_by_text("KEGG GSEA", exact=True).click()
    gui.settle(page, 2500)
    gui.shot(page, shots, f"{prefix}6_kegg_gsea")
    section(page, "Provenance")
    gui.shot(page, shots, f"{prefix}7_provenance")


def _download(page: Page, button, target_dir: Path) -> Path:
    with page.expect_download(timeout=300_000) as info:
        button.click()
    download = info.value
    path = target_dir / download.suggested_filename
    download.save_as(path)
    return path


def step_downloads(page: Page, shots: Path, project: Path) -> None:
    run = results.result_runs(project)[0]
    out = shots.parent / "downloads"
    out.mkdir(parents=True, exist_ok=True)
    gui.goto(page, "reports", str(project))
    gui.shot(page, shots, "G1_reports_page")
    checks = []
    report = results.report_path(run)
    checks.append(("report.html", report, _download(page, page.get_by_role("button", name=f"Download {report.name}"), out)))
    items = results.downloads(run)
    groups = list(dict.fromkeys(item.group for item in items if item.group != "Report"))
    wanted = {
        "Counts": ["l2/normalized_counts.tsv"],
        "DEG tables": ["trt_vs_untrt/all_genes.tsv", "trt_vs_untrt/significant.tsv"],
        "Enrichment tables and members": ["kegg/trt_vs_untrt/significant/kegg_all_terms.tsv", "kegg/trt_vs_untrt/significant/kegg_members.tsv"],
        "Figures (PNG / PDF)": ["l2/contrasts/trt_vs_untrt/ma_plot.png", "l2/contrasts/trt_vs_untrt/ma_plot.pdf"],
    }
    for group, labels in wanted.items():
        members = [item for item in items if item.group == group]
        for label in labels:
            gui.choose(page, f"{group} ({len(members)})", label)
            button = page.locator("[data-testid=stDownloadButton]").filter(has_text="Download").nth(1 + groups.index(group))
            item = next(item for item in members if item.label == label)
            checks.append((label, item.path, _download(page, button, out)))
    gui.click(page, "Prepare delivery package (.zip)")
    archive = results.delivery_archive_path(run)
    zipped = _download(page, page.get_by_role("button", name=f"Download {archive.name}"), out)
    checks.append(("delivery zip", archive, zipped))
    gui.shot(page, shots, "G2_delivery_package_ready")
    rows = []
    for label, server, local in checks:
        rows.append({"file": label, "server": str(server), "sha256_server": sha256(server), "sha256_download": sha256(local)})
        assert rows[-1]["sha256_server"] == rows[-1]["sha256_download"], rows[-1]
    with zipfile.ZipFile(zipped) as handle:
        names = handle.namelist()
        mismatched = [name for name in names if hashlib.sha256(handle.read(name)).hexdigest() != sha256(run.delivery / name.split("/", 1)[1])]
    assert not mismatched, mismatched
    EVIDENCE["downloads"] = rows
    EVIDENCE["delivery_zip"] = {"members": len(names), "members_matching_delivery_sha256": len(names) - len(mismatched)}


def step_report(browser, shots: Path, project: Path) -> None:
    report = results.report_path(results.result_runs(project)[0])
    context = browser.new_context(offline=True, viewport={"width": 1400, "height": 1000})
    requested: list[str] = []
    context.on("request", lambda request: requested.append(request.url))
    page = context.new_page()
    page.goto(report.as_uri())
    page.wait_for_load_state("load")
    images = page.evaluate("() => Array.from(document.images).map(i => [i.alt, i.complete && i.naturalWidth > 0])")
    page.screenshot(path=str(shots / "H1_report_offline.png"))
    external = [url for url in requested if not url.startswith(("file:", "data:"))]
    assert images and all(ok for _alt, ok in images) and not external, (images, external)
    headings = page.evaluate("() => Array.from(document.querySelectorAll('h2')).map(h => h.textContent)")
    EVIDENCE["report_offline"] = {"report": str(report), "images": len(images), "all_rendered": True, "external_requests": external, "h2": headings}
    context.close()


def step_unavailable(page: Page, shots: Path, missing: Path, l1_only: Path) -> None:
    gui.goto(page, "results", str(missing))
    section(page, "Differential expression")
    body = gui.text_of(page)
    assert "MA plot: trt_vs_untrt: unavailable" in body and "Shrinkage: not recorded" in body, body[:3000]
    gui.shot(page, shots, "I1_missing_deg_outputs")
    section(page, "Enrichment")
    page.get_by_text("KEGG ORA", exact=True).click()
    gui.settle(page, 1500)
    body = gui.text_of(page)
    assert "KEGG dotplot (all significant DEGs): unavailable" in body, body[:3000]
    gui.shot(page, shots, "I2_missing_enrichment_figure")
    page.get_by_text("GO GSEA", exact=True).click()
    gui.settle(page, 1500)
    body = gui.text_of(page)
    assert "GO GSEA: unavailable" in body, body[:3000]
    gui.shot(page, shots, "I3_missing_enrichment_module")
    gui.goto(page, "results", str(l1_only))
    body = gui.text_of(page)
    assert "Differential expression" not in body and "Enrichment" not in body, body[:3000]
    gui.shot(page, shots, "I4_l1_only_results")
    gui.goto(page, "reports", str(l1_only))
    body = gui.text_of(page)
    assert "DEG tables" not in body and "Enrichment tables" not in body, body[:3000]
    gui.shot(page, shots, "I5_l1_only_downloads")
    EVIDENCE["unavailable_and_l1_only"] = "passed"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("steps", nargs="+")
    parser.add_argument("--shots", type=Path, required=True)
    parser.add_argument("--project", type=Path, default=Path(gui.PROJECT))
    parser.add_argument("--prefix", default="F")
    parser.add_argument("--base", default=gui.BASE)
    parser.add_argument("--missing", type=Path)
    parser.add_argument("--l1-only", type=Path)
    parser.add_argument("--evidence", type=Path)
    args = parser.parse_args()
    gui.BASE = args.base
    args.shots.mkdir(parents=True, exist_ok=True)
    with sync_playwright() as playwright:
        browser = playwright.chromium.launch(executable_path=gui.CHROME, headless=True)
        context = browser.new_context(viewport={"width": 1600, "height": 2200}, accept_downloads=True)
        page = context.new_page()
        for step in args.steps:
            if step == "browse":
                step_browse(page, args.shots, args.project, args.prefix)
            elif step == "downloads":
                step_downloads(page, args.shots, args.project)
            elif step == "report":
                step_report(browser, args.shots, args.project)
            elif step == "unavailable":
                step_unavailable(page, args.shots, args.missing, args.l1_only)
            else:
                raise SystemExit(f"unknown step {step}")
            print(f"step {step}: ok")
        browser.close()
    print(json.dumps(EVIDENCE, indent=1))
    if args.evidence:
        args.evidence.write_text(json.dumps(EVIDENCE, indent=1) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
