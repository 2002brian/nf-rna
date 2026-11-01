"""Drive the nf-rna GUI in headless Chrome for the v1.7.0 pre-release GUI hardening.

Usage (the GUI must already be running with a gui.yaml whose first allowed root is DATA; it is not started here):

    python scripts/dev/gui_hardening_acceptance.py STEP [STEP ...] --base URL --data DIR --shots DIR --evidence FILE [--out DIR]

Steps (fixtures under DATA are made by the caller; formal runs are only read):
  picker     "Use this folder" / "Use this file" in the path picker; a symlink out of the roots is not offered
  url        ?project= outside the allowed roots (directly and through a symlink) is refused
  metadata   grid edits survive bulk assignment, a page switch and a paste; replacing needs a tick; saved file
  contrasts  disabled reasons for Add / Save contrasts; the saved table is redrawn after saving
  run        disabled reason for Run; the full validation report is dropped after the project changes
  live       a RUNNING run refreshes on its own and stops refreshing after it ends (needs --live-pid-file)
  sorting    DEG and enrichment tables of GUI-6B-L2 sort numerically by column header
  reference  New Project > Other reference: fields visible, no premature error
  figures    the Phase 7b PDF / TIFF download checks (Results and Reports)
  pages      Projects, Datasets, Project overview, Results, Reports of airway-gui-demo open without an exception

Needs the ``gui-dev`` extra. This is a development tool, not part of the test suite.
"""

from __future__ import annotations

import argparse
import json
import re
import sys
import time
from pathlib import Path

from playwright.sync_api import Page, sync_playwright

sys.path.insert(0, str(Path(__file__).parent))
import gui_browser_acceptance as gui  # noqa: E402
import gui_workspace_acceptance as ws  # noqa: E402

from rnaseq.gui import results  # noqa: E402

DEMO = "/media/duobook02/brian/nf-rna/projects/airway-gui-demo"
EVIDENCE = ws.EVIDENCE
expect = ws.expect


def main_text(page: Page) -> str:
    return ws.main_text(page)


def browse_open(page: Page, folder: str) -> None:
    ws.pick(page, "Subfolders", folder)
    gui.click(page, "Open folder")


def step_picker(page: Page, data: Path, shots: Path) -> None:
    gui.goto(page, "datasets", project=None)
    gui.click(page, "Go to root")
    page.get_by_label("Subfolders", exact=True).click()
    offered = [item.inner_text() for item in page.get_by_role("option").all()]
    page.keyboard.press("Escape")
    expect("picker: a symlink leaving the allowed roots is not offered", "escape" not in offered and "fastq" in offered, offered)
    expect("picker: hidden links are reported", "not shown" in main_text(page), "caption present")
    browse_open(page, "fastq")
    browse_open(page, "hard-fq")
    gui.click(page, "Use this folder")  # settle() fails on any Streamlit exception
    field = page.get_by_label("Dataset folder", exact=True).input_value()
    text = main_text(page)
    expect("picker: Use this folder selects the folder, no exception", field == str(data / "fastq" / "hard-fq") and "FASTQ files" in text, field)
    gui.shot(page, shots, "H1_use_this_folder")

    gui.goto(page, "inputs", project=str(data / "projects" / "hard-meta"))
    ws.pick(page, "Files", "counts.csv")  # the Browse expander is open while the field is empty
    gui.click(page, "Use this file")
    field = page.get_by_label("Count matrix file (CSV, gene_id first column, integer counts)", exact=True).input_value()
    expect("picker: Use this file selects the file, no exception", field == str(data / "counts.csv") and page.get_by_role("button", name="Import count matrix").count() == 1, field)
    gui.shot(page, shots, "H2_use_this_file")


def step_url(page: Page, data: Path, shots: Path) -> None:
    outside = data.parent / "outside" / "outsider"
    for name, value in (("direct", str(outside)), ("symlink", str(data / "escape" / "outsider")), ("dotdot", str(data / ".." / "outside" / "outsider"))):
        gui.goto(page, "project", project=value)
        text = main_text(page) + ws.sidebar_text(page)
        refused = "was not opened" in text and "outsider" not in page.locator("h1").all_inner_texts()
        expect(f"url: ?project= outside the roots is refused ({name})", refused and "project=" not in page.url, page.url)
    gui.shot(page, shots, "H3_url_refused")
    gui.goto(page, "project", project=str(data / "projects" / "hard-meta"))
    expect("url: ?project= inside the roots opens", "hard-meta" in page.locator("h1").all_inner_texts(), page.url)


def grid_column(page: Page, values: list[str], column: int) -> None:
    gui.grid_paste(page, "\n".join(values), row=0, column=column, editor_index=0)
    gui.settle(page, 600)


def step_metadata(page: Page, data: Path, shots: Path) -> None:
    project = data / "projects" / "hard-meta"
    gui.goto(page, "inputs", project=str(project))
    editor = page.locator("[data-testid=stDataFrame]").first
    editor.scroll_into_view_if_needed()
    # Typed into the grid (unsaved): donor and batch columns.
    grid_column(page, ["P1", "P2", "P3", "P1", "P2", "P3"], 1)
    grid_column(page, ["b1", "b1", "b2", "b2", "b2", "b1"], 2)
    text = main_text(page)
    expect("metadata: grid edits are marked unsaved", "Unsaved changes" in text, re.findall(r"Unsaved changes[^\n]*", text)[:1])
    gui.shot(page, shots, "H4_grid_edits_unsaved")

    # Bulk assignment of condition must keep the typed donor and batch values.
    page.get_by_role("tab", name="Assign values to samples").click()
    gui.settle(page)
    for value, samples in (("Control", ["C1", "C2", "C3"]), ("Treatment", ["T1", "T2", "T3"])):
        box = page.get_by_label("Samples", exact=True)
        for sample in samples:
            box.click()
            box.fill(sample)
            page.get_by_role("option", name=sample, exact=True).click()
            gui.settle(page, 200)
        page.keyboard.press("Escape")
        ws.pick(page, "Column", "condition")
        gui.fill(page, "Value", value)
        gui.click(page, "Apply")
    grid = ws.grid_text(page)
    kept = all(token in grid for token in ("P1", "P2", "P3", "b1", "b2", "Control", "Treatment"))
    expect("metadata: typed donor/batch values survive bulk assignment", kept and "condition = 'Treatment' applied" in main_text(page), grid[:300])
    gui.shot(page, shots, "H5_bulk_kept_edits")

    # Typed edit, then a page switch and back: the edit is still there.
    grid_column(page, ["P1b"], 1)
    # In-app navigation (sidebar links); a URL reload would start a new session.
    page.locator("[data-testid=stSidebarNav]").get_by_text("Contrasts", exact=True).click()
    gui.settle(page, 1500)
    assert "Contrasts" in page.locator("h1").all_inner_texts(), page.locator("h1").all_inner_texts()
    page.locator("[data-testid=stSidebarNav]").get_by_text("Input & Metadata", exact=True).click()
    gui.settle(page, 1500)
    grid = ws.grid_text(page)
    expect("metadata: an unsaved grid edit survives a page switch", "P1b" in grid and "Unsaved changes" in main_text(page), grid[:300])

    # Paste a table that would replace a typed value: needs a tick.
    page.get_by_role("tab", name="Paste table (TSV / Excel)").click()
    area = page.get_by_label("Paste rows copied from Excel (tab-separated) or CSV, header row first (must include sample_id)")
    area.fill("sample_id\tdonor\tRIN\nC1\tP1\t9.1\nC2\t\t8.7\n")
    area.press("Control+Enter")
    gui.settle(page)
    apply = page.get_by_role("button", name="Apply pasted table", exact=True)
    text = main_text(page)
    expect("metadata: a paste that replaces a typed value needs a tick", apply.is_disabled() and "would be replaced" in text and "P1b" in ws.grid_text(page), re.findall(r"\d+ existing value\(s\) would be replaced", text))
    gui.shot(page, shots, "H6_paste_needs_confirmation")
    page.get_by_text("Replace these 1 existing value(s)", exact=True).click()
    gui.settle(page)
    gui.click(page, "Apply pasted table")
    gui.click(page, "Save metadata")
    saved = (project / "metadata.csv").read_text(encoding="utf-8")
    expected = ("sample_id,donor,batch,condition,RIN\nC1,P1,b1,Control,9.1\nC2,P2,b1,Control,8.7\nC3,P3,b2,Control,\n"
                "T1,P1,b2,Treatment,\nT2,P2,b2,Treatment,\nT3,P3,b1,Treatment,\n")
    expect("metadata: saved metadata.csv has every typed, assigned and pasted value", saved == expected, saved)
    expect("metadata: no unsaved changes after saving", "No unsaved changes" in main_text(page), "caption")
    gui.shot(page, shots, "H7_metadata_saved")


def step_contrasts(page: Page, data: Path, shots: Path) -> None:
    gui.goto(page, "contrasts", project=str(data / "projects" / "hard-meta"))
    text = main_text(page)
    expect("contrasts: Save contrasts says why it is disabled", "Save contrasts is disabled: Add at least one contrast" in text, "caption")
    gui.choose(page, "Numerator (test group)", "Treatment")
    gui.choose(page, "Denominator (reference group)", "Treatment")
    expect("contrasts: Add says why it is disabled", "must be different levels" in main_text(page), "caption")
    gui.choose(page, "Denominator (reference group)", "Control")
    gui.click(page, "Add")
    gui.shot(page, shots, "H8_contrast_composed")
    gui.click(page, "Save contrasts")
    text = main_text(page)
    expect("contrasts: saved, and the saved table is redrawn", "Wrote 1 contrast(s)" in text and "Treatment_vs_Control" in ws.grid_text(page), "success + table")
    gui.shot(page, shots, "H9_contrast_saved")


def step_run(page: Page, data: Path, shots: Path) -> None:
    project = data / "projects" / "hard-meta"
    gui.goto(page, "run", project=str(project))
    expect("run: Run says why it is disabled (no case ID)", "is disabled: enter a case ID" in main_text(page), "caption")
    gui.click(page, "Validate and refresh the plan")
    expect("run: full validation report shown", "Full validation report" in main_text(page), "expander")
    metadata = project / "metadata.csv"
    metadata.write_text(metadata.read_text(encoding="utf-8").replace("8.7", "8.8"), encoding="utf-8")
    # Same session (sidebar navigation), so the stored report would still be there.
    page.locator("[data-testid=stSidebarNav]").get_by_text("Contrasts", exact=True).click()
    gui.settle(page, 1500)
    page.locator("[data-testid=stSidebarNav]").get_by_text("Validate & Run", exact=True).click()
    gui.settle(page, 1500)
    text = main_text(page)
    expect("run: the report of an older project state is removed", "Full validation report" not in text and "full validation report was removed" in text, "caption")
    gui.shot(page, shots, "H10_stale_report_removed")


def _status_stamp(page: Page) -> str | None:
    found = re.search(r"(?:last read|Read from the run files at) (\d\d:\d\d:\d\d)", main_text(page))
    return found.group(1) if found else None


def step_live(page: Page, data: Path, shots: Path, run_dir: Path, pid_file: Path) -> None:
    project = data / "projects" / "hard-live"
    gui.goto(page, "run", project=str(project))
    first = _status_stamp(page)
    text = main_text(page)
    expect("live: a RUNNING run announces automatic updates", "updates automatically every 10 s" in text and "How to stop this run" in text, first)
    page.get_by_label("Case ID", exact=True).fill("UNSUBMITTED")  # typed, not submitted: must survive the refreshes
    page.wait_for_timeout(23_000)
    second = _status_stamp(page)
    expect("live: the status section re-reads the run files without a click", first is not None and second is not None and second != first, [first, second])
    expect("live: other widgets are not reset by the refresh", page.get_by_label("Case ID", exact=True).input_value() == "UNSUBMITTED", "case ID kept")
    gui.shot(page, shots, "H11_live_refresh")
    state = json.loads((run_dir / "run_state.json").read_text())
    (run_dir / "run_state.json").write_text(json.dumps({**state, "status": "SUCCESS", "completed_at": "2026-11-01T12:30:00+08:00"}))
    page.wait_for_timeout(14_000)
    text = main_text(page)
    expect("live: the page shows the final state on its own", "SUCCESS" in text and "press Refresh to update" in text, _status_stamp(page))
    stopped = _status_stamp(page)
    page.wait_for_timeout(23_000)
    expect("live: refreshing stops after the run ended", _status_stamp(page) == stopped, [stopped, _status_stamp(page)])
    gui.shot(page, shots, "H12_live_stopped")


def _grid_rows(page: Page, grid) -> tuple[list[str], list[list[str]]]:
    """Header and the rendered rows of a data grid, from its accessible table (cell values as text)."""

    rows = page.evaluate("(el) => Array.from(el.querySelectorAll('tr')).map(r => Array.from(r.children).map(c => c.textContent))", grid.element_handle())
    return rows[0], rows[1:]


def _click_header(page: Page, grid, header: str, offset: int = 0) -> None:
    """Click a column header of the canvas grid; columns are found from the vertical separator lines of the header row."""

    from io import BytesIO

    from PIL import Image

    grid.scroll_into_view_if_needed()
    page.wait_for_timeout(300)
    box = grid.bounding_box()
    names, _rows = _grid_rows(page, grid)
    strip = Image.open(BytesIO(page.screenshot(clip={"x": box["x"], "y": box["y"] + 4, "width": box["width"], "height": 28}))).convert("L")
    width, height = strip.size
    background = strip.getpixel((width // 2, 2))
    lines = [x for x in range(2, width - 2) if all(abs(strip.getpixel((x, y)) - background) > 6 for y in range(4, height - 4))]
    edges = [x for index, x in enumerate(lines) if index == 0 or x - lines[index - 1] > 3]
    bounds = [0, *edges, width]
    if header not in names:
        raise SystemExit(f"{header!r} not among the rendered columns {names}")
    column = names.index(header) + offset  # offset: a row-selection checkbox column drawn before the data
    page.mouse.click(box["x"] + (bounds[column] + bounds[column + 1]) / 2, box["y"] + 18)
    page.wait_for_timeout(800)


def _numbers(rows: list[list[str]], index: int) -> list[float]:
    return [float(row[index]) for row in rows if row[index] not in ("", "None", "NA")]


def _check_sorted(page: Page, grid, header: str, label: str, shots: Path, shot: str, source: Path, offset: int = 0) -> None:
    """Sort by clicking the header; the rows must be in numeric order and start at the file's numeric minimum."""

    _click_header(page, grid, header, offset)
    names, rows = _grid_rows(page, grid)
    index = next(position for position, name in enumerate(names) if name.endswith(header))
    values = _numbers(rows, index)
    column = [row[header] for row in results.read_table(source).rows if row[header] not in results.NA_TOKENS]
    numeric_min = min(column, key=float)
    lexical_min = min(column)
    ok = len(values) > 5 and values == sorted(values) and values[0] == float(numeric_min)
    expect(f"sorting: {label} sorts numerically by its header", ok,
           {"header": names[index], "first rows": [row[index] for row in rows[:5]], "file numeric min": numeric_min, "file lexical min": lexical_min})
    gui.shot(page, shots, shot)


def step_sorting(page: Page, data: Path, shots: Path) -> None:
    gui.goto(page, "results", DEMO)
    ws.pick(page, "Run", "GUI-6B-L2")
    page.get_by_text("Differential expression", exact=True).click()
    gui.settle(page, 2000)
    grid = page.locator("[data-testid=stMain] [data-testid=stDataFrame]").last
    deg = Path(DEMO) / "runs" / "GUI-6B-L2" / "20261030-210152+0800" / "downstream" / "l2" / "contrasts" / "trt_vs_untrt" / "all_genes.tsv"
    _check_sorted(page, grid, "baseMean", "DEG table baseMean", shots, "H13_deg_sorted_baseMean", deg)
    _check_sorted(page, grid, "log2FoldChange", "DEG table log2FoldChange (negative values)", shots, "H14_deg_sorted_lfc", deg)
    page.get_by_text("Enrichment", exact=True).click()
    gui.settle(page, 2000)
    page.get_by_text("KEGG GSEA", exact=True).click()
    gui.settle(page, 2000)
    page.get_by_text("Show all evaluated terms (default: significant only)", exact=True).click()
    gui.settle(page, 2000)
    grid = page.locator("[data-testid=stMain] [data-testid=stDataFrame]").filter(has_text="enrichmentScore").last
    terms = Path(DEMO) / "runs" / "GUI-6B-L2" / "20261030-210152+0800" / "downstream" / "l2" / "enrichment" / "gsea_kegg" / "trt_vs_untrt" / "all_terms.tsv"
    _check_sorted(page, grid, "NES", "KEGG GSEA terms NES", shots, "H15_gsea_sorted_nes", terms, offset=1)


def step_reference(page: Page, data: Path, shots: Path) -> None:
    gui.goto(page, "new", project=None)
    gui.fill(page, "Project name", "never-created")
    gui.choose(page, "Reference", "Other (advanced)")
    text = main_text(page)
    errors = page.locator("[data-testid=stMain] [data-testid=stAlertContentError]").all_inner_texts()
    expect("reference: Other shows its fields at once, without an error", "Managed reference folder" in text and not errors and "Enter the managed reference folder" in text, errors)
    gui.shot(page, shots, "H16_other_reference")
    expect("reference: no project was created", not (data / "projects" / "never-created").exists(), "")


def step_pages(page: Page, data: Path, shots: Path) -> None:
    for path in ("", "datasets", "project", "results", "reports", "contrasts", "inputs", "run"):
        gui.goto(page, path, DEMO)
        expect(f"pages: {path or 'projects'} opens without an exception", page.locator("[data-testid=stException]").count() == 0, page.locator("h1").all_inner_texts())


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("steps", nargs="+")
    parser.add_argument("--base", required=True)
    parser.add_argument("--data", type=Path, required=True)
    parser.add_argument("--shots", type=Path, required=True)
    parser.add_argument("--evidence", type=Path, required=True)
    parser.add_argument("--out", type=Path)
    parser.add_argument("--live-run", type=Path)
    parser.add_argument("--live-pid-file", type=Path)
    args = parser.parse_args()
    gui.BASE = args.base
    args.shots.mkdir(parents=True, exist_ok=True)
    with sync_playwright() as playwright:
        browser = playwright.chromium.launch(executable_path=gui.CHROME, headless=True)
        context = browser.new_context(viewport={"width": 1500, "height": 1100}, permissions=["clipboard-read", "clipboard-write"], accept_downloads=True)
        page = context.new_page()
        for step in args.steps:
            print(f"== {step}")
            started = time.time()
            if step == "live":
                step_live(page, args.data, args.shots, args.live_run, args.live_pid_file)
            elif step == "figures":
                import gui_figure_download_acceptance as figures

                args.out.mkdir(parents=True, exist_ok=True)
                figures.step_results(page, args.shots, args.out)
                figures.step_reports(page, args.shots, args.out)
            else:
                globals()[f"step_{step}"](page, args.data, args.shots)
            print(f"   {time.time() - started:.0f} s")
        browser.close()
    previous = json.loads(args.evidence.read_text()) if args.evidence.is_file() else {}
    args.evidence.write_text(json.dumps({**previous, **EVIDENCE}, indent=2) + "\n")
    failed = [name for name, item in EVIDENCE.items() if not item["pass"]]
    print(f"{len(EVIDENCE) - len(failed)}/{len(EVIDENCE)} passed")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
