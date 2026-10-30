"""Drive the nf-rna GUI in headless Chrome for the v1.6.0 Phase 6a acceptance.

Usage (the GUI must already be running on 127.0.0.1:8501):

    python scripts/dev/gui_browser_acceptance.py STEP [STEP ...] --shots DIR

Steps: create, inputs, metadata, contrasts, validate, run, status, existing,
fastq, pages.  Needs the ``gui-dev`` extra (Playwright); it uses the
installed Google Chrome, so no Playwright browser download is needed.
This is a development tool, not part of the test suite.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from urllib.parse import quote

from playwright.sync_api import Page, sync_playwright

BASE = "http://127.0.0.1:8501"
CHROME = "/usr/bin/google-chrome"
PROJECT = "/media/duobook02/brian/nf-rna/projects/airway-gui-demo"
AIRWAY = "/home/brian/nf-rna-work/testdata/airway"
EXISTING = "/home/brian/nf-rna-work/testdata/airway-e2e"
FASTQ_DIR = "/media/duobook02/brian/nf-rna/fastq/gui-6a-recognition-test"
FASTQ_PROJECT = "gui-6a-fastq-test"


def settle(page: Page, extra: int = 900) -> None:
    """Wait until Streamlit has finished the current rerun."""

    page.wait_for_timeout(extra)
    page.wait_for_function("() => !document.querySelector('[data-testid=stStatusWidget]')", timeout=180_000)
    page.wait_for_timeout(300)
    errors = page.locator("[data-testid=stException]")
    if errors.count():
        raise SystemExit("Streamlit exception: " + errors.first.inner_text())


def goto(page: Page, path: str, project: str | None = PROJECT) -> None:
    page.goto(f"{BASE}/{path}" + (f"?project={quote(project)}" if project else ""))
    page.wait_for_selector("[data-testid=stSidebar]")
    settle(page, 1500)


def fill(page: Page, label: str, value: str) -> None:
    field = page.get_by_label(label, exact=True)
    field.fill(value)
    field.press("Enter")
    settle(page)


def choose(page: Page, label: str, option: str) -> None:
    """Pick a selectbox option by typing it (Streamlit selectboxes are searchable)."""

    field = page.get_by_label(label, exact=True)
    field.click()
    field.fill(option)
    page.get_by_role("option", name=option, exact=True).click()
    settle(page)


def choose_containing(page: Page, label: str, text: str) -> None:
    page.get_by_label(label, exact=True).click()
    page.get_by_role("option").filter(has_text=text).first.click()
    settle(page)


def click(page: Page, name: str) -> None:
    page.get_by_role("button", name=name, exact=True).click()
    settle(page)


def shot(page: Page, shots: Path, name: str) -> None:
    page.screenshot(path=str(shots / f"{name}.png"), full_page=True)
    print(f"screenshot {name}.png")


def text_of(page: Page) -> str:
    return page.locator("[data-testid=stMain]").inner_text()


def grid_paste(page: Page, text: str, *, row: int, column: int, editor_index: int = 0) -> None:
    """Paste tab/newline separated cells into the data editor, as from an Excel copy.

    The grid is a canvas: the first data cell is focused with a click and the
    target reached with the arrow keys; the text is put on the browser
    clipboard and pasted with Ctrl+V.
    """

    grid = page.locator("[data-testid=stDataFrame]").nth(editor_index)
    grid.click(position={"x": 40, "y": 35 + 17})
    page.wait_for_timeout(200)
    for _ in range(column):
        page.keyboard.press("ArrowRight")
    for _ in range(row):
        page.keyboard.press("ArrowDown")
    page.evaluate("text => navigator.clipboard.writeText(text)", text)
    page.keyboard.press("Control+V")
    page.wait_for_timeout(800)


def step_create(page: Page, shots: Path) -> None:
    goto(page, "new", project=None)
    shot(page, shots, "B1_new_project_empty")
    fill(page, "Project name", "airway-gui-demo")
    page.get_by_text("raw_counts — Integer gene-level raw-count matrix").click()
    settle(page)
    page.get_by_text("L2 — L1 plus differential expression").click()
    settle(page)
    choose_containing(page, "Experimental design", "paired_two_group")
    fill(page, "Pairing metadata column", "cell")
    choose_containing(page, "Gene identifier type in the count matrix", "ENSEMBL")
    shot(page, shots, "B2_new_project_filled")
    click(page, "Create project")
    assert "Created project" in text_of(page), text_of(page)
    shot(page, shots, "B3_new_project_created")


def step_inputs(page: Page, shots: Path) -> None:
    goto(page, "inputs")
    fill(page, "Count matrix file (CSV, gene_id first column, integer counts)", f"{AIRWAY}/counts.csv")
    click(page, "Import count matrix")
    body = text_of(page)
    assert "63677 genes × 8 samples" in body, body
    shot(page, shots, "C1_counts_imported")


def step_metadata(page: Page, shots: Path) -> None:
    goto(page, "inputs")
    # Excel-style paste of the cell column straight into the grid (rows in count-matrix order).
    cells = ["N61311", "N61311", "N052611", "N052611", "N080611", "N080611", "N061011", "N061011"]
    grid_paste(page, "\n".join(cells), row=0, column=1)
    shot(page, shots, "C2_grid_paste_cell")
    click(page, "Save metadata")
    body = text_of(page)
    assert "Saved" in body and "missing" in body, body  # incomplete metadata may be saved
    shot(page, shots, "C3_saved_incomplete")
    # Save and reload: a fresh page must show the saved cell values.
    goto(page, "inputs")
    saved = Path(PROJECT, "metadata.csv").read_text(encoding="utf-8")
    assert "SRR1039508,N61311," in saved, saved
    # Batch-assign the condition (entered by the user, never derived from names).
    page.get_by_role("tab", name="Assign values to samples").click()
    for value, samples in (("untrt", ["SRR1039508", "SRR1039512", "SRR1039516", "SRR1039520"]), ("trt", ["SRR1039509", "SRR1039513", "SRR1039517", "SRR1039521"])):
        box = page.get_by_label("Samples", exact=True)
        for sample in samples:
            box.click()
            box.fill(sample)
            page.get_by_role("option", name=sample, exact=True).click()
            settle(page, 200)
        page.keyboard.press("Escape")
        assert all(page.get_by_text(sample, exact=True).count() for sample in samples)
        choose(page, "Column", "condition")
        fill(page, "Value", value)
        click(page, "Apply")
    shot(page, shots, "C4_batch_assign")
    click(page, "Save metadata")
    assert Path(PROJECT, "metadata.csv").read_bytes() == Path(AIRWAY, "metadata.csv").read_bytes()
    shot(page, shots, "C5_metadata_saved")
    # TSV fallback: paste the whole table as copied from Excel.
    tsv = Path(AIRWAY, "metadata.csv").read_text(encoding="utf-8").replace(",", "\t")
    page.get_by_role("tab", name="Paste table (TSV / Excel)").click()
    page.get_by_label("Paste rows copied from Excel (tab-separated) or CSV, header row first").fill(tsv)
    click(page, "Apply pasted table")
    click(page, "Save metadata")
    assert Path(PROJECT, "metadata.csv").read_bytes() == Path(AIRWAY, "metadata.csv").read_bytes()
    shot(page, shots, "C6_tsv_paste_saved")


def step_contrasts(page: Page, shots: Path) -> None:
    goto(page, "contrasts")
    shot(page, shots, "D1_contrasts_initial")
    choose(page, "Numerator (test group)", "trt")
    choose(page, "Denominator (reference group)", "untrt")
    click(page, "Add")
    body = text_of(page)
    assert "Positive log2FC = higher in trt" in body, body
    shot(page, shots, "D2_contrast_preview")
    replace = page.get_by_text("Replace the saved contrasts", exact=True)
    if replace.count():
        replace.click()
        settle(page)
    click(page, "Save contrasts")
    assert "Wrote 1 contrast(s)" in text_of(page), text_of(page)
    shot(page, shots, "D3_contrasts_saved")
    pin = page.get_by_role("button", name="Pin to the current snapshot")
    if pin.count():
        pin.first.click()
        settle(page, 1500)
    assert "Pinned to" in text_of(page), text_of(page)
    shot(page, shots, "D4_kegg_pinned")


def step_validate(page: Page, shots: Path) -> None:
    goto(page, "run")
    click(page, "Validate")
    body = text_of(page)
    assert "Validation result:" in body and "PASS" in body, body
    shot(page, shots, "E1_validate")


def step_run(page: Page, shots: Path, case_id: str) -> None:
    goto(page, "run")
    fill(page, "Case ID", case_id)
    button = page.get_by_role("button", name="Run L2", exact=True)
    # Two quick clicks: the second lands while the first submission is running.
    button.click()
    page.wait_for_timeout(150)
    try:
        button.click(timeout=2000)
    except Exception as exc:  # disabled or detached after the first submission
        print(f"second click not delivered: {type(exc).__name__}")
    settle(page, 3000)
    shot(page, shots, "E2_run_submitted")
    print(text_of(page)[:2000])


def step_status(page: Page, shots: Path, name: str, project: str = PROJECT) -> None:
    goto(page, "run", project=project)
    click(page, "Refresh")
    shot(page, shots, name)
    print(text_of(page)[-3000:])


def step_fastq(page: Page, shots: Path) -> None:
    goto(page, "new", project=None)
    fill(page, "Project name", FASTQ_PROJECT)
    page.get_by_text("L1 — Expression-level QC").click()
    settle(page)
    click(page, "Create project")
    assert "Created project" in text_of(page), text_of(page)
    goto(page, "inputs", project=f"/media/duobook02/brian/nf-rna/projects/{FASTQ_PROJECT}")
    fill(page, "FASTQ folder", FASTQ_DIR)
    body = text_of(page)
    assert "sample,fastq_1,fastq_2,strandedness" in body, body
    shot(page, shots, "C7_fastq_detected")
    page.get_by_role("button", name="Copy these").click()
    settle(page)
    assert "Imported" in text_of(page) or "FASTQ input" in text_of(page), text_of(page)
    shot(page, shots, "C8_fastq_imported")


def step_pages(page: Page, shots: Path) -> None:
    goto(page, "", project=PROJECT)
    shot(page, shots, "A1_dashboard")
    for path, name in (("new", "B0_new_project"), ("inputs", "C0_inputs"), ("contrasts", "D0_contrasts"), ("run", "E0_run")):
        goto(page, path)
        shot(page, shots, name)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("steps", nargs="+")
    parser.add_argument("--shots", type=Path, required=True)
    parser.add_argument("--case-id", default="GUI-L2-001")
    parser.add_argument("--name", default="E9_status")
    parser.add_argument("--project", default=PROJECT)
    args = parser.parse_args()
    args.shots.mkdir(parents=True, exist_ok=True)
    with sync_playwright() as playwright:
        browser = playwright.chromium.launch(executable_path=CHROME, headless=True)
        context = browser.new_context(viewport={"width": 1400, "height": 1000}, permissions=["clipboard-read", "clipboard-write"])
        page = context.new_page()
        for step in args.steps:
            print(f"== {step}")
            if step == "run":
                step_run(page, args.shots, args.case_id)
            elif step == "status":
                step_status(page, args.shots, args.name, args.project)
            else:
                globals()[f"step_{step}"](page, args.shots)
        browser.close()


if __name__ == "__main__":
    sys.exit(main())
