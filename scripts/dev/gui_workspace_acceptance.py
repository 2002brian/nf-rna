"""Drive the nf-rna Local Research Workspace in headless Chrome (v1.7.0 Phase 7a acceptance).

Usage (the GUI must already be running; nothing is started, imported or written):

    python scripts/dev/gui_workspace_acceptance.py --shots DIR --evidence FILE [--base URL]

Checks, on existing projects and runs only:
  projects   the project list with species, input, backend, samples, contrasts, last run, execution, QC
  datasets   a FASTQ folder (D5), the airway count/metadata folder, an empty folder
  workspace  GUI project (airway-gui-demo) and CLI project (airway-e2e) overviews and their run history
  switching  switch projects in the sidebar and back; the per-project run choice is kept
  qc         F-1 strandedness WARNING (REL_D5_NEG) and F-3 R warnings (REL_D3_COUNTS), project, run and results pages
  states     SUCCESS, FAILED (fixture) and INCOMPLETE (GUI-KILL-001) on the Validate & Run page, log filter
  results    Open results / Open report from the project overview go to the chosen run; figures render

Reuses the helpers of gui_browser_acceptance.py; needs the ``gui-dev`` extra.
This is a development tool, not part of the test suite.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from playwright.sync_api import Page, sync_playwright

sys.path.insert(0, str(Path(__file__).parent))
import gui_browser_acceptance as gui  # noqa: E402

DEMO = "/media/duobook02/brian/nf-rna/projects/airway-gui-demo"
E2E = "/home/brian/nf-rna-work/testdata/airway-e2e"
D5_NEG = "/home/brian/nf-rna-validation/v161-reliability/nfrna/rel_d5_neg"
D3_COUNTS = "/home/brian/nf-rna-validation/v161-reliability/nfrna/rel_d3_counts"
FAILED = "/home/brian/nf-rna-work/testdata/gui7a-fixtures/failed-demo"
D5_FASTQ = "/home/brian/nf-rna-validation/v161-phase3/data/D5/fastq"
AIRWAY = "/home/brian/nf-rna-work/testdata/airway"
EMPTY = "/media/duobook02/brian/nf-rna/fastq"
EVIDENCE: dict[str, object] = {}


def main_text(page: Page) -> str:
    return page.locator("[data-testid=stMain]").inner_text()


def sidebar_text(page: Page) -> str:
    return page.locator("[data-testid=stSidebar]").inner_text()


def expect(name: str, condition: bool, detail: object = "") -> None:
    EVIDENCE[name] = {"pass": bool(condition), "detail": detail}
    print(("PASS " if condition else "FAIL ") + name + (f": {detail}" if detail else ""))


def grid_text(page: Page) -> str:
    """Text of every data grid (their accessible table is hidden, so innerText would be empty)."""

    return page.evaluate("() => Array.from(document.querySelectorAll('[data-testid=stMain] [data-testid=stDataFrame]')).map(g => g.textContent).join('\\n')")


def pick(page: Page, label: str, text: str) -> None:
    """Choose a selectbox option by typing part of it (the list is filtered) and pressing Enter."""

    field = page.get_by_label(label, exact=True)
    field.click()
    field.fill(text)
    page.wait_for_timeout(300)
    page.keyboard.press("Enter")
    gui.settle(page)


def selected(page: Page, label: str) -> str:
    """The shown value of the main-area selectbox with this label."""

    return page.evaluate(
        "(label) => { const box = Array.from(document.querySelectorAll('[data-testid=stMain] [data-testid=stSelectbox]'))"
        ".find(b => (b.querySelector('label') || {}).textContent?.trim() === label);"
        " return box ? (box.querySelector('input') || {}).value || '' : ''; }", label)


def step_projects(page: Page, shots: Path) -> None:
    gui.goto(page, "", None)
    gui.shot(page, shots, "W1_projects")
    text = main_text(page) + grid_text(page)
    expect("projects: list shows GUI, CLI, FASTQ and fixture projects", all(name in text for name in ("airway-gui-demo", "airway_dex_e2e", "rel_d5_neg", "rel_d3_counts", "failed-demo")))
    columns = ("species", "quantification", "reference", "samples", "contrasts", "last run", "execution", "scientific QC")
    expect("projects: columns", all(column in text for column in columns), [column for column in columns if column not in text])


def step_datasets(page: Page, shots: Path) -> None:
    gui.goto(page, "datasets", None)
    gui.fill(page, "Dataset folder", D5_FASTQ)
    gui.shot(page, shots, "W2_datasets_fastq_d5")
    text = main_text(page)
    expect("datasets: D5 FASTQ folder counted", "FASTQ files\n14" in text or "FASTQ files\n\n14" in text, text[:400])
    expect("datasets: paired-end recognition lists samples", "D5_R_A" in grid_text(page) and "D5_F_B" in grid_text(page))
    expect("datasets: single-end files reported by the validator rule", "D5_SE_A" in text)
    gui.fill(page, "Dataset folder", AIRWAY)
    page.get_by_role("tab", name="Count matrix").click()
    gui.settle(page)
    gui.shot(page, shots, "W3_datasets_airway_counts")
    text = main_text(page)
    expect("datasets: airway count matrix recognised (validator)", "63677 genes × 8 samples" in text, [line for line in text.splitlines() if "genes ×" in line])
    page.get_by_role("tab", name="Metadata").click()
    gui.settle(page)
    expect("datasets: metadata preview", "SRR1039508" in grid_text(page))
    gui.fill(page, "Dataset folder", EMPTY)
    expect("datasets: empty mounted folder", "No FASTQ files, count matrices, metadata or samplesheets" in main_text(page))


def step_workspace(page: Page, shots: Path) -> None:
    gui.goto(page, "project", DEMO)
    gui.shot(page, shots, "W4_project_gui_demo")
    text = main_text(page)
    expect("workspace: GUI project summary", all(item in text + grid_text(page) for item in ("airway-gui-demo", "Homo sapiens", "trt_vs_untrt", "Analysis history")))
    grid = grid_text(page)
    expect("workspace: GUI project run history (4 runs, INCOMPLETE shown)", all(case in grid for case in ("GUI-6B-L2", "GUI-L2-002", "GUI-KILL-001", "GUI-L2-001")) and "INCOMPLETE" in grid)
    gui.goto(page, "project", E2E)
    gui.shot(page, shots, "W5_project_cli_e2e")
    text = main_text(page)
    expect("workspace: CLI project opens", "airway_dex_e2e" in text and "AIRWAY-E2E" in grid_text(page))


def step_switching(page: Page, shots: Path) -> None:
    gui.goto(page, "project", DEMO)
    pick(page, "Run", "GUI-KILL-001")
    pick(page, "Project", "rel_d5_neg")
    expect("switching: sidebar switch to rel_d5_neg", "rel_d5_neg" in page.locator("h1").first.inner_text() and "REL_D5_NEG" in grid_text(page))
    gui.shot(page, shots, "W6_switched_to_rel_d5_neg")
    pick(page, "Project", "airway-gui-demo")
    chosen = selected(page, "Run")
    expect("switching: back to airway-gui-demo keeps its run choice", "GUI-KILL-001" in chosen, chosen)
    expect("switching: URL follows the project", "airway-gui-demo" in page.url, page.url)
    gui.shot(page, shots, "W7_switched_back_choice_kept")


def step_qc(page: Page, shots: Path) -> None:
    gui.goto(page, "project", D5_NEG)
    text = main_text(page)
    expect("qc: F-1 WARNING on project overview", "Scientific QC: WARNING" in text and "Reads contradict the declared `forward` for: D5_R_A, D5_R_B" in grid_text(page))
    expect("qc: execution stays SUCCESS", "SUCCESS" in grid_text(page))
    gui.shot(page, shots, "W8_qc_f1_rel_d5_neg")
    gui.goto(page, "results", D5_NEG)
    text = main_text(page)
    expect("qc: F-1 strandedness on Results (v1.6.1 view unchanged)", "Strandedness check: WARNING" in text)
    gui.shot(page, shots, "W9_results_f1_strandedness")
    gui.goto(page, "project", D3_COUNTS)
    text = main_text(page)
    expect("qc: F-3 R warnings on project overview", "WARNINGS: 6" in grid_text(page) and "shrinkage: lfcShrink(type = apeglm)" in grid_text(page))
    gui.shot(page, shots, "W10_qc_f3_rel_d3_counts")
    gui.goto(page, "results", D3_COUNTS)
    page.get_by_text("Differential expression", exact=True).first.click()
    gui.settle(page, 1500)
    text = main_text(page)
    expect("qc: F-3 R warnings expander on Results (v1.6.1 view unchanged)", "R warnings recorded during the fit and this contrast (6)" in text)
    gui.shot(page, shots, "W11_results_f3_warnings")


def step_states(page: Page, shots: Path) -> None:
    gui.goto(page, "run", FAILED)
    text = main_text(page)
    expect("states: FAILED read from run_state.json", "FAILED: Phase 7a acceptance fixture" in text)
    expect("states: validation overview", "The project is valid and can be run." in text or "cannot run yet" in text)
    gui.shot(page, shots, "W12_run_failed_fixture")
    page.get_by_text("Only warnings and errors").click()
    gui.settle(page)
    log = page.locator("[data-testid=stCode]").last.inner_text()
    expect("states: log filter keeps WARNING/ERROR lines only", "ERROR downstream L1" in log and "WARNING synthetic" in log and "INFO case run started" not in log, log)
    gui.shot(page, shots, "W13_run_log_filter")
    gui.goto(page, "run", DEMO)
    pick(page, "Run", "GUI-KILL-001")
    text = main_text(page)
    expect("states: INCOMPLETE (GUI-KILL-001)", "INCOMPLETE — the run did not finish" in text)
    gui.shot(page, shots, "W14_run_incomplete")
    pick(page, "Run", "GUI-6B-L2")
    expect("states: SUCCESS (GUI-6B-L2)", "SUCCESS. Results:" in main_text(page))
    gui.shot(page, shots, "W15_run_success_settings")


def step_results(page: Page, shots: Path) -> None:
    gui.goto(page, "project", DEMO)
    pick(page, "Run", "GUI-L2-001")
    gui.click(page, "Open results")
    gui.settle(page, 2000)
    chosen = selected(page, "Run")
    expect("results: Open results preselects the chosen run", "results" in page.url and "GUI-L2-001" in chosen, chosen)
    widths = page.evaluate("() => Array.from(document.querySelectorAll('[data-testid=stMain] img')).map(i => i.complete ? i.naturalWidth : 0)")
    expect("results: QC figures render", bool(widths) and all(width > 0 for width in widths), widths)
    gui.shot(page, shots, "W16_results_from_project")
    gui.goto(page, "project", DEMO)
    pick(page, "Run", "GUI-6B-L2")
    gui.click(page, "Open report & downloads")
    gui.settle(page, 2000)
    text = main_text(page)
    expect("results: Open report preselects the chosen run", "reports" in page.url and "GUI-6B-L2" in selected(page, "Run") and "Download report_" in text, page.url)
    gui.shot(page, shots, "W17_reports_from_project")


STEPS = ("projects", "datasets", "workspace", "switching", "qc", "states", "results")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--shots", type=Path, required=True)
    parser.add_argument("--evidence", type=Path, required=True)
    parser.add_argument("--base", default="http://127.0.0.1:8701")
    parser.add_argument("steps", nargs="*", default=list(STEPS))
    args = parser.parse_args()
    gui.BASE = args.base
    args.shots.mkdir(parents=True, exist_ok=True)
    with sync_playwright() as playwright:
        browser = playwright.chromium.launch(executable_path=gui.CHROME, headless=True)
        page = browser.new_context(viewport={"width": 1440, "height": 1000}).new_page()
        for step in args.steps:
            print(f"== {step}")
            try:
                globals()[f"step_{step}"](page, args.shots)
            except Exception as exc:  # record and continue with the next step
                expect(f"{step}: completed without error", False, repr(exc)[:500])
                gui.shot(page, args.shots, f"ERROR_{step}")
        browser.close()
    args.evidence.write_text(json.dumps(EVIDENCE, indent=2) + "\n", encoding="utf-8")
    failed = [name for name, item in EVIDENCE.items() if not item["pass"]]
    print(f"{len(EVIDENCE) - len(failed)}/{len(EVIDENCE)} checks passed")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
