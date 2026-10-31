"""Download figures as PDF / TIFF through the GUI in headless Chrome (v1.7.0 Phase 7b acceptance).

Usage (the GUI must already be running; it is not started here):

    python scripts/dev/gui_figure_download_acceptance.py --shots DIR --out DIR --evidence FILE [--base URL]

On run GUI-6B-L2 of airway-gui-demo (read-only):
  results   Results > Differential expression: the volcano's PDF, TIFF 300 dpi and TIFF 600 dpi buttons
  reports   Reports & downloads > Figures: a GSEA running-score plot (no pipeline TIFF) as TIFF 300,
            a GO dotplot as PDF and TIFF 600, and the PCA as TIFF 300 (the default format)

Every downloaded file is compared with its source: PDF and pipeline TIFF byte for byte, a converted
running-score TIFF pixel for pixel with the pipeline's 300 dpi PNG, and every TIFF for format, RGB,
LZW, dpi tags and pixel size = round(page_pt / 72 * dpi).  Needs the ``gui-dev`` extra.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import subprocess
import sys
from pathlib import Path

from PIL import Image, ImageChops
from playwright.sync_api import Page, sync_playwright

sys.path.insert(0, str(Path(__file__).parent))
import gui_browser_acceptance as gui  # noqa: E402
import gui_workspace_acceptance as ws  # noqa: E402

from rnaseq.gui import figure_export  # noqa: E402

DEMO = "/media/duobook02/brian/nf-rna/projects/airway-gui-demo"
RUN = Path(DEMO) / "runs" / "GUI-6B-L2" / "20261030-210152+0800"
DOWN = RUN / "downstream"
EVIDENCE = ws.EVIDENCE


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def save_download(page: Page, button, out: Path, prefix: str) -> Path:
    with page.expect_download(timeout=120_000) as info:
        button.click()
    download = info.value
    target = out / f"{prefix}__{download.suggested_filename}"
    download.save_as(target)
    return target


def tiff_facts(path: Path, stem: Path, dpi: int) -> tuple[bool, dict]:
    with Image.open(path) as image:
        image.load()
        expected = figure_export.pixels_for(figure_export.pdf_page_points(stem.with_suffix(".pdf")), dpi)
        facts = {"format": image.format, "mode": image.mode, "compression": image.info.get("compression"),
                 "dpi": [round(float(value)) for value in image.info.get("dpi", (0, 0))], "size": list(image.size), "expected_size": list(expected)}
    ok = facts["format"] == "TIFF" and facts["mode"] == "RGB" and facts["compression"] == "tiff_lzw" and facts["dpi"] == [dpi, dpi] and tuple(facts["size"]) == expected
    return ok, facts


def step_results(page: Page, shots: Path, out: Path) -> None:
    gui.goto(page, "results", DEMO)
    ws.pick(page, "Run", "GUI-6B-L2")
    page.get_by_text("Differential expression", exact=True).click()
    gui.settle(page, 2000)
    stem = DOWN / "l2" / "contrasts" / "trt_vs_untrt" / "volcano"
    main = page.locator("[data-testid=stMain]")
    main.get_by_role("button", name="TIFF 600 dpi", exact=True).first.scroll_into_view_if_needed()
    page.mouse.wheel(0, 250)
    page.wait_for_timeout(500)
    page.screenshot(path=str(shots / "F1_results_deg_figure_buttons.png"))
    pdf = save_download(page, main.get_by_role("button", name="PDF", exact=True).first, out, "results")
    ws.expect("results: volcano PDF download equals the pipeline PDF", pdf.read_bytes() == stem.with_suffix(".pdf").read_bytes(), pdf.name)
    opened = subprocess.run(["pdfinfo", str(pdf)], capture_output=True, text=True)
    ws.expect("results: downloaded PDF opens", opened.returncode == 0 and re.search(r"^Pages:\s+1$", opened.stdout, re.MULTILINE) is not None,
              [line for line in opened.stdout.splitlines() if line.startswith(("Pages", "Page size", "Producer"))] or opened.stderr)
    t300 = save_download(page, main.get_by_role("button", name="TIFF 300 dpi", exact=True).first, out, "results")
    ok, facts = tiff_facts(t300, stem, 300)
    ws.expect("results: volcano TIFF 300 equals the pipeline TIFF and is 300 dpi", ok and t300.read_bytes() == stem.with_suffix(".tiff").read_bytes(), {"file": t300.name, **facts})
    t600 = save_download(page, main.get_by_role("button", name="TIFF 600 dpi", exact=True).first, out, "results")
    ok, facts = tiff_facts(t600, stem, 600)
    ws.expect("results: volcano TIFF 600 is 600 dpi at the figure's physical size", ok and t600.name.endswith("volcano_600dpi.tiff"), {"file": t600.name, **facts})


def report_download(page: Page, out: Path, figure: str, form: str, prefix: str) -> Path:
    ws.pick(page, "Figure (73)", figure)
    page.get_by_text(form, exact=True).click()
    gui.settle(page)
    section = page.locator("[data-testid=stMain]")
    return save_download(page, section.get_by_role("button", name="Download", exact=True).last, out, prefix)


def step_reports(page: Page, shots: Path, out: Path) -> None:
    gui.goto(page, "reports", DEMO)
    ws.pick(page, "Run", "GUI-6B-L2")
    gui.settle(page, 1500)
    text = ws.main_text(page)
    ws.expect("reports: Figures (PDF / TIFF) section with 73 figures, TIFF 300 dpi by default",
              "Figures (PDF / TIFF)" in text and "Figure (73)" in text and "TIFF 300 dpi: pipeline TIFF, unchanged" in text, "caption shown")
    pca = report_download(page, out, "l1/pca", "TIFF 300 dpi", "reports")
    ws.expect("reports: PCA TIFF 300 equals the pipeline TIFF", pca.read_bytes() == (DOWN / "l1" / "pca.tiff").read_bytes(), pca.name)
    running = "l2/enrichment/gsea_kegg/trt_vs_untrt/running_score/negative_1_hsa04115"
    tif = report_download(page, out, running, "TIFF 300 dpi", "reports")
    page.get_by_text("a raster is never enlarged", exact=False).scroll_into_view_if_needed()
    page.screenshot(path=str(shots / "F2_reports_running_score_tiff300.png"))
    stem = DOWN / running
    ok, facts = tiff_facts(tif, stem, 300)
    with Image.open(tif) as a, Image.open(stem.with_name(stem.name + "_300dpi.png")) as b:
        same = ImageChops.difference(a.convert("RGB"), b.convert("RGB")).getbbox() is None
    ws.expect("reports: running-score TIFF 300 is the pipeline 300 dpi PNG, pixel for pixel", ok and same and "lossless conversion" in ws.main_text(page), {"file": tif.name, **facts})
    dot = "l2/enrichment/go/trt_vs_untrt/up/dotplot_BP"
    pdf = report_download(page, out, dot + "_simplified", "PDF", "reports")
    ws.expect("reports: dotplot PDF equals the pipeline PDF", pdf.read_bytes() == (DOWN / (dot + "_simplified.pdf")).read_bytes(), pdf.name)
    t600 = report_download(page, out, dot, "TIFF 600 dpi", "reports")
    page.get_by_text("a raster is never enlarged", exact=False).scroll_into_view_if_needed()
    page.screenshot(path=str(shots / "F3_reports_dotplot_tiff600.png"))
    ok, facts = tiff_facts(t600, DOWN / dot, 600)
    ws.expect("reports: dotplot TIFF 600 is 600 dpi at the figure's physical size", ok and "rendered from the pipeline's vector PDF" in ws.main_text(page), {"file": t600.name, **facts})


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--shots", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--evidence", type=Path, required=True)
    parser.add_argument("--base", default="http://127.0.0.1:8701")
    args = parser.parse_args()
    gui.BASE = args.base
    for folder in (args.shots, args.out):
        folder.mkdir(parents=True, exist_ok=True)
    before = {str(path): sha(path) for path in sorted(RUN.rglob("*")) if path.is_file()}
    with sync_playwright() as playwright:
        browser = playwright.chromium.launch(executable_path=gui.CHROME, headless=True)
        page = browser.new_context(viewport={"width": 1440, "height": 1000}, accept_downloads=True).new_page()
        for step in (step_results, step_reports):
            print(f"== {step.__name__}")
            try:
                step(page, args.shots, args.out)
            except Exception as exc:  # record and continue
                ws.expect(f"{step.__name__}: completed without error", False, repr(exc)[:500])
                gui.shot(page, args.shots, f"ERROR_{step.__name__}")
        browser.close()
    after = {str(path): sha(path) for path in sorted(RUN.rglob("*")) if path.is_file()}
    ws.expect("run files unchanged (SHA-256)", before == after, f"{len(before)} files")
    args.evidence.write_text(json.dumps(EVIDENCE, indent=2) + "\n", encoding="utf-8")
    failed = [name for name, item in EVIDENCE.items() if not item["pass"]]
    print(f"{len(EVIDENCE) - len(failed)}/{len(EVIDENCE)} checks passed")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
