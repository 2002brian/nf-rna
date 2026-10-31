# nf-rna v1.7.0 Phase 7b — PDF/TIFF figure downloads

Branch `feature/v1.7.0`, after Phase 7a (`1385318`). Date 2026-10-31.
Scope: download options only. No R code, figure, threshold, statistic, manifest,
run or delivery package was changed; no pipeline was run.

## 1. Inventory

Every formal figure is drawn by `nf_rna_save_figure` (`src/rnaseq/r/figures.R`)
directly at its final physical size on each device. Run GUI-6B-L2 of
airway-gui-demo (profile `general`, raster 300 dpi) holds 73 figures:

| Figure | Count | PNG preview (150 dpi) | `_300dpi.png` | TIFF (300 dpi, LZW, RGB) | Cairo PDF |
|---|---|---|---|---|---|
| PCA, sample correlation, boxplots (L1) | 4 | yes | – | yes | yes |
| Volcano, MA plot, DEG heatmap (L2) | 3 | yes | yes | yes | yes |
| GO / KEGG ORA dotplots (incl. GO simplified) | 21 | yes | yes | yes | yes |
| GO / KEGG GSEA dotplots (incl. simplified) | 7 | yes | yes | yes | yes |
| GSEA running-score plots | 38 | yes | yes | **no** | yes |

Before 7b the GUI offered PNG and PDF on *Results* and PNG/PDF files on
*Reports & Downloads*; no TIFF. The delivery package already contains
`figures/tiff_300dpi/` (unchanged).

## 2. What changed

* `src/rnaseq/gui/figure_export.py` (new, no Streamlit): chooses the source of a
  TIFF at 300 or 600 dpi and produces the bytes.
  1. pipeline TIFF whose recorded dpi equals the request → offered byte for byte;
  2. pipeline `<name>_<dpi>dpi.png` → lossless re-encoding to TIFF (pixels unchanged);
  3. pipeline PDF → poppler `pdftocairo` at exactly that dpi, cropped to
     `round(page_pt / 72 × dpi)` px (R's half-to-even `round`, as in `nf_rna_save_figure`);
  4. otherwise unavailable, with the reason. A raster is never enlarged.
  Output: TIFF, RGB on white, LZW, dpi in the resolution tags; encoded in memory on click.
* *Results*: each figure has **PNG · PDF** and **TIFF 300 dpi · TIFF 600 dpi** buttons;
  a missing source gives a disabled button whose tooltip says why.
* *Reports & Downloads*: new **Figures (PDF / TIFF)** section — figure selector (all
  figures of the run), format PDF / TIFF 300 dpi (default) / TIFF 600 dpi / PNG, and a
  caption naming the source. The old flat "Figures (PNG / PDF)" file group is replaced
  by this section.
* `results.py`: `Figure.stem`, `figure_stems(run)`; running-score listing recognises any
  `_<n>dpi.png` (the `elsevier` profile writes `_500dpi.png`).
* Tests: `tests/test_gui_figure_export.py` (6 small tests), two assertions in
  `tests/test_gui_results.py`.
* Dev tools: `scripts/dev/figure_export_check.py` (all figures of one run),
  `scripts/dev/gui_figure_download_acceptance.py` (real downloads in headless Chrome).
* Docs: `docs/gui.md` (Figure downloads section, limitation), `DECISIONS.md`,
  `README_latest.md`, `README_zh-TW.md`.

Result per source on GUI-6B-L2:

| Option | Source | Figures |
|---|---|---|
| TIFF 300 dpi | pipeline TIFF, unchanged | 35 |
| TIFF 300 dpi | lossless conversion of the pipeline's 300 dpi PNG | 38 (running-score) |
| TIFF 600 dpi | pipeline PDF rendered at 600 dpi | 73 |
| PDF | pipeline PDF, unchanged | 73 |

## 3. Verification

**Every figure of GUI-6B-L2** (`figure_export_check.py`, evidence
`/home/brian/nf-rna-validation/gui-7b/figure_export_check_GUI-6B-L2.json`): 73/73 pass.

* every PDF opens (`pdfinfo`);
* every TIFF (146): format TIFF, mode RGB, LZW, dpi tag = requested dpi, pixel size =
  `round(page_pt / 72 × dpi)` (e.g. heatmap 306 × 577 pt → 1275 × 2404 px at 300, 2550 × 4808 at 600);
* 300 dpi from pipeline TIFF: byte-identical; from running-score PNG: pixel-identical;
* 600 dpi: compared with an independent Ghostscript rendering of the same PDF at
  600 dpi, mean |difference| ≤ 1.48/255 (anti-aliasing only); every colour covering
  ≥ 0.2 % of the pipeline's own raster is present (0 missing in 73 figures);
* the run's 853 files: SHA-256 identical before and after.

**Browser** (`gui_figure_download_acceptance.py`, GUI on 127.0.0.1:8701 with the 7a
acceptance `gui.yaml`; the user's GUI on 8501 untouched): 10/10 pass
(`docs/dev/screenshots/7b/browser_checks.json`).

| Check | Result |
|---|---|
| Results › volcano PDF = pipeline PDF; opens (1 page, cairo 1.18.6) | PASS |
| Results › volcano TIFF 300 dpi = pipeline TIFF, 1004 × 1004 px, 300 dpi | PASS |
| Results › volcano TIFF 600 dpi: 2008 × 2008 px, 600 dpi, RGB, LZW | PASS |
| Reports › Figures section, 73 figures, TIFF 300 dpi default | PASS |
| Reports › PCA TIFF 300 dpi = pipeline TIFF | PASS |
| Reports › KEGG GSEA running-score TIFF 300 dpi = pipeline `_300dpi.png` pixel for pixel, 2125 × 1488 px | PASS |
| Reports › GO dotplot (simplified) PDF = pipeline PDF | PASS |
| Reports › GO dotplot TIFF 600 dpi: 2725 × 2658 px, 600 dpi | PASS |
| Run files unchanged (853 files, SHA-256) | PASS |

Screenshots: `docs/dev/screenshots/7b/F1`–`F3` (buttons, captions);
`F4` pipeline TIFF 300 dpi (left) beside the 600 dpi TIFF from the PDF (right, resampled).

**Tests**: GUI tests (`tests/test_gui*.py`) pass; full suite — see §6.
`git diff 1385318 -- src/rnaseq/r src/rnaseq/workflows workflow tools src/rnaseq/*.py` is empty.

## 4. Finding: PDF and raster layouts differ slightly (pre-existing, not scientific)

While verifying the 600 dpi option, the pipeline's own PDF and its PNG/TIFF of the
same figure turned out not to share exactly one layout: text rows and the label
column are measured slightly differently by R's PDF device (72 dpi metrics) and by
its raster devices, so the panel edges move. poppler and Ghostscript render the PDF
identically, so this is in the PDF, not in the conversion. Axis-line displacement
PDF vs pipeline raster on GUI-6B-L2:

| Figure | Max displacement |
|---|---|
| PCA | 0.96 pt |
| boxplots | 1.44 pt |
| running-score plots | 1.32 pt |
| MA plot, volcano | 2.76, 3.00 pt |
| ORA / GSEA dotplots (long term names) | up to 10.08 pt (≈ 3.6 mm) |
| sample correlation, heatmap | not measurable (no axis lines) |

Points, colours, labels, term order and legends are the same (F4); no label is
clipped. It is therefore a figure-layout inconsistency in v1.6.1 output, not a
change of any scientific value, and figure code is outside the 7b scope.
Consequence: the 600 dpi TIFF has the PDF's layout; the 300 dpi TIFF has the raster
layout. This is stated in the GUI caption and in `docs/gui.md`.

Proposed follow-up (not done): if PDF and raster must match exactly, either (a) let
`nf_rna_save_figure` also write a 600 dpi TIFF with R's raster device (new runs only),
or (b) investigate the metric difference with `tools/figure_review/compare_runs.py`
and record the tolerance in `docs/figure-standard.md`. Both need a figure-standard
decision by the owner.

## 5. Not done / limits

* 600 dpi needs poppler (`pdftocairo`, `pdfinfo`) on the server; it is a system package
  here and not part of `environment.yml`. Without it the option is shown as unavailable.
* No TIFF above 600 dpi, no SVG, no presets, no figure editing, no legends (as specified).
* Under the `elsevier` profile both TIFF options come from the PDF (its pipeline TIFF is
  500 dpi and is not relabelled); tested with a fixture, no real elsevier run exists here.
* Downloaded names are the figure's own (`volcano.tiff`, `volcano_600dpi.tiff`); two
  contrasts give the same file name, as with PNG/PDF before.
* Version stays 1.6.1; no CHANGELOG entry (release prep after 7b).

## 6. Full suite

1133 passed, 2 skipped, 0 failed (26 min; locked downstream R prefix `0e721ff3…-200d98f8759eaf17`). The two skips are the same as in Phase 7a (real SAMtools/featureCounts, real SJ project). New in 7b: the 6 tests of `test_gui_figure_export.py` (1135 collected in total).
