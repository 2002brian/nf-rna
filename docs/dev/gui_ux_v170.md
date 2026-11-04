# v1.7.0 GUI plotting, performance and UX

**When:** 2026-11-04, on `feature/v1.7.0` after `6ee4fe7`.
**Why:** the first real FASTQ GUI end-to-end run (NGS1130901, Salmon, L1) found three problems:
- inconsistent L1 figure formatting;
- slow button interactions;
- a workflow that was not self-explanatory.

**Scope:** figure presentation, GUI responsiveness and guidance. No analysis method, threshold, default or run semantics changed.
**Evidence:** `/home/brian/nf-rna-validation/gui-ux-170/`, containing the benchmarks, the L1 comparison runs, screenshots and test logs.

## 1. Baseline

### 1.1 Architecture

| Part | What it is |
|---|---|
| GUI | Streamlit 1.66.1. `src/rnaseq/gui/pages.py` holds every page, and `st.navigation` gives one script per page. Every widget interaction reruns the page function. State is read from the project files on every rerun; `st.session_state` holds only unsaved input. |
| Figures | `src/rnaseq/r/figures.R` is the owner-approved figure standard (`docs/figure-standard.md`): presets, profiles, bundled Arimo font, stable group colours, content-aware sizes. L1 figures are drawn by `l1_analysis.R`: PCA and the boxplots with ggplot2, the correlation heatmap with pheatmap. |
| GUI caching already present | doctor checks (10 min), the dataset scan (60 s), count-matrix recognition, reference-manifest checks and result tables (all keyed by mtime). The run-status section is a fragment that polls every 10 s only while a run is live, and logs are read from their end. |

### 1.2 Performance

Full-script rerun cost of each page (Streamlit AppTest, in-process, isolated `gui.yaml`):
- **Data:** a copy of the 8-sample Human Salmon project, with the real Ensembl 116 reference read-only and stub FASTQs.
- **Method:** five reruns after one cold run.

| Page | Rerun median before |
|---|---|
| New Project | 9,272 ms |
| Contrasts | 5,196 ms |
| Projects | 4,701 ms |
| Validate & Run | 2,657 ms |
| Project overview | 2,642 ms |
| Input & Metadata | 2,616 ms |
| Results | 31 ms |
| Reports | 9 ms |
| Datasets | 3 ms |

Root cause, from cProfile:
- `validate_project` runs on every rerun of Project overview, Input & Metadata, Contrasts (twice) and Validate & Run. It calls `load_local_reference`, which recomputes the **SHA-256 of the multi-GB reference FASTA, GTF and transcriptome**: 2.66 s of a 2.68 s call, 2.08 s of it in `hashlib.update`.
- `registered_references()` does the same for every registered reference: 4.8 s on the Projects page each rerun, and twice inside `reference_availability` (9.5 s) on New Project.

Everything else measured at most 11 ms.

### 1.3 L1 figures

The run used the standard preset with the verified Arimo font. Colours were already stable by group (`nf_rna_group_colours`). The remaining inconsistencies were:
- **Correlation heatmap text:** pheatmap draws its title at 1.3 × the base size (10.4 pt bold), and bold annotation names and legend titles, while the ggplot figures use an 8 pt bold title and plain legend titles.
- **No 300 dpi PNG for L1:** L2 figures have had `<name>_300dpi.png` since v1.5.0; L1 had only the 150 dpi preview, TIFF and PDF.
- **Results page layout:** figures were shown two per row at their preview pixel size (354–719 px), so their on-screen scale depended on the column. An optional `pca_by_batch` shifted the pairing (PCA + batch PCA, correlation + raw boxplot, …). Each figure had five separate download buttons.

### 1.4 UX

- No page except Project overview showed where the user is in the workflow.
- FASTQ preprocessing (raw vs pretrimmed) was hidden under *Advanced settings* and missing from the Run summary.
- New Project did not show the design formula it would create.
- Paired designs got feedback only from validation after saving.
- Detecting FASTQs vs copying them, and read-only validation vs writing `planning/`, were not distinguished.

## 2. Changes

| File | Change |
|---|---|
| `references.py` | Opt-in SHA-256 memo: `enable_checksum_memo()` / `disable_checksum_memo()`. It is **off by default**: the CLI, `rnaseq run` and every pipeline process hash in full as before. When on, a digest is reused only while the file identity (realpath, device, inode, size, mtime_ns, ctime_ns) is unchanged. It applies only to files of at least 16 MiB, so project metadata, contrasts and configs are always hashed. A digest is stored only if the identity is the same before and after reading. |
| `gui/pages.py` | `gui_config()` enables the memo. `validated()` wraps validation in a spinner that appears after 0.5 s and explains the one-time first check. Other additions: the workflow bar (Input & Metadata, Contrasts, Validate & Run, Results); New Project explanations, formula preview and missing-field line, with preprocessing visible; FASTQ preview labels; the pairing check; Run-page case ID, write/read-only and L1/L2 notes, plus preprocessing in the summary; grouped L1 figures at one scale; one *Download* menu per figure with PNG 300 dpi. |
| `gui/services.py` | `pairing_preview()`: the validator's paired-design rule, applied to the unsaved table (display only). |
| `gui/results.py` | `Figure.size_mm` / `group` and `figure_sizes()`, both from `figure_manifest.json`. |
| `gui/figure_export.py` | `png_option()` / `png_bytes()`: PNG 300 dpi from the pipeline's `_300dpi.png` (byte for byte), or else from the vector PDF via poppler. As for TIFF, a raster is never enlarged. |
| `r/figures.R` | `nf_rna_pheatmap_theme_text()`, applied only to the standard correlation heatmap. It restyles text only: the title row height follows pheatmap's own 1.5 × title-height rule. |
| `r/l1_analysis.R` | `formats = c("png", "png_raster", "tiff", "pdf")` on all five L1 figures. |
| Docs | `figure-standard.md`, `gui.md`, `CHANGELOG.md`; scope corrections in `release_notes_v1.7.0.md` and `release_readiness_v170.md`. |
| Tests | `tests/test_gui_ux_v170.py` (17 new tests), plus an additive L1 raster assertion in `tests/test_downstream.py`. |

Not changed:
- profiles and the pipeline raster dpi (300; Elsevier 500);
- the Arimo font policy: the standard preset never substitutes another font, and Arimo is metric-compatible with Arial;
- colours, group-to-colour mapping, sizes of the other figures, legend positions (pheatmap keeps its right-hand legends);
- the DEG heatmap, file names, the manifest schema, QC thresholds, run IDs, isolation and resume behaviour.

A 600 dpi LZW TIFF is provided by the GUI download (from the vector PDF), not written by the pipeline.

## 3. Performance after

### AppTest

Same method as §1.2, final code; median rerun in ms.

| Page | Before | After |
|---|---|---|
| New Project | 9,272 | 56 |
| Contrasts | 5,196 | 33 |
| Projects | 4,701 | 36 |
| Validate & Run | 2,657 | 71 |
| Project overview | 2,642 | 25 |
| Input & Metadata | 2,616 | 24 |
| Results | 31 | 46 (the page now validates for the workflow bar and reads the figure manifest) |

### Browser

Method:
- headless Chrome against two isolated servers on the same data: 8706 = `git archive HEAD`, 8705 = the new code;
- time from the action to the end of the rerun it triggers, from a MutationObserver on `data-test-script-state`, with a 600 ms quiet window;
- n = 5 per action; values are median (SD) in ms.

| Action | Before | After |
|---|---|---|
| Initial load of the Projects page (warm server) | 4,863 (9) | 254 (6) |
| Page switch (Project overview ↔ Validate & Run) | 2,700 (27) | 115 (761)¹ |
| Select the Human FASTQ project in the sidebar | 2,765 (7) | 179 (19) |
| Enter a FASTQ folder (Input & Metadata) | 2,710 (40) | 138 (30) |
| Metadata tool input (Input & Metadata) | 2,726 (10) | 153 (13) |
| Save metadata | 7,950 (1,298)² | 216 (15) |
| Validate and refresh the plan | 2,764 (8) | 197 (17) |
| Run-history *Refresh* | 130 (9) | 123 (10) |
| Results page load | 259 (12) | 265 (8) |
| Results section switch (n = 4) | 131 (17) | 135 (10) |
| New Project: change the backend | 9,455 (11) | 175 (12) |
| First load after the server starts (n = 1) | 5,376 | 6,858 |

Notes:
1. One of the five switches took 2.0 s: the first visit to Validate & Run in the session fills the cached runtime-prefix check.
2. *Save* reruns twice (save, then `st.rerun`), so each rerun paid the validation.

Interpretation:
- Every action that validates a project got 15–54 × faster, and every interaction now finishes in about 0.1–0.3 s.
- Actions that did not validate (Refresh, Results) are unchanged, as expected.
- The first load after a GUI start still hashes the references once (5–7 s, n = 1 each). This is not improved and is announced by the spinner.
- The remaining time is Streamlit's own rerun and render. Backend I/O is now 1–11 ms per call (§1.2); legitimately expensive work is the FASTQ copy and the run itself, both unchanged.

## 4. Scientific invariance

The L1 backend was run on the real 8-sample NGS1130901 inputs, read-only from the protected run's work directory, into `gui-ux-170/l1_compare/`:
- `before` = `git archive HEAD` scripts;
- `after` = the working tree;
- the same R prefix the protected run used.

| Check | Result |
|---|---|
| `normalized_counts.tsv`, `vst.tsv`/`.csv`, `source_counts.csv`, `pca_scores.tsv`, `pca_variance.tsv`, `sample_correlation.tsv`, `sample_distribution_qc.tsv`, `library_size_qc.tsv`, `backend_summary.json` (filtering, normalization, VST method, flags), `r_session_info.txt` | byte-identical before vs after, **and** byte-identical to the protected run's own outputs |
| Correlation heatmap clustering (`display_order`), annotation colours, cell size | identical |
| Standard preset: PCA and both boxplots (PNG, TIFF) | byte-identical; PDFs pixel-identical at 300 dpi |
| Standard preset: correlation heatmap | changed as intended (typography; height 92.4 → 91.4 mm) |
| Legacy preset: every existing PNG/TIFF | byte-identical; PDFs pixel-identical at 300 dpi; only the four `_300dpi.png` files are added |
| New `_300dpi.png` | 300 dpi, the same pixel size as the TIFF |
| PDFs | embed only Arimo (vector text) |
| TIFFs | RGB, LZW, 300 dpi |

Contrast definitions, the DESeq2 model, `l2_analysis.R` and QC policies were not touched.

## 5. Tests

| Suite | Result |
|---|---|
| Figure suites (`test_figures`, `_presets`, `_settings`, `_review_tools`, `_device_metrics`, `_migration`, `_remediation`, `_dotplots`) | 218 passed |
| GUI suites (`test_gui*.py`, including the new `test_gui_ux_v170.py`) | 105 passed (final code) |
| Full suite (locked runtime `0e721ff3…-200d98f8759eaf17`, `NF_RNA_REQUIRE_R=1`), final code | **1,195 passed, 2 skipped, 0 failed** (26 min 51 s). The 2 skips are the known environment skips: real SAMtools/featureCounts, real SJ project. An earlier run before the last edits had 2 failures: `test_l1_runs_for_two_samples_where_pc2_does_not_exist[legacy/standard]` asserted the exact L1 format set `{png, tiff, pdf}`, and the expected set was updated to include `png_raster` (still an exact set). |

One existing test (`test_results_page_works_without_poppler`) looked up the button label `TIFF 600 dpi`. The label was restored rather than the test changed.

**Browser checks, isolated GUI on port 8705.** Screenshots are in `gui-ux-170/shots/`.

The HISAT2 configuration path, without running anything (`hisat2_config_check.json`):
- New Project: HISAT2 + featureCounts, reverse, pretrimmed, Ensembl 116, paired design with `patient`.
- The formula preview reads `~ patient + condition`.
- `project.yaml` has the expected input, upstream, design, reference and scope.
- FASTQ preview, then copy.
- Metadata pasted, applied and saved; the pairing check showed 3/3 complete pairs.
- Validate & Run: the project is valid; *Run L1* is disabled with its reason until a case ID is entered, then enabled. It was not pressed.

The Salmon path:
- the Results page of the copied NGS1130901 L1 run, with grouped figures at a common scale and the Download menu (screenshots);
- the pairing check flags a missing Tumor sample and a blank `patient`.

The protected project, run, FASTQ, reference registry, indexes and the user's GUI (port 8501) were not touched.

## 6. Remaining issues

1. **First load after a GUI start** hashes the references once (5–7 s). A persistent digest cache would remove it, but it would need a trust decision (cross-process staleness), so it was not done.
2. **Heatmap legends sit on the right** (pheatmap), while PCA and the boxplots have them below. Moving them would need a different heatmap implementation, so it was left.
3. **The correlation heatmap has a title; PCA and the boxplots do not.** The typography is now the same, but whether L1 QC figures carry embedded titles is a figure-standard decision for the owner.
4. **The DEG heatmap (L2)** still uses pheatmap's 1.3 × title. It was not in scope, and a change would be an L2 figure change.
5. **PDF vs raster text layout** (Phase 7b finding) is unchanged: TIFF 600 dpi and PNG 300 dpi rendered from the PDF follow the PDF layout.
6. **Copying FASTQ still blocks the page:** the importer has no progress callback.
7. **Real-world UX is not yet validated by the user.**

## 7. HISAT2 E2E readiness (NGS1130901_hisat2)

The settings planned for the next manual test were checked in the browser on a disposable project (§5). The real project was not created:
- Human, paired-end, pretrimmed, reverse, HISAT2 + featureCounts, Ensembl 116 / GRCh38.p14;
- paired design on `patient`, formula `~ patient + condition`.

Every setting is selectable and shown before *Run*, the pairing check applies, and validation accepts the project. With the automated suites passing (§5), the GUI is ready for that manual test. The real-world UX itself is validated only by that test.
