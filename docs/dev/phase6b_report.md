# nf-rna v1.6.0 Phase 6b — results, reports, downloads and release preparation

Branch `feature/v1.6.0`, on top of the Phase 6a commits `7cdb329` and
`ed4c4ca`. Phase 6b code commit: `80f1d08`. Nothing was pushed, tagged or
released.

## Current state before Phase 6b (code as of `ed4c4ca`)

Recorded from the code and from the airway runs on this machine.

### What a run writes

`runs/<case>/<run>/` of a raw-count L2 run (`airway-gui-demo`, GUI-L2-002):

| Path | Content |
|---|---|
| `run_state.json` | state, phase, start/end, `delivery` path |
| `frozen/` | `downstream_contract.json` (`analysis_level`, `analysis.enrichment`, annotation, design, source), `project.yaml`, `metadata.csv`, `contrasts.csv`, input and execution manifests |
| `downstream_inputs/` | the staged copies the Nextflow processes read (`execution_inputs.json`, `execution_manifest.yaml` with the runtime identity: lock, wheel, nf-rna version, source revision) |
| `downstream/l1/` | `backend_summary.json` (filter, normalization, VST, package versions), `library_size_qc.tsv`, `sample_distribution_qc.tsv`, `pca_scores.tsv`, `pca_variance.tsv`, `sample_correlation.tsv`, `normalized_counts.tsv`, `vst.tsv/.csv`; figures `pca`, `pca_by_batch` (with a batch), `sample_correlation`, `boxplot_log2_raw_counts`, `boxplot_vst` as `.png` + `.pdf` + `.tiff`; `figure_manifest.json`, `scientific_provenance.json`, `r_session_info.txt` |
| `downstream/l2/` | `normalized_counts.tsv`; per contrast `contrasts/<id>/`: `all_genes.tsv`, `significant.tsv`, `upregulated.tsv`, `downregulated.tsv`, `backend_summary.json`, `contrast_metadata.json` (DEG definition, shrinkage method/status/error, annotation, software), `volcano`, `ma_plot`, `heatmap` (`.png`, `_300dpi.png`, `.pdf`, `.tiff`), heatmap tables |
| `downstream/l2/enrichment/go/` | `go_backend_summary.json` (significance semantics, ORA universe per contrast × foreground × ontology, figure records); `<contrast>/<significant|up|down>/GO_<ont>.tsv` (significant), `GO_<ont>_all_terms.tsv`, `GO_<ont>_members.tsv` (v1.5, significant pathways), `GO_<ont>_gene_members.tsv` (legacy layout), `GO_<ont>_simplified.tsv`, `dotplot_<ont>[_simplified].png/.pdf` |
| `…/enrichment/kegg/` | `kegg_backend_summary.json` (snapshot ID, release, SHA-256, selection), `kegg_snapshot_manifest.json`; `<contrast>/<foreground>/kegg.tsv`, `kegg_all_terms.tsv`, `kegg_members.tsv`, `kegg_gene_members.tsv`, `dotplot.png/.pdf` |
| `…/enrichment/gsea_go/` | `gsea_backend_summary.json` (ranking metric, seed, nPermSimple, eps; per ontology counts and figure records); `<contrast>/<BP|MF|CC>/all_terms.tsv`, `significant.tsv`, `members.tsv`, `gsea_core_members.tsv`, `dotplot[_simplified]`, `running_score/<sign>_<rank>_<ID>.png/.pdf` |
| `…/enrichment/gsea_kegg/` | as GO GSEA, one level less (no ontology), `core_members.tsv` |
| `downstream/report/` | `report.html`, `enrichment_status.tsv` |
| `delivery/` | built by `service.assemble_delivery`: `report_<date>.html`, `README_<date>.md`, `counts/` (raw counts and VST as CSV, `artifact_manifest.json`), `tables/` (every TSV/CSV of `downstream/` at the same relative path), `figures/png|pdf|tiff_300dpi/` (every figure), `methods_and_versions/` (provenance JSON, R session info, figure manifests, contrast metadata, frozen manifests, run state), `multiqc/` for FASTQ runs, and `delivery_manifest.yaml` (SHA-256 and size of every file). 414 files, 90 MB for airway L2. **No archive is written**. |
| `provenance/run_provenance.yaml` | nf-rna version, git commit, Nextflow and Python versions, lock and runtime prefix, upstream implementation |

An L1 run has `downstream/l1`, `downstream/report` and the L1 part of
`delivery/`; no `l2/`.

### How report.html is produced

The Nextflow process `TECHNICAL_REPORT` (`TECHNICAL_REPORT_NO_ENRICHMENT`,
`TECHNICAL_REPORT_L1`) runs `python -m rnaseq.workflow_support report` with the
frozen contract, the staged inputs, the L1/L2 directories and the enrichment
directories. `workflow_support.report` (L2) and `_report_l1_only` write one
HTML file. Figures are embedded as base64 PNG (`_image_html`), there is no
external stylesheet, script or font, so the file is already self-contained.
Missing figures are written as "unavailable — reason" from the backend's
`figures` record. `service.assemble_delivery` copies it as
`delivery/report_<date>.html`.

Sections present before 6b: overview (case, run, input type, organism, number
of samples, design, contrasts), experimental design, L1 QC with PCA, sample
correlation and boxplots, L2 DEG per contrast (counts, thresholds, volcano, MA,
heatmap), enrichment status table, ORA (significance rule, universe, GO
simplified dotplots), GSEA (ranking, dotplots, running-score plots, top
terms), methods and reproducibility (R and DESeq2 versions, source revision).
Missing from the required list: a sample table, reference and quantification,
the LFC shrinkage method/status, and the software versions/provenance beyond R
and DESeq2 (lock hash, apeglm/ashr, clusterProfiler, OrgDb, GO.db, KEGG
snapshot).

### GUI before 6b

The Run page showed the path of a finished run's `delivery/` folder; there was
no result view or download.

## Changed files

| File | Change |
|---|---|
| `src/rnaseq/gui/results.py` (new) | Read-only access to a run's outputs (no Streamlit, no pandas): `result_runs`, `Table`/`read_table` (TSV as strings), `qc_view`, `contrast_ids`, `deg_view`, `significant_ids`, `filter_deg`, `enrichment_modules` (facts, legacy flag), `enrichment_results`, `running_score_figures`, `pathway_members`, `provenance`, `downloads`, `report_path`, `delivery_size`, `build_delivery_archive`, `scp_command`, `too_large` |
| `src/rnaseq/gui/pages.py` | Results page (QC & exploratory, Differential expression, Enrichment, Provenance; L2 sections only for L2 runs), Reports & Downloads page, "Open in Results" link on a successful run |
| `src/rnaseq/gui/config.py` | `download_limit_mb` (default 500) |
| `src/rnaseq/gui/__init__.py` | install hint uses the Git URL (nf-rna is not on PyPI) |
| `src/rnaseq/workflow_support.py` | report sections: Samples, Reference and quantification, DEG rule and LFC shrinkage per contrast, Software versions and provenance (both L1-only and L2 reports) |
| `src/rnaseq/__init__.py`, `CITATION.cff` | version 1.6.0, no `date-released` |
| `pyproject.toml` | `gui = ["streamlit>=1.66,<2"]` |
| `CHANGELOG.md` | `1.6.0 — Unreleased`, opening with the v1.5.1 inclusion and the v1.5.0 semantics warning |
| `docs/dev/release_notes_v1.6.0.md` (new) | GitHub Release text |
| `docs/dev/release_checklist.md` | rewritten for v1.6.0 |
| `README.md`, `README_zh-TW.md`, `docs/quickstart.md`, `docs/runtime.md`, `docs/gui.md` | v1.6.0, GUI results/downloads, Git-URL install of the `gui` extra |
| `docs/dev/DECISIONS.md` | Phase 6b table |
| `tests/test_gui_results.py` (new) | 8 tests on a synthetic run (5 service, 3 AppTest) |
| `tests/test_workflow_support.py` | new report sections; FASTQ reference-section branch |
| `scripts/dev/gui_results_acceptance.py` (new) | Playwright acceptance of the new pages |
| `docs/dev/screenshots/6b/` | 27 screenshots |

No R script, Nextflow workflow, lock, reference, existing run or MRPan file
was changed. `git diff ed4c4ca -- src/rnaseq/r workflow` is empty.

## Acceptance evidence

Evidence files are in `/home/brian/nf-rna-validation/gui-6b/` (outside the
repository). The GUI ran on 127.0.0.1:8501 from the committed tree `80f1d08`;
browser checks used headless Google Chrome via Playwright.

### Both projects can be browsed

`airway-gui-demo` (GUI-built, run GUI-6B-L2) and `testdata/airway-e2e`
(CLI-built, AIRWAY-E2E) were browsed through QC, DEG, all four enrichment
analyses and provenance (screenshots `F1`–`F7` and `J1`–`J7`). Each step checks
that no Streamlit exception occurred.

### DEG table: DUSP1

Searching `DUSP1` on the DEG page; the row read from the rendered grid
(accessibility table of the grid) against `all_genes.tsv`:

| | grid | `all_genes.tsv` |
|---|---|---|
| gene_id | ENSG00000120129 | ENSG00000120129 |
| log2FoldChange | 2.94784960588663 | 2.94784960588663 |
| padj | 1.82328987063955e-125 | 1.82328987063955e-125 |
| log2FoldChange_shrunk | 2.93659125491481 | 2.93659125491481 |

Identical strings for both projects (`evidence-airway-gui-demo.json`,
`evidence-airway-e2e.json`). The page states "Shrinkage method: apeglm …
status success" and that DEGs are called on the raw log2FoldChange (screenshot
`F2`).

### KEGG pathway members

KEGG ORA, all significant DEGs: clicking the first row of the significant
terms table (hsa04060, Cytokine-cytokine receptor interaction) shows its member
table. All 29 rendered rows equal, cell by cell and in order, the 29 rows of
`kegg_members.tsv` for hsa04060 (both projects; screenshot `F3`).

### Figures render

Checked in the browser (`naturalWidth > 0`): volcano, MA plot and heatmap;
simplified GO ORA dotplots; 10 GO BP GSEA running-score plots (5 per NES sign);
KEGG GSEA dotplot and running-score plots (screenshots `F2`, `F4`–`F6`).

### Unavailable outputs and L1-only

Truncated copies of the GUI-L2-002 run, inside the allowed root
(`~/nf-rna-work/testdata/gui6b-fixtures/`, deleted after the check):

- `missing-outputs`: `ma_plot.*`, `contrast_metadata.json`, the KEGG
  significant `dotplot.png` and `gsea_backend_summary.json` removed. Shown: "MA
  plot: trt_vs_untrt: unavailable — ma_plot.png was not produced (L2 results
  from before nf-rna v1.5)", "Shrinkage: not recorded …", "KEGG dotplot (all
  significant DEGs): unavailable — …", "GO GSEA: unavailable —
  gsea_backend_summary.json was not produced" (`I1`–`I3`).
- `l1-only`: contract level L1, no `l2/`. Results offers only *QC &
  exploratory* and *Provenance*; Reports & Downloads lists no DEG or
  enrichment group, and "Delivery package: unavailable" (no delivery in the
  copy) (`I4`, `I5`).

The same cases, plus a pre-v1.5 (legacy) enrichment summary and a failed
shrinkage, are covered by `tests/test_gui_results.py` without R.

### Downloads: SHA-256

Each file was downloaded through the browser and compared with the server file:

| File | Equal |
|---|---|
| `report_20261030.html` (delivered report) | yes |
| `l2/normalized_counts.tsv` | yes |
| `trt_vs_untrt/all_genes.tsv` | yes |
| `trt_vs_untrt/significant.tsv` | yes |
| `kegg/…/significant/kegg_all_terms.tsv` | yes |
| `kegg/…/significant/kegg_members.tsv` | yes |
| `ma_plot.png` | yes |
| `ma_plot.pdf` | yes |
| delivery package zip | yes (download = `~/.cache/nf-rna/gui/packages/…zip`) |

All 414 members of the zip have the SHA-256 of the corresponding
`delivery/` file (`G1`, `G2`; `evidence-airway-gui-demo.json`).

### report.html offline

The delivered report of GUI-6B-L2 opened from `file://` in headless Chrome with
the browser context offline: 33 images, all rendered, no request other than
`file:`/`data:` (`H1`). Sections (h2): Run and analysis overview, **Samples**,
**Reference and quantification**, Experimental design, L1 — quality control,
L2 — differential expression (with **DEG rule and LFC shrinkage**), L2 —
functional enrichment status, L2 — over-representation analysis, L2 —
preranked GSEA, Methods and reproducibility, **Software versions and
provenance**.

### Airway rerun: results unchanged

`airway-gui-demo` was validated and run again from the GUI Run page
(case GUI-6B-L2, process launcher; 21:01:52–21:10:45 including the creation of
the new runtime prefix `nf-rna-downstream-linux-64-0e721ff3cf1983f6-8b8770dc88641846`
for the 1.6.0 wheel). Recorded commit `80f1d08`, nf-rna 1.6.0, same lock
`0e721ff3…`, pinned KEGG snapshot `hsa_kegg20261028_20261030T061908Z`.

- Tables: 135 of 135 `.tsv`/`.csv` under `downstream/` byte-identical to
  `airway-e2e` (`compare_tables_e2e_vs_6b.txt`).
- Figures against the pre-6b GUI run GUI-L2-002: 142/142 PNG and 35/35 TIFF
  byte-identical; 73/73 PDF identical in every content and object stream once
  the Cairo `CreationDate` and the cross-reference byte offsets that shift with
  it are excluded (`compare_figures_6a_vs_6b.txt`). The report change did not
  change any table or figure; `enrichment_status.tsv` is identical.

### Test suites

- With the locked R runtime of this build on PATH (`…-8b8770dc88641846/bin`):
  **1080 passed, 2 skipped** in 25 min (`pytest-full/pytest-full.log`; the two
  skips need real SAMtools/featureCounts and a real SJ project, as in 6a). The
  two FASTQ reference-section report tests were added after this run started
  and pass on their own (`tests/test_workflow_support.py -k reference_section`).
- Without the GUI extra (Python venv without Streamlit or pandas, system R):
  1018 passed, 50 skipped, 4 failed, 10 errors. The 6 GUI page tests skip;
  the 5 service tests of `test_gui_results.py` run. The 14 failures/errors are
  the font/figure tests that need R packages missing from the system R
  (`test_figure_presets`, `test_figure_migration`, `test_figure_remediation`),
  the same set as the Phase 6a run of this environment (which was identical on
  `main`); the 6a local-path failure is gone.
- `tests/test_release_metadata.py` passes with `1.6.0 — Unreleased` and no
  `date-released`.

### Phase 6a operations

The rerun itself went through the 6a Run page (validate, plan refresh,
background launch, single active run). The 6a page walk (`pages` step of
`gui_browser_acceptance.py`: Dashboard, New Project, Input & Metadata,
Contrasts, Run) passed on the 6b build (`K_6a_*` screenshots), and
`tests/test_gui.py` passes unchanged.

## Not done, not verified, known limitations

- No real FASTQ run was browsed (no FASTQ pipeline was run). The MultiQC
  download and the FASTQ branch of the report's Reference and quantification
  section are covered by unit tests only.
- Legacy (pre-v1.5) enrichment display was verified only with the synthetic
  fixture; no pre-v1.5 run exists on this machine.
- The grid comparison reads the cells rendered in the browser. The DEG grid
  shows the key columns first (gene_id, symbol, log2FoldChange,
  log2FoldChange_shrunk, padj, pvalue, baseMean) so they are on screen.
- A downloaded file is read into the Streamlit server's memory when its button
  is clicked (bounded by `download_limit_mb`). Larger files show a path and an
  `scp` command instead; the command uses the server's hostname and current
  user, which may need adjusting for an SSH alias.
- Delivery zips stay in `~/.cache/nf-rna/gui/packages/` (one per run, about the
  size of `delivery/`); nothing removes them automatically.
- `report.html` is shown in the page in an iframe (`components.html`); a large
  report (2–3 MB) takes a moment to appear. Download it to read it fully.
- The release itself (finalize date, merge, push, CI, tag, GitHub Release,
  clean-environment install with and without `.[gui]`) is left to the
  maintainer: `docs/dev/release_checklist.md`.

## Commands

GUI on the server (as in 6a):

```bash
tmux new -s nf-rna-gui
cd /home/brian/projects/nf-rna
/home/brian/miniforge3/envs/nf-rna/bin/rnaseq gui
```

From the Mac: `ssh -N -L 8501:127.0.0.1:8501 brian@<server>`, then open
http://127.0.0.1:8501 → Results / Reports & Downloads.

Acceptance (GUI running on 8501):

```bash
python scripts/dev/gui_results_acceptance.py browse downloads report --shots DIR
python scripts/dev/gui_results_acceptance.py browse --project /home/brian/nf-rna-work/testdata/airway-e2e --shots DIR --prefix J
python scripts/dev/gui_results_acceptance.py unavailable --shots DIR --missing <copy> --l1-only <copy>
```
