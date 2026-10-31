# nf-rna v1.7.0 Phase 7a — Local Research Workspace

Branch `feature/v1.7.0`, created from `main` at `7a92bdd` (= tag `v1.6.1`,
annotated, tagger Brian Wu, 2026-10-31).

## 1. Starting point: what v1.6.1 already provides

Recorded before any change, so that Phase 7a extends the existing GUI and does
not duplicate it.

### GUI modules (`src/rnaseq/gui/`)

| Module | Role in v1.6.1 | Streamlit? |
|---|---|---|
| `config.py` | `gui.yaml` (`allowed_roots`, `project_parent`, `condarc`, `log_tail_lines`, `launcher`, `download_limit_mb`); `GuiConfig.resolve_allowed` (resolves symlinks, must lie inside an *available* root); `root_status` (roots under `/media`, `/mnt`, `/run/media` must be on a mounted filesystem); `ensure_data_layout` (creates `projects/`, `fastq/` only when mounted) | no |
| `services.py` | `discover_projects` (two levels below each root), `project_row`, `recent_runs`, `registered_references`, metadata read/merge/save through the project validator, `discover_fastq` (validator FASTQ naming rule), design view / contrast definition / `design_proposal` (the `rnaseq design` transaction), KEGG pin status, `reuse_candidates` (accepted by `preflight_reuse_upstream`), `run_preflight`, `scope_label` | no |
| `jobs.py` | single active GUI run (`active_run.lock`, `O_CREAT|O_EXCL`), launch records and launcher logs in `~/.cache/nf-rna/gui/`, systemd or detached-process launcher, `run_display_state` (RUNNING without process → INCOMPLETE, nothing written), `tail` | no |
| `results.py` | read-only run outputs: `result_runs`, QC view (incl. F-1 strandedness from the frozen contract), DEG view (incl. F-3 `r_warnings`, shrinkage), enrichment modules/result sets/members, provenance (incl. F-2 threads), report path, downloads, verified delivery zip | no |
| `pages.py` | 7 pages: Dashboard, New Project, Input & Metadata, Contrasts, Run, Results, Reports & Downloads; sidebar project selector (kept in the URL) and "Open project by path" | yes |

### Project, run and result data the GUI reads

| Source | Read through |
|---|---|
| `project.yaml` | `rnaseq.project.load_project` (pydantic `ProjectConfig`) |
| `metadata.csv`, `contrasts.csv`, input | `rnaseq.validators.validate_project` (`ValidationReport`) |
| `runs/<case>/<run>/run_state.json`, `logs/rnaseq.process.json` | `rnaseq.run_status.discover_runs`, `effective_state`, `build_status`, `render_status` |
| `frozen/downstream_contract.json` (analysis level, `upstream_qc`), `frozen/project.yaml`, `frozen/contrasts.csv` | `results.ResultRun` |
| `provenance/run_provenance.yaml` | `results.provenance` |
| `downstream/l1`, `downstream/l2`, `delivery/` (report, manifest) | `results.*` |
| Reference registry | `rnaseq.references.registered_local_references` |

### Gaps against the Phase 7a brief

| Brief item | v1.6.1 state |
|---|---|
| Dataset browser | Only the `path_picker` inside Input & Metadata (one folder at a time, no sizes, no file kinds, no readability, no symlink reporting). FASTQ recognition exists (`discover_fastq`). |
| Assign a dataset to a project | Exists per kind on Input & Metadata (`import_count_matrix`, `import_fastq_samplesheet`, metadata paste/import). No size/free-space information before a FASTQ copy. Origin of imported FASTQ only in `planning/imported_fastq_samplesheet.csv`, not shown. |
| Project overview | Dashboard table: id, scope, input type, latest run, state. No species, backend, reference, sample count, contrasts, QC summary, next step. |
| Project switching | Sidebar selector. The run selectors of Run, Results and Reports use global widget keys (`run_selected`, `results_run`, `reports_run`, `run_log`), so a selection can leak between projects. |
| Run history | Per project on the Run page; no QC column; no link from a run to its results or report. |
| Execution vs scientific QC | Execution state shown; F-1/F-3 shown only deep inside Results. No per-run QC roll-up. |
| Upstream reuse hint | Only accepted candidates listed; refused runs and the refusal reason are not shown. |
| Validation | Raw `rnaseq validate` text only after pressing *Validate*. |
| Log viewer | Last N lines; no warning/error filter. |

No new analysis engine, project schema, database, job queue, login, REST API
or service is needed for any of the gaps.

## 2. What Phase 7a adds or changes

All changes are in the GUI package (`src/rnaseq/gui/`), its tests, a
development script and documentation. No CLI command, project schema, run
format, workflow, R script or lock changed (§7).

| Area | Change | Files |
|---|---|---|
| Navigation | Sidebar groups *Workspace* (Projects, Datasets, Project overview), *Set up* (New Project, Input & Metadata, Contrasts), *Analysis* (Validate & Run, Results, Reports & Downloads); under the project selector, a one-line summary (scope, input, samples, latest execution and QC). | `pages.py` |
| Projects (was Dashboard) | One row per project: scope, samples, contrasts, last run, execution, scientific QC, runs, species, input, quantification, reference, path. *Open* goes to Project overview. Unavailable data folders are reported first; storage, doctor checks and references moved below. | `pages.py`, `workspace.project_summary` |
| Datasets (new) | Browse a folder in the allowed roots; files, total size, files per kind, unreadable files, ignored out-of-root symlinks; recognition with the validator's FASTQ and count-matrix rules and the samplesheet reader; *use in project* through the existing importers, with a copy-size / free-space check and an allowed-root check of every FASTQ source. | `datasets.py`, `pages.py` |
| Project overview (new) | Summary, analysis settings, workflow row (Dataset → … → Results) with the next step and its findings, data locations (project input, original FASTQ folders), analysis history with execution and QC per run, QC checks of the chosen run, *Open results* / *Open report & downloads* / *Run details and logs*. | `workspace.py`, `pages.py` |
| Validate & Run | Validator findings always shown, grouped by workflow step; settings as a table with what the scope produces; refused reuse candidates with the preflight's reason; *Run* disabled while invalid; run table with execution and QC; QC checks of the selected run; links to Results/Reports; status dashboard in an expander; log viewer starts at `rnaseq.log` and can show only warning/error lines. | `pages.py`, `workspace.reuse_options`, `workspace.alert_lines` |
| Results | Scientific QC summary at the top of *QC & exploratory*. Everything else unchanged. | `pages.py` |
| Input & Metadata | FASTQ imports go through the same guarded helper as Datasets (size / free space, allowed roots); the original FASTQ folder is shown. | `pages.py` |
| Per-project state | Run selectors, case ID, result sections/contrasts and report selectors are keyed per project or run directory; the run choice and case ID survive switching projects. A project opened by URL survives page switches. | `pages.py` |

### Fixed along the way (GUI only)

- **Project lost after a page switch** when the project had been opened by URL
  (bookmark, reload): `st.switch_page` and page links drop the query string,
  and the project lived only there. `current_project` now copies it into the
  session and restores the URL. Found by the browser acceptance.
- **Selections leaking between projects**: global widget keys (`run_selected`,
  `results_run`, `reports_run`, `run_log`, `rep::<group>`) are now per project
  or per run directory.

## 3. Existing services used

| Need | Existing service |
|---|---|
| Path guard, mount check | `GuiConfig.resolve_allowed`, `root_status`, `available_roots` |
| Project files | `project.load_project`, `validators.validate_project` |
| FASTQ recognition | `validators.validate_fastq` (via `datasets.recognise_fastq`) |
| Count matrix recognition | `validators.validate_counts` |
| Samplesheet recognition | `project._read_fastq_samplesheet` (the `rnaseq new --fastq-samplesheet` reader) |
| Import | `project.import_fastq_samplesheet`, `project.import_count_matrix`, `project.fastq_samplesheet_text`; metadata through `services.save_metadata` |
| Runs and state | `run_status.discover_runs`, `jobs.run_display_state` (→ `effective_state`), `build_status`, `render_status` |
| QC records | frozen `downstream_contract.json` `upstream_qc` (F-1), provenance `strandedness_qc`, `l1/backend_summary.json` `flags`, `contrast_metadata.json` `r_warnings` and `shrinkage` (F-3) |
| Reuse | `service.preflight_reuse_upstream` |
| Launch | `services.run_preflight`, `jobs.launch_run` (single active run lock unchanged) |
| Results / reports | `results.*` unchanged |

## 4. Data flow

```
Datasets page                       Project files (authority)             Runs (authority)
-------------                       -------------------------             ----------------
folder in allowed root ─scan──▶ kinds, sizes, permissions
   │ recognise (validator rules)
   ▼
FASTQ / samplesheet ──import_fastq_samplesheet (copy)──▶ input/fastq/ + planning/imported_fastq_samplesheet.csv
count matrix ───────import_count_matrix (copy)────────▶ input/counts.csv
metadata ──────────editor (unsaved) ─save_metadata───▶ metadata.csv
                                    contrasts page ──▶ contrasts.csv, project.yaml (rnaseq design)
                                    project.yaml, metadata.csv,
                                    contrasts.csv ──validate_project──▶ workflow steps, findings
                                                     Run ──rnaseq run──▶ runs/<case>/<run>/
                                                                         run_state.json ─▶ execution status
                                                                         frozen contract, l1/l2 records ─▶ scientific QC
                                                                         downstream/, delivery/ ─▶ Results, Reports
```

The GUI writes only through the existing importers, `save_metadata`, the
`rnaseq design` transaction, the KEGG pin and `rnaseq run`, exactly as in
v1.6.1. The new pages and services only read.

## 5. Acceptance

### Automated tests

| Suite | Result |
|---|---|
| `tests/test_gui_workspace.py` (new, 13 tests: dataset kinds/sizes/permissions, symlink escape, copy plan, execution states, F-1/F-3 QC roll-up with no file written, run history, project summary, workflow steps, Datasets page, Datasets → project FASTQ import by copy with originals unchanged, project switching keeps state, Run page validation and per-project case ID, project opened by URL) | 13 passed |
| `tests/test_gui.py`, `tests/test_gui_results.py` (v1.6.x; two results tests updated for the per-project widget keys) | 22 passed |
| Full suite (`pytest`, locked downstream R prefix on PATH, on the tree with every change except the final project-table column order and the Datasets import AppTest, both GUI-only and covered by the 35 GUI tests afterwards) | **1126 passed, 2 skipped, 0 failed** in 26 min (v1.6.1: 1114 passed, 2 skipped); the skips are the same two (real SAMtools/featureCounts, real SJ project) |

### Browser acceptance on the server

The GUI ran from this checkout with `rnaseq gui --port 8701` on the server
(Streamlit 1.66.1, bound to `127.0.0.1` only, verified with `ss -ltnp`),
with an acceptance `gui.yaml` that adds two v1.6.1 validation folders as
data roots: `/home/brian/nf-rna-validation/v161-reliability/nfrna` and
`/home/brian/nf-rna-validation/v161-phase3/data`. The user's
`~/.config/nf-rna/gui.yaml` and the GUI already running on port 8501 were not
touched. `scripts/dev/gui_workspace_acceptance.py` drove headless Google
Chrome: **27/27 checks passed**
(`docs/dev/screenshots/7a/browser_checks.json`).

| # | Brief criterion | Evidence | Result |
|---|---|---|---|
| 1 | Datasets in allowed roots are shown | D5 FASTQ folder: 14 files, 85.7 MB, 14 FASTQ, 0 unreadable (W2); airway folder (W3); the empty `fastq/` on the mounted data disk | PASS |
| 2 | FASTQ / raw counts / metadata recognised | D5 paired-end samples listed by the validator rule, the single-end `D5_SE_*` files reported as unpaired (W2); airway `counts.csv` 63677 genes × 8 samples by `validate_counts` (W3); airway metadata preview | PASS |
| 3 | Existing CLI and GUI projects open | `airway-gui-demo` (GUI, W4), `airway-e2e` (CLI, W5), `rel_d5_neg`, `rel_d3_counts` (CLI) | PASS |
| 4 | Switch projects, correct UI state | Run GUI-KILL-001 chosen in airway-gui-demo → sidebar switch to rel_d5_neg (W6) → back: GUI-KILL-001 still chosen, URL follows (W7); AppTest also covers the case ID | PASS |
| 5 | Sample count, contrasts, history | Projects table (W1), Project overview (W4, W5, W8) | PASS |
| 6 | SUCCESS / FAILED / INCOMPLETE read | GUI-6B-L2 SUCCESS (W15), GUI-KILL-001 INCOMPLETE (W14), FAILED fixture with its recorded error (W12); log filter (W13) | PASS |
| 7 | F-1 and F-3 warnings shown | REL_D5_NEG: Scientific QC WARNING, "Reads contradict the declared `forward` for: D5_R_A, D5_R_B", execution SUCCESS (W8); Results strandedness WARNING (W9). REL_D3_COUNTS: "WARNINGS: 6", apeglm line-search step (W10); Results R-warnings expander with 6 (W11) | PASS |
| 8 | Results / Reports on airway runs | *Open results* from the overview opens GUI-L2-001 with all figures rendered (W16); *Open report & downloads* opens GUI-6B-L2 with its report download (W17) | PASS |
| 9 | No scientific table, figure or provenance modified | SHA-256 of all 4048 files of airway-gui-demo, airway-e2e, rel_d5_neg, rel_d3_counts and rel_d4_sal before and after the browser acceptance: identical. No new folder on the data disk. | PASS |
| 10 | GUI starts on the server; Mac/Windows via SSH tunnel | Started on the server and driven over HTTP on 127.0.0.1:8701; listens on loopback only, so the documented `ssh -N -L 8501:127.0.0.1:8501 user@server` is the access path. **The tunnel itself was not exercised**: no key-based SSH login was available in this session (also not exercised in Phase 6). | PARTIAL |

The FAILED case uses one synthetic fixture
(`/home/brian/nf-rna-work/testdata/gui7a-fixtures/failed-demo`, a copy of the
airway-gui-demo project files and a FAILED `run_state.json` with a short log);
no real FAILED run exists and no pipeline was run. The fixture was deleted
after acceptance, because `testdata/` is one of the default allowed roots and
it would otherwise appear in the everyday GUI. Evidence (logs, hashes, all
screenshots, full-suite log) is in `/home/brian/nf-rna-validation/gui-7a/`.

### Screenshots (`docs/dev/screenshots/7a/`)

| File | Shows |
|---|---|
| `W1_projects.png` | Projects list and recent runs with execution and QC |
| `W2_datasets_fastq_d5.png`, `W3_datasets_airway_counts.png` | Datasets: FASTQ folder, count matrix recognition |
| `W4_project_gui_demo.png`, `W5_project_cli_e2e.png` | Project overview of a GUI and a CLI project |
| `W6_switched_to_rel_d5_neg.png`, `W7_switched_back_choice_kept.png` | Project switching |
| `W8_qc_f1_rel_d5_neg.png`, `W9_results_f1_strandedness.png` | F-1 strandedness WARNING |
| `W10_qc_f3_rel_d3_counts.png`, `W11_results_f3_warnings.png` | F-3 R warnings |
| `W12_run_failed_fixture.png`, `W13_run_log_filter.png`, `W14_run_incomplete.png`, `W15_run_success_settings.png` | Validate & Run: states, settings, logs |
| `W16_results_from_project.png`, `W17_reports_from_project.png` | Results and Reports opened from the overview |

## 6. Scientific Engine unchanged

- `git diff --name-only v1.6.1 -- src/rnaseq/r src/rnaseq/workflows workflow 'src/rnaseq/*.py'`
  is empty: no R script, Nextflow file, lock, model, validator, planner,
  service, execution or CLI module changed. Source changes are limited to
  `src/rnaseq/gui/` (`pages.py` changed; `datasets.py`, `workspace.py` new).
- The GUI recomputes no statistic: the QC roll-up only reads recorded
  verdicts; tables and figures are still read by the unchanged `results.py`.
- F-1/F-2/F-3 records are read as written; the strandedness, R-warning and
  provenance views of v1.6.1 are unchanged (W9, W11).
- `run_state.json` is never written by the new code (tests compare file
  digests before and after; the browser acceptance compares 4048 SHA-256).
- No new run was started; the downstream runtime and lock are untouched.

No scientific problem was found in this phase.

## 7. Not done or not verified

- **SSH tunnel from a Mac/Windows client** not exercised (see criterion 10).
- **No "use in place" FASTQ mode.** Assignment copies into the project
  (existing importer; `input.path` must lie inside the project). The GUI shows
  size and free space first. A reference-in-place mode would change the input
  contract and is a separate decision.
- **Origin of an imported count matrix** is not recorded by the importer, so
  it is not shown.
- **A real FASTQ import into a project in the browser** was not done (it
  would copy data onto the data disk); it is covered by an AppTest with small
  files (originals byte-identical afterwards) and by the v1.6.x importer tests.
- **No run was started from the GUI in this phase** (the brief forbids
  re-running FASTQ pipelines; the launch path is unchanged from v1.6.1).
- **`/dashboard` is not a URL of the Projects page** (Streamlit serves the
  default page at `/`); unchanged from v1.6.x, documented.
- **Scope switching** (running an L2 project as L1) remains a project.yaml
  edit.
- Version stays 1.6.1; CHANGELOG and release notes for v1.7.0 belong to the
  release step after Phase 7b.
