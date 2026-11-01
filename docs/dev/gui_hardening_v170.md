# nf-rna v1.7.0 — pre-release GUI hardening

Branch `feature/v1.7.0`, starting from `ac6471d` (v1.7.0 finalization, not yet
tagged). Version stays **1.7.0**; release dates are unchanged. Only
`src/rnaseq/gui/` (pages, services, results), GUI tests, a development
acceptance script and documentation change. The scientific engine (DESeq2,
Salmon, HISAT2/featureCounts, tximport, GO/KEGG, GSEA, thresholds, R scripts,
reference management, Nextflow workflows, Conda locks) is untouched:

```
git diff ac6471d -- src ':!src/rnaseq/gui'   → empty
git diff ac6471d -- workflow pyproject.toml environment.yml   → empty
```

The usability review behind this work was made on v1.6.1. Every reported
item was therefore checked again on the v1.7.0 code (`ac6471d`) before
anything was changed.

## 1. Finding matrix

| # | Reported item | v1.7.0 verdict | How it was established |
|---|---|---|---|
| P0-1 | *Use this folder* / *Use this file* crash | **CONFIRMED** | AppTest on `ac6471d`: `StreamlitWidgetAlreadyInstantiatedError` for `ds_text` and `counts_text`; browser on `ac6471d` (port 8704): exception shown (`original_behaviour.json`, `O1_original_use_this_folder.png`) |
| P0-1b | Path picker follows a symlinked folder out of the allowed roots | **CONFIRMED** (code) | *Open folder* set the folder to `cwd/<name>` without resolving it, and *Subfolders* listed symlinked folders, so a link to an outside folder could be opened and its contents listed. Typed paths were already checked (`resolve_allowed`). |
| P0-2a | Bulk assignment discards unsaved grid edits | **CONFIRMED** | Browser on `ac6471d`: donors P1–P3 typed into the grid, then *Assign values to samples* → donors empty (`original_behaviour.json`, `O2_original_bulk_lost_edits.png`) |
| P0-2b | *Paste table* / *Import CSV* / Datasets *Load into the metadata editor* replace the table | **CONFIRMED** (code) | Each built a new table from the pasted rows only: values of samples or columns not in the paste, and every unsaved edit, were dropped without a warning |
| P0-2c | Unsaved grid edits lost when switching page | **CONFIRMED** (code) | Grid edits lived only in the data editor's widget state, which Streamlit drops when the widget is not drawn; the docs claimed they were kept |
| P0-2d | No indication of unsaved metadata | **CONFIRMED** | Only a "Saved" message, which stayed visible after later edits |
| P0-2e | *Reload* left the old editor's edits applied to the reloaded table | **CONFIRMED** (code), not reproduced | The editor key was not renewed on reload |
| R-1 | Run status needs manual *Refresh* | **CONFIRMED** | Documented as "no automatic polling" |
| R-2 | No way to cancel a run | **CONFIRMED** — **not implemented** (§3.4) | Process ownership inspected; cancellation not safe in every phase |
| D-1 | *Run* disabled without a reason | **ALREADY_FIXED** in Phase 7a for active run, invalid project and missing case ID; an invalid case ID was still unexplained (**CONFIRMED** for that case) | `runner_page` captions |
| D-2 | *Save contrasts* / *Add* disabled without a reason | **CONFIRMED** | No caption or help |
| D-3 | Case IDs should be invented automatically | **NOT CHANGED** by decision | A case ID names the analysis; a generated name would be a silent choice. The placeholder suggests a form. |
| F-1 | Reference *Other*: fields hidden, error before input, CLI flags in the message | **CONFIRMED** (code) | Fields were inside the collapsed *Advanced settings*; with the default *local* source and an empty root, `reference_options` raised "Managed local reference requires --reference-root." on first display |
| T-1 | Results tables sort lexically | **CONFIRMED** | `pages.frame` built every column with `dtype=str`, so a column-header click sorted text (`-0.0001…` before `-5.3`) |
| I-1 | FASTQ import has no progress | **CONFIRMED**; the completion message was never visible either (an immediate `st.rerun()` replaced it) | `_import_samplesheet` |
| S-1 | `?project=` opens any local folder with a project.yaml | **CONFIRMED** | AppTest on `ac6471d`: a project outside the roots opened (with only a warning), and every page, including saving metadata and starting a run, could act on it |
| S-2 | Stale *Full validation report* | **CONFIRMED** (code) | The report was kept in the session and shown after project.yaml, metadata or contrasts changed |
| S-3 | Contrasts page not refreshed after saving | **CONFIRMED** (code) | The saved-contrasts table was drawn before the save; composed contrasts of levels removed from the metadata were kept |
| N-1 | Project navigation | **NOT_REPRODUCED** | Sidebar switching, *Open* buttons and URL reloads keep the right project (existing tests and browser page checks); no defect found beyond S-1 |

## 2. Root causes and fixes

### P0-1 Path picker

*Root cause.* The buttons wrote the text input's widget key
(`st.session_state[text_key] = …`) after `st.text_input` had been drawn in
the same run, which Streamlit forbids.

*Fix* (`pages.path_picker`). All picker buttons now use `on_click`
callbacks, which run before the next rerun draws the widgets. Every target
(*Use this folder*, *Use this file*, *Open folder*, *Up*, *Go to root*) is
resolved with `GuiConfig.resolve_allowed` inside the callback; a refused
target is reported and nothing changes. The browsed folder is re-resolved on
every run and falls back to the first root when it is outside (a stale or
forged session value). Folder and file entries that resolve outside the roots
are not listed ("N link(s) pointing outside the allowed data folders are not
shown").

### S-1 Project URL guard

*Fix* (`pages.current_project`, `pages.allowed_project`). The project from
the session or `?project=` must resolve (symlinks and `..` included) inside
an available allowed root and contain project.yaml; otherwise it is dropped
from the session and the URL, and the page says it was not opened. The same
rule as *Open project by path*. The Project overview's "outside the available
data folders, can be viewed" branch is gone.

### P0-2 Metadata

*Root cause.* The editor's data (`metadata::<project>`) did not include the
user's grid edits, which exist only in the widget state of `st.data_editor`.
Every tool built its result from that data or from the pasted rows alone, then
renewed the editor key, which discards the widget state.

*Fix* (`pages._metadata_section`, `services.merge_metadata_table`,
`services.assign_metadata_value`, `services.metadata_difference`).

- After each rerun the editor's output (data plus edits) is kept as
  `metadata::<project>::working`. The tools are drawn above the grid but run
  after it, so they act on the current edits. When Streamlit has dropped the
  editor's state (another page was shown), the kept table becomes the editor's
  data, so unsaved edits survive page switches.
- Paste, CSV import and the Datasets page **merge by `sample_id`**: samples
  and columns not in the new table keep their values, an empty new cell
  changes nothing, new samples are added as rows (with a warning when they are
  not input samples), new columns are added. A new table with a missing or
  repeated `sample_id`, or a `sample_id` that is repeated in the table being
  edited, is refused with the reason.
- Before applying, the page states what would change. Replacing any existing
  non-empty value lists every replaced cell (sample, column, current, new)
  and needs a tick; bulk assignment does the same.
- An **Unsaved changes** banner compares the table with metadata.csv
  (changed values, rows not in the file, rows of the file removed, new
  columns); *Reload metadata.csv* needs a tick while there are unsaved
  changes and uses a fresh editor.
- *Save metadata* still goes through `services.save_metadata` (the project's
  metadata validator; structurally invalid tables are not written). After a
  save the editor is renewed from the saved table.

Sample identity is never inferred or reordered: rows keep their order, new
samples are appended, and nothing is derived from file names.

### R-1 Run status refresh

*Fix* (`pages._run_status_section`, `pages._run_status_body`,
`pages.live_run_in_project`). The runs section of *Validate & Run* is a
Streamlit fragment. When a run of the project is starting or running (GUI
lock of this project, or a run whose `run_state.json` is CREATED/PLANNED/
RUNNING and whose process is alive or cannot be checked), the fragment reruns
itself every 10 s: it re-reads `run_state.json`, the process record, the
`rnaseq status` view and the log tail only. The project is not re-validated
and the rest of the page is not redrawn, so typed input is kept. When the run
reaches a final state, the fragment redraws the page once and the timer stops.

### R-2 Cancel — not implemented (release finding)

Process ownership, from the code:

- The GUI starts `rnaseq run` either as a transient systemd user unit or as a
  process in its own session (`jobs.launch_run`); the launch record holds its
  PID, boot ID and start time.
- `rnaseq run` starts each Nextflow stage in a further new session
  (`service._run_command`) and, **only while a Nextflow stage runs**, turns
  SIGTERM/SIGHUP into a clean stop: it sends SIGTERM to Nextflow, waits up to
  300 s, kills it if needed, and records the run as INTERRUPTED.
- Outside a Nextflow stage, SIGTERM has its default effect. The `rnaseq`
  process ends without recording anything, so the GUI shows the run as
  INCOMPLETE. While the run is CREATED, `ensure_downstream_runtime` may be
  running `conda env create` **directly into the final runtime prefix**. That
  child is in rnaseq's process group, not in Nextflow's session, so it is not
  forwarded a signal. Stopping there leaves a half-built prefix without its
  marker. Later runs refuse it ("lacks its wheel identity marker … Remove it
  manually").

So a GUI button could stop the process tree cleanly only in some phases.
Making it safe everywhere would need changes outside the GUI: signal handling
for the whole run, and atomic runtime provisioning. Both are out of scope for
this release. Following the brief, no stop button was added. Instead, a
RUNNING run on *Validate & Run* has a **How to stop this run** panel. It
shows the recorded PID, a `ps` command to check the PID first, and
`kill -TERM <pid>`. The panel explains the INTERRUPTED / INCOMPLETE outcomes,
warns against `kill -9` and against stopping during preparation, and adds a
warning while the run is CREATED. Nothing is signalled by the GUI.

### D-1/D-2 Disabled controls

*Run*, *Save contrasts* and *Add* show the reason as a caption and as the
button's tooltip: active run, validation problems, missing or invalid case
ID (`service.validate_case_id`'s message), no contrast, the *Replace* tick,
equal levels, a contrast already in the list.

### F-1 Other reference

*Other (advanced)* shows its fields at once, in a bordered box under the
reference choice: a managed reference folder with its manifest name, or a
genome FASTA, an annotation GTF and an optional HISAT2 index. Until the
required fields are filled, the page shows a note naming them, not an error,
and *Create project* says to fill them in. Messages name GUI fields, not CLI
flags. Validation still goes through `wizard.reference_options`.

### T-1 Numeric tables

`results.numeric_columns` marks a column numeric when every non-missing cell
parses as a number (`NA`, `NaN`, empty and `NULL` are missing; `Inf`/`-Inf`
are numbers). Columns whose name ends in "id" (gene_id, ID, entrez_id,
geneID …), ratio columns (`5/100`), TRUE/FALSE and text columns stay text.
`pages.frame` gives numeric columns a float or integer dtype, so the grid
sorts them numerically, with missing values last. Floats are shown with 6
significant digits (`%.6g`), so p-values such as `7.40803e-132` are not
rounded to 0. The DEG, QC, enrichment-term and member tables all use it. The
values are parsed from the file's text for display only. Filtering, the
*Sort by* selector, *Only significant* and every download still use the TSV
text and the original files.

### I-1 FASTQ import

The existing importer has no progress hook, and adding one would mean
changing `project.import_fastq_samplesheet` (shared with the CLI). Per the
brief, the GUI now shows a blocking-operation warning while copying, with
the number of files and bytes. On completion it shows a message that
survives the rerun: files, size, elapsed time, originals unchanged, and the
kept samplesheet. A failure says whether files may already have been copied.
No percentage is shown.

### S-2/S-3 Stale messages and contrasts

The *Full validation report* is stored with a fingerprint (mtime and size of
project.yaml, metadata.csv, contrasts.csv and the input) and removed, with a
note, when the project has changed. After *Save contrasts* the page redraws
from the written files (saved table, validation), and the success message is
kept. Composed contrasts whose levels are no longer in the metadata are
removed with a warning and never saved.

## 3. Test evidence

### 3.1 Unit and AppTest (`tests/test_gui_hardening.py`, 20 tests)

| Case | Tests |
|---|---|
| 1 Folder selection | `test_use_this_folder_selects_the_browsed_folder_without_an_exception` |
| 2 File selection | `test_use_this_file_selects_the_browsed_file_without_an_exception`, `test_the_path_picker_never_browses_or_selects_outside_the_allowed_roots` |
| 3 URL guard | `test_a_url_project_outside_the_allowed_roots_is_not_opened` (direct, symlink, `..`; Project, Input, Run pages), `test_a_url_project_inside_the_allowed_roots_still_opens` |
| 4 Edits survive bulk assignment | `test_manual_metadata_edits_survive_bulk_assignment_and_are_saved` (condition, donor = pair ID, batch; saved file checked) |
| 5 Edits and imports | `test_pasting_merges_into_unsaved_edits_and_replacing_needs_a_tick`, `test_datasets_metadata_load_merges_with_unsaved_edits`, `test_reload_discards_unsaved_edits_only_after_a_tick`, four service tests of merge / assign / difference / identity errors |
| 6 Numeric sorting | `test_numeric_columns_sort_numerically_with_na_last`, `test_enrichment_columns_are_numeric_and_identifiers_and_ratios_stay_text`, `test_numeric_display_leaves_the_tsv_files_unchanged` |
| 7 Refresh without validation | `test_an_active_run_refreshes_its_status_section_only` (the fragment body runs with `validate_project` replaced by a function that fails) |
| others | `test_disabled_run_and_contrast_controls_say_why`, `test_a_full_validation_report_of_an_older_project_state_is_not_shown`, `test_other_reference_fields_are_visible_and_not_an_error_before_they_are_filled` |

The grid edits in the page tests use a test double of `st.data_editor`. It
applies edits stored per editor key, the way Streamlit applies the edits kept
in the editor's own widget state. AppTest cannot type into the grid, so the
real grid was exercised in the browser (§3.2). Against `ac6471d`, 18 of the
20 tests fail. The two that pass check behaviour that was already correct: an
allowed project opens by URL, and the TSV is unchanged.

`tests/test_gui_results.py`: one assertion updated (the DEG `padj` column is
now a float with the file's exact value, and `symbol` stays text).

`tests/test_gui*.py` and `tests/test_release_metadata.py`: 66 passed.

### 3.2 Browser (headless Chrome, `scripts/dev/gui_hardening_acceptance.py`)

The test GUI was started on port 8703 with its own `gui.yaml`
(`/home/brian/nf-rna-validation/gui-hardening-170/`). Its roots were a
fixture data folder in that directory and `/media/duobook02/brian/nf-rna`,
the latter for reading GUI-6B-L2. The user's GUI on 8501 and the GUI on 8701
were not touched. Results: `browser_checks.json` (copied to
`docs/dev/screenshots/hardening-170/`), **47/47 passed**.

| Mandatory case | Browser checks |
|---|---|
| 1, 2 | *Use this folder* and *Use this file* select the path with no exception; the `escape` symlink to an outside folder is not offered |
| 3 | `?project=` outside the roots refused directly, through a symlink and through `..`; inside opens |
| 4, 5 | donor and batch typed into the grid, then condition assigned in bulk: all kept; an unsaved edit survives a sidebar page switch; a paste replacing a typed value is disabled until ticked; the saved metadata.csv equals the expected file byte for byte; "No unsaved changes" afterwards |
| 6 | header clicks on GUI-6B-L2: DEG `baseMean` and `log2FoldChange` ascending start at the file's numeric minimum (`-5.325…`; lexical minimum `-0.000119…`); KEGG GSEA `NES` starts at `-1.864…` (lexical `-0.707…`) |
| 7 | a fixture RUNNING run (process record of a `sleep` process): status time advanced from 11:51:14 to 11:51:34 with no click, a typed but unsubmitted case ID stayed; after `run_state.json` became SUCCESS the page showed SUCCESS on its own and the time stayed the same for the next 23 s |
| 8 | no cancel button (§2, R-2); the *How to stop this run* panel with `kill -TERM <pid>` is shown for the RUNNING fixture |
| 9 | Projects, Datasets, Project overview, Input & Metadata, Contrasts, Validate & Run, Results, Reports of airway-gui-demo open without an exception; disabled reasons, stale report and Other reference checked |
| 10 | the Phase 7b download checks repeated: volcano PDF byte-identical and opens; TIFF 300 byte-identical; TIFF 600 at 2008×2008 px, 600 dpi, LZW; PCA TIFF 300, running-score TIFF 300 pixel-identical to the R PNG, dotplot PDF and TIFF 600 |

Original behaviour (port 8704, `git archive ac6471d`):
`original_behaviour.json`; screenshots `O1`, `O2`.

Formal runs: the 2,640 files of `airway-gui-demo` (including run
GUI-6B-L2) had the same SHA-256 before and after all browser sessions
(`airway-gui-demo.{before,after}.sha256`). No pipeline was run, and no
FASTQ was processed.

### 3.3 Full suite

Locked runtime `0e721ff3…-200d98f8759eaf17`, on the hardening tree:
**1,156 passed, 2 skipped** (real SAMtools/featureCounts and the real SJ
project are not present), 0 failed, 26 min
(`/home/brian/nf-rna-validation/gui-hardening-170/fullsuite.log`).

## 4. Outstanding limitations

- **No GUI stop button** (R-2). Stopping a run uses `kill -TERM` from a
  terminal, as shown in the run panel. Stopping during runtime provisioning
  leaves a prefix that must be removed by hand.
- **FASTQ import has no progress bar** (I-1); the page is blocked until the
  copy ends.
- Numeric table cells are shown with 6 significant digits; the full values
  are in the TSV downloads (unchanged).
- Grid edits typed in the last second before leaving a page are kept only if
  the grid had sent them (Streamlit sends an edit when the cell is left).
- A browser reload starts a new Streamlit session: unsaved metadata edits are
  then lost (as for every unsaved widget value). The *Unsaved changes* banner
  says that the table is not saved.
- Merging never deletes rows; to remove a sample, delete its row in the grid.

## 5. Release blockers

None from this review. R-2 is a documented limitation (manual stop with a
safe command), not a blocker. The unchanged release-level conditions stand: a
green `Release gate` CI run on the pushed `main`, and the maintainer's
decision (`release_readiness_v170.md`).

## 6. Real FASTQ acceptance checklist (manual, for the maintainer)

Not run in this review: no real FASTQ end-to-end test was executed. Use a
small real paired-end dataset with a known kit and a known expected result,
in a new project on an allowed root.

1. **Datasets.** Browse to the FASTQ folder with the picker and press *Use
   this folder* (no error). Check the file count, the total size and that
   every sample is recognised with R1/R2 pairs. Check that the samplesheet's
   absolute paths are the files you meant.
2. **New Project.** FASTQ; choose the backend. For HISAT2 + featureCounts,
   set the strandedness explicitly from the kit (dUTP/TruSeq Stranded =
   reverse); Salmon may use `auto`. Choose the registered reference; with
   *Other*, the fields show at once. Set scope L2, the design, the pairing
   column and batch if any.
3. **Input.** On Datasets or Input & Metadata, check the copy plan (files,
   bytes, free space), then copy; wait for the completion message. Verify
   **input-file identity**: `sha256sum` of a few source files equals the
   copies in `<project>/input/fastq/`, and the sample ↔ file assignment in
   `planning/imported_fastq_samplesheet.csv` is correct.
4. **Metadata.** Type or paste condition, pair and batch values; use *Assign
   values to samples*; switch pages and come back. Check that nothing was
   lost and that *Unsaved changes* lists what you expect. Save. Open
   metadata.csv and verify the **study design** sample by sample against the
   lab record.
5. **Contrasts.** Add numerator = treated, denominator = control. Check the
   sentence "Positive log2FC = higher in <treated>" (**contrast direction**)
   and save.
6. **Validate.** *Validate and refresh the plan*; no errors; read the
   warnings (KEGG snapshot pin).
7. **L1 run.** With scope L1 (project.yaml) or as the first run: enter a case
   ID, *Run*, and watch the status refresh on its own until SUCCESS. Check
   the scientific QC panel: strandedness CONSISTENT for the declared value
   (**strandedness verification**), L1 flags, MultiQC report; inspect every
   warning (**QC and warnings**).
8. **L2 run.** Reuse the upstream results of the L1 run if offered. After
   SUCCESS check the DEG counts and the R warnings expander. Sort the DEG table
   by `padj` and `log2FoldChange` with the column headers (numeric order).
   Check that a known marker gene has the expected sign (**contrast
   direction**), then look at the enrichment results.
9. **Results/Reports.** Download all_genes.tsv and compare it with
   `delivery/` (`sha256sum`). *Prepare delivery package* and check the zip
   against `delivery_manifest.yaml` (**output file integrity**). Download a
   volcano as PDF, TIFF 300 dpi and TIFF 600 dpi and open them.
10. **Stop/recovery (optional, on a throw-away run).** While the run is past
    preparation and in a Nextflow stage, follow *How to stop this run*
    (`kill -TERM <pid>`). Check that the run becomes INTERRUPTED, that `ps`
    shows no remaining `nextflow` process of that run, and that
    `rnaseq retry --retry-of CASE/RUN` starts a new attempt. Do not stop a
    run that is preparing (CREATED).

## 7. Final state

The fixes, tests, acceptance script and documents are one commit on
`feature/v1.7.0` on top of `ac6471d` (its SHA is reported with the results).
Nothing was pushed, tagged or merged.
