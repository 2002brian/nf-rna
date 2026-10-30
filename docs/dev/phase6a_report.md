# nf-rna v1.6.0 Phase 6a — GUI core (project → run)

Branch `feature/v1.6.0`, created from `main` at `9a7352b` (v1.5.1, not tagged).

## Current state before Phase 6a (code as of `9a7352b`)

Recorded from the code, not from documentation.

### Functions the GUI can call as they are

| Need | Existing function |
|---|---|
| Load and validate a project | `project.load_project`, `validators.validate_project` (returns a `ValidationReport` with counts/FASTQ/metadata/contrast summaries, issues, `state`, `execution_ready`) |
| Validation report text | `cli.render_validation_report(report)` |
| Create a project | `project.create_project(..., scaffold=True)`; the interactive wizard always creates a scaffold |
| Metadata template / design | `project.derive_design`, `project.metadata_template_header`, `project.resolve_enrichment_defaults` |
| Metadata and contrast checks | `validators.validate_metadata`, `validate_contrasts`, sample agreement, model-matrix rank, confounding, pairing (`ValidationReport.pairing`) |
| FASTQ recognition | `validators.validate_fastq(dir, layout, report)` with `FASTQ_NAME_PATTERN` (`<sample>[_Lnnn]_R1/_R2[_nnn].fastq.gz`) |
| Contrast proposal | `design_setup.build_proposal`, `generate_contrasts` (ID `<numerator>_vs_<denominator>`), `write_design`, `rollback_design` |
| Planning artifacts | `planner.generate_plan(report)`; `rnaseq run` refuses a stale plan (`execution.require_fresh_plan`) |
| Run preflight | `service.validate_case_id`, `service.preflight_reuse_upstream`, `service.prepare_service_run` (resources, fresh plan, reference readiness, Nextflow, Conda, downstream-runtime preflight, KEGG snapshot/REST) |
| Run execution | `rnaseq run PROJECT --case-id CASE [--reuse-upstream CASE/RUN] --yes` (`service.execute_service_run`) |
| Run state | `run_status.discover_runs`, `effective_state`, `build_status`, `render_status`, `process_liveness` |
| Doctor / runtime | `execution.doctor_checks`, `doctor_readiness`; `downstream_runtime.ensure_downstream_runtime` |
| References | `references.registered_local_references`, `compatible_registered_references` |
| KEGG snapshot | `enrichment_config.latest_kegg_snapshot`, `KEGG_CODES`; the validator warns `kegg_snapshot_not_pinned` |

### Logic that existed only in the interactive CLI

- The wizard's choice lists, descriptions and defaults (species, input type,
  layout, preprocessing, backend, strandedness per backend, reference source
  and default per backend, scope per input type, design types, gene ID types,
  default pairing column `patient`), the conversion of a registered reference
  into `reference.root`/`manifest`, the `auto` resource preview and
  `_reference_options` — all inside `cli.new_project`.
- The `rnaseq design` transaction inside `cli.design_command`: blocking errors
  except `missing_contrast_rows`, the residual-degrees-of-freedom check, the
  enrichment selection and annotation decision (`_design_enrichment`), write →
  re-validate → roll back → refresh the plan.
- Importing FASTQ files or a count matrix into an existing project: only
  `create_project` copied them, at creation time.
- No function pinned a KEGG snapshot in project.yaml (users edited the file).

### Existing run-state mechanism

- `create_case_run` allocates `runs/<case>/<run-id>/` (Asia/Taipei timestamp,
  `-NN` suffix on collision) and writes `run_state.json` as `CREATED`.
- `execute_service_run` records the owning process in
  `logs/rnaseq.process.json` (pid, hostname, boot ID, kernel start ticks),
  appends to `logs/rnaseq.log`, provisions the downstream runtime while the
  state is still `CREATED`, then moves through `RUNNING` with phases
  `freeze` → `upstream` → `downstream` → `delivery`, ending `SUCCESS`, `FAILED`
  (with `error`) or `INTERRUPTED` (SIGINT/SIGTERM/SIGHUP). Nextflow stdout and
  stderr go to `logs/<stage>.stdout.log` / `.stderr.log`.
- `run_status.effective_state` reports a recorded `CREATED`/`RUNNING` run whose
  process is gone as `INTERRUPTED` with a "stale" note, without writing.
- Runs refuse a dirty or unidentified source checkout
  (`require_identified_source`), and each nf-rna wheel gets its own runtime
  prefix, so the first run of a new build provisions a new Conda prefix.

## Changed and added files

| File | Change |
|---|---|
| `src/rnaseq/gui/__init__.py` | `rnaseq gui` launcher: Streamlit command bound to 127.0.0.1, extra detection, install hint. No Streamlit import. |
| `src/rnaseq/gui/config.py` | `~/.config/nf-rna/gui.yaml`, allowed-roots guard (`resolve_allowed`), mount check, `projects/` and `fastq/` creation on a mounted data root. |
| `src/rnaseq/gui/jobs.py` | Background `rnaseq run`: atomic active-run lock, systemd transient unit or detached process, launch records, launcher log, launch → run matching, INCOMPLETE label, log tail. |
| `src/rnaseq/gui/services.py` | GUI helpers composed from existing services: project discovery, runs, metadata save/paste, FASTQ discovery, design view, contrast proposal, KEGG pin state, reuse candidates, run preflight. |
| `src/rnaseq/gui/pages.py`, `app.py` | The five Streamlit pages and the entry point. |
| `src/rnaseq/wizard.py` | New: wizard choices/defaults, `reference_options`, `managed_reference_settings`, `resource_suggestion` (moved out of `cli.py`). |
| `src/rnaseq/design_setup.py` | New shared functions: `design_blockers`, `residual_df_blocker`, `design_enrichment_selection`, `design_needs_gene_id_type`, `resolve_design_enrichment`, `apply_design`, `DesignRejected`, `pin_kegg_snapshot`. |
| `src/rnaseq/project.py` | New: `import_fastq_samplesheet`, `import_count_matrix`, `fastq_samplesheet_text`; `create_project` uses the shared FASTQ copy. |
| `src/rnaseq/downstream_runtime.py` | New: `planned_runtime_prefix` (which prefix the next run uses; no environment is created). |
| `src/rnaseq/cli.py` | Wizard and `rnaseq design` call the shared functions (same prompts and messages); new `rnaseq gui [--port]`. |
| `pyproject.toml` | Extras `gui` (Streamlit) and `gui-dev` (+ Playwright). |
| `tests/test_gui.py` | 14 tests: allowed roots, mount check, metadata save, paste parsing, GUI contrasts = `rnaseq design`, FASTQ recognition/import, single active run, stale lock, INCOMPLETE, argument-list command / loopback, `rnaseq gui` without the extra, and AppTest pages (dashboard, new project, run button disabled while active). |
| `tests/test_design_lifecycle.py` | Two rollback tests patch the shared module the rollback moved to. |
| `scripts/dev/gui_browser_acceptance.py` | Playwright/Chrome driver for the browser acceptance below. |
| `docs/gui.md`, `docs/systemd/nf-rna-gui.service`, `README.md`, `README_zh-TW.md` | GUI documentation, user-service example (not installed), README sections. |
| `docs/dev/DECISIONS.md`, this report, `docs/dev/screenshots/6a/` | Decisions, evidence, screenshots. |

## Background execution design

1. *Run* first performs `rnaseq run`'s own preflight in the GUI
   (`validate_project`, `validate_case_id`, `preflight_reuse_upstream`,
   `prepare_service_run`), so a refusal is shown before anything starts.
2. The lock `~/.cache/nf-rna/gui/active_run.lock` is created with
   `O_CREAT|O_EXCL`. If it exists and its process is alive, the submission is
   refused; if the process is gone it is removed and creation is retried once.
3. `rnaseq run PROJECT --case-id CASE --yes [--reuse-upstream CASE/RUN]` is
   started as an argument list: as a transient unit `nf-rna-run-<id>` via
   `systemd-run --user` (under a systemd GUI service, or with lingering), or
   otherwise as a detached process in a new session. Its stdout/stderr go to
   `~/.cache/nf-rna/gui/logs/<id>.log`; CONDARC from gui.yaml is set when the
   environment has none.
4. The process identity (pid, hostname, boot ID, kernel start ticks) is stored
   in the lock and in `launches/<id>.json`.
5. Everything shown afterwards is read from files: the run directory is found
   by matching that identity with `logs/rnaseq.process.json`, its state comes
   from `run_state.json` through `run_status.effective_state`/`build_status`,
   and liveness from `run_status.process_liveness`. After a GUI restart the
   same files give the same view. The lock is removed by the first page
   render that finds its process gone.

## Acceptance evidence

Environment: `admin-MS-7E34`, Streamlit 1.66.1, Google Chrome headless via
Playwright 1.64.1 (`scripts/dev/gui_browser_acceptance.py`), GUI started
detached as from tmux (`rnaseq gui`, 127.0.0.1:8501). `/media/duobook02` is the
mounted ntfs3 disk `/dev/sde1`; the GUI created `projects/` and `fastq/` there
on first start (both were missing). Lingering is off, so `launcher: auto`
selected the **process** launcher; the kill test used a second GUI on port
8502 with `launcher: systemd`. Screenshots: `docs/dev/screenshots/6a/`.

| Item | Result | Evidence |
|---|---|---|
| `pip install -e '.[gui]'`, `rnaseq gui` starts | PASS | Streamlit server on `127.0.0.1:8501` only (`ss -ltnp`: `LISTEN 127.0.0.1:8501`); every page renders with no exception. |
| Without the extra the CLI is unchanged | PASS | Separate venv with `.[dev]` only (no Streamlit): `rnaseq gui` → exit 1 with the install command; full suite 1011 passed, 47 skipped, 5 failed + 10 errors — the same figure tests fail identically on `main` (`9a7352b`) in that venv (system R 4.3.3 lacks the locked R packages/fonts), plus one real failure, fixed: `gui/config.py` contained the literal `/home/brian`, which the local-path guard forbids (defaults now use `~`). |
| Full suite with the locked R runtime | PASS | `1071 passed, 2 skipped in 1677.93s` (the two machine-specific skips of earlier releases). `tests/test_gui.py`: 14 passed, also re-run after the config fix together with the path guard (15 passed). |
| New project `airway-gui-demo` created in the browser | PASS | Raw counts, L2, paired_two_group, pairing column `cell`, enrichment on, ENSEMBL (B1–B3). |
| Counts imported | PASS | `/home/brian/nf-rna-work/testdata/airway/counts.csv` → `input/counts.csv`, "63677 genes × 8 samples" (C1). |
| Metadata: Excel paste into the grid | PASS | The `cell` column was pasted from the clipboard with Ctrl+V into the data editor (C2). |
| Metadata: incomplete save, then reload | PASS | Saved with the conditions still empty; the missing values were listed under `condition` (C3); a fresh page showed the saved `cell` values. |
| Metadata: assign values, TSV paste | PASS | `untrt`/`trt` assigned to the selected samples (C4); then the whole table pasted as TSV (C6). Both times metadata.csv was byte-identical to `testdata/airway/metadata.csv`. |
| Contrast trt vs untrt (pair_id = cell) | PASS | Preview "Positive log2FC = higher in trt (than in untrt)" (D2); saved by the shared `rnaseq design` transaction (D3): `trt_vs_untrt,condition,trt,untrt`, `condition_order: [untrt, trt]`. |
| KEGG snapshot pinned | PASS | Validator warning shown, *Pin to the current snapshot* wrote `annotation.enrichment.kegg.snapshot: hsa_kegg20261028_20261030T061908Z` (D4) — the snapshot `airway-e2e` used. |
| GUI Validate and CLI `rnaseq validate` | PASS | GUI: `PASS WITH WARNINGS`, plan refreshed (E1). CLI on the project: exit 0, `PASS WITH WARNINGS` (only warning: 30208 all-zero genes, same as airway-e2e), `READY`. |
| Same structure as the CLI | PASS | The same project made with `rnaseq new --scaffold … --pair-id cell --gene-id-type ENSEMBL` + `rnaseq design --yes` (`/home/brian/nf-rna-validation/gui-6a/cli-equivalent/`): identical file tree, metadata.csv and contrasts.csv; project.yaml differs only by the GUI's KEGG pin (3 lines). |
| Run L2 + GO/KEGG/GSEA from the GUI to SUCCESS | PASS | Case `GUI-L2-001`, run `20261030-181614+0800`: first run showed "creates the downstream runtime" (E2a) and RUNNING/preparing while prefix `…-6a0d6922b98826fb` was built (E4); SUCCESS in 9m 14s (E6). Provenance `git_commit: 7cdb329…`; KEGG `snapshot_id: hsa_kegg20261028_20261030T061908Z`. The first attempt was correctly refused by the preflight because screenshots written into the repo made the checkout dirty (E2x); screenshots were then written outside the repo. |
| Result tables equal airway-e2e | PASS | SHA-256 of every `.tsv`/`.csv` under `downstream/` vs `airway-e2e/runs/AIRWAY-E2E/20261030-153417+0800`: **135 of 135 identical**, none missing or extra (`/home/brian/nf-rna-validation/gui-6a/compare_tables.py`). |
| Close and restart Streamlit during the run | PASS | Streamlit killed 30 s into the run (0 Streamlit processes; `rnaseq run` pid 1708887 still running); restarted 20 s later; the Run page found the run again as RUNNING, elapsed 0m 59s (E5), and later SUCCESS. |
| Two quick clicks on Run | PASS | One launch record, one run directory (`runs/GUI-L2-001/20261030-181614+0800`), one `rnaseq run` process. A second browser tab showed *Run L2* disabled with "A GUI run is active" (E3). |
| Kill a run's process | PASS | `GUI-KILL-001` (systemd unit `nf-rna-run-20261030-182628-6196d9`) killed with `kill -9` in the downstream phase. systemd stopped the unit's remaining processes (no Nextflow/Java/Rscript left). GUI: INCOMPLETE, "process no longer exists" (E7); `run_state.json` unchanged (SHA-256 checked; still RUNNING; `rnaseq status`: INTERRUPTED (stale)); lock removed. |
| New run after the kill | PASS | `GUI-L2-002` started from the same GUI (systemd launcher), SUCCESS in 6m 18s (E9); 135/135 tables identical to airway-e2e. |
| Existing CLI project `testdata/airway-e2e` | PASS | Shown on the dashboard and Run page: SUCCESS 8m 26s, all stages done, delivery path (E8). Nothing in it was written. |
| FASTQ recognition and samplesheet | PASS | Ten empty files in `/media/duobook02/brian/nf-rna/fastq/gui-6a-recognition-test/` (4 samples, `Trt2_S4` in two lanes): table of 5 R1/R2 assignments and samplesheet (C7); imported into a scaffold FASTQ project, which `rnaseq validate` read as 4 paired-end samples (C8). Nothing was executed. Afterwards the test folder and the test project were deleted. |
| Page screenshots | Done | A1 dashboard, B0 new project, C0 input & metadata, D0 contrasts, E0 run, plus the step screenshots above. |

## Not done / limitations

- The GUI service was not installed or enabled, and lingering was not
  enabled (as instructed). The systemd launcher was exercised only through a
  second GUI forced to `launcher: systemd`; a GUI running as a systemd service
  was not tested.
- With `launcher: process`, killing the `rnaseq` process with SIGKILL leaves
  its Nextflow child running (it is not in a unit). SIGTERM is handled by
  `rnaseq run` itself, which stops Nextflow and records INTERRUPTED.
- FASTQ import copies files inside the Streamlit request; very large copies
  block that page until done.
- Covariates other than `batch` and the pairing column, custom formulas and
  imported contrasts stay CLI-only (`rnaseq new --metadata/--covariate …`).
- `rnaseq retry` and `rnaseq clean` have no GUI controls.
- Result browsing and downloads are Phase 6b.

## Start commands

Server (inside tmux):

```bash
tmux new -s nf-rna-gui
cd /home/brian/projects/nf-rna
/home/brian/miniforge3/envs/nf-rna/bin/rnaseq gui
```

Mac:

```bash
ssh -N -L 8501:127.0.0.1:8501 brian@<server>
open http://127.0.0.1:8501
```
