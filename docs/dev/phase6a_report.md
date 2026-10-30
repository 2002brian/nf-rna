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
