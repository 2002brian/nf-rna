# Browser GUI

`rnaseq gui` serves a Streamlit GUI on `127.0.0.1` (default port 8501). It is a
front end to the same Python services as the CLI: projects it creates are
`rnaseq new` projects, its validation is `rnaseq validate`, its contrast step
is `rnaseq design`, and its runs are `rnaseq run`. It keeps no database and no
analysis state of its own.

Result browsing and downloads are not part of v1.6.0 Phase 6a; a finished run
shows the path of its `delivery/` folder.

## Install and start

```bash
pip install -e '.[gui]'       # source checkout; or pip install 'nf-rna[gui]'
tmux new -s nf-rna-gui        # or attach: tmux attach -t nf-rna-gui
rnaseq gui                    # --port 8600 to use another port
```

The server listens only on the loopback interface; there is no option to
listen on other addresses. From another machine, forward the port over SSH:

```bash
ssh -N -L 8501:127.0.0.1:8501 user@server    # then open http://127.0.0.1:8501
```

Without the `gui` extra, every CLI command works as before and `rnaseq gui`
prints the install command.

## Settings: `~/.config/nf-rna/gui.yaml`

Written with defaults the first time the GUI opens:

```yaml
allowed_roots:                       # the only folders the GUI may browse, read or write
- /media/duobook02/brian/nf-rna
- /home/brian/nf-rna-work/testdata
project_parent: /media/duobook02/brian/nf-rna/projects   # default location of new projects
condarc: /home/brian/nf-rna-work/conda/condarc.yml       # set as CONDARC for runs (unless already set)
log_tail_lines: 200
launcher: auto                       # auto, systemd or process (see "How runs are started")
```

Every path typed or browsed in the GUI is resolved (symlinks included) and must
lie inside an available allowed root. A root under `/media`, `/mnt` or
`/run/media` is available only when it is on a mounted filesystem, so an
unmounted data disk is reported as "not mounted" instead of being mistaken for
an empty local folder. When the data root (the parent of `project_parent`) is
mounted, `projects/` and `fastq/` are created in it if missing; nothing is
created when it is not mounted.

## Pages

The sidebar selects the page and the current project (also kept in the URL, so
a reload or a bookmark returns to it).

### Dashboard

nf-rna version, allowed roots and their mount state, the projects found in the
allowed roots (two folder levels deep) with their latest run, recent runs of
all projects, `rnaseq doctor` checks (cached for 10 minutes; *Re-run checks*
refreshes), and the registered references. An active GUI run is announced at
the top.

### New Project

The questions of the interactive `rnaseq new` wizard, with the same choices
and defaults: name, location, species, FASTQ or raw counts; for FASTQ the
quantification backend (Salmon or HISAT2 + featureCounts), layout,
strandedness and reference (registered references first; iGenomes for Salmon;
local or custom roots under *Advanced settings*); the analysis scope (qc for
FASTQ only, L1, L2), the experimental design and pairing column, a known
technical batch, L2 enrichment, the gene identifier type, and CPU/memory
limits (`auto` by default).

HISAT2 + featureCounts has no default strandedness: choose `reverse`,
`forward` or `unstranded`. The page explains the common kits (dUTP/TruSeq
Stranded = reverse) and how to find out when unsure. Salmon offers `auto`.

Like the wizard, the GUI creates a scaffold: project.yaml, a header-only
metadata.csv for the chosen design, an empty contrasts.csv and the input
folder. Data is added on the next page.

### Input & Metadata

**Input.** Raw counts: choose a count-matrix file inside the allowed roots; it
is copied to the project's `input/counts.csv` (as `rnaseq new --counts`).
FASTQ: choose a folder; sample IDs and R1/R2 pairs are recognised by the
validator's file-name rule (`<sample>[_Lnnn]_R1|R2[_nnn].fastq.gz`), shown as
a table and as a samplesheet (`sample,fastq_1,fastq_2,strandedness`) you can
download. *Copy into the project* copies the files into `input/fastq/` (as
`rnaseq new --fastq-samplesheet`); an existing samplesheet can be imported
instead. FASTQ files cannot be uploaded through the browser.

**Metadata.** An editable table with `sample_id`, the design variables (for
example `cell` and `condition` for `~ cell + condition`) and any other columns
already in metadata.csv. Sample IDs are prefilled from the count-matrix header
or the FASTQ names; conditions are never derived from file names. Ways to fill
it:

- type or paste into the table (a block copied from Excel pastes into the
  selected cell and those below/right of it);
- *Paste table*: paste rows copied from Excel (tab-separated) or CSV text,
  header first;
- *Import CSV*: load a metadata file;
- *Assign values to samples*: set one column to one value for selected samples.

*Save metadata* checks the table with the project's metadata validator. A
structural problem (no or duplicate `sample_id`, duplicate column) is shown
beside the column and nothing is written. Empty design values are allowed:
the table is saved, and the missing values are listed so the work can be
finished later. Disagreement with the count matrix or FASTQ samples is listed
too.

### Contrasts

Shows the design type, DESeq2 formula, pairing column and the levels of the
comparison factor with their sample counts. Contrasts are composed from those
levels (numerator = test group, denominator = reference group); *Propose all
pairs* applies the `rnaseq design` rule. Each contrast shows its ID, formula,
numerator and denominator, and the sentence "Positive log2FC = higher in
<numerator>".

*Save contrasts* runs the `rnaseq design` transaction: it writes contrasts.csv,
freezes `design.condition_order`, sets the scope to L2, re-validates (rolling
back on failure) and refreshes the planning artifacts. Full-rank, confounding
and pairing checks come from the project validator. Saved contrasts are
replaced only after ticking *Replace the saved contrasts*.

When KEGG ORA or GSEA is selected and the project does not pin a KEGG
snapshot, the validator's warning is shown with *Pin to the current snapshot*,
which writes `annotation.enrichment.kegg.snapshot` to project.yaml.

### Run

- *Validate* shows the `rnaseq validate` report and, when the project is
  valid, refreshes the planning artifacts (as `rnaseq plan`).
- Settings summary: scope, input, design, every contrast with its direction,
  enrichment and KEGG snapshot, resources, and the exact command.
- Downstream runtime: whether the locked runtime for this nf-rna build exists,
  or that the first run will create it (several minutes).
- Case ID (required) and, for FASTQ projects with a successful earlier run,
  the upstream results to reuse (offered only when `--reuse-upstream`'s own
  preflight accepts them; a single candidate is preselected). Raw-count
  projects have no upstream stage, so there is nothing to reuse.
- *Run L1* / *Run L2* (named after the project's scope) runs the `rnaseq run`
  preflight, then starts the run in the background.
- Runs: every recorded run of the project with its state, the `rnaseq status`
  dashboard of the selected run (stages, tasks, elapsed time, errors), and a
  log viewer (last N lines of the run's logs or the GUI launcher log).
  *Refresh* re-reads the files; there is no automatic polling.

States are read from each run's `run_state.json`. A run recorded as RUNNING
whose `rnaseq` process no longer exists is shown as **INCOMPLETE** (`rnaseq
status` calls it INTERRUPTED (stale)); the run's files are not changed. A run
that is preparing (freezing inputs, creating the runtime) is shown as RUNNING
with that explanation.

## How runs are started

*Run* starts `rnaseq run PROJECT --case-id CASE --yes [--reuse-upstream ...]`
as an argument list, never through a shell:

- **systemd** — a transient user unit `nf-rna-run-<id>.service`
  (`systemd-run --user`). It belongs to the user's service manager, so
  stopping or restarting the GUI (or its service) does not stop the run.
- **process** — a child in a new session (`start_new_session`), which outlives
  the GUI process.

`launcher: auto` uses systemd when the GUI itself runs as a systemd service,
or when lingering is enabled for the user; otherwise process. (Without
lingering, the user's service manager stops at the last logout, and its units
with it.)

Only one GUI run can be active. Submitting creates
`~/.cache/nf-rna/gui/active_run.lock` atomically; a second click or a second
browser tab is refused while the lock's process is alive. The lock is removed
when that process has ended. Each submission also leaves
`~/.cache/nf-rna/gui/launches/<id>.json` (command, launcher, process identity)
and `~/.cache/nf-rna/gui/logs/<id>.log` (the command's own output, which
exists before the run's log folder does). These files are what the GUI uses to
find its runs again after a restart; the analysis state stays in the run.

## Running the GUI as a user service

`docs/systemd/nf-rna-gui.service` is an example; adjust the paths, then:

```bash
mkdir -p ~/.config/systemd/user
cp docs/systemd/nf-rna-gui.service ~/.config/systemd/user/
systemctl --user daemon-reload
systemctl --user enable --now nf-rna-gui.service
sudo loginctl enable-linger "$USER"     # keep user services (and runs) running after logout
systemctl --user status nf-rna-gui.service
journalctl --user -u nf-rna-gui.service -f
systemctl --user list-units 'nf-rna-run-*'   # runs started from the GUI
```

Nothing in nf-rna installs or enables this unit.
