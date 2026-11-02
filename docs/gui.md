# Browser GUI

`rnaseq gui` serves a Streamlit GUI on `127.0.0.1` (default port 8501). It is a
front end to the same Python services as the CLI: projects it creates are
`rnaseq new` projects, its validation is `rnaseq validate`, its contrast step
is `rnaseq design`, and its runs are `rnaseq run`. It keeps no database and no
analysis state of its own. Results are read from the run's own files; the GUI
never recomputes a statistic or redraws a figure.

## Install and start

```bash
python -m pip install "nf-rna[gui] @ git+https://github.com/2002brian/nf-rna.git@v1.7.0"
# or, in a source checkout: pip install -e '.[gui]'
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

**Optional: poppler for TIFF 600 dpi.** The *TIFF 600 dpi* figure download
renders the pipeline's PDF with poppler (`pdftocairo`, `pdfinfo`), which nf-rna
does not install:

```bash
sudo apt install poppler-utils          # Debian/Ubuntu/WSL2
# or, into the GUI's Conda environment: conda install -c conda-forge poppler
```

Without poppler only that option is unavailable (its button says why and the
*Projects* page shows a notice); PDF, TIFF 300 dpi and everything else work.

## Settings: `~/.config/nf-rna/gui.yaml`

Written with defaults the first time the GUI opens. The default workspace is
`~/nf-rna-work` in the home directory of the user running the GUI (paths are
written to the file expanded, for example `/home/<user>/nf-rna-work`):

```yaml
allowed_roots:                       # the only folders the GUI may browse, read or write
- ~/nf-rna-work
project_parent: ~/nf-rna-work/projects        # default location of new projects
condarc: ~/nf-rna-work/conda/condarc.yml      # set as CONDARC for runs (unless already set); written only if the file exists, otherwise null
log_tail_lines: 200
launcher: auto                       # auto, systemd or process (see "How runs are started")
download_limit_mb: 500               # larger files: server path and scp command instead of a download
```

The values in an existing `gui.yaml` always take precedence; the file is never
rewritten once it exists. A server whose data lives elsewhere, for example on a
separate data disk, sets its own folders. This is a custom example, not a
default:

```yaml
allowed_roots:
- /media/<disk>/<user>/nf-rna       # data disk: projects/ and fastq/ live here
- /home/<user>/nf-rna-work/testdata
project_parent: /media/<disk>/<user>/nf-rna/projects
condarc: /home/<user>/nf-rna-work/conda/condarc.yml
```

`project_parent` must lie inside an allowed root: the GUI creates nothing for a
`project_parent` outside them, and the *New project* location is checked
against the roots like every other path.

Every path typed or browsed in the GUI is resolved (symlinks included) and must
lie inside an available allowed root. The folder browser never opens or lists
a link that resolves outside the roots (it says how many it hides), and a
project named in the URL (`?project=`) is opened only when it is inside an
available root; otherwise the page says it was not opened. A root under `/media`, `/mnt` or
`/run/media` is available only when it is on a mounted filesystem, so an
unmounted data disk is reported as "not mounted" instead of being mistaken for
an empty local folder. When the data root (the parent of `project_parent`) is
available, `projects/` and `fastq/` are created in it if missing; nothing is
created when it is not mounted. The default workspace `~/nf-rna-work` itself is
created on first use when it is the configured data root; any other missing
root is reported and never created. Existing projects, FASTQ folders, runs and
the reference registry are never moved or changed by these settings.

## Pages

The sidebar groups the pages in the order of the work — **Workspace**
(Projects, Datasets, Project overview), **Set up** (New Project, Input &
Metadata, Contrasts) and **Analysis** (Validate & Run, Results, Reports &
Downloads) — and selects the current project. The current project is kept in
the URL, so a reload or a bookmark returns to it; under the selector the
sidebar shows its scope, input, sample count and the status of its latest run.

Selections that belong to a project (the run chosen on Project overview,
Validate & Run, Results and Reports; the case ID; the metadata being edited;
the contrasts being composed) are kept per project, so switching to another
project and back returns to the same run and the same unsaved edits.

### Execution status and scientific QC

Every run is shown with two separate statuses:

| | Values | Source |
|---|---|---|
| **Execution** | RUNNING, SUCCESS, FAILED, INCOMPLETE | the run's `run_state.json` and process liveness: CREATED/RUNNING with a live process is RUNNING; RUNNING without its process, a recorded INTERRUPTED and an unknown state are INCOMPLETE |
| **Scientific QC** | PASS, WARNING, INCONCLUSIVE, NOT_EVALUATED | a display roll-up of verdicts the run itself recorded (below) |

The scientific QC roll-up reads, for each run:

- the strandedness check of FASTQ runs on the HISAT2 route (v1.6.1, frozen
  contract `upstream_qc.strandedness`): CONSISTENT → PASS, WARNING → WARNING,
  ADVISORY (declared unstranded, reads look stranded; unstranded counting is
  valid) → WARNING, INCONCLUSIVE → INCONCLUSIVE, Salmon route → NOT_EVALUATED;
- the L1 sample flags (`downstream/l1/backend_summary.json`, for example
  `LOW_LIBRARY_SIZE`): none → PASS, any → WARNING;
- per contrast, the R warnings of the DESeq2 fit, `results()` and shrinkage
  (v1.6.1, `contrast_metadata.json` `r_warnings`): none → PASS, any → WARNING;
  runs before v1.6.1 → NOT_EVALUATED; and a shrinkage that did not succeed →
  WARNING.

The roll-up is WARNING if any check is WARNING, otherwise INCONCLUSIVE if any
is INCONCLUSIVE, otherwise PASS if any check passed, otherwise NOT_EVALUATED.
Each check is listed with the verdict exactly as recorded. The roll-up never
changes a recorded verdict, `run_state.json` or the execution status: a run
with a QC warning is still SUCCESS.

### Projects

The research projects found in the allowed roots (two folder levels deep), one
row each: project, species, input type, quantification (Salmon, HISAT2 +
featureCounts, or none for imported counts), reference, sample count (rows of
metadata.csv), contrasts (contrasts.csv), scope, number of runs, the last run,
its execution status and scientific QC. *Open* goes to the project's overview.
Below: recent runs of all projects, the data folders and their mount state,
`rnaseq doctor` checks (cached for 10 minutes; *Re-run checks* refreshes) and
the registered references. An unavailable data folder (for example an
unmounted disk) is reported at the top. An active GUI run is announced at the
top.

### Datasets

Browse a folder inside the allowed roots (unmounted roots are not offered) and
see what it holds, three folder levels deep (run, work and output folders are
skipped): the number of files, total size, files per kind (FASTQ, count
matrix, sample metadata, FASTQ samplesheet, contrasts, other tables, other),
and files without read permission. Kinds come from the file name and, for
tables, the header line only; FASTQ contents are never read. A symbolic link
that resolves outside the allowed roots is listed as ignored and never
followed, and symlinked folders are not entered.

Recognition uses the existing rules: FASTQ samples and R1/R2 pairs by the
validator's file-name rule (paired-end or single-end; other files in the
folder are ignored); count matrices by the validator's count-matrix checks
(genes × samples); samplesheets by the `rnaseq new --fastq-samplesheet`
reader, which also reports any FASTQ path outside the allowed roots;
metadata as a preview.

*Use in a project* acts on the current project (sidebar), with the existing
importers only:

- FASTQ folder or samplesheet → a FASTQ project of the same layout: shows how
  many files and bytes would be copied into `input/fastq/` and the free space
  there, then copies with `import_fastq_samplesheet` (refused when a file is
  outside the allowed roots, the project is not on an available root, or the
  space is insufficient). The originals are not moved or changed; where they
  came from is shown on Input & Metadata and Project overview (from
  `planning/imported_fastq_samplesheet.csv`).
- Count matrix → a raw-count project: `import_count_matrix` (replacing an
  existing matrix needs a tick).
- Metadata → merged by `sample_id` into the project's metadata editor
  (unsaved edits included, replacements need a tick, as on Input & Metadata),
  unsaved; *Save metadata* on Input & Metadata checks and writes it.

### Project overview

The current project: species, input, sample count, contrasts, quantification,
reference, layout and strandedness, and the analysis settings (scope, design
formula, DEG thresholds, enrichment). The workflow row shows Dataset →
Project → Metadata → Contrasts → Validate → Run → Results, each done, to do,
needing attention or not needed, from `rnaseq validate` and the runs; the
first unfinished step is named with its findings and a link to its page.
*Data locations* gives the project input, the original FASTQ folder(s), the
metadata and contrasts files and the runs folder.

*Analysis history* lists every run (newest first) with level, execution
status, scientific QC, start and end time, and whether results and a report
exist. For the chosen run it shows each recorded QC check, and *Open
results*, *Open report & downloads* and *Run details and logs* go to that run
on the corresponding page.

### New Project

The questions of the interactive `rnaseq new` wizard, with the same choices
and defaults: name, location, species, FASTQ or raw counts; for FASTQ the
quantification backend (*Salmon (pseudoalignment)*, the default, or
*HISAT2 + featureCounts (genome alignment)*; `project.yaml` stores `salmon` or
`hisat2_featurecounts`), layout,
strandedness and reference (registered references first; iGenomes for Salmon;
local or custom roots under *Advanced settings*); the analysis scope (qc for
FASTQ only, L1, L2), the experimental design and pairing column, a known
technical batch, L2 enrichment, the gene identifier type, and CPU/memory
limits (`auto` by default).

Managed references come only from this workstation's registry
(`~/.config/nf-rna/references.yaml`, filled with `rnaseq reference register
ROOT`); the list follows the species and the quantification backend, and a
choice made for one combination is not carried over to another. When no
registered reference fits, the page says "<species> <backend> reference not
prepared", lists registered references of that species that cannot be used
with the reason and the preparation command (`rnaseq reference prepare-hisat2
ROOT`, `rnaseq reference prepare ROOT`), and fills in no path. Before the
project is created, a managed reference is loaded again and checked against
the species and backend (a Salmon index never counts as a HISAT2 index); the
validator checks it once more at run time.

*Other (advanced)* shows its fields under the reference choice as soon as it
is chosen: a managed reference folder (with its manifest file name), or a
genome FASTA and an annotation GTF (and optionally a prepared HISAT2 index).
Until the required fields are filled the page names them; it does not show an
error before you have had the chance to fill them.

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
instead. FASTQ files cannot be uploaded through the browser. The copy blocks
the page until every file is copied (the importer reports no per-file
progress, so none is shown); the page says so while copying and, at the end,
how many files and bytes were copied in how long.

**Metadata.** An editable table with `sample_id`, the design variables (for
example `cell` and `condition` for `~ cell + condition`) and any other columns
already in metadata.csv. Sample IDs are prefilled from the count-matrix header
or the FASTQ names; conditions are never derived from file names. Ways to fill
it:

- type or paste into the table (a block copied from Excel pastes into the
  selected cell and those below/right of it);
- *Paste table*: paste rows copied from Excel (tab-separated) or CSV text,
  header first, with a `sample_id` column;
- *Import CSV*: load a metadata file;
- *Assign values to samples*: set one column to one value for selected samples.

Every tool works on the table as it is on screen, unsaved edits included.
*Paste table*, *Import CSV* and Datasets' *Load into the metadata editor*
**merge by `sample_id`**. Samples and columns that are not in the new table
keep their values, and an empty cell in the new table changes nothing. New
samples are added as rows, with a warning when they are not samples of the
input, and new columns are added. Before anything is applied, the page says
what would change. Replacing a value that is not empty lists each replaced
cell (sample, column, current and new value) and needs a tick. The same holds
for *Assign values to samples*. A new table without `sample_id`, or with a
repeated or empty `sample_id`, is refused. Merging never removes a row; delete
a row in the table to remove a sample.

An **Unsaved changes** banner under the table says how it differs from
metadata.csv (changed values, rows not yet in the file, removed rows, new
columns). Unsaved edits are kept when you switch to another page or project
and come back; they are lost when the browser tab is reloaded or closed.
*Reload metadata.csv* discards them only after a tick.

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
<numerator>". *Add* and *Save contrasts* say why they are disabled (same
level twice, contrast already listed, no contrast, *Replace* not ticked).
After saving, the page is redrawn from the written files. Composed contrasts
whose levels are no longer in the metadata are removed, with a warning.

*Save contrasts* runs the `rnaseq design` transaction: it writes contrasts.csv,
freezes `design.condition_order`, sets the scope to L2, re-validates (rolling
back on failure) and refreshes the planning artifacts. Full-rank, confounding
and pairing checks come from the project validator. Saved contrasts are
replaced only after ticking *Replace the saved contrasts*.

When KEGG ORA or GSEA is selected and the project does not pin a KEGG
snapshot, the validator's warning is shown with *Pin to the current snapshot*,
which writes `annotation.enrichment.kegg.snapshot` to project.yaml.

### Validate & Run

- The validator's findings are always shown, grouped by workflow step
  (Dataset, Project, Metadata, Contrasts) with the page that fixes them; the
  page says whether the project can run. *Validate and refresh the plan*
  shows the full `rnaseq validate` report and, when the project is valid,
  refreshes the planning artifacts (as `rnaseq plan`). The report is removed,
  with a note, when project.yaml, the metadata, the contrasts or the input
  change afterwards, so an outdated "PASS" is never shown.
- Settings summary: scope (with what that scope produces), species, input,
  quantification, strandedness, reference, design, DEG thresholds, every
  contrast with its direction, enrichment and KEGG snapshot, resources; the
  exact command is in an expander.
- Downstream runtime: whether the locked runtime for this nf-rna build exists,
  or that the first run will create it (several minutes).
- Case ID (required) and, for FASTQ projects with a successful earlier run,
  the upstream results to reuse (offered only when `--reuse-upstream`'s own
  preflight accepts them; a single candidate is preselected). Raw-count
  projects have no upstream stage, so there is nothing to reuse.
  Earlier SUCCESS runs that the preflight refuses are listed with the
  preflight's own reason (for example a missing upstream handoff), and runs
  made by another nf-rna version are listed as not reusable, because the run's
  upstream identity check would refuse them after the run had started.
- *Run L1* / *Run L2* (named after the project's scope) runs the `rnaseq run`
  preflight, then starts the run in the background. It is disabled while the
  project is invalid, without a valid case ID, or while a GUI run is active;
  the reason is shown under the button. Case IDs are never generated.
- Runs: every recorded run of the project with execution status and
  scientific QC, the selected run's state and QC checks, *Open results* /
  *Open report & downloads*, the `rnaseq status` dashboard (stages, tasks,
  elapsed time, errors), and a log viewer (last N lines of the run's logs or
  the GUI launcher log; `rnaseq.log` first; *Only warnings and errors* keeps
  the lines with ERROR, WARN/WARNING, FAILED or Exception).
  While a run of the project is starting or running, this section re-reads
  the run files every 10 seconds by itself. Only this section is redrawn: the
  project is not re-validated and typed input elsewhere on the page is kept.
  When the run ends, the page is redrawn once and the automatic updates stop.
  *Refresh* re-reads the files at any time.
- **Stopping a run.** The GUI has no stop button (see
  [Known limitations](#known-limitations-v170)). A running run shows *How to
  stop this run*: its recorded `rnaseq` PID, a `ps` command to check the PID
  first, and `kill -TERM <pid>`. During a Nextflow stage, `rnaseq` stops
  Nextflow and its tasks and records the run as INTERRUPTED
  (`rnaseq retry` can start a new attempt). Outside a stage, `rnaseq` stops
  without recording anything and the run is shown as INCOMPLETE. Do not use
  `kill -9`, which leaves Nextflow running. Do not stop a run that is
  preparing (CREATED) while it creates the downstream runtime.

States are read from each run's `run_state.json`. A run recorded as RUNNING
whose `rnaseq` process no longer exists is shown as **INCOMPLETE** (`rnaseq
status` calls it INTERRUPTED (stale)); the run's files are not changed. A run
that is preparing (freezing inputs, creating the runtime) is shown as RUNNING
with that explanation.

### Results

Choose a run of the current project that has downstream outputs (any state; a
run that did not finish shows what it wrote). Sections follow the run's
analysis level: an L1 run shows only *QC & exploratory* and *Provenance*.
*QC & exploratory* starts with the run's scientific QC summary (open when it
is WARNING).

- **QC & exploratory**: the L1 summary (genes filtered and retained,
  normalization, VST), PCA, sample correlation, raw and VST boxplots, the QC
  tables, and the MultiQC report of FASTQ runs. FASTQ runs also show the
  strandedness check (v1.6.1). On the HISAT2 route it shows, per sample,
  the featureCounts `-s 1`/`-s 2` assignments, their ratio and the status.
  A warning appears when the reads contradict the declared strandedness, and a
  note when an unstranded declaration meets stranded reads. Samples without
  clear evidence are shown as inconclusive. The check never changes the
  declaration or the counts. The Salmon route shows it as not evaluated.
- **Differential expression** (one contrast at a time): the numbers of
  significant, up and down genes and the DEG thresholds and rule from
  `contrast_metadata.json`; the LFC shrinkage method and status (a failed
  shrinkage is shown as an error); the R warnings recorded during the fit,
  `results()` and shrinkage, verbatim with their step and count, in an
  expander (v1.6.1; the analysis completed, so these are not errors);
  volcano, MA plot and heatmap; and the
  `all_genes.tsv` table, searchable by symbol or gene_id, filterable by padj and
  |log2FC|, sortable, with *Only significant* (the genes of `significant.tsv`).
  DEGs are called on the raw `log2FoldChange`; `log2FoldChange_shrunk` is for
  display and ranking only.
- **Enrichment** (GO ORA, KEGG ORA, GO GSEA, KEGG GSEA): significance
  semantics and rule, BH adjustment, ORA universe definition and size, GSEA
  ranking metric, seed and nPermSimple, the KEGG snapshot and whether the
  project pinned it. Per result set (ORA foreground × ontology, GSEA ontology):
  status and reason, dotplots, simplified GO dotplots, running-score plots, the
  significant-terms table (or all evaluated terms), and the member genes of the
  pathway selected in the table (from `*_members.tsv`, significant pathways
  only). Output written before v1.5 has no significance-semantics record and
  is labelled **legacy**.
- **Provenance**: nf-rna version and commit, Nextflow, the downstream lock
  SHA-256 and runtime prefix, input and reference, the R/Bioconductor package
  versions and the KEGG snapshot recorded by the backends. For FASTQ runs
  (v1.6.1) it also shows the strandedness-check status and the CPUs (threads)
  that each relevant upstream process actually used, read from the run's own
  Nextflow trace, and Salmon's recorded `--threads`.

Tables are read from the TSV text. For display, columns whose every value is
a number become numeric columns, so clicking a column header sorts numerically
(missing values, `NA`, last). Floating-point values are shown with 6
significant digits, so very small p-values are not rounded to 0.
Identifiers (gene_id, ID, Entrez IDs), ratios such as `GeneRatio` and text
stay as written. Filters and *Sort by* use the file's values, and every
download is the original file, unchanged.
Figures are the pipeline's PNG files, each with PNG, PDF, *TIFF 300 dpi* and
*TIFF 600 dpi* downloads (see [Figure downloads](#figure-downloads-pdf--tiff)). An
output that does not exist is shown as *unavailable* with the backend's
recorded reason (for example "no significant pathway"); nothing is filled in.

### Reports & Downloads

- The run's `report.html` (the delivered copy), shown in the page or
  downloaded. It is self-contained: figures are embedded and nothing is loaded
  from the network.
- Single files by group: counts (normalized, VST), QC tables, DEG tables
  (all / significant / up / down per contrast), enrichment tables and member
  tables. L1 runs offer no L2 file.
- *Figures (PDF / TIFF)*: every figure of the run (one per pipeline PDF), as
  PDF, TIFF 300 dpi (the default), TIFF 600 dpi or the PNG report preview. The
  caption says where the chosen file comes from.
- The delivery package: *Prepare delivery package* checks every file of
  `delivery/` against `delivery_manifest.yaml` and writes a zip to
  `~/.cache/nf-rna/gui/packages/` (the run directory is not written), then
  offers it for download.

Files and packages larger than `download_limit_mb` are not offered for
download; the page shows their server path and an `scp` command instead.

### Figure downloads (PDF / TIFF)

Nothing is redrawn from data and nothing is written into the run; generated
TIFF files are made in memory when the button is clicked.

| Download | Source |
|---|---|
| PDF | the pipeline's Cairo PDF, byte for byte |
| TIFF 300 dpi | the pipeline's TIFF, byte for byte, when it was drawn at 300 dpi (every figure except GSEA running-score plots, under the default profiles) |
| TIFF 300 dpi, running-score plots | the pipeline's `<name>_300dpi.png`, drawn by R at 300 dpi, converted losslessly (every pixel unchanged) |
| TIFF 600 dpi | the pipeline's PDF rasterized at 600 dpi with poppler (`pdftocairo`) |

A TIFF at a given dpi always has `round(page_pt / 72 × dpi)` pixels, the same
rule as the pipeline's own devices, RGB on white, LZW, with the dpi in its
resolution tags. A raster is never enlarged to reach a higher dpi: without a
raster drawn at that dpi and without a PDF (or without poppler for the PDF),
the option is shown as unavailable with the reason. Under the `elsevier`
profile (500 dpi pipeline TIFF) both TIFF options are rendered from the PDF.

The 600 dpi TIFF follows the PDF's text layout. R's PDF device measures text
slightly differently from its raster devices, so in the pipeline's own output
the panel edge of a PDF can sit up to about 3 mm (10 pt, on enrichment
dotplots with long term names) from that of its PNG/TIFF; data, colours, labels
and order are the same. The PDF and TIFF 300 dpi downloads are the pipeline's
files as they are.

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

## Known limitations (v1.7.0)

- Assigning FASTQ files to a project copies them into the project (the
  existing importer; the project schema requires `input.path` inside the
  project). The Datasets page shows the size and free space before copying,
  but there is no "use in place" mode.
- The original location of an imported count matrix is not recorded by the
  importer, so Project overview cannot show it.
- `/dashboard` is not a URL of the Projects page; use the server root.
- Scope changes (for example running an L2 project as L1) remain project.yaml
  edits; the Run page runs the project's scope.
- TIFF 600 dpi needs poppler (`pdftocairo`, `pdfinfo`) on the server; without
  it the option is unavailable. There is no TIFF beyond 600 dpi, no SVG and no
  figure editing.
- **No stop button.** `rnaseq run` stops cleanly on SIGTERM only while a
  Nextflow stage runs. Before that, the downstream runtime may be being
  created directly in its final folder, and stopping would leave it half-built
  (later runs then refuse it until it is removed by hand). The GUI therefore
  only shows how to stop a run from a terminal (*How to stop this run*).
- FASTQ copying shows no per-file progress (the importer has none); the page
  is blocked until it ends and then reports what was copied.
- Numbers in the result tables are displayed with 6 significant digits; the
  downloaded TSV files have the full values.
- Unsaved metadata edits live in the browser session: they survive page and
  project switches, not a reload of the tab.

## Known limitations (v1.6.0)

- One user, one machine: the server listens on 127.0.0.1, there are no accounts, and one GUI run can be active at a time.
- Verified end to end on raw-count projects. FASTQ recognition and import were verified, but a FASTQ run has not been started or browsed from the GUI; the MultiQC download is tested only with fixtures.
- The `systemd` launcher was verified from a GUI started in a terminal with `launcher: systemd`; running the GUI itself as the example service was not tested.
- With the `process` launcher, SIGKILL of `rnaseq` leaves its Nextflow child running (SIGTERM stops both).
- Legacy (pre-v1.5) enrichment display was verified with test fixtures only.
- FASTQ copying blocks the Input page until it finishes; other covariates, `rnaseq retry` and `rnaseq clean` are CLI-only; delivery zips accumulate in `~/.cache/nf-rna/gui/packages/`.

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
