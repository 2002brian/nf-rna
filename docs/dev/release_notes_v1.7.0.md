# nf-rna v1.7.0 — Local Research Workspace & Publication Figure Export

> **No analysis value changes.** DESeq2, tximport, LFC shrinkage, GO/KEGG ORA,
> GSEA, Salmon, HISAT2 + featureCounts commands and parameters, scientific
> thresholds, references and annotation and the Conda locks are identical to
> v1.6.1, and every numerical table is unchanged for projects that pass
> validation; v1.6.1 results do not need to be recomputed. What changes: the
> optional browser GUI, how L1 and L2 figures are drawn, which HISAT2-route
> files are published and copied on reuse, and stricter input validation —
> **FASTQ sample IDs and contrast IDs must match `^[A-Za-z0-9][A-Za-z0-9._-]*$`**
> (breaking for names with spaces, non-ASCII letters or shell characters).

Full list of changes: [CHANGELOG](https://github.com/2002brian/nf-rna/blob/v1.7.0/CHANGELOG.md#170--2026-11-06).
Release readiness: [`docs/dev/release_readiness_v170.md`](https://github.com/2002brian/nf-rna/blob/v1.7.0/docs/dev/release_readiness_v170.md).

## Local Research Workspace

- **Datasets** — browse the server's data folders (only the allowed roots of
  `~/.config/nf-rna/gui.yaml`; unmounted disks are shown as unavailable). File
  kinds, sizes and unreadable files; FASTQ and count-matrix samples recognised
  with the validator's rules; symlinks leaving the allowed roots are ignored. A
  dataset is assigned to a project with the existing importers (FASTQ is copied
  after a path and free-space check; source files are never moved or changed).
- **Project overview** — species, input, quantification, reference, samples,
  contrasts, L1/L2 settings, the workflow steps Dataset → Project → Metadata →
  Contrasts → Validate → Run → Results with the next step, and every run of the
  project with buttons to its results and report. Projects are switched from
  the sidebar and keep their own run choice.
- **Run history and QC display** — execution status (RUNNING, SUCCESS, FAILED,
  INCOMPLETE) is shown separately from scientific QC (PASS, WARNING,
  INCONCLUSIVE, NOT_EVALUATED). QC is a roll-up of verdicts the run recorded
  (v1.6.1 strandedness check, L1 sample flags, DESeq2/apeglm/ashr warnings,
  shrinkage status); it never changes the execution status or a run file.
- **Validate & Run** — settings summary, validation findings grouped by step,
  reuse hints with the reason a run cannot be reused, log filter for warnings
  and errors; *Run* is disabled while the project is invalid.

## PDF/TIFF figure downloads

Every figure (PCA, sample correlation, boxplots, volcano, MA plot, DEG
heatmap, GO/KEGG ORA and GSEA dotplots, GSEA running-score plots) can be
downloaded from *Results* and *Reports & Downloads* as **PDF**, **TIFF 300 dpi**
(default) or **TIFF 600 dpi**:

| Download | Source |
|---|---|
| PDF | the pipeline's Cairo PDF, byte for byte |
| TIFF 300 dpi | the pipeline's TIFF, byte for byte; for GSEA running-score plots a lossless TIFF of the R-drawn 300 dpi PNG |
| TIFF 600 dpi | the pipeline's vector PDF rendered at 600 dpi with poppler |

Nothing is redrawn from data, nothing is written into runs, and a raster is
never enlarged. Every TIFF is RGB, LZW, with its dpi in the resolution tags and
`round(page_pt / 72 × dpi)` pixels.

## Safer data entry and run monitoring

Found in a usability review and fixed before release
([`docs/dev/gui_hardening_v170.md`](https://github.com/2002brian/nf-rna/blob/v1.7.0/docs/dev/gui_hardening_v170.md)):

- **Metadata is never discarded silently.** Assigning values, pasting a table,
  importing a CSV or loading metadata from *Datasets* now works on the table
  as it is on screen, unsaved edits included, and merges by `sample_id`.
  Replacing an existing value lists each cell and needs a tick. An *Unsaved
  changes* banner shows what is not yet in metadata.csv, and unsaved edits are
  kept when you switch pages or projects.
- **Folder and file browser** — *Use this folder* / *Use this file* no longer
  fail with a Streamlit error, and links that point outside the allowed data
  folders are not shown or opened.
- **Project links** — a `?project=` URL outside the allowed data folders is
  refused.
- **Run status updates itself** every 10 s while a run is active (only the
  status section; the project is not re-validated) and stops when the run
  ends. A running run shows how to stop it safely from a terminal.
- **Result tables sort numerically** when you click a column header
  (baseMean, log2FoldChange, p-values, NES, …); identifiers stay text, and
  downloads are the original files.
- Disabled *Run*, *Save contrasts* and *Add* say why; the *Other* reference
  fields appear at once without an early error; an outdated full validation
  report is no longer shown; FASTQ copying says that it blocks the page and
  reports what was copied.
- **Clearer backend names** — New Project, Project overview and the project
  list show *Salmon (pseudoalignment)* (the default) and *HISAT2 +
  featureCounts (genome alignment)*. Project files are unchanged.

## Figures (standard preset)

- **Volcano and MA plots** have concise axis titles; the statistics, the
  display-only shrinkage and the DEG rule move to an external caption
  (`figure_caption`, shown under the figure in the HTML report).
- **DEG heatmap samples are no longer clustered**: they are ordered by pair
  (paired designs) or by group, then sample-sheet order. Genes stay clustered;
  `heatmap_zscores.tsv` is unchanged.
- **Enrichment dotplots** (GO/KEGG ORA and GSEA) use one blue scale on
  −log₁₀(adjusted P) with per-figure limits, so significant terms no longer
  share one colour. NES is never coloured with it.
- **L1 figures** gain a 300 dpi PNG; with a custom formula such as
  `~ condition + batch` they are coloured by `condition`, as the control plane
  groups samples. The legacy preset is unchanged.

## HISAT2 storage and upstream reuse

- Lane SAMs (≈242 GB for an 8-sample, ~30 M-pair project) and fastp
  prepared FASTQs are no longer published; only the fastp and HISAT2 reports
  are. BAMs, counts, strandedness evidence and MultiQC are published as before.
- `--reuse-upstream` and retries copy only what downstream and delivery read
  and reference the source run's alignment files in
  `referenced_alignment_artifacts.json`. No links are created; the copy is
  staged and verified before it is moved into place.
- L1 → L2 reuse no longer fails on checksum bookkeeping (`sha256_source`); the
  SHA-256 itself must still match.

## Stricter validation and safer run files (2026-11-04 audit)

- Unsafe FASTQ sample IDs and contrast IDs, literal `NA` identifiers and
  design values that R would read as missing, non-ASCII or underscore-grouped
  numbers, and a merged FASTQ next to the per-lane files of the same sample
  are rejected before anything runs.
- `run_state.json` and control-plane files are written atomically; the
  delivered run state is the final `SUCCESS` state.

## Install

```bash
python -m pip install "nf-rna[gui] @ git+https://github.com/2002brian/nf-rna.git@v1.7.0"
rnaseq --version          # 1.7.0
rnaseq gui                # http://127.0.0.1:8501; from a laptop: ssh -N -L 8501:127.0.0.1:8501 user@server
```

Optional, for TIFF 600 dpi only:

```bash
sudo apt install poppler-utils           # or: conda install -c conda-forge poppler
```

Without poppler, TIFF 600 dpi is shown as unavailable with this command; PDF,
TIFF 300 dpi and every other page work. Runtime requirements are unchanged
(Linux x86-64 or WSL2, Nextflow, Conda; see `docs/runtime.md`). The first
v1.7.0 run provisions a fresh downstream environment from the unchanged lock.

## Compatibility

- Result tables and statistics are unchanged. Figure images and their records
  change (above); the HISAT2 route publishes fewer files and its handoff
  manifest records `alignment.storage` and `reuse.source_pipeline`. Projects
  that use IDs outside the safe character set must rename the FASTQ files
  together with the metadata `sample_id`, or the contrast IDs. New runs record
  `nf_rna_version` 1.7.0.
- **`--reuse-upstream` does not accept a v1.6.1 (or earlier) source run**: the
  frozen upstream identity includes the nf-rna version, as in every release.
  Run the upstream stage once with v1.7.0. The GUI lists such runs as not
  reusable, with this reason.
- GUI pages renamed: *Dashboard* → *Projects*, *Run* → *Validate & Run*. The
  CLI is unchanged.

## Known limitations

- **PDF and raster layouts differ slightly** (already so in v1.6.1 output; not
  changed): R's PDF device measures text slightly differently from its PNG/TIFF
  devices, so panel edges sit 1–3 pt apart, up to ≈10 pt (3.6 mm) on enrichment
  dotplots with long term names. Data, colours, labels and order are identical.
  TIFF 600 dpi follows the PDF layout; TIFF 300 dpi follows the raster layout.
- FASTQ assigned to a project is copied into it (no "use in place" mode); the
  origin of an imported count matrix is not recorded.
- No TIFF above 600 dpi, no SVG, no figure editing. Under the `elsevier`
  profile both TIFF options come from the PDF (tested with fixtures only).
- Single user on 127.0.0.1, one active GUI run. SSH-tunnel access is documented
  but was not exercised in the v1.7.0 acceptance. The FAILED state was verified
  with a synthetic fixture.
- **No stop button** in the GUI: `rnaseq run` stops cleanly on SIGTERM only
  during a Nextflow stage, and stopping while the downstream runtime is being
  created would leave it half-built. A running run shows the `kill -TERM`
  command and when not to use it.
- FASTQ copying shows no progress bar (the importer has none). Result-table
  numbers are displayed with 6 significant digits; the TSV files are
  unchanged.
- **Interrupted upstream reuse**: a reuse interrupted while copying ends
  `INTERRUPTED`/`FAILED` with an empty `upstream/`. Do not `rnaseq retry` a run
  interrupted before its upstream handoff was frozen (it would recompute the
  whole upstream); start a new `--reuse-upstream` run.
- **KEGG snapshot**: an unpinned project uses the latest local snapshot, so a
  newer snapshot changes KEGG results (`rnaseq validate` warns). Pin
  `annotation.enrichment.kegg.snapshot` for published results.
- **FASTQ identity caching**: a FASTQ's SHA-256 is reused while its resolved
  path, size and mtime are unchanged; an in-place rewrite that keeps size and
  mtime keeps the old digest. Replace FASTQs by new files.
- **Low-depth VST**: a project with ≥ 1,000 retained genes but fewer than
  1,000 with mean normalized count > 5 stops in L1 (DESeq2 `vst()` `nsub`
  error); no result is produced.
- The scientific limitations of v1.6.1 are unchanged (CHANGELOG 1.6.1, Known
  limitations).
