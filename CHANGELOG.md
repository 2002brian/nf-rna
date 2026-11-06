# Changelog

All notable changes to this project are documented here.

## 1.7.0 — 2026-11-06

**Local Research Workspace & Publication Figure Export.** v1.7.0 changes the optional browser GUI, how the L1 and L2 figures are drawn, which HISAT2-route files are published and copied on reuse, and makes input validation stricter (2026-11-04 audit remediation). DESeq2, tximport, LFC shrinkage, GO/KEGG ORA, GSEA, Salmon, HISAT2 + featureCounts commands and parameters, the scientific thresholds, the reference and annotation handling and both downstream Conda locks are unchanged from v1.6.1 (`git diff v1.6.1 -- workflow/envs` is empty). In `src/rnaseq/r` only figure code, figure records and export formats changed (`figures.R`, the figure sections of `l1_analysis.R`, `l2_analysis.R`, `go_analysis.R`, `kegg_analysis.R`); `workflow/hisat2_featurecounts.nf` changes only its `publishDir` patterns and a sample-ID guard. Every numerical L1/L2 table is byte-identical to v1.6.1 for projects that pass validation; v1.6.1 results do not need to be recomputed. Figure images, figure records (`figure_caption`, `axes`, `colour`), the HISAT2 handoff manifest and the set of published HISAT2 files do change (sections below). **Breaking:** FASTQ sample IDs and contrast IDs are restricted to a safe character set. Release readiness: `docs/dev/release_readiness_v170.md`.

### Changed — L2 and enrichment figures, standard preset (`docs/figure-standard.md`)

Figure-only; no table, statistic, threshold, gene selection, clustering of genes or enrichment value changes, and the legacy preset is unchanged.

- **Volcano / MA: concise axes, details in an external caption.** Volcano x `log₂ fold change`, y `−log₁₀(adjusted P value)`; MA x `log₁₀(mean normalized count)`, y `log₂ fold change`. No caption, method name, BH abbreviation, contrast sentence or DEG rule is drawn in the image. The plotted data are unchanged (`log2FoldChange_shrunk`, or the unshrunken value when shrinkage failed; DESeq2 `padj`; log₁₀ `baseMean`). Each record gains `figure_caption` (plotted statistics, shrinkage being display-only, the unshrunken DEG rule, contrast direction) and `axes` with the source columns; the HTML report shows the caption under the figure.
- **DEG heatmap samples are no longer clustered.** Genes keep their clustering; samples are drawn by pair (first appearance in the sample sheet) then reference/numerator level for paired designs, and by level then sample-sheet order otherwise. `heatmap_zscores.tsv` is unchanged. Heatmap text follows the theme like the correlation heatmap; the order rule is also in `figure_caption`.
- **Enrichment dotplots use a single-hue blue adjusted-P scale** (`nf-rna Standard Enrichment Scale v3 (Blues)`; darker = smaller adjusted P, lightest `#9ECAE1` stays visible on white) in GO/KEGG ORA and GO/KEGG GSEA. Colour is mapped on −log₁₀(adjusted P) with per-figure limits (the smallest drawn adjusted P rounded down to a decade, to the cutoff) and decade ticks, so significant terms spread over the ramp instead of sharing one colour; colourbar `Adjusted P value`. ORA dotplots now also colour by adjusted P (redundant with their unchanged x axis). ORA x is `−log₁₀(adjusted P value)` with the database in the subtitle; GSEA x is `Normalized enrichment score (Tumor vs Normal)`. NES is never coloured with this scale. Transformation, limits, breaks and clipping are recorded per figure (`colour`, `figure_caption`).

### Fixed — HISAT2 artifact lifecycle and upstream reuse (`docs/dev/hisat2_artifact_reuse_v170.md`)

- **Lane SAM files are no longer published.** `HISAT2_ALIGN` published every uncompressed lane SAM to `upstream/hisat2_featurecounts/alignment/lane_summaries/` (about 242 GB for an 8-sample, ~30 M-pair project), although nothing reads them after `SORT_LANE_BAM`. Its `publishDir` now has `pattern: '*.hisat2.summary'`. The SAM is still the channel input of `SORT_LANE_BAM` and stays in the task work directory, which `rnaseq clean` removes. The BAMs, BAI files, flagstat, count-only BAMs, featureCounts tables and summaries, strandedness evidence, count matrix and MultiQC report are published as before. No alignment or counting command, parameter or output value changes. The SHA-256 of `workflow/hisat2_featurecounts.nf`, recorded as provenance, changes; it is not part of the `--reuse-upstream` comparison.
- **Prepared FASTQs are no longer published.** `FASTP_PREPARE` had the same unfiltered `publishDir`, so `qc/fastp/prepared/*.fastq.gz` was published as well. For pretrimmed input these are byte copies of the source FASTQs (26 GB for NGS1130901); for raw input they are fastp's trimmed reads. Only `fastp_*.html` and `fastp_*.json` are published now. The prepared reads remain the input of `HISAT2_ALIGN` and `FASTQC_PROCESSED`, and a reuse references them in the source run like the alignment artifacts.
- **`--reuse-upstream` and retries no longer copy alignment artifacts.** Both copied the source run's whole `upstream/` tree with `shutil.copytree`, so an L2 freeze re-read and re-wrote every SAM and BAM (observed: 38+ minutes and ~198 GB read plus ~198 GB written before any analysis). The new run now gets its own copy of every file downstream and delivery read: the count matrix, sample map, featureCounts tables and summaries, strandedness evidence, MultiQC, FastQC/fastp reports, HISAT2 summaries and flagstat.
  - The HISAT2 route's `*.sam`, `*.bam`, `*.bai`, `*.cram`, `*.crai` and `*.csi`, plus `qc/fastp/prepared/` FASTQs, stay in the run that produced them. They are listed with path, size and owning run in `upstream/hisat2_featurecounts/referenced_alignment_artifacts.json`.
  - The handoff manifest records `alignment.storage: referenced_in_source_run` and `reuse.source_pipeline`, the workflow identity that produced the outputs.
  - No hard links or symbolic links are created. The copy is staged and moved into place only after every file's size has been verified, so a failed or interrupted copy leaves an empty `upstream/` and a `FAILED` run.
  - Reuse compatibility checks are unchanged and run before anything is copied.
  - Salmon (`nfcore_rnaseq`) upstreams are copied completely, as before.
  - L1 runs made before this fix, including those with lane SAMs and prepared FASTQs, are reused the same way. Those files are referenced and never copied.

### Fixed — 2026-11-04 audit remediation (`docs/dev/audit/2026-11-04-cloud-audit.zh-TW.md`, `docs/dev/audit/2026-11-04-v170-remediation.zh-TW.md`)

Validation is stricter in the cases below; projects that pass were never affected. Regression tests: `tests/audit/`.

- **FASTQ sample IDs are restricted (D-06, security).** A sample ID taken from a FASTQ file name must match `^[A-Za-z0-9][A-Za-z0-9._-]*$` (new error `unsafe_fastq_sample_id`). Sample IDs are interpolated into the HISAT2 + featureCounts task scripts and output names, so a file such as `a$(id)_R1.fastq.gz` could run a command and `a b_R1.fastq.gz` failed only after alignment. **Breaking:** FASTQ names with spaces, non-ASCII letters or shell characters must be renamed (together with the metadata `sample_id`). The HISAT2 workflow also refuses such IDs in its samplesheet before any task starts.
- **Contrast IDs are restricted (D-02, path traversal).** `contrast_id` becomes the L2 output directory `contrasts/<contrast_id>`; hand-written or imported IDs such as `../../escape`, `a/b` or `with space` are now rejected (`unsafe_contrast_id`) with the rule `rnaseq design` already used for the IDs it generates.
- **L1 → L2 `--reuse-upstream` no longer fails on checksum bookkeeping (D-01).** The upstream identity compared the per-file `sha256_source` (`computed` in the run that hashed a FASTQ, `cache` in later runs), so the usual L2 reuse of an L1 upstream was refused although every path, size, mtime and SHA-256 matched. That field is now left out of the comparison; everything else, including the SHA-256, must still match. The compatibility rules are otherwise unchanged. The source must be a `SUCCESS` run of the same project with an upstream handoff manifest, and these frozen files must match: `nfcore.params.json`, `nfcore.local.config` and `samplesheet.csv` byte for byte, and in `input_manifest.yaml` the nf-rna pipeline name and version plus the `input` (per FASTQ: path, resolved path, size, mtime and SHA-256), `upstream` and `reference` sections. Design, metadata, contrasts and analysis scope may differ. The workflow SHA-256s (`workflow/main.nf`, `workflow/hisat2_featurecounts.nf`) are provenance in `run_provenance.yaml`, `upstream_run.yaml` and the HISAT2 handoff manifest; they are not part of this comparison. Salmon and HISAT2 v1.7.0 L1 runs made before this fix, including HISAT2 runs made before the workflow's sample-ID guard, can therefore be reused when the rules above hold. The guard only refuses sample IDs that validation never writes, so it does not change any output.
- **Merged and per-lane FASTQs of one sample are refused (R-03).** A lane-less FASTQ (`S1_R1.fastq.gz`) next to lane files of the same sample (`S1_L001_R1_001.fastq.gz`) is almost always their concatenation and was aligned and merged as well, doubling the sample's reads (new error `mixed_merged_and_lane_fastq`).
- **Values R reads differently are rejected (D-03, D-04, D-05).** R's `read.csv` reads a literal `NA` as missing: `NA` as a metadata `sample_id` (`r_missing_sample_id`), as a categorical design value (`r_missing_design_value`), in `pair_id` (`r_missing_pair_id`) or as a count-matrix `gene_id` (`r_missing_gene_id`) is now an error. Counts and continuous design values must use plain ASCII decimal syntax (e.g. `40`, `4e1`, `40.0`); PEP 515 underscores (`1_000`), non-ASCII digits and non-breaking spaces, which Python accepted and R turns into `NA`, are rejected (`non_numeric_count`, `invalid_continuous_value`).
- **Run state and control-plane files are written atomically (R-01).** `run_state.json`, frozen contracts and manifests are written to a temporary file, flushed and renamed, so an interruption or a polling reader (`rnaseq status --watch`, GUI) never sees a truncated file. An unreadable `run_state.json` is kept as `run_state.json.corrupt-<time>` and logged instead of being silently discarded.
- **The delivered run state is final (D-07).** `delivery/methods_and_versions/run_state_<date>.json` was copied while the run was still `RUNNING` (`completed_at: null`); it is now the `SUCCESS` state recorded once the package verifies, with the same `completed_at`.
- **L1 figures group by the same factor as the control plane (D-08).** With a custom formula such as `~ condition + batch`, the L1 PCA, boxplots and correlation annotation coloured samples by the last term (`batch`); they now use `condition` whenever the formula models it, as `rnaseq.validators.condition_factor` does. Default formulas (condition last) are unchanged. No statistic depends on this choice.

### Added — Local Research Workspace (GUI)

- **Local Dataset Browser** (*Datasets* page). Browses folders inside the allowed roots of `gui.yaml` only; roots on an unmounted disk are shown as unavailable and never browsed. Shows file count, total size, unreadable files and file kinds (FASTQ, count matrix, metadata, samplesheet, contrasts, other), recognised from the file name and, for tables, the header line; FASTQ contents are not read. FASTQ and count-matrix samples are recognised with the validator's existing rules. Symlinks that resolve outside the allowed roots are listed as ignored and never followed. A dataset is assigned to a project with the existing importers: FASTQ files are copied into the project after a path and free-space check; the source files are never moved or changed. There is no central FASTQ database.
- **Project Workspace** (*Project overview* page; the dashboard becomes *Projects*). Species, input type, quantification backend, reference, sample count, contrasts, L1/L2 settings, the original FASTQ folder of imported FASTQ, the workflow steps Dataset → Project → Metadata → Contrasts → Validate → Run → Results with their state and the next step, and buttons that open a run's results or report. Projects can be switched from the sidebar; each project keeps its own run selection and case ID.
- **Run History.** Every run of a project with level, start and end time, execution status, scientific QC and whether results and a report exist; the *Validate & Run* page (formerly *Run*) adds the same history, a settings summary table, validation findings grouped by workflow step, the runs that cannot be reused with the preflight's own reason, and a log filter for warnings and errors. *Run* is disabled while the project is invalid.
- **QC display.** Execution status (RUNNING, SUCCESS, FAILED, INCOMPLETE) is shown separately from scientific QC (PASS, WARNING, INCONCLUSIVE, NOT_EVALUATED). Scientific QC is a display roll-up of verdicts the run itself recorded: the v1.6.1 strandedness check, the L1 sample flags, the DESeq2/apeglm/ashr R warnings and the shrinkage status, each with its recorded value. A QC warning never changes the execution status or `run_state.json`, and no verdict is re-judged. Runs that predate a check are NOT_EVALUATED, not PASS.

### Added — PDF/TIFF figure downloads (GUI)

- *Results* offers every figure as **PNG, PDF, TIFF 300 dpi and TIFF 600 dpi**; *Reports & Downloads* gains a **Figures (PDF / TIFF)** section covering every figure of a run (PCA, sample correlation, boxplots, volcano, MA plot, DEG heatmap, GO/KEGG ORA and GSEA dotplots, GSEA running-score plots), with TIFF 300 dpi as the default.
- Nothing is redrawn from data and nothing is written into a run. PDF and the 300 dpi TIFF are the pipeline's own files, byte for byte. GSEA running-score plots, which have no pipeline TIFF, get a lossless TIFF of their R-drawn 300 dpi PNG (every pixel unchanged). TIFF 600 dpi is the pipeline's vector PDF rasterized with poppler at exactly 600 dpi. Every TIFF is RGB on white, LZW, with the dpi in its resolution tags and `round(page_pt / 72 × dpi)` pixels, the rule of the pipeline's own devices. A raster is never enlarged; a figure or format that cannot be produced is shown as unavailable with the reason.
- **Poppler** (`pdftocairo`, `pdfinfo`) is needed only for TIFF 600 dpi and is not installed by nf-rna. Without it the option is shown as unavailable with the install command, the *Projects* page says so, and every other download and page works.

### Changed — L1 figures (`docs/dev/gui_ux_v170.md`)

- Every L1 figure (PCA, PCA by batch, sample correlation, both sample boxplots) also writes `<name>_300dpi.png` (the profile's raster dpi) beside the 150 dpi preview, as the L2 figures have since v1.5.0. Existing file names, the TIFF (LZW, profile dpi) and the PDF are unchanged.
- Standard preset, sample correlation heatmap: the text now follows the theme of the ggplot figures (PCA, boxplots): title at the base size in bold instead of pheatmap's 1.3 × base, annotation names and legend titles plain, legend labels at the small size. Values, colours, cell size and the clustering order are unchanged. Legacy-preset figures are unchanged (PNG/TIFF byte-identical, PDF pixel-identical).
- The correlation matrix, PCA scores and variance, VST, normalized counts, distribution statistics, filtering and QC flags are byte-identical before and after (checked on the 8-sample NGS1130901 data).

### Changed — GUI responsiveness and workflow clarity (`docs/dev/gui_ux_v170.md`)

- **Much faster interactions.** Each Streamlit interaction re-ran the project validation, which re-computed the SHA-256 of every multi-GB reference file (2.7 s per validation for the Human reference; 9.3 s per interaction on New Project). The GUI process now reuses the digest of an unchanged large file (≥ 16 MiB; identity = device, inode, size, mtime, ctime); a changed, replaced or touched file is hashed again. The memo exists only inside the GUI process: `rnaseq` commands and every pipeline run still hash every file in full. Median rerun on the 8-sample Human project: New Project 9.3 s → 0.06 s, Contrasts 5.2 s → 0.03 s, Projects 4.7 s → 0.03 s, Project overview / Input & Metadata / Validate & Run 2.6 s → 0.02–0.07 s. The first validation after the GUI starts still hashes once and says so.
- **Workflow bar** on Input & Metadata, Contrasts, Validate & Run and Results: the steps Dataset → Project → Metadata → Contrasts → Validate → Run → Results with their state, the steps done on the current page, and the next step with a link.
- **New Project**: FASTQ preprocessing (raw / pretrimmed) is shown with the other FASTQ settings instead of under *Advanced settings* (default unchanged); one-line explanations of Salmon versus HISAT2 + featureCounts and of unpaired / paired / multi-group designs; the DESeq2 formula and metadata columns the project will get (from the same function that creates it); and, before *Create project*, what is still missing or which folder will be written.
- **Input & Metadata**: FASTQ detection and samplesheet selection are labelled as previews (nothing is copied until *Copy* / *Import*); a live **pairing check** for paired designs shows, for the table as edited (saved or not), each pair with its samples per condition and what is missing or duplicated, by the validator's rule.
- **Validate & Run**: what the case ID is and that a re-used case ID never overwrites a run, the case IDs used before, which actions write files (*Validate and refresh the plan* writes `planning/`; validation is otherwise read-only), and how L1 relates to L2.
- **Results**: L1 figures are grouped (*Sample structure*, *Expression distributions*) and shown at one common on-screen scale, so text has the same size in every figure; captions give the physical size; the downloads of a figure are in one *Download* menu: PDF (vector), **PNG 300 dpi** (new; the pipeline's file, or the PDF rendered at 300 dpi for runs made before v1.7.0), TIFF 300/600 dpi and the PNG preview. L2 figures use the same menu.

### Fixed (GUI)

- A project opened by URL (bookmark, reload) was lost when changing page.
- Run selectors used global widget keys, so a run selection leaked between projects.
- The *Validate & Run* reuse hint offered SUCCESS runs made by another nf-rna version, which the run then refused after it had started; they are now listed as not reusable with the reason.
- GSEA running-score listings treated a `_<n>dpi.png` raster of a non-300 dpi profile (`elsevier`, 500 dpi) as a separate plot.
- Pre-release hardening (`docs/dev/gui_hardening_v170.md`):
  - *Use this folder* / *Use this file* raised a Streamlit error. The browser no longer lists or opens links that resolve outside the allowed roots.
  - `?project=` opened any local project; it must now be inside an available allowed root.
  - Metadata: bulk assignment, *Paste table*, *Import CSV* and the Datasets metadata load discarded unsaved grid edits, and paste/import replaced the whole table. They now act on the table on screen and merge by `sample_id`. Replacing non-empty values needs a tick, an *Unsaved changes* banner compares the table with metadata.csv, unsaved edits survive page switches, and *Reload* needs a tick.
  - Result tables sorted lexically by column header; numeric columns are now numeric for display, with `NA` last. Identifiers and ratios stay text, and downloads are unchanged.
  - The run section of *Validate & Run* refreshes itself every 10 s while a run is active, without re-validating the project, and stops at a final state. It shows how to stop a running run from a terminal.
  - New Project: the reference choice is kept per species and backend, so "Other", chosen while Human had no registered reference, no longer hid the registered Mouse reference after a switch. A species/backend without a registered reference says "reference not prepared" with the registration or preparation command. A managed reference is re-checked for species and backend before the project is created (`docs/dev/human_hisat2_reference_v170.md`; the missing Human entry was a workstation registry gap, fixed with `rnaseq reference register`).
  - Quantification backends are labelled *Salmon (pseudoalignment)* and *HISAT2 + featureCounts (genome alignment)* in New Project, Project overview and the project list (previously the wizard's descriptions in New Project and *Salmon (+ tximport)* / *HISAT2 + featureCounts* elsewhere). Salmon stays the default; `project.yaml` keeps `salmon` / `hisat2_featurecounts`.
  - The GUI's default settings named one workstation's data disk as the first allowed root and the project location. A first start now writes the portable default workspace `~/nf-rna-work` (allowed root, with `projects/` and `fastq/`, created on first use) for the current user. An existing `gui.yaml` is never rewritten and keeps precedence, so installations with their own folders need no change. `projects/` and `fastq/` are no longer created for a `project_parent` outside the allowed roots.
  - Disabled *Run*, *Save contrasts* and *Add* state their reason. *Other (advanced)* reference fields show at once without a premature error or CLI flag names. An outdated full validation report is removed. The Contrasts page redraws after saving. FASTQ copying warns that it blocks the page and reports the result.

### Compatibility

- No analysis value or result table changes. Figure images and their records change (L1/L2 figure sections above); the HISAT2 route no longer publishes lane SAMs or prepared FASTQs, and its handoff manifest records `alignment.storage` and `reuse.source_pipeline`. Projects with FASTQ sample IDs or contrast IDs outside the new safe character set, literal `NA` identifiers, or merged plus per-lane FASTQs of one sample now fail validation. New runs record `nf_rna_version` `1.7.0`.
- Local AI-agent development assets (`.claude/skills/`, `.agents/skills/`) are no longer tracked or shipped: the sdist no longer includes `.claude/skills`, and `.claude/` and `.agents/` are git-ignored. The figure engine, `docs/figure-standard.md` and `tools/figure_review/` are unchanged and still shipped.
- **Upstream reuse across versions.** The frozen upstream identity that `rnaseq run --reuse-upstream` compares includes the nf-rna version (`input_manifest.yaml` `pipeline.version`), so a v1.7.0 run refuses a source run made by v1.6.1 or earlier (`frozen input_manifest.yaml differs in its input, upstream, or reference identity`), as every earlier release did across versions. Rerun the upstream stage with v1.7.0 once; reuse between v1.7.0 runs works as before. The GUI's reuse hint now lists such runs as not reusable with this reason instead of offering them.
- The runtime environment is keyed by lock and nf-rna build, so the first v1.7.0 run provisions a fresh downstream environment from the unchanged lock.
- GUI page names changed: *Dashboard* → *Projects*, *Run* → *Validate & Run*; new pages *Datasets* and *Project overview*. The CLI is unchanged.

### Known limitations

- **PDF and raster figure layouts differ slightly** (pre-existing in v1.6.1 output, not changed). R's PDF device measures text slightly differently from its PNG/TIFF devices, so a figure's panel edge can sit 1–3 pt apart between its PDF and its PNG/TIFF, and up to about 10 pt (≈3.6 mm) on enrichment dotplots with long term names. Data, colours, labels and order are identical. The TIFF 600 dpi download follows the PDF's layout; the TIFF 300 dpi download follows the raster layout. The R plotting engine is not changed in this release.
- FASTQ assigned to a project is copied into it (the project requires `input.path` inside the project); there is no "use in place" mode. The original location of an imported count matrix is not recorded, so it is not shown.
- No TIFF above 600 dpi, no SVG, no figure editing or presets. Under the `elsevier` profile (500 dpi pipeline TIFF) both TIFF options are rendered from the PDF; this was tested with fixtures only.
- The GUI remains single-user on 127.0.0.1 with one active GUI run. Access from another computer through an SSH tunnel is documented but was not exercised in the v1.7.0 acceptance. `/dashboard` is not a URL; the *Projects* page is the server root.
- The FAILED execution state was verified with a synthetic fixture project; no real failed run was available.
- The GUI has no stop button: `rnaseq run` stops cleanly on SIGTERM only during a Nextflow stage, and stopping while the downstream runtime is being created leaves it half-built. A running run shows `kill -TERM <pid>` with these caveats. FASTQ copying shows no progress bar.
- GUI limitations of v1.6.0 not listed above are unchanged (`docs/gui.md`).
- **Interrupted upstream reuse.** A `--reuse-upstream` run interrupted while it copies the source upstream ends `INTERRUPTED` (or `FAILED`) with an empty `upstream/`; nothing partial is used. Do not `rnaseq retry` a run interrupted before its upstream handoff was frozen (its contract has `upstream_handoff: null`): the retry is treated as a run without reuse and recomputes the whole upstream. Start a new run with `--reuse-upstream` instead. HISAT2 runs on a custom (unmanaged) reference whose index is a directory cannot be retried (unchanged).
- **KEGG snapshot.** KEGG gene sets come from a local snapshot under `~/.cache/nf-rna/genesets/kegg/`. An unpinned project uses the latest local snapshot (or downloads one on the first KEGG run), so a newer snapshot changes KEGG results; `rnaseq validate` warns `kegg_snapshot_not_pinned`. Pin `annotation.enrichment.kegg.snapshot` for published results (unchanged since v1.5.0).
- **FASTQ identity caching.** A FASTQ's SHA-256 is computed once and reused from the user cache (`fastq-sha256/checksums.json`) while its resolved path, size and mtime (ns) are unchanged. A file rewritten in place with the same size and a restored mtime would keep its old digest; replace FASTQs by new files, not by in-place edits that preserve mtime.
- **Low-depth VST.** L1 uses DESeq2 `vst()` when 1,000 or more genes are retained and the exact `varianceStabilizingTransformation()` below that. `vst()` also needs 1,000 genes with mean normalized count > 5, so a low-depth project with ≥ 1,000 retained genes but fewer than 1,000 above that level stops in L1 with `less than 'nsub' rows with mean normalized count > 5` (`docs/dev/scientific_ci_mutations.md`). The run fails; no result is produced. The fallback condition is not changed in this release.

## 1.6.1 — 2026-10-31

**Result-changing.** Three corrections found by the independent scientific validation of v1.6.0 (`docs/dev/scientific_corrections_v161.md`) change GSEA leading edges and GSEA Monte Carlo results, and DESeq2 independent filtering (padj, DEG sets, ORA). Existing v1.6.0 runs are not rewritten; recompute a project to obtain corrected results. Three further changes add QC, provenance and warning records and change no analysis value (`docs/dev/scientific_reliability_v161.md`). Release readiness, evidence and limitations: `docs/dev/release_readiness_v161.md`.

### Fixed

- **GSEA leading edge (GO and KEGG).** DOSE 4.4.0 `leading_edge()` locates the running-sum minimum of an ES < 0 gene set among post-step scores at hits only and can report the wrong boundary. nf-rna now recomputes `core_enrichment`, `rank` and the `leading_edge` tag string by the standard running-sum definition (Subramanian et al. 2005) and checks the recomputed extremum against the fgsea ES. Affects `core_enrichment`, the `leading_edge` tags/signal, `gsea_core_members.tsv`/`core_members.tsv` and the `leading_edge` flag of `members.tsv` for ES < 0 sets only; ES, NES, p-values and FDR are untouched by this fix.
- **GSEA seed.** GSEA was called with `seed = TRUE`, which makes DOSE 4.4.0 re-seed with `set.seed(.Random.seed)`, i.e. the fixed state 10403, whatever `annotation.enrichment.gsea.seed` / `kegg.gsea.seed` said. nf-rna now calls `set.seed(seed)` and then GSEA with `seed = FALSE`, so the configured seed determines the fgsea Monte Carlo draw. NES, p-values, FDR and borderline significant sets change within Monte Carlo error; the recorded seed is now the seed actually used, and the provenance records the RNG kind.
- **DESeq2 independent filtering alpha.** `results()` used DESeq2's default `alpha = 0.1` while DEGs are called at `padj < thresholds.padj` (default 0.05). `alpha` is now `thresholds.padj`, as the DESeq2 vignette recommends. The filter threshold, which genes have `padj = NA`, padj values, DEG sets and the ORA universe/foreground can change; `pvalue`, `stat` and LFCs do not. `backend_summary.json` and `contrast_metadata.json` record the alpha and the chosen baseMean filter threshold. `thresholds.padj` must now lie strictly between 0 and 1.

### Added (QC and provenance; analysis values unchanged)

Findings of the Phase 3 upstream FASTQ validation (`docs/dev/scientific_reliability_v161.md`).

- **Strandedness consistency check (HISAT2 + featureCounts).** A wrongly declared strandedness used to give a successful run with wrong counts and no warning. Each count-only BAM is now also counted with featureCounts `-s 1` and `-s 2`, and the declaration is judged per sample:
  - `CONSISTENT`;
  - `CONTRADICTS_DECLARATION` or `APPEARS_UNSTRANDED` (warning);
  - `APPEARS_STRANDED` (advisory);
  - `INCONCLUSIVE` (fewer than 1,000 strand-assignable fragments, or no clear direction).

  The result is shown in the run log, the report and the GUI, and is recorded in provenance. The check never changes the declared strandedness or the counts, and never fails a run.
- **Actual upstream thread counts.** Nextflow traces of FASTQ runs gain `process`, `tag`, `cpus` and `memory` columns. `run_provenance.yaml` (`execution_threads`), the report and the GUI record the CPUs each relevant upstream process actually used (HISAT2, featureCounts, fastp, Trim Galore, Salmon), Salmon's `cmd_info.json` threads and the tool versions per process, with a note that alignment and quantification can depend on the thread count. CPU allocation is unchanged.
- **Recorded R warnings in L2.** DESeq2, apeglm and ashr warnings are recorded verbatim, by step and with counts, in `contrast_metadata.json` (`r_warnings`), `backend_summary.json` and the L2 provenance, and shown in the report and the GUI. They are not muffled, not escalated, and do not change any result.

### Compatibility

- **Results.** GO/KEGG GSEA (leading edge of ES < 0 sets; NES, p-values, FDR and borderline significance through the effective seed) and DESeq2 `padj`, DEG sets and ORA differ from v1.6.0. Do not mix v1.6.0 and v1.6.1 results in one analysis; recompute.
- **Upstream reuse across versions.** The frozen `nfcore.local.config` of a FASTQ run now ends with the `trace.fields` line, so `rnaseq run --reuse-upstream` refuses a source run made by v1.6.0 or earlier (`frozen nfcore.local.config differs`). The error's hint names only CPU/memory limits; the cause here is the version change. Rerun the upstream stage with v1.6.1. Reuse between v1.6.1 runs is unchanged.
- New outputs only: `qc/strandedness/` on the HISAT2 route, and new keys in the handoff manifest, frozen contract, `run_provenance.yaml`, `contrast_metadata.json`, `backend_summary.json` and the L2 provenance. No existing key, table or column changes meaning.

### Runtime

- The downstream Conda lock and the nf-core/rnaseq 3.26.0 pin are unchanged. The runtime environment is keyed by lock and nf-rna build, so the first v1.6.1 run provisions a fresh downstream environment from the same lock.
- The HISAT2 route runs one extra process per sample (`STRAND_EVIDENCE`: two featureCounts passes, 2 CPUs).

### Known limitations

- **Thread-dependent upstream results.** HISAT2 2.2.3 places multimapping reads differently at different `-p`; Salmon's EM gives different `NumReads` at different thread counts. Thread counts are recorded, not fixed; `execution.max_cpus` can lower them.
- **Salmon is not bitwise reproducible at the same thread count.** Rerunning Salmon 1.10.3 on identical input with identical arguments (`--threads 6`) gave identical `NumReads` but different `EffectiveLength` (105–117 transcripts, through the fragment-length estimate) and therefore `TPM`. Through the tximport length offsets, L1 normalized counts differed by at most 1.03e-3 (relative) and VST by 6.4e-5 on the yeast test data; the effect on DE was not measured. Byte-identical Salmon-route reruns cannot be expected.
- The strandedness check covers the HISAT2 route only (Salmon route: `NOT_EVALUATED`), warns but never fails, and its limits were calibrated on one real stranded library (yeast, reverse) plus unstranded human and simulated data.
- Not validated with natural data: forward-stranded and single-end libraries (simulation only), mouse FASTQ, multi-lane merging, preprocessed/`--skip_trimming` input, other aligners, quantifiers and annotation releases. `auto` strandedness near nf-core's 0.8 threshold was not explored. L1 R warnings (VST, PCA) are not recorded.
- The independent validations of v1.6.0 (Phase 2, downstream) and v1.6.1 (Phase 3, upstream) do not PASS under their frozen criteria; their FAIL verdicts stand. Their root causes are the defects fixed here (D-1, D-2, M-1), tool behaviour (thread dependence, Salmon's library-type classifier) and acceptance-criterion calibration; no nf-rna computation error remains open. See `docs/dev/release_readiness_v161.md`.

## 1.6.0 — 2026-10-30

**This release includes the v1.5.1 test fixes**: v1.5.1 was prepared but not tagged separately, so its changes ship in v1.6.0 (see [1.5.1](#151--2026-10-30)). **Upgrading from v1.4.x: the enrichment significance semantics changed in v1.5.0.** GO/KEGG ORA and GSEA call a term significant when the Benjamini–Hochberg adjusted p-value is ≤ 0.05, the ORA universe changed, and KEGG gene sets come from a local snapshot that you should pin for published results; read the [1.5.0 notes](#150--2026-10-30) before comparing results across versions. v1.4.x results are not recomputed and keep their v1.4 semantics.

v1.6.0 adds an optional local browser GUI. DESeq2, GO/KEGG ORA, GSEA, the significance semantics ("v1.5"), the downstream Conda lock, the Nextflow workflows, the R scripts and every result table and figure are unchanged from v1.5.0/v1.5.1 (an airway L2 rerun gives byte-identical tables, PNG and TIFF files; PDFs differ only in their embedded creation date); v1.5.x results do not need to be recomputed.

### Added — browser GUI (`rnaseq gui`, optional `gui` extra)

- `rnaseq gui` serves a Streamlit GUI on 127.0.0.1 (default port 8501, `--port`); use an SSH port forward from another machine. Install with the `gui` extra (`pip install -e '.[gui]'`, or `"nf-rna[gui] @ git+https://github.com/2002brian/nf-rna.git@v1.6.0"`). Without the extra the CLI is unchanged and `rnaseq gui` prints the install command.
- Pages: Dashboard (projects, recent runs, `rnaseq doctor`, registered references, mount state of the data roots), New Project (the `rnaseq new` wizard), Input & Metadata (count matrix or FASTQ folder import, an Excel-paste metadata grid that may be saved incomplete), Contrasts (the `rnaseq design` transaction, KEGG snapshot pinning), Run (validate, plan, background `rnaseq run` with one active GUI run, status and logs), Results and Reports & Downloads.
- The GUI calls the same services as the CLI: wizard choices, the `rnaseq design` transaction and the FASTQ/count import were moved into shared functions, with unchanged CLI prompts and messages. It keeps no database; run state is read from the run files.
- Runs started from the GUI are ordinary `rnaseq run --yes` runs, started as a transient `systemd --user` unit when the GUI runs under systemd or lingering is enabled, otherwise as a detached process; they survive closing the browser or restarting the GUI. A run whose process has died is shown as INCOMPLETE without changing its files.
- Results Explorer: QC and PCA, per-contrast DEG summary (counts and thresholds from `contrast_metadata.json`), volcano, MA plot and heatmap, a searchable DEG table showing `log2FoldChange` and `log2FoldChange_shrunk` with the shrinkage method and status, GO/KEGG ORA and GSEA tables (significant or all terms) with the member genes of a selected pathway, dotplots, simplified GO dotplots and running-score plots, the significance semantics, ORA universe sizes, GSEA ranking metric and seed, the KEGG snapshot and whether it was pinned, and provenance. Results produced before v1.5 are labelled legacy. Figures are the pipeline's own PNG files; nothing is recomputed. Missing outputs are shown as unavailable with the backend's reason; L1 runs show no L2 sections.
- Reports & Downloads: the run's `report.html` (viewed in the page or downloaded), DEG/count/QC/enrichment TSVs, PNG/PDF figures and the delivery package as a zip verified against `delivery_manifest.yaml`. Files above `download_limit_mb` (default 500 MB) are offered as a server path and `scp` command.
- `docs/gui.md`; an example `systemd --user` unit in `docs/systemd/` (never installed or enabled by nf-rna).

### Changed

- `report.html` gains a Samples table, a Reference and quantification section, the DEG rule and LFC shrinkage method/status per contrast, and a Software versions and provenance table (nf-rna version, source revision, wheel and lock SHA-256, R, DESeq2, apeglm, ashr, clusterProfiler, OrgDb, GO.db, KEGG snapshot). The report remains a single self-contained file with embedded figures; no table or figure changes.

### Runtime

- The downstream lock is unchanged. The runtime environment is keyed by lock and nf-rna build, so the first v1.6.0 run provisions a fresh downstream environment from the same lock.

### Known limitations

- The GUI is a single-user, local tool: one active GUI run at a time, loopback only, no accounts.
- Verified end to end on raw-count projects (airway, L2 with GO/KEGG ORA and GSEA). FASTQ projects: file recognition and import were verified in the GUI, but no FASTQ run was started or browsed from the GUI for this release; the MultiQC download and the FASTQ wording of the report's *Reference and quantification* section are covered by unit tests only.
- The `systemd` launcher was verified with `launcher: systemd` from a GUI started in a terminal; the GUI running as the example `systemd --user` service was not tested. With the detached-process launcher, killing `rnaseq` with SIGKILL leaves its Nextflow child running (a normal stop is clean).
- Display of pre-v1.5 (legacy) enrichment output was verified with test fixtures only.
- Copying FASTQ files into a project blocks that page until the copy finishes. Covariates other than `batch` and the pairing column, `rnaseq retry` and `rnaseq clean` remain CLI-only.
- Delivery zips prepared by the GUI stay in `~/.cache/nf-rna/gui/packages/` until removed by hand.

## 1.5.1 — 2026-10-30

**Not released separately**: v1.5.1 was prepared but never tagged; its changes ship in v1.6.0.

Maintenance release. It fixes how three R-dependent tests decide whether they can run; nothing that nf-rna installs or executes has changed. DESeq2, GO/KEGG ORA, GSEA, the enrichment significance semantics ("v1.5"), the downstream Conda lock, the Nextflow workflows and every result table are identical to v1.5.0. Results produced with v1.5.0 do not need to be recomputed.

### Fixed

- Tests: the KEGG snapshot-resolution test and the offline KEGG backend tests (`tests/test_kegg.py`, `tests/test_enrichment_mapping_qc.py`, `tests/test_figure_dotplots.py`) now require the R packages they load (`jsonlite`, `openssl`) through the shared `require_r_packages` gate, so they skip on a system R without `openssl` (or fail under `NF_RNA_REQUIRE_R=1`) instead of failing in CI. The `v1.5.0` tag predates this fix; it is unchanged and remains valid.

### Runtime

- The downstream lock is unchanged (`nf-rna-downstream-linux-64.lock.yml`). Because the runtime environment is keyed by lock and installed nf-rna build, the first v1.5.1 run provisions a fresh downstream environment from the same lock.

## 1.5.0 — 2026-10-30

Scientific-engine release: publication-ready DEG outputs and a new enrichment significance contract. **GO/KEGG ORA and GSEA now call a term significant when the Benjamini–Hochberg adjusted p-value (p.adjust) is ≤ 0.05.** KEGG gene sets come from a SHA-256-verified local snapshot; **pin `annotation.enrichment.kegg.snapshot` for results you will publish.** DESeq2 statistics and the DEG definition are unchanged from v1.4.1. Results produced with v1.4.x are not recomputed and keep their v1.4 semantics.

### Changed — enrichment significance semantics

- **Enrichment significance semantics have changed (semantics "v1.5").** GO ORA, KEGG ORA, GO GSEA and KEGG GSEA now all call a term significant when `p.adjust` (Benjamini–Hochberg) `<= 0.05`. clusterProfiler is called with `pvalueCutoff = 1` and `qvalueCutoff = 1`, and nf-rna adds a `significant` column to every all-terms table. In v1.4.x, ORA used `pvalue <= 0.05 AND qvalue <= 0.2` (p.adjust not used) and GSEA used `p.adjust <= 0.05 AND pvalue <= 0.05`. Existing results are not recomputed; new summaries and provenance record `significance_semantics: "v1.5"`. The `pvalue_cutoff`/`qvalue_cutoff`/`padj_cutoff` fields of `annotation.enrichment` are still accepted but no longer select terms.
- **ORA universe.** The ORA universe is now the genes with a non-NA DESeq2 padj, mapped to Entrez and intersected with the gene-set database (each GO ontology, or the KEGG snapshot). It was every gene with a finite p-value. Universe sizes, the Entrez mapping rate and clusterProfiler's BgRatio N are recorded per foreground.
- **KEGG gene sets come from a local snapshot.** KEGG REST is downloaded once (KEGGREST) into `~/.cache/nf-rna/genesets/kegg/<org>/<snapshot_id>/` (`TERM2GENE.tsv`, `TERM2NAME.tsv`, `manifest.json` with URLs, retrieval time, KEGG release information and SHA-256). KEGG ORA and GSEA use `enricher()`/`GSEA()` on that snapshot instead of `enrichKEGG`/`gseKEGG`. Without `annotation.enrichment.kegg.snapshot` the latest existing snapshot is used. Snapshot ID and SHA-256 are verified and recorded in provenance. The run preflight needs KEGG REST only when no snapshot exists. `resource_provider` defaults to `kegg_rest_snapshot`; the old value is accepted.
- **GSEA reproducibility.** Both GO and KEGG GSEA now use `set.seed(seed)` immediately before the call, `seed = TRUE`, `nPermSimple = 10000` and `eps = 0`, all recorded. KEGG GSEA previously had no explicit seed. The ranking metric is unchanged (DESeq2 Wald `stat`).
- **GO provenance** records the OrgDb, GO.db and clusterProfiler versions.

### Added

- Long-format member tables for the significant pathways of every analysis (`GO_<ont>_members.tsv`, `kegg_members.tsv`, GSEA `members.tsv`), keyed by `gene_id`. Columns: symbol, Entrez ID, log2FoldChange, padj and the pathway's p.adjust; GSEA tables add NES and leading-edge membership.
- GO redundancy reduction: `simplify(cutoff = 0.7)` tables and dotplots (`*_simplified`) beside the unchanged full results.
- GSEA running-score plots (`enrichplot::gseaplot2`) for up to 5 significant gene sets per NES sign, titled with NES and adjusted P.
- Every enrichment figure is also written as a PNG at the profile's raster dpi (`_300dpi.png`) beside the vector PDF. When nothing is significant, no figure is drawn and the summary says why.
- L2 (Phase 1): DEG tables gain `symbol`, `description`, group means of normalized counts and `log2FoldChange_shrunk` (apeglm, ashr fallback, never `normal`). Also new: `normalized_counts.tsv`, per-contrast `contrast_metadata.json`, an MA plot, symbol-labelled volcano and heatmap, and 300 dpi PNGs. The DEG definition and every DESeq2 statistic are unchanged.

- `rnaseq validate` and the run log (`logs/rnaseq.log`) warn when `project.yaml` sets the legacy `pvalue_cutoff`/`qvalue_cutoff`/`padj_cutoff` enrichment fields, which no longer select terms (warning, not error). For an unpinned KEGG analysis they name the snapshot that will be used and mark it "not pinned"; a pinned snapshot that does not exist is a validation error. Pin `annotation.enrichment.kegg.snapshot` for published results.
- `report.html` shows the MA plot, the simplified GO dotplots (GO ORA of the significant DE genes; GO GSEA) and up to two running-score plots per NES sign for GO and KEGG GSEA. A figure that does not exist is shown as "unavailable" with its reason.

### Fixed

- GO GSEA no longer emits "object not found" warnings from its cleanup code (results unchanged).
- `rnaseq.go/gsea/kegg.execute_*` (Python entry points) now build the R backend configuration with the same function as the Nextflow path, including the contrast numerator/denominator. They no longer read summary fields the R backends stopped writing.

### Documentation and tooling

- `examples/gsea_two_group` no longer sets the retired cutoff fields and uses `resource_provider: kegg_rest_snapshot`.
- `tools/figure_review/compare_runs.py --expected-changes FILE` gates a release-to-release comparison on a declared, reasoned list of intended changes (changed figure kinds, new figures by ID pattern and kind, added formats); any undeclared difference or removed figure fails. `tools/figure_review/expectations/v1.4.1-to-v1.5.0.json` declares the v1.5.0 changes.

### Runtime

- The supported runtime is unchanged: Linux x86-64, including Windows WSL2. The osx-arm64 downstream lock is maintained as an unqualified candidate and is not a supported platform.
- The downstream lock adds `bioconductor-apeglm` and `r-ashr` (Phase 1); the first run provisions the new downstream environment. Enrichment additionally uses enrichplot, KEGGREST, GOSemSim, GO.db and openssl, which were already locked as clusterProfiler dependencies and are now probed by the runtime check.

## 1.4.1 — 2026-10-29

Patch release from the `v1.4.0` post-release user-workflow audit. DESeq2, tximport, nf-core/rnaseq 3.26.0, QC and enrichment calculations, Nextflow workflows, reference validation, the upstream-reuse compatibility contract, the provenance schema, and scientific result tables are unchanged.

### Compatibility

- `rnaseq run --reuse-upstream` is now refused for raw-count projects, which have no upstream stage; previously the option was silently ignored. Remove the option from such commands.

### Fixed

- **Wizard next steps.** `rnaseq new` printed the same "run L1 for QC" instructions for L1 and L2 projects, so following them in an L2 project failed validation. L1 projects now list validate → plan → run L1 → review QC → `rnaseq design`; L2 projects list `rnaseq design` → validate → run L2 and state that the L2 run also performs the L1 QC. The commands show the actual project path, the FASTQ steps state the sample-ID rule, and the review and validation report say that L2 needs `rnaseq design` before it can run.
- **Validation blockers.** When validation failed for metadata, FASTQ, or contrast reasons, `rnaseq validate` also reported `reference genome not configured` for a correctly configured iGenomes (or custom) reference. It now reports `Project validation failed; fix the errors listed above before execution.` and names a reference problem only when the reference is actually not ready.
- **Raw-count upstream reuse.** `rnaseq design` no longer suggests `--reuse-upstream` for raw-count projects, and `rnaseq run` rejects the option for them (see Compatibility).
- **Early reuse preflight.** `--reuse-upstream` is checked before the confirmation prompt and before any run directory, input freezing, or runtime provisioning: the `CASE-ID/RUN-ID` format, the case and run IDs, the existence and SUCCESS state of the source run, and its upstream handoff manifest. Such mistakes no longer leave a FAILED run behind. The byte-level compatibility checks still run once the new run's inputs are frozen and are unchanged.
- **Reuse compatibility errors.** `frozen nfcore.local.config differs` now explains that the file records the local CPU/memory ceiling (changed `execution.max_cpus` / `execution.max_memory_gb`, or `auto` resolving to a different capacity on another host) and how to proceed; format and source errors include an example and point to `rnaseq status PROJECT --all`.
- **`rnaseq status --case-id`** is accepted as an alias of `--case`, matching `rnaseq run` and `rnaseq clean`.

### Documentation

- README (English and Traditional Chinese) and Quick Start describe both supported workflows, staged L1 → L2 and direct L2, and that `rnaseq design` writes the contrasts and refreshes the planning artifacts without needing a previous run. The obsolete statement that the L1 → L2 journey requires the development version was removed.
- FASTQ naming and sample-ID rules, with the example `Ctrl1_S1_L001_R1_001.fastq.gz` → sample `Ctrl1_S1`.
- Case ID versus run ID, and the exact `--reuse-upstream CASE-ID/RUN-ID` format.
- Foreground execution, monitoring from a second terminal, the effect of closing the terminal, shutdown, sleep, and `wsl --shutdown`, and optional `tmux`/`nohup` use.
- `rnaseq clean --dry-run` before cleanup.
- Upstream reuse requires identical resource limits; `docs/runtime.md` distinguishes `auto` limits from the 8 CPU / 12 GiB floor and fallback, and its supported-runtime section was updated from `v1.3.0`.

### Tests

- Tests that run the real L1/L2 backends are gated on exactly the R components the runtime requires, so an incomplete R skips them (or fails them under `NF_RNA_REQUIRE_R=1`) instead of failing part-way.
- Regression tests for the wizard instructions, validation blockers, upstream-reuse preflight and messages, the status alias, and a check that every `rnaseq` command and option shown in the user documentation exists.

## 1.4.0 — 2026-10-19

Changes since `v1.3.0`. Unless stated otherwise, DESeq2, tximport, nf-core/rnaseq 3.26.0, contrast semantics, enrichment methods, and scientific result tables are unchanged.

### Template-first project lifecycle

- `rnaseq new` now writes templates only: a `metadata.csv` header derived from the design choices (`sample_id,condition`, plus `batch` for a known technical batch and the pairing column for paired designs) and a header-only `contrasts.csv`. The provisional model is derived from the same choices. The wizard's analysis scope defaults to L1.
- New `rnaseq design PROJECT` reads the completed metadata, validates the design, and proposes contrasts for confirmation: later-vs-earlier for two groups (with `swap`), and every pairwise contrast in first-appearance order for multi-group designs (`--condition-order` overrides the order). On confirmation it writes `contrasts.csv`, freezes `design.condition_order`, sets the scope to L2, and regenerates the planning artifacts in the same transaction, so `rnaseq run` does not report a stale plan. Existing non-empty contrasts need `--force`; QC-only and invalid projects are refused.
- L1 validates and runs with an empty `contrasts.csv`; contrasts are required only for L2.
- New design validation: residual degrees of freedom from the expanded model matrix (a warning at validation, refused by `rnaseq design` and the L2 launch), clearer explanations of confounded batch, and warnings for partially confounded or single-sample batches, unused batch columns, and wrong-case column names. `rnaseq validate` prints batch-by-condition counts.
- Design flags are never silently ignored: `--covariate`, `--continuous-covariate`, and `--condition-column` without `--metadata`, `--pair-id` without a paired design, and design flags with the QC scope are errors.
- Reference flags are never silently ignored either. Non-interactive `rnaseq new` still uses iGenomes for Salmon when no reference flag is given and never consults the registry of local references (only the wizard offers them). `--reference-root` or `--reference-manifest` without `--reference-source` now selects that local reference, where they were previously dropped in favour of iGenomes for Salmon; the custom flags select `custom` the same way. Flags that do not belong to the selected source are errors. The project review names the selected reference and says when the iGenomes default applied, and the analysis plan records the iGenomes genome.

### Functional enrichment by default for L2

- Configuring a project for L2 (`rnaseq new --preset L2`, or `rnaseq design` moving an L1 project to L2) now selects every supported method: `analysis.enrichment: [go, kegg, gsea]`, i.e. GO ORA, KEGG ORA, GO GSEA, and KEGG GSEA. This applies equally to batch-adjusted and multi-group designs. `--no-enrichment` records an explicit opt-out (`[]`); a project already at L2 keeps its recorded selection. Statistical methods, thresholds, and the GSEA ranking are unchanged.
- The required annotation contract is collected at creation: the wizard asks for the count-matrix gene ID type (ENSEMBL, ENTREZID, SYMBOL) and shows the OrgDb; non-interactive L2 creation and `rnaseq design --yes` require `--gene-id-type` unless enrichment is opted out. Organisms without a supported OrgDb are refused with the opt-out as the remedy; enrichment is never dropped silently.
- `rnaseq run` and `rnaseq retry` check that the online KEGG REST resource is reachable before any run exists when KEGG ORA or GSEA is selected.
- The technical report lists an explicit status for every method and contrast (`SUCCESS`, `NO_SIGNIFICANT_TERMS`, `NOT_APPLICABLE`, `BLOCKED`, `FAILED`, `NETWORK_UNAVAILABLE`, `NOT_REQUESTED`), also delivered as `tables/report/enrichment_status.tsv` and summarized in the delivery README. Methods text describes GSEA only when it was selected. Design warnings from validation (e.g. partial batch confounding) are frozen into the downstream contract and stated in the report and README.

### GO ORA gene-member tables

- Fixed `GO_{BP,MF,CC}_gene_members.tsv` padding: an unmapped tested gene (NA Entrez ID) added one all-NA row for every term-gene pair, inflating member tables by orders of magnitude (11 GB per table and 43 GB of GO output for a 14-sample mouse project) and GO ORA runtime/memory. Members are now exactly the mapped source genes of each significant term's `geneID`. GO term tables, statistics, mapping tables, and figures are unchanged.

### L1-to-L2 upstream reuse

- `rnaseq run --reuse-upstream CASE-ID/RUN-ID` now works across the L1 → `rnaseq design` → L2 change for FASTQ projects (Salmon/tximport and HISAT2 + featureCounts). The frozen input manifest is compared only on what determines the upstream result: input file identities, upstream configuration, reference, and pipeline name and version. Scope, design, metadata, and contrasts may differ. The frozen nf-core parameters and local config remain byte-compared, the frozen samplesheet is now byte-compared too, and an unreadable manifest is never treated as compatible.

### L1 QC and figures

- L1 adds per-sample distribution boxplots of log2(raw count + 1) (`boxplot_log2_raw_counts`) and of blind VST values (`boxplot_vst`), with Tukey box statistics in `sample_distribution_qc.tsv`. Boxes follow the L1/PCA sample order and are filled by condition; a complete `batch` column is shown as an annotation row, never as a correction. The figures are included in the figure manifest, the delivery package, and the L1 and combined HTML reports. The PCA panel no longer collapses for two-sample projects: when PC2 has no variance (constant at 0), the y axis gets a fixed extent of ±20% of the PC1 span and a single break at 0 instead of a thin strip with overlapping tick labels, and the figure manifest records this. Scores, explained variance and equal axis scaling are unchanged.
- All R figures go through a central figure policy (fonts, typography, Okabe–Ito group colours, appropriate sequential/diverging scales) and are exported as PNG, 300 dpi TIFF, and PDF; delivery adds `figures/pdf` and per-module `figure_manifest.json`. Row z-scores keep undefined values as NA.
- PCA sample labels are placed with ggrepel so tight clusters stay legible; ggrepel is now part of the downstream runtime package probe. DEG heatmap titles are compact and no longer run past the page edge.
- GO and KEGG ORA dotplots carry a subtitle naming the comparison and the foreground gene set (significant, or higher in the numerator/denominator), and the report shows the correct KEGG ORA term counts.

### Figure standard

- New project option `figures: {preset, profile}` selects the figure style. `preset` is `standard` (canonical), `bioinfo` (a backward-compatible alias, resolved to `standard` everywhere) or `legacy`; `profile` is `general` (default), `nature`, `plos` or `elsevier`. The selection is validated, frozen with the run, and recorded in the figure manifests.
- **New projects** created by `rnaseq new` now write `figures: {preset: standard, profile: general}`, so they render with the standard preset. **Projects without a `figures` section keep the `legacy` preset on the `general` profile** and render exactly as before; the legacy preset is unchanged and every legacy figure, source-data table and manifest key reproduces the 1.3.0 output. Existing projects change style only if `figures:` is added.
- The `standard` preset covers every figure type: volcano, PCA (including the batch PCA), sample boxplots, sample correlation heatmap, DEG heatmaps, and the GO, KEGG and GSEA dotplots. It draws with the bundled Arimo font (SIL OFL 1.1, verified by hash before the first figure; never substituted by Arial, Helvetica or DejaVu Sans and never downloaded) and sizes figures to their content within the selected profile's limits (columns, maximum size, minimum text size, raster dpi). PLOS, which defines no single column, uses an 85 mm compact width for volcano and PCA.
- Approved colours of the standard preset: group colours 1-8 (more than eight levels keep slots 1-8 and add deterministic, colour-vision-aware extension colours), `nf-rna Standard Batch Palette v1` for batch annotations and batch PCA, volcano up/down colours, a diverging heatmap scale, and `nf-rna Standard Enrichment Scale v1` for adjusted P in GSEA dotplots (the ORA mark colour is its first stop). Scales, limits, transformations and all scientific mappings are unchanged. Colour-vision limitations are documented in `docs/figure-standard.md`.
- Figure manifests (`nf-rna.figure-manifest.v1`, schema unchanged) describe the applied standard preset: the font actually used and its SHA-256, the profile, and the palette names; additive keys only. Legacy manifests keep their original keys.
- `docs/figure-standard.md` is the authoritative description of presets, profiles, typography, the output contract and verification. Review tools for figure comparison, palette and profile audits and package checks are shipped under `tools/figure_review/`, and the sdist now includes `tools/` and `.claude/skills`.

### Fixed: stale text metrics across graphics devices

- R can reuse the text metric it last computed for a font when the next graphics device has a different resolution (PDF 72 dpi, PNG 150 dpi, TIFF 300/500 dpi). The first measurement on the new device then occasionally returned the previous device's units, mis-sizing axis and legend rows by a few points; in a controlled stress test about 3% of outputs were affected in the legacy drawing order. Every output device is now primed with a throwaway measurement before drawing, which draws nothing and changes no correct output. In whole-pipeline stress runs 0 of 6,120 figure files deviated with the fix and 3 of 6,120 without it. The cause was reproduced in plain R and grid; details and evidence are in `docs/figure-standard.md`.

### Execution and operations

- `rnaseq status` is an execution dashboard: route-specific pipeline stages, per-sample progress, running tasks with elapsed time, and the configured resource ceiling. `--watch` follows one run until it finishes, and `--all` lists one line per run. It remains read-only.
- New `rnaseq clean PROJECT --case-id CASE --run RUN` removes one finished run's disposable Nextflow work directory, located only from the run's recorded provenance. SUCCESS runs are eligible; FAILED and INTERRUPTED runs need `--force`; running runs are always refused. `--dry-run` shows the target. Cleanup is never automatic.
- `RNASEQ_WORK_ROOT` (an absolute path) relocates only the Nextflow task work directory to `<RNASEQ_WORK_ROOT>/<case>/<run>/work`; the launch directory and shared caches stay in place. `rnaseq doctor` reports the work root and its free space separately.

### Runtime

- Retired the remaining Docker runtime support: `rnaseq new` no longer writes `runtime.execution_image` (existing values still load and are ignored), macOS is reported as unsupported, retry of runs frozen for the Docker backend is refused, and `Dockerfile`, `environment.docker.yml`, and `.dockerignore` are removed from the tree (they remain at the `v1.2.1` tag). Linux Nextflow commands are unchanged.

### Tests

- Status timing tests are timezone-independent, and the publication figure delivery contract is enforced by tests.

## 1.3.0 — 2026-09-26

Runtime release. The supported runtime is now **Linux x86-64 or Windows WSL2 with Nextflow and Conda**; Docker is no longer required or used. nf-core/rnaseq 3.26.0, DESeq2, tximport, contrast semantics, and the enrichment methods are unchanged apart from the ORA mapping-QC correction and the Salmon/tximport sample-identity fix below.

Scientific qualification: the final real-data qualification (4T1, GRCm39, Salmon/tximport) was run at revision `9753e5d89ea56e0afe936733f6903368423348be`. The later runtime and usability changes listed under *Runtime and usability* were not part of that qualification run; they do not change DESeq2, tximport, sample identity, enrichment, references, or delivered scientific results, and are covered by the regression suite.

- FASTQ upstream runs in Conda: nf-core/rnaseq 3.26.0 (revision `e7ca46272c8f9d5ceee3f71759f4ba551d3217a4`) with `-profile conda`, and the HISAT2/featureCounts graph in an exactly pinned Conda environment with the same tool versions and builds as the former containers.
- Downstream L1, L2, enrichment, and report tasks run in a locked linux-64 Conda environment (`workflow/envs/locks/nf-rna-downstream-linux-64.lock.yml`, verified against `SHA256SUMS`) into which the non-editable nf-rna wheel is installed. One immutable prefix is created per (lock, wheel) pair, so an upgrade provisions a new prefix instead of colliding with an earlier run's environment.
- Every run records the exact nf-rna source commit. A source checkout reports its `HEAD`; an installed distribution reports the commit pip recorded for a `git+…@<ref>` install. Runs are refused when no commit is available (directory, sdist, or wheel installs) or the checkout has uncommitted changes. Provenance also records the lock and wheel SHA-256 and the bundled R-script inventory.
- Includes every correctness fix of 1.2.1 (see below): categorical label identity (`1` ≠ `01`), explicit schema 1.3 formula-variable types, frozen variable types passed to L1/L2, inferred-continuous warnings, the typed experimental-design report table, blank-line raw-count delivery validation, and the ORA significance-rule disclosure. The Docker image-identity changes of 1.2.1 do not apply to the Conda runtime.
- Fixed GO and KEGG ORA treating `annotation.mapping_warning_rate` as a blocking threshold. All four enrichment backends (GO/KEGG ORA and GO/KEGG preranked GSEA) now share one mapping QC: at or above `mapping_warning_rate` is `PASS`; from `minimum_mapping_rate` up to the warning threshold is `WARNING` and the backend runs; below `minimum_mapping_rate` is `BLOCKED`. **ORA results can change for projects whose tested-gene mapping rate lies between the two thresholds (0.50–0.70 by default)**, which 1.2.x blocked. ORA summaries now record the QC status and both thresholds.
- Raw-count delivery keeps the valid source column order of `input/counts.csv` instead of requiring it to match metadata order; analysis still uses the frozen metadata order.
- `rnaseq doctor`, `rnaseq run`, and `rnaseq retry` verify that Conda's effective channel configuration lists `conda-forge` before `bioconda`, as nf-core/rnaseq `-profile conda` requires; the Salmon route refuses to start otherwise. `rnaseq doctor` now ends with `Overall: READY` or `Overall: NOT READY` and exits with code 1 when any check fails.
- Fixed the Salmon reference contract. `rnaseq reference prepare` had built the decoy-aware index from the registered Ensembl `cdna.all` FASTA, whose versioned transcript IDs never equal the GTF `transcript_id`, so nf-core/rnaseq `CUSTOM_TX2GENE` failed, and which omits non-coding transcripts such as lncRNA. The builder now derives the transcriptome as nf-core does (nf-core GTF filter plus `rsem-prepare-reference --gtf`, bioconda RSEM 1.3.3), requires a PASS transcript-ID contract, and publishes the index under `salmon/gtf_derived/`. Every built Salmon declaration must now reference a PASS contract artifact, so validation, planning, `rnaseq doctor`, and run preflight reject an incompatible Salmon reference before nf-core runs. `--rebuild-salmon` replaces such a declaration while archiving the previous manifest and leaving the previous index in place. **Salmon results change for references rebuilt this way**, because the index now covers every GTF transcript on the genome.
- `rnaseq run` and `rnaseq retry` start Nextflow with plain logging (`NXF_ANSI_LOG=false`), no stdin, and in its own session without a controlling terminal. A run launched as a background job of an interactive shell (`nohup rnaseq run ... &`) was previously stopped by SIGTTOU when Nextflow's ANSI log ran `stty` on the terminal. Interrupting `rnaseq` (Ctrl-C, SIGTERM, or SIGHUP when not ignored) is forwarded to Nextflow as SIGTERM for a clean shutdown.
- `rnaseq run` no longer copies raw FASTQs into the run or hashes them repeatedly before launching Nextflow. Previously a run hashed every FASTQ twice (plan-freshness check and frozen input manifest) and copied it into `frozen/input/fastq`, about four times the input size in I/O before Nextflow started, which on a WSL `/mnt` drive meant well over 100 GB for a 37 GB project. Planning, freshness checks, validation and `rnaseq doctor` now identify FASTQs by path, size and mtime without reading them; the frozen input manifest records each file's resolved path, size, mtime and SHA256, computed at most once per file and reused while path, size and mtime are unchanged. FASTQs are referenced in place by the frozen samplesheet. Re-run `rnaseq plan` once after upgrading; planning manifests from earlier versions no longer match.
- Fixed Salmon/tximport sample identity. L1 and L2 imported `quant.sf` files by position, pairing the sample-sorted `quant_sf` mapping written by the control plane with the metadata sample order, so when metadata was not in sorted sample order each sample received another sample's quantification. Both now resolve `quant.sf` by sample name, require an exact one-to-one `sample_id` → path mapping, and check that tximport's columns equal the metadata samples. Earlier releases' Salmon/tximport route used the same positional import. **Salmon/tximport results change for projects whose metadata sample order differs from sorted order.** Raw-count and HISAT2/featureCounts routes are unaffected.
- macOS and native Windows are not supported by 1.3.0. v1.2.1 remains the last Docker-based release.

### Runtime and usability

- `execution.max_cpus` / `execution.max_memory_gb` accept `auto`, now the default for new projects and for projects without an `execution` section: the local Nextflow ceiling is sized from the running machine (CPU affinity and cgroup limits respected) minus an OS reserve, never below the 8 CPU / 12 GiB contract floor when the machine has it, and falling back to 8 CPUs / 12 GiB if capacity cannot be detected. Explicit integers are used unchanged. The decision is recorded as `resource_policy` in provenance and the frozen execution manifest.
- With a prebuilt Salmon index, nf-core/rnaseq `SALMON_QUANT`'s memory *request* is sized from the index (frozen as `frozen/nfcore.tuning.config`) so samples can quantify concurrently; its CPUs and all Salmon arguments are unchanged. The HISAT2 route now writes a Nextflow trace to `provenance/upstream.trace.txt`.
- `rnaseq status PROJECT` shows the latest run's state, phase, per-stage task progress (completed/cached/failed, running task names), delivery state and effective resources; `--watch` refreshes until the run finishes; `--case`/`--run` select a run; `--all` keeps the previous one-block-per-run listing. It is read-only, and reports a RUNNING run whose `rnaseq` process has gone as INTERRUPTED rather than guessing success.
- Each run writes its own execution log, `logs/rnaseq.log`, and process record; two launches with the same case ID no longer share a log. Ctrl-C/SIGTERM/SIGHUP during a run now records `INTERRUPTED` (previously the run stayed `RUNNING`), and unexpected errors record `FAILED` with a traceback in the run log. The delivery package still excludes logs.
- `rnaseq retry` accepts INTERRUPTED source runs as well as FAILED ones, with the same frozen-contract, identity and provenance checks. SUCCESS runs and runs still recorded RUNNING/CREATED (including a stale RUNNING run whose process has gone) remain non-retryable; the source run's INTERRUPTED state is kept as is.

## 1.2.1 — 2026-09-24

Correctness release for 1.2.0. DESeq2, tximport, nf-core/rnaseq 3.26.0, contrast semantics, and enrichment methods are unchanged.

- Fixed R metadata parsing that merged distinct categorical labels such as `1` and `01` (or `1`/`1.0`, `1e2`/`100`, and very long numeric IDs) into one factor level before DESeq2. Metadata is now read as text and the frozen variable types are applied afterwards, preserving each label; reference-level order is unchanged. **Results can change for projects whose categorical covariates or `pair_id` contained such labels**, because 1.2.0 fitted a model with merged levels. Other projects produce identical DESeq2 results.
- Schema 1.3 projects must declare a type in `design.variables` for every formula variable; incomplete projects are rejected before analysis.
- L1/L2 now receive design-variable types from the frozen, validated contract instead of re-reading `project.yaml`.
- Legacy schema 1.0–1.2 projects keep their inference, but each formula variable inferred as continuous now produces a validation warning.
- The report's Experimental design section lists the resolved variable types from the frozen contract and marks undeclared ones as inferred.
- For legacy schema 1.0–1.2 projects, `pca_scores.tsv` now includes an annotation column for each resolved design variable; PC coordinates and explained variance are unchanged.
- Fixed post-analysis delivery validation rejecting raw-count inputs that contain physically blank CSV lines (such as a trailing blank line), which input validation and R already ignore. The delivered `raw_counts.csv` remains byte-identical to the input and its manifest counts only real gene rows. This is a delivery-validation fix; DESeq2 inputs and results are unchanged.
- The GO/KEGG ORA significance rule (raw p-value ≤ `pvalue_cutoff` and q-value ≤ `qvalue_cutoff`; `p.adjust` is reported but not used for selection) is now documented and stated in the report. The rule itself is unchanged.
- The first-party Docker image must report the same nf-rna version as the CLI, is resolved to its local RepoDigest, and runs by that frozen `repository@sha256:…` reference. Production acceptance also requires a clean release revision label matching a clean source checkout. After upgrading, set existing projects' `runtime.execution_image` to the 1.2.1 image and re-run `rnaseq plan`.

## 1.2.0 — 2026-09-21

- Added explicitly typed categorical and continuous metadata covariates for generalized additive fixed-effect DESeq2 designs.
- Added covariate-adjusted differential expression for designs such as `~ age + condition`, `~ batch + condition`, and `~ batch + age + condition`, while preserving paired designs and directional categorical contrasts.
- Added known-batch adjustment through the DESeq2 design without modifying raw counts, together with batch-aware PCA outputs.
- Added preflight full-rank model-matrix validation and schema 1.3 while retaining read compatibility with project schemas 1.0–1.2.
- Extended scientific provenance with fitted categorical levels, continuous-variable summaries, model-matrix rank and columns, and the actual DESeq2 fit method.
- Qualified the generalized-design milestone on both raw-count and Salmon/tximport Docker/Nextflow production routes, including preservation of tximport normalization offsets.

## 1.1.2 — 2026-09-21

- Replaced the custom GHCR registry preflight with standard GitHub and Docker publishing actions for stable release tags.
- `v1.1.0` and `v1.1.1` remain valid source releases, but neither has an official GHCR image because publication stopped before image build or push. Official GHCR distribution begins with `v1.1.2`.
- This release contains no scientific-analysis, R, or Nextflow changes.

## 1.1.1 — 2026-09-21

- Fixed the fail-closed GHCR registry check to normalize case-insensitive HTTP challenge field names before validating the exact trusted GHCR Bearer challenge.
- `v1.1.0` remains a valid source release, but no official GHCR image was published for it because publication stopped before image build or push. Official GHCR distribution is deferred to `v1.1.1`.
- This hotfix contains no scientific-analysis, R, or Nextflow changes.

## 1.1.0 — 2026-09-20

- Added production GO and KEGG over-representation analysis (ORA) through the supported Nextflow graph, including stable GO terminal-state artifacts and term accounting.
- Corrected ORA to use the finite-p-value, statistically tested-gene universe for each contrast and deduplicated enrichment selections.
- Added structured scientific provenance, complete R session records, and OCI source-revision propagation through downstream outputs.
- Added guarded, release-only GHCR publication for the first official image, `ghcr.io/2002brian/nf-rna:1.1.0`; the official image is qualified for `linux/amd64` only. `linux/arm64` is not yet independently qualified.

## 1.0.0 — 2026-09-16

- Finalized Nextflow-owned downstream execution with a local first-party Docker runtime, immutable runs/retry/upstream reuse, managed-reference provenance, and curated report/delivery artifacts.
- Released the HISAT2 + featureCounts and nf-core/rnaseq Salmon/tximport routes, L1 expression QC, L2 DESeq2 differential expression, and optional GO/KEGG preranked GSEA.
- Qualified the paired-end Mouse HISAT2 + featureCounts path through reuse, metadata-ordered staging, L1, L2, no-enrichment technical reporting, delivery, and SUCCESS.
- Added a delivery-wide SHA-256 integrity manifest for final curated delivery packages.

## 1.0.0rc1 — 2026-09-12

- Added an active downstream L2 gate that rejects contrasts with fewer than two biological samples in either group before DESeq2 is launched; L1 remains available.
- Made the package version the single source for wheel/sdist metadata, runtime provenance, and `rnaseq --version`.
- Added immutable failed-run retry, installed-artifact workflow smoke coverage, and the routine release CI gate.

## 0.5.1 — 2026-09-06

- Corrected the nf-core/rnaseq 3.26.0 Salmon handoff to require and checksum the augmented tx2gene mapping used by tximport.
- Added opt-in production acceptance for checksum-bound managed references, explicit biological pairing, and preflight model-matrix rank validation.
- Added immutable requested/observed runtime identity, source commit and workflow hashes to run provenance.
- Connected raw and processed FastQC outputs to HISAT2 MultiQC and corrected downstream source labels.

## 0.5.0 — 2026-09-05

- Added the explicitly configured HISAT2 + featureCounts FASTQ backend, including a source-aware count-matrix handoff to existing L1/L2 analysis.
- Added local-reference HISAT2 index preparation, explicit count/strand contracts, and backend-aware readiness/planning.
- Added the reviewed interactive and explicit non-interactive project-creation interface, including strict FASTQ samplesheet import and safe scaffolding.
- Kept Salmon/tximport and external raw-count routes compatible with their existing contracts.

## 0.4.3 — 2026-09-03

Initial public release of nf-rna.

- Versioned FASTQ and raw-count project contracts with read-only validation.
- Pinned nf-core/rnaseq 3.26.0 FASTQ route with Salmon and immutable case/run provenance.
- First-party L1 QC, L2 DESeq2, and preranked GO/KEGG GSEA execution paths.
- Deterministic planning, scoped delivery assembly, and runtime/resource preflight checks.
- Public packaging controls that exclude biological inputs, references, indices, runs, work directories, and machine-local metadata by default.
