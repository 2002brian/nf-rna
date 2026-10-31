# nf-rna v1.7.0 — Local Research Workspace & Publication Figure Export

> **No analysis changes.** v1.7.0 changes only the optional browser GUI.
> DESeq2, tximport, LFC shrinkage, GO/KEGG ORA, GSEA, Salmon, HISAT2 +
> featureCounts, scientific thresholds, references and annotation, the R
> analysis and plotting scripts, the Nextflow workflows and the Conda locks are
> identical to v1.6.1. v1.6.1 results do not need to be recomputed.

Full list of changes: [CHANGELOG](https://github.com/2002brian/nf-rna/blob/v1.7.0/CHANGELOG.md#170--2026-10-31).
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

- Results, tables, figure files, manifests, reports and delivery packages are
  unchanged. New runs record `nf_rna_version` 1.7.0.
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
- The scientific limitations of v1.6.1 are unchanged (CHANGELOG 1.6.1, Known
  limitations).
