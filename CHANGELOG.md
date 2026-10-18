# Changelog

All notable changes to this project are documented here.

## 1.4.0 — Unreleased

The release date is set when the release is finalized. Changes since `v1.3.0`. Unless stated otherwise, DESeq2, tximport, nf-core/rnaseq 3.26.0, contrast semantics, enrichment methods, and scientific result tables are unchanged.

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
