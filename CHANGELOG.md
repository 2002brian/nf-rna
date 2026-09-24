# Changelog

All notable changes to this project are documented here.

## 1.3.0 — 2026-09-24

Runtime release. The supported runtime is now **Linux x86-64 or Windows WSL2 with Nextflow and Conda**; Docker is no longer required or used. nf-core/rnaseq 3.26.0, DESeq2, tximport, contrast semantics, and the enrichment methods are unchanged apart from the ORA mapping-QC correction below.

- FASTQ upstream runs in Conda: nf-core/rnaseq 3.26.0 (revision `e7ca46272c8f9d5ceee3f71759f4ba551d3217a4`) with `-profile conda`, and the HISAT2/featureCounts graph in an exactly pinned Conda environment with the same tool versions and builds as the former containers.
- Downstream L1, L2, enrichment, and report tasks run in a locked linux-64 Conda environment (`workflow/envs/locks/nf-rna-downstream-linux-64.lock.yml`, verified against `SHA256SUMS`) into which the non-editable nf-rna wheel is installed. One immutable prefix is created per (lock, wheel) pair, so an upgrade provisions a new prefix instead of colliding with an earlier run's environment.
- Every run records the exact nf-rna source commit. A source checkout reports its `HEAD`; an installed distribution reports the commit pip recorded for a `git+…@<ref>` install. Runs are refused when no commit is available (directory, sdist, or wheel installs) or the checkout has uncommitted changes. Provenance also records the lock and wheel SHA-256 and the bundled R-script inventory.
- Includes every correctness fix of 1.2.1 (see below): categorical label identity (`1` ≠ `01`), explicit schema 1.3 formula-variable types, frozen variable types passed to L1/L2, inferred-continuous warnings, the typed experimental-design report table, blank-line raw-count delivery validation, and the ORA significance-rule disclosure. The Docker image-identity changes of 1.2.1 do not apply to the Conda runtime.
- Fixed GO and KEGG ORA treating `annotation.mapping_warning_rate` as a blocking threshold. All four enrichment backends (GO/KEGG ORA and GO/KEGG preranked GSEA) now share one mapping QC: at or above `mapping_warning_rate` is `PASS`; from `minimum_mapping_rate` up to the warning threshold is `WARNING` and the backend runs; below `minimum_mapping_rate` is `BLOCKED`. **ORA results can change for projects whose tested-gene mapping rate lies between the two thresholds (0.50–0.70 by default)**, which 1.2.x blocked. ORA summaries now record the QC status and both thresholds.
- Raw-count delivery keeps the valid source column order of `input/counts.csv` instead of requiring it to match metadata order; analysis still uses the frozen metadata order.
- `rnaseq doctor`, `rnaseq run`, and `rnaseq retry` verify that Conda's effective channel configuration lists `conda-forge` before `bioconda`, as nf-core/rnaseq `-profile conda` requires; the Salmon route refuses to start otherwise. `rnaseq doctor` now ends with `Overall: READY` or `Overall: NOT READY` and exits with code 1 when any check fails.
- macOS and native Windows are not supported by 1.3.0. v1.2.1 remains the last Docker-based release.

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
