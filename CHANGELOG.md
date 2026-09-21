# Changelog

All notable changes to this project are documented here.

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
