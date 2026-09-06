# Changelog

All notable changes to this project are documented here.

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
