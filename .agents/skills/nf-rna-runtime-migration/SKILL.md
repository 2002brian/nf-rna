---
name: nf-rna-runtime-migration
description: Guide Docker-to-native Nextflow and Conda migration work in nf-rna without changing scientific semantics.
---

# Purpose

Guide safe migration of nf-rna from Docker-based execution to reproducible native Nextflow + Conda execution without changing scientific semantics.

# Target architecture

```text
rnaseq CLI
-> Nextflow
-> pinned nf-core/rnaseq
-> Nextflow-managed upstream Conda environments
-> nf-rna downstream workflow
-> immutable pinned downstream Conda runtime
```

# Non-negotiable scientific invariants

- Do not change DESeq2 statistical semantics during runtime migration.
- Preserve explicit biological contrasts and contrast direction.
- Preserve paired/covariate-aware design behavior.
- Preserve tximport offsets.
- Preserve GO/KEGG ORA and GSEA semantics.
- Preserve raw-count staging/alignment behavior.
- Preserve provenance sufficient to reproduce scientific results.
- Do not introduce automatic pairwise contrasts.
- Keep runtime migration separate from scientific feature development.

# Runtime invariants

- Never replace Docker with uncontrolled host Python/R dependencies.
- Nextflow remains the workflow/process orchestrator.
- Runtime dependencies must be pinned.
- The exact nf-rna code executed downstream must be identifiable.
- Do not assume the caller's activated Conda environment is the production runtime.
- Do not claim platform support without real execution evidence.

# Supported-platform policy

Candidate targets are `linux-64` and `osx-arm64`. A platform becomes supported only after a real native smoke execution succeeds.

# Migration order

1. Establish the immutable downstream Conda runtime contract.
2. Resolve and lock `linux-64`.
3. Resolve and lock `osx-arm64`.
4. Qualify the downstream runtime independently.
5. Convert first-party downstream Nextflow execution.
6. Convert nf-core/rnaseq from Docker to Conda.
7. Migrate `rnaseq doctor`.
8. Migrate runtime provenance.
9. Run real end-to-end native smoke tests.
10. Remove Docker/GHCR assets only after native qualification.

Do not skip gates merely because unit tests pass.

# Milestone behavior

Before editing:

- Inspect the current implementation.
- Identify affected runtime assumptions.
- Identify scientific code that must remain untouched.

After editing:

- Run targeted tests and the full pytest suite.
- Run real runtime/integration tests when the milestone affects execution.
- Run `git diff --check` and inspect the diff for unintended scientific changes.
- Report unresolved blockers separately from optional improvements.

# Classification

Use `BLOCKER`, `REQUIRED`, `OPTIONAL`, and `OUT_OF_SCOPE`.

Do not spend implementation effort on `OPTIONAL` items while `BLOCKER` or `REQUIRED` migration work remains.

# Docker removal gate

Docker/GHCR artifacts may be removed only when:

- The downstream native runtime is reproducibly locked.
- Upstream nf-core native execution is qualified.
- Required provenance replacement exists.
- Linux native smoke passes.
- macOS ARM64 native smoke passes if macOS ARM64 is claimed supported.
- The full regression suite passes.

# Cache policy

Runtime migration must not silently delete resumable Nextflow state. Shared Conda caches must remain separate from scientific deliverables. Cache/work cleanup is a separate operational milestone unless required for runtime correctness.

# Release rule

Do not commit, push, tag, publish, or create a release unless explicitly requested.
