English | [繁體中文](README_zh-TW.md)

# Upstream execution boundary

The local profile invokes the external, pinned `nf-core/rnaseq 3.26.0` pipeline with `-profile conda` and then this first-party DSL2 downstream workflow:

```text
rnaseq CLI → nf-core/rnaseq → standardized upstream outputs → downstream Nextflow workflow → R modules
```

The Python control plane freezes a version-pinned input contract and records the stable nf-core handoff boundary. Its frozen `analysis_level` is the graph selector: an `L1` project runs L1 and the L1 technical report only; an `L2` project runs L1, L2, then any independently selected GO ORA, KEGG ORA, GO preranked GSEA (BP/MF/CC), and KEGG preranked GSEA before the technical HTML report. The graph never infers L2 from contrasts, metadata, or available workflow modules. Server executor settings remain intentionally deferred.

Each downstream computational process explicitly declares `conda params.downstream_runtime_prefix`. After the user confirms a run, the control plane creates or verifies a deterministic cached Conda prefix from the platform lock, installs a non-editable nf-rna wheel with `pip --no-deps`, and freezes the lock, wheel, and bundled R-script identities into the run contract. Nextflow activates that prefix for every first-party task; R scripts are located with `importlib.resources` from the installed wheel, never `/opt/nf-rna/r` or the caller's editable environment. The FASTQ upstream also runs in Conda: nf-core/rnaseq with `-profile conda`, and the HISAT2/featureCounts graph in its exactly pinned Conda environment.

Independent `ENRICHMENT_ANALYSIS` tasks may run concurrently when their summed
CPU and memory requests fit the effective aggregate local budget. Nextflow's
local executor provides this scheduling guard without a fixed `maxForks` cap.

Each enrichment task publishes its own backend directory below the immutable
L2 output: `downstream/l2/enrichment/go/`, `kegg/`, `gsea_go/`, and
`gsea_kegg/`. The process never publishes their shared
`enrichment/` directory as a task output. Instead, each task publishes its
module-specific `enrichment/*` output relative to the L2 publish root, so one
successful backend cannot mask the other when `overwrite: false` is enforced.

When GSEA is enabled, the final report receives the collected output channels
from both GSEA backends. It therefore starts only after GO GSEA and KEGG GSEA
have both completed; a missing or failed enabled backend prevents report and
delivery finalization.

## Local resource policy

The checked-in local configuration uses explicit, non-scientific runtime
classes: SMALL (1 CPU, 2 GiB, 2 h), MEDIUM (4 CPUs, 8 GiB, 8 h), and LARGE
(8 CPUs, 12 GiB, 12 h). L1 and the technical report use SMALL; L2 and each
GSEA backend use MEDIUM. Sample-level upstream and independent enrichment tasks
are not artificially serialized; Nextflow admits ready tasks while their summed
requests fit the frozen effective aggregate ceiling. nf-core process-specific
label requests remain intact.

These limits are intended to avoid local memory oversubscription and do
not alter any inputs, model, filtering, thresholds, ranking, enrichment
calculation, or published artifact. Future server profiles may override the
resource directives explicitly. Use `rnaseq doctor [PROJECT]` before a local
run to inspect host capacity and the Conda runtime.
