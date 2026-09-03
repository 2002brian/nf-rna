# Upstream execution boundary

The local Docker profile invokes the external, pinned `nf-core/rnaseq 3.26.0` pipeline and then this first-party DSL2 downstream workflow:

```text
rnaseq CLI → nf-core/rnaseq → standardized upstream outputs → downstream Nextflow workflow → R modules
```

The Python control plane freezes a version-pinned input contract and records the stable nf-core handoff boundary. Its frozen `analysis_level` is the graph selector: an `L1` project runs L1 and the L1 technical report only; an `L2` project runs L1, L2, optional GO preranked GSEA (BP/MF/CC) and KEGG preranked GSEA, then the technical HTML report. The graph never infers L2 from contrasts, metadata, or available workflow modules. It does not invoke GO or KEGG ORA in production. Server executor settings remain intentionally deferred.

The default downstream workflow runs at most one `ENRICHMENT_ANALYSIS` task at a
time (`maxForks 1`) to avoid concurrent R enrichment jobs exhausting a local
Docker host. This does not alter module selection or calculation; a future
server-specific Nextflow configuration can override that process directive.

Each enrichment task publishes its own backend directory below the immutable
L2 output: `downstream/l2/enrichment/gsea_go/` and
`downstream/l2/enrichment/gsea_kegg/`. The process never publishes their shared
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
(6 CPUs, 12 GiB, 12 h). L1 and the technical report use SMALL; L2 and each
GSEA backend use MEDIUM. `ENRICHMENT_ANALYSIS` retains `maxForks 1`, and the
frozen upstream local configuration also limits nf-core `SALMON_QUANT` to one
MEDIUM task at a time. The global local resource ceiling is LARGE.

These limits are intended to avoid local Docker memory oversubscription and do
not alter any inputs, model, filtering, thresholds, ranking, enrichment
calculation, or published artifact. Future server profiles may override the
resource directives explicitly. Use `rnaseq doctor [PROJECT]` before a local
run to inspect host/Docker capacity and architecture; dynamic upstream
nf-core image architectures must be checked in the frozen Nextflow trace.
