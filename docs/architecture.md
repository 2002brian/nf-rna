# Architecture

## Scope

`nf-rna` separates configuration validation, deterministic planning, upstream quantification, downstream analysis, and client delivery. It is not a workflow engine replacement and it does not mutate source metadata or inputs during validation. Its stable CLI and Python namespace remain `rnaseq`.

```text
project.yaml + metadata.csv + contrasts.csv + input
                         |
                         v
              validate / deterministic plan
                         |
                         v
                 immutable case/run freeze
                         |
      +------------------+------------------+
      |                                     |
      v                                     v
FASTQ: nf-core/rnaseq 3.26.0          raw counts: staged matrix
      |                                     |
      +---------- frozen execution boundary-+
                         |
                         v
               first-party Nextflow workflow
                         |
              L1 -----------> L1 report
               |
               +-- L2 (only for preset L2) --> optional GSEA --> L2 report
                         |
                         v
                    curated delivery
```

## Project contracts

`project.yaml` is strict and versioned. It specifies an `L1` or `L2` preset, input route, design formula, thresholds, and optional GSEA selection. Unknown configuration fields are rejected. `metadata.csv` permits additional columns but requires `sample_id` and all formula variables. `contrasts.csv` has an exact four-column contract.

The project preset is frozen as `analysis_level` in the downstream contract. That immutable field selects the workflow graph:

- **L1:** L1 analysis and L1-only report; no control-plane L2 or enrichment process is constructed.
- **L2:** L1, L2, then only configured GSEA backends, followed by the regular report.

The graph does not infer L2 from the presence of contrasts, condition metadata, or an installed R module.

`paired_two_group` is a specialization of the existing additive-design path, not a separate analysis pipeline. Its explicit `pair_id` column, formula, contrast-specific pair membership, complete-pair counts, and analyzed-sample counts are validated before execution and frozen into the input manifest and downstream contract. L1, the single DESeq2 L2 fit, explicit contrast extraction, and optional preranked GSEA then follow the same graph as unpaired projects.

## FASTQ boundary

For `upstream.quantification.method: salmon`, the control plane invokes pinned nf-core/rnaseq 3.26.0 with Salmon pseudoalignment and a frozen params file. Its stable downstream handoff contains validated `quant.sf` and tx2gene material needed by tximport. For `hisat2_featurecounts`, a first-party Nextflow graph performs lane alignment, per-sample BAM merge/index, featureCounts and deterministic matrix assembly; it hands off a staged canonical integer matrix to `DESeqDataSetFromMatrix`. The downstream workflow receives a narrow task-staged bundle rather than arbitrary host paths from provenance.

FASTQ preprocessing is declared as `raw` or `pretrimmed`. The Salmon route maps `pretrimmed` to nf-core's native `skip_trimming: true`; the HISAT2 route bypasses fastp while retaining QC. Local references are checksum-bound by a manifest and require the index appropriate to the selected backend. HISAT2 does not consume a Salmon index.

## Immutable case/run lifecycle

Every authorized execution creates `runs/<case-id>/<run-id>/` with frozen configuration/contracts, logs, upstream outputs, staged downstream inputs, downstream artifacts, provenance, and a curated delivery tree. The run ID is an internal identity; client-facing delivery filenames use the date only.

Nextflow launch/cache/work state belongs in a local execution root outside the persistent immutable run. Reuse of upstream data is explicit and requires a matching frozen contract. Delivery finalization uses an allowlist and a scoped, symlink-safe AppleDouble cleanup pass restricted to the new delivery root. It then writes `delivery_manifest.yaml`: a deterministic relative-path SHA-256 and byte-size inventory of delivered files. The manifest is deliberately excluded from its own inventory.

## Reproducibility boundaries

Planning is deterministic for unchanged source inputs. Execution provenance captures selected versions, contracts, references, and runtime facts. Runtime timestamps, container scheduling, external KEGG availability, and external resource evolution are not deterministic scientific inputs and are recorded as such rather than hidden.
