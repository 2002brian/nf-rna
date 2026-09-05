# Scientific contract

## Intent and boundaries

This repository produces technical bulk RNA-seq analysis artifacts. It does not provide diagnosis, treatment recommendations, clinical interpretation, or automated biological conclusions. Operators remain responsible for experimental design, sample identity, reference selection, and interpretation.

## Input integrity

- Count matrices must have `gene_id` first, unique nonblank sample headers, at least one gene, and non-negative mathematically integral values.
- Metadata is preserved as supplied. It must contain nonblank unique `sample_id` values and all variables in the additive design formula.
- The count/metadata sample sets must agree exactly; their ordering may differ.
- Contrasts are case-sensitive, directional, and must reference a factor and existing levels present in the design.
- FASTQ preprocessing is user-declared. The pipeline does not infer whether data are trimmed from names, paths, or contents.
- The Salmon route preserves tximport average-transcript-length offsets. The HISAT2 + featureCounts route uses only nonnegative integral gene counts with `DESeqDataSetFromMatrix`; the two sources are never combined.
- HISAT2 + featureCounts requires explicit `unstranded`, `forward`, or `reverse` library strandedness. It counts exons grouped by `gene_id`, excludes multimappers and ambiguous overlaps, and counts from a separate BAM that excludes secondary/supplementary records (`samtools -F 0x900`) while retaining the original tagged diagnostic BAM; paired-end data count fragments with both mates mapped while excluding chimeras.

## L1

L1 performs technical expression QC. It removes all-zero genes, then genes whose total count is below the configured L1 rule; records the retained gene universe; normalizes; produces blind VST, PCA, sample correlation, and QC tables/figures. It never automatically removes samples.

L1 does not perform differential-expression inference. An L1 preset never launches the control-plane L2 or enrichment graph.

## L2 differential expression

L2 fits the validated additive design using DESeq2 and extracts each configured contrast in `numerator / denominator` direction. The reported sign convention is `log2FC > 0` higher in the numerator and `log2FC < 0` higher in the denominator. The exact L1 gene universe is used. DESeq2 independent filtering remains enabled and native adjusted p-values and unshrunken log2 fold changes are retained.

Significance is `padj < thresholds.padj` and `abs(log2FoldChange) >= thresholds.abs_log2fc`. The pipeline warns about limited replication and blocks unsuitable unreplicated two-group DE inference; it does not alter DESeq2 behavior to force a result.

## Preranked GSEA

When `analysis.enrichment: gsea` is selected for L2, the public setting expands only to GO preranked GSEA (BP, MF, CC) and KEGG preranked GSEA. GO/KEGG ORA is not in the v1.0 production path.

The ranking uses all finite tested-gene DESeq2 Wald statistics, without a DEG, adjusted-p-value, p-value, or effect-size prefilter. Mapping/collapse and tie ordering are deterministic and audited. GO uses supported local organism annotation packages. KEGG uses its declared online provider; network failure is reported rather than silently replaced by a different resource.

## Reference and provenance

Reference identity is explicit: species, assembly/release, checksums, selected-backend index status, and index-build provenance are frozen for a run. HISAT2 preparation records annotation-derived splice sites and its genome-index inputs; it never reuses a Salmon index. No reference is downloaded, auto-discovered, or silently migrated during validation/planning. External runtime behavior and versions are captured in provenance without changing the scientific configuration.
