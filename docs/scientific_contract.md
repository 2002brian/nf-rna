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
- New deliveries make source semantics explicit. Imported and featureCounts sources deliver nonnegative integral `delivery/counts/raw_counts.csv`; Salmon/tximport delivers possibly non-integral `delivery/counts/estimated_counts.csv`. `delivery/counts/vst.csv` is normalized/transformed expression for visualization, never a model-fitting input. Its artifact manifest records source type, construction method, checksum, dimensions and ordered samples.
- New nf-core/rnaseq 3.26.0 runs require `salmon.merged.tx2gene_augmented.tsv`, including its self-mappings for quantified transcripts absent from the ordinary GTF mapping. Its path, role and SHA-256 are frozen. Historical ordinary mappings remain readable only with an explicit `historical_ordinary` label.
- HISAT2 + featureCounts requires explicit `unstranded`, `forward`, or `reverse` library strandedness. It counts exons grouped by `gene_id`, excludes multimappers and ambiguous overlaps, and counts from a separate BAM that excludes secondary/supplementary records (`samtools -F 0x900`) while retaining the original tagged diagnostic BAM; paired-end data count fragments with both mates mapped while excluding chimeras.

## L1

L1 performs technical expression QC. It removes all-zero genes, then genes whose total count is below the configured L1 rule; records the retained gene universe; normalizes; produces blind VST, PCA, sample correlation, and QC tables/figures. It never automatically removes samples.

L1 does not perform differential-expression inference. An L1 preset never launches the control-plane L2 or enrichment graph.

## L2 differential expression

L2 fits the validated additive design using DESeq2 and extracts each configured contrast in `numerator / denominator` direction. The reported sign convention is `log2FC > 0` higher in the numerator and `log2FC < 0` higher in the denominator. The exact L1 gene universe is used. DESeq2 independent filtering remains enabled and native adjusted p-values and unshrunken log2 fold changes are retained.

Biological pairing is declared by `design.pairing_column`; it is never inferred from formula order and is independent of paired-end read layout. Each block must contain exactly one observation from each requested contrast level. The additive model matrix is constructed during validation and rank-deficient designs are rejected before execution.

Significance is `padj < thresholds.padj` and `abs(log2FoldChange) >= thresholds.abs_log2fc`. The pipeline warns about limited replication and blocks unsuitable unreplicated two-group DE inference; it does not alter DESeq2 behavior to force a result.

## Preranked GSEA

When `analysis.enrichment: gsea` is selected for L2, the public setting expands only to GO preranked GSEA (BP, MF, CC) and KEGG preranked GSEA. GO/KEGG ORA is not in the v1.0 production path.

The ranking uses all finite tested-gene DESeq2 Wald statistics, without a DEG, adjusted-p-value, p-value, or effect-size prefilter. Mapping/collapse and tie ordering are deterministic and audited. GO uses supported local organism annotation packages. KEGG uses its declared online provider; network failure is reported rather than silently replaced by a different resource.

## Reference and provenance

Reference identity is explicit: species, assembly/release, checksums, selected-backend index status, and index-build provenance are frozen for a run. HISAT2 preparation records annotation-derived splice sites and its genome-index inputs; it never reuses a Salmon index. No reference is downloaded, auto-discovered, or silently migrated during validation/planning. External runtime behavior and versions are captured in provenance without changing the scientific configuration.

Managed Salmon and HISAT2 reference preparation is host-native only. It validates the controlled builder executable/version set before staging, records non-scientific host-builder facts separately from biological/index identity, and atomically publishes only validated artifacts. Existing valid Docker-built manifests remain readable and are not relabelled or rebuilt merely because their historical builder mode differs.

Schema 1.2 adds an explicit HISAT2 strategy discriminator. Legacy
graph-embedded splice-site manifests remain readable. A
`genome_only_runtime_splices` reference instead requires a non-empty,
checksum-bound GTF-derived splice-site file and passes it at alignment runtime;
the index-builder and runtime-aligner versions are distinct provenance fields.
An unvalidated version combination is not execution-ready. Its tiny Human
FASTQ smoke fixture is technical software validation only, never biological
evidence.

Production acceptance is opt-in with `reference.acceptance: production`. It requires `reference.source: local`, manifest schema 1.1, `purpose: production`, verified local asset hashes, and a complete index bound to the same FASTA/GTF identity. Synthetic/fixture identities and legacy manifests are rejected. To migrate schema 1.0 safely, review the provider/release/assembly identity and every local checksum, change the schema to 1.1, and deliberately set `purpose`; never add `purpose: production` merely to satisfy validation. Optional source URLs/accessions are provenance only. An optional upstream checksum is reported as verified only when it is supplied and matches the local asset.

The first intended human production identity is `Homo sapiens`, Ensembl release 116, GRCh38 patch p14. This repository configures the identity contract only; it does not contain or download those assets.
