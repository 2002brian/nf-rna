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
- **Strandedness check (HISAT2 + featureCounts, v1.6.1).** Each count-only BAM is counted again with featureCounts `-s 1` and `-s 2`, with all other options unchanged (`STRAND_EVIDENCE`). The check compares `r = assigned(-s 1) / assigned(-s 2)` with the declaration. It decides nothing with fewer than 1,000 strand-assignable fragments, calls a strand when one direction has at least 3× the other, and calls the reads unstranded when `r` lies within 1/1.5..1.5; anything else is `INCONCLUSIVE`. Per sample, the status is one of:
  - `CONSISTENT`;
  - `CONTRADICTS_DECLARATION` (warning: the opposite direction dominates);
  - `APPEARS_UNSTRANDED` (warning: a stranded declaration on balanced reads);
  - `APPEARS_STRANDED` (advisory: an unstranded declaration on stranded reads; unstranded counting stays valid);
  - `INCONCLUSIVE`.

  Opposite-strand gene overlaps move `r` towards 1 and cannot make the opposite direction dominate. The check reports only: it never changes the declared strandedness, never recounts, and never fails a run. Its result is in `upstream/hisat2_featurecounts/qc/strandedness/strandedness_qc.json`, the handoff manifest, the frozen contract (`upstream_qc.strandedness`), the run log, the report and the GUI. The Salmon route records the library type Salmon used and is `NOT_EVALUATED` by this check: with a declared library type Salmon records only compatible orientations, and under `auto` nf-core infers the type itself.

## L1

L1 performs technical expression QC. It removes all-zero genes, then genes whose total count is below the configured L1 rule; records the retained gene universe; normalizes; produces blind VST, PCA, sample correlation, and QC tables/figures. Per-sample distribution boxplots of log2(raw count + 1) and of blind VST values are drawn in the L1/PCA sample order, with their Tukey box statistics written to `sample_distribution_qc.tsv`; a missing or duplicated sample stops the plot. A `batch` column is shown only as an annotation; plotted values are never batch-corrected. L1 never automatically removes samples.

L1 does not perform differential-expression inference. An L1 preset never constructs the Nextflow L2 or enrichment graph.

## L2 differential expression

L2 fits the validated additive design using DESeq2 and extracts each configured contrast in `numerator / denominator` direction. The reported sign convention is `log2FC > 0` higher in the numerator and `log2FC < 0` higher in the denominator. The exact L1 gene universe is used. DESeq2 independent filtering remains enabled and native adjusted p-values and unshrunken log2 fold changes are retained.

Biological pairing is declared by `design.type: paired_two_group` and `design.pair_id`; it is never inferred from formula order, sample names, metadata order, or paired-end read layout. For each explicit contrast, every pair must contain exactly one numerator and one denominator observation. Incomplete or duplicated pair structures are rejected, not silently excluded. L2 requires at least two complete biological pairs. The declared additive formula is used directly by one DESeq2 object, and the model matrix is constructed during validation so rank-deficient or confounded designs (including a batch that is completely confounded with condition, or constant within every biological pair) are rejected before execution. A model without residual degrees of freedom, counted from the expanded model matrix, is reported as a validation warning (L1 QC stays available) and is refused where L2 is configured (`rnaseq design`) and launched. A `batch` term in the provisional L1 design does not batch-correct any plotted value: L1 creates no `removeBatchEffect` or corrected matrix, and formal inference is the L2 fit.

Significance is `padj < thresholds.padj` and `abs(log2FoldChange) >= thresholds.abs_log2fc`. The pipeline warns about limited replication and blocks unsuitable unreplicated two-group DE inference; it does not alter DESeq2 behavior to force a result.

### DEG definition (unchanged in v1.5.0)

The DEG definition is the v1.4.1 definition and is not changed by the v1.5.0 outputs below: `padj < thresholds.padj` AND `abs(log2FoldChange) >= thresholds.abs_log2fc`, where `log2FoldChange` is DESeq2's raw (unshrunken) maximum-likelihood estimate and `padj` is DESeq2's Benjamini–Hochberg adjusted Wald p-value with independent filtering, optimised at `alpha = thresholds.padj` (v1.6.1; the chosen baseMean filter threshold is recorded). `baseMean`, `log2FoldChange`, `lfcSE`, `stat`, `pvalue` and `padj` are produced by the same single fit and `results()` call as in v1.4.1. `log2FoldChange_shrunk` never enters the definition, the up/down split, the GO/KEGG ORA foregrounds, or the GSEA ranking (Wald `stat`).

**Effect-size filtering, not thresholded testing.** nf-rna tests H0: LFC = 0 (Wald) and then keeps genes whose raw |LFC| reaches the configured threshold. This is effect-size filtering: the reported p-values do not test |LFC| > threshold, and a gene passing both filters is not shown to have |LFC| above the threshold with the stated error rate. Thresholded hypothesis testing (`results(..., lfcThreshold = t, altHypothesis = "greaterAbs")`) tests H0: |LFC| ≤ t directly and yields different, generally larger p-values; nf-rna does not use it. LRT, interaction models and `lfcThreshold` testing are not part of this release.

### Raw and shrunken LFC

- `log2FoldChange` (raw MLE) defines significance and direction and is the value to report as the tested effect.
- `log2FoldChange_shrunk` is an empirical-Bayes shrunken estimate for display and ranking (volcano/MA x/y axis, sorting by effect size). It reduces the large, noisy LFCs of low-count or high-dispersion genes. It is computed per contrast from the same fitted object with `lfcShrink`, in this order:
  1. The contrast is a coefficient of `resultsNames(dds)` (`<factor>_<numerator>_vs_<denominator>`): `lfcShrink(type = "apeglm", coef = ...)`.
  2. Otherwise (both levels belong to the factor but the numerator/denominator pair is not a coefficient): on a copy of the fitted object, the factor is releveled so the denominator is the reference and `nbinomWaldTest` is rerun with the original size factors and dispersions; apeglm is applied to that copy's coefficient only if the copy's raw `log2FoldChange` and `pvalue` reproduce the original results (identical NA pattern, `all.equal` mean relative difference ≤ 1e-4; a reparametrized IRLS refit differs from the original by about 1e-6). Otherwise the copy is discarded.
  3. Otherwise: `lfcShrink(type = "ashr", contrast = c(factor, numerator, denominator))`.
  `type = "normal"` is never used. Method, coefficient or contrast, whether a relevel was used and its verification result are recorded in `contrast_metadata.json`.
- **R warnings (v1.6.1).** Warnings raised by the shared fit (`DESeq()`, or the gene-wise fallback), by `results()` and by shrinkage (including the relevel refit) are recorded in `contrast_metadata.json` `r_warnings`, by step, with the original message, the calling function and the number of occurrences. For example, apeglm's "the line search routine failed" from `nbinomGLM`. They are also recorded in `backend_summary.json` (`r_warning_status`, `r_warning_count`), the L2 `scientific_provenance.json`, the report and the GUI. Recording does not muffle, change or escalate a warning; a warning is never an analysis failure. Errors keep their existing handling: the documented dispersion fallback, a shrinkage failure (`shrinkage.status: failed`), or a failed analysis.
- If shrinkage fails, `log2FoldChange_shrunk` is `NA` for every gene, `shrinkage.status` is `failed` and the error message is recorded. Raw LFC is never written into the shrunken column. The volcano and MA plots then draw raw LFC, labelled on the axis, caption and manifest as "raw log2 fold change (shrinkage failed)"; DEG calls are unaffected.

### Group means

`mean_normcount_numerator` / `mean_normcount_denominator` are the arithmetic means of `counts(dds, normalized = TRUE)` over the analysed samples in the numerator / denominator level. The normalized counts use the size factors (or tximport normalization factors) of the complete fitted object. They are descriptive statistics, not model estimates: with covariates (for example batch) or a paired design, `log2(mean_num / mean_den)` can differ from `log2FoldChange`, and can even have the opposite sign, because the model adjusts for the other terms and the means do not. Report the model LFC as the effect.

### Gene annotation

DEG annotation is independent of the enrichment selection. The OrgDb is selected from `organism.species` (`org.Hs.eg.db` for Homo sapiens, `org.Mm.eg.db` for Mus musculus; other species are not annotated) and the identifier namespace from `annotation.input_id_type`. Without an `annotation` block, the namespace is inferred only when every `gene_id` matches one unambiguous pattern (Ensembl `ENS…G…`, or all-numeric Entrez); otherwise `symbol` and `description` are `NA` and the status is `NOT_AVAILABLE`. The OrgDb package version is recorded in `contrast_metadata.json` and `scientific_provenance.json`.

- `gene_id` is the original identifier, including any version suffix, and is the unique key of every DEG table. For ENSEMBL, a trailing `.N` version is removed for the OrgDb query only (`gene_annotation.R`, the same rule used by GO/KEGG ORA and GSEA mapping).
- `symbol` = `mapIds(..., column = "SYMBOL", multiVals = "first")`; `description` = `mapIds(..., column = "GENENAME", multiVals = "first")`. Missing mappings are `NA`. The number of gene IDs with more than one distinct symbol is recorded (`one_to_many_symbol_gene_ids`). Annotation never adds or removes rows: each table has exactly one row per `gene_id`.

### L2 output schema

Per contrast (`contrasts/<contrast_id>/`):

| File | Columns / content |
|---|---|
| `all_genes.tsv` | `gene_id`, `symbol`, `description`, `baseMean`, `log2FoldChange`, `lfcSE`, `stat`, `pvalue`, `padj`, `mean_normcount_numerator`, `mean_normcount_denominator`, `log2FoldChange_shrunk` — every retained gene |
| `significant.tsv`, `upregulated.tsv`, `downregulated.tsv` | same columns; DEG rows, split by the sign of raw `log2FoldChange` |
| `contrast_metadata.json` | `numerator`/`denominator` `{level, n_samples}`, `design_formula`, `fit_method`, `deg_definition` (thresholds, rule, testing), `shrinkage` (`method`, `coef`, `contrast`, `releveled`, `relevel_verification`, `status`, `error`), `group_means` definition, `annotation` (OrgDb package/version, namespace and its source, mapping counts and rates), `software` (R, DESeq2, apeglm, ashr versions), `r_warnings` (`status` `NONE`/`WARNINGS`; `fit` and `contrast` lists of `{step, message, call, count}`) |
| `volcano.*`, `ma_plot.*`, `heatmap.*` | PNG preview (150 dpi), `<name>_300dpi.png` (profile raster dpi), TIFF (raster dpi), Cairo PDF |
| `heatmap_genes.tsv` | `gene_id`, `symbol` |

Per L2 run: `normalized_counts.tsv` — `gene_id`, `symbol`, then one column per sample (`counts(dds, normalized = TRUE)`, retained genes). Missing values are written as `NA`.

## Preranked GSEA

`analysis.enrichment` accepts an explicit list containing any of `go`, `kegg`, and `gsea`. `gsea` expands to GO preranked GSEA (BP, MF, CC) and KEGG preranked GSEA; GO and KEGG ORA are independent production modules. Public selections are deduplicated and stored in canonical `go`, `kegg`, `gsea` order, so each backend runs at most once. The historical complete four-backend spelling is normalized to this full public selection. Configuring a project for L2 (`rnaseq new --preset L2` or `rnaseq design`) selects all three by default; `[]` is the explicit opt-out (`--no-enrichment`), and any selection requires an `annotation` block with `organism` and `input_id_type`. The final report gives every method and contrast an explicit status, with `NOT_REQUESTED` for unselected methods.

## Enrichment significance (semantics v1.5)

All four enrichment analyses — GO ORA, KEGG ORA, GO GSEA and KEGG GSEA — call a gene set **significant when `p.adjust <= 0.05`**, where `p.adjust` is the Benjamini–Hochberg adjusted p-value. From v1.5.0 this is the primary, and only, enrichment significance criterion; it is fixed and not configurable.

- **clusterProfiler's return condition vs nf-rna's decision.** clusterProfiler is always called with `pAdjustMethod = "BH"`, `pvalueCutoff = 1` and `qvalueCutoff = 1` (GSEA has no q-value cutoff), so it returns every term it evaluated. nf-rna then adds a `significant` column (`TRUE` exactly when `p.adjust` is not NA and `p.adjust <= 0.05`). `pvalue`, `p.adjust` and `qvalue` are all kept. Evaluated terms are those whose gene set falls within `minGSSize`..`maxGSSize`; for ORA, a term also needs at least one foreground gene.
- **Tables.** Every all-terms table carries the `significant` column and keeps the non-significant terms: `GO_<ont>_all_terms.tsv`, `kegg_all_terms.tsv` and GSEA `all_terms.tsv`. The significant tables are the rows where it is `TRUE`: `GO_<ont>.tsv`, `kegg.tsv`, GSEA `significant.tsv`, and the `positive_enrichment.tsv` / `negative_enrichment.tsv` split by NES sign.
- **Metadata.** Every enrichment summary and its `scientific_provenance.json` record `significance_semantics: "v1.5"`, the rule and the clusterProfiler call cutoffs.
- **Legacy configuration fields.** `annotation.enrichment.go.pvalue_cutoff/qvalue_cutoff`, `kegg.ora.pvalue_cutoff/qvalue_cutoff` and `gsea`/`kegg.gsea` `pvalue_cutoff/padj_cutoff` are still accepted for schema compatibility but no longer select terms. When a project.yaml sets any of them, `rnaseq validate` and the run log emit a `legacy_enrichment_cutoff` warning naming the fields; it is never an error.

**v1.4 vs v1.5.** Up to v1.4.x, ORA significance was `pvalue <= pvalue_cutoff AND qvalue <= qvalue_cutoff` (defaults 0.05 / 0.2), with `p.adjust` reported but not used. GSEA significance was `p.adjust <= padj_cutoff AND pvalue <= pvalue_cutoff` (defaults 0.05 / 0.05). Results produced before v1.5 are not recomputed or rewritten. Their summaries have no `significance_semantics` field, which identifies them as v1.4 semantics.

## Over-representation analysis

- **Foregrounds.** These come from the L2 DEG definition, which is unchanged: `significant` is `significant.tsv`, `up` is `upregulated.tsv` (log2FoldChange > 0, higher in the numerator) and `down` is `downregulated.tsv` (log2FoldChange < 0). Gene IDs are mapped to Entrez with the OrgDb. A one-to-many source ID keeps every Entrez target, and repeated targets are deduplicated. A foreground with fewer than `annotation.minimum_mapped_foreground` mapped genes is `NOT_APPLICABLE`.
- **Universe.** The universe is the genes whose DESeq2 `padj` is not NA, that is, the genes that passed independent filtering and outlier handling. They are mapped to Entrez and then intersected with the genes of the gene-set database: the Entrez IDs annotated in that GO ontology (OrgDb), or the genes of the KEGG snapshot. Each foreground must be a subset of the padj-non-NA genes. For every foreground and ontology/database the summary records `padj_non_na_genes`, `entrez_mapped_genes`, `entrez_mapping_rate`, `database_genes`, `universe_genes`, `foreground_in_universe` and `clusterprofiler_bgratio_n`. The last is the N of clusterProfiler's `BgRatio`, which equals `universe_genes`. (Up to v1.4 the universe was every gene with a finite p-value, including genes whose padj was NA.)
- **Calculation.** GO uses `enrichGO(ont = BP/MF/CC, keyType = "ENTREZID", minGSSize = 10, maxGSSize = 500)`; these are clusterProfiler's defaults, now passed explicitly. KEGG uses `enricher(TERM2GENE, TERM2NAME)` from the KEGG snapshot with `annotation.enrichment.kegg.ora.min_gs_size/max_gs_size` (default 10/500).

## Preranked GSEA calculation

- **Ranking metric.** The ranking metric is the DESeq2 Wald statistic `stat` of the contrast, unmodified; `log2FoldChange_shrunk` is never used. There is no DEG, p-value, adjusted-p-value or effect-size prefilter. Genes with a non-finite `stat` are excluded.
- **Mapping and ties.** Gene IDs are mapped to Entrez. When several genes map to the same Entrez ID, the one with the largest `|stat|` is kept (ties by gene_id). The list is sorted by decreasing `stat`, then Entrez ID, so ties are ordered deterministically.
- **Calls.** GO uses `gseGO(ont = BP/MF/CC)`; KEGG uses `GSEA(TERM2GENE, TERM2NAME)` from the same KEGG snapshot as ORA. Gene-set size is `min_gs_size`..`max_gs_size` (default 10..500).
- **Reproducibility.** `set.seed(seed)` runs immediately before each call (`annotation.enrichment.gsea.seed` / `kegg.gsea.seed`, default 1), and the call uses `seed = FALSE` so that DOSE does not re-seed (v1.6.1), with `nPermSimple = 10000` and `eps = 0`. All of these, and the ranking metric name, are written to the summary (`gsea`) and provenance. With the same snapshot, packages and input, repeated runs give identical tables.

ORA and GSEA summaries retain their scientific outcome: `SUCCESS`, `NO_SIGNIFICANT_TERMS`, `NOT_APPLICABLE`, `BLOCKED`, `FAILED`, or for KEGG `NETWORK_UNAVAILABLE`. The corresponding structured provenance records `COMPLETED` for successful execution.

## Gene-set sources and versions

- **GO.** GO annotations come from the OrgDb (`org.Hs.eg.db`, `org.Mm.eg.db`) through clusterProfiler. No GO snapshot is made: GO results depend on the package versions, which every GO summary and provenance records (OrgDb, GO.db, clusterProfiler).
- **Where KEGG snapshots live.** KEGG gene sets come from a local snapshot of KEGG REST, stored outside the repository in `${XDG_CACHE_HOME:-~/.cache}/nf-rna/genesets/kegg/<organism code>/<snapshot_id>/`:
  - `TERM2GENE.tsv`: from `https://rest.kegg.jp/link/pathway/<org>`, pathway ID and Entrez gene ID.
  - `TERM2NAME.tsv`: from `https://rest.kegg.jp/list/pathway/<org>`, without the trailing ` - <organism>` suffix.
  - `manifest.json`: source URLs, UTC retrieval time, KEGG release information (`/info/pathway` "Last update", plus the raw `/info/<org>` and `/info/pathway` text), the SHA-256 and size of both files, and the pathway and gene counts.
- **Pin the snapshot for publication.** An unpinned project uses the latest local snapshot, so creating a newer snapshot changes later KEGG results. `rnaseq validate` and the run log (`logs/rnaseq.log`) name the snapshot an unpinned project uses and mark it "not pinned". For results that will be published, set `annotation.enrichment.kegg.snapshot: <snapshot ID>`; a pinned snapshot that does not exist locally is a validation error.
- **Which snapshot is used.** It is downloaded with KEGGREST and published atomically. A project can pin one with `annotation.enrichment.kegg.snapshot`. Otherwise the latest existing snapshot (by retrieval time) is used, and one is created only when none exists. A lock directory makes KEGG ORA and KEGG GSEA, which may start together, share the same snapshot.
- **Verification and recording.** Before use, both SHA-256 values are recomputed and must match the manifest. The snapshot ID, selection, release and SHA-256 values are written to every KEGG summary and provenance, and the manifest is copied as `kegg_snapshot_manifest.json`.
- **No live KEGG queries.** `enrichKEGG`/`gseKEGG` and `use_internal_data = TRUE` are not used. If no verified snapshot is available, for example because the download fails, the analysis is `NETWORK_UNAVAILABLE`, records `reproducible: false` ("NOT_FULLY_REPRODUCIBLE") and produces no KEGG result.
- **Provider name.** `resource_provider: kegg_rest_snapshot` is the default; the pre-v1.5 value `online_kegg_rest_via_clusterprofiler` is accepted and behaves identically.

## Enrichment output schema (v1.5 additions)

- **Long member tables.** There is one row per pathway × `gene_id` for the significant pathways. Columns: `pathway_id`, `pathway_description`, `gene_id`, `symbol`, `entrez_id`, `log2FoldChange`, `padj`, `pathway_p.adjust`; GSEA adds `NES` and `leading_edge` (TRUE when the gene is in `core_enrichment`).
  - Files: GO ORA `GO_<ont>_members.tsv`, KEGG ORA `kegg_members.tsv`, GSEA `members.tsv`.
  - `gene_id` is the key; `symbol` is the L2 annotation, falling back to the OrgDb symbol of the Entrez mapping.
  - ORA rows are the foreground genes in the pathway. GSEA rows are every ranked gene of the gene set, with the gene_id kept by the ranking.
  - The earlier member tables (`GO_<ont>_gene_members.tsv`, `kegg_gene_members.tsv`, `gsea_core_members.tsv`, `core_members.tsv`) are unchanged.
- **GO redundancy reduction.** `clusterProfiler::simplify(cutoff = 0.7, by = "p.adjust", select_fun = min, measure = "Wang")` is applied to the significant terms only. It writes `GO_<ont>_simplified.tsv` (ORA) or `significant_simplified.tsv` (GSEA), plus `dotplot_<ont>_simplified` / `dotplot_simplified`. The full results are unchanged.
- **GSEA running-score plots.** `running_score/{positive,negative}_<n>_<ID>` are drawn with `enrichplot::gseaplot2` for up to 5 significant gene sets per NES sign, chosen by smallest `p.adjust` (ties by ID). They use the analysis' own gseaResult and therefore its ranked list. Each title states the ID, description, NES and adjusted P.
- **Figure formats.** Every enrichment figure is a vector PDF plus a PNG at the profile's raster dpi (`_300dpi.png`) and the 150 dpi report preview; dotplots also keep the TIFF.
- **No empty figures.** When nothing is significant, no figure is drawn and the summary's `figures` record says why.

Annotation mapping QC is shared by all four enrichment backends: GO and KEGG ORA (evaluated on the Entrez mapping of the padj-non-NA genes) and GO and KEGG preranked GSEA (evaluated on the finite-statistic ranked source genes). `annotation.mapping_warning_rate` defaults to 0.70 and `annotation.minimum_mapping_rate` defaults to 0.50. A mapping rate at or above the warning threshold is `PASS`; a rate from the blocking threshold up to (but not including) the warning threshold is `WARNING` and the backend still executes; a rate below the blocking threshold is `BLOCKED` before formal ORA or GSEA. Every backend records finite source genes, mapped and unmapped source genes, mapping rate, final ranked genes, one-to-many mappings, duplicate-target collapse, both thresholds, identifier namespaces, annotation database/version, and resulting QC status. A warning never invalidates the technical report. Existing projects that contain only the legacy `minimum_mapping_rate` are normalized as `mapping_warning_rate` with a blocking threshold of `min(0.50, legacy value)` when validated and frozen, so they never become stricter merely by migration.

## Figures

Figures are views of the result tables; no figure changes a calculation, a threshold or a result set. Every R backend draws through one policy (`src/rnaseq/r/figures.R`) and writes a `figure_manifest.json` beside its outputs.

**Output.** Each figure is drawn directly, at its final physical size, on three devices: a PNG preview (150 dpi) for the HTML report, a raster TIFF (300 dpi, RGB, LZW, opaque white background, no alpha channel) and a Cairo PDF (vector text and lines, one page). Nothing is resampled from another output. Sizes are snapped to whole points because Cairo PDF pages are whole points; the manifest records the exact millimetres, pixels and dpi. The L2 volcano, MA and DEG heatmap additionally write `<name>_300dpi.png` at the profile's raster dpi. Delivery places them under `figures/png/`, `figures/tiff_300dpi/` and `figures/pdf/`; manifests go to `methods_and_versions/`. A `.pdf` is not a guarantee that every element is vector.

**Style.** The project option `figures: {preset, profile}` selects the figure style: `standard` (alias `bioinfo`) or `legacy`. New projects write `standard`; a project without the section renders with `legacy`. A style changes appearance only (colours, fonts, sizes, layout), never a calculation, a threshold, a result table or a file name; see `docs/figure-standard.md`.

**Profiles.** The layout targets live in one table: `general` (used by production; 85/180 mm columns, at most 180 × 220 mm, 8/7/6 pt text, 300 dpi), and `nature`, `plos` and `elsevier`, which record those publishers' published artwork targets (for example 89/183 mm columns and 5–7 pt for Nature research figures; at most 190.5 × 222.3 mm and 8–12 pt for PLOS ONE; 300/500/1000 dpi halftone/combination/line art for Elsevier). They are layout targets, not proof that a journal will accept a figure. Production figures fit the `general` profile; a figure that would exceed it is reported as a warning rather than resized after drawing.

**Typography.** The preferred family is Arial. The supported Linux/WSL2 runtime (the locked Conda environment) does not contain Arial, so nf-rna uses the explicit fallback **DejaVu Sans**, which that environment ships, and records the requested family, the family used, the font file and the reason in the manifest. Font availability is checked with `systemfonts` (present in the locked runtime); it refuses to draw if neither family is installed, and if `systemfonts` itself is unavailable it uses DejaVu Sans and records the choice as unverified. Theme and heatmap sizes are set in points, ggplot2 text geoms use `size.unit = "pt"`, and the family is set on every device so heatmap text uses it too.

**Colour.** Categorical groups get a deterministic named mapping in sorted level order: Okabe–Ito for up to eight levels, otherwise the HCL qualitative "Dark 3" palette — never an interpolation between two colours. PCA encodes the group by both colour and shape (shape up to 12 levels; beyond that colour only, as recorded); the batch PCA colours by batch and keeps the group as the shape. Group colours follow the design's level order, so a group has the same colour in the PCA and the DEG heatmap annotation. Volcano directions use vermillion (higher in the numerator), blue (higher in the denominator) and grey. Diverging scales (ColorBrewer RdBu, midpoint 0) are used only where zero is meaningful: row z-scores, and sample correlation when a negative correlation occurs. An all-positive correlation matrix uses a sequential scale from the observed minimum to 1. Sequential scales for magnitudes start at zero, and a degenerate range (such as an all-zero matrix) is widened so colour breaks stay valid.

**DEG heatmap.** Rows are the top `heatmap_top_n` (50) significant genes by padj (ties by gene ID); the complete significant set stays in `significant.tsv`. Values are blind VST row z-scores over the numerator and denominator samples. Missing values stay missing: they are drawn in grey, never as 0. A gene whose z-score is undefined (constant, or fewer than two finite values) is not drawn and is listed in the manifest; `heatmap_zscores.tsv` keeps every selected gene's unclipped z-scores with `NA` preserved. Colours saturate at |z| ≥ 2 (legend ≤−2 / ≥2) for display only; clustering (complete linkage, Euclidean, as before) uses the unclipped values and is switched off, and recorded, if missing values leave a distance undefined. If no selected gene has a defined z-score, no heatmap is drawn: `heatmap_status` is `NOT_APPLICABLE` and `heatmap_reason` says why. Row labels are gene symbols (gene ID when the symbol is NA; a repeated symbol gets its gene ID appended) and shrink from 7 pt to no less than 6 pt; if the rows still do not fit, labels are hidden (and recorded as hidden) rather than enlarging the figure.

**Volcano.** The x-axis is `log2FoldChange_shrunk` (symmetric, labelled with the shrinkage method and `numerator / denominator`); the y-axis is −log10 of DESeq2's BH-adjusted P. Colour follows the unchanged DEG definition on the raw LFC. A dashed line marks the padj threshold; because the |log2FC| threshold applies to the raw LFC, it is stated in the caption rather than drawn on the shrunken axis (the dashed |log2FC| lines are drawn only in the raw-LFC fallback after a shrinkage failure, which the axis and caption state). A padj that underflows to 0 is drawn at the smallest positive double and marked with a triangle; the table keeps the reported value. Genes without a padj or a finite shown LFC are not drawn and are not counted in the legend. The top 10 drawn DEGs by padj (ties by gene ID) are labelled with their symbol, or gene ID when the symbol is NA.

**MA plot.** x = log10(`baseMean`), y = `log2FoldChange_shrunk` (raw LFC fallback as for the volcano), coloured by the DEG definition, with a line at 0.

**Enrichment dotplots.** ORA plots show the top 15 significant terms by adjusted P (an adjusted P of 0 is drawn at the same display floor as the volcano and marked); GSEA plots show the top 15 significant terms by |NES|, ordered by NES, with a zero line and an axis stating which contrast group positive NES favours (ranking: DESeq2 Wald statistic). GSEA colour encodes adjusted P on fixed limits [0, 0.05] (the v1.5 significance cutoff), so all GSEA dotplots are comparable. GO dotplots have a `_simplified` companion, and GSEA results add running-score plots (see the enrichment output schema). Long descriptions are wrapped and, beyond two lines, truncated; each plot's `*_source_data.tsv` keeps the full ID and description of every plotted term.

## Reference and provenance

Reference identity is explicit: species, assembly/release, checksums, selected-backend index status, and optional index-build provenance are frozen for a run. The frozen execution contract carries the git commit of the executed nf-rna source as `execution.downstream_runtime.source_revision`, together with the locked Conda runtime identity (lock and wheel SHA-256, bundled R-script inventory); every R module receives those values and records them in its structured provenance. A run is refused when no source commit is available or the source checkout has uncommitted changes. (The historical Docker runtime of v1.2.1 and earlier carried the OCI image revision label as `execution.source_revision` instead.) Every production R stage additionally emits `scientific_provenance.json` and `r_session_info.txt`, recording script/runtime identity, R and package versions, identifier namespace, design/contrasts, normalized effective configuration, output schema, UTC timestamp, and its scientific result summary. Prebuilt Salmon and HISAT2 indexes are supported directly; their checksums and declared strategies are validated without requiring nf-rna to have created them. HISAT2 records annotation-derived splice sites and its genome-index inputs; it never reuses a Salmon index. Interactive project creation may discover a compatible production reference from the user's machine-local registry, but it validates the manifest before offering it and writes the resolved root/manifest into the project. No reference is downloaded or silently migrated during validation/planning. External runtime behavior and versions are captured in provenance without changing the scientific configuration.

**Actual upstream threads (v1.6.1).** HISAT2's placement of multimapping reads and Salmon's multithreaded EM depend on the thread count, so a rerun with other thread counts reproduces counts only up to these differences (Phase 3 validation, N-1 and N-3). Every Nextflow trace of a FASTQ run therefore adds the `process`, `tag`, `cpus` and `memory` columns; trace fields never change scheduling.

From the run's own trace, `run_provenance.yaml` `execution_threads` records, per process, the CPUs every task used, after any `resourceLimits` cap: HISAT2 `-p`, featureCounts `-T` and Salmon `--threads` are all `task.cpus`. The processes are fastp, HISAT2, featureCounts and the strand evidence on the HISAT2 route, and Trim Galore/fastp and Salmon on the Salmon route.

The record also includes Salmon's `--threads` from each sample's `cmd_info.json` and the tool versions per process. It is copied to the frozen contract (`upstream_qc.execution_threads`), the report and the GUI. Runs without a `cpus` trace column, or with a reused upstream, are recorded as `NOT_RECORDED` with the reason; no thread count is inferred from `execution.max_cpus`.

A Salmon index is usable only with a PASS transcript-ID contract: its transcriptome is derived from the genome FASTA and GTF as nf-core/rnaseq does, every transcript name equals exactly one GTF `transcript_id` with one `gene_id` (no identifier normalization), and every GTF transcript on the genome is indexed, because nf-core builds tx2gene by exact name matching. Managed Salmon and HISAT2 reference preparation is host-native only. It validates the controlled builder executable/version set before staging, records non-scientific host-builder facts separately from biological/index identity, and atomically publishes only validated artifacts. Existing valid Docker-built manifests remain readable and are not relabelled or rebuilt merely because their historical builder mode differs.

Schema 1.2 adds an explicit HISAT2 strategy discriminator. Legacy
graph-embedded splice-site manifests remain readable. A
`genome_only_runtime_splicesites` reference instead requires a non-empty,
checksum-bound GTF-derived splice-site file and passes it at alignment runtime;
the index-builder and runtime-aligner versions are distinct provenance fields.
An unvalidated version combination is not execution-ready. Its tiny Human
FASTQ smoke fixture is technical software validation only, never biological
evidence.

Production acceptance is opt-in with `reference.acceptance: production`. It requires `reference.source: local`, manifest schema 1.1, `purpose: production`, verified local asset hashes, and a complete index bound to the same FASTA/GTF identity. Synthetic/fixture identities and legacy manifests are rejected. To migrate schema 1.0 safely, review the provider/release/assembly identity and every local checksum, change the schema to 1.1, and deliberately set `purpose`; never add `purpose: production` merely to satisfy validation. Optional source URLs/accessions are provenance only. An optional upstream checksum is reported as verified only when it is supplied and matches the local asset.

The first intended human production identity is `Homo sapiens`, Ensembl release 116, GRCh38 patch p14. This repository configures the identity contract only; it does not contain or download those assets.
