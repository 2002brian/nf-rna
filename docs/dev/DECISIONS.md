# Development decisions

Non-blocking implementation choices, recorded as they are made.

## v1.5.0 Phase 1 — DEG scientific upgrade (2026-10-30)

| Decision | Reason |
|---|---|
| The downstream locks were regenerated with the documented `conda env create --dry-run --json` method, from a spec that pins every v1.4.1 lock entry exactly and adds `bioconda::bioconductor-apeglm` and `conda-forge::r-ashr` (project `CONDARC=/home/brian/nf-rna-work/conda/condarc.yml`, strict conda-forge > bioconda). | A free re-resolution would have moved DESeq2 and other packages. The pinned resolution is purely additive (15 new packages per platform, no existing entry changed), so v1.4.1 statistics run on the identical DESeq2 build. `nf-rna-downstream.yml`, both locks and `SHA256SUMS` are updated together. |
| The osx-arm64 lock was resolved on linux with `CONDA_SUBDIR=osx-arm64 CONDA_OVERRIDE_OSX=13.0`. | Same method as the existing lock; it is a resolution, not a native qualification. |
| Relevel verification uses `all.equal` (mean relative difference ≤ 1e-4) on raw `log2FoldChange` and `pvalue`, with an identical NA pattern. | A reparametrized IRLS refit differs from the original by ~1e-6 relative (measured on the paired example: 5.9e-7 LFC, 1.2e-6 p-value); 1e-6 rejected valid copies. 1e-4 is still far below any reportable difference. |
| The contrast interface always names two levels of one factor, so ashr is reached only when the releveled copy fails verification (for example after DESeq2 replaced outlier counts for ≥7 replicates, which a plain `nbinomWaldTest` on the copy does not reproduce). | Follows the specified order; the copy is never used unverified. |
| Without an `annotation` block, DEG annotation infers ENSEMBL (`^ENS…G\d+(\.\d+)?$`) or ENTREZID (all numeric) only when every gene ID matches; otherwise symbol/description stay NA. | DEG annotation must not depend on the enrichment selection, but SYMBOL or mixed namespaces cannot be inferred safely. |
| `description` is OrgDb `GENENAME` (`mapIds`, `multiVals = "first"`). | The only description field in both supported OrgDbs. |
| The ENSEMBL version-stripping rule moved to `gene_annotation.R` (`nf_rna_annotation_key`) and GO/KEGG/GSEA call it. | One mapping rule shared by DEG annotation and enrichment; the regex and behaviour are unchanged. |
| The 300 dpi PNG is an extra file `<name>_300dpi.png` (format `png_raster`, the profile's raster dpi) for the volcano, MA plot and DEG heatmap; the 150 dpi `<name>.png` preview is unchanged. | Keeps existing file names and the HTML report unchanged; delivery copies it to `figures/png/` automatically. |
| On the shrunken volcano axis the dashed |log2FC| lines are omitted; the rule (on raw LFC) is in the caption. They are drawn in the raw-LFC fallback. | A raw-LFC threshold drawn on a shrunken axis would misplace the cut-off. |
| Volcano labels: top 10 drawn DEGs by padj (ties by gene_id), `ggrepel` seed 1. | Simple, deterministic; ggrepel is already locked. |
| The MA plot is drawn with the volcano figure scope and colours. | Same encoding; no new style values. |
| The DEG heatmap keeps the L1 blind VST row z-scores; rows are labelled by symbol, with `symbol (gene_id)` for repeated symbols. | The existing transformation is retained as specified; repeated labels would be ambiguous. |
| `normalized_counts.tsv` is written once per L2 run (`downstream/l2/`); the L1 `normalized_counts.tsv` is unchanged. | Normalized counts do not depend on the contrast. |
| `contrast_metadata.json` is delivered under `methods_and_versions/` with the other scientific JSON records. | Delivery copied only allow-listed JSON names. |
| Acceptance ran `prepare_l2`/`execute_l2` directly in a v1.4.1 worktree and in the feature branch, on `examples/gsea_two_group` (plus a reversed contrast) and `examples/paired_two_group`. | `rnaseq run` refuses uncommitted sources and would provision runtimes; the direct call is the same L2 code path for raw-count input, without FASTQ processing. |
| The package version stays 1.4.1 and CHANGELOG.md is not edited during Phase 1. | No release in this phase; the release-metadata test requires an `X.Y.Z — Unreleased` heading to match a bumped version, which belongs to the release step. |
