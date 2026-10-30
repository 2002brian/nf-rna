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

## v1.5.0 Phase 2 — acceptance prerequisites (2026-10-30)

| Decision | Reason |
|---|---|
| The enrichment baseline runs the four R backends with configs shaped like `workflow_support.enrichment_config`, not through `rnaseq.go/gsea/kegg.execute_*`. | The legacy Python entry points fail identically in v1.4.1 and HEAD (no contrast numerator/denominator for the ORA subtitle; `summary["mapping"]` no longer written). The Nextflow config is the supported path. The legacy path is left unchanged. |
| airway lives in a separate env, `/home/brian/nf-rna-work/conda/airway-testdata-env` (project `CONDARC`). | The spec forbids adding it to the downstream runtime; a conda env is reproducible without touching a global R library. |
| The airway project uses `design.type: paired_two_group` with `pair_id: cell` (formula `~ cell + condition`). | Each cell line has exactly one trt and one untrt sample, so the cell line is the pairing unit; the formula is the one specified. |
| airway counts are written unchanged (all 63,677 genes, Ensembl IDs without version); the existing L2 all-zero/low-total filter applies at analysis time. | Test data must not be modified; the filter is part of the pipeline under test. |
| The airway `project.yaml` declares an ENSEMBL → org.Hs.eg.db annotation block with `analysis.enrichment: []`. | DEG annotation uses the declared namespace; enrichment selection is left to Phase 2 acceptance. |

## v1.5.0 Phase 2 — enrichment scientific upgrade (2026-10-30)

| Decision | Reason |
|---|---|
| The significance cutoff (0.05, BH) and `nPermSimple = 10000` are constants in `enrichment_common.R`, not configuration. The `pvalue_cutoff`/`qvalue_cutoff`/`padj_cutoff` fields stay in the schema but are ignored; the summary says so. | The rule is fixed by the specification; removing the fields would make existing project.yaml files invalid. |
| GO ORA passes `minGSSize = 10, maxGSSize = 500` explicitly. | These are the clusterProfiler defaults it used implicitly before, so results do not change, and now they are recorded. |
| GO database genes per ontology are the Entrez IDs with a direct GO annotation of that ontology in the OrgDb (`select(keytype = "ONTOLOGY")`). The universe passed to `enrichGO` is the padj-non-NA mapped genes ∩ these genes. The summary also records clusterProfiler's BgRatio N. | Every gene annotated anywhere in an ontology has a direct annotation there. On airway, BgRatio N equals the computed universe for every ontology and foreground, so the definition is verified rather than assumed. |
| Member tables (`*_members.tsv`) cover significant pathways only. For ORA they list the foreground genes in each pathway; for GSEA, every ranked gene of the set with a `leading_edge` flag. The old member tables are kept. | All-terms member tables would be very large (GO BP × up to 500 genes). Keeping the old files keeps existing delivery paths. |
| The member-table `symbol` comes from the L2 `all_genes.tsv` `symbol` column (Phase 1 annotation), with the OrgDb symbol of the Entrez mapping as fallback. | One symbol per gene_id across DEG and enrichment tables. |
| Running-score plots: one figure per gene set, per GO ontology and for KEGG, at most 5 per NES sign, chosen by smallest p.adjust (ties by ID). They get PNG preview, 300 dpi PNG and PDF, but no TIFF. | Per-ontology plots match the directory layout. The preview PNG is what the figure review tools (`compare_runs.py`) and the report expect. |
| Running-score plots are not a standard-preset figure type. | Migrating a new figure type needs owner review under the figure standard; they use the run's context and are identical across presets (verified by `compare_runs` in `test_figure_dotplots`). |
| simplify runs on the significant terms only (`by = "p.adjust", select_fun = min, measure = "Wang"`). A result with fewer than 2 significant terms is copied unchanged, with the reason recorded. | The full result stays untouched; Wang needs no information content. |
| KEGG snapshot ID: `<org>_kegg<pathway last-update yyyymmdd>_<UTC retrieval stamp>`. "Latest" means the latest `retrieved_at_utc` in a valid manifest. A `.lock` directory serialises resolution, and the snapshot is staged then renamed into place. | KEGG REST no longer reports a "Release" number; `/info/pathway` "Last update" is the release information it gives, and the raw `/info` texts are stored too. KEGG ORA and GSEA run as parallel Nextflow tasks and must share one snapshot. |
| A pinned snapshot is `annotation.enrichment.kegg.snapshot` (validated as a plain directory name). `resource_provider` accepts the old `online_kegg_rest_via_clusterprofiler` value as an alias. | Existing projects stay valid without edits. |
| The service preflight probes KEGG REST only when no usable snapshot exists. A pinned snapshot that is missing is a preflight error. | With a snapshot, KEGG analyses no longer need the network. |
| enrichplot, KEGGREST, GOSemSim, GO.db and openssl are added to the runtime package probes, not to `nf-rna-downstream.yml`. msigdbr (optional item E) was not added. | They are already in both locks as clusterProfiler dependencies, so the lock and checksums stay unchanged; `ggrepel` set the precedent. msigdbr would need a lock regeneration and brings `babelgene` with it (not "msigdbr only"). |
| The Python `execute_*` entry points and `workflow_support.enrichment_config` both use `rnaseq.enrichment_config.backend_config`. | Fixes the missing contrast levels; the two paths can no longer drift. |
| The version is now the development version `1.5.0` (`__init__.py`, CITATION.cff without `date-released`) and CHANGELOG has `1.5.0 — Unreleased`. | The specification requires a v1.5.0 development CHANGELOG entry, and the release-metadata guard accepts an `Unreleased` heading only with the matching development version. No tag or release. |
| The offline KEGG figure-review fixture writes a synthetic snapshot and replaces `enricher()`/`GSEA()` (GSEA returns a synthetic `gseaResult`). | The network touch points moved from `url()`/`enrichKEGG`/`gseKEGG` to the snapshot download; the fixture must still need no network. |
