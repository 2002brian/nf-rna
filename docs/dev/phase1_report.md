# nf-rna v1.5.0 Phase 1 report — DEG scientific upgrade

Date: 2026-10-30. Branch `feature/v1.5.0`, created from `v1.4.1` (`f69058e`). Phase 2 (enrichment upgrade), GUI work and the v1.5.0 release are out of scope.

## Scientific definitions

**The DEG statistical definition did not change.** Same single DESeq2 fit, same `results(contrast = c(factor, numerator, denominator), independentFiltering = TRUE)`, same predicate `padj < thresholds.padj AND abs(raw log2FoldChange) >= thresholds.abs_log2fc`. `log2FoldChange_shrunk` is display and ranking only. On the acceptance data, `baseMean`, `log2FoldChange`, `lfcSE`, `stat`, `pvalue` and `padj` are byte-identical to v1.4.1, and so are the significant, up and down gene sets (see below). The heatmap transformation (L1 blind VST, row z-score) is unchanged.

## Modified files

| File | Change |
|---|---|
| `src/rnaseq/r/l2_analysis.R` | Normalized counts, annotation, group means, LFC shrinkage (apeglm → verified relevel + apeglm → ashr; failure → NA), new columns, `normalized_counts.tsv`, `contrast_metadata.json`, shrunken volcano with labels, MA plot, symbol heatmap labels, extra 300 dpi PNG |
| `src/rnaseq/r/gene_annotation.R` (new) | Shared OrgDb key rule (ENSEMBL version removed for the query only) and DEG annotation (`mapIds`, `multiVals = "first"`) |
| `src/rnaseq/r/go_analysis.R`, `kegg_analysis.R`, `gsea_analysis.R` | Use the shared key rule (identical behaviour) |
| `src/rnaseq/r/figures.R` | `nf_rna_save_figure` accepts the opt-in `png_raster` format |
| `src/rnaseq/l2.py` | New result header; `ma_plot.png` and `contrast_metadata.json` required; passes `gene_annotation` settings |
| `src/rnaseq/workflow_support.py` | Nextflow `l2-config` passes the same `gene_annotation` settings |
| `src/rnaseq/models.py` | `deg_gene_annotation()` (OrgDb from `organism.species`, namespace from `annotation.input_id_type`) |
| `src/rnaseq/downstream.py`, `downstream_runtime.py` | apeglm and ashr in the R runtime checks |
| `src/rnaseq/service.py` | Delivers `contrast_metadata.json` under `methods_and_versions/` |
| `workflow/envs/nf-rna-downstream.yml`, `workflow/envs/locks/*.lock.yml`, `SHA256SUMS` | apeglm 1.32.0 and ashr 2.2-63 (additive, see DECISIONS) |
| `tests/test_l2.py`, `tests/test_doctor_native.py`, `tests/test_installed_runtime.py` | New outputs; new linux-64 lock SHA-256 |
| `docs/scientific_contract.md`, `docs/figure-standard.md`, `docs/dev/DECISIONS.md` | Contract, figure notes, decisions |

## Final output schema

Per contrast, `downstream/l2/contrasts/<contrast_id>/`:

- `all_genes.tsv`, `significant.tsv`, `upregulated.tsv`, `downregulated.tsv`:
  `gene_id, symbol, description, baseMean, log2FoldChange, lfcSE, stat, pvalue, padj, mean_normcount_numerator, mean_normcount_denominator, log2FoldChange_shrunk`
- `contrast_metadata.json`: `contrast_id, factor, numerator{level,n_samples}, denominator{level,n_samples}, design_formula, fit_method, deg_definition{padj_below, abs_log2fc_at_least, rule, testing}, shrinkage{method, coef, contrast, releveled, relevel_verification, status, error, column, use}, group_means{columns, definition}, annotation{status, orgdb_package, orgdb_version, input_id_type, input_id_type_source, version_suffix, one_to_many_rule, genes, symbol_mapped, description_mapped, symbol_mapping_rate, description_mapping_rate, one_to_many_symbol_gene_ids}, software{R, DESeq2, apeglm, ashr}`
- `volcano`, `ma_plot`, `heatmap`: `.png` (150 dpi preview), `_300dpi.png`, `.tiff` (300 dpi), `.pdf` (Cairo vector)
- `heatmap_genes.tsv`: `gene_id, symbol`; other existing files unchanged

Per run: `downstream/l2/normalized_counts.tsv`: `gene_id, symbol, <sample columns>`.

## Runtime

The linux-64 and osx-arm64 locks were re-resolved with every v1.4.1 entry pinned plus the two new packages; both diffs are insertion-only (15 packages each). New linux-64 lock SHA-256 `0e721ff3cf1983f6f4411a64279dd1cd67e5f9d9ff266dd46f873b05c6540213`. A prefix was created from the new lock with the production command (`CONDA_NO_PLUGINS=true CONDA_SOLVER=classic conda env create --prefix <scratch>/r-prefix --file workflow/envs/locks/nf-rna-downstream-linux-64.lock.yml`, `CONDARC=/home/brian/nf-rna-work/conda/condarc.yml`): DESeq2 1.50.2 (same build as v1.4.1), apeglm 1.32.0, ashr 2.2.63, org.Mm.eg.db 3.22.0. The production prefix under `~/.cache/nf-rna/runtime/prefixes/` is provisioned by the next `rnaseq run` from this commit (the name changes because the lock SHA-256 changed). No reference was rebuilt.

## Acceptance

Test data (no FASTQ pipeline): `examples/gsea_two_group` (111 Entrez mouse genes, 3 vs 3, `~ condition`, contrast `Treatment_vs_Control`, plus an added reversed contrast `Control_vs_Treatment` to exercise the relevel path) and `examples/paired_two_group` (`~ patient + condition`, `Metastasis_vs_Primary`, gene-wise dispersion fallback). Copies of the same project directories were run by both versions.

```bash
git worktree add --detach $S/wt-v1.4.1 v1.4.1
# v1.4.1 code + v1.4.1 locked prefix
PATH=$OLD_PREFIX/bin:$PATH PYTHONPATH=$S/wt-v1.4.1/src python run_l2.py accept/old/gsea_two_group accept/old/paired_two_group
# feature branch + new locked prefix
PATH=$S/r-prefix/bin:$PATH PYTHONPATH=src python run_l2.py accept/new/gsea_two_group accept/new/paired_two_group
```

(`run_l2.py` calls `execute_l2(prepare_l2(validate_project(project), run_id=None))`.)

```
rnaseq from .../wt-v1.4.1/src/rnaseq/__init__.py 1.4.1
accept/old/gsea_two_group SUCCESS [('Treatment_vs_Control', 60), ('Control_vs_Treatment', 60)]
accept/old/paired_two_group SUCCESS [('Metastasis_vs_Primary', 1)]
rnaseq from /home/brian/projects/nf-rna/src/rnaseq/__init__.py 1.4.1
accept/new/gsea_two_group SUCCESS [('Treatment_vs_Control', 60), ('Control_vs_Treatment', 60)]
accept/new/paired_two_group SUCCESS [('Metastasis_vs_Primary', 1)]
```

**1–3, 5. Equivalence, uniqueness, mapping, shrinkage** (`compare.py`: compares the six statistics per gene as strings and floats, NA patterns, significant/up/down gene sets; counts rows per `gene_id`):

```
gsea_two_group/Control_vs_Treatment: genes=111 unique_gene_id=True stats_byte_identical=True max_abs_diff=0 significant_old=60 new=60 same=True up_same=True down_same=True symbol=111/111 description=111/111 shrink=success/apeglm coef=condition_Control_vs_Treatment releveled=True shrunk!=raw=111/111
gsea_two_group/Treatment_vs_Control: genes=111 unique_gene_id=True stats_byte_identical=True max_abs_diff=0 significant_old=60 new=60 same=True up_same=True down_same=True symbol=111/111 description=111/111 shrink=success/apeglm coef=condition_Treatment_vs_Control releveled=False shrunk!=raw=111/111
paired_two_group/Metastasis_vs_Primary: genes=3 unique_gene_id=True stats_byte_identical=True max_abs_diff=0 significant_old=1 new=1 same=True up_same=True down_same=True symbol=0/3 description=0/3 shrink=success/apeglm coef=condition_Metastasis_vs_Primary releveled=True shrunk!=raw=3/3
```

Mapping success: gsea example (ENTREZID, org.Mm.eg.db 3.22.0) symbol 111/111 (100%), description 111/111 (100%). The paired example uses synthetic IDs (`GeneA`…), so no namespace can be inferred: 0/3, status `NOT_AVAILABLE` with the reason recorded. An ENSEMBL check of `gene_annotation.R` (versioned IDs, one ID with two symbols, one unknown ID):

```
                gene_id symbol
1  ENSMUSG00000030359.7   Pzp2
2    ENSMUSG00000020804  Aanat
3 ENSMUSG00000025375.12   Aatk
4  ENSMUSG00000000562.1 Adora3        # OrgDb symbols: Adora3,Tmigd3 -> first
5  ENSMUSG99999999999.1   <NA>
status SUCCESS, symbol_mapped 4/5, one_to_many_symbol_gene_ids 1; inferred namespace without annotation block: ENSEMBL
```

All three shrinkage branches ran: direct coefficient (`Treatment_vs_Control`), verified relevel (`Control_vs_Treatment`, paired contrast), and ashr (the paired contrast before the relevel tolerance was widened from 1e-6 to 1e-4: `"method": "ashr", "status": "success"`).

**Shrinkage failure** (scratch copy of the R scripts with `lfcShrink` overridden to `stop()`, same backend config):

```
shrunk all NA: True significant: 60
{'method': 'apeglm', 'coef': 'condition_Treatment_vs_Control', ..., 'status': 'failed', 'error': 'forced shrinkage failure (validation)', ...}
volcano statistics.fold_change: DESeq2 unshrunken log2FoldChange (fallback: shrinkage failed)
```

The fallback volcano axis reads "raw log2 fold change (shrinkage failed)" and the caption "LFC shrinkage: failed, raw LFC shown".

**4. Independent recomputation** (`estimateSizeFactorsForMatrix` on the raw counts after the L1 filter, then means per level, compared with `normalized_counts.tsv` and `all_genes.tsv` of `Treatment_vs_Control`):

```
11287 Pzp2: norm max|diff|=2.84e-13  mean_num out=211.145592 recomputed=211.145592  mean_den out=19.947374 recomputed=19.947374
11298 Aanat: norm max|diff|=3.41e-13  mean_num out=213.495675 recomputed=213.495675  mean_den out=20.961769 recomputed=20.961769
11471 Actl7b: norm max|diff|=4.55e-13  mean_num out=39.212464 recomputed=39.212464  mean_den out=247.820340 recomputed=247.820340
```

**6. Figures:**

```
volcano.png:        PNG image data, 1062 x 708, 8-bit/color RGB
volcano_300dpi.png: PNG image data, 2125 x 1417, 8-bit/color RGB
ma_plot_300dpi.png: PNG image data, 2125 x 1417, 8-bit/color RGB
heatmap_300dpi.png: PNG image data, 2125 x 2325, 8-bit/color RGB
pdfinfo volcano.pdf / ma_plot.pdf / heatmap.pdf: Pages 1; 510 x 340 / 510 x 340 / 510 x 558 pts; header %PDF-
```

The PNGs were also inspected visually (volcano labels by symbol, MA plot, symbol-labelled heatmap, fallback volcano).

**7. CLI and delivery:** `rnaseq validate` → `PASS WITH WARNINGS`, `READY`; `rnaseq plan` → `READY` with planning artifacts on the acceptance copy. Focused tests (locked R prefix on `PATH`, `NF_RNA_REQUIRE_R=1` for the R-dependent set):

```
tests/test_l2.py test_downstream_runtime_contract.py test_workflow_support.py test_doctor_native.py
      test_go.py test_kegg.py test_gsea.py test_enrichment_mapping_qc.py           144 passed (after the lock-hash test update)
tests/test_doctor_native.py test_installed_runtime.py test_figures.py test_figure_presets.py
      test_figure_settings.py test_figure_review_tools.py                          220 passed
tests/test_service.py tests/test_cli.py -k "deliver or l2 or report"               32 passed
tests/test_release_metadata.py test_packaging_contract.py test_runtime_dispatch.py 26 passed
```

The full test suite, `tools/figure_review` legacy/candidate comparisons and a Nextflow `rnaseq run` were not executed in this phase.

## Unresolved issues

- The production downstream prefix for the new lock has not been provisioned through `rnaseq run`; it is created on the first run from a clean commit. osx-arm64 was resolved but not qualified natively.
- On real multi-replicate data (≥7 per group) DESeq2 may replace outlier counts; a releveled copy then fails verification and the contrast uses ashr. This is by design and recorded, but has not been exercised on real data.
- `tools/figure_review/compare_runs.py` was not run; the volcano/MA/heatmap changes were checked by eye only.
- The HTML report (`report.html`) still embeds only the volcano and heatmap previews; it does not yet show the MA plot or the new columns.
- CHANGELOG and version are left for the release step.
