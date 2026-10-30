# nf-rna v1.5.0 Phase 2 — report

Branch `feature/v1.5.0`, based on the Phase 1 commit `0a3382a`.

## Before Phase 2 (code as of `0a3382a`)

Recorded from the code, not from documentation.

| | GO ORA | KEGG ORA | GO GSEA | KEGG GSEA |
|---|---|---|---|---|
| Function | `enrichGO(ont = BP/MF/CC, keyType = "ENTREZID", readable = FALSE)` | `enrichKEGG(keyType = "ncbi-geneid", use_internal_data = FALSE)` | `gseGO(ont = BP/MF/CC, keyType = "ENTREZID")` | `gseKEGG(keyType = "ncbi-geneid")` |
| pAdjustMethod | `BH` (config `go.p_adjust_method`) | `BH` (config) | `BH` (config) | `BH` (config) |
| pvalueCutoff / qvalueCutoff in the call | 1 / 1 | 1 / 1 | 1 / — | 1 / — |
| minGSSize / maxGSSize | not passed (clusterProfiler defaults 10 / 500) | 10 / 500 (config) | 10 / 500 (config) | 10 / 500 (config) |
| nf-rna significance | `pvalue <= 0.05 AND qvalue <= 0.2` (config); p.adjust not used | `pvalue <= 0.05 AND qvalue <= 0.2` | `p.adjust <= 0.05 AND pvalue <= 0.05` | `p.adjust <= 0.05 AND pvalue <= 0.05` |
| Seed / permutations | — | — | `set.seed(1)` (config `gsea.seed`), `seed = TRUE`, `eps = 0`, fgsea default `nPermSimple` (1000) | no `set.seed`, `seed = TRUE`, `eps = 0`, default `nPermSimple` |

- **ORA gene input:** three foregrounds per contrast, read from Phase 1 L2 tables:
  `significant` (`significant.tsv`), `up` (`upregulated.tsv`) and `down`
  (`downregulated.tsv`). These use the unchanged DEG definition `padj <
  threshold AND |raw log2FC| >= threshold`, with up meaning log2FC > 0
  (higher in the numerator). Gene IDs are mapped to Entrez with
  `AnnotationDbi::select` (one-to-many: every Entrez target is kept; targets
  are deduplicated). A foreground with fewer than `minimum_mapped_foreground`
  (5) mapped genes is `NOT_APPLICABLE`.
- **ORA universe:** genes in `all_genes.tsv` with a **finite p-value**
  (`ora_helpers.R`), mapped to Entrez. Genes whose padj is NA after
  independent filtering stayed in the universe. The universe passed to
  `enrichGO`/`enrichKEGG` was intersected with the gene-set database inside
  clusterProfiler, but nf-rna did not record that size.
- **GSEA ranking metric:** the DESeq2 Wald statistic `stat` from
  `all_genes.tsv`. Genes with a non-finite stat are excluded. Gene IDs are
  mapped to Entrez; a duplicate Entrez keeps the largest |stat| (ties go to
  the lexically first `gene_id`). The list is sorted by decreasing stat, then
  Entrez ID.
- **KEGG source:** clusterProfiler downloads KEGG live from KEGG REST on every
  call (`provider: online_kegg_rest_via_clusterprofiler`). R first probes
  `https://rest.kegg.jp/list/pathway/<org>`, and the service runs a network
  preflight before a run. No snapshot is kept and the KEGG release is not
  recorded, so results can change when KEGG changes.
- **GO source:** the OrgDb (`org.Hs.eg.db` / `org.Mm.eg.db`) GO annotations
  via clusterProfiler. The OrgDb and clusterProfiler versions are recorded;
  GO.db is not.
- **Member tables:** `GO_<ont>_gene_members.tsv` (ORA, significant terms, keyed
  by Entrez with `original_gene_id`), `gsea_core_members.tsv` and
  `core_members.tsv` (GSEA leading edge only, with `stat`), and
  `kegg_gene_members.tsv`. None of them has log2FoldChange, padj or the
  pathway's p.adjust.
- **Figures:** one dotplot per result set with significant terms (top 15; PNG
  at 150 dpi, TIFF, PDF). There are no GSEA running-score plots and no GO
  redundancy reduction.
- **Python entry points** (`rnaseq.go/gsea/kegg.execute_*`): the backend config
  lacks contrast numerator/denominator, so the ORA figure subtitle stops. The
  GO post-processing reads `summary["mapping"]`, which the R backend no longer
  writes, and the KEGG ORA report reads `pathway_count`/`unique_target_genes`,
  which can be missing. The Nextflow path (`workflow_support.enrichment_config`)
  is not affected.
- **Delivery and report:** every TSV and figure under `downstream/` is copied;
  JSON is copied only from an allow-list. The HTML report states the ORA
  predicate `pvalue <= cutoff AND qvalue <= cutoff`.

## Acceptance prerequisites (2026-10-30)

### 1. Phase 1 enrichment baseline

Phase 1 moved the ENSEMBL version-suffix rule into `gene_annotation.R`
(`nf_rna_annotation_key`), which `go_analysis.R`, `gsea_analysis.R` and
`kegg_analysis.R` now call. Before enrichment is changed, v1.4.1 and the
Phase 1 HEAD were run on `examples/gsea_two_group` (mouse, 111 Entrez genes).

| | v1.4.1 | Phase 1 HEAD |
|---|---|---|
| Source | `git worktree` at `v1.4.1` (`f69058e`) | `0a3382a` |
| R runtime | old-lock prefix `nf-rna-downstream-linux-64-18197d49…-61fa2616…` | prefix built from the Phase 1 lock (`0e721ff3…`) |
| L2 input to enrichment | v1.4.1 `all_genes.tsv` (9 columns) | Phase 1 `all_genes.tsv` (12 columns) |

Procedure: fresh copies of the example; `prepare_l2`/`execute_l2`; then the
four R backends (`go`, `gsea-go`, `kegg` ORA, `gsea-kegg`) with configs built
the way `workflow_support.enrichment_config` builds them for the Nextflow path
(output dirs `enrichment/{go,gsea_go,kegg,gsea_kegg}`, contrast
numerator/denominator included). Both versions ran within the same minute, so
the online KEGG resource was the same.

Result: **consistent**.

- All four backends: return code 0, status `SUCCESS` in both versions.
- Same file set (140 files: 96 TSV/JSON/TXT, 36 figures, 8 R logs).
- 84 of the 96 TSV/JSON/TXT files are byte-identical. All of the
  significant, all-terms, gene-member, ranked-list, mapping, summary and
  provenance files match. The other 12 differ only in absolute paths:
  `backend_config.json` (scratch paths), `figure_manifest.json` (font file path
  inside the R prefix) and `r_session_info.txt` (BLAS path inside the R
  prefix). They are identical once the paths are normalized.
- Figures: 34 of 36 PNG/PDF/TIFF files are byte-identical. The two GO
  `dotplot_CC.pdf` files (up, down) differ only in the PDF `/CreationDate`,
  which falls in a different second. Their decompressed drawing streams and
  300 dpi renderings (`pdftoppm`) are identical.

The three extra columns in the Phase 1 `all_genes.tsv` (`symbol`,
`description` and the group means/shrunken LFC) do not change any enrichment
input or output.

Note (pre-existing, not changed): the legacy Python entry points
`rnaseq.go.execute_go` and the GSEA/KEGG equivalents cannot complete on this
example in either version. Their backend config omits the contrast
numerator/denominator that the ORA figure subtitle requires, and after that
the Python post-processing reads a `summary["mapping"]` key the R backend no
longer writes. The supported Nextflow path (`enrichment_config`) supplies both,
so the baseline used that config shape.

### 2. Realistic test data: airway

- Generator: `scripts/dev/make_airway_testdata.R`, run with the separate env
  `/home/brian/nf-rna-work/conda/airway-testdata-env`, created with
  `CONDARC=/home/brian/nf-rna-work/conda/condarc.yml` and containing
  `bioconductor-airway`, `bioconductor-summarizedexperiment` and `r-jsonlite`.
  airway is not added to the nf-rna downstream runtime.
- Output (not in the repository): `/home/brian/nf-rna-work/testdata/airway/`
  holds `counts.csv`, `metadata.csv`, `contrasts.csv`, `project.yaml` and
  `testdata_provenance.json`.
- airway **1.30.0**, SummarizedExperiment 1.40.0, R 4.5.3. There are 63,677
  Ensembl genes (no version suffix) × 8 samples, written unchanged. The
  metadata maps `dex` → `condition` (trt/untrt) and keeps `cell` (4 cell lines,
  each with one trt and one untrt sample).
- Project: `paired_two_group`, `pair_id: cell`, formula `~ cell + condition`,
  contrast `trt_vs_untrt` (condition, numerator `trt`, denominator `untrt`),
  `padj < 0.05`, `|log2FC| ≥ 1`. The annotation block is ENSEMBL → org.Hs.eg.db,
  with `analysis.enrichment: []`.
- `rnaseq validate`: PASS WITH WARNINGS (30,208 all-zero genes; no genes
  filtered at validation), READY.

Phase 1 L2 (`execute_l2`, prefix built from the Phase 1 lock, about 41 s):
status SUCCESS. 22,369 genes remain after the existing all-zero/low-total
filter, with 1,018 significant (531 up, 487 down) and 4,337 NA padj from
independent filtering.

| Gene | gene_id | baseMean | log2FoldChange | padj | log2FoldChange_shrunk | significant |
|---|---|---|---|---|---|---|
| DUSP1 | ENSG00000120129 | 3409.0 | +2.948 | 1.8e-125 | +2.937 | yes |
| KLF15 | ENSG00000163884 | 561.1 | +4.459 | 2.2e-77 | +4.433 | yes |
| PER1 | ENSG00000179094 | 776.6 | +3.192 | 3.3e-81 | +3.177 | yes |
| ZBTB16 | ENSG00000109906 | 385.1 | +7.353 | 2.3e-40 | +7.273 | yes |

All four known dexamethasone-induced genes are significant with positive
log2FC, so the contrast direction is correct and no data were changed.

- Shrinkage: apeglm, coef `condition_trt_vs_untrt` on a releveled copy
  (reference `untrt`; the original fit had `trt` as reference). The copy passed
  verification, status `success`, and no NA shrunken values.
- Annotation: org.Hs.eg.db 3.22.0. 17,209/22,369 (76.9%) genes have a
  symbol and description; 242 gene_ids map to more than one symbol (first
  taken).
- Figures (volcano, MA plot, heatmap) were inspected visually. The top-10
  volcano labels include DUSP1, and the heatmap rows include KLF15, PER1 and
  DUSP1.

All Phase 2 acceptance uses these airway results
(`/home/brian/nf-rna-work/testdata/airway/downstream/l2/`).

## Changes

### Modified and new files

| File | Change |
|---|---|
| `src/rnaseq/r/enrichment_common.R` (new) | v1.5 constants (p.adjust ≤ 0.05 BH, nPermSimple 10000, top 5 running-score sets, simplify 0.7, dotplot top 15, figure formats); the `significant` column; ORA universe record; long member tables; running-score selection; GO simplify |
| `src/rnaseq/r/kegg_snapshot.R` (new) | KEGG snapshot download (KEGGREST), manifest, SHA-256 verification, latest/pinned selection, lock |
| `src/rnaseq/r/ora_helpers.R` | ORA universe source = padj non-NA |
| `src/rnaseq/r/gsea_term_filtering.R` | significance = p.adjust ≤ 0.05; `significant` column |
| `src/rnaseq/r/go_analysis.R` | universe ∩ GO ontology genes, explicit min/max GS size, members, simplify, 300 dpi, metadata |
| `src/rnaseq/r/gsea_analysis.R` | `nPermSimple`, members, simplify, running-score plots, 300 dpi, metadata |
| `src/rnaseq/r/kegg_analysis.R` | snapshot + `enricher()`/`GSEA()` (no live KEGG); seed; members; running-score plots; metadata |
| `src/rnaseq/r/figures.R` | `formats` passed through the dotplots; `gglist` drawing; `nf_rna_gsea_running_score` |
| `src/rnaseq/r/provenance.R` | `snapshot_root` omitted from the normalized config (machine path) |
| `src/rnaseq/enrichment_config.py` (new) | one backend-config builder for Nextflow and the Python entry points; KEGG cache root |
| `src/rnaseq/workflow_support.py` | uses the builder; HTML report text (universe, v1.5 rule, snapshot, methods incl. v1.4 rule) |
| `src/rnaseq/go.py`, `gsea.py`, `kegg.py` | bug fix: shared config (contrast levels), post-processing reads existing fields only; reports/provenance record v1.5 semantics and gene-set source |
| `src/rnaseq/models.py` | `annotation.enrichment.kegg.snapshot`; `resource_provider` default `kegg_rest_snapshot` (old value accepted) |
| `src/rnaseq/service.py` | preflight probes KEGG REST only without a snapshot; delivers `kegg_snapshot_manifest.json` |
| `src/rnaseq/downstream.py`, `downstream_runtime.py` | runtime probes: enrichplot, KEGGREST, GOSemSim, GO.db, openssl (already locked) |
| `src/rnaseq/__init__.py`, `CITATION.cff`, `CHANGELOG.md` | development version 1.5.0, `1.5.0 — Unreleased` entry (Phase 1 + 2, semantics change) |
| `docs/scientific_contract.md`, `docs/figure-standard.md`, `docs/dev/DECISIONS.md`, `.gitignore` | contract, figure notes, decisions, KEGG data guard |
| `tools/figure_review/kegg_offline_profile.R`, `generate_runs.py` | offline fixture: synthetic snapshot, `enricher()`/`GSEA()` replaced |
| `scripts/dev/make_airway_testdata.R` (new) | airway test data (prerequisite) |
| tests | updated `test_go`, `test_gsea`, `test_kegg`, `test_enrichment_mapping_qc`, `test_figure_dotplots`, `test_workflow_support`; new: Python entry point uses the Nextflow config (`test_go`), snapshot latest/pinned/tamper/missing (`test_kegg`) |

No lock, environment YAML or SHA256SUMS change in Phase 2.

### New enrichment contract (summary; full text in `docs/scientific_contract.md`)

- **Significance:** significant = `p.adjust` (BH) ≤ 0.05 for GO ORA, KEGG ORA, GO GSEA and KEGG GSEA. clusterProfiler is called with `pvalueCutoff = 1, qvalueCutoff = 1`, so the all-terms tables keep the non-significant terms and add a `significant` column; `pvalue`, `p.adjust` and `qvalue` are kept. Every summary and provenance records `significance_semantics: "v1.5"`; the v1.4 rules are documented, and old results are not rewritten.
- **ORA universe:** padj non-NA genes → Entrez → ∩ gene-set database genes (GO ontology / KEGG snapshot). min/max GS size 10/500 (GO now explicit; KEGG configured).
- **ORA input:** `significant.tsv` / `upregulated.tsv` / `downregulated.tsv` from L2, unchanged.
- **GSEA ranking:** DESeq2 Wald `stat`, unchanged and never the shrunken LFC. NA excluded; for a duplicate Entrez the largest |stat| is kept. Seed `set.seed(1)` + `seed = TRUE`, `nPermSimple = 10000`, `eps = 0`; all recorded.
- **Gene sets:** GO from the OrgDb (OrgDb, GO.db and clusterProfiler versions recorded; no snapshot). KEGG from a SHA-256-verified local snapshot.

### Outputs actually produced on airway (`trt_vs_untrt`)

Per GO ORA foreground (`significant`, `up`, `down`) × ontology (BP, MF, CC):
- `GO_<ont>_all_terms.tsv` (with `significant`), `GO_<ont>.tsv`, `GO_<ont>_simplified.tsv`
- `GO_<ont>_members.tsv` (new long format), `GO_<ont>_gene_members.tsv` (unchanged)
- `dotplot_<ont>` and `dotplot_<ont>_simplified`, each as `.png`, `_300dpi.png`, `.tiff`, `.pdf`, with `_source_data.tsv`

Per KEGG ORA foreground:
- `kegg_all_terms.tsv`, `kegg.tsv`, `kegg_members.tsv`, `kegg_gene_members.tsv`
- `dotplot.{png,_300dpi.png,tiff,pdf}`
- `kegg_snapshot_manifest.json` at the backend root

Per GO GSEA ontology:
- `all_terms.tsv`, `significant.tsv`, `positive_enrichment.tsv`, `negative_enrichment.tsv`, `significant_simplified.tsv`
- `members.tsv`, `gsea_core_members.tsv`
- `dotplot` and `dotplot_simplified` (4 formats each)
- `running_score/{positive,negative}_<n>_<GO ID>.{png,_300dpi.png,pdf}`: 10 for BP, 10 for MF, 8 for CC (CC has only 3 significant negative sets)

KEGG GSEA:
- the same tables without simplify
- 10 running-score plots

Counts: 392 files in `enrichment/`, 118 of them TSV tables; 52 tables carry `significant` + `p.adjust`.

| Result set | evaluated | significant (p.adjust ≤ 0.05) | simplified |
|---|---|---|---|
| GO ORA significant BP / MF / CC | 4531 / 808 / 445 | 558 / 53 / 35 | 214 / 24 / 20 |
| GO ORA up BP / MF / CC | 3715 / 665 / 390 | 98 / 20 / 9 | 59 / 12 / 6 |
| GO ORA down BP / MF / CC | 3721 / 588 / 347 | 361 / 28 / 34 | 167 / 12 / 23 |
| KEGG ORA significant / up / down | 322 / 275 / 291 | 16 / 5 / 26 | — |
| GO GSEA BP / MF / CC | 5684 / 1130 / 718 | 245 / 17 / 42 | 145 / 12 / 21 |
| KEGG GSEA | 350 | 29 (23 positive, 6 negative) | — |

### KEGG snapshot scheme

`~/.cache/nf-rna/genesets/kegg/<org>/<snapshot_id>/` (honours `XDG_CACHE_HOME`; outside the repository; `.gitignore` also guards `genesets/`, `TERM2GENE.tsv`, `TERM2NAME.tsv`):
- `TERM2GENE.tsv` comes from `https://rest.kegg.jp/link/pathway/<org>` and `TERM2NAME.tsv` from `https://rest.kegg.jp/list/pathway/<org>`, both via KEGGREST.
- `manifest.json` holds the URLs, UTC retrieval time, release information, raw `/info/<org>` and `/info/pathway` text, SHA-256 and size of both files, and the pathway/gene counts.
- **Snapshot ID:** `<org>_kegg<pathway last update>_<UTC stamp>`, e.g. `hsa_kegg20261028_20261030T061908Z` (372 pathways, 9451 genes, KEGG pathway database "Last update 2026/10/28").
- **Selection:** `annotation.enrichment.kegg.snapshot` if set; otherwise the latest existing snapshot; a new download only when none exists, under a lock so parallel KEGG ORA/GSEA tasks share one. On airway, ORA created the snapshot and GSEA used it as the latest.
- **Verification:** SHA-256 values are verified before every use and recorded in the summary and provenance.
- **Failure:** no snapshot means `NETWORK_UNAVAILABLE`, `reproducible: false` ("NOT_FULLY_REPRODUCIBLE"), and no result. There is no fallback to live KEGG or to `use_internal_data`.
- The download succeeded here; a mouse snapshot (`mmu`) was also created by the cross-check below.

## Acceptance (airway, Phase 1 L2 results)

Runtime: R prefix built from the Phase 1 lock (`0e721ff3…`), Python `/home/brian/miniforge3/envs/nf-rna/bin/python`. Every analysis ran through the fixed Python entry points (`rnaseq.go.execute_go`, `rnaseq.gsea.execute_gsea`, `rnaseq.kegg.execute_kegg(mode = "ora" | "gsea")`) with scratch driver `run_enrich_py.py`:

```
python run_enrich_py.py /home/brian/nf-rna-work/testdata/airway go,gsea,kegg,gsea-kegg      # run 1
python run_enrich_py.py <copy of airway with the same downstream/l2> go,gsea,kegg,gsea-kegg  # run 2
go SUCCESS 112s / gsea SUCCESS 97s / kegg SUCCESS 24s / gsea-kegg SUCCESS 107s           (run 1)
go SUCCESS 105s / gsea SUCCESS 91s / kegg SUCCESS 24s / gsea-kegg SUCCESS 108s           (run 2)
```

1. **`significant` ≡ `p.adjust ≤ 0.05`.** Across all 52 result tables with both columns, 26,604 rows were checked with 0 mismatches. The full tables keep 22,404 non-significant rows. Every significant/positive/negative/simplified table contains only `TRUE`.
2. **Universe, mapping, DEGs.**
   - 18,032 padj-non-NA genes (of 22,369; 4,337 NA); 15,049 mapped to Entrez, a mapping rate of **83.46%** (recorded; lower than 100% as expected for Ensembl → Entrez); 15,218 unique Entrez IDs.
   - DEGs: 1,018 significant, 531 up and 487 down; mapped foregrounds 949 / 490 / 459.
   - **Final universe:** GO BP 12,385, MF 12,848, CC 12,992; KEGG 6,057.
   - For every foreground and database, the universe equals clusterProfiler's BgRatio N.
   - An independent R recomputation from `all_genes.tsv` with AnnotationDbi (no nf-rna code) gives the same 18,032 / 15,049 (0.8346) / 15,218 / 12,385 / 12,848 / 12,992 / 6,057.
   - GSEA ranking: 22,369 finite stats, 17,209 mapped (76.93%), 17,425 ranked Entrez; 6,542 of them are in the KEGG snapshot.
3. **Snapshot SHA-256.** Recomputed with Python `hashlib`: TERM2GENE `d228cd3e…2a6c` and TERM2NAME `ad28c97b…bd6c` match `manifest.json` and the provenance. ORA and GSEA used the same snapshot ID.
4. **Members.** 16 member tables, 0 unmapped symbols (0.00%), 0 missing gene_ids. Largest: GO GSEA BP with 35,438 rows; KEGG ORA significant with 321 rows.
5. **Running-score plots.** 38 drawn, each as PNG 2125×1488 at 300 dpi plus a one-page PDF. The PDF titles (`pdftotext`) give NES and adjusted P matching `significant.tsv`, e.g. hsa04140 Autophagy – animal: NES = 2.01, adjusted P = 1.72e-05. Inspected visually.
6. **Two runs, same snapshot and input.** Same 392-file set; all 118 TSV tables are byte-identical (sha256 lists identical); all four backend summary JSONs are identical.
7. **DEG outputs unchanged.** The 13 L2 TSV/JSON/CSV files (`contrasts/trt_vs_untrt/*.tsv`, `contrast_metadata.json`, `normalized_counts.tsv`, …) have the same SHA-256 before and after Phase 2 (`sha256sum -c`: 13 OK).

Cross-check: `examples/gsea_two_group` (mouse, Entrez): L2 and all four analyses `SUCCESS` through the Python entry points, with a new `mmu` snapshot.

Tests (R prefix from the Phase 1 lock):
- enrichment set (`test_go`, `test_gsea`, `test_kegg`, `test_enrichment_mapping_qc`, `test_figure_dotplots`): 57 passed, plus the two new tests.
- broad set (`test_workflow_support`, `test_enrichment_defaults`, `test_config`, `test_service`, `test_cli`, `test_figure_settings`, `test_downstream*`, `test_doctor_native`, `test_installed_runtime`, `test_runtime_dispatch`, `test_workflow`, `test_user_workflow`, `test_planner`, `test_run_status`, design tests, `test_figure_review_tools`, `test_figures`): 609 passed, 1 skipped, 1 failed. The failure was a stale provider assertion, fixed and rerun.
- release metadata, packaging, CLI, workflow support after the version bump: 77 passed.
- The full suite was not run.

## Known limitations

- **Nextflow path not run end to end.** `rnaseq run` was not executed: it needs a commit and provisions the runtime. The Nextflow enrichment config is built by the same `backend_config` as the Python entry points, and the unit tests cover it.
- **KEGG release information.** KEGG REST no longer publishes a "Release" number; the snapshot records the pathway database "Last update" date and the raw `/info` texts instead.
- **GO is not snapshotted.** Results are reproducible only for the recorded OrgDb, GO.db and clusterProfiler versions, which the lock pins.
- **Tied ranking statistics.** fgsea warns that 1.43% of the ranked statistics are tied; nf-rna orders ties by Entrez ID, and two runs were byte-identical.
- **Harmless warnings.** `gsea_analysis.R` prints about 30 "object not found" warnings from its existing `on.exit(rm(...))`, which predates this phase and was left unchanged.
- **Member tables** cover significant pathways only.
- **HTML report gaps.** `report.html` states the new rules but does not show the running-score plots, the simplified dotplots or the MA plot; they are delivered under `figures/`.
- **Optional item E (MSigDB Hallmark) was not done.** msigdbr is not in the lock and would bring its dependency `babelgene` (not "msigdbr only"), which needs a lock regeneration.
- **Development version.** The package version is now the development version 1.5.0. Existing provisioned runtime prefixes (old lock) do not contain apeglm/ashr, so the first `rnaseq run` provisions the new lock.
