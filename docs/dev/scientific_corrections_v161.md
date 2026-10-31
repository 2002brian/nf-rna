# Scientific corrections for v1.6.1: D-1, D-2, M-1

Branch `fix/v1.6.1-scientific` was created from the `v1.6.0` tag (ff4fe5c).

| Commit | Content |
|---|---|
| `3cf0386` | The corrections, the regression tests and the CHANGELOG entry |
| Separate later commit | This document and the evidence in [`validation/corrections_v161/`](../../validation/corrections_v161/) |

The findings come from the Phase 2 independent validation of v1.6.0, recorded on branch `validation/v1.6.0`, commit `b886f11`:

- `docs/validation/report.md`
- `validation/results/phase2/triage.md`

That evidence, and the Phase 2 FAIL verdict, are unchanged. Nothing on `validation/v1.6.0`, the `v1.6.0` tag or any existing v1.6.0 run was modified.

**All three corrections can change results.** v1.6.0 results are not rewritten. A project must be rerun to obtain corrected results, and the rerun gets a new run directory.

## Summary

| ID | v1.6.0 behaviour | Corrected behaviour | Outputs that change | Independent check (Phase 2 reference) |
|---|---|---|---|---|
| D-1 | Leading edge from DOSE 4.4.0 `leading_edge()`, which is wrong for some ES < 0 sets | nf-rna recomputes the leading edge by the standard running-sum definition | `core_enrichment`, `leading_edge` tag string, core-member tables, `members.tsv` `leading_edge` flag (ES < 0 sets only) | 33,544 / 33,544 sets equal fgsea `leadingEdge` (v1.6.0: 29,990) |
| D-2 | `seed = TRUE`: DOSE re-seeds with a fixed state, so the configured seed is ignored | `set.seed(seed)` then `seed = FALSE`: the configured seed drives fgsea | GSEA NES, p-value, p.adjust and borderline significance (Monte Carlo level); ES unchanged | Significant sets identical to independent fgsea with the same seed in 20/20 collections (seeds 1 and 42) |
| M-1 | DESeq2 `results()` used `alpha = 0.1` while DEGs are called at `padj < 0.05` | `alpha = thresholds.padj` | Filter threshold, `padj` / `padj = NA` pattern, DEG sets, ORA universe, foreground and significant terms | Identical to the reference's alpha = 0.05 analysis in 6/6 contrasts and 48/48 ORA analyses |

## D-1: GSEA leading edge

### Locked DOSE behaviour

The locked runtime is `0e721ff3…`, with DOSE 4.4.0, clusterProfiler 4.18.4 and fgsea 1.36.2. In it, `DOSE:::leading_edge()`:

1. keeps only the running-score rows at hits (`position == 1`), after the hit's step;
2. for ES < 0, takes `which.min()` over those rows and reports the hits from that row to the end.

The running-sum minimum of an ES < 0 set lies on a miss just *before* a hit. When that hit carries a large weight, an earlier hit can have a lower post-step value, so DOSE reports a boundary that is too early. Its leading edge then includes extra genes.

Which fields are affected:

- `tags`, and therefore `signal`, inherit the same error.
- `rank` and `list` are computed by DOSE on the full running sum and were correct.
- The positive-ES branch is correct, because the maximum lies on a hit.

The counterexample is in `tests/test_scientific_corrections_v161.py`. It is a 14-gene list and a 5-gene set:

| Source | Leading edge |
|---|---|
| DOSE | g9/g11/g14 |
| Standard definition and fgsea | g14 |

### Correction

`src/rnaseq/r/gsea_leading_edge.R` recomputes the leading edge for both GO GSEA (`gseGO`) and KEGG GSEA (`GSEA`). It runs on the `gseaResult` before any table, member list or figure is written.

**Positions.** These are the positions in the `gseaResult` `geneList`, the exact vector passed to fgsea: DESeq2 `stat` descending, with ties broken by ascending Entrez ID. Gene identifiers are the Entrez IDs of that list, mapped back to the source `gene_id` exactly as in v1.6.0. Members of a gene set that are absent from the list are ignored, as in fgsea.

**Running sum.** Weight exponent 1 (fgsea `gseaParam` = 1); other exponents are rejected. For hits h₁ < … < h_k:

- after hit j: top_j = Σ_{i≤j} |r_i| / N_R − (h_j − j)/(N − k)
- just before hit j: bottom_j = top_j − |r_j| / N_R

**Leading edge.**

- ES > 0: hits 1..i, where top_i is the maximum.
- ES < 0: hits i..k, where bottom_i is the minimum.
- A tied extremum resolves to its first occurrence along the list. This matches fgsea and GSEA desktop.

**Derived fields.** `rank` keeps DOSE's definition: for ES < 0 it is counted from the minimum to the end. `tags`, `list` and `signal` are recomputed with DOSE's formulas from the corrected boundary.

**Guards.**

- The recomputed extremum must equal the fgsea ES within 1e-9.
- The number of members in the list must equal `setSize`.

Otherwise the analysis fails rather than writing a silent leading edge.

**Provenance.**

- `gsea.leading_edge` states the definition.
- Each GO ontology and each KEGG contrast summary records `leading_edge` (an audit): the number of terms whose `core_enrichment` differs from DOSE's, split by ES sign, and the maximum ES difference.

ES, NES, p-values, p.adjust and significance are not touched by this fix.

### Evidence: correction alone, applied to the exact v1.6.0 GSEA

Script: `validation/corrections_v161/scripts/d1_isolation.R`. Result: `results/d1_isolation.tsv`.

The GSEA of all 16 collections of VAL-D1 and VAL-D2 was recomputed exactly as v1.6.0 did it (seed 1, `seed = TRUE`) from each run's own `ranked_gene_list.tsv`. The correction was then applied.

| Check | Result |
|---|---|
| Recomputation reproduces v1.6.0 `all_terms.tsv` (ES, NES, p, p.adjust, q, core, rank, tag) | 16/16 collections |
| ES, NES, pvalue, p.adjust, qvalue, setSize unchanged by the correction (in-memory `identical`) | 16/16 |
| `core_enrichment` changed | 3,554 of 33,544 terms. All 3,554 are among the 14,908 terms with ES < 0, and **0** have ES > 0 |
| Changed among significant terms (p.adjust ≤ 0.05) | 162 |
| `rank` changed | 0 |
| `leading_edge` tag string changed | 3,523. The other 31 differ only below the tag's rounding to whole percentages |
| Recomputed extremum − fgsea ES | 0 in every term |
| Leading edge equal to the independent reference fgsea `leadingEdge` (as sets) | v1.6.0: 29,990 / 33,544. Corrected: **33,544 / 33,544** |

The 3,554 and 162 are exactly the Phase 2 D-1 figures. No unexpected difference was found.

## D-2: GSEA seed

### Locked DOSE behaviour

In DOSE 4.4.0, `GSEA_fgsea()` runs `if (seed) set.seed(.Random.seed)` immediately before `fgsea()`. `set.seed()` uses only the first element of `.Random.seed`, which is the RNG-kind code (10403 for Mersenne-Twister/Inversion/Rejection). So with `seed = TRUE`:

- every call is re-seeded to the same state, `set.seed(10403)`;
- the `set.seed(seed)` that nf-rna ran just before has no effect.

This is confirmed in the locked runtime by `test_configured_seed_controls_fgsea_and_is_reproducible`: seeds 7 and 8 give identical results, and both equal `set.seed(10403); fgsea(...)`.

### Correction

GO and KEGG GSEA keep `set.seed(as.integer(gsea_cfg$seed))` immediately before `gseGO`/`GSEA`, and now pass `seed = FALSE`. DOSE then leaves the RNG state alone.

The following are unchanged:

- the fgsea call: multilevel, `eps = 0`, `nPermSimple = 10000`, `nproc = 0`, exponent 1;
- the gene sets, the ranking and minGSSize/maxGSSize;
- the BH adjustment;
- the significance rule `p.adjust ≤ 0.05`.

Provenance `gsea.seed_call` describes the call, and `gsea.rng_kind` records the RNG kind. The recorded `seed` is now the seed actually used.

### Evidence

Script: `scripts/compare_v161.py`. Results:

- `results/gsea_fixed_vs_reference_and_v160.tsv`
- `results/reproducibility.tsv`
- `results/seed_1_vs_42.tsv`

**Agreement with the independent reference.** The corrected nf-rna (FIX-D1 and FIX-D2 with the project seed 1, FIX-D1-S42 with seed 42) was compared with the Phase 2 reference's direct `set.seed(seed); fgsea(...)` with the same seed. The reference has its own gene sets and ranking, a serial backend, and no clusterProfiler or DOSE.

| Check | Result |
|---|---|
| Significant sets (p.adjust ≤ 0.05) identical | **20/20 collections**: 16 for seed 1, 4 for seed 42 |
| NES / p-value / p.adjust, maximum relative difference | ≤ 6.5e-15 in 18/20 collections. ≤ 2.1e-4 in D2 LP_vs_ML GO_BP (571 sets) and GO_MF (1 set); see open finding U-1 |
| ES | ≤ 6e-15 in all 20 |
| Leading edge equal to the reference | 100% in all 20 collections, both ES signs |

**Reproducibility.**

| Check | Result |
|---|---|
| Same seed: FIX-D1 vs FIX-D1-REP | 125/125 TSV tables byte-identical |
| Different seed: FIX-D1 (1) vs FIX-D1-S42 (42) | All 87 non-GSEA tables identical. GSEA: ES identical; NES differs in every set, and p-value in all but 1–2 sets per collection. Significance differs for 28 (BP), 3 (MF), 5 (CC) and 1 (KEGG) sets |

The test suite does not require every pathway to differ between seeds.

**Change from v1.6.0.** These figures are for the 16 collections of D1 and D2 and include D-1. The ranking is unchanged because DESeq2 `stat` does not depend on alpha.

| Quantity | Change |
|---|---|
| ES | Identical in all 16 collections |
| NES | Changes in every set. Median \|ΔNES\| per collection 0.002–0.007; maximum 0.027 |
| p-value | Changes in nearly every set. Median \|log10 ratio\| per collection 0.007–0.015; maximum 0.96 (a factor of about 9, for a multilevel p-value near 1e-11: D2 LP_vs_ML KEGG mmu00190) |
| Significant sets | 3,016 → 2,869 in total: 69 gained, 216 lost, all near the threshold |

The asymmetry in significant sets is expected. Phase 2 (triage T-4) showed that the fixed v1.6.0 draw (state 10403) lies in the upper tail of the Monte Carlo distribution for GO BP. For example, D2 LP_vs_Basal had 714 significant sets in v1.6.0, against 656–699 over eight independent seeds; with seed 1 it now has 660.

Per-collection figures are in `gsea_fixed_vs_reference_and_v160.tsv`:

| Columns | Meaning |
|---|---|
| `significant_gained` / `significant_lost` | Sets entering or leaving significance relative to v1.6.0 |
| `nes_abs_change_*` | Absolute NES change |
| `pvalue_abs_log10_ratio_*` | Absolute log10 ratio of the p-values |

## M-1: DESeq2 independent-filtering alpha

### Basis

**DESeq2 documentation** (locked DESeq2 1.50.2):

- `?results`, `alpha`: "the significance cutoff used for optimizing the independent filtering (by default 0.1). If the adjusted p-value cutoff (FDR) will be a value other than 0.1, alpha should be set to that value."
- The vignette, section *p-values and adjusted p-values*, says the same and shows `results(dds, alpha = 0.05)`.

**nf-rna configuration.** The DEG FDR threshold is the project setting `thresholds.padj` (default 0.05, configurable). The DEG rule is `padj < thresholds.padj` and `|log2FoldChange| ≥ thresholds.abs_log2fc`.

### Correction

`l2_analysis.R` calls `results(..., independentFiltering = TRUE, alpha = cfg$thresholds$padj)`. Independent filtering itself is unchanged: the filter statistic is baseMean, the threshold is chosen to maximise rejections at `padj < alpha`, and BH is used.

What is recorded:

- `backend_summary.json`: `independent_filtering_alpha` and `independent_filtering_threshold`, the latter from `metadata(res)$filterThreshold`.
- `contrast_metadata.json`, `deg_definition.independent_filtering`: alpha, filter statistic, threshold and rule.
- The L2 run summary: `independent_filtering_alpha`.

`thresholds.padj` must now satisfy 0 < padj < 1, because DESeq2 rejects alpha outside this range. A project with `padj: 0` or `padj: 1` now fails validation instead of failing inside R. Neither value gave a meaningful DEG call before.

pvalue, stat, log2FoldChange, lfcSE and the shrunken LFC do not depend on alpha. All are byte-identical to v1.6.0 in every contrast.

### Evidence: corrected nf-rna vs v1.6.0 and vs the reference's prespecified alpha = 0.05 sensitivity analysis

Results: `results/m1_deseq2_alpha.tsv` and `results/m1_ora.tsv`.

| Dataset / contrast | baseMean filter threshold (0.1 → 0.05) | padj NA (v1.6.0 → fixed) | DEGs (v1.6.0 → fixed: +gained / −lost) | Equals reference alpha 0.05 |
|---|---|---|---|---|
| D1 airway (paired) | 4.50 → 6.14 | 4,337 → 5,204 | 1,018 → 999 (+6/−25) | NA pattern, padj (≤5e-15), DEG set |
| D1a (two-group) | 3.86 → 4.50 | 3,954 → 4,387 | 872 → 878 (+7/−1) | yes |
| D1b (batch) | 4.50 → 6.14 | 4,337 → 5,204 | 1,018 → 999 (+6/−25) | yes |
| D2 LP_vs_Basal | 1.70 → 1.70 | 743 → 743 | 6,173 → 6,173 (0/0) | yes |
| D2 ML_vs_Basal | 1.70 → 1.29 | 743 → 372 | 6,319 → 6,314 (0/−5) | yes |
| D2 LP_vs_ML | 2.94 → 3.79 | 1,486 → 1,858 | 2,931 → 2,936 (+6/−1) | yes |

D1b equals D1 because both have the same counts and the same effective model. Note that a lower alpha does not always filter more: in D2 ML_vs_Basal, alpha 0.05 chose a lower baseMean threshold.

The threshold recorded in `backend_summary.json` is rounded to 4 decimals by `jsonlite`. The relative difference from the reference (≤ 1.7e-5) is that rounding; see U-2.

**ORA.** There are 48 analyses: D1 and D2 × GO BP/MF/CC and KEGG × up/down/all DEGs.

- The universe and foreground sizes and the significant term sets equal the reference's alpha 0.05 analysis in 48/48 analyses, with p.adjust ≤ 5e-15.
- Changes from v1.6.0:
  - D1 universe shrinks by about 2%, for example GO BP 12,385 → 12,138.
  - D2 LP_vs_Basal is unchanged.
  - D2 ML_vs_Basal universe grows, for example GO BP 15,974 → 16,248.
  - D2 LP_vs_ML universe shrinks.

Significant ORA terms, v1.6.0 → fixed:

| Analysis | Before → after (+gained/−lost) |
|---|---|
| D1 GO BP, all DEGs | 558 → 579 (+34/−13) |
| D1 GO BP, up | 98 → 98 (+4/−4) |
| D1 GO BP, down | 361 → 364 (+18/−15) |
| D1 GO MF, all | 53 → 60 |
| D1 GO CC, all | 35 → 32 |
| D1 KEGG, all | 16 → 15 |
| D2 LP_vs_Basal | No change in any collection |
| D2 LP_vs_ML GO BP, all | 794 → 823 (+34/−5) |
| D2 ML_vs_Basal GO BP, all | 879 → 890 (+39/−28) |
| Largest change in any analysis | 39 terms |

The gained and lost term IDs are listed per analysis in `m1_ora.tsv`.

GSEA is not affected by M-1, because the ranking uses `stat`.

## Reproducibility and migration

**Existing results.** v1.6.0 runs stay valid records of what v1.6.0 computed, and nothing rewrites them. To obtain corrected results, rerun the project (`rnaseq plan`, then `rnaseq run`). A corrected run can be recognised by:

- `source_revision` (3cf0386 or later);
- the provenance fields `gsea.seed_call` (mentions `seed = FALSE`), `gsea.leading_edge` and `independent_filtering_alpha`.

v1.6.0 runs lack all three fields.

**Comparing across versions.**

- **DE:** pvalue, stat and LFCs compare directly. padj, NA patterns and DEGs differ as listed above.
- **ORA:** results differ wherever the DEG set or universe changed.
- **GSEA:** ES compares directly. NES, p and padj differ at the Monte Carlo level, and v1.6.0 leading edges of ES < 0 sets should not be used.

**Reproducing v1.6.0 GSEA numbers.** They are a single draw, equivalent to `set.seed(10403); fgsea(...)` with the same arguments.

**Seeds.** With the correction, the seed controls the GSEA draw. Different seeds give different but equally valid p-values. Report the seed, which provenance records, and expect borderline sets near p.adjust 0.05 to move between seeds.

**Projects.** No project-file change is needed, unless a project sets `thresholds.padj` to 0 or 1.

**Version.** The package version on this branch is `1.6.1`, a development version: CHANGELOG `1.6.1 — Unreleased`, and CITATION.cff without a release date. No tag was created.

The corrected validation runs were made from commit `3cf0386`, before the version change, so they still report `1.6.0`. `source_revision` identifies them.

The release date is set when the release is finalised.

## Tests

`tests/test_scientific_corrections_v161.py` (9 tests) covers:

- the DOSE counterexample;
- agreement with `fgsea::calcGseaStat(returnLeadingEdge = TRUE)` over 3,000 random sets: both ES signs, ties in the running-sum extremum, members absent from the list, and rank-ordered output;
- that the corrected `gseaResult` changes only the leading-edge fields, and that misaligned ES or setSize fails;
- the seed: reproducible, different seeds reach fgsea, equality with direct `set.seed(seed); fgsea()`, and the characterisation of DOSE `seed = TRUE`;
- static checks of both backends;
- the DESeq2 alpha semantics;
- the padj bounds.

`tests/test_l2.py` additionally checks the recorded alpha and filter threshold of a real L2 run. One existing lifecycle assertion in `tests/test_gsea.py` was updated for the new `corrected` object.

**Full-suite result.** The full suite on `3cf0386` (locked runtime `0e721ff3…`) gave 1,087 passed, 2 skipped and 4 failed. Both causes were corrected in the follow-up commit:

- **Three `test_figure_dotplots.py` failures.** The offline KEGG figure fixture (`tools/figure_review/kegg_offline_profile.R`) mocked `GSEA()` with an ES unrelated to its gene sets, which the new ES guard rejects. The mock now reports the real running-sum ES of its synthetic sets.
- **`test_release_metadata.py`.** The repository's development-version convention needed the version set to `1.6.1 — Unreleased`.

`test_figure_dotplots.py`, `test_figure_review_tools.py`, `test_release_metadata.py`, `test_packaging_contract.py` and `test_cli.py` then pass: 70 tests.

## Open findings

- **U-1. Monte Carlo-level GSEA differences in 2/20 collections.** For D2 LP_vs_ML, GO BP (571 sets, NES and p up to 2.1e-4 relative) and GO MF (1 set) differ from the reference fgsea run with the same seed, although ES and the significant sets are identical.
  - Phase 2's fixed-state check (seed 10403) showed the same pattern in two *other* collections (D1 GO MF, D2 ML_vs_Basal GO BP), so it moves with the seed.
  - The likely cause is a difference in how the random stream is consumed between `gseGO`'s fgsea call and the reference's call. Candidates are the order of the pathway list or the parallel and serial backends. This was not isolated.
  - It does not affect the correctness of D-2. The configured seed is used and reproducible, and the statistics are a valid fgsea draw.
- **U-2. Precision of recorded values.** `backend_summary.json` rounds `independent_filtering_threshold` to 4 decimals (jsonlite's default). This is the same specification gap as Phase 2 T-2 for mapping rates. The computation uses full precision.
- **U-3. Upstream reporting.** The DOSE `leading_edge()` ES < 0 defect and the `seed = TRUE` behaviour should be reported to DOSE upstream. Neither has been reported.
- **Not addressed here.**
  - The Phase 2 S-GS-6/S-GS-9 rule calibration (triage T-4) is a criteria question, not a software defect.
  - Phase 3 (upstream FASTQ) remains pending.
  - `rank` for ES < 0 keeps DOSE's convention (counted from the minimum, which is a miss) rather than any other convention.
