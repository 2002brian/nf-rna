# Scientific corrections for v1.6.1: D-1, D-2, M-1

Branch `fix/v1.6.1-scientific` was created from the `v1.6.0` tag (ff4fe5c).

| Commit | Content |
|---|---|
| `3cf0386` | The corrections, the regression tests and the CHANGELOG entry |
| `7d85046` | Version `1.6.1 — Unreleased`; consistent ES in the offline KEGG GSEA figure fixture |
| `1cbe8fb` | This document and the evidence in [`validation/corrections_v161/`](../../validation/corrections_v161/) |
| Later commit | The U-1 investigation, the full-suite result on HEAD and the final acceptance |

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
| NES / p-value / p.adjust, maximum relative difference | ≤ 6.5e-15 in 18/20 collections. ≤ 2.1e-4 in D2 LP_vs_ML GO_BP (571 sets) and GO_MF (1 set); explained in [U-1 investigation](#u-1-investigation) |
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

The corrected validation runs (FIX-*) were made from commit `3cf0386`, before the version change. Their provenance therefore reports `nf_rna_version` `1.6.0`, but `source_revision` `3cf038630ff56301a2d48e672a640121e70a0497` identifies the corrected code unambiguously. The runs have not been rewritten, and their provenance was not edited after the fact.

Checks on HEAD:

- `rnaseq --version` prints `1.6.1`.
- `CITATION.cff` has `version: "1.6.1"` and no `date-released`.
- `CHANGELOG.md` has `## 1.6.1 — Unreleased`.
- A new analysis records `nf_rna_version` from `rnaseq.__version__`, through `PIPELINE_VERSION` and the runtime identity, so it records `1.6.1`.
- A new analysis records `source_revision` from the git HEAD of a source checkout, or from the commit pip recorded for a git install. A checkout with uncommitted changes reports `+dirty`, and `rnaseq run` refuses it (`tests/test_source_identity.py`).

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

**Full suite on `1cbe8fb`** (locked runtime `0e721ff3…`, `NF_RNA_REQUIRE_R=1`, log in `/home/brian/nf-rna-validation/v161/fullsuite/pytest_1cbe8fb.log`): **1,091 passed, 2 skipped, 0 failed** in 26 min 43 s. All four earlier failures are resolved.

Both skips are environment-gated and were also skipped on `3cf0386`:

- `test_featurecounts_real_semantics.py:55` needs real SAMtools and featureCounts.
- `test_run_status.py:664` needs a real SJ project that is not present.

Commits after `1cbe8fb` change only this document and `validation/` evidence, so the result applies to HEAD's code and tests.

## U-1 investigation

**Question.** With seed 1, D2 LP_vs_ML GO BP (571 sets) and GO MF (1 set) differ from the reference fgsea in NES and p-value by up to 2.1e-4 relative. ES and the significant sets are identical. The investigation asked four things:

- Is the cause RNG consumption order, or RNG state shared across gene sets or contrasts?
- Can the difference be explained once the seed, RNG kind, ranking and gene sets are fixed?
- Is there an algorithmic difference?
- Is the difference within Monte Carlo error?

**Method.** The investigation used minimal reproductions in the locked runtime (fgsea 1.36.2) and the fgsea 1.36.2 source package from Bioconductor 3.22 (`fgsea_1.36.2.tar.gz`, SHA-256 `d6996ff5…5f5a`).

- Scripts: `validation/corrections_v161/u1/` (`run_u1.sh` runs all of them).
- Outputs: `u1/results/`.

**Acceptance criterion.** This was fixed before the comparison was evaluated, but after the discrepancy had been observed:

- A same-seed difference is acceptable as Monte Carlo-level when |log2(p_nf-rna / p_reference)| ≤ the reference's own `log2err` (fgsea's estimate of the standard deviation of log2 p) for every set.
- The significant sets must also be identical.
- Bitwise equality of independent runs is not required.

### Findings

**1. Not RNG consumption, shared state or backend.** The reference's gene sets and fgsea call were held fixed (seed 1, `SerialParam`, the reference parameters). Only the input statistic was varied:

| Input to fgsea | vs reference output | vs nf-rna output |
|---|---|---|
| A: reference ranked list (the DESeq2 `stat` at full double precision) | identical | 571 NES / 333 p differ (GO BP) |
| B: nf-rna ranked list (the same `stat`, as nf-rna reads it from `all_genes.tsv`, written with 15 significant digits) | 571 NES / 333 p differ (GO BP) | **identical** |
| A, with all collection sets passed unfiltered, as DOSE does | identical | as A |
| A, with `MulticoreParam(4)` instead of serial | identical | as A |

GO MF behaves the same way: 1 set differs.

- Each fgsea call seeds its own C++ engine from a single `sample.int()` draw. No RNG state is carried across gene sets, collections or contrasts.
- The backend does not change the result.
- Passing the unfiltered set list, which is what DOSE does, gives the same result.
- **Caveat.** Shuffling the *order* of the gene sets does give a different (equally valid) draw for every set. nf-rna and the reference present the sets in the same order, so order is not the cause here. It does mean that a same-seed result also depends on the order of the gene sets, for example after an annotation update.

**2. Mechanism.**

- **Input.** nf-rna builds the GSEA ranking from `all_genes.tsv`, which R's `write.table` writes with 15 significant digits. The reference used the in-memory statistic. Both lists have the same gene order and differ only in the last digits (about 1e-15 relative).
- **The flip.** That is enough to flip a floating-point comparison in exactly one permutation of the simple stage. The permutation is reached through the same random stream (C++ seed 66608964 from `set.seed(1)`).
  - **GO BP.** In permutation 1089 of 10,000, the random sets of size 10 and 11 have maximum positive and maximum negative running-sum deviations equal to at least 13 digits (±0.2006017, ±0.2006125).
    - In `fastGSEA.cpp`, `calcGseaStatCumulative()` sets the `std` ES to 0 when the two are exactly equal (`res == resDown`).
    - With input B they are exactly equal, so the permutation scores ES = 0 and is counted as both ≥ 0 and ≤ 0.
    - With input A the negative deviation is larger by a few ulp, so ES = −0.2006.
    - So for every set of size 10 or 11 (571 sets), `nGeZero` differs by 1 and the zero-side means differ by one permutation's contribution. NES changes by about 1/5,000, and p by 1/(1 + nGeZero) where p depends on that count.
    - Five sets whose observed ES lies between −0.2006 and 0 also shift their ES-side count by 1.
  - **GO MF.** In permutation 1079 the random ES of a size-10 set equals the observed ES of GO:0051880 (0.2247770495326) to within rounding. Under A it counts as ≥, under B as ≤. One p-value changes by one count, and NES is unchanged.
- **Same algorithm.** Both runs use the same algorithm, on the same random draws. The difference is a discontinuity of the permutation statistic at exact ties, exposed by a 1e-15 change in the input.
- **Phase 2 cases.** The two seed-10403 cases from Phase 2 (D1 GO MF, D2 ML_vs_Basal GO BP) fit the same mechanism with a different random stream. They were not re-examined.

**3. Size against Monte Carlo error.** From `u1/results/monte_carlo_error.txt`:

| Collection | Sets differing | max \|log2 p ratio\| / log2err | Significance changes | Seed 1 vs seed 42 in the reference, median \|log2 p ratio\| for the same sets |
|---|---|---|---|---|
| GO BP | 571 | 0.02 | 0 | 0.029 |
| GO MF | 1 | 0.02 | 0 | 0.011 |

The largest difference is 2% of one Monte Carlo standard deviation. It is about 100 times smaller than the variation between two seeds. The acceptance criterion is met.

**Assessment.**

- U-1 is not a defect in D-2, fgsea, DOSE or nf-rna's GSEA call. It needs no code change.
- nf-rna results remain a valid fgsea draw for the recorded seed.
- Exact agreement with an independent run requires both the same seed and a bit-identical input statistic.
- Writing `all_genes.tsv` at full precision would remove this source of difference. That is a precision-specification question, together with U-2, rather than a scientific correction, and it was not changed here.

**Limits.**

- The investigation did not establish why near-exact ties between permutation deviations (or with an observed ES) arise as often as two in about 20 collections × 10,000 permutations. The ranked list has 75 genes whose statistic is exactly 0, which carry zero weight in the running sum, and 75 duplicated values in total. These may contribute, but this was not shown.
- Multilevel p-values change only through their denominator `(modeFraction + 1)/(nPermSimple + 1)`, which takes the same ±1 count. The multilevel sampling itself was not involved.

## Acceptance of the v1.6.1 scientific corrections

| Item | Status |
|---|---|
| D-1 leading edge | **Accepted.** 33,544 / 33,544 sets equal the independent fgsea `leadingEdge`. ES, NES, p, p.adjust and significance are unchanged by the correction itself |
| D-2 seed | **Accepted.** The configured seed is effective, recorded and reproducible. Significant sets equal the independent same-seed fgsea in 20/20 collections. The remaining same-seed differences (U-1) are explained and lie within the Monte Carlo criterion |
| M-1 DESeq2 alpha | **Accepted.** NA pattern, padj, DEG sets and 48/48 ORA analyses equal the reference's alpha 0.05 analysis |
| U-1 | **Explained**, no code change needed. A residual limit is recorded in [U-1 investigation](#u-1-investigation) |
| Full test suite on HEAD's code | 1,091 passed, 2 skipped (environment), 0 failed |
| Version and provenance | `1.6.1` / `1.6.1 — Unreleased` / no release date. New runs record `1.6.1` and the source revision. Existing runs are unchanged |

What this acceptance does **not** cover:

- **Downstream only.** It covers the three corrections, on the downstream analyses (DE, ORA, GSEA) of the Phase 2 datasets.
- **Upstream FASTQ processing has not been validated.** Phase 3 has not started. Upstream must not be read as having passed.
- **No method is fully validated.** It does not claim that every scientific method in nf-rna is fully validated.
- **The Phase 2 v1.6.0 verdict stays FAIL.**

**Release blockers for the scientific corrections:** none known.

**Open before a release:**

- the release decision itself (date, tag), which is left to the maintainer;
- U-2 and U-3, which are tracked but are not scientific blockers;
- Phase 3, if upstream validation is a release requirement.

## Open findings

- **U-1. Explained; see [U-1 investigation](#u-1-investigation).** The differences in 2/20 collections come from the last-digit precision of the gene-level statistic that is passed to fgsea. They do not come from RNG consumption or from an algorithmic difference, and they lie far inside the Monte Carlo error. One part is not explained: why near-exact ties in permutation ES occur this often. That does not change the assessment.
- **U-2. Precision of recorded values.** `backend_summary.json` rounds `independent_filtering_threshold` to 4 decimals (jsonlite's default). This is the same specification gap as Phase 2 T-2 for mapping rates. The computation uses full precision.
- **U-3. Upstream reporting.** The DOSE `leading_edge()` ES < 0 defect and the `seed = TRUE` behaviour should be reported to DOSE upstream. Neither has been reported.
- **Not addressed here.**
  - The Phase 2 S-GS-6/S-GS-9 rule calibration (triage T-4) is a criteria question, not a software defect.
  - Phase 3 (upstream FASTQ) remains pending.
  - `rank` for ES < 0 keeps DOSE's convention (counted from the minimum, which is a miss) rather than any other convention.
