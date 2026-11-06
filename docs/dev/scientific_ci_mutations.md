# Scientific gate: mutation evidence

- **Date:** 2026-11-05.
- **Base:** `ci/scientific-release-gates` working tree on `feature/v1.7.0` @ `1de489e`.
- **Runtime:** the locked downstream prefix built with `conda env create --file workflow/envs/locks/nf-rna-downstream-linux-64.lock.yml`, conda 26.9.4, `NF_RNA_REQUIRE_R=1`.

## Method

- **Copy.** Each run used a disposable copy built from `git ls-files -co --exclude-standard`, which is the set a clean checkout plus the uncommitted non-ignored work contains.
- **Precheck.** The copy was checked for the oracle scripts and fixtures before testing.
- **Mutation.** Exactly one string replacement was applied, and only the smallest relevant test subset was run.
- **Order.** Mutations ran one after another, never concurrently.
- **Control.** The same subsets ran first on an unmodified copy, and all of them passed.
- **Verdict rule.** A mutation counts as detected only when the intended test *fails* (not errors) with an assertion message that names the violated invariant. Missing files, import errors, missing R packages or fixture errors are classified INFRA and do not count.

## Results

| # | Injected defect (file: old → new) | Invariant violated | Test subset | Observed failure | Verdict |
|---|---|---|---|---|---|
| 1 | `src/rnaseq/r/l2_analysis.R`: `results(dds, contrast = c(factor, numerator, denominator)` → `c(factor, denominator, numerator)` | log2FoldChange > 0 means higher in the numerator; a planted Treatment-high gene is never called down | `test_deseq2_designs.py -k contrast_direction` | `test_contrast_direction_and_recovery_of_planted_effects`: `assert (not ({planted up genes} & {called down genes}))` | DETECTED |
| 2 | `l2_analysis.R`: `alpha = cfg$thresholds$padj` → `alpha = 0.1` | independent filtering is optimised at α = the declared padj threshold (0.05) | `test_independent_filtering.py -k optimised_at_the_declared` | `independent_filtering_threshold`: `assert 1.45 == 2.86449292762733 ± 5.0e-05` (the oracle's α = 0.05 cut-off) | DETECTED |
| 3 | `src/rnaseq/downstream.py`: `"minimum_total_count": 10` → `5` | genes with total count < 10 are removed, ≥ 10 retained (L1 and L2) | `test_low_count_filter.py -k "retains_exactly or l1_normalized"` | `retained boundary genes: ['ABOVE_11', 'ABOVE_25', 'EXACT_10', 'MID_05', 'MID_07', 'MID_09', 'SPARSE_12']`; L1 `normalized_counts.tsv` likewise | DETECTED |
| 4 | `l2_analysis.R`: `DESeqDataSetFromTximport(txi, …)` → `DESeqDataSetFromMatrix(round(txi$counts), …)` | Salmon estimates enter DESeq2 with tximport average-transcript-length offsets | `test_tximport_input.py` | `('Gene000', 'baseMean')`: L2 differs from the tximport-offset oracle | DETECTED |
| 5 | `src/rnaseq/r/gsea_analysis.R`: `order(-abs(candidates$stat), …` → `order(abs(candidates$stat), …` | an Entrez ID reached from several genes keeps the largest \|stat\| (ties by gene_id) | `test_enrichment_identifiers.py -k "ranked_list and gsea_go"` | `test_gsea_ranked_list_follows_the_mapping_policy_exactly[gsea_go]`: ranked list ≠ the independently built list | DETECTED |
| 6 | `src/rnaseq/r/go_analysis.R`: `nf_rna_ora_universe(tested, tested_mapping, genes)` → `nf_rna_ora_universe(retained, mapping, genes)` | the ORA universe is the padj-non-NA genes → Entrez → genes annotated in the ontology | `test_enrichment_identifiers.py -k go_ora_universe` | all 6 `test_go_ora_universe_and_hypergeometric_statistics[*]`: `universe_genes` `assert 654 == 449` (BP; MF 603 ≠ 405; CC 640 ≠ 436) | DETECTED |

Controls on unmodified copies:

- Mutations 1, 3, 4 and 5: 15/15 passed in 1,178 s.
- Mutations 2, 5 and 6 on the revised fixtures: passed in 759 s.

## Coverage gaps the mutations exposed (now closed)

A first round, run before the fixes, produced misleading results. They are recorded here because they show why mutation testing and fixture preconditions are both needed.

- **Missing oracles (INFRA, all reference-based results invalid).** `.gitignore` excludes `**/reference/` and `**/references/`, so the oracle scripts were absent from the copies. They now live in `tests/science/oracles/`, and `tests/test_ci_gate.py::test_every_file_the_scientific_gate_needs_is_visible_to_git` fails if any gate file becomes ignored.
- **Filter threshold not detectable.** No simulated gene had a total count between 5 and 9, so a 10 → 5 change was invisible. Added `test_low_count_filter.py`, with hand-written boundary genes.
- **Alpha not detectable.** On the design fixtures DESeq2 chose the same baseMean cut-off at α 0.05 and 0.1. Added `test_independent_filtering.py`, whose first test proves the fixture discriminates.
- **ORA universe not detectable.** No gene of the enrichment fixture had NA padj, so "retained" and "padj-non-NA" were the same set. The enrichment fixture now carries 450 low-expression annotated genes (205 NA padj in the oracle), plus a precondition test.

## Production defect found while building the fixtures (not fixed here)

- **Where:** `src/rnaseq/r/l1_analysis.R:41`.
- **What it does:** it uses `vst(dds, blind = TRUE)` whenever `nrow(dds) >= 1000`.
- **Why that fails:** DESeq2's `vst()` additionally requires at least 1000 genes with mean normalized count > 5. Otherwise it stops with `less than 'nsub' rows with mean normalized count > 5`.
- **Consequence:** a low-depth project with ≥ 1000 retained genes, fewer than 1000 of them expressed above 5, therefore fails L1.
- **Status:** reproduced with a 1,200-gene fixture. The scientific fixture was resized so it stays inside the supported domain. The production fallback condition is a separate change for the L1 owner.
