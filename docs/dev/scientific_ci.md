# Scientific CI and release gates

This document describes how nf-rna's automated checks protect the statistical
analysis. It covers the workflows, the independent regression fixtures, the
skip policy, and what a release needs on top of CI.

## Workflows

| Workflow | Trigger | What it proves | Typical time |
|---|---|---|---|
| `release-gate.yml` (fast) | every PR, push to `main` | Python tests, figure R packages, wheel/sdist build and inventories, fresh-install smoke | ~5 min |
| `science.yml` → `scientific-gate` | every PR; push to `main`, `feature/**`, `release/**`, `fix/**`; manual | Runs the **full suite against the exact locked downstream R runtime** and the pinned HISAT2/featureCounts environment, with `NF_RNA_REQUIRE_R=1`. Then `scripts/ci/pytest_gate.py` enforces the skip policy and the required scientific modules | measured locally (4 cores, 2026-11-05): 92 min for the 1,420-test suite (the slowest are the R fixture builds: the enrichment fixture alone takes ~10 min); +5–10 min to build uncached prefixes. Job timeout: 180 min |
| `science.yml` (weekly `schedule`, or manual with `rebuild_environments`) | Monday 02:41 UTC | Same as above, but without the prefix cache. The prefixes are rebuilt from the locks, so channel or Bioconductor download drift is caught before a release needs the lock | ~100–110 min |

### Runtime identity

- **Conda:** Miniforge 26.9.4-1, installed by `conda-incubator/setup-miniconda` (pinned to a commit SHA) with the channel order nf-rna documents (`conda-forge`, `bioconda`, strict priority, no `defaults`). This is the `conda` executable that nf-rna requires on `PATH` (`rnaseq.execution.check_upstream_conda`); the `tests/test_execution.py` service-lifecycle tests use it unmocked.
- **Downstream R:** built with `conda env create --yes --prefix PREFIX --file workflow/envs/locks/nf-rna-downstream-linux-64.lock.yml`, the same operation that `rnaseq.downstream_runtime` performs. No second package manager is involved. The lock is checked against `workflow/envs/locks/SHA256SUMS` first.
- **HISAT2 / SAMtools / featureCounts:** built from `workflow/envs/hisat2-featurecounts-linux-64.yml`. That file is checked against `HISAT2_LINUX_CONDA_ENV_SHA256` in `src/rnaseq/execution.py`, and the prefix is passed as `NF_RNA_HISAT2_CONDA_PREFIX`.
- **Lock verification:** `scripts/ci/verify_conda_lock.py` compares each prefix's `conda-meta` records with its lock in both directions (name, version and build). For the R prefix it also requires an R library directory for every locked `r-*` / `bioconductor-*` package.
  - This matters because Bioconductor *data* packages (`org.Hs.eg.db`, `org.Mm.eg.db`, `GO.db`) download their tarballs from bioconductor.org in a post-link script. When that download fails, the environment creation **still exits 0**, leaving an incomplete runtime. This was observed while building this gate in a network-restricted container.
- **Actions:** every `uses:` is pinned to a full commit SHA, with the release tag in a comment.
- **Python:** `ci/requirements-ci.txt` is used as a pip constraints file. It pins the exact versions of every Python package in the test environment.
- **Approved platform tools:** Nextflow 26.04.6 (same as `release-gate.yml`), Temurin Java 21, and Ubuntu's `poppler-utils`. poppler is a test-only PDF rasteriser; its version comes from the runner image.
- **Cache key:** Miniforge version plus the SHA-256 of both lock files. A lock change always rebuilds the prefixes.

### Skip policy

- With `NF_RNA_REQUIRE_R=1`, a missing `Rscript` or R package **fails** the test instead of skipping it (`tests/conftest.py`).
- `pytest_gate.py` fails the run on:
  - any failure or error;
  - fewer collected tests than the threshold;
  - a required scientific module that did not pass in full;
  - **any skip that is not listed in `ci/approved_skips.json`**.
- Every approved entry has a category and a justification. `tests/test_ci_gate.py` forbids approving R, runtime, Nextflow, SAMtools, featureCounts or poppler skips.
- Approved today:

| Category | Skip | Why |
|---|---|---|
| local-real-data | maintainer's V130-HISAT2 run | Real-data evidence is a release-checklist item |
| local-developer-tool | untracked figure-review skill | Not in Git by design |
| platform | root-only permission tests | Cannot fail as root; hosted runners are not root |
| platform | procfs start times | Restricted sandboxes only |

- The step summary and `reports/gate.json` record the counts (collected, passed, failed, errors, skipped, xfailed) and every approved skip.

## Independent scientific regression fixtures (`tests/science`)

The fixtures are simulated with the Python standard library. nf-rna is driven only through its public entry points (`validate_project` → `prepare_*` → `execute_*`).

Expectations never come from nf-rna's own R scripts. They come from three sources:

- **Simulation truth.** These checks are approximate, with generous margins.
- **Closed-form statistics in `tests/science/_science.py`.** Median-of-ratios size factors, the hypergeometric upper tail, Benjamini–Hochberg, and the GSEA running sum.
- **Reference R scripts in `tests/science/oracles/`.** They call DESeq2, tximport and AnnotationDbi directly and do not source `src/rnaseq/r`.

| Module | Fixture | Checks | Exact or tolerance |
|---|---|---|---|
| `test_deseq2_designs.py` | two-group (795 human Ensembl genes, 3 v 3); paired (4 patients, strong patient effect); batch-adjusted (2 batches × 2 conditions) | gene filter and median-of-ratios size factors; every DESeq2 column versus a direct DESeq2 call; IF threshold; `pvalue_na` / `padj_na`; Cook's applied for 3 v 3 but **not** for `~ patient + condition`; the paired and batch model matrices are actually fitted (they differ from `~ condition` and recover more planted effects); contrast direction; independent paired log-ratio signs; DEG tables follow the raw-LFC rule and not the shrunken LFC; apeglm shrinks; volcano counts, heatmap gene selection and figure manifests agree with the DEG tables | DESeq2 columns: rel 1e-8 (the same locked DESeq2 on both sides; only the 15-digit text round-trip differs). apeglm: abs 1e-6. Size factors: rel 1e-9. Direction: exact. Recovery: ≥ 0.8. Empirical false discoveries: ≤ 10%. Paired signs: ≥ 98% of DEGs and every planted DEG |
| `test_enrichment_identifiers.py` | the two-group panel, with planted GO BP terms (DNA replication up, inflammatory response down), version suffixes, unmapped, one-to-many and many-to-one IDs, 450 low-expression annotated null genes (so independent filtering gives NA padj), and an offline KEGG snapshot | precondition: the padj-non-NA universe is a strict subset of the retained-gene universe; GO and KEGG ORA universe (padj-non-NA → Entrez → ontology or snapshot genes), foreground size, gene-set sizes and hit lists; **every** ORA p-value against the hypergeometric tail and BH; the planted terms are significant in their own direction only; the GSEA ranked list against an independent implementation of the documented mapping, collapse and tie policy; ranking-summary identities; **every** GSEA setSize, enrichment score and leading edge recomputed; NES signs; fixed-seed reproducibility | ORA: rel 1e-9. Universe, ranking and leading edge: exact. ES: abs 1e-9. fgsea p-values and NES are permutation based: only seed reproducibility (byte-identical tables) and the planted signs and significance are required |
| `test_low_count_filter.py` | 200 simulated genes plus hand-written boundary genes with totals 0, 3, 4, 5, 7, 9, 10, 11, 25 and a sparse 12 | **which** genes survive the total-count ≥ 10 rule in L1 (normalized counts, VST, library totals) and in L2; every execution path receives threshold 10 | exact |
| `test_independent_filtering.py` | 1,700 genes with mean counts from 1 to 2,000 | precondition: the oracle selects different baseMean cut-offs at α = 0.05 and 0.1 on this fixture. Then nf-rna's cut-off, its NA-padj set and its padj values equal the α = 0.05 oracle; NA padj occurs exactly below the cut-off | cut-off: abs 5e-5 (jsonlite's 4 decimals); padj: rel 1e-8; NA set: exact |
| `test_tximport_input.py` | 6 Salmon `quant.sf` files with an isoform-length shift in Treatment | per-gene NumReads aggregation; L2 versus direct `tximport` + `DESeqDataSetFromTximport`; the result differs from a no-offset analysis, so the length offsets are actually used | rel 1e-9 / 1e-8 |

The committed panel `tests/science/fixtures/gene_panel.tsv` was drawn from org.Hs.eg.db 3.22.0 by `build_gene_panel.R`. The tests do not call that script.

### Mutation evidence

Each mutation was applied to a disposable copy of the source tree, and only the smallest relevant subset was run. Defects, invariants, the observed assertions and the unmodified control are recorded in `docs/dev/scientific_ci_mutations.md`.

Every fixture that a mutation depends on carries an explicit **precondition test** proving that the fixture can tell the correct behaviour from the defect:

- the oracle selects different filter cut-offs at α 0.05 and 0.1;
- independent filtering produces NA padj for annotated genes;
- boundary genes sit on both sides of the total-count threshold of 10.

Without these preconditions, a fixture that happens not to exercise the code path passes silently. That happened twice while building this gate (alpha, ORA universe).

## Release policy

| Level | When | Required evidence |
|---|---|---|
| 1. Fast PR checks | every PR | `release-gate.yml` green |
| 2. Scientific integration | every PR that touches `src/`, `workflow/`, the locks or `tests/`; and the release commit | `science.yml` `scientific-gate` green: 0 failed, 0 unapproved skips, all required modules passed. Keep the `scientific-gate-reports` artifact |
| 3. Lock reproducibility | weekly; and before tagging | a green scheduled or `rebuild_environments` run within 7 days of the tag |
| 4. Real-data end-to-end | before every release; not in CI | One real FASTQ project per supported quantification route that the release changes (HISAT2 + featureCounts; Salmon), run with the release commit. Record in `docs/dev/release_readiness_<version>.md`: commit, workflow SHA-256, lock SHA-256, run ID, `run_state` SUCCESS, delivery manifest verification, DEG and enrichment counts, and their differences from the previous release's run on the same data, with an explanation of every difference |

A full production HISAT2 run is **not** required per PR. A release without level-4 evidence is not releasable.
