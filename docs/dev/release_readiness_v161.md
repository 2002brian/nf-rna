# nf-rna v1.6.1 release readiness review

Date: 2026-10-31. Reviewed tree: branch `fix/v1.6.1-reliability`, the commit
that adds this file (its parent is `89a4510`). This review adds no analysis
method, changes no threshold, reruns no Phase 3 run and modifies no existing
run. It relies on existing evidence, plus the version and packaging checks
listed in §10.

**Verdict: ready for release, subject to the maintainer's decision and to a
green `Release gate` CI run on the pushed `main` (the v1.6.1 commits have never
run in CI).** No release blocker is known (§8). The frozen FAIL verdicts of the
Phase 2 and Phase 3 validations stand and are disclosed as known findings (§6).

## 1. Git state and ancestry

The working tree is clean. `v1.6.0` (`ff4fe5c`) is `main` and `origin/main`.
`fix/v1.6.1-reliability` is a linear descendant of `v1.6.0`, with no merge
commits:

| Commit | Branch where it was made | Content |
|---|---|---|
| `3cf0386` | `fix/v1.6.1-scientific` | D-1, D-2, M-1 corrections, regression tests, CHANGELOG |
| `7d85046` | `fix/v1.6.1-scientific` | Version `1.6.1 — Unreleased`; offline KEGG GSEA fixture |
| `1cbe8fb` | `fix/v1.6.1-scientific` | Correction evidence (`validation/corrections_v161/`) |
| `033ce62` | `fix/v1.6.1-scientific` | U-1 investigation and acceptance of the corrections |
| `9c37d84` | `fix/v1.6.1-reliability` | F-1, F-2, F-3 code and tests |
| `89a4510` | `fix/v1.6.1-reliability` | F-1/F-2/F-3 report, contract, GUI docs, acceptance evidence |
| this commit | `fix/v1.6.1-reliability` | Release finalization: date, release notes, checklist, this review |

All commits since `v1.6.0` are authored and committed by Brian Wu and carry no
trailers. `fix/v1.6.1-scientific` (`033ce62`) is fully contained.

Evidence branches, **not** merged and not needed in `main`:

| Branch | Head | Base | Adds |
|---|---|---|---|
| `validation/v1.6.0` | `b886f11` | `v1.6.0` | Phase 1 and Phase 2 validation (`validation/`, `docs/validation/report.md`) |
| `validation/v1.6.1-phase3` | `64c6977` | `033ce62` | Phase 3 upstream validation (`validation/phase3/` only) |

Neither touches `src/` or `workflow/`, and neither overlaps a file changed on
`fix/v1.6.1-reliability`.

**Merge recommendation.** Fast-forward `main` to `fix/v1.6.1-reliability`
(`git merge --ff-only`). Keep the two validation branches separate. The
v1.6.1 documents cite them by branch and commit, so push them as branches if
those citations should resolve for readers outside this machine. Do not merge
them into `main` and do not tag them.

## 2. Version and metadata

| Item | Value |
|---|---|
| `src/rnaseq/__init__.py` (`hatch` version source) | `1.6.1` |
| `CITATION.cff` | `version: "1.6.1"`, `date-released: 2026-10-31` |
| `CHANGELOG.md` first heading | `## 1.6.1 — 2026-10-31` |
| `tests/test_release_metadata.py` | passes (§10) |
| Wheel and sdist metadata | `Version: 1.6.1`, extra `gui` → `streamlit<2,>=1.66` (§10) |
| README_latest.md, README_zh-TW.md, docs/quickstart.md, docs/gui.md, docs/runtime.md | name `v1.6.1` as the latest release and install from `@v1.6.1` |
| Downstream Conda locks, nf-core/rnaseq pin | unchanged since `v1.6.0` |
| New runs | record `nf_rna_version` `1.6.1` and the source revision |

The release date is 2026-10-31. If the tag is made on a later day, change the
date in `CHANGELOG.md` and `CITATION.cff` together (the metadata test requires
both to agree) and the CHANGELOG anchors `#161--2026-10-31` in the READMEs and
quickstart.

The FIX-* validation runs were made from `3cf0386`, before the version change,
and record `nf_rna_version` `1.6.0` with `source_revision` `3cf0386…`. They are
not rewritten. `docs/dev/scientific_corrections_v161.md` says "Phase 3 has not
started" and "1.6.1 — Unreleased"; both were true at `033ce62` and that record
is kept as written.

## 3. What v1.6.1 changes

### 3.1 Changes that alter analysis values

From the Phase 2 validation of v1.6.0
(`docs/dev/scientific_corrections_v161.md`):

| ID | Defect in v1.6.0 | Correction | Values that change |
|---|---|---|---|
| D-1 | GSEA leading edge from DOSE 4.4.0 `leading_edge()`, wrong for some ES < 0 sets | nf-rna recomputes the leading edge by the standard running-sum definition, checked against the fgsea ES | `core_enrichment`, `leading_edge` tags, core-member tables and the `leading_edge` flag of `members.tsv`, for ES < 0 sets only (3,554 of 33,544 sets in the Phase 2 data, 162 of them significant). ES, NES, p, FDR are not changed by this fix |
| D-2 | `seed = TRUE` made DOSE re-seed to a fixed state; the configured seed had no effect | `set.seed(seed)`, then GSEA with `seed = FALSE` | GSEA NES, p-value, p.adjust and borderline significance, at the Monte Carlo level; ES unchanged |
| M-1 | DESeq2 independent filtering used `alpha = 0.1` while DEGs are called at `padj < 0.05` | `alpha = thresholds.padj` | Filter threshold, `padj`, `padj = NA` pattern, DEG sets, ORA universe, foregrounds and significant terms (Phase 2 data: up to 25 DEGs and 39 ORA terms per analysis). `pvalue`, `stat`, LFCs and GSEA unchanged |

Users must recompute to obtain corrected results. v1.6.0 runs are not
rewritten.

### 3.2 Changes that do not alter analysis values

From the Phase 3 upstream validation
(`docs/dev/scientific_reliability_v161.md`):

| ID | Finding | Change | Evidence that values are unchanged |
|---|---|---|---|
| F-1 | A wrong strandedness declaration on the HISAT2 route ran to SUCCESS with wrong counts and no warning | Per-sample check from featureCounts `-s 1` / `-s 2` on the count-only BAM; `CONSISTENT`, warning (`CONTRADICTS_DECLARATION`, `APPEARS_UNSTRANDED`), advisory (`APPEARS_STRANDED`) or `INCONCLUSIVE`; report only | Counts, featureCounts summaries and L1 tables byte-identical to the Phase 3 runs (REL_D5_NEG, REL_D4_HF) |
| F-2 | Thread counts change upstream results but were not recorded | Trace columns `process,tag,cpus,memory`; actual CPUs per process, Salmon `cmd_info.json` threads and tool versions in provenance, contract, report and GUI | CPU allocation and tool arguments unchanged; trace CPUs equal the executed `hisat2 -p` and Salmon `--threads` |
| F-3 | DESeq2/apeglm warnings stayed in the Nextflow work directory | Warnings recorded verbatim, by step and with counts, in `contrast_metadata.json`, `backend_summary.json`, L2 provenance, report and GUI; not muffled or escalated | D3 DE tables byte-identical to P3_D3_HF; D1 7 L1 and 125 L2 tables (DE, ORA, GSEA) byte-identical to FIX-D1 |

New outputs only: `qc/strandedness/` (HISAT2 route) and new keys in existing
JSON/YAML records. The HISAT2 route runs one extra process per sample.

## 4. Scientific correctness evidence

| Area | Evidence | Source |
|---|---|---|
| Downstream DE, normalization, VST/PCA, shrinkage, annotation, ORA, GSEA statistics (v1.6.0) | Independent reference agreement to ≤ 5.1e-15 relative for every deterministic quantity; Phase 2 implementation concordance PASS for all areas except the GSEA leading edge (D-1) and the stochastic GSEA rules | `validation/v1.6.0` `b886f11`, `docs/validation/report.md` |
| D-1 | 33,544 / 33,544 leading edges equal the independent fgsea `leadingEdge` (v1.6.0: 29,990); the correction alone changes no ES/NES/p/FDR | `validation/corrections_v161/` |
| D-2 | Significant sets equal independent fgsea with the same seed in 20/20 collections (seeds 1 and 42); residual same-seed differences (U-1) ≤ 2% of one Monte Carlo SD, from input-statistic precision | `docs/dev/scientific_corrections_v161.md` |
| M-1 | padj NA pattern, padj (≤ 5e-15) and DEG sets equal the reference alpha 0.05 analysis in 6/6 contrasts; ORA in 48/48 analyses | same |
| Upstream FASTQ (v1.6.1 code at `033ce62`) | For identical inputs and execution settings, counts, import, normalization, DE and annotation equal the independent reference, mostly to ≤ 5e-15 and by exact sets. HISAT2 re-alignment at nf-rna's `-p 4` reproduces nf-rna's counts and summaries in 16/16 samples; the reference's DE script on nf-rna's counts gives ≤ 4.9e-15 and the same DEG set | `validation/v1.6.1-phase3` `64c6977`, `validation/phase3/report.md` |
| F-1/F-2/F-3 | 61 acceptance checks: 60 PASS, 1 FAIL kept and explained (§5) | `validation/reliability_v161/results/acceptance.tsv` |
| Test suite | 1,114 passed, 2 skipped (environment-gated), 0 failed, locked runtime `0e721ff3…`, `NF_RNA_REQUIRE_R=1`, code `9c37d84`. The commits after it change documentation, validation evidence and release metadata only | `docs/dev/scientific_reliability_v161.md` |

Every result is conditional on the locked versions (R 4.5.3, Bioconductor
3.22, DOSE 4.4.0, clusterProfiler 4.18.4, fgsea 1.36.2, nf-core/rnaseq 3.26.0,
Salmon 1.10.3, HISAT2 2.2.3, Subread 2.0.6), the KEGG 2026-10-28 snapshots and
the datasets tested. The independent references were written by an AI session
of the same model family; human review of the reference scripts is
recommended before external citation (Phase 2 §8).

## 5. Reproducibility limitations

### 5.1 Thread count

- **HISAT2 2.2.3** places multimapping reads differently at different `-p`
  (`-p 4` vs `-p 8`: 9,419 differing records on the same reads), and the
  difference reaches counts (D3: 3–8 genes per sample) and borderline DESeq2
  values. At the same `-p` the alignments were identical (Phase 3 N-1).
- **Salmon's EM** gives different `NumReads` at different thread counts
  (Phase 3 N-3: 57,165 of 2,579,656 transcript rows between `-p 1` and `-p 8`;
  gene-level p99 |log2| difference 0.001). In a 2 vs 2 design this moved the
  independent-filtering threshold from 46.9 to 12.2 and left only 24 of 37 DEGs
  in common between `-p 8` and nf-rna's `-p 6`; `-p 1` and `-p 6` gave the same
  37 DEGs.
- nf-rna fixes HISAT2 at `cpus 4` and nf-core runs Salmon with 6, but
  `execution.max_cpus` (Nextflow `resourceLimits`) can lower both. v1.6.1
  records the actual counts (F-2); it does not fix them.

### 5.2 Salmon `EffectiveLength` at the same thread count

Established during the F-1/F-2/F-3 acceptance; no new analysis was run for
this review.

- **Observation.** The v1.6.1 Salmon-route rerun of D4 (REL_D4_SAL) and the
  Phase 3 run P3_D4_SAL used the same Salmon command, `cmd_info.json`,
  parameters, inputs and `--threads 6`. `NumReads` and `estimated_counts.csv`
  were identical. `EffectiveLength` and `TPM` were not (max relative difference
  0.33% and 0.19%), because the estimated fragment-length mean differed
  (165.897 vs 165.840).
- **Reproduction.** Three standalone Salmon 1.10.3 runs on identical trimmed
  reads, index and arguments at `--threads 6` differed pairwise in
  `EffectiveLength` for 105–117 of the 124 transcripts of the D4 test
  reference, never in `NumReads`
  (`validation/reliability_v161/results/salmon_same_thread_repro.tsv`). This is
  Salmon's own run-to-run variation, not an nf-rna change.
- **Downstream path.** nf-rna imports Salmon with `tximport(type = "salmon")`
  and `DESeqDataSetFromTximport`. Gene counts are sums of `NumReads` and stay
  identical; the abundance-weighted gene lengths come from `EffectiveLength`
  and enter DESeq2 as per-gene normalization offsets.
- **Measured effect.** D4 L1: normalized counts differ by at most 1.03e-3
  relative, VST by at most 6.4e-5.
- **Not measured.** The effect on DESeq2 statistics, padj, DEG sets and
  enrichment of a same-thread Salmon rerun. The thread-count case (§5.1) shows
  that borderline DEGs of very small designs can move when the
  independent-filtering threshold jumps; a same-thread rerun is expected to
  perturb the offsets far less, but this was not shown.
- **Consequence.** Byte-identical Salmon-route reruns cannot be expected, even
  with identical settings. The acceptance FAIL for `quant.sf` (REL_D4_SAL) is
  kept as FAIL with the triage row `EXPLAINED`.

### 5.3 Other

- GSEA p-values are a Monte Carlo draw for the recorded seed; other seeds give
  equally valid, different p-values, and sets near p.adjust 0.05 can change
  (D-2, U-1).
- `backend_summary.json` rounds `independent_filtering_threshold` to 4 decimals
  (U-2); the computation uses full precision.

## 6. Known findings

### 6.1 Frozen validation verdicts (unchanged)

Both independent validations **do not PASS under their frozen criteria**, and
those verdicts are not rewritten.

| Validation | Frozen FAIL | Root cause | Status in v1.6.1 |
|---|---|---|---|
| Phase 2 (downstream, v1.6.0) | GSEA leading edge E-GS-4 | D-1 | Corrected |
| Phase 2 | GSEA stochastic rules S-GS-6 (GO 1/12, KEGG 2/4), S-GS-9 (GO 9/12) | Rules stricter than Monte Carlo variability (triage T-4); related D-2 | D-2 corrected; criterion issue remains |
| Phase 3 (upstream, v1.6.1 code) | E-UP-4, E-UP-5, SE-4 (counts and summaries not identical) | HISAT2 thread dependence (N-1); exact match at matched `-p` in 16/16 | Tool behaviour; thread counts now recorded (F-2) |
| Phase 3 | N-UP-8, FC-3 (HISAT2 route DE values, retained genes) | Propagation of N-1; ≤ 4.9e-15 on matched counts | Tool behaviour |
| Phase 3 | SA-1 D4 standalone `-l A` (Salmon `IU` vs protocol `ISR`) | Salmon's classifier vs nf-core's 0.8 rule on the same statistic (N-2); nf-rna's call (`reverse`) is the correct one here | Tool behaviour; near nf-core's threshold |
| Phase 3 | HF-1 / ST-3 / SE-3 gene-level accuracy vs `unique` truth | Criterion denominator too narrow (T-1) | Criterion issue |
| Phase 3 | ST-2 D4 wrong-direction share 0.174–0.179 vs ≤ 0.1 | Antisense overlaps in the compact yeast reference (T-2) | Criterion issue |
| F-1/F-2/F-3 acceptance | REL_D4_SAL `quant.sf` not byte-identical | Salmon same-thread variation (§5.2) | Explained, kept as FAIL |

No FAIL was traced to an nf-rna computation error that remains open.

### 6.2 Other findings

- **Upstream reuse across the version boundary.** v1.6.1 appends a
  `trace.fields` line to the frozen `nfcore.local.config`, so
  `rnaseq run --reuse-upstream` refuses a source run made by v1.6.0 or earlier
  ("frozen nfcore.local.config differs"). Refusing is the fail-closed
  behaviour, and such an upstream has neither the strand evidence nor the
  `cpus` trace column. The error hint attributes the difference only to
  CPU/memory limits, which is misleading in this case. Documented in the
  CHANGELOG *Compatibility* section; reuse between v1.6.1 runs is unaffected.
- **N-4 apeglm line-search warnings** (D3 featureCounts route: 6) are now
  recorded (F-3). No shrunken LFC was NA or non-finite; shrinkage status
  `success`.
- **U-1** explained (precision of the statistic passed to fgsea); **U-2**
  (rounded threshold in JSON) and **U-3** (DOSE defects not yet reported
  upstream) open, not scientific blockers.

## 7. Scope of validation and applicability limits

Validated with natural data: human (Ensembl 116 / GRCh38) and mouse count
matrices downstream; human (unstranded) and yeast (reverse-stranded, compact
test reference) FASTQ upstream; paired-end; one sample per FASTQ pair;
nf-core/rnaseq 3.26.0 with Salmon, and HISAT2 + featureCounts.

Not validated, or validated only conditionally:

| Area | Status |
|---|---|
| Forward-stranded libraries | Simulation only (D5-F) |
| Single-end libraries | Simulation only (D5 SE) |
| Mouse FASTQ (and other species upstream) | Not run; mouse validated as counts only |
| Multi-lane merging (several FASTQ pairs per sample) | Not run |
| Preprocessed or `--skip_trimming` input | Not run |
| Strandedness `auto` near nf-core's 0.8 threshold, antisense-rich protocols | Not explored; D4 sits at 0.81–0.85 |
| F-1 limits (3×, 1.5×, 1,000 fragments) | Calibrated on one real stranded library (yeast reverse), unstranded human, simulated F/R/U/SE; low-depth, natural forward and other organisms not shown |
| F-1 on the Salmon route | Not evaluated (`NOT_EVALUATED`) |
| Other aligners, quantifiers, Ensembl/GENCODE releases | Not run |
| DE stability under re-execution for small designs | One dataset (2 vs 2), thread-count case only |
| Exception paths (gene-wise dispersion fallback, shrinkage failure, `NETWORK_UNAVAILABLE`) | Not forced in Phase 2 or 3; the ashr fallback was reached naturally (Phase 2 SH-4 PASS) |
| Continuous covariates, SYMBOL input, ORA NOT_APPLICABLE, mapping-QC WARNING/BLOCKED | Did not arise in Phase 2 |
| L1 R warnings (VST, PCA) | Not recorded |
| GUI | As in v1.6.0: verified end to end on raw-count projects; FASTQ runs not started or browsed from the GUI. The new strandedness, thread and R-warning views were checked on the acceptance runs |
| Platform | Linux x86-64 / WSL2 only (unchanged) |

## 8. Release blockers

**None known.**

- No scientific-correctness blocker: every frozen FAIL is either corrected
  (D-1, D-2, M-1) or traced to tool behaviour or criterion calibration.
- F-1, the one silent failure mode found upstream, is now reported as a
  warning. Whether a warning (rather than a hard failure) is sufficient is a
  maintainer decision; the evidence does not support hard failure as the
  default.
- Pre-tag gate, not a blocker: `Release gate` CI must pass on the pushed
  `main`. No v1.6.1 commit has run in CI.

## 9. Future nonblocking improvements

1. Optional strict mode for F-1 (fail on `CONTRADICTS_DECLARATION`), after
   calibration on natural forward-stranded, single-end, low-depth and
   antisense-rich libraries and other species.
2. A strand check for the Salmon route (for example from `lib_format_counts.json`
   under `-l A`, or a featureCounts pass on a pseudo-BAM).
3. An option to pin upstream thread counts, or a warning when
   `execution.max_cpus` lowers HISAT2's or Salmon's fixed CPUs.
4. Measure the DE and enrichment effect of same-thread Salmon reruns
   (§5.2), and document Salmon's `EffectiveLength` variation in the methods text.
5. A specific `--reuse-upstream` error when the source run was made by another
   nf-rna version.
6. Record L1 R warnings (VST, PCA) as F-3 does for L2.
7. Full-precision values in `backend_summary.json` and `all_genes.tsv` (U-1, U-2).
8. Report the DOSE `leading_edge()` and `seed = TRUE` behaviour upstream (U-3).
9. Validation protocol: fix thread counts in exact-match criteria, calibrate
   stochastic GSEA criteria with more seeds, widen the accuracy denominators
   (T-1) and the wrong-direction threshold for compact references (T-2).
10. Validate the untested scopes of §7 (multi-lane, mouse FASTQ, natural SE and
    forward libraries, preprocessed input, exception paths).

## 10. Checks run for this review

Only version and packaging checks; no analysis, no full suite, no FASTQ run.

| Check | Result |
|---|---|
| `git status`, ancestry, authorship (`git log v1.6.0..HEAD`) | clean; linear on `v1.6.0`; 6 + 1 commits by Brian Wu, no trailers |
| `git diff --stat v1.6.0 HEAD -- workflow/envs` | empty (locks unchanged) |
| `tests/test_release_metadata.py` | 2 passed |
| `python -m build` (wheel and sdist) | `nf_rna-1.6.1-py3-none-any.whl` (1.5 MB), `nf_rna-1.6.1.tar.gz` (19.8 MB, includes the GUI screenshots as in v1.6.0); METADATA `Version: 1.6.1`, extra `gui` → `streamlit<2,>=1.66`. The wheel contains `rnaseq/workflows/{main.nf,hisat2_featurecounts.nf,nextflow.config}` (with `STRAND_EVIDENCE`) and `rnaseq/r/{l2_analysis.R,gsea_leading_edge.R}` |
| `tests/release_smoke.py` on the wheel and the sdist, each in a fresh venv | pass (both) |
| `rnaseq --version` from each installed artifact | `1.6.1` (both) |
