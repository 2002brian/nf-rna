# nf-rna v1.6.1

> **This release changes results.** It corrects three issues found by an
> independent scientific validation of v1.6.0. GSEA leading edges (negative-ES
> gene sets), GSEA NES/p-values/FDR (the configured seed is now effective) and
> DESeq2 `padj`, DEG sets and ORA (independent-filtering `alpha` now equals the
> padj threshold) can differ from v1.6.0. v1.6.0 results are not rewritten;
> **recompute a project to obtain corrected results**, and do not mix v1.6.0
> and v1.6.1 results in one analysis.

v1.6.1 also adds QC and provenance records that change no analysis value: a
strandedness consistency check on the HISAT2 + featureCounts route, the actual
upstream thread counts, and the R warnings of the L2 analysis. Full list of
changes: [CHANGELOG](https://github.com/2002brian/nf-rna/blob/v1.6.1/CHANGELOG.md#161--2026-10-31).

## Corrections (result-changing)

- **GSEA leading edge (GO and KEGG).** DOSE 4.4.0 reports a wrong leading-edge
  boundary for some ES < 0 gene sets. nf-rna now recomputes `core_enrichment`,
  the `leading_edge` tags and the core-member tables by the standard running-sum
  definition, checked against the fgsea ES. In the validation data, 3,554 of
  33,544 gene sets changed (all ES < 0); all 33,544 now equal an independent
  fgsea `leadingEdge`. ES, NES, p-values and FDR are not changed by this fix.
- **GSEA seed.** v1.6.0 called GSEA with `seed = TRUE`, which made DOSE re-seed
  to a fixed state, so `annotation.enrichment.gsea.seed` / `kegg.gsea.seed` had
  no effect. nf-rna now calls `set.seed(seed)` and GSEA with `seed = FALSE`.
  NES, p-values, FDR and borderline significant sets change within Monte Carlo
  error; the recorded seed is the seed actually used.
- **DESeq2 independent filtering.** `results()` used `alpha = 0.1` while DEGs
  are called at `padj < thresholds.padj` (default 0.05). `alpha` is now
  `thresholds.padj`. The filter threshold, `padj`, `padj = NA`, DEG sets and the
  ORA universe and foregrounds can change; `pvalue`, `stat` and LFCs do not.
  `thresholds.padj` must now lie strictly between 0 and 1.

## QC and provenance (analysis values unchanged)

- **Strandedness check (HISAT2 + featureCounts).** A wrongly declared
  strandedness used to give a successful run with wrong counts and no warning.
  Each sample is now also counted with featureCounts `-s 1` and `-s 2`, and the
  declaration is judged `CONSISTENT`, contradicted or apparently unstranded
  (warning), apparently stranded (advisory) or `INCONCLUSIVE`. The result
  appears in the run log, `report.html`, the GUI and the provenance. **The check
  only reports:** it never changes the declaration or the counts and never fails
  a run. The Salmon route is not checked (`NOT_EVALUATED`).
- **Actual upstream thread counts.** HISAT2 and Salmon results depend on the
  thread count. The CPUs each upstream process actually used, Salmon's
  `--threads` and the tool versions are recorded in `run_provenance.yaml`, the
  report and the GUI. CPU allocation is unchanged.
- **R warnings.** DESeq2, apeglm and ashr warnings are recorded verbatim, by
  step and with counts, in `contrast_metadata.json`, `backend_summary.json`,
  the L2 provenance, the report and the GUI. They are neither muffled nor
  treated as failures.

## Upgrade notes

- No project-file change is needed, unless a project sets `thresholds.padj` to
  0 or 1.
- `rnaseq run --reuse-upstream` cannot reuse an upstream run made by v1.6.0 or
  earlier ("frozen nfcore.local.config differs", because v1.6.1 adds trace
  columns to that file). Run the upstream stage again with v1.6.1.
- The downstream Conda lock and nf-core/rnaseq 3.26.0 are unchanged. The first
  v1.6.1 run provisions a fresh downstream environment from the same lock.

## Install

Supported runtime: Linux x86-64 or Windows WSL2 with Conda, Java 17+ and
Nextflow (see [docs/runtime.md](https://github.com/2002brian/nf-rna/blob/v1.6.1/docs/runtime.md)).
nf-rna is not on PyPI; install it from the Git tag, which lets every run record
its source commit:

```bash
conda create --yes --name nf-rna python=3.11
conda activate nf-rna
python -m pip install "nf-rna[gui] @ git+https://github.com/2002brian/nf-rna.git@v1.6.1"   # without the GUI: "git+https://github.com/2002brian/nf-rna.git@v1.6.1"
rnaseq --version      # 1.6.1
rnaseq doctor
```

## Known limitations

- **Reproducibility of upstream results.** HISAT2 places multimapping reads
  differently at different thread counts, and Salmon's EM gives different
  counts at different thread counts. Salmon 1.10.3 is also not bitwise
  reproducible at the same thread count: in the validation, reruns with
  identical input and arguments gave identical `NumReads` but different
  `EffectiveLength` and `TPM`, which reach DESeq2 through the tximport length
  offsets (L1 normalized counts differed by up to about 0.1%). Borderline DEGs of
  very small designs can change between executions. Thread counts are recorded,
  not fixed.
- **Strandedness check limits.** HISAT2 route only; warning, not failure. Its
  limits were set from one real reverse-stranded library (yeast) plus
  unstranded human and simulated data.
- **Not validated with natural data:** forward-stranded and single-end libraries
  (simulation only), mouse and other non-human FASTQ, multi-lane merging,
  preprocessed or `--skip_trimming` input, other aligners, quantifiers and
  annotation releases. L1 R warnings are not recorded.
- **Validation verdicts.** The independent validations of v1.6.0 (downstream)
  and of v1.6.1 (upstream FASTQ) did not pass under their pre-registered
  criteria, and those verdicts stand. Their failures are explained by the
  defects corrected in this release, by tool behaviour (thread dependence, the
  library-type classifier) and by criterion calibration. For identical inputs
  and settings, nf-rna's counts, DESeq2 statistics, ORA and GSEA agree with
  independent reference implementations. Details:
  [release readiness review](https://github.com/2002brian/nf-rna/blob/v1.6.1/docs/dev/release_readiness_v161.md).
- The GUI limitations of v1.6.0 still apply (single user, local only; FASTQ runs
  not verified end to end from the GUI).
