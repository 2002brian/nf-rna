# Scientific reliability hardening for v1.6.1: F-1, F-2, F-3

**Branch.** `fix/v1.6.1-reliability`, created from `fix/v1.6.1-scientific` at `033ce62`, which is v1.6.0 plus the D-1/D-2/M-1 corrections. The Git ancestry was checked before branching:

- `v1.6.0` is an ancestor of `fix/v1.6.1-scientific`.
- `validation/v1.6.1-phase3` descends from `033ce62` but is **not** merged; it adds only `validation/`, with no `src/` or `workflow/` changes.

**Source of the findings.** The Phase 3 upstream FASTQ validation, `validation/phase3/report.md` on branch `validation/v1.6.1-phase3` (commit `64c6977`). Its verdicts are unchanged by this work: Phase 3 still does not PASS under its frozen criteria, and every FAIL listed there stands.

| Commit | Content |
|---|---|
| `9c37d84` | Code and tests for F-1, F-2 and F-3 |
| Later commit | This document, contract/GUI/CHANGELOG updates, acceptance evidence (`validation/reliability_v161/`) |

## Summary

| ID | Root cause (Phase 3) | Correction | Changes analysis values? |
|---|---|---|---|
| F-1 | The HISAT2 route never checks the declared strandedness. A reverse library declared `forward` ran to SUCCESS with 7–9% of fragments assigned, all to the wrong genes, and no warning | A per-sample strandedness consistency check from the same count-only BAM counted with featureCounts `-s 1` and `-s 2`. It reports only | No: counts, declaration and run status are unchanged |
| F-2 | HISAT2 multimapper placement and Salmon EM depend on the thread count (Phase 3 N-1, N-3), but the effective thread counts were not in the provenance; the trace had no `cpus` column | Trace columns `process,tag,cpus,memory`. The actual CPUs per upstream process, Salmon's `cmd_info.json` threads and the tool versions are recorded in provenance, the contract, the report and the GUI, with a reproducibility note | No: CPU allocation and tool arguments are unchanged |
| F-3 | DESeq2/apeglm warnings (for example apeglm's line-search warnings on D3) stayed in the Nextflow work directory | Warnings are recorded verbatim, by step, with call and count, in the L2 outputs, provenance, report and GUI | No: warnings are neither muffled nor escalated |

## F-1: strandedness consistency check

**Evidence.**

- A new workflow process, `STRAND_EVIDENCE` (`workflow/hisat2_featurecounts.nf`), takes the same count-only BAM, GTF and options as `FEATURECOUNTS` (`-t exon -g gene_id -Q 0 --primary`, plus the pair options), and counts with `-s 1` and with `-s 2`.
- It publishes only the two `.summary` files, to `upstream/hisat2_featurecounts/qc/strandedness/`.
- Its output feeds neither `ASSEMBLE_COUNTS` nor MultiQC. The counts used downstream still come from `FEATURECOUNTS` with the declared strandedness.

**Rule.** The rule is in `rnaseq.hisat2_featurecounts.strandedness_verdict`; the run-level summary is in `evaluate_strandedness`. With `r = assigned(-s 1) / assigned(-s 2)` per sample:

| Condition | Observed |
|---|---|
| assigned(-s 1) + assigned(-s 2) < 1,000 | insufficient evidence → `INCONCLUSIVE` |
| r ≥ 3 | forward |
| r ≤ 1/3 | reverse |
| 1/1.5 ≤ r ≤ 1.5 | unstranded |
| otherwise | ambiguous → `INCONCLUSIVE` |

| Declared | Observed | Status | Severity |
|---|---|---|---|
| any | same as declared | `CONSISTENT` | pass |
| forward or reverse | the opposite strand | `CONTRADICTS_DECLARATION` | warning |
| forward or reverse | unstranded | `APPEARS_UNSTRANDED` | warning |
| unstranded | forward or reverse | `APPEARS_STRANDED` | advisory: unstranded counting is valid; check the protocol |
| any | insufficient evidence or ambiguous | `INCONCLUSIVE` | no judgement made |

The run-level status is the most severe sample: `WARNING`, then `ADVISORY`, then `INCONCLUSIVE`, then `CONSISTENT`.

**Why a ratio and these limits.**

- The ratio compares two strand hypotheses on the same reads, so it is independent of species, depth, mapping rate and rRNA content. A fixed assigned-read percentage is not.
- Opposite-strand gene overlaps assign correctly sensed fragments to the antisense neighbour, which moves `r` towards 1. In Phase 3:
  - a simulated, perfectly reverse yeast library gave r = 0.075–0.096, and 99.8% of its wrong-direction reads lay on annotated opposite-strand overlaps;
  - the real D4 library gave 0.17–0.18;
  - unstranded data gave 0.99–1.00.
- Overlaps cannot make the opposite direction dominate. A contradiction therefore requires a 3× dominance of the wrong direction; the correctly declared libraries were 5.6× to 13× away from that. A library that is truly stranded but overlap-dense can at worst fall into the ambiguous band and become `INCONCLUSIVE`.
- The 1,000-fragment minimum avoids judging samples with too few reads (D5-R subsampled to 0.2%: 394 fragments → `INCONCLUSIVE`).

**What happens with the result.**

- Recorded in `strandedness_qc.json`, the handoff manifest, the frozen contract (`upstream_qc.strandedness`) and `run_provenance.yaml`.
- Each warning or advisory is written to the run log (`logs/rnaseq.log`) as `WARNING:` or `NOTE:`.
- Shown in `report.html` ("Strandedness check", in the Reference and quantification section) and in the GUI (*QC & exploratory*).
- The declared strandedness is never changed, nothing is recounted, and the run is never failed.

**Hard failure: not enabled.** Phase 3 provides one real stranded library (yeast, reverse) plus simulations. That is not enough to rule out real libraries that sit near the limits, such as antisense-rich protocols, low depth or other organisms. A warning is the conservative default. An opt-in strict mode can be considered once more stranded datasets have been validated.

**Salmon route.** The status is `NOT_EVALUATED`, with the reason. With a declared library type, Salmon's `lib_format_counts.json` records only compatible orientations, so it is not evidence against the declaration. Under `auto`, nf-core infers the type itself (Phase 3 N-2). The library type Salmon used per sample is recorded.

## F-2: actual thread provenance

- **Trace columns.** `render_trace_fields_config()` appends `trace.fields = 'task_id,…,cpus,memory'` (Nextflow's defaults plus `process`, `tag`, `cpus` and `memory`) to the frozen `nfcore.local.config`. That file is passed to nf-core, to the HISAT2 route and to the downstream run. Trace fields do not alter scheduling, `task.cpus` or any tool argument.
- **Actual thread counts.** `upstream_execution_threads()` reads the run's own trace: the latest `pipeline_info/execution_trace_*.txt` for nf-core, or `provenance/upstream.trace.txt` for the HISAT2 route. For each process that uses threads in counting or quantification it records the CPUs of every completed or cached task:
  - HISAT2 route: `FASTP_PREPARE`, `HISAT2_ALIGN`, `FEATURECOUNTS`, `STRAND_EVIDENCE`;
  - Salmon route: `TRIMGALORE`, `FASTP`, `SALMON_QUANT`.

  These tools receive `task.cpus` as their thread argument (`hisat2 -p`, `featureCounts -T`, `salmon --threads`), so the trace value is the thread count actually used, after any `resourceLimits` cap. It is not taken from `execution.max_cpus`.
- **Other records.** Salmon's `--threads` from each sample's `cmd_info.json`, and the tool versions per process (from nf-core's versions file; for the HISAT2 route, the pinned versions).
- **Where it appears.** `run_provenance.yaml` (`execution_threads`), the frozen contract (`upstream_qc.execution_threads`), `report.html` ("Upstream execution threads", with the reproducibility note) and the GUI *Provenance* table.
- **Older or reused runs.** Without a `cpus` column, or with an upstream reused from another run, the record is `NOT_RECORDED` with the reason.

## F-3: R warnings in L2

**Where warnings are captured.** `l2_analysis.R` wraps the following calls in `withCallingHandlers`:

- the shared fit `DESeq()` and the gene-wise fallback;
- each `results()`;
- each `lfcShrink()` (apeglm, releveled apeglm, ashr) and the relevel refit.

The handler records each warning's step, message and call, and returns without `muffleWarning`. The warning therefore still reaches R's own output, and the computation is unchanged. Identical warnings are collapsed with a count.

**Where the record appears.**

- `contrast_metadata.json` `r_warnings`: `status` (`NONE`/`WARNINGS`), `meaning`, `fit` and `contrast` lists. It is delivered under `methods_and_versions/l2/contrasts/<id>/`.
- `backend_summary.json`: `r_warning_status`, `r_warning_count`.
- The L2 `scientific_provenance.json` result summary.
- `report.html`: under each contrast's shrinkage line.
- The GUI: an expander on the *Differential expression* page.

**What stays the same.**

- An error keeps its existing handling: the documented dispersion fallback, `shrinkage.status: failed` with the error message, or a failed analysis.
- A warning is never turned into a failure.

## Files changed

| File | Change |
|---|---|
| `workflow/hisat2_featurecounts.nf` | `STRAND_EVIDENCE` process and its call |
| `src/rnaseq/hisat2_featurecounts.py` | Strandedness QC constants, `featurecounts_assigned`, `strandedness_verdict`, `evaluate_strandedness` |
| `src/rnaseq/execution.py` | Strandedness QC in both handoff manifests; `TRACE_FIELDS`, `render_trace_fields_config`, `upstream_execution_threads` |
| `src/rnaseq/service.py` | Trace fields in the frozen config; `upstream_qc` in the contract and run provenance; run-log warnings |
| `src/rnaseq/run_status.py` | `STRAND_EVIDENCE` in the featureCounts progress stage |
| `src/rnaseq/r/l2_analysis.R` | Warning capture and records |
| `src/rnaseq/workflow_support.py` | Report: strandedness check, upstream threads, R warnings |
| `src/rnaseq/gui/results.py`, `src/rnaseq/gui/pages.py` | GUI: strandedness on QC, R warnings on DE, threads on Provenance |
| `tests/test_scientific_reliability_v161.py` (new, 23 tests); `tests/test_l2.py` | Regression tests |
| `docs/scientific_contract.md` | Strandedness check, R warnings, actual threads. Also two corrections of statements that v1.6.1 had made stale: GSEA now uses `seed = FALSE` after `set.seed(seed)` (D-2), and independent filtering is optimised at `alpha = thresholds.padj` (M-1) |
| `docs/gui.md`, `CHANGELOG.md` | User-facing description |

## Evidence

### Tests

- **New tests:** `tests/test_scientific_reliability_v161.py` (23 tests). They cover:
  - the verdict on all Phase 3 evidence cases (D3, D4, D5-F/R/SE consistent; NEG and D4-as-forward contradicting; D5-U-as-reverse appears unstranded; D5-R-as-unstranded appears stranded);
  - inconclusive cases (minimum evidence; ratio 2);
  - one-sided evidence and the run-level priority;
  - that the workflow's strand evidence uses exactly the counting options of `FEATURECOUNTS` apart from `-s` and `-o`, and feeds nothing downstream;
  - trace fields and thread extraction for both routes, plus runs without `cpus`;
  - report rendering;
  - that the R script records warnings without `muffleWarning`.
- **Extended test:** the real L2 run in `tests/test_l2.py` now checks the `r_warnings` record and its count in `backend_summary.json`.
- **Focused existing modules, all passing:** `test_hisat2_featurecounts`, `test_workflow`, `test_execution`, `test_workflow_support`, `test_run_status`, `test_gui_results`, `test_gui` and `test_workflow_assets` (170 passed, 1 environment skip); `test_service`, `test_resource_policy`, `test_l2` and `test_featurecounts_pinned_conda_semantics` (116 passed).
- **Full suite**, run once at integration on code `9c37d84` with the documentation of this commit, locked runtime `0e721ff3…`, `NF_RNA_REQUIRE_R=1`: **1,114 passed, 2 skipped, 0 failed** (26 min 10 s). The two skips are the environment-gated tests also skipped before: real SAMtools/featureCounts, and a real SJ project.

### Acceptance runs

The acceptance runs used nf-rna `9c37d84`, installed from the commit. They are small runs with no large FASTQ pipeline. Data are under `/home/brian/nf-rna-validation/v161-reliability/`; the script is `validation/reliability_v161/scripts/compare_acceptance.py`; the results are `validation/reliability_v161/results/acceptance.tsv`, with 61 checks and 1 triage row.

| Run | Input | Result |
|---|---|---|
| REL_D5_NEG | D5 reverse library declared `forward` (Phase 3 negative control) | **`WARNING`**: both samples `CONTRADICTS_DECLARATION` (r 0.096, 0.075). Declaration kept; 2 run-log warnings; shown in the report, the frozen contract and the GUI. Counts, featureCounts summaries and L1 tables are byte-identical to P3_D5_NEG |
| REL_D4_HF | D4 real yeast, reverse | **`CONSISTENT`** (r 0.174, 0.179): no false warning despite antisense overlaps. Counts, summaries and L1 tables are byte-identical to P3_D4_HF |
| (subsample) | D5-R count-only BAM at 0.2% (394 strand-assignable fragments) | **`INCONCLUSIVE`**: not forced |
| REL_D4_SAL | D4, Salmon route | Strandedness `NOT_EVALUATED`, library types ISR recorded. `SALMON_QUANT` trace CPUs 6 equal `cmd_info.json` `--threads` 6; tool versions recorded |
| REL_D5_NEG, REL_D4_HF | — | Trace has `cpus`. `HISAT2_ALIGN` 4 CPUs equals the `hisat2 -p 4` of the executed `.command.sh`; `FEATURECOUNTS` 4, `STRAND_EVIDENCE` 2. Shown in the GUI provenance and the report |
| REL_D3_COUNTS | Phase 3 D3 HISAT2-route counts as raw-count L2 | **6 apeglm "line search" warnings recorded**, the same number Phase 3 found in nf-rna's work-directory stderr. Status `WARNINGS`, shrinkage `success`. Present in the delivered `contrast_metadata.json`, the provenance, the report (verbatim) and the GUI. `all_genes`, `significant`, `up`/`down` and the normalized counts are byte-identical to P3_D3_HF |
| REL_D1 | airway D1, L2 with GO/KEGG ORA and GO/KEGG GSEA | **All 7 L1 and 125 L2 TSV tables (DE, ORA, GSEA) are byte-identical to FIX-D1** (v1.6.1 corrections, `3cf0386`). `r_warnings` recorded as `NONE` |

**One FAIL in the acceptance table, kept and explained.** REL_D4_SAL `quant.sf` is not byte-identical to P3_D4_SAL.

- `NumReads` is identical, and so are `estimated_counts.csv` and every count.
- `EffectiveLength` and `TPM` differ, by at most 0.33% relative, through Salmon's fragment-length estimate (mean 165.897 vs 165.840).
- The Salmon command, `cmd_info.json`, parameters and inputs are identical. The only configuration difference is the new trace-fields line, which cannot affect a task.
- Three standalone Salmon 1.10.3 runs on identical trimmed reads at `--threads 6` also differ pairwise in `EffectiveLength` (105–117 rows), never in `NumReads` (`results/salmon_same_thread_repro.tsv`). This is Salmon's own run-to-run variation, which answers the Phase 3 open question: **Salmon is not bitwise reproducible even at the same thread count**.
- Effect: through the tximport length offsets, D4 L1 normalized counts differ by ≤ 1.03e-3 relative and VST by ≤ 6.4e-5. This is the F-2/N-3 limitation, not a consequence of this change.

## Effect on scientific results

**No analysis value is changed by F-1, F-2 or F-3.** Shown by:

- byte-identical counts, featureCounts summaries and L1 tables (HISAT2 route);
- identical Salmon `NumReads`/estimated counts;
- byte-identical DE tables on identical counts;
- byte-identical DE, ORA and GSEA tables for D1 (132 tables).

New outputs:

- `qc/strandedness/` (two summaries per sample, plus the verdict file);
- new keys in the handoff manifest, the frozen contract, `run_provenance.yaml`, `contrast_metadata.json`, `backend_summary.json` and the L2 provenance;
- new report and GUI sections;
- extra trace columns.

`STRAND_EVIDENCE` costs two extra featureCounts passes per sample (2 CPUs each).

## Not supported or not validated

- **Hard failure** on a strandedness contradiction is not enabled (see F-1).
- **The Salmon route** is not strand-checked by nf-rna (`NOT_EVALUATED`).
- **Validation data.** The strand limits were set from one real stranded library (yeast, reverse), unstranded human data and simulated F/R/U/SE data. Antisense-rich protocols, very low-depth samples, natural forward-stranded data and other organisms are not validated; such samples are expected to fall into `INCONCLUSIVE` rather than receive a false warning, but this is not demonstrated.
- **Thread counts are recorded, not controlled.** Results reproduce exactly only at the same thread counts for HISAT2. For Salmon, even the same thread count gives run-to-run variation in effective lengths, while its counts were identical here.
- **Runs before v1.6.1** show "not evaluated" / "not recorded" / "R warnings: not recorded"; their results are not changed.
- **L1 R warnings** (`vst`, PCA) are not captured; F-3 covers the L2 DESeq2/apeglm/ashr steps.

## Effect on v1.6.1 release readiness

- **F-1** has gone from a silent failure mode to a reported warning, so it is no longer a hidden risk.
- **F-2** and **F-3** are closed for the scope above.
- **No new release blocker was found.**
- **Phase 3 itself is not re-verdicted.** Its frozen FAILs (E-UP-4/5, SE-4, N-UP-8, FC-3 on the HISAT2 route, SA-1 D4 standalone, T-1, T-2) remain FAIL. Their root causes are tool behaviour and criterion issues, as the Phase 3 report states.
- The release decision remains with the maintainer.
