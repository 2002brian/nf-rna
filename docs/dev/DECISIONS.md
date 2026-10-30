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

## v1.5.0 release prep (2026-10-30)

| Decision | Reason |
|---|---|
| A legacy cutoff field is reported only when it is written in project.yaml (pydantic `model_fields_set`), in one `legacy_enrichment_cutoff` warning naming every such field. | Defaults are not the user's choice, and the wizard never writes these fields, so new projects do not warn. |
| `kegg_snapshot_not_pinned` is a warning only when KEGG ORA or GSEA is selected. It names the snapshot that `kegg_snapshot.R` would choose (`latest_kegg_snapshot` mirrors `nf_rna_latest_snapshot`), or says that none exists yet. A pinned snapshot missing locally is the error `kegg_snapshot_missing`. | Validation can only predict the choice; the backend records the snapshot actually used. A missing pinned snapshot would fail the run preflight anyway. |
| These two warnings are repeated in `logs/rnaseq.log` at run start (`ENRICHMENT_WARNING_CODES`). | They concern the reproducibility of that specific run. |
| Report: simplified GO ORA dotplots only for the `significant` foreground, and simplified GO GSEA dotplots per ontology. Running-score plots: ranks 1–2 per NES sign, per GO ontology and for KEGG. Missing figures are shown as "unavailable — reason", using the backend's recorded reason when there is one. | Keeps the report readable; every figure remains in the delivery `figures/`. |
| The GO GSEA `on.exit` cleanup removes only objects that still exist (`rm(list = intersect(..., ls()))`). | The success path already removed them; the warnings were noise. Tables are byte-identical. |
| Two preflight tests set `XDG_CACHE_HOME` to an empty temporary directory. | They assumed no KEGG snapshot exists, which now depends on the machine. |

## v1.5.1 maintenance release (2026-10-30)

| Decision | Reason |
|---|---|
| Ship the CI test-gating fix `f66f96e` as v1.5.1 instead of moving the published `v1.5.0` tag. | `v1.5.0` (`c283d31`) is public; moving it would change what existing installs from the tag resolve to. The fix touches tests only, so v1.5.0 stays a valid release. |
| v1.5.1 changes only the version, release metadata and docs; no source, workflow, lock or R script. | A maintenance release must not change scientific results; CHANGELOG states that v1.5.0 results need not be recomputed. |

## v1.5.0 finalization (2026-10-30)

| Decision | Reason |
|---|---|
| Release date 2026-10-30 in CHANGELOG and CITATION.cff, set in the finalization commit before tagging. | Same practice as `f69058e` (v1.4.1); the release-metadata guard requires the CHANGELOG heading and `date-released` to agree. |
| osx-arm64 is documented as an unqualified candidate lock, not a release blocker. | The only supported runtime is Linux x86-64/WSL2 (`docs/runtime.md`, `execution_backend()`); macOS would need its own native qualification in a release that targets it. |
| `compare_runs.py --expected-changes` with a versioned JSON (`expectations/v1.4.1-to-v1.5.0.json`), instead of widening `--require-identical-except`. | Every allowed difference is listed with a reason (changed kinds, new figures by ID pattern and kind, added formats). Unlisted kinds must stay identical (or differ only by a declared added format with identical pixels, PDF, TIFF and source data), and removed figures always fail. |
| `examples/gsea_two_group` drops the retired cutoff fields and uses `kegg_rest_snapshot`, with a commented pin line instead of a real snapshot ID. | Snapshot IDs are machine-local; a hard-coded ID would make the example fail validation on other machines. |

## v1.6.0 Phase 6a — GUI core (2026-10-30)

| Decision | Reason |
|---|---|
| The GUI is the package `rnaseq.gui`; only `pages.py` and `app.py` import Streamlit. The `gui` extra is `streamlit>=1.40,<2`; Playwright is in a separate `gui-dev` extra, not in `dev`. | The service modules and their tests run without the extra; CI installs `.[dev]` and needs neither Streamlit nor Playwright. |
| Wizard choices, defaults and `reference_options` moved to `rnaseq.wizard`; the `rnaseq design` transaction to `design_setup` (`design_blockers`, `residual_df_blocker`, `resolve_design_enrichment`, `apply_design`, `DesignRejected`); FASTQ/count import to `project` (`import_fastq_samplesheet`, `import_count_matrix`, `fastq_samplesheet_text`). The CLI calls the same functions; its prompts and messages are unchanged. | One engine: the GUI must not copy wizard or design logic. |
| `wizard.resource_suggestion` takes the capacity as an argument. | Existing CLI tests patch `rnaseq.cli.detect_local_resource_capacity`; the CLI passes its own value through. |
| The two design-rollback tests now patch `rnaseq.validators.validate_project` and `rnaseq.planner.generate_plan`. | The rollback moved from `cli.design_command` to `design_setup.apply_design`; the tests check the same behaviour at its new location. |
| The GUI creates projects as scaffolds (`create_project(..., scaffold=True)`), then imports data on the Input page. | The interactive wizard is scaffold-first; strict validation happens at Validate and Run. |
| FASTQ files are copied into `input/fastq/` (with the samplesheet kept as `planning/imported_fastq_samplesheet.csv`), never linked. Importing into a non-empty `input/fastq/` is refused. | Same as `rnaseq new --fastq-samplesheet`; the validator requires `input.path` inside the project. |
| *Save metadata* refuses only structural errors (missing/duplicate `sample_id`, blank/duplicate column, malformed row); empty design values are saved and listed. Sample agreement with the input is reported after saving. | Clarification §2: incomplete work can be saved; full validation is Validate. |
| Contrasts in the GUI are any chosen numerator/denominator pairs of the comparison factor, with the `rnaseq design` ID rule `<numerator>_vs_<denominator>`; the condition order frozen in project.yaml is the first-appearance order (or the already frozen order). | Same files as `rnaseq design`; contrasts.csv is the authority on direction. |
| Only the designs the wizard creates are offered (`two_group`, `paired_two_group` with a pairing column, `multi_group`, optional `batch`). Other covariates remain the CLI's `--covariate`/`--continuous-covariate` import route. | The GUI must not invent formula logic; the wizard covers the same set. |
| *Pin to the current snapshot* writes `annotation.enrichment.kegg.snapshot` with `design_setup.pin_kegg_snapshot`, which replaces only the `annotation` block and validates the result. | Keeps the rest of project.yaml as written. |
| *Validate* on the Run page also refreshes the planning artifacts when the project is valid. | `rnaseq run` refuses a stale plan, and metadata or pin edits make it stale; this is what `rnaseq plan` does. |
| A recorded `CREATED`/`RUNNING` run whose process is gone is labelled INCOMPLETE in the GUI, from `run_status.effective_state` (which says INTERRUPTED, stale). Nothing is written to the run. | The existing mechanism detects it; the label distinguishes it from a run that recorded its own interruption. |
| Runs start as a transient `systemd --user` unit when the GUI runs under systemd (`INVOCATION_ID`) or lingering is enabled; otherwise as a detached process (`start_new_session=True`). `gui.yaml` `launcher:` can force either. | Clarification §1. Under a service, children share its cgroup and stop with it. Without lingering the user manager stops at logout, while a detached process does not. On this server lingering is off, so `auto` = process. |
| The transient unit receives PATH, HOME, locale, CONDARC, JAVA_* and RNASEQ_/NXF_/NF_RNA_/CONDA_/XDG_ variables via `--setenv`; stdout/stderr go to the GUI launcher log (`append:`). | A transient unit starts from the user manager's environment, not the GUI's. |
| The active-run lock is `~/.cache/nf-rna/gui/active_run.lock`, created with `O_CREAT|O_EXCL`. It is stale when its process identity is not alive (`run_status.process_liveness`: pid, boot ID, start ticks) or, without a process yet, after 120 s. | Atomic single submission; process identity guards against PID reuse. |
| A thread waits on each detached child. | An unreaped child is a zombie, and `kill(pid, 0)` would report it alive. |
| Launch records (`launches/<id>.json`) and launcher logs (`logs/<id>.log`) live in `~/.cache/nf-rna/gui/`. A launch is matched to its run by the pid and start ticks in `logs/rnaseq.process.json`. | The run directory exists only after `rnaseq run` passes its preflight; analysis state stays in the run. |
| The runtime check (`downstream_runtime.planned_runtime_prefix`) builds or materialises the wheel exactly as a run does and checks whether that prefix exists; it is cached per source revision and skipped while a run is active. | The prefix name includes the wheel SHA-256 (deterministic: two builds gave the same hash); a run building the same wheel at the same time is avoided. |
| Allowed roots under `/media`, `/mnt` and `/run/media` must be on a mounted filesystem (nearest mount point is not `/`). | Clarification §4: an unmounted disk must not look like an empty local folder. |
| The default allowed roots and CONDARC are the paths given in the specification, with the home-directory ones written as `~/nf-rna-work/...`; `gui.yaml` is written once with them expanded. | They are this server's data locations; other installations edit gui.yaml. The shipped package must not contain `/home/brian` (local-path guard in `test_figure_review_tools`). |
| `rnaseq gui` replaces itself with `python -m streamlit run` (`os.execv`), bound to 127.0.0.1, with the file watcher off and the toolbar in viewer mode. | A service manager then supervises the server directly; no "Deploy" button, no reloads from source edits. |
| No automatic polling; *Refresh* re-reads the run files. | Specification: manual refresh only. |
| Browser acceptance is `scripts/dev/gui_browser_acceptance.py` (Playwright driving the installed Google Chrome). The grid paste uses the real clipboard and Ctrl+V. | AppTest cannot drive `st.data_editor`; a synthetic paste event is ignored by the grid. |

## v1.6.0 Phase 6b — results, reports, release prep (2026-10-30)

| Decision | Reason |
|---|---|
| The Results and Reports pages read `downstream/` of the selected run (the R outputs) and the delivered `report_*.html`; the reading code is `rnaseq.gui.results` (no Streamlit). | `downstream/` is what the delivery copies; one place for both pages, testable without the GUI extra. |
| Tables are read with pandas as text (`dtype=str`, no NA conversion). Filtering and sorting convert a column with `to_numeric` only to choose and order rows. | The displayed values are the TSV strings, so they equal the files exactly (no float formatting). |
| *Only significant* in the DEG table is membership of `significant.tsv`; the enrichment default view is the backend's significant table (`GO_<ont>.tsv`, `kegg.tsv`, `significant.tsv`), *Show all* is the all-terms table. | Significance is never recomputed in the GUI. |
| Enrichment output whose backend summary lacks `significance_semantics: "v1.5"` is labelled legacy; its rule is stated as the pre-v1.5 rule and its member table falls back to the legacy `*_gene_members.tsv`/`*core_members.tsv`. | Spec: legacy results must not be presented with v1.5 definitions. |
| Missing outputs are shown as "unavailable — reason", the reason taken from the backend summary's `figures` record (`workflow_support._figure_reason`) or a fixed explanation. | Same rule and helper as `report.html`. |
| Pathway members: clicking a row of the terms table (`st.dataframe` row selection) selects the pathway; a selectbox lists the pathways of the member table. Members come from `*_members.tsv`, which covers significant pathways only (a non-significant pathway says so). | v1.5 member tables are written for significant pathways only (Phase 2 decision). |
| The run's existing delivery has no archive, so *Prepare delivery package* zips the delivery folder after `service.verify_delivery_manifest` passes, into `~/.cache/nf-rna/gui/packages/` (written file by file, reused while newer than the manifest). | No existing delivery function writes an archive, and existing runs must not be modified; the zip contains the verified delivery files unchanged. |
| `download_limit_mb` (gui.yaml, default 500): larger single files and delivery packages show the server path and an `scp` command instead of a download. Downloads use a callable `data`, so a file is read only when its button is clicked. | Spec; avoids reading every offered file on each rerun. |
| The `gui` extra requires `streamlit>=1.66,<2`. | The pages use callable download data, dataframe row selection and `width="stretch"`, which older 1.x releases lack; 1.66.1 is the tested version. |
| Reports page lists single files per group in a selectbox with one download button, not one button per file. | Hundreds of figure and enrichment files per L2 run. |
| `report.html` adds Samples (staged metadata.csv), Reference and quantification (from the frozen project.yaml and source type), the DEG rule and shrinkage per contrast (from `contrast_metadata.json`, "not recorded" for pre-v1.5 L2), and Software versions and provenance (execution manifest runtime identity; package versions from the L1/L2/enrichment records). Both the L1-only and L2 reports. | Spec list of report sections; every value comes from a file the run already has, and no table or figure changes (`enrichment_status.tsv` unchanged). |
| The report's raw-count "Reference and quantification" says that nf-rna performed no quantification and used no reference, with the provider/method recorded in project.yaml. | Raw-count projects have no reference; the section must not invent one. |
| Version 1.6.0 with `1.6.0 — Unreleased` and no `date-released`; v1.5.1 is not tagged and its CHANGELOG entry stays as a historical section that 1.6.0 links to. | Release-metadata guard; the finalize commit sets the date. |
| GUI install instructions use `"nf-rna[gui] @ git+https://github.com/2002brian/nf-rna.git@vX.Y.Z"`, not `pip install 'nf-rna[gui]'`. | nf-rna is not published on PyPI; the 6a text named a package index install that does not exist. |
| Phase 6b browser acceptance is `scripts/dev/gui_results_acceptance.py`; the unavailable and L1-only views were checked on truncated copies of the GUI-L2-002 run in `~/nf-rna-work/testdata/gui6b-fixtures/`, deleted afterwards. | Spec allows a fixture or a truncated project copy; no real run was modified. |
