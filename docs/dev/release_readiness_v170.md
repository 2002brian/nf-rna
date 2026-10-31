# nf-rna v1.7.0 release readiness

Release: **v1.7.0 — Local Research Workspace & Publication Figure Export**,
release date 2026-10-31. Branch `feature/v1.7.0`. This review reuses the
Phase 7a/7b acceptance evidence and adds only version, packaging and GUI smoke
checks (§5); no analysis and no FASTQ pipeline was run.

**Verdict: ready for release, subject to the maintainer's decision and to a
green `Release gate` CI run on the pushed `main` (no v1.7.0 commit has run in
CI).** No release blocker is known (§7).

## 1. Git state and ancestry

| Item | State |
|---|---|
| `main`, `origin/main` (local tracking ref), `v1.6.1` | `7a92bdd` |
| `feature/v1.7.0` | `1385318` (7a) → `f04dce8` (7b) → finalization commit; linear on `v1.6.1`, no merges |
| Authorship `v1.6.1..HEAD` | Brian Wu as author and committer, no trailers (CLAUDE.md) |
| Working tree | clean after the finalization commit |
| Remote | not reachable from this session (no SSH key, no `gh`); `origin/main` is the last fetched value |

`main` fast-forwards to `feature/v1.7.0`.

## 2. Scope and what is unchanged

Changed since `v1.6.1`: `src/rnaseq/gui/` (new `datasets.py`, `workspace.py`,
`figure_export.py`; `pages.py`, `results.py`), `src/rnaseq/__init__.py`
(version), tests, docs, dev scripts.

```
git diff --stat v1.6.1 HEAD -- src ':!src/rnaseq/gui'   → src/rnaseq/__init__.py only (version)
git diff --stat v1.6.1 HEAD -- src/rnaseq/r workflow workflow/envs pyproject.toml environment.yml   → empty
```

Unchanged: DESeq2, tximport, shrinkage, GO/KEGG ORA, GSEA, Salmon, HISAT2,
featureCounts, thresholds, references/annotation, R analysis and plotting
scripts, Nextflow workflows, both Conda locks, nf-core/rnaseq 3.26.0, package
dependencies (`gui` extra: `streamlit<2,>=1.66`; Pillow and pandas come with
Streamlit). Existing formal runs were not modified (SHA-256 checks in §4, §5).

## 3. Version consistency

| Where | Value |
|---|---|
| `src/rnaseq/__init__.py` (single version source; Hatch metadata and `PIPELINE_VERSION`) | `1.7.0` |
| `CITATION.cff` | `version: "1.7.0"`, `date-released: 2026-10-31` |
| `CHANGELOG.md` | `## 1.7.0 — 2026-10-31` (first heading) |
| Wheel / sdist METADATA | `Version: 1.7.0` |
| `rnaseq --version` from each installed artifact; GUI sidebar | `1.7.0` |
| README_latest.md, README_zh-TW.md, docs/quickstart.md, docs/gui.md, docs/runtime.md | name `v1.7.0` as the latest release, install from `@v1.7.0`, CHANGELOG anchor `#170--2026-10-31` |
| README.md (historical v1.3.0 README) | pointer line names v1.7.0 |

If the tag is made on a later day, change the date in `CHANGELOG.md` and
`CITATION.cff` together (`tests/test_release_metadata.py`) and the
`#170--2026-10-31` anchors in the READMEs, quickstart and release notes.

## 4. Feature evidence (existing, Phase 7a/7b)

| Area | Evidence | Result |
|---|---|---|
| Local Dataset Browser, Project Workspace, Run History, QC display | `docs/dev/phase7a_report.md`, `docs/dev/screenshots/7a/browser_checks.json` | 27/27 browser checks; 4,048 run/project files SHA-256 identical before/after |
| PDF/TIFF downloads, all 73 figures of airway-gui-demo GUI-6B-L2 | `docs/dev/phase7b_report.md`, `/home/brian/nf-rna-validation/gui-7b/` | 73/73: PDF byte-identical and opens; 146 TIFFs correct format/RGB/LZW/dpi/pixel size; 300 dpi byte- or pixel-identical to pipeline files; 600 dpi vs Ghostscript ≤ 1.48/255, no dominant colour missing; 853 run files unchanged |
| PDF/TIFF in the browser | `docs/dev/screenshots/7b/browser_checks.json` | 10/10 |
| Full suite on `f04dce8` (locked runtime `0e721ff3…-200d98f8759eaf17`) | `docs/dev/phase7b_report.md` §6 | 1,133 passed, 2 skipped (environment), 0 failed |

## 5. Release checks run for this review

| Check | Result |
|---|---|
| `tests/test_release_metadata.py` | 2 passed |
| `tests/test_gui*.py`, `test_release_metadata.py`, `test_figure_settings.py` (finalization tree, incl. new poppler and reuse-hint tests) | 129 passed |
| Version-sensitive subset on the finalization tree (`-k "reuse or version or provenance or release or smoke"`, locked runtime) | 82 passed |
| `python -m build` (code of the finalization tree; release docs were completed afterwards, CI rebuilds) | `nf_rna-1.7.0-py3-none-any.whl` (1.5 MB), `nf_rna-1.7.0.tar.gz` (22.5 MB) |
| Wheel inventory | `rnaseq/gui/{__init__,app,config,datasets,figure_export,jobs,pages,results,services,workspace}.py`, `rnaseq/workflows/{main.nf,hisat2_featurecounts.nf,nextflow.config}`, `rnaseq/r/*.R`, bundled Arimo font and manifest; no `__pycache__` |
| sdist inventory | GUI sources, tests (incl. `test_gui_figure_export.py`), docs, workflow, Phase 7a/7b screenshots |
| `tests/release_smoke.py` on the wheel and on the sdist, each in a fresh venv | pass (both); `rnaseq --version` 1.7.0; `rnaseq gui --help` works without the `gui` extra |
| GUI from the installed wheel with `[gui]` (fresh venv, Streamlit 1.66.1, port 8702, acceptance `gui.yaml`) — figure downloads | 10/10 (`/home/brian/nf-rna-validation/release-170/figure_downloads_wheel.json`) |
| Same server — workspace smoke (projects, datasets, workspace, switching, qc, results) | 21/22; the one miss is the synthetic `failed-demo` fixture project, deliberately removed after Phase 7a |
| Same wheel, server started with `PATH` = venv only (no poppler) | 6/6: *Projects* notice with install command; TIFF 600 dpi disabled on Results and Reports with the reason; PDF, TIFF 300 dpi and the running-score TIFF still byte/pixel-identical; no traceback in the server log (`figure_downloads_wheel_nopoppler.json`, screenshots `docs/dev/screenshots/7b/F5`, `F6`) |
| GUI-6B-L2 run files during every browser check | 853 files, SHA-256 unchanged |

The user's own GUI (port 8501) and `~/.config/nf-rna/gui.yaml` were not touched.

## 6. Changes made during finalization

* Version 1.7.0, CHANGELOG, CITATION, READMEs, quickstart, gui, runtime docs,
  release notes (`docs/dev/release_notes_v1.7.0.md`), this review, checklist.
* **Poppler**: `figure_export.poppler_status()`; the *Projects* page shows whether
  TIFF 600 dpi is available and, if not, the install command; the unavailable
  reason names the command. Install documented in `docs/gui.md`, README_latest,
  README_zh-TW, quickstart and the release notes. Tests: unit test with poppler
  hidden; AppTest of *Results* with poppler hidden.
* **Reuse hint (GUI, display only)**: found while checking compatibility. The
  frozen upstream identity that `--reuse-upstream` compares at run time
  includes the nf-rna version (`service._upstream_identity`,
  `input_manifest.yaml` `pipeline.version`); verified with a REL_D5_NEG manifest
  (1.6.1 vs 1.7.0 → not equal; same version → equal). The 7a hint only ran the
  cheap preflight and would have offered v1.6.1 runs, which the run then refuses
  after it has started. The hint now lists runs of another nf-rna version as not
  reusable with the reason. The reuse rule itself is unchanged. Test added.

## 7. Release blockers

None known. Before tagging, the `Release gate` CI must pass on the pushed
`main` (checklist step 4).

## 8. Known limitations (CHANGELOG, release notes)

* PDF vs raster text layout differs (1–3 pt; up to ≈10 pt on dotplots); kept,
  R plotting engine not changed (`phase7b_report.md` §4).
* `--reuse-upstream` refuses v1.6.1 source runs (version in the upstream
  identity, as in every release); rerun the upstream stage once.
* TIFF 600 dpi needs poppler (system package, not in `environment.yml`).
* No "use in place" FASTQ mode; count-matrix origin not recorded.
* No TIFF > 600 dpi, SVG or figure editing; `elsevier` profile TIFF paths tested
  with fixtures only.
* SSH-tunnel access not exercised; FAILED state verified with a fixture; single
  user on 127.0.0.1.
* All scientific limitations and frozen validation verdicts of v1.6.1 stand
  (`release_readiness_v161.md`).
