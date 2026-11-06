# v1.7.0 release cleanup and repository hygiene (2026-11-06)

Inventory taken on `feature/v1.7.0` at `9a983db`. Scope: repository hygiene
only; no algorithm, default, schema, CLI, workflow, lock or test change.

## 1. Tracked files

445 tracked files. `git ls-files -ci --exclude-standard` is empty (nothing is
tracked against `.gitignore`). No tracked cache, build output, log, backup or
editor file was found.

| Area | Files | Verdict |
|---|---|---|
| `src/`, `workflow/`, `tests/` (incl. `tests/science/` oracles and fixtures) | 187 | keep: production code, Nextflow, scientific gate |
| `ci/`, `scripts/ci/`, `.github/workflows/` | 6 | keep: required by `test_ci_gate.py` and the release/science gates |
| `workflow/envs/locks/` (two locks + `SHA256SUMS`), `environment*.yml` | 6 | keep: runtime identity |
| `src/rnaseq/fonts/arimo/` (4 TTF, `FONT_MANIFEST.json`, `OFL.txt`) | 6 | keep: bundled fonts and integrity manifest |
| `examples/` (smoke-test FASTQ/reference, ~9.5 MB) | 42 | keep: used by tests, whitelisted in `.gitignore` |
| `validation/corrections_v161/`, `validation/reliability_v161/` | 33 | keep: v1.6.1 scientific evidence |
| `tools/figure_review/`, `scripts/dev/` | 29 | keep: figure-review tooling (tested) and acceptance scripts cited by reports |
| `reference_tools/validate_gtf_derived_transcriptome.py` | 1 | keep: maintained user utility (imports `rnaseq.references`); not referenced elsewhere — review whether to document it |
| `README.md` (frozen v1.3.0) / `README_latest.md` | 2 | keep both: deliberate (`DECISIONS.md`); `pyproject.toml` `readme` is `README_latest.md` |
| `docs/dev/` phase reports, release notes, readiness, audit (without screenshots) | 24 | keep: historical record; excluded from the sdist |
| `docs/dev/screenshots/` (104 files, 14.7 MB) | 104 | human review: acceptance evidence cited by the reports by name (some only by directory); excluded from the sdist |
| `CLAUDE.md` | 1 | keep: repository policy; not in the sdist |

No file met the "confirmed unnecessary and unused" bar, so no tracked file was
removed.

## 2. `.gitignore`

Added common local tool directories: `.venv/`, `venv/`, `.mypy_cache/`,
`.ruff_cache/`, `.ipynb_checkpoints/`, `.idea/`, `.vscode/`. None matches a
tracked file; `test_ci_gate.py` (gate files not ignored) still passes.

## 3. Branches and worktrees

Merged into `feature/v1.7.0` (ancestor check):

| Branch | Worktree | Worktree state |
|---|---|---|
| `fix/v1.7.0-audit-finalize` (`5cba156`) | `../nf-rna-audit-check` | clean; ignored only `__pycache__/`, regenerated `examples/simple_two_group/planning/` |
| `fix/v1.7.0-hisat2-artifact-reuse` (`1de489e`) | `../nf-rna-hisat2-reuse` | same |
| `fix/v1.7.0-l2-figure-standardization` (`ddb6807`) | `../nf-rna-l2fig` | clean; ignored only `__pycache__/` |
| `ci/scientific-release-gates` (`0246d15`, remote only) | — | merged by PR #10 |

Other merged local branches: `develop`, `feature/l2-enrichment-default`,
`feature/plotting-standard`, `feature/v1.5.0`, `feature/v1.6.0`,
`fix/v1.6.1-reliability`, `fix/v1.6.1-scientific`.

Not merged (keep): `validation/v1.6.0` (3 commits, Phase 2 evidence),
`validation/v1.6.1-phase3` (6 commits, Phase 3 report and scripts).

Remote (GitHub, no branch is protected and no ruleset is visible publicly):

| Remote branch | State |
|---|---|
| `ci/scientific-release-gates`, `fix/v1.7.0-audit-remediation`, `develop`, `feature/l2-enrichment-default`, `feature/native-reference-delivery`, `feature/v1.5.0`, `migration/nfcore-conda` | merged (ancestor) |
| `feature/generalized-design-covariates` | patch-equivalent (`git cherry`: 0 unmerged; PR #6) |
| `hotfix/v1.2.1` | not merged; holds tag `v1.2.1` (`62e2ce4`), which keeps the commits if the branch goes |
| `release/v1.2.0` | 2 commits not merged (version bump; superseded) |
| `audit/cloud-2026-11-04` | 1 commit not merged (audit report and `tests/audit/` tests) |

Deletion of any branch, remote branch or worktree is left to the maintainer.

## 4. Remaining release items

* `docs/dev/release_checklist.md` "State" still describes the 2026-10-31
  finalize (`ac6471d`, "nothing pushed"); `feature/v1.7.0` has since gained the
  audit, HISAT2 reuse, L2 figure and CI-gate merges and is on `origin`.
* Release date is `2026-10-31` (`CHANGELOG.md`, `CITATION.cff`, anchors). The
  checklist requires changing it in one commit if the tag is made later.
* The checklist's evidence table (1,133 tests on `f04dce8`) predates the later
  merges; the release gate must pass on the final `main`.
