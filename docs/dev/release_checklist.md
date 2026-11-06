# nf-rna v1.7.0 release checklist

State (2026-11-06): v1.7.0 is finalized on `chore/v1.7.0-release-finalization`,
branched from `feature/v1.7.0` at `c81d789`, with release date **2026-11-06**
(CHANGELOG `## 1.7.0 — 2026-11-06`, CITATION.cff `date-released: 2026-11-06`,
`#170--2026-11-06` anchors in README.md, README_zh-TW.md,
docs/quickstart.md and the release notes). It replaces the 2026-10-31 date of
the first finalization (`ac6471d`), which was never tagged. `main` and
`origin/main` are at `v1.6.1` (`7a92bdd`), an ancestor, so `main`
fast-forwards. Nothing for v1.7.0 has been tagged or released; `origin/feature/v1.7.0`
is at `9a983db`.

Readiness, evidence and limitations: `docs/dev/release_readiness_v170.md`
(§10 is the final state). The steps are run by the maintainer, in order; stop
at the first failure. If the tag is made on a later day than 2026-11-06, first
change the date in `CHANGELOG.md`, `CITATION.cff`, the `#170--2026-11-06`
anchors and this file, in one commit (`tests/test_release_metadata.py`
requires CHANGELOG and CITATION to agree).

## Scope

- Since `v1.6.1`: Phase 7a workspace (`1385318`), Phase 7b PDF/TIFF downloads
  (`f04dce8`), GUI hardening and UX, L1 figure standardization, 2026-11-04 audit
  remediation (`5cba156`), HISAT2 artifact lifecycle and reuse (`1de489e`), L2
  figure standardization (`ddb6807`, merge `45787df`), scientific CI gates
  (PR #10, `9a983db`), repository cleanup (`e1343f6`, merge `c81d789`).
- Changed outside the GUI: figure code and figure records in `src/rnaseq/r`;
  `publishDir` patterns and a sample-ID guard in `workflow/hisat2_featurecounts.nf`;
  validation, atomic state writes and upstream reuse in the control plane.
- Unchanged: statistics and thresholds, `workflow/main.nf`, both Conda locks
  (`workflow/envs/locks/SHA256SUMS`), package dependencies, `environment*.yml`,
  nf-core/rnaseq 3.26.0.

## Gates

| Gate | State | Evidence |
|---|---|---|
| Version 1.7.0 in `__init__.py`, CITATION, CHANGELOG, wheel/sdist METADATA | done | `tests/test_release_metadata.py`; §10 of the readiness report |
| Documentation matches the integrated changes (CHANGELOG, release notes, READMEs, quickstart) | done | this finalization commit |
| Known limitations recorded (interrupted reuse, KEGG snapshot, FASTQ identity cache, low-depth VST, GUI) | done | CHANGELOG 1.7.0 *Known limitations*, release notes |
| Scientific CI fixtures, oracles, gate scripts, locks, fonts tracked and unchanged | done | `tests/test_ci_gate.py`; `git diff 9a983db -- tests/science ci scripts/ci workflow/envs src/rnaseq/fonts` empty |
| Local release tests (metadata, packaging, CLI, CI gate, workflow assets) | done | readiness §10 |
| Wheel and sdist build reproducibly; `release_smoke.py` in fresh venvs | done | readiness §10 |
| `Scientific validation` on PR #10 (`0246d15`) | **pass** | GitHub run 42043814911 |
| `Release gate` on PR #10 (`0246d15`) | **pass** | GitHub run 42043815085 |
| `Scientific validation` on the pushed `feature/v1.7.0` head | **outstanding** | push of `9a983db` (run 42145251344) was in progress at finalization; the final head must be pushed and pass |
| `Release gate` on the pushed `main` (final commit) | **outstanding** | step 4 |
| Tag, GitHub Release, tagged-install smoke | **outstanding** | steps 5–7 |

Release blockers: none known. Outstanding gates are the CI runs on the final
commit and the publication steps below.

## 1. Review

```bash
cd /home/brian/projects/nf-rna
git status                                                     # clean
git switch feature/v1.7.0
git merge --ff-only chore/v1.7.0-release-finalization
git merge-base --is-ancestor v1.6.1 feature/v1.7.0 && echo ancestor-ok
git diff --stat v1.6.1 feature/v1.7.0 -- workflow/envs workflow/main.nf   # expect no output
git log --format='%an|%cn' v1.6.1..feature/v1.7.0 | sort | uniq -c
#   Brian Wu|Brian Wu for every commit except the GitHub merge of PR #10 (2002brian|GitHub)
git log --format='%(trailers:key=Co-Authored-By)' v1.6.1..feature/v1.7.0 | grep -c . # 0
```

## 2. Push the release branch and wait for Scientific validation

```bash
git push origin feature/v1.7.0
```

Wait until `Scientific validation` for the pushed head is green
(GitHub → Actions, or `gh run watch` where `gh` is installed).

## 3. Fast-forward and push main

```bash
git switch main
git merge --ff-only feature/v1.7.0
python -m pytest -q tests/test_release_metadata.py
git push origin main
```

## 4. Wait for CI on main

```bash
gh run list --workflow "Release gate" --branch main --limit 1 --json headSha,conclusion
gh run list --workflow "Scientific validation" --branch main --limit 1 --json headSha,conclusion
git rev-parse main                                             # both headSha must equal this
```

Do not tag until both runs are green on this commit.

## 5. Tag

```bash
git tag -a v1.7.0 -m "nf-rna v1.7.0" main
git show v1.7.0 --no-patch
git push origin v1.7.0
```

## 6. GitHub Release

```bash
gh release create v1.7.0 --title "nf-rna v1.7.0 — Local Research Workspace & Publication Figure Export" \
  --notes-file docs/dev/release_notes_v1.7.0.md --verify-tag --latest
```

Without `gh`: GitHub → Releases → *Draft a new release*, tag `v1.7.0`, the same
title, body = `docs/dev/release_notes_v1.7.0.md`, *Set as the latest release*.

## 7. Tagged-install smoke test

```bash
conda create --yes --name nf-rna-170 python=3.11
conda run -n nf-rna-170 python -m pip install "nf-rna[gui] @ git+https://github.com/2002brian/nf-rna.git@v1.7.0"
conda run -n nf-rna-170 rnaseq --version      # 1.7.0
conda run -n nf-rna-170 python -c "import importlib.metadata as m; print(m.distribution('nf-rna').read_text('direct_url.json'))"   # vcs_info.commit_id = git rev-parse v1.7.0^{commit}
conda run -n nf-rna-170 rnaseq doctor         # Overall: READY
conda env remove --yes --name nf-rna-170
```

## Not to be done

- Do not move, recreate or delete existing tags (`v1.3.0` … `v1.6.1`).
- Do not rewrite history to change authorship; commits are authored by Brian Wu with no AI attribution (CLAUDE.md).
- Do not relabel the frozen Phase 2 or Phase 3 FAIL verdicts of v1.6.x as PASS.
