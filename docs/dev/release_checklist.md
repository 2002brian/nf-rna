# nf-rna v1.7.0 release checklist

State: v1.7.0 is finalized on `feature/v1.7.0` with release date
**2026-10-31** (CHANGELOG `## 1.7.0 — 2026-10-31`, CITATION.cff
`date-released: 2026-10-31`). The finalize commit is the head of
`feature/v1.7.0`; its parent is `f04dce8`. `main` and `origin/main` are at
`v1.6.1` (`7a92bdd`), an ancestor, so `main` fast-forwards. Nothing for v1.7.0
has been pushed, tagged or released. The v1.6.1 checklist is in the `v1.6.1`
tag.

Readiness, evidence and limitations: `docs/dev/release_readiness_v170.md`.
The steps are run by the maintainer, in order; stop at the first failure. If
the tag is made on a later day than 2026-10-31, first change the date in
`CHANGELOG.md`, `CITATION.cff`, the `#170--2026-10-31` anchors (README_latest.md,
README_zh-TW.md, docs/quickstart.md, release notes) and this file, in one commit.

## Scope

- Since `v1.6.1`: Phase 7a local research workspace (`1385318`), Phase 7b
  PDF/TIFF figure downloads (`f04dce8`), release finalization. GUI only.
- Unchanged: `src/rnaseq/r`, `workflow`, both Conda locks, `pyproject.toml`
  dependencies, `environment.yml`, nf-core/rnaseq 3.26.0, thresholds.

## Release evidence

| Check | Result | Where |
|---|---|---|
| Phase 7a browser acceptance | 27/27; 4,048 files unchanged | `docs/dev/phase7a_report.md` |
| Phase 7b figure export, all 73 figures | 73/73; browser 10/10; run unchanged | `docs/dev/phase7b_report.md` |
| Full suite on `f04dce8`, locked runtime `0e721ff3…` | 1,133 passed, 2 skipped, 0 failed | `docs/dev/phase7b_report.md` §6 |
| Build, inventories, `release_smoke.py` (wheel, sdist, fresh venvs) | pass; `rnaseq --version` 1.7.0 | `release_readiness_v170.md` §5 |
| GUI from installed wheel: downloads / workspace / no poppler | 10/10, 21/22 (removed fixture), 6/6 | `release_readiness_v170.md` §5 |

## Release blockers

None known. CI on the pushed `main` must pass before tagging (step 4).

## 1. Review

```bash
cd /home/brian/projects/nf-rna
git status                                          # clean
git log --oneline v1.6.1..feature/v1.7.0            # 3 commits, linear
git merge-base --is-ancestor v1.6.1 feature/v1.7.0 && echo ancestor-ok
git diff --stat v1.6.1 feature/v1.7.0 -- src/rnaseq/r workflow   # expect no output
git log --format='%an|%cn|%(trailers:only)' v1.6.1..feature/v1.7.0 | sort -u   # Brian Wu|Brian Wu|
```

## 2. Fast-forward main

```bash
git switch main
git merge --ff-only feature/v1.7.0
git log --oneline -1                                # the finalize commit
python -m pytest -q tests/test_release_metadata.py
```

## 3. Push main

```bash
git push origin main
```

## 4. Wait for CI

```bash
gh run list --workflow "Release gate" --branch main --limit 3
gh run watch "$(gh run list --workflow 'Release gate' --branch main --limit 1 --json databaseId --jq '.[0].databaseId')" --exit-status
gh run list --workflow "Release gate" --branch main --limit 1 --json headSha --jq '.[0].headSha'   # must equal: git rev-parse main
```

Do not tag until this run is green on the finalize commit.

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
gh release view v1.7.0 --web
```

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
