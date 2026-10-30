# nf-rna v1.6.0 release checklist

State: v1.6.0 is finalized on `feature/v1.6.0` with release date **2026-10-30**
(CHANGELOG `## 1.6.0 — 2026-10-30`, CITATION.cff `date-released: 2026-10-30`).
The finalize commit is the head of `feature/v1.6.0`; its parent is `41d394c`.
Local `main` is at `9a7352b` (v1.5.1 maintenance commit, never pushed),
`origin/main` at `f66f96e`; both are ancestors of `feature/v1.6.0`, so `main`
fast-forwards. **v1.5.1 is not tagged**: its changes ship in v1.6.0 and its
CHANGELOG section says so. Nothing for v1.6.0 has been pushed, tagged or
released.

The remaining steps are run by the maintainer, in order; stop at the first
failure.

## Scope

- Since `v1.5.0` (`c283d31`): the v1.5.1 test-gating fix (`f66f96e`) and its
  release metadata (`9a7352b`); the browser GUI (`rnaseq.gui`, `rnaseq gui`,
  extras `gui` = `streamlit>=1.66,<2` and `gui-dev` = `gui` + Playwright);
  wizard, `rnaseq design` and FASTQ/count import moved into shared functions
  used by the CLI and the GUI; new `report.html` sections; 1.6.0 metadata.
- Unchanged since `v1.5.0`: `src/rnaseq/r/`, `workflow/`, `src/rnaseq/workflows/`,
  `src/rnaseq/templates/`, both Conda locks (`git diff --stat v1.5.0 HEAD --
  src/rnaseq/r workflow src/rnaseq/workflows src/rnaseq/templates` is empty).
  DESeq2, GO/KEGG ORA, GSEA and the "v1.5" significance semantics are unchanged.

## Release evidence

| Check | Result | Where |
|---|---|---|
| Airway L2 rerun from the GUI, 1.6.0 build `80f1d08` | SUCCESS; 135/135 tables byte-identical to `airway-e2e`; 142/142 PNG and 35/35 TIFF identical, 73/73 PDF identical except the creation date | `docs/dev/phase6b_report.md` |
| GUI acceptance (Phase 6a, 6b) | all items passed: project → metadata → contrasts → run (`airway-gui-demo`), restart, single active run, kill → INCOMPLETE; DEG table and KEGG members equal to the TSVs; downloads SHA-256 equal; report offline with all images | `docs/dev/phase6a_report.md`, `docs/dev/phase6b_report.md` |
| Full suite with the locked R runtime | 1080 passed, 2 skipped (`80f1d08`); the 2 report tests added in `41d394c` pass | `phase6b_report.md` |
| Suite without the GUI extra (no Streamlit/pandas, system R 4.3.3) | 1018 passed, 50 skipped, 4 failed, 10 errors. The 14 are `test_figure_presets` (4) and `test_figure_migration`/`test_figure_remediation` setup errors (10); `tests/test_figure_*.py` on `main` (`9a7352b`) in the same venv: 4 failed, 67 passed, 10 errors, the **identical** 14 tests as HEAD. Cause: the system R lacks the locked figure packages/fonts. Not a regression. | `/home/brian/nf-rna-validation/gui-6b/pytest-main-baseline.log`, `pytest-head-figs.log` |
| Build (finalize tree) | `nf_rna-1.6.0-py3-none-any.whl`, `nf_rna-1.6.0.tar.gz`; METADATA `Version: 1.6.0`, `Provides-Extra: gui` → `streamlit<2,>=1.66`; wheel contains `rnaseq/gui/*` and the bundled workflows | local `python -m build` |
| `tests/release_smoke.py` on the wheel and the sdist (fresh venvs) | pass | local |
| Wheel without extra | `rnaseq --version` 1.6.0; `rnaseq gui` exits 1 with the install command | local |
| Wheel with `[gui]` | Streamlit 1.66.1; `rnaseq gui --port 8613` serves on 127.0.0.1 only, `/_stcore/health` ok, pages 200 | local |
| Release metadata guard | `tests/test_release_metadata.py` passes | finalize commit |

## Release blockers

None known. Before tagging, CI on the pushed `main` must pass (step 4). CI
never ran on the v1.6.0 commits (they were never pushed).

## Known limitations (documented in CHANGELOG, release notes, docs/gui.md)

- GUI is single-user and local (127.0.0.1), one active GUI run.
- Verified end to end on raw-count projects only; FASTQ recognition/import was
  verified, but a FASTQ run was not started or browsed from the GUI.
- The GUI as the example `systemd --user` service was not tested; the
  `systemd` launcher was tested from a terminal-started GUI. Lingering is not
  enabled on this server.
- `process` launcher: SIGKILL of `rnaseq` leaves its Nextflow child running.
- Legacy (pre-v1.5) enrichment display tested with fixtures only.
- FASTQ copy blocks the page; other covariates, `retry` and `clean` are CLI-only;
  GUI delivery zips accumulate in `~/.cache/nf-rna/gui/packages/`.
- The sdist is ~19.7 MB because it includes the 6a/6b screenshots under
  `docs/dev/screenshots/` (nf-rna is not distributed through PyPI).

## 1. Review

```bash
cd /home/brian/projects/nf-rna
git status                                    # clean
git log --oneline -8 feature/v1.6.0           # finalize commit, 41d394c, 80f1d08, ed4c4ca, 7cdb329, 9a7352b, f66f96e
git diff --stat v1.5.0 feature/v1.6.0 -- src/rnaseq/r workflow src/rnaseq/workflows src/rnaseq/templates   # expect no output
git show feature/v1.6.0 --stat
```

## 2. Fast-forward main

```bash
git switch main
git merge --ff-only feature/v1.6.0
git log --oneline -1                          # the finalize commit
python -m pytest -q tests/test_release_metadata.py
```

## 3. Push main

```bash
git push origin main
```

(`origin/main` is at `f66f96e`, an ancestor; the push is a fast-forward.)

## 4. Wait for CI

The `Release gate` workflow runs on push to `main` (pytest with `.[dev]`, build,
distribution inventory, fresh-install smoke of wheel and sdist).

```bash
gh run list --workflow "Release gate" --branch main --limit 3
gh run watch "$(gh run list --workflow 'Release gate' --branch main --limit 1 --json databaseId --jq '.[0].databaseId')" --exit-status
```

Do not tag until this run is green on the finalize commit
(`gh run list ... --json headSha` must show `git rev-parse main`).

## 5. Tag

```bash
git tag -a v1.6.0 -m "nf-rna v1.6.0" main
git show v1.6.0 --no-patch
git push origin v1.6.0
```

Do not create a `v1.5.1` tag.

## 6. GitHub Release

```bash
gh release create v1.6.0 --title "nf-rna v1.6.0" \
  --notes-file docs/dev/release_notes_v1.6.0.md --verify-tag --latest
gh release view v1.6.0 --web
```

The notes open with the v1.5.0 enrichment-semantics warning for v1.4.x users
and link to the tag's CHANGELOG, `docs/runtime.md` and `docs/gui.md`.

## 7. Tagged-install smoke test

CLI only:

```bash
conda create --yes --name nf-rna-160 python=3.11
conda run -n nf-rna-160 python -m pip install "git+https://github.com/2002brian/nf-rna.git@v1.6.0"
conda run -n nf-rna-160 rnaseq --version      # 1.6.0
conda run -n nf-rna-160 python -c "import importlib.metadata as m; print(m.distribution('nf-rna').read_text('direct_url.json'))"   # vcs_info.commit_id = git rev-parse v1.6.0^{commit}
conda run -n nf-rna-160 rnaseq doctor         # Overall: READY
conda run -n nf-rna-160 rnaseq gui            # exit 1, prints the gui-extra install command
```

With the GUI extra:

```bash
conda create --yes --name nf-rna-160-gui python=3.11
conda run -n nf-rna-160-gui python -m pip install "nf-rna[gui] @ git+https://github.com/2002brian/nf-rna.git@v1.6.0"
conda run -n nf-rna-160-gui python -c "import streamlit, rnaseq; print(streamlit.__version__, rnaseq.__version__)"   # >= 1.66, 1.6.0
tmux new -s nf-rna-160-gui 'conda run -n nf-rna-160-gui rnaseq gui --port 8610'
curl -sf http://127.0.0.1:8610/_stcore/health      # ok
```

In a browser (`ssh -N -L 8610:127.0.0.1:8610 brian@<server>`, then
http://127.0.0.1:8610): Dashboard lists the projects; Results of `airway-e2e`
shows the DEG table and figures; Reports & Downloads downloads `report.html`.
Then stop the tmux session and remove the environments:

```bash
tmux kill-session -t nf-rna-160-gui
conda env remove --yes --name nf-rna-160
conda env remove --yes --name nf-rna-160-gui
```

## Not to be done

- Do not move, recreate or delete existing tags (`v1.3.0`, `v1.4.0`, `v1.4.1`, `v1.5.0`).
- Do not rewrite history to change authorship; commits are authored by Brian Wu with no AI attribution (CLAUDE.md).
- Do not enable the example `systemd --user` unit or lingering as part of the release.
- MRPan is not part of this release.
