# nf-rna v1.6.0 release checklist

State: v1.6.0 is prepared on `feature/v1.6.0` (Phase 6a GUI core, Phase 6b
results, reports and release preparation). The version is `1.6.0`, CHANGELOG
says `1.6.0 — Unreleased` and CITATION.cff has no `date-released`; the release
date is set by the finalize commit. Local `main` is at `9a7352b` (the v1.5.1
maintenance commit, never pushed or tagged); `origin/main` is at `f66f96e`.
v1.5.1 is **not** tagged separately: its changes ship in v1.6.0. Nothing for
v1.6.0 has been merged, pushed, tagged or released.

These are the steps for the maintainer to run by hand. Run them in order; stop
at the first failure.

## Scope

- Since `v1.5.0`: the v1.5.1 test-gating fix (`f66f96e`), the v1.5.1 release
  metadata (`9a7352b`), the browser GUI (`rnaseq.gui`, `rnaseq gui`, the `gui`
  and `gui-dev` extras), shared wizard/design/import functions used by CLI and
  GUI, and new `report.html` sections.
- Unchanged: R scripts, Nextflow workflows, the downstream locks, DESeq2/
  enrichment semantics ("v1.5"), every result table and figure (airway rerun:
  135/135 tables byte-identical, `docs/dev/phase6b_report.md`).

## 1. Review

```bash
cd /home/brian/projects/nf-rna
git status                                   # clean
git log --oneline v1.5.0..feature/v1.6.0     # f66f96e, 9a7352b, the 6a and 6b commits
git diff --stat v1.5.0 feature/v1.6.0
git diff v1.5.0 feature/v1.6.0 -- src/rnaseq/r workflow   # expect no output
```

## 2. Finalize (on main)

```bash
git switch main
git merge --ff-only feature/v1.6.0
```

Edit, in one commit:

- `CHANGELOG.md`: `## 1.6.0 — Unreleased` → `## 1.6.0 — YYYY-MM-DD`.
- `CITATION.cff`: add `date-released: YYYY-MM-DD` after `version: "1.6.0"`.
- `docs/dev/release_notes_v1.6.0.md`: no change needed unless the date is mentioned.

```bash
python -m pytest -q tests/test_release_metadata.py
git add CHANGELOG.md CITATION.cff
git commit -m "chore(release): finalize nf-rna v1.6.0"
```

The commit is authored by Brian Wu with no AI attribution trailer (CLAUDE.md).

## 3. Full tests before pushing

With the locked R runtime on PATH (`NF_RNA_REQUIRE_R=1` makes missing R packages fail instead of skip):

```bash
PATH=$HOME/.cache/nf-rna/runtime/prefixes/<current prefix>/bin:$PATH python -m pytest -q
```

And in an environment with only `.[dev]` (no GUI extra, as CI):

```bash
conda create --yes -p /tmp/nf-rna-ci python=3.11
/tmp/nf-rna-ci/bin/python -m pip install -e '.[dev]'
/tmp/nf-rna-ci/bin/python -m pytest -q          # GUI page tests skip; service tests run
```

## 4. Push main and wait for CI

```bash
git push origin main
gh run list --branch main --limit 3
gh run watch <run-id> --exit-status
```

Do not tag until GitHub Actions on the finalize commit is green.

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
gh release view v1.6.0
```

The notes start with the v1.5.0 enrichment-semantics warning for v1.4.x users.

## 7. Clean-environment install test

CLI only:

```bash
conda create --yes --name nf-rna-160 python=3.11
conda run -n nf-rna-160 python -m pip install "git+https://github.com/2002brian/nf-rna.git@v1.6.0"
conda run -n nf-rna-160 rnaseq --version     # 1.6.0
conda run -n nf-rna-160 rnaseq doctor        # Overall: READY
conda run -n nf-rna-160 rnaseq gui           # prints the gui-extra install command, exit 1
```

With the GUI extra:

```bash
conda create --yes --name nf-rna-160-gui python=3.11
conda run -n nf-rna-160-gui python -m pip install "nf-rna[gui] @ git+https://github.com/2002brian/nf-rna.git@v1.6.0"
conda run -n nf-rna-160-gui python -c "import streamlit; print(streamlit.__version__)"   # >= 1.66
conda run -n nf-rna-160-gui rnaseq gui --port 8610   # in tmux; open via ssh -L 8610:127.0.0.1:8610
```

In the browser: Dashboard lists the projects; Results of `airway-e2e` shows the
DEG table and figures; Reports & Downloads downloads `report.html`. Stop the
server (Ctrl-C) and remove the test environments:

```bash
conda env remove --yes --name nf-rna-160
conda env remove --yes --name nf-rna-160-gui
```

## Not to be done

- Do not move, recreate or delete existing tags (`v1.3.0`, `v1.4.0`, `v1.4.1`, `v1.5.0`).
- Do not rewrite history to change authorship.
- Do not enable the example `systemd --user` unit or lingering as part of the release.
- MRPan is not part of this release.
