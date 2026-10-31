# nf-rna v1.6.1 release checklist

State: v1.6.1 is finalized on `fix/v1.6.1-reliability` with release date
**2026-10-31** (CHANGELOG `## 1.6.1 — 2026-10-31`, CITATION.cff
`date-released: 2026-10-31`). The finalize commit is the head of
`fix/v1.6.1-reliability`; its parent is `89a4510`. `main` and `origin/main`
are at `v1.6.0` (`ff4fe5c`), an ancestor, so `main` fast-forwards. Nothing for
v1.6.1 has been pushed, tagged or released. The v1.6.0 checklist is in the
`v1.6.0` tag.

Readiness, evidence, known findings and limitations:
`docs/dev/release_readiness_v161.md`. The remaining steps are run by the
maintainer, in order; stop at the first failure. If the tag is made on a later
day than 2026-10-31, first change the date in `CHANGELOG.md`, `CITATION.cff`,
the `#161--2026-10-31` anchors (README_latest.md, README_zh-TW.md,
docs/quickstart.md, release notes) and this file, in one commit.

## Scope

- Since `v1.6.0`: D-1, D-2, M-1 (result-changing; `3cf0386`, `7d85046`,
  evidence `1cbe8fb`, `033ce62`); F-1, F-2, F-3 (QC/provenance, no value
  change; `9c37d84`, `89a4510`); release finalization.
- Changed analysis code: `src/rnaseq/r/{gsea_analysis,kegg_analysis,enrichment_common,l2_analysis}.R`,
  new `gsea_leading_edge.R`, `workflow/hisat2_featurecounts.nf` (`STRAND_EVIDENCE`).
- Unchanged: both Conda locks (`git diff --stat v1.6.0 HEAD -- workflow/envs`
  is empty), nf-core/rnaseq 3.26.0, Salmon/HISAT2/featureCounts arguments, CPU
  allocation, scientific thresholds.
- Not in `main`: `validation/v1.6.0` (`b886f11`) and
  `validation/v1.6.1-phase3` (`64c6977`), cited as evidence.

## Release evidence

| Check | Result | Where |
|---|---|---|
| D-1/D-2/M-1 against the independent reference | 33,544/33,544 leading edges; significant GSEA sets 20/20 collections; DEG sets 6/6, ORA 48/48 | `docs/dev/scientific_corrections_v161.md` |
| F-1/F-2/F-3 acceptance (nf-rna `9c37d84`) | 61 checks: 60 PASS, 1 FAIL kept and explained (Salmon `EffectiveLength`); D1 132 tables byte-identical to FIX-D1 | `docs/dev/scientific_reliability_v161.md`, `validation/reliability_v161/results/acceptance.tsv` |
| Phase 2 / Phase 3 frozen verdicts | FAIL under frozen criteria, kept; no open nf-rna computation error | `docs/dev/release_readiness_v161.md` §6 |
| Full suite, locked runtime `0e721ff3…`, `NF_RNA_REQUIRE_R=1` | 1,114 passed, 2 skipped (environment), 0 failed, on code `9c37d84`; later commits change docs, evidence and metadata only | `docs/dev/scientific_reliability_v161.md` |
| `tests/test_release_metadata.py` (finalize tree) | 2 passed | local |
| Build (finalize tree) | `nf_rna-1.6.1-py3-none-any.whl`, `nf_rna-1.6.1.tar.gz`; METADATA `Version: 1.6.1`, extra `gui` → `streamlit<2,>=1.66`; wheel contains the workflows and R scripts | local `python -m build` |
| `tests/release_smoke.py` on the wheel and the sdist (fresh venvs) | pass; `rnaseq --version` 1.6.1 | local |

## Release blockers

None known. Before tagging, CI on the pushed `main` must pass (step 4). No
v1.6.1 commit has run in CI.

## Known limitations (CHANGELOG, release notes, readiness review)

- Upstream results depend on HISAT2/Salmon thread counts (recorded, not fixed);
  Salmon is not bitwise reproducible at the same thread count
  (`EffectiveLength`, `TPM`; `NumReads` identical).
- Strandedness check: HISAT2 route only, warning only, calibrated on one real
  stranded library.
- Not validated with natural data: forward-stranded and single-end libraries,
  mouse/other FASTQ, multi-lane merging, preprocessed input.
- `--reuse-upstream` of a v1.6.0 run is refused (frozen `nfcore.local.config`
  differs); the error hint names only CPU/memory limits.
- GUI limitations of v1.6.0 unchanged.

## 1. Review

```bash
cd /home/brian/projects/nf-rna
git status                                         # clean
git log --oneline v1.6.0..fix/v1.6.1-reliability   # 7 commits, linear
git merge-base --is-ancestor v1.6.0 fix/v1.6.1-reliability && echo ancestor-ok
git diff --stat v1.6.0 fix/v1.6.1-reliability -- workflow/envs   # expect no output
git log --format='%an|%cn|%(trailers:only)' v1.6.0..fix/v1.6.1-reliability | sort -u   # Brian Wu|Brian Wu|
```

## 2. Fast-forward main

```bash
git switch main
git merge --ff-only fix/v1.6.1-reliability
git log --oneline -1                               # the finalize commit
python -m pytest -q tests/test_release_metadata.py
```

## 3. Push main (and, optionally, the evidence branches)

```bash
git push origin main
# Optional, so the evidence cited by the v1.6.1 documents is reachable:
git push origin validation/v1.6.0 validation/v1.6.1-phase3
```

`origin/main` is at `ff4fe5c`, an ancestor; the push is a fast-forward. Do not
merge the validation branches into `main`.

## 4. Wait for CI

```bash
gh run list --workflow "Release gate" --branch main --limit 3
gh run watch "$(gh run list --workflow 'Release gate' --branch main --limit 1 --json databaseId --jq '.[0].databaseId')" --exit-status
gh run list --workflow "Release gate" --branch main --limit 1 --json headSha --jq '.[0].headSha'   # must equal: git rev-parse main
```

Do not tag until this run is green on the finalize commit.

## 5. Tag

```bash
git tag -a v1.6.1 -m "nf-rna v1.6.1" main
git show v1.6.1 --no-patch
git push origin v1.6.1
```

## 6. GitHub Release

```bash
gh release create v1.6.1 --title "nf-rna v1.6.1" \
  --notes-file docs/dev/release_notes_v1.6.1.md --verify-tag --latest
gh release view v1.6.1 --web
```

The notes open with the result-changing warning and link to the tag's
CHANGELOG, `docs/runtime.md` and the readiness review.

## 7. Tagged-install smoke test

```bash
conda create --yes --name nf-rna-161 python=3.11
conda run -n nf-rna-161 python -m pip install "git+https://github.com/2002brian/nf-rna.git@v1.6.1"
conda run -n nf-rna-161 rnaseq --version      # 1.6.1
conda run -n nf-rna-161 python -c "import importlib.metadata as m; print(m.distribution('nf-rna').read_text('direct_url.json'))"   # vcs_info.commit_id = git rev-parse v1.6.1^{commit}
conda run -n nf-rna-161 rnaseq doctor         # Overall: READY
conda env remove --yes --name nf-rna-161
```

## Not to be done

- Do not move, recreate or delete existing tags (`v1.3.0` … `v1.6.0`).
- Do not rewrite history to change authorship; commits are authored by Brian Wu with no AI attribution (CLAUDE.md).
- Do not relabel the frozen Phase 2 or Phase 3 FAIL verdicts as PASS.
- Do not merge or tag the validation branches.
