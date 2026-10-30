# nf-rna v1.5.1 release checklist

State: v1.5.1 is a maintenance release prepared on `main`, on top of `f66f96e` (the R-package test-gating fix). The published `v1.5.0` annotated tag stays on `c283d31` and is not moved. No GitHub Release exists for v1.5.0. The release metadata is dated 2026-10-30. Nothing for v1.5.1 has been tagged, pushed or released yet.

## Scope

- Changed since `v1.5.0`: three test files (`f66f96e`) and the v1.5.1 version and release metadata. No source, workflow, R script, lock, reference or example change.
- Scientific analysis semantics are identical to v1.5.0 (DESeq2, GO/KEGG ORA, GSEA, significance semantics "v1.5"). v1.5.0 results do not need to be recomputed.

## Platform policy

The supported runtime is **Linux x86-64, including Windows WSL2** (`docs/runtime.md`; `execution_backend()` refuses every other platform). The repository keeps `workflow/envs/locks/nf-rna-downstream-osx-arm64.lock.yml` as an **unqualified candidate**: resolved on linux, never installed or smoke-tested on a Mac. This is not a release blocker.

## Release evidence

| Check | Result | Where |
|---|---|---|
| GitHub Actions on `f66f96e` | SUCCESS | GitHub Actions, `main` |
| v1.5.0 release evidence (full suite, airway end-to-end, Phase 1/2 acceptance, figure gate) | still applies: no analysis code changed | `v1.5.0` checklist at tag `v1.5.0` |
| Release metadata, package build, KEGG-gating tests | maintenance release commit | maintenance release commit |

## Remaining manual steps

1. Review: `git log --oneline v1.5.0..main` (expect `f66f96e` and the v1.5.1 release commit) and `git diff --stat v1.5.0 main`.
2. Push `main`:

   ```
   cd /home/brian/projects/nf-rna
   git push origin main
   ```

   Wait for GitHub Actions on the release commit to pass.
3. Annotated tag on the release commit, in the style of earlier tags:

   ```
   git tag -a v1.5.1 -m "nf-rna v1.5.1" main
   git show v1.5.1 --no-patch
   git push origin v1.5.1
   ```

4. GitHub release from the tag, with the CHANGELOG 1.5.1 section as notes:

   ```
   awk '/^## 1.5.1 /{flag=1; next} /^## 1.5.0 /{flag=0} flag' CHANGELOG.md > /tmp/nf-rna-v1.5.1-notes.md
   gh release create v1.5.1 --title "nf-rna v1.5.1" --notes-file /tmp/nf-rna-v1.5.1-notes.md --verify-tag --latest
   ```

   The 1.5.1 notes describe only the maintenance fix; the scientific changes are in the 1.5.0 section of CHANGELOG.md.
5. Check the release:
   - Install from the tag in a fresh env: `python -m pip install "git+https://github.com/2002brian/nf-rna.git@v1.5.1"`.
   - `rnaseq --version` must print 1.5.1, and `rnaseq doctor` must report READY.

## Not to be done

- Do not move, recreate or delete existing tags (`v1.3.0`, `v1.4.0`, `v1.4.1`, `v1.5.0`).
- Do not rewrite history to change authorship. Commits are authored by Brian Wu without AI attribution trailers (CLAUDE.md).
- MRPan is not part of this release.
