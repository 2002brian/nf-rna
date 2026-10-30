# nf-rna v1.5.0 release checklist

State: the release is finalized on `feature/v1.5.0`. The release metadata is dated 2026-10-30. Nothing has been merged, tagged, pushed or released yet. `feature/v1.5.0` was branched from `main` at `f69058e` (v1.4.1), and `main` has not moved since, so `main` can be fast-forwarded.

## Platform policy

The supported runtime for v1.5.0 is **Linux x86-64, including Windows WSL2** (`docs/runtime.md`; `execution_backend()` refuses every other platform). macOS, and therefore osx-arm64, is not a supported runtime. A future GUI would reach a Linux server from a browser and does not change this.

- The repository keeps `workflow/envs/locks/nf-rna-downstream-osx-arm64.lock.yml` as an **unqualified candidate**. It was resolved on linux, has never been installed or smoke-tested on a Mac, and the docs say so.
- This is not a release blocker. A Mac validation belongs to a future release that decides to support macOS (see `.agents/skills/nf-rna-runtime-migration/SKILL.md`: a platform becomes supported only after a real native smoke execution).

## Release evidence (no rerun needed)

| Check | Result | Where |
|---|---|---|
| Full test suite (on `aabea4a`) | 1056 passed, 2 skipped (real SAMtools/featureCounts; absent local SJ project) | `/home/brian/nf-rna-validation/v150-release-prep-full-suite.log` |
| End-to-end `rnaseq run` on the airway counts project (L2 + GO/KEGG/GSEA, Nextflow) | SUCCESS. All 125 tables byte-identical to the Phase 2 acceptance results; summaries, contract metadata and snapshot manifest identical apart from paths and timestamps; `report.html` inspected | `airway-e2e/runs/AIRWAY-E2E/20261030-153417+0800` |
| Phase 1 DEG acceptance | DESeq2 statistics and significant sets identical to v1.4.1 | `docs/dev/phase1_report.md` |
| Phase 2 enrichment acceptance | significance ≡ p.adjust ≤ 0.05; universe = definition; snapshot SHA-256 verified; two runs byte-identical; DEG outputs unchanged | `docs/dev/phase2_report.md` |
| Figures, v1.4.1 vs v1.5.0 (`compare_runs.py --expected-changes tools/figure_review/expectations/v1.4.1-to-v1.5.0.json`) | gate PASS: 59 declared changes, 0 undeclared differences, 0 removed figures; L1 figures identical | figure-review runs from the release-prep session |
| Finalization checks | release metadata, figure-review tools, example validation | finalization commit |

## Remaining manual steps

1. Review: `git log --oneline v1.4.1..feature/v1.5.0`.
2. Merge into `main` (fast-forward only):

   ```
   cd /home/brian/projects/nf-rna
   git checkout main
   git merge --ff-only feature/v1.5.0
   git log --oneline -1   # must be the finalization commit
   ```

3. Annotated tag, in the style of earlier tags:

   ```
   git tag -a v1.5.0 -m "nf-rna v1.5.0"
   git show v1.5.0 --no-patch
   ```

4. Push the branch and the tag:

   ```
   git push origin main
   git push origin feature/v1.5.0   # optional: keep the feature branch on the remote
   git push origin v1.5.0
   ```

5. GitHub release from the tag, with the CHANGELOG 1.5.0 section as notes:

   ```
   awk '/^## 1.5.0 /{flag=1; next} /^## 1.4.1 /{flag=0} flag' CHANGELOG.md > /tmp/nf-rna-v1.5.0-notes.md
   gh release create v1.5.0 --title "nf-rna v1.5.0" --notes-file /tmp/nf-rna-v1.5.0-notes.md --verify-tag
   ```

   The notes open with the enrichment semantics change (BH-adjusted p ≤ 0.05) and the advice to pin the KEGG snapshot.
6. Check the release:
   - Install from the tag in a fresh env: `python -m pip install "git+https://github.com/2002brian/nf-rna.git@v1.5.0"`.
   - `rnaseq --version` must print 1.5.0, and `rnaseq doctor` must report READY.
   - The first run provisions the new downstream runtime (lock `0e721ff3…`).

## Not to be done

- Do not move or recreate existing tags (`v1.3.0`, `v1.4.0`, `v1.4.1`).
- Do not rewrite history to change authorship. Commits are authored by Brian Wu without AI attribution trailers (CLAUDE.md).
- MRPan is not part of this release.
