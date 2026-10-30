# nf-rna v1.5.0 release checklist

State: release prep on `feature/v1.5.0`. These steps are not automated and remain to be done by hand, in order. Nothing here has been tagged, pushed or released.

`feature/v1.5.0` was branched from `main` at `f69058e` (v1.4.1) and `main` has not moved since, so the branch can be fast-forwarded.

## Remaining manual steps

1. **Review** the branch: `git log --oneline v1.4.1..feature/v1.5.0` and `docs/dev/phase1_report.md`, `docs/dev/phase2_report.md`, `docs/dev/DECISIONS.md`.
2. **osx-arm64 on a real Mac.** The osx-arm64 downstream lock was resolved on linux
   (`CONDA_SUBDIR=osx-arm64`, `CONDA_OVERRIDE_OSX=13.0`) and has never been installed natively.
   - On an Apple-silicon Mac, create the env:
     `conda env create --prefix <dir> --file workflow/envs/locks/nf-rna-downstream-osx-arm64.lock.yml`
   - Probe the packages with `Rscript -e` and `requireNamespace(...)`: DESeq2, apeglm, ashr, clusterProfiler, enrichplot, KEGGREST, GOSemSim, GO.db, openssl, org.Hs.eg.db, org.Mm.eg.db.
   - Run `pytest tests/test_l2.py tests/test_go.py tests/test_gsea.py tests/test_kegg.py` with that R on PATH.
   - Run the airway L2 + enrichment once and compare its tables with the linux run.
   - If osx-arm64 is not qualified, say so in the release notes; it is not a supported runtime platform in `docs/runtime.md` anyway.
3. **Pin the KEGG snapshot for any result that will be published.** This is a recommendation for users; nothing needs to be done for the release itself.
4. **Finalize the release metadata** in one commit, as in `f69058e`:
   - `CHANGELOG.md`: `## 1.5.0 — Unreleased` becomes `## 1.5.0 — <YYYY-MM-DD>`, and the "Development version … Not released." line is removed.
   - `CITATION.cff`: add `date-released: <YYYY-MM-DD>` (version is already `1.5.0`).
   - README (en, zh-TW) and `docs/quickstart.md`: latest stable release becomes `v1.5.0`, with the dated CHANGELOG anchor.
   - `docs/runtime.md`: "Supported runtime (v1.4.1)" becomes v1.5.0.
   - Run `pytest tests/test_release_metadata.py` (the dates must agree) and the full suite.
5. **Merge into `main`** (fast-forward):

   ```
   git checkout main
   git merge --ff-only feature/v1.5.0
   ```

6. **Tag**, as an annotated tag in the style of the earlier tags:

   ```
   git tag -a v1.5.0 -m "nf-rna v1.5.0"
   git show v1.5.0 --stat | head
   ```

7. **Push** the branch and the tag:

   ```
   git push origin main
   git push origin v1.5.0
   ```

8. **GitHub release.** Create it from tag `v1.5.0` with the CHANGELOG 1.5.0 section as the notes. State prominently that **enrichment significance semantics changed** (v1.5: p.adjust (BH) ≤ 0.05; ORA universe = padj-non-NA genes; KEGG from a local snapshot), and say that v1.4 results are not recomputed. For example:
   `gh release create v1.5.0 --title "nf-rna v1.5.0" --notes-file <CHANGELOG 1.5.0 section>`
9. **After release:**
   - The first `rnaseq run` on each machine provisions the new downstream runtime (lock `0e721ff3…`, which adds apeglm and ashr). Older prefixes from the `18197d49…` lock stay valid for v1.4.1 only.
   - MRPan is not part of this release and was not changed.

## Not to be done on this branch

- Do not move or recreate existing tags (`v1.3.0`, `v1.4.0`, `v1.4.1`).
- Do not rewrite history to change authorship. Commits are authored by Brian Wu without AI attribution trailers (CLAUDE.md).
