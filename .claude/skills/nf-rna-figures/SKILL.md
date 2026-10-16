---
name: nf-rna-figures
description: Use when adding, changing, migrating or reviewing a figure in nf-rna (volcano, PCA, heatmap, boxplot, dotplot, figure presets/profiles, fonts, export). Routes to nf-rna's own figure standard, R implementation and review tools.
---

# nf-rna figures

nf-rna owns its whole figure system. Everything you need is inside this repository; do
not consult, import, copy from or depend on any file outside it.

## Read first

1. `docs/figure-standard.md` — the single authoritative standard (status, presets,
   profiles, typography, output contract, verification, contributor rules). Read it
   before changing anything; update it in the same change as the code.
2. `src/rnaseq/r/figures.R` — the implementation and the only place style values live:
   `NF_RNA_STYLE_PRESETS`, `NF_RNA_FIGURE_PROFILES`, `NF_RNA_PRESET_ALIASES`,
   `NF_RNA_STANDARD_FIGURES`, the accessors (`nf_rna_*`), `nf_rna_figure_scope`,
   `nf_rna_save_figure`, the font verification functions.
3. `src/rnaseq/figure_settings.py` and `FiguresConfig` in `src/rnaseq/models.py` — how
   `figures: {preset, profile}` is resolved, aliased, frozen and passed to R.
4. `src/rnaseq/fonts/arimo/FONT_MANIFEST.json` — the bundled font, its license and pinned hashes.

Do not copy palette values, sizes or layout constants into this skill, into scripts or
into other documents; read them from the files above.

## Rules

- Values come from the preset/profile tables through accessors. A script never
  hard-codes a colour, shape, size, alpha or font family.
- A figure type uses the standard preset only if it is listed in `NF_RNA_STANDARD_FIGURES`
  and is drawn through `nf_rna_figure_scope(context, type)`. Other figure types keep the
  legacy appearance until they are migrated; do not restyle unrelated figures.
- Never change scientific calculations, thresholds, normalization, clustering, gene
  selection, source-data table contents, output file names/paths/formats or the
  manifest schema (`nf-rna.figure-manifest.v1`).
- Legacy output must stay unchanged. The standard preset's font is the bundled Arimo with
  no substitution; never fall back to Arial, Helvetica, Liberation Sans or DejaVu Sans
  for it, and never install or download fonts at run time.
- Colours that the standard marks provisional stay provisional until the owner approves
  them; do not finalize them.
- Presets and profiles are independent; a profile's size limits and raster dpi always apply.

## Workflow for a figure change

1. Make the change in `figures.R` / the analysis script (`l1_analysis.R`, `l2_analysis.R`,
   `go_analysis.R`, `gsea_analysis.R`, `kegg_analysis.R`).
2. Run the plotting tests: `tests/test_figures.py`, `tests/test_figure_presets.py`,
   `tests/test_figure_settings.py`, `tests/test_figure_review_tools.py`.
3. Compare actual output with the review tools in `tools/figure_review/`
   (usage in each script's `--help`):
   - `make_demo_project.py` for a realistic project (the shipped example is degenerate),
   - `generate_runs.py` once with the legacy selection and once with the candidate one,
   - `compare_runs.py --require-identical-except <figure kinds>` to prove that only the
     intended figure kinds changed and to inspect the legacy / candidate / side-by-side
     PNGs, physical sizes, declared vs embedded PDF fonts and source-data identity.
   Look at the images; automated checks are not a substitute.
4. Update `docs/figure-standard.md` and the tests with the code.
5. For packaging-affecting changes (new assets, scripts, fonts), build the wheel and run
   a clean installation outside the repository to confirm every asset is shipped.
