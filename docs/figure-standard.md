# nf-rna figure standard

This is the authoritative description of how nf-rna draws, sizes, names and records
its figures. nf-rna owns its whole figure system: every plotting function, colour,
layout rule, font policy and export step lives in this repository
(`src/rnaseq/r/figures.R` and the analysis scripts that call it), and the font it
draws with ships inside the package (`src/rnaseq/fonts/arimo`). No figure code,
configuration or font is read from outside nf-rna at install, test or run time, and
nothing outside nf-rna is a dependency. Visual conventions adopted from other
published guidance are implemented here natively and specified in this document, not
by reference to another file.

## Status

| Item | State |
|---|---|
| `legacy` preset | Active. Reproduces the appearance nf-rna has always produced. |
| `standard` preset | **Active for volcano, PCA, sample boxplots, the sample correlation heatmap, DEG heatmaps and the enrichment dotplots** (`volcano`, `pca`, `pca_by_batch`, `boxplot_*`, `sample_correlation`, `heatmap`, and the GO, GSEA and KEGG ORA/GSEA dotplots `dotplot_*`, `dotplot`). Every figure type is now migrated; the legacy preset remains for projects that select it. |
| `bioinfo` | Backward-compatible alias of `standard`; resolved to `standard` everywhere. |
| Profiles (`general`, `nature`, `plos`, `elsevier`) | Active with either preset. |
| Arimo font | Bundled (SIL OFL 1.1), verified before the first standard figure is drawn. |
| Enrichment colour scale | **Approved** as `nf-rna Standard Enrichment Scale v1` (see Continuous scales). |
| Group colours 1–8 | **Approved** (1–3 earlier, 4–8 in Phase 4). Beyond eight levels slots 1–8 are kept and a deterministic maximin extension supplies the rest. |
| Batch palette | **Approved** as `nf-rna Standard Batch Palette v1`. |

Migrated figure types are listed once, in `NF_RNA_STANDARD_FIGURES` (`figures.R`): `volcano`, `pca`, `boxplot`, `correlation`, `heatmap`, `dotplot`.

## Selecting a style

```yaml
figures:
  preset: standard     # standard | legacy   (bioinfo is accepted as an alias of standard)
  profile: general     # general | nature | plos | elsevier   (default: general)
```

* The **preset** controls the visual style. The **profile** controls the publication
  specification (column widths, size limits, text sizes, raster dpi). They are
  independent: every preset works with every profile, and a selected profile is
  never ignored. A profile's raster dpi always applies (300; Elsevier 500).
* `rnaseq new` writes `preset: standard, profile: general`.
* A `project.yaml` without `figures` resolves to `legacy` / `general`, and so does a
  frozen downstream contract without `figures`.
* `bioinfo` is normalised to `standard` when a project or a frozen contract is read;
  frozen files themselves are never rewritten. A contract and a project that spell the
  same preset differently are not a "mixed" run.
* The resolved setting is frozen with the run (`downstream_contract.json`, key
  `figures`). Retries clone it unchanged. The contract and the frozen project must
  agree; a mismatch stops the run ("mixed figure settings").
* Changing the figure selection does not affect upstream reuse
  (`--reuse-upstream`): reuse compares the nf-core inputs and the input, upstream and
  reference identity only.
* The R backends receive the resolved setting as `figure_style`
  (`{preset, profile, source}`) in their JSON config. The older `figures.profile` key
  is still read if `figure_style` is absent; if both are present they must agree.
  An unknown preset or profile stops the run.

Single sources of truth: Python resolves defaults and aliases only in
`rnaseq/figure_settings.py` (and `PRESET_ALIASES` in `rnaseq/models.py`); R holds all
style values only in `NF_RNA_STYLE_PRESETS`, `NF_RNA_PRESET_ALIASES` and
`NF_RNA_FIGURE_PROFILES` in `figures.R`. Scripts read them through accessors
(`nf_rna_categorical_colors`, `nf_rna_direction_colours`, `nf_rna_point_alpha`,
`nf_rna_pca_label_pt`, `nf_rna_layout_width_mm`, …) and never hard-code a colour, shape,
size or font. A figure type is drawn with the standard preset by asking for its scoped
context, `nf_rna_figure_scope(context, "volcano")`; any other type receives the legacy
context unchanged.

## Group colours

The colour of a group is decided by its position in the **universe** of the design variable's levels
(`levels(metadata[[variable]])`, the complete ordered set), never by the levels a particular figure
contains. `nf_rna_group_colours(universe, present, context)` returns the colours of `present`, so
a subset, a reordering or a missing group never moves a colour, and PCA, boxplots, the correlation
heatmap annotation and the DEG heatmap annotation always agree. Batch uses its own palette
(`nf_rna_batch_colours`), so group and batch colours never coincide. Beyond eight levels the
standard preset keeps slots 1–8 and appends colours chosen by `nf_rna_extend_palette`: one at a
time from a fixed HCL pool, each maximising the smallest CIE Lab distance to all earlier colours
under normal vision and simulated deuteranopia, protanopia and tritanopia (no randomness; a longer
palette always starts with the shorter one). The legacy preset keeps its historic HCL “Dark 3” set.
The manifest names the strategy per figure (`encoding.palette`). Groups 1–8 are approved; the batch
palette is approved as `nf-rna Standard Batch Palette v1` (its manifest name; beyond eight levels the name carries the suffix “+ deterministic maximin extension from slot 9”). Distinguishability of the implemented palettes is reported
by `tools/figure_review/palette_extension_review.R`: slots 1–8 are at least 18.6 CIEDE2000 apart
under normal vision (6.5 under deuteranopia, limited by the approved groups 2 and 3); the 9th colour
is at least 12.8 from every earlier one under every simulation, and beyond ten levels colour
vision deficiencies make some pairs weak, so shape encodes identity as well.

## Continuous scales

Three kinds of colour are kept separate: the **discrete group/batch palettes** (categories), the **diverging
expression scale** of the DEG heatmap and the correlation heatmap (a value with a meaningful zero), and the
**sequential enrichment scale** (adjusted P of the GSEA dotplots; its first, most-significant stop is also the
single mark colour of the ORA dotplots). Low adjusted P is dark in both presets, and the scale limits (0 to the
GSEA adjusted-P cutoff) are unchanged. The standard enrichment scale is a Viridis-hue ramp of its own and shares
no colour with the group, batch or diverging colours. It is approved as `nf-rna Standard Enrichment Scale v1`
(stops `#4B0055 #353E7C #007094 #009B95 #00BE7D`, recorded with that name in the manifest `colour_policy.enrichment`);
its colour-vision limitations are those documented for the palettes above.

Review finding for the expression scale: the standard diverging ends (`standard$diverging`) are the same
colours as group colours 1 and 2, so an annotation bar and the expression extremes share colours. This follows the
approved values and is unchanged; `tools/figure_review/scale_review.R` reports the contrast and colour-vision
metrics of the current scale and alternatives without changing anything.

## Profiles (publication specification)

| Profile | Single / double column (mm) | Max (mm) | Base / small / min text (pt) | Raster dpi |
|---|---|---|---|---|
| general | 85 / 180 | 180 × 220 | 8 / 7 / 6 | 300 |
| nature | 89 / 183 | 183 × 170 | 7 / 6 / 5 | 300 |
| plos | – / 190.5 | 190.5 × 222.3 | 9 / 8 / 8 | 300 |
| elsevier | 90 / 190 | 190 × 240 | 8 / 7 / 7 | 500 |

These are layout targets taken from publisher artwork guidance, not a guarantee that a
journal will accept a figure. PNG previews are drawn at 150 dpi.

## Presets (visual style)

The values below are defined in `NF_RNA_STYLE_PRESETS` (`figures.R`); this table
describes them and must be updated together with the code.

| Property | `legacy` | `standard` (migrated figure types) |
|---|---|---|
| Group colours 1–3 | `#0072B2` `#D55E00` `#009E73` | `#4E78A7` `#CC5A6A` `#59A14F` (approved) |
| Group colours 4–8 | Okabe–Ito order | `#843C39 #56B4E9 #BCBD22 #6A3D9A #8C6D31` (approved) |
| Batch colours | same as groups | independent palette `#393B79 #585858 #17BECF #E6AB02 #CC79A7 #A6761D #7F7F7F #E08214` (`nf-rna Standard Batch Palette v1`, approved) |
| Volcano up / down / not significant | `#D55E00` / `#0072B2` / `#BDBDBD` | `#B2182B` / `#2166AC` / `#BDBDBD` |
| Heatmap diverging low / mid / high | `#2166AC` / `#F7F7F7` / `#B2182B` | `#4E78A7` / `#F7F7F7` / `#CC5A6A` |
| Sequential (correlation, GSEA) | ColorBrewer Blues | unchanged |
| Missing values | `#A0A0A0` | `#A0A0A0` |
| PCA label size | profile minimum text size | 6 pt, raised to the profile minimum where that is larger (PLOS 8 pt, Elsevier 7 pt) |
| Volcano point alpha | opaque | 0.85 |
| Volcano legend counts | yes | yes (drawn genes only) |
| Shapes | 12 fixed shapes, then colour only | unchanged (PCA keeps its shape distinctions) |
| Volcano size | profile double column × 120 mm | profile single column, square (aspect 1.0), legend stacked below |
| Sample boxplot size | double column × (65 mm + label height) | content-aware: width from the number of samples (minimum 60 mm, at most the profile double column), height from the measured labels and a 42 mm panel; the legend wraps into more rows until it fits |
| Sample correlation heatmap | double column × 150 mm, no annotation | content-aware square cells (6 mm, shrinking only when the profile maximum would be exceeded; sample labels are hidden only when a cell is smaller than the label height); group (and batch) annotation in the PCA colours |
| Enrichment dotplots (GO, GSEA, KEGG) | double column, height from the number of label lines, legend on the right | content-aware: width = measured label/axis width + a 55 mm panel (minimum 85 mm), height = measured fixed parts + one row per term and label line; legend below, stacked when it is wider than the figure; subtitle and x-axis title wrapped to the figure width with the real font metrics; term text, order, point size and colour mappings, axes and numeric labels unchanged |
| Enrichment colours | ORA marks use the primary colour; GSEA adjusted P uses reversed Blues | ORA marks use the first stop of the enrichment scale; GSEA adjusted P uses the enrichment scale (low P dark), limits 0 to the cutoff as before |
| DEG heatmap | double column, height from rows | content-aware: 5 mm sample columns, row height from the label font (7 pt down to the profile minimum, then labels are hidden and rows are 1.5 mm); width and height are measured from labels, trees, legend and title, never above the profile maximum |
| PCA size | profile double column, height from the data, legend right | profile single column, height from the data (minimum 40 mm), legend below |
| PCA with constant PC2 | y range ±20% of the PC1 span, one y break at the constant (see below) | same |
| Font | Arial, else DejaVu Sans (recorded fallback) | Arimo only, no substitution |

“Single column” is the profile's single column; for a profile without one (PLOS) the standard preset uses a compact 85 mm (`compact_width_mm`, never above the profile maximum). Legacy figures are unchanged.

### PCA with a constant PC2

Both presets draw PCA with equal axis scaling, so the panel's aspect ratio is that of the scores. When PC2 has no variance (two samples always; more samples only if they lie on one line), its scores are constant (0 up to floating-point rounding, as `prcomp` centres them) and an equal-scale panel would be a thin strip whose default y breaks overlap. `nf_rna_pca_y_axis()` (`figures.R`) treats PC2 as constant when its range is no larger than `sqrt(.Machine$double.eps)` times the PC1 range, and then, in both the legacy and the standard preset:

* displays a y range of ±20% of the PC1 span around the constant (0, or the constant itself if it is not 0 within that tolerance); the figure height is sized for this range;
* draws a single y break, at the constant, so no tick label can overlap or repeat;
* records it in the figure manifest: the PCA entry's `axes` ends with `PC2 is constant (no variance): y range fixed at +/-20% of the PC1 span with a single break at 0`.

Only the displayed y range and breaks change. `pca_scores.tsv` and `pca_variance.tsv` (for example PC1 = 100%, PC2 = 0%), the axis titles, equal axis scaling, sample positions, labels, colours and shapes, and the PNG/TIFF/PDF exports are unchanged; no separation is added. A PC2 with real spread, however small, keeps the previous display (at least ±1% of the PC1 span), so multi-sample PCA figures are unaffected.

## Typography and the bundled font

* **legacy**: prefers Arial; if it is not installed it uses DejaVu Sans and records
  the fallback and its reason in the manifest. Unchanged.
* **standard**: Arimo only. Arial, Helvetica, Liberation Sans and DejaVu Sans are never
  substituted.
  * Arimo 1.33 is bundled in `src/rnaseq/fonts/arimo`: four static TrueType files
    (regular, bold, italic, bold italic), the SIL Open Font License 1.1 text, and
    `FONT_MANIFEST.json` (source repository and commit, license, SHA-256 and size of
    every file). Static rather than variable files are bundled because Cairo embeds
    them under their own names (`Arimo-Regular`, …), whereas a variable font's bold
    instance is embedded under a generic `CairoFont-*` name that cannot be verified.
  * Source: https://github.com/googlefonts/arimo (commit pinned in the manifest).
    Arimo is not available as a conda-forge or bioconda package, so it is vendored.
  * The OFL permits redistribution with the license text; the font is not sold by
    itself and carries no Reserved Font Name.
* Before the first standard figure is drawn nf-rna (1) checks the bundled files against
  the pinned SHA-256, (2) points fontconfig at the bundled directory with a generated
  configuration (any other Arimo build is hidden), (3) checks that `systemfonts` lists
  exactly the bundled files, and (4) draws a probe PDF in all four faces and reads the
  embedded font names back, which must be exactly the four Arimo names. fontconfig
  would otherwise map an unknown family to another font without any message. Any
  failure stops the run with an actionable error before a figure is drawn.
* fontconfig reads its configuration once per R process, so activation happens when the
  figure context is created. Nothing is installed into the system or the Conda
  environment.

## Output contract (unchanged by presets)

* Every figure is drawn directly at its final physical size on each device, with no
  resampling: PNG preview, raster TIFF (LZW, RGB, opaque white) and Cairo PDF. Sizes are
  snapped to whole points so the three outputs share one size.
* `nf_rna_save_figure(..., formats = ...)` also accepts `png_raster`: an additional PNG,
  `<name>_<raster dpi>dpi.png`, drawn at the profile's raster dpi. Since v1.5.0 the L2 volcano,
  MA plot (`ma_plot`, drawn with the volcano scope) and DEG heatmap use it; the 150 dpi
  preview keeps its name and remains the report image.
* Since v1.5.0 the volcano and MA plot draw `log2FoldChange_shrunk` (raw LFC only as a labelled
  fallback when shrinkage failed), label the top 10 DEGs by symbol, and DEG heatmap rows are labelled
  by symbol; see `docs/scientific_contract.md`. Colours and significance are unchanged.
* Since v1.5.0 every enrichment dotplot (ORA and GSEA, including the GO `*_simplified` dotplots) also
  writes `png_raster`. Dotplot selection and the GSEA colour limits follow the v1.5 significance rule
  (`p.adjust <= 0.05`; colour limits [0, 0.05]).
* Since v1.5.0, GSEA running-score plots (`running_score/<sign>_<n>_<ID>`, kind
  "GSEA running-score plot") are drawn with `enrichplot::gseaplot2` through
  `nf_rna_gsea_running_score` in PNG preview, `png_raster` and PDF. They are not a standard-preset
  figure type (not in `NF_RNA_STANDARD_FIGURES`): they use the run's context, so they are identical
  under the legacy and standard presets. The colour comes from `nf_rna_enrichment_mark`, the base size
  from the profile, and the size is double column × 0.7.
* File names and paths, the three formats and the source-data tables
  (`*_source_data.tsv`, `pca_scores.tsv`, `all_genes.tsv`, …) do not depend on the
  preset or profile.
* Standard heatmaps record the drawn row/column order (`display_order`) and their cell sizes in
  their figure record; the order is pheatmap's / `nf_rna_cluster`'s, never a new clustering.
* `figure_manifest.json` keeps schema `nf-rna.figure-manifest.v1`. A `figure_style`
  object (requested name, canonical preset, profile, migrated types, types drawn, the
  verified font with its SHA-256, embedded names and bundled-font provenance) is added
  only when a non-legacy preset is requested, and each figure drawn with the standard
  preset carries its own `figure_style` entry. Legacy manifests are unchanged.
  When standard figures were drawn, the manifest-level `font` and `colour_policy` describe the standard
  preset actually used (Arimo with its SHA-256; the standard direction, diverging, sequential and group
  colours), and `colour_policy` gains the additive keys `applied_preset`, `categorical_colours`, `batch`
  and `enrichment`. Batch figures name the batch palette in `encoding.palette`. Legacy manifests keep
  exactly their original keys and values.
* Figures never alter a scientific result: z-scores, saturation, clustering and gene
  selection are display-only and documented in the manifest.

## Verification

* `tests/test_figures.py`, `tests/test_figure_presets.py`,
  `tests/test_figure_settings.py`: policy, presets/aliases/profiles, propagation,
  bundled font, retry and frozen-contract compatibility.
* `tools/figure_review/make_demo_project.py`: a deterministic synthetic project (two
  groups, three batches, optional multi-level grouping variable and long sample labels) that
  gives a realistic volcano, batch PCA, boxplots and heatmaps.
* `tools/figure_review/review_matrix.py`: legacy-versus-standard reports over 6, 12 and 26
  samples, four groups and long labels, gated on "only migrated kinds change" and on the drawn
  clustering order equalling the legacy code path's order (`display_order.R`).
* `tools/figure_review/review_enrichment.py`: legacy-versus-standard reports for GO, GSEA, KEGG ORA and
  KEGG GSEA over many/few terms, long/short names, non-significant and empty results. `compare_runs.py` reports for
  every dotplot the dimensions, embedded PDF fonts, source-data identity, term-order identity (and drawn order =
  source order), legend/label completeness (words found in the PDF) and a clipping check (ink on the outermost
  pixel frame of the PNG); the gate fails on any of them.
* Offline KEGG: `tools/figure_review/kegg_offline_profile.R` (loaded through `R_PROFILE_USER` by
  `generate_runs.py --kegg-scenario ...`) replaces only the KEGG network calls (the resource probe, `enrichKEGG`,
  `gseKEGG`) so the real KEGG backend runs offline. Pathway names are real; **the statistics are synthetic**, every
  output directory is marked `SYNTHETIC_KEGG_FIXTURE.txt`, and the fixture is never used in production.
* `tools/figure_review/scale_review.R`: review of the continuous scales.
* `tools/figure_review/heatmap_sizes.R`: standard heatmaps for 20/50/100 genes, 6–60 samples and
  short/long labels, with a size table.
* `tools/figure_review/generate_runs.py`: reproducible figure sets from a project
  (L1, L2, and synthetic GO/GSEA) using only the installed package.
* `tools/figure_review/compare_runs.py`: for each matched figure, the legacy and
  candidate PNGs, a side-by-side PNG, a pixel-difference image, actual physical and
  pixel dimensions, declared font versus the font names embedded in the PDF, and
  source-data identity. `--require-identical` makes it a parity gate;
  `--require-identical-except volcano PCA ...` (manifest kinds) allows only the listed figure kinds to
  change (they must still carry identical source data and a matching embedded font).
* Release audit tools: `tools/figure_review/build_gallery.py` (one self-contained HTML gallery of
  legacy/standard pairs with dimensions, fonts, palettes, source-data identity and visual warnings),
  `palette_audit.R` (current palettes plus unapproved proposals; groups 1-3 fixed), `profile_audit.py`
  (limits, dpi, exact PDF font sizes per profile), `invariance_audit.py` (byte comparison of every
  non-figure file of two runs), `manifest_review.py` (manifest backward compatibility and provenance
  fields) and `package_audit.py` (wheel/sdist font hashes, licence, no external references). A live KEGG
  read-only check is `generate_runs.py --kegg-scenario live`; it is never part of the test suite.
  `pdf_layout_repeat.py` repeats one offline KEGG render and reports distinct PDF layouts (see the known limitation).
* `tools/figure_review/palette_review.R`: contrast and colour-vision-deficiency review
  of candidate palettes.

## Rules for contributors

1. Add or change a colour, shape, size or font only in the preset or profile table.
2. Migrating a figure type: add it to `NF_RNA_STANDARD_FIGURES`, obtain its scoped
   context with `nf_rna_figure_scope`, read every style value through an accessor, add a
   `layout` entry for its size, and prove with `compare_runs.py --require-identical-except`
   that nothing else changed.
3. Keep legacy output byte-stable (modulo PDF creation dates) when changing shared code;
   prove it with `compare_runs.py --require-identical` against a legacy run.
4. Do not read files, fonts or configuration from outside the repository or the
   installed package, and do not add a dependency for plotting.
5. Do not change statistics, thresholds, normalization, clustering or gene selection in
   figure code.

## Stale text metrics across graphics devices (fixed in Phase 4B)

**Symptom.** Drawn on several Cairo devices in one R process (PNG, TIFF, PDF), a figure's TIFF or PDF was
intermittently laid out differently from its PNG: rows holding text (the x-axis row, the bottom legend)
were off by a few points and the panel correspondingly shorter or taller. Pixels of the correct draws,
source data, manifests and page sizes were not affected. In the png -> tiff -> pdf sequence after the sizing
probe it hit 3.07% of outputs (276 of 9000: TIFF 218, PDF 58).

**Cause (established by intervention, not by correlation).** R remembers the metric it last computed for a
font specification and reuses it when the next text height request names the same font, without regard to
which device it was computed on. When the device changes resolution (PDF 72, PNG 150, TIFF 300/500 dpi), the
first measurement on the new device can therefore come back in the previous device's units: exactly 1/2, 1/4 or
2 times the right height. The row-height differences measured in drawn layouts equal those factors (2.4, 3.6
and 5 pt for a 7 pt label). A plain R + grid script with no ggplot2 and no nf-rna code reproduces it; it appears
only when the first font measured on a device equals the last one measured on the previous device (0 of 3200
runs otherwise), it is more frequent after a `gc()` or a short busy period, and a throw-away measurement in
another font at the start of each device removed it entirely (0 of 759 runs, against 262 of 363 and 323 of 381
without). It is not caused by fontconfig, `systemfonts`, the shared font cache, CPU load itself, garbage
collection timing or the figure's gtable (each ruled out by experiment). Why the hit is intermittent
(roughly 0-50% depending on incidental state) is not explained further.

**Fix.** `nf_rna_prime_device_metrics()` measures one string in a font no figure uses (1.37 pt bold italic)
immediately after a device is opened, before anything is measured or drawn: in `nf_rna_save_figure` for every
output device and in `nf_rna_with_probe_device`. It draws nothing and changes no figure. Legacy figures go
through the same helpers and are protected the same way; their output is unchanged when it was correct.

**Evidence for the chosen fix.** Natural sequence (probe device, then png, tiff, pdf), 3000 runs under CPU
load per arm: control 276 of 9000 outputs wrong; priming 0 of 9000 (95% upper bound 0.043%); isolated process
per output 0 of 9000 (0.128% per format); draw-verify-redraw 0 of 9000 and, with a fault injected after every
draw, 1,106 detections and 0 misses in 4500 outputs. Production code: 0 of 9000 natural and 0 of 4500
fault-injected outputs wrong. Whole pipeline, 120 complete runs under load (6,120 figure files): 0 deviating
with the fix (upper bound 0.063%), 3 deviating without it (all the KEGG GSEA PDF). Cost per figure: priming
+0.03 s, verify-redraw +0.34 s, isolated processes +1.4 s.

Reproduce and measure with `tools/figure_review/pdf_layout_repeat.py` (a correct installation yields a single
layout); `tests/test_figure_device_metrics.py` guards that every device is primed.
