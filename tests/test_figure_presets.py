"""R-side figure presets: resolution, aliasing, independence from profiles, the bundled Arimo font and legacy parity.

``standard`` is the canonical publication preset and ``bioinfo`` its backward-compatible alias.
The standard preset is applied per figure type (volcano and PCA in this release); every other
figure keeps the legacy appearance and the legacy font.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import pytest

from test_figures import run_r

PROFILES = ["general", "nature", "plos", "elsevier"]
FONT_DIR = Path(__file__).resolve().parents[1] / "src" / "rnaseq" / "fonts" / "arimo"


def ctx(tmp_path: Path, style: str, extra: str = "") -> dict:
    return run_r(f'''
cfg <- {style}
ctx <- nf_rna_figure_context(cfg)
{extra}
''', tmp_path)


def test_a_config_without_figure_style_is_legacy_general(tmp_path):
    result = ctx(tmp_path, "list()", 'cat(toJSON(list(req=ctx$requested_preset, app=ctx$preset, profile=ctx$profile$name, scoped=length(ctx$standard_figures)), auto_unbox=TRUE))')
    assert result == {"req": "legacy", "app": "legacy", "profile": "general", "scoped": 0}


@pytest.mark.parametrize("preset,canonical", [("legacy", "legacy"), ("standard", "standard"), ("bioinfo", "standard")])
@pytest.mark.parametrize("profile", PROFILES)
def test_preset_and_profile_are_independent(tmp_path, preset, canonical, profile):
    style = f'list(figure_style=list(preset="{preset}", profile="{profile}", source="explicit"))'
    result = ctx(tmp_path, style, 'cat(toJSON(list(req=ctx$requested_preset, given=ctx$requested_name, profile=ctx$profile$name, base=ctx$profile$base_pt, dpi=ctx$profile$raster_dpi), auto_unbox=TRUE))')
    assert (result["req"], result["given"], result["profile"]) == (canonical, preset, profile)
    expected = {"general": (8, 300), "nature": (7, 300), "plos": (9, 300), "elsevier": (8, 500)}[profile]
    assert (result["base"], result["dpi"]) == expected  # a selected profile is never ignored; Elsevier keeps 500 dpi


def test_legacy_with_a_journal_profile_is_applied_not_ignored(tmp_path):
    result = ctx(tmp_path, 'list(figure_style=list(preset="legacy", profile="nature"))', 'cat(toJSON(list(max=ctx$profile$max_width_mm, base=ctx$profile$base_pt), auto_unbox=TRUE))')
    assert result == {"max": 183, "base": 7}


@pytest.mark.parametrize("style", [
    'list(figure_style=list(preset="modern"))',
    'list(figure_style=list(preset="standard", profile="science"))',
    'list(figure_style=list(preset=c("legacy", "standard")))',
    'list(figure_style=list(preset=1))',
    'list(figure_style="standard")',
    'list(figure_style=list(profile="cell"))',
])
def test_unknown_preset_or_profile_stops_instead_of_defaulting(tmp_path, style):
    with pytest.raises(AssertionError, match="Execution halted|unknown figure|must be an object"):
        ctx(tmp_path, style)


def test_the_older_figures_profile_hook_still_works_and_must_agree(tmp_path):
    ok = ctx(tmp_path, 'list(figures=list(profile="plos"))', 'cat(toJSON(list(p=ctx$profile$name, r=ctx$requested_preset), auto_unbox=TRUE))')
    assert ok == {"p": "plos", "r": "legacy"}
    agree = ctx(tmp_path, 'list(figures=list(profile="plos"), figure_style=list(preset="legacy", profile="plos"))', 'cat(toJSON(list(p=ctx$profile$name), auto_unbox=TRUE))')
    assert agree == {"p": "plos"}
    with pytest.raises(AssertionError, match="conflicting figure profiles"):
        ctx(tmp_path, 'list(figures=list(profile="nature"), figure_style=list(preset="legacy", profile="plos"))')


def test_bioinfo_is_an_alias_that_normalises_to_standard(tmp_path):
    result = run_r('''
a <- nf_rna_read_figure_style(list(figure_style=list(preset="bioinfo", profile="nature")))
b <- nf_rna_read_figure_style(list(figure_style=list(preset="standard", profile="nature")))
cat(toJSON(list(a=a$preset, a_given=a$given, b=b$preset, b_given=b$given, aliases=as.list(NF_RNA_PRESET_ALIASES), names=names(NF_RNA_STYLE_PRESETS)), auto_unbox=TRUE))
''', tmp_path)
    assert (result["a"], result["a_given"], result["b"], result["b_given"]) == ("standard", "bioinfo", "standard", "standard")
    assert result["aliases"] == {"bioinfo": "standard"} and result["names"] == ["legacy", "standard"]  # the alias has no table of its own


def test_legacy_accessors_return_exactly_the_historic_values(tmp_path):
    result = ctx(tmp_path, "list()", '''
cat(toJSON(list(
  direction=identical(nf_rna_direction_colours(ctx), NF_RNA_DIRECTION_COLOURS),
  diverging=identical(nf_rna_diverging_ends(ctx), NF_RNA_DIVERGING_ENDS),
  sequential=identical(nf_rna_sequential_stops(ctx), NF_RNA_SEQUENTIAL_STOPS),
  missing=identical(nf_rna_na_colour(ctx), NF_RNA_NA_COLOUR),
  primary=identical(nf_rna_primary_colour(ctx), NF_RNA_OKABE_ITO[["blue"]]),
  cat3=identical(nf_rna_categorical_colors(c("a", "b", "c"), ctx), nf_rna_categorical_colors(c("a", "b", "c"))),
  batch=identical(nf_rna_batch_colors(c("a", "b"), ctx), nf_rna_categorical_colors(c("a", "b"), ctx)),
  shapes=identical(nf_rna_categorical_shapes(c("a", "b"), ctx), nf_rna_categorical_shapes(c("a", "b"))),
  palette_name=nf_rna_categorical_palette_name(3, ctx),
  alpha_is_na=is.na(nf_rna_point_alpha(ctx)),
  pca_pt=nf_rna_pca_label_pt(ctx) == ctx$profile$min_pt,
  width=nf_rna_layout_width_mm(ctx, "volcano", 180),
  legacy_hex=unname(nf_rna_categorical_colors(c("a", "b"), ctx))
), auto_unbox=TRUE))''')
    assert all(v is True for k, v in result.items() if k not in {"palette_name", "legacy_hex", "width"})
    assert result["palette_name"] == "Okabe-Ito" and result["legacy_hex"] == ["#0072B2", "#D55E00"] and result["width"] == 180


def test_standard_is_scoped_to_the_migrated_types_and_everything_else_stays_legacy(tmp_path):
    result = ctx(tmp_path, 'list(figure_style=list(preset="standard", profile="general"))', '''
volcano <- nf_rna_figure_scope(ctx, "volcano"); pca <- nf_rna_figure_scope(ctx, "pca")
heat <- nf_rna_figure_scope(ctx, "unmigrated_type")
cat(toJSON(list(
  figures=NF_RNA_STANDARD_FIGURES,
  unscoped_cat=unname(nf_rna_categorical_colors(c("a", "b", "c"), ctx)),
  heat_is_base=identical(heat, ctx), heat_font=heat$font$requested,
  pca_cat=unname(nf_rna_categorical_colors(c("a", "b", "c"), pca)),
  direction=as.list(nf_rna_direction_colours(volcano)),
  alpha=nf_rna_point_alpha(volcano), pca_pt=nf_rna_pca_label_pt(pca),
  font=pca$font$family, verified=pca$font$verified, fallback=pca$font$fallback_used, embedded=unlist(pca$font$embedded_fonts),
  preset=pca$preset, scoped=ctx$scoped, diverging=as.list(nf_rna_diverging_ends(pca)),
  counts=NF_RNA_STYLE_PRESETS$standard$volcano_legend_counts, shapes=nf_rna_categorical_shapes(c("a", "b"), pca)
), auto_unbox=TRUE))''')
    assert result["figures"] == ["volcano", "pca", "boxplot", "correlation", "heatmap", "dotplot"]
    assert result["unscoped_cat"] == ["#0072B2", "#D55E00", "#009E73"] and result["heat_is_base"] and result["heat_font"] == "Arial"
    assert result["pca_cat"] == ["#4E78A7", "#CC5A6A", "#59A14F"]
    assert result["direction"] == {"Up": "#B2182B", "Down": "#2166AC", "Not significant": "#BDBDBD"}
    assert (result["alpha"], result["pca_pt"], result["counts"]) == (0.85, 6, True)
    assert result["font"] == "Arimo" and result["verified"] is True and result["fallback"] is False and result["preset"] == "standard"
    assert sorted(result["embedded"]) == ["Arimo-Bold", "Arimo-BoldItalic", "Arimo-Italic", "Arimo-Regular"]
    assert result["scoped"] == ["volcano", "pca"] and result["diverging"]["low"] == "#4E78A7"
    assert result["shapes"] == [16, 17]  # PCA shape distinctions are preserved


def test_standard_values_are_declared_in_the_single_preset_table(tmp_path):
    result = run_r('''
b <- NF_RNA_STYLE_PRESETS$standard; l <- NF_RNA_STYLE_PRESETS$legacy
cat(toJSON(list(cat=unname(b$categorical), batch=unname(b$batch), direction=as.list(b$direction), diverging=as.list(b$diverging),
  alpha=b$point_alpha, pca=b$pca_label_pt, counts=b$volcano_legend_counts, font=b$font, layout=b$layout, has_tiff_dpi="tiff_dpi" %in% names(b),
  legacy_font=l$font, legacy_layout=is.null(l$layout), names=names(NF_RNA_STYLE_PRESETS)), auto_unbox=TRUE, null="null"))
''', tmp_path)
    assert result["cat"][:3] == ["#4E78A7", "#CC5A6A", "#59A14F"] and len(result["cat"]) == 8 and len(set(result["cat"])) == 8
    assert len(result["batch"]) >= 8 and not set(result["batch"]) & set(result["cat"][:3])
    assert result["direction"] == {"Up": "#B2182B", "Down": "#2166AC", "Not significant": "#BDBDBD"}
    assert result["diverging"] == {"low": "#4E78A7", "mid": "#F7F7F7", "high": "#CC5A6A"}
    assert (result["alpha"], result["pca"], result["counts"]) == (0.85, 6, True)
    assert result["has_tiff_dpi"] is False  # the profile always decides the raster dpi
    assert result["font"] == {"policy": "strict", "preferred": "Arimo", "fallback": None, "bundled": "arimo"}
    assert set(result["layout"]) == {"volcano", "pca", "boxplot", "correlation", "heatmap", "dotplot"} and result["legacy_layout"]
    assert result["legacy_font"] == {"policy": "preferred_with_recorded_fallback", "preferred": "Arial", "fallback": "DejaVu Sans"}


@pytest.mark.parametrize("profile,width", [("general", 85), ("nature", 89), ("plos", 85), ("elsevier", 90)])
def test_standard_figures_use_the_profile_single_column_and_never_exceed_it(tmp_path, profile, width):
    result = ctx(tmp_path, f'list(figure_style=list(preset="standard", profile="{profile}"))', '''
v <- nf_rna_figure_scope(ctx, "volcano")
cat(toJSON(list(volcano=nf_rna_layout_width_mm(v, "volcano", 999), pca=nf_rna_layout_width_mm(v, "pca", 999), aspect=nf_rna_layout_aspect(v, "volcano")), auto_unbox=TRUE))''')
    assert result["volcano"] == width and result["pca"] == width  # plos has no single column: the standard preset's compact width (never above the profile maximum)
    assert result["aspect"] == 1.0


def test_legacy_font_resolution_is_unchanged(tmp_path):
    result = run_r('f <- nf_rna_resolve_font(); cat(toJSON(list(req=f$requested, fam=f$family, fb=f$fallback_used, pref=NF_RNA_FONT_PREFERRED, fallback=NF_RNA_FONT_FALLBACK), auto_unbox=TRUE))', tmp_path)
    assert result["req"] == "Arial" and result["pref"] == "Arial" and result["fallback"] == "DejaVu Sans"
    assert result["fam"] in {"Arial", "DejaVu Sans"} and result["fb"] == (result["fam"] != "Arial")


def test_a_legacy_run_never_touches_fontconfig_and_standard_does_not_change_legacy_fonts(tmp_path):
    result = run_r('''
before <- Sys.getenv("FONTCONFIG_FILE", unset="")
legacy <- nf_rna_figure_context(list(figure_style=list(preset="legacy", profile="general")))
after_legacy <- Sys.getenv("FONTCONFIG_FILE", unset="")
standard <- nf_rna_figure_context(list(figure_style=list(preset="standard", profile="general")))
cat(toJSON(list(untouched=identical(before, after_legacy), activated=Sys.getenv("FONTCONFIG_FILE") != before,
  same_legacy_font=identical(legacy$font, standard$font), legacy_family=legacy$font$family), auto_unbox=TRUE))
''', tmp_path)
    assert result["untouched"] and result["activated"] and result["same_legacy_font"] and result["legacy_family"] in {"Arial", "DejaVu Sans"}


# ------------------------------------------------------------------ bundled Arimo


def test_bundled_font_files_match_the_pinned_manifest_and_license_is_present():
    manifest = json.loads((FONT_DIR / "FONT_MANIFEST.json").read_text(encoding="utf-8"))
    assert manifest["family"] == "Arimo" and manifest["license"] == "SIL Open Font License 1.1" and manifest["reserved_font_name"] is None
    assert manifest["source"]["repository"] == "https://github.com/googlefonts/arimo" and len(manifest["source"]["commit"]) == 40
    assert set(manifest["faces"]) == {"regular", "bold", "italic", "bold_italic"}
    for name, entry in manifest["files"].items():
        data = (FONT_DIR / name).read_bytes()
        assert hashlib.sha256(data).hexdigest() == entry["sha256"] and len(data) == entry["bytes"]
    license_text = (FONT_DIR / "OFL.txt").read_text(encoding="utf-8")
    assert "SIL OPEN FONT LICENSE Version 1.1" in license_text and "Arimo Project Authors" in license_text
    assert sorted(path.name for path in FONT_DIR.iterdir()) == sorted([*manifest["files"], "FONT_MANIFEST.json"])  # only what is needed
    assert not any("[" in path.name for path in FONT_DIR.iterdir())  # static fonts, not variable


def test_the_bundled_font_is_verified_activated_and_embedded_under_its_own_name(tmp_path):
    result = run_r('''
bundle <- nf_rna_verify_bundled_font("arimo")
nf_rna_activate_bundled_font("arimo")
font <- nf_rna_require_strict_font("Arimo", "arimo")
cat(toJSON(list(dir=bundle$dir, font=font), auto_unbox=TRUE, null="null"))
''', tmp_path)
    font = result["font"]
    assert font["family"] == "Arimo" and font["verified"] is True and font["fallback_used"] is False and font["variable"] is False
    assert sorted(font["embedded_fonts"]) == ["Arimo-Bold", "Arimo-BoldItalic", "Arimo-Italic", "Arimo-Regular"]
    assert Path(font["file"]).resolve() == (FONT_DIR / "Arimo-Regular.ttf").resolve()
    manifest = json.loads((FONT_DIR / "FONT_MANIFEST.json").read_text(encoding="utf-8"))
    assert font["sha256"] == manifest["files"]["Arimo-Regular.ttf"]["sha256"] == hashlib.sha256((FONT_DIR / "Arimo-Regular.ttf").read_bytes()).hexdigest()
    assert font["bundled"]["license"] == "SIL Open Font License 1.1" and font["bundled"]["source_commit"] == manifest["source"]["commit"]


def test_a_modified_or_missing_bundled_file_stops_the_run(tmp_path):
    import shutil

    copy = tmp_path / "fonts" / "arimo"
    shutil.copytree(FONT_DIR, copy)
    (copy / "Arimo-Bold.ttf").write_bytes((copy / "Arimo-Bold.ttf").read_bytes() + b"x")
    result = run_r(f'''
nf_rna_bundled_font_dir <- function(name) "{copy}"
tampered <- tryCatch({{ nf_rna_verify_bundled_font("arimo"); "no error" }}, error=function(e) conditionMessage(e))
file.copy(file.path("{FONT_DIR}", "Arimo-Bold.ttf"), file.path("{copy}", "Arimo-Bold.ttf"), overwrite=TRUE)
file.remove(file.path("{copy}", "Arimo-Italic.ttf"))
missing <- tryCatch({{ nf_rna_verify_bundled_font("arimo"); "no error" }}, error=function(e) conditionMessage(e))
cat(toJSON(list(tampered=tampered, missing=missing), auto_unbox=TRUE))
''', tmp_path)
    assert "does not match its pinned SHA-256" in result["tampered"] and "Arimo-Bold.ttf" in result["tampered"]
    assert "bundled font file is missing" in result["missing"] and "never substituted" in result["missing"]


def test_another_installed_arimo_build_is_hidden_so_only_the_bundled_files_are_drawn(tmp_path):
    other = tmp_path / "other"
    other.mkdir()
    (other / "Arimo-Regular.ttf").write_bytes((FONT_DIR / "Arimo-Regular.ttf").read_bytes())
    result = run_r(f'''
conf <- file.path(getwd(), "extra.conf")
writeLines(c('<?xml version="1.0"?>', '<fontconfig>', '<include ignore_missing="yes">/etc/fonts/fonts.conf</include>', '<dir>{other}</dir>', '</fontconfig>'), conf)
Sys.setenv(FONTCONFIG_FILE=conf)
font <- nf_rna_require_strict_font("Arimo", "arimo")
paths <- unique(normalizePath(systemfonts::system_fonts()$path[systemfonts::system_fonts()$family == "Arimo"]))
cat(toJSON(list(file=font$file, paths=paths), auto_unbox=FALSE))
''', tmp_path)
    assert all(str(FONT_DIR.resolve()) in path for path in result["paths"]) and len(result["paths"]) == 4


def test_fontconfig_initialised_before_activation_is_detected_not_ignored(tmp_path):
    result = run_r('''
invisible(systemfonts::system_fonts())  # fontconfig is now initialised without the bundled directory
msg <- tryCatch({ nf_rna_require_strict_font("Arimo", "arimo"); "no error" }, error=function(e) conditionMessage(e))
cat(toJSON(list(msg=msg), auto_unbox=TRUE))
''', tmp_path)
    assert result["msg"] != "no error" and "Arimo" in result["msg"] and "never substituted" in result["msg"]


def test_strict_font_report_never_substitutes_another_family(tmp_path):
    result = run_r('''
missing <- nf_rna_strict_font_report("NoSuchFamilyForNfRnaTests")
installed <- systemfonts::system_fonts()$family[[1]]
present <- nf_rna_strict_font_report(installed)
cat(toJSON(list(missing=missing, present_family=present$family, requested=installed, present_ok=present$available && file.exists(present$file) && grepl("^[0-9a-f]{64}$", present$sha256), file=present$file, sha=present$sha256), auto_unbox=TRUE, null="null"))
''', tmp_path)
    assert result["missing"]["available"] is False and result["missing"]["family"] is None
    assert "NoSuchFamilyForNfRnaTests" in result["missing"]["requirement"] and "never substituted" in result["missing"]["requirement"]
    assert result["present_family"] == result["requested"] and result["present_ok"]
    assert result["sha"] == hashlib.sha256(Path(result["file"]).read_bytes()).hexdigest()  # SHA-256, not md5


def test_drawn_font_check_rejects_fontconfig_substitution(tmp_path):
    """An uninstalled family is mapped by fontconfig to some other font; only the embedded name reveals it."""

    result = run_r('''
fake <- nf_rna_verify_drawn_font("NoSuchFontNfRna")
cat(toJSON(list(fake=fake, have_dejavu="DejaVu Sans" %in% systemfonts::system_fonts()$family,
  dejavu=if ("DejaVu Sans" %in% systemfonts::system_fonts()$family) nf_rna_verify_drawn_font("DejaVu Sans") else NULL), auto_unbox=TRUE, null="null"))
''', tmp_path)
    assert result["fake"]["matched"] is False and result["fake"]["embedded_fonts"]
    if result["have_dejavu"]:
        assert result["dejavu"]["matched"] is True and set(result["dejavu"]["embedded_fonts"]) <= {"DejaVuSans", "DejaVuSans-Bold", "DejaVuSans-Oblique", "DejaVuSans-BoldOblique"}


def test_pdf_font_name_reader_strips_the_subset_tag(tmp_path):
    result = run_r('''
f <- file.path(getwd(), "probe.pdf")
grDevices::cairo_pdf(f, width=2, height=1, family="sans"); grid::grid.text("abc"); invisible(dev.off())
cat(toJSON(list(names=nf_rna_pdf_font_names(f)), auto_unbox=TRUE))
''', tmp_path)
    assert result["names"] and all("+" not in name for name in (result["names"] if isinstance(result["names"], list) else [result["names"]]))


def test_strict_font_is_required_with_an_actionable_error_when_missing(tmp_path):
    result = run_r('''
msg <- tryCatch({ nf_rna_require_strict_font("NoSuchFontNfRna"); "no error" }, error=function(e) conditionMessage(e))
cat(toJSON(list(msg=msg), auto_unbox=TRUE))
''', tmp_path)
    assert "Required figure font 'NoSuchFontNfRna' is not available" in result["msg"]
    assert "never substituted" in result["msg"]


# ------------------------------------------------------------------ manifest


def _manifest(tmp_path: Path, style: str, extra: str = "") -> dict:
    ctx(tmp_path, style, f'{extra}\ndir.create("out"); nf_rna_write_figure_manifest(ctx, "out"); cat(toJSON(list(ok=TRUE), auto_unbox=TRUE))')
    return json.loads((tmp_path / "out" / "figure_manifest.json").read_text(encoding="utf-8"))


def test_legacy_manifest_has_no_new_keys(tmp_path):
    for style in ("list()", 'list(figure_style=list(preset="legacy", profile="general", source="legacy_default"))'):
        manifest = _manifest(tmp_path, style)
        assert "figure_style" not in manifest
        assert manifest["schema_version"] == "nf-rna.figure-manifest.v1"
        assert set(manifest) == {"schema_version", "profile", "font", "typography", "colour_policy", "devices", "r_version", "package_versions", "figures"}


def test_a_standard_request_is_disclosed_in_the_manifest_with_the_verified_font(tmp_path):
    manifest = _manifest(tmp_path, 'list(figure_style=list(preset="bioinfo", profile="elsevier"))', 'invisible(nf_rna_figure_scope(ctx, "pca"))')
    style = manifest["figure_style"]
    assert style["requested_preset"] == "bioinfo" and style["preset"] == "standard" and style["profile"] == "elsevier"
    assert style["standard_figure_types"] == ["volcano", "pca", "boxplot", "correlation", "heatmap", "dotplot"] and style["standard_figures_drawn"] == ["pca"]
    assert style["font"]["family"] == "Arimo" and len(style["font"]["sha256"]) == 64 and style["font"]["verified"] is True
    assert manifest["profile"]["raster_dpi"] == 500  # the Elsevier requirement survives the preset
    assert manifest["font"]["family"] == "Arimo" and manifest["font"]["sha256"] == style["font"]["sha256"]  # the manifest-level font is the one actually drawn
    assert manifest["colour_policy"]["applied_preset"] == "standard" and manifest["colour_policy"]["direction"]["Up"] == "#B2182B"
    assert manifest["schema_version"] == "nf-rna.figure-manifest.v1"


def test_a_scoped_figure_is_recorded_on_the_run_and_embeds_arimo(tmp_path):
    result = ctx(tmp_path, 'list(figure_style=list(preset="standard", profile="general"))', '''
library(ggplot2)
v <- nf_rna_figure_scope(ctx, "volcano")
p <- ggplot(data.frame(x=1:3, y=1:3), aes(x, y)) + geom_point() + labs(x="Higher in Treatment", y="-log10(adjusted P, BH)", title="Bold title") + nf_rna_theme(v)
dir.create("out")
w <- nf_rna_layout_width_mm(v, "volcano", 180)
rec <- nf_rna_save_figure(v, p, file.path("out", "volcano"), w, w, "out")
cat(toJSON(list(n=length(ctx$records), own=length(v$records), id=ctx$records[[1]]$id, style=ctx$records[[1]]$figure_style, mm=rec$width_mm,
  pdf=nf_rna_pdf_font_names("out/volcano.pdf")), auto_unbox=TRUE, null="null"))''')
    assert result["n"] == 1 and result["own"] == 0 and result["id"] == "volcano"
    assert result["style"]["preset"] == "standard" and result["style"]["font"] == "Arimo"
    assert 84.5 < result["mm"] < 85.5  # snapped to whole points; a single column, not 180 mm
    assert all(name.startswith("Arimo") for name in result["pdf"]) and "Arimo-Regular" in result["pdf"] and "Arimo-Bold" in result["pdf"]
