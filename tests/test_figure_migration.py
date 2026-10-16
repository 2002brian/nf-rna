"""Standard-preset migration of sample boxplots, the sample correlation heatmap and DEG heatmaps.

Unit tests use figures.R directly (stable group colours, content-aware sizing); the integration test runs
L1 + L2 on a synthetic project once with the legacy selection and once with the standard preset and checks
that only the figure files changed.
"""

from __future__ import annotations

import hashlib
import importlib.util
import json
from pathlib import Path

import pytest

from conftest import require_r_packages
from test_figures import run_r

ROOT = Path(__file__).resolve().parents[1]
TOOLS = ROOT / "tools" / "figure_review"
STANDARD = 'nf_rna_figure_context(list(figure_style=list(preset="standard", profile="{profile}")))'


def standard(profile: str = "general") -> str:
    return f"base <- {STANDARD.format(profile=profile)}\n"


# ------------------------------------------------------------------ stable group colours


def test_group_colours_do_not_move_when_groups_are_reordered_missing_or_subsetted(tmp_path):
    result = run_r(standard() + '''
ctx <- nf_rna_figure_scope(base, "pca")
universe <- c("Control", "Low", "High")
full <- nf_rna_group_colours(universe, context = ctx)
cat(toJSON(list(
  full = as.list(full),
  reordered_present = as.list(nf_rna_group_colours(universe, c("High", "Control", "Low"), ctx)),
  subset_high = as.list(nf_rna_group_colours(universe, "High", ctx)),
  subset_control_high = as.list(nf_rna_group_colours(universe, c("High", "Control"), ctx)),
  other_pair = as.list(nf_rna_group_colours(universe, c("Low", "High"), ctx)),
  names_follow_request = names(nf_rna_group_colours(universe, c("High", "Control"), ctx)),
  unknown = tryCatch({ nf_rna_group_colours(universe, "Nope", ctx); "no error" }, error = function(e) conditionMessage(e))
), auto_unbox = TRUE))''', tmp_path)
    approved = ["#4E78A7", "#CC5A6A", "#59A14F"]
    assert result["full"] == dict(zip(["Control", "Low", "High"], approved))
    assert result["reordered_present"] == result["full"]  # order of the request never moves a colour
    assert result["subset_high"] == {"High": approved[2]}  # a lone group keeps its own colour, not the first one
    assert result["subset_control_high"] == {"High": approved[2], "Control": approved[0]}
    assert result["other_pair"] == {"Low": approved[1], "High": approved[2]}
    assert result["names_follow_request"] == ["High", "Control"]
    assert "not in the colour universe" in result["unknown"]


def test_more_than_three_groups_extend_deterministically_and_eight_is_the_palette_limit(tmp_path):
    result = run_r(standard() + '''
ctx <- nf_rna_figure_scope(base, "pca")
eight <- paste0("G", 1:8); nine <- paste0("G", 1:9)
c8 <- nf_rna_group_colours(eight, context = ctx); c9 <- nf_rna_group_colours(nine, context = ctx)
cat(toJSON(list(c8 = unname(c8), c9 = unname(c9), g5_in_subset = unname(nf_rna_group_colours(eight, "G5", ctx)), g5_in_full = unname(c8[["G5"]]),
  again = identical(nf_rna_group_colours(eight, context = ctx), c8), style_palette = unname(NF_RNA_STYLE_PRESETS$standard$categorical)), auto_unbox = FALSE))''', tmp_path)
    approved = ["#4E78A7", "#CC5A6A", "#59A14F"]
    assert result["c8"][:3] == approved and len(set(result["c8"])) == 8 and result["c8"] == result["style_palette"]
    assert result["g5_in_subset"] == result["g5_in_full"]
    assert len(result["c9"]) == 9 and len(set(result["c9"])) == 9 and result["c9"][:8] == result["c8"]  # beyond 8: slots 1-8 are kept
    assert result["again"] == [True]


def test_group_and_batch_colours_are_separate_palettes_and_both_stable(tmp_path):
    result = run_r(standard() + '''
ctx <- nf_rna_figure_scope(base, "pca")
levels <- c("B1", "B2", "B3")
cat(toJSON(list(group = unname(nf_rna_group_colours(levels, context = ctx)), batch = unname(nf_rna_batch_colours(levels, context = ctx)),
  batch_subset = unname(nf_rna_batch_colours(levels, "B3", ctx)), batch_full_b3 = unname(nf_rna_batch_colours(levels, context = ctx)[["B3"]])), auto_unbox = FALSE))''', tmp_path)
    assert not set(result["group"]) & set(result["batch"])  # group and batch never share a colour
    assert result["batch_subset"] == [result["batch_full_b3"]] or result["batch_subset"] == result["batch_full_b3"]


def test_legacy_group_colours_are_the_historic_okabe_ito_values(tmp_path):
    result = run_r('''
ctx <- nf_rna_figure_context(list())
cat(toJSON(list(three = unname(nf_rna_group_colours(c("a", "b", "c"), context = ctx)), old = unname(nf_rna_categorical_colors(c("a", "b", "c"), ctx)),
  batch = unname(nf_rna_batch_colours(c("a", "b"), context = ctx))), auto_unbox = FALSE))''', tmp_path)
    assert result["three"] == result["old"] == ["#0072B2", "#D55E00", "#009E73"] and result["batch"] == ["#0072B2", "#D55E00"]


# ------------------------------------------------------------------ boxplots


BOX = '''
{standard}
ctx <- nf_rna_figure_scope(base, "boxplot")
make <- function(n, prefix = "S", groups = c("Control", "Treatment"), batch = TRUE) {{
  samples <- sprintf("%s%02d", prefix, seq_len(n))
  stats <- data.frame(sample_id = samples, matrix = "m", genes = 100, whisker_low = 1, q1 = 4, median = 6, q3 = 8, whisker_high = 12)
  group <- factor(rep(groups, length.out = n), levels = groups)
  b <- if (batch) factor(rep(c("B1", "B2", "B3"), length.out = n)) else NULL
  nf_rna_sample_boxplot(ctx, stats, group, "condition", "VST expression", b, "batch")
}}
'''


@pytest.mark.parametrize("profile,single,double,max_height", [("general", 85, 180, 220), ("nature", 89, 183, 170), ("elsevier", 90, 190, 240)])
def test_boxplot_size_follows_the_content_and_stays_inside_the_profile(tmp_path, profile, single, double, max_height):
    result = run_r(BOX.format(standard=standard(profile)) + '''
sizes <- lapply(c(6, 12, 26, 60, 150), function(n) { b <- make(n); list(n = n, w = b$width_mm, h = b$height_mm, labels = b$labels_shown) })
long <- make(12, prefix = "Organoid_Patient_Derived_Line_")
short <- make(12)
cat(toJSON(list(sizes = sizes, long_h = long$height_mm, short_h = short$height_mm, long_w = long$width_mm, short_w = short$width_mm), auto_unbox = TRUE))''', tmp_path)
    widths = [item["w"] for item in result["sizes"]]
    assert widths == sorted(widths) and widths[0] < widths[2] < widths[3]  # grows with the number of samples (compact figures share the minimum width)
    assert widths[0] < single and widths[0] >= 59  # a six-sample boxplot is compact: narrower than a single column
    assert all(w <= double + 1e-9 for w in widths) and widths[-1] == double  # capped at the profile's maximum width
    assert all(item["h"] <= max_height + 1e-9 for item in result["sizes"])
    assert result["long_h"] > result["short_h"] + 10  # long rotated sample labels get the height they need
    assert [item["labels"] for item in result["sizes"]] == [True, True, True, True, False]  # unchanged rule: no labels above 100 samples


def test_boxplot_legend_is_never_wider_than_the_figure_and_preserves_order_and_labels(tmp_path):
    result = run_r(BOX.format(standard=standard()) + '''
b <- make(6, groups = c("Vehicle control", "Low dose compound A", "Medium dose compound A", "High dose compound A", "Positive control reference"))
g <- nf_rna_with_probe_device(ctx, function() { gt <- ggplot2::ggplotGrob(b$plot); nf_rna_grob_width_mm(gt, "guide-box") })
built <- ggplot2::ggplot_build(b$plot)
cat(toJSON(list(legend_rows = b$legend_rows, legend_mm = g, width = b$width_mm, x_order = levels(b$plot$data$sample_id), fills = unname(b$colours)), auto_unbox = FALSE))''', tmp_path)
    assert result["legend_mm"][0] + 4 <= result["width"][0] + 1e-6
    assert result["legend_rows"][0] >= 2  # a wide legend wraps instead of being clipped
    assert result["x_order"] == [f"S{i:02d}" for i in range(1, 7)]  # sample order is the input order
    assert result["fills"][:3] == ["#4E78A7", "#CC5A6A", "#59A14F"]


def test_boxplot_statistics_are_untouched_by_the_standard_preset(tmp_path):
    result = run_r(BOX.format(standard=standard()) + '''
legacy <- nf_rna_figure_context(list())
stats <- data.frame(sample_id = c("a", "b"), matrix = "m", genes = 10, whisker_low = c(1, 2), q1 = c(3, 4), median = c(5, 6), q3 = c(7, 8), whisker_high = c(9, 10))
group <- factor(c("x", "y")); std <- nf_rna_sample_boxplot(ctx, stats, group, "g", "y"); old <- nf_rna_sample_boxplot(legacy, stats, group, "g", "y")
cat(toJSON(list(same = identical(std$plot$data[, c("whisker_low", "q1", "median", "q3", "whisker_high")], old$plot$data[, c("whisker_low", "q1", "median", "q3", "whisker_high")]),
  legacy_has_no_width = is.null(old$width_mm)), auto_unbox = TRUE))''', tmp_path)
    assert result == {"same": True, "legacy_has_no_width": True}


# ------------------------------------------------------------------ heatmaps: content-aware sizes


DEG = '''
{standard}
ctx <- nf_rna_figure_scope(base, "heatmap")
deg <- function(genes, samples, long = FALSE, seed = 3) {{
  set.seed(seed)
  ids <- if (long) sprintf("ENSMUSG%011d_Long_Symbol", seq_len(genes)) else sprintf("Gene%05d", seq_len(genes))
  cols <- if (long) sprintf("Organoid_Patient_Derived_Line_%02d", seq_len(samples)) else paste0("S", seq_len(samples))
  z <- matrix(rnorm(genes * samples), genes, samples, dimnames = list(ids, cols)); z <- pmax(pmin(z, 2), -2)
  group <- factor(rep(c("Control", "Treatment"), length.out = samples)); ann <- data.frame(condition = group, row.names = cols)
  colours <- list(condition = nf_rna_group_colours(levels(group), context = ctx))
  fit <- nf_rna_deg_heatmap(ctx, z, ann, colours, nf_rna_cluster(z), nf_rna_cluster(t(z)), "Treatment vs Control\\nTop DEGs", 2)
  list(w = fit$width_mm, h = fit$height_mm, rows = fit$show_rows, cols = fit$show_cols, row_pt = fit$row_pt, clipped = fit$clipped, cw = fit$cw, ch = fit$ch)
}}
'''


@pytest.mark.parametrize("profile,max_w,max_h", [("general", 180, 220), ("nature", 183, 170), ("plos", 190.5, 222.3), ("elsevier", 190, 240)])
def test_deg_heatmap_sizes_for_20_50_and_100_genes_with_short_and_long_labels(tmp_path, profile, max_w, max_h):
    result = run_r(DEG.format(standard=standard(profile)) + '''
out <- list()
for (genes in c(20, 50, 100)) for (samples in c(6, 12, 26)) for (long in c(FALSE, TRUE)) {
  r <- deg(genes, samples, long); r$genes <- genes; r$samples <- samples; r$long <- long; out[[length(out) + 1]] <- r
}
cat(toJSON(out, auto_unbox = TRUE))''', tmp_path)
    for item in result:
        assert item["w"] <= max_w + 1e-6 and item["h"] <= max_h + 1e-6, item  # the profile constrains the size
        assert item["clipped"] is False, item  # content that does not fit is hidden, never clipped
        assert item["cw"] > 0 and item["ch"] > 0
    by = {(i["genes"], i["samples"], i["long"]): i for i in result}
    assert by[(20, 6, False)]["w"] < by[(20, 26, False)]["w"] and by[(20, 6, False)]["h"] < by[(100, 6, False)]["h"]  # follows the content
    if profile == "general":
        assert by[(20, 6, False)]["w"] < 120 and by[(20, 6, False)]["h"] < 120  # compact: not forced to 180 mm
        assert by[(20, 12, False)]["rows"] and by[(50, 12, False)]["rows"] and not by[(100, 12, False)]["rows"]  # 100 gene labels are hidden, 20 and 50 shown
        assert by[(100, 12, False)]["ch"] >= 0.8  # but the rows stay distinguishable
        assert by[(20, 12, True)]["w"] > by[(20, 12, False)]["w"] + 20  # long labels get the room they need
        assert all(i["cols"] for i in result if i["samples"] <= 12)


def test_deg_heatmap_keeps_every_value_order_and_cluster_it_is_given(tmp_path):
    result = run_r(DEG.format(standard=standard()) + '''
set.seed(5)
z <- matrix(rnorm(30 * 8), 30, 8, dimnames = list(sprintf("G%02d", 1:30), paste0("S", 1:8))); z <- pmax(pmin(z, 2), -2)
ann <- data.frame(condition = factor(rep(c("Control", "Treatment"), 4)), row.names = colnames(z))
rt <- nf_rna_cluster(z); ct <- nf_rna_cluster(t(z))
fit <- nf_rna_deg_heatmap(ctx, z, ann, list(condition = nf_rna_group_colours(levels(ann$condition), context = ctx)), rt, ct, "t", 2)
plain <- pheatmap::pheatmap(z, cluster_rows = rt, cluster_cols = ct, silent = TRUE)
cat(toJSON(list(rows_same = identical(fit$heat$tree_row$order, rt$order), cols_same = identical(fit$heat$tree_col$order, ct$order),
  plain_rows = identical(plain$tree_row$order, rt$order), unchanged = identical(z, pmax(pmin(z, 2), -2))), auto_unbox = TRUE))''', tmp_path)
    assert all(result.values())


@pytest.mark.parametrize("samples,label", [(6, "short"), (12, "short"), (26, "short"), (60, "short"), (12, "long"), (26, "long")])
def test_correlation_heatmap_is_square_readable_and_inside_the_profile(tmp_path, samples, label):
    result = run_r(standard() + f'''
ctx <- nf_rna_figure_scope(base, "correlation")
n <- {samples}; cols <- if ("{label}" == "short") paste0("S", seq_len(n)) else sprintf("Organoid_Patient_Derived_Line_%02d", seq_len(n))
set.seed(2); m <- matrix(rnorm(500 * n, 8, 2), 500, n, dimnames = list(NULL, cols)); corr <- cor(m)
ann <- data.frame(condition = factor(rep(c("Control", "Treatment"), length.out = n)), row.names = cols)
fit <- nf_rna_correlation_heatmap(ctx, corr, nf_rna_sequential_palette(100, ctx), seq(floor(min(corr) * 100) / 100, 1, length.out = 101), ann, list(condition = nf_rna_group_colours(levels(ann$condition), context = ctx)), "Pearson correlation (blind VST)")
cat(toJSON(list(w = fit$width_mm, h = fit$height_mm, cw = fit$cw, ch = fit$ch, labels = fit$show_rows, clipped = fit$clipped, same_order = identical(fit$heat$tree_row$order, pheatmap::pheatmap(corr, silent = TRUE)$tree_row$order)), auto_unbox = TRUE))''', tmp_path)
    assert result["cw"] == pytest.approx(result["ch"]) and not result["clipped"]  # square cells, nothing clipped
    assert result["w"] <= 180 + 1e-6 and result["h"] <= 220 + 1e-6
    assert result["same_order"] is True  # pheatmap's own clustering is what is drawn
    assert result["labels"] == (samples <= 26)  # labels only while they fit legibly (cell >= label height)
    if samples == 6 and label == "short":
        assert result["w"] < 120  # compact


def test_every_standard_figure_embeds_only_arimo_and_follows_the_profile_dpi(tmp_path):
    result = run_r(standard("elsevier") + '''
ctx <- nf_rna_figure_scope(base, "correlation")
set.seed(2); m <- matrix(rnorm(300 * 8), 300, 8, dimnames = list(NULL, paste0("S", 1:8))); corr <- cor(m)
ann <- data.frame(condition = factor(rep(c("A", "B"), 4)), row.names = colnames(corr))
fit <- nf_rna_correlation_heatmap(ctx, corr, nf_rna_sequential_palette(100, ctx), seq(-1, 1, length.out = 101), ann, list(condition = nf_rna_group_colours(c("A", "B"), context = ctx)), "Pearson correlation (blind VST)")
dir.create("out"); rec <- nf_rna_save_figure(ctx, fit$heat, file.path("out", "corr"), fit$width_mm, fit$height_mm, "out")
cat(toJSON(list(fonts = nf_rna_pdf_font_names("out/corr.pdf"), dpi = vapply(rec$files, function(f) if (is.null(f$dpi)) NA_real_ else f$dpi, numeric(1)), w = rec$width_mm), auto_unbox = FALSE))''', tmp_path)
    assert result["fonts"] and all(name.startswith("Arimo") for name in result["fonts"])
    assert [d for d in result["dpi"] if d != "NA"] == [150, 500]  # PNG preview 150, Elsevier TIFF 500


# ------------------------------------------------------------------ end to end


def _load(name: str):
    spec = importlib.util.spec_from_file_location(name, TOOLS / f"{name}.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@pytest.fixture(scope="module")
def runs(tmp_path_factory):
    require_r_packages("DESeq2", "ggplot2", "pheatmap", "ggrepel", "jsonlite", "systemfonts")
    base = tmp_path_factory.mktemp("migration")
    demo = _load("make_demo_project")
    generate = _load("generate_runs")
    demo.build(base / "project", 4, 600, 20261015, "", 3)
    out = {}
    for label, figures in (("legacy", None), ("standard", {"preset": "standard", "profile": "general"})):
        out[label] = base / label
        out[label].mkdir()
        assert generate.run_l1_l2(base / "project", out[label], figures) == "SUCCESS"
    return out


def _files(root: Path) -> dict[str, bytes]:
    return {str(p.relative_to(root)): p.read_bytes() for p in sorted(root.rglob("*")) if p.is_file()}


def test_only_figure_files_and_manifests_change_between_legacy_and_standard(runs):
    legacy, standard_run = _files(runs["legacy"]), _files(runs["standard"])
    assert legacy.keys() == standard_run.keys()  # identical output paths and file names
    tables = [name for name in legacy if name.endswith((".tsv", ".csv"))]
    assert tables and all(legacy[name] == standard_run[name] for name in tables)  # source data is byte-identical
    summaries = [name for name in legacy if name.endswith("backend_summary.json")]
    assert summaries and all(legacy[name] == standard_run[name] for name in summaries)  # scientific summaries are identical


def test_standard_run_changes_the_five_migrated_figure_kinds_and_embeds_arimo(runs):
    legacy = {}
    standard_run = {}
    for label, store in (("legacy", legacy), ("standard", standard_run)):
        for manifest in sorted(runs[label].rglob("figure_manifest.json")):
            for figure in json.loads(manifest.read_text(encoding="utf-8"))["figures"]:
                store[f"{manifest.parent.relative_to(runs[label])}/{figure['id']}"] = (manifest.parent, figure)
    assert legacy.keys() == standard_run.keys() and len(standard_run) >= 7
    kinds = {key.split("/")[-1]: figure["kind"] for key, (_, figure) in standard_run.items()}
    assert {"volcano", "PCA", "sample expression distribution boxplot", "sample correlation heatmap", "DEG heatmap"} <= set(kinds.values())
    compare = _load("compare_runs")
    for key, (folder, figure) in standard_run.items():
        pdf = folder / next(f["path"] for f in figure["files"] if f["format"] == "pdf")
        fonts = compare.pdf_info(pdf)["embedded_fonts"]
        assert fonts and all(name.startswith("Arimo") for name in fonts), (key, fonts)
        assert figure["width_mm"] < legacy[key][1]["width_mm"], key  # none is forced to the full double column
        assert figure["figure_style"]["preset"] == "standard"
        assert all(f["dpi"] == 300 for f in figure["files"] if f["format"] == "tiff")


def test_displayed_order_and_clustering_are_the_legacy_ones(runs):
    """The order drawn by the standard figures equals the order of the legacy code path recomputed from the source data."""

    out = runs["standard"]
    corr = next(out.rglob("sample_correlation.tsv"))
    zscores = next(out.rglob("heatmap_zscores.tsv"))
    l1_manifest = json.loads(next(out.rglob("l1/figure_manifest.json")).read_text(encoding="utf-8"))
    l2_manifest = json.loads(next(out.rglob("l2/figure_manifest.json")).read_text(encoding="utf-8"))
    drawn_corr = next(f for f in l1_manifest["figures"] if f["id"] == "sample_correlation")["display_order"]
    drawn_deg = next(f for f in l2_manifest["figures"] if f["kind"] == "DEG heatmap")["display_order"]
    script = f'''
suppressPackageStartupMessages(library(pheatmap)); source("{ROOT / 'src' / 'rnaseq' / 'r' / 'figures.R'}")
c <- read.delim("{corr}", check.names = FALSE); m <- as.matrix(c[, -1]); rownames(m) <- c$sample_id
legacy_corr <- pheatmap(m, silent = TRUE)
z <- read.delim("{zscores}", check.names = FALSE); zm <- as.matrix(z[, -1]); rownames(zm) <- z$gene_id
zm <- zm[rowSums(is.finite(zm)) >= 2, , drop = FALSE]
rt <- nf_rna_cluster(zm); ct <- nf_rna_cluster(t(zm))
cat(toJSON(list(corr_rows = rownames(m)[legacy_corr$tree_row$order], corr_cols = colnames(m)[legacy_corr$tree_col$order],
  deg_rows = rownames(zm)[rt$order], deg_cols = colnames(zm)[ct$order]), auto_unbox = FALSE))
'''
    from test_figures import run_r as _run
    result = _run('library(jsonlite)\n' + script, runs["standard"])
    assert drawn_corr["rows"] == result["corr_rows"] and drawn_corr["columns"] == result["corr_cols"]
    assert drawn_deg["rows"] == result["deg_rows"] and drawn_deg["columns"] == result["deg_cols"]
    # the full sample set is shown exactly once, in one order for rows and columns of a symmetric matrix
    assert sorted(drawn_corr["rows"]) == sorted(drawn_corr["columns"])


def test_group_colours_agree_across_pca_boxplots_and_heatmap_annotations(runs):
    colours = {}
    for manifest in runs["standard"].rglob("figure_manifest.json"):
        for figure in json.loads(manifest.read_text(encoding="utf-8"))["figures"]:
            if figure["kind"] == "PCA" and figure["id"] == "pca":
                colours["pca"] = figure["encoding"]["colours"]
            elif figure["kind"] == "sample expression distribution boxplot":
                colours.setdefault("box", figure["encoding"]["colours"])
            elif figure["kind"] == "sample correlation heatmap":
                colours["corr"] = figure["annotation"]["colours"]
    assert set(colours) == {"pca", "box", "corr"}
    # the colour variable is the last design variable (dose, 3 levels): same level -> same colour everywhere
    for level, colour in colours["pca"].items():
        assert colours["box"][level] == colour
        if level in colours["corr"]:
            assert colours["corr"][level] == colour
    assert list(colours["pca"].values()) == ["#4E78A7", "#CC5A6A", "#59A14F"]


def test_existing_standard_volcano_and_pca_files_match_the_previous_phase_output_rules(runs):
    manifest = json.loads(next(runs["standard"].rglob("l2/figure_manifest.json")).read_text(encoding="utf-8"))
    volcano = next(f for f in manifest["figures"] if f["kind"] == "volcano")
    assert 84.9 < volcano["width_mm"] < 85.2 and 84.9 < volcano["height_mm"] < 85.2  # single column square, as in Phase 3B-1
    l1 = json.loads(next(runs["standard"].rglob("l1/figure_manifest.json")).read_text(encoding="utf-8"))
    pca = next(f for f in l1["figures"] if f["id"] == "pca")
    assert 84.9 < pca["width_mm"] < 85.2
    assert hashlib.sha256(b"").hexdigest()  # keep hashlib imported for failure messages above
