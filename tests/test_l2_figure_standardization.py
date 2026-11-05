"""v1.7.0 L2 publication-figure refinement (standard preset).

Volcano/MA axis text names exactly the plotted statistics, the DEG heatmap keeps gene clustering
but orders samples by a documented deterministic rule, and every enrichment dotplot colours
adjusted P with the red (most significant) to blue scale.  Scientific tables are unchanged.
"""

from __future__ import annotations

import json
import subprocess
from pathlib import Path

import pytest

from conftest import require_r_packages
from test_figures import run_r

ROOT = Path(__file__).resolve().parents[1]
R_ROOT = ROOT / "src" / "rnaseq" / "r"

STD = 'std <- nf_rna_figure_context(list(figure_style = list(preset = "standard", profile = "general")))\nleg <- nf_rna_figure_context(list())\n'


# ------------------------------------------------------------------ 1-2. volcano / MA axis text


def test_volcano_and_ma_text_names_the_plotted_statistics(tmp_path):
    result = run_r(STD + '''
v <- nf_rna_figure_scope(std, "volcano")
cat(toJSON(list(
  shrunk = nf_rna_deg_plot_text(v, "Tumor", "Normal", TRUE, "apeglm", 0.05, 1),
  raw = nf_rna_deg_plot_text(v, "Tumor", "Normal", FALSE, NULL, 0.05, 1),
  legacy = nf_rna_deg_plot_text(leg, "Tumor", "Normal", TRUE, "apeglm", 0.05, 1)), auto_unbox = TRUE))''', tmp_path)
    shrunk, raw, legacy = result["shrunk"], result["raw"], result["legacy"]
    assert shrunk["volcano_x"] == "Shrunken log₂ fold change (Tumor vs Normal)" and shrunk["lfc_column"] == "log2FoldChange_shrunk"
    assert shrunk["ma_y"] == shrunk["volcano_x"]
    assert raw["volcano_x"] == "log₂ fold change (Tumor vs Normal)" and raw["lfc_column"] == "log2FoldChange"
    assert "Shrunken" not in raw["volcano_x"] and "shrinkage failed" in raw["caption"]
    for text in (shrunk, raw):
        # y is -log10(DESeq2 padj): labelled as the adjusted P value, never as a raw p-value
        assert text["p_column"] == "padj" and text["volcano_y"] == "−log₁₀(adjusted P value)"
        assert text["ma_x"] == "log₁₀(mean normalized count)" and text["mean_column"] == "baseMean"
        assert "apeglm" in shrunk["caption"] and "unshrunken" in text["caption"]
        for label in (text["volcano_x"], text["volcano_y"], text["ma_x"], text["ma_y"]):
            assert "padj" not in label and "baseMean" not in label and "/" not in label  # no column names or ratio slashes
    # legacy keeps its historical text byte for byte
    assert legacy["volcano_x"] == "shrunken log2 fold change (apeglm) (Tumor / Normal)"
    assert legacy["volcano_y"] == "-log10(adjusted P, BH)" and legacy["ma_x"] == "log10(baseMean)"
    assert legacy["caption"] == "DEG: padj < 0.05 and |raw log2FC| ≥ 1; LFC shrinkage: apeglm"


def test_l2_script_plots_padj_and_baseMean_with_the_text_helper():
    script = (R_ROOT / "l2_analysis.R").read_text(encoding="utf-8")
    # the volcano y value is derived from result_df$padj only (never pvalue)
    assert "plot_padj = ifelse(is.na(result_df$padj), NA_real_, pmax(result_df$padj, padj_floor))" in script
    assert "aes(x = lfc, y = -log10(plot_padj)" in script
    assert "ma_df <- data.frame(x = log10(result_df$baseMean), y = shown_lfc" in script
    assert "labs(x = deg_text$volcano_x, y = deg_text$volcano_y" in script and "labs(x = deg_text$ma_x, y = deg_text$ma_y" in script
    assert "pvalue)" not in script.split("plot_df <- data.frame(")[1].split("\n")[0]


# ------------------------------------------------------------------ 3-6. DEG heatmap sample order


def test_paired_sample_order_is_by_pair_then_condition_and_deterministic(tmp_path):
    result = run_r(STD + '''
samples <- c("N_1201020_P2", "N_17997_P2", "N_3031393_P2", "N_3038215_P2", "T_1201020_P1", "T_17997_P1", "T_3031393_P1", "T_3038215_P1")
condition <- rep(c("Normal", "Tumor"), each = 4); pair <- sub("^[NT]_([0-9]+)_P[12]$", "\\\\1", samples)
a <- nf_rna_heatmap_sample_order(samples, condition, c("Normal", "Tumor"), pair)
shuffled <- c(5, 2, 8, 1, 3, 7, 4, 6)
b <- nf_rna_heatmap_sample_order(samples[shuffled], condition[shuffled], c("Normal", "Tumor"), pair[shuffled])
u <- nf_rna_heatmap_sample_order(c("t2", "c1", "t1", "c2", "c3"), c("T", "C", "T", "C", "C"), c("C", "T"))
r <- nf_rna_heatmap_sample_order(c("p1_a", "p2_a", "p1_b", "p1_c"), c("T", "T", "N", "T"), c("N", "T"), c("p1", "p2", "p1", "p1"))
cat(toJSON(list(a = a, again = nf_rna_heatmap_sample_order(samples, condition, c("Normal", "Tumor"), pair), b = b, u = u, r = r)))''', tmp_path)
    assert result["a"] == ["N_1201020_P2", "T_1201020_P1", "N_17997_P2", "T_17997_P1", "N_3031393_P2", "T_3031393_P1", "N_3038215_P2", "T_3038215_P1"]
    assert result["again"] == result["a"]
    # pair order follows first appearance in the given sample sheet; within a pair the reference level comes first
    # shuffled sheet: T_1201020, N_17997, T_3038215, N_1201020, N_3031393, ... -> pairs 1201020, 17997, 3038215, 3031393
    assert result["b"] == ["N_1201020_P2", "T_1201020_P1", "N_17997_P2", "T_17997_P1", "N_3038215_P2", "T_3038215_P1", "N_3031393_P2", "T_3031393_P1"]
    assert result["u"] == ["c1", "c2", "c3", "t2", "t1"]  # unpaired: reference level, then sample-sheet order
    assert result["r"] == ["p1_b", "p1_a", "p1_c", "p2_a"]  # replicates within a pair keep sample-sheet order


PAIRED_SAMPLES = ["T_p3", "N_p1", "T_p1", "N_p3", "N_p2", "T_p2", "N_p4", "T_p4"]  # deliberately interleaved sample sheet


def _write_l2_fixture(root: Path, preset: str) -> Path:
    import numpy as np

    rng = np.random.default_rng(7)
    genes = [f"g{i:03d}" for i in range(1, 301)]
    patient_effect = {"p1": 0.0, "p2": 0.4, "p3": -0.3, "p4": 0.2}
    rows, vst = [], []
    for index, gene in enumerate(genes):
        base = rng.uniform(4, 9)
        lfc = 2.0 if index < 25 else (-2.0 if index < 50 else 0.0)
        values = []
        for sample in PAIRED_SAMPLES:
            condition, patient = sample.split("_")
            mu = np.exp(base + patient_effect[patient] + (lfc * np.log(2) if condition == "T" else 0.0))
            values.append(int(rng.negative_binomial(20, 20 / (20 + mu))))
        rows.append([gene, *values])
        vst.append([gene, *[f"{np.log2(v + 1):.6f}" for v in values]])
    work = root / preset
    work.mkdir(parents=True)
    (work / "counts.csv").write_text("gene_id," + ",".join(PAIRED_SAMPLES) + "\n" + "\n".join(",".join(map(str, r)) for r in rows) + "\n")
    (work / "vst.tsv").write_text("gene_id\t" + "\t".join(PAIRED_SAMPLES) + "\n" + "\n".join("\t".join(r) for r in vst) + "\n")
    (work / "metadata.csv").write_text("sample_id,patient,condition\n" + "\n".join(f"{s},{s.split('_')[1]},{'Tumor' if s[0] == 'T' else 'Normal'}" for s in PAIRED_SAMPLES) + "\n")
    config = {
        "contrasts": [{"contrast_id": "Tumor_vs_Normal", "factor": "condition", "numerator": "Tumor", "denominator": "Normal"}],
        "counts": str(work / "counts.csv"), "metadata": str(work / "metadata.csv"), "l1_vst": str(work / "vst.tsv"),
        "design_variable_types": {"condition": "categorical", "patient": "categorical"}, "formula": "~ patient + condition", "pair_id": "patient",
        "figure_style": {"preset": preset, "profile": "general", "source": "explicit"},
        "filter": {"minimum_total_count": 10}, "heatmap_top_n": 30, "output_dir": str(work / "out"),
        "samples": PAIRED_SAMPLES, "source_type": "raw_counts", "thresholds": {"abs_log2fc": 1.0, "padj": 0.05},
    }
    path = work / "config.json"
    path.write_text(json.dumps(config))
    return path


@pytest.fixture(scope="module")
def paired_l2_runs(tmp_path_factory):
    executable = require_r_packages("DESeq2", "apeglm", "ashr", "ggplot2", "ggrepel", "pheatmap", "jsonlite")
    pytest.importorskip("numpy")
    root = tmp_path_factory.mktemp("l2fig")
    outputs = {}
    for preset in ("standard", "legacy"):
        config = _write_l2_fixture(root, preset)
        run = subprocess.run([executable, str(R_ROOT / "l2_analysis.R"), "--config", str(config)], cwd=config.parent, capture_output=True, text=True, check=False)
        assert run.returncode == 0, run.stderr[-3000:]
        outputs[preset] = config.parent / "out"
    return outputs


def _records(out: Path) -> dict:
    manifest = json.loads((out / "figure_manifest.json").read_text())
    return {record["kind"]: record for record in manifest["figures"]}


def test_standard_deg_heatmap_orders_samples_without_clustering_and_keeps_gene_clustering(paired_l2_runs):
    standard, legacy = _records(paired_l2_runs["standard"]), _records(paired_l2_runs["legacy"])
    heat = standard["DEG heatmap"]
    assert heat["display_order"]["columns"] == ["N_p3", "T_p3", "N_p1", "T_p1", "N_p2", "T_p2", "N_p4", "T_p4"]
    assert heat["clustering"]["columns"].startswith("off (not clustered); display order: pairs in order of first appearance")
    assert heat["clustering"]["rows"] == "complete linkage, Euclidean distance on unclipped z-scores"
    assert legacy["DEG heatmap"]["clustering"]["columns"] == "complete linkage, Euclidean distance on unclipped z-scores"  # legacy unchanged
    # gene (row) clustering is the same hclust in both presets
    assert heat["display_order"]["rows"] and len(heat["display_order"]["rows"]) == heat["display_subset"]["drawn"]


def test_heatmap_values_and_all_scientific_tables_are_byte_identical_across_presets(paired_l2_runs):
    contrast = Path("contrasts") / "Tumor_vs_Normal"
    names = ["all_genes.tsv", "significant.tsv", "upregulated.tsv", "downregulated.tsv", "heatmap_genes.tsv", "heatmap_zscores.tsv", "backend_summary.json", "contrast_metadata.json"]
    for name in names:
        assert (paired_l2_runs["standard"] / contrast / name).read_bytes() == (paired_l2_runs["legacy"] / contrast / name).read_bytes(), name
    assert (paired_l2_runs["standard"] / "normalized_counts.tsv").read_bytes() == (paired_l2_runs["legacy"] / "normalized_counts.tsv").read_bytes()
    header = (paired_l2_runs["standard"] / contrast / "heatmap_zscores.tsv").read_text().splitlines()[0].split("\t")
    assert header == ["gene_id", *PAIRED_SAMPLES]  # the z-score table keeps the sample-sheet columns; only the drawing is reordered


def test_volcano_and_ma_manifest_axes_match_the_plotted_columns(paired_l2_runs):
    records = _records(paired_l2_runs["standard"])
    volcano, ma = records["volcano"], records["MA plot"]
    assert volcano["axes"]["y_column"] == "padj" and volcano["axes"]["y"] == "−log₁₀(adjusted P value)"
    shrunk = "log2FoldChange_shrunk" in volcano["statistics"]["fold_change"]
    assert volcano["axes"]["x_column"] == ("log2FoldChange_shrunk" if shrunk else "log2FoldChange")
    assert volcano["axes"]["x"].startswith("Shrunken log₂ fold change" if shrunk else "log₂ fold change")
    assert ma["axes"]["x_column"] == "baseMean" and ma["axes"]["x"] == "log₁₀(mean normalized count)" and ma["axes"]["y_column"] == volcano["axes"]["x_column"]
    for record in (volcano, ma, records["DEG heatmap"]):
        assert any(f["format"] == "png_raster" and f["dpi"] == 300 for f in record["files"]) and any(f["format"] == "pdf" for f in record["files"])


# ------------------------------------------------------------------ 7-8, 10. dotplot colour scale and meaning


def test_enrichment_scale_runs_from_red_most_significant_to_blue_and_signed_scale_is_separate(tmp_path):
    result = run_r(STD + '''
d <- nf_rna_figure_scope(std, "dotplot")
rgb <- function(x) unname(t(grDevices::col2rgb(x))[, 1:3, drop = FALSE])
cat(toJSON(list(scale = nf_rna_enrichment_colours(d), rgb = rgb(nf_rna_enrichment_colours(d)), signed = unname(nf_rna_signed_colours(d)),
  signed_rgb = rgb(nf_rna_signed_colours(d)), name = NF_RNA_STYLE_PRESETS$standard$enrichment_name, legacy = nf_rna_enrichment_colours(leg)), auto_unbox = TRUE))''', tmp_path)
    rgb = result["rgb"]
    assert rgb[0][0] > rgb[0][2] and rgb[-1][2] > rgb[-1][0]  # first (adjusted P = 0) red, last (cutoff) blue
    redness = [r - b for r, _, b in rgb]
    assert redness == sorted(redness, reverse=True)  # monotonic red -> blue
    assert result["name"] == "nf-rna Standard Enrichment Scale v2"
    low, mid, high = result["signed_rgb"]
    assert low[2] > low[0] and high[0] > high[2] and min(mid) > 230  # blue - white - red, white centre
    assert result["signed"][1] not in result["scale"]  # p.adjust never shares the signed (white-centred) palette
    assert result["legacy"] == ["#08306B", "#2171B5", "#6BAED6", "#C6DBEF"]  # legacy unchanged


DOTPLOT_SETUP = STD + '''
ora <- data.frame(ID = sprintf("GO:%07d", 1:6), Description = sprintf("Process %d", 1:6), GeneRatio = "5/100", BgRatio = "50/5000",
  pvalue = c(1e-9, 1e-6, 1e-4, 1e-3, 5e-3, 1e-2), p.adjust = c(1e-8, 1e-5, 1e-3, 0.01, 0.03, 0.049), qvalue = 0.01, Count = 6:1, stringsAsFactors = FALSE)
gsea <- data.frame(ID = sprintf("hsa%05d", 1:6), Description = sprintf("Pathway %d", 1:6), setSize = 20L + 1:6, enrichmentScore = 0.5, NES = c(2.1, 1.8, 1.2, -1.3, -1.9, -2.4),
  pvalue = 1e-4, p.adjust = c(1e-6, 0.002, 0.04, 0.03, 0.01, 1e-4), qvalue = 0.01, rank = 100L, leading_edge = "tags=40%", stringsAsFactors = FALSE)
ora_before <- ora; gsea_before <- gsea
dir.create("out", showWarnings = FALSE)
inspect <- function(p) { b <- ggplot2::ggplot_build(p); layer <- which(vapply(p$layers, function(l) inherits(l$geom, "GeomPoint"), logical(1)))[[1]]
  sc <- b$plot$scales$get_scales("colour")
  list(colour_var = rlang::as_label(p$mapping$colour), x_var = rlang::as_label(p$mapping$x), size_var = rlang::as_label(p$mapping$size),
       title = sc$name, limits = sc$get_limits(), colours = b$data[[layer]]$colour) }
'''


def test_ora_and_gsea_dotplots_map_adjusted_p_to_the_red_blue_scale(tmp_path):
    result = run_r(DOTPLOT_SETUP + '''
ora_rec <- nf_rna_ora_dotplot(std, ora, nf_rna_ora_x_label(std, "GO BP", "BH"), file.path("out", "ora"), "out", padj_cutoff = 0.05, p_adjust_method = "BH",
  subtitle = nf_rna_ora_subtitle(std, list(numerator = "Tumor", denominator = "Normal"), "significant", database = "GO BP"))
ora_plot <- inspect(ggplot2::last_plot())
gsea_rec <- nf_rna_gsea_dotplot(std, gsea, list(contrast_id = "x", numerator = "Tumor", denominator = "Normal"), 0.05, "BH", file.path("out", "gsea"), "out")
gsea_plot <- inspect(ggplot2::last_plot()); gp <- ggplot2::last_plot()
cat(toJSON(list(ora = ora_plot, gsea = gsea_plot, ora_rec = ora_rec$colour, gsea_rec = gsea_rec$colour_meaning, ora_x = ora_rec$x, gsea_x = gp$labels$x,
  ora_subtitle = ora_rec$subtitle, ora_xlabel = nf_rna_ora_x_label(std, "GO BP", "BH"), legacy_xlabel = nf_rna_ora_x_label(leg, "GO BP", "BH"),
  unchanged = identical(ora, ora_before) && identical(gsea, gsea_before), scale = nf_rna_enrichment_colours(nf_rna_figure_scope(std, "dotplot"))), auto_unbox = TRUE))''', tmp_path)
    for kind in ("ora", "gsea"):
        plot = result[kind]
        assert plot["colour_var"] == "p.adjust" and plot["title"] == "Adjusted P value (BH)" and plot["limits"] == [0, 0.05]
    assert result["ora"]["x_var"] == "-log10(plot_padj)" and result["ora"]["size_var"] == "Count"  # axes and size unchanged
    assert result["gsea"]["x_var"] == "NES" and result["gsea"]["size_var"] == "setSize"  # NES stays on x, never on colour
    # the most significant point is drawn in the first (red) stop, the least significant is the bluest
    ora_colours = result["ora"]["colours"]
    assert ora_colours[0].upper().startswith(result["scale"][0].upper())
    def red_minus_blue(h): return int(h[1:3], 16) - int(h[5:7], 16)
    assert red_minus_blue(ora_colours[0]) > 0 > red_minus_blue(ora_colours[-1])
    assert [red_minus_blue(c) for c in ora_colours] == sorted((red_minus_blue(c) for c in ora_colours), reverse=True)
    assert result["ora_rec"]["column"] == "p.adjust" and result["ora_rec"]["scale"] == "nf-rna Standard Enrichment Scale v2"
    assert "NES is drawn on the x axis" in result["gsea_rec"]["signed_statistic"]
    assert result["gsea_x"] == "Normalized enrichment score (Tumor vs Normal)"
    assert result["ora_xlabel"] == "−log₁₀(adjusted P value)" and result["legacy_xlabel"] == "-log10(GO BP adjusted P, BH)"
    assert result["ora_subtitle"].startswith("GO BP: Tumor vs Normal")
    assert result["unchanged"] is True  # the enrichment tables handed to the plots are not modified


def test_dotplot_source_tables_are_byte_identical_across_presets(tmp_path):
    run_r(DOTPLOT_SETUP + '''
for (name in c("legacy", "standard")) { c <- if (name == "legacy") leg else std
  nf_rna_ora_dotplot(c, ora, nf_rna_ora_x_label(c, "KEGG", "BH"), file.path("out", paste0(name, "_ora")), "out")
  nf_rna_gsea_dotplot(c, gsea, list(contrast_id = "x", numerator = "Tumor", denominator = "Normal"), 0.05, "BH", file.path("out", paste0(name, "_gsea")), "out") }
cat(toJSON(list(ok = TRUE), auto_unbox = TRUE))''', tmp_path)
    out = tmp_path / "out"
    for kind in ("ora", "gsea"):
        assert (out / f"legacy_{kind}_source_data.tsv").read_bytes() == (out / f"standard_{kind}_source_data.tsv").read_bytes()
