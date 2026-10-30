"""Standard-preset enrichment dotplots (GO, GSEA, KEGG ORA, KEGG GSEA) and the offline KEGG fixture."""

from __future__ import annotations

import hashlib
import importlib.util
import json
from pathlib import Path

import pytest

from conftest import require_r_packages, require_rscript
from test_figures import run_r

ROOT = Path(__file__).resolve().parents[1]
TOOLS = ROOT / "tools" / "figure_review"


def _load(name: str):
    spec = importlib.util.spec_from_file_location(name, TOOLS / f"{name}.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


compare = _load("compare_runs")

SETUP = '''
{ctx}
terms <- function(n, long = FALSE) {{
  desc <- if (long) sprintf("Regulation of a very long biological process description number %d involving many interacting components", seq_len(n)) else sprintf("Short term %d", seq_len(n))
  data.frame(ID = sprintf("T%04d", seq_len(n)), Description = desc, GeneRatio = sprintf("%d/200", 40 - seq_len(n)), BgRatio = sprintf("%d/5000", 100 + seq_len(n)),
    pvalue = 10^(-2 - seq_len(n)), p.adjust = 10^(-1.5 - seq_len(n)), qvalue = 10^(-1.5 - seq_len(n)), Count = as.integer(40 - seq_len(n)), stringsAsFactors = FALSE)
}}
gterms <- function(n, long = FALSE) {{
  t <- terms(n, long)
  data.frame(ID = t$ID, Description = t$Description, setSize = as.integer(30 + 7 * seq_len(n)), enrichmentScore = 0.5, NES = rev(seq(-2.4, 2.6, length.out = n)),
    pvalue = 10^(-3 - seq_len(n) / 4), p.adjust = 0.04 * seq_len(n) / n / 2, qvalue = 0.01, rank = as.integer(100 * seq_len(n)), leading_edge = "tags=40%", stringsAsFactors = FALSE)
}}
'''


def preamble(preset: str, profile: str = "general") -> str:
    ctx = f'ctx <- nf_rna_figure_context(list(figure_style=list(preset="{preset}", profile="{profile}")))'
    return SETUP.format(ctx=ctx) + "dir.create('out', showWarnings = FALSE)\n"


# ------------------------------------------------------------------ continuous enrichment scale


def test_the_enrichment_scale_is_its_own_sequential_scale_separate_from_group_and_expression_colours(tmp_path):
    result = run_r('''
std <- NF_RNA_STYLE_PRESETS$standard; leg <- NF_RNA_STYLE_PRESETS$legacy
ctx_std <- nf_rna_figure_scope(nf_rna_figure_context(list(figure_style = list(preset = "standard"))), "dotplot"); ctx_leg <- nf_rna_figure_context(list())
L <- function(h) { rgb <- t(grDevices::col2rgb(h)) / 255; 0.2126 * rgb[, 1] + 0.7152 * rgb[, 2] + 0.0722 * rgb[, 3] }
cat(toJSON(list(standard = nf_rna_enrichment_colours(ctx_std), legacy = nf_rna_enrichment_colours(ctx_leg), legacy_expected = rev(leg$sequential[-1]),
  categorical = unname(std$categorical), diverging = unname(std$diverging), batch = unname(std$batch), lightness = L(nf_rna_enrichment_colours(ctx_std)),
  mark_std = nf_rna_enrichment_mark(ctx_std), mark_legacy = nf_rna_enrichment_mark(ctx_leg), primary_legacy = nf_rna_primary_colour(ctx_leg)), auto_unbox = FALSE))''', tmp_path)
    assert result["legacy"] == result["legacy_expected"]  # legacy keeps its scale: low P dark
    assert result["mark_legacy"] == result["primary_legacy"]
    scale = result["standard"]
    assert len(scale) >= 4 and not set(scale) & set(result["categorical"] + result["diverging"] + result["batch"])  # nothing shared with the discrete or diverging palettes
    assert result["lightness"] == sorted(result["lightness"])  # most significant (first) is darkest, lightness rises monotonically
    assert result["mark_std"] == [scale[0]] or result["mark_std"] == scale[0]


# ------------------------------------------------------------------ unit tests of the two dotplot functions


def _tsv(path: Path) -> bytes:
    return path.read_bytes()


@pytest.mark.parametrize("kind", ["ora", "gsea"])
def test_standard_and_legacy_dotplots_share_source_data_term_order_and_mappings(tmp_path, kind):
    call = ('nf_rna_ora_dotplot(c, top, "-log10(adjusted P, BH)", file.path("out", stem), "out", list(), "133Gy vs 22Gy \\u2014 Significant DE genes")' if kind == "ora" else
            'nf_rna_gsea_dotplot(c, top, list(contrast_id = "x", numerator = "133Gy", denominator = "22Gy"), 0.05, "BH", file.path("out", stem), "out")')
    run_r(preamble("standard") + f'''
legacy <- nf_rna_figure_context(list()); std <- ctx
top <- {"terms(12, TRUE)" if kind == "ora" else "gterms(12, TRUE)"}
for (name in c("legacy", "standard")) {{ c <- if (name == "legacy") legacy else std; stem <- name; {call} }}
cat(toJSON(list(ok = TRUE), auto_unbox = TRUE))''', tmp_path)
    out = tmp_path / "out"
    assert _tsv(out / "legacy_source_data.tsv") == _tsv(out / "standard_source_data.tsv")  # the plotted values and their order are identical
    ids = [line.split("\t")[1] for line in (out / "standard_source_data.tsv").read_text().splitlines()[1:]]
    assert ids == [f"T{i:04d}" for i in range(1, 13)]
    for name in ("legacy", "standard"):
        assert all((out / f"{name}.{ext}").is_file() for ext in ("png", "tiff", "pdf"))  # every format kept
    assert compare.pdf_info(out / "standard.pdf")["embedded_fonts"] and all(f.startswith("Arimo") for f in compare.pdf_info(out / "standard.pdf")["embedded_fonts"])
    assert not all(f.startswith("Arimo") for f in compare.pdf_info(out / "legacy.pdf")["embedded_fonts"])  # legacy keeps its own font
    assert compare.edge_ink(out / "standard.png") == 0  # nothing runs past the figure edge


def test_gsea_colour_limits_size_mapping_and_direction_are_preserved(tmp_path):
    result = run_r(preamble("standard") + '''
legacy <- nf_rna_figure_context(list()); top <- gterms(8)
std_ctx <- nf_rna_figure_scope(ctx, "dotplot")
for (name in c("legacy", "standard")) { c <- if (name == "legacy") legacy else ctx
  rec <- nf_rna_gsea_dotplot(c, top, list(contrast_id = "x", numerator = "A", denominator = "B"), 0.05, "BH", file.path("out", name), "out") }
recs <- lapply(list(legacy, ctx), function(x) x$records[[1]])
cat(toJSON(list(colour = lapply(recs, function(r) r$colour), size = lapply(recs, function(r) r$size), x = lapply(recs, function(r) r$x), order = lapply(recs, function(r) r$display_order),
  standard_scale = recs[[2]]$enrichment_colours, terms = unlist(recs[[2]]$term_order)), auto_unbox = TRUE))''', tmp_path)
    assert result["colour"][0] == result["colour"][1]  # same quantity, same limits (0 .. cutoff), same kind
    assert result["colour"][0]["limits"] == [0, 0.05] and result["size"][0] == result["size"][1] and result["x"][0] == result["x"][1]
    assert result["order"][0] == result["order"][1] == "NES descending"
    assert result["terms"] == [f"T{i:04d}" for i in range(1, 9)]  # ordered by NES descending, as in legacy
    assert len(result["standard_scale"]) >= 4


@pytest.mark.parametrize("profile,max_w,max_h", [("general", 180, 220), ("nature", 183, 170), ("plos", 190.5, 222.3), ("elsevier", 190, 240)])
def test_dotplot_size_follows_the_measured_labels_and_term_count_inside_the_profile(tmp_path, profile, max_w, max_h):
    result = run_r(preamble("standard", profile) + '''
sizes <- list()
for (n in c(1, 3, 8, 15)) for (long in c(FALSE, TRUE)) {
  rec <- nf_rna_ora_dotplot(ctx, terms(n, long), "-log10(adjusted P, BH)", file.path("out", sprintf("o_%d_%s", n, long)), "out")
  g <- nf_rna_gsea_dotplot(ctx, gterms(n, long), list(contrast_id = "x", numerator = "A", denominator = "B"), 0.05, "BH", file.path("out", sprintf("g_%d_%s", n, long)), "out")
  sizes[[length(sizes) + 1]] <- list(n = n, long = long, ow = rec$width_mm, oh = rec$height_mm, gw = g$width_mm, gh = g$height_mm)
}
cat(toJSON(sizes, auto_unbox = TRUE))''', tmp_path)
    by = {(i["n"], i["long"]): i for i in result}
    for item in result:
        assert item["ow"] <= max_w + 1e-6 and item["gw"] <= max_w + 1e-6 and item["oh"] <= max_h + 1e-6 and item["gh"] <= max_h + 1e-6
    assert by[(1, False)]["oh"] < by[(8, False)]["oh"] < by[(15, False)]["oh"]  # more terms, more height
    assert by[(3, False)]["ow"] <= by[(3, True)]["ow"]  # labels are measured: wider labels never get a narrower figure
    assert by[(3, False)]["ow"] < max_w - 20  # compact: a three-term figure is not forced to the full profile width
    for stem in (tmp_path / "out").glob("*.png"):
        assert compare.edge_ink(stem) == 0, stem.name  # no clipping at any term count or label length


def test_empty_term_sets_are_a_caller_decision_and_the_dotplot_helpers_are_not_reached(tmp_path):
    """The backends only call the dotplot helpers when there is at least one significant term; the legacy rule is unchanged."""

    for script in ("go_analysis.R", "gsea_analysis.R", "kegg_analysis.R"):
        text = (ROOT / "src" / "rnaseq" / "r" / script).read_text(encoding="utf-8")
        assert "nrow(" in text and "dotplot(" in text
    figures_r = (ROOT / "src" / "rnaseq" / "r" / "figures.R").read_text(encoding="utf-8")
    assert 'NF_RNA_STANDARD_FIGURES <- c("volcano", "pca", "boxplot", "correlation", "heatmap", "dotplot")' in figures_r


# ------------------------------------------------------------------ the offline KEGG fixture


@pytest.fixture(scope="module")
def generate():
    # openssl: the KEGG backend verifies the fixture snapshot's SHA-256 manifest
    require_r_packages("DESeq2", "clusterProfiler", "org.Mm.eg.db", "ggplot2", "jsonlite", "systemfonts", "openssl")
    return _load("generate_runs")


def _kegg_run(generate, base: Path, scenario: str, figures: dict | None, name: str, monkeypatch) -> Path:
    monkeypatch.setenv("HTTPS_PROXY", "http://127.0.0.1:9")  # any real KEGG access would fail: the fixture must not need it
    monkeypatch.setenv("https_proxy", "http://127.0.0.1:9")
    out = base / name
    out.mkdir(parents=True)
    codes = generate.run_enrichment(out, figures, scenario, only_kegg=True, empty=scenario == "empty")
    assert codes == {"kegg": 0, "gsea_kegg": 0}, (name, codes, list((out / "enrichment").glob("*.stderr.log")))
    return out / "enrichment"


def test_the_offline_kegg_fixture_exercises_both_plotting_paths_without_network(generate, tmp_path, monkeypatch):
    run = _kegg_run(generate, tmp_path, "many", {"preset": "standard", "profile": "general"}, "standard", monkeypatch)
    for backend in ("kegg", "gsea_kegg"):
        assert (run / backend / "SYNTHETIC_KEGG_FIXTURE.txt").read_text().startswith("SYNTHETIC KEGG FIXTURE")  # outputs are marked as synthetic
    summary = json.loads((run / "kegg" / "kegg_backend_summary.json").read_text())
    assert summary["status"] == "SUCCESS"
    dotplots = sorted(path.relative_to(run).as_posix() for path in run.rglob("dotplot.png"))
    assert dotplots == ["gsea_kegg/Treatment_vs_Control/dotplot.png", "kegg/Treatment_vs_Control/down/dotplot.png", "kegg/Treatment_vs_Control/significant/dotplot.png", "kegg/Treatment_vs_Control/up/dotplot.png"]
    for backend in ("kegg", "gsea_kegg"):
        manifest = json.loads((run / backend / "figure_manifest.json").read_text())
        assert manifest["figure_style"]["standard_figures_drawn"] == ["dotplot"]
        for figure in manifest["figures"]:
            if figure["kind"] == "GSEA running-score plot":  # not a standard-preset figure type
                assert [item["format"] for item in figure["files"]] == ["png", "png_raster", "pdf"]
                continue
            assert figure["figure_style"]["font"] == "Arimo"
            assert [item["format"] for item in figure["files"]] == ["png", "png_raster", "tiff", "pdf"]
            assert figure["term_order"] and figure["layout"]["label_lines"] >= 1


def test_the_offline_fixture_is_deterministic_and_covers_empty_and_nonsignificant_results(generate, tmp_path, monkeypatch):
    first = _kegg_run(generate, tmp_path, "few", None, "a", monkeypatch)
    second = _kegg_run(generate, tmp_path, "few", None, "b", monkeypatch)
    tables = sorted(p.relative_to(first) for p in first.rglob("*.tsv"))
    assert tables and all((first / t).read_bytes() == (second / t).read_bytes() for t in tables)  # byte-identical on repeat
    significant = (first / "kegg" / "Treatment_vs_Control" / "significant" / "kegg.tsv").read_text().splitlines()
    assert len(significant) == 4  # header + the three "few" pathways
    for scenario in ("nonsig", "empty"):
        run = _kegg_run(generate, tmp_path, scenario, None, scenario, monkeypatch)
        assert not list(run.rglob("dotplot.png"))  # nothing significant: no figure is drawn, in either backend
        assert json.loads((run / "kegg" / "kegg_backend_summary.json").read_text())["status"] in {"NO_SIGNIFICANT_TERMS", "NOT_APPLICABLE"}
        assert json.loads((run / "gsea_kegg" / "gsea_kegg_backend_summary.json").read_text())["status"] in {"NO_SIGNIFICANT_TERMS", "NOT_APPLICABLE"}


def test_kegg_legacy_and_standard_share_source_data_and_term_order_and_standard_is_complete(generate, tmp_path, monkeypatch):
    legacy = _kegg_run(generate, tmp_path, "long", None, "legacy", monkeypatch)
    standard = _kegg_run(generate, tmp_path, "long", {"preset": "standard", "profile": "general"}, "standard", monkeypatch)
    report = compare.compare(legacy, standard, tmp_path / "report")
    # 4 dotplots plus the running-score plots, which stay identical across presets (not a standard figure type)
    assert report["summary"]["matched_figures"] == 4 + len(list((legacy / "gsea_kegg").rglob("running_score/*.pdf"))) and report["summary"]["unmatched"] == 0
    assert compare.gate_failures(report, ["ORA dotplot", "GSEA dotplot"]) == []
    for record in (r for r in report["figures"] if "running_score" not in r["figure"]):
        assert record["source_data"]["identical"] is True
        assert record["terms"]["order_identical"] is True and record["terms"]["drawn_matches_source"] is True
        assert record["clipping"]["candidate_edge_ink_px"] == 0
        assert record["font"]["candidate"]["embedded_in_pdf"] and all(f.startswith("Arimo") for f in record["font"]["candidate"]["embedded_in_pdf"])
        if record["legend"]["candidate"]["verified"]:
            assert record["legend"]["candidate"]["complete"] is True
        assert record["candidate"]["width_mm"] < record["legacy"]["width_mm"]  # content-aware, not the full double column
    digest = hashlib.sha256(b"".join((legacy / name).read_bytes() for name in sorted(p.relative_to(legacy).as_posix() for p in legacy.rglob("*.tsv")))).hexdigest()
    assert digest == hashlib.sha256(b"".join((standard / name).read_bytes() for name in sorted(p.relative_to(standard).as_posix() for p in standard.rglob("*.tsv")))).hexdigest()


def test_the_offline_fixture_files_are_self_contained_and_marked_synthetic():
    profile = (TOOLS / "kegg_offline_profile.R").read_text(encoding="utf-8")
    assert "NOT scientific results" in profile and "SYNTHETIC_KEGG_FIXTURE.txt" in profile
    listing = (TOOLS / "fixtures" / "kegg_list_pathway_mmu.txt").read_text(encoding="utf-8").splitlines()
    assert len(listing) >= 40 and all(line.startswith("path:mmu") and line.endswith(" - Mus musculus (house mouse)") for line in (l.split("\t")[0] + "\t" + l.split("\t")[1] for l in listing))
    assert any(len(line.split("\t")[1]) > 60 for line in listing)  # long KEGG names are represented
    assert require_rscript()
