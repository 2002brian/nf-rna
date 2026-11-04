"""v1.7.0 GUI responsiveness, workflow clarity and L1 figure consistency.

* references.sha256_file: the opt-in per-process memo used by the GUI (off for the CLI and runs).
* services.pairing_preview: the live paired-design check of the metadata table being edited.
* figure_export.png_option: PNG at the submission dpi from the pipeline file or the vector PDF.
* results.qc_view: physical sizes and display groups of the L1 figures.
* figures.R: the standard correlation heatmap's text follows nf_rna_theme (display only).
* pages: workflow bar, New Project explanations, pairing check, Results download menus.
"""

from __future__ import annotations

import json
import os
import shutil
from pathlib import Path

import pytest

from rnaseq import references
from rnaseq.gui import figure_export, results, services


# --------------------------------------------------------------------- checksum memo

@pytest.fixture
def memo(monkeypatch):
    monkeypatch.setattr(references, "CHECKSUM_MEMO_MIN_BYTES", 1024)
    references.disable_checksum_memo()
    yield
    references.disable_checksum_memo()


def _count_hashes(monkeypatch) -> list[int]:
    calls = []
    real = references.hashlib.sha256

    def counting():
        calls.append(1)
        return real()

    monkeypatch.setattr(references.hashlib, "sha256", counting)
    return calls


def test_checksum_memo_is_off_by_default_so_the_cli_always_hashes(tmp_path, memo, monkeypatch):
    big = tmp_path / "big.fa"
    big.write_bytes(b"A" * 4096)
    calls = _count_hashes(monkeypatch)
    first = references.sha256_file(big)
    second = references.sha256_file(big)
    assert first == second and len(calls) == 2


def test_checksum_memo_reuses_an_unchanged_large_file_and_rehashes_any_change(tmp_path, memo, monkeypatch):
    big = tmp_path / "big.fa"
    big.write_bytes(b"A" * 4096)
    references.enable_checksum_memo()
    calls = _count_hashes(monkeypatch)
    digest = references.sha256_file(big)
    assert references.sha256_file(big) == digest and len(calls) == 1
    # same size, different content, explicit later mtime: re-hashed and the new digest returned
    big.write_bytes(b"C" * 4096)
    stat = big.stat()
    os.utime(big, ns=(stat.st_atime_ns, stat.st_mtime_ns + 10_000_000))
    changed = references.sha256_file(big)
    assert changed != digest and len(calls) == 2
    # a replaced file (new inode) is re-hashed too
    replacement = tmp_path / "new.fa"
    replacement.write_bytes(b"A" * 4096)
    replacement.replace(big)
    assert references.sha256_file(big) == digest and len(calls) == 3


def test_checksum_memo_never_applies_to_small_project_files(tmp_path, memo, monkeypatch):
    small = tmp_path / "metadata.csv"
    small.write_text("sample_id,condition\nA,x\n")
    references.enable_checksum_memo()
    calls = _count_hashes(monkeypatch)
    references.sha256_file(small)
    references.sha256_file(small)
    assert len(calls) == 2


def test_disabling_the_memo_forgets_every_digest(tmp_path, memo, monkeypatch):
    big = tmp_path / "big.fa"
    big.write_bytes(b"A" * 4096)
    references.enable_checksum_memo()
    references.sha256_file(big)
    references.disable_checksum_memo()
    calls = _count_hashes(monkeypatch)
    references.sha256_file(big)
    assert len(calls) == 1


def test_the_gui_enables_the_memo_and_a_checksum_mismatch_is_still_reported(tmp_path, memo, monkeypatch):
    pytest.importorskip("streamlit")
    from rnaseq.gui import pages

    monkeypatch.setenv("NF_RNA_GUI_CONFIG", str(tmp_path / "gui.yaml"))
    (tmp_path / "gui.yaml").write_text(f"allowed_roots:\n- {tmp_path}\nproject_parent: {tmp_path}\n")
    pages.gui_config()
    assert references._CHECKSUM_MEMO is not None
    big = tmp_path / "genome.fa"
    big.write_bytes(b"A" * 4096)
    good = references.sha256_file(big)
    big.write_bytes(b"T" * 4096)
    stat = big.stat()
    os.utime(big, ns=(stat.st_atime_ns, stat.st_mtime_ns + 10_000_000))
    assert references.sha256_file(big) != good  # an asset modified in place no longer matches its manifest digest


# --------------------------------------------------------------------- pairing preview

COLUMNS = ["sample_id", "patient", "condition"]


def _rows(*items):
    return [dict(zip(COLUMNS, item)) for item in items]


def test_pairing_preview_counts_complete_pairs():
    preview = services.pairing_preview(COLUMNS, _rows(("N1", "p1", "Normal"), ("T1", "p1", "Tumor"), ("N2", "p2", "Normal"), ("T2", "p2", "Tumor")), "patient", "condition")
    assert preview.levels == ("Normal", "Tumor") and preview.complete == 2 and preview.problems == 0
    assert preview.pairs[0].samples == {"Normal": ("N1",), "Tumor": ("T1",)}


def test_pairing_preview_reports_missing_duplicate_and_blank_members():
    preview = services.pairing_preview(COLUMNS, _rows(
        ("N1", "p1", "Normal"), ("T1", "p1", "Tumor"),
        ("N2", "p2", "Normal"),                         # no Tumor
        ("N3", "p3", "Normal"), ("N3b", "p3", "Normal"), ("T3", "p3", "Tumor"),  # two Normal
        ("T4", "", "Tumor"),                            # no pair ID
    ), "patient", "condition")
    status = {row.pair: row.status for row in preview.pairs}
    assert status["p1"] == "complete" and status["p2"] == "no Tumor sample"
    assert status["p3"].startswith("2 Normal samples") and preview.blank_pair == ("T4",)
    assert preview.complete == 1 and preview.problems == 3


def test_pairing_preview_needs_two_levels_and_both_columns():
    preview = services.pairing_preview(COLUMNS, _rows(("A", "p1", "x"), ("B", "p1", "y"), ("C", "p1", "z")), "patient", "condition")
    assert preview.complete == 0 and "exactly 2 condition levels" in preview.pairs[0].status
    assert services.pairing_preview(["sample_id", "condition"], [], "patient", "condition") is None


# --------------------------------------------------------------------- PNG export

def _png(path: Path, dpi: int, size=(30, 20)) -> None:
    from PIL import Image

    Image.new("RGB", size, "white").save(path, dpi=(dpi, dpi))


def test_png_option_prefers_the_pipeline_raster_byte_for_byte(tmp_path):
    stem = tmp_path / "pca"
    _png(tmp_path / "pca_300dpi.png", 300)
    (tmp_path / "pca.pdf").write_bytes(b"%PDF-1.4\n")
    option = figure_export.png_option(stem)
    assert option.available and option.method == figure_export.PIPELINE_PNG.format(dpi=300) and option.file_name == "pca_300dpi.png"
    assert figure_export.png_bytes(option) == (tmp_path / "pca_300dpi.png").read_bytes()


def test_png_option_never_enlarges_the_preview_and_says_why_when_unavailable(tmp_path, monkeypatch):
    stem = tmp_path / "pca"
    _png(tmp_path / "pca.png", 150)
    option = figure_export.png_option(stem)
    assert not option.available and "no PDF" in option.reason
    (tmp_path / "pca.pdf").write_bytes(b"%PDF-1.4\n")
    monkeypatch.setattr(shutil, "which", lambda name: None)
    option = figure_export.png_option(stem)
    assert not option.available and "poppler" in option.reason


def test_png_option_renders_the_vector_pdf_at_the_exact_dpi(tmp_path):
    if not figure_export.poppler_available():
        pytest.skip("poppler (pdftocairo, pdfinfo) is not installed")
    from PIL import Image

    # A 144 x 72 pt (2 x 1 in) one-page PDF.
    Image.new("RGB", (144, 72), "white").save(tmp_path / "pca.pdf", resolution=72.0)
    option = figure_export.png_option(tmp_path / "pca")
    assert option.available and option.method == figure_export.PNG_PDF_RENDER.format(dpi=300)
    from io import BytesIO
    from PIL import Image

    with Image.open(BytesIO(figure_export.png_bytes(option))) as image:
        assert image.size == (600, 300) and image.mode == "RGB" and round(image.info["dpi"][0]) == 300


# --------------------------------------------------------------------- L1 figure sizes and groups

def test_qc_view_records_physical_sizes_and_display_groups(tmp_path):
    run_dir = tmp_path / "run"
    l1 = run_dir / "downstream" / "l1"
    l1.mkdir(parents=True)
    for stem in ("pca", "sample_correlation", "boxplot_log2_raw_counts", "boxplot_vst"):
        _png(l1 / f"{stem}.png", 150)
    (l1 / "figure_manifest.json").write_text(json.dumps({"figures": [
        {"id": "pca", "width_mm": 85.0, "height_mm": 97.7}, {"id": "sample_correlation", "width_mm": 121.7, "height_mm": 91.4},
        {"id": "boxplot_vst", "width_mm": 60.0, "height_mm": 73.4}, {"id": "broken", "width_mm": "x"},
    ]}))
    view = results.qc_view(results.ResultRun(run_dir, "CASE", "RUN", "SUCCESS", "L1"))
    by_label = {item.label: item for item in view.figures}
    assert by_label["PCA (blind VST)"].size_mm == (85.0, 97.7) and by_label["PCA (blind VST)"].group == results.L1_STRUCTURE
    assert by_label["Per-sample blind VST expression"].group == results.L1_DISTRIBUTION
    assert by_label["Per-sample log2(raw count + 1), before normalization"].size_mm is None  # not in the manifest: no invented size


# --------------------------------------------------------------------- figures.R: correlation heatmap text

def test_standard_correlation_heatmap_text_matches_the_ggplot_theme(tmp_path):
    from test_figures import run_r

    result = run_r('''
ctx <- nf_rna_figure_context(list(figure_style=list(preset="standard", profile="general")))
c <- nf_rna_figure_scope(ctx, "correlation")
set.seed(1); m <- matrix(rnorm(60), 10); corr <- cor(m); colnames(corr) <- rownames(corr) <- paste0("S", 1:6)
ann <- data.frame(condition=factor(rep(c("A", "B"), 3)), row.names=colnames(corr))
cols <- list(condition=nf_rna_group_colours(c("A", "B"), context=c))
plain <- pheatmap::pheatmap(corr, annotation_col=ann, annotation_colors=cols, main="Pearson correlation (blind VST)", fontsize=8, silent=TRUE)
fit <- nf_rna_correlation_heatmap(c, corr, nf_rna_sequential_palette(100, c), seq(min(corr), 1, length.out=101), ann, cols, "Pearson correlation (blind VST)")
g <- fit$heat$gtable
grob <- function(name) g$grobs[[which(g$layout$name == name)]]
texts <- function(x) if (inherits(x, "text")) list(x) else if (inherits(x, "gTree")) unlist(lapply(x$children, texts), recursive=FALSE) else list()
fonts <- function(name) vapply(texts(grob(name)), function(t) as.integer(t$gp[["font"]] %||% 1L), integer(1))
sizes <- function(name) vapply(texts(grob(name)), function(t) as.numeric(t$gp[["fontsize"]]), numeric(1))
cat(toJSON(list(main_size=grob("main")$gp[["fontsize"]], main_font=as.integer(grob("main")$gp[["font"]]),
  names_font=fonts("col_annotation_names"), legend_fonts=fonts("annotation_legend"), legend_sizes=sizes("annotation_legend"),
  key_sizes=sizes("legend"), order_same=identical(fit$heat$tree_col$order, plain$tree_col$order),
  base=c$profile$base_pt, small=c$profile$small_pt), auto_unbox=TRUE))
''', tmp_path)
    values = lambda key: set(result[key]) if isinstance(result[key], list) else {result[key]}  # noqa: E731 (auto_unbox makes 1-vectors scalars)
    assert result["main_size"] == result["base"] and result["main_font"] == 2  # theme plot.title: base_pt, bold (pheatmap: 1.3 x)
    assert values("names_font") == {1} and values("legend_fonts") == {1}  # plain, like ggplot legend titles
    assert values("legend_sizes") == {result["base"], result["small"]}  # title base_pt, labels small_pt
    assert values("key_sizes") == {result["small"]}
    assert result["order_same"]  # display only: the clustering order is pheatmap's


# --------------------------------------------------------------------- pages (AppTest)

def _app(name: str, project: Path | None = None):
    pytest.importorskip("streamlit")
    from streamlit.testing.v1 import AppTest

    page = AppTest.from_string(f"from rnaseq.gui import pages\npages.{name}()\n", default_timeout=60)
    if project is not None:
        page.session_state["project"] = str(project)
    return page


@pytest.fixture
def ux_home(tmp_path, monkeypatch):
    import yaml

    data = tmp_path / "data"
    (data / "projects").mkdir(parents=True)
    settings = tmp_path / "gui.yaml"
    settings.write_text(yaml.safe_dump({"allowed_roots": [str(data)], "project_parent": str(data / "projects"), "condarc": None, "launcher": "process"}), encoding="utf-8")
    monkeypatch.setenv("NF_RNA_GUI_CONFIG", str(settings))
    monkeypatch.setenv("XDG_CACHE_HOME", str(tmp_path / "cache"))
    yield data
    references.disable_checksum_memo()


def _paired_project(parent: Path) -> Path:
    from conftest import BASE_COUNTS
    from rnaseq.models import DesignType, InputType, Preset, Species
    from rnaseq.project import create_project, import_count_matrix

    project = create_project(project_name="paired", destination=parent, species=Species.HUMAN, preset=Preset.L2, design_type=DesignType.PAIRED_TWO_GROUP,
                             input_type=InputType.RAW_COUNTS, layout=None, scaffold=True, pair_id="patient", gene_id_type="ENSEMBL")
    counts = parent / "paired-counts.csv"
    counts.write_text(BASE_COUNTS, encoding="utf-8")
    import_count_matrix(project, counts)
    return project


def test_new_project_shows_the_formula_preprocessing_and_what_is_missing(ux_home):
    page = _app("new_project_page").run()
    assert not page.exception
    assert page.radio(key="np_input").value == "fastq"
    assert page.radio(key="np_pre").value == "raw"  # visible, default unchanged
    assert any("Still needed before" in item.value and "a project name" in item.value for item in page.caption)
    page.selectbox(key="np_design").set_value("paired_two_group").run()
    assert any("`~ patient + condition`" in item.value for item in page.markdown)
    assert any("Paired:" in item.value for item in page.caption)
    page.text_input(key="np_pair").set_value("donor").run()
    assert any("`~ donor + condition`" in item.value for item in page.markdown)


def test_inputs_page_shows_the_workflow_and_a_live_pairing_check(ux_home):
    project = _paired_project(ux_home / "projects")
    (project / "metadata.csv").write_text("sample_id,patient,condition\nC1,p1,A\nT1,p1,B\nC2,p2,A\nT2,p2,B\nC3,p3,A\nT3,,B\n", encoding="utf-8")
    page = _app("inputs_page", project).run()
    assert not page.exception, page.exception
    assert any("Dataset" in item.value and "→" in item.value and "Metadata" in item.value for item in page.markdown)
    warnings = [item.value for item in page.warning]
    assert any("Pairing check" in text and "2 complete pair(s) of 3" in text and "1 sample(s) without a `patient` value" in text for text in warnings)


def test_results_page_offers_every_figure_format_in_one_menu(ux_home, monkeypatch):
    from test_gui_results import _run, _write

    project = ux_home / "projects" / "demo"
    _write(project / "project.yaml", "project: {id: demo}\n")
    _run(project)
    page = _app("results_page", project).run()
    assert not page.exception
    keys = [button.key for button in page.button] + [button.key for button in page.get("download_button")]
    figure_keys = [key for key in keys if key and key.startswith("qc::")]
    assert {key.rsplit("::", 1)[1] for key in figure_keys} >= {"pdf", "png300", "tiff300", "tiff600", "png"}
    assert any("one common scale" in item.value for item in page.caption)


def test_run_summary_shows_fastq_preprocessing_and_explains_validation(ux_home, monkeypatch):
    import yaml
    from test_local_reference import _local_fastq_project

    project, _reference = _local_fastq_project(ux_home)
    config = yaml.safe_load((project / "project.yaml").read_text())
    config["input"]["preprocessing"] = "pretrimmed"
    (project / "project.yaml").write_text(yaml.safe_dump(config, sort_keys=False))
    monkeypatch.setattr("rnaseq.downstream_runtime.planned_runtime_prefix", lambda: (ux_home / "prefix", True))
    page = _app("runner_page", project).run()
    assert not page.exception, page.exception
    summary = next(frame.value for frame in page.dataframe if "item" in frame.value.columns and "strandedness" in set(frame.value["item"]))
    rows = dict(zip(summary["item"], summary["value"]))
    assert rows["FASTQ preprocessing"].startswith("pretrimmed")
    assert "`" not in rows.get("design", "")
    assert any("planning files" in item.value for item in page.caption)
    run = next(button for button in page.button if button.key == "run_start")
    assert run.disabled  # no case ID yet; the reason is stated
    assert any("is disabled: enter a case ID" in item.value for item in page.caption)
