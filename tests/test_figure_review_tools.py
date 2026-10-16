"""Figure-review tooling (tools/figure_review): comparison workflow and palette review."""

from __future__ import annotations

import csv
import importlib.util
import json
import subprocess
from pathlib import Path

import pytest

from conftest import require_rscript

TOOLS = Path(__file__).resolve().parents[1] / "tools" / "figure_review"
spec = importlib.util.spec_from_file_location("compare_runs", TOOLS / "compare_runs.py")
compare_runs = importlib.util.module_from_spec(spec)
spec.loader.exec_module(compare_runs)


def png(path: Path, pixels: dict[tuple[int, int], tuple[int, int, int]] | None = None, size=(6, 4)) -> None:
    pixels = pixels or {}
    rows = [b"".join(bytes(pixels.get((x, y), (255, 255, 255))) for x in range(size[0])) for y in range(size[1])]
    compare_runs.write_png(path, size[0], size[1], rows)


def make_run(root: Path, *, pixel=None, width_mm=179.92, font="DejaVu Sans", source="a\n1\n", style=None, extra_figure=False) -> Path:
    folder = root / "downstream" / "l1"
    folder.mkdir(parents=True)
    (folder / "pca_source_data.tsv").write_text(source, encoding="utf-8")
    png(folder / "pca.png", pixel)
    figures = [{"id": "pca", "width_mm": width_mm, "height_mm": 50.0, "kind": "PCA",
                "files": [{"path": "pca.png", "format": "png", "dpi": 150, "width_px": 6, "height_px": 4}],
                "source_data": [{"path": "pca_source_data.tsv", "md5": "x", "description": "d"}]}]
    if extra_figure:
        png(folder / "extra.png")
        figures.append({"id": "extra", "width_mm": 10, "height_mm": 10, "files": [{"path": "extra.png", "format": "png", "dpi": 150}], "source_data": []})
    manifest = {"schema_version": "nf-rna.figure-manifest.v1", "profile": {"name": "general"}, "r_version": "R test",
                "font": {"requested": "Arial", "family": font, "fallback_used": True, "file": "/f.ttf"}, "figures": figures}
    if style:
        manifest["figure_style"] = style
    (folder / "figure_manifest.json").write_text(json.dumps(manifest), encoding="utf-8")
    return root


def run(tmp_path, **candidate):
    legacy = make_run(tmp_path / "legacy")
    cand = make_run(tmp_path / "candidate", **candidate)
    return compare_runs.compare(legacy, cand, tmp_path / "out")


def test_identical_runs_are_reported_identical(tmp_path):
    (tmp_path / "out").mkdir()
    report = run(tmp_path)
    assert report["summary"] == {"matched_figures": 1, "by_status": {"identical": 1}, "unmatched": 0}
    assert report["figures"][0]["diff_png"] is None and report["figures"][0]["pixels"]["identical"]


def test_one_changed_pixel_is_found_and_highlighted(tmp_path):
    (tmp_path / "out").mkdir()
    report = run(tmp_path, pixel={(2, 1): (0, 0, 0)})
    figure = report["figures"][0]
    assert figure["status"] == "appearance differs"
    assert figure["pixels"]["differing_pixels"] == 1 and figure["pixels"]["max_channel_delta"] == 255
    diff = compare_runs.read_png(tmp_path / "out" / "assets" / figure["diff_png"])
    assert diff[3][1][2 * 3:2 * 3 + 3] == b"\xe0\x00\x00"


def test_a_size_change_is_a_layout_difference_and_unmatched_figures_are_listed(tmp_path):
    (tmp_path / "out").mkdir()
    report = run(tmp_path, width_mm=183.1, extra_figure=True)
    assert report["figures"][0]["status"] == "layout differs"
    assert report["only_in_candidate"] == ["downstream/l1/extra"] and report["summary"]["unmatched"] == 1


def test_changed_source_data_is_never_called_identical(tmp_path):
    (tmp_path / "out").mkdir()
    report = run(tmp_path, source="a\n2\n")
    assert report["figures"][0]["source_data"]["identical"] is False and report["figures"][0]["status"] != "identical"


def test_font_and_style_are_reported_from_each_manifest(tmp_path):
    (tmp_path / "out").mkdir()
    report = run(tmp_path, font="Arimo", style={"requested_preset": "bioinfo", "preset": "standard", "profile": "general", "standard_figures_drawn": ["pca"],
                                                 "font": {"family": "Arimo", "sha256": "ab" * 32, "verified": True}})
    assert report["legacy_context"]["font_family"] == "DejaVu Sans" and report["candidate_context"]["font_family"] == "Arimo"
    assert report["candidate_context"]["requested_preset"] == "bioinfo" and report["candidate_context"]["preset"] == "standard"
    assert report["candidate_context"]["standard_font_family"] == "Arimo"
    assert "Arimo" in compare_runs.render_html(report) and "standard" in compare_runs.render_html(report)


def test_every_figure_gets_a_side_by_side_png_and_a_font_verdict(tmp_path):
    (tmp_path / "out").mkdir()
    report = run(tmp_path, pixel={(1, 1): (0, 0, 0)})
    record = report["figures"][0]
    width, height, _channels, _rows = compare_runs.read_png(tmp_path / "out" / "assets" / record["side_by_side_png"])
    assert (width, height) == (6 + 24 + 6, 4)  # both previews at their own pixel size
    assert set(record["font"]) == {"legacy", "candidate"} and record["font"]["candidate"]["pdf_matches_declared"] is False  # no PDF in the fixture


def test_the_scoped_gate_allows_only_the_listed_kinds_to_change(tmp_path):
    (tmp_path / "out").mkdir()
    changed = run(tmp_path, pixel={(1, 1): (0, 0, 0)})
    assert compare_runs.gate_failures(changed, ["volcano"])[0].endswith("unexpected change (appearance differs)")
    assert not [f for f in compare_runs.gate_failures(changed, ["PCA"]) if "unexpected change" in f]  # PCA may change ...
    assert any("embedded PDF font" in f for f in compare_runs.gate_failures(changed, ["pca"]))  # ... but must show a matching font
    changed["figures"][0]["source_data"]["identical"] = False
    assert any("source data differs" in f for f in compare_runs.gate_failures(changed, ["pca"]))


def test_cli_writes_a_deterministic_report_and_can_require_identity(tmp_path):
    make_run(tmp_path / "a")
    make_run(tmp_path / "b")
    args = ["--legacy", str(tmp_path / "a"), "--candidate", str(tmp_path / "b")]
    assert compare_runs.main([*args, "--out", str(tmp_path / "r1"), "--require-identical"]) == 0
    assert compare_runs.main([*args, "--out", str(tmp_path / "r2")]) == 0
    assert (tmp_path / "r1" / "report.json").read_bytes() == (tmp_path / "r2" / "report.json").read_bytes()
    assert (tmp_path / "r1" / "report.html").read_bytes() == (tmp_path / "r2" / "report.html").read_bytes()
    make_run(tmp_path / "c", width_mm=1.0)
    assert compare_runs.main(["--legacy", str(tmp_path / "a"), "--candidate", str(tmp_path / "c"), "--out", str(tmp_path / "r3"), "--require-identical"]) == 1


def test_png_reader_handles_every_cairo_filter(tmp_path):
    """Round-trip a real cairo PNG (which uses the adaptive filters) through the decoder."""

    executable = require_rscript()
    target = tmp_path / "cairo.png"
    code = f'png("{target}", width=200, height=150, type="cairo", bg="white"); plot(1:5, pch=19, col="red"); text(3, 3, "a"); invisible(dev.off())'
    result = subprocess.run([executable, "-e", code], capture_output=True, text=True, check=False)
    assert result.returncode == 0, result.stderr
    width, height, channels, rows = compare_runs.read_png(target)
    assert (width, height) == (200, 150) and channels in (3, 4) and len(rows) == 150 and all(len(row) == 200 * channels for row in rows)
    assert any(row[i:i + 3] != b"\xff\xff\xff" for row in rows for i in range(0, len(row), channels))  # something was drawn


def test_palette_review_keeps_the_approved_colours_and_writes_the_review_files(tmp_path):
    executable = require_rscript()
    result = subprocess.run([executable, str(TOOLS / "palette_review.R"), str(tmp_path)], capture_output=True, text=True, check=False)
    if result.returncode != 0 and "colorspace" in result.stderr:
        pytest.skip("colorspace unavailable")
    assert result.returncode == 0, result.stderr
    for name in ("palette_groups", "palette_batch", "palette_batch_demo"):
        assert (tmp_path / f"{name}.png").stat().st_size > 1000 and (tmp_path / f"{name}.pdf").stat().st_size > 1000
    data = json.loads((tmp_path / "palette_metrics.json").read_text(encoding="utf-8"))
    assert data["approved_groups_1_3"] == ["#4E78A7", "#CC5A6A", "#59A14F"]
    assert data["status"].startswith("CANDIDATES ONLY")
    assert len(data["batch_candidates"]) >= 2 and all(len(v) >= 8 for v in data["batch_candidates"].values())
    assert all(len(v) == 5 for v in data["group_candidates_4_8"].values())
    rows = list(csv.DictReader((tmp_path / "palette_metrics.tsv").open(encoding="utf-8"), delimiter="\t"))
    assert {r["vision"] for r in rows} == {"normal", "deuteranopia", "protanopia", "tritanopia"}
    assert all(float(r["min_de2000_within"]) > 0 for r in rows)
    # The review never edits the preset table: candidates stay out of figures.R.
    figures_r = (TOOLS.parents[1] / "src" / "rnaseq" / "r" / "figures.R").read_text(encoding="utf-8")
    assert "NF_RNA_BIOINFO_ACTIVE" not in figures_r and "categorical=c(" in figures_r


# ------------------------------------------------------------ self-containment


ROOT = Path(__file__).resolve().parents[1]
EXTERNAL = ("bioinfo_toolkit", "bioinfo-plotting", "/home/brian", "04_plots", "PLOTTING_GUIDE")


def test_nf_rna_does_not_reference_an_external_plotting_toolkit_or_local_paths():
    """The figure system is owned by nf-rna: nothing shipped, built or documented points outside it."""

    files = [*(ROOT / "src").rglob("*"), *(ROOT / "tools").rglob("*"), *(ROOT / "workflow").rglob("*"), ROOT / "pyproject.toml", ROOT / "docs" / "figure-standard.md", ROOT / "environment.yml", *(ROOT / ".claude").rglob("*")]
    offenders = []
    for path in files:
        if path.is_file() and "__pycache__" not in path.parts and path.suffix not in {".pyc", ".gz", ".png", ".pdf"}:
            text = path.read_text(encoding="utf-8", errors="ignore")
            offenders += [f"{path.relative_to(ROOT)}: {needle}" for needle in EXTERNAL if needle in text]
    assert offenders == []


def test_r_scripts_source_only_files_shipped_next_to_themselves():
    import re

    for script in (ROOT / "src" / "rnaseq" / "r").glob("*.R"):
        for line in script.read_text(encoding="utf-8").splitlines():
            if re.match(r"\s*source\(", line):
                assert "dirname(normalizePath(script_file))" in line, f"{script.name}: {line.strip()}"


def test_figure_standard_document_exists_and_matches_the_code_flags():
    document = (ROOT / "docs" / "figure-standard.md").read_text(encoding="utf-8")
    figures_r = (ROOT / "src" / "rnaseq" / "r" / "figures.R").read_text(encoding="utf-8")
    assert 'NF_RNA_STANDARD_FIGURES <- c("volcano", "pca", "boxplot", "correlation", "heatmap", "dotplot")' in figures_r and "volcano" in document and "PCA" in document
    assert "NF_RNA_BIOINFO_ACTIVE" not in figures_r and "Declared but not active" not in document
    for colour in ("#4E78A7", "#CC5A6A", "#59A14F", "#B2182B", "#2166AC"):
        assert colour in document and colour in figures_r
    assert "docs/figure-standard.md" in figures_r


def test_the_project_skill_points_only_at_nf_rna_files_and_does_not_duplicate_constants():
    import re

    skill = ROOT / ".claude" / "skills" / "nf-rna-figures" / "SKILL.md"
    text = skill.read_text(encoding="utf-8")
    assert text.startswith("---\nname: nf-rna-figures\n")
    paths = set(re.findall(r"`((?:docs|src|tools|tests)/[^`\s]+)`", text))
    assert paths, "the skill must name the repository files it relies on"
    for path in paths:
        assert any(ROOT.glob(path.split(" ")[0])) or (ROOT / path).exists(), f"skill references a missing file: {path}"
    assert not re.search(r"#[0-9A-Fa-f]{6}\b", text)  # no palette constants
    assert not re.search(r"\b\d+(\.\d+)?\s*(mm|dpi)\b", text)  # no layout constants


# ------------------------------------------------------------ Phase 4 audit tools


def _load(name: str):
    import sys

    sys.path.insert(0, str(TOOLS))
    module_spec = importlib.util.spec_from_file_location(name, TOOLS / f"{name}.py")
    module = importlib.util.module_from_spec(module_spec)
    module_spec.loader.exec_module(module)
    return module


def test_invariance_audit_requires_identical_tables_and_ignores_figure_style(tmp_path):
    audit = _load("invariance_audit")
    for name, preset in (("a", "legacy"), ("b", "standard")):
        run = tmp_path / name
        (run / "l1").mkdir(parents=True)
        (run / "l1" / "normalized_counts.tsv").write_text("g\t1\n", encoding="utf-8")
        (run / "l1" / "pca.png").write_bytes(name.encode())
        (run / "backend_config.json").write_text(json.dumps({"out": str(run), "figure_style": {"preset": preset}}), encoding="utf-8")
    out = tmp_path / "r.json"
    assert audit.main(["--legacy", str(tmp_path / "a"), "--candidate", str(tmp_path / "b"), "--out", str(out)]) == 0
    report = json.loads(out.read_text(encoding="utf-8"))
    assert report["excluded_figure_artifacts"] == 1 and not report["metadata_differences_beyond_figure_style"]
    (tmp_path / "b" / "l1" / "normalized_counts.tsv").write_text("g\t2\n", encoding="utf-8")
    assert audit.main(["--legacy", str(tmp_path / "a"), "--candidate", str(tmp_path / "b")]) == 1


def test_pdf_text_sizes_are_read_from_the_content_stream(tmp_path):
    import zlib

    audit = _load("profile_audit")
    content = zlib.compress(b"BT\n6 0 0 -6 10 10 Tm\n/f-0-0 1 Tf\n(a)Tj\nET BT\n0 -8 -8 0 5 5 Tm\n/f-0-0 1 Tf\n(b)Tj\nET")
    pdf = tmp_path / "x.pdf"
    pdf.write_bytes(b"%PDF-1.5\n1 0 obj\n<<>>\nstream\n" + content + b"\nendstream\nendobj\n")
    assert audit.text_sizes(pdf) == {"runs": 2, "implied_min_pt": 6.0, "implied_max_pt": 8.0}


def test_manifest_review_reports_added_and_removed_keys(tmp_path):
    review = _load("manifest_review")
    for name, record in (("a", {"id": "pca", "width_mm": 1}), ("b", {"id": "pca", "width_mm": 1, "figure_style": {}})):
        (tmp_path / name).mkdir()
        (tmp_path / name / "figure_manifest.json").write_text(json.dumps({"schema_version": "nf-rna.figure-manifest.v1", "figures": [record]}), encoding="utf-8")
    assert review.compat(tmp_path / "a", tmp_path / "b")["added_keys"] == {"./pca": ["figure_style"]}
    assert review.main(["--compat", str(tmp_path / "a"), str(tmp_path / "b")]) == 0
    assert review.main(["--compat", str(tmp_path / "b"), str(tmp_path / "a")]) == 1


def test_pdf_layout_repeat_groups_pdfs_by_content_not_by_date(tmp_path):
    import zlib

    repeat = _load("pdf_layout_repeat")

    def pdf(name: str, height: str, date: str) -> Path:
        path = tmp_path / name
        stream = zlib.compress(f"0 0 100 {height} re f".encode())
        path.write_bytes(b"%PDF-1.5\n1 0 obj\n<< /CreationDate (D:" + date.encode() + b") >>\nstream\n" + stream + b"\nendstream\nendobj\n")
        return path

    groups = repeat.group_by_content([pdf("a.pdf", "200", "20261015120000"), pdf("b.pdf", "200", "20261015120001"), pdf("c.pdf", "170", "20261015120000")])
    assert sorted(len(files) for files in groups.values()) == [1, 2]
