"""Phase 4 remediation: PLOS compact width and L2 crash, PCA label floor, approved palettes, manifest accuracy."""

from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import pytest

from conftest import require_r_packages
from test_figures import run_r

ROOT = Path(__file__).resolve().parents[1]
TOOLS = ROOT / "tools" / "figure_review"
GROUPS_1_8 = ["#4E78A7", "#CC5A6A", "#59A14F", "#843C39", "#56B4E9", "#BCBD22", "#6A3D9A", "#8C6D31"]
BATCH = ["#393B79", "#585858", "#17BECF", "#E6AB02", "#CC79A7", "#A6761D", "#7F7F7F", "#E08214"]


def standard(profile: str = "general") -> str:
    return f'base <- nf_rna_figure_context(list(figure_style=list(preset="standard", profile="{profile}")))\nctx <- nf_rna_figure_scope(base, "pca")\n'


def _load(name: str):
    import sys

    sys.path.insert(0, str(TOOLS))
    spec = importlib.util.spec_from_file_location(name, TOOLS / f"{name}.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


# ------------------------------------------------------------------ palettes


def test_group_slots_one_to_eight_are_the_approved_colours_and_the_extension_is_deterministic(tmp_path):
    result = run_r(standard() + '''
n <- function(k) unname(nf_rna_group_colours(paste0("G", seq_len(k)), context = ctx))
lab <- function(h) grDevices::convertColor(t(grDevices::col2rgb(h)) / 255, from = "sRGB", to = "Lab")
c12 <- n(12); L <- lab(c12)
dist <- as.matrix(stats::dist(L))
cat(toJSON(list(c8 = n(8), c9 = n(9), c12 = c12, again = identical(c12, n(12)), min_new_to_old = sapply(9:12, function(i) min(dist[i, 1:8])),
  min_all = min(dist[upper.tri(dist)]), name9 = nf_rna_categorical_palette_name(9, ctx), name8 = nf_rna_categorical_palette_name(8, ctx),
  legacy9 = unname(nf_rna_categorical_colors(paste0("G", 1:9))), legacy_name = nf_rna_categorical_palette_name(9)), auto_unbox = FALSE))''', tmp_path)
    assert result["c8"] == GROUPS_1_8
    assert result["c9"][:8] == GROUPS_1_8 and result["c12"][:9] == result["c9"]  # prefix-stable: adding a level never recolours another
    assert result["again"] == [True] and len(set(result["c12"])) == 12
    assert min(result["min_new_to_old"]) >= 15 and result["min_all"][0] >= 15  # CIE Lab distance of the extension colours to every earlier one
    assert result["name8"] == ["standard groups 1-8 (approved)"] and "maximin extension from slot 9" in result["name9"][0]
    assert result["legacy_name"] == ["HCL qualitative 'Dark 3'"] and len(set(result["legacy9"])) == 9  # legacy behaviour is untouched


def test_the_batch_palette_is_the_reported_candidate_with_a_darker_neutral_and_carries_its_approved_name(tmp_path):
    result = run_r(standard() + '''
cat(toJSON(list(batch = unname(nf_rna_batch_colours(paste0("B", 1:8), context = ctx)), name = nf_rna_categorical_palette_name(3, ctx, batch = TRUE),
  name9 = nf_rna_categorical_palette_name(9, ctx, batch = TRUE), b9 = unname(nf_rna_batch_colours(paste0("B", 1:9), context = ctx)),
  legacy_name = nf_rna_categorical_palette_name(3, batch = TRUE)), auto_unbox = FALSE))''', tmp_path)
    assert result["batch"] == BATCH and "#BAB0AC" not in result["batch"]
    assert result["name"][0] == "nf-rna Standard Batch Palette v1" and result["legacy_name"] == ["Okabe-Ito"]
    assert result["b9"][:8] == BATCH and len(set(result["b9"])) == 9 and "maximin extension" in result["name9"][0]


@pytest.mark.parametrize("profile,expected", [("general", 6), ("nature", 6), ("plos", 8), ("elsevier", 7)])
def test_pca_labels_are_six_points_but_never_below_the_profile_minimum(tmp_path, profile, expected):
    result = run_r(standard(profile) + 'cat(toJSON(list(standard = nf_rna_pca_label_pt(ctx), legacy = nf_rna_pca_label_pt(nf_rna_figure_context(list())), legacy_min = nf_rna_figure_profile("' + profile + '")$min_pt), auto_unbox = TRUE))', tmp_path)
    assert result["standard"] == expected


# ------------------------------------------------------------------ end to end through L1 and L2


@pytest.fixture(scope="module")
def runs(tmp_path_factory):
    require_r_packages("DESeq2", "ggplot2", "pheatmap", "ggrepel", "jsonlite", "systemfonts")
    base = tmp_path_factory.mktemp("remediation")
    _load("make_demo_project").build(base / "project", 4, 600, 20261015, "", 3)
    generate = _load("generate_runs")
    out = {}
    for profile in ("general", "nature", "plos", "elsevier"):
        out[profile] = base / profile
        out[profile].mkdir()
        # the L2 DEG heatmap crashed here for the plos profile (no single column) before the fix
        assert generate.run_l1_l2(base / "project", out[profile], {"preset": "standard", "profile": profile}) == "SUCCESS"
    return out


PROFILES = {"general": (85, 180, 6), "nature": (89, 183, 5), "plos": (85, 190.5, 8), "elsevier": (90, 190, 7)}


@pytest.mark.parametrize("profile", PROFILES)
def test_standard_l2_runs_in_every_profile_with_the_expected_compact_widths_and_text_floor(runs, profile):
    single, maximum, floor = PROFILES[profile]
    audit = _load("profile_audit")
    seen = {}
    for manifest in sorted(runs[profile].rglob("figure_manifest.json")):
        data = json.loads(manifest.read_text(encoding="utf-8"))
        for figure in data["figures"]:
            pdf = next(manifest.parent / item["path"] for item in figure["files"] if item["format"] == "pdf")
            sizes = audit.text_sizes(pdf)
            assert sizes["implied_min_pt"] >= floor - 1e-6, (figure["id"], sizes)  # every text run, PCA labels included
            assert figure["width_mm"] <= maximum + 0.05 and figure["height_mm"] <= data["profile"]["max_height_mm"] + 0.05
            seen[figure["kind"]] = figure["width_mm"]
    assert abs(seen["volcano"] - single) < 0.5 and abs(seen["PCA"] - single) < 0.5
    assert "DEG heatmap" in seen  # the figure that failed for plos


def test_manifest_level_font_and_colour_policy_describe_what_was_drawn(runs):
    for profile, root in runs.items():
        for manifest in root.rglob("figure_manifest.json"):
            data = json.loads(manifest.read_text(encoding="utf-8"))
            assert data["schema_version"] == "nf-rna.figure-manifest.v1"
            assert data["font"]["family"] == "Arimo" and len(data["font"]["sha256"]) == 64 and data["font"]["fallback_used"] is False
            assert data["figure_style"]["font"]["sha256"] == data["font"]["sha256"] and data["profile"]["name"] == profile
            policy = data["colour_policy"]
            assert policy["applied_preset"] == "standard" and policy["direction"]["Up"] == "#B2182B" and policy["categorical_colours"] == GROUPS_1_8
            assert policy["batch"]["colours"] == BATCH and policy["batch"]["name"] == "nf-rna Standard Batch Palette v1"
            assert policy["diverging"]["colours"] == {"low": "#4E78A7", "mid": "#F7F7F7", "high": "#CC5A6A"}
            for figure in data["figures"]:
                assert figure["figure_style"]["font_sha256"] == data["font"]["sha256"]
                if figure["kind"] == "PCA" and figure["encoding"]["colour"] == "batch":
                    assert "nf-rna Standard Batch Palette v1" in figure["encoding"]["palette"]  # the batch palette is named as such, not as the group palette


def test_legacy_manifests_keep_their_original_keys(tmp_path):
    result = run_r('''
ctx <- nf_rna_figure_context(list(figure_style=list(preset="legacy", profile="general")))
dir.create("o"); nf_rna_write_figure_manifest(ctx, "o")
m <- jsonlite::fromJSON("o/figure_manifest.json", simplifyVector = FALSE)
cat(toJSON(list(keys = names(m), policy = names(m$colour_policy), font = m$font$family), auto_unbox = TRUE))''', tmp_path)
    assert result["keys"] == ["schema_version", "profile", "font", "typography", "colour_policy", "devices", "r_version", "package_versions", "figures"]
    assert result["policy"] == ["categorical", "direction", "diverging", "sequential", "missing"]
