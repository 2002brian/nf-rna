"""Independent DESeq2 regression gate: two-group, paired and batch-adjusted designs.

Each design is simulated with a planted condition effect.  nf-rna's L2 output is
then checked against three independent sources:

* the simulation truth (direction and recovery of planted effects; approximate);
* closed-form calculations implemented in ``_science`` (filtering, size factors,
  paired log-ratio signs);
* ``oracles/deseq2_reference.R``, which calls DESeq2 directly from the
  documented contract without sourcing nf-rna's R scripts (numerically exact up to
  text round-trip, because both run the same locked DESeq2).
"""

from __future__ import annotations

import json
import math
from pathlib import Path

import pytest

from _science import (
    Design, by_gene, close, gene_panel, median_of_ratios, read_tsv, require_downstream_r_runtime, run_l2,
    run_reference, simulate, write_project,
)

CONTRAST = "Treatment_vs_Control"
STAT_COLUMNS = ("baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")
# Same locked DESeq2 on both sides: differences can only come from the 15-significant-digit text round-trip.
EXACT_REL = 1e-8
SHRINK_ABS = 1e-6  # apeglm is a numerical optimiser; identical inputs give identical output, allow last-digit noise


def _synthetic_genes(n: int) -> list[str]:
    return [f"G{index:04d}" for index in range(1, n + 1)]


def _effects(genes: list[str], up: int, down: int, size: float) -> dict[str, float]:
    return {**{gene: size for gene in genes[:up]}, **{gene: -size for gene in genes[up:up + down]}}


def _two_group(low_expression: int = 0) -> tuple[Design, dict, dict[str, str]]:
    """``low_expression`` null genes get mean counts of 1-4, so independent filtering gives them NA padj."""

    panel = gene_panel()
    genes = [row["ensembl_id"] for row in panel]
    roles = {row["ensembl_id"]: row["role"] for row in panel}
    effects = {gene: 2.0 for gene in genes if roles[gene] == "up_term"}
    effects |= {gene: -2.0 for gene in genes if roles[gene] == "down_term"}
    samples = ["C1", "C2", "C3", "T1", "T2", "T3"]
    design = Design(samples, {"condition": ["Control"] * 3 + ["Treatment"] * 3}, [0.8, 1.1, 1.3, 0.9, 1.2, 0.7],
                    "~ condition", {"condition": "categorical"})
    nulls = [gene for gene in genes if roles[gene] == "null"]
    low = {gene: 1.0 + 3.0 * index / max(low_expression - 1, 1) for index, gene in enumerate(nulls[:low_expression])}
    return design, simulate(genes, design, seed=11, effects=effects, roles=roles, base_override=low), roles


def _paired() -> tuple[Design, object]:
    genes = _synthetic_genes(600)
    samples = ["N1", "N2", "N3", "N4", "T1", "T2", "T3", "T4"]
    design = Design(samples, {"patient": ["P1", "P2", "P3", "P4"] * 2, "condition": ["Control"] * 4 + ["Treatment"] * 4},
                    [1.0, 0.8, 1.2, 0.9, 1.1, 0.7, 1.3, 1.0], "~ patient + condition",
                    {"patient": "categorical", "condition": "categorical"}, pair_id="patient", design_type="paired_two_group")
    return design, simulate(genes, design, seed=23, effects=_effects(genes, 30, 30, 1.3), block_sd=1.5, block_column="patient")


def _batch() -> tuple[Design, object]:
    genes = _synthetic_genes(600)
    samples = ["C1", "C2", "C3", "C4", "T1", "T2", "T3", "T4"]
    design = Design(samples, {"batch": ["B1", "B1", "B2", "B2", "B1", "B1", "B2", "B2"], "condition": ["Control"] * 4 + ["Treatment"] * 4},
                    [1.0, 0.8, 1.2, 0.9, 1.1, 0.7, 1.3, 1.0], "~ batch + condition",
                    {"batch": "categorical", "condition": "categorical"})
    return design, simulate(genes, design, seed=37, effects=_effects(genes, 30, 30, 1.3), block_sd=1.5, block_column="batch")


class Analysed:
    def __init__(self, root: Path, design: Design, sim, roles: dict[str, str] | None = None):
        self.root, self.design, self.sim, self.roles = root, design, sim, roles or {}
        self.l2 = run_l2(root)
        self.contrast = self.l2 / "contrasts" / CONTRAST
        self.all = by_gene(read_tsv(self.contrast / "all_genes.tsv"))
        self.significant = by_gene(read_tsv(self.contrast / "significant.tsv"))
        self.up = by_gene(read_tsv(self.contrast / "upregulated.tsv"))
        self.down = by_gene(read_tsv(self.contrast / "downregulated.tsv"))
        self.summary = json.loads((self.contrast / "backend_summary.json").read_text(encoding="utf-8"))
        self.references: dict[str, tuple[dict, dict]] = {}

    def reference(self, formula: str) -> tuple[dict[str, dict[str, str]], dict]:
        if formula not in self.references:
            prefix = self.root / f"reference_{len(self.references)}"
            categorical = ",".join(name for name in self.design.columns if formula.count(name))
            run_reference("deseq2_reference.R", self.root / "input" / "counts.csv", self.root / "metadata.csv", formula, categorical, 0.05, prefix)
            self.references[formula] = (by_gene(read_tsv(Path(f"{prefix}_results.tsv"))),
                                        json.loads(Path(f"{prefix}_info.json").read_text(encoding="utf-8")))
        return self.references[formula]

    def recall(self, rows: dict[str, dict[str, str]]) -> float:
        planted = [gene for gene, lfc in self.sim.truth_lfc.items() if lfc]
        found = [gene for gene in planted if gene in rows and _deg(rows[gene])]
        return len(found) / len(planted)


def _deg(row: dict[str, str]) -> bool:
    padj, lfc = row["padj"], row["log2FoldChange"]
    return padj != "NA" and lfc != "NA" and float(padj) < 0.05 and abs(float(lfc)) >= 1


@pytest.fixture(scope="module")
def designs(tmp_path_factory) -> dict[str, Analysed]:
    require_downstream_r_runtime("jsonlite", "DESeq2", "apeglm", "ggplot2", "pheatmap")
    built: dict[str, Analysed] = {}
    design, sim, roles = _two_group()
    built["two_group"] = Analysed(write_project(tmp_path_factory.mktemp("two") / "two_group", sim, design, enrichment=[], human_annotation=True), design, sim, roles)
    for name, factory in (("paired", _paired), ("batch", _batch)):
        design, sim = factory()
        built[name] = Analysed(write_project(tmp_path_factory.mktemp(name) / name, sim, design, enrichment=[], human_annotation=False), design, sim)
    return built


@pytest.mark.parametrize("name", ["two_group", "paired", "batch"])
def test_filtering_and_size_factors_match_an_independent_median_of_ratios(designs, name):
    analysed = designs[name]
    retained = {gene: row for gene, row in analysed.sim.counts.items() if sum(row) >= 10}
    assert set(analysed.all) == set(retained)
    assert analysed.summary["input_genes"] == len(analysed.sim.counts)
    assert analysed.summary["retained_genes"] == len(retained)
    expected = median_of_ratios(retained)
    library = read_tsv(analysed.l2.parent / "l1" / "library_size_qc.tsv")
    assert [row["sample_id"] for row in library] == analysed.design.samples
    for row, factor in zip(library, expected):
        assert math.isclose(float(row["size_factor"]), factor, rel_tol=1e-9), row


@pytest.mark.parametrize("name", ["two_group", "paired", "batch"])
def test_l2_statistics_equal_a_direct_deseq2_reference(designs, name):
    analysed = designs[name]
    reference, info = analysed.reference(analysed.design.formula)
    assert set(reference) == set(analysed.all)
    for gene, expected in reference.items():
        observed = analysed.all[gene]
        for column in STAT_COLUMNS:
            assert close(_n(observed[column]), _n(expected[column]), rel=EXACT_REL), (name, gene, column, observed[column], expected[column])
        assert close(_n(observed["log2FoldChange_shrunk"]), _n(expected["log2FoldChange_shrunk"]), rel=0, abs_=SHRINK_ABS), (name, gene)
    # backend_summary.json is written by jsonlite with its default 4 decimal digits.
    assert analysed.summary["independent_filtering_threshold"] == pytest.approx(info["filter_threshold"], abs=5e-5)
    assert analysed.summary["padj_na"] == sum(row["padj"] == "NA" for row in reference.values())
    assert analysed.summary["pvalue_na"] == sum(row["pvalue"] == "NA" for row in reference.values())


def test_cooks_filtering_applies_to_replicated_groups_but_not_to_the_paired_design(designs):
    # DESeq2 computes Cook's distance only for samples with >= 3 replicates of the
    # same model-matrix row; ~ patient + condition has none, so nothing is flagged.
    assert designs["two_group"].reference("~ condition")[1]["cooks_applied"] is True
    paired_reference, paired_info = designs["paired"].reference("~ patient + condition")
    assert paired_info["cooks_applied"] is False
    assert designs["paired"].summary["pvalue_na"] == 0 == sum(row["pvalue"] == "NA" for row in paired_reference.values())


def test_contrast_direction_and_recovery_of_planted_effects(designs):
    analysed = designs["two_group"]
    up_genes = {gene for gene, lfc in analysed.sim.truth_lfc.items() if lfc > 0}
    down_genes = {gene for gene, lfc in analysed.sim.truth_lfc.items() if lfc < 0}
    # Direction is exact: a planted Treatment-high gene can never be called down.
    assert not (up_genes & set(analysed.down)) and not (down_genes & set(analysed.up))
    # Recovery and false discoveries are statistical (approximate) guarantees with generous margins.
    assert len(up_genes & set(analysed.up)) / len(up_genes) >= 0.8
    assert len(down_genes & set(analysed.down)) / len(down_genes) >= 0.8
    false = set(analysed.significant) - up_genes - down_genes
    assert len(false) <= 0.1 * len(analysed.significant)
    for gene in (up_genes | down_genes) & set(analysed.significant):
        assert math.copysign(1, float(analysed.all[gene]["log2FoldChange"])) == math.copysign(1, analysed.sim.truth_lfc[gene])


def test_paired_design_is_fitted_and_its_direction_matches_manual_paired_differences(designs):
    analysed = designs["paired"]
    provenance = json.loads((analysed.l2 / "scientific_provenance.json").read_text(encoding="utf-8"))["design"]
    columns = provenance["model_matrix"]["column_names"]
    assert columns == ["(Intercept)", "patientP2", "patientP3", "patientP4", "conditionTreatment"]
    assert provenance["model_matrix"]["full_rank"] is True
    # The fitted model is the paired one, not an unpaired fit of the same data ...
    unpaired, _ = analysed.reference("~ condition")
    differing = sum(not close(_n(analysed.all[g]["stat"]), _n(unpaired[g]["stat"]), rel=1e-6) for g in analysed.all)
    assert differing > 0.9 * len(analysed.all)
    # ... and, with a strong patient effect, pairing recovers more planted effects.
    assert analysed.recall(analysed.all) > analysed.recall(unpaired)
    # Independent direction check: mean over patients of log2((T + 0.5) / (N + 0.5))
    # on median-of-ratios normalised counts.  GLM estimates and empirical paired
    # differences are not identical, so agreement is required for >= 98% of DEGs
    # and for every planted gene that is called.
    retained = {g: row for g, row in analysed.sim.counts.items() if sum(row) >= 10}
    factors = median_of_ratios(retained)
    agree = 0
    for gene, row in analysed.significant.items():
        values = [count / factor for count, factor in zip(analysed.sim.counts[gene], factors)]
        paired = sum(math.log2((values[4 + i] + 0.5) / (values[i] + 0.5)) for i in range(4)) / 4
        same = math.copysign(1, paired) == math.copysign(1, float(row["log2FoldChange"]))
        agree += same
        if analysed.sim.truth_lfc[gene]:
            assert same, gene
    assert agree >= 0.98 * len(analysed.significant)


def test_batch_term_is_fitted_and_improves_recovery(designs):
    analysed = designs["batch"]
    provenance = json.loads((analysed.l2 / "scientific_provenance.json").read_text(encoding="utf-8"))["design"]
    assert provenance["model_matrix"]["column_names"] == ["(Intercept)", "batchB2", "conditionTreatment"]
    unadjusted, _ = analysed.reference("~ condition")
    assert analysed.recall(analysed.all) >= analysed.recall(unadjusted)
    assert any(not close(_n(analysed.all[g]["stat"]), _n(unadjusted[g]["stat"]), rel=1e-6) for g in analysed.all)


@pytest.mark.parametrize("name", ["two_group", "paired", "batch"])
def test_deg_tables_follow_the_raw_lfc_rule_not_the_shrunken_lfc(designs, name):
    analysed = designs[name]
    expected = {gene for gene, row in analysed.all.items() if _deg(row)}
    assert set(analysed.significant) == expected
    assert set(analysed.up) == {g for g in expected if float(analysed.all[g]["log2FoldChange"]) > 0}
    assert set(analysed.down) == {g for g in expected if float(analysed.all[g]["log2FoldChange"]) < 0}
    shrunk = [(float(row["log2FoldChange"]), float(row["log2FoldChange_shrunk"])) for row in analysed.all.values()
              if row["log2FoldChange_shrunk"] != "NA"]
    assert len(shrunk) == len(analysed.all)  # apeglm succeeded for every gene
    # apeglm shrinks towards zero; allow numerical noise and rare non-monotone genes.
    assert sum(abs(s) <= abs(r) + 1e-6 for r, s in shrunk) >= 0.99 * len(shrunk)
    metadata = json.loads((analysed.contrast / "contrast_metadata.json").read_text(encoding="utf-8"))
    assert metadata["shrinkage"]["status"] == "success" and metadata["shrinkage"]["method"] == "apeglm"


@pytest.mark.parametrize("name", ["two_group", "paired"])
def test_figure_manifests_agree_with_the_deg_tables(designs, name):
    analysed = designs[name]
    manifest = json.loads((analysed.l2 / "figure_manifest.json").read_text(encoding="utf-8"))
    figures = {item["id"].rsplit("/", 1)[-1]: item for item in manifest["figures"] if CONTRAST in item["id"]}
    volcano = figures["volcano"]["metadata"] if "metadata" in figures["volcano"] else figures["volcano"]
    counts = volcano["counts"]
    assert counts["Up"] == len(analysed.up) and counts["Down"] == len(analysed.down)
    assert counts["Up"] + counts["Down"] + counts["Not significant"] + volcano["display"]["genes_not_drawn"] == len(analysed.all)
    heatmap_genes = [row["gene_id"] for row in read_tsv(analysed.contrast / "heatmap_genes.tsv")]
    ranked = sorted(analysed.significant, key=lambda gene: (float(analysed.significant[gene]["padj"]), gene))
    assert heatmap_genes == ranked[: len(heatmap_genes)]
    assert len(heatmap_genes) == min(len(ranked), analysed.summary["heatmap_top_n"])


def _n(value: str) -> float | None:
    return None if value in {"NA", ""} else float(value)
