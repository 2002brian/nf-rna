"""Independent-filtering alpha gate.

nf-rna calls ``results(alpha = thresholds.padj)`` so that DESeq2's independent
filtering is optimised at the declared DEG FDR (v1.6.1, docs/scientific_contract.md).
On the broader design fixtures DESeq2 happens to choose the same baseMean
cut-off for alpha 0.05 and 0.1, so an alpha regression would go unnoticed there.
This fixture has many low-expression genes; the oracle first proves that alpha
0.05 and 0.1 select different cut-offs on it, and nf-rna must match the 0.05 one.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from _science import Design, by_gene, close, median_of_ratios, read_tsv, require_downstream_r_runtime, run_l2, run_reference, simulate, write_project

CONTRAST = "Treatment_vs_Control"


@pytest.fixture(scope="module")
def run(tmp_path_factory):
    require_downstream_r_runtime("jsonlite", "DESeq2", "apeglm", "ggplot2", "pheatmap")
    design = Design(["C1", "C2", "C3", "T1", "T2", "T3"], {"condition": ["Control"] * 3 + ["Treatment"] * 3},
                    [1.0, 0.9, 1.1, 1.0, 1.2, 0.8], "~ condition", {"condition": "categorical"})
    # 1700 genes keep >= 1000 genes above a mean normalized count of 5, which DESeq2::vst() in L1
    # requires (l1_analysis.R only falls back to varianceStabilizingTransformation below 1000 rows).
    genes = [f"G{index:04d}" for index in range(1, 1701)]
    effects = {gene: (1.5 if i % 8 == 0 else -1.5) for i, gene in enumerate(genes) if i % 8 in (0, 1)}
    # Expression from 1 to 2000 mean counts: many genes sit where independent filtering decides.
    sim = simulate(genes, design, seed=8, effects=effects, base_range=(1.0, 2000.0))
    kept = {gene: row for gene, row in sim.counts.items() if sum(row) >= 10}
    factors = median_of_ratios(kept)
    assert sum(sum(c / f for c, f in zip(row, factors)) / len(row) > 5 for row in kept.values()) >= 1000
    root = write_project(tmp_path_factory.mktemp("independent_filtering") / "independent_filtering", sim, design, enrichment=[], human_annotation=False)
    l2 = run_l2(root)
    oracle = {}
    for alpha in ("0.05", "0.1"):
        prefix = root / f"oracle_alpha_{alpha}"
        run_reference("deseq2_reference.R", root / "input" / "counts.csv", root / "metadata.csv", "~ condition", "condition", alpha, prefix)
        oracle[alpha] = (by_gene(read_tsv(Path(f"{prefix}_results.tsv"))), json.loads(Path(f"{prefix}_info.json").read_text(encoding="utf-8")))
    return l2, oracle


def _na(rows):
    return {gene for gene, row in rows.items() if row["padj"] == "NA"}


def test_fixture_discriminates_between_alpha_005_and_01(run):
    _l2, oracle = run
    (rows_005, info_005), (rows_01, info_01) = oracle["0.05"], oracle["0.1"]
    assert abs(info_005["filter_threshold"] - info_01["filter_threshold"]) > 0.1
    assert _na(rows_005) != _na(rows_01)


def test_independent_filtering_is_optimised_at_the_declared_padj_threshold(run):
    l2, oracle = run
    contrast = l2 / "contrasts" / CONTRAST
    observed = by_gene(read_tsv(contrast / "all_genes.tsv"))
    expected, info = oracle["0.05"]
    summary = json.loads((contrast / "backend_summary.json").read_text(encoding="utf-8"))
    # backend_summary.json carries jsonlite's default 4 decimals.
    assert summary["independent_filtering_threshold"] == pytest.approx(info["filter_threshold"], abs=5e-5)
    assert summary["independent_filtering_threshold"] != pytest.approx(oracle["0.1"][1]["filter_threshold"], abs=0.05)
    assert _na(observed) == _na(expected)
    for gene, row in expected.items():
        assert close(None if observed[gene]["padj"] == "NA" else float(observed[gene]["padj"]),
                     None if row["padj"] == "NA" else float(row["padj"]), rel=1e-8), gene
    # DESeq2's rule: padj is NA exactly for genes whose filter statistic (baseMean) is below the chosen cut-off.
    threshold = info["filter_threshold"]
    assert all(float(observed[gene]["baseMean"]) < threshold for gene in _na(observed))
    assert all(float(row["baseMean"]) >= threshold for gene, row in observed.items() if gene not in _na(observed))
    metadata = json.loads((contrast / "contrast_metadata.json").read_text(encoding="utf-8"))
    assert metadata["deg_definition"]["independent_filtering"]["alpha"] == 0.05
