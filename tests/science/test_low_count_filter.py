"""Low-count filtering boundary gate (L1 and L2).

The documented contract (docs/scientific_contract.md) removes all-zero genes and
genes whose total count over all analysed samples is < 10.  The DESeq2 design
fixtures contain no gene near that boundary, so a changed threshold would pass
them; this focused fixture puts hand-written genes on both sides of every
boundary and checks *which* genes survive, in L1 (normalized counts, VST) and in
L2 (all_genes.tsv), via the public entry points.  The expected set is derived
from the hand-written totals here, not from nf-rna.
"""

from __future__ import annotations

import csv
import json

import pytest

from _science import Design, by_gene, read_tsv, require_downstream_r_runtime, run_l2, simulate, write_project

SAMPLES = ["C1", "C2", "C3", "T1", "T2", "T3"]
# name -> per-sample counts; the name records the total.
BOUNDARY = {
    "ZERO_00": [0, 0, 0, 0, 0, 0],
    "LOW_03": [1, 0, 1, 0, 1, 0],
    "LOW_04": [1, 1, 0, 1, 1, 0],
    "MID_05": [1, 1, 1, 1, 1, 0],
    "MID_07": [2, 1, 1, 1, 1, 1],
    "MID_09": [2, 2, 1, 2, 1, 1],
    "EXACT_10": [2, 2, 1, 2, 2, 1],
    "ABOVE_11": [2, 2, 2, 2, 2, 1],
    "ABOVE_25": [5, 4, 4, 4, 4, 4],
    # Total >= 10 although most samples are zero: the rule is the total, not per-sample presence.
    "SPARSE_12": [12, 0, 0, 0, 0, 0],
}
EXPECTED_KEPT = {"EXACT_10", "ABOVE_11", "ABOVE_25", "SPARSE_12"}
EXPECTED_REMOVED = {"ZERO_00", "LOW_03", "LOW_04", "MID_05", "MID_07", "MID_09"}


@pytest.fixture(scope="module")
def filtered(tmp_path_factory):
    require_downstream_r_runtime("jsonlite", "DESeq2", "apeglm", "ggplot2", "pheatmap")
    assert {name for name, row in BOUNDARY.items() if sum(row) >= 10} == EXPECTED_KEPT
    assert all(int(name.rsplit("_", 1)[1]) == sum(row) for name, row in BOUNDARY.items())
    design = Design(SAMPLES, {"condition": ["Control"] * 3 + ["Treatment"] * 3}, [1.0] * 6, "~ condition", {"condition": "categorical"})
    background = [f"G{index:04d}" for index in range(1, 201)]
    sim = simulate(background, design, seed=101, effects={gene: 2.0 for gene in background[:20]})
    assert min(sum(row) for row in sim.counts.values()) >= 10
    for name, row in BOUNDARY.items():
        sim.genes.append(name)
        sim.counts[name] = row
        sim.truth_lfc[name] = 0.0
    root = write_project(tmp_path_factory.mktemp("filter") / "low_count_filter", sim, design, enrichment=[], human_annotation=False)
    return root, run_l2(root), sim


def test_l2_retains_exactly_the_genes_with_total_count_of_at_least_ten(filtered):
    _root, l2, sim = filtered
    genes = set(by_gene(read_tsv(l2 / "contrasts" / "Treatment_vs_Control" / "all_genes.tsv")))
    assert genes & set(BOUNDARY) == EXPECTED_KEPT, f"retained boundary genes: {sorted(genes & set(BOUNDARY))}"
    assert not genes & EXPECTED_REMOVED
    assert genes == {gene for gene, row in sim.counts.items() if sum(row) >= 10}
    summary = json.loads((l2 / "contrasts" / "Treatment_vs_Control" / "backend_summary.json").read_text(encoding="utf-8"))
    assert (summary["input_genes"], summary["filtered_genes"], summary["retained_genes"]) == (len(sim.counts), len(EXPECTED_REMOVED), len(genes))


def test_l1_normalized_counts_and_vst_use_the_same_retained_genes(filtered):
    _root, l2, _sim = filtered
    l1 = l2.parent / "l1"
    for name in ("normalized_counts.tsv", "vst.tsv"):
        with (l1 / name).open(encoding="utf-8", newline="") as handle:
            genes = {row[0] for row in list(csv.reader(handle, delimiter="\t"))[1:]}
        assert genes & set(BOUNDARY) == EXPECTED_KEPT, (name, sorted(genes & set(BOUNDARY)))
    library = read_tsv(l1 / "library_size_qc.tsv")
    removed_total = {sample: sum(BOUNDARY[name][index] for name in EXPECTED_REMOVED) for index, sample in enumerate(SAMPLES)}
    for row in library:
        assert int(float(row["input_total"])) - int(float(row["retained_total"])) == removed_total[row["sample_id"]], row


def test_the_documented_threshold_is_the_one_every_execution_path_receives():
    from rnaseq.downstream import L1_FILTER
    from rnaseq.workflow_support import FILTER

    assert L1_FILTER["minimum_total_count"] == FILTER["minimum_total_count"] == 10
