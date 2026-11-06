"""Shared helpers for the independent scientific regression gate (tests/science).

Everything here is deliberately independent of nf-rna's R scripts: data are
simulated with the standard library, expectations are computed either from the
simulation truth, from closed-form statistics implemented here (median-of-ratios
size factors, hypergeometric tail, Benjamini-Hochberg, the GSEA running sum), or
from reference R scripts in ``tests/science/oracles`` that call DESeq2 and
AnnotationDbi directly and never source ``src/rnaseq/r``.

nf-rna itself is driven only through its public Python entry points
(``validate_project`` -> ``prepare_*`` -> ``execute_*``), so the control plane
(project validation, frozen configuration, backend configs) is exercised too.
"""

from __future__ import annotations

import csv
import hashlib
import json
import math
import os
import random
import subprocess
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable

import yaml

TESTS = Path(__file__).resolve().parents[1]
if str(TESTS) not in sys.path:
    sys.path.insert(0, str(TESTS))

from conftest import base_config, require_downstream_r_runtime  # noqa: E402

SCIENCE = Path(__file__).resolve().parent
FIXTURES = SCIENCE / "fixtures"
# Not "reference(s)/": .gitignore excludes those names to keep biological reference data out of Git.
ORACLES = SCIENCE / "oracles"
UP_TERM = "GO:0006260"    # DNA replication; 40 panel genes planted 4x up in Treatment
DOWN_TERM = "GO:0006954"  # inflammatory response; 40 panel genes planted 4x down
KEGG_UP, KEGG_DOWN = "hsa99001", "hsa99002"


# ----------------------------------------------------------------- simulation
def poisson(rng: random.Random, lam: float) -> int:
    """Exact Poisson draw: Knuth's method on chunks of at most 500 (no underflow)."""

    total = 0
    while lam > 0:
        chunk = min(lam, 500.0)
        lam -= chunk
        limit, k, product = math.exp(-chunk), 0, rng.random()
        while product > limit:
            k += 1
            product *= rng.random()
        total += k
    return total


def negative_binomial(rng: random.Random, mean: float, dispersion: float) -> int:
    """Gamma-Poisson mixture: Var = mean + dispersion * mean^2 (DESeq2's parametrisation)."""

    shape = 1.0 / dispersion
    return poisson(rng, rng.gammavariate(shape, mean / shape))


@dataclass
class Design:
    samples: list[str]
    columns: dict[str, list[str]]           # metadata columns other than sample_id
    size_factors: list[float]
    formula: str
    variables: dict[str, str]
    pair_id: str | None = None
    design_type: str = "two_group"


@dataclass
class Simulation:
    genes: list[str]
    counts: dict[str, list[int]]
    truth_lfc: dict[str, float]             # log2(Treatment/Control) planted per gene
    roles: dict[str, str] = field(default_factory=dict)


def simulate(genes: list[str], design: Design, *, seed: int, effects: dict[str, float],
             block_sd: float = 0.0, block_column: str | None = None, roles: dict[str, str] | None = None,
             base_range: tuple[float, float] = (30.0, 3000.0), base_override: dict[str, float] | None = None) -> Simulation:
    """Counts with a planted condition effect and an optional per-gene block effect.

    ``effects`` maps gene -> log2 fold change (Treatment / Control).  With
    ``block_column`` every gene gets an independent log-normal multiplier per block
    level (a patient or batch effect of standard deviation ``block_sd`` on log2).
    """

    rng = random.Random(seed)
    condition = design.columns["condition"]
    blocks = design.columns.get(block_column, []) if block_column else []
    counts: dict[str, list[int]] = {}
    for gene in genes:
        base = math.exp(rng.uniform(math.log(base_range[0]), math.log(base_range[1])))
        base = (base_override or {}).get(gene, base)
        dispersion = 0.02 + 2.0 / base
        block_effect = {level: 2 ** rng.gauss(0.0, block_sd) for level in sorted(set(blocks))}
        row = []
        for index, _sample in enumerate(design.samples):
            mean = base * design.size_factors[index]
            if condition[index] == "Treatment":
                mean *= 2 ** effects.get(gene, 0.0)
            if blocks:
                mean *= block_effect[blocks[index]]
            row.append(negative_binomial(rng, mean, dispersion))
        counts[gene] = row
    return Simulation(genes, counts, {gene: effects.get(gene, 0.0) for gene in genes}, roles or {})


def counts_csv(sim: Simulation, samples: list[str]) -> str:
    lines = ["gene_id," + ",".join(samples)]
    lines += [gene + "," + ",".join(str(value) for value in sim.counts[gene]) for gene in sim.genes]
    return "\n".join(lines) + "\n"


def metadata_csv(design: Design) -> str:
    names = list(design.columns)
    lines = ["sample_id," + ",".join(names)]
    lines += [sample + "," + ",".join(design.columns[name][index] for name in names) for index, sample in enumerate(design.samples)]
    return "\n".join(lines) + "\n"


def sha256_text(text: str) -> str:
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


# ------------------------------------------------------------------ projects
def gene_panel() -> list[dict[str, str]]:
    with (FIXTURES / "gene_panel.tsv").open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_project(root: Path, sim: Simulation, design: Design, *, enrichment: list[str], human_annotation: bool,
                  kegg_snapshot: str | None = None) -> Path:
    config = base_config()
    config["schema_version"] = "1.3"
    config["project"]["id"] = root.name
    design_block: dict[str, object] = {"type": design.design_type, "formula": design.formula, "variables": design.variables}
    if design.pair_id:
        design_block["pair_id"] = design.pair_id
    config["design"] = design_block
    config["analysis"] = {"enrichment": enrichment}
    if human_annotation:
        config["organism"] = {"species": "Homo sapiens"}
        kegg: dict[str, object] = {
            "resource_provider": "kegg_rest_snapshot",
            "ora": {"p_adjust_method": "BH", "min_gs_size": 10, "max_gs_size": 500},
            "gsea": {"minimum_ranked_genes": 50, "p_adjust_method": "BH", "min_gs_size": 10, "max_gs_size": 500, "seed": 1},
        }
        if kegg_snapshot:
            kegg["snapshot"] = kegg_snapshot
        config["annotation"] = {
            "organism": "Homo sapiens", "input_id_type": "ENSEMBL", "target_id_type": "ENTREZID",
            "mapping_warning_rate": 0.70, "minimum_mapping_rate": 0.50, "minimum_mapped_foreground": 5,
            "enrichment": {
                "go": {"p_adjust_method": "BH"},
                "gsea": {"minimum_ranked_genes": 50, "min_gs_size": 10, "max_gs_size": 500, "p_adjust_method": "BH", "seed": 1},
                "kegg": kegg,
            },
        }
    (root / "input").mkdir(parents=True)
    (root / "input" / "counts.csv").write_text(counts_csv(sim, design.samples), encoding="utf-8")
    (root / "metadata.csv").write_text(metadata_csv(design), encoding="utf-8")
    (root / "contrasts.csv").write_text("contrast_id,factor,numerator,denominator\nTreatment_vs_Control,condition,Treatment,Control\n", encoding="utf-8")
    (root / "project.yaml").write_text(yaml.safe_dump(config, sort_keys=False), encoding="utf-8")
    return root


def run_l2(root: Path) -> Path:
    from rnaseq.l2 import execute_l2, prepare_l2
    from rnaseq.validators import validate_project

    report = validate_project(root)
    assert report.is_valid, [issue.message for issue in report.errors]
    result = execute_l2(prepare_l2(report, run_id=None))
    assert result.status == "SUCCESS"
    return result.output_dir


def run_enrichment(root: Path) -> None:
    from rnaseq.go import execute_go, prepare_go
    from rnaseq.gsea import execute_gsea, prepare_gsea
    from rnaseq.kegg import execute_kegg, prepare_kegg
    from rnaseq.validators import validate_project

    report = validate_project(root)
    assert execute_go(prepare_go(report, run_id=None)).status in {"SUCCESS", "NO_SIGNIFICANT_TERMS"}
    assert execute_gsea(prepare_gsea(report, run_id=None)).status in {"SUCCESS", "NO_SIGNIFICANT_TERMS"}
    for mode in ("ora", "gsea"):
        assert execute_kegg(prepare_kegg(report, run_id=None, mode=mode)).status in {"SUCCESS", "NO_SIGNIFICANT_TERMS"}


def write_kegg_snapshot(cache_home: Path, snapshot_id: str, pathways: dict[str, list[str]], names: dict[str, str]) -> None:
    """An offline KEGG snapshot in nf-rna's documented layout (TERM2GENE/TERM2NAME + SHA-256 manifest)."""

    directory = cache_home / "nf-rna" / "genesets" / "kegg" / "hsa" / snapshot_id
    directory.mkdir(parents=True)
    term2gene = "term\tgene\n" + "".join(f"{term}\t{gene}\n" for term in sorted(pathways) for gene in sorted(set(pathways[term]), key=int))
    term2name = "term\tname\n" + "".join(f"{term}\t{names[term]}\n" for term in sorted(names))
    (directory / "TERM2GENE.tsv").write_text(term2gene, encoding="utf-8")
    (directory / "TERM2NAME.tsv").write_text(term2name, encoding="utf-8")
    files = {name: {"sha256": hashlib.sha256((directory / name).read_bytes()).hexdigest(), "bytes": (directory / name).stat().st_size}
             for name in ("TERM2GENE.tsv", "TERM2NAME.tsv")}
    manifest = {"schema": "nf-rna.kegg-snapshot.v1", "snapshot_id": snapshot_id, "organism": "hsa",
                "retrieved_at_utc": "2026-01-01T00:00:00Z", "kegg_release": "synthetic test fixture", "files": files,
                "pathways": len(pathways), "genes": len({gene for genes in pathways.values() for gene in genes})}
    (directory / "manifest.json").write_text(json.dumps(manifest, indent=1), encoding="utf-8")


# ------------------------------------------------------------- reference R
def run_reference(script: str, *arguments: str | Path) -> None:
    require_downstream_r_runtime()
    result = subprocess.run(["Rscript", str(ORACLES / script), *map(str, arguments)], capture_output=True, text=True, check=False)
    assert result.returncode == 0, f"{script} failed:\n{result.stdout}\n{result.stderr}"


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def by_gene(rows: Iterable[dict[str, str]], key: str = "gene_id") -> dict[str, dict[str, str]]:
    return {row[key]: row for row in rows}


def number(value: str) -> float | None:
    return None if value in {"", "NA", "NaN"} else float(value)


def close(a: float | None, b: float | None, *, rel: float, abs_: float = 0.0) -> bool:
    if a is None or b is None:
        return a is None and b is None
    return math.isclose(a, b, rel_tol=rel, abs_tol=abs_)


# ------------------------------------------------- closed-form statistics
def median_of_ratios(counts: dict[str, list[int]]) -> list[float]:
    """DESeq2 estimateSizeFactors (Anders & Huber 2010): genes with no zero, geometric-mean reference."""

    rows = [row for row in counts.values() if all(value > 0 for value in row)]
    logs = [[math.log(value) for value in row] for row in rows]
    reference = [sum(row) / len(row) for row in logs]
    factors = []
    for column in range(len(rows[0])):
        ratios = sorted(row[column] - ref for row, ref in zip(logs, reference))
        middle = len(ratios) // 2
        median = ratios[middle] if len(ratios) % 2 else (ratios[middle - 1] + ratios[middle]) / 2
        factors.append(math.exp(median))
    return factors


def log_comb(n: int, k: int) -> float:
    return math.lgamma(n + 1) - math.lgamma(k + 1) - math.lgamma(n - k + 1)


def hypergeometric_upper(k: int, universe: int, set_size: int, drawn: int) -> float:
    """P(X >= k) for X ~ Hypergeometric(universe, set_size, drawn); the one-sided ORA test."""

    top = min(set_size, drawn)
    logs = [log_comb(set_size, i) + log_comb(universe - set_size, drawn - i) - log_comb(universe, drawn) for i in range(k, top + 1)]
    peak = max(logs)
    return min(1.0, math.exp(peak) * sum(math.exp(value - peak) for value in logs))


def benjamini_hochberg(pvalues: list[float]) -> list[float]:
    n = len(pvalues)
    order = sorted(range(n), key=lambda i: pvalues[i], reverse=True)
    adjusted, running = [0.0] * n, 1.0
    for rank_from_top, index in enumerate(order):
        rank = n - rank_from_top
        running = min(running, pvalues[index] * n / rank)
        adjusted[index] = min(running, 1.0)
    return adjusted


def gsea_enrichment_score(ranked: list[tuple[str, float]], members: set[str]) -> float:
    """Subramanian et al. (2005) running-sum ES with weight p = 1 over a decreasing ranked list."""

    hits = [(gene, abs(stat)) for gene, stat in ranked if gene in members]
    total = sum(weight for _gene, weight in hits)
    miss_step = 1.0 / (len(ranked) - len(hits))
    running, best = 0.0, 0.0
    for gene, stat in ranked:
        running = running + abs(stat) / total if gene in members else running - miss_step
        if abs(running) > abs(best):
            best = running
    return best
