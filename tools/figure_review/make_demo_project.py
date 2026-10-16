#!/usr/bin/env python3
"""Write a deterministic synthetic two-group project with a batch for figure review.

The shipped example is tiny and nearly degenerate (a few dozen genes with perfectly
separated groups).  This generator builds a project with a realistic volcano shape and
a batch-coloured PCA, using only the standard library (a seeded negative-binomial
simulation).  It is a review fixture only: the counts are not real data.

    make_demo_project.py OUT_DIR [--samples-per-group 6] [--genes 4000] [--seed 20261015]
"""

from __future__ import annotations

import argparse
import math
import random
import sys
from pathlib import Path


def poisson(rng: random.Random, mean: float) -> int:
    if mean < 30:
        limit, k, p = math.exp(-mean), 0, 1.0
        while True:
            p *= rng.random()
            if p <= limit:
                return k
            k += 1
    return max(0, round(rng.gauss(mean, math.sqrt(mean))))


def build(out: Path, per_group: int, genes: int, seed: int, prefix: str = "", dose_levels: int = 0) -> None:
    rng = random.Random(seed)
    out.mkdir(parents=True, exist_ok=True)
    groups = ["Control", "Treatment"]
    batches = ["B1", "B2", "B3"]
    samples = []
    for g, group in enumerate(groups):
        for i in range(per_group):
            samples.append((f"{prefix}{group[0]}{i + 1}", group, batches[(i + g) % len(batches)]))
    base = [math.exp(rng.uniform(math.log(5), math.log(3000))) for _ in range(genes)]
    effect = [0.0] * genes
    batch_effect = [{b: 1.0 for b in batches} for _ in range(genes)]
    for gene in range(genes):
        draw = rng.random()
        if draw < 0.12:
            effect[gene] = rng.choice([-1, 1]) * rng.uniform(0.7, 3.0)
        if rng.random() < 0.25:
            batch_effect[gene]["B2"] = rng.uniform(0.6, 1.6)
            batch_effect[gene]["B3"] = rng.uniform(0.6, 1.6)
    size = {s[0]: rng.uniform(0.8, 1.25) for s in samples}
    header = ["gene_id", *[s[0] for s in samples]]
    lines = [",".join(header)]
    for gene in range(genes):
        row = [f"Gene{gene + 1:05d}"]
        for sample, group, batch in samples:
            mean = base[gene] * size[sample] * batch_effect[gene][batch] * (2 ** (effect[gene] if group == "Treatment" else 0.0))
            dispersion = 0.08
            lam = rng.gammavariate(1 / dispersion, mean * dispersion) if mean > 0 else 0.0
            row.append(str(poisson(rng, lam)))
        lines.append(",".join(row))
    (out / "counts.csv").write_text("\n".join(lines) + "\n", encoding="utf-8")
    # An optional multi-level grouping variable (listed last in the formula, so it is the PCA/boxplot colour variable).
    doses = [f"D{(index * 7 // 3) % dose_levels + 1}" for index in range(len(samples))] if dose_levels else []
    header_extra, rows_extra = (",dose", doses) if dose_levels else ("", [""] * len(samples))
    (out / "metadata.csv").write_text(f"sample_id,condition,batch{header_extra}\n" + "".join(f"{s},{g},{b}{(',' + d) if d else ''}\n" for (s, g, b), d in zip(samples, rows_extra)), encoding="utf-8")
    (out / "contrasts.csv").write_text("contrast_id,factor,numerator,denominator\nTreatment_vs_Control,condition,Treatment,Control\n", encoding="utf-8")
    (out / "project.yaml").write_text(
        'schema_version: "1.0"\nproject:\n  id: demo_two_group\n  pipeline: bulk_rnaseq\n  preset: L2\norganism:\n  species: Mus musculus\n'
        "input:\n  type: raw_counts\n  path: counts.csv\ndesign:\n  type: two_group\n  formula: \"~ batch + condition" + (" + dose" if dose_levels else "") + "\"\n"
        "metadata_file: metadata.csv\ncontrasts_file: contrasts.csv\nupstream:\n  engine: external\n  provider: external_provider\n  quantification_method: unknown\n"
        "reference:\n  source: igenomes\n  genome: null\nthresholds:\n  padj: 0.05\n  abs_log2fc: 1.0\n",
        encoding="utf-8",
    )


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("out", type=Path)
    parser.add_argument("--samples-per-group", type=int, default=6)
    parser.add_argument("--genes", type=int, default=4000)
    parser.add_argument("--seed", type=int, default=20261015)
    parser.add_argument("--sample-prefix", default="", help="prefix for every sample ID (long-label review)")
    parser.add_argument("--dose-levels", type=int, default=0, help="add a multi-level grouping variable with this many levels (colour variable)")
    args = parser.parse_args(argv)
    build(args.out, args.samples_per_group, args.genes, args.seed, args.sample_prefix, args.dose_levels)
    print(f"wrote {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
