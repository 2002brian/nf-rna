"""Salmon/tximport input regression gate.

Synthetic quant.sf files whose transcripts have sample-specific effective lengths
(an isoform-length shift in Treatment), so the average-transcript-length offsets
of DESeqDataSetFromTximport change the statistics.  nf-rna's L1/L2 R backends
(driven through the frozen run contract) must:

* aggregate NumReads per gene exactly as tximport does (checked against sums
  computed here from the quant.sf text);
* equal ``oracles/tximport_reference.R`` with offsets (rel 1e-8), and
* differ from the no-offset analysis of the same counts, proving the offsets are used.
"""

from __future__ import annotations

import csv
import math
import random
from pathlib import Path

from _science import by_gene, close, negative_binomial, read_tsv, require_downstream_r_runtime, run_reference
from test_r_design_metadata import _frozen_configs, _run_r

SAMPLES = ("C1", "C2", "C3", "T1", "T2", "T3")
GENES = 150


def _quant_inputs(root: Path) -> tuple[dict[str, str], Path, dict[str, dict[str, float]]]:
    rng = random.Random(4242)
    transcripts = {f"Gene{g:03d}": [f"Tx{g:03d}_{i}" for i in range(1 + g % 3)] for g in range(GENES)}
    base = {tx: math.exp(rng.uniform(math.log(40), math.log(2000))) for txs in transcripts.values() for tx in txs}
    length = {tx: rng.uniform(600, 4000) for tx in base}
    gene_sums: dict[str, dict[str, float]] = {}
    quant: dict[str, str] = {}
    for index, sample in enumerate(SAMPLES):
        treated = sample.startswith("T")
        records = []
        sums: dict[str, float] = {}
        for gene, txs in transcripts.items():
            for position, tx in enumerate(txs):
                effect = 3.0 if treated and int(gene[4:]) < 20 else 1.0
                # Isoform-length shift: in Treatment the first isoform of every third gene is 40% longer,
                # so it yields more reads at unchanged abundance.  Only the length offsets can remove that.
                effective = length[tx] * (1.4 if treated and position == 0 and int(gene[4:]) % 3 == 0 else 1.0) * (0.97 + 0.01 * index)
                reads = float(negative_binomial(rng, base[tx] * effect * effective / 2000.0, 0.05))
                records.append((tx, effective, reads))
                sums[gene] = sums.get(gene, 0.0) + reads
        # tximport weights transcript lengths by abundance, so TPM must be consistent with NumReads / EffectiveLength.
        rate = sum(reads / effective for _tx, effective, reads in records)
        lines = ["Name\tLength\tEffectiveLength\tTPM\tNumReads"]
        lines += [f"{tx}\t{length[tx]:.0f}\t{effective:.3f}\t{1e6 * reads / effective / rate:.6f}\t{reads:.3f}" for tx, effective, reads in records]
        path = root / sample / "quant.sf"
        path.parent.mkdir(parents=True)
        path.write_text("\n".join(lines) + "\n", encoding="utf-8")
        quant[sample] = str(path)
        gene_sums[sample] = sums
    tx2gene = root / "tx2gene.tsv"
    tx2gene.write_text("transcript_id\tgene_id\n" + "".join(f"{tx}\t{gene}\n" for gene, txs in transcripts.items() for tx in txs), encoding="utf-8")
    return quant, tx2gene, gene_sums


def test_tximport_offsets_reach_deseq2_and_match_an_independent_reference(project_factory, tmp_path):
    require_downstream_r_runtime("tximport", "DESeq2", "jsonlite", "yaml", "ggplot2", "pheatmap")
    metadata = "sample_id,condition\n" + "".join(f"{s},{'Treatment' if s.startswith('T') else 'Control'}\n" for s in SAMPLES)
    work, l1, l2 = _frozen_configs(project_factory, tmp_path, formula="~ condition", variables={"condition": "categorical"}, metadata=metadata)
    quant, tx2gene, gene_sums = _quant_inputs(tmp_path / "salmon")
    for config in (l1, l2):
        config.pop("counts")
        config.update(source_type="salmon_tximport", quant_sf=dict(reversed(list(quant.items()))), tx2gene=str(tx2gene), tx2gene_mapping="test")
    for script, config in (("l1_analysis", l1), ("l2_analysis", l2)):
        result = _run_r(script, config, work)
        assert result.returncode == 0, result.stderr

    with (work / "l1" / "source_counts.csv").open(encoding="utf-8", newline="") as handle:
        source = {row["gene_id"]: row for row in csv.DictReader(handle)}
    for sample in SAMPLES:
        for gene, expected in gene_sums[sample].items():
            assert math.isclose(float(source[gene][sample]), expected, rel_tol=1e-9), (sample, gene)

    files = tmp_path / "files.tsv"
    files.write_text("".join(f"{s}\t{quant[s]}\t{'Treatment' if s.startswith('T') else 'Control'}\n" for s in SAMPLES), encoding="utf-8")
    run_reference("tximport_reference.R", files, tx2gene, tmp_path / "txref")
    observed = by_gene(read_tsv(work / "l2" / "contrasts" / l2["contrasts"][0]["contrast_id"] / "all_genes.tsv"))
    offset = by_gene(read_tsv(tmp_path / "txref_offset.tsv"))
    no_offset = by_gene(read_tsv(tmp_path / "txref_no_offset.tsv"))
    assert observed.keys() == offset.keys()
    number = lambda value: None if value == "NA" else float(value)  # noqa: E731
    for gene, row in offset.items():
        for column in ("baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj"):
            assert close(number(observed[gene][column]), number(row[column]), rel=1e-8), (gene, column)
    differing = sum(not close(number(observed[g]["log2FoldChange"]), number(no_offset[g]["log2FoldChange"]), rel=1e-6) for g in observed)
    assert differing > 0.5 * len(observed), "length offsets must change the estimates"
