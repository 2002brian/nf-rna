"""Real-R checks that Salmon/tximport keeps every quant.sf with its own sample.

Regression for 4T1 run grcm39-v130-salmon-r5: L1 and L2 flattened the
``quant_sf`` mapping and renamed the files positionally with the metadata
sample order.  The control plane writes that mapping with sorted keys, so with
non-alphabetical metadata every sample was analysed with another sample's
counts (a control library was modelled as treated and vice versa).
"""

from __future__ import annotations

import csv
import json
import random
import subprocess
from pathlib import Path

import pytest

from conftest import require_r_packages
from test_r_design_metadata import R_SCRIPTS, _frozen_configs, _run_r


R_PACKAGES = ("jsonlite", "yaml", "DESeq2", "ggplot2", "pheatmap", "tximport")
# Metadata order is deliberately neither alphabetical nor grouped by condition.
METADATA_ORDER = ("T2", "C1", "T3", "C2", "T1", "C3")
CONDITION = {"C1": "Control", "C2": "Control", "C3": "Control", "T1": "Treatment", "T2": "Treatment", "T3": "Treatment"}
# The quant_sf mapping order differs from both metadata and alphabetical order.
MAPPING_ORDER = ("C3", "T1", "C1", "T3", "T2", "C2")
GENES = 120


def _metadata() -> str:
    return "sample_id,condition\n" + "".join(f"{sample},{CONDITION[sample]}\n" for sample in METADATA_ORDER)


def _write_quant(path: Path, reads: dict[str, int]) -> None:
    path.parent.mkdir(parents=True)
    lines = ["Name\tLength\tEffectiveLength\tTPM\tNumReads"]
    total = sum(reads.values()) or 1
    lines += [f"{tx}\t1000\t800.000\t{1e6 * value / total:.4f}\t{value}.000" for tx, value in reads.items()]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _salmon_inputs(root: Path) -> tuple[dict[str, str], Path, dict[str, dict[str, int]]]:
    """Distinguishable quant.sf per sample; genes 0-19 are 8x up in Treatment."""

    rng = random.Random(20260925)
    base = [rng.randint(40, 2000) for _ in range(GENES)]
    expected: dict[str, dict[str, int]] = {}
    quant_sf: dict[str, str] = {}
    for index, sample in enumerate(sorted(CONDITION)):
        effect = 8 if CONDITION[sample] == "Treatment" else 1
        # A per-sample offset makes every sample's column unique.
        reads = {
            f"Tx{gene:03d}": base[gene] * (effect if gene < 20 else 1) + rng.randint(0, 60) + 7 * index
            for gene in range(GENES)
        }
        path = root / sample / "quant.sf"
        _write_quant(path, reads)
        expected[sample] = {f"Gene{gene:03d}": reads[f"Tx{gene:03d}"] for gene in range(GENES)}
    for sample in MAPPING_ORDER:
        quant_sf[sample] = str(root / sample / "quant.sf")
    tx2gene = root / "tx2gene.tsv"
    tx2gene.write_text("transcript_id\tgene_id\n" + "".join(f"Tx{gene:03d}\tGene{gene:03d}\n" for gene in range(GENES)), encoding="utf-8")
    return quant_sf, tx2gene, expected


def _salmon_configs(project_factory, tmp_path: Path):
    work, l1, l2 = _frozen_configs(
        project_factory, tmp_path, formula="~ condition", variables={"condition": "categorical"}, metadata=_metadata(),
    )
    quant_sf, tx2gene, expected = _salmon_inputs(tmp_path / "salmon")
    for config in (l1, l2):
        assert config["samples"] == list(METADATA_ORDER)
        config.pop("counts")
        config.update(source_type="salmon_tximport", quant_sf=dict(quant_sf), tx2gene=str(tx2gene), tx2gene_mapping="test")
    assert list(l1["quant_sf"]) != l1["samples"] and list(l1["quant_sf"]) != sorted(l1["samples"])
    return work, l1, l2, expected


def _independent_l2(work: Path, l1: dict, contrast: tuple[str, str, str]) -> Path:
    """DESeq2 on the same quant.sf files, with each file named explicitly."""

    output = work / "reference.tsv"
    files = ", ".join(f"{json.dumps(sample)} = {json.dumps(l1['quant_sf'][sample])}" for sample in METADATA_ORDER)
    conditions = ", ".join(json.dumps(CONDITION[sample]) for sample in METADATA_ORDER)
    script = f"""
suppressPackageStartupMessages({{ library(DESeq2); library(tximport) }})
files <- c({files})
tx2gene <- read.delim({json.dumps(l1['tx2gene'])}, check.names = FALSE, stringsAsFactors = FALSE)[, 1:2]
txi <- tximport(files, type = "salmon", tx2gene = tx2gene)
metadata <- data.frame(condition = factor(c({conditions})), row.names = names(files))
dds <- DESeq(DESeqDataSetFromTximport(txi, colData = metadata, design = ~ condition), quiet = TRUE)
result <- results(dds, contrast = c({", ".join(json.dumps(item) for item in contrast)}), independentFiltering = TRUE)
write.table(data.frame(gene_id = rownames(result), as.data.frame(result))[, c("gene_id", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")],
  {json.dumps(str(output))}, sep = "\\t", quote = FALSE, row.names = FALSE, na = "NA")
"""
    result = subprocess.run(["Rscript", "-e", script], capture_output=True, text=True, check=False)
    assert result.returncode == 0, result.stderr
    return output


def _read_tsv(path: Path) -> dict[str, dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return {row["gene_id"]: row for row in csv.DictReader(handle, delimiter="\t")}


def _number(value: str) -> float | None:
    return None if value == "NA" else float(value)


def test_l1_and_l2_keep_each_quant_sf_with_its_own_sample(project_factory, tmp_path):
    require_r_packages(*R_PACKAGES)
    work, l1, l2, expected = _salmon_configs(project_factory, tmp_path)

    result = _run_r("l1_analysis", l1, work)
    assert result.returncode == 0, result.stderr
    with (work / "l1" / "source_counts.csv").open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle))
    assert list(rows[0]) == ["gene_id", *METADATA_ORDER]
    for sample in METADATA_ORDER:
        observed = {row["gene_id"]: float(row[sample]) for row in rows}
        assert observed == {gene: float(value) for gene, value in expected[sample].items()}, sample
    # Normalized counts and VST are per-sample transforms of the same columns.
    normalized = _read_tsv(work / "l1" / "normalized_counts.tsv")
    for sample in METADATA_ORDER:
        ratios = [float(normalized[gene][sample]) / expected[sample][gene] for gene in expected[sample]]
        assert max(ratios) == pytest.approx(min(ratios), rel=1e-6), sample

    result = _run_r("l2_analysis", l2, work)
    assert result.returncode == 0, result.stderr
    contrast = l2["contrasts"][0]
    triple = (contrast["factor"], contrast["numerator"], contrast["denominator"])
    assert triple == ("condition", "Treatment", "Control")
    observed = _read_tsv(work / "l2" / "contrasts" / contrast["contrast_id"] / "all_genes.tsv")
    reference = _read_tsv(_independent_l2(work, l1, triple))
    assert observed.keys() == reference.keys()
    for gene, row in reference.items():
        for column in ("baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj"):
            assert _number(observed[gene][column]) == pytest.approx(_number(row[column]), rel=1e-9, abs=1e-12), (gene, column)
    # The designed 8x Treatment effect is recovered with the right sign.
    for gene in (f"Gene{index:03d}" for index in range(20)):
        assert float(observed[gene]["log2FoldChange"]) == pytest.approx(3, abs=0.25)


def _write_config_text(work: Path, script: str, text: str) -> subprocess.CompletedProcess[str]:
    path = work / f"{script}.raw.json"
    path.write_text(text, encoding="utf-8")
    return subprocess.run(
        ["Rscript", str(R_SCRIPTS / f"{script}.R"), "--config", str(path)], cwd=work, capture_output=True, text=True, check=False,
    )


def _quant_sf_json(entries: list[tuple[str, str]]) -> str:
    # Written by hand so that duplicate keys reach R exactly as given.
    return "{" + ", ".join(f"{json.dumps(key)}: {json.dumps(value)}" for key, value in entries) + "}"


def _invalid_mapping_cases(quant_sf: dict[str, str]) -> dict[str, tuple[str, str]]:
    entries = list(quant_sf.items())
    extra_path = entries[0][1]
    return {
        "missing_sample": (_quant_sf_json(entries[:-1]), f"quant_sf is missing metadata samples: {entries[-1][0]}"),
        "extra_sample": (_quant_sf_json([*entries, ("X9", extra_path)]), "quant_sf lists samples absent from metadata: X9"),
        "duplicate_sample": (_quant_sf_json([*entries[:-1], (entries[0][0], entries[-1][1])]), f"quant_sf lists a sample_id more than once: {entries[0][0]}"),
        "empty_sample_id": (_quant_sf_json([*entries, ("", extra_path)]), "quant_sf contains an entry without a sample_id"),
        "unnamed_mapping": (json.dumps([value for _key, value in entries]), "quant_sf must map each sample_id to its quant.sf path"),
        "empty_path": (_quant_sf_json([*entries[:-1], (entries[-1][0], "")]), f"quant_sf path is missing or invalid for sample: {entries[-1][0]}"),
        "null_path": (_quant_sf_json(entries[:-1])[:-1] + f", {json.dumps(entries[-1][0])}: null}}", f"quant_sf path is missing or invalid for sample: {entries[-1][0]}"),
    }


@pytest.mark.parametrize("script", ["l1_analysis", "l2_analysis"])
@pytest.mark.parametrize(
    "case", ["missing_sample", "extra_sample", "duplicate_sample", "empty_sample_id", "unnamed_mapping", "empty_path", "null_path"],
)
def test_invalid_quant_sf_identity_fails_before_tximport(project_factory, tmp_path, script, case):
    require_r_packages(*R_PACKAGES)
    work, l1, l2, _expected = _salmon_configs(project_factory, tmp_path)
    config = l1 if script == "l1_analysis" else l2
    mapping, message = _invalid_mapping_cases(config["quant_sf"])[case]
    text = json.dumps({key: value for key, value in config.items() if key != "quant_sf"})
    result = _write_config_text(work, script, text[:-1] + f', "quant_sf": {mapping}}}')
    assert result.returncode != 0
    assert message in result.stderr, result.stderr
    assert "reading in files" not in result.stdout + result.stderr
    assert not (work / script[:2] / "source_counts.csv").exists()
