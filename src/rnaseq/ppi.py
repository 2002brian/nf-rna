"""Experimental, offline PPI candidate preparation; not a production run entry point.

The L2 DEG definition remains authoritative. This module checks that the
delivered significant set agrees with its recorded definition, retains original
gene IDs, and prepares query IDs without resolving them to proteins.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import math
import re
import shutil
import tempfile
from pathlib import Path


TAXA = {"Homo sapiens": 9606, "Mus musculus": 10090}
ID_TYPES = ("ENSEMBL", "ENTREZID", "SYMBOL")
SOURCE_FILES = ("all_genes.tsv", "significant.tsv", "contrast_metadata.json")
FIELDS = ("gene_id", "query_id", "input_id_type", "symbol", "log2FoldChange", "padj")


def _number(value: str) -> float | None:
    if value in {"", "NA"}:
        return None
    number = float(value)
    if not math.isfinite(number):
        raise ValueError("Non-finite statistics are not supported.")
    return number


def _table(data: bytes) -> dict[str, dict[str, str]]:
    reader = csv.DictReader(io.StringIO(data.decode("utf-8")), delimiter="\t")
    fields = reader.fieldnames or []
    if len(fields) != len(set(fields)) or not {"gene_id", "log2FoldChange", "padj"} <= set(fields):
        raise ValueError("DEG table requires unique headers and gene_id, log2FoldChange, padj.")
    records = {}
    for row in reader:
        if None in row or any(value is None for value in row.values()):
            raise ValueError("Malformed DEG row.")
        gene = row["gene_id"]
        if not gene.strip() or gene == "NA" or gene != gene.strip() or gene in records:
            raise ValueError("DEG gene IDs must be nonblank, unique and without surrounding whitespace.")
        padj = _number(row["padj"])
        _number(row["log2FoldChange"])
        if padj is not None and not 0 <= padj <= 1:
            raise ValueError("padj must be in [0, 1].")
        records[gene] = row
    return records


def prepare_ppi_inputs(contrast_dir: Path, output_dir: Path, *, organism: str, input_id_type: str) -> Path:
    """Create a new candidate bundle from an existing L2 contrast directory.

    Reads/hashes the same byte snapshots and validates before writing outputs.
    ID type/species must be declared explicitly; no symbol or species inference.
    A valid empty significant set produces EMPTY rather than an execution error.
    """
    if organism not in TAXA or input_id_type not in ID_TYPES:
        raise ValueError("Unsupported organism or input ID type.")
    contrast_dir, output_dir = Path(contrast_dir), Path(output_dir)
    if output_dir.exists() or output_dir.is_symlink():
        raise FileExistsError(output_dir)
    sources = {name: (contrast_dir / name).read_bytes() for name in SOURCE_FILES}
    metadata = json.loads(sources["contrast_metadata.json"])
    contrast_id = metadata["contrast_id"]
    if not isinstance(contrast_id, str) or not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9._-]*", contrast_id):
        raise ValueError("Invalid contrast ID.")
    for direction in ("numerator", "denominator"):
        if not isinstance(metadata[direction]["level"], str) or not metadata[direction]["level"].strip():
            raise ValueError("Contrast direction must contain nonblank level labels.")
    definition = metadata["deg_definition"]
    alpha, effect = definition["padj_below"], definition["abs_log2fc_at_least"]
    if (type(alpha) not in (int, float) or type(effect) not in (int, float)
            or not 0 < alpha < 1 or not math.isfinite(effect) or effect < 0):
        raise ValueError("Invalid recorded DEG thresholds.")
    all_genes = _table(sources["all_genes.tsv"])
    significant = _table(sources["significant.tsv"])
    if not all_genes:
        raise ValueError("all_genes.tsv must contain genes.")
    expected = {
        gene for gene, row in all_genes.items()
        if (p := _number(row["padj"])) is not None and p < alpha
        and (lfc := _number(row["log2FoldChange"])) is not None and abs(lfc) >= effect
    }
    if set(significant) != expected:
        raise ValueError("Significant gene set disagrees with the recorded DEG definition.")
    for gene, row in significant.items():
        for field in ("log2FoldChange", "padj"):
            if _number(row[field]) != _number(all_genes[gene][field]):
                raise ValueError("Significant statistics disagree with all_genes.tsv.")
    candidates = []
    query_sources: dict[str, list[str]] = {}
    for gene in sorted(significant):
        query = re.sub(r"\.\d+$", "", gene) if input_id_type == "ENSEMBL" else gene
        query_sources.setdefault(query, []).append(gene)
        candidates.append({**{key: all_genes[gene].get(key, "NA") for key in FIELDS},
                           "query_id": query, "input_id_type": input_id_type})
    sets = {
        "all": candidates,
        "up": [row for row in candidates if float(row["log2FoldChange"]) > 0],
        "down": [row for row in candidates if float(row["log2FoldChange"]) < 0],
    }
    summary = {
        "schema": "nf-rna.ppi-candidates.v1", "experimental": True,
        "status": "READY_FOR_MAPPING" if candidates else "EMPTY",
        "mapping_status": "NOT_RUN", "organism": organism, "taxon_id": TAXA[organism],
        "input_id_type": input_id_type, "contrast_id": contrast_id,
        "numerator": metadata["numerator"], "denominator": metadata["denominator"],
        "deg_definition": definition, "selection": "existing L2 significant set, validated against recorded thresholds",
        "counts": {name: len(rows) for name, rows in sets.items()},
        "zero_lfc_candidates": sum(float(row["log2FoldChange"]) == 0 for row in candidates),
        "query_id_collisions": {key: value for key, value in sorted(query_sources.items()) if len(value) > 1},
        "sources": {name: {"sha256": hashlib.sha256(data).hexdigest(), "size_bytes": len(data)}
                    for name, data in sources.items()},
        "limitations": ["No protein mapping, network retrieval, PPI inference or enrichment has run.",
                        "Source bytes are recorded; source-run authenticity is not verified by this standalone helper."],
    }
    output_dir.parent.mkdir(parents=True, exist_ok=True)
    staging = Path(tempfile.mkdtemp(prefix=".ppi-candidates-", dir=output_dir.parent))
    try:
        for name, rows in sets.items():
            with (staging / f"{name}.tsv").open("w", encoding="utf-8", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=FIELDS, delimiter="\t", lineterminator="\n")
                writer.writeheader()
                writer.writerows(rows)
        for name, data in sources.items():
            (staging / f"source_{name}").write_bytes(data)
        summary["outputs"] = {path.name: {"sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
                                          "size_bytes": path.stat().st_size}
                              for path in sorted(staging.iterdir())}
        (staging / "summary.json").write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
        # Publish without replacing a pre-existing destination, including an empty one.
        output_dir.mkdir()
        try:
            # The summary is the completion marker and is always published last.
            for path in sorted(staging.iterdir(), key=lambda item: (item.name == "summary.json", item.name)):
                shutil.move(str(path), output_dir / path.name)
        except BaseException:
            shutil.rmtree(output_dir)
            raise
    finally:
        shutil.rmtree(staging)
    return output_dir


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--contrast-dir", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    parser.add_argument("--organism", choices=TAXA, required=True)
    parser.add_argument("--input-id-type", choices=ID_TYPES, required=True)
    args = parser.parse_args()
    prepare_ppi_inputs(args.contrast_dir, args.outdir, organism=args.organism, input_id_type=args.input_id_type)


if __name__ == "__main__":
    main()
