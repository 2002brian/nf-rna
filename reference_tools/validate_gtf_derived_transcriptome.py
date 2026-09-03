"""Validate a one-time GTF-derived transcriptome before Salmon indexing.

This utility is intentionally outside the production control-plane package.  It
creates an auditable JSON artifact for a reference-library build without
changing the reference manifest or any project configuration.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
from collections import defaultdict
from pathlib import Path

import yaml


ATTRIBUTE_PATTERN = re.compile(r'(\S+) "(.*?)(?<!\\)";')
VALID_DNA = frozenset("ACGTUNRYKMSWBDHV")


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def parse_gtf(path: Path) -> tuple[set[str], dict[str, set[str]]]:
    transcript_ids: set[str] = set()
    transcript_genes: dict[str, set[str]] = defaultdict(set)
    with path.open("r", encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9:
                raise ValueError(f"Malformed GTF line {line_number}: expected nine tab-separated columns.")
            attributes = dict(ATTRIBUTE_PATTERN.findall(fields[8]))
            transcript_id = attributes.get("transcript_id")
            gene_id = attributes.get("gene_id")
            if transcript_id:
                transcript_ids.add(transcript_id)
                if gene_id:
                    transcript_genes[transcript_id].add(gene_id)
    return transcript_ids, transcript_genes


def parse_fasta(path: Path) -> tuple[set[str], int, int, list[str], list[str]]:
    identifiers: set[str] = set()
    duplicate_ids: list[str] = []
    invalid_records: list[str] = []
    sequence_length = 0
    current_id: str | None = None
    current_length = 0
    with path.open("r", encoding="utf-8") as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            line = raw_line.rstrip("\n\r")
            if not line:
                continue
            if line.startswith(">"):
                if current_id is not None and current_length == 0:
                    invalid_records.append(f"{current_id}: empty sequence")
                current_id = line[1:].split(maxsplit=1)[0]
                current_length = 0
                if not current_id:
                    invalid_records.append(f"line {line_number}: blank FASTA identifier")
                    continue
                if current_id in identifiers:
                    duplicate_ids.append(current_id)
                identifiers.add(current_id)
                continue
            if current_id is None:
                invalid_records.append(f"line {line_number}: sequence before FASTA header")
                continue
            sequence = line.upper()
            unsupported = set(sequence) - VALID_DNA
            if unsupported:
                invalid_records.append(
                    f"{current_id}: unsupported sequence symbols {''.join(sorted(unsupported))}"
                )
            current_length += len(sequence)
            sequence_length += len(sequence)
    if current_id is None:
        invalid_records.append("no FASTA records")
    elif current_length == 0:
        invalid_records.append(f"{current_id}: empty sequence")
    return identifiers, len(identifiers), sequence_length, duplicate_ids, invalid_records


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference-root", type=Path, required=True)
    parser.add_argument("--transcriptome", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    root = args.reference_root.resolve()
    transcriptome = args.transcriptome.resolve()
    output = args.output.resolve()
    manifest_path = root / "reference_manifest.yaml"
    manifest = yaml.safe_load(manifest_path.read_text(encoding="utf-8"))
    files = manifest["files"]
    genome = (root / files["genome_fasta"]["path"]).resolve()
    gtf = (root / files["annotation_gtf"]["path"]).resolve()

    gtf_ids, transcript_genes = parse_gtf(gtf)
    fasta_ids, fasta_count, total_length, duplicates, invalid_records = parse_fasta(transcriptome)
    fasta_only = sorted(fasta_ids - gtf_ids)
    gtf_only = sorted(gtf_ids - fasta_ids)
    one_gene = sorted(identifier for identifier in fasta_ids if len(transcript_genes.get(identifier, set())) == 1)
    zero_gene = sorted(identifier for identifier in fasta_ids if len(transcript_genes.get(identifier, set())) == 0)
    multi_gene = sorted(identifier for identifier in fasta_ids if len(transcript_genes.get(identifier, set())) > 1)
    passed = not (fasta_only or zero_gene or multi_gene or duplicates or invalid_records)
    payload = {
        "reference": {
            key: manifest["reference"][key]
            for key in ("species", "provider", "release", "assembly", "assembly_patch")
        },
        "inputs": {
            "genome": {"path": str(genome.relative_to(root)), "sha256": sha256_file(genome)},
            "gtf": {"path": str(gtf.relative_to(root)), "sha256": sha256_file(gtf)},
        },
        "generated_transcriptome": {
            "path": str(transcriptome.relative_to(root)),
            "sha256": sha256_file(transcriptome),
            "transcript_count": fasta_count,
            "total_sequence_length": total_length,
            "file_size_bytes": transcriptome.stat().st_size,
        },
        "identifier_compatibility": {
            "fasta_unique_ids": fasta_count,
            "gtf_unique_transcript_ids": len(gtf_ids),
            "intersection": len(fasta_ids & gtf_ids),
            "fasta_only": len(fasta_only),
            "gtf_only": len(gtf_only),
            "fasta_mapping_rate": len(fasta_ids & gtf_ids) / fasta_count if fasta_count else 0.0,
            "gtf_mapping_rate": len(fasta_ids & gtf_ids) / len(gtf_ids) if gtf_ids else 0.0,
            "fasta_only_examples": fasta_only[:20],
            "gtf_only_examples": gtf_only[:20],
        },
        "tx2gene": {
            "one_gene_mappings": len(one_gene),
            "zero_gene_mappings": len(zero_gene),
            "multi_gene_mappings": len(multi_gene),
            "zero_gene_examples": zero_gene[:20],
            "multi_gene_examples": multi_gene[:20],
        },
        "fasta_validation": {
            "duplicate_transcript_ids": len(duplicates),
            "duplicate_examples": sorted(set(duplicates))[:20],
            "invalid_records": invalid_records[:20],
        },
        "status": "PASS" if passed else "BLOCKED",
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.with_name(f".{output.name}.{os.getpid()}.tmp")
    temporary.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    temporary.replace(output)
    print(json.dumps(payload, sort_keys=True))
    return 0 if passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
