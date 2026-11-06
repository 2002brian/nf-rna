"""Experimental offline STRING snapshot import, mapping and induced networks.

No HTTP requests, neighbor expansion or enrichment tests are performed.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import io
import json
import re
import shutil
import tempfile
from contextlib import contextmanager
from datetime import datetime, timezone
from pathlib import Path

from rnaseq.ppi import FIELDS, ID_TYPES, TAXA


SNAPSHOT_SCHEMA = "nf-rna.string-snapshot.v1"
MAPPING_FIELDS = ("query_id", "input_id_type", "string_id", "taxon_id")


def _identity(path: Path) -> dict:
    digest = hashlib.sha256()
    size = 0
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
            size += len(chunk)
    return {"sha256": digest.hexdigest(), "size_bytes": size}


@contextmanager
def _bundle(output: Path):
    if output.exists() or output.is_symlink():
        raise FileExistsError(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    stage = Path(tempfile.mkdtemp(prefix=".ppi-", dir=output.parent))
    try:
        yield stage
        output.mkdir()
        try:
            for path in sorted(stage.iterdir(), key=lambda p: (p.name in {"summary.json", "manifest.json"}, p.name)):
                shutil.move(str(path), output / path.name)
        except BaseException:
            shutil.rmtree(output)
            raise
    finally:
        shutil.rmtree(stage)


def _json(path: Path, data: dict) -> None:
    path.write_text(json.dumps(data, indent=2, ensure_ascii=False, allow_nan=False) + "\n", encoding="utf-8")


def _write_table(path: Path, fields, rows) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def _protein(value: str, taxon: int) -> None:
    if not re.fullmatch(rf"{taxon}\.[A-Za-z0-9_.-]+", value):
        raise ValueError(f"Invalid STRING protein ID or wrong taxon: {value!r}")


def _mappings(path: Path, taxon: int):
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames != list(MAPPING_FIELDS):
            raise ValueError(f"Mapping headers must be {MAPPING_FIELDS} in order.")
        for row in reader:
            if (None in row or any(value is None or not value or value != value.strip() for value in row.values())
                    or row["input_id_type"] not in ID_TYPES or row["taxon_id"] != str(taxon)
                    or row["query_id"] == "NA"):
                raise ValueError("Invalid mapping row, namespace or taxon.")
            _protein(row["string_id"], taxon)
            yield row


def _edges(path: Path, taxon: int):
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8") as handle:
        fields = handle.readline().split()
        if (len(fields) != len(set(fields)) or len(fields) < 3
                or fields[:2] != ["protein1", "protein2"] or "combined_score" not in fields):
            raise ValueError("Links require protein1, protein2 and combined_score headers.")
        yield fields, None
        for line in handle:
            values = line.split()
            if len(values) != len(fields):
                raise ValueError("Malformed STRING links row.")
            row = dict(zip(fields, values))
            _protein(row["protein1"], taxon)
            _protein(row["protein2"], taxon)
            for field in fields[2:]:
                if not row[field].isdigit() or not 0 <= int(row[field]) <= 1000:
                    raise ValueError("STRING download scores must be integer values in [0, 1000].")
            yield fields, row


def create_snapshot(mapping: Path, links: Path, output: Path, *, version: str, taxon_id: int,
                    network_type: str, source_urls: list[str]) -> Path:
    """Import explicit namespace mappings and STRING download-format links.

    Source URLs and version are operator declarations, not independently verified.
    A snapshot may contain a subset; absence of a mapping/edge is not absence in STRING.
    """
    if (taxon_id not in TAXA.values() or network_type not in {"functional", "physical"}
            or not re.fullmatch(r"\d+\.\d+", version)
            or not source_urls or any(not url.startswith("https://") for url in source_urls)):
        raise ValueError("Declare a supported taxon, network type, numeric STRING version and HTTPS source URLs.")
    with _bundle(Path(output)) as stage:
        mapping_name = "mapping.tsv"
        links_name = "links.txt.gz" if Path(links).suffix == ".gz" else "links.txt"
        shutil.copyfile(mapping, stage / mapping_name)
        shutil.copyfile(links, stage / links_name)
        # Validate copied bytes. Keep memory bounded even for full-species downloads.
        for _ in _mappings(stage / mapping_name, taxon_id):
            pass
        for _ in _edges(stage / links_name, taxon_id):
            pass
        _json(stage / "manifest.json", {
            "schema": SNAPSHOT_SCHEMA, "provider": "STRING", "version": version,
            "taxon_id": taxon_id, "network_type": network_type,
            "scope": "provided_tables", "source_urls": source_urls,
            "imported_at_utc": datetime.now(timezone.utc).isoformat(),
            "files": {"mapping": {"path": mapping_name, **_identity(stage / mapping_name)},
                      "links": {"path": links_name, **_identity(stage / links_name)}},
        })
    return Path(output)


def _verified_file(root: Path, record: dict) -> Path:
    name = record["path"]
    if not isinstance(name, str) or not name or Path(name).name != name:
        raise ValueError("Resource paths must be single filenames inside the snapshot.")
    path = root / name
    if path.is_symlink() or not path.is_file() or _identity(path) != {key: record[key] for key in ("sha256", "size_bytes")}:
        raise ValueError(f"Snapshot checksum/size mismatch or invalid file: {name}")
    return path


def extract_network(candidates: Path, snapshot: Path, output: Path, *, selection: str = "all",
                    required_score: int = 700, minimum_mapping_rate: float = 0.0) -> Path:
    """Map only unambiguous queries and extract their undirected induced graph.

    Mapping rate is the proportion of selected original genes with exactly one
    protein match. All alternatives are audited; ambiguous genes are excluded.
    Multiple unambiguous genes may point to the same protein, without combining LFCs.
    """
    if (selection not in {"all", "up", "down"} or type(required_score) is not int
            or not 0 <= required_score <= 1000 or not 0 <= minimum_mapping_rate <= 1):
        raise ValueError("Invalid selection, score or mapping-rate threshold.")
    candidates, snapshot, output = Path(candidates), Path(snapshot), Path(output)
    summary_bytes = (candidates / "summary.json").read_bytes()
    candidate_summary = json.loads(summary_bytes)
    if candidate_summary["schema"] != "nf-rna.ppi-candidates.v1":
        raise ValueError("Unsupported candidate schema.")
    data_path = _verified_file(candidates, {**candidate_summary["outputs"][f"{selection}.tsv"], "path": f"{selection}.tsv"})
    data_bytes = data_path.read_bytes()
    if hashlib.sha256(data_bytes).hexdigest() != candidate_summary["outputs"][f"{selection}.tsv"]["sha256"]:
        raise ValueError("Candidate changed while reading.")
    with io.StringIO(data_bytes.decode("utf-8")) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames != list(FIELDS):
            raise ValueError("Invalid candidate table headers.")
        genes = list(reader)
    id_type = candidate_summary["input_id_type"]
    if (id_type not in ID_TYPES or candidate_summary["counts"][selection] != len(genes)
            or len({row["gene_id"] for row in genes}) != len(genes)):
        raise ValueError("Invalid candidate namespace or duplicate genes.")
    if any(None in row or any(value is None for value in row.values()) or row["input_id_type"] != id_type
           or not row["gene_id"].strip() or not row["query_id"].strip() for row in genes):
        raise ValueError("Invalid candidate row.")
    manifest_bytes = (snapshot / "manifest.json").read_bytes()
    manifest = json.loads(manifest_bytes)
    taxon = candidate_summary["taxon_id"]
    if (manifest["schema"] != SNAPSHOT_SCHEMA or manifest["provider"] != "STRING"
            or manifest["taxon_id"] != taxon or TAXA.get(candidate_summary["organism"]) != taxon
            or manifest["network_type"] not in {"physical", "functional"}
            or not re.fullmatch(r"\d+\.\d+", manifest["version"])):
        raise ValueError("Snapshot schema, provider, version, network type or taxon mismatch.")
    paths = {role: _verified_file(snapshot, manifest["files"][role]) for role in ("mapping", "links")}
    queries = {row["query_id"] for row in genes}
    matches = {query: set() for query in queries}
    for row in _mappings(paths["mapping"], taxon):
        if row["input_id_type"] == id_type and row["query_id"] in queries:
            matches[row["query_id"]].add(row["string_id"])
    audit, nodes = [], {}
    counts = {"MAPPED": 0, "UNMAPPED": 0, "AMBIGUOUS_EXCLUDED": 0}
    for gene in sorted(genes, key=lambda row: row["gene_id"]):
        proteins = sorted(matches[gene["query_id"]])
        status = "MAPPED" if len(proteins) == 1 else "AMBIGUOUS_EXCLUDED" if proteins else "UNMAPPED"
        counts[status] += 1
        for protein in proteins or [""]:
            audit.append({**gene, "string_id": protein, "mapping_status": status})
        if status == "MAPPED":
            nodes.setdefault(proteins[0], []).append(gene["gene_id"])
    rate = counts["MAPPED"] / len(genes) if genes else None
    if rate is not None and rate < minimum_mapping_rate:
        raise ValueError(f"Unambiguous gene mapping rate {rate:.3f} is below {minimum_mapping_rate:.3f}.")
    edges = {}
    edge_fields = ["protein1", "protein2", "combined_score"]
    self_loops = 0
    for fields, row in _edges(paths["links"], taxon):
        edge_fields = fields
        if row is None:
            continue
        a, b = sorted((row["protein1"], row["protein2"]))
        if a not in nodes or b not in nodes or int(row["combined_score"]) < required_score:
            continue
        if a == b:
            self_loops += 1
            continue
        normalized = {**row, "protein1": a, "protein2": b}
        if (a, b) in edges and edges[a, b] != normalized:
            raise ValueError("Conflicting duplicate edge scores in snapshot.")
        edges[a, b] = normalized
    # Detect resource changes during streaming before publishing any results.
    for role in paths:
        _verified_file(snapshot, manifest["files"][role])
    degree = {protein: 0 for protein in nodes}
    for a, b in edges:
        degree[a] += 1
        degree[b] += 1
    status = "NO_CANDIDATES" if not genes else "NO_MAPPED_PROTEINS" if not nodes else "NO_EDGES" if not edges else "COMPLETE"
    with _bundle(output) as stage:
        _write_table(stage / "mapping.tsv", (*FIELDS, "string_id", "mapping_status"), audit)
        _write_table(stage / "nodes.tsv", ("string_id", "gene_ids", "degree"),
                     ({"string_id": protein, "gene_ids": json.dumps(sorted(nodes[protein])), "degree": degree[protein]}
                      for protein in sorted(nodes)))
        _write_table(stage / "edges.tsv", edge_fields, (edges[key] for key in sorted(edges)))
        (stage / "candidate_summary.json").write_bytes(summary_bytes)
        (stage / "snapshot_manifest.json").write_bytes(manifest_bytes)
        (stage / "selected_candidates.tsv").write_bytes(data_bytes)
        _json(stage / "summary.json", {
            "schema": "nf-rna.ppi-network.v1", "experimental": True, "status": status,
            "contrast_id": candidate_summary["contrast_id"], "selection": selection,
            "numerator": candidate_summary["numerator"], "denominator": candidate_summary["denominator"],
            "taxon_id": taxon, "string_version": manifest["version"], "network_type": manifest["network_type"],
            "required_score": required_score, "minimum_mapping_rate": minimum_mapping_rate,
            "mapping_counts": counts, "mapping_rate": rate, "node_count": len(nodes), "edge_count": len(edges),
            "isolated_node_count": sum(value == 0 for value in degree.values()), "self_loop_rows_excluded": self_loops,
            "mapping_policy": "Exclude ambiguous queries; retain all mapping alternatives in audit table.",
            "limitations": ["Only the provided snapshot tables are searched; they may be incomplete.",
                            "No neighbor expansion or statistical enrichment test.",
                            "Functional associations do not establish direct physical binding.",
                            "Manifest metadata and source-run authenticity are not independently verified."],
            "outputs": {path.name: _identity(path) for path in sorted(stage.iterdir())},
        })
    return output


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    snapshot = commands.add_parser("snapshot")
    for name in ("mapping", "links", "outdir"):
        snapshot.add_argument(f"--{name}", type=Path, required=True)
    snapshot.add_argument("--version", required=True)
    snapshot.add_argument("--taxon-id", type=int, choices=TAXA.values(), required=True)
    snapshot.add_argument("--network-type", choices=("functional", "physical"), required=True)
    snapshot.add_argument("--source-url", action="append", required=True)
    network = commands.add_parser("network")
    for name in ("candidates", "snapshot", "outdir"):
        network.add_argument(f"--{name}", type=Path, required=True)
    network.add_argument("--selection", choices=("all", "up", "down"), default="all")
    network.add_argument("--required-score", type=int, default=700)
    network.add_argument("--minimum-mapping-rate", type=float, default=0.0)
    args = parser.parse_args()
    if args.command == "snapshot":
        create_snapshot(args.mapping, args.links, args.outdir, version=args.version, taxon_id=args.taxon_id,
                        network_type=args.network_type, source_urls=args.source_url)
    else:
        extract_network(args.candidates, args.snapshot, args.outdir, selection=args.selection,
                        required_score=args.required_score, minimum_mapping_rate=args.minimum_mapping_rate)


if __name__ == "__main__":
    main()
