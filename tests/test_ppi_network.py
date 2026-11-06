"""Scientific contract tests using small synthetic, explicitly declared snapshots."""

import csv
import gzip
import hashlib
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

from rnaseq.ppi import prepare_ppi_inputs
from rnaseq.ppi_network import create_snapshot, extract_network


class PpiNetworkTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        contrast = self.root / "contrast"
        contrast.mkdir()
        table = "gene_id\tsymbol\tlog2FoldChange\tpadj\n" + "".join(
            f"{gene}\t{gene}\t{2 if gene != 'B' else -2}\t0.01\n" for gene in "ABCDE")
        for name in ("all_genes.tsv", "significant.tsv"):
            (contrast / name).write_text(table)
        (contrast / "contrast_metadata.json").write_text(json.dumps({
            "contrast_id": "B_vs_A", "numerator": {"level": "B"}, "denominator": {"level": "A"},
            "deg_definition": {"padj_below": 0.05, "abs_log2fc_at_least": 1},
        }))
        self.candidates = self.root / "candidates"
        prepare_ppi_inputs(contrast, self.candidates, organism="Homo sapiens", input_id_type="SYMBOL")
        self.mapping = self.root / "mapping.tsv"
        self.mapping.write_text("query_id\tinput_id_type\tstring_id\ttaxon_id\n" + "".join(
            f"{gene}\tSYMBOL\t9606.P{protein}\t9606\n" for gene, protein in
            [("A", 1), ("B", 2), ("C", 3), ("D", 4), ("D", 5), ("A", 1)]))
        self.links = self.root / "links.txt"
        self.links.write_text(
            "protein1 protein2 experimental combined_score\n"
            "9606.P1 9606.P2 600 700\n9606.P2 9606.P1 600 700\n"
            "9606.P2 9606.P3 500 699\n9606.P1 9606.P99 900 1000\n"
            "9606.P3 9606.P3 900 900\n")
        self.snapshot = self.root / "snapshot"
        self.output = self.root / "network"

    def snapshot_import(self, links=None, network_type="functional"):
        create_snapshot(self.mapping, links or self.links, self.snapshot, version="12.5", taxon_id=9606,
                        network_type=network_type, source_urls=["https://example.org/synthetic-test-fixture"])

    def extract(self, **kwargs):
        return extract_network(self.candidates, self.snapshot, self.output, **kwargs)

    def table(self, name):
        with (self.output / name).open() as handle:
            return list(csv.DictReader(handle, delimiter="\t"))

    def summary(self):
        return json.loads((self.output / "summary.json").read_text())

    def test_induced_network_ambiguity_evidence_isolates_and_hashes(self):
        self.snapshot_import()
        self.extract()
        summary = self.summary()
        self.assertEqual(summary["mapping_counts"], {"MAPPED": 3, "UNMAPPED": 1, "AMBIGUOUS_EXCLUDED": 1})
        self.assertEqual(summary["mapping_rate"], 0.6)
        self.assertEqual((summary["node_count"], summary["edge_count"], summary["isolated_node_count"]), (3, 1, 1))
        self.assertEqual(summary["self_loop_rows_excluded"], 1)
        edge = self.table("edges.tsv")[0]
        self.assertEqual(edge, {"protein1": "9606.P1", "protein2": "9606.P2", "experimental": "600", "combined_score": "700"})
        self.assertEqual(len([row for row in self.table("mapping.tsv") if row["gene_id"] == "D"]), 2)
        self.assertEqual(self.table("nodes.tsv")[-1]["degree"], "0")
        for name, identity in summary["outputs"].items():
            self.assertEqual(hashlib.sha256((self.output / name).read_bytes()).hexdigest(), identity["sha256"])

    def test_score_threshold_direction_and_no_expansion(self):
        self.snapshot_import()
        self.extract(selection="down")
        self.assertEqual(self.summary()["status"], "NO_EDGES")
        self.assertEqual([row["string_id"] for row in self.table("nodes.tsv")], ["9606.P2"])
        self.output = self.root / "higher"
        self.extract(required_score=701)
        self.assertEqual(self.summary()["edge_count"], 0)

    def test_gzip_and_physical_declaration(self):
        compressed = self.root / "links.gz"
        with gzip.open(compressed, "wb") as handle:
            handle.write(self.links.read_bytes())
        self.snapshot_import(compressed, network_type="physical")
        self.extract()
        self.assertEqual(self.summary()["network_type"], "physical")

    def test_reject_corrupt_snapshot_and_candidates(self):
        self.snapshot_import()
        path = self.snapshot / "links.txt"
        original = path.read_bytes()
        path.write_bytes(original + b"\n")
        with self.assertRaisesRegex(ValueError, "checksum"):
            self.extract()
        path.write_bytes(original)
        (self.candidates / "all.tsv").write_text("tampered")
        with self.assertRaisesRegex(ValueError, "checksum"):
            self.extract()
        self.assertFalse(self.output.exists())

    def test_reject_wrong_taxon_path_escape_and_mapping_rate(self):
        self.snapshot_import()
        with self.assertRaisesRegex(ValueError, "mapping rate"):
            self.extract(minimum_mapping_rate=0.8)
        manifest = self.snapshot / "manifest.json"
        data = json.loads(manifest.read_text())
        data["taxon_id"] = 10090
        manifest.write_text(json.dumps(data))
        with self.assertRaisesRegex(ValueError, "taxon mismatch"):
            self.extract()
        data["taxon_id"] = 9606
        data["files"]["links"]["path"] = "../links.txt"
        manifest.write_text(json.dumps(data))
        with self.assertRaisesRegex(ValueError, "single filenames"):
            self.extract()
        self.assertFalse(self.output.exists())

    def test_import_rejects_species_mismatch_and_api_score_scale(self):
        self.mapping.write_text(self.mapping.read_text().replace("9606.P1", "10090.P1"))
        with self.assertRaisesRegex(ValueError, "wrong taxon"):
            self.snapshot_import()
        self.assertFalse(self.snapshot.exists())
        self.mapping.write_text("query_id\tinput_id_type\tstring_id\ttaxon_id\n")
        self.links.write_text("protein1 protein2 combined_score\n9606.P1 9606.P2 0.9\n")
        with self.assertRaisesRegex(ValueError, "integer values"):
            self.snapshot_import()

    def test_no_mappings_and_empty_candidates_are_distinct(self):
        self.mapping.write_text("query_id\tinput_id_type\tstring_id\ttaxon_id\n")
        self.snapshot_import()
        self.extract()
        self.assertEqual(self.summary()["status"], "NO_MAPPED_PROTEINS")
        selected = self.candidates / "down.tsv"
        selected.write_text(selected.read_text().splitlines()[0] + "\n")
        summary_path = self.candidates / "summary.json"
        data = json.loads(summary_path.read_text())
        data["counts"]["down"] = 0
        data["outputs"]["down.tsv"] = {"sha256": hashlib.sha256(selected.read_bytes()).hexdigest(), "size_bytes": selected.stat().st_size}
        summary_path.write_text(json.dumps(data))
        self.output = self.root / "empty"
        self.extract(selection="down")
        self.assertEqual(self.summary()["status"], "NO_CANDIDATES")
        self.assertIsNone(self.summary()["mapping_rate"])

    def test_conflicting_reverse_edges_are_rejected(self):
        self.links.write_text(self.links.read_text().replace("9606.P2 9606.P1 600 700", "9606.P2 9606.P1 601 700"))
        self.snapshot_import()
        with self.assertRaisesRegex(ValueError, "Conflicting duplicate"):
            self.extract()

    def test_header_only_links_preserve_evidence_columns(self):
        self.links.write_text("protein1 protein2 experimental combined_score\n")
        self.snapshot_import()
        self.extract()
        self.assertEqual(self.summary()["status"], "NO_EDGES")
        self.assertEqual((self.output / "edges.tsv").read_text(), "protein1\tprotein2\texperimental\tcombined_score\n")

    def test_mapping_namespace_is_not_inferred_from_symbol(self):
        self.mapping.write_text(self.mapping.read_text().replace("\tSYMBOL\t", "\tENTREZID\t"))
        self.snapshot_import()
        self.extract()
        self.assertEqual(self.summary()["status"], "NO_MAPPED_PROTEINS")

    def test_invalid_parameters_are_rejected(self):
        for kwargs in ({"required_score": True}, {"required_score": 1001},
                       {"minimum_mapping_rate": float("nan")}, {"selection": "neighbors"}):
            with self.subTest(kwargs=kwargs), self.assertRaises(ValueError):
                self.extract(**kwargs)

    def test_deterministic_network_and_no_overwrite(self):
        self.snapshot_import()
        self.extract()
        first = {path.name: path.read_bytes() for path in self.output.iterdir()}
        with self.assertRaises(FileExistsError):
            self.extract()
        self.output = self.root / "second"
        self.extract()
        self.assertEqual(first, {path.name: path.read_bytes() for path in self.output.iterdir()})

    def test_many_genes_per_protein_retains_original_genes(self):
        self.mapping.write_text(self.mapping.read_text() + "E\tSYMBOL\t9606.P1\t9606\n")
        self.snapshot_import()
        self.extract()
        self.assertEqual(json.loads(self.table("nodes.tsv")[0]["gene_ids"]), ["A", "E"])
        self.assertEqual(self.summary()["mapping_rate"], 0.8)

    def test_cli_offline_end_to_end(self):
        subprocess.run([sys.executable, "-m", "rnaseq.ppi_network", "snapshot", "--mapping", str(self.mapping),
                        "--links", str(self.links), "--outdir", str(self.snapshot), "--version", "12.5",
                        "--taxon-id", "9606", "--network-type", "functional", "--source-url", "https://example.org/test"], check=True)
        subprocess.run([sys.executable, "-m", "rnaseq.ppi_network", "network", "--candidates", str(self.candidates),
                        "--snapshot", str(self.snapshot), "--outdir", str(self.output)], check=True)
        self.assertEqual(self.summary()["status"], "COMPLETE")


if __name__ == "__main__":
    unittest.main()
