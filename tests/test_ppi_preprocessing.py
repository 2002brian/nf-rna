"""Offline contract tests; no R, STRING, network or runtime provisioning."""

import csv
import hashlib
import json
import tempfile
import unittest
from pathlib import Path

from rnaseq.ppi import prepare_ppi_inputs


class PpiPreprocessingTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.source = self.root / "contrast"
        self.source.mkdir()
        self.out = self.root / "prepared"
        self.rows = [
            ["ENSG1.1", "A", "2", "0.01", "-9"],
            ["ENSG1.2", "A", "-3", "0.02", "9"],
            ["ENSG2.1", "B", "0", "0.01", "0"],
            ["ENSG3.1", "C", "4", "0.05", "4"],
            ["ENSG4.1", "D", "5", "NA", "5"],
        ]
        self.meta = {
            "contrast_id": "B_vs_A", "numerator": {"level": "B", "n_samples": 3},
            "denominator": {"level": "A", "n_samples": 3},
            "deg_definition": {"padj_below": 0.05, "abs_log2fc_at_least": 1},
        }
        self.write()

    def write(self, significant=None):
        for name, rows in (("all_genes.tsv", self.rows),
                           ("significant.tsv", self.rows[:2] if significant is None else significant)):
            with (self.source / name).open("w", newline="") as handle:
                writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
                writer.writerow(["gene_id", "symbol", "log2FoldChange", "padj", "log2FoldChange_shrunk"])
                writer.writerows(rows)
        (self.source / "contrast_metadata.json").write_text(json.dumps(self.meta))

    def prepare(self, id_type="ENSEMBL"):
        return prepare_ppi_inputs(self.source, self.out, organism="Homo sapiens", input_id_type=id_type)

    def test_selection_direction_collisions_and_hashes(self):
        before = {path.name: path.read_bytes() for path in self.source.iterdir()}
        self.prepare()
        summary = json.loads((self.out / "summary.json").read_text())
        self.assertEqual(summary["counts"], {"all": 2, "up": 1, "down": 1})
        self.assertEqual(summary["taxon_id"], 9606)
        self.assertEqual(summary["numerator"]["level"], "B")
        self.assertEqual(summary["mapping_status"], "NOT_RUN")
        self.assertEqual(summary["query_id_collisions"], {"ENSG1": ["ENSG1.1", "ENSG1.2"]})
        with (self.out / "up.tsv").open() as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual(rows[0]["gene_id"], "ENSG1.1")  # raw LFC, not shrunken LFC
        self.assertEqual(rows[0]["query_id"], "ENSG1")
        for name, data in before.items():
            self.assertEqual((self.source / name).read_bytes(), data)
            self.assertEqual((self.out / f"source_{name}").read_bytes(), data)
            self.assertEqual(summary["sources"][name]["sha256"], hashlib.sha256(data).hexdigest())
        for name, record in summary["outputs"].items():
            self.assertEqual(hashlib.sha256((self.out / name).read_bytes()).hexdigest(), record["sha256"])

    def test_deterministic_bundle(self):
        self.prepare()
        first = {path.name: path.read_bytes() for path in self.out.iterdir()}
        self.out = self.root / "other"
        self.prepare()
        self.assertEqual(first, {path.name: path.read_bytes() for path in self.out.iterdir()})

    def test_no_overwrite_even_empty_destination(self):
        self.out.mkdir()
        with self.assertRaises(FileExistsError):
            self.prepare()

    def test_empty_is_valid_and_zero_lfc_is_retained_only_in_all(self):
        self.rows = [self.rows[2]]
        self.write(significant=[])
        self.prepare()
        self.assertEqual(json.loads((self.out / "summary.json").read_text())["status"], "EMPTY")
        self.out = self.root / "zero"
        self.meta["deg_definition"]["abs_log2fc_at_least"] = 0
        self.write(significant=self.rows)
        self.prepare()
        self.assertEqual(json.loads((self.out / "summary.json").read_text())["counts"], {"all": 1, "up": 0, "down": 0})

    def test_reject_inconsistent_set_or_statistics(self):
        for selected in ([], [self.rows[0]], [self.rows[0], [*self.rows[1][:2], "-4", "0.02", "9"]]):
            with self.subTest(selected=selected):
                self.write(significant=selected)
                with self.assertRaises(ValueError):
                    self.prepare()
                self.assertFalse(self.out.exists())

    def test_reject_duplicate_gene_and_nonfinite_statistic(self):
        self.rows.append(self.rows[0])
        self.write()
        with self.assertRaises(ValueError):
            self.prepare()
        self.rows.pop()
        self.rows[0][2] = "inf"
        self.write()
        with self.assertRaises(ValueError):
            self.prepare()

    def test_symbols_are_preserved_without_version_stripping(self):
        self.prepare(id_type="SYMBOL")
        with (self.out / "all.tsv").open() as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual(rows[0]["query_id"], rows[0]["gene_id"])

    def test_reject_invalid_threshold_and_organism(self):
        self.meta["deg_definition"]["abs_log2fc_at_least"] = float("nan")
        self.write()
        with self.assertRaises(ValueError):
            self.prepare()
        with self.assertRaises(ValueError):
            prepare_ppi_inputs(self.source, self.out, organism="unknown", input_id_type="SYMBOL")


if __name__ == "__main__":
    unittest.main()
