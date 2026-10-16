#!/usr/bin/env python3
"""Scientific-invariance audit between two nf-rna runs of the same project (stdlib only; review aid).

Every file of the two run trees is compared byte for byte except the figure artefacts themselves (PNG / PDF / TIFF and
``figure_manifest.json``, which legitimately differ between presets).  Files are grouped into the scientific areas the
figure standard must not touch (normalisation, differential expression, gene selection, clustering/ordering,
correlation, enrichment, figure source data).  Tables (``.tsv`` / ``.csv``) must be byte-identical.  Metadata files
(JSON / YAML / Markdown / logs) legitimately mention the run directory, timestamps and the selected figure style, so
they are compared after normalising the run root and ISO timestamps, and every remaining difference is listed.

    invariance_audit.py --legacy RUN_A --candidate RUN_B [--out report.json]
"""

from __future__ import annotations

import argparse
import difflib
import hashlib
import json
import re
import sys
from pathlib import Path

FIGURE_SUFFIXES = {".png", ".pdf", ".tiff", ".tif"}
AREAS = [
    ("normalisation", re.compile(r"(normalized_counts|vst\.(tsv|csv)|source_counts|library_size_qc|filtering_summary|sample_distribution_qc)")),
    ("differential expression", re.compile(r"l2/contrasts/[^/]+/(all_genes|significant|upregulated|downregulated)\.tsv$")),
    ("gene selection", re.compile(r"(heatmap_genes|heatmap_zscores)\.tsv$")),
    ("clustering / sample ordering / PCA", re.compile(r"(pca_scores|pca_variance)\.tsv$")),
    ("correlation values", re.compile(r"sample_correlation\.tsv$")),
    ("enrichment statistics and term tables", re.compile(r"enrichment/.*(GO_|all_terms|significant|positive_enrichment|negative_enrichment|kegg\.tsv|kegg_all_terms|gene_members|core_members|members|ranked_gene_list|rank_mapping|kegg_mapping|gene_mapping)")),
    ("figure source-data tables", re.compile(r"(source_data|_source)\.tsv$|dotplot_source_data\.tsv$")),
]
TIMESTAMP = re.compile(r"\d{4}-\d{2}-\d{2}[T ]\d{2}:\d{2}:\d{2}(\.\d+)?(Z|[+-]\d{2}:?\d{2})?")


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def is_figure_artifact(path: Path) -> bool:
    return path.suffix.lower() in FIGURE_SUFFIXES or path.name == "figure_manifest.json"


def area_of(rel: str) -> str:
    for name, pattern in AREAS:
        if pattern.search(rel):
            return name
    return "other (metadata, logs, inputs)"


def normalised(path: Path, root: Path, other_root: Path) -> list[str]:
    text = path.read_text(encoding="utf-8", errors="replace")
    for token in {str(root), str(other_root), str(root.resolve()), str(other_root.resolve())}:
        text = text.replace(token, "<RUN>")
    return TIMESTAMP.sub("<TIME>", text).splitlines()


def leaf_differences(x, y, path="") -> list[str]:
    """Paths of the leaves at which two parsed JSON/YAML documents differ."""

    if isinstance(x, dict) and isinstance(y, dict):
        out = []
        for key in sorted(set(x) | set(y)):
            if key not in x or key not in y:
                out.append(f"{path}/{key} (only in {'candidate' if key in y else 'legacy'})")
            else:
                out += leaf_differences(x[key], y[key], f"{path}/{key}")
        return out
    if isinstance(x, list) and isinstance(y, list) and len(x) == len(y):
        return [d for i, (u, v) in enumerate(zip(x, y)) for d in leaf_differences(u, v, f"{path}[{i}]")]
    return [] if x == y else [f"{path}: {str(x)[:40]!r} -> {str(y)[:40]!r}"]


def structured(path: Path, root: Path, other_root: Path):
    import yaml

    text = "\n".join(normalised(path, root, other_root))
    return json.loads(text) if path.suffix == ".json" else yaml.safe_load(text)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--legacy", type=Path, required=True)
    parser.add_argument("--candidate", type=Path, required=True)
    parser.add_argument("--out", type=Path)
    args = parser.parse_args(argv)
    a, b = args.legacy.resolve(), args.candidate.resolve()
    files_a = {p.relative_to(a).as_posix() for p in a.rglob("*") if p.is_file()}
    files_b = {p.relative_to(b).as_posix() for p in b.rglob("*") if p.is_file()}
    report = {"only_in_legacy": sorted(files_a - files_b), "only_in_candidate": sorted(files_b - files_a), "areas": {}, "excluded_figure_artifacts": 0, "metadata_differences": {}}
    for rel in sorted(files_a & files_b):
        pa, pb = a / rel, b / rel
        if is_figure_artifact(pa):
            report["excluded_figure_artifacts"] += 1
            continue
        area = area_of(rel)
        entry = report["areas"].setdefault(area, {"files": 0, "byte_identical": 0, "different": []})
        entry["files"] += 1
        if sha(pa) == sha(pb):
            entry["byte_identical"] += 1
            continue
        if pa.suffix in {".tsv", ".csv"}:
            entry["different"].append({"file": rel, "kind": "TABLE DIFFERS"})
            continue
        left, right = normalised(pa, a, b), normalised(pb, b, a)
        if pa.suffix in {".json", ".yaml", ".yml"}:
            try:
                diff = leaf_differences(structured(pa, a, b), structured(pb, b, a))
            except Exception:  # not parseable: fall back to a line diff
                diff = [line for line in difflib.unified_diff(left, right, lineterm="", n=0) if line[:1] in "+-" and line[:3] not in ("+++", "---")]
        else:
            diff = [line for line in difflib.unified_diff(left, right, lineterm="", n=0) if line[:1] in "+-" and line[:3] not in ("+++", "---")]
        if diff:
            report["metadata_differences"][rel] = diff[:12]
            entry["different"].append({"file": rel, "kind": "metadata differs after normalisation", "lines": len(diff)})
        else:
            entry["byte_identical"] += 1  # differs only in run root / timestamps
            entry.setdefault("equal_after_normalisation", []).append(rel)
    tables_ok = all(not [d for d in v["different"] if d["kind"] == "TABLE DIFFERS"] for v in report["areas"].values())
    style_only = re.compile(r"figure_style|figures|preset|profile|/source|style")
    unexpected = {rel: [d for d in lines if not style_only.search(d)] for rel, lines in report["metadata_differences"].items()}
    report["metadata_differences_beyond_figure_style"] = {rel: lines for rel, lines in unexpected.items() if lines}
    report["tables_byte_identical"] = tables_ok
    report["structure_identical"] = not report["only_in_legacy"] and not report["only_in_candidate"]
    if args.out:
        args.out.write_text(json.dumps(report, indent=1, sort_keys=True) + "\n", encoding="utf-8")
    for name, v in sorted(report["areas"].items()):
        flagged = [d for d in v["different"]]
        print(f"{name}: {v['files']} files, {v['byte_identical']} identical" + (f", DIFFERENT: {[d['file'] for d in flagged]}" if flagged else ""))
    print(f"figure artefacts excluded: {report['excluded_figure_artifacts']}; only in legacy: {len(report['only_in_legacy'])}; only in candidate: {len(report['only_in_candidate'])}; tables byte-identical: {tables_ok}; metadata files differing: {len(report['metadata_differences'])}, of which beyond the figure-style selection: {len(report['metadata_differences_beyond_figure_style'])}")
    return 0 if tables_ok and report["structure_identical"] and not report["metadata_differences_beyond_figure_style"] else 1


if __name__ == "__main__":
    sys.exit(main())
