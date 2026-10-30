"""The one enrichment R-backend configuration.

Both the Nextflow ENRICHMENT_ANALYSIS task (``workflow_support.enrichment_config``)
and the Python entry points (``rnaseq.go/gsea/kegg.execute_*``) build their
``--config`` JSON here, so the two paths cannot drift apart.
"""

from __future__ import annotations

import os
from pathlib import Path
from typing import Any, Iterable

from rnaseq.models import ANNOTATION_DATABASES

KEGG_CODES = {"Homo sapiens": "hsa", "Mus musculus": "mmu"}
BACKEND_DIRECTORIES = {"go": "go", "kegg": "kegg", "gsea-go": "gsea_go", "gsea-kegg": "gsea_kegg"}
BACKEND_SCRIPTS = {"go": "go_analysis.R", "kegg": "kegg_analysis.R", "gsea-go": "gsea_analysis.R", "gsea-kegg": "kegg_analysis.R"}


def kegg_snapshot_root(organism_code: str) -> Path:
    """User-cache directory of KEGG gene-set snapshots (never the repository)."""

    return Path(os.environ.get("XDG_CACHE_HOME", Path.home() / ".cache")) / "nf-rna" / "genesets" / "kegg" / organism_code


def kegg_snapshot_available(organism_code: str, snapshot: str | None) -> bool:
    """True when the configured snapshot, or any snapshot, exists locally."""

    root = kegg_snapshot_root(organism_code)
    if snapshot:
        return (root / snapshot / "manifest.json").is_file()
    return root.is_dir() and any((child / "manifest.json").is_file() for child in root.iterdir() if not child.name.startswith("."))


def backend_config(
    kind: str, *, annotation: dict[str, Any], contrasts: Iterable[dict[str, Any]], l2_dir: Path, output_dir: Path,
    figure_style: dict[str, Any], runtime: dict[str, Any] | None = None,
) -> dict[str, Any]:
    """Return the R backend config of one enrichment module (go, kegg, gsea-go, gsea-kegg)."""

    if kind not in BACKEND_DIRECTORIES:
        raise ValueError(f"unsupported enrichment module: {kind}")
    organism = annotation["organism"]
    orgdb = ANNOTATION_DATABASES.get(organism)
    if orgdb is None:
        raise ValueError(f"unsupported enrichment organism: {organism}")
    items = []
    for item in contrasts:
        root = l2_dir / "contrasts" / item["contrast_id"]
        items.append({
            "contrast_id": item["contrast_id"],
            # Figure labelling and direction only: the sign convention of the contrast.
            "numerator": item.get("numerator"),
            "denominator": item.get("denominator"),
            "all_genes": str(root / "all_genes.tsv"),
            "significant": str(root / "significant.tsv"),
            "up": str(root / "upregulated.tsv"),
            "down": str(root / "downregulated.tsv"),
        })
    cfg: dict[str, Any] = {
        "output_dir": str(output_dir), "annotation": annotation, "orgdb_package": orgdb, "contrasts": items,
        "figure_style": figure_style,
    }
    if runtime is not None:
        cfg["runtime"] = runtime
    if kind in {"kegg", "gsea-kegg"}:
        code = KEGG_CODES[organism]
        kegg = annotation.get("enrichment", {}).get("kegg", {})
        cfg.update({"mode": "ora" if kind == "kegg" else "gsea", "kegg": {
            "organism_code": code, "provider": "kegg_rest_snapshot",
            "snapshot": kegg.get("snapshot"), "snapshot_root": str(kegg_snapshot_root(code)),
        }})
    return cfg
