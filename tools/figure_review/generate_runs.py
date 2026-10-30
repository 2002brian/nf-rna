#!/usr/bin/env python3
"""Generate a reproducible set of nf-rna figures for review or comparison.

Runs the shipped example project through L1 + L2 (its own heatmap, volcano, PCA,
correlation and boxplot figures) and a small synthetic mouse GO ORA + GSEA-GO
fixture (dotplots), using only the installed ``rnaseq`` package and the R scripts
shipped inside it.  Nothing outside nf-rna is read.  KEGG runs only with --kegg-scenario: an
OFFLINE synthetic fixture, or ``live`` for the real KEGG REST service (read-only integration check).

    generate_runs.py --example PROJECT_DIR --out RUN_DIR [--figures "{preset: legacy, profile: nature}"]

``--figures`` is written into the copied project's ``project.yaml`` (and into the
enrichment configs as ``figure_style``); omit it to reproduce a project that
predates the figures section.  ``Rscript`` (with ggplot2, pheatmap, DESeq2,
clusterProfiler and org.Mm.eg.db) must be on PATH.
"""

from __future__ import annotations

import argparse
import json
import shutil
import subprocess
import sys
from importlib import resources
from pathlib import Path

import yaml

GO_CONFIG = {"pvalue_cutoff": 0.05, "qvalue_cutoff": 0.2, "p_adjust_method": "BH"}
GSEA_CONFIG = {"minimum_ranked_genes": 50, "min_gs_size": 10, "max_gs_size": 500, "pvalue_cutoff": 0.05, "padj_cutoff": 0.05, "p_adjust_method": "BH", "seed": 1}

FIXTURE = '''
suppressMessages(library(org.Mm.eg.db)); set.seed(1)
terms <- function(go) unique(AnnotationDbi::select(org.Mm.eg.db, keys=go, keytype="GOALL", columns="ENTREZID")$ENTREZID)
immune <- head(sort(terms("GO:0006955")), 150)
cycle <- head(sort(setdiff(terms("GO:0007049"), immune)), 150)
bg <- head(setdiff(sort(AnnotationDbi::keys(org.Mm.eg.db, "ENTREZID")), c(immune, cycle)), 2500)
genes <- c(immune, cycle, bg)
stat <- c(rnorm(length(immune), 3, 1), rnorm(length(cycle), -3, 1), rnorm(length(bg), 0, 1))
write.table(data.frame(gene_id=genes, pvalue=0.5, padj=0.9, stat=stat, log2FoldChange=stat/3), "{base}/all_genes.tsv", sep="\\t", quote=FALSE, row.names=FALSE)
write.table(data.frame(gene_id=immune), "{base}/significant.tsv", sep="\\t", quote=FALSE, row.names=FALSE)
write.table(data.frame(gene_id=character()), "{base}/empty.tsv", sep="\\t", quote=FALSE, row.names=FALSE)
'''


def run_l1_l2(example: Path, out: Path, figures: dict | None) -> str:
    from rnaseq.l2 import execute_l2, prepare_l2
    from rnaseq.validators import validate_project

    project = out / "simple_two_group"
    shutil.copytree(example, project, ignore=shutil.ignore_patterns("downstream", "planning", "runs", "delivery"))
    if figures is not None:
        config = yaml.safe_load((project / "project.yaml").read_text(encoding="utf-8"))
        config["figures"] = figures
        (project / "project.yaml").write_text(yaml.safe_dump(config, sort_keys=False), encoding="utf-8")
    result = execute_l2(prepare_l2(validate_project(project), run_id=None))
    return str(result.status)


KEGG_ORA = {"p_adjust_method": "BH", "pvalue_cutoff": 0.05, "qvalue_cutoff": 0.2, "min_gs_size": 10, "max_gs_size": 500}
KEGG_GSEA = {**GSEA_CONFIG, "padj_cutoff": 0.05}
HERE = Path(__file__).resolve().parent


def _synthetic_kegg_snapshot(base: Path) -> Path:
    """A SYNTHETIC KEGG snapshot for the offline fixture: real pathway names, synthetic gene links."""

    import hashlib

    root = base / "kegg_snapshots"
    target = root / "synthetic_fixture"
    target.mkdir(parents=True)
    pathways = [line.split("\t") for line in (HERE / "fixtures" / "kegg_list_pathway_mmu.txt").read_text(encoding="utf-8").splitlines() if line]
    genes = [line.split("\t")[0] for line in (base / "all_genes.tsv").read_text(encoding="utf-8").splitlines()[1:]]
    names = [(pid.removeprefix("path:"), name.removesuffix(" - Mus musculus (house mouse)")) for pid, name in pathways]
    pairs = [(pid, gene) for index, (pid, _name) in enumerate(names) for gene in genes[index::len(names)]]
    (target / "TERM2GENE.tsv").write_text("term\tgene\n" + "".join(f"{a}\t{b}\n" for a, b in pairs), encoding="utf-8")
    (target / "TERM2NAME.tsv").write_text("term\tname\n" + "".join(f"{a}\t{b}\n" for a, b in names), encoding="utf-8")
    files = {name: {"sha256": hashlib.sha256((target / name).read_bytes()).hexdigest()} for name in ("TERM2GENE.tsv", "TERM2NAME.tsv")}
    manifest = {"schema": "nf-rna.kegg-snapshot.v1", "snapshot_id": "synthetic_fixture", "organism": "mmu", "retrieved_at_utc": "1970-01-01T00:00:00Z",
                "kegg_release": "SYNTHETIC offline fixture; not KEGG data", "source_urls": {}, "files": files,
                "pathways": len(names), "genes": len(set(genes))}
    (target / "manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    return root


def run_enrichment(out: Path, figures: dict | None, kegg_scenario: str | None = None, only_kegg: bool = False, empty: bool = False) -> dict[str, int]:
    base = out / "enrichment"
    base.mkdir(parents=True)
    subprocess.run(["Rscript", "-e", FIXTURE.replace("{base}", str(base))], check=True)
    if empty:
        # No foreground genes and a ranking without signal: every enrichment backend reports "nothing significant".
        import random

        (base / "significant.tsv").write_text("gene_id\n", encoding="utf-8")
        rng = random.Random(7)
        lines = (base / "all_genes.tsv").read_text(encoding="utf-8").splitlines()
        header = lines[0].split("\t")
        rewritten = [lines[0]]
        for line in lines[1:]:
            fields = line.split("\t")
            noise = rng.gauss(0, 1)
            fields[header.index("stat")] = f"{noise:.6g}"
            fields[header.index("log2FoldChange")] = f"{noise / 3:.6g}"
            rewritten.append("\t".join(fields))
        (base / "all_genes.tsv").write_text("\n".join(rewritten) + "\n", encoding="utf-8")
    r_dir = Path(str(resources.files("rnaseq.r")))
    runtime = lambda script: {"script": script, "source_revision": "figure-review", "container_image": "none", "output_schema": "nf-rna.scientific-provenance.v1"}
    contrast = {"contrast_id": "Treatment_vs_Control", "numerator": "133Gy", "denominator": "22Gy", "all_genes": f"{base}/all_genes.tsv",
                "significant": f"{base}/significant.tsv", "up": f"{base}/empty.tsv", "down": f"{base}/empty.tsv"}
    annotation = {"organism": "Mus musculus", "input_id_type": "ENTREZID", "target_id_type": "ENTREZID", "mapping_warning_rate": 0, "minimum_mapping_rate": 0, "minimum_mapped_foreground": 5}
    configs = {
        "go": ("go_analysis.R", {"output_dir": f"{base}/go", "orgdb_package": "org.Mm.eg.db", "annotation": {**annotation, "enrichment": {"go": GO_CONFIG}}, "contrasts": [contrast], "runtime": runtime("go_analysis.R")}),
        "gsea_go": ("gsea_analysis.R", {"output_dir": f"{base}/gsea_go", "orgdb_package": "org.Mm.eg.db", "annotation": {**annotation, "gene_symbol_output": True, "enrichment": {"gsea": GSEA_CONFIG}}, "contrasts": [contrast], "runtime": runtime("gsea_analysis.R")}),
    }
    if only_kegg:
        configs = {}
    env = None
    if kegg_scenario:
        # OFFLINE KEGG fixture (synthetic statistics over real pathway names): the real backend runs with the
        # network touch points replaced through R_PROFILE_USER; see kegg_offline_profile.R.
        rows = (base / "significant.tsv").read_text(encoding="utf-8").splitlines()
        half = 1 + (len(rows) - 1) // 2
        (base / "kegg_up.tsv").write_text("\n".join(rows[:half]) + "\n", encoding="utf-8")
        (base / "kegg_down.tsv").write_text("\n".join([rows[0], *rows[half:]]) + "\n", encoding="utf-8")
        kegg_contrast = {**contrast, "up": f"{base}/kegg_up.tsv", "down": f"{base}/kegg_down.tsv"}
        kegg_annotation = {**annotation, "gene_symbol_output": True, "enrichment": {"kegg": {"ora": KEGG_ORA, "gsea": KEGG_GSEA}}}
        if kegg_scenario == "live":
            from rnaseq.enrichment_config import kegg_snapshot_root

            kegg_resource = {"organism_code": "mmu", "provider": "kegg_rest_snapshot", "snapshot": None, "snapshot_root": str(kegg_snapshot_root("mmu"))}
        else:
            kegg_resource = {"organism_code": "mmu", "provider": "kegg_rest_snapshot", "snapshot": "synthetic_fixture", "snapshot_root": str(_synthetic_kegg_snapshot(base))}
        for name, mode in (("kegg", "ora"), ("gsea_kegg", "gsea")):
            configs[name] = ("kegg_analysis.R", {"mode": mode, "output_dir": f"{base}/{name}", "orgdb_package": "org.Mm.eg.db", "annotation": kegg_annotation, "kegg": kegg_resource,
                                                 "contrasts": [kegg_contrast], "runtime": runtime("kegg_analysis.R")})
        import os
        env = dict(os.environ) if kegg_scenario == "live" else {**os.environ, "R_PROFILE_USER": str(HERE / "kegg_offline_profile.R"), "NF_RNA_KEGG_FIXTURE_DIR": str(HERE / "fixtures"), "NF_RNA_KEGG_SCENARIO": kegg_scenario}
    codes = {}
    for name, (script, config) in configs.items():
        if figures is not None:
            config = {**config, "figure_style": {**figures, "source": "explicit"}}
        path = base / f"{name}_config.json"
        path.write_text(json.dumps(config), encoding="utf-8")
        result = subprocess.run(["Rscript", str(r_dir / script), "--config", str(path)], capture_output=True, text=True, env=env if name in ("kegg", "gsea_kegg") else None)
        codes[name] = result.returncode
        if result.returncode:
            (base / f"{name}.stderr.log").write_text(result.stderr, encoding="utf-8")
    return codes


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--example", type=Path, required=True, help="the examples/simple_two_group project (read only; it is copied)")
    parser.add_argument("--out", type=Path, required=True, help="new, empty or absent output directory")
    parser.add_argument("--kegg-scenario", choices=["many", "few", "long", "short", "nonsig", "empty", "live"], help="also run the KEGG ORA and GSEA backends against the OFFLINE synthetic fixture (kegg_offline_profile.R); 'live' instead uses the real KEGG REST service (read-only GET requests; separate integration check, never part of the test suite)")
    parser.add_argument("--skip-analysis", action="store_true", help="skip L1 + L2 (enrichment only)")
    parser.add_argument("--only-kegg", action="store_true", help="in the enrichment step run only the KEGG backends (needs --kegg-scenario)")
    parser.add_argument("--empty-foreground", action="store_true", help="enrichment fixture with no foreground genes and a signal-free ranking")
    parser.add_argument("--skip-enrichment", action="store_true", help="only run L1 + L2 (no synthetic GO/GSEA fixture)")
    parser.add_argument("--figures", help="YAML mapping for the project's figures section, e.g. '{preset: standard, profile: general}'")
    args = parser.parse_args(argv)
    if args.out.exists() and any(args.out.iterdir()):
        parser.error(f"--out must be empty or absent: {args.out}")
    args.out.mkdir(parents=True, exist_ok=True)
    figures = yaml.safe_load(args.figures) if args.figures else None
    import rnaseq

    print(f"rnaseq imported from {Path(rnaseq.__file__).parent}")
    status = "SUCCESS" if args.skip_analysis else run_l1_l2(args.example.resolve(), args.out, figures)
    codes = {} if args.skip_enrichment else run_enrichment(args.out, figures, args.kegg_scenario, args.only_kegg, args.empty_foreground)
    print(f"L1/L2 status: {status}; enrichment exit codes: {codes}")
    return 0 if status == "SUCCESS" and not any(codes.values()) else 1


if __name__ == "__main__":
    sys.exit(main())
