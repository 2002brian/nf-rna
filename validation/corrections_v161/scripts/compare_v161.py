"""Before/after and independent-reference comparison for the v1.6.1 scientific corrections.

Inputs (outside Git):
  v1.6.0 runs   $V/nfrna/<dataset>/runs/VAL-<case>/<stamp>/downstream/l2   (Phase 2, frozen)
  fixed runs    $W/nfrna/<dataset>/runs/FIX-<case>/<stamp>/downstream/l2
  reference     $V/reference/outputs/<dataset>/contrasts/<contrast>        (Phase 2 blind reference)

Writes TSV tables to the output directory given on the command line.
Usage: python -I compare_v161.py <V> <W> <outdir>
"""
from __future__ import annotations

import csv
import hashlib
import json
import math
import statistics
import sys
from pathlib import Path

V, W, OUT = (Path(arg) for arg in sys.argv[1:4])
OUT.mkdir(parents=True, exist_ok=True)
csv.field_size_limit(1 << 30)

RUNS = {  # case -> (dataset, v1.6.0 case, fixed case, reference GSEA seed)
    "D1": ("D1", "VAL-D1", "FIX-D1", 1), "D1a": ("D1a", "VAL-D1a", "FIX-D1a", 1), "D1b": ("D1b", "VAL-D1b", "FIX-D1b", 1),
    "D2": ("D2", "VAL-D2", "FIX-D2", 1), "D1-S42": ("D1s42", None, "FIX-D1-S42", 42),
}
ORA_COLLECTIONS = {"GO_BP": ("go", "GO_BP_all_terms.tsv"), "GO_MF": ("go", "GO_MF_all_terms.tsv"), "GO_CC": ("go", "GO_CC_all_terms.tsv"), "KEGG": ("kegg", "kegg_all_terms.tsv")}
GSEA_COLLECTIONS = {"GO_BP": ("gsea_go", "BP"), "GO_MF": ("gsea_go", "MF"), "GO_CC": ("gsea_go", "CC"), "KEGG": ("gsea_kegg", None)}
FOREGROUNDS = ("up", "down", "significant")


def l2_dir(base: Path, project: str, case: str) -> Path:
    stamps = sorted((base / "nfrna" / project / "runs" / case).glob("*/downstream/l2"))
    if len(stamps) != 1:
        raise SystemExit(f"expected one run for {project}/{case}, found {stamps}")
    return stamps[0]


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def num(value: str | None) -> float:
    return math.nan if value in (None, "", "NA") else float(value)


def rel(a: float, b: float) -> float:
    if math.isnan(a) and math.isnan(b):
        return 0.0
    if math.isnan(a) or math.isnan(b):
        return math.inf
    return 0.0 if a == b else abs(a - b) / max(abs(a), abs(b))


def max_rel(pairs) -> float:
    values = [rel(a, b) for a, b in pairs]
    return max(values) if values else 0.0


def write(name: str, rows: list[dict]) -> None:
    if not rows:
        return
    with (OUT / name).open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows({k: (f"{v:.3g}" if isinstance(v, float) else v) for k, v in row.items()} for row in rows)


def deg_set(rows: list[dict[str, str]]) -> set[str]:
    return {row["gene_id"] for row in rows}


def contrasts(l2: Path) -> list[str]:
    return sorted(path.name for path in (l2 / "contrasts").iterdir() if path.is_dir())


# ---------------------------------------------------------------- M-1: DESeq2 alpha
de_rows, ora_rows = [], []
for case in ("D1", "D1a", "D1b", "D2"):
    project, old_case, new_case, _ = RUNS[case]
    old_l2, new_l2 = l2_dir(V, project, old_case), l2_dir(W, project, new_case)
    for cid in contrasts(new_l2):
        old_dir, new_dir = old_l2 / "contrasts" / cid, new_l2 / "contrasts" / cid
        ref_dir = V / "reference" / "outputs" / project / "contrasts" / cid
        old = {r["gene_id"]: r for r in read_tsv(old_dir / "all_genes.tsv")}
        new = {r["gene_id"]: r for r in read_tsv(new_dir / "all_genes.tsv")}
        ref = {r["gene_id"]: r for r in read_tsv(ref_dir / "sensitivity_alpha005" / "all_genes.tsv")}
        summary = json.loads((new_dir / "backend_summary.json").read_text())
        thresholds = json.loads((ref_dir / "filter_threshold.json").read_text())
        old_na = {g for g, r in old.items() if r["padj"] == "NA"}
        new_na = {g for g, r in new.items() if r["padj"] == "NA"}
        ref_na = {g for g, r in ref.items() if r["padj"] == "NA"}
        old_deg, new_deg = deg_set(read_tsv(old_dir / "significant.tsv")), deg_set(read_tsv(new_dir / "significant.tsv"))
        ref_deg = deg_set(read_tsv(ref_dir / "sensitivity_alpha005" / "significant.tsv"))
        new_up, new_down = deg_set(read_tsv(new_dir / "upregulated.tsv")), deg_set(read_tsv(new_dir / "downregulated.tsv"))
        old_up, old_down = deg_set(read_tsv(old_dir / "upregulated.tsv")), deg_set(read_tsv(old_dir / "downregulated.tsv"))
        de_rows.append({
            "dataset": case, "contrast": cid, "alpha_recorded": summary.get("independent_filtering_alpha"),
            "filter_threshold_v160_alpha0.1_ref": thresholds["alpha_0.1"]["filterThreshold"],
            "filter_threshold_fixed": summary.get("independent_filtering_threshold"),
            "filter_threshold_ref_alpha0.05": thresholds["alpha_0.05"]["filterThreshold"],
            "filter_threshold_rel_diff_fixed_vs_ref": rel(float(summary["independent_filtering_threshold"]), thresholds["alpha_0.05"]["filterThreshold"]),
            "padj_na_v160": len(old_na), "padj_na_fixed": len(new_na), "padj_na_ref_alpha0.05": len(ref_na),
            "padj_na_pattern_fixed_equals_ref": new_na == ref_na, "padj_na_added": len(new_na - old_na), "padj_na_removed": len(old_na - new_na),
            "pvalue_identical_v160_fixed": all(old[g]["pvalue"] == new[g]["pvalue"] and old[g]["stat"] == new[g]["stat"] and old[g]["log2FoldChange"] == new[g]["log2FoldChange"] for g in old),
            "padj_max_rel_diff_fixed_vs_ref": max_rel((num(new[g]["padj"]), num(ref[g]["padj"])) for g in new),
            "deg_v160": len(old_deg), "deg_fixed": len(new_deg), "deg_gained": len(new_deg - old_deg), "deg_lost": len(old_deg - new_deg),
            "up_v160": len(old_up), "up_fixed": len(new_up), "down_v160": len(old_down), "down_fixed": len(new_down),
            "deg_ref_alpha0.05": len(ref_deg), "deg_fixed_equals_ref": new_deg == ref_deg,
        })
        ref_ora = ref_dir / "sensitivity_alpha005" / "ora"
        if not (new_l2 / "enrichment" / "go" / cid).is_dir():
            continue
        for coll, (module, filename) in ORA_COLLECTIONS.items():
            for fg in FOREGROUNDS:
                def load(l2: Path) -> dict[str, dict[str, str]]:
                    path = l2 / "enrichment" / module / cid / fg / filename
                    return {r["ID"]: r for r in read_tsv(path)} if path.is_file() else {}
                o, n = load(old_l2), load(new_l2)
                r_terms = {r["term_id"]: r for r in read_tsv(ref_ora / coll / fg / "terms.tsv")} if (ref_ora / coll / fg / "terms.tsv").is_file() else {}
                r_summary = json.loads((ref_ora / coll / fg / "summary.json").read_text()) if (ref_ora / coll / fg / "summary.json").is_file() else {}
                sig = lambda table, key="significant": {t for t, r in table.items() if r[key] == "TRUE"}  # noqa: E731
                o_sig, n_sig, r_sig = sig(o), sig(n), sig(r_terms)
                bg = lambda table: {r["BgRatio"].split("/")[1] for r in table.values()}  # noqa: E731
                fgn = lambda table: {r["GeneRatio"].split("/")[1] for r in table.values()}  # noqa: E731
                common = set(n) & set(r_terms)
                ora_rows.append({
                    "dataset": case, "contrast": cid, "collection": coll, "foreground": fg,
                    "universe_v160": "/".join(sorted(bg(o))), "universe_fixed": "/".join(sorted(bg(n))), "universe_ref_alpha0.05": r_summary.get("N", ""),
                    "foreground_v160": "/".join(sorted(fgn(o))), "foreground_fixed": "/".join(sorted(fgn(n))), "foreground_ref_alpha0.05": r_summary.get("n", ""),
                    "significant_v160": len(o_sig), "significant_fixed": len(n_sig), "gained": len(n_sig - o_sig), "lost": len(o_sig - n_sig),
                    "significant_ref_alpha0.05": len(r_sig), "fixed_equals_ref": n_sig == r_sig,
                    "terms_fixed": len(n), "terms_ref": len(r_terms), "padj_max_rel_diff_fixed_vs_ref": max_rel((num(n[t]["p.adjust"]), num(r_terms[t]["p.adjust"])) for t in common),
                    "gained_terms": ";".join(sorted(n_sig - o_sig)), "lost_terms": ";".join(sorted(o_sig - n_sig)),
                })
write("m1_deseq2_alpha.tsv", de_rows)
write("m1_ora.tsv", ora_rows)


# ---------------------------------------------------------------- D-1/D-2: GSEA
def gsea_table(l2: Path, cid: str, coll: str) -> dict[str, dict[str, str]]:
    module, sub = GSEA_COLLECTIONS[coll]
    root = l2 / "enrichment" / module / cid
    path = (root / sub if sub else root) / "all_terms.tsv"
    return {r["ID"]: r for r in read_tsv(path)}


def as_set(value: str) -> frozenset[str]:
    return frozenset(x for x in value.split("/") if x)


gsea_rows = []
for case in ("D1", "D1-S42", "D2"):
    project, old_case, new_case, seed = RUNS[case]
    new_l2 = l2_dir(W, project, new_case)
    old_l2 = l2_dir(V, project, old_case) if old_case else None
    provenance = json.loads((new_l2 / "enrichment" / "gsea_go" / "gsea_backend_summary.json").read_text())["gsea"]
    for cid in contrasts(new_l2):
        for coll in GSEA_COLLECTIONS:
            n = gsea_table(new_l2, cid, coll)
            ref = {r["set_id"]: r for r in read_tsv(V / "reference" / "outputs" / project.removesuffix("s42") / "contrasts" / cid / "gsea" / coll / f"seed_{seed}.tsv")}
            common = sorted(set(n) & set(ref))
            n_sig = {t for t, r in n.items() if r["significant"] == "TRUE"}
            r_sig = {t for t, r in ref.items() if num(r["padj"]) <= 0.05}
            neg = [t for t in common if num(n[t]["enrichmentScore"]) < 0]
            row = {
                "case": case, "contrast": cid, "collection": coll, "seed_recorded": provenance["seed"], "reference_seed": seed,
                "terms_fixed": len(n), "terms_ref": len(ref), "terms_common": len(common),
                "es_max_rel_vs_ref": max_rel((num(n[t]["enrichmentScore"]), num(ref[t]["ES"])) for t in common),
                "nes_max_rel_vs_ref": max_rel((num(n[t]["NES"]), num(ref[t]["NES"])) for t in common),
                "pvalue_max_rel_vs_ref": max_rel((num(n[t]["pvalue"]), num(ref[t]["pval"])) for t in common),
                "padj_max_rel_vs_ref": max_rel((num(n[t]["p.adjust"]), num(ref[t]["padj"])) for t in common),
                "significant_fixed": len(n_sig), "significant_ref": len(r_sig), "significant_equal_ref": n_sig == r_sig,
                "significant_only_fixed": len(n_sig - r_sig), "significant_only_ref": len(r_sig - n_sig),
                "leading_edge_equal_ref": sum(as_set(n[t]["core_enrichment"]) == as_set(ref[t]["leadingEdge"]) for t in common),
                "leading_edge_equal_ref_negative_es": sum(as_set(n[t]["core_enrichment"]) == as_set(ref[t]["leadingEdge"]) for t in neg), "negative_es_terms": len(neg),
            }
            if old_l2 is not None:
                o = gsea_table(old_l2, cid, coll)
                both = sorted(set(o) & set(n))
                o_sig = {t for t, r in o.items() if r["significant"] == "TRUE"}
                nes_diff = [abs(num(n[t]["NES"]) - num(o[t]["NES"])) for t in both]
                p_ratio = [abs(math.log10(num(n[t]["pvalue"]) / num(o[t]["pvalue"]))) for t in both if num(o[t]["pvalue"]) > 0 and num(n[t]["pvalue"]) > 0]
                row.update({
                    "terms_v160": len(o), "es_identical_v160": all(n[t]["enrichmentScore"] == o[t]["enrichmentScore"] for t in both),
                    "nes_changed_v160": sum(n[t]["NES"] != o[t]["NES"] for t in both), "nes_abs_change_median": statistics.median(nes_diff), "nes_abs_change_max": max(nes_diff),
                    "pvalue_changed_v160": sum(n[t]["pvalue"] != o[t]["pvalue"] for t in both), "pvalue_abs_log10_ratio_median": statistics.median(p_ratio), "pvalue_abs_log10_ratio_max": max(p_ratio),
                    "padj_changed_v160": sum(n[t]["p.adjust"] != o[t]["p.adjust"] for t in both),
                    "significant_v160": len(o_sig), "significant_gained": len(n_sig - o_sig), "significant_lost": len(o_sig - n_sig),
                    "core_changed_v160": sum(as_set(n[t]["core_enrichment"]) != as_set(o[t]["core_enrichment"]) for t in both),
                })
            gsea_rows.append(row)
write("gsea_fixed_vs_reference_and_v160.tsv", gsea_rows)


# ---------------------------------------------------------------- reproducibility
def digests(l2: Path) -> dict[str, str]:
    return {str(p.relative_to(l2)): hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(l2.rglob("*.tsv"))}


rep_rows = []
a, b = digests(l2_dir(W, "D1", "FIX-D1")), digests(l2_dir(W, "D1", "FIX-D1-REP"))
rep_rows.append({"comparison": "FIX-D1 vs FIX-D1-REP (same seed 1)", "tsv_files": len(a), "identical": sum(a[k] == b.get(k) for k in a),
                 "differing": ";".join(k for k in a if a[k] != b.get(k))})
s42 = l2_dir(W, "D1s42", "FIX-D1-S42")
c = digests(s42)
gsea_keys = [k for k in a if k.startswith("enrichment/gsea_")]
rep_rows.append({"comparison": "FIX-D1 (seed 1) vs FIX-D1-S42 (seed 42), GSEA tables", "tsv_files": len(gsea_keys), "identical": sum(a[k] == c.get(k) for k in gsea_keys),
                 "differing": ";".join(k for k in gsea_keys if a[k] != c.get(k))})
non_gsea = [k for k in a if not k.startswith("enrichment/gsea_")]
rep_rows.append({"comparison": "FIX-D1 vs FIX-D1-S42, non-GSEA tables (seed must not matter)", "tsv_files": len(non_gsea), "identical": sum(a[k] == c.get(k) for k in non_gsea),
                 "differing": ";".join(k for k in non_gsea if a[k] != c.get(k))})
seed_rows = []
for coll in GSEA_COLLECTIONS:
    one, other = gsea_table(l2_dir(W, "D1", "FIX-D1"), "trt_vs_untrt", coll), gsea_table(s42, "trt_vs_untrt", coll)
    both = sorted(set(one) & set(other))
    seed_rows.append({"collection": coll, "terms": len(both), "es_identical": all(one[t]["enrichmentScore"] == other[t]["enrichmentScore"] for t in both),
                      "pvalue_differs": sum(one[t]["pvalue"] != other[t]["pvalue"] for t in both), "nes_differs": sum(one[t]["NES"] != other[t]["NES"] for t in both),
                      "significant_seed1": sum(one[t]["significant"] == "TRUE" for t in both), "significant_seed42": sum(other[t]["significant"] == "TRUE" for t in both),
                      "significance_differs": sum(one[t]["significant"] != other[t]["significant"] for t in both)})
write("reproducibility.tsv", rep_rows)
write("seed_1_vs_42.tsv", seed_rows)
print("written to", OUT)
