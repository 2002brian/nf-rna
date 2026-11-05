"""Independent regression gate for GO/KEGG ORA and preranked GSEA.

One simulated human (Ensembl) two-group project with planted effects on two GO
BP terms (up: DNA replication, down: inflammatory response) and a matching
offline KEGG snapshot.  The panel deliberately contains version-suffixed IDs,
unmapped IDs, Ensembl IDs with several Entrez targets (one-to-many) and Entrez
IDs reached from several Ensembl IDs (many-to-one).

Expectations come from ``oracles/annotation_reference.R`` (AnnotationDbi
queried directly) and from closed-form statistics in ``_science`` -- not from
nf-rna's R scripts:

* ORA universe = padj-non-NA genes -> Entrez -> genes annotated in the ontology
  (or the KEGG snapshot), exact;
* every ORA p-value is the hypergeometric upper tail of its own GeneRatio/BgRatio,
  and p.adjust is Benjamini-Hochberg over the evaluated terms (rel 1e-9);
* the ranked list follows the documented mapping policy, exact;
* every GSEA setSize and enrichment score is recomputed from the ranked list and
  the independently queried gene sets (abs 1e-9), as is the leading edge;
* fgsea p-values/NES are permutation-based: they are only required to be
  reproducible for a fixed seed, and the planted sets to be significant with the
  right sign.
"""

from __future__ import annotations

import hashlib
import json
import math
import random
from pathlib import Path

import pytest

from _science import (
    DOWN_TERM, KEGG_DOWN, KEGG_UP, UP_TERM, benjamini_hochberg, by_gene, gsea_enrichment_score, hypergeometric_upper,
    read_tsv, require_downstream_r_runtime, run_enrichment, run_l2, run_reference, write_kegg_snapshot, write_project,
)
from test_deseq2_designs import CONTRAST, _two_group

SNAPSHOT = "hsa_synthetic_20260101T000000Z"
ONTOLOGIES = ("BP", "MF", "CC")


class Enriched:
    pass


def _strip(gene_id: str) -> str:
    head, dot, tail = gene_id.rpartition(".")
    return head if dot and tail.isdigit() else gene_id


@pytest.fixture(scope="module")
def run(tmp_path_factory) -> Enriched:
    require_downstream_r_runtime("jsonlite", "DESeq2", "apeglm", "ggplot2", "pheatmap", "clusterProfiler", "org.Hs.eg.db", "GO.db", "fgsea", "openssl")
    work = tmp_path_factory.mktemp("enrichment")
    # 450 low-expression annotated null genes (mean 1-4): independent filtering gives many of them
    # NA padj, so the padj-non-NA ORA universe is a strict subset of the retained genes (checked below).
    design, sim, roles = _two_group(low_expression=450)
    ids = work / "ids.txt"
    ids.write_text("\n".join(sim.genes) + "\n", encoding="utf-8")
    (work / "none.txt").write_text("", encoding="utf-8")
    run_reference("annotation_reference.R", ids, work / "none.txt", work / "ref")
    mapping: dict[str, list[str]] = {}
    for row in read_tsv(work / "ref_mapping.tsv"):
        mapping.setdefault(row["ENSEMBL"], []).append(row["ENTREZID"])
    entrez_of = lambda role: sorted({e for g in sim.genes if roles[g] == role for e in mapping.get(_strip(g), [])}, key=int)  # noqa: E731
    rng = random.Random(5)
    null = entrez_of("null")
    pathways = {KEGG_UP: entrez_of("up_term"), KEGG_DOWN: entrez_of("down_term")}
    pathways |= {f"hsa990{index:02d}": rng.sample(null, rng.randint(12, 40)) for index in range(10, 22)}
    names = {term: f"Synthetic pathway {term}" for term in pathways}
    monkeypatch = pytest.MonkeyPatch()
    monkeypatch.setenv("XDG_CACHE_HOME", str(work / "cache"))
    write_kegg_snapshot(work / "cache", SNAPSHOT, pathways, names)
    root = write_project(work / "two_group_enrichment", sim, design, enrichment=["go", "kegg", "gsea"], human_annotation=True, kegg_snapshot=SNAPSHOT)
    try:
        l2 = run_l2(root)
        run_enrichment(root)
        state = Enriched()
        state.root, state.l2, state.sim, state.roles, state.mapping, state.pathways = root, l2, sim, roles, mapping, pathways
        state.all = by_gene(read_tsv(l2 / "contrasts" / CONTRAST / "all_genes.tsv"))
        state.enrichment = l2 / "enrichment"
        terms = set()
        for ontology in ONTOLOGIES:
            terms |= {row["ID"] for row in read_tsv(state.enrichment / "gsea_go" / CONTRAST / ontology / "all_terms.tsv")}
            for foreground in ("significant", "up", "down"):
                path = state.enrichment / "go" / CONTRAST / foreground / f"GO_{ontology}_all_terms.tsv"
                if path.is_file():
                    terms |= {row["ID"] for row in read_tsv(path)}
        (work / "terms.txt").write_text("\n".join(sorted(terms)) + "\n", encoding="utf-8")
        run_reference("annotation_reference.R", ids, work / "terms.txt", work / "ref")
        state.members: dict[str, set[str]] = {}
        for row in read_tsv(work / "ref_term_members.tsv"):
            state.members.setdefault(row["GOALL"], set()).add(row["ENTREZID"])
        state.annotated = {o: set((work / f"ref_annotated_{o}.txt").read_text(encoding="utf-8").split()) for o in ONTOLOGIES}
        state.monkeypatch = monkeypatch
        yield state
    finally:
        monkeypatch.undo()


def _tested_entrez(run) -> set[str]:
    return {e for gene, row in run.all.items() if row["padj"] != "NA" for e in run.mapping.get(_strip(gene), [])}


def _foreground_entrez(run, name: str) -> set[str]:
    rows = read_tsv(run.l2 / "contrasts" / CONTRAST / {"significant": "significant.tsv", "up": "upregulated.tsv", "down": "downregulated.tsv"}[name])
    return {e for row in rows for e in run.mapping.get(_strip(row["gene_id"]), [])}


def test_fixture_exercises_every_identifier_case(run):
    roles = set(run.roles.values())
    assert {"one_to_many", "many_to_one", "unmapped"} <= roles
    assert any("." in gene for gene in run.sim.genes)
    assert all(len(run.mapping[_strip(g)]) >= 2 for g in run.sim.genes if run.roles[g] == "one_to_many")
    assert all(_strip(g) not in run.mapping for g in run.sim.genes if run.roles[g] == "unmapped")


def _retained_entrez(run) -> set[str]:
    return {e for gene in run.all for e in run.mapping.get(_strip(gene), [])}


@pytest.mark.parametrize("ontology", ONTOLOGIES)
def test_fixture_separates_the_padj_non_na_universe_from_all_retained_genes(run, ontology):
    # Precondition for the universe checks: genes removed by independent filtering are
    # annotated, so a universe built from all retained genes would differ.
    assert sum(row["padj"] == "NA" for row in run.all.values()) >= 50
    assert _tested_entrez(run) & run.annotated[ontology] < _retained_entrez(run) & run.annotated[ontology]


@pytest.mark.parametrize("ontology", ONTOLOGIES)
@pytest.mark.parametrize("foreground", ["up", "down"])
def test_go_ora_universe_and_hypergeometric_statistics(run, ontology, foreground):
    summary = json.loads((run.enrichment / "go" / "go_backend_summary.json").read_text(encoding="utf-8"))
    record = summary["contrasts"][0]["foregrounds"][foreground]["ontologies"][ontology]
    universe = _tested_entrez(run) & run.annotated[ontology]
    assert record["universe"]["universe_genes"] == len(universe)
    foreground_entrez = _foreground_entrez(run, foreground)
    assert record["universe"]["foreground_in_universe"] == len(foreground_entrez & universe)
    terms = read_tsv(run.enrichment / "go" / CONTRAST / foreground / f"GO_{ontology}_all_terms.tsv")
    assert terms, "the planted effects must give evaluated terms"
    _check_ora_terms(terms, universe, foreground_entrez, run.members, max_size=500)


def test_go_ora_finds_each_planted_term_only_in_its_own_direction(run):
    up = by_gene(read_tsv(run.enrichment / "go" / CONTRAST / "up" / "GO_BP_all_terms.tsv"), "ID")
    down = by_gene(read_tsv(run.enrichment / "go" / CONTRAST / "down" / "GO_BP_all_terms.tsv"), "ID")
    assert up[UP_TERM]["significant"] == "TRUE" and down[DOWN_TERM]["significant"] == "TRUE"
    assert up.get(DOWN_TERM, {}).get("significant", "FALSE") == "FALSE"
    assert down.get(UP_TERM, {}).get("significant", "FALSE") == "FALSE"


@pytest.mark.parametrize("foreground", ["up", "down"])
def test_kegg_ora_universe_and_statistics_from_the_offline_snapshot(run, foreground):
    summary = json.loads((run.enrichment / "kegg" / "kegg_backend_summary.json").read_text(encoding="utf-8"))
    assert summary["resource"]["snapshot_id"] == SNAPSHOT and summary["resource"]["sha256_verified"] is True
    # Offline: the pinned local snapshot was used and nothing was downloaded from KEGG REST.
    assert summary["resource"]["selection"].startswith("configured") and summary["resource"]["created_by_this_run"] is False
    snapshot_genes = {gene for genes in run.pathways.values() for gene in genes}
    universe = _tested_entrez(run) & snapshot_genes
    assert summary["contrasts"][0]["universe"]["kegg"]["universe_genes"] == len(universe)
    terms = read_tsv(run.enrichment / "kegg" / CONTRAST / foreground / "kegg_all_terms.tsv")
    _check_ora_terms(terms, universe, _foreground_entrez(run, foreground), {k: set(v) for k, v in run.pathways.items()}, max_size=500)
    planted = KEGG_UP if foreground == "up" else KEGG_DOWN
    assert by_gene(terms, "ID")[planted]["significant"] == "TRUE"


def _check_ora_terms(terms, universe, foreground, members, *, max_size):
    pvalues = []
    for row in terms:
        k, n = map(int, row["GeneRatio"].split("/"))
        m, total = map(int, row["BgRatio"].split("/"))
        assert total == len(universe), row["ID"]
        assert n == len(foreground & universe), row["ID"]
        assert m == len(members[row["ID"]] & universe), row["ID"]
        hits = set(row["geneID"].split("/"))
        assert len(hits) == k == int(row["Count"]) and hits == members[row["ID"]] & universe & foreground, row["ID"]
        assert 10 <= m <= max_size
        expected = hypergeometric_upper(k, total, m, n)
        assert math.isclose(float(row["pvalue"]), expected, rel_tol=1e-9, abs_tol=1e-300), (row["ID"], row["pvalue"], expected)
        pvalues.append(float(row["pvalue"]))
    for row, adjusted in zip(terms, benjamini_hochberg(pvalues)):
        assert math.isclose(float(row["p.adjust"]), adjusted, rel_tol=1e-9, abs_tol=1e-300), row["ID"]
        assert (row["significant"] == "TRUE") == (float(row["p.adjust"]) <= 0.05)


def _expected_ranking(run) -> list[tuple[str, float, str]]:
    """nf-rna's documented GSEA ranking policy, implemented independently.

    Every finite-stat source gene contributes one row per Entrez target
    (one-to-many targets propagate the same statistic); a target reached from several
    sources keeps the largest |stat| (ties: original gene_id, then normalized ID);
    the list is sorted by decreasing stat, ties by Entrez ID as text.
    """

    candidates: dict[str, list[tuple[float, str]]] = {}
    for gene, row in run.all.items():
        if row["stat"] in {"NA", ""}:
            continue
        for entrez in run.mapping.get(_strip(gene), []):
            candidates.setdefault(entrez, []).append((float(row["stat"]), gene))
    retained = []
    for entrez, rows in candidates.items():
        stat, gene = min(rows, key=lambda item: (-abs(item[0]), item[1], _strip(item[1])))
        retained.append((entrez, stat, gene))
    return sorted(retained, key=lambda item: (-item[1], item[0]))


@pytest.mark.parametrize("backend", ["gsea_go", "gsea_kegg"])
def test_gsea_ranked_list_follows_the_mapping_policy_exactly(run, backend):
    ranked = read_tsv(run.enrichment / backend / CONTRAST / "ranked_gene_list.tsv")
    observed = [(row["entrez_id"], float(row["stat"]), row["original_gene_id"]) for row in ranked]
    assert observed == _expected_ranking(run)
    assert [int(row["rank"]) for row in ranked] == list(range(1, len(ranked) + 1))


def test_gsea_ranking_summary_counts_are_reconstructible(run):
    summary = json.loads((run.enrichment / "gsea_go" / CONTRAST / "gsea_ranking_summary.json").read_text(encoding="utf-8"))
    finite = [gene for gene, row in run.all.items() if row["stat"] not in {"NA", ""}]
    mapped = [gene for gene in finite if _strip(gene) in run.mapping]
    rows = sum(len(run.mapping[_strip(gene)]) for gene in mapped)
    expected = _expected_ranking(run)
    assert summary["finite_stat_source_genes"] == len(finite)
    assert summary["mapped_source_genes"] == len(mapped)
    assert summary["final_ranked_genes"] == len(expected) == rows - summary["duplicate_target_rows_collapsed"]
    assert summary["one_to_many_source_ids"] == sum(len(run.mapping[_strip(g)]) > 1 for g in mapped) > 0
    assert summary["duplicate_target_rows_collapsed"] > 0


def _running_extremes(ranked: list[tuple[str, float]], members: set[str]) -> tuple[list[str], list[str]]:
    """Leading edge by the running-sum definition: hits up to the maximum (ES > 0) or from the minimum (ES < 0)."""

    hits = [index for index, (gene, _stat) in enumerate(ranked) if gene in members]
    total = sum(abs(ranked[i][1]) for i in hits)
    miss = 1.0 / (len(ranked) - len(hits))
    running, top, bottom, seen = 0.0, [], [], 0
    for index, (gene, stat) in enumerate(ranked):
        if gene in members:
            bottom.append(running - miss * 0)  # value just before this hit (after preceding misses)
            running += abs(stat) / total
            top.append(running)
            seen += 1
        else:
            running -= miss
    peak = max(range(len(top)), key=lambda i: (top[i], -i))
    trough = min(range(len(bottom)), key=lambda i: (bottom[i], i))
    names = [ranked[i][0] for i in hits]
    return names[: peak + 1], names[trough:]


@pytest.mark.parametrize("ontology", ONTOLOGIES)
def test_go_gsea_set_sizes_scores_and_leading_edges_are_recomputed_exactly(run, ontology):
    ranked = [(row["entrez_id"], float(row["stat"])) for row in read_tsv(run.enrichment / "gsea_go" / CONTRAST / "ranked_gene_list.tsv")]
    ranked_ids = {gene for gene, _ in ranked}
    terms = read_tsv(run.enrichment / "gsea_go" / CONTRAST / ontology / "all_terms.tsv")
    assert terms
    for row in terms:
        members = run.members[row["ID"]] & ranked_ids
        _check_gsea_term(row, ranked, members)


def test_kegg_gsea_set_sizes_scores_and_leading_edges_are_recomputed_exactly(run):
    ranked = [(row["entrez_id"], float(row["stat"])) for row in read_tsv(run.enrichment / "gsea_kegg" / CONTRAST / "ranked_gene_list.tsv")]
    ranked_ids = {gene for gene, _ in ranked}
    terms = read_tsv(run.enrichment / "gsea_kegg" / CONTRAST / "all_terms.tsv")
    assert {row["ID"] for row in terms} == {term for term, genes in run.pathways.items() if 10 <= len(set(genes) & ranked_ids) <= 500}
    for row in terms:
        _check_gsea_term(row, ranked, set(run.pathways[row["ID"]]) & ranked_ids)


def _check_gsea_term(row, ranked, members):
    assert int(row["setSize"]) == len(members), row["ID"]
    assert math.isclose(float(row["enrichmentScore"]), gsea_enrichment_score(ranked, members), abs_tol=1e-9), row["ID"]
    positive, negative = _running_extremes(ranked, members)
    core = row["core_enrichment"].split("/")
    assert core == (positive if float(row["enrichmentScore"]) > 0 else negative), row["ID"]
    assert (row["significant"] == "TRUE") == (float(row["p.adjust"]) <= 0.05)


def test_planted_gene_sets_have_the_correct_nes_sign_and_significance(run):
    go = by_gene(read_tsv(run.enrichment / "gsea_go" / CONTRAST / "BP" / "all_terms.tsv"), "ID")
    kegg = by_gene(read_tsv(run.enrichment / "gsea_kegg" / CONTRAST / "all_terms.tsv"), "ID")
    for table, up, down in ((go, UP_TERM, DOWN_TERM), (kegg, KEGG_UP, KEGG_DOWN)):
        assert float(table[up]["NES"]) > 0 and table[up]["significant"] == "TRUE"
        assert float(table[down]["NES"]) < 0 and table[down]["significant"] == "TRUE"
    for directory in (run.enrichment / "gsea_go" / CONTRAST / "BP", run.enrichment / "gsea_kegg" / CONTRAST):
        significant = {row["ID"] for row in read_tsv(directory / "significant.tsv")}
        positive = {row["ID"] for row in read_tsv(directory / "positive_enrichment.tsv")}
        negative = {row["ID"] for row in read_tsv(directory / "negative_enrichment.tsv")}
        assert positive | negative == significant and not positive & negative


def test_gsea_is_reproducible_for_a_fixed_seed(run):
    from rnaseq.gsea import execute_gsea, prepare_gsea
    from rnaseq.kegg import execute_kegg, prepare_kegg
    from rnaseq.validators import validate_project

    tables = [run.enrichment / "gsea_go" / CONTRAST / o / "all_terms.tsv" for o in ONTOLOGIES] + [run.enrichment / "gsea_kegg" / CONTRAST / "all_terms.tsv"]
    before = [hashlib.sha256(path.read_bytes()).hexdigest() for path in tables]
    report = validate_project(run.root)
    execute_gsea(prepare_gsea(report, run_id=None))
    execute_kegg(prepare_kegg(report, run_id=None, mode="gsea"))
    assert [hashlib.sha256(path.read_bytes()).hexdigest() for path in tables] == before
