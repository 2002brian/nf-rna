"""Regression tests for the v1.6.1 scientific corrections (docs/dev/scientific_corrections_v161.md).

D-1: GSEA leading edge by the standard running-sum definition (DOSE 4.4.0 is wrong for ES < 0).
D-2: the configured GSEA seed controls fgsea (DOSE 4.4.0 seed = TRUE ignores it).
M-1: DESeq2 results(alpha) equals the declared DEG FDR threshold.
"""
from __future__ import annotations

import subprocess
from pathlib import Path

import pytest

from conftest import base_config, require_r_packages
from rnaseq.validators import validate_project

R_ROOT = Path(__file__).parents[1] / "src" / "rnaseq" / "r"
HELPER = R_ROOT / "gsea_leading_edge.R"

# A 14-gene list where the running-sum minimum (-0.472, just before g14) lies after
# hit g9, whose post-step value (-0.307) is the lowest of the hits: DOSE 4.4.0
# reports g9/g11/g14, the standard leading edge is g14 alone.
DOSE_COUNTEREXAMPLE = '''
stats <- c(2.8, 1.7, 1.6, 1.3, -0.3, -0.5, -0.6, -0.7, -1.2, -1.4, -1.5, -2.5, -3.1, -4.2); names(stats) <- paste0("g", 1:14)
members <- c("g2", "g5", "g9", "g11", "g14")
'''

# 400 genes, 40 gene sets of 15-60 genes; the first ten are shifted so some sets are enriched.
SYNTHETIC_GSEA = '''
set.seed(11)
genes <- as.character(1:400)
stats <- sort(setNames(rnorm(400), genes), decreasing=TRUE)
sets <- lapply(1:40, function(i) sample(genes, sample(15:60, 1)))
for (i in 1:5) sets[[i]] <- unique(c(sets[[i]], names(stats)[1:20 + i]))
for (i in 6:10) sets[[i]] <- unique(c(sets[[i]], names(stats)[380:400 - i]))
names(sets) <- sprintf("SET%02d", 1:40)
term2gene <- data.frame(term=rep(names(sets), lengths(sets)), gene=unlist(sets, use.names=FALSE))
run_gsea <- function(seed_arg) suppressMessages(clusterProfiler::GSEA(geneList=stats, TERM2GENE=term2gene, minGSSize=10, maxGSSize=500, eps=0,
  pvalueCutoff=1, pAdjustMethod="BH", verbose=FALSE, seed=seed_arg, nPermSimple=1000))
'''


def _rscript(code: str) -> None:
    result = subprocess.run(["Rscript", "-e", code], capture_output=True, text=True, check=False)
    assert result.returncode == 0, result.stdout + result.stderr


def test_leading_edge_corrects_the_dose_negative_es_counterexample():
    require_r_packages("DOSE", "fgsea")
    _rscript(f'''suppressPackageStartupMessages({{ library(DOSE); library(fgsea) }}); source("{HELPER}")
{DOSE_COUNTEREXAMPLE}
es <- calcGseaStat(stats, match(members, names(stats)))
standard <- nf_rna_leading_edge(stats, members, es)
fg <- names(stats)[calcGseaStat(stats, match(members, names(stats)), returnLeadingEdge=TRUE)$leadingEdge]
dose <- DOSE:::leading_edge(list(DOSE:::gseaScores(stats, members, exponent=1)))
stopifnot(es < 0, identical(standard$core, "g14"), setequal(standard$core, fg), abs(standard$extremum - es) < 1e-12,
          identical(dose$core_enrichment[[1]], c("g9", "g11", "g14")),  # the locked DOSE defect this replaces
          standard$rank == dose$rank[[1]])  # DOSE's rank already uses the full running sum
''')


def test_leading_edge_matches_fgsea_for_both_signs_ties_and_absent_members():
    require_r_packages("fgsea")
    # Integer statistics with many ties give tied running-sum extrema; fgsea resolves
    # them to the first extremum, as nf-rna does.  Members absent from the list are ignored.
    _rscript(f'''suppressPackageStartupMessages(library(fgsea)); source("{HELPER}")
set.seed(5); counts <- c(positive=0L, negative=0L, tied=0L)
for (trial in 1:3000) {{
  stats <- sort(sample(c(-4:-1, 1:4), 30, replace=TRUE), decreasing=TRUE); names(stats) <- sprintf("g%02d", 1:30)
  members <- c(sample(names(stats), sample(3:12, 1)), "absent")
  hits <- sort(match(members[members %in% names(stats)], names(stats)))
  es <- calcGseaStat(stats, hits)
  if (es == 0) next
  fg <- names(stats)[calcGseaStat(stats, hits, returnLeadingEdge=TRUE)$leadingEdge]
  le <- nf_rna_leading_edge(stats, members, es)
  if (!setequal(le$core, fg) || abs(le$extremum - es) > 1e-12 || le$set_size != length(hits)) stop(sprintf("trial %d disagrees with fgsea", trial))
  # Positions are those of the ranked list: the core is reported in list order.
  if (!identical(le$core, names(stats)[sort(match(le$core, names(stats)))])) stop("core not in rank order")
  w <- abs(stats[hits]); top <- cumsum(w) / sum(w) - (hits - seq_along(hits)) / (30 - length(hits)); bottom <- top - w / sum(w)
  extreme <- if (es > 0) top else bottom
  if (sum(abs(extreme - es) < 1e-12) > 1) counts[["tied"]] <- counts[["tied"]] + 1L
  counts[[if (es > 0) "positive" else "negative"]] <- counts[[if (es > 0) "positive" else "negative"]] + 1L
}}
stopifnot(counts[["positive"]] > 500, counts[["negative"]] > 500, counts[["tied"]] > 20)
''')


def test_corrected_gsea_result_keeps_statistics_and_replaces_only_leading_edge_fields():
    require_r_packages("clusterProfiler", "DOSE", "fgsea")
    _rscript(f'''suppressPackageStartupMessages({{ library(clusterProfiler); library(fgsea) }}); source("{HELPER}")
{SYNTHETIC_GSEA}
set.seed(1); before <- run_gsea(FALSE)
corrected <- nf_rna_correct_leading_edge(before); after <- corrected$result
a <- before@result; b <- after@result
stats_cols <- c("ID", "Description", "setSize", "enrichmentScore", "NES", "pvalue", "p.adjust", "qvalue")
stopifnot(identical(a[, stats_cols], b[, stats_cols]), identical(names(a), names(b)), identical(before@geneList, after@geneList))
set.seed(1); direct <- fgsea(pathways=sets, stats=stats, minSize=10, maxSize=500, eps=0, nPermSimple=1000, gseaParam=1, nproc=0)
direct <- direct[match(b$ID, direct$pathway), ]
same <- mapply(function(x, y) setequal(strsplit(x, "/", fixed=TRUE)[[1]], y), b$core_enrichment, direct$leadingEdge)
stopifnot(all(same), corrected$audit$terms == nrow(b), corrected$audit$core_enrichment_changed == sum(a$core_enrichment != b$core_enrichment),
          corrected$audit$core_enrichment_changed_positive_es == 0L, corrected$audit$max_abs_es_difference < 1e-12,
          all(b$core_enrichment[b$enrichmentScore >= 0] == a$core_enrichment[a$enrichmentScore >= 0]),
          all(b$rank == a$rank), all(grepl("^tags=[0-9]+%, list=[0-9]+%, signal=[0-9]+%$", b$leading_edge)))
# Misaligned statistics or membership must fail loudly, never produce a silent core set.
broken <- before; broken@result$enrichmentScore[[1]] <- broken@result$enrichmentScore[[1]] + 0.01
stopifnot(inherits(try(nf_rna_correct_leading_edge(broken), silent=TRUE), "try-error"))
broken <- before; broken@result$setSize[[1]] <- broken@result$setSize[[1]] + 1L
stopifnot(inherits(try(nf_rna_correct_leading_edge(broken), silent=TRUE), "try-error"))
''')


def test_configured_seed_controls_fgsea_and_is_reproducible():
    require_r_packages("clusterProfiler", "DOSE", "fgsea")
    _rscript(f'''suppressPackageStartupMessages({{ library(clusterProfiler); library(fgsea) }})
{SYNTHETIC_GSEA}
numbers <- function(x) x@result[order(x@result$ID), c("ID", "enrichmentScore", "NES", "pvalue", "p.adjust")]
set.seed(7); first <- run_gsea(FALSE); set.seed(7); again <- run_gsea(FALSE); set.seed(8); other <- run_gsea(FALSE)
stopifnot(identical(numbers(first), numbers(again)))  # same seed, same result
stopifnot(!identical(numbers(first)$pvalue, numbers(other)$pvalue))  # another seed reaches fgsea (not every set need differ)
stopifnot(identical(numbers(first)$enrichmentScore, numbers(other)$enrichmentScore))  # ES is deterministic
set.seed(7); direct <- fgsea(pathways=sets, stats=stats, minSize=10, maxSize=500, eps=0, nPermSimple=1000, gseaParam=1, nproc=0)
direct <- direct[match(numbers(first)$ID, direct$pathway), ]
stopifnot(identical(numbers(first)$pvalue, direct$pval), identical(numbers(first)$NES, direct$NES))  # exactly fgsea after set.seed(seed)
# Locked DOSE 4.4.0 with seed = TRUE re-seeds with set.seed(.Random.seed): every configured seed gives one fixed draw.
set.seed(7); legacy7 <- run_gsea(TRUE); set.seed(8); legacy8 <- run_gsea(TRUE)
set.seed(10403); fixed <- fgsea(pathways=sets, stats=stats, minSize=10, maxSize=500, eps=0, nPermSimple=1000, gseaParam=1, nproc=0)
stopifnot(identical(numbers(legacy7), numbers(legacy8)), identical(numbers(legacy7)$pvalue, fixed[match(numbers(legacy7)$ID, fixed$pathway), ]$pval))
''')


def test_gsea_backends_pass_seed_false_and_recompute_the_leading_edge():
    common = (R_ROOT / "enrichment_common.R").read_text(encoding="utf-8")
    for name, call in (("gsea_analysis.R", "gseGO("), ("kegg_analysis.R", "GSEA(")):
        script = (R_ROOT / name).read_text(encoding="utf-8")
        assert 'source(file.path(dirname(normalizePath(script_file)), "gsea_leading_edge.R"))' in script
        assert call in script and "seed=FALSE" in script and "seed=TRUE" not in script
        assert script.index("set.seed(as.integer(gsea_cfg$seed))") < script.index(call)
        # The correction runs on the gseaResult before any table, member list or figure is written.
        assert script.index("nf_rna_correct_leading_edge(result)") < script.index("gsea_term_tables(raw")
        assert "leading_edge=corrected$audit" in script
    assert "seed = FALSE" in common and "leading_edge=NF_RNA_LEADING_EDGE_DEFINITION" in common


def test_deseq2_alpha_is_the_declared_padj_threshold():
    script = (R_ROOT / "l2_analysis.R").read_text(encoding="utf-8")
    assert "results(dds, contrast = c(factor, numerator, denominator), independentFiltering = TRUE, alpha = cfg$thresholds$padj)" in script
    assert "independent_filtering_alpha = cfg$thresholds$padj" in script
    assert "filterThreshold" in script


def test_deseq2_alpha_changes_independent_filtering_as_documented():
    require_r_packages("DESeq2")
    # DESeq2 chooses the filter that maximises rejections at padj < alpha; the chosen
    # alpha is recorded in the result metadata and only padj (not pvalue) may change.
    _rscript('''suppressPackageStartupMessages(library(DESeq2))
set.seed(2); dds <- makeExampleDESeqDataSet(n=3000, m=8, betaSD=1); dds <- DESeq(dds, quiet=TRUE)
r10 <- results(dds, independentFiltering=TRUE); r05 <- results(dds, independentFiltering=TRUE, alpha=0.05)
stopifnot(metadata(r10)$alpha == 0.1, metadata(r05)$alpha == 0.05, identical(r10$pvalue, r05$pvalue), identical(r10$stat, r05$stat),
          is.numeric(metadata(r05)$filterThreshold), all(is.na(r05$padj) == (r05$baseMean <= metadata(r05)$filterThreshold | is.na(r05$pvalue))))
''')


@pytest.mark.parametrize("padj", [0.0, 1.0])
def test_padj_threshold_must_be_a_valid_deseq2_alpha(project_factory, padj):
    config = base_config()
    config["thresholds"]["padj"] = padj
    report = validate_project(project_factory(config=config))
    assert not report.is_valid
