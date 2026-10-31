# U-1 minimal reproduction, D2 LP_vs_ML.  Runs fgsea (seed 1, the reference parameters) on
# A = the reference ranked list (17 significant digits) and B = the nf-rna ranked list (the
# DESeq2 stat as written to all_genes.tsv, 15 significant digits), and on A with the pathway
# order shuffled, with all collection sets, and with MulticoreParam.
# Usage: Rscript u1_repro.R <GO_BP|GO_MF> <BP|MF>
suppressPackageStartupMessages({library(fgsea); library(BiocParallel); library(data.table)})
V <- "/home/brian/nf-rna-validation"
args <- commandArgs(TRUE); coll <- args[1]; sub <- args[2]
l2 <- Sys.glob(file.path(V, "v161/nfrna/D2/runs/FIX-D2/*/downstream/l2"))
rd <- read.delim(file.path(V, "v1.6.0/reference/outputs/D2/contrasts/LP_vs_ML/gsea/ranked_list.tsv"), colClasses=c(entrez_id="character"))
statsA <- setNames(rd$stat, rd$entrez_id)                       # reference, full precision
nd <- read.delim(file.path(l2, "enrichment/gsea_go/LP_vs_ML/ranked_gene_list.tsv"), colClasses=c(entrez_id="character"))
statsB <- setNames(nd$stat, nd$entrez_id)                       # nf-rna, 15 significant digits
stopifnot(identical(names(statsA), names(statsB)))
cat("stats: n =", length(statsA), " max rel diff A vs B =", max(abs(statsA-statsB)/abs(statsA)), "\n")
cat("exact ties in A:", sum(duplicated(statsA)), " in B:", sum(duplicated(statsB)), "\n")
cl <- readRDS(file.path(V, "v1.6.0/reference/work/collections/mmu.rds"))[[coll]]
sets <- lapply(cl$sets, function(s) s[s %in% names(statsA)])
pw <- sets[lengths(sets) >= 10 & lengths(sets) <= 500]
ref <- fread(file.path(V, "v1.6.0/reference/outputs/D2/contrasts/LP_vs_ML/gsea", coll, "seed_1.tsv"), colClasses=c(set_id="character"))
nf <- fread(file.path(l2, "enrichment/gsea_go/LP_vs_ML", sub, "all_terms.tsv"), colClasses=c(ID="character"))
run <- function(stats, pathways, bp=SerialParam()) { set.seed(1); r <- fgsea(pathways=pathways, stats=stats, minSize=10, maxSize=500, eps=0, nPermSimple=10000, gseaParam=1, scoreType="std", BPPARAM=bp); r[match(names(pw), r$pathway)] }
cmp <- function(r, tgt_nes, tgt_p) c(nes_diff=sum(abs(r$NES/tgt_nes-1) > 1e-12), p_diff=sum(abs(r$pval/tgt_p-1) > 1e-12))
o_ref <- list(nes=ref$NES[match(names(pw), ref$set_id)], p=ref$pval[match(names(pw), ref$set_id)])
o_nf <- list(nes=nf$NES[match(names(pw), nf$ID)], p=nf$pvalue[match(names(pw), nf$ID)])
report <- function(label, r) cat(sprintf("%-42s vs reference: NES %4d p %4d | vs nf-rna: NES %4d p %4d\n", label,
  cmp(r, o_ref$nes, o_ref$p)[1], cmp(r, o_ref$nes, o_ref$p)[2], cmp(r, o_nf$nes, o_nf$p)[1], cmp(r, o_nf$nes, o_nf$p)[2]))
report("A: reference stats, serial", rA <- run(statsA, pw))
report("B: nf-rna stats (15 digits), serial", run(statsB, pw))
set.seed(99); report("A, pathway order shuffled", run(statsA, pw[sample(length(pw))]))
report("A, all collection sets (unfiltered list)", run(statsA, sets))
report("A, MulticoreParam(4)", run(statsA, pw, MulticoreParam(4)))
saveRDS(list(statsA=statsA, statsB=statsB, pw=pw), paste0("inputs_", coll, ".rds"))
