# GO_MF: locate the permutation whose random ES is compared with the observed ES of GO:0051880
suppressPackageStartupMessages(library(fgsea))
inp <- readRDS(commandArgs(TRUE)[1])
prep <- function(stats) fgsea:::preparePathwaysAndStats(inp$pw, stats, 10, 500, 1, "std")
pa <- prep(inp$statsA); pb <- prep(inp$statsB)
set.seed(1); seed <- sample.int(10^9, 1)
es <- function(pp) calcGseaStat(pp$stats, pp$filtered[["GO:0051880"]], gseaParam = 1)
ea <- es(pa); eb <- es(pb)
sizes <- c(length(pa$filtered[["GO:0051880"]]), max(pa$sizes))
batch <- function(pp, e, iters) fgsea:::calcGseaStatCumulativeBatch(pp$stats, 1, c(e, 0), sizes, iters, seed, "std")
cnt <- function(pp, e, iters) { b <- batch(pp, e, iters); c(b$leEs[1], b$geEs[1]) }
lo <- 0; hi <- 10000
while (hi - lo > 1) { mid <- (lo + hi) %/% 2; if (identical(cnt(pa, ea, mid), cnt(pb, eb, mid))) lo <- mid else hi <- mid }
ra <- batch(pa, ea, hi); ra0 <- batch(pa, ea, hi - 1); rb <- batch(pb, eb, hi); rb0 <- batch(pb, eb, hi - 1)
cat(sprintf("permutation %d; observed ES A %.17g B %.17g (diff %.3g)\n", hi, ea, eb, ea - eb))
cat(sprintf("A: random ES ~%.15g  le %d ge %d | B: random ES ~%.15g  le %d ge %d\n",
  (ra$geZeroSum - ra0$geZeroSum)[1], ra$leEs[1] - ra0$leEs[1], ra$geEs[1] - ra0$geEs[1],
  (rb$geZeroSum - rb0$geZeroSum)[1], rb$leEs[1] - rb0$leEs[1], rb$geEs[1] - rb0$geEs[1]))
