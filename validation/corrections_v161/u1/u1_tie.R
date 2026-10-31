# U-1, GO_BP: locate the one permutation whose zero-side counts differ between the two
# inputs, and read its random ES under scoreType "std", "pos" and "neg".
# fgsea 1.36.2 src/fastGSEA.cpp, calcGseaStatCumulative(): for "std", if the maximum
# positive deviation equals the maximum negative deviation exactly (res == resDown),
# the ES is set to 0, which then counts as both >= 0 and <= 0.
# Values are read as differences of cumulative sums (resolution about 1e-13), so
# "equal" below means equal at that resolution.
# Usage: Rscript u1_tie.R inputs_GO_BP.rds
suppressPackageStartupMessages(library(fgsea))
inp <- readRDS(commandArgs(TRUE)[1])
prep <- function(stats) fgsea:::preparePathwaysAndStats(inp$pw, stats, 10, 500, 1, "std")
pa <- prep(inp$statsA); pb <- prep(inp$statsB)
sizes <- c(10L, 11L, max(pa$sizes))  # K = max size keeps the random stream of the real run
set.seed(1); seed <- sample.int(10^9, 1)
batch <- function(pp, iters, type) fgsea:::calcGseaStatCumulativeBatch(pp$stats, 1, rep(0, 3), sizes, iters, seed, type)
lo <- 0; hi <- 10000
while (hi - lo > 1) {
  mid <- (lo + hi) %/% 2
  if (identical(batch(pa, mid, "std")$geZero[1:2], batch(pb, mid, "std")$geZero[1:2])) lo <- mid else hi <- mid
}
cat("first differing permutation:", hi, "of 10000 (C++ seed", seed, ")\n")
one <- function(pp, type, field) (batch(pp, hi, type)[[field]] - batch(pp, hi - 1, type)[[field]])[1:2]
for (lab in c("A reference stats (full precision)", "B nf-rna stats (15 significant digits)")) {
  pp <- if (startsWith(lab, "A")) pa else pb
  pos <- one(pp, "pos", "geZeroSum"); neg <- one(pp, "neg", "leZeroSum")
  std <- one(pp, "std", "geZeroSum") + one(pp, "std", "leZeroSum")
  ge <- one(pp, "std", "geZero"); le <- one(pp, "std", "leZero")
  for (j in 1:2) cat(sprintf("%s, prefix size %d: max positive deviation %.15g, max negative deviation %.15g, std ES %.15g, counted >=0: %d, <=0: %d\n",
    lab, sizes[j], pos[j], neg[j], std[j], ge[j], le[j]))
}
