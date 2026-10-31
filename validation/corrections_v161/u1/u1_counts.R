# U-1: per-set permutation counts of the fgsea simple stage for the two inputs (same seed,
# same pathways, same order).  Usage: Rscript u1_counts.R inputs_<collection>.rds
suppressPackageStartupMessages({library(fgsea); library(BiocParallel); library(data.table)})
inp <- readRDS(commandArgs(TRUE)[1])
simple <- function(stats) {
  pp <- fgsea:::preparePathwaysAndStats(inp$pw, stats, 10, 500, 1, "std")
  g <- do.call(rbind, lapply(pp$filtered, calcGseaStat, stats=pp$stats, returnLeadingEdge=TRUE, scoreType="std"))
  set.seed(1); seeds <- sample.int(10^9, 1)
  r <- fgsea:::fgseaSimpleImpl(pathwayScores=unlist(g[, "res"]), pathwaysSizes=pp$sizes, pathwaysFiltered=pp$filtered,
    leadingEdges=rep(list(NULL), length(pp$filtered)), permPerProc=10000, seeds=seeds, toKeepLength=length(pp$filtered),
    stats=pp$stats, BPPARAM=SerialParam(), scoreType="std")
  list(r=r, order=names(pp$stats), stats=pp$stats)
}
a <- simple(inp$statsA); b <- simple(inp$statsB)
cat("sorted gene order identical:", identical(a$order, b$order), "\n")
d <- which(a$order != b$order); cat("positions whose gene differs:", length(d), "\n")
if (length(d)) print(data.frame(pos=d, geneA=a$order[d], statA=sprintf("%.17g", a$stats[d]), geneB=b$order[d], statB=sprintf("%.17g", b$stats[d]))[1:min(10, length(d)), ])
cols <- c("nLeZero","nGeZero","nLeEs","nGeEs")
diff <- a$r[, ..cols] != b$r[, ..cols]
cat("pathways with any count change:", sum(rowSums(diff) > 0), "\n")
print(colSums(diff))
ch <- which(rowSums(diff) > 0)
x <- cbind(a$r[ch, .(pathway, size, ES)], a$r[ch, ..cols] - b$r[ch, ..cols])
print(x[, .N, by=.(size, dLeZero=nLeZero, dGeZero=nGeZero, dLeEs=nLeEs, dGeEs=nGeEs)][order(size)])
cat("geZeroMean/leZeroMean changed w/o count change:", sum(abs(a$r$geZeroMean - b$r$geZeroMean) > 1e-12 * abs(a$r$geZeroMean) | abs(a$r$leZeroMean - b$r$leZeroMean) > 1e-12 * abs(a$r$leZeroMean)), "\n")
k <- c(ch[a$r$size[ch] == 10][1], ch[a$r$size[ch] == 11][1])
for (i in k) cat(sprintf("size %d: A nGeZero %d nLeZero %d geZeroSum %.17g leZeroSum %.17g | B nGeZero %d nLeZero %d geZeroSum %.17g leZeroSum %.17g\n",
  a$r$size[i], a$r$nGeZero[i], a$r$nLeZero[i], a$r$geZeroMean[i]*a$r$nGeZero[i], a$r$leZeroMean[i]*a$r$nLeZero[i],
  b$r$nGeZero[i], b$r$nLeZero[i], b$r$geZeroMean[i]*b$r$nGeZero[i], b$r$leZeroMean[i]*b$r$nLeZero[i]))
cat("nGeZero+nLeZero A (size10):", a$r$nGeZero[k[1]] + a$r$nLeZero[k[1]], " B:", b$r$nGeZero[k[1]] + b$r$nLeZero[k[1]], "\n")
neg <- ch[a$r$ES[ch] < 0]; cat("changed sets with ES<0:", length(neg), " NES-relevant leZeroMean rel change max:", max(abs(a$r$leZeroMean[neg]/b$r$leZeroMean[neg] - 1)), "\n")
ra <- a$r$ES/ifelse(a$r$ES > 0, a$r$geZeroMean, abs(a$r$leZeroMean)); rb <- b$r$ES/ifelse(b$r$ES > 0, b$r$geZeroMean, abs(b$r$leZeroMean))
cat("NES changed (>1e-12 rel): pos", sum(abs(ra/rb-1)[a$r$ES>0] > 1e-12), " neg", sum(abs(ra/rb-1)[a$r$ES<=0] > 1e-12), "\n")
cat("all sets: sizes 10/11 count", sum(a$r$size %in% 10:11), "\n")
