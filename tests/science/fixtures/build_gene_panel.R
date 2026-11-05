# Regenerates tests/science/fixtures/gene_panel.tsv from org.Hs.eg.db / GO.db.
#
# The panel is committed; this script documents how it was drawn and is run only
# when the fixture must change (Rscript build_gene_panel.R OUTPUT.tsv).  It uses
# the locked downstream runtime (org.Hs.eg.db 3.22.0, GO.db 3.22.0).  The
# regression tests never call it: they compute every expectation from the
# committed panel with code independent of nf-rna's R scripts.
suppressPackageStartupMessages({ library(AnnotationDbi); library(org.Hs.eg.db) })
args <- commandArgs(trailingOnly = TRUE)
out <- if (length(args)) args[[1]] else "gene_panel.tsv"
set.seed(20261105)
UP_TERM <- "GO:0006260"    # DNA replication (BP)
DOWN_TERM <- "GO:0006954"  # inflammatory response (BP)
ens <- suppressMessages(AnnotationDbi::select(org.Hs.eg.db, keys = keys(org.Hs.eg.db, "ENSEMBL"), keytype = "ENSEMBL", columns = "ENTREZID"))
ens <- ens[!is.na(ens$ENTREZID), ]
per_ens <- tapply(ens$ENTREZID, ens$ENSEMBL, function(x) length(unique(x)))
per_ent <- tapply(ens$ENSEMBL, ens$ENTREZID, function(x) length(unique(x)))
bp <- suppressMessages(AnnotationDbi::select(org.Hs.eg.db, keys = "BP", keytype = "ONTOLOGYALL", columns = "ENTREZID"))
bp_genes <- unique(bp$ENTREZID)
term_genes <- function(term) unique(suppressMessages(AnnotationDbi::select(org.Hs.eg.db, keys = term, keytype = "GOALL", columns = "ENTREZID"))$ENTREZID)
unique_ens <- ens[ens$ENSEMBL %in% names(per_ens)[per_ens == 1] & ens$ENTREZID %in% names(per_ent)[per_ent == 1], ]
pick_term <- function(term, n, exclude) {
  pool <- unique_ens[unique_ens$ENTREZID %in% setdiff(term_genes(term), exclude), ]
  pool[sort(sample(nrow(pool), n)), ]
}
up <- pick_term(UP_TERM, 40, term_genes(DOWN_TERM))
down <- pick_term(DOWN_TERM, 40, c(term_genes(UP_TERM), up$ENTREZID))
used <- c(up$ENTREZID, down$ENTREZID)
null_pool <- unique_ens[unique_ens$ENTREZID %in% setdiff(bp_genes, c(term_genes(UP_TERM), term_genes(DOWN_TERM))), ]
null <- null_pool[sort(sample(nrow(null_pool), 700)), ]
# Ensembl IDs with several Entrez targets, and Entrez IDs with several Ensembl IDs.
o2m_ids <- head(sort(names(per_ens)[per_ens >= 2]), 4)
m2o_entrez <- head(sort(names(per_ent)[per_ent == 2]), 3)
m2o_ids <- sort(unique(ens$ENSEMBL[ens$ENTREZID %in% m2o_entrez]))
panel <- rbind(
  data.frame(ensembl_id = up$ENSEMBL, role = "up_term"),
  data.frame(ensembl_id = down$ENSEMBL, role = "down_term"),
  data.frame(ensembl_id = null$ENSEMBL, role = "null"),
  data.frame(ensembl_id = o2m_ids, role = "one_to_many"),
  data.frame(ensembl_id = m2o_ids, role = "many_to_one"),
  data.frame(ensembl_id = paste0("ENSG9999999999", 1:5), role = "unmapped")
)
panel <- panel[!duplicated(panel$ensembl_id), ]
# Every fifth ID carries a version suffix: nf-rna must strip it for the OrgDb query only.
versioned <- seq_len(nrow(panel)) %% 5 == 0
panel$ensembl_id[versioned] <- paste0(panel$ensembl_id[versioned], ".", 1 + seq_len(nrow(panel))[versioned] %% 9)
write.table(panel, out, sep = "\t", quote = FALSE, row.names = FALSE)
cat("terms:", UP_TERM, DOWN_TERM, "rows:", nrow(panel), "\n")
