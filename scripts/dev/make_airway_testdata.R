#!/usr/bin/env Rscript
# Development test data: Bioconductor `airway` (Himes et al. 2014) as an nf-rna
# raw-counts project.  Human airway smooth muscle, 4 cell lines x
# dexamethasone treated/untreated, Ensembl gene IDs, 8 samples.
#
# airway is NOT part of the nf-rna downstream runtime.  Run this script with an
# R that has airway installed separately, for example:
#
#   CONDARC=/home/brian/nf-rna-work/conda/condarc.yml conda create -y \
#     --prefix /home/brian/nf-rna-work/conda/airway-testdata-env \
#     bioconductor-airway bioconductor-summarizedexperiment r-jsonlite
#   /home/brian/nf-rna-work/conda/airway-testdata-env/bin/Rscript \
#     scripts/dev/make_airway_testdata.R [OUTPUT_DIR]
#
# OUTPUT_DIR defaults to /home/brian/nf-rna-work/testdata/airway and is kept
# outside the repository.  The counts are written unchanged (all genes, no
# filtering); only the metadata column names are mapped to the nf-rna schema:
# dex -> condition (trt/untrt), cell -> cell.

suppressPackageStartupMessages({
  library(airway)
  library(SummarizedExperiment)
  library(jsonlite)
})

args <- commandArgs(trailingOnly = TRUE)
out <- if (length(args)) args[[1]] else "/home/brian/nf-rna-work/testdata/airway"
dir.create(out, recursive = TRUE, showWarnings = FALSE)

data("airway", package = "airway", envir = environment())
counts <- assay(airway, "counts")
stopifnot(is.numeric(counts), all(counts == round(counts)), !anyDuplicated(rownames(counts)))
storage.mode(counts) <- "integer"

samples <- colnames(airway)
metadata <- data.frame(sample_id = samples,
                       cell = as.character(colData(airway)$cell),
                       condition = as.character(colData(airway)$dex),
                       stringsAsFactors = FALSE)
stopifnot(setequal(metadata$condition, c("trt", "untrt")), all(table(metadata$cell, metadata$condition) == 1))

count_table <- data.frame(gene_id = rownames(counts), counts, check.names = FALSE)
write.csv(count_table, file.path(out, "counts.csv"), row.names = FALSE, quote = FALSE)
write.csv(metadata, file.path(out, "metadata.csv"), row.names = FALSE, quote = FALSE)
write.csv(data.frame(contrast_id = "trt_vs_untrt", factor = "condition", numerator = "trt", denominator = "untrt"),
          file.path(out, "contrasts.csv"), row.names = FALSE, quote = FALSE)

writeLines(c(
  'schema_version: "1.2"',
  "project:",
  "  id: airway_dex",
  "  pipeline: bulk_rnaseq",
  "  preset: L2",
  "organism:",
  "  species: Homo sapiens",
  "input:",
  "  type: raw_counts",
  "  path: counts.csv",
  "design:",
  "  type: paired_two_group",
  '  formula: "~ cell + condition"',
  "  pair_id: cell",
  "metadata_file: metadata.csv",
  "contrasts_file: contrasts.csv",
  "upstream:",
  "  engine: external",
  "  provider: bioconductor_airway",
  "  quantification_method: integer_gene_counts",
  "reference:",
  "  source: igenomes",
  "  genome: null",
  "thresholds:",
  "  padj: 0.05",
  "  abs_log2fc: 1.0",
  "annotation:",
  "  organism: Homo sapiens",
  "  input_id_type: ENSEMBL",
  "  target_id_type: ENTREZID",
  "  gene_symbol_output: true",
  "analysis:",
  "  enrichment: []"
), file.path(out, "project.yaml"))

files <- c("counts.csv", "metadata.csv", "contrasts.csv", "project.yaml")
provenance <- list(
  dataset = "Bioconductor airway (Himes et al. 2014, GEO GSE52778)",
  airway_version = as.character(packageVersion("airway")),
  summarizedexperiment_version = as.character(packageVersion("SummarizedExperiment")),
  r_version = R.version.string,
  generated_at_utc = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
  genes = nrow(counts), samples = ncol(counts),
  gene_id_namespace = "Ensembl gene ID without version suffix (as shipped in airway)",
  metadata_mapping = list(condition = "colData(airway)$dex", cell = "colData(airway)$cell"),
  contrast = "trt vs untrt (condition), design ~ cell + condition",
  md5 = as.list(setNames(unname(tools::md5sum(file.path(out, files))), files))
)
write_json(provenance, file.path(out, "testdata_provenance.json"), auto_unbox = TRUE, pretty = TRUE)
cat(sprintf("airway %s: %d genes x %d samples -> %s\n", provenance$airway_version, nrow(counts), ncol(counts), out))
