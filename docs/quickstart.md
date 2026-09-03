# Quick start

## 1. Install and inspect the runtime

```bash
conda env create -f environment.yml
conda activate nf-rna
docker build -t rnaseq-control-plane:latest .
rnaseq doctor
```

`rnaseq new` is interactive. Choose a project name, destination, species, input type, layout (for FASTQ), L1/L2 preset, and supported design.

```bash
rnaseq new
```

It creates `project.yaml`, `metadata.csv`, `contrasts.csv`, `input/`, and `planning/`. It never creates a placeholder biological matrix or FASTQ.

## 2. Configure a raw-count project

Place a CSV count matrix at `input/counts.csv`, with `gene_id` as its first column. Populate metadata and strict contrasts, then use a project definition such as:

```yaml
schema_version: "1.1"
project: {id: demo_counts, pipeline: bulk_rnaseq, preset: L2}
organism: {species: Mus musculus}
input: {type: raw_counts, path: input/counts.csv}
design: {type: two_group, formula: "~ condition"}
metadata_file: metadata.csv
contrasts_file: contrasts.csv
upstream: {engine: external, provider: external_provider, quantification_method: unknown}
reference: {source: igenomes, genome: null}
thresholds: {padj: 0.05, abs_log2fc: 1.0}
analysis: {enrichment: []}
```

## 3. Configure a FASTQ project

Place matching reads under `input/fastq/` using supported `*_R1` / `*_R2` names. Set `preprocessing` accurately:

```yaml
input:
  type: fastq
  path: input/fastq
  layout: paired_end
  preprocessing: pretrimmed
```

`pretrimmed` is appropriate only when the reads were already adapter/quality trimmed; it causes nf-core trimming to be skipped. Use `raw` otherwise. Configure an execution-ready iGenomes or managed local reference; do not commit a real reference root, index, or FASTQs.

## 4. Validate and plan

```bash
rnaseq validate path/to/project
rnaseq plan path/to/project
```

Validation is read-only. Planning emits deterministic artifacts only after a valid project. Review warnings, contrast direction, readiness, and planning artifacts before executing.

## 5. Run deliberately

```bash
rnaseq doctor path/to/project
rnaseq run path/to/project --case-id CASE-001 --profile local --yes
rnaseq status path/to/project
```

This creates a new immutable case/run directory. L1 projects stop after L1 and the L1 report. L2 projects may run L2 and explicitly selected GSEA. Never use a small fixture or an unreplicated comparison to draw biological conclusions.

## Public smoke fixture

The bundled smoke fixture is intentionally small and public:

```bash
rnaseq validate examples/nfcore_smoke_test
rnaseq plan examples/nfcore_smoke_test
rnaseq run examples/nfcore_smoke_test --case-id SMOKE-001 --profile local --yes
```

It is an integration check, not a biological dataset. It may pull containers and consume local resources. Do not run it automatically in CI without making those costs explicit.
