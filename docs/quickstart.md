# Quick start

## 1. Install and inspect the runtime

```bash
conda env create -f environment.yml
conda activate nf-rna
docker build -t rnaseq-control-plane:latest .
rnaseq doctor
```

`rnaseq new` is an interactive reviewed wizard. It limits species to Human or Mouse, offers FASTQ or raw counts, Salmon or HISAT2 + featureCounts when FASTQ is selected, managed/custom reference choices, and QC-only/L1/L2 scopes. It checks imported files before writing and never overwrites a populated target.

```bash
rnaseq new
```

It creates `project.yaml`, `metadata.csv`, `contrasts.csv`, `input/`, and `planning/`. Choose import mode to copy source data, or `--scaffold` for templates only; a scaffold deliberately remains invalid until real inputs and analytical metadata are supplied.

The exact same interface is scriptable without prompts:

```bash
rnaseq new --name demo --destination projects --species mouse \
  --input-type raw_counts --counts source/counts.csv \
  --metadata source/metadata.csv --contrasts source/contrasts.csv \
  --preset L2 --design-type two_group --condition-column condition --yes
```

For FASTQ import, pass `--fastq-samplesheet` with exact columns `sample,fastq_1,fastq_2,strandedness`. The wizard retains source lane filenames, infers layout, rejects mixed layouts or strandedness, and requires explicit strandedness for HISAT2 + featureCounts.

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

For a biological paired/repeated design, declare the block explicitly, for example `design: {type: paired, formula: "~ subject_id + condition", pairing_column: subject_id}`. This is unrelated to `input.layout: paired_end`.

For a final production-intended run, use a reviewed managed reference and immutable runtime selection:

```yaml
reference:
  source: local
  root: /absolute/reference-root
  manifest: reference_manifest.yaml
  acceptance: production
runtime:
  control_plane_image: rnaseq-control-plane:0.5.1
```

The managed manifest must use schema 1.1 with deliberate `purpose: production`. The initial human identity contract is Ensembl 116, GRCh38.p14; assets are not downloaded by this project.

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

This creates a new immutable case/run directory. QC-only FASTQ projects stop after the upstream backend and MultiQC delivery. L1 projects stop after L1 and the L1 report. L2 projects may run L2 and explicitly selected GSEA. Never use a small fixture or an unreplicated comparison to draw biological conclusions.

## Public smoke fixture

The bundled smoke fixture is intentionally small and public:

```bash
rnaseq validate examples/nfcore_smoke_test
rnaseq plan examples/nfcore_smoke_test
rnaseq run examples/nfcore_smoke_test --case-id SMOKE-001 --profile local --yes
```

It is an integration check, not a biological dataset. It may pull containers and consume local resources. Do not run it automatically in CI without making those costs explicit.
