# nf-rna

`nf-rna` is a production-oriented bulk RNA-seq workflow built around Nextflow, nf-core/rnaseq, Salmon, tximport, DESeq2, and clusterProfiler. It validates project contracts, creates deterministic plans, runs a pinned FASTQ workflow when explicitly authorized, and assembles curated technical delivery artifacts.

The stable command-line interface is `rnaseq`; the Python import namespace is `rnaseq`. The validated production runtime image intentionally remains `rnaseq-control-plane:latest` for compatibility with existing immutable run provenance.

This repository contains source code and deliberately small public test fixtures only. It does **not** publish client FASTQs, references, indices, immutable run directories, delivery packages, or Nextflow work state.

## What it supports

- **FASTQ route:** paired- or single-end FASTQs through pinned [nf-core/rnaseq 3.26.0](https://nf-co.re/rnaseq/3.26.0/) with Salmon pseudoalignment, then a frozen tximport handoff.
- **Raw-count route:** validated integer count matrices with `gene_id`, metadata, and explicit contrasts.
- **L1:** expression QC, filtering, normalization, VST, PCA, sample correlation, and a technical report.
- **L2:** L1 plus DESeq2 contrasts and, when configured, preranked GSEA using GO BP/MF/CC and KEGG.

L1 is an execution boundary: an L1 project never invokes the control plane's L2/DESeq2 or enrichment tasks merely because metadata or a contrast file is present. See the [architecture](docs/architecture.md).

## Scientific contract

The project validates only supported designs and explicit contrast directions. It does not infer biological groups, add covariates, alter metadata, silently change reference strategy, or reinterpret already-trimmed FASTQs. L2 uses DESeq2; GSEA ranks all finite tested-gene Wald statistics without a DEG or effect-size prefilter. Results are technical outputs, not medical or biological interpretation. Full details are in the [scientific contract](docs/scientific_contract.md).

## Installation

Python 3.11+ is required. The supplied Conda environment includes the R and Bioconductor runtime used by the first-party downstream image.

```bash
git clone <YOUR-REPOSITORY-URL>
cd RNA-seq
conda env create -f environment.yml
conda activate nf-rna
docker build -t rnaseq-control-plane:latest .
```

For FASTQ execution, also install Nextflow and start Docker Desktop (or another compatible Docker daemon). Run `rnaseq doctor` before a real run. See [runtime requirements](docs/runtime.md).

## Quick start

Create a project interactively, populate its input files, then validate before planning or executing:

```bash
rnaseq new
rnaseq validate path/to/project
rnaseq plan path/to/project
rnaseq doctor path/to/project
```

An immutable execution is explicitly authorized with a filesystem-safe client case ID:

```bash
rnaseq run path/to/project --case-id CASE-001 --profile local --yes
rnaseq status path/to/project
```

`rnaseq run` is not a dry run: it can invoke containers, process FASTQs, and consume substantial compute and storage. It creates a new immutable run; it does not overwrite a previous case/run. The [quickstart](docs/quickstart.md) shows minimal FASTQ and raw-count configurations.

## FASTQ configuration

FASTQ projects specify layout and preprocessing explicitly. The declaration is authoritative; preprocessing is never inferred from a filename or directory.

```yaml
input:
  type: fastq
  path: input/fastq
  layout: paired_end
  preprocessing: raw       # or: pretrimmed
```

- `raw` keeps nf-core trimming enabled.
- `pretrimmed` freezes a native boolean `skip_trimming: true` for nf-core, so reads are not trimmed again.

FASTQ execution uses an execution-ready reference. A managed local reference is declared with a portable placeholder rather than a machine-specific path:

```yaml
reference:
  source: local
  root: <REFERENCE_ROOT>
  manifest: reference_manifest.yaml
```

`<REFERENCE_ROOT>` is replaced by the operator's absolute local reference directory at deployment time; it is not a path to commit. The manifest binds asset checksums and the selected Salmon strategy. Reference preparation/adoption is explicit and is never performed implicitly by validation or planning.

## Raw-count configuration

Raw-count projects set `input.type: raw_counts`, provide a non-negative integral matrix whose first column is `gene_id`, and use strict contrast columns:

```yaml
input:
  type: raw_counts
  path: input/counts.csv
upstream:
  engine: external
  provider: external_provider
  quantification_method: unknown
```

`metadata.csv` is extensible but must contain `sample_id` and each variable in the configured formula. `contrasts.csv` is strict and has exactly `contrast_id,factor,numerator,denominator`.

## Small public smoke fixture

`examples/simple_fastq_two_group` contains tiny structural FASTQ fixtures for validation tests only. `examples/nfcore_smoke_test` contains a small, public nf-core test-dataset fixture suitable for an explicitly authorized local smoke run; it is not a biological study and must not be interpreted as one.

```bash
rnaseq validate examples/nfcore_smoke_test
rnaseq plan examples/nfcore_smoke_test
rnaseq run examples/nfcore_smoke_test --case-id SMOKE-001 --profile local --yes
```

The smoke route downloads or uses local container assets as needed and is intentionally not run by `pytest`.

## Validation status

At the public-packaging review, the non-expensive regression suite completed with **187 passed, 1 deselected** (the deselected test is an explicitly expensive L2 determinism test). A fresh L1 FASTQ smoke run completed upstream nf-core/Salmon, staged the immutable handoff, ran L1 and an L1-only report, assembled delivery, and did not invoke control-plane L2 or GSEA.

These checks demonstrate the repository state at review time; repeat them in your own runtime before relying on a release.

## Runtime requirements and limitations

- Local execution requires Docker, Nextflow, adequate local disk, and a compatible host/container architecture. `rnaseq doctor` reports prerequisites and capacity warnings.
- The default local resource ceiling is 6 CPUs, 12 GiB memory, and 12 hours; it is a scheduler bound and not a scientific parameter.
- L2 requires appropriate biological replication. The pipeline intentionally does not make an unreplicated 1-vs-1 design suitable for DESeq2 inference.
- KEGG GSEA uses the configured online KEGG route and can report network unavailability; no silent alternate provider is substituted.
- This software generates technical analysis artifacts only. It does not make clinical decisions or provide biological interpretation.

See [architecture](docs/architecture.md), [scientific contract](docs/scientific_contract.md), [runtime](docs/runtime.md), and [quickstart](docs/quickstart.md).

## License and citation

Released under the [MIT License](LICENSE). If you use this work, cite the specific release; a machine-readable record is available in [CITATION.cff](CITATION.cff).
