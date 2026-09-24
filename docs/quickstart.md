# nf-rna Quick Start

This guide is for production users on Linux x86-64 or Windows WSL2. Install Python 3.11+, Git, Conda (for example Miniforge) with the `conda-forge` and `bioconda` channels configured in that order (`conda config --add channels bioconda && conda config --add channels conda-forge && conda config --set channel_priority strict`), Bash, Java 17+, and Nextflow before starting. Nextflow runs on the host and runs every analysis task in a Conda environment; Docker is not required. See the [official Nextflow installation guide](https://docs.seqera.io/nextflow/install) for Java and Nextflow.

## Install a released nf-rna version

nf-rna is installed from Git tags; it is not published to PyPI. Installing from
Git lets pip record the exact source commit, which every run records; nf-rna
refuses to run from an install without a commit. Replace `X.Y.Z` with the
released version you intend to use:

```bash
RELEASE_VERSION=X.Y.Z
conda create --yes --name nf-rna python=3.11
conda activate nf-rna
python -m pip install --upgrade pip
python -m pip install "git+https://github.com/2002brian/nf-rna.git@v${RELEASE_VERSION}"
rnaseq --version
rnaseq doctor
```

The installed package contains required workflow assets. `rnaseq doctor` is read-only: it checks Java, Nextflow, Conda and its channel order, the nf-core/rnaseq pin, the downstream Conda lock, and the nf-rna source revision without creating environments or downloading anything, and ends with `Overall: READY` or `Overall: NOT READY`. The first run creates the locked downstream Conda environment; Nextflow creates the nf-core/rnaseq and HISAT2/featureCounts process environments in a shared Conda cache.

Releases `v1.1.2` through `v1.2.1` used a Docker execution image published as
`ghcr.io/2002brian/nf-rna:X.Y.Z`; `v1.2.1` is the last Docker-based release.
`v1.0.0` predates GHCR publication, and `v1.1.0`/`v1.1.1` have no official image.

## Create and run a project

```bash
mkdir -p ~/projects/rnaseq-projects
cd ~/projects/rnaseq-projects
rnaseq new
cd <new-project>
# Add FASTQs to input/fastq/, or a count matrix to input/counts.csv.
# Complete metadata.csv and contrasts.csv.
rnaseq validate .
rnaseq plan .
rnaseq doctor .
rnaseq run . --case-id CASE-001 --profile local --yes
rnaseq status .
```

The wizard writes `project.yaml`, `metadata.csv`, `contrasts.csv`, `input/`, and `planning/`; it does not import data unless explicit import flags are supplied. FASTQ projects require an execution-ready reference.

## References and inputs

Register an existing checksum-bound managed reference once:

```bash
rnaseq reference register /absolute/reference-root
```

Registration validates and records a local pointer to the existing manifest; it does not download, copy, or build assets. If you need to create an index rather than adopt one, use the optional builder workflow in [Runtime](runtime.md#managed-reference-builder-host-native).

For `input/counts.csv`, the first column is `gene_id` and remaining columns are sample IDs. For FASTQ, accurately declare `input.layout`, `input.preprocessing`, `upstream.strandedness`, the selected quantification backend, and a matching reference. HISAT2 + featureCounts requires explicit `unstranded`, `forward`, or `reverse` strandedness.

Only `--profile local` is implemented. HPC, SLURM, and Apptainer are not currently supported.

## Reproducible runtime

The release contract is:

```text
CLI X.Y.Z  ↔  Git tag vX.Y.Z  ↔  one source commit, recorded in every run
```

Install an explicit release tag for formal analyses. `rnaseq doctor PROJECT`
reports the source revision, the downstream Conda lock, and any provisioned
runtime. Run provenance records the nf-rna version and source commit, the
downstream Conda lock filename and SHA-256, the installed wheel SHA-256, the
bundled R-script inventory, the nf-core/rnaseq version and revision, the
reference manifest and file checksums, and executed workflow hashes.

## Historical v1.0.0 compatibility

The immutable `v1.0.0` release predates both GHCR publication and the
version-matched project default. It has no official GHCR image. Its `rnaseq
new` wizard writes `execution_image: nf-rna:latest`; users working from that
historical source release should build that local development image before
planning or running:

```yaml
runtime:
  execution_image: nf-rna:latest
```

For example, update the generated `project.yaml` in an editor after creating the project. This is a v1.0.0 compatibility note only; it is not part of the normal Quick Start.
