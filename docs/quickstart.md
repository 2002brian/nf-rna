# nf-rna Quick Start

This document distinguishes the intended next-release installation UX from the immutable historical `v1.0.0` release. It does not create or announce a new release.

## Canonical path for the next published patch release

On Linux or WSL, install Python 3.11+, Git, Docker with a running daemon, Bash, Java 17+, and Nextflow before starting. Nextflow executes on the host; Docker executes analysis tasks. Follow the [official Nextflow installation guide](https://docs.seqera.io/nextflow/install) for Java and Nextflow.

When a patch release containing this change is published, substitute its actual version for `RELEASE_VERSION`:

```bash
RELEASE_VERSION=<published-version>
python3.11 -m venv ~/.venvs/nf-rna-${RELEASE_VERSION}
source ~/.venvs/nf-rna-${RELEASE_VERSION}/bin/activate
python -m pip install --upgrade pip
python -m pip install "git+https://github.com/2002brian/nf-rna.git@v${RELEASE_VERSION}"
docker pull ghcr.io/2002brian/nf-rna:${RELEASE_VERSION}
rnaseq --version
rnaseq doctor
```

No PyPI distribution exists. This is a non-editable Git installation and the installed package contains the required workflow assets. `doctor` is read-only: it does not pull or build an image.

Create, configure, inspect, and run a project:

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

The CLI maps its package version directly to the initial project image: CLI `X.Y.Z` writes `ghcr.io/2002brian/nf-rna:X.Y.Z`. A prerelease CLI such as `1.0.1rc1` writes its equally explicit prerelease image tag; publish and qualify that image before offering the prerelease to users. No manual YAML replacement is part of this path.

## References and inputs

The wizard creates `project.yaml`, `metadata.csv`, `contrasts.csv`, `input/`, and `planning/`; it does not import data unless explicit import flags are used. FASTQ projects require an execution-ready reference. Register an existing checksum-bound managed reference once:

```bash
rnaseq reference register /absolute/reference-root
```

Registration validates and records a local pointer to the existing manifest; it does not download, copy, or build assets. If you need to create an index rather than adopt one, use the optional builder workflow in [Runtime](runtime.md#managed-reference-builder-host-native).

For `input/counts.csv`, the first column is `gene_id` and remaining columns are sample IDs. For FASTQ, accurately declare `input.layout`, `input.preprocessing`, `upstream.strandedness`, the selected quantification backend, and a matching reference. HISAT2 + featureCounts requires explicit `unstranded`, `forward`, or `reverse` strandedness.

Only `--profile local` is implemented. HPC, SLURM, and Apptainer are not currently supported.

## Historical v1.0.0 compatibility

The currently published `v1.0.0` release is installable and qualified with `ghcr.io/2002brian/nf-rna:1.0.0`, but its project wizard predates the version-matched default. Use the following only for that released version:

```bash
python3.11 -m venv ~/.venvs/nf-rna-1.0.0
source ~/.venvs/nf-rna-1.0.0/bin/activate
python -m pip install "git+https://github.com/2002brian/nf-rna.git@v1.0.0"
docker pull ghcr.io/2002brian/nf-rna:1.0.0
rnaseq new
cd <new-project>
sed -i 's|execution_image: nf-rna:latest|execution_image: ghcr.io/2002brian/nf-rna:1.0.0|' project.yaml
```

The `sed` command is a v1.0.0-only compatibility step. It is not part of the canonical next-patch Quick Start and does not modify the v1.0.0 tag.

## Runtime identity and reproducibility

`rnaseq doctor PROJECT` reports requested and Docker-observed image identity. Use a versioned tag for normal reproducible work, or the qualified v1.0.0 digest for an immutable identity:

```text
ghcr.io/2002brian/nf-rna@sha256:ee60405181783ff075a1f4a9f452990c651e91152f5a9837c2a5df44a838decb
```

`ghcr.io/2002brian/nf-rna:latest` is only a convenience tag. The run provenance records the requested reference, Docker-observed image ID/repository digest, OCI revision label when available, CLI version, and executed workflow hashes.
