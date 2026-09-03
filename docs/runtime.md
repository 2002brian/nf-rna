# Runtime and reproducibility

## Required local components

- Python 3.11+ and the supplied Conda environment.
- Docker Desktop or a compatible running Docker daemon.
- Nextflow for FASTQ execution.
- The first-party image built locally as `rnaseq-control-plane:latest`. This is the validated production runtime-image name; it is intentionally distinct from the public `nf-rna` project name.
- Sufficient local disk for Nextflow cache/work and the selected reference.

```bash
conda env create -f environment.yml
conda activate nf-rna
docker build -t rnaseq-control-plane:latest .
rnaseq doctor
```

Run `rnaseq doctor <PROJECT>` before an authorized FASTQ run. It checks the runtime image and its required executables/packages, reports host/Docker facts, and evaluates project/reference readiness. It does not execute a workflow.

## Versioned execution

The FASTQ route is pinned to nf-core/rnaseq 3.26.0. The first-party downstream image includes Python, R, DESeq2, tximport, plotting, and enrichment packages. Container image availability, host architecture, and registry access remain operator/runtime responsibilities.

## Local resource policy

The default local contracts are intentionally conservative:

| class | CPUs | memory | time |
| --- | ---: | ---: | ---: |
| SMALL | 1 | 2 GiB | 2 h |
| MEDIUM | 4 | 8 GiB | 8 h |
| LARGE | 6 | 12 GiB | 12 h |

The global local ceiling is 6 CPUs, 12 GiB, and 12 hours. L1/report processes are SMALL; L2 and each GSEA task are MEDIUM; enrichment and Salmon quantification are constrained to avoid uncontrolled concurrent memory use. These are scheduling bounds only and do not alter quantification, filtering, DESeq2, thresholds, rankings, or GSEA calculations.

## Filesystem layout

Persistent immutable artifacts live beneath the project `runs/` directory. Operational Nextflow state lives in a local execution root outside that tree. On macOS the default is under the user's cache directory; on Linux it follows the standard cache location. Set `RNASEQ_EXECUTION_ROOT` to an absolute operator-owned scratch location if required.

Do not commit project run directories, workflow work directories, logs, delivery packages, references, indices, or production inputs. The repository's `.gitignore` and narrow `.dockerignore` enforce this default for new Git users.

## Determinism and external dependencies

Validation and planning are deterministic for unchanged inputs. Execution captures immutable contracts, checksums, versions, and selected resource profile. Timestamps, Docker scheduling, network conditions, and the external KEGG service are inherently runtime-dependent and are recorded rather than treated as deterministic results.
