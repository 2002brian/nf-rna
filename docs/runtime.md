# Runtime and reproducibility

## Required local components

- Python 3.11+ and the supplied Conda environment.
- Docker Desktop or a compatible running Docker daemon.
- Nextflow for FASTQ execution.
- A first-party control-plane image selected by `runtime.control_plane_image`.
- Sufficient local disk for Nextflow cache/work and the selected reference.

```bash
conda env create -f environment.yml
conda activate nf-rna
docker build -t rnaseq-control-plane:latest .
rnaseq doctor
```

## Managed-reference builder (host-native)

Both managed-reference commands run native executables on the host; they never
invoke Docker or fall back to a container:

```bash
mamba env create -f environment.reference-builder.yml
mamba activate nf-rna-reference-builder
rnaseq reference prepare /absolute/reference-root --threads 4
rnaseq reference prepare-hisat2 /absolute/reference-root --threads 4
```

The separately pinned builder environment contains Salmon 1.10.3, HISAT2
2.2.1 (including `hisat2_extract_splice_sites.py`), and RSEM 1.3.3. Each
command resolves absolute executable paths and validates versions before it
creates staging output. Salmon retains its decoy-aware gentrome strategy with
`k=31`; HISAT2 retains the annotation-aware genome index plus generated splice
sites. Docker remains required for FASTQ workflow execution and the downstream
control-plane contract, not for index construction.

Preparation records host resource facts, selected threads, source-relative
paths/checksums, tool version output, tokenized build arguments, OS and
architecture in the managed reference manifest. It builds in a sibling staging
directory, validates index artifacts, then atomically publishes the index and
manifest. A failed staging directory is deliberately retained for inspection;
the previous published index and manifest remain unchanged. Human
annotation-aware HISAT2 construction can require substantially more memory
than a 64 GiB machine. Swap availability is not evidence that this build is
ready.

### Human Ensembl 116 genome-only HISAT2 compatibility

The production Human Ensembl 116/GRCh38.p14 bundle may declare
`genome_only_runtime_splices`: a HISAT2 2.2.3 genome-only index and a
registered splice-site file derived from the same GTF. This strategy never
claims graph-embedded splice sites. The first-party 2.2.1 runtime receives
exactly one `--known-splicesite-infile` argument only for this strategy; the
manifest separately records `index_builder_version` and
`runtime_aligner_version`. Until its dedicated smoke acceptance is recorded,
the 2.2.3-builder/2.2.1-runtime pair is `requires_smoke_validation` and is
not execution-ready. The narrow builder environment is
`environment.reference-builder-hisat2-2.2.3.yml`; it does not silently change
the validated 2.2.1 graph-embedded builder/runtime contract.

## Delivery count semantics

New runs place source-labelled matrices in `delivery/counts/` with an
`artifact_manifest.json` containing the checksum, dimensions, sample order and
DESeq2 construction method. `raw_counts.csv` exists only for imported integer
counts and featureCounts integer gene counts. Salmon/tximport emits
`estimated_counts.csv`, which can be non-integer and retains tximport offset
semantics. `vst.csv` is normalized/transformed expression for visualization;
it is never used for DESeq2 fitting and must not replace the raw or estimated
matrix. QC-only runs may deliver the canonical upstream source matrix but do
not fabricate VST. Do not compare raw-count magnitudes across samples without
an appropriate normalization.

Run `rnaseq doctor <PROJECT>` before an authorized FASTQ run. It checks the runtime image and its required executables/packages, reports host and Docker CPU/RAM, the selected local ceiling, free space at the configured Nextflow work location, and project/reference readiness. It warns when Docker has fewer CPUs or less memory than the selected local ceiling. It does not execute a workflow.

Only the local execution profile is implemented. Workstation/HPC and SLURM execution remain deferred to the resource-profile milestone.

`rnaseq-control-plane:latest` is allowed for non-production development and is labelled accordingly. Production acceptance requires an image digest or versioned tag plus a Docker-observed image ID/digest. `rnaseq doctor <PROJECT>` reports requested and observed identities. Each run also freezes those identities, the Git commit resolved from the installed source checkout, and SHA-256 values for the executed first-party workflow files.

On Linux and WSL, downstream Nextflow task containers run with the invoking host user's numeric UID:GID. nf-rna freezes a per-run Nextflow override equivalent to Docker `--user $(id -u):$(id -g)`, so host-created task directories remain writable without `sudo`, `chmod 777`, or changing the fixed image user. macOS keeps its existing Docker Desktop behavior and does not receive this override. `rnaseq doctor` reports the selected policy and Linux/WSL mapping before execution.

## Versioned execution

The Salmon FASTQ route is pinned to nf-core/rnaseq 3.26.0. The HISAT2 route uses Biocontainers build tags for HISAT2 2.2.1, SAMtools 1.21, Subread/featureCounts 2.0.6, FastQC 0.12.1 and fastp 0.24.0, plus Seqera Wave MultiQC 1.33; each resolved Docker digest is recorded only when the runtime can inspect it. The first-party downstream image includes Python, R, DESeq2, tximport, plotting, and enrichment packages. Container image availability, host architecture, and registry access remain operator/runtime responsibilities.

## Milestone A validation status

On 2026-09-05, the pinned production semantics fixture passed under Docker Desktop 29.7.2 on an arm64 host (the pinned amd64 process images ran under Docker emulation): SAMtools 1.21 `sha256:783c6646029a306ec5e4162009dc1a20d8f6c528f7c380e5b4affbf12d9112e5` and Subread/featureCounts 2.0.6 `sha256:114390a783c77f7739d86e474bedfa5a4e65309a2f71d4db430803fb04601f5d`. It verified single-end counting, paired fragments counted once, forward/reverse strands, ambiguous-overlap exclusion, both-mates and chimeric-fragment policy, and technical-lane merging. A primary `NH:i:2` alignment remained excluded after the production `samtools view -bh -F 0x900` transformation removed its secondary record; primary-only filtering therefore did not make a multimapper appear unique.

The successful raw L2 smoke run was `MILESTONE-A-RAW-L2-FINAL/20260905-115310+0800`; its first-party workflow identity was `nf-rna/hisat2_featurecounts` with source SHA-256 `e5408f74c1afd4c41bcd51e7b06868d7b7df9da4b77aff2ddb5483492810d05c`. It completed raw fastp/FastQC, HISAT2, original and published count-only BAM lineage, featureCounts, MultiQC, canonical matrix/sample map, `featurecounts_raw_counts` → `DESeqDataSetFromMatrix`, L1, L2, report, and delivery. The independent pretrimmed L1 branch also completed as `MILESTONE-A-PRETRIMMED-L1/20260905-115200+0800`. These immutable validation records are intentionally outside the source distribution.

The smoke intentionally selected no enrichment backend, so the documented no-enrichment L2 route was exercised; GO/KEGG enrichment was not claimed for synthetic gene identifiers. This is a synthetic execution check, not biological evidence.

## Single-end HISAT2 QC acceptance

On 2026-09-05, the separate interactive-wizard single-end QC run `WIZARD-HISAT2-QC-R3/20260905-144237+0800` passed in a local validation workspace that is not distributed with the source. It exercised raw preprocessing, HISAT2, featureCounts, FastQC/MultiQC, and QC delivery for the single-end synthetic fixture only. It did not exercise paired fragment accounting or downstream L1/L2/enrichment; those claims belong respectively to the paired acceptance record below and the Milestone A smoke records above.

## Paired-end HISAT2 QC acceptance

On 2026-09-05, a public interactive-wizard paired-end QC run passed as `PAIRED-HISAT2-QC-UAT-R3/20260905-150643+0800`. Its immutable validation record remains outside the source distribution. Its durable synthetic fixture is `tests/fixtures/hisat2_paired_raw_v2`; regenerate a new empty fixture root with `python tests/fixtures/hisat2_paired_raw/generate_fixture.py --root PATH/TO/EMPTY/hisat2_paired_raw_v2`, then prepare its managed index with `rnaseq reference prepare-hisat2 PATH/TO/EMPTY/hisat2_paired_raw_v2/reference`. The fixture is forward stranded (HISAT2 `FR`; featureCounts `-s 1`), has PairAlpha lanes 001/002 and separate PairBeta, and defines counts before execution: PairAlpha `GeneA=3, GeneB=1`; PairBeta `GeneA=0, GeneB=2`.

The successful run retained synchronized R1/R2 pairs after raw fastp processing, trimmed the three short-insert adapter pairs (two reads and 66 bases per lane), mapped every lane at 100%, and exactly reproduced the fixture counts. featureCounts executed `-p --countReadPairs -B -C`; PairAlpha's eight alignment records represented four fragments and PairBeta's four records represented two fragments. Original and count-only BAMs passed `samtools quickcheck`; count-only BAMs had no `0x900` records and retained `NH:i:1` on every mapped record. Its MultiQC contains FastQC, fastp, HISAT2, and featureCounts, and the QC scope created no downstream L1/L2/enrichment artifacts. The run freezes workflow SHA-256 `ead01d1a6c84f20ad8d7f2c4943c97be20d748cc6adff0fe060520e852706d33` and the resolved container identities in provenance. This remains a software-contract fixture, not biological evidence.

## Local resource policy

`rnaseq new` records an optional aggregate local execution ceiling under
`execution.max_cpus` and `execution.max_memory_gb`. It is an execution
preference, not a biological or statistical setting. The interactive wizard
detects host CPU/memory on Linux/WSL and macOS, reserves approximately 20%, and
rounds memory conservatively; a 20-CPU / 64-GiB machine normally suggests 16
CPUs / 48 GiB. Existing projects without this section retain 8 CPUs / 12 GiB.
`rnaseq doctor PROJECT` compares the requested ceiling with host and Docker
capacity on the machine actually running the project. Each run freezes the
effective local Nextflow configuration for upstream and downstream workflows.

The ceiling is total executor capacity, not a per-task request. HISAT2,
Salmon, featureCounts, process directives and `maxForks` keep their own
declared limits. Reference preparation remains independent: use its explicit
`--threads` option.

The default local contracts are intentionally conservative:

| class | CPUs | memory | time |
| --- | ---: | ---: | ---: |
| SMALL | 1 | 2 GiB | 2 h |
| MEDIUM | 4 | 8 GiB | 8 h |
| LARGE | 8 | 12 GiB | 12 h |

The global local ceiling is 8 CPUs, 12 GiB, and 12 hours, with one project run at a time. Salmon and HISAT2 alignment each use `maxForks 1`; BAM sorting/processing and featureCounts use `maxForks 2`; cohort aggregation, MultiQC, L1/L2, enrichment, and reporting use `maxForks 1`. The first-party HISAT2 workflow binds fastp, HISAT2, supported SAMtools operations, and featureCounts threads to `task.cpus`. Each FASTQ run freezes one path-free local Nextflow configuration and its SHA-256, and supplies it to nf-core Salmon, HISAT2/featureCounts, and downstream. These are scheduling bounds only and do not alter quantification, filtering, DESeq2, thresholds, rankings, or GSEA calculations.

## Filesystem layout

Persistent immutable artifacts live beneath the project `runs/` directory. Operational Nextflow state lives in a local execution root outside that tree. On macOS the default is under the user's cache directory; on Linux it follows the standard cache location. Set `RNASEQ_EXECUTION_ROOT` to an absolute operator-owned scratch location if required.

Do not commit project run directories, workflow work directories, logs, delivery packages, references, indices, or production inputs. The repository's `.gitignore` and narrow `.dockerignore` enforce this default for new Git users.

## Determinism and external dependencies

Validation and planning are deterministic for unchanged inputs. Execution captures immutable contracts, checksums, versions, and selected resource profile. Timestamps, Docker scheduling, network conditions, and the external KEGG service are inherently runtime-dependent and are recorded rather than treated as deterministic results.
