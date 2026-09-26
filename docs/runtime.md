# Runtime and reproducibility

## Supported runtime (v1.3.0)

nf-rna v1.3.0 supports one runtime: **Linux x86-64, including Windows WSL2,
with Nextflow and Conda**. Docker is not required and is never invoked. Native
Windows and macOS are not supported by v1.3.0.

Required local components:

- Python 3.11+ for the `rnaseq` control plane, installed from a Git release tag
  so that pip records the source commit (see the [README](../README.md)).
- Conda on `PATH` (for example Miniforge), with the `conda-forge` and `bioconda`
  channels in that order in its effective configuration, as nf-core/rnaseq
  `-profile conda` requires. Other channels may also be present.
- Bash, Java 17+, and Nextflow for host-side workflow execution.
- Sufficient local disk for Conda environments, Nextflow cache/work, and the
  selected reference.

Execution path for a FASTQ project:

```text
rnaseq → Nextflow → nf-core/rnaseq 3.26.0 (-profile conda) → Salmon → quant.sf
       → nf-rna staging → locked downstream Conda env → tximport → DESeq2
       → enrichment → report → delivery
```

The HISAT2 + featureCounts route replaces the nf-core step with nf-rna's
first-party Nextflow graph, run in the exactly pinned
`workflow/envs/hisat2-featurecounts-linux-64.yml` Conda environment.

Runtime identity:

- **Source revision.** Every run records the git commit of the nf-rna source it
  executed. A source checkout reports its `HEAD`; an installed distribution
  reports the commit pip recorded for a `git+…@<ref>` install. nf-rna refuses to
  start a run when no commit is available (a plain directory, sdist, or wheel
  install) or when the source checkout has uncommitted changes.
- **Downstream environment.** R, DESeq2, tximport, plotting, and enrichment
  packages come from the reviewed lock
  `workflow/envs/locks/nf-rna-downstream-linux-64.lock.yml`, verified against
  `SHA256SUMS`. nf-rna creates one prefix per (lock, nf-rna wheel) pair under
  its runtime cache, installs the non-editable nf-rna wheel with
  `pip --no-deps`, and never updates a prefix in place, so an upgrade gets a new
  prefix while earlier runs' prefixes stay intact.
- **Upstream environments.** nf-core/rnaseq is pinned to 3.26.0 (revision
  `e7ca46272c8f9d5ceee3f71759f4ba551d3217a4`) and runs with `-profile conda`;
  Nextflow creates its process environments in a shared upstream Conda cache.
- **Provenance.** The downstream contract, `execution_manifest.yaml`, and
  `run_provenance.yaml` record the nf-rna version and source revision, the lock
  filename and SHA-256, the wheel filename, SHA-256 and origin, the bundled
  R-script inventory and digest, the upstream runtime, and the reference
  manifest and file checksums. Each R module's `scientific_provenance.json`
  records the same runtime identity and its package versions.

## Historical Docker releases (v1.1.2–v1.2.1)

The sections below that mention execution images, GHCR, or Docker describe
releases up to v1.2.1, which ran processes in Docker. v1.2.1 is the last
Docker-based release.

### Official container images (historical)

The historical `v1.0.0` release predates GHCR and has no official image.
`v1.1.0` and `v1.1.1` remain valid source releases, but neither has an official
GHCR image because their publication workflows stopped before image build or
push. Beginning with the `v1.1.2` release, official images are published
through GitHub Container Registry as:

```bash
docker pull ghcr.io/2002brian/nf-rna:X.Y.Z
```

Use an explicit version tag (for example, `1.1.2`) or a digest for a
reproducible analysis. `ghcr.io/2002brian/nf-rna:latest` is a convenience tag
only.
Official release images are qualified for `linux/amd64`. `linux/arm64` is not
yet independently qualified.

Automatic GHCR publication ended with v1.2.1: the publication workflow was
removed for v1.3.0, so v1.3.0 and later releases publish no Docker image and
leave the existing tags, including `latest`, unchanged. Through v1.2.1, normal
GHCR publication started automatically when a GitHub Release was
published. Official publishing accepted stable `vX.Y.Z` release tags only,
checked that the checked-out release commit and package version matched that tag,
then built and pushed both the explicit version tag and `latest`. It ran
`rnaseq --version` against the published version-tagged image afterward.
Per-release-tag workflow concurrency prevented ordinary duplicate runs.

The workflow did not create, edit, move, or republish a Git tag or GitHub
Release. The immutable historical `v1.0.0` GitHub Release predates this
workflow and has no official GHCR image. Official GHCR distribution begins with
the `v1.1.2` release. The `v1.1.2` hotfix changes no scientific analysis, R,
or Nextflow behavior.

Building an image from source is a developer/offline workflow; see
[Development installation](development.md). It is not part of normal user
installation.

The image runtime is intentionally verified with a non-login shell, matching
the environment Nextflow task containers inherit:

```bash
docker run --rm ghcr.io/2002brian/nf-rna:X.Y.Z \
  sh -c 'command -v ps && ps --version'
```

Do not use `sh -lc` for this check: login-shell startup files may replace the
image `PATH` and do not represent Nextflow task execution.

## Managed-reference builder (host-native)

Prebuilt, checksum-bound indexes are first-class managed-reference inputs.
Register them in `reference_manifest.yaml`; validation and planning do not
require nf-rna to have built them. The two builder commands are optional
host-native conveniences and never invoke Docker or fall back to a container:

```bash
mamba env create -f environment.reference-builder.yml
mamba activate nf-rna-reference-builder
rnaseq reference prepare /absolute/reference-root --threads 4
# replace an existing Salmon declaration (for example a rejected cDNA index):
rnaseq reference prepare /absolute/reference-root --threads 4 --rebuild-salmon
rnaseq reference prepare-hisat2 /absolute/reference-root --threads 4
rnaseq reference register /absolute/reference-root
```

`register` performs manifest validation only; it does not prepare or rebuild
assets. It writes a machine-local reference pointer registry used by the
interactive `rnaseq new` wizard after species and backend selection. A
candidate must be a production manifest and have the selected backend's fully
validated asset; stale registrations are ignored so manual and custom routes
remain usable.

The separately pinned builder environment contains Salmon 1.10.3, bioconda
RSEM 1.3.3 (the version nf-core/rnaseq 3.26.0 uses; its binaries report
`v1.3.1`, so the Conda package record is the checked identity), and exactly
HISAT2 2.2.3 (including `hisat2_extract_splice_sites.py`).
Other 2.2.x releases are rejected so that reference construction remains
reproducible.
`hisat2-build --version` is the authoritative version check. The splice-site
helper has no semantic-version interface (`-v` is verbose mode), so preflight
instead requires the helper beside `hisat2-build` and verifies its supported
`-h` help contract. Provenance records the helper path and its association with
the checked HISAT2 binary rather than inventing a helper version.
Each command resolves absolute executable paths and validates versions before it
creates staging output. Salmon preparation derives its transcriptome exactly as
nf-core/rnaseq 3.26.0 does (`MAKE_TRANSCRIPTS_FASTA`): the GTF is filtered like
`CUSTOM_GTFFILTER` (records on genome FASTA sequences that carry a
`transcript_id`) and `rsem-prepare-reference --gtf` extracts transcripts from
the genome FASTA, so transcript names equal GTF `transcript_id`. It then
requires a PASS transcript-ID contract and builds the decoy-aware gentrome
index with `k=31` under `salmon/gtf_derived/`. The registered
`files.transcript_fasta` (for example Ensembl cDNA) is not used: its versioned
IDs (`ENSMUST00000200568.2`) never equal an Ensembl GTF `transcript_id`
(`ENSMUST00000200568` with a separate `transcript_version`), and `cdna.all`
omits non-coding transcripts such as lncRNA.

Every built Salmon declaration must reference a checksum-bound transcript-ID
contract artifact (`salmon.validation.artifact`) with status PASS: every
transcriptome ID equals exactly one GTF `transcript_id` with one `gene_id`, no
ID is duplicated, and every GTF transcript on a genome FASTA sequence is in the
transcriptome. No identifier is normalized. nf-core/rnaseq builds tx2gene by
exact equality between Salmon transcript names and GTF attributes, so
validation, planning, `rnaseq doctor`, and run preflight reject a Salmon
reference without such a contract. `--rebuild-salmon` replaces such a
declaration: the previous manifest is archived as
`reference_manifest.<sha256-prefix>.yaml`, the previous index directory is left
untouched, and the new provenance records what it replaced and why. HISAT2 builds a
genome-only index and registers GTF-derived splice sites for the runtime
`--known-splicesite-infile` argument.

For a prebuilt Salmon index, declare `index`, `version`, `strategy`,
`source_transcriptome_sha256`, and `validation.artifact` (a PASS transcript-ID
contract; `rnaseq reference adopt-salmon-index` computes and writes it);
`decoy_aware` additionally requires the matching `source_genome_sha256`. For a prebuilt HISAT2 index, declare `index_prefix`,
`version`, `strategy: genome_only_runtime_splicesites`, `genome_fasta_sha256`,
`source_gtf_sha256`, and a checksum-bound `splice_sites` asset. All paths are
reference-root-relative. The manifest validates asset/index consistency; native
versus external builder mode is provenance, but the declared HISAT2 version
must be exactly 2.2.3.
An HISAT2-only reference need not include `files.transcript_fasta`; that asset
is required only when a Salmon index is declared or Salmon preparation is used.

For example, a manifest can register external assets without copying or
rebuilding them:

```yaml
salmon:
  status: built
  index: vendor/salmon/index
  version: "1.10.3"
  strategy: decoy_aware
  source_transcriptome_sha256: "<files.transcript_fasta SHA-256>"
  source_genome_sha256: "<files.genome_fasta SHA-256>"
  validation: {artifact: vendor/salmon/index.transcript_id_contract.json, transcript_id_contract: exact}
hisat2:
  status: built
  index_prefix: vendor/hisat2/genome
  version: "2.2.3"
  strategy: genome_only_runtime_splicesites
  genome_fasta_sha256: "<files.genome_fasta SHA-256>"
  source_gtf_sha256: "<files.annotation_gtf SHA-256>"
  splice_sites_gtf_sha256: "<files.annotation_gtf SHA-256>"
  splice_sites: {path: vendor/hisat2/splice_sites.txt, sha256: "<SHA-256>"}
```

Preparation records host resource facts, selected threads, source-relative
paths/checksums, tool version output, tokenized build arguments, OS and
architecture in the managed reference manifest. It builds in a sibling staging
directory, validates index artifacts, then atomically publishes the index and
manifest. A failed staging directory is deliberately retained for inspection;
the previous published index and manifest remain unchanged. Human HISAT2
construction can require substantially more memory than a 64 GiB machine.
Swap availability is not evidence that this build is ready.

### Human Ensembl 116 genome-only HISAT2 compatibility

The production Human Ensembl 116/GRCh38.p14 bundle may declare
`genome_only_runtime_splicesites`: a HISAT2 genome-only index and a
registered splice-site file derived from the same GTF. This strategy never
claims graph-embedded splice sites. The first-party 2.2.3 runtime receives
exactly one `--known-splicesite-infile` argument only for this strategy; the
manifest separately records `index_builder_version` and
`runtime_aligner_version`. A prebuilt index without an explicit
`runtime_compatibility: validated` remains `requires_smoke_validation`, even
when its builder version is 2.2.3. This prevents a version-pin update from
silently substituting for real FASTQ acceptance evidence.

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

Run `rnaseq doctor <PROJECT>` before an authorized FASTQ run. On Linux/WSL2 it checks Java, Nextflow, Conda, the effective Conda channel order (`conda-forge` before `bioconda`, read from `conda config --show channels --json`), the nf-core/rnaseq pin and cache, the pinned HISAT2/featureCounts Conda environment, the downstream Conda lock and any provisioned prefix, the nf-rna source revision, host CPU/RAM, the selected local ceiling, free space at the configured Nextflow work location, and project/reference readiness. It ends with `Overall: READY`, or `Overall: NOT READY` and exit code 1 when any check fails. It does not create environments or execute a workflow. `rnaseq run` and `rnaseq retry` repeat the channel check for the nf-core Salmon route before creating a run.

Only the local execution profile is implemented. Workstation/HPC and SLURM execution remain deferred to the resource-profile milestone.

Historical Docker runtime (v1.2.1 and earlier): new projects used a version-matched `runtime.execution_image`; a local `nf-rna:latest` image remains allowed only for non-production development. Production acceptance requires a versioned tag or digest plus a Docker-observed image ID/digest. `rnaseq doctor <PROJECT>` reports requested and observed identities. Each run freezes those identities and the OCI build revision label as the canonical execution `source_revision`; it is propagated unchanged into every downstream module configuration, scientific provenance document, technical report, and delivery copy. Official release builds supply the exact release commit SHA, while dirty development builds supply an explicit dirty identity. An unlabeled container is explicitly recorded as `unlabeled-container-image`; nf-rna never substitutes the control-plane checkout, a tag, or `unknown`. `runtime.control_plane_image` remains a read-only compatibility alias for existing projects; newly created and serialized configuration uses `runtime.execution_image`. See [v1.0.0 compatibility](quickstart.md#historical-v100-compatibility) for the historical project-default limitation.

In that historical Docker runtime, `rnaseq` validated, froze, and launched Nextflow; Nextflow owned Docker container launch through explicit process-level `container params.first_party_image` directives. The one first-party execution image contains the installed `rnaseq.workflow_support` package and `/opt/nf-rna/r` scripts used by L1, L2, GSEA, and technical-report tasks. Docker is the current local process runtime; an Apptainer/Singularity profile can be added later without moving scientific execution into Python.

In that historical Docker runtime on Linux and WSL, downstream Nextflow task containers ran with the invoking host user's numeric UID:GID. nf-rna freezes a per-run Nextflow override equivalent to Docker `--user $(id -u):$(id -g)`, so host-created task directories remain writable without `sudo`, `chmod 777`, or changing the fixed image user. macOS keeps its existing Docker Desktop behavior and does not receive this override. `rnaseq doctor` reports the selected policy and Linux/WSL mapping before execution.

## Versioned execution

The Salmon FASTQ route is pinned to nf-core/rnaseq 3.26.0 (revision `e7ca462`) with `-profile conda`. The HISAT2 route uses an exactly pinned Conda environment with HISAT2 2.2.3, SAMtools 1.21, Subread/featureCounts 2.0.6, FastQC 0.12.1, fastp 0.24.0, and MultiQC 1.33 (the same versions and builds as the historical Docker images). The downstream Conda lock pins Python, R, DESeq2, tximport, plotting, and enrichment packages exactly. Conda channel availability and network access for first-time environment creation remain operator responsibilities.

## Milestone A validation status (historical Docker runtime)

On 2026-09-05, the pinned production semantics fixture passed under Docker Desktop 29.7.2 on an arm64 host (the pinned amd64 process images ran under Docker emulation): SAMtools 1.21 `sha256:783c6646029a306ec5e4162009dc1a20d8f6c528f7c380e5b4affbf12d9112e5` and Subread/featureCounts 2.0.6 `sha256:114390a783c77f7739d86e474bedfa5a4e65309a2f71d4db430803fb04601f5d`. It verified single-end counting, paired fragments counted once, forward/reverse strands, ambiguous-overlap exclusion, both-mates and chimeric-fragment policy, and technical-lane merging. A primary `NH:i:2` alignment remained excluded after the production `samtools view -bh -F 0x900` transformation removed its secondary record; primary-only filtering therefore did not make a multimapper appear unique.

The successful raw L2 smoke run was `MILESTONE-A-RAW-L2-FINAL/20260905-115310+0800`; its first-party workflow identity was `nf-rna/hisat2_featurecounts` with source SHA-256 `e5408f74c1afd4c41bcd51e7b06868d7b7df9da4b77aff2ddb5483492810d05c`. It completed raw fastp/FastQC, HISAT2, original and published count-only BAM lineage, featureCounts, MultiQC, canonical matrix/sample map, `featurecounts_raw_counts` → `DESeqDataSetFromMatrix`, L1, L2, report, and delivery. The independent pretrimmed L1 branch also completed as `MILESTONE-A-PRETRIMMED-L1/20260905-115200+0800`. These immutable validation records are intentionally outside the source distribution.

The smoke intentionally selected no enrichment backend, so the documented no-enrichment L2 route was exercised; GO/KEGG enrichment was not claimed for synthetic gene identifiers. This is a synthetic execution check, not biological evidence.

## Single-end HISAT2 QC acceptance

On 2026-09-05, the separate interactive-wizard single-end QC run `WIZARD-HISAT2-QC-R3/20260905-144237+0800` passed in a local validation workspace that is not distributed with the source. It exercised raw preprocessing, HISAT2, featureCounts, FastQC/MultiQC, and QC delivery for the single-end synthetic fixture only. It did not exercise paired fragment accounting or downstream L1/L2/enrichment; those claims belong respectively to the paired acceptance record below and the Milestone A smoke records above.

## Paired-end HISAT2 QC acceptance

On 2026-09-05, a public interactive-wizard paired-end QC run passed as `PAIRED-HISAT2-QC-UAT-R3/20260905-150643+0800`. Its immutable validation record remains outside the source distribution. Its durable synthetic fixture is `tests/fixtures/hisat2_paired_raw_v2`; regenerate a new empty fixture root with `python tests/fixtures/hisat2_paired_raw/generate_fixture.py --root PATH/TO/EMPTY/hisat2_paired_raw_v2`, then register a matching prebuilt index or optionally run `rnaseq reference prepare-hisat2 PATH/TO/EMPTY/hisat2_paired_raw_v2/reference`. The fixture is forward stranded (HISAT2 `FR`; featureCounts `-s 1`), has PairAlpha lanes 001/002 and separate PairBeta, and defines counts before execution: PairAlpha `GeneA=3, GeneB=1`; PairBeta `GeneA=0, GeneB=2`.

The successful run retained synchronized R1/R2 pairs after raw fastp processing, trimmed the three short-insert adapter pairs (two reads and 66 bases per lane), mapped every lane at 100%, and exactly reproduced the fixture counts. featureCounts executed `-p --countReadPairs -B -C`; PairAlpha's eight alignment records represented four fragments and PairBeta's four records represented two fragments. Original and count-only BAMs passed `samtools quickcheck`; count-only BAMs had no `0x900` records and retained `NH:i:1` on every mapped record. Its MultiQC contains FastQC, fastp, HISAT2, and featureCounts, and the QC scope created no downstream L1/L2/enrichment artifacts. The run freezes workflow SHA-256 `ead01d1a6c84f20ad8d7f2c4943c97be20d748cc6adff0fe060520e852706d33` and the resolved container identities in provenance. This remains a software-contract fixture, not biological evidence.

## Local resource policy

`execution.max_cpus` and `execution.max_memory_gb` set the aggregate local
execution ceiling. They are execution preferences, not biological or
statistical settings. Each accepts a positive integer or `auto`:

- `auto` (the default for new projects and for projects without an
  `execution` section) sizes the ceiling from the machine that runs the case:
  the CPUs and memory available to the process (CPU affinity and cgroup v1/v2
  CPU/memory limits are respected) minus an OS reserve of
  `max(1, ceil(CPUs / 8))` CPUs and `max(4, ceil(15% of memory))` GiB. It never
  selects less than the 8 CPU / 12 GiB contract floor when the machine itself
  has that much, so small machines behave as the former fixed default did. On a
  28-CPU / 94-GiB WSL host, `auto` selects 24 CPUs / 79 GiB.
- An explicit integer is a deliberate limit and is used unchanged. It is only
  clamped, with a visible warning, when it exceeds the machine's capacity.
- If capacity cannot be detected, `auto` falls back to 8 CPUs / 12 GiB instead
  of blocking the run.

`rnaseq new` writes `auto` unless `--cpus`/`--memory-gb` are given; the
interactive wizard shows what `auto` selects on the current machine.
`rnaseq doctor PROJECT` reports host capacity, container-runtime capacity, the
requested project budget and how it was chosen, and the effective budget on the
machine actually running the project. Docker Desktop and WSL allocations constrain the effective
budget; native Linux uses the host ceiling rather than double-counting Docker's
repeated host values. Each run freezes that effective local Nextflow
configuration for upstream and downstream workflows.

The ceiling is total executor capacity, not a per-task request. HISAT2,
featureCounts, and other process directives keep their own declared requests,
and no task's CPU count (tool thread count) is changed. The single per-process
adjustment is a scheduling request for nf-core/rnaseq `SALMON_QUANT` when a
prebuilt Salmon index is used: its memory request is sized from the index on
disk (`ceil(index GiB x 1.15 + 2)`, capped at nf-core's own 36 GB) and frozen as
`frozen/nfcore.tuning.config`, which is passed only to nf-core/rnaseq. Salmon's
peak memory is dominated by the in-memory index, so this lets several samples
quantify concurrently instead of each request claiming the whole ceiling; Salmon
arguments, including `--threads`, are unchanged. Reference preparation remains
independent: use its explicit `--threads` option.

Every run records the decision in `provenance/run_provenance.yaml` and
`frozen/execution_manifest.yaml` (`resource_policy`): detected usable CPUs and
memory (with affinity and cgroup observations), the OS reserve, the selected and
effective ceilings, whether each was `auto` or explicit, any fallback, and the
Salmon scheduling request. A retry inherits the source run's frozen resource
configuration.

## Monitoring a run

```bash
rnaseq status PROJECT            # latest run: state, phase, task progress, resources
rnaseq status PROJECT --watch    # refresh every 10 s (--interval N, minimum 2) until it finishes; Ctrl-C stops watching only
rnaseq status PROJECT --case CASE-ID [--run RUN-ID]
rnaseq status PROJECT --all      # one summary block per recorded run
```

`rnaseq status` is read-only and does not need Nextflow to be running. It reads
the run's `run_state.json`, its own logs, the nf-core/rnaseq execution trace
(`upstream/nfcore_rnaseq/pipeline_info/`), the HISAT2 route's
`provenance/upstream.trace.txt`, the downstream trace, and frozen provenance.
Running task names come from Nextflow's plain log minus the tasks the trace has
already finished. A run is SUCCESS only when its durable state says so. A run
recorded as RUNNING whose `rnaseq` process no longer exists (checked by PID,
process start time and boot ID) is shown as INTERRUPTED (stale); one launched on
another host, or before per-run process records existed, is shown as RUNNING
(unverified).

Each run owns its execution log, `runs/CASE/RUN/logs/rnaseq.log`: start (PID,
host, command), selected resources, every state/phase transition, each Nextflow
launch and exit code, and the final outcome, including a traceback for an
unexpected error. Nextflow's output stays in `logs/upstream.*.log` and
`logs/downstream.*.log`. Two launches with the same case ID therefore never share
a log, and redirecting `rnaseq run` output to a file is not needed for
diagnostics. Interrupting `rnaseq run` (Ctrl-C, SIGTERM, or SIGHUP when not
ignored) stops Nextflow and records the run as INTERRUPTED; an unexpected error
records FAILED. None of these logs is copied into the delivery package.

`rnaseq retry PROJECT --retry-of CASE/RUN` accepts a source run recorded as
FAILED or INTERRUPTED; it refuses SUCCESS and any run still recorded as CREATED
or RUNNING. Retry reads only the durable `run_state.json`: a RUNNING run that
`rnaseq status` shows as INTERRUPTED (stale) has never recorded an end state and
is not retryable, because status never rewrites it; start a new run instead. The
source run keeps its recorded state, and the new attempt records it as
`retry_of.status`.

The default local contracts are intentionally conservative:

| class | CPUs | memory | time |
| --- | ---: | ---: | ---: |
| SMALL | 1 | 2 GiB | 2 h |
| MEDIUM | 4 | 8 GiB | 8 h |
| LARGE | 8 | 12 GiB | 12 h |

The default project ceiling is 8 CPUs, 12 GiB, and 12 hours. The frozen config sets Nextflow's aggregate local executor ceiling and `resourceLimits`, while preserving nf-core labels and first-party per-process requests. Independent sample-level and enrichment tasks may run concurrently when their summed requests fit that ceiling; dependency edges still sequence dependent stages. The first-party HISAT2 workflow binds fastp, HISAT2, supported SAMtools operations, and featureCounts threads to `task.cpus`. Each FASTQ run freezes one path-free local Nextflow configuration and its SHA-256, and supplies it to nf-core Salmon, HISAT2/featureCounts, and downstream. These are scheduling bounds only and do not alter quantification, filtering, DESeq2, thresholds, rankings, or GSEA calculations.

## Filesystem layout

Persistent immutable artifacts live beneath the project `runs/` directory. Operational Nextflow state lives in a local execution root outside that tree. On macOS the default is under the user's cache directory; on Linux it follows the standard cache location. Set `RNASEQ_EXECUTION_ROOT` to an absolute operator-owned scratch location if required.

Do not commit project run directories, workflow work directories, logs, delivery packages, references, indices, or production inputs. The repository's `.gitignore` and narrow `.dockerignore` enforce this default for new Git users.

## Determinism and external dependencies

Validation and planning are deterministic for unchanged inputs. Execution captures immutable contracts, checksums, versions, and selected resource profile. Timestamps, process scheduling, network conditions, and the external KEGG service are inherently runtime-dependent and are recorded rather than treated as deterministic results.
