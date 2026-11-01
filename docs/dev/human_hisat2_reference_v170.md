# nf-rna v1.7.0 — Human HISAT2 reference resolution

Branch `feature/v1.7.0`, from `a290564`. Reported problem: in *New Project*,
Mouse + HISAT2 showed a reference path, Human + HISAT2 an empty one.

## 1. Root cause

**Primary (configuration): the Human reference was prepared and validated but
never registered.** *New Project* offers managed references only from this
workstation's registry, `/home/brian/.config/nf-rna/references.yaml`. Both the
GUI on port 8501 and the CLI read that file (`HOME=/home/brian`, no
`XDG_CONFIG_HOME`). It is read through
`references.compatible_registered_references(species, backend)`, the same
function `rnaseq new` uses. Until 2026-11-01 the registry had a single entry,
`/home/brian/reference/mus_musculus/ensembl_116_GRCm39`, so Human had no
registered reference for either backend. Human + HISAT2 then had only
*Other (advanced)*, with empty fields. Human + Salmon also had no managed
reference; only iGenomes and *Other* were offered there.

A reference enters the registry only through `rnaseq reference register`
(`cli.py`, `register_local_reference`). Building or runtime-validating a
reference does not register it. The Human manifest was last revised on
2026-10-18 (HISAT2 runtime validation PASS); the registry was last written
on 2026-10-30, with Mouse only.

**Secondary (GUI): a stale selection hid the Mouse reference after a species
switch.** The reference selectbox had a single key. With Human selected,
*Other (advanced)* was the only option and became the value. After switching
to Mouse, Streamlit kept *Other*, because that label is also a Mouse option,
so the registered Mouse path was hidden. This was reproduced in the browser
before the fix (`browser_checks_before_registration.keyfix-pre.json`: Mouse +
HISAT2 showed "Other (advanced)").

Ruled out:
- the GUI reading another registry or root;
- the Human reference being absent, incomplete or incompatible (§2);
- a defect in `compatible_registered_references`. Once Human is registered,
  it returns `/home/brian/reference/homo_sapiens/ensembl_116_GRCh38` for both
  backends, and switching species or backend selects the right reference.
  These are 4 of the new tests, which also pass on `a290564`.

## 2. Human reference inventory (actual files)

Root `/home/brian/reference/homo_sapiens/ensembl_116_GRCh38`, manifest
`reference_manifest.yaml` (schema 1.1, purpose production, Ensembl 116,
GRCh38.p14, SHA-256 `a5dcebd8…3200`). Older manifests kept beside it:
`reference_manifest.2cd2257c1ff9.yaml`, `reference_manifest.c4ce70179855.yaml`,
`reference_manifest.yaml.bak`.

| Asset | Path | Check |
|---|---|---|
| Genome FASTA | `fasta/Homo_sapiens.GRCh38.dna.primary_assembly.fa` (3,151,425,851 B) | SHA-256 `1e74081a…ef3b` = manifest |
| Annotation GTF | `annotation/Homo_sapiens.GRCh38.116.gtf` (4,659,761,880 B) | SHA-256 `474626b4…d9e7` = manifest |
| HISAT2 index | `hisat2/index/genome.{1..8}.ht2` (complete numbered family) | manifest `hisat2.status: built`, HISAT2 2.2.3, strategy `genome_only_runtime_splicesites`, `genome_fasta_sha256` = FASTA, `source_gtf_sha256` = GTF |
| Splice sites | `hisat2/annotation/splice_sites.txt` | SHA-256 `2d18305d…4118` = manifest |
| Runtime compatibility | manifest `hisat2.runtime_compatibility: validated` | validated 2026-10-18 by direct HISAT2 2.2.3 alignment of PRJEB3366 (4 runs × 1 M pairs): overall alignment 95.1–98.1 %, featureCounts assigned 43–76 % (ERR188428 low from rDNA), annotated junction fraction 0.89–0.99; evidence bundle `/home/brian/nf-rna-work/human-runtime-validation/evidence/hisat2_runtime_validation.tar.gz` |
| Salmon | `salmon/gtf_derived/index` (decoy-aware, Salmon 1.10.3) | separate from the HISAT2 index; never used for HISAT2 |

Index content (`hisat2-inspect` from `nf-rna-reference-builder`):
194 sequences, 3,099,750,718 bp, headers `chromosome:GRCh38:…`. The index's
sequence names equal the 194 FASTA headers exactly, and all 70 GTF seqnames
are among them. Name lists are in
`/home/brian/nf-rna-validation/human-hisat2-ref/`.

Loaded with the registry's own loader (`_load_local_reference_root`):
species Homo sapiens, production, `hisat2_index` set,
`hisat2_runtime_ready = True`.

**Verdict:** a complete, compatible, production Human HISAT2 reference
exists. Nothing was missing, so nothing was downloaded or rebuilt.

## 3. Mouse comparison

| | Human | Mouse |
|---|---|---|
| Root | `/home/brian/reference/homo_sapiens/ensembl_116_GRCh38` | `/home/brian/reference/mus_musculus/ensembl_116_GRCm39` |
| Identity | Ensembl 116 GRCh38.p14, production | Ensembl 116 GRCm39, production |
| HISAT2 index | `hisat2/index/genome.*.ht2`, runtime-ready | same layout, runtime-ready |
| Index names = FASTA names | yes (194) | yes (61) |
| Registered before the fix | **no** | yes (manifest `8bce2065…4a7b`) |
| `compatible_registered_references(…, "hisat2_featurecounts")` before | `[]` | Mouse root |

(`ensembl_116_GRCm39_broken_chr1` beside the Mouse root is a deliberate
fault fixture; it is not registered.)

## 4. Fix

1. **Registration (workstation configuration, existing command).**
   ```bash
   rnaseq reference register /home/brian/reference/homo_sapiens/ensembl_116_GRCh38
   ```
   Run on 2026-11-01. The registry gained only the Human entry (identity and
   `manifest_sha256 a5dcebd8…`); the Mouse entry is unchanged. The previous
   file is kept as
   `/home/brian/nf-rna-validation/human-hisat2-ref/references.yaml.before`.
   No file under `/home/brian/reference` changed (manifest SHA-256 identical,
   and nothing newer than the backup). To undo, copy that file back.
2. **GUI (`src/rnaseq/gui/pages.py`, `services.py`).**
   - The reference selectbox is keyed per species and backend, so a choice
     made for one combination is never carried over to another.
   - When no registered reference fits, the page says **"<species> <backend>
     reference not prepared"**. It lists the registered references of that
     species that cannot be used and why (no HISAT2 index, with
     `rnaseq reference prepare-hisat2 ROOT`; HISAT2 runtime compatibility not
     validated; no Salmon index, with `rnaseq reference prepare ROOT`), and
     shows the registry path and `rnaseq reference register`. No path is
     filled in.
   - Before *Create project*, a managed reference (registered or typed under
     *Other*) is loaded again with the same loader and checked against the
     project's species and backend (`services.check_project_reference`). The
     check is the same rule as the registry: same species, and the backend's
     own index, so a Salmon index never stands in for a HISAT2 index. It also
     runs live under *Other*, so a Mouse root typed for a Human project shows
     an error and is refused.
   - At run time the project validator still re-checks the species
     (`expected_species`) and HISAT2 runtime readiness
     (`production_hisat2_runtime_unvalidated`), unchanged.

Nothing in HISAT2, reference preparation, annotation mapping, assembly
identities, the scientific engine or any run was changed.

## 5. Test evidence

`tests/test_gui_references.py` (10 tests; synthetic references and registry
in a temporary `XDG_CONFIG_HOME`):

| Case | Test |
|---|---|
| 1–4 Human/Mouse × Salmon/HISAT2 | `test_each_species_and_backend_offers_only_its_own_registered_reference` |
| 5–7 Human→Mouse, Mouse→Human, Salmon→HISAT2 | `test_switching_species_and_backend_never_keeps_an_incompatible_selection`, `test_other_chosen_for_lack_of_a_human_reference_is_not_carried_over_to_mouse` |
| 8 Missing reference | `test_a_missing_human_hisat2_reference_is_reported_not_filled_in` (also: a Mouse root under *Other* refused for Human) |
| 9 Incomplete reference | `test_an_incomplete_human_reference_names_what_is_missing` (Salmon-only; HISAT2 not runtime-validated), `test_a_salmon_index_is_never_taken_for_a_hisat2_index` |
| 10 Valid custom reference | `test_registered_and_valid_custom_references_create_valid_human_hisat2_projects` (registered and typed under *Other*; project.yaml and validator checked) |
| creation / preflight | `test_the_creation_check_refuses_another_species_and_a_broken_root`, `test_the_validator_still_refuses_a_mismatched_reference_at_run_time` |

On `a290564`, 5 of these fail: the "not prepared" report, the incomplete
reference, the creation checks and the carried-over *Other*. The 4 switching
and valid-reference tests pass there, which shows the registry rule itself was
correct. GUI, reference and release-metadata tests: 134 passed.

Browser (headless Chrome, `scripts/dev/gui_reference_acceptance.py`, GUI on
port 8703 with the real registry; evidence in
`/home/brian/nf-rna-validation/human-hisat2-ref/`): **11/11**. All four
species/backend combinations show their own root. Mouse→Human and
Human→Mouse switch correctly, and Salmon→HISAT2 drops iGenomes. A Human +
HISAT2 FASTQ scaffold was created with
`reference.root = /home/brian/reference/homo_sapiens/ensembl_116_GRCh38`, and
the project validator loads it as GRCh38.p14 with HISAT2 runtime-ready and no
reference error. The only findings are the expected empty-scaffold ones:
`missing_fastq_files`, `missing_metadata_rows`. With a copy of the old
registry (port 8704), Human + HISAT2 shows "Homo sapiens HISAT2 reference not
prepared" with no path, and Mouse + HISAT2 keeps the Mouse reference
(`browser_checks_before_registration.json`). Screenshots: `shots/R0`, `R1`,
`R2`.

## 6. Remaining limitations

- The registry is per workstation. Another machine or user account needs its
  own `rnaseq reference register`. The GUI does not scan for unregistered
  reference folders.
- *Other → Genome FASTA + annotation GTF* (custom) references are checked only
  by the project validator, not before creation.
- The live *Other* check caches its result per manifest mtime; a manifest
  changed within the same second is re-checked only at creation.
- This review did not re-run the HISAT2 runtime validation of 2026-10-18; it
  relies on the recorded evidence and on the checksums verified here.

## 7. Ready for a real Human HISAT2 FASTQ run?

**Yes, as far as the reference is concerned.** The registered Human reference
is complete (checksums verified), GRCh38-consistent across FASTA, GTF and
index, runtime-validated for HISAT2 2.2.3, and accepted by the project
validator. No real FASTQ run was made in this review. The first real run
should follow the checklist in `gui_hardening_v170.md` §6, with the
strandedness chosen explicitly for HISAT2 + featureCounts.
