# nf-rna downstream runtime contract

`nf-rna-downstream.yml` is the reviewed, human-readable dependency source for
first-party downstream analysis only. It intentionally excludes nf-core
alignment/QC tools, HISAT2/featureCounts tools, development dependencies, and
Docker-only `procps-ng`.

It is not a production runtime by itself. Each supported platform must have a
reviewed exact Conda lock in `locks/`, with its SHA-256 recorded in the run's
runtime provenance. Locks retain Conda's resolved channel, subdir, version,
and build for every package. A lock is produced only from a successful platform
resolution of this source; never hand-edit or synthesize a lock.

The runtime identity has two independent immutable parts:

1. The platform lock SHA-256 identifies Conda packages.
2. A non-editable nf-rna wheel SHA-256 identifies the exact Python package and
   its bundled `rnaseq/r/*.R` scripts. Install that wheel with `pip --no-deps`
   after creating the locked environment; do not use an editable install or
   the caller's active environment.

The wheel's origin is recorded as `wheel.origin`:

- `installed-distribution` (production): the normally installed nf-rna is
  verified against its pip `RECORD` hashes and re-packed, with the Python
  standard library only, into a deterministic uncompressed wheel. Identical
  installed content always yields the identical wheel SHA-256; no Git checkout,
  build tool, or network access is needed.
- `source-checkout-build` (development): a source checkout builds the wheel with
  `python -m build` from the `dev` extra.

A cached prefix whose recorded wheel SHA-256 differs from the requested wheel
is rejected, never updated in place.

For each qualified release, provenance must record the platform, lock filename
and SHA-256, wheel filename and SHA-256, nf-rna version/source revision, and a
digest of the bundled R-script inventory. The wheel digest remains authoritative
because the R scripts are package data shipped with the wheel.

No platform is supported merely because this source resolves. Native smoke
execution remains the platform-qualification gate.
