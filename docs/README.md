# nf-rna documentation

The project homepage is the [README](../README.md) ([繁體中文](../README_zh-TW.md)). This page indexes the rest of the documentation for the current release, v1.7.0.

## User documentation

| Page | Contents |
|---|---|
| [Quick Start](quickstart.md) | installation, project creation, FASTQ naming, metadata templates, `rnaseq design`, upstream reuse, monitoring, cleanup, references |
| [Browser GUI](gui.md) | installing and starting `rnaseq gui`, `gui.yaml`, every page, figure downloads, GUI limitations, user service |
| [Runtime and reproducibility](runtime.md) | supported runtime, managed references, local resource policy, monitoring, filesystem layout, published HISAT2 outputs |
| [Scientific contract](scientific_contract.md) | input integrity, L1, L2 differential expression, DEG definition, GO/KEGG ORA and GSEA, output schemas, provenance |
| [Figure standard](figure-standard.md) | figure presets, profiles, colour scales, typography and the output contract |
| [CHANGELOG](../CHANGELOG.md) | changes and known limitations of every release |
| [User service example](systemd/nf-rna-gui.service) | `systemd --user` unit for the GUI (not installed automatically) |

## Developer documentation

| Page | Contents |
|---|---|
| [Development installation](development.md) | source checkout, editable install, development runs, tests |
| [Architecture](architecture.md) | project contracts, execution ownership, run lifecycle, reproducibility boundaries |
| [Scientific CI](dev/scientific_ci.md) | the scientific release gate and its [mutation evidence](dev/scientific_ci_mutations.md) |
| [v1.6.1 scientific corrections](dev/scientific_corrections_v161.md) and [reliability fixes](dev/scientific_reliability_v161.md) | the result-changing corrections of v1.6.1 and their evidence |
| Release notes: [v1.7.0](dev/release_notes_v1.7.0.md), [v1.6.1](dev/release_notes_v1.6.1.md), [v1.6.0](dev/release_notes_v1.6.0.md) | the GitHub Release texts |

Scientific validation evidence of earlier releases is in [`validation/`](../validation/).

## Historical archive

| Page | Contents |
|---|---|
| [v1.3.0 README](archive/README_v1.3.0.md) | the README of v1.3.0, archived unchanged |

The documentation of any release is at its Git tag, for example
<https://github.com/2002brian/nf-rna/tree/v1.6.1>.
