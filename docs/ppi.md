# Experimental offline PPI analysis

The `v1.8.0` development branch provides candidate preparation, explicit-namespace
STRING mapping and induced-network extraction. These developer entry points are
not yet wired into `rnaseq run`, its immutable run lifecycle, Nextflow, GUI or
delivery reports. The package version remains 1.7.0 until release preparation.

## Inputs

Prepare candidates from an existing L2 contrast, outside its immutable run:

```bash
PYTHONPATH=src python -m rnaseq.ppi \
  --contrast-dir /data/project/runs/CASE/RUN/downstream/l2/contrasts/B_vs_A \
  --outdir /data/ppi/candidates \
  --organism 'Homo sapiens' --input-id-type ENSEMBL
```

Use original L2 files, not a GUI-filtered export. The selected genes must agree
with the recorded DEG definition. Raw LFC determines up/down direction.

## Import a local STRING snapshot

The mapping input is an explicit, normalized TSV with these exact columns:

```text
query_id	input_id_type	string_id	taxon_id
ENSG00000141510	ENSEMBL	9606.ENSP00000269305	9606
```

This is **not** the native STRING aliases file or the API response format.
The operator must prepare the mapping using the matching STRING version and
namespace, retaining all matches. Use ENSEMBL gene IDs without version suffixes,
ENTREZID or SYMBOL exactly as declared. Do not label mixed aliases as one ID type
or select the first API match. Species and ID type are never inferred.

Links use the STRING download format: whitespace-separated `protein1 protein2`
followed by integer evidence-score columns, including `combined_score` (0–1000).
Gzipped links are supported and streamed. API scores on the 0–1 scale must first
be explicitly converted; this importer rejects fractional scores.

```bash
PYTHONPATH=src python -m rnaseq.ppi_network snapshot \
  --mapping /data/string/mapping.tsv \
  --links /data/string/9606.protein.links.full.v12.5.txt.gz \
  --outdir /data/ppi/string-snapshot \
  --version 12.5 --taxon-id 9606 --network-type functional \
  --source-url https://version-12-5.string-db.org/
```

The version above is an example, not an automatically chosen release. Declare
the actual downloaded release, URLs and network type. `physical` is also
supported as a declaration, with the corresponding physical links as input.
Version, network type, mapping correctness and table completeness are operator
declarations, not independently authenticated. `scope: provided_tables` records
that the data may be a subset. A missing mapping or edge means absent from these
provided tables, not necessarily absent from STRING.

The importer copies and validates both resources, then writes a manifest with
version, taxon, network type, source URLs, import time, checksums and sizes.

## Extract a network

```bash
PYTHONPATH=src python -m rnaseq.ppi_network network \
  --candidates /data/ppi/candidates \
  --snapshot /data/ppi/string-snapshot \
  --outdir /data/ppi/network-all \
  --selection all --required-score 700 --minimum-mapping-rate 0.5
```

`selection` is `all`, `up` or `down`. Edges with combined score **greater than or
equal to** `required-score` are retained. The default score is 700. The default
minimum mapping rate is 0 (report without blocking); the example explicitly
requires 0.5. These are configurable policies, not validated universal cutoffs.
Mapping rate counts original selected genes with exactly one protein match.
Empty input has a null mapping rate and is distinct from an unmapped input.

Ambiguous queries are excluded from network construction; all alternatives are
retained in the mapping audit. Multiple unambiguous genes mapping to one protein
are preserved as a list of original gene IDs. No LFC averaging is performed.
No neighbors are added. Isolated nodes are retained. Identical undirected
reverse edges are deduplicated; conflicting scores fail. Self-loops are excluded
and counted. Evidence columns are copied without recomputing combined scores.

Outputs:

| File | Meaning |
| --- | --- |
| `mapping.tsv` | Original gene statistics, protein matches and mapping status |
| `nodes.tsv` | Unique protein IDs, original gene IDs as JSON lists, unweighted degree |
| `edges.tsv` | Induced undirected edges and available evidence scores |
| `selected_candidates.tsv` | Exact selected input bytes |
| `candidate_summary.json` | Candidate provenance |
| `snapshot_manifest.json` | Resource provenance |
| `summary.json` | Parameters, mapping rate, counts, status and output hashes |

Statuses are `NO_CANDIDATES`, `NO_MAPPED_PROTEINS`, `NO_EDGES`, `COMPLETE`.
The summary is published last as the completion marker. Existing outputs are
never overwritten; after abrupt termination, an incomplete destination may need
operator inspection. This helper does not provide the production run lifecycle.

The network does not constitute an enrichment test or causal evidence. Functional
associations do not establish physical binding; a declared physical network also
does not establish experimental validation for every edge. There are no hub
significance scores, pathway tests or figures in this stage.

Reference: [STRING API documentation](https://string-db.org/help/api/) recommends
mapping identifiers before querying networks and pinning the release for
production use. For large analyses, use downloaded data rather than repeated API
queries.
