"""Validate a GTF-derived transcriptome before Salmon indexing.

This utility writes the auditable transcript-ID contract artifact that a
Salmon declaration in reference_manifest.yaml references, without changing the
manifest or any project configuration.  The contract itself lives in
``rnaseq.references.transcript_id_contract`` and is the one the reference
builder and loader enforce: every transcriptome ID must exactly equal one GTF
transcript_id with one gene_id; no identifier is normalized.
"""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path

from rnaseq.references import (
    _load_local_reference_root,
    transcript_id_contract,
    transcript_id_contract_summary,
)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference-root", type=Path, required=True)
    parser.add_argument("--transcriptome", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    root = args.reference_root.resolve()
    transcriptome = args.transcriptome.resolve()
    output = args.output.resolve()
    # The Salmon declaration is what this artifact validates, so do not require it here.
    reference, _manifest = _load_local_reference_root(root, "reference_manifest.yaml", None, validate_salmon=False)
    payload = transcript_id_contract(
        identity={
            "species": reference.species, "provider": reference.provider, "release": reference.release,
            "assembly": reference.assembly, "assembly_patch": reference.assembly_patch,
        },
        genome_fasta=reference.genome_fasta,
        annotation_gtf=reference.annotation_gtf,
        transcriptome=transcriptome,
        transcriptome_relative_path=transcriptome.relative_to(root).as_posix(),
    )
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.with_name(f".{output.name}.{os.getpid()}.tmp")
    temporary.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    temporary.replace(output)
    print(transcript_id_contract_summary(payload))
    return 0 if payload["status"] == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())
