"""Create the minimal transcript FASTA matched to the pinned smoke-test GTF."""

from __future__ import annotations

import gzip
import re
from pathlib import Path


ROOT = Path(__file__).parent / "reference"
GTF = ROOT / "genes.gtf.gz"
SOURCE = ROOT / "transcriptome.fasta"
OUTPUT = ROOT / "transcriptome.gtf_matched.fasta"


def main() -> None:
    with gzip.open(GTF, "rt", encoding="utf-8") as handle:
        identifiers = set(re.findall(r'transcript_id "([^"]+)"', handle.read()))
    retained = 0
    total = 0
    keep = False
    with SOURCE.open(encoding="utf-8") as source, OUTPUT.open("w", encoding="utf-8", newline="\n") as output:
        for line in source:
            if line.startswith(">"):
                total += 1
                keep = line[1:].split()[0] in identifiers
                retained += int(keep)
            if keep:
                output.write(line)
    if retained != len(identifiers):
        raise RuntimeError(f"Retained {retained} of {len(identifiers)} GTF transcript records.")
    print(f"Retained {retained} of {total} source transcript records.")


if __name__ == "__main__":
    main()
