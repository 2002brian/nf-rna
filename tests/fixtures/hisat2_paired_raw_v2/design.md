# Paired-end HISAT2 UAT fixture

The library is forward stranded: R1 is on the transcript strand, so HISAT2 uses `--rna-strandness FR` and featureCounts uses `-s 1`.

`PairAlpha` has two technical lanes from one library: lane 001 has two GeneA fragments and lane 002 has one GeneA plus one GeneB fragment. `PairBeta` is a separate sample with two GeneB fragments. All six fragments are uniquely mappable. Three ordinary fragments have 180 bp inserts; three adapter-bearing fragments have 60 bp inserts, allowing paired-end overlap detection to trim standard Illumina adapters while retaining 60 high-quality biological bases.

Expected counts are independently specified by this construction, before execution: PairAlpha = GeneA 3 / GeneB 1; PairBeta = GeneA 0 / GeneB 2. A fragment contributes once under `-p --countReadPairs -B -C`; BAM alignment-record counts are therefore twice the eligible fragment count.
