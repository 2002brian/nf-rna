nextflow.enable.dsl=2

/* First-party conventional bulk-RNA-seq alignment/counting route.
 *
 * Each source row is one technical lane.  The lane identity is preserved to
 * alignment, lanes are grouped only by the declared sample ID, then exactly
 * one coordinate-sorted BAM is counted per biological sample.
 */

params.input = null
params.outdir = null
params.fasta = null
params.gtf = null
params.hisat2_index = null
params.hisat2_index_basename = 'genome'
params.hisat2_splice_sites = null
params.hisat2_use_runtime_splices = false
params.layout = null
params.strandedness = null
params.pretrimmed = false
params.assembly_script = null
params.hisat2_container = 'quay.io/biocontainers/hisat2:2.2.3--h8471819_0'
params.samtools_container = 'quay.io/biocontainers/samtools:1.21--h50ea8bc_0'
params.subread_container = 'quay.io/biocontainers/subread:2.0.6--he4a0461_2'
params.fastp_container = 'quay.io/biocontainers/fastp:0.24.0--h125f33a_0'
params.fastqc_container = 'quay.io/biocontainers/fastqc:0.12.1--hdfd78af_0'
params.multiqc_container = 'community.wave.seqera.io/library/multiqc:1.33--ee7739d47738383b'

def hisatStrand(value, layout) {
    if (value == 'unstranded') return ''
    if (layout == 'paired_end') return value == 'forward' ? '--rna-strandness FR' : '--rna-strandness RF'
    return value == 'forward' ? '--rna-strandness F' : '--rna-strandness R'
}

def featureCountsStrand(value) { ['unstranded':'0', 'forward':'1', 'reverse':'2'][value] }

process FASTQC_RAW {
    tag { "raw ${sample}" }
    container params.fastqc_container
    cpus 2
    memory '2 GB'
    publishDir "${params.outdir}/qc/fastqc/raw", mode: 'copy', overwrite: false
    input:
    tuple val(sample), val(strandedness), path(reads)
    output:
    path 'raw_*_fastqc.html', emit: html
    path 'raw_*_fastqc.zip', emit: zip
    script:
    """
    fastqc --threads ${task.cpus} ${reads}
    for f in *_fastqc.html *_fastqc.zip; do mv "\$f" "raw_\$f"; done
    """
}

process FASTP_PREPARE {
    tag { sample }
    container params.fastp_container
    cpus 4
    memory '6 GB'
    publishDir "${params.outdir}/qc/fastp", mode: 'copy', overwrite: false
    input:
    tuple val(sample), val(strandedness), path(reads)
    output:
    tuple val(sample), val(strandedness), path('prepared/*.fastq.gz'), emit: prepared
    path 'fastp_*.html', optional: true, emit: fastp_html
    path 'fastp_*.json', optional: true, emit: fastp_json
    script:
    def paired = reads.size() == 2
    def pretrimmed = params.pretrimmed == true || params.pretrimmed == 'true'
    def r1 = reads[0].baseName.replaceFirst(/\\.fastq$/, '')
    def r2 = paired ? reads[1].baseName.replaceFirst(/\\.fastq$/, '') : null
    if (paired) {
        """
        mkdir -p prepared
        if ${pretrimmed}; then cp ${reads[0]} prepared/${r1}.fastq.gz; cp ${reads[1]} prepared/${r2}.fastq.gz
        else fastp --thread ${task.cpus} -i ${reads[0]} -I ${reads[1]} --detect_adapter_for_pe -o prepared/${r1}.fastq.gz -O prepared/${r2}.fastq.gz -h fastp_${r1}.html -j fastp_${r1}.json
        fi
        """
    } else {
        """
        mkdir -p prepared
        if ${pretrimmed}; then cp ${reads[0]} prepared/${r1}.fastq.gz
        else fastp --thread ${task.cpus} -i ${reads[0]} -o prepared/${r1}.fastq.gz -h fastp_${r1}.html -j fastp_${r1}.json
        fi
        """
    }
}

process HISAT2_ALIGN {
    tag { sample }
    container params.hisat2_container
    cpus 4
    memory '6 GB'
    publishDir "${params.outdir}/alignment/lane_summaries", mode: 'copy', overwrite: false
    input:
    tuple val(sample), val(strandedness), path(reads)
    path index
    path known_splices
    output:
    tuple val(sample), val(strandedness), path("*.lane.sam"), path("*.hisat2.summary"), emit: aligned
    script:
    def paired = reads.size() == 2
    def strand = hisatStrand(strandedness, params.layout)
    def inputs = paired ? "-1 ${reads[0]} -2 ${reads[1]}" : "-U ${reads[0]}"
    def lane = reads[0].baseName.replaceFirst(/\\.fastq$/, '')
    def spliceArgument = (params.hisat2_use_runtime_splices == true || params.hisat2_use_runtime_splices == 'true') ? "--known-splicesite-infile ${known_splices}" : ''
    """
    hisat2 -p ${task.cpus} -x ${index}/${params.hisat2_index_basename} ${inputs} ${strand} ${spliceArgument} --summary-file ${lane}.hisat2.summary -S ${lane}.lane.sam
    """
}

process FASTQC_PROCESSED {
    tag { "processed ${sample}" }
    container params.fastqc_container
    cpus 2
    memory '2 GB'
    publishDir "${params.outdir}/qc/fastqc/processed", mode: 'copy', overwrite: false
    input:
    tuple val(sample), val(strandedness), path(reads)
    output:
    path 'processed_*_fastqc.html', emit: html
    path 'processed_*_fastqc.zip', emit: zip
    script:
    """
    fastqc --threads ${task.cpus} ${reads}
    for f in *_fastqc.html *_fastqc.zip; do mv "\$f" "processed_\$f"; done
    """
}

process SORT_LANE_BAM {
    tag { sample }
    container params.samtools_container
    cpus 4
    memory '6 GB'
    input:
    tuple val(sample), val(strandedness), path(sam), path(summary)
    output:
    tuple val(sample), val(strandedness), path("*.lane.bam"), path(summary), emit: sorted
    script:
    """
    samtools sort -@ ${task.cpus} -o ${sam.baseName}.bam ${sam}
    samtools quickcheck -v ${sam.baseName}.bam
    """
}

process MERGE_AND_INDEX {
    tag { sample }
    container params.samtools_container
    cpus 4
    memory '6 GB'
    publishDir "${params.outdir}/bam", mode: 'copy', overwrite: false
    input:
    tuple val(sample), val(strandedness_values), path(bams), path(summaries)
    output:
    tuple val(sample), path("${sample}.bam"), path("${sample}.bam.bai"), path("${sample}.flagstat.txt"), emit: merged
    script:
    def uniqueStrandedness = strandedness_values.toSet()
    if (uniqueStrandedness.size() != 1) {
        error "All technical lanes for ${sample} must declare the same strandedness; observed: ${uniqueStrandedness}"
    }
    """
    samtools merge -@ ${task.cpus} -o ${sample}.merged.bam ${bams}
    samtools sort -@ ${task.cpus} -o ${sample}.bam ${sample}.merged.bam
    samtools index -@ ${task.cpus} ${sample}.bam
    samtools quickcheck -v ${sample}.bam
    samtools flagstat ${sample}.bam > ${sample}.flagstat.txt
    """
}

process FEATURECOUNTS {
    tag { sample }
    container params.subread_container
    cpus 4
    memory '6 GB'
    publishDir "${params.outdir}/counts/per_sample", mode: 'copy', overwrite: false
    input:
    tuple val(sample), path(bam), path(bai), path(flagstat), path(gtf)
    output:
    tuple val(sample), path("${sample}.counts.txt"), path("${sample}.counts.txt.summary"), emit: counted
    script:
    def paired = params.layout == 'paired_end'
    def pairArgs = paired ? '-p --countReadPairs -B -C' : ''
    def strand = featureCountsStrand(params.strandedness)
    """
    featureCounts -T ${task.cpus} -a ${gtf} -o ${sample}.counts.txt -t exon -g gene_id -s ${strand} -Q 0 --primary ${pairArgs} ${bam}
    """
}

process PREPARE_COUNT_BAM {
    tag { sample }
    container params.samtools_container
    cpus 2
    memory '3 GB'
    publishDir "${params.outdir}/bam/count_only", mode: 'copy', overwrite: false
    input:
    tuple val(sample), path(bam), path(bai), path(flagstat)
    output:
    tuple val(sample), path("${sample}.countable.bam"), path("${sample}.countable.bam.bai"), path(flagstat), emit: countable
    script:
    """
    # Retain the original coordinate-sorted diagnostic BAM.  This separate
    # count input excludes secondary (0x100) and supplementary (0x800)
    # records without modifying NH or other alignment tags.
    samtools view -@ ${task.cpus} -bh -F 0x900 -o ${sample}.countable.bam ${bam}
    samtools index -@ ${task.cpus} ${sample}.countable.bam
    samtools quickcheck -v ${sample}.countable.bam
    """
}

process ASSEMBLE_COUNTS {
    tag 'canonical featureCounts matrix'
    container 'quay.io/biocontainers/python:3.10.4'
    cpus 1
    memory '2 GB'
    publishDir "${params.outdir}/counts", mode: 'copy', overwrite: false
    input:
    path counts
    path assembly_script
    output:
    path 'canonical_counts.csv'
    path 'sample_map.csv'
    script:
    def args = counts.collect { "--sample ${it.baseName.replace('.counts','')}=${it.name}" }.join(' ')
    """
    python ${assembly_script} ${args} --out canonical_counts.csv
    printf 'sample_id,featurecounts_file\\n' > sample_map.csv
    for f in ${counts}; do printf '%s,%s\\n' "\${f%.counts.txt}" "\$f" >> sample_map.csv; done
    """
}

process MULTIQC {
    tag 'MultiQC'
    container params.multiqc_container
    cpus 1
    memory '2 GB'
    publishDir "${params.outdir}/multiqc", mode: 'copy', overwrite: false
    input:
    path reports
    output:
    path 'multiqc_report.html'
    script:
    """
    multiqc --force --filename multiqc_report.html .
    """
}

workflow {
    if (!params.input || !params.outdir || !params.fasta || !params.gtf || !params.hisat2_index || !params.assembly_script || !(params.layout in ['paired_end','single_end']) || !(params.strandedness in ['unstranded','forward','reverse']) || ((params.hisat2_use_runtime_splices == true || params.hisat2_use_runtime_splices == 'true') && !params.hisat2_splice_sites)) {
        error 'Specify --input, --outdir, --fasta, --gtf, --hisat2_index, --assembly_script, supported --layout and explicit --strandedness; runtime splice mode also requires --hisat2_splice_sites.'
    }
    reads = Channel.fromPath(params.input).splitCsv(header: true).map { row ->
        def files = row.fastq_2 ? [file(row.fastq_1), file(row.fastq_2)] : [file(row.fastq_1)]
        tuple(row.sample, row.strandedness, files)
    }
    FASTQC_RAW(reads)
    FASTP_PREPARE(reads)
    FASTQC_PROCESSED(FASTP_PREPARE.out.prepared)
    knownSplices = params.hisat2_splice_sites ? file(params.hisat2_splice_sites) : file(params.gtf)
    HISAT2_ALIGN(FASTP_PREPARE.out.prepared, Channel.value(file(params.hisat2_index)), Channel.value(knownSplices))
    SORT_LANE_BAM(HISAT2_ALIGN.out.aligned)
    MERGE_AND_INDEX(SORT_LANE_BAM.out.sorted.groupTuple())
    PREPARE_COUNT_BAM(MERGE_AND_INDEX.out.merged)
    FEATURECOUNTS(PREPARE_COUNT_BAM.out.countable.combine(Channel.value(file(params.gtf))).map { sample, bam, bai, flagstat, gtf -> tuple(sample, bam, bai, flagstat, gtf) })
    ASSEMBLE_COUNTS(FEATURECOUNTS.out.counted.map { sample, counts, summary -> counts }.collect(), Channel.value(file(params.assembly_script)))
    multiqc_reports = FEATURECOUNTS.out.counted
        .map { sample, counts, summary -> summary }
        .mix(FASTP_PREPARE.out.fastp_json)
        .mix(HISAT2_ALIGN.out.aligned.map { sample, strandedness, sam, summary -> summary })
        .mix(FASTQC_RAW.out.html)
        .mix(FASTQC_RAW.out.zip)
        .mix(FASTQC_PROCESSED.out.html)
        .mix(FASTQC_PROCESSED.out.zip)
    MULTIQC(multiqc_reports.collect())
}
