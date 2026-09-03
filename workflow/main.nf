nextflow.enable.dsl=2

/*
 * First-party downstream orchestration. Statistical work remains in the
 * versioned R backends; Python only materializes task-local JSON contracts.
 */

params.contract = null
params.inputs = null
params.outdir = null
params.r_scripts = '/opt/rnaseq/src/rnaseq/r'
params.enrichment = ''
params.analysis_level = null

process L1_ANALYSIS {
    tag 'L1 expression QC'
    publishDir params.outdir, mode: 'copy', overwrite: false
    stageInMode 'copy'
    input:
    path contract
    path inputs
    output:
    path 'l1'
    script:
    """
    python -m rnaseq.workflow_support l1-config --contract $contract --inputs $inputs --out l1-config.json
    Rscript ${params.r_scripts}/l1_analysis.R --config l1-config.json
    """
}

process L2_ANALYSIS {
    tag 'L2 DESeq2'
    publishDir params.outdir, mode: 'copy', overwrite: false
    stageInMode 'copy'
    input:
    path l1
    path contract
    path inputs
    output:
    path 'l2'
    script:
    """
    python -m rnaseq.workflow_support l2-config --contract $contract --inputs $inputs --l1 $l1 --out l2-config.json
    Rscript ${params.r_scripts}/l2_analysis.R --config l2-config.json
    """
}

process TECHNICAL_REPORT {
    tag 'HTML technical report'
    publishDir params.outdir, mode: 'copy', overwrite: false
    stageInMode 'copy'
    input:
    path l1
    path l2
    path enrichment_dirs
    path contract
    path inputs
    output:
    path 'report'
    script:
    """
    python -m rnaseq.workflow_support report --contract $contract --inputs $inputs --l1 $l1 --l2 $l2 --enrichment $enrichment_dirs --out report
    """
}

process TECHNICAL_REPORT_NO_ENRICHMENT {
    tag 'HTML technical report (no enrichment selected)'
    publishDir params.outdir, mode: 'copy', overwrite: false
    stageInMode 'copy'
    input:
    path l1
    path l2
    path contract
    path inputs
    output:
    path 'report'
    script:
    """
    python -m rnaseq.workflow_support report --contract $contract --inputs $inputs --l1 $l1 --l2 $l2 --out report
    """
}

process TECHNICAL_REPORT_L1 {
    tag 'HTML technical report (L1 only)'
    publishDir params.outdir, mode: 'copy', overwrite: false
    stageInMode 'copy'
    input:
    path l1
    path contract
    path inputs
    output:
    path 'report'
    script:
    """
    python -m rnaseq.workflow_support report --contract $contract --inputs $inputs --l1 $l1 --out report
    """
}

process ENRICHMENT_ANALYSIS {
    tag { module }
    maxForks 1
    // Each task emits exactly one backend-specific directory.  Publishing the
    // shared parent directory would make gsea-go and gsea-kegg collide.
    publishDir "${params.outdir}/l2", mode: 'copy', overwrite: false
    stageInMode 'copy'
    input:
    tuple val(module), path(l2)
    path contract
    path inputs
    output:
    tuple val(module), path('enrichment/*')
    script:
    def rScript = module == 'gsea-go' ? 'gsea_analysis.R' : 'kegg_analysis.R'
    """
    python -m rnaseq.workflow_support enrichment-config --kind ${module} --contract $contract --inputs $inputs --l2 $l2 --out enrichment-config.json
    Rscript ${params.r_scripts}/${rScript} --config enrichment-config.json
    """
}

workflow {
    if( !params.contract || !params.inputs || !params.outdir || !params.analysis_level ) error 'Specify --contract, --inputs, --outdir and --analysis_level'
    if( !(params.analysis_level in ['L1', 'L2']) ) error 'analysis_level must be L1 or L2'
    contract = Channel.value(file(params.contract))
    inputs = Channel.value(file(params.inputs))
    l1 = L1_ANALYSIS(contract, inputs)
    modules = params.enrichment ? params.enrichment.split(',').findAll { it } : []
    allowedEnrichment = ['gsea-go', 'gsea-kegg']
    if( modules.any { !(it in allowedEnrichment) } ) error 'Only internal GSEA backends gsea-go and gsea-kegg are supported'
    if( params.analysis_level == 'L1' ) {
        if( modules ) error 'L1 analysis_level cannot enable enrichment'
        TECHNICAL_REPORT_L1(l1, contract, inputs)
    } else {
        l2 = L2_ANALYSIS(l1, contract, inputs)
        if( modules ) {
            enrichment_dirs = ENRICHMENT_ANALYSIS(Channel.fromList(modules).combine(l2), contract, inputs)
                .map { module, directory -> directory }
                .collect()
            TECHNICAL_REPORT(l1, l2, enrichment_dirs, contract, inputs)
        } else {
            TECHNICAL_REPORT_NO_ENRICHMENT(l1, l2, contract, inputs)
        }
    }
}
