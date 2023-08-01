#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { checkInput_data; checkInput_outdir } from './modules-function-definitions.nf'
data_dir         = checkInput_data(params.data)
outdir           = checkInput_outdir(params.outdir)
qc_dir           = file("${outdir}/1_Read_QC", type: 'dir')
outdir.mkdir()

process run_fastqc {
    label 'fastqc'
    tag { "fastqc: ${sample}" }
    
    input:
    tuple val(sample), path(reads)
    
    output:
    tuple val(sample), path("${sample}*.{html,zip}"), emit: qc_html
    
    """
    fastqc ${reads.findAll().join(' ') } \
        --threads ${task.cpus} \
        --noextract
    """
}

process run_multiqc {
    label 'multiqc'
    tag { 'multiqc: all' }
    publishDir "${qc_dir}", mode: 'copy', overwrite: true

    input:
    path(dir)
    
    output:
    tuple path("multiqc_report.html"), path("multiqc_data"), emit: multiQC
    
    """
    multiqc . --force
    """
}
