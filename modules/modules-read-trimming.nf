#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { checkInput_data; checkInput_outdir; checkInput_stranded; checkInput_trim } from './modules-function-definitions.nf'
data_dir         = checkInput_data(params.data)
outdir           = checkInput_outdir(params.outdir)
trim_dir         = file("${outdir}/2_Read_Trimming", type: 'dir')
stranded         = checkInput_stranded(params.singleEnd,params.pairedEnd)
trim_params      = checkInput_trim(params.trim,stranded)
outdir.mkdir()

process run_trimmomatic {
    label 'trimmomatic'
    tag { "trimmomatic: ${sample}" }
    publishDir "${trim_dir}/${sample}", mode: 'copy', overwrite: false
    
    input:
    tuple val(sample), path(reads)
    
    output:
    tuple val(sample), path("${sample}_trimmed{_1P.fastq.gz,_2P.fastq.gz,.fastq.gz}"), emit: read_pairs_trimmed
    tuple val(sample), path("${sample}_trim_out.log"), emit: trimmomatic_reports
    
    """
    ln -s /opt/Trimmomatic-0.39/adapters/*.fa .

    if [[ ${stranded} == "paired-end" ]]
    then
        java -jar /opt/Trimmomatic-0.39/trimmomatic-0.39.jar PE \
            -threads ${task.cpus} \
            ${reads.findAll().join(' ')} \
            -baseout ${sample}_trimmed.fastq.gz \
            ${trim_params} 2> ${sample}_trim_out.log
    elif [[ ${stranded} == "single-end" ]]
    then
        java -jar /opt/Trimmomatic-0.39/trimmomatic-0.39.jar SE \
            -threads ${task.cpus} \
            ${reads.findAll().join(' ')} \
            ${sample}_trimmed.fastq.gz \
            ${trim_params} 2> ${sample}_trim_out.log
    else :
    fi
    """
}

process run_fastqc_trimming {
    label 'fastqc'
    tag { "fastqc: ${sample}" }
    
    input:
    tuple val(sample), path(reads)
    
    output:
    tuple val(sample), path("${sample}*.{html,zip}"), emit: fastqc_reports
    
    """
    fastqc ${reads.findAll().join(' ') } \
        --threads ${task.cpus} \
        --noextract
    """
}

process run_multiqc_trimming {
    label 'multiqc'
    tag { 'multiqc: trimmomatic' }
    publishDir "${trim_dir}", mode: 'copy', overwrite: false

    input:
    path(dir)
    
    output:
    path("qc_trimmomatic"), emit: multiqc_report
    
    """
    multiqc . --force --outdir qc_trimmomatic
    """
}
