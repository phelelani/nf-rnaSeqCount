#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { checkInput_genome; checkInput_genes; checkInput_data; checkInput_outdir } from './modules-function-definitions.nf'
genome           = checkInput_genome(params.genome)
genes            = checkInput_genes(params.genes)
index_dir        = file(genome, type: 'file', checkIfExists: true).getParent()
data_dir         = checkInput_data(params.data)
outdir           = checkInput_outdir(params.outdir)
align_dir        = file("${outdir}/3_Read_Alignment", type: 'dir')        
outdir.mkdir()

process run_star_align {
    label 'star'
    tag { "star: ${sample}" }
    publishDir "${align_dir}/${sample}", mode: 'copy', overwrite: false
    
    input:
    tuple val(sample), path(reads)
    
    output:
    tuple val(sample), path("${sample}_Aligned.out.bam"), emit: star_alignments
    tuple val(sample), path("${sample}_Log.final.out"), emit: star_reports
    
    """
    STAR --runMode alignReads \
        --genomeDir ${index_dir} \
        --readFilesCommand gunzip -c \
        --readFilesIn ${reads.findAll().join(' ')} \
        --runThreadN ${task.cpus} \
        --outSAMtype BAM Unsorted \
        --outSAMunmapped Within \
        --outSAMattributes Standard \
        --outFileNamePrefix ${sample}_ 
    """
}

process run_multiqc_star {
    label 'multiqc'
    tag { 'multiqc: star' }
    publishDir "${align_dir}", mode: 'copy', overwrite: false

    input:
    path(dir) 
    
    output:
    path("qc_star"), emit: multiqc_report
    
    """
    multiqc . --force --outdir qc_star
    """
}
