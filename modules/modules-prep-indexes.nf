#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { checkInput_genome; checkInput_genes } from './modules-function-definitions.nf'
genome           = checkInput_genome(params.genome)
genes            = checkInput_genes(params.genes)
index_dir        = file(genome, type: 'file', checkIfExists: true).getParent()

process run_star_index {
    label 'star'
    tag { "star: index" }
    publishDir "${index_dir}", mode: 'copy', overwrite: true
    
    output:
    tuple val("starIndex"), path("*"), emit: star_index
    
    """
    STAR --runThreadN ${task.cpus} \
        --runMode genomeGenerate \
        --genomeDir . \
        --genomeFastaFiles ${genome} \
        --sjdbGTFfile ${genes} \
        --sjdbOverhang 99
     """
}

process run_bowtie_index {
    label 'bowtie'
    tag { "bowtie: index" }
    publishDir "${index_dir}", mode: 'copy', overwrite: true  

    output:
    tuple val("bowtieIndex"), path("*"), emit: bowtie_index
    """
    bowtie2-build --threads ${task.cpus} ${genome} genome
    """
}
