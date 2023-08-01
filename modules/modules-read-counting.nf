#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { checkInput_genome; checkInput_genes; checkInput_data; checkInput_outdir } from './modules-function-definitions.nf'
genome           = checkInput_genome(params.genome)
genes            = checkInput_genes(params.genes)
index_dir        = file(genome, type: 'file', checkIfExists: true).getParent()
data_dir         = checkInput_data(params.data)
outdir           = checkInput_outdir(params.outdir)
counts_dir       = file("${outdir}/4_Read_Counts", type: 'dir')
outdir.mkdir()

// USE HTSEQCOUNTS TO GET RAW READ COUNTS
process run_htseqcount {
    label 'htseqcount'
    tag { "htseqcount: ${sample}" }
    publishDir "${counts_dir}/htseqCounts", mode: 'copy', overwrite: true
    
    input:
    tuple val(sample), path(bam)
    
    output:
    tuple val(sample), path("${sample}.txt"), emit: htseqcounts_raw
    
    """
    htseq-count --format=bam \
        --order=name \
        --idattr=gene_id \
        --stranded=reverse \
        --mode=union \
        --type=exon \
        ${bam} ${genes} > ${sample}.txt
    """
}

// // COLLECT FILES FOR HTSEQ CLEANUP
// CLEANUP HTSEQCOUNTS
process run_clean_htseqcounts  {
    tag { 'htseqcounts: cleanup' }
    publishDir "${counts_dir}/htseqCounts", mode: 'copy', overwrite: false
    
    input:
    path(htseqcounts_list)
    
    output:
    path("gene_counts_final.txt"), emit: htseqcounts_cleaned

    """
    clean_htseqcounts.sh ${htseqcounts_list} gene_counts_final.txt
    """    
}

// USE FEATURECOUNTS TO GET RAW GENE COUNTS
process run_featurecounts {
    label 'featurecounts'
    tag { 'featurecounts: all' }
    publishDir "${counts_dir}/featureCounts", mode: 'copy', overwrite: false
    
    input:
    path(bam_list)
    
    output:
    path('gene_counts.txt'), emit: featurecounts_raw
    path("gene_counts.txt.summary"), emit: featurecounts_reports
    
    """
    featureCounts -p -B -C -P -J -s 2 \
        -G ${genome} -J \
        -t exon \
        -d 40 \
        -g gene_id \
        -a ${genes} \
        -T ${task.cpus} \
        -o gene_counts.txt \
        `< ${bam_list}`
    """
}

// CLEANUP FEATURECOUNTS GENE COUNTS
process run_clean_featurecounts {
    tag { 'featurecounts: cleanup' }
    publishDir "${counts_dir}/featureCounts", mode: 'copy', overwrite: false
    
    input:
    path(raw_counts)
    
    output:
    path("gene_counts_final.txt"), emit: featurecounts_cleaned
    
    """
    clean_featurecounts.sh ${raw_counts} gene_counts_final.txt
    """
}

process run_multiqc_counts {
    label 'multiqc'
    tag { 'multiqc: counts' }
    publishDir "${counts_dir}", mode: 'copy', overwrite: false

    input:
    path(dir) 
    
    output:
    path("qc_read-counts"), emit: multiqc_report
    
    """
    multiqc . --force --outdir qc_read-counts
    """
}
