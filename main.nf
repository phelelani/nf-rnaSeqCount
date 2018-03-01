#!/usr/bin/env nextflow

params.data = "/spaces/phelelani/ssc_data/data_trimmed/inflated" 
params.out = "/spaces/phelelani/ssc_data/assembly/star_results"
params.genes = "/global/blast/reference_genomes/Homo_sapiens/Ensembl/GRCh38/Annotation/Genes/genes.gtf"
params.refSeq = "/global/blast/reference_genomes/Homo_sapiens/Ensembl/GRCh38/Sequence/WholeGenomeFasta/genome.fa"
params.genome = "/global/blast/reference_genomes/Homo_sapiens/Ensembl/GRCh38/Sequence/STARIndex"

data_path = params.data
out_path = file(params.out)
genes = params.genes
genome = params.genome
refSeq = params.refSeq

out_path.mkdir()

read_pair = Channel.fromFilePairs("${data_path}/*R[1,2].fastq", type: 'file')

// 1. Align reads to reference genome
process runSTAR_process {
    cache = true
    executor 'pbs'
    queue 'batch'
    cpus 9
    memory '50 GB'
    time '100h'
    scratch '$HOME/tmp'
    tag { sample }
    stageInMode 'symlink'
    stageOutMode 'rsync'
    publishDir "$out_path/${sample}", mode: 'copy', overwrite: false

    input:
    set sample, file(reads) from read_pair

    output:
    set sample, "${sample}*" into star_results
    set sample, file("${sample}_Aligned.sortedByCoord.out.bam") into bams_htseqCounts, bams_featureCounts
    
    """
    ~/applications/STAR-2.5.3a/source/STAR --runMode alignReads \
        --genomeDir $genome \
        --readFilesIn ${reads.get(0)} ${reads.get(1)} \
        --runThreadN 18 \
        --outSAMtype BAM SortedByCoordinate \
        --outFileNamePrefix ${sample}_
    """
}

// 2. Get raw counts using HTSeq-count
process runHTSeqCount_process {
    cache = true
    executor 'pbs'
    queue 'batch'
    cpus 3
    memory '5 GB'
    time '100h'
    scratch '$HOME/tmp'
    tag { sample }
    stageInMode 'symlink'
    stageOutMode 'rsync'
    publishDir "$out_path/htseqCounts", mode: 'copy', overwrite: false

    input:
    set sample, file(bam) from bams_htseqCounts

    output:
    set sample, "${sample}.txt" into htseqCounts
    
    """
    htseq-count -f bam \
        -r pos \
        -i gene_id \
        -a 10 \
        -s reverse \
        -m union \
        -t exon \
        $bam $genes > ${sample}.txt
    """
}

// 3a. Get all the bam file locations to process with featureCounts
bams_featureCounts
.collectFile () { item -> [ 'sample_bams.txt', "${item.get(1)}" + ' ' ] }
.set { sample_bams }

// 3b. Get raw counts using featureCounts
process runFeatureCounts_process {
    cache = true
    executor 'pbs'
    queue 'batch'
    cpus 9
    memory '50 GB'
    time '100h'
    scratch '$HOME/tmp'
    tag { sample }
    stageInMode 'symlink'
    stageOutMode 'rsync'
    publishDir "$out_path/featureCounts", mode: 'copy', overwrite: false

    input:
    file(samples) from sample_bams

    output:
    file('gene_counts*') into featureCounts
    
    """
    featureCounts -p -B -C -P -J -s 2 \
        -G $refSeq -J \
        -t exon \
        -d 40 \
        -g gene_id \
        -a $genes \
        -T 18 \
        -o gene_counts.txt \
        `< ${samples}`
    """
}

workflow.onComplete {
    """
    Pipeline execution summary
    ---------------------------
    Completed at: ${workflow.complete}
    Duration    : ${workflow.duration}
    Success     : ${workflow.success}
    workDir     : ${workflow.workDir}
    exit status : ${workflow.exitStatus}
    Error report: ${workflow.errorReport ?: '-'}
    """
}

workflow.onError {
    println "Oops... Pipeline execution stopped with the following message: ${workflow.errorMessage}"
}
