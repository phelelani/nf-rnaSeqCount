#!/usr/bin/env nextflow

// PIPELINE PARAMETERS - Edit if brave... Else, specify options on command line
params.data    = "/spaces/phelelani/ssc_data/data_trimmed/inflated"                                                   // Path to where the input data is located (where fastq files are located).
params.out     = "/spaces/phelelani/ssc_data/nf-rnaSeqCount"                                                          // Path to where the output should be directed.
params.genome  = "/global/blast/reference_genomes/Homo_sapiens/Ensembl/GRCh38/Sequence/WholeGenomeFasta/genome.fa"    // The whole genome sequence
params.index   = "/global/blast/reference_genomes/Homo_sapiens/Ensembl/GRCh38/Sequence/STARIndex"                     // Path to where the STAR index files are locaded
params.genes   = "/global/blast/reference_genomes/Homo_sapiens/Ensembl/GRCh38/Annotation/Genes/genes.gtf"             // The genome annotation file
params.bind    = '/global/;/spaces/'                                                                                  // Paths to be passed onto the singularity image

// DO NOT EDIT FROM HERE
data_path      = file(params.data, type: 'dir')                                                                       // Path to where the input data is located (where fastq files are located). 
out_path       = file(params.out, type: 'dir')                                                                        // Path to where the output should be directed.
genome         = file(params.genome, type: 'file')                                                                    // The whole genome sequence
index          = file(params.index, type: 'dir')                                                                      // Path to where the STAR index files are locaded 
genes          = file(params.genes, type: 'file')                                                                     // The genome annotation file 
bind           = params.bind.split(';')                                                                               // Paths to be passed onto the singularity image
//======================================================================================================
//
//
//
//======================================================================================================
// HELP MENU
if (params.help) {
    log.info ''
    log.info "===================================="
    log.info "         nf-rnaSeqCount v0.1        "
    log.info "===================================="
    log.info ''
    log.info 'USAGE: '
    log.info 'nextflow run main.nf --data /path/to/data --out /path/to/output --genome /path/to/genome.fa --genes /path/to/genes.gtf --index /path/to/STARIndex'
    log.info ''
    log.info 'HELP: '
    log.info 'nextflow run main.nf --help'
    log.info ''
    log.info 'MANDATORY ARGUEMENTS:'
    log.info '    --data     FOLDER    Path to where the input data is located (fastq | fq)'
    log.info '    --out      FOLDER    Path to where the output should be directed (will be created if it does not exist).'
    log.info '    --genome   FILE      The whole genome sequence (fasta | fa | fna)'
    log.info '    --index    FOLDER    Path to where the STAR index files are locaded'
    log.info '    --genes    FILE      The genome annotation file (gtf)'
    log.info '    --bind     FOLDER(S) Paths to be passed onto the singularity image'
    log.info ''
    log.info "====================================\n"
    exit 1
}

// RUN INFO
log.info "===================================="
log.info "           nf-rnaSeqCount           "
log.info "===================================="
log.info "Input data          : ${data_path}"
log.info "Output path         : ${out_path}"
log.info "Genome              : ${genome}"
log.info "Genome Index (STAR) : ${index}"
log.info "Genome annotation   : ${genes}"
log.info "Paths to bind       : ${bind}"
log.info "====================================\n"
//======================================================================================================
//
//
//
//======================================================================================================
// PIPELINE START
//Create output directory
out_path.mkdir()

// Get input reads
read_pair = Channel.fromFilePairs("${data_path}/*R[1,2].fastq", type: 'file') 
                   .ifEmpty { error "ERROR - Data input: \nOooops... Cannot find any '.fastq' or '.fq' files in ${data_path}. Please specify a folder with '.fastq' or '.fq' files."}


// 1. Align reads to reference genome
process runSTAR_process {
    cpus 6
    memory '40 GB'
    time '10h'
    scratch '$HOME/tmp'
    tag { sample }
    publishDir "$out_path/${sample}", mode: 'copy', overwrite: false

    input:
    set sample, file(reads) from read_pair
    
    output:
    set sample, "${sample}*" into star_results
    set sample, file("${sample}_Aligned.sortedByCoord.out.bam") into bams_htseqCounts, bams_featureCounts
    
    """
    STAR --runMode alignReads \
        --genomeDir ${index} \
        --readFilesIn ${reads.get(0)} ${reads.get(1)} \
        --runThreadN 5 \
        --outSAMtype BAM SortedByCoordinate \
        --outFileNamePrefix ${sample}_
    """
}

// 2. Get raw counts using HTSeq-count
process runHTSeqCount_process {
    cpus 4
    memory '5 GB'
    time '10h'
    scratch '$HOME/tmp'
    tag { sample }
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

// 3. Get raw counts using featureCounts
process runFeatureCounts_process {
    cpus 6
    memory '5 GB'
    time '10h'
    scratch '$HOME/tmp'
    tag { 'featureCounts - ALL' }
    publishDir "$out_path/featureCounts", mode: 'copy', overwrite: false

    input:
    file(samples) from sample_bams

    output:
    file('gene_counts*') into featureCounts
    
    """
    featureCounts -p -B -C -P -J -s 2 \
        -G $genome -J \
        -t exon \
        -d 40 \
        -g gene_id \
        -a $genes \
        -T 5 \
        -o gene_counts.txt \
        `< ${samples}`
    """
}

// 4a. Collect files for STAR QC
star_results.collectFile () { item -> [ 'qc_star.txt', "${item.get(1).find { it =~ 'Log.final.out' } }" + ' ' ] }
.set { qc_star }

// 4b. Collect files for HTSeq QC
htseqCounts
.collectFile () { item -> [ 'qc_htseqcounts.txt', "${item.get(1)}" + ' ' ] }
.set { qc_htseqcounts }

// 4c. Collect files for featureCounts QC
featureCounts
.collectFile () { item -> [ 'qc_featurecounts.txt', "${item.find { it =~ 'txt.summary' } }" + ' ' ] }
.set { qc_featurecounts }

// 4. Get QC for STAR, HTSeqCounts and featureCounts
process runMultiQC_process {
    cpus 1
    memory '5 GB'
    time '10h'
    scratch '$HOME/tmp'
    tag { 'MultiQC - ALL' }
    publishDir "$out_path/report_QC", mode: 'copy', overwrite: false

    input:
    file(star) from qc_star
    file(htseqcounts) from qc_htseqcounts
    file(featurecounts) from qc_featurecounts

    output:
    file('*') into multiQC
    
    """
    multiqc `< ${star}` `< ${htseqcounts}` `< ${featurecounts}` --force
    """
}
//======================================================================================================
//
//
//
//======================================================================================================
// WORKFLOW SUMMARY
workflow.onComplete {
    println "===================================="
    println "Pipeline execution summary:"
    println "===================================="
    println "Execution command   : ${workflow.commandLine}"
    println "Execution name      : ${workflow.runName}"
    println "Workflow start      : ${workflow.start}"
    println "Workflow end        : ${workflow.complete}"
    println "Workflow duration   : ${workflow.duration}"
    println "Workflow completed? : ${workflow.success}"
    println "Work directory      : ${workflow.workDir}"
    println "Project directory   : ${workflow.projectDir}"
    println "Execution directory : ${workflow.launchDir}"
    println "Configuration files : ${workflow.configFiles}"
    println "Workflow containers : ${workflow.container}"
    println "exit status : ${workflow.exitStatus}"
    println "Error report: ${workflow.errorReport ?: '-'}"
    println "===================================="
}

workflow.onError {
    println "Oohhh DANG IT!!... Pipeline execution stopped with the following message: ${workflow.errorMessage}"
}
//======================================================================================================
