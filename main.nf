#!/usr/bin/env nextflow
//
//  DO NOT EDIT FROM HERE!! - Unless you brave like King Shaka of course! 
/*  ======================================================================================================
 *  HELP MENU
 *  ======================================================================================================
 */
if (params.help) {
    log.info ''
    log.info "===================================="
    log.info "         nf-rnaSeqCount v0.1        "
    log.info "===================================="
    log.info ''
    log.info 'USAGE: '
    log.info 'nextflow run main.nf --data "/path/to/data" --filetype "type" --out "/path/to/output" --genome "/path/to/genome.fa" --index "/path/to/STARIndex" --genes "/path/to/genes.gtf" --bind "/path/to/bind_1;/path/to/bind_2" -profile "profile" '
    log.info ''
    log.info 'HELP: '
    log.info 'nextflow run main.nf --help'
    log.info ''
    log.info 'MANDATORY ARGUEMENTS:'
    log.info '    --data      FOLDER     Path to where the input data is located (where fastq files are located)'
    log.info '    --filetype  STRING     Extension of the input FASTQ files (fastq | fq | fastq.gz | fq.gz | fastq.bz2 | fq.bz2)'
    log.info '    --out       FOLDER     Path to where the output should be directed (will be created if it does not exist).'
    log.info '    --genome    FILE       The whole genome sequence (fasta | fa | fna)'
    log.info '    --index     FOLDER     Path to where the STAR index files are locaded'
    log.info '    --genes     FILE       The genome annotation file (gtf)'
    log.info '    --bind      FOLDER(S)  Paths to be passed onto the singularity image (Semi-colon separated)'
    log.info '     -profile   SRTING     Executor to be used'
    log.info ''
    log.info "====================================\n"
    exit 1
}
//
//
/*  ======================================================================================================
 *  CHECK ALL USER INPUTS
 *  ======================================================================================================
 */
if (params.data == null) {
    exit 1, "\nPlease enter a directory with input FASTQ/FASTQ.GZ files."
} else{
    data_path = file(params.data, type: 'dir')  // Path to where the input data is located (where fastq files are located).
}

switch ( params.filetype ) {
case ['fastq','fq']:
    ext = params.filetype
    read_file_cmd = ''
    break
case ['fastq.gz','fq.gz']:
    ext = params.filetype
    read_file_cmd = '--readFilesCommand gunzip -c'
    break
case ['fastq.bz2','fq.bz2']:
    ext = params.filetype
    read_file_cmd = '--readFilesCommand bunzip2 -c'
    break
case null:
    ext = 'fastq.gz'
    read_file_cmd = '--readFilesCommand gunzip -c'
    break
}

if(params.out == null) {
    params.out = "${baseDir}/results_nf-rnaSeqCount"
} else{
    out_path = file(params.out, type: 'dir')   // Path to where the output should be directed.
}

if(params.genome == null) {
    exit 1, "Please provide a FASTA sequence of the reference genome."
} else{
    genome = file(params.genome, type: 'file')  // The whole genome sequence.
}

if(params.index == null) {
    exit 1, "Please provide a STAR index."
} else{
    index = file(params.index, type: 'dir')  // Path to where the STAR index files are locaded.
}

if(params.genes == null) {
    exit 1, "Please provide an annotation GTF file."
} else{
    genes = file(params.genes, type: 'file')  // The genome annotation file.
}

if(params.bind == null) {
    bind = "No paths specified for binding to Singularity images! I will assume the data lies somewhere in the '$HOME' directory."
} else{
    bind = params.bind.split(';')  // Paths to be passed onto the singularity image
}
//
//
/*  ======================================================================================================
 *  RUN INFO
 *  ======================================================================================================
 */
log.info "===================================="
log.info "           nf-rnaSeqCount           "
log.info "===================================="
log.info "Input data            : $data_path"
log.info "Input file extension  : $ext"
log.info "STAR readFile command : $read_file_cmd"
log.info "Output path           : $out_path"
log.info "Genome                : $genome"
log.info "Genome Index (STAR)   : $index"
log.info "Genome annotation     : $genes"
log.info "Paths to bind         : $bind"
log.info "====================================\n"
//
//
/*  ======================================================================================================
 *  PIPELINE START
 *  ======================================================================================================
 */
// CREATE OUTPUT DIRECTORY
out_path.mkdir()

// GET INPUT READS
read_pair = Channel.fromFilePairs("${data_path}/*{R,read}[1,2].${ext}", type: 'file') 
.ifEmpty { error "ERROR - Data input: \nOooops... Cannot find any '.fastq' or '.fq' files in ${data_path}. Please specify a folder with '.fastq' or '.fq' files."}


// 1. ALIGN READS TO REFERENCE GENOME
process runSTAR_process {
    cpus 6
    memory '60 GB'
    time '20h'
    scratch '$HOME/tmp'
    tag { sample }
    publishDir "$out_path/${sample}", mode: 'copy', overwrite: false, pattern: "${sample}*.{out,tab}"

    input:
    set sample, file(reads) from read_pair
    
    output:
    set sample, file("${sample}*.{out,tab}") into star_results
    set sample, file("${sample}_Aligned.sortedByCoord.out.bam") into bams_htseqCounts, bams_featureCounts
    
    """
    STAR --runMode alignReads \
        --genomeDir ${index} ${read_file_cmd} \
        --readFilesIn ${reads.get(0)} ${reads.get(1)} \
        --runThreadN 5 \
        --outSAMtype BAM SortedByCoordinate \
        --outFileNamePrefix ${sample}_
    """
}

// 2. GET RAW COUNTS USING HTSEQ-COUNT
// 2a - Use HTSeqCount to get the raw gene counts
process runHTSeqCount_process {
    cpus 6
    memory '60 GB'
    time '20h'
    scratch '$HOME/tmp'
    tag { sample }
    publishDir "$out_path/htseqCounts", mode: 'copy', overwrite: false

    input:
    set sample, file(bam) from bams_htseqCounts

    output:
    set sample, "${sample}.txt" into htseqCounts
    set sample, "${sample}.txt" into htseqCounts_cleanup
    
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

// 2b - Collect files for HTSeq cleanup
htseqCounts_cleanup
.collectFile () { item -> [ 'cleanup_htseqcounts.txt', "${item.get(1)}" + '\n' ] }
.set { list_htseqcounts }

// 2c - Cleanup HTSeqCounts
process runCleanHTSeqCounts_process {
    cpus 1
    memory '10 GB'
    time '5h'
    scratch '$HOME/tmp'
    tag { 'HTSEqCounts Cleanup' }
    publishDir "$out_path/htseqCounts", mode: 'copy', overwrite: false
    
    input:
    file(file_list) from list_htseqcounts

    output:
    file(out_file) into htseqCounts_cleaned

    shell:
    in_file = "${file_list}"
    out_file = "gene_counts_final.txt"
    template "clean_htseqCounts.sh"
}


// 3. GET RAW COUNTS USING FEATURECOUNTS
// 3a - Get all the bam file locations to process with featureCounts
bams_featureCounts
.collectFile () { item -> [ 'sample_bams.txt', "${item.get(1)}" + ' ' ] }
.set { sample_bams }

// 3b - Use featureCounts to get raw gene counts
process runFeatureCounts_process {
    cpus 6
    memory '60 GB'
    time '20h'
    scratch '$HOME/tmp'
    tag { 'featureCounts - ALL' }
    publishDir "$out_path/featureCounts", mode: 'copy', overwrite: false

    input:
    file(samples) from sample_bams

    output:
    file('gene_counts*') into featureCounts
    file('gene_counts.txt') into featureCounts_raw
    
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

// 3c - Cleanup featureCounts gene counts
process runCleanFeatureCounts_process {
    cpus 1
    memory '10 GB'
    time '5h'
    scratch '$HOME/tmp'
    tag { 'featureCounts Cleanup' }
    publishDir "$out_path/featureCounts", mode: 'copy', overwrite: false
    
    input:
    file(raw_counts) from featureCounts_raw

    output:
    file(out_file) into featureCounts_cleaned

    shell:
    in_file = "${raw_counts}"
    out_file = "gene_counts_final.txt"
    template "clean_featureCounts.sh"
}

// 4. DO QC FOR ALL THE OUTPUT FILES 
// 4a - Collect files for STAR QC
star_results.collectFile () { item -> [ 'qc_star.txt', "${item.get(1).find { it =~ 'Log.final.out' } }" + ' ' ] }
.set { qc_star }

// 4b - Collect files for HTSeq QC
htseqCounts
.collectFile () { item -> [ 'qc_htseqcounts.txt', "${item.get(1)}" + ' ' ] }
.set { qc_htseqcounts }

// 4c - Collect files for featureCounts QC
featureCounts
.collectFile () { item -> [ 'qc_featurecounts.txt', "${item.find { it =~ 'txt.summary' } }" + ' ' ] }
.set { qc_featurecounts }

// 4d Get QC for STAR, HTSeqCounts and featureCounts
process runMultiQC_process {
    cpus 1
    memory '10 GB'
    time '5h'
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
//
//
/*  ======================================================================================================
 *  WORKFLOW SUMMARY
 *  ======================================================================================================
 */
workflow.onComplete {
    println "\n\n===================================="
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
