#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/*  ======================================================================================================
 *  HELP MENU
 *  ======================================================================================================
 */
line="=".multiply(100)
ver="nf-rnaSeqCount v0.2"
if (params.help) {
    println "clear".execute().text    
    println "\n${line}"
    println "#".multiply(48 - ("${ver}".size() / 2 )) + "  ${ver}   " + "#".multiply(48 - ("${ver}".size() / 2 ))
    println "${line}\n"
    println "USAGE:"
    println "nextflow run nf-rnaSeqCount -profile \"slurm\" --data \"/path/to/data\" --genome \"/path/to/genome.fa\" --genes \"/path/to/genes.gtf\"\n" 
    println "HELP:"
    println "nextflow run nf-rnaSeqCount --help\n"
    println "MANDATORY ARGUEMENTS:"
    println "-profile     STRING    Executor to be used. Available options:"
    println "\t\t\t\t\"standard\"           : Local execution (no job scheduler)."
    println "\t\t\t\t\"slurm\"              : SLURM scheduler."
    println "--workflow   STRING    To specify which step of the workflow you are running (see https://github.com/phelelani/nf-rnaSeqCount)."
    println "                       Availeble options:"
    println "\t\t\t\t\"genome-indexing\"    : For indexing your reference genome using STAR and Bowtie2."
    println "\t\t\t\t\"read-qc\"            : For performing general QC on your reads using FastQC. "
    println "\t\t\t\t\"read-trimming\"      : For trimming low quality bases and removing adapters from your reads using Trimmmomatic."
    println "\t\t\t\t\"read-alignment\"     : For aligning your reads to your reference genome using STAR."
    println "\t\t\t\t\"read-counting\"      : For counting features in your reads using HTSeq-count and featureCounts."
    println "--data       FOLDER    Path to where the input data (FASTQ files) is located. Supported FASTQ files:"
    println "\t\t\t\t[ fastq | fastq.gz | fastq.bz2 | fq | fq.gz | fq.bz2 ]"
    println "--genome     FILE      The whole genome FASTA sequence. Supported FASTA files:"
    println "\t\t\t\t[ fasta | fa | fna ]"
    println "--genes      FILE      The genome annotation GFT file. Supported GTF file:"
    println "\t\t\t\t[ gtf ]"
    println "OPTIONAL ARGUEMENTS:"
    println "--help                 To show this menu."
    println "--outdir     FOLDER    Path to where the output should be directed."
    println "                       Default: \$PWD/results_nf-rnaSeqCount)."
    println "--pairedEnd            If working with paired-end FASTQ files (default)."
    println "--singleEnd            If working with single-end FASTQ files."
    println "--trim       STRING    Parameters for Trimmomatic. See http://www.usadellab.org/cms/index.php?page=trimmomatic for a more detailed use."
    println "                       The default parameters for Trimmomatic I have given you here (for both paird- and single-end sequences) are:"
    println "\t\t\t\tFor paired-end: \"ILLUMINACLIP:TruSeq3-PE-2.fa:2:30:10:8:true TRAILING:28 MINLEN:40\""
    println "\t\t\t\tFor single-end: \"ILLUMINACLIP:TruSeq3-SE.fa:2:30:10:8:true TRAILING:28 MINLEN:40\""
    // println "--max_memory STRING    Maximum memory you have access to."
    // println "                       Default: \"200.GB\""
    // println "--max_cpus   STRING    Maximum CPUs you have access to."
    // println "                       Default: \"24\""
    // println "--max_time   STRING    Maximum time you have access to."
    // println "                       Default: \"24.h\""
    println "${line}\n"
    exit 1
}

params.workflow  = null
params.genome    = null
params.genes     = null
params.data      = null
params.outdir    = null
params.singleEnd = null
params.pairedEnd = null
params.trim      = null

include { workflow_error; checkInput_data; main_data_error; trim_data_error;
         bams_error; checkInput_genome; checkInput_genes; checkInput_outdir;
         checkInput_stranded; checkInput_trim } from './modules/modules-function-definitions.nf'

workflow         = params.workflow
stranded         = checkInput_stranded(params.singleEnd,params.pairedEnd)
trim_params      = checkInput_trim(params.trim,stranded)
ext              = "fastq,fastq.gz,fastq.bz2,fq,fq.gz,fq.bz2"

// USER INPUT WORKFLOW: WHICH ANALYSIS TO RUN!
switch (workflow) {
    case [null]:
        exit 1, workflow_error(workflow)
    case ["genome-indexing"]:
        include { run_star_index; run_bowtie_index } from './modules/modules-prep-indexes.nf'
        genome = checkInput_genome(params.genome)
        genes = checkInput_genes(params.genes)
        index_dir = file(genome, type: 'file', checkIfExists: true).getParent()
        break
    case ["read-qc"]:
        include { run_fastqc; run_multiqc } from './modules/modules-read-qc.nf'
        data_dir = checkInput_data(params.data)
        outdir = checkInput_outdir(params.outdir)
        switch (stranded) {
            case ["paired-end"]:
                samples = Channel.fromFilePairs("${data_dir}/**{R,read,_}[1,2]*.{${ext}}", type: 'file')
                    .ifEmpty { exit 1, main_data_error(data_dir) }
                break
            case["single-end"]:
                samples = Channel.fromFilePairs("${data_dir}/**.{${ext}}", type: 'file', size:1)
                    .ifEmpty { exit 1, main_data_error(data_dir) }
                break
        }
        break
    case ["read-trimming"]:
        include { run_trimmomatic; run_fastqc_trimming; run_multiqc_trimming } from './modules/modules-read-trimming.nf'
        data_dir = checkInput_data(params.data)
        outdir = checkInput_outdir(params.outdir)
        switch (stranded) {
            case ["paired-end"]:
                samples = Channel.fromFilePairs("${data_dir}/**{R,read,_}[1,2]*.{${ext}}", type: 'file')
                    .ifEmpty { exit 1, main_data_error(data_dir) }
                break
            case["single-end"]:
                samples = Channel.fromFilePairs("${data_dir}/**.{${ext}}", type: 'file', size:1)
                    .ifEmpty { exit 1, main_data_error(data_dir) }
                break
        }
        break
    case ["read-alignment"]:
        include { run_star_align; run_multiqc_star } from './modules/modules-read-alignment.nf'
        genome = checkInput_genome(params.genome)
        genes = checkInput_genes(params.genes)
        index_dir = file(genome, type: 'file', checkIfExists: true).getParent()
        data_dir = checkInput_data(params.data)
        outdir = checkInput_outdir(params.outdir)        
        switch (stranded) {
            case ["paired-end"]:
                samples = Channel.fromFilePairs("${data_dir}/**{R,read,_}[1,2]*.{${ext}}", type: 'file')
                    .ifEmpty { exit 1, main_data_error(data_dir) }
                break
            case["single-end"]:
                samples = Channel.fromFilePairs("${data_dir}/**.{${ext}}", type: 'file', size:1)
                    .ifEmpty { exit 1, main_data_error(data_dir) }
                break
        }
        break
    case ["read-counting"]:
        include { run_htseqcount; run_clean_htseqcounts; run_featurecounts;
                 run_clean_featurecounts; run_multiqc_counts } from './modules/modules-read-counting.nf'
        genome = checkInput_genome(params.genome)
        genes = checkInput_genes(params.genes)
        index_dir = file(genome, type: 'file', checkIfExists: true).getParent()
        data_dir = checkInput_data(params.data)
        outdir = checkInput_outdir(params.outdir)
        align_dir = file("${outdir}/3_Read_Alignment", type: 'dir')        
        bams = Channel.fromFilePairs("${align_dir}/**_Aligned.out.bam", size:-1) { 
            it -> "${it.baseName.replace(/_Aligned.out/, "")}" }
            .ifEmpty { exit 1, bams_error(align_dir) }
        break
    default:
        exit 1, workflow_error(workflow)
}

//  ======================================================================================================
//  RUN INFO
//  ======================================================================================================
// options="nf-rnaSeqCount v0.2 - Input | Output | Parameters:"
// println "\n" + "=".multiply(100)
// println "#".multiply(48 - ("${options}".size() / 2 )) + "  ${options}  " + "#".multiply(48 - ("${options}".size() / 2 ))
// println "=".multiply(100)
// println "Input data              : ${data_dir}"
// println "Input data type         : ${stranded}"
// println "Output directories      : ${outdir}"
// println ' '.multiply(26) + "- ${qc_dir.baseName}"
// println ' '.multiply(26) + "- ${trim_dir.baseName}"
// println ' '.multiply(26) + "- ${align_dir.baseName}"
// println ' '.multiply(26) + "- ${counts_dir.baseName}"
// println "Genome                  : ${genome}"
// println "Genome annotation       : ${genes}"
// println "Trimmomatic parameters  : ${trim_params}"
// println "=".multiply(100)
// println " "

// ======================================================================================================
//  PIPELINE START
// ======================================================================================================

// PREPARE GENOME INDEXES
workflow GENOME_INDEXING {
    main:
    run_star_index()
    run_bowtie_index()
}

// 
workflow READ_QC {
    take:
    samples
    
    main:
    run_fastqc(samples)
    run_fastqc.out.qc_html
        .map { it -> it[1] }
        .collect()
        .set { fastqc_out }
    run_multiqc(fastqc_out)
}

// 
workflow READ_TRIMMING {
    take:
    sample

    main:
    run_trimmomatic(sample)
    run_fastqc_trimming(run_trimmomatic.out.read_pairs_trimmed)
    run_trimmomatic.out.trimmomatic_reports
        .join(run_fastqc_trimming.out.fastqc_reports)
        .map { it -> it[1,-1] }
        .flatten()
        .collect()
        .set{ reports }
    run_multiqc_trimming(reports)
}

// 
workflow READ_ALIGNMENT {
    take:
    samples

    main:
    run_star_align(samples)
    run_star_align.out.star_reports
        .map { it -> it[1]}
        .flatten()
        .collect()
        .set { reports }
    run_multiqc_star(reports)
}

//
workflow READ_COUNTING {
    take:
    bams

    main:
    run_htseqcount(bams)
    run_htseqcount.out.htseqcounts_raw
        .collectFile() { item -> [ 'cleanup_htseqcounts.txt', "${item.get(1)}" + '\n' ] }
        .set { htseqcounts_counts }
    run_clean_htseqcounts(htseqcounts_counts)
    bams
        .collectFile() { item -> [ 'sample_bams.txt', "${item.get(1).join().toString()}" + ' ' ] }
        .set { featurecounts_bams }
    run_featurecounts(featurecounts_bams)
    run_clean_featurecounts(run_featurecounts.out.featurecounts_raw)
    run_htseqcount.out.htseqcounts_raw
        .map { it -> it[1] }
        .collect()
        .concat(run_featurecounts.out.featurecounts_reports)
        .flatten()
        .toList()
        .set { reports } 
    run_multiqc_counts(reports)
}

workflow {
    switch (workflow) {
        case [null]:
            exit 1, workflow_error(workflow)
        case ['genome-indexing']:
            GENOME_INDEXING()
            break
        case['read-qc']:
            READ_QC(samples)
            break
        case['read-trimming']:
            READ_TRIMMING(samples)
            break
        case['read-alignment']:
            READ_ALIGNMENT(samples)
            break
        case['read-counting']:
            READ_COUNTING(bams)
            break
        default:
            exit 1, workflow_error(workflow)
            break
    }
}
