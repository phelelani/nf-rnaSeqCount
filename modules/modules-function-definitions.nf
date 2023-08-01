#!/usr/bin/env nextflow
nextflow.enable.dsl=2

line="=".multiply(100)
/*  ======================================================================================================
 *  CHECK ALL USER INPUTS
 *  ======================================================================================================
 */


// MAIN USER INPUT ERRORS
def workflow_error(input_workflow) {
    """
${line}
I do not recognise the \'--workflow ${input_workflow}\' option you have given me, or you have not given me any \'--workflow\' option at all!
The allowed options for \'--workflow\' are:
\tgenome-indexing   : For indexing your reference genome using STAR and Bowtie2.
\tread-qc           : For performing general QC on your reads using FastQC. 
\tread-trimming     : For trimming low quality bases and removing adapters from your reads using Trimmmomatic.
\tread-alignment    : For aligning your reads to your reference genome using STAR.
\tread-counting     : For counting features in your reads using HTSeq-count and featureCounts.
\tcomplete-workflow : For getting a summary of QC through the analysis using MultiQC.
\nPlease use one of the above options with \'--workflow\' to run the nf-rnaSeqCount workflow!
${line}
"""
}

// USER PARAMETER INPUT: DATA DIRECTORY
def checkInput_data(input) {
    switch (input) {
        case [null]:
            exit 1, """
${line}
I do not recognise the \'--data ${input}\' option you have given me, or you have not given me any \'--data\' option at all!
Please provide a valid directory with you input FASTQ reads with the \'--data\' option to run the nf-rnaSeqCount workflow! 
${line}
"""
            break
        default:
            data_dir = file(input, type: 'dir', checkIfExists: true)
            break
    }
    return data_dir
}


// DATA ABSENT ERRORS
def main_data_error(input_dir) {
    """
${line}
There are no FASTQ file in the directory:
\t${input_dir}
Please ensure that you have given me the correct directory for you FASTQ input reads using the \'--data\' option!
${line}
"""
}

def trim_data_error(input_dir) {
    """
${line}
There are no FASTQ file in the directory:
\t${input_dir}
Are you sure you ran the READ TRIMMING STEP using \'--workflow run.ReadTrimming\' ??
Please ensure that you have ran the READ TRIMMING STEP successfully and try again!
${line}
"""
}

def bams_error(input) {
    """
There are no BAM files in the directory:
\t${input_dir}
Are you sure you ran the READ ALIGNMENT STEP using \'--workflow run.ReadAlignment\' ??
Please ensure that you have ran the READ ALIGNMENT STEP successfully and try again!
"""
}

// USER PARAMETER INPUT: GENOME FASTA FILE
def checkInput_genome(input) {
    switch (input) {
        case [null]:
            exit 1, """
${line}
I do not recognise the \'--genome ${input}\' option you have given me, or you have not given me any \'--genome\' option at all!
Please provide a valid FASTA file (.fasta or .fa) for your reference genome with the \'--genome\' option to run the nf-rnaSeqCount workflow! 
${line}
"""
            break
        default:
            genome = file(input, type: 'file', checkIfExists: true)
            break
    }
    return genome
}

// USER PARAMETER INPUT: GENOME ANNOTATION FILE (GFT/GFF)
def checkInput_genes(input) {
    switch (input) {
        case [null]:
            exit 1, """
${line}
I do not recognise the \'--genes ${input}\' option you have given me, or you have not given me any \'--genes\' option at all!
Please provide a valid GTF annotation file (.gtf) for your reference genome with the \'--genes\' option to run the nf-rnaSeqCount workflow! 
${line}
"""
            break
        default:
            genes = file(input, type: 'file', checkIfExists: true)
            break
    }
    return genes
}

// USER PARAMETER INPUT: OUTPUT DIRECTORY
def checkInput_outdir(input_dir) {
    switch (input_dir) {
        case [null]:
            outdir = file("${launchDir}/results_nf-rnaSeqCount", type: 'dir')
            break
        default:
            outdir = file(params.outdir, type: 'dir')
            break
    }
    return outdir
}

// USER STRANDED MODE: ARE WE DOING PAIRED- OR SINGLE-END?
def checkInput_stranded(input_single,input_paired){
    if (input_single == null && input_paired == null) {
        stranded = "paired-end"
    } else if(input_single) {
        stranded = "single-end"
    } else if(input_paired){
        stranded = "paired-end"
    } else {}
    return stranded
}

// USER PARAMETER INPUT: TRIMMOMATIC OPTIONS
def checkInput_trim(input_trim,input_stranded) {
    switch (input_trim) {
        case [null]:
            switch (input_stranded) {
                case ["paired-end"]:
                    trim_params = "ILLUMINACLIP:TruSeq3-PE-2.fa:2:30:10:8:keepBothReads TRAILING:28 MINLEN:40"
                    break
                case ["single-end"]:
                    trim_params = "ILLUMINACLIP:TruSeq3-SE.fa:2:30:10 TRAILING:28 MINLEN:40"
                    break
            }
            break
        default:
            trim_params = input_trim
            break
    }
    return trim_params
}

// // TRIMMED
// def getTrimmedFASTQ(input_dir,input_stranded) {
//     ext = "fastq,fastq.gz,fastq.bz2,fq,fq.gz,fq.bz2"
//     switch (input_stranded) {
//         case ["paired-end"]:
//             samples = Channel.fromFilePairs("${input_dir}/**_trimmed_{1P,2P}.{${ext}}", type: 'file')
//                 .ifEmpty { exit 1, trim_data_error }
//             break
//         case["single-end"]:
//             samples = Channel.fromFilePairs("${input_dir}/**_trimmed.{${ext}}", type: 'file', size:1)
//                 .ifEmpty { exit 1, trim_data_error(input_dir) }
//             break
//     }
// }

// summary="nf-rnaSeqCount v0.2 - Execution Summary:"
// workflow.onComplete {
//     println "\n${line}"
//     println "#".multiply(48 - ("${summary}".size() / 2 )) + "  ${summary}  " + "#".multiply(48 - ("${summary}".size() / 2 ))    
//     println "${line}"
//     println "Execution command   : ${workflow.commandLine}"
//     println "Execution name      : ${workflow.runName}"
//     println "Workflow start      : ${workflow.start}"
//     println "Workflow end        : ${workflow.complete}"
//     println "Workflow duration   : ${workflow.duration}"
//     println "Workflow completed? : ${workflow.success}"
//     println "Work directory      : ${workflow.workDir}"
//     println "Project directory   : ${workflow.projectDir}"
//     println "Execution directory : ${workflow.launchDir}"
//     println "Configuration files : ${workflow.configFiles}"
//     println "Workflow containers : ${workflow.container}"
//     println "exit status         : ${workflow.exitStatus}"
//     println "Error report        : ${workflow.errorReport ?: '-'}"
//     println "${line}\n"
//     println "\n"
// }

// workflow.onError {
//     println "Oohhh DANG IT!!... Pipeline execution stopped with the following message: ${workflow.errorMessage}"
// }

