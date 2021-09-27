#!/usr/bin/env nextflow
println "clear".execute().text

/*  ======================================================================================================
 *  HELP MENU
 *  ======================================================================================================
 */
line="=".multiply(100)
ver="nf-rnaSeqCount v0.2"
if (params.help) {
    println "\n${line}"
    println "#".multiply(48 - ("${ver}".size() / 2 )) + "  ${ver}   " + "#".multiply(48 - ("${ver}".size() / 2 ))
    println "${line}\n"
    println "USAGE:"
    println "nextflow run nf-rnaSeqCount -profile \"slurm\" --data \"/path/to/data\" --genome \"/path/to/genome.fa\" --genes \"/path/to/genes.gtf\"\n" 
    println "HELP:"
    println "nextflow run nf-rnaSeqCount --help\n"
    println "MANDATORY ARGUEMENTS:"
    println "-profile     STRING    Executor to be used. Available options:"
    println "\t\t\t\t\"standard\"          : Local execution (no job scheduler)."
    println "\t\t\t\t\"slurm\"             : SLURM scheduler."
    println "--mode       STRING    To specify which step of the workflow you are running (see https://github.com/phelelani/nf-rnaSeqCount)."
    println "                       Availeble options:"
    println "\t\t\t\t\"prep.Containers\"   : For downloading Singularity containers used in this workflow."
    println "\t\t\t\t\"prep.Indexes\"      : For indexing your reference genome using STAR and Bowtie2."
    println "\t\t\t\t\"run.ReadQC\"        : For performing general QC on your reads using FastQC. "
    println "\t\t\t\t\"run.ReadTrimming\"  : For trimming low quality bases and removing adapters from your reads using Trimmmomatic."
    println "\t\t\t\t\"run.ReadAlignment\" : For aligning your reads to your reference genome using STAR."
    println "\t\t\t\t\"run.ReadCounting\"  : For counting features in your reads using HTSeq-count and featureCounts."
    println "\t\t\t\t\"run.MultiQC\"       : For getting a summary of QC through the analysis using MultiQC.\n"
    println "--data       FOLDER    Path to where the input data (FASTQ files) is located. Supported FASTQ files:"
    println "\t\t\t\t[ fastq | fastq.gz | fastq.bz2 | fq | fq.gz | fq.bz2 ]"
    println "--genome     FILE      The whole genome FASTA sequence. Supported FASTA files:"
    println "\t\t\t\t[ fasta | fa | fna ]"
    println "--genes      FILE      The genome annotation GFT file. Supported GTF file:"
    println "\t\t\t\t[ gtf ]"
    println "OPTIONAL ARGUEMENTS:"
    println "--help                 To show this menu."
    println "--out        FOLDER    Path to where the output should be directed."
    println "                       Default: \$PWD/results_nf-rnaSeqCount)."
    println "--from       STRING    Specify to resume workflow from the QC or trimming step. Options:"
    println "\t\t\t\t\"run.ReadQC\"        : To resume from the QC step (default)."
    println "\t\t\t\t\"run.ReadTrimming\"  : To resume from the trimming step."
    println "--pairedEnd            If working with paired-end FASTQ files (default)."
    println "--singleEnd            If working with single-end FASTQ files."
    println "--trim       STRING    Parameters for Trimmomatic. See http://www.usadellab.org/cms/index.php?page=trimmomatic for a more detailed use."
    println "                       The default parameters for Trimmomatic I have given you here (for both paird- and single-end sequences) are:"
    println "\t\t\t\tFor paired-end: \"ILLUMINACLIP:TruSeq3-PE-2.fa:2:30:10:8:true TRAILING:28 MINLEN:40\""
    println "\t\t\t\tFor single-end: \"ILLUMINACLIP:TruSeq3-SE.fa:2:30:10:8:true TRAILING:28 MINLEN:40\""
    println "--max_memory STRING    Maximum memory you have access to."
    println "                       Default: \"200.GB\""
    println "--max_cpus   STRING    Maximum CPUs you have access to."
    println "                       Default: \"24\""
    println "--max_time   STRING    Maximum time you have access to."
    println "                       Default: \"24.h\""
    println "${line}\n"
    exit 1
}

/*  ======================================================================================================
 *  CHECK ALL USER INPUTS
 *  ======================================================================================================
 */

// MAIN USER INPUT ERRORS
data_error = """
${line}
Oooh no!! Looks like there's a serious issue in your command! 
I do not recognise the \'--data ${params.data}\' option you have given me, or you have not given me any \'--data\' option at all!
Please provide a valid directory with you input FASTQ reads with the \'--data\' option to run the nf-rnaSeqCount workflow! 
${line}
"""

genome_error = """
${line}
Oooh no!! Looks like there's a serious issue in your command! 
I do not recognise the \'--genome ${params.genome}\' option you have given me, or you have not given me any \'--genome\' option at all!
Please provide a valid FASTA file (.fasta or .fa) for your reference genome with the \'--genome\' option to run the nf-rnaSeqCount workflow! 
${line}
"""

genes_error = """
${line}
Oooh no!! Looks like there's a serious issue in your command! 
I do not recognise the \'--genes ${params.genes}\' option you have given me, or you have not given me any \'--genes\' option at all!
Please provide a valid GTF annotation file (.gtf) for your reference genome with the \'--genes\' option to run the nf-rnaSeqCount workflow! 
${line}
"""

mode_error = """
${line}
Oooh no!! Looks like there's an serious issue in your command! 
I do not recognise the \'--mode ${params.mode}\' option you have given me, or you have not given me any \'--mode\' option at all!
The allowed options for \'--mode\' are:
\tprep.Containers\t\t: For downloading Singularity containers used in this workflow.
\tprep.Indexes\t\t: For indexing your reference genome using STAR and Bowtie2.
\trun.ReadQC\t\t: For performing general QC on your reads using FastQC. 
\trun.ReadTrimming\t: For trimming low quality bases and removing adapters from your reads using Trimmmomatic.
\trun.ReadAlignment\t: For aligning your reads to your reference genome using STAR.
\trun.ReadCounting\t: For counting features in your reads using HTSeq-count and featureCounts.
\trun.MultiQC\t\t: For getting a summary of QC through the analysis using MultiQC.
\nPlease use one of the above options with \'--mode\' to run the nf-rnaSeqCount workflow!
${line}
"""
from_error = """
${line}
Oooh no!! Looks like there's an serious issue in your command! 
I do not recognise the \'--from ${params.from}\' option you have given me!
The allowed options for \'--from\' are:
\trun.ReadQC\t\t: To resume from the QC step.
\trun.ReadTrimming\t: To resume from the trimming step.
\nPlease use one of the above options with \'--from\' to run the nf-rnaSeqCount workflow!
${line}
"""

// EMPTY LIST FOR COLLECTING ALL THE PATHS TO BIND TO SINGULARITY IMAGE
bind_dirs = []

// USER PARAMETER INPUT: DATA DIRECTORY
switch (params.data) {
    case [null]:
        data_dir = "NOT SPECIFIED!"
        break
    default:
        data_dir = file(params.data, type: 'dir')
        bind_dirs.add(data_dir)
        break
}

// USER PARAMETER INPUT: GENOME FASTA FILE
switch (params.genome) {
    case [null]:
        genome = "NOT SPECIFIED!"
        break
    default:
        genome = file(params.genome, type: 'file', checkIfExists: true)
        index_dir = genome.getParent()
        bind_dirs.add(genome.getParent())
        break
}

// USER PARAMETER INPUT: GENOME ANNOTATION FILE (GFT/GFF)
switch (params.genes) {
    case [null]:
        genes = "NOT SPECIFIED!"
        break
    default:
        genes = file(params.genes, type: 'file', checkIfExists: true)
        bind_dirs.add(genes.getParent())
        break
}

// USER PARAMETER INPUT: OUTPUT DIRECTORY
switch (params.out) {
    case [null]:
        out_dir = file("${PWD}/results_nf-rnaSeqCount", type: 'dir')
        break
    default:
        out_dir = file(params.out, type: 'dir')
        bind_dirs.add(out_dir)
        break
}

// USER INPUT RESUME FROM: WHERE TO PICK UP??
if(params.from == null) {
    resume_from = params.from
} else if(params.from in ["run.ReadTrimming", "run.ReadQC"]) {
    resume_from = params.from
} else {
    exit 1, "$from_error"
}


// USER STRANDED MODE: ARE WE DOING PAIRED- OR SINGLE-END?
if (params.singleEnd == null && params.pairedEnd == null) {
    stranded = "paired-end"
} else if(params.singleEnd) {
    stranded = "single-end"
} else if(params.pairedEnd){
    stranded = "paired-end"
} else {}

def breakIfNull(parameter,error) {
    if (parameter == null) {
        exit 1, error
    } else {}
}

// USER PARAMETER INPUT: TRIMMOMATIC OPTIONS
switch (params.trim) {
    case [null]:
        switch (stranded) {
            case ["paired-end"]:
                trim_params = "ILLUMINACLIP:TruSeq3-PE-2.fa:2:30:10:8:keepBothReads TRAILING:28 MINLEN:40"
                break
            case ["single-end"]:
                trim_params = "ILLUMINACLIP:TruSeq3-SE.fa:2:30:10 TRAILING:28 MINLEN:40"
                break
        }
        break
    default:
        trim_params = params.trim
        break
}

// OUTPUT DIRECTORIES
qc_dir        = file("${out_dir}/1_Read_QC", type: 'dir')
trim_dir      = file("${out_dir}/2_Read_Trimming", type: 'dir')
align_dir     = file("${out_dir}/3_Read_Alignment", type: 'dir')
counts_dir    = file("${out_dir}/4_Read_Counts", type: 'dir')
multiqc_dir   = file("${out_dir}/5_MultiQC", type: 'dir')

ext = "fastq,fastq.gz,fastq.bz2,fq,fq.gz,fq.bz2"

// DATA ABSENT ERRORS
main_data_error = """
${line}
Oooh no!! Looks like there's a serious issue in your input data! There are no FASTQ file in the directory:
\t${data_dir}
Please ensure that you have given me the correct directory for you FASTQ input reads using the \'--data\' option!
${line}
"""

trim_data_error = """
${line}
Oooh no!! Looks like there's a serious issue in your input data! There are no FASTQ file in the directory:
\t${trim_dir}
Are you sure you ran the READ TRIMMING STEP using \'--mode run.ReadTrimming\' ??
Please ensure that you have ran the READ TRIMMING STEP successfully and try again!
${line}
"""

bams_error = """
Ooops!! Looks like there's a serious issue in your input data! There are no BAM files in the directory:
\t${align_dir}
Are you sure you ran the READ ALIGNMENT STEP using \'--mode run.ReadAlignment\' ??
Please ensure that you have ran the READ ALIGNMENT STEP successfully and try again!
"""

// USER INPUT MODE: WHICH ANALYSIS TO RUN!
switch (params.mode) {
    case [null]:
        exit 1, "$mode_error"
    
    case ["prep.Containers", "prep.Indexes"]:
        mode = params.mode
        switch (mode) {
            case ["prep.Containers"]:
                break
            case ["prep.Indexes"]:
                breakIfNull(params.genome,"$genome_error")
                breakIfNull(params.genes,"$genes_error")
                break
        }
        break

    case ["run.ReadQC", "run.ReadTrimming", "run.ReadAlignment", "run.ReadCounting", "run.MultiQC"]:
        mode = params.mode
        switch (mode) {
            case ["run.ReadQC"]:
                breakIfNull(params.data,"$data_error")
                // GET DATA BASED ON THE STRANDEDNESS
                switch (stranded) {
                    case ["paired-end"]:
                        read_pairs = Channel.fromFilePairs("${data_dir}/*.{${ext}}", type: 'file', size:-1) { "all_reads" }
                            .ifEmpty { exit 1, "$main_data_error" }
                        break
                    case ["single-end"]:
                        read_pairs = Channel.fromFilePairs("${data_dir}/*.{${ext}}", type: 'file', size:-1) { "all_reads" }
                            .ifEmpty { exit 1, "$main_data_error" }
                        break
                }
                break
            case ["run.ReadTrimming"]:
                breakIfNull(params.data,"$data_error")
                // GET DATA BASED ON THE STRANDEDNESS
                switch (stranded) {
                    case ["paired-end"]:
                        read_pairs = Channel.fromFilePairs("${data_dir}/*{R,read,_}[1,2]*.{${ext}}", type: 'file')
                            .ifEmpty { exit 1, "$main_data_error" }
                        break
                    case["single-end"]:
                        read_pairs = Channel.fromFilePairs("${data_dir}/*.{${ext}}", type: 'file', size:1)
                            .ifEmpty { exit 1, "$main_data_error" }
                        break
                }
                break
            case ["run.ReadAlignment"]:
                breakIfNull(params.genome,"$genome_error")
                breakIfNull(params.genes,"$genes_error")
                switch (resume_from) {
                    case [null]:
                        breakIfNull(params.data,"$data_error")
                        // GET DATA BASED ON THE STRANDEDNESS
                        switch (stranded) {
                            case ["paired-end"]:
                                read_pairs = Channel.fromFilePairs("${data_dir}/*{R,read,_}[1,2]*.{${ext}}", type: 'file')
                                    .ifEmpty { exit 1, "$main_data_error" }
                                break
                            case ["single-end"]:
                                read_pairs = Channel.fromFilePairs("${data_dir}/*.{${ext}}", type: 'file', size:1)
                                    .ifEmpty { exit 1, "$main_data_error" }
                                break
                        }
                        break
                        
                    case ["run.ReadTrimming"]:
                        // GET DATA BASED ON THE STRANDEDNESS
                        switch (stranded) {
                            case ["paired-end"]:
                                read_pairs = Channel.fromFilePairs("${trim_dir}/*[1,2]P*.{${ext}}", type: 'file')
                                    .ifEmpty { exit 1, "$trim_data_error" }
                                break
                            case ["single-end"]:
                                read_pairs = Channel.fromFilePairs("${trim_dir}/*.{${ext}}", type: 'file', size:1)
                                    .ifEmpty { exit 1, "$trim_data_error" }
                                break
                        }
                        break
                }
                break
            case ["run.ReadCounting"]:                
                bams = Channel.fromFilePairs("${align_dir}/**_Aligned.out.bam", size:-1) { 
                    file -> "${file.baseName.replace(/_Aligned.out/, "")}" }
                    .ifEmpty { exit 1, "$bams_error" }
                bams.into { bams_htseqCounts; bams_featureCounts}
                break
            case ["run.MultiQC"]:
                break
        }
        out_dir.mkdir()
        break
    default:
        exit 1, "$mode_error"
        break
}

// USER PARAMETER INPUT: PATHS TO BE BINDED TO THE IMAGE
bind_dirs = bind_dirs
    .unique()
    .collect { it -> "-B ${it}"}
    .join("\n" + ' '.multiply(26))
    .toString()

//  ======================================================================================================
//  RUN INFO
//  ======================================================================================================
options="nf-rnaSeqCount v0.2 - Input/Output and Parameters:"
println "\n" + "=".multiply(100)
println "#".multiply(48 - ("${options}".size() / 2 )) + "  ${options}  " + "#".multiply(48 - ("${options}".size() / 2 ))
println "=".multiply(100)
println "Input data              : $data_dir"
println "Input data type         : $stranded"
println "Output directory        : $out_dir"
println ' '.multiply(26) + "- ${qc_dir.baseName}"
println ' '.multiply(26) + "- ${trim_dir.baseName}"
println ' '.multiply(26) + "- ${align_dir.baseName}"
println ' '.multiply(26) + "- ${counts_dir.baseName}"
println ' '.multiply(26) + "- ${multiqc_dir.baseName}"
println "Genome                  : $genome"
println "Genome annotation       : $genes"
println "Trimmomatic parameters  : $trim_params"
println "Paths to bind           : $bind_dirs"
println "=".multiply(100)
println " "
// ======================================================================================================
//  PIPELINE START
// ======================================================================================================


// ========== THIS SECTION IS FOR PREPPING DATA (SINGULARITY IMAGES, STAR INDEXES AND BOWTIE INDEXES)
switch (mode) {
        //
    case ['prep.Containers']: 
        base = "docker://phelelani/nf-rnaseqcount:"
        images = ["star", "htseqcount", "featurecounts", "multiqc", "bowtie2", "fastqc", "trimmomatic"]
        
        process run_DownloadContainers {
            label 'mini'
            tag { "Downloading: ${base}${image}" }
            maxForks 1

            input:
            each image from images
            
            """
            singularity pull --force --dir \$HOME/.singularity/cache/ ${base}${image}
            """
        }
        break
        // ==========
        
        //
    case ['prep.Indexes']:
        process run_GenerateSTARIndex {
            label 'maxi'
            tag { "Generate Star Index" }
            publishDir "$index_dir", mode: 'copy', overwrite: true
            
            output:
            set val("starIndex"), file("*") into star_index
            
            """
            STAR --runThreadN ${task.cpus} \
                --runMode genomeGenerate \
                --genomeDir . \
                --genomeFastaFiles ${genome} \
                --sjdbGTFfile ${genes} \
                --sjdbOverhang 99
            """
        }

        process run_GenerateBowtie2Index {
            label 'maxi'
            tag { "Generate Bowtie2 Index" }
            publishDir "$index_dir", mode: 'copy', overwrite: true
        
            output:
            set val("bowtieIndex"), file("*") into bowtie_index
            
            """
            bowtie2-build --threads ${task.cpus} ${genome} genome
            """
        }       
        break
        // ========== PREPPING STEPS/OPTIONS END HERE!


        // MAIN WORKFLOW - STEP 1 (OPTIONAL): PERFORM QC ON INPUT FASTQ FILES!
    case['run.ReadQC']: 
        process run_QualityChecks {
            label 'midi'
            tag { samples }
            publishDir "${qc_dir}", mode: 'copy', overwrite: true
            
            input:
            set val(samples), file(reads) from read_pairs
            
            output:
            set val(samples), file("*.{html,zip}") into qc_html

            """
            fastqc ${reads.findAll().join(' ') } \
                --threads ${task.cpus} \
                --noextract
            """
        }
        break
        // --------------------

        // MAIN WORKFLOW - STEP 2 (OPTIONAL): TRIMMING OF INPUT FASTQ FILES
    case['run.ReadTrimming']:
        process run_ReadTrimming {
            label 'maxi'
            tag { sample }
            publishDir "${trim_dir}", mode: 'copy', overwrite: true
            
            input:
            set sample, file(reads) from read_pairs
            
            output:
            set sample, file("${sample}_trimmed*.fastq.gz") into read_pairs_trimmed

            """
            ln -s /opt/Trimmomatic-0.39/adapters/*.fa .

            if [[ ${stranded} == "paired-end" ]]
            then
                java -jar /opt/Trimmomatic-0.39/trimmomatic-0.39.jar PE \
                    -threads ${task.cpus} \
                    ${reads.findAll().join(' ')} \
                    -baseout ${sample}_trimmed.fastq.gz \
                    ${trim_params}
            elif [[ ${stranded} == "single-end" ]]
            then
                java -jar /opt/Trimmomatic-0.39/trimmomatic-0.39.jar SE \
                    -threads ${task.cpus} \
                    ${reads.findAll().join(' ')} \
                    ${sample}_trimmed.fastq.gz \
                    ${trim_params}
            else :
            fi
            """
        }
        break
        // --------------------        

    case['run.ReadAlignment']:
        process run_STAR {
            label 'maxi'
            tag { sample }
            publishDir "${align_dir}", mode: 'copy', overwrite: true

            input:
            set sample, file(reads) from read_pairs
    
            output:
            set sample, file("${sample}*.{out,tab}"), file("${sample}_Aligned.out.bam") into star_alignments
            
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
        break
        // ==========
        
        case['run.ReadCounting']:
        // USE HTSEQCOUNTS TO GET RAW READ COUNTS
        process run_HTSeqCount {
            label 'mini'
            tag { sample }
            publishDir "${counts_dir}/htseqCounts", mode: 'copy', overwrite: true

            input:
            set sample, file(bam) from bams_htseqCounts

            output:
            set sample, "${sample}.txt" into htseqCounts
            set sample, "${sample}.txt" into htseqCounts_cleanup
    
            """
            htseq-count -f bam \
                -r name \
                -i gene_id \
                -a 10 \
                -s reverse \
                -m union \
                -t exon \
                $bam $genes > ${sample}.txt
            """
        }

        // COLLECT FILES FOR HTSEQ CLEANUP
        htseqCounts_cleanup
            .collectFile () { item -> [ 'cleanup_htseqcounts.txt', "${item.get(1)}" + '\n' ] }
            .set { list_htseqcounts }

        // CLEANUP HTSEQCOUNTS
        process run_CleanHTSeqCounts {
            label 'mini'
            tag { 'HTSEqCounts Cleanup' }
            publishDir "${counts_dir}/htseqCounts", mode: 'copy', overwrite: false
    
            input:
            file(file_list) from list_htseqcounts

            output:
            file(out_file) into htseqCounts_cleaned

            shell:
            in_file = "${file_list}"
            out_file = "gene_counts_final.txt"
            template "clean_htseqCounts.sh"
        }

        // GET RAW COUNTS USING FEATURECOUNTS
        // GET ALL THE BAM FILE LOCATIONS TO PROCESS WITH FEATURECOUNTS
        bams_featureCounts
            .collectFile () { item -> [ 'sample_bams.txt', "${item.get(1).join().toString()}" + ' ' ] }
            .set { sample_bams }
        

        // USE FEATURECOUNTS TO GET RAW GENE COUNTS
        process run_FeatureCounts {
            label 'maxi'
            tag { 'featureCounts - ALL' }
            publishDir "${counts_dir}/featureCounts", mode: 'copy', overwrite: false

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
                -T ${task.cpus} \
                -o gene_counts.txt \
                `< ${samples}`
            """
        }

        // CLEANUP FEATURECOUNTS GENE COUNTS
        process run_CleanFeatureCounts {
            label 'mini'
            tag { 'featureCounts Cleanup' }
            publishDir "${counts_dir}/featureCounts", mode: 'copy', overwrite: false
            
            input:
            file(raw_counts) from featureCounts_raw

            output:
            file(out_file) into featureCounts_cleaned
            
            shell:
            in_file = "${raw_counts}"
            out_file = "gene_counts_final.txt"
            template "clean_featureCounts.sh"
        }
        break
        // ==========

    case['run.MultiQC']:
        process run_MultiQC {
            label 'mini'
            tag { 'MultiQC - ALL' }
            publishDir "${multiqc_dir}", mode: 'copy', overwrite: false
            
            output:
            file('*') into multiQC
            
            """
            multiqc ${out_dir} --force
            """
        }
        break
        // ==========
}


// ======================================================================================================
//  WORKFLOW SUMMARY
//  ======================================================================================================
summary="nf-rnaSeqCount v0.2 - Execution Summary:"
workflow.onComplete {
    println "\n${line}"
    println "#".multiply(48 - ("${summary}".size() / 2 )) + "  ${summary}  " + "#".multiply(48 - ("${summary}".size() / 2 ))    
    println "${line}"
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
    println "exit status         : ${workflow.exitStatus}"
    println "Error report        : ${workflow.errorReport ?: '-'}"
    println "${line}\n"
    println "\n"
}

workflow.onError {
    println "Oohhh DANG IT!!... Pipeline execution stopped with the following message: ${workflow.errorMessage}"
}
//======================================================================================================
