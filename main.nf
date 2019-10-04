#!/usr/bin/env nextflow

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
    log.info '    --out       FOLDER     Path to where the output should be directed (will be created if it does not exist).'
    log.info '    --genome    FILE       The whole genome sequence (fasta | fa | fna)'
    log.info '    --genes     FILE       The genome annotation file (gtf)'
    log.info '     -profile   SRTING     Executor to be used'
    log.info ''
    log.info "====================================\n"
    exit 1
}


/*  ======================================================================================================
 *  CHECK ALL USER INPUTS
 *  ======================================================================================================
 */

// MAIN USER INPUT ERRORS
line="=".multiply(100)

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
\tprep.STARIndex\t\t: For indexing your reference genome using STAR.
\tprep.BowtieIndex\t: For indexing your reference genome using Bowtie2.
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

// USER PARAMETER INPUT: DATA DIRECTORY
if (params.data == null) {
    exit 1, "$data_error"
} else{
    data_dir = file(params.data, type: 'dir')
}

// USER PARAMETER INPUT: OUTPUT DIRECTORY
if(params.out == null) {
    out_dir = file("${PWD}/resultsf-rnaSeqCount", type: 'dir')
} else{
    out_dir = file(params.out, type: 'dir')
}

// USER PARAMETER INPUT: GENOME FASTA FILE
if(params.genome == null) {
    exit 1, "$genome_error"
} else{
    genome = file(params.genome, type: 'file')
    index = new File(params.genome).getParent()
}

// USER PARAMETER INPUT: GENOME ANNOTATION FILE (GFT/GFF)
if(params.genes == null) {
    exit 1, "$genes_error"
} else{
    genes = file(params.genes, type: 'file') 
}

// USER INPUT MODE: WHICH ANALYSIS TO RUN!
if(params.mode == null) {
    exit 1, "$mode_error"
} else {
    mode = params.mode
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
if(params.singleEnd == null && params.pairedEnd == null) {
    stranded = "paired-end"
} else if(params.singleEnd) {
    stranded = "single-end"
} else if(params.pairedEnd){
    stranded = "paired-end"
} else {}

// USER PARAMETER INPUT: PATHS TO BE BINDED TO THE IMAGE
bind_dir      = [params.data, params.out, params.genome.split('/')[0..-2].join('/'), params.genes.split('/')[0..-2].join('/')]
    .unique()
    .collect { it -> "-B ${it}"}
    .join("\n" + ' '.multiply(26))
    .toString()
trim_params   = params.trim

// OUTPUT DIRECTORIES
out_dir.mkdir()
qc_dir        = file("${out_dir}/1_Read_QC", type: 'dir')
trim_dir      = file("${out_dir}/2_Read_Trimming", type: 'dir')
align_dir     = file("${out_dir}/3_Read_Alignment", type: 'dir')
counts_dir    = file("${out_dir}/4_Read_Counts", type: 'dir')
multiqc_dir   = file("${out_dir}/5_MultiQC", type: 'dir')
ext           = "fastq,fastq.gz,fastq.bz2,fq,fq.gz,fq.bz2"

//  ======================================================================================================
//  RUN INFO
//  ======================================================================================================

println "\n" + "=".multiply(100)
println " ".multiply(45) + "nf-rnaSeqCount v0.2"
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
println "Paths to bind           : $bind_dir"
println "=".multiply(100)

// ======================================================================================================
//  PIPELINE START
// ======================================================================================================

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

// GET INPUT DATA DEPENDING ON USER "MODE"
if(mode in ["prep.Containers", "prep.STARIndex", "prep.BowtieIndex"]) {
    // OPTIONS FOR PREPARING DATA
    switch (mode) {
        case ["prep.Containers"]:
            
            break
        case ["prep.STARIndex","prep.BowtieIndex"]:
            index_dir = genome.getParent()
            println index_dir
            break
    }
} else if(mode in ["run.ReadQC", "run.ReadTrimming", "run.ReadAlignment", "run.ReadCounting", "run.MultiQC"]) {
    // OPTIONS FOR PERFORMING THE ANALYSES
    switch (mode) {
        case["run.ReadQC"]:
            if(stranded == "paired-end") {
                read_pairs = Channel.fromFilePairs("${data_dir}/*{R,read}[1,2]*.{${ext}}", type: 'file')
                    .ifEmpty { exit 1, "$main_data_error" }
            } else if(stranded == "single-end") {
                read_pairs = Channel.fromFilePairs("${data_dir}/*.{${ext}}", type: 'file', size:1)
                    .ifEmpty { exit 1, "$main_data_error" }
            }
            break
        case["run.ReadTrimming"]:
            if(stranded == "paired-end") {
                read_pairs = Channel.fromFilePairs("${data_dir}/*{R,read}[1,2]*.{${ext}}", type: 'file')
                    .ifEmpty { exit 1, "$main_data_error" }
            } else if(stranded == "single-end") {
                read_pairs = Channel.fromFilePairs("${data_dir}/*.{${ext}}", type: 'file', size:1)
                    .ifEmpty { exit 1, "$main_data_error" }
            }
            break
        case["run.ReadAlignment"]:
            if(resume_from == null) {
                if(stranded == "paired-end") {
                    read_pairs = Channel.fromFilePairs("${data_dir}/*{R,read}[1,2]*.{${ext}}", type: 'file')
                        .ifEmpty { exit 1, "$main_data_error" }
                } else if(stranded == "single-end") {
                    read_pairs = Channel.fromFilePairs("${data_dir}/*.{${ext}}", type: 'file', size:1)
                        .ifEmpty { exit 1, "$main_data_error" }
                }
            } else if(resume_from == "run.ReadTrimming") {
                if(stranded == "paired-end") {
                    read_pairs = Channel.fromFilePairs("${trim_dir}/*[1,2]P*.{${ext}}", type: 'file')
                        .ifEmpty { exit 1, "$trim_data_error" }
                } else if(stranded == "single-end") {
                    read_pairs = Channel.fromFilePairs("${trim_dir}/*.{${ext}}", type: 'file', size:1)
                        .ifEmpty { exit 1, "$trim_data_error" }
                }
            }
            break
        case["run.ReadCounting"]:
            Channel.fromFilePairs("${align_dir}/**/*_Aligned.out.bam", size:-1) { file -> "${file.baseName.replace(/_Aligned.out/, "")}" }
                .into { bams_htseqCounts; bams_featureCounts; check }
            //.ifEmpty { exit 1, "$bams_error" }
            check.view()
            break
        case["run.MultiQC"]:
            // CHECK IF THE OUTPUT DIRECTORY IS EMPTY! GIVE ERROR IF ITS EMPTY
            
            break
    }
} else {
    exit 1, "$mode_error"
}

// ========== THIS SECTION IS FOR PREPPING DATA (SINGULARITY IMAGES, STAR INDEXES AND BOWTIE INDEXES)
switch (mode) {
        //
    case ['prep.Containers']: // <<<<<< WORKS ALL GOOD HERE!
        base = "shub://phelelani/nf-rnaSeqCount:"
        images = Channel.from( ["${base}star", "${base}htseqcount", "${base}featurecounts", "${base}multiqc", "${base}trinity", "${base}fastqc", "${base}trimmomatic"] )
        
        process run_DownloadContainers {
            label 'mini'
            scratch '$HOME/tmp'
            tag { "Downloading: ${link}" }
            publishDir "$PWD/containers", mode: 'copy', overwrite: true
            
            input:
            each link from images
            
            output:
            file("*.sif") into containers
            
            """
            singularity pull nf-rnaSeqCount-${link.substring(32,)}.sif ${link}
            """
        }

        containers.subscribe { println "${it}" }
        break
        // ==========
        
        //
    case ['prep.STARIndex']: /// <<<<<<< WORKS ALL GOOD HERE!
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

        star_index.subscribe { println "${it}" }
        break
        // ==========
        
        //
    case ['prep.BowtieIndex']: // <<<<<< WORKS! ALL GOOD
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
    
        bowtie_index.subscribe { println "${it}" }
        break
        // ========== PREPPING STEPS/OPTIONS END HERE!


        // MAIN WORKFLOW - STEP 1 (OPTIONAL): PERFORM QC ON INPUT FASTQ FILES!
    case['run.ReadQC']: // wWORKS FINE!
        // process run_QualityChecks {
        //     label 'midi'
        //     tag { sample }
        //     publishDir "${qc_dir}", mode: 'copy', overwrite: true
            
        //     input:
        //     set sample, file(reads) from read_pairs
            
        //     output:
        //     set sample, file("${sample}*.html") into qc_html
        //     set sample, file("${sample}*.zip") into qc_multiqc

        //     """
        //     fastqc ${reads.get(0)} ${reads.get(1)} \
        //         --threads ${task.cpus} \
        //         --noextract
        //     """
        // }

        // qc_html.subscribe { println "${it}" }
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
            set sample, file("${sample}*{1,2}P*") into read_pairs_trimmed

            """
            java -jar /opt/Trimmomatic-0.39/trimmomatic-0.39.jar PE \
                -threads ${task.cpus} \
                -trimlog trimlog_${sample}.log \
                ${reads.get(0)} ${reads.get(1)} \
                -baseout ${sample}_trimmed.fastq.gz \
                \$(sed 's|ILLUMINACLIP:|ILLUMINACLIP:/opt/Trimmomatic-0.39/adapters/|' <<< "${trim_params}")
            """
        }
        break
        // --------------------        

// ${read_file_cmd} \
// def checkFileExtension(){
//     switch() {
//         case ['fastq','fq']:
//             read_file_cmd = ''
//             break
//         case ['fastq.gz','fq.gz']:
//             read_file_cmd = '--readFilesCommand gunzip -c'
//             break
//         case ['fastq.bz2','fq.bz2']:
//             read_file_cmd = '--readFilesCommand bunzip2 -c'
//             break
//         case null:
//             read_file_cmd = '--readFilesCommand gunzip -c'
//             break
//     }
// }

    case['run.ReadAlignment']: // <<<<< ALL GOOD - MUST FIX THE READFILE COMMAND
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
                --genomeDir ${index} \
                --readFilesCommand gunzip -c \
                --readFilesIn ${reads.get(0)} ${reads.get(1)} \
                --runThreadN ${task.cpus} \
                --outSAMtype BAM Unsorted \
                --outFileNamePrefix ${sample}_
            """
        }
        
        star_alignments.subscribe { println "${it}" }
        break
        // ==========
        
        case['run.ReadCounting']:
        // GET INPUT DATA
        // Channel.fromFilePairs("$align_dir/**_Aligned.out.bam", size:-1) { 
        //     file -> "${file.baseName.replace(/_Aligned.out/, "")}" 
        // }.into { bams_htseqCounts; bams_featureCounts }

        // USE HTSEQCOUNTS TO GET RAW READ COUNTS
        process run_HTSeqCount {
            label 'mini'
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
            label 'midi'
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

        featureCounts_cleaned.view()
        break
        // ==========

    case['run.MultiQC']:
        process run_MultiQC {
            label 'mini'
            tag { 'MultiQC - ALL' }
            publishDir "${multiqc_dir}/report_QC", mode: 'copy', overwrite: false
            
            output:
            file('*') into multiQC
            
            """
            multiqc ${out_dir} --force
            """
        }

        multiQC.view()
        break
        // ==========
}


// ======================================================================================================
//  WORKFLOW SUMMARY
//  ======================================================================================================
workflow.onComplete {
    println "\n${line}"
    println "Pipeline execution summary:"
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
    println "exit status : ${workflow.exitStatus}"
    println "Error report: ${workflow.errorReport ?: '-'}"
    println "${line}"
}

workflow.onError {
    println "Oohhh DANG IT!!... Pipeline execution stopped with the following message: ${workflow.errorMessage}"
}
//======================================================================================================
