#!/usr/bin/env nextflow
println "clear".execute().text
//  DO NOT EDIT FROM HERE!! - Unless you brave like King Shaka of course! 
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
    println "--data       FOLDER    Path to where the input data (FASTQ files) is located. Supported FASTQ files:"
    println "\t\t\t\t[ fastq | fastq.gz | fastq.bz2 | fq | fq.gz | fq.bz2 ]"
    println "--genome     FILE      The whole genome FASTA sequence. Supported FASTA files:"
    println "\t\t\t\t[ fasta | fa | fna ]"
    println "--genes      FILE      The genome annotation GFT file. Supported GTF file:"
    println "\t\t\t\t[ gtf ]"
    println "--mode       STRING    To specify which step of the workflow you are running (see https://github.com/phelelani/nf-rnaSeqCount)."
    println "                       Availeble options:"
    println "\t\t\t\t\"prep.Containers\"   : For downloading Singularity containers used in this workflow."
    println "\t\t\t\t\"prep.STARIndex\"    : For indexing your reference genome using STAR."
    println "\t\t\t\t\"prep.BowtieIndex\"  : For indexing your reference genome using Bowtie2."
    println "\t\t\t\t\"run.ReadQC\"        : For performing general QC on your reads using FastQC. "
    println "\t\t\t\t\"run.ReadTrimming\"  : For trimming low quality bases and removing adapters from your reads using Trimmmomatic."
    println "\t\t\t\t\"run.ReadAlignment\" : For aligning your reads to your reference genome using STAR."
    println "\t\t\t\t\"run.ReadCounting\"  : For counting features in your reads using HTSeq-count and featureCounts."
    println "\t\t\t\t\"run.MultiQC\"       : For getting a summary of QC through the analysis using MultiQC.\n"
    println "OPTIONAL ARGUEMENTS:"
    println "--help                 To show this menu."
    println "--out        FOLDER    Path to where the output should be directed (default: \$PWD/results_nf-rnaSeqCount)."
    println "--from       STRING    Specify to resume workflow from the QC or trimming step. Options:"
    println "\t\t\t\t\"run.ReadQC\"        : To resume from the QC step (default)."
    println "\t\t\t\t\"run.ReadTrimming\"  : To resume from the trimming step."
    println "--pairedEnd            If working with paired-end FASTQ files (default)."
    println "--singleEnd            If working with single-end FASTQ files."
    println "--trim       STRING    Parameters for Trimmomatic. See http://www.usadellab.org/cms/index.php?page=trimmomatic for a more detailed use."
    println "                       The default parameters for Trimmomatic I have given you here (for both paird- and single-end sequences) are:"
    println "\t\t\t\tFor paired-end: \"ILLUMINACLIP:TruSeq3-PE-2.fa:2:30:10:8:true TRAILING:28 MINLEN:40\""
    println "\t\t\t\tFor single-end: \"ILLUMINACLIP:TruSeq3-SE.fa:2:30:10:8:true TRAILING:28 MINLEN:40\""
    println "--max_memory STRING    Maximum memory you have access to (default: \"200.GB\")"
    println "--max_cpus   STRING    Maximum CPUs you have access to (default: \"24\")"
    println "--max_time   STRING    Maximum time you have access to(default: \"24.h\")"
    println "${line}\n"
    exit 1
}

from       = null
pairedEnd  = null
singleEnd  = null
trim       = null
max_memory = 200.GB
max_cpus   = 24
max_time   = 24.h

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
// ---- THESE DO NOT REQUIRE DATA!!
if (params.mode in [ "prep.Containers", "prep.STARIndex", "prep.BowtieIndex", "run.MultiQC" ]) {
    if (params.data == null) {
        data_dir = "YOU HAVEN'T SPECIFIED THE DITA DIRECTORY YET! PLEASE SPECIFY BEFORE RUNNING THE WORKFLOW"
    } else{
        data_dir = file(params.data, type: 'dir')
    }
} else if (params.data == null && prarams.mode in [ "run.ReadQC", "run.ReadTrimming", "run.ReadAlignment", "run.ReadCounting"]) {
    exit 1, "$data_error"
} else{
    data_dir = file(params.data, type: 'dir')
}

// USER PARAMETER INPUT: OUTPUT DIRECTORY
if(params.out == null) {
    out_dir = file("${PWD}/results_nf-rnaSeqCount", type: 'dir')
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
bind_dir      = [ params.data, out_dir, new File("${params.genome}").getParent(), new File("${params.genes}").getParent() ]
    .unique()
    .collect { it -> "-B ${it}"}
    .join("\n" + ' '.multiply(26))
    .toString()

if(params.trim == null) {
    if(stranded == "paired-end") {
        trim_params = "ILLUMINACLIP:TruSeq3-PE-2.fa:2:30:10:8:keepBothReads TRAILING:28 MINLEN:40"
    } else if(stranded == "single-end") {
        trim_params = "ILLUMINACLIP:TruSeq3-SE.fa:2:30:10 TRAILING:28 MINLEN:40"
    }
} else{
    trim_params = params.trim
}

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
println "Paths to bind           : $bind_dir"
println "=".multiply(100)
println " "
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
            break
    }
} else if(mode in ["run.ReadQC", "run.ReadTrimming", "run.ReadAlignment", "run.ReadCounting", "run.MultiQC"]) {
    // OPTIONS FOR PERFORMING THE ANALYSES
    switch (mode) {
        case["run.ReadQC"]:
            if(stranded == "paired-end") {
                read_pairs = Channel.fromFilePairs("${data_dir}/*.{${ext}}", type: 'file', size:-1) { "all_reads" }
                    .ifEmpty { exit 1, "$main_data_error" }
            } else if(stranded == "single-end") {
                read_pairs = Channel.fromFilePairs("${data_dir}/*.{${ext}}", type: 'file', size:-1) { "all_reads" }
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
            bams = Channel.fromFilePairs("${align_dir}/**_Aligned.out.bam", size:-1) { 
                file -> "${file.baseName.replace(/_Aligned.out/, "")}" }
                .ifEmpty { exit 1, "$bams_error" }
            bams.into { bams_htseqCounts; bams_featureCounts}
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
            maxForks 2
            publishDir "$PWD/containers", mode: 'copy', overwrite: true
            
            input:
            each link from images
            
            output:
            file("*.sif") into containers
            
            """
            singularity pull nf-rnaSeqCount-${link.substring(32,)}.sif ${link}
            """
        }
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

        star_index.subscribe { 
            println "\nSTAR index files generated:"
            it[1].each { 
                item -> println "\t${item}" 
            }
            println " "
        }
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
    
        bowtie_index.subscribe { 
            println "\nBowtie2 index files generated:"
            it[1].each { 
                item -> println "\t${item}" 
            }
            println " "
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
        qc_html.subscribe {
            println "\nFastQC files generated for all the samples:"
            it[1].each {
                item -> println "\t${item}" 
            }
            println " "
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
            set sample, file("${sample}*{1,2}P*") into read_pairs_trimmed

            """
            ln -s /opt/Trimmomatic-0.39/adapters/*.fa .

            if [[ ${stranded} === "paired-end" ]]
            then
                java -jar /opt/Trimmomatic-0.39/trimmomatic-0.39.jar PE \
                    ${reads.findAll().join(' ')} \
                    -threads ${task.cpus} \
                    -trimlog trimlog_${sample}.log \
                    -baseout ${sample}_trimmed.fastq.gz \
                    ${trim_params}
            elif [[ ${stranded} == "singl-end" ]]
            then
                java -jar /opt/Trimmomatic-0.39/trimmomatic-0.39.jar SE \
                    ${reads.findAll().join(' ')} \
                    -threads ${task.cpus} \
                    -trimlog trimlog_${sample}.log \
                    ${trim_params}
            fi
            """
        }
        //${reads.get(0)} ${reads.get(1)} \
        read_pairs_trimmed.subscribe {
            println "\nTrimmed FASTQ files generated for ${it[0]}:"
            it[1].each {
                item -> println "\t${item}" 
            }
            println " "
        }
        break
        // --------------------        

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
                --readFilesIn ${reads.findAll().join(' ')} \
                --runThreadN ${task.cpus} \
                --outSAMtype BAM Unsorted \
                --outFileNamePrefix ${sample}_
            """
        }
        star_alignments.subscribe {
            println "\nAlignment files generated for ${it[0]}:"
            it[1].each {
                item -> println "\t${item}" 
            }
            println " "
        }        
        break
        // ==========
        
        case['run.ReadCounting']:
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
        htseqCounts_cleaned.subscribe {
            println "\nRead Counts generated by HTSeqCounts for all sample:"
            println "\t${it}" 
            println " "
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
        featureCounts_cleaned.subscribe {
            println "\nRead Counts generated by featureCounts for all sample:"
            println "\t${it}" 
            println " "
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
