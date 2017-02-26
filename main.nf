#!/usr/bin/env nextflow

params.data = "/spaces/phelelani/ssc_data/data_trimmed/inflated" 
params.out = "/spaces/phelelani/ssc_data/assembly/assembly_results"
params.genes = "/global/blast/reference_genomes/Homo_sapiens/UCSC/hg19/Annotation/Genes/genes.gtf"
params.genome = "/global/blast/reference_genomes/Homo_sapiens/UCSC/hg19/Sequence/Bowtie2Index/genome"

data_path = params.data
out_path = file(params.out)
genes = params.genes
genome = params.genome

out_path.mkdir()

//read_pair = Channel.fromFilePairs("${data_path}/caskiSubset_0*R{1,2}.fq", type: 'file')
read_pair = Channel.fromFilePairs("${data_path}/*R{1,2}.fastq", type: 'file')
//read_pair.subscribe { println it }

// 1. Align reads to reference genome
process runTophat_process {
    cache = true
    executor 'pbs'
    queue 'WitsLong'
    cpus 11
    memory '200 GB'
    time '20h'
    scratch '$HOME/tmp'
    tag { sample }
    stageInMode 'symlink'
    stageOutMode 'rsync'
    publishDir "$out_path/${sample}", mode: 'copy', overwrite: false

    input:
    set sample, file(reads) from read_pair

    output:
    set sample, "tophat_out_${sample}/*" into tophat_results_cufflinks, tophat_results_cuffquant

    """
    tophat -p 10 -G $genes -o tophat_out_${sample} --library-type=fr-firststrand $genome ${reads.get(0)} ${reads.get(1)}
    """
}

// 2. Assemble the reads.
process runCufflinks_process {
    cache = true
    executor 'pbs'
    queue 'WitsLong'
    cpus 11
    memory '50 GB'
    time '20h'
    scratch '$HOME/tmp'
    tag { sample }
    stageInMode 'symlink'
    stageOutMode 'rsync'
    publishDir "$out_path/${sample}", mode: 'copy', overwrite: true
    
    input:
    set sample, items from tophat_results_cufflinks

    output:
    set sample, file("cufflinks_out_${sample}/*") into cufflinks_results
    file("cufflinks_out_${sample}/transcripts.gtf") into gtf_files

    """
    cufflinks -p 10 -o cufflinks_out_${sample} ${items.find {item -> item.endsWith('accepted_hits.bam')} } --library-type=fr-firststrand
    """
}

// 3. Collect all the "transcripts.gtf" locations and put them into a single file that is needed by cuffmerge
gtf_files
.collectFile () { item -> [ 'assemblies.txt', "${item}" + '\n' ] }
.set { assembly_names }

// 4. Merge all the assemblies
process runCuffmerge_process {
    cache = true
    echo = true
    executor 'pbs'
    queue 'WitsLong'
    cpus 11
    memory '50 GB'
    time '20h'
    scratch '$HOME/tmp'
    tag { 'merge' }
    stageInMode 'symlink'
    stageOutMode 'rsync'
    publishDir "$out_path", mode: 'copy', overwrite: true
    
    input:
    file(gtf) from assembly_names
    
    output:
    file("merged_assemblies/*") into merged_results
    file("merged_assemblies/merged.gtf") into merged_assembly
    
    """
    cuffmerge -g $genes -s $genome".fa" -o merged_assemblies -p 10 ${gtf}
    """
}

// 5. 
process runCuffquant_process {
    cache = true
    executor 'pbs'
    queue 'WitsLong'
    cpus 11
    memory '50 GB'
    time '20h'
    scratch '$HOME/tmp'
    tag { sample }
    stageInMode 'symlink'
    stageOutMode 'rsync'
    publishDir "$out_path/${sample}", mode: 'copy', overwrite: true
    
    input:
    set sample, items from tophat_results_cuffquant
    file(gtf) from merged_assembly

    output:
    set sample, file("cuffquant_out_${sample}/*") into cuffquant_results

    """
    cuffquant -p 10 -o cuffquant_out_${sample} -b $genome".fa" -u ${gtf} --library-type=fr-firststrand ${items.find {item -> item.endsWith('accepted_hits.bam')} }
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

