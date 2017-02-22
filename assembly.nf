#!/usr/bin/env nextflow

params.data = "/spaces/phelelani/ssc_data/data_trimmed/inflated" 
params.out = "/home/phelelani/scleroderma_analysis/jobs/assembly/nextflow_pipeline/the_results"
params.genes = "/global/blast/reference_genomes/Homo_sapiens/UCSC/hg19/Annotation/Genes/genes.gtf"
params.genome = "/global/blast/reference_genomes/Homo_sapiens/UCSC/hg19/Sequence/Bowtie2Index/genome"

data_path = params.data
out_path = file(params.out)
genes = params.genes
genome = params.genome

out_path.mkdir()

read_pair = Channel.fromFilePairs("${data_path}/caskiSubset_0*R{1,2}.fq", type: 'file')

process runTophat {
    cache = true
    //echo = true
    executor 'pbs'
    queue 'WitsLong'
    cpus 11
    memory '200 GB'
    time '20h'
    scratch true
    tag { sample }
    publishDir "$out_path/${sample}", mode: 'symlink', overwrite: true

    input:
    set sample, file(reads) from read_pair

    output:
    set sample, "tophat_out_${sample}/*" into tophat_results 

    """
    tophat -p 10 -G $genes -o tophat_out_${sample} --library-type=fr-firststrand $genome ${reads.get(0)} ${reads.get(1)}
    """
}


process runCufflinks {
    cache = true
    //echo = true
    executor 'pbs'
    queue 'WitsLong'
    cpus 11
    memory '50 GB'
    time '20h'
    scratch true
    tag { sample }
    publishDir "$out_path/${sample}", mode: 'symlink', overwrite: true
    
    input:
    set sample, items from tophat_results

    output:
    set sample, file("cuff_out_${sample}/*") into cufflinks_results
    file("cuff_out_${sample}/transcripts.gtf") into gtf_files

    """
    cufflinks -p 10 -o cuff_out_${sample} ${items.find {item -> item.endsWith('accepted_hits.bam')} }
    """
}

gtf_files
.collectFile () { item -> [ 'assemblies.txt', "${item}" + '\n' ] }
.set { assembly_names }

process runMerge {
    cache = true
    echo = true
    executor 'pbs'
    queue 'WitsLong'
    cpus 11
    memory '50 GB'
    time '20h'
    scratch true
    tag { 'merge' }
    publishDir "$out_path", mode: 'symlink', overwrite: true
    
    input:
    file(gtf) from assembly_names
    
    output:
    file("merged_assemblies/*") into merged_results
    
    """
    cuffmerge -g $genes -s $genome".fa" -o merged_assemblies -p 10 ${gtf}
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




/*
all_gtfs.subscribe{ println "${it.text}" }
merged_results.subscribe{ println it }
*/

//.flatMap()
//.flatten()
//.collectFile(name: 'assemblies.txt', newLine: true)
//gtf_files.subscribe { println it }

//allChannels = Channel.create()
//allChannels.bind(cufflink_results)
//.subscribe { 
//     println it.text
// }
//.subscribe { println it }
//.subscribe {
//    println it //"${it.name} contains:"
    //    println it.text
// }

//all_gtfs.subscribe { println it }