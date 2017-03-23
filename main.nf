#!/usr/bin/env nextflow

params.data = "/spaces/phelelani/ssc_data/data_trimmed/inflated" 
params.out = "/spaces/phelelani/ssc_data/assembly/assembly_results"
params.genes = "/global/blast/reference_genomes/Homo_sapiens/UCSC/hg19/Annotation/Genes/genes.gtf"
params.genome = "/global/blast/reference_genomes/Homo_sapiens/UCSC/hg19/Sequence/Bowtie2Index/genome"
params.ref = "/global/blast/reference_genomes/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa"
params.rg_script="$PWD/get_rg_info.sh"

data_path = params.data
out_path = file(params.out)
genes = params.genes
genome = params.genome
ref = params.ref
rg_script = params.rg_script

out_path.mkdir()

//read_pair = Channel.fromFilePairs("${data_path}/caskiSubset_0*R{1,2}.fq", type: 'file')
read_pair = Channel.fromFilePairs("${data_path}/*R[1,2].fastq", type: 'file')
//read_pair = Channel.fromFilePairs("${data_path}/CASE_03*R[1,2].fastq", type: 'file')

// 1. Align reads to reference genome
process runTophat_process {
    cache = true
    executor 'pbs'
    queue 'WitsLong'
    cpus 9
    memory '100 GB'
    time '100h'
    scratch '$HOME/tmp'
    tag { sample }
    stageInMode 'symlink'
    stageOutMode 'rsync'
    publishDir "$out_path/${sample}", mode: 'copy', overwrite: false

    input:
    set sample, file(reads) from read_pair

    output:
    set sample, "tophat_out_${sample}/*" into tophat_results
    set sample, file("tophat_out_${sample}/accepted_hits.bam"), file(reads) into tophat_raw_bams

    """
    tophat -p 8 -G $genes -o tophat_out_${sample} --library-type=fr-firststrand $genome ${reads.get(0)} ${reads.get(1)}
    """
}


// 2. Reorder to match reference genome
process runReorderSam_process {
    cache = true
    executor 'pbs'
    queue 'WitsLong'
    cpus 9
    memory '100 GB'
    time '100h'
    scratch '$HOME/tmp'
    tag { sample }
    stageInMode 'symlink'
    stageOutMode 'rsync'
    publishDir "$out_path/${sample}", mode: 'copy', overwrite: false

    input:
    set sample, file(bam), file(reads) from tophat_raw_bams

    output:
    set sample, file("tophat_out_${sample}/${sample}_accepted_hits.bam"), file(reads) into reordered_bams

    """
    java -jar /opt/exp_soft/bioinf/picard-tools/ReorderSam.jar \
        I=${bam} \
        O=tophat_out_${sample}/${sample}_accepted_hits.bam \
        R=$ref
    """
}

// 3. Sort reads
process runSamtoolsSort_process {
    cache = true
    executor 'pbs'
    queue 'WitsLong'
    cpus 9
    memory '100 GB'
    time '100h'
    scratch '$HOME/tmp'
    tag { sample }
    stageInMode 'symlink'
    stageOutMode 'rsync'
    publishDir "$out_path/${sample}", mode: 'copy', overwrite: false

    input:
    set sample, file(bam), file(reads) from reordered_bams

    output:
    set sample, file("tophat_out_${sample}/${sample}_accepted_hits_sorted.bam"), file(reads) into reordered_sorted_bams

    """
    samtools sort -O bam \
        -T tophat_out_${sample}/${sample}_accepted_hits_sorted \
        -o tophat_out_${sample}/${sample}_accepted_hits_sorted.bam \
        -@ 8 ${bam}
    """
}

// 4. Add read group information
process runAddReadGroups_process {
    cache = true
    executor 'pbs'
    queue 'WitsLong'
    cpus 9
    memory '100 GB'
    time '100h'
    scratch '$HOME/tmp'
    tag { sample }
    stageInMode 'symlink'
    stageOutMode 'rsync'
    publishDir "$out_path/${sample}", mode: 'copy', overwrite: false

    input:
    set sample, file(bam), file(reads) from reordered_sorted_bams
    
    output:
    set sample, file("tophat_out_${sample}/${sample}_accepted_hits_final.bam"), file("tophat_out_${sample}/${sample}_accepted_hits_final.bam.bai") into bams_cufflinks, bams_cuffquant
    
    """
    java -jar /opt/exp_soft/bioinf/picard-tools/AddOrReplaceReadGroups.jar \
        I=${bam} \
        O=tophat_out_${sample}/${sample}_accepted_hits_final.bam `${rg_script} ${reads.get(0)}`
    samtools index tophat_out_${sample}/${sample}_accepted_hits_final.bam
    """
}

// 5. Assemble the reads.
process runCufflinks_process {
    cache = true
    executor 'pbs'
    queue 'WitsLong'
    cpus 9
    memory '100 GB'
    time '100h'
    scratch '$HOME/tmp'
    tag { sample }
    stageInMode 'symlink'
    stageOutMode 'rsync'
    publishDir "$out_path/${sample}", mode: 'copy', overwrite: true
    
    input:
    set sample, file(bam_file), file(index) from bams_cufflinks

    output:
    set sample, file("cufflinks_out_${sample}/*") into cufflinks_results
    file("cufflinks_out_${sample}/transcripts.gtf") into gtf_files

    """
    cufflinks -p 8 -o cufflinks_out_${sample} ${bam_file} --library-type=fr-firststrand
    """
}

// 6. Collect all the "transcripts.gtf" locations and put them into a single file that is needed by cuffmerge
gtf_files
.collectFile () { item -> [ 'assemblies.txt', "${item}" + '\n' ] }
.set { assembly_names }

// 7. Merge all the assemblies
process runCuffmerge_process {
    cache = true
    echo = true
    executor 'pbs'
    queue 'WitsLong'
    cpus 9
    memory '50 GB'
    time '100h'
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
    cuffmerge -g $genes -s $ref -o merged_assemblies -p 8 ${gtf}
    """
}

// 8. 
process runCuffquant_process {
    cache = true
    executor 'pbs'
    queue 'WitsLong'
    cpus 9
    memory '100 GB'
    time '100h'
    scratch '$HOME/tmp'
    tag { sample }
    stageInMode 'symlink'
    stageOutMode 'rsync'
    publishDir "$out_path/${sample}", mode: 'copy', overwrite: true
    
    input:
    set sample, file(bam_file), file(index) from bams_cuffquant
    file(gtf) from merged_assembly.first()

    output:
    set sample, file("cuffquant_out_${sample}/*") into cuffquant_results
    file("cuffquant_out_${sample}/abundances.cxb") into cuffquant_cbx

    """
    cuffquant -p 8 -o cuffquant_out_${sample} -b ${ref} -u ${gtf} --library-type=fr-firststrand ${bam_file}
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

/*    """
    cd $out_path
    find . -maxdepth 3 -iname accepted_hits.bam -exec bash -c 'cp \$1 `sed 's%accepted_hits.bam%\${1:2:15}_accepted_hits.bam%' <<< \$1`' - {} \\;
    cuffquant -p 7 -o cuffquant_out_${sample} -b $genome".fa" -u ${gtf} --library-type=fr-firststrand `find . -maxdepth 3 -iname *_accepted_hits.bam -printf "%p," | sed '//s/,\$//'`
    """
    """
    cufflinks -p 10 -o cufflinks_out_${sample} ${items.find {item -> item.endsWith('accepted_hits.bam')} } --library-type=fr-firststrand
    """
*/
