# nf-rnaSeqCount

[![https://www.singularity-hub.org/static/img/hosted-singularity--hub-%23e32929.svg](https://www.singularity-hub.org/static/img/hosted-singularity--hub-%23e32929.svg)](https://singularity-hub.org/collections/770)

`nf-rnaSeqCount` is a [`Nextflow`](http://nextflow.io/) pipeline for obtaining raw read counts for RNA-seq data using a given reference genome and annotation. To use the `nf-rnaSeqCount` pipeline, the following dependencies are required:
   1. Installed softwares:
      - [`Nextflow`](https://www.nextflow.io/)
      - [`Singularity`](http://singularity.lbl.gov/)
   2. `Singularity` [containers](https://www.singularity-hub.org/collections/770) with the required applications/programs for executing the workflow:
      - `nf-rnaSeqCount-fastqc.sif`
      - `nf-rnaSeqCount-featurecounts.sif`
      - `nf-rnaSeqCount-htseqcount.sif`
      - `nf-rnaSeqCount-multiqc.sif`
      - `nf-rnaSeqCount-star.sif`
      - `nf-rnaSeqCount-trimmomatic.sif`
      - `nf-rnaSeqCount-trinity.sif`
   3. Reference genome, annotation and indexes
      - Reference genome (`.fa`/`.fasta`) and genome annotation (`.gtf`) files.
      - Reference genome indexes (`bowtie2` & `STAR` - see *1.3.* below on how to generate the indexes).

---

<p align="center">
  <img width="1000" src="nf-rnaSeqCount.png">
</p>

---

## 1. Obtaining the `nf-rnaSeqCount` pipeline and preparing data
First, you need to clone the `nf-rnaSeqCount` repository onto you machine. You can eisther use `git` or `nextflow` (see the two methods below). I recommend using `nextflow` and creating you own `config` file (will explain later) for executing the workflow in the directory of your choosing. The rest of this documentation assumes that you have used `nextflow` to clone this workflow - If your're an expert and have used `git` to clone the workflow - you know what to do :)
```bash
## Using nextflow
nextflow pull https://github.com/phelelani/nf-rnaSeqCount.git
```
Content of the repository:
```bash
nf-rnaSeqCount
  |--containers                       ## Folder for Singularity images and recipes (in case you want to build yourself). All downloaded images go here!
  |  |--Singularity.fastqc            ## Singularity recipe file for 
  |  |--Singularity.featureCounts     ## Singularity recipe file for 
  |  |--Singularity.htseqCount        ## Singularity recipe file for 
  |  |--Singularity.multiQC           ## Singularity recipe file for 
  |  |--Singularity.star              ## Singularity recipe file for 
  |  |--Singularity.trimmomatic       ## Singularity recipe file for 
  |  |--Singularity.trinity           ## Singularity recipe file for 
  |--templates                        ## Folder for extra scripts for the pipeline.
  |  |--clean_featureCounts.sh        ## Script for 
  |  |--clean_htseqCounts.sh          ## Script for 
  |--LICENSE                          ## Duh!
  |--main.config                      ## User configuration file! All inputs, outputs and options GO HERE!! ONLY file that SHOULD be modified by user!
  |--main.nf                          ## Main h3avarcall nextflow scripts.
  |--nextflow.config                  ## Pipeline configuration file! DO NOT EDIT!!!
  |--nextfnf-rnaSeqCount.png          ## Pipeline 
  |--README.md                        ## Duh!
```

### 1.1. Download test datasets (optional)
*__NB__: Skip this section if you have your own data to analyse using this workflow! This section is only for getting data to practice using the `nf-rnaSeqCount` workflow!* 

- [x] Download and decompress the mouse reference genome along with its annotation:
```
wget -c -O genome.fa.gz ftp://ftp.ensembl.org/pub/release-68/fasta/mus_musculus/dna/Mus_musculus.GRCm38.68.dna.toplevel.fa.gz
wget -c -O genes.gtd.gz ftp://ftp.ensembl.org/pub/release-68/gtf/mus_musculus/Mus_musculus.GRCm38.68.gtf.gz
gunzip genome.fa.gz
gunzip genes.gtf.gz
```

- [x] Download RNA-seq test dataset from H3ABioNet:
```
lftp -e 'mirror -c --use-pget-n=5 http://h3data.cbio.uct.ac.za/assessments/RNASeq/practice/dataset; bye'
for sample in sample{37..42}_R{1,2}.fastq.gz; do wget -c http://h3data.cbio.uct.ac.za/assessments/RNASeq/practice/dataset/$sample; done
```

### 1.2. Download the `Singularity` containers (required to execute the pipeline):
```bash
nextflow run nf-rnaSeqCount -profile slurm --mode prep.Containers
```

### 1.3. Generating genome indexes.
The only file that you should edit is the `main.config` file. I recommend creating you own `config` file using the template below:
```groovy
/*==================================================================================================
 * THIS FILE IS USED TO SPECIFY INPUT, OUTPUTS AND PARAMETERS. THE FOLLOWING OPTIONS ARE THE ALLOWED:
 * ==================================================================================================
 * data         : Path to where the data is (FASTQ files).
 * out          : Path to store output results.
 * bundle       : GATK-b37-bundle list file.
 * mode         : Worflow step to perform. Can be any of [ do.GetContainers | do.GenomeIndexing | do.QC | do.ReadTrimming | do.ReadAlignment | do.VarianCalling | do.VariantFiltering | do.MultiQC].
 * trim         : Trimming options for Trimmomatic.
 * resources    : Location of the GATK-b37-bundle folder.
 * from         : pipeline step to resume pipeline from. Can be any of [ do.QC | do.ReadTrimming | do.ReadAlignment | do.VarianCalling | do.VariantFiltering ].
 * params.help  : Print help menu.
 * ==================================================================================================
 * BELOW ARE THE DEFAULT PARAMETERS! YOU'RE MORE THAN WELCOME TO CHANGE AS DESIRED!
 * ==================================================================================================
 */

params {
    data         = "$baseDir/data"
    out          = "$baseDir/results"
    mode         = "do.QC"
    trim         = "ILLUMINACLIP:TruSeq3-PE-2.fa:2:30:10:8:true TRAILING:28 MINLEN:40"
    from         = null
}
```

- [x] `prep.STARIndex`   : for generating `STAR` indexes.
- [x] `prep.BowtieIndex` : for generating `bowtie2` indexes.

To generate the genome indexes, run the following commands:

- [x] Generate ```STAR``` indexes
```
nextflow run nf-rnaSeqCount -profile slurm --mode prep.STARIndex 
```

- [x] Generate ```bowtie2``` indexes
```
nextflow run nf-rnaSeqCount -profile slurm --mode prep.BowtieIndex
```

---

## 2. Executing the main `nf-rnaSeqCount` pipeline

### 2.1. Read QC (optional):

To perform the QC of your fastq files, use this command:
```bash
nextflow run nf-rnaSeqCount -profile slurm --mode run.ReadQC
```

### 2.2. Read Trimming (optional):

To run the trimming step of the `nf-rnaSeqCount` pipeline, use this command:
```bash
nextflow run nf-rnaSeqCount -profile slurm --mode run.ReadTrimming
```

### 2.3. Read Alignment:

To run the read alignment step of the `nf-rnaSeqCount` pipeline, use this comman (NB: can be run with `--from run.ReadTrimming` if you would like to use your trimmed reads):
```bash
nextflow run nf-rnaSeqCount -profile slurm --mode run.ReadAlignment
```

### 2.4. Read Counting:
This step uses the `BAM` file outputs generated by the read alignment step! You **MUST** run STEP 2.3 (`--mode run.ReadAlignment`) before running this step:
```bash
nextflow run nf-rnaSeqCount -profile slurm --mode run.ReadCounting
```

### 2.6. Workflow QC (optional):
This step performs a Quality Check of the different pipeline steps that have been ran. You need to run at least ONE step of the pipeline to be able to run this MultiQC step!
```bash
nextflow run nf-rnaSeqCount -profile slurm --mode run.MultiQC 
```

---


## 3. Explore `nf-rnaSeqCount` results

```
- [1] Read QC (optional)         =>    `<output_directory>/1_RQC`
- [2] Read Trimming (optional)   =>    `<output_directory>/2_Read_Trimming`
- [3] Read Alignment             =>    `<output_directory>/3_Read_Alignment`
- [4] Read Counting              =>    `<output_directory>/4_Read_Counts`
- [5] MultiQC                    =>    `<output_directory>/5_MultiQC
- [6] Workflow tracing           =>    `<output_directory>/workflow-tracing
```
In addition to the 5 directories created for each step in the results directory, a directory `workflow-tracing` is created to monitor the resources used in each step. This directory will contain 4 files for each step (--mode) of the workflow:
- `nf-rnaSeqCount_<mode>_report.html`
- `nf-rnaSeqCount_<mode>_timeline.html`
- `nf-rnaSeqCount_<mode>_trace.txt`

These files contain detailed information on the resources (CPU, MEMORY and TIME) usage of each of the process in the different pipeline steps. The `<output_directory>` directory structure is summarized below:

```bash
<output_directory>
  |--1_Read_QC
  |  |--<sample_1>_R1.fastqc.html .. <sample_N>_R1.fastqc.html
  |  |--<sample_1>_R2.fastqc.html .. <sample_N>_R1.fastqc.html
  |--2_Read_Trimming
  |  |--<sample_1>.1P.fastq.gz .. <sample_N>.1P.fastq.gz
  |  |--<sample_1>.2P.fastq.gz .. <sample_N>.2P.fastq.gz
  |--3_Read_Alignment
  |  |--<sample_1>_Aligned.out.bam .. <sample_N>_Aligned.out.bam
  |  |--<sample_1>_Log.final.out .. <sample_N>_Log.final.out
  |  |--<sample_1>_Log.out .. <sample_N>_Log.out
  |  |--<sample_1>_Log.progress.out .. <sample_N>_Log.progress.out
  |  |--<samplle_1>_SJ.out.tab .. <sample>_SJ.out.tab
  |--4_Read_Counts
  |  |--featureCounts
  |  |  |--gene_counts_final.txt
  |  |  |--gene_counts.txt
  |  |  |--gene_counts.txt.jcounts
  |  |  |--gene_counts.txt.summary
  |  |--htseqCounts
  |  |  |--gene_counts_final.txt
  |  |  |--<sample>.txt .. <sample>.txt
  |--5_MultiQC
  |  |--multiqc_data
  |  |--multiqc_report.html
  |--workflow-tracing
  |  |--nf-rnaSeqCount_run.MultiQC_{report.html,timeline.html,trace.txt}
  |  |--nf-rnaSeqCount_run.ReadAlignment_{report.html,timeline.html,trace.txt}
  |  |--nf-rnaSeqCount_run.ReadCounting_{report.html,timeline.html,trace.txt}
  |  |--nf-rnaSeqCount_run.ReadTrimming_{report.html,timeline.html,trace.txt}   
  |  |--nf-rnaSeqCount_run.ReadQC_{report.html,timeline.html,trace.txt}
```
##We're working on further improving the pipleine and the associated documentation, feel free to share comments and suggestions!

---
