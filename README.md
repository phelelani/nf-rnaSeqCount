# nf-rnaSeqCount

[![https://www.singularity-hub.org/static/img/hosted-singularity--hub-%23e32929.svg](https://www.singularity-hub.org/static/img/hosted-singularity--hub-%23e32929.svg)](https://singularity-hub.org/collections/770)

`nf-rnaSeqCount` is a [Nextflow](http://nextflow.io/) pipeline for obtaining raw read counts for RNA-seq data using a given reference genome and annotation. To use the `nf-rnaSeqCount` pipeline, the following dependencies are required:

   1. Softwares
      - [x] [Nextflow](https://www.nextflow.io/)
      - [x] [Singularity](http://singularity.lbl.gov/)

   2. Singularity Containers
      - [x] https://www.singularity-hub.org/collections/770

   3. Reference Genome and Indexes
      - Reference genome (`.fa`/`.fasta`) and genome annotation (`.gtf`) files.
      - Reference genome indexes (`bowtie2` & `STAR` - see *1.3.* below on how to generate the indexes).

<p align="center">
  <img width="1000" src="nf-rnaSeqCount.png">
</p>

## 1. Obtaining the `nf-rnaSeqCount` pipeline and preparing data
First, you need to clone the `nf-rnaSeqCount` repository onto you machine. You can eisther use `git` or `nextflow` (see the two methods below). I recommend using `nextflow` and creating you own `config` file (will explain later) for executing the workflow in the directory of your choosing.
```bash
## Using nextflow
nextflow pull https://github.com/h3abionet/h3avarcall.git

## Using git
git clone https://github.com/h3abionet/h3avarcall.git
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

### 1.1. Download test datasets 
NB: Skip this section if you have your own data to analyse using this workflow! 
SAY SOMETHING ABOUT THE PRACTICE DATA!
- [x] Download the mouse reference genome along with its annotation:
```
lftp -e 'pget -n10 ftp://ftp.ensembl.org/pub/release-68/fasta/mus_musculus/dna/Mus_musculus.GRCm38.68.dna.toplevel.fa.gz; bye'
lftp -e 'pget -n10 ftp://ftp.ensembl.org/pub/release-68/gtf/mus_musculus/Mus_musculus.GRCm38.68.gtf.gz; bye'
```

- [x] Download RNA-seq test dataset from H3ABioNet:
```
for sample in sample{37..42}_R{1,2}.fastq.gz; do wget -c http://h3data.cbio.uct.ac.za/assessments/RNASeq/practice/dataset/$sample; done
```

### 1.2. Download the `Singularity` containers (required to execute the pipeline):
```bash
nextflow run main.nf -profile slurm --mode prep.Containers
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

- [x] Generate ```STAR``` index
```
nextflow run main.nf -profile slurm --mode prep.STARIndex 
```

- [x] Generate ```bowtie2``` index
```
nextflow run main.nf -profile slurm --mode prep.BowtieIndex
```

## 2. Executing the main `nf-rnaSeqCount` pipeline

### 2.1. Read QC (optional): <br/>

To perform the QC of your fastq files, use this command:
```bash
nextflow run main.nf -profile slurm --mode run.ReadQC
```

### 2.2. Read Trimming (optional):<br/>

To run the trimming step of the `nf-rnaSeqCount` pipeline, use this command:
```bash
nextflow run main.nf -profile slurm --mode run.ReadTrimming
```

### 2.3. Read Alignment <br/>

To run the read alignment step of the `nf-rnaSeqCount` pipeline, use this comman (NB: can be run with `--from run.ReadTrimming` if you would like to use your trimmed reads):
```bash
nextflow run main.nf -profile slurm --mode run.ReadAlignment
```

### 2.4. Read Counting <br/>
This step uses the `BAM` file outputs generated by the read alignment step! You **MUST** run STEP 2.3 (`--mode run.ReadAlignment`) before running this step:
```bash
nextflow run main.nf -profile slurm --mode run.ReadCounting
```

### 2.6. Workflow QC (MultiQC - Optional)
This step performs a Quality Check of the different pipeline steps that have been ran. You need to run at least ONE step of the pipeline to be able to run this MultiQC step!
```bash
nextflow run main.nf -profile slurm --mode run.MultiQC 
```

## 3. Explore `nf-rnaSeqCount` results

```
- [1] Read QC (optional)         =>    `<results>/1_RQC`
- [2] Read Trimming (optional)   =>    `<results>/2_Read_Trimming`
- [3] Read Alignment             =>    `<results>/3_Read_Alignment`
- [4] Variant Calling            =>    `<results>/4_Read_Counts`
- [5] MultiQC                    =>    `<results>/5_MultiQC
- [6] Workflow tracing           =>    `<results>/workflow-tracing
```
In each of these folders, a sub-folder "`workflow_report`"  is created. It  contains 4 different files (`h3avarcall_report.html`, `h3avarcall_timeline.html`, `h3avarcall_workflow.dot` and `h3avarcall_trace.txt`) containing detailed information on the resources (CPU, MEMORY and TIME) usage of each process in the different pipeline steps. <br/> 
The `results` directory structure within `h3avarcall` repository can be summarized as below:

```bash
h3avarcall
  |--results
  |  |--1_QC
  |  |  |--workflow_report
  |  |  |  |--h3avarcall_report.html
  |  |  |  |--h3avarcall_timeline.html
  |  |  |  |--h3avarcall_workflow.dot
  |  |  |  |--h3avarcall_trace.txt
  |  |  |--<sample_1>_R1.fastqc.html .. <sample_N>_R1.fastqc.html
  |  |  |--<sample_1>_R2.fastqc.html .. <sample_N>_R1.fastqc.html
  |  |--2_Read_Trimming
  |  |  |--workflow_report
  |  |  |  |--h3avarcall_report.html
  |  |  |  |--h3avarcall_timeline.html
  |  |  |  |--h3avarcall_workflow.dot
  |  |  |  |--h3avarcall_trace.txt
  |  |  |--<sample_1>.1P.fastq.gz .. <sample_N>.1P.fastq.gz
  |  |  |--<sample_1>.2P.fastq.gz .. <sample_N>.2P.fastq.gz
  |  |--3_Read_Alignment
  |  |  |--workflow_report
  |  |  |  |--h3avarcall_report.html
  |  |  |  |--h3avarcall_timeline.html
  |  |  |  |--h3avarcall_workflow.dot
  |  |  |  |--h3avarcall_trace.txt
  |  |  |--<sample_1>_md.recal.bam .. <sample_N>_md.recal.bam
  |  |  |--<sample_1>_md.recal.bai .. <sample_N>_md.recal.bai
  |  |--4_Variant_Calling
  |  |  |--workflow_report
  |  |  |  |--h3avarcall_report.html
  |  |  |  |--h3avarcall_timeline.html
  |  |  |  |--h3avarcall_workflow.dot
  |  |  |  |--h3avarcall_trace.txt
  |  |  |--chr_1_genotyped.vcf.gz .. chr_22_genotyped.vcf.gz
  |  |  |--chr_1_genotyped.vcf.gz.tbi .. chr_22_genotyped.vcf.gz.tbi
  |  |--5_Variant_Filtering
  |  |  |--workflow_report
  |  |  |  |--h3avarcall_report.html
  |  |  |  |--h3avarcall_timeline.html
  |  |  |  |--h3avarcall_workflow.dot
  |  |  |  |--h3avarcall_trace.txt
  |  |  |--genome.SNP-recal.vcf.gz
  |  |  |--genome.SNP-recal.vcf.gz.tbi
  |  |--MultiQC
  |  |  |--multiqc_data
  |  |  |--multiqc_report.html
  |--work
  |  |--<There's a lot of folders here! Lets not worry about them for today!>
```
##We're working on further improving the pipleine and the associated documentation, feel free to share comments and suggestions!
