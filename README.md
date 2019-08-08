# nf-rnaSeqCount
```nf-rnaSeqCount``` is a [Nextflow](http://nextflow.io/) pipeline for obtaining raw read counts for RNA-seq data using a given reference genome and annotation. The 

<p align="center">
  <img width="1000" src="nf-rnaSeqCount.png">
</p>

# 1. Pipeline Dependencies
To use the rnaSeqCount pipeline, the following dependencies are required:
### The data
- [x] Download the mouse reference genome along with its annotation:
```
lftp -e 'pget -n10 ftp://ftp.ensembl.org/pub/release-68/fasta/mus_musculus/dna/Mus_musculus.GRCm38.68.dna.toplevel.fa.gz; bye'
lftp -e 'pget -n10 ftp://ftp.ensembl.org/pub/release-68/gtf/mus_musculus/Mus_musculus.GRCm38.68.gtf.gz; bye'
```

- [x] Download RNA-seq test dataset from H3ABioNet:
```
for sample in sample{37..42}_R{1,2}.fastq.gz; do wget -c http://h3data.cbio.uct.ac.za/assessments/RNASeq/practice/dataset/$sample; done
```

### 1.1. Softwares
- [x] [Nextflow](https://www.nextflow.io/)
- [x] [Singularity](http://singularity.lbl.gov/)
- [x] [R](https://www.r-project.org/)

### 1.2.  Singularity Containers
- [x] https://www.singularity-hub.org/collections/770

### 1.3. Reference Genome and Indexes
- [x] Reference Genome (.fa) and Genome Annotation (.gtf) files
- [x] Reference Genome Indexes (```bowtie2``` & ```STAR``` - see *3.* below on how to generate)

# 2. Optaining the ```nf-rnaSeqCount``` pipeline
The ```nf-rnaSeqCount``` pipeline can be obtain using any of the following methods:

### 2.1. Using the ```git``` command (recommended):
- [x] ```git clone https://github.com/phelelani/nf-rnaSeqCount.git```

### 2.2. Using the ```nextflow``` command:
- [x] ```nextflow pull phelelani/nf-rnaSeqCount```
- [x] ```nextflow pull https://github.com/phelelani/nf-rnaSeqCount.git```
- [x] ```nextflow clone phelelani/nf-rnaSeqCount <target-dir>```

# 3. Generating genome indexes.
To generate the ```STAR``` and ```bowtie2``` indexes for the reference genome, the ```Singularity``` containers first need to be downloaded from [Singularity Hub](ttps://www.singularity-hub.org). The ```prepareData.nf``` script can be used to download and prepare data (generate indexes) to be be used with the ```nf-rnaSeqCount``` pipeline. The ```prepareData.nf``` can be run in three different modes:
- [x] ```getContainers```: for downloading the required ```Singularity``` containers.
- [x] ```generateStarIndex```: for generating ```STAR``` indexes.
- [x] ```generateBowtieIndex```: for generating ```bowtie2``` indexes.

To generate the genome indexes, run the following commands in the pipeline directory:

### 3.1 Download ```Singularity``` containers:
```
nextflow run prepareData.nf --mode getContainers -profile pbsPrepare
```

### 3.2. Generate ```STAR``` index
```
nextflow run prepareData.nf --mode generateStarIndex -profile pbsPrepare
```

### 3.3. Generate ```bowtie2``` index
```
nextflow run prepareData.nf --mode generateBowtieIndex -profile pbsPrepare
```

NB: The ```-profile``` option can either be one of depending on the scheduler.

# 4. Pipeline Execution
The ```nf-rnaSeqCount``` pipeline can be run in one of two ways:

### 4.1. By editing the ```parameters.config``` file and specifying the parameters (recommended)
Edit ```main.config```:
```
/*
 *  USE THIS FILE TO SPECIFY YOUR PARAMETERS. ALLOWED PARAMETERS ARE AS FOLLOWS:
 *  ============================================================================
 *  data      : Path to where the input data is located (where fastq files are located).
 *  filetype  : Extension of the input FASTQ files (fastq | fq | fastq.gz | fq.gz | fastq.bz2 | fq.bz2).
 *  out       : Path to where the output should be directed (will be created if it does not exist).
 *  genome    : The whole genome sequence (fasta | fa | fna).
 *  index     : Path to where the STAR index files are locaded.
 *  genes     : The genome annotation file (gtf).
 *  bind      : Paths to be passed onto the singularity image (Semi-colon separated).
 *  help      : Print out help menu. Passed as "--help" to the "main.nf" script for detailed information
 */
params {
    data      = "/path/to/data"
    out       = "/path/to/output"
    filetype  = "fastq.gz"
    genome    = "/path/to/genome.fa"  
    index     = "/path/to/STARIndex"
    genes     = "/path/to/genes.gtf" 
    bind      = "/path/to/bind_1;/path/to/bind_2"
    help      = null
}
```

Then run the pipeline:
```
nextflow run main.nf
```

### 4.2. Directly from the command line by supplying the required parameters
```
nextflow run main.nf \
    --data "/path/to/data" \
    --filetype "filetype"
    --out "/path/to/output" \
    --genome "/path/to/genome.fa" \
    --index "/path/to/STARIndex" \
    --genes "/path/to/genes.gtf" \
    --bind "/path/to/bind_1;/another/path/to/bind_2"
```

# References
