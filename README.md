# nf-rnaSeqCount
*nf-rnaSeqCount* is a [Nextflow](http://nextflow.io/) pipeline for obtaining raw read counts for RNA-seq data using a given reference genome and annotation. This pipeline 

<p align="center">
  <img width="700" src="nf-rnaSeqCount.png">
</p>

# 1. Pipeline Dependencies
To use the rnaSeqCount pipeline, the following dependencies are required:
### 1.1. Softwares
- [x] [Nextflow](https://www.nextflow.io/)
- [x] [Singularity](http://singularity.lbl.gov/)
- [x] [R](https://www.r-project.org/)

### 1.2.  Singularity Containers
- [x] [STAR](https://github.com/alexdobin/STAR) - ```shub://phelelani/nf-rnaSeqCount:star```
- [x] [HTSeq-Counts](https://htseq.readthedocs.io/en/release_0.9.1/overview.html) - ```shub://phelelani/nf-rnaSeqCount:htseqcount```
- [x] [featureCounts](http://subread.sourceforge.net/) - ```shub://phelelani/nf-rnaSeqCount:featurecounts```
- [x] [MultiQC](http://multiqc.info/) - ```shub://phelelani/nf-rnaSeqCount:multiqc```

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
To generate the ```STAR``` and ```bowtie2``` indexes for the reference genome, run the following commands inside the downloaded nf-rnaSeqCount repository:
### 3.1. ```STAR``` index
```
sh scripts/generate_star_index.sh "/path/to/genome.fa"
```

### 3.2. ```bowtie2``` index
```
sh scripts/generate_bowtie_index.sh "/path/to/genome.fa"
```

# 4. Pipeline Execution
The ```nf-rnaSeqCount``` pipeline can be run in one of two ways:

### 4.1. By editing the ```parameters.config``` file and specifying the parameters (recommended)
Edit parameters.config:
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
    filetype  = "fastq.gz"
    out       = "/path/to/output"
    genome    = "/path/to/genome.fa'"  
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
