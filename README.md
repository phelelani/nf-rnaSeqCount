# nf-rnaSeqCount
*nf-rnaSeqCount* is a [Nextflow](http://nextflow.io/) pipeline for obtaining raw read counts for RNA-seq data using a given reference genome and annotation. This pipeline 

<p align="center">
  <img height="480" src="nf-rnaSeqCount.png">
</p>

# Pipeline Dependencies
To use the rnaSeqCount pipeline, the following dependencies are required:
## _*Softwares*_
- [x] [Nextflow](https://www.nextflow.io/)
- [x] [Singularity](http://singularity.lbl.gov/)

## _*Singularity Containers*_
- [x] [STAR](https://github.com/alexdobin/STAR) - ```shub://phelelani/nf-rnaSeqCount:star```
- [x] [HTSeq-Counts](https://htseq.readthedocs.io/en/release_0.9.1/overview.html) - ```shub://phelelani/nf-rnaSeqCount:htseqcount```
- [x] [featureCounts](http://subread.sourceforge.net/) - ```shub://phelelani/nf-rnaSeqCount:featurecounts```
- [x] [MultiQC](http://multiqc.info/) - ```shub://phelelani/nf-rnaSeqCount:multiqc```

## _*Reference Genome and Indexes*_
- [x] Reference Genome (.fa) and Genome Annotation (.gtf) files
- [x] Reference Genome Indexes (```bowtie2``` & ```STAR``` - see below on how to generate)

To generate the ```STAR``` and ```bowtie2``` indexes for the reference genome, run the following commands:
```
singularity exec --cleanenv containers/phelelani-rnaSeqCount-master-star.simg STAR --runThreadN 4 --runMode genomeGenerate --genomeDir <> --genomeFastaFiles <>
```

# Pipeline Execution

Edit main.nf:
```
params.data = "/path/to/data"
params.out = "/path/to/output"
params.genes = "/path/to/genes.gtf"
params.refSeq = "/path/to/genome.fa"
params.genome = "/path/to/STARIndex"
```

To run the pipeline:
```
nextflow run main.nf
```

# References
