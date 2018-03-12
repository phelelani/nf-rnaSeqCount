# rnaSeqCount - " something catchy goes in here :) "
Some fancy text
----

# Pipeline Dependencies
## _*Softwares*_
- [x] [Nextflow](https://www.nextflow.io/)
- [x] [Singularity](http://singularity.lbl.gov/)

## _*Singularity Containers*_
- [x] [STAR](https://github.com/alexdobin/STAR) - shub://phelelani/rnaSeqCount:star
- [x] [HTSeq-Counts](https://htseq.readthedocs.io/en/release_0.9.1/overview.html) - shub://phelelani/rnaSeqCount:htseqcount
- [x] [featureCounts](http://subread.sourceforge.net/) - shub://phelelani/rnaSeqCount:featurecounts
- [x] [MultiQC](http://multiqc.info/) - shub://phelelani/rnaSeqCount:multiqc

## _*Reference Genome and Indexes*_
- [x] Reference Genome & Annotation
- [x] STAR Index

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
