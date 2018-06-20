#!/usr/bin/bash

input=$1
directory=$(dirname $input)
echo $directory
echo "$directory:$directory"

singularity exec --cleanenv --bind /global:/global --bind "$directory:$directory" ../containers/phelelani-nf-rnaSeqCount-master-star.simg \
    STAR --runThreadN 4 \
    --runMode genomeGenerate \
    --genomeDir $directory \
    --genomeFastaFiles $input
    
