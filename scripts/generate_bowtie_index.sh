#!/usr/bin/bash

input=$1
directory=$(dirname $input)
echo $directory
echo "$directory:$directory"

singularity exec --cleanenv --bind /global:/global --bind "$directory:$directory" ../containers/phelelani-nf-rnaSeqCount-master-trinity.simg \
    bowtie2-build $input $directory
    
