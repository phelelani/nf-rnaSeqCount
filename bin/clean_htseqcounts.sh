#!/usr/bin/env bash

header='Geneid'
counter=0

while read line
do
    sample=$(basename ${line:0:${#line}-4})
    header="$header\t$sample"
    
    if [[ $counter == 0 ]]
    then 
        file_a=$line
        (( counter ++ ))
    else
        file_b=$line
        join -t $'\t' $file_a $file_b > tmp
        mv tmp the_matrix
        file_a=the_matrix
    fi
done < $1

sed -i '1s/^/'"$header\n"'/' the_matrix
mv the_matrix $2
