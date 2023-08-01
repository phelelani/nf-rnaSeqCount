#!/usr/bin/env bash

## This script is for tidying up the gene counts output from featureCounts and create a matrix ready for use.
while read line
do 
    if [ ${line:0:6} == "Geneid" ]
    then
        read -a arr <<< $line
        for i in "${!arr[@]}"
        do 
            fixed_name=$(basename  $(sed 's|_Aligned.sortedByCoord.out.bam||' <<< ${arr[$i]}))
            arr[$i]="$fixed_name"
        done
        
        echo -e ${arr[@]} | sed 's/ /\t/g'
    else
        echo -e "$line"
    fi 
done < <(sed '1d' $1 | cut -f 1,7-) > $2
