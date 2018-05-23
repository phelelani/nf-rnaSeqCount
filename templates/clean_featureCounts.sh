#!/bin/bash

while read line
do 
    if [ ${line:0:1} == "#" ]
    then
        :
    elif [ ${line:0:6} == "Geneid" ]
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
done < tmp_genes.txt > gene_counts_final.txt

