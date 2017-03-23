#!/bin/bash

file_1=$1

sample_name=${file_1:0:11}

while IFS=: read machine run_no flowcell flowcell_lane barcode
do
    echo -e 'RGID='$flowcell'.'$barcode'.'$flowcell_lane' RGLB='$flowcell'.'$flowcell_lane' RGPL=ILLUMINA RGPU='$flowcell'.'$barcode'.'$flowcell_lane' RGSM='$sample_name
done <<< "$(LC_ALL=C grep -m 10 --no-filename ^'@HWI' $file_1 | cut -d: -f 1-4,10 | sort -u)"

