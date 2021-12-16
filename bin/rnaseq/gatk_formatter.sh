#!/bin/bash

# A script to format gatk coverage file

export file=$1
export output1=$2
export output2=$3
export target=$4


cat  $file|awk ' { print $1,$2 } ' |tr ':' '\t'|awk ' { print $1,$2,$2+1,$3}' | tr ' ' '\t'|sed '/^Locus/d'  > $output1
cat $target|awk ' { print $1,$2,$3,$4}' |tr ' ' '\t' > Target.bed

bedtools intersect -loj -a Target.bed  -b $output1|awk ' { print $1,$2,$3,$4,0,$NF } ' |tr ' ' '\t' |sed 's/\./0/g' > $output2
rm Target.bed