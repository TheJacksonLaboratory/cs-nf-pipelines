#!/usr/bin/env bash

name_string=${1}
merged_vcf=${2}

perl -ne 'print "$1\n" if /SUPP_VEC=([^,;]+)/' ${merged_vcf} | \
    sed -e 's/\(.\)/\1 /g' > ${name_string}.merged.overlap.txt

grep -v "#" ${merged_vcf} | \
    cut -f 1,2,3 > ${name_string}.merged.SVs.txt

printf "chr\tpos\tSV\tpbsv\tsniffels\n" > ${name_string}.merged.overlap.annotated.txt && \
    paste ${name_string}.merged.SVs.txt ${name_string}.merged.overlap.txt | \
    sed -e 's/ /\t/g' >> ${name_string}.merged.overlap.annotated.txt