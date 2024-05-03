#!/usr/env/bin bash

sampleID=${1}

cat ${sampleID}.sniffles_calls.vcf | awk '$1 ~ /^#/ {print $0;next} {print $0 | "sort -k1,1 -k2,2n"}' > ${sampleID}.sniffles_calls_sorted.vcf

grep "#" ${sampleID}.sniffles_calls_sorted.vcf > ${sampleID}.sniffles_sorted_prefix.vcf

grep -v "#" ${sampleID}.sniffles_calls_sorted.vcf | \
    sed -e 's/None/0/g' >> ${sampleID}.sniffles_sorted_prefix.vcf;