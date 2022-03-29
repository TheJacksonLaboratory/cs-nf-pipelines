#!/bin/bash

## vcf cli params
export in_vcf=$1
export out_vcf=$2

grep '^#' ${in_vcf} > ${out_vcf} && \
grep '^chr' ${in_vcf} | \
LC_ALL=C sort -t $'\t' -k1,1 -k2,2n >> ${out_vcf} 
