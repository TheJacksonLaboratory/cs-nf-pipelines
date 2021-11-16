#!/usr/bin/env nextflow

println params.contact

process trim {
  println "Trimmomatic"

"""
trimmomatic PE \
'${params.fq_path_test}'/SRR1641108.1_R1.fastq.gz \
'${params.fq_path_test}'/SRR1641108.1_R2.fastq.gz \
output/output_forward_paired.fq.gz \
output/output_forward_unpaired.fq.gz \
output/output_reverse_paired.fq.gz \
output/output_reverse_unpaired.fq.gz \
LEADING:3 \
TRAILING:3 \
MINLEN:36
"""
}
