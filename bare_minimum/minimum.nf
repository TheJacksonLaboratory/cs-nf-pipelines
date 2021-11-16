#!/usr/bin/env nextflow

params.fq_path_test = "$HOME/bitbucket/nf_dsl2_testing/bare_minimum/cache/rnaseqs" // default
params.fq_path
params.outdir = "."
println params.fq_path
println params.fq_path

process bar {
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
