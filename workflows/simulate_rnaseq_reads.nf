#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include {TRIM} from '../modules/trimmomatic'
include {RSEM_EXPRESSION, RSEM_SIMULATE_READS} from '../modules/rsem'
read_ch = Channel.fromFilePairs("${params.sample_folder}/*_R{1,2}_*${params.extension}",checkExists:true )

workflow SIMULATE_RNASEQ {
 ref_files = file("${params.ref_files}/*")
 TRIM(read_ch)
 RSEM_EXPRESSION(TRIM.output, ref_files)
 RSEM_SIMULATE_READS(RSEM_EXPRESSION.output)
}
