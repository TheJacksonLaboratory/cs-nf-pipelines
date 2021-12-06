#!/usr/bin/env nextflow
nextflow.enable.dsl=2
include {TRIM} from './modules/trim'
include {RSEM_EXPRESSION; RSEM_SIMULATE_READS} from './modules/rsem'

read_ch = Channel.fromFilePairs("${params.fq_path}/*_R{1,2}${params.extension}",checkExists:true )

workflow{
 ref_files = file("${params.ref_files}/*")
 TRIM(read_ch)
 RSEM_EXPRESSION(TRIM.output, ref_files)
 // RSEM_SIMULATE_READS(TRIM.output, ref_files)
}
