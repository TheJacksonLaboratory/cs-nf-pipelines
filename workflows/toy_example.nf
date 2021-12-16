#!/usr/bin/env nextflow

// 1.) enable DSL2
nextflow.enable.dsl=2

// 2.) bring in modules
include {TRIM} from '../modules/trimmomatic'
include {RSEM_REF_PULL;RSEM_REF_BUILD;RSEM_EXPRESSION} from '../modules/rsem'

// 3.) reads channel and reference files
read_ch = Channel.fromFilePairs("${params.fq_path}/*_R{1,2}_*${params.extension}",checkExists:true )
ref_files = file("${params.ref_files}/*")

// 4.) The main workflow
workflow TOY_EXAMPLE {
 TRIM(read_ch)
 RSEM_EXPRESSION(TRIM.output, ref_files)
}
