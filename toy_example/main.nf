#!/usr/bin/env nextflow

// 1.) enable DSL2
nextflow.enable.dsl=2

// 2.) bring in modules
include 'modules/toy'

// 3.) reads channel
read_ch = Channel.fromFilePairs("${params.fq_path}/*_R{1,2}${params.extension}",checkExists:true )

// 4.) reference channel
// RSEM_REF_BUILD(RSEM_REF_PULL())
TRIM(read_ch)
// RSEM_EXPRESSION(TRIM.output, RSEM_REF_BUILD.output)
