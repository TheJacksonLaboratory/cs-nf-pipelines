#!/usr/bin/env nextflow

// 1.) enable DSL2
nextflow.enable.dsl=2

// 2.) bring in modules
include {RSEM_REF_PULL;RSEM_REF_BUILD;RSEM_EXPRESSION} from '../modules/rsem'

// 3.) reads channel
read_ch = Channel.fromFilePairs("${params.fq_path}/*_R{1,2}_*${params.extension}",checkExists:true )

// 4.) The main workflow
workflow{
 if( params.ref_build == 'true')
   ref_files = RSEM_REF_BUILD(RSEM_REF_PULL())
 else if ( params.ref_build == 'false' )
   ref_files = file("${params.ref_files}/*")

 TRIM(read_ch)
 RSEM_EXPRESSION(TRIM.output, ref_files)
}
