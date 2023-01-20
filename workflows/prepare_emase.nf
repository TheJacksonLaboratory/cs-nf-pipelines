#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// import modules
include {help} from "${projectDir}/bin/help/prepare_emase.nf"
include {param_log} from "${projectDir}/bin/log/prepare_emase.nf"
include {EMASE_PREPARE_EMASE} from "${projectDir}/modules/emase/emase_prepare_emase"
include {CLEAN_TRANSCRIPT_LISTS} from "${projectDir}/modules/utility_modules/clean_prepEmase_transcriptList"


// help if needed
if (params.help){
    help()
    exit 0
}

// log params
param_log()

// main workflow
workflow PREPARE_EMASE {
    EMASE_PREPARE_EMASE()
    CLEAN_TRANSCRIPT_LISTS(EMASE_PREPARE_EMASE.out.pooled_transcript_info)
}