#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// import modules
include {help} from "${projectDir}/bin/help/generate_transcriptome.nf"
include {param_log} from "${projectDir}/bin/log/generate_transcriptome.nf"
include {EMASE_PREPARE_EMASE} from "${projectDir}/modules/emase/emase_prepare_emase"


// help if needed
if (params.help){
    help()
    exit 0
}

// log params
param_log()

// main workflow
workflow GENERATE_TRANSCRIPTOME {
    EMASE_PREPARE_EMASE()    
}