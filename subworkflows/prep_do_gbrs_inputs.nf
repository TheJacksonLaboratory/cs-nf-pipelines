#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// import modules
include {help} from "${projectDir}/bin/help/prep_do_gbrs_input.nf"
include {param_log} from "${projectDir}/bin/log/prep_do_gbrs_input.nf"
include {DO_TRANSITION_PROBABILITIES} from "${projectDir}/modules/R/do_transition_probablities"
include {PARSE_TRANSITION_PROBABILITIES} from "${projectDir}/modules/python/parse_transprobs"

// help if needed
if (params.help){
    help()
    exit 0
}

// log params
param_log()

// main workflow
workflow PREP_DO_GBRS_INPUT {
    DO_TRANSITION_PROBABILITIES()
    PARSE_TRANSITION_PROBABILITIES(DO_TRANSITION_PROBABILITIES.out.h5_files)
}
