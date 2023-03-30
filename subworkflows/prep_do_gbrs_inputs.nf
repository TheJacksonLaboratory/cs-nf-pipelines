#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// import modules
include {help} from "${projectDir}/bin/help/prep_do_gbrs_input.nf"
include {param_log} from "${projectDir}/bin/log/prep_do_gbrs_input.nf"
include {DO_TRANSITION_PROBABILITIES} from "${projectDir}/modules/R/do_transition_probablities"
include {PARSE_TRANSITION_PROBABILITIES} from "${projectDir}/modules/python/parse_transprobs"
include {PARSE_GENE_POSITONS} from "${projectDir}/modules/python/parse_gene_positions"
include {GENERATE_GRID_FILE} from "${projectDir}/modules/R/generate_grid_file"

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
    PARSE_GENE_POSITONS(DO_TRANSITION_PROBABILITIES.out.gene_list_tsv)
    GENERATE_GRID_FILE()
}
