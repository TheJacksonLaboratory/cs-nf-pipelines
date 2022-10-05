#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// import modules
include {CHECK_DESIGN} from '../modules/utility_modules/chipseq_check_design'



// main workflow
workflow CHIPSEQ {

  if (params.input)     { ch_input = file(params.input, checkIfExists: true) } else { exit 1, 'Samples design file not specified!' }

  // Step 1: CHECK_DESIGN
  CHECK_DESIGN(ch_input)


}
