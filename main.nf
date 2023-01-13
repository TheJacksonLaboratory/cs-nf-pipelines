#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// import workflow of interest
if (params.workflow == "prepare_emase"){
  include {prepare_emase} from './workflows/prepare_emase'
}

// conditional to kick off appropriate workflow
workflow{
  if (params.workflow == "prepare_emase"){
    PREPARE_EMASE()
    }
}
