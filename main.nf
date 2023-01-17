#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// import workflow of interest
if (params.workflow == "prepare_emase"){
  include {PREPARE_EMASE} from './workflows/prepare_emase'
}
if (params.workflow == "emase"){
  include {EMASE} from './workflows/emase'
}

// conditional to kick off appropriate workflow
workflow{
  if (params.workflow == "prepare_emase"){
    PREPARE_EMASE()
    }
  if (params.workflow == "emase"){
    EMASE()
    }

}
