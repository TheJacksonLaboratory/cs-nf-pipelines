#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// import workflow of interest
if (params.workflow == "generate_pseudoreference"){
  include {GENERATE_PSEUDOREFERENCE} from './workflows/generate_pseudoreference'
}
if (params.workflow == "prepare_emase"){
  include {PREPARE_EMASE} from './workflows/prepare_emase'
}
if (params.workflow == "emase"){
  include {EMASE} from './workflows/emase'
}
if (params.workflow == "gbrs"){
  include {GBRS} from './workflows/gbrs'
}
if (params.workflow == "prep_do_gbrs_inputs"){
  include {PREP_DO_GBRS_INPUT} from './subworkflows/prep_do_gbrs_inputs'
}

// conditional to kick off appropriate workflow
workflow{
  if (params.workflow == "generate_pseudoreference") {
    GENERATE_PSEUDOREFERENCE()
  }

  if (params.workflow == "prepare_emase"){
    PREPARE_EMASE()
  }
 
  if (params.workflow == "emase"){
    EMASE()
  }

  if (params.workflow == "gbrs"){
    GBRS()
  }

  if (params.workflow == "prep_do_gbrs_inputs"){
    PREP_DO_GBRS_INPUT()
  }

}
