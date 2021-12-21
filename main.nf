#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// import workflow of interest
if (params.workflow == "toy_example"){
  include {TOY_EXAMPLE} from './workflows/toy_example'
}
else if (params.workflow == "rnaseq"){
  include {RNASEQ} from './workflows/rnaseq'
}

// conditional to kick off appropriate workflow
workflow{
if (params.workflow == "toy_example"){
  TOY_EXAMPLE()
  }
else if (params.workflow == "rnaseq"){
  RNASEQ()
  }
}
