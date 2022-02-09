#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// import workflow of interest
if (params.workflow == "toy_example"){
  include {TOY_EXAMPLE} from './workflows/toy_example'
}
if (params.workflow == "rnaseq"){
  include {RNASEQ} from './workflows/rnaseq'
}
if (params.workflow == "wes"){
  include {WES} from './workflows/wes'
}
if (params.workflow == "wgs"){
  include {WGS} from './workflows/wgs'
}
// conditional to kick off appropriate workflow
workflow{
  if (params.workflow == "toy_example"){
    TOY_EXAMPLE()
    }
  if (params.workflow == "rnaseq"){
    RNASEQ()
    }
  if (params.workflow == "wes"){
    WES()
    }
  if (params.workflow == "wgs"){
    WGS()
    }
}
