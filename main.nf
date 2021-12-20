#!/usr/bin/env nextflow

nextflow.enable.dsl=2


if (params.workflow == "toy_example"){
  include {TOY_EXAMPLE} from './workflows/toy_example'
  workflow{
    TOY_EXAMPLE()
  }
}

if (params.workflow == "RNASEQ"){
  include {RNASEQ} from './workflows/rnaseq'
  workflow{
    RNASEQ()
  }
}
