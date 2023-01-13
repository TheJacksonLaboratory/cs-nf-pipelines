#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// import workflow of interest
if (params.workflow == "generate_transcriptome"){
  include {GENERATE_TRANSCRIPTOME} from './workflows/generate_transcriptome'
}

// conditional to kick off appropriate workflow
workflow{
  if (params.workflow == "generate_transcriptome"){
    GENERATE_TRANSCRIPTOME()
    }
}
