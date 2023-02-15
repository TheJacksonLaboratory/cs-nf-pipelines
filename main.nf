#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// import workflow of interest
if (params.workflow == "rnaseq"){
  include {RNASEQ} from './workflows/rnaseq'
}
if (params.workflow == "wes"){
  include {WES} from './workflows/wes'
}
if (params.workflow == "pdx_wes"){
  include {PDX_WES} from './workflows/pdx_wes'
}
if (params.workflow == "wgs"){
  include {WGS} from './workflows/wgs'
}
if (params.workflow == "rrbs"){
  include {RRBS} from './workflows/rrbs'
}
if (params.workflow == "atac"){
  include {ATAC} from './workflows/atac'
}
if (params.workflow == "sv"){
  include {SV} from './workflows/sv'
}
// conditional to kick off appropriate workflow
workflow{
  if (params.workflow == "rnaseq"){
    RNASEQ()
    }
  if (params.workflow == "wes"){
    WES()
    }
  if (params.workflow == "pdx_wes"){
    PDX_WES()
  }
  if (params.workflow == "wgs"){
    WGS()
    }
  if (params.workflow == "rrbs"){
    RRBS()
    }
  if (params.workflow == "atac"){
    ATAC()
    }
  if (params.workflow == "sv"){
    SV()
  } 
}
