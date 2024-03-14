#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// import modules
if (params.data_type == 'pacbio') {
  include {PACBIO} from "${projectDir}/subworkflows/pacbio"
}
else if (params.data_type == 'illumina') {
  include {ILLUMINA} from "${projectDir}/subworkflows/illumina"
}
else if (params.data_type == 'ont') {
  include {ONT} from "${projectDir}/subworkflows/ont"
}
else {
  // if it is not a validate data type, an acceptable string
  exit 1, "'--data_type': \"${params.data_type}\" is not valid, supported options are 'pacbio' or 'illumina' or 'ont'" 
}

// main workflow
workflow GERMLINE_SV {
  if (params.data_type == 'pacbio') {
    PACBIO()
  } else if (params.data_type == 'illumina') {
    ILLUMINA()
  } else if (params.data_type == 'ont') {
    ONT()
  } 
}
