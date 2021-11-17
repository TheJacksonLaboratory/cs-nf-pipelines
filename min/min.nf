#!/usr/bin/env nextflow

// print a param from config
println params.description

// create a channel of pair reads
reads_ch = Channel.fromFilePairs("${params.fq_path}/*_R{1,2}${params.extension}")
println reads_ch.view()

// define process: keep it one tool and one container per process
process trim {
    
  input:
  
  output: 
  script:  
  """
  echo ouch
  """
}
