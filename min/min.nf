#!/usr/bin/env nextflow

// print a param from config
println params.description

// create a channel of pair reads
read_ch = Channel.fromFilePairs("${params.fq_path}/*_R{1,2}${params.extension}")

// just a test
println reads_ch.view()

/*
  define process: keep it one tool and one container per process
  the process container and publishDir (output directory) are defined in nextflow.config
*/

process trim {
  // required: this is where you define the channel to be used and variable names
  input:
  file sampleId, reads from read_ch

/*
  required: this is where you define the channel to be created from variables
  output will be specific to the container so you must know how the container functions
*/

  output:
  file "*.txt" into file_ch

  /*
  required: the script/command entered here will be run by the container
  note that you can surround the script in """echo here""" or you can point
  at a script in the templates folder (folder must be named templates)
  */

  script:
  template 'trim.sh'

}
