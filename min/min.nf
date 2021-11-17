#!/usr/bin/env nextflow

// enable DSL2
nextflow.enable.dsl=2

// log important info
log.info """\
 R N A S E Q - N F   P I P E L I N E
 ===================================
 fq_path        : ${params.fq_path}
 outdir         : ${params.outdir}
 """

// create a channel of pair reads
read_ch = Channel.fromFilePairs("${params.fq_path}/*_R{1,2}${params.extension}",checkExists:true )

/*
  define process: keep it one tool and one container per process
  the process container and publishDir (output directory) are defined in nextflow.config
*/

process trim {
  
  // required: this is where you define the channel to be used and variable names
  input:
  tuple val(sampleId), file(reads)

/*
  required: this is where you define the channel to be created from variables
  output will be specific to the container so you must know how the container functions
*/

  output:
  file "*.txt"
  /*
  required: the script/command entered here will be run by the container
  note that you can surround the script in """echo here""" or you can point
  at a script in the templates folder (folder must be named templates)
  */

  script:
  template "trim.sh"

}
workflow{
 data=read_ch
 trim(data)
}

workflow.onComplete {
	log.info ( workflow.success ? "\nDone! Output in : $params.outdir\n" : "Oops .. something went wrong" )
}
