#!/usr/bin/env nextflow

// 1.) enable DSL2
nextflow.enable.dsl=2

// 2.) log important info

log.info """\
 TOY R N A S E Q - N F   P I P E L I N E
 ===================================
 fq_path        : ${params.fq_path}
 outdir         : ${params.outdir}
 """

// 3.) create a channel of pair reads

read_ch = Channel.fromFilePairs("${params.fq_path}/*_R{1,2}${params.extension}",checkExists:true )

/*
   4.) define trimming process: keep it one tool and one container per process
   the process container and publishDir (output directory) are defined in nextflow.config
*/

process trim {

  // 4.a.) required: this is where you define the channel to be used and variable names

  input:
  tuple val(sampleId), file(reads)

  /*
     4.b.) required: this is where you define the channel to be created from variables
     output will be specific to the container so you must know how the container functions
  */

  output:
  file "*.fastq.gz"

  /*
     4.c.) required: the script/command entered here will be run by the container
     note that you can surround the script in """echo here""" or you can point
     at a script in the templates folder (folder must be named templates)


  */

  script:
  """
  echo ${params.t_lead} > ${sampleId}.txt

  trimmomatic \
  PE \
  /home/guglib/rnaseqs/PE/${reads[0]} \
  /home/guglib/rnaseqs/PE/${reads[1]} \
  /home/guglib/test/output_forward_paired.fastq.gz \
  /home/guglib/test/output_forward_unpaired.fastq.gz \
  /home/guglib/test/output_reverse_paired.fastq.gz \
  /home/guglib/test/output_reverse_unpaired.fastq.gz \
  LEADING:${t_lead} \
  TRAILING:${t_trail} \
  MINLEN:${min_len}
  """

}

// 5. Use RSEM for quantification

process quant{
input:
tuple val(sampleId), file(reads)
output:
file "*.txt"
}

workflow{
 data=read_ch
 trim(data)
}

workflow.onComplete {
	log.info ( workflow.success ? "\nDone! Output in : $params.outdir\n" : "Oops .. something went wrong" )
}
