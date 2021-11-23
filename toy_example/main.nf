#!/usr/bin/env nextflow

// 1.) enable DSL2
nextflow.enable.dsl=2

// 2.) log important info

log.info """\
 TOY EXAMPLE   P I P E L I N E
 ===================================
 fq_path        : ${params.fq_path}
 outdir         : ${params.outdir}
 """

// 3.) create a channel of pair reads

read_ch = Channel.fromFilePairs("${params.fq_path}/*_R{1,2}${params.extension}",checkExists:true )

/*
   4.) define trimming process: keep it one tool and one container per process
   the process container is defined in nextflow.config
*/

process trim {

  // 4.a.) required: this is where you define the channel to be used and variable names

  publishDir "${params.outdir}/trimmed"

  input:
  tuple val(sampleId), file(reads)

  /*
     4.b.) required: this is where you define the channel to be created from variables
     output will be specific to the container so you must know how the container functions
  */

  output:
  tuple val(sampleId),file('*.fastq.gz')


  /*
     4.c.) required: the script/command entered here will be run by the container
     note that you can surround the script in """echo here""" or you can point
     at a script in the templates folder (folder must be named templates)
  */

  script:
  """
  trimmomatic \
  PE \
  ${params.fq_path}/${reads[0]} \
  ${params.fq_path}/${reads[1]} \
  ${sampleId}_R1_paired.fastq.gz \
  ${sampleId}_R1_unpaired.fastq.gz \
  ${sampleId}_R2_paired.fastq.gz \
  ${sampleId}_R2_unpaired.fastq.gz \
  LEADING:${params.t_lead} \
  TRAILING:${params.t_trail} \
  MINLEN:${params.min_len}
  """

}

// 5. Use RSEM for quantification

process quant{
  publishDir "${params.outdir}/rsem"

  input:
  tuple val(sampleId), file(trimmed)


  output:
  file "*.txt"

  script:

  '''
  echo hello > world.txt
  '''
}


workflow{
 trim_ch=trim(read_ch)
 quant(trim_ch)
}

workflow.onComplete {
	log.info ( workflow.success ? "\nDone! Output in : $params.outdir\n" : "Oops .. something went wrong" )
}
