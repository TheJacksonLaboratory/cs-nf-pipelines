#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// import modules
include {CHECK_DESIGN} from '../modules/utility_modules/chipseq_check_design'



// main workflow
workflow CHIPSEQ {

  if (params.input)     { ch_input = file(params.input, checkIfExists: true) } else { exit 1, 'Samples design file not specified!' }

  // Step 1: CHECK_DESIGN
  CHECK_DESIGN(ch_input)

  /*
  * Create channels for input fastq files
  */

  if (params.read_type == 'SE'){
      read_ch = CHECK_DESIGN.out[0]
                   .splitCsv(header:true, sep:',')
                   .map { row -> [ row.sample_id, [ file(row.fastq_1, checkIfExists: true) ] ] }
  } else {
      read_ch = CHECK_DESIGN.out[0]
                   .splitCsv(header:true, sep:',')
                   .map { row -> [ row.sample_id, [ file(row.fastq_1, checkIfExists: true), file(row.fastq_2, checkIfExists: true) ] ] }
  }

  /*
   * Create a channel with [sample_id, control id, antibody, replicatesExist, multipleGroups]
   */
  control_ch = CHECK_DESIGN.out[1]
      .splitCsv(header:true, sep:',')
      .map { row -> [ row.sample_id, row.control_id, row.antibody, row.replicatesExist.toBoolean(), row.multipleGroups.toBoolean() ] }


  // Test to see if sample sheet inputs are OK.

  println "\nChecking Sample Sheet Inputs : \n"

  read_ch.view { "value: $it" }
  control_ch.view { "value: $it" }


}
