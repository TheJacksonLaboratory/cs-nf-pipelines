#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// import modules
include {CHECK_DESIGN} from '../modules/utility_modules/chipseq_check_design'
include {SAMTOOLS_FAIDX} from '../modules/samtools/samtools_faidx'
include {MAKE_GENOME_FILTER} from '../modules/utility_modules/chipseq_make_genome_filter'
include {FASTQC} from '../modules/fastqc/fastqc'
include {TRIM_GALORE} from '../modules/trim_galore/trim_galore'
include {READ_GROUPS} from '../modules/utility_modules/read_groups'
include {BWA_MEM} from '../modules/bwa/bwa_mem'
include {SAMTOOLS_FILTER} from '../modules/samtools/samtools_filter'

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


  // Reference genome
  ch_fasta = file(params.fasta, checkIfExists: true)


  // Step 2: Make genome filter
  SAMTOOLS_FAIDX(ch_fasta)
  MAKE_GENOME_FILTER(SAMTOOLS_FAIDX.out, params.blacklist)

  // Step 3: Fastqc
  FASTQC(read_ch)
  
  // Step 4: Trim Galore
  TRIM_GALORE(read_ch)

  // Step 5: Get Read Group Information
  READ_GROUPS(TRIM_GALORE.out.trimmed_fastq, "gatk")

  // Step 6: BWA-MEM
  bwa_mem_mapping = TRIM_GALORE.out.trimmed_fastq.join(READ_GROUPS.out.read_groups)
  BWA_MEM(bwa_mem_mapping)

  // Step 7: Samtools Removing Unmapped
  SAMTOOLS_FILTER(BWA_MEM.out, '-F 0x0100')

}
