process BISMARK_ALIGNMENT {
  tag "$sampleID"

  cpus 8
  memory {60.GB * task.attempt}
  time {30.hour * task.attempt}
  errorStrategy 'retry' 
  maxRetries 1

  container 'CONTAINER_TBD'

  publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID : 'bwa_mem' }", pattern: "*.sam", mode:'copy'

  input:
  tuple val(sampleID), file(fq_reads)

  output:
  tuple val(sampleID), file("*.sam"), emit: sam

  script:
  log.info "----- Bismark Alignment Running on: ${sampleID} -----"

  if (params.read_type == "SE"){
    inputfq="-1 ${fq_reads[0]}"
    }
  if (params.read_type == "PE"){
    inputfq="-1 ${fq_reads[0]} -2 ${fq_reads[1]}"
    }

  if ${params.non_directional} {
    directionality = '--non_directional'
  } 

  """
  bismark --bowtie2 -p ${task.cpus} ${directionality} -L ${params.seedlength} -N ${params.seed_mismatch} -minins ${params.MinInsert} -maxins ${params.MaxInsert}  --output_dir {in_4} --unmapped --ambiguous ${params.ref_fa_index} ${inputfq}
  """ // NOTE: OUTPUT DIR IS CALLED....IS THIS A NEEDED THING? HOW IS OUTPUT NAMED? 
}


