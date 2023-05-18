process FASTQ_SORT {

  tag "$sampleID"

  cpus 1
  memory 50.GB
  time 2.h
  errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.mem} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

  container 'quay.io/biocontainers/fastq-tools:0.8.3--hbd632db_2'

  publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID + '/deconvoluted_reads': 'deconvoluted_reads' }", pattern: "*.fastq", mode:'copy', enabled: params.keep_intermediate

  input:
  tuple val(sampleID), file(reads)
  val(suffix)

  output:
  tuple val(sampleID), file("*sorted*{1,2}.fastq"), emit: sorted_fastq

  script:
  command_two = params.read_type == 'PE' ? "fastq-sort --id ${reads[1]} > ${sampleID}_sorted_${suffix}_2.fastq" : ''

  """
  fastq-sort --id ${reads[0]} > ${sampleID}_sorted_${suffix}_1.fastq
  ${command_two}
  """
}
