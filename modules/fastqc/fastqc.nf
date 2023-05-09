
process FASTQC {
  tag "$sampleID"

  cpus 8
  memory 4.GB
  time '10:00:00'
  errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.mem} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

  container 'quay.io/biocontainers/fastqc:0.11.9--hdfd78af_1'
  
  publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID+'/stats' : 'fastqc' }", pattern: "*_fastqc.{zip,html}", mode:'copy'

  input:
  tuple val(sampleID), file(fq_reads)

  output:
  tuple val(sampleID), file("*_fastqc.{zip,html}"), emit: quality_stats


  script:

  """
    fastqc --quiet -t ${task.cpus} ${fq_reads}
  """
}
