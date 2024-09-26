process CHAIN_CONVERT {
  tag "$sampleID"

  cpus 1
  memory 10.GB
  time '10:00:00'
  errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

  container 'quay.io/jaxcompsci/g2gtools:0.1.31'

  publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID+'/stats' : 'g2gtools' }", pattern: "*.log", mode:'copy'

  input:
  tuple val(sampleID), file(bam_shifted)

  output:
  tuple val(sampleID), file("*.tmp.mm10.bam"), emit: converted_bam
  tuple val(sampleID), file("*g2gconvert.log"), emit: log

  when: params.chain != null

  script:
  """
  g2gtools convert \
  -r -f bam -c ${params.chain} \
  -i ${bam_shifted[0]} \
  -o ${sampleID}.tmp.mm10.bam 2> ${sampleID}_g2gconvert.log
  """
}
