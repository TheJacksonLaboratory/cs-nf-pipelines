process CHAIN_FILTER_READS {
  tag "$sampleID"

  cpus 2
  memory 4.GB
  time = '10:00:00'
  errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

  publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID+'/stats' : 'gatk' }", pattern: "*.log", mode: 'copy'
  container 'broadinstitute/gatk:4.2.4.1'

  input:
  tuple val(sampleID), file(bam_sort_mm10), file(ReadName_unique)
  

  output:
  tuple val(sampleID), path("*.tmp2.mm10.bam"), emit: bam
  tuple val(sampleID), file("*_FilterSamReads.log"), emit: filterReads_log

  when: params.chain != null

  script:
  """
  gatk FilterSamReads \
  -I ${bam_sort_mm10[0]} \
  -RLF ReadName_unique \
  --FILTER excludeReadList \
  --VALIDATION_STRINGENCY LENIENT \
  -O ${sampleID}.tmp2.mm10.bam \
  > ${sampleID}_FilterSamReads.log 2>&1
  """
}
