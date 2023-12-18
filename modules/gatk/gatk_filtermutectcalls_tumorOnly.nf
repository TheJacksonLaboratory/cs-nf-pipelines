process GATK_FILTERMUECTCALLS {
  tag "$sampleID"

  cpus = 1
  memory = 15.GB
  time '05:00:00'
  errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

  container 'broadinstitute/gatk:4.4.0.0'

  publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID : 'gatk' }", pattern: "*_mutect2_somatic.filtered.vcf.gz", mode:'copy'
  publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID+'/stats'  : 'gatk' }", pattern: "*.filteringStats.tsv", mode:'copy'

  input:
  tuple val(sampleID), path(vcf), path(tbi), path(stats)

  output:
  tuple val(sampleID), file("*_mutect2_somatic.filtered.vcf.gz"), file("*_mutect2_somatic.filtered.vcf.gz.tbi"), emit: mutect2_vcf_tbi
  tuple val(sampleID), file("*.filteringStats.tsv"), emit: stats

  script:
  //Estimate somatic variants using Mutect2
  String my_mem = (task.memory-1.GB).toString()
  my_mem =  my_mem[0..-4]
  """
  gatk --java-options "-Xmx${my_mem}G" FilterMutectCalls \
    -R ${params.ref_fa} \
    -V ${vcf} \
    --stats ${stats} \
    -O ${sampleID}_mutect2_somatic.filtered.vcf.gz
  """
}