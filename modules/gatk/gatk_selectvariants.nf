process GATK_SELECTVARIANTS {
  tag "$sampleID"

  cpus = 1
  memory = 6.GB
  time = '03:00:00'
  errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

  container 'broadinstitute/gatk:4.2.4.1'

  publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID : 'gatk' }", pattern: "*.vcf", mode:'copy', enabled: params.keep_intermediate
  publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID : 'gatk' }", pattern: "*filtered_dbsnpID.vcf", mode:'copy'

  input:
  tuple val(sampleID), file(vcf), file(idx)
  val(indel_snp)
  val(suffix)

  output:
  tuple val(sampleID), file("*.vcf"), emit: vcf
  tuple val(sampleID), file("*.idx"), emit: idx

  script:
  String my_mem = (task.memory-1.GB).toString()
  my_mem =  my_mem[0..-4]
  """
  gatk --java-options "-Xmx${my_mem}G" SelectVariants \
  -R ${params.ref_fa} \
  -V ${vcf} \
  -select-type ${indel_snp} \
  -O ${sampleID}_${suffix}.vcf
  """
}