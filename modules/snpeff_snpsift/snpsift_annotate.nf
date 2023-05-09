process SNPSIFT_ANNOTATE {
  tag "$sampleID"

  cpus = 1
  memory = 6.GB
  time = '06:00:00'
  errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.mem} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

  container 'quay.io/jaxcompsci/snpeff_snpsift_5.1:v5.1d'
  
  publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID : 'snpsift' }", pattern:"*dbsnpID.vcf", mode:'copy'
  publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID : 'snpeff' }", pattern:"*.vcf", mode:'copy', enabled: params.workflow == 'amplicon' ? true : false

  input:
  tuple val(sampleID), file(vcf)
  path(annot_source)
  path(annot_index)
  val(output_suffix)

  output:
  tuple val(sampleID), file("*.vcf"), emit: vcf

  script:

  String my_mem = (task.memory-1.GB).toString()
  my_mem =  my_mem[0..-4]

  """
  java -Xmx${my_mem}G -jar /opt/snpEff/SnpSift.jar \
  annotate -noDownload -id ${annot_source} ${vcf} > ${vcf.baseName}_${output_suffix}.vcf
  """
}
