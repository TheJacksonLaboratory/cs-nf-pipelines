process BCFTOOLS_MERGEDEEPVAR {

  cpus 2
  memory 50.GB
  time {5.hour * task.attempt}
  errorStrategy 'retry' 
  maxRetries 2

  tag "$sampleID"
  errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

  container 'quay.io/biocontainers/bcftools:1.15--h0ea216a_2'

  input:
  tuple val(sampleID), file(vcf), file(gvcf), file(vcf_index), file(gvcf_index)

  output:
  tuple val(sampleID), file("*_sorted_deepvariant.vcf.gz"), file("*_sorted_deepvariant.vcf.gz.tbi"), emit: vcf
  tuple val(sampleID), file("*_sorted_deepvariant.gvcf.gz"), file("*_sorted_deepvariant.gvcf.gz"), emit: gvcf

  script:

  """
  ## GVCF
  # concatenate  gvcfs
  bcftools concat ${gvcf} --allow-overlaps --remove-duplicates -Oz -o ${sampleID}_deepvariant.gvcf.gz 

  # sort
  bcftools sort ${sampleID}_deepvariant.gvcf.gz -Oz -o ${sampleID}_sorted_deepvariant.gvcf.gz

  # index vcfs
  bcftools index --tbi ${sampleID}_sorted_deepvariant.gvcf.gz

  ## VCF
  # concatenate  vcfs
  bcftools concat ${vcf} --allow-overlaps --remove-duplicates -Oz -o ${sampleID}_deepvariant.vcf.gz 

  # sort
  bcftools sort ${sampleID}_deepvariant.vcf.gz -Oz -o ${sampleID}_sorted_deepvariant.vcf.gz

  # index vcfs
  bcftools index --tbi ${sampleID}_sorted_deepvariant.vcf.gz

  """
}
