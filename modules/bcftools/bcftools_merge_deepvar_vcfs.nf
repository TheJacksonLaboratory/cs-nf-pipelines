process BCFTOOLS_MERGEDEEPVAR {

  cpus 2
  memory 50.GB
  time {5.hour * task.attempt}
  errorStrategy 'retry' 
  maxRetries 2

  tag "$sampleID"
  errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

  container 'quay.io/biocontainers/bcftools:1.15--h0ea216a_2'

  publishDir "${params.pubdir}/${sampleID}", pattern: "*.*vcf*", mode:'copy'

  input:
  tuple val(sampleID), path(vcf), path(gvcf), path(vcf_index), path(gvcf_index)

  output:
  tuple val(sampleID), path("*_sorted_deepvariant.vcf.gz"), path("*_sorted_deepvariant.vcf.gz.tbi"), emit: vcf_idx
  tuple val(sampleID), path("*_sorted_deepvariant.gvcf.gz"), path("*_sorted_deepvariant.gvcf.gz.tbi"), emit: gvcf_idx

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
