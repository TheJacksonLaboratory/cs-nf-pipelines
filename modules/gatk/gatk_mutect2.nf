process GATK_MUTECT2 {
  tag "$sampleID"

  cpus = 1
  memory = 15.GB
  time = '10:00:00'

  container 'broadinstitute/gatk:4.0.5.1'

  publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID : 'gatk' }", pattern: "*_somatic.vcf.gz", mode:'copy'

  input:
  tuple val(sampleID), file(tumor_bam)
  tuple val(sampleID), file(normal_bam)

  output:
  tuple val(sampleID), file("*_somatic.vcf.gz"), emit: somatic_vcf

  script:
  log.info "----- GATK Mutect2 Running on: ${sampleID} -----" 
  
  """
  gatk --java-options "-Xmx${my_mem}G" Mutect2 \\
    -R ${params.ref} \\
    -I ${tumor_bam} \\
    -tumor ${sampleID} \\
    -I ${normal_bam} \\
    -normal ${sampleID} \\
    -O ${sampleID}_somatic.vcf.gz
  """
