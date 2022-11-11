process GATK_MUTECT2 {
  tag "$meta.patient"

  cpus = 1
  memory = 15.GB
  time = '10:00:00'

  container 'broadinstitute/gatk:4.0.5.1'

  publishDir "${params.pubdir}/${ params.organize_by=='sample' ? "$meta.patient" : 'gatk' }", pattern: "*_somatic.vcf.gz", mode:'copy'

  input:
  tuple val(sampleID), file(normal_bam), file(normal_bai), val(meta)
  tuple val(sampleID), file(tumor_bam), file(normal_bai), val(meta)

  output:
  tuple val(meta), file("*_somatic.vcf.gz"), emit: mutect_vcf
  tuple val(meta), file("*_somatic.vcf.gz.tbi"), emit: mutect_vcf_tbi

  script:
  log.info "----- GATK Mutect2 Running on: ${sampleID} -----" 
  //Estimate somatic variants using Mutect2
  
  """
  gatk --java-options "-Xmx${my_mem}G" Mutect2 \\
    -R ${params.ref_fa} \\
    -I ${tumor_bam} \\
    -tumor ${sampleID} \\
    -I ${normal_bam} \\
    -normal ${sampleID} \\
    -O ${meta.patient}_somatic.vcf.gz
  """
