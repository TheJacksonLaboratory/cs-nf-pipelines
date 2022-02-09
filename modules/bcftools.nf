process BCF_SORT {
  tag "sampleID"

  cpus = 1
  memory = 6.GB
  time = '06:00:00'
  clusterOptions = '-q batch'


  container 'gatk-4.1.9.0_samtools-1.3.1_snpEff_4.3_vcftools_bcftools.sif'

  input:
  tuple val(sampleID), file(vcf)

  output:
  tuple val(sampleID), file("*.vcf"), emit: vcf

  script:
  log.info "----- BCFTools Sort Running on: ${sampleID} -----"

  """
  bcftools sort -o ${sampleID}_only_${indel_snp}.vcf ${vcf}
  """

}
