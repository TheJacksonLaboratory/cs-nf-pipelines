process CAT_ANNOTATE {
  
  cpus = 1
  memory = 10.GB
  time = '23:00:00'
  clusterOptions = 'q batch'
  
  input:
  tuple val(sampleID), file(snp_vcf)

  output:
  tuple val(sampleID), file("*.vcf"), emit: vcf

  container 'gatk-4.1.9.0_samtools-1.3.1_snpEff_4.3_vcftools_bcftools.sif'

  script:
  """
  cat ${snp_vcf} | vcf-annotate -a ${params.dbSNP} -c CHROM,FROM,TO,ID > ${sampleID}_variants_filtered_dbsnp.vcf
  """
}
