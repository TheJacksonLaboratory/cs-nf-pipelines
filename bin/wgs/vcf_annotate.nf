process VCF_ANNOTATE {
  tag "$sampleID"

  cpus = 1
  memory = 10.GB
  time = '23:00:00'
  clusterOptions = 'q batch'

  input:
  tuple val(sampleID), file(snp_vcf)
  val(indel_snp)

  output:
  tuple val(sampleID), file("*.vcf"), emit: vcf

  // vcftools container needed
  container 'gatk-4.1.9.0_samtools-1.3.1_snpEff_4.3_vcftools_bcftools.sif'

  script:
  log.info "----- CAT VCF-ANNOTATE Running on: ${sampleID} -----"

  if (params.gen_org=='mouse'){
    // make sure it does not break
    delta="CHROM,POS,ID,REF,ALT"
  }
  else if (params.gen_org=='human'){
    delta="CHROM,POS,ID,REF,ALT"
  }

  if (indel_snp == 'INDEL'){
    output_suffix = 'indel_filtered_dbsnp.vcf'
  }
  if (indel_snp =='SNP'){
    output_suffix = 'snp_filtered_dbsnp.vcf'
  }
  if (indel_snp == 'BOTH'){
    output_suffix = 'snp_indel_filtered_dbsnp.vcf'
  }

  """
  cat ${snp_vcf} | vcf-annotate -a ${params.dbSNP} -c ${delta} > ${sampleID}_${output_suffix}
  """
}
