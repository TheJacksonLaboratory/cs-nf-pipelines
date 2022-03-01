process CAT_HUMAN{
  tag "$sampleID"

  cpus 1
  memory 2.GB
  time '00:10:00'
  clusterOptions '-q batch'

  container 'gatk-4.1.6.0_samtools-1.3.1_snpEff_4.3_vcftools_bcftools.sif'

  input:
  tuple val(sampleID), file(vcf)
  val(indel_snp)

  output:
  tuple	val(sampleID), file("*.vcf"), emit: vcf

  script:
  // the pl in here needs to be discovered. this will happen when making the container cook book
  """
  cat ${vcf} | /snpEff_v4_3/snpEff/scripts/vcfEffOnePerLine.pl > ${sampleID}_oneperline_${indel_snp}.vcf
  """
}
