process MERGE_COLUMNS_TUMOR {
  tag "$sampleID"

  cpus 1
  memory 4.GB
  time '04:00:00'

  container 'quay.io/jaxcompsci/bedtools-python2:2.26.0'

  input:
  tuple val(sampleID), val(meta), file(vcf)

  output:
  tuple val(sampleID), file("*.vcf"), emit: single_column_chrom_tumor_vcf

  script:
  """
   python \
  ${projectDir}/bin/sv/merge_columns_pon.py \
  ${supportedChromVcf} \
  ${tumor}_single_column_${chrom}.vcf \
  ${tumor} 
  """
}