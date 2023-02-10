process MERGE_COLUMNS_PAIR {
  tag "$sampleID"

  cpus 1
  memory 4.GB
  time '04:00:00'

  container 'quay.io/jaxcompsci/bedtools-python2:2.26.0'

  input:
  tuple val(sampleID), file(vcf)

  output:
  tuple val(sampleID), file("*.vcf"), emit: single_column_chrom_pair_vcf

  script:
  """
   python \
  ${projectDir}/bin/sv/merge_columns.py \
  ${supportedChromVcf} \
  ${sampleID}_single_column_${chrom}.vcf \
  ${tumor} \
  ${normal}
  """
}