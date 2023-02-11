process MERGE_COLUMNS {
  tag "$sampleID"

  cpus 1
  memory 4.GB
  time '04:00:00'

  container 'quay.io/jaxcompsci/bedtools-python2:2.26.0'

  input:
  tuple val(sampleID), file(vcf), val(meta), TBD

  output:
  tuple val(sampleID), file("*.vcf"), emit: single_column_chrom_pair_vcf

  script:

  normal = meta.normal_id
  tumor = meta.tumor_id

  """
   python \
  ${projectDir}/bin/sv/merge_columns.py \
  ${vcf} \
  ${sampleID}_single_column_${chrom}.vcf \
  ${normal} \
  ${tumor}
  """
}
