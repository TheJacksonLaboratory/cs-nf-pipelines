process ADD_FINAL_ALLELE_COUNTS {
  tag "$sampleID"

  cpus 1
  memory 4.GB
  time '04:00:00'

  container 'quay.io/jaxcompsci/bedtools-python2:2.26.0'

  input:
  tuple val(sampleID), file(vcf)

  output:
  tuple val(sampleID), file("*.vcf"), emit: merge_final_vcf

  script:
  """
   python \
  ${projectDir}/bin/sv/add_final_allele_counts_to_vcf.py \
  -v ${sampleID}_pre_count_${chrom}.vcf \
  -o ${sampleID}_final_${chrom}.vcf \
  """
}