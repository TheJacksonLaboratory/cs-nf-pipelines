process VCF_TO_BED {
  tag "$sampleID"

  cpus 1
  memory 4.GB
  time '04:00:00'

  container 'quay.io/jaxcompsci/bedtools-python2:2.26.0'

  input:
  tuple val(sampleID), file(vcf)

  output:
  tuple val(sampleID), file("*.vcf"), emit: candidate_merged_bed

  script:
  """
   python \
  ${projectDir}/bin/sv/vcf_to_bed.py \
  ${sampleID}_candidate_merged_${chrom}.vcf \
  | bedtools \
  merge \
  > ${sampleID}_candidate_merged_${chrom}.bed 
  """
}