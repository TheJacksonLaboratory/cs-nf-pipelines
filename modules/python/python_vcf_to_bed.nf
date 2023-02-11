process VCF_TO_BED {
  tag "$sampleID"

  cpus 1
  memory 4.GB
  time '04:00:00'

  container 'quay.io/jaxcompsci/bedtools-python2:2.26.0'

  input:
  tuple val(sampleID), file(vcf), val(meta), val(chrom)

  output:
  tuple val(sampleID), file("*.bed"), val(meta), val(chrom), emit: bed

  script:
  """
   python \
  ${projectDir}/bin/sv/vcf_to_bed.py \
  ${vcf} \
  | bedtools \
  merge \
  > ${sampleID}_candidate_merged_${chrom}.bed 
  """
}