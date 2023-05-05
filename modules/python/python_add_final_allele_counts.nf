process ADD_FINAL_ALLELE_COUNTS {
  tag "$sampleID"

  cpus 1
  memory 4.GB
  time '04:00:00'

  container 'quay.io/jaxcompsci/bedtools-python3:2.26.0'

  input:
  tuple val(sampleID), file(vcf), val(meta), val(chrom)

  output:
  tuple val(sampleID), file("*.vcf"), val(meta), val(chrom), emit: vcf

  script:
  """
   python \
  ${projectDir}/bin/pta/add_final_allele_counts_to_vcf.py \
  -v ${vcf} \
  -o ${sampleID}_final_${chrom}.vcf \
  """
}