process ADD_NYGC_ALLELE_COUNTS_PAIR {
  tag "$sampleID"

  cpus 1
  memory 4.GB
  time '04:00:00'

  container 'quay.io/jaxcompsci/bedtools-python2:2.26.0'

  input:
  tuple val(sampleID), val(meta), file(vcf)

  output:
  tuple val(sampleID), file("*.vcf"), emit: add_nygc_allele_counts_pair_vcf

  script:
  """
   python \
  ${projectDir}/bin/sv/add_nygc_allele_counts_to_vcf.py \
  -t ${tumor_bam} \
  -n ${normal_bam} \
  -v ${sampleID}_single_column_${chrom}.vcf \
  -b 10 \
  -m 10 \
  -o ${sampleID}_pre_count_${chrom}.vcf
  """
}