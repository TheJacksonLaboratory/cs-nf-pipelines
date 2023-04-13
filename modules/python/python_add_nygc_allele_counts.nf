process ADD_NYGC_ALLELE_COUNTS {
  tag "$sampleID"

  cpus 1
  memory 120.GB
  time '24:00:00'

  container 'quay.io/jaxcompsci/bedtools-python3:2.26.0'

  input:
  tuple val(sampleID), file(vcf), val(meta), file(normal_bam), file(normal_bai), file(tumor_bam), file(tumor_bai), val(chrom)

  output:
  tuple val(sampleID), file("*.vcf"), val(meta), val(chrom), emit: vcf

  script:
  """
   python \
  ${projectDir}/bin/pta/add_nygc_allele_counts_to_vcf.py \
  -t ${tumor_bam} \
  -n ${normal_bam} \
  -v ${vcf} \
  -b 10 \
  -m 10 \
  -o ${sampleID}_pre_count_${chrom}.vcf
  """
}
