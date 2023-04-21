process FILTER_PON {
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
  ${projectDir}/bin/pta/filter_pon.py \
        --bed ${params.pon_bed} \
        --chrom ${chrom} \
        --vcf ${vcf} \
        --out ${sampleID}_pon_final_${chrom}.vcf
  """
}
