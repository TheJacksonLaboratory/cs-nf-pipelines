process BCF_SPLITMULTIALLELICREGIONS {
  tag "$sampleID"

  cpus = 1
  memory = 6.GB
  time = '06:00:00'

  container 'quay.io/biocontainers/bcftools:1.15--h0ea216a_2'

  input:
  tuple val(sampleID), path(interval), val(index)

  output:
  tuple val(sampleID), file("*.vcf"), emit: vcf

  script:

  """
  bcftools \
  norm \
  -m \
  -any \
  --threads ${task.cpus} \
  --regions ${index} \
  --no-version \
  -f ${params.ref_fa} \
  -o ${sampleID}_split.vcf \
  ${vcf}
"""
}