process BCF_MERGECALLERS {
  tag "$sampleID"

  cpus = 1
  memory = 6.GB
  time = '06:00:00'

  container 'quay.io/biocontainers/bcftools:1.15--h0ea216a_2'

  input:
  tuple val(sampleID), file(vcf)

  output:
  tuple val(sampleID), file("*.vcf"), emit: vcf

  script:

  """
  bcftools \
  merge \
  -r ${chrom} \
  --force-samples \
  --no-version \
  --threads ${task.cpus} \
  -f PASS, SUPPORT
  -F x \
  -m none \
  -o ${sampleID}_${chrom}.vcf \
  -i called_by:join,num_callers:sum,MNV_ID:join,supported_by:join \
  ${vcf}
"""
}