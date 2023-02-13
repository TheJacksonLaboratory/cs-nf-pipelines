process REMOVE_CONTIG {
  tag "$sampleID"

  cpus 1
  memory 4.GB
  time '04:00:00'

  container 'quay.io/jaxcompsci/bedtools-python3:2.26.0'

  input:
  tuple val(sampleID), file(vcf), val(meta), val(normal_name), val(tumor_name), val(tool)

  output:
  tuple val(sampleID), file("*.vcf"), val(meta), val(normal_name), val(tumor_name), val(tool), emit: remove_contig_vcf

  script:
  """
  python \
  ${projectDir}/bin/sv/remove_contig.py \
  ${vcf} \
  ${vcf.baseName}_removeContig.vcf
  """
}