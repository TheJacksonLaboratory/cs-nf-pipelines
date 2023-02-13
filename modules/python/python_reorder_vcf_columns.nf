process REORDER_VCF_COLUMNS {
  tag "$sampleID"

  cpus 1
  memory 4.GB
  time '04:00:00'

  container 'quay.io/jaxcompsci/bedtools-python2:2.26.0'

  input:
  tuple val(sampleID), path(vcf), path(idx), val(meta)
  
  output:
  tuple val(sampleID), path("*_mnv_final_filtered_merged_reordered.vcf"), val(meta), emit: vcf

  script:
  
  normal = meta.normal_id
  tumor = meta.tumor_id
  
  """
  python \
  ${projectDir}/bin/sv/reorder_vcf.py \
  ${vcf} \
  ${vcf.baseName}_mnv_final_filtered_merged_reordered.vcf \
  ${normal} ${tumor}
  """
}