process MERGE_COLUMNS {
  tag "$sampleID"

  cpus 1
  memory 4.GB
  time '04:00:00'

  container 'quay.io/jaxcompsci/bedtools-python3:2.26.0'

  input:
  tuple val(sampleID), file(vcf), file(tbi), val(meta), val('empty_name'), val('empty_name'), val(chrom)

  output:
  tuple val(sampleID), file("*.vcf"), val(meta), val(chrom), emit: mergeColumn_vcf

  script:

  normal = meta.normal_id
  tumor = meta.tumor_id

  """
   python \
  ${projectDir}/bin/pta/merge_columns.py \
  ${vcf} \
  ${sampleID}_single_column_${chrom}.vcf \
  ${normal} \
  ${tumor}
  """
}
  /*
  NOTE: This script will take 'tumor' and 'normal' names and match string based on a simple split on '_'. 
        Sample names are currently <tool>_<tumor_name> or <tool>_<normal_name>. This script merged based on the index[1] of split('_'). 
  */
