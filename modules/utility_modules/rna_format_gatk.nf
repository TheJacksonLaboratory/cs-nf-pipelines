process FORMAT_GATK {
  tag "$sampleID"

  cpus 1
  memory 15.GB
  time '24:00:00'

  container 'quay.io/biocontainers/bedtools:2.23.0--h5b5514e_6'
  file(params.ref_fai)

  input:
  tuple val(sampleID), file(txt)
  val(L)

  output:
  tuple val(sampleID), file("*_gatk_formatter.txt"), emit: txt

  script:
  """
  chmod +x ${projectDir}/bin/rnaseq/gatk_formatter.sh
  ${projectDir}/bin/rnaseq/gatk_formatter.sh ${txt} ${sampleID}_gatk_temp2.txt ${sampleID}_gatk_formatter.txt ${L}
  """
  // This is a script to format gatk coverage file for subsequent use in log aggregation 
}