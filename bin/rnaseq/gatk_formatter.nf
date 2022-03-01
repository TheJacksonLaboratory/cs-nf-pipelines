process FORMAT_GATK {
  tag "$sampleID"

  cpus 1
  memory 15.GB
  time '24:00:00'
  clusterOptions '-q batch'

  // change to bedtools container (make sure cat statment in container)
  container 'broadinstitute/gatk:4.2.4.1'
  file(params.ref_fai)

  input:
  tuple val(sampleID), file(txt)
  val(L)

  output:
  tuple val(sampleID), file("*_gatk_formatter.txt"), emit: txt

  script:
  log.info "----- GATK Formatter Running on: ${sampleID} -----"
  """
  chmod +x ${params.gatk_form}
  ${params.gatk_form} ${txt} ${sampleID}_gatk_temp2.txt ${sampleID}_gatk_formatter.txt ${L}
  """
}
process COVCALC_GATK {
  tag "$sampleID"

  cpus 1
  memory 15.GB
  time '24:00:00'
  clusterOptions '-q batch'

  container 'python_2.7.sif'

  // store in /stats
  publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID : 'gatk' }", pattern: "*.bed", mode:'copy'

  input:
  tuple val(sampleID), file(txt)
  val(filename)

  output:
  tuple val(sampleID), file("*.bed"), emit: bed

  script:
  log.info "----- GATK COVCALC Running on: ${sampleID} -----"

  """
  python ${params.cov_calc} ${txt} ${sampleID}_${filename}_avg_median_coverage.bed
  """
}
