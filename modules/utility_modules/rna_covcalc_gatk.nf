process COVCALC_GATK {
  tag "$sampleID"

  cpus 1
  memory 15.GB
  time '24:00:00'

  container 'quay.io/jaxcompsci/python-bz2file:np_2.7.18'

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
  python ${projectDir}/bin/rnaseq/coveragecalculator.py ${txt} ${sampleID}_${filename}_avg_median_coverage.bed
  """
}
