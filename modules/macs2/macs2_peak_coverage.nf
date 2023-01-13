process PEAK_COVERAGE {
  tag "$sampleID"

  cpus = 1
  memory 1.GB
  time '01:00:00'

  container 'quay.io/biocontainers/macs2:2.2.7.1--py39hbf8eff0_4'  

  input:
  tuple val(sampleID), file(narrow_peaks)

  output:
  tuple val(sampleID), file("*_peaks.narrowPeak.saf")

  shell:
  log.info "----- Get coverage in each peak on ${sampleID} -----"
  '''
  awk 'OFS="\\t" {print $1"."$2"."$3, $1, $2, $3, "."}' !{narrow_peaks} \
  > !{sampleID}_peaks.narrowPeak.saf
  '''
}
