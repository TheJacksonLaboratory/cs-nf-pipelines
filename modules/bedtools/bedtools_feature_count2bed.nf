process FEATURE_COUNT2BED {
  tag "$sampleID"

  cpus = 1

  publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID : 'bash-utils' }", pattern: "*_peaks_countMatrix.mm10.bed", mode: 'copy'
  container 'quay.io/biocontainers/bedtools:2.23.0--h5b5514e_6'

  input:
  tuple val(sampleID), file(peak_cnt_matrx)

  output:
  tuple val(sampleID), file("*_peaks_countMatrix.mm10.bed")

  shell:
  log.info "----- Feature Count to Bed on ${sampleID} -----"
  '''
  tail -n +3 !{peak_cnt_matrx} \
  | awk -F $'\\t' 'BEGIN {OFS = FS} { print $2, $3, $4, $7, $6 }' \
  > !{sampleID}_peaks_countMatrix.mm10.bed
  '''
}
