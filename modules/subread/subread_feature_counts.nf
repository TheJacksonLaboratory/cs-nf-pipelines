process FEATURE_COUNTS {
  tag "$sampleID"

  cpus = 1

  publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID : 'subread' }", pattern: "*_peaks_countMatrix.txt", mode: 'copy'
  container 'quay.io/biocontainers/subread:1.6.4--h84994c4_1'

  input:
  tuple val(sampleID), file(processed_bams)
  tuple val(sampleID), file(peak_cvg_saf)

  output:
  tuple val(sampleID), file("*_peaks_countMatrix.txt")

  script:
  log.info "----- Feature Counts on ${sampleID} -----"
  """
  featureCounts \
  -a ${peak_cvg_saf} \
  -F SAF -p \
  -T $task.cpus \
  -o ${sampleID}_peaks_countMatrix.txt \
  ${processed_bams[0]}
  """
}
