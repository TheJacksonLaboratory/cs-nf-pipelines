process FRIP_READS_IN_PEAKS {
  tag "$sampleID"

  cpus 2
  memory 4.GB 
  time '04:00:00'

  container 'quay.io/biocontainers/bedtools:2.23.0--h5b5514e_6'

  input:
  tuple val(sampleID), file(processed_bams)
  tuple val(sampleID), file(narrow_peaks)

  output:
  tuple val(sampleID), file("reads_in_peaks.tmp.ba*")

  script:
  log.info "----- Fraction of reads in peaks (FRiP) on ${sampleID} -----"
  """
  bedtools sort \
  -i ${narrow_peaks} \
  | bedtools merge -i stdin \
  | bedtools intersect -u -nonamecheck \
  -a ${processed_bams[0]} \
  -b stdin \
  -ubam \
  > reads_in_peaks.tmp.bam
  """
}
