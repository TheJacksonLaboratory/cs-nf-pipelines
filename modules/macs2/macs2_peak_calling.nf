process PEAK_CALLING {
  tag "$sampleID"

  cpus = 1

  publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID : 'macs2' }", pattern: "*_peaks.narrowPeak", mode: 'copy'
  publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID : 'macs2' }", pattern: "*_summits.bed", mode: 'copy'
  publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID : 'macs2' }", pattern: "*.log", mode: 'copy'
  container 'quay.io/biocontainers/macs2:2.2.7.1--py39hbf8eff0_4'  


  input:
  tuple val(sampleID), file(processed_bams)

  output:
  tuple val(sampleID), file("*_peaks.narrowPeak"), emit: np
  tuple val(sampleID), file("*_summits.bed")
  tuple val(sampleID), file("*_macs2.log")


  script:
  log.info "----- Performing Peak Calling on on ${sampleID} -----"
  """
  macs2 callpeak \
  -f BAMPE \
  --nomodel \
  -g ${params.g} \
  --keep-dup all \
  --cutoff-analysis \
  --tempdir ${params.tmpdir} \
  -n ${sampleID} \
  -t ${processed_bams[0]} \
  --outdir . \
  > ${sampleID}_macs2.log 2>&1
  """
}
