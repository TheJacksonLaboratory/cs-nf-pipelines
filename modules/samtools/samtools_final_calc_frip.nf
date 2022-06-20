process FINAL_CALC_FRIP {
  tag "$sampleID"

  cpus = 1

  publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID : 'samtools' }", pattern: "*_Fraction_reads_in_peak.txt", mode: 'copy'
  container 'library://taihpw/collection/samtools-atac:1.3.1'

  input:
  tuple val(sampleID), file(processed_bams)
  tuple val(sampleID), file(reads_peaks_bams)

  output:
  tuple val(sampleID), file("*_Fraction_reads_in_peak.txt")

  shell:
  log.info "----- Final Calculate (FRiP) on ${sampleID} -----"
  '''
  total_reads=$(samtools view -c !{processed_bams[0]})
  reads_in_peaks=$(samtools view -c !{reads_peaks_bams[0]})
  FRiP=$(awk "BEGIN {print "${reads_in_peaks}"/"${total_reads}"}")
  echo -e ${FRiP}"\\t"${total_reads} \
  > !{sampleID}_Fraction_reads_in_peak.txt
  '''
}
