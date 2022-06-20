process FILTER_RMMULTI_SORT {
  tag "$sampleID"

  cpus = 1

  container 'library://taihpw/collection/samtools-atac:1.3.1'

  input:
  tuple val(sampleID), file(shift_bams)

  output:
  tuple val(sampleID), file("*.sorted.rmDup.rmChrM.rmMulti.filtered.shifted.ba*")

  script:
  log.info "----- Re-sorting shifted bam on ${sampleID} -----"
  """
  samtools sort \
  -@ $task.cpus -O bam \
  -o ${sampleID}.sorted.rmDup.rmChrM.rmMulti.filtered.shifted.bam \
  ${shift_bams[0]}
  """

}
