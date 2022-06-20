process RM_DUPES_READS {
  tag "$sampleID"

  cpus = 1

  container 'library://taihpw/collection/samtools-atac:1.3.1'

  input:
  tuple val(sampleID), file(marked_bam_file)
  tuple val(sampleID), file(marked_bai_file)

  output:
  tuple val(sampleID), file("*.sorted.rmDup.bam"), emit: rmDup_bam
  tuple val(sampleID), file("*.sorted.rmDup.bam.bai"), emit: rmDup_bai

  script:
  log.info "----- Samtools Removing PCR Duplicates on: ${sampleID} -----"
  """
  samtools view -h -b -F 1024 ${marked_bam_file} > ${sampleID}.sorted.rmDup.bam

  samtools index ${sampleID}.sorted.rmDup.bam
  """

}
