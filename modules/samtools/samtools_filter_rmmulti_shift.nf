process FILTER_RMMULTI_SHIFT {
  tag "$sampleID"

  cpus = 1

  publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID : 'samtools' }", pattern: "*.sorted.rmDup.rmChrM.rmMulti.filtered.ba*", mode: 'copy' 
  container 'library://taihpw/collection/samtools-atac:1.3.1'

  input:
  tuple val(sampleID), file(mtdna_bam_file)
  tuple val(sampleID), file(mtdna_bai_file)

  output:
  tuple val(sampleID), file("*.shift.tmp0.ba*")
  tuple val(sampleID), file("*.sorted.rmDup.rmChrM.rmMulti.filtered.ba*"), emit: srf_bam

  script:
  log.info "----- Filter Non-Unique and Include Only 'properly mapped reads' Alignments on ${sampleID} -----"
  """
  samtools view -@ $task.cpus -h -q 30 ${mtdna_bam_file} \
  > ${sampleID}.sorted.rmDup.rmChrM.rmMulti.bam

  samtools view -@ ${params.threads} -h -b -F 1804 -f 2 \
  ${sampleID}.sorted.rmDup.rmChrM.rmMulti.bam \
  > ${sampleID}.sorted.rmDup.rmChrM.rmMulti.filtered.bam

  samtools index \
  ${sampleID}.sorted.rmDup.rmChrM.rmMulti.filtered.bam

  samtools sort \
  -@ $task.cpus -O bam \
  -o ${sampleID}.shift.tmp0.bam \
  ${sampleID}.sorted.rmDup.rmChrM.rmMulti.filtered.bam

  samtools index \
  ${sampleID}.shift.tmp0.bam
  """

}
