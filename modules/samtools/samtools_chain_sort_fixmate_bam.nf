process CHAIN_SORT_FIXMATE_BAM {
  tag "$sampleID"

  cpus = 1

  publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID : 'samtools' }", pattern: "*.${params.processed_bam_extension}.ba*", mode: 'copy'
  container 'library://taihpw/collection/samtools-atac:1.3.1'

  input:
  tuple val(sampleID), file(bam_mm10)

  output:
  tuple val(sampleID), file("*.${params.processed_bam_extension}.ba*")

  when: params.chain != null

  script:
  log.info "----- Performing name sort the bam on ${sampleID} -----"
  """
  samtools sort \
  -n \
  -@ $task.cpus -O bam \
  -o ${sampleID}.tmp3.mm10.bam ${bam_mm10[0]}

  samtools fixmate \
  -O bam ${sampleID}.tmp3.mm10.bam ${sampleID}.tmp4.mm10.bam

  samtools sort \
  -@ $task.cpus -O bam \
  -o ${sampleID}.tmp5.mm10.bam ${sampleID}.tmp4.mm10.bam

  samtools index ${sampleID}.tmp5.mm10.bam

  samtools view ${sampleID}.tmp5.mm10.bam \
  -h 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 X Y \
  | grep -ve 'SN:MT*' > tmp.sam

  samtools view -b tmp.sam \
  > ${sampleID}.${params.processed_bam_extension}.bam

  samtools index \
  ${sampleID}.${params.processed_bam_extension}.bam
  """
}
